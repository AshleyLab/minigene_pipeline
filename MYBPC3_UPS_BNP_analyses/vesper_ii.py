from __future__ import annotations
import re
from pathlib import Path
import pandas as pd
import warnings
import os
from itertools import combinations
import yaml

# silence future warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)

### Determine the full path of which pipeline is in ###
vesper_dir = os.path.dirname(os.path.realpath(__file__))

### Load config.yaml ###
config_file = os.path.join(vesper_dir, "vesper_config.yaml")
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

### Assign variables from vesper_config.yaml ###
ROLE_MAP = config["ROLE_MAP"]
PRIMER_RANGES = config["PRIMER_RANGES"]

maf_cutoff = config["maf_cutoff"]
cx_count_files_dir = config["cx_count_files_dir"]
ipsc_list_path = config["ipsc_list_path"]
mm2_sample_rep_raw_variants_dir = config["vesper_i_raw_variants_dir"]
vesperii_output_dir = config["vesperii_output_dir"]
os.makedirs(vesperii_output_dir, exist_ok=True)

combined_reps_count_files_dir = os.path.join(vesperii_output_dir, "combined_reps_hisat2_count_files")
combined_t1t2_count_files_dir = os.path.join(vesperii_output_dir, "combined_t1t2_hisat2_count_files")

mm2_sample_variants_dir = os.path.join(vesperii_output_dir, "sample_mm2_mnv_counts")
combined_t1t2_mm2_vars_dir = os.path.join(vesperii_output_dir, "combined_t1t2_mm2_mnv_counts")
mm2_hisat_paired_counts_dir = os.path.join(vesperii_output_dir, "mm2_hisat2_paired_counts")
mm2_snv_counts_dir = os.path.join(vesperii_output_dir, "mm2_snv_counts")

def build_blacklist(ranges):
    """Convert [(start, end), ...] to a set of integers."""
    out = set()
    for start, end in ranges:
        out.update(range(start, end + 1))
    return out

# precompute blacklists (fast lookup)
PRIMER_BLACKLISTS = {
    tile: build_blacklist(ranges)
    for tile, ranges in PRIMER_RANGES.items()
}


###############################################################################################################################
                                ### OBTAIN VARIANT STATISTICS FROM CX PIPELINE COUNT FILES ###
###############################################################################################################################

def combine_reps_by_protein_tile(cx_count_files_dir, combined_reps_count_files_dir):
    """
    Combine hisat2 count files across replicates, per protein and tile by summing counts and averaging enrichment, and removing primer regions.
    """
    cx_count_files_dir = Path(cx_count_files_dir)
    combined_reps_count_files_dir = Path(combined_reps_count_files_dir)
    combined_reps_count_files_dir.mkdir(parents=True, exist_ok=True)

    # match: rep1_iPSC_t1.csv, rep12_BNP_t2.csv, etc.
    fn_re = re.compile(r"^(rep\d+)_([^_]+)_(t\d+)\.csv$", re.IGNORECASE)

    keep_cols = ["cDNA", "reference", "observed_nucleotide", "exp_enrichment", "num_obs_exp", "num_obs_con"]

    def load_one_csv(fp: Path) -> pd.DataFrame:
        df = pd.read_csv(fp)

        missing = [c for c in keep_cols if c not in df.columns]
        if missing:
            raise ValueError(f"{fp} is missing columns: {missing}")

        df = df[keep_cols].copy()

        # Variant = cDNA + reference + ">" + observed_nucleotide
        df["Variant"] = (
            df["cDNA"].astype(str).str.strip()
            + df["reference"].astype(str).str.strip()
            + ">"
            + df["observed_nucleotide"].astype(str).str.strip()
        )

        df = df.drop(columns=["cDNA", "reference", "observed_nucleotide"])

        # Make sure numeric cols are numeric (coerce bad values to NaN then fill)
        for col in ["exp_enrichment", "num_obs_exp", "num_obs_con"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        return df

    # Discover protein directories automatically
    protein_dirs = [p for p in cx_count_files_dir.iterdir() if p.is_dir()]

    if not protein_dirs:
        raise FileNotFoundError(f"No protein directories found under {cx_count_files_dir}")

    for prot_dir in sorted(protein_dirs):
        protein = prot_dir.name

        # bucket files by tile
        tile_buckets: dict[str, list[Path]] = {}

        for fp in prot_dir.glob("*.csv"):
            m = fn_re.match(fp.name)
            if not m:
                continue

            _rep, protein_from_file, tile = m.groups()

            # safety check: directory name must match filename protein
            if protein_from_file != protein:
                raise ValueError(
                    f"Protein mismatch: directory={protein}, file={fp.name}"
                )

            tile_buckets.setdefault(tile, []).append(fp)

        for tile, fps in sorted(tile_buckets.items()):
            if len(fps) == 0:
                continue

            per_rep_dfs = [load_one_csv(fp) for fp in fps]
            combined = pd.concat(per_rep_dfs, ignore_index=True)

            out = (
                combined.groupby("Variant", as_index=False)
                        .agg(
                            num_obs_exp=("num_obs_exp", "sum"),
                            num_obs_con=("num_obs_con", "sum"),
                            sum_exp_enrichment=("exp_enrichment", "sum"),
                        )
            )

            # get number of per_rep_dfs files for averaging
            num_reps = len(per_rep_dfs)
            out['mean_exp_enrichment'] = out['sum_exp_enrichment'] / num_reps

            # position for sorting + filtering
            out["Variant_position"] = (out["Variant"].str.extract(r"(-?\d+)", expand=False).astype(float))

            # remove primer regions
            # --- tile-specific position filters (from global PRIMER_BLACKLISTS) ---
            blacklist = PRIMER_BLACKLISTS.get(tile)
            if blacklist is not None:
                out = out[~out["Variant_position"].isin(blacklist)]

            # sort, then drop helper column
            out = (out.sort_values(by="Variant_position", na_position="last").drop(columns=["Variant_position"]))

            out_path = combined_reps_count_files_dir / f"{protein}_{tile}.csv"
            out.to_csv(out_path, index=False)

def merge_tiles_per_protein(combined_reps_count_files_dir, combined_t1t2_count_files_dir):
    """
    Combines T1 and T2 count files per protein into single T1T2 files, then perform MAF filtering.
    """
    in_dir = Path(combined_reps_count_files_dir)
    out_dir = Path(combined_t1t2_count_files_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # matches: iPSC_t1.csv, BNP_t2.csv, MYBPC3_t1.csv, etc.
    fn_re = re.compile(r"^(?P<protein>.+)_(?P<tile>t\d+)\.csv$", re.IGNORECASE)

    # Load all files, bucket by protein
    buckets: dict[str, dict[str, Path]] = {}
    for fp in sorted(in_dir.glob("*.csv")):
        m = fn_re.match(fp.name)
        if not m:
            continue
        protein = m.group("protein")
        tile = m.group("tile")
        buckets.setdefault(protein, {})[tile] = fp

    if not buckets:
        raise FileNotFoundError(f"No files like {{protein}}_t{{n}}.csv found in {in_dir}")

    def load_and_suffix(fp: Path, tile: str) -> pd.DataFrame:
        df = pd.read_csv(fp)
        if "Variant" not in df.columns:
            raise ValueError(f"{fp} missing required column 'Variant'")

        rename = {c: f"{c}_{tile}" for c in df.columns if c != "Variant"}
        df = df.rename(columns=rename)

        # Ensure Variant is clean
        df["Variant"] = df["Variant"].astype(str).str.strip()
        return df

    # Columns you expect / want special fill rules for
    count_cols_base = ["num_obs_exp", "num_obs_con"]
    enrich_cols_base = ["sum_exp_enrichment", "mean_exp_enrichment"]

    for protein, tile_map in sorted(buckets.items()):
        # Load whichever tiles exist
        df_t1 = load_and_suffix(tile_map["t1"], "t1") if "t1" in tile_map else None
        df_t2 = load_and_suffix(tile_map["t2"], "t2") if "t2" in tile_map else None

        if df_t1 is None and df_t2 is None:
            continue

        # Outer merge (keeps Variants present in either tile)
        if df_t1 is None:
            merged = df_t2.copy()
        elif df_t2 is None:
            merged = df_t1.copy()
        else:
            merged = df_t1.merge(df_t2, on="Variant", how="outer")

        # --- Apply your fill rules ---
        # For counts: missing -> 0
        for base in count_cols_base:
            for tile in ("t1", "t2"):
                col = f"{base}_{tile}"
                if col in merged.columns:
                    merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0).astype(int)
                else:
                    merged[col] = 0

        # For enrichment columns: missing -> '-'
        for base in enrich_cols_base:
            for tile in ("t1", "t2"):
                col = f"{base}_{tile}"
                if col in merged.columns:
                    # keep numeric if present, but represent missing as '-'
                    merged[col] = pd.to_numeric(merged[col], errors="coerce")
                    merged[col] = merged[col].where(merged[col].notna(), "-")
                else:
                    merged[col] = "-"

        # --- Build t1t2 columns ---
        merged["num_obs_exp_t1t2"] = merged["num_obs_exp_t1"] + merged["num_obs_exp_t2"]
        merged["num_obs_con_t1t2"] = merged["num_obs_con_t1"] + merged["num_obs_con_t2"]

        # sum_exp_enrichment_t1t2: treat '-' as 0
        sum1 = pd.to_numeric(merged["sum_exp_enrichment_t1"].replace("-", 0), errors="coerce").fillna(0)
        sum2 = pd.to_numeric(merged["sum_exp_enrichment_t2"].replace("-", 0), errors="coerce").fillna(0)
        merged["sum_exp_enrichment_t1t2"] = sum1 + sum2

        # mean_exp_enrichment_t1t2 rules:
        # - adopt the non '-' value if only one exists
        # - if both exist, average them
        m1 = pd.to_numeric(merged["mean_exp_enrichment_t1"].replace("-", pd.NA), errors="coerce")
        m2 = pd.to_numeric(merged["mean_exp_enrichment_t2"].replace("-", pd.NA), errors="coerce")

        both = m1.notna() & m2.notna()
        only1 = m1.notna() & m2.isna()
        only2 = m1.isna() & m2.notna()

        merged["mean_exp_enrichment_t1t2"] = "-"
        merged.loc[only1, "mean_exp_enrichment_t1t2"] = m1[only1]
        merged.loc[only2, "mean_exp_enrichment_t1t2"] = m2[only2]
        merged.loc[both,  "mean_exp_enrichment_t1t2"] = (m1[both] + m2[both]) / 2

        #  sort by signed position extracted from Variant
        merged["Variant_position"] = merged["Variant"].str.extract(r"(-?\d+)", expand=False).astype(int)
        merged = merged.sort_values("Variant_position", na_position="last")

        # create and append depth_obs_exp and depth_obs_con columns (depth is specific to variant position)
        merged["depth_obs_con_t1t2"] = (merged.groupby("Variant_position")["num_obs_con_t1t2"].transform("sum"))
        merged["depth_obs_exp_t1t2"] = (merged.groupby("Variant_position")["num_obs_exp_t1t2"].transform("sum"))

        # calculate MAF
        # total counts per variant
        merged["total_counts"] = merged["num_obs_exp_t1t2"] + merged["num_obs_con_t1t2"]
        # denominator: sum of total_counts per position
        pos_totals = merged.groupby("Variant_position")["total_counts"].transform("sum")
        # MAF
        merged["MAF"] = merged["total_counts"] / pos_totals
        merged = merged.drop(columns=["total_counts", "Variant_position"])

        # --- Rename exp/con columns based on protein role map (if defined) ---
        if protein in ROLE_MAP:
            exp_label = ROLE_MAP[protein]["exp"]
            con_label = ROLE_MAP[protein]["control"]

            merged = merged.rename(columns={
                "num_obs_exp_t1t2": f"count_{protein}_{exp_label}_T1T2",
                "num_obs_con_t1t2": f"count_{protein}_{con_label}_T1T2",
                "depth_obs_exp_t1t2": f"depth_{protein}_{exp_label}_T1T2",
                "depth_obs_con_t1t2": f"depth_{protein}_{con_label}_T1T2"
            })

        out_path = out_dir / f"{protein}_T1T2.csv"
        merged.to_csv(out_path, index=False)

        # remove rows where MAF < maf_cutoff
        merged = merged[merged["MAF"] >= maf_cutoff].reset_index(drop=True)
        out_MAF_path = out_dir / f"{protein}_T1T2_MAF_filtered.csv"
        merged.to_csv(out_MAF_path, index=False)


###############################################################################################################################
                                                ###  PROCESS VESPER_I FILES ###
###############################################################################################################################

def cdna_variants_per_sample(mm2_sample_rep_raw_variants_dir, col="cdna_variants"):
    """
    Combine minimap2 MNV variant calls across replicates, per sample and tile, counting occurrences of each unique variant, and removing primer regions.
    """
    mm2_sample_rep_raw_variants_dir = Path(mm2_sample_rep_raw_variants_dir)
    pat = re.compile(r"(?P<rep>\d+)_(?P<sample>.+?)_(?P<tile>T\d)\.csv$", re.IGNORECASE)
    buckets = {}

    for p in sorted(mm2_sample_rep_raw_variants_dir.glob("*.csv")):
        m = pat.match(p.name)
        if not m:
            continue

        sample = m["sample"]
        tile = m["tile"].lower()   # normalize T1 / t1
        key = f"{sample}_{tile}"
        df = pd.read_csv(p)

        if col not in df.columns:
            raise ValueError(f"{p.name} is missing required column '{col}'")

        buckets.setdefault(key, []).append(df[[col]].copy())

    outputs = []
    for key, parts in buckets.items():
        if not parts:
            continue
        cat = pd.concat(parts, ignore_index=True)

        # drop NAs; remove this line if you want to count NaNs as a category
        cat = cat.dropna(subset=[col])

        counts = (
            cat.groupby(col, as_index=False)
               .size()
               .rename(columns={"size": "count"})
               .sort_values("count", ascending=False)
        )


        # --- position + primer filtering ---
        counts["Variant_position"] = (counts[col].astype(str).str.extract(r"(-?\d+)", expand=False).astype(float))

        tile = key.split("_")[-1]  # t1 / t2
        blacklist = PRIMER_BLACKLISTS.get(tile)
        if blacklist is not None:
            counts = counts[~counts["Variant_position"].isin(blacklist)]

        counts = counts.drop(columns=["Variant_position"])

        mm2_sample_variants_out_dir = Path(mm2_sample_variants_dir)
        mm2_sample_variants_out_dir.mkdir(parents=True, exist_ok=True)

        sample, tile = key.rsplit("_", 1)
        tile_out = tile.upper()   # t1 -> T1
        out_path = mm2_sample_variants_out_dir / f"{sample}_{tile_out}.csv"
        counts.to_csv(out_path, index=False)
        outputs.append(str(out_path))

    if not outputs:
        raise FileNotFoundError("No files matched pattern {rep}_{sample}_{tile}.csv")

    return outputs


def merge_mm2_vars_T1T2_counts(mm2_sample_variants_dir, combined_t1t2_mm2_vars_dir):
    """
    Merge {sample}_{condition}_T1.csv and {sample}_{condition}_T2.csv by cdna_variants and summing counts, producing {sample}_{condition}_T1T2.csv files.
    Also produces filtered versions keeping only variants present in the provided iPSC list.
    """
    mm2_sample_variants_dir = Path(mm2_sample_variants_dir)
    combined_t1t2_mm2_vars_dir = Path(combined_t1t2_mm2_vars_dir)
    combined_t1t2_mm2_vars_dir.mkdir(parents=True, exist_ok=True)

    pat = re.compile(r"^(?P<prefix>.+)_(?P<tile>T[12])\.csv$", re.IGNORECASE)

    # bucket files by prefix (sample_condition)
    buckets: dict[str, dict[str, Path]] = {}

    for fp in sorted(mm2_sample_variants_dir.glob("*.csv")):
        m = pat.match(fp.name)
        if not m:
            continue
        prefix = m.group("prefix")              # "{sample}_{condition}"
        tile = m.group("tile").upper()          # "T1" or "T2"
        buckets.setdefault(prefix, {})[tile] = fp

    if not buckets:
        raise FileNotFoundError(f"No *_T1.csv / *_T2.csv files found in {mm2_sample_variants_dir}")

    def load_one(fp: Path, prefix: str, tile: str) -> pd.DataFrame:
        df = pd.read_csv(fp)
        required = {"cdna_variants", "count"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{fp.name} missing columns: {sorted(missing)}")

        df = df[["cdna_variants", "count"]].copy()
        df["cdna_variants"] = df["cdna_variants"].astype(str).str.strip()

        # Rename count -> count_{sample}_{condition}_T1/T2  (prefix already "{sample}_{condition}")
        df = df.rename(columns={"count": f"count_{prefix}_{tile}"})

        # If the same cdna_variants appears multiple times in a file, sum it (defensive)
        df = df.groupby("cdna_variants", as_index=False).sum(numeric_only=True)
        return df

    for prefix, tile_map in sorted(buckets.items()):
        df_t1 = load_one(tile_map["T1"], prefix, "T1") if "T1" in tile_map else None
        df_t2 = load_one(tile_map["T2"], prefix, "T2") if "T2" in tile_map else None

        # Outer merge so we keep variants present in either tile
        if df_t1 is None and df_t2 is None:
            continue
        elif df_t1 is None:
            merged = df_t2.copy()
            # ensure missing T1 column exists
            merged[f"count_{prefix}_T1"] = 0
        elif df_t2 is None:
            merged = df_t1.copy()
            merged[f"count_{prefix}_T2"] = 0
        else:
            merged = df_t1.merge(df_t2, on="cdna_variants", how="outer")

        # Fill missing counts with 0 and force ints
        for tile in ("T1", "T2"):
            col = f"count_{prefix}_{tile}"
            if col not in merged.columns:
                merged[col] = 0
            merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0).astype(int)

        # Keep columns in the exact order you want
        merged = merged[["cdna_variants", f"count_{prefix}_T1", f"count_{prefix}_T2"]]
        merged[f"count_{prefix}_T1T2"] = merged[f"count_{prefix}_T1"] + merged[f"count_{prefix}_T2"]

        out_path = combined_t1t2_mm2_vars_dir / f"{prefix}_T1T2.csv"
        merged.to_csv(out_path, index=False)

        ### filter by iPSC-allowed variants ###
        ipsc_list = pd.read_csv(ipsc_list_path)
        allowed = set(ipsc_list["Variant"].astype(str).str.strip())

        def keep_allowed(s: str) -> str:
            if pd.isna(s):
                return ""
            # split on commas with optional whitespace (handles ", " and "," etc.)
            parts = [p.strip() for p in re.split(r"\s*,\s*", str(s).strip()) if p.strip()]
            kept = [p for p in parts if p in allowed]
            return ", ".join(kept)

        mm2_t1t2_ipsc_filtered = merged.copy()
        mm2_t1t2_ipsc_filtered["cdna_variants"] = mm2_t1t2_ipsc_filtered["cdna_variants"].apply(keep_allowed)

        # drop rows where nothing remains
        mm2_t1t2_ipsc_filtered = mm2_t1t2_ipsc_filtered[mm2_t1t2_ipsc_filtered["cdna_variants"].str.len() > 0].reset_index(drop=True)

        mm2_t1t2_ipsc_filtered = (mm2_t1t2_ipsc_filtered.groupby("cdna_variants", as_index=False)[[f"count_{prefix}_T1", f"count_{prefix}_T2", f"count_{prefix}_T1T2"]].sum())

        # save filtered output (recommended)
        out_path = combined_t1t2_mm2_vars_dir / f"{prefix}_T1T2_iPSC_list_vars.csv"
        mm2_t1t2_ipsc_filtered.to_csv(out_path, index=False)

def merge_T1T2_mm2_mnv_by_conditions(combined_t1t2_mm2_vars_dir):
    """
    Merge {protein}_{condition}_T1T2_iPSC_list_vars.csv files across conditions to get {protein}_T1T2_iPSC_list_vars.csv files
    """
    in_dir = Path(combined_t1t2_mm2_vars_dir)

    # filename example: BNP_low_T1T2_iPSC_list_vars.csv
    pat = re.compile(r"^(?P<protein>.+?)_(?P<condition>.+?)_T1T2_iPSC_list_vars\.csv$", re.IGNORECASE)

    # bucket filepaths by protein
    buckets: dict[str, list[tuple[str, Path]]] = {}
    for fp in sorted(in_dir.glob("*_T1T2_iPSC_list_vars.csv")):
        m = pat.match(fp.name)
        if not m:
            continue
        protein = m.group("protein")
        condition = m.group("condition")
        buckets.setdefault(protein, []).append((condition, fp))

    if not buckets:
        raise FileNotFoundError(f"No *_T1T2_iPSC_list_vars.csv files found in {in_dir}")

    outputs: list[str] = []

    for protein, cond_files in sorted(buckets.items()):
        merged: pd.DataFrame | None = None

        for condition, fp in sorted(cond_files, key=lambda x: x[0]):
            df = pd.read_csv(fp)

            # expected columns
            c1 = f"count_{protein}_{condition}_T1"
            c2 = f"count_{protein}_{condition}_T2"
            c12 = f"count_{protein}_{condition}_T1T2"

            required = {"cdna_variants", c1, c2, c12}
            missing = required - set(df.columns)
            if missing:
                raise ValueError(f"{fp.name} missing columns: {sorted(missing)}")

            df = df[["cdna_variants", c1, c2, c12]].copy()
            df["cdna_variants"] = df["cdna_variants"].astype(str).str.strip()

            # ensure numeric + fill NA with 0
            for c in (c1, c2, c12):
                df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

            if merged is None:
                merged = df
            else:
                merged = merged.merge(df, on="cdna_variants", how="outer")

        assert merged is not None

        # fill any missing counts created by outer merges
        count_cols = [c for c in merged.columns if c.startswith("count_")]
        for c in count_cols:
            merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0).astype(int)

        # optional: sort by first numeric position in cdna_variants
        merged["_pos"] = merged["cdna_variants"].str.extract(r"(-?\d+)", expand=False).astype(float)
        merged = merged.sort_values("_pos", na_position="last").drop(columns=["_pos"]).reset_index(drop=True)

        out_path = in_dir / f"{protein}_T1T2_iPSC_list_vars.csv"
        merged.to_csv(out_path, index=False)
        outputs.append(str(out_path))

    return outputs

def filter_mm2_t1t2_by_maf(combined_t1t2_mm2_vars_dir, combined_t1t2_count_files_dir):
    """
    Filter by MAF i.e. keep only variants in {protein}_T1T2_iPSC_list_vars.csv files found in the corresponding {protein}_T1T2_MAF_filtered.csv file
    """
    mm2_dir = Path(combined_t1t2_mm2_vars_dir)
    maf_dir = Path(combined_t1t2_count_files_dir)

    # Match ONLY: {protein}_T1T2_iPSC_list_vars.csv
    # This will NOT match: {protein}_{condition}_T1T2_iPSC_list_vars.csv
    fn_re = re.compile(r"^(?P<protein>[^_]+)_T1T2_iPSC_list_vars\.csv$", re.IGNORECASE)

    outputs: list[str] = []

    for fp in sorted(mm2_dir.glob("*.csv")):
        m = fn_re.match(fp.name)
        if not m:
            continue

        protein = m.group("protein")

        # Load mm2 t1t2 file
        df = pd.read_csv(fp)

        if "cdna_variants" not in df.columns:
            raise ValueError(f"{fp.name} missing required column 'cdna_variants'")

        # Keep only columns ending in T1T2 (plus cdna_variants)
        t1t2_cols = [c for c in df.columns if c.endswith("T1T2")]
        keep_cols = ["cdna_variants"] + t1t2_cols
        df = df[keep_cols].copy()

        # Load matching MAF file
        maf_fp = maf_dir / f"{protein}_T1T2_MAF_filtered.csv"
        if not maf_fp.exists():
            raise FileNotFoundError(f"Missing MAF file for {protein}: {maf_fp}")

        maf_df = pd.read_csv(maf_fp)
        if "Variant" not in maf_df.columns:
            raise ValueError(f"{maf_fp.name} missing required column 'Variant'")

        allowed = set(maf_df["Variant"].astype(str).str.strip())

        # Filter each row's cdna_variants list
        splitter = re.compile(r"\s*,\s*")

        def keep_allowed_list(s: str) -> str:
            if pd.isna(s):
                return ""
            parts = [p.strip() for p in splitter.split(str(s).strip()) if p.strip()]
            kept = [p for p in parts if p in allowed]
            return ", ".join(kept)

        df["cdna_variants"] = df["cdna_variants"].apply(keep_allowed_list)

        # Drop rows that become empty after filtering
        df = df[df["cdna_variants"].str.len() > 0].reset_index(drop=True)

        # Ensure numeric counts; sum duplicates after filtering (recommended)
        for c in t1t2_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

        df = df.groupby("cdna_variants", as_index=False)[t1t2_cols].sum()

        # Optional: sort by first position in string for nicer output
        df["_pos"] = df["cdna_variants"].str.extract(r"(-?\d+)", expand=False).astype(float)
        df = df.sort_values("_pos", na_position="last").drop(columns=["_pos"])

        out_fp = mm2_dir / f"{protein}_T1T2_iPSC_list_MAF_filtered.csv"
        df.to_csv(out_fp, index=False)
        outputs.append(str(out_fp))

    if not outputs:
        raise FileNotFoundError(
            f"No files matched {{protein}}_T1T2_iPSC_list_vars.csv in {mm2_dir}"
        )

    return outputs

def get_variant_pairs(combined_t1t2_mm2_vars_dir, mm2_hisat_paired_counts_dir):
    """
    Processes mm2 MNVs after iPSC list and MAF filtering to get paired variant counts.
    """
    combined_t1t2_mm2_vars_dir = Path(combined_t1t2_mm2_vars_dir)
    mm2_hisat_paired_counts_dir = Path(mm2_hisat_paired_counts_dir)
    mm2_hisat_paired_counts_dir.mkdir(parents=True, exist_ok=True)

    ### list all CSV files in the input directory that has filenames ending with MAF_filtered.csv ###
    csv_files = sorted([p for p in combined_t1t2_mm2_vars_dir.glob("*.csv") if p.name.endswith("MAF_filtered.csv")])
    if not csv_files:
        raise FileNotFoundError(f"No *MAF_filtered.csv files found in {combined_t1t2_mm2_vars_dir}")
    
    ### protein parser: grabs everything before "_T1T2" ###
    prot_re = re.compile(r"^(?P<protein>.+?)_T1T2_.*MAF_filtered\.csv$", re.IGNORECASE)

    outputs = []

    for csv_file in csv_files:
        m = prot_re.match(csv_file.name)
        protein = m.group("protein") if m else csv_file.stem  # fallback: filename without .csv

        ### read csv into df ###
        df = pd.read_csv(csv_file)                              
        if "cdna_variants" not in df.columns:
            raise ValueError(f"{csv_file.name} missing required column 'cdna_variants'")

        ### pick all count columns (or restrict to T1T2 counts only) ###
        count_cols = [c for c in df.columns if c.startswith("count_")]
        if not count_cols:
            raise ValueError(f"{csv_file.name} has no columns starting with 'count_'")

        ### initialise pair_counts dictionary to store variant pairs and their counts ###                    
        pair_counts = {}

        for _, row in df.iterrows():
            variants = [
                v.strip()
                for v in re.split(r"\s*,\s*", str(row["cdna_variants"]))
                if v.strip()
            ]

            ### determine pairs for EACH ROW ONLY ###
            if len(variants) == 1:
                pairs = [(variants[0], "-")]
            else:
                pairs = [tuple(sorted(p)) for p in combinations(variants, 2)]

            ### accumulate counts ###
            for v1, v2 in pairs:
                key = (v1, v2)
                if key not in pair_counts:
                    pair_counts[key] = {c: 0 for c in count_cols}

                for c in count_cols:
                    pair_counts[key][c] += int(pd.to_numeric(row[c], errors="coerce") or 0)

        ### convert to dataframe ###
        out = pd.DataFrame([{"Variant1": k[0], "Variant2": k[1], **v} for k, v in pair_counts.items()])

        ### sort by numeric position for readability ###
        out["_p1"] = out["Variant1"].str.extract(r"(-?\d+)", expand=False).astype(float)
        out["_p2"] = out["Variant2"].str.extract(r"(-?\d+)", expand=False).astype(float)
        out = out.sort_values(["_p1", "_p2", "Variant1", "Variant2"]).drop(columns=["_p1", "_p2"])

        out_path = mm2_hisat_paired_counts_dir / f"{protein}_paired_counts.csv"
        out.to_csv(out_path, index=False)
        outputs.append(str(out_path))

    return outputs

def sum_individual_snvs(df):
    # identify count columns (keeps this generic)
    count_cols = [c for c in df.columns if c.startswith("count_")]

    # split comma-separated variants
    df = df.copy()
    df["snv_variant"] = df["cdna_variants"].str.split(r"\s*,\s*")

    # explode into one variant per row
    flat = df.explode("snv_variant")

    # clean
    flat["snv_variant"] = flat["snv_variant"].str.strip()
    flat = flat[flat["snv_variant"] != ""]

    # ensure numeric counts
    for c in count_cols:
        flat[c] = pd.to_numeric(flat[c], errors="coerce").fillna(0)

    # sum per variant
    out = (
        flat.groupby("snv_variant", as_index=False)[count_cols]
            .sum()
    )

    # optional: sort by numeric position for readability
    out["_pos"] = out["snv_variant"].str.extract(r"(-?\d+)", expand=False).astype(int)
    out = out.sort_values("_pos").drop(columns="_pos").reset_index(drop=True)

    return out


def count_snvs_from_mm2_mnvs(combined_t1t2_mm2_vars_dir, mm2_snv_counts_dir):
    """
    calculate count of each variant regardless of whether it's travelling alone or with others from mm2 MNV calls ({protein}_T1T2_iPSC_list_MAF_filtered.csv files) after T1T2 combination, iPSC list variant retention and MAF filtering.
    """
    combined_t1t2_mm2_vars_dir = Path(combined_t1t2_mm2_vars_dir)
    mm2_snv_counts_dir = Path(mm2_snv_counts_dir)
    mm2_snv_counts_dir.mkdir(parents=True, exist_ok=True)

    pattern = f"*_T1T2_iPSC_list_MAF_filtered.csv"
    files = sorted(combined_t1t2_mm2_vars_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {combined_t1t2_mm2_vars_dir}")

    written = []
    for p in files:
        df = pd.read_csv(p)
        count_cols = [c for c in df.columns if c.startswith("count_")]
        if "cdna_variants" not in df.columns or not count_cols:
            raise ValueError(f"{p.name} missing required columns: {'cdna_variants'} and/or {count_cols!r}")

        ### identify each individual variant from comma-separated cdna_variants ###
        df = df.copy()
        df["snv_variant"] = df["cdna_variants"].str.split(r"\s*,\s*")

        ### explode into one variant per row and clean ###
        flat = df.explode("snv_variant")
        flat["snv_variant"] = flat["snv_variant"].str.strip()
        flat = flat[flat["snv_variant"] != ""]

        # Ensure numeric counts
        for count_col in count_cols:
            flat[count_col] = pd.to_numeric(flat[count_col], errors="coerce").fillna(0)

        ### sum per variant ###
        out = flat.groupby("snv_variant", as_index=False)[count_cols].sum()

        ### sort by numeric position for readability ###
        out["_pos"] = out["snv_variant"].str.extract(r"(-?\d+)", expand=False).astype(int)
        out = out.sort_values("_pos").drop(columns="_pos").reset_index(drop=True)

        protein = p.stem.split("_T1T2")[0]
        out_path = mm2_snv_counts_dir / f"{protein}_snv_counts.csv"
        out.to_csv(out_path, index=False)
        written.append(str(out_path))

    return written


def annotate_paired_counts_perc(mm2_snv_counts_dir, mm2_hisat_paired_counts_dir):
    """
    calculate percentage of times Variant1 and Variant2 appear in paired counts relative to their total mm2 snv counts, producing {protein}_paired_counts_pct.csv files.
    """
    mm2_snv_counts_dir = Path(mm2_snv_counts_dir)
    mm2_hisat_paired_counts_dir = Path(mm2_hisat_paired_counts_dir)

    # match protein from filenames like "BNP_snvs_counts.csv" and "BNP_paired_counts.csv"
    snv_pat = re.compile(r"^(?P<protein>.+)_snv_counts\.csv$", re.IGNORECASE)
    pair_pat = re.compile(r"^(?P<protein>.+)_paired_counts\.csv$", re.IGNORECASE)

    # index SNV count files by protein
    snv_files: dict[str, Path] = {}
    for p in mm2_snv_counts_dir.glob("*.csv"):
        m = snv_pat.match(p.name)
        if m:
            snv_files[m.group("protein")] = p

    written: list[str] = []

    for pair_fp in sorted(mm2_hisat_paired_counts_dir.glob("*_paired_counts.csv")):
        m = pair_pat.match(pair_fp.name)
        if not m:
            continue

        protein = m.group("protein")
        snv_fp = snv_files.get(protein)
        if snv_fp is None:
            raise FileNotFoundError(
                f"No matching SNV totals file for protein={protein}. "
                f"Expected {protein}_snv_counts.csv in {mm2_snv_counts_dir}"
            )

        # read inputs
        df_pairs = pd.read_csv(pair_fp)
        df_snvs = pd.read_csv(snv_fp)

        # basic checks
        for req in ("Variant1", "Variant2"):
            if req not in df_pairs.columns:
                raise ValueError(f"{pair_fp.name} missing required column '{req}'")

        if "snv_variant" not in df_snvs.columns:
            raise ValueError(f"{snv_fp.name} missing required column 'snv_variant'")

        # count_* columns must exist in both files; we only use the intersection
        pair_count_cols = [c for c in df_pairs.columns if c.startswith("count_")]
        snv_count_cols = [c for c in df_snvs.columns if c.startswith("count_")]
        count_cols = [c for c in pair_count_cols if c in set(snv_count_cols)]

        if not count_cols:
            raise ValueError(
                f"No shared count_* columns between:\n"
                f"  pairs: {pair_fp.name} ({pair_count_cols})\n"
                f"  snvs : {snv_fp.name} ({snv_count_cols})"
            )

        # normalize strings
        df_pairs = df_pairs.copy()
        df_pairs["Variant1"] = df_pairs["Variant1"].astype(str).str.strip()
        df_pairs["Variant2"] = df_pairs["Variant2"].astype(str).str.strip()

        df_snvs = df_snvs.copy()
        df_snvs["snv_variant"] = df_snvs["snv_variant"].astype(str).str.strip()

        # make numeric
        for c in count_cols:
            df_pairs[c] = pd.to_numeric(df_pairs[c], errors="coerce").fillna(0)
            df_snvs[c] = pd.to_numeric(df_snvs[c], errors="coerce")

        # build lookup maps for totals: totals_map[col][variant] = total
        totals_map = {
            c: dict(df_snvs[["snv_variant", c]].dropna(subset=["snv_variant"]).values)
            for c in count_cols
        }


        def pct(variant: str, paired_val: float, totals_for_col: dict) -> float:
            v = variant.strip()
            if not v or v == "-":
                return float("nan")
            total = totals_for_col.get(v)
            if total is None or pd.isna(total) or float(total) == 0.0:
                return float("nan")
            return float(paired_val) / float(total)

        # add percent columns for each count column
        for c in count_cols:
            suffix = c.replace("count_", "", 1)  # e.g. "BNP_low_T1T2"
            v1_col = f"Var1%_{suffix}"
            v2_col = f"Var2%_{suffix}"

            totals_for_col = totals_map[c]

            df_pairs[v1_col] = [
                pct(v1, pv, totals_for_col) for v1, pv in zip(df_pairs["Variant1"], df_pairs[c])
            ]
            df_pairs[v2_col] = [
                pct(v2, pv, totals_for_col) for v2, pv in zip(df_pairs["Variant2"], df_pairs[c])
            ]

        # reorder columns exactly like you asked (counts first, then % columns)
        pct_cols = []
        for c in count_cols:
            suffix = c.replace("count_", "", 1)
            pct_cols.extend([f"Var1%_{suffix}", f"Var2%_{suffix}"])

        ordered = ["Variant1", "Variant2"] + count_cols + pct_cols
        df_pairs = df_pairs[ordered]

        out_fp = mm2_hisat_paired_counts_dir / f"{protein}_paired_counts_pct.csv"
        df_pairs.to_csv(out_fp, index=False)
        written.append(str(out_fp))

    if not written:
        raise FileNotFoundError(f"No *_paired_counts.csv files found in {mm2_hisat_paired_counts_dir}")

    return written


def recalc_pct_to_hisat2_counts(mm2_hisat_paired_counts_dir, combined_t1t2_count_files_dir):
    """
    Multiple hisat2 counts with mm2's paired counts' Var1% and Var2% to get recalculated hisat2 counts for paired vars.
    """
    mm2_hisat_paired_counts_dir = Path(mm2_hisat_paired_counts_dir)
    combined_t1t2_count_files_dir = Path(combined_t1t2_count_files_dir)

    pat_pct = re.compile(r"^(?P<protein>.+)_paired_counts_pct\.csv$", re.IGNORECASE)
    pat_maf = re.compile(r"^(?P<protein>.+)_T1T2_MAF_filtered\.csv$", re.IGNORECASE)

    # index maf files by protein
    maf_map: dict[str, Path] = {}
    for p in combined_t1t2_count_files_dir.glob("*_T1T2_MAF_filtered.csv"):
        m = pat_maf.match(p.name)
        if m:
            maf_map[m.group("protein")] = p

    written: list[str] = []

    for pct_fp in sorted(mm2_hisat_paired_counts_dir.glob("*_paired_counts_pct.csv")):
        m = pat_pct.match(pct_fp.name)
        if not m:
            continue
        protein = m.group("protein")

        maf_fp = maf_map.get(protein)
        if maf_fp is None:
            raise FileNotFoundError(
                f"No matching {protein}_T1T2_MAF_filtered.csv found in {combined_t1t2_count_files_dir}"
            )

        df = pd.read_csv(pct_fp)
        maf = pd.read_csv(maf_fp)

        # Required columns
        if not {"Variant1", "Variant2"}.issubset(df.columns):
            raise ValueError(f"{pct_fp.name} must contain Variant1 and Variant2")

        if "Variant" not in maf.columns:
            raise ValueError(f"{maf_fp.name} must contain column 'Variant'")

        # Identify % columns dynamically
        var1_pct_cols = [c for c in df.columns if c.startswith("Var1%_")]
        var2_pct_cols = [c for c in df.columns if c.startswith("Var2%_")]
        if not var1_pct_cols or not var2_pct_cols:
            raise ValueError(f"{pct_fp.name} missing Var1%_/Var2%_ columns")

        # Sanity: keep only suffixes that exist in BOTH Var1% and Var2%
        suffixes1 = {c.replace("Var1%_", "", 1) for c in var1_pct_cols}
        suffixes2 = {c.replace("Var2%_", "", 1) for c in var2_pct_cols}
        suffixes = sorted(suffixes1 & suffixes2)
        if not suffixes:
            raise ValueError(f"{pct_fp.name}: no matching Var1%_/Var2%_ suffix pairs found")

        # Normalize variant strings
        df["Variant1"] = df["Variant1"].astype(str).str.strip()
        df["Variant2"] = df["Variant2"].astype(str).str.strip()
        maf["Variant"] = maf["Variant"].astype(str).str.strip()

        # For each suffix, we need a hisat2 count column named count_{suffix}
        # Example suffix: "BNP_low_T1T2" -> hisat2 column "count_BNP_low_T1T2"
        needed_count_cols = [f"count_{suf}" for suf in suffixes]
        missing_counts = [c for c in needed_count_cols if c not in maf.columns]
        if missing_counts:
            raise ValueError(f"{maf_fp.name} missing required hisat2 count columns: {missing_counts}")
        
        # build lookup maps: count_col -> {Variant: total}
        totals_map: dict[str, dict[str, float]] = {}
        for count_col in needed_count_cols:
            maf[count_col] = pd.to_numeric(maf[count_col], errors="coerce")
            sub = maf[["Variant", count_col]].dropna(subset=["Variant", count_col]).copy()
            totals_map[count_col] = dict(zip(sub["Variant"], sub[count_col]))

        for suf in suffixes:
            v1_pct_col = f"Var1%_{suf}"
            v2_pct_col = f"Var2%_{suf}"
            count_col = f"count_{suf}"

            p1 = pd.to_numeric(df[v1_pct_col], errors="coerce")
            p2 = pd.to_numeric(df[v2_pct_col], errors="coerce")

            t1 = pd.to_numeric(df["Variant1"].map(totals_map[count_col]), errors="coerce")
            t2 = pd.to_numeric(df["Variant2"].map(totals_map[count_col]), errors="coerce")

            v1_new = f"Var1_recalc_hisat2_{suf}"
            v2_new = f"Var2_recalc_hisat2_{suf}"

            df[v1_new] = (p1 * t1).round().fillna(0).astype(int)
            df[v2_new] = (p2 * t2).round().fillna(0).astype(int)


        out_fp = mm2_hisat_paired_counts_dir / f"{protein}_paired_counts_w_recalc_hisat2.csv"
        df.to_csv(out_fp, index=False)
        written.append(str(out_fp))

        ### save a condensed version for proptest ###
        keep_cols = [c for c in df.columns if c in {"Variant1", "Variant2"} or "recalc" in c]
        df_condensed = df[keep_cols].copy()
        
        # rename column names from Var{num}_recalc_hisat2_{protein}_{condition}_T1T2 to calculated_Variant{num}_mpileup_count_{protein}_{condition}
        df_condensed = df_condensed.rename(columns=lambda c: c.replace("_T1T2", ""))
        df_condensed = df_condensed.rename(columns=lambda c: c.replace("Var1_recalc_hisat2", "calculated_Variant1_mpileup_count"))
        df_condensed = df_condensed.rename(columns=lambda c: c.replace("Var2_recalc_hisat2", "calculated_Variant2_mpileup_count"))

        # keep rows where neither "Variant1" nor "Variant2" is "-" to get co_travelling SNV pairs only
        df_co_travelling_vars = df_condensed[(df_condensed["Variant1"] != "-") & (df_condensed["Variant2"] != "-")].reset_index(drop=True)
        df_co_travelling_vars_path = (mm2_hisat_paired_counts_dir / f"{protein}_co_travelling_snv_pairs.csv")
        df_co_travelling_vars.to_csv(df_co_travelling_vars_path, index=False)
        written.append(str(df_co_travelling_vars_path))

        # keep rows where either "Variant1" or "Variant2" is "-" to get single SNV counts only
        df_single_vars = df_condensed[(df_condensed["Variant1"] == "-") | (df_condensed["Variant2"] == "-")].reset_index(drop=True)
        df_single_vars_path = (mm2_hisat_paired_counts_dir / f"{protein}_single_var_reads.csv")
        df_single_vars.to_csv(df_single_vars_path, index=False)
        written.append(str(df_single_vars_path))


    if not written:
        raise FileNotFoundError(f"No *_paired_counts_pct.csv files found in {mm2_hisat_paired_counts_dir}")

    return written


def main():
    combine_reps_by_protein_tile(cx_count_files_dir, combined_reps_count_files_dir)
    merge_tiles_per_protein(combined_reps_count_files_dir, combined_t1t2_count_files_dir)
    cdna_variants_per_sample(mm2_sample_rep_raw_variants_dir, col="cdna_variants")
    merge_mm2_vars_T1T2_counts(mm2_sample_variants_dir, combined_t1t2_mm2_vars_dir)
    merge_T1T2_mm2_mnv_by_conditions(combined_t1t2_mm2_vars_dir)
    filter_mm2_t1t2_by_maf(combined_t1t2_mm2_vars_dir, combined_t1t2_count_files_dir)
    get_variant_pairs(combined_t1t2_mm2_vars_dir, mm2_hisat_paired_counts_dir)
    count_snvs_from_mm2_mnvs(combined_t1t2_mm2_vars_dir, mm2_snv_counts_dir)
    annotate_paired_counts_perc(mm2_snv_counts_dir, mm2_hisat_paired_counts_dir)
    recalc_pct_to_hisat2_counts(mm2_hisat_paired_counts_dir, combined_t1t2_count_files_dir)


if __name__ == "__main__":
    main()
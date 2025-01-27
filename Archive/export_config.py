import yaml
import sys
import getopt

def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def export_variables(config):
    print(f'ROOT_OUTPUT_DIR="{config["root_output_dir"]}"')
    print(f'REFERENCE_FASTA="{config["reference_fasta"]}"')
    print(f'ARR_NUMBER="{config["arr_number"]}"')
    print(f'DEEPVARIANT_SIF="{config["deepvariant_sif"]}"')
    print(f'MODEL_TYPE="{config["model_type"]}"')

def main(argv):
    config_file = ''
    try:
        opts, args = getopt.getopt(argv, "c:", ["config="])
    except getopt.GetoptError:
        print('export_config.py -c <configfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-c", "--config"):
            config_file = arg
    if not config_file:
        print('No config file provided. Use -c to specify the config file.')
        sys.exit(2)

    config = load_config(config_file)
    export_variables(config)

if __name__ == "__main__":
    main(sys.argv[1:])
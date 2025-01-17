import os
import subprocess
import sys

def run_command(command, check=True):
    """
    Runs a shell command and handles errors gracefully.
    """
    try:
        print(f"Running: {command}")
        subprocess.run(command, shell=True, check=check, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command '{command}' failed with exit code {e.returncode}")
        sys.exit(1)

def ensure_plotly_installed():
    """Ensure Plotly is installed in the virtual environment."""
    try:
        import plotly
    except ImportError:
        print("Plotly is not installed. Installing it now...")
        run_command("pip install plotly==5.24.1")

def main():
    # Ensure script is running in a virtual environment
    if not hasattr(sys, 'real_prefix') and not (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
        print("Please activate your virtual environment before running this script.")
        sys.exit(1)

    # Define installation paths
    env_bin = os.path.join(sys.prefix, "bin")

    print(f"Virtual environment detected: {sys.prefix}")
    print(f"Installing tools into: {env_bin}")

    # Upgrade pip and install Python dependencies
    run_command("pip install --upgrade pip")
    run_command("pip install -r requirements.txt")

    # Ensure Plotly is installed
    ensure_plotly_installed()

    # Install Chopper
    run_command("wget https://github.com/wdecoster/chopper/releases/download/v0.7.0/chopper-musl.zip")
    run_command("unzip chopper-musl.zip -d chopper_bin")
    run_command(f"mv chopper_bin/chopper {env_bin}")
    run_command(f"chmod +x {env_bin}/chopper")

    # Install Minimap2
    run_command("wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2")
    run_command("tar -xvf minimap2-2.28_x64-linux.tar.bz2")
    run_command(f"mv minimap2-2.28_x64-linux/minimap2 {env_bin}")

    # Install GMAP
    run_command("wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2024-11-20.tar.gz")
    run_command("tar -xvzf gmap-gsnap-2024-11-20.tar.gz")
    run_command("cd gmap-2024-11-20 && ./configure --prefix=$VIRTUAL_ENV && make && make install && cd ..")

    # Install Samtools
    run_command("wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2")
    run_command("tar -xvjf samtools-1.16.1.tar.bz2")
    run_command("cd samtools-1.16.1 && ./configure --prefix=$VIRTUAL_ENV --disable-lzma --without-curses && make && make install && cd ..")

    # Clean up temporary files
    run_command("rm -f gmap-gsnap-2024-11-20.tar.gz samtools-1.16.1.tar.bz2 chopper-musl.zip minimap2-2.28_x64-linux.tar.bz2")
    run_command("rm -rf gmap-2024-11-20 samtools-1.16.1 chopper_bin minimap2-2.28_x64-linux")

    # Verify installations
    run_command(f"{env_bin}/samtools --version")
    run_command(f"{env_bin}/gmap --version")
    run_command(f"{env_bin}/chopper --version")
    run_command(f"{env_bin}/minimap2 --version")

    print("\nAll tools installed successfully!")

if __name__ == "__main__":
    main()

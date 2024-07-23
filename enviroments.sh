#!/bin/bash
# Dependencies Installation

# Function to check the OS and set the appropriate environment YAML file
detect_os_and_create_env() {
    case "$(uname -s)" in
        Darwin)
            # macOS
            echo "Detected macOS"
            env_file="enviroments_macos.yaml"
            ;;
        Linux)
            # Linux
            echo "Detected Linux"
            env_file="enviroments_linux.yaml"
            ;;
        CYGWIN*|MINGW32*|MINGW64*|MSYS*|MINGW*)
            # Windows
            echo "Detected Windows"
            env_file="enviroments_windows.yaml"
            ;;
        *)
            echo "Unsupported OS"
            exit 1
            ;;
    esac

    # Print the current working directory
    echo "Current working directory: $(pwd)"

    # Add channels to Conda configuration
    echo "Adding channels to Conda config"
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # Create the Conda environment from the YAML file
    echo "Creating Conda environment from ${env_file}"
    conda env create --file "$env_file"

    echo "Environment creation process completed."
}

# Execute the function
detect_os_and_create_env

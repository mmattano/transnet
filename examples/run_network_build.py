"""Example of how to run the build_networks.py script to create a network for a model organism.

This script demonstrates using the maintenance/build_networks.py script to create a network
for the mouse organism in debug mode (which creates a smaller network for testing purposes).
"""

# Import necessary modules
import os
import subprocess
import sys

def main():
    """Run the build_networks.py script with the desired parameters."""
    
    # Get the absolute path to the build_networks.py script
    script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    build_script = os.path.join(script_dir, "maintenance", "build_networks.py")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.join(script_dir, "test_data")
    os.makedirs(output_dir, exist_ok=True)
    
    # Set the command to run
    cmd = [
        sys.executable,  # Use the current Python interpreter
        build_script,
        "--output-dir", output_dir,
        "--organisms", "mouse",
        "--debug"  # Use debug mode for faster execution
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    
    # Execute the command
    try:
        subprocess.run(cmd, check=True)
        print(f"Network successfully built! Check the output in: {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Error building network: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
    
import os
import sys
import numpy as np

# Import the Python vegetation model modules
from src import GAPpyModel, CODENAME, VERSION_ID

def main():
    """Main entry point for the vegetation model."""
    
    print(f"Starting {CODENAME} {VERSION_ID} - Python Implementation")
    print("=" * 80)
    
    # Create and run the model
    model = GAPpyModel()
    
    # You can pass a filelist if needed
    filelist = ""
    if len(sys.argv) > 1:
        filelist = sys.argv[1]
    
    try:
        model.run(filelist)
        print("Model execution completed successfully")
    except Exception as e:
        print(f"Error during model execution: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
import sys
import os
import streamlit.web.cli as stcli

def resolve_path(path):
    """
    Returns the absolute path to a resource, handling PyInstaller's temp paths.
    """
    if getattr(sys, "frozen", False):
        # Running as compiled executable
        basedir = sys._MEIPASS
    else:
        # Running as normal python script
        basedir = os.path.dirname(__file__)
    return os.path.join(basedir, path)

if __name__ == "__main__":
    # Path to the streamlit app script
    app_path = resolve_path("app.py")
    
    # Configure system arguments to mimic "streamlit run app.py"
    sys.argv = [
        "streamlit",
        "run",
        app_path,
        "--global.developmentMode=false",
    ]
    
    # Run the streamlit cli
    sys.exit(stcli.main())

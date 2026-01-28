import sys
import os
import traceback
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
    try:
        print("🚀 Starting PadDx_miR Designer...")
        
        # Path to the streamlit app script
        app_path = resolve_path("app.py")
        
        if not os.path.exists(app_path):
            print(f"❌ Error: App file not found at {app_path}")
            input("Press Enter to exit...")
            sys.exit(1)

        print(f"📂 App path resolved to: {app_path}")
        
        # Configure system arguments to mimic "streamlit run app.py"
        # We also create a .streamlit folder in the temp dir to avoid config warnings
        sys.argv = [
            "streamlit",
            "run",
            app_path,
            "--global.developmentMode=false",
            "--server.headless=true",  # Automatically open browser without asking
        ]
        
        print("⏳ Launching Streamlit server...")
        sys.exit(stcli.main())

    except Exception as e:
        print("\n\n❌ An unexpected error occurred:")
        traceback.print_exc()
        print("\n")
        input("Press Enter to exit...")
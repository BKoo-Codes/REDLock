import sys
import os
import traceback
import time
import webbrowser
import threading
import streamlit.web.cli as stcli

def resolve_path(path):
    if getattr(sys, "frozen", False):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(__file__)
    return os.path.join(basedir, path)

def open_browser():
    """Waits for the server to start and then opens the browser."""
    print("⏳ Waiting for server to initialize...")
    time.sleep(2)  # Give Streamlit a moment to start
    url = "http://localhost:8501"
    print(f"🚀 Opening browser at {url}...")
    webbrowser.open(url)

if __name__ == "__main__":
    try:
        print("🚀 Starting REDLock Designer...")
        app_path = resolve_path("app.py")
        
        if not os.path.exists(app_path):
            print(f"❌ Error: App file not found at {app_path}")
            input("Press Enter to exit...")
            sys.exit(1)

        # Start browser opener in a separate thread
        threading.Thread(target=open_browser, daemon=True).start()
        
        sys.argv = [
            "streamlit",
            "run",
            app_path,
            "--global.developmentMode=false",
            "--server.headless=true",
            "--server.port=8501",
            "--browser.gatherUsageStats=false" # Disable stats collection for cleaner log
        ]
        
        print("💻 Starting Streamlit server...")
        sys.exit(stcli.main())

    except Exception as e:
        print("\n\n❌ An unexpected error occurred:")
        traceback.print_exc()
        input("Press Enter to exit...")

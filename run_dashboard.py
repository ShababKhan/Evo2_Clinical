#!/usr/bin/env python
"""
Dashboard Launcher for Evo2_Clinical

This script launches the Evo2_Clinical visualization dashboard.
Run this script to start the Streamlit-based interactive dashboard.
"""

import os
import sys
from pathlib import Path

# Add the project root to Python path to ensure imports work correctly
project_root = Path(__file__).parent.absolute()
sys.path.insert(0, str(project_root))

# Check for required packages and install if needed
try:
    import streamlit
    import plotly
except ImportError:
    print("Installing required packages...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "streamlit", "plotly"])

from evo2_clinical.visualization.dashboard import run_dashboard

if __name__ == "__main__":
    print("Launching Evo2_Clinical Dashboard...")
    # Use Streamlit's CLI to run the dashboard module
    os.system(f"{sys.executable} -m streamlit run {__file__} -- --server.headless=true")
else:
    # When this file is imported by streamlit, run the actual dashboard
    run_dashboard()
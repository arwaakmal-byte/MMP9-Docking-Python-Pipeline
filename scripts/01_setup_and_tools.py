"""
Stage 1 — Environment setup and tool installation.

Installs RDKit, gemmi, and Meeko, then downloads a prebuilt AutoDock Vina
binary directly from GitHub Releases. The `pip install vina` package is
deliberately NOT used here: on Colab it consistently fails to build
(`egg_info did not run successfully`), most likely due to a missing
system-level Boost/SWIG dependency for its C++ extension. Downloading the
official prebuilt binary sidesteps that entirely and is the approach this
whole pipeline was validated against.
"""

import os
import stat
import subprocess
import sys

import requests

from config import DOCKING_DIR, VINA_BIN_PATH


def mount_drive_if_colab():
    if os.path.isdir("/content/drive/MyDrive"):
        print("Drive already mounted.")
        return
    try:
        from google.colab import drive  # type: ignore

        drive.mount("/content/drive")
    except ImportError:
        print("Not running in Colab — skipping Drive mount. "
              "Make sure MMP9_PROJECT_DIR points to a real, persistent path.")


def install_python_packages():
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "-q",
         "rdkit", "gemmi", "meeko", "scipy"],
        check=True,
    )
    import meeko
    from rdkit import Chem

    print(f"Meeko {meeko.__version__} / RDKit {Chem.rdBase.rdkitVersion} imported successfully.")


def install_vina_binary():
    if os.path.exists(VINA_BIN_PATH):
        print("Vina binary already present.")
    else:
        release_data = requests.get(
            "https://api.github.com/repos/ccsb-scripps/AutoDock-Vina/releases/latest"
        ).json()
        target_asset = next(
            a for a in release_data["assets"]
            if a["name"].endswith("linux_x86_64") and a["name"].startswith("vina_")
        )
        r = requests.get(target_asset["browser_download_url"])
        r.raise_for_status()
        with open(VINA_BIN_PATH, "wb") as f:
            f.write(r.content)
        print(f"Downloaded Vina binary: {target_asset['name']}")

    # Always (re-)set the executable bit — a binary that already existed
    # from a prior Drive session may not have it set in a fresh runtime.
    st = os.stat(VINA_BIN_PATH)
    os.chmod(VINA_BIN_PATH, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

    result = subprocess.run([VINA_BIN_PATH, "--version"], capture_output=True, text=True)
    print(result.stdout.strip())


if __name__ == "__main__":
    mount_drive_if_colab()
    install_python_packages()
    install_vina_binary()

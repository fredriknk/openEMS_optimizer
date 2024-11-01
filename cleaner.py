import os
import shutil
from pathlib import Path
import time

def cleanup_runs_directory(directory_path, keep=3):
    """
    Deletes all but the 'keep' newest folders in the specified directory.

    :param directory_path: Path to the directory containing the folders.
    :param keep: Number of newest folders to keep.
    """
    # Convert to absolute path
    directory = Path(directory_path).resolve()

    # Ensure the directory exists
    if not directory.is_dir():
        print(f"The directory {directory} does not exist.")
        return

    # List all items in the directory and filter out only directories
    items = [item for item in directory.iterdir() if item.is_dir()]

    # Sort the directories by modification time (newest first)
    items.sort(key=lambda x: x.stat().st_mtime, reverse=True)

    # Keep the 'keep' newest folders
    folders_to_keep = items[:keep]
    folders_to_delete = items[keep:]

    if not folders_to_delete:
        print("No folders to delete.")
        return

    # Delete the folders to delete
    for folder in folders_to_delete:
        try:
            shutil.rmtree(folder)
            print(f"Deleted folder: {folder}")
        except Exception as e:
            print(f"Could not delete folder {folder}: {e}")

if __name__ == "__main__":
    # Specify the path to the 'runs' directory
    runs_directory = "./runs"
    # Call the cleanup function
    while True:
        #print the current time HH:MM:SS - DD/MM/YYYY
        print(time.strftime("%H:%M:%S - %d/%m/%Y", time.localtime(time.time())))
        cleanup_runs_directory(runs_directory, keep=5)
        time.sleep(600)
        
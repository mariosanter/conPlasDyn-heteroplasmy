import os
import warnings
from git import Repo

def _get_git_commit():
    try:
        repo_dir = os.path.dirname(os.path.abspath(__file__))
        repo = Repo(repo_dir, search_parent_directories=True)
        commit = repo.head.object.hexsha[:7]

        # Only check untracked files inside "module"
        module_path = "module/"
        untracked_in_module = [
            f for f in repo.untracked_files if f.startswith(module_path)
        ]

        if untracked_in_module:
            warnings.warn("Repository has untracked files in module!", UserWarning)
            commit += "-dirty"

        return commit
    except Exception:
        return "unknown"

__version__ = "0.1"
__commit__ = _get_git_commit()

print(f"Loaded ConPlasDyn package v{__version__} (commit {__commit__})")
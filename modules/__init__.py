import os
import warnings
from git import Repo

def _get_git_commit():
    try:
        repo_dir = os.path.dirname(os.path.abspath(__file__))
        repo = Repo(repo_dir, search_parent_directories=True)
        commit = repo.head.object.hexsha[:7]

        # Check if repo has uncommitted changes
        if repo.is_dirty(untracked_files=True):
            warnings.warn("Repository has uncommitted changes!", UserWarning)
            commit += "-dirty"

        return commit
    except Exception:
        return "unknown"

__version__ = "0.1"
__commit__ = _get_git_commit()

print(f"Loaded ConPlasDyn package v{__version__} (commit {__commit__})")
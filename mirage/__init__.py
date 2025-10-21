import re
import sys

try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # For Python < 3.8, fallback to backport
    from importlib_metadata import version, PackageNotFoundError

__version_commit__ = ''
_regex_git_hash = re.compile(r'.*\+g(\w+)')

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = 'dev'

if '+' in __version__:
    match = _regex_git_hash.match(__version__)
    if match:
        __version_commit__ = match.groups()[0]
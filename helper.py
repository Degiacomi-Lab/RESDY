

import shutil

def _is_tool(name):
    """Check whether `name` is on PATH and marked as executable.
    Parameters:
    -----------
    name : string
        name of tool, e.g. wget or curl

    Returns:
    --------
    is_tool : bool
        True if tool exists

    """

    return shutil.which(name) is not None

def get_download_tool():

    if _is_tool('wget'):
        return 'wget'

    elif _is_tool('curl'):
        return 'curl'

    else:
        return None
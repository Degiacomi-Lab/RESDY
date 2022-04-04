

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
    """Get the name of the tool used for commandline download
    Returns:
    --------
        tool_name : string
            Name of the tool currently wget and curl supported if fails None returned

    """

    if _is_tool('wget'):
        return 'wget'

    elif _is_tool('curl'):
        return 'curl'

    else:
        return None
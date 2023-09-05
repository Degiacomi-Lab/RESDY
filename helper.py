import shutil
import sys, os


def get_download_tool():
    """Get the name of the tool used for commandline download
    Returns:
    --------
        tool_name : string
            Name of the tool currently wget and curl supported if fails None returned

    """

    if shutil.which('wget') is not None:
        return 'wget'

    elif shutil.which('curl') is not None:
        return 'curl'

    else:
        return None
    
    
class ShutUp(object):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, *args):
        sys.stdout.close()
        sys.stdout =  self._stdout
import os
import subprocess
end_url = '4HHB.A'
web_url = "https://www.rcsb.org/fasta/chain/" + end_url + '/download'
subprocess.check_call("wget -O newfile.txt " + web_url, shell=True)
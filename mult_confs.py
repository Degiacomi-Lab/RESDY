import re
import glob
import os
import numpy as np

def split_confs(pdb):
    print('Checking for alternate conformations')
    try:
        f = open(pdb + '.pdb')
        number = 1
        name = pdb + '_alt_' + str(number) + '.pdb'
    except Exception as e:
        print('Error %s'%e)
        print('Print failure opening file.')

    f_write = open(name, 'w')
    for line in f:
        try:
            line = str(line)

            if re.search('MASTER', line):
                os.remove(name)
                break

            else:

                if re.search('ENDMDL', line):
                    f_write.write(line)
                    f_write.close()
                    number = number + 1
                    name = pdb + '_alt_' + str(number) + '.pdb'
                    f_write = open(name, 'w')

                else:
                    f_write.write(line)
        except:
            continue
    return()
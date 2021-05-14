#! /usr/bin/python3
# ! coding=utf-8


import os
import sys
from libs import *


if __name__ == '__main__':
    local_path = os.path.split(os.path.realpath(__file__))[0]
    cutadapt_path = '{}/bin/cutadapt'.format(local_path)
    sys.path.append(cutadapt_path)
    main(local_path)

import os
from utils import *

def GenerateResult(tablefile, fastqcfile, fig):
    os.system('mkdir RESULT/qc')
    os.system('cp Fastqc/*.html RESULT/qc')


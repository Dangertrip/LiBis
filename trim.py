'''
Tested
'''
from utils import *
import os
class Trim():
    
    def check(self):
        if not toolcheck('trim_galore -v'):
            return False,'trim_galore not found!'
        if os.path.exists('Trim'):
            return False,'"Trim" file/dic exist! Please delete it.'
        os.mkdir('Trim')
        return True,''

    def setpath(self,path):
        self.path = path+'Trim'

    def run(self,filename):
        p = Pshell('trim_galore --gz -o '+self.path+' '+filename)
        p.process()

if __name__=='__main__':
    tr = Trim()
    #tr.check()
    tr.setpath('./')
    tr.run('../trimtest/SRR1248444_1.1.1.1.fastq')



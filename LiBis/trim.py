'''
Tested
'''
from .utils import *
import os
class Trim():
    
    def check(self,nocheck=False):
        if not toolcheck('trim_galore -v'):
            return False,'trim_galore not found!'
        if os.path.exists('Trim'):
            if nocheck:
                print('"Trim" file/dic exist! But --nocheck enabled, so continue running')
            else:
                return False,'"Trim" file/dic exist! Please delete it.'
        else:
            os.mkdir('Trim')
        return True,''

    def setpath(self,path):
        self.path = path+'Trim'

    def run(self,filename,pair):
        if pair==1:
            p = Pshell('trim_galore --gz -o '+self.path+' '+filename)
        else:
            p = Pshell('trim_galore --gz --paired -o '+self.path+' '+filename)
        p.process()

if __name__=='__main__':
    tr = Trim()
    #tr.check()
    tr.setpath('./')
    tr.run('../trimtest/SRR1248444_1.1.1.1.fastq')



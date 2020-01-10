'''
Tested
'''
from .utils import *
import os

class Fastqc():


    def check(self,nockeck=False):
        #return True,''
        if not toolcheck('fastqc --help'):
            return False,'Fastqc command not found'
        if os.path.exists('Fastqc'):
            if nockeck:
                print('Fastqc file or dir exists! But --nockeck enabled, so continue running')
            else:
                return False,'Fastqc file or dir exists'
        else:
            os.mkdir('Fastqc')
        return True,''

    def setpath(self,path):
        self.path = path+'Fastqc'

    def run(self,filename):
        pshell=Pshell('fastqc -o '+self.path+' '+filename)
        pshell.process()

    

if __name__=="__main__":
    a = Fastqc()
    #print(a.check())
    a.setpath('./')
    a.run('../trimtest/SRR1248444_1.1.1.1.fastq')

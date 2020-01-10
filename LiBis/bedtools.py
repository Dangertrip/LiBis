'''
Tested
'''
import os
from .utils import *

class Bedtools:

    def check(self,nockeck=False):
        if not toolcheck('bedtools --version'):
            return False,'Bedtools not found!'
        return True,''

    def setparam(self,param):
        self.genome = param['genome']
        if 'bin' in param and param['bin']!=None:
            self.bin = param['bin']
        else:
            self.bin=1000000

    def makewindow(self):
        path=os.path.abspath(__file__)
        path = path[:path.rfind('/')+1]
        filename = path+'chromsize/'+self.genome + '.chrom.sizes'
        print(filename)
        outputname = self.genome+'_'+str(self.bin)+'.bed'
        os.system('bedtools makewindows -g '+filename+' -w '+str(self.bin)+' > '+outputname)
        self.binfile = outputname

    def intersect(self,names):
        sample=0
        result=[]
        for name in names:
            temp = str(sample)+'.intersect'
            output = str(sample)+'.bed'
            os.system('bedtools intersect -loj -a '+self.binfile+' -b '+ name +" | awk -v OFS='\t' '{print $1,$2,$3,$7}' > BED_FILE/"+temp)
            os.system('bedtools groupby -i BED_FILE/'+temp+' -g 1,2,3 -c 4 -o mean > '+'BED_FILE/'+output)
            result.append('BED_FILE/'+output)
            sample=sample+1
        return result

    '''
    Two plots from segmented genome average: TSNE/PCA and heatmap
    '''
if __name__=="__main__":
    b = Bedtools()
    param={'genome':'hg19'}
    b.setparam(param)
    print(b.check())
    b.makewindow()
    name = ['FWAC.bed']
    print(b.intersect(name))
    print(b.binfile)

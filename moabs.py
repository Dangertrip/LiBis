'''
Tested
This is a wapper for MOABS in this pipeline. It contains commands for BSMAP, bam sort, bam deduplicate, methylation calling.
'''

from utils import *
from clipmode import clipmode
import os

class Mcall():

    def check(self):
        #return True,''
        if not toolcheck('mcall'):
            return False,'Mcall not found!'
        if os.path.exists('BED_FILE'):
            return False,'"BED_FILE" exists! Please delete "BED_FILE"'
        os.mkdir("BED_FILE")
        return True,''

    def setpath(self,path):
        self.path = path+'BED_FILE/'

    def setparam(self,param):
        mcallparam=None
        if 'mcall' in param:
            mcallparam=param['mcall']
        self.refpath = param['ref']
        self.extraparam=''

    def run(self,name):
        '''
        4 muted 37 41 46 49
        '''
        cmd='mcall -m '+name+' -r '+self.refpath+' -p 8 1>>'+self.path+'log 2>>'+self.path+'err'
        p = Pshell(cmd)
        p.process()
        generatedfile=['.G.bed','.HG.bed','_stat.txt']
        newname=[]
        for g in generatedfile:
            os.rename(name+g,self.path+name[name.rfind('/')+1:]+g)
            newname.append(self.path+name[name.rfind('/')+1:]+g)
       
        cmd="awk -v OFS='\t' '{if (NR!=1) print $1,$2,$3,$4}' "+newname[0]+'> '+newname[0]+'.short.bed'
        p.change(cmd)
        p.process()
        cmd="rm mSuite.G.bed"
        p.change(cmd)
        p.process()

        return newname[0],newname[2]

        



class Bsmap():

    def check(self):
        #return True,''
        if not toolcheck('bsmap -h'):
            return False,'BSMAP not found!'
        if os.path.exists('BAM_FILE'):
            return False,'"BAM_FILE" exists! Please delete "BAM_FILE"'
        os.mkdir("BAM_FILE")
        return True,''
   
    def setpath(self,path):
        self.path = path+'BAM_FILE/'

    def setparam(self,param):
        self.extraparam=''
        bsmapparam=None
        if 'bsmap' in param:
            bsmapparam = param['bsmap']
        '''
        Once I'm ready to add all parameters for bsmap and mcall (After finish the dictionary 
        to parameter function, I will delete variable named ***param )
        '''
        self.refpath = param['ref']

    def normalmode(self,file,param={}):
        '''
        3 muted 96 98 100
        '''
        #f = file.strip().split()
        f = file
        purename = RemoveFastqExtension(f[0][f[0].rfind('/')+1:])
        name = self.path+purename+'.bam'
        logname = self.path+purename+'.record'
        if len(f)==1:
            cmd = 'bsmap -a '+f[0]+' -d '+self.refpath+' -o '+name+' -n 0 1>>BAM_FILE/bsmap_log 2>'+logname
        else:
            cmd = 'bsmap -a '+f[0]+' -b '+f[1]+' -d '+self.refpath+' -o '+name+' -n 0 1>>BAM_FILE/bsmap_log 2>'+logname
        p = Pshell(cmd)
        p.process()
        p.change('samtools sort -@ 4 '+name+' -o '+name+'.sorted.bam')
        p.process()
        p.change('mv '+name+'.sorted.bam '+name)
        p.process()
        return name,logname

    def clipping(self,filenames,param={}):
        '''
        I should return a bam file name and a log file name here
        '''
        newname,log = clipmode(filenames,param)
        return newname,log

if __name__=="__main__":
    #mcall = Mcall()
    #print(mcall.check())
    #bsmap = Bsmap()
    #bsmap.setparam({'ref':'/data/dsun/ref/humanigenome/hg19.fa'})
    #print(bsmap.check())
    #bsmap.setpath('./')
    #bsmap.normalmode(['Trim/head.fq'])
    mcall = Mcall()
    mcall.setpath('./')
    mcall.setparam({'ref':'/data/dsun/ref/humanigenome/hg19.fa'})
    mcall.run('BAM_FILE/head_combine.bam')

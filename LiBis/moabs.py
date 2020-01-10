'''
Tested
This is a wapper for MOABS in this pipeline. It contains commands for BSMAP, bam sort, bam deduplicate, methylation calling.
'''

from .utils import *
from .clipmode import clipmode
import os

class Mcall():

    def check(self,nocheck=False):
        #return True,''
        if not toolcheck('mcall'):
            return False,'Mcall not found!'
        if os.path.exists('BED_FILE'):
            if nocheck:
                print('"BED_FILE" exists! But --nocheck enabled, so continue running.')
            else:
                return False,'"BED_FILE" exists! Please delete "BED_FILE"'
        else:
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

    def __init__(self):
        self.samtools_version = 1.1#samtoolsversion()
    
    def samtools_sort(self,p, inputfile, outputfile):
        if self.samtools_version<=1.3:
            p.change('samtools sort -f -@ 4 '+inputfile+' '+outputfile)
            p.process()
        else:
            p.change('')
            p.process()
    
    def check(self,nocheck=False):
        #return True,''
        if not toolcheck('bsmap -h'):
            return False,'BSMAP not found!'
        if os.path.exists('BAM_FILE'):
            if nocheck:
                print('"BAM_FILE" exists! But --nocheck enabled, so continue running.')
            else:
                return False,'"BAM_FILE" exists! Please delete "BAM_FILE"'
        else:
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

    def normalmode(self,file,given_label,threads,param={}):
        '''
        3 muted 96 98 100
        '''
        #f = file.strip().split()
        f = file
        purename = given_label#RemoveFastqExtension(f[0][f[0].rfind('/')+1:])
        name = self.path+purename+'.bam'
        logname = self.path+purename+'.record'
        if len(f)==1:
            cmd = 'bsmap -a '+f[0]+' -d '+self.refpath+' -o '+name+' -S 123 -n 1 -p ' + threads + ' 1>>BAM_FILE/bsmap_log 2>'+logname
        else:
            cmd = 'bsmap -a '+f[0]+' -b '+f[1]+' -d '+self.refpath+' -o '+name+' -S 123 -n 1 -p ' + threads + ' 1>>BAM_FILE/bsmap_log 2>'+logname
        p = Pshell(cmd)
        p.process()
        p.change('samtools sort -f -@ 4 '+name+' '+name+'.sorted.bam')
        p.process()
        p.change('mv '+name+'.sorted.bam '+name)
        p.process()
        return name,logname

    def clipping(self,filenames, param, given_bam_file,given_label):
        '''
        I should return a bam file name and a log file name here
        '''
        newname,log = clipmode(filenames,param, given_bam_file,given_label)
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

import sys
import multiprocessing
import subprocess
import os
import time
from datetime import datetime

timer = str(int(time.time()))


class Pshell():

    def __init__(self,command):
        self.command = command
        self.out = ""

    def change(self,command):
        self.command = command

    def getcommand(self,command):
        return command

    def process(self):
        t = subprocess.Popen(self.command,shell=True,stdout = subprocess.PIPE,stderr = subprocess.PIPE)
        self.out = ''#t.stdout.read().decode()
        self.err = t.stderr.read().decode()
        with open('LiBis_command_log_'+timer+'.txt','a') as f: 
            f.write(str(datetime.now())+'\tcommand: '+self.command+'\n')
            f.write('Output:\n')
            f.write(self.err+'\n\n')

    def get_result(self):
        return self.out

class reads():

    def __init__(self,chromosome,startpos,order,mismatch):
       # self.fqname = fqname
        self.chromosome = chromosome
        self.startpos = startpos
        self.sum = 0
        self.order = order
        self.chain = 0
        self.mismatch = mismatch

    def __str__(self):
        print(self.__dict__)

    def equal(self,rreads):
        if (rreads.getChrom()==self.chromosome)and(rreads.getStartPos()==self.startpos): return True
        return False

    def getChrom(self):
        return self.chromosome

    def getOrder(self):
        return self.order

    def getSum(self):
        return self.sum

    def getStartPos(self):
        return self.startpos

    def join(self,tail,mis):
        self.mismatch = self.mismatch + mis
        self.sum = tail.getOrder()-self.getOrder()

    def getChain(self):
        return self.chain

    def getMismatch(self):
        return self.mismatch

    def canjoin(self,tail,step,xx,binsize):
        if (tail.getChrom()==self.chromosome):
            if (tail.getStartPos() == self.startpos + step*(self.order-tail.getOrder()) or 
            tail.getStartPos() == self.startpos - step*(self.order-tail.getOrder())):
            #if tail.getStartPos()>self.startpos: self.chain = 0
            #else: self.chain = 1
                return True
            if xx!=binsize:
                if xx % step == abs(self.startpos-tail.getStartPos())%step:
                    return True
        return False


def readsjoin(head,tail,step,xx,binsize):
    
    if (head[0]==tail[0]):
        if not overlap(head,tail,step,binsize): return False
        if (tail[1]==head[1]+step*(head[2]-tail[2]) or  tail[1]==head[1]-step*(head[2]-tail[2])):
            return True
        if xx!=binsize:
            if xx % step == abs(head[1]-tail[1])%step:
                return True
    return False
#[s[2][3:],int(s[3]),file_order,mismatch,0] [chr,startpos,fileorder,mismatch,0]
#chromosome,startpos,order,mismatch,sum
def overlap(a,b,step,length_bin):
#Detect 2 kinds of overlapping
# 1. fragments combination overlap: fragment 1 and fragment 2,3,4.. can not appear in two different combinations
# 2. combinations can not overlap in reference genome.
    ordera = a[2]
    orderb = b[2]
    if (orderb<ordera):
        k=a
        a=b
        b=k
    ordera = a[2]
    orderb = b[2]
    sa = a[1]
    sb = b[1]
    suma = a[4]
    sumb = b[4]
    if (a[0]==b[0]) and ((sa<sb and sa+suma*step+length_bin>sb) or (sb<sa and sb+sumb*step+length_bin>sa)): 
        return True
    if orderb-ordera-suma<=(length_bin/step)-1: return True
    return False

def RemoveFastqExtension(name):
    '''
    Tested
    This function is used for remove extensions like .gz .fq .fastq.gz
    '''
    newname = name
    if newname[-3:].lower()=='.gz':
        newname=newname[:-3]
    if newname[-3:].lower()=='.fq':
        newname=newname[:-3]
    if newname[-6:].lower()=='.fastq':
        newname=newname[:-6]
    return newname

def toolcheck(command):
    '''
    Tested
    '''
    t = subprocess.Popen(command,shell=True,stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out = t.stdout.read().decode()
    err = t.stdout.read().decode()
    if ' command not found' in out or ' command not found' in err:
        return False
    return True

def samtoolsversion():
    t = subprocess.Popen('samtools --version',shell=True,stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out = t.stdout.read().decode()
    err = t.stdout.read().decode()
    version = out.split('\n')[0].split()[1].strip()
    version = float(version)
    return version

def readf(filename):
    '''
    Tested
    '''
    with open(filename) as f:
        lines = f.readlines()
    return lines

def union(files):
    '''
    Tested
    '''
    dic={}
    for f in files:
        lines = readf(f)
        for line in lines:
            temp = line.strip().split()
            if temp[3]=='.': continue
            key=temp[0]+'\t'+temp[1]+'\t'+temp[2]
            if key in dic:
                dic[key].append(float(temp[3]))
            else:
                dic[key]=[float(temp[3])]
    del_arr=[]
    for key in dic:
        if len(dic[key])!=len(files):
            del_arr.append(key)
    for k in del_arr:
        del(dic[k])
    return dic
            
def exist(file):
    for f in file:
        if isinstance(f,list):
            for ff in f:
                if not os.path.exists(ff):
                    return ff,False
        else:
            if not os.path.exists(f):
                return f,False
    return '',True


def removeFileIfExist(name):
    if os.path.exists(name):
        os.remove(name)


def bam_file_name_format(s):
    return s.strip().split(',')


def fastq_file_name_format(s):
    #split samples using blank. split double end files using comma(,);
    ans=[]
    for f in s:
        if ',' in f:
            ff=f.strip().split(',')
        else:
            ff = [f]
        ans.append(ff)
    return ans



if __name__=="__main__":
    #dic=union(['BED_FILE/head_combine.bam.G.bed.short.bed','BED_FILE/head_combine.bam.G.bed.short.bed'])
    #print(dic)
    #print(RemoveFastqExtension('SRR1248444_2.fastq.gz'))
    #files=[]
    #for i in range(15):
    #    files.append('../example/scWGBS/BED_FILE/'+str(i)+'.bed')
    #data=[]
    #dic=union(files)
    #columns = ['chrom','start','end']
    #print(methdic)
    #methdata=list(map(lambda x:x.split()+dic[x],dic))
    #columns.extend(list(map(lambda x:'F'+str(x),list(range(1,len(methdata[0])-2)))))
    #columns.extend(['MII','MII','MII','MII','MII','2iESC','2iESC','2iESC','2iESC','2iESC','Ser_ESC','Ser_ESC','Ser_ESC','Ser_ESC','Ser_ESC'])
    #print(methdata)
    #print(columns)
    #from pandas import DataFrame
    #from bsplot import *
    #df = DataFrame(methdata,columns=columns)
    #point_cluster(df,'point_clutser.png')
    print(samtoolsversion())

from .moabs import *
from .fastqc import Fastqc
from .trim import Trim
from .utils import *
from .NumExtractor import *
from .bedtools import Bedtools
from pandas import DataFrame
from .bsplot import *
import os

def computeQC(name,fastqc,param={}):
    '''
    Computing function for qc. Return two arrays: The first one contains filenames for next step computing.
    The second one contains 
    '''
    fn=[]
    filename=''
    for nn in name:
        temp=[]
        for n in nn:
            filename = filename+n+' '
            newname = RemoveFastqExtension(n)
            temp.append(newname+'_fastqc.html')
        fn.append(temp)
    fastqc.run(filename)
    return name,fn

def computeTrim(name,trim,param={}):
    fn=''
    trimmedname=[]
    triminfo=[]
    print(name)
    for n in name:
        temp=[]
        infotemp=[]
        fn=''
        fileorder=1
        for nn in n:
            newname = RemoveFastqExtension(nn)
            fn = fn+nn+' '
            if len(n)==1:
                temp.append('Trim/'+newname+'_trimmed.fq.gz')
            else:
                temp.append('Trim/'+newname+'_val_'+str(fileorder)+'.fq.gz')
                fileorder+=1
            infotemp.append('Trim/'+nn+'_trimming_report.txt')
        trimmedname.append(temp)
        triminfo.append(infotemp)
        fn.strip()
        if len(n)==1:
            trim.run(fn,1)
        else:
            trim.run(fn,2)
    name=trimmedname
    return name,triminfo

def computeBsmap(name,bsmap,param):
    #return ['BAM_FILE/6P1_notrim_val_1_val_1_combine.bam','BAM_FILE/6P1_notrim_val_1_val_1_combine.bam'],[['BAM_FILE/6P1_notrim_val_1_val_1_originallog.record','BAM_FILE/6P1_notrim_val_1_val_1_split_log.record'],['BAM_FILE/6P1_notrim_val_1_val_1_originallog.record','BAM_FILE/6P1_notrim_val_1_val_1_split_log.record']]
    if not exist(name):
        raise "Fastq not found"
    clip = param['clip']
    given_bam_files = param['bamfiles']
    labels = param['label']
    bamname=[]
    bamlogname=[]
    threads = param['threads']
    for n,given_bam,given_label in zip(name,given_bam_files,labels):
        if clip:
            newname,logname = bsmap.clipping(n,param,given_bam,given_label)
        else:
            if given_bam!='':
                logname = ''
                newname = given_bam
            else:
                newname,logname = bsmap.normalmode(n,given_label,threads)
        bamname.append(newname)
        bamlogname.append(logname)
    return bamname,bamlogname

def computeMcall(name,mcall,param):
    if not exist(name):
        raise "BAM not found"
    bedname=[]
    logname=[]
    # print(name)
    for n in name:
        # print(n)
        #bn="BED_FILE/6P1_notrim_val_1_val_1_combine.bam.G.bed"
        #ln="BED_FILE/6P1_notrim_val_1_val_1_combine.bam_stat.txt"
        bn,ln=mcall.run(n)
        bedname.append(bn)
        logname.append(ln)
    return bedname,logname

def computeProcess(param):
    fastqc,trim,bsmap,mcall,bedtools = Fastqc(), Trim(), Bsmap(), Mcall(), Bedtools()
    ProcessObject=[fastqc,trim,bsmap,mcall,bedtools]

    if not param['moabs']:
        for obj in ProcessObject:
            status,word=obj.check(param['nocheck'])
        #if not status:
        #    raise Exception(word)
    
    #check all software have been installed successfully.
    for obj in ProcessObject[:-1]:
        obj.setpath('./')
    
    bsmap.setparam(param)
    mcall.setparam(param)

    name = param['name']
    resultfilename={}
    if param['qc']:
        name,qcresult=computeQC(name,fastqc,{}) 
        resultfilename['qc']=qcresult
        
        if os.path.exists('RESULT') and not os.path.exists('RESULT/qc'):
            os.system('mkdir RESULT/qc')
        os.system('cp Fastqc/*.html RESULT/qc')
    
    if param['trim']:
        name,trimresult = computeTrim(name,trim,{})
        resultfilename['trim']=trimresult

    name,bsmapresult = computeBsmap(name,bsmap,param) 
    resultfilename['bsmap'] = bsmapresult
    #print(resultfilename) 
    #resultfilename={}
    #resultfilename['trim']=
    #resultfilename['bsmap']
    #name=['BAM_FILE/6P1_notrim_val_1_val_1_combine.bam']
    
    if 'mcall' not in param or not param['mcall']:
        return
    
    name,mcallresult = computeMcall(name,mcall,param)
    resultfilename['mcall'] = mcallresult
    meth_cpg_bed_name=name
    '''
    All computing job is done here. Plotting is behind.
    Table:
    Filename, Label, trim ratio(TrimOutputExtractor),input reads,mapped reads, unique mapped reads,
    clipped reads, unique clipped reads, all mapped reads, unique all, mapping ratio, unique mapping ratio
    
    Plots:
    1. Extract trim ratio from trim report files. Generate a bar plot
    2. Extract bsmap mapping ratio 
    3. 
    
    '''
    if 'mcall' not in resultfilename or 'bsmap' not in resultfilename:
        return
    if 'plot' not in param or not param['plot']:
        return

    originalfilename = param['name']
    filelabel = param['label']
    marker=['Filename','Label']
    marker.append('QC')
    if param['qc']:
        qcresult = []
        prepath = "qc/"
        for q in resultfilename['qc']:
            if len(q)==1:
                qcresult.append('<a href="'+prepath+q[0]+'">QC'+'</a>')
            else:
                qcresult.append('<a href="'+prepath+q[0]+'">QC1'+'</a>,'+'<a href="'+prepath+q[1]+'">QC2'+'</a>')
    else:
        qcresult=["/" for i in range(len(originalfilename))]
    marker.append('Trim')
    if param['trim']:
        trimresult =TrimOutputExtractor(resultfilename['trim'])#,param['clip'])
    else:
        trimresult=["/" for i in range(len(originalfilename))]

    '''
    Result from BsmapResult:
    [[total reads,mapped reads,uniquely mapped reads, clipped reads, unique clipped reads,
    all mapped reads, all uniquely mapped reads, mapping ratio, uniquely mapping ratio],...]
    '''
    marker.extend(['Input reads','mapped reads','uniquely mapped reads','clipped reads', \
            'uniquely clipped reads','all mapped reads','all uniquely mapped reads',\
            'mapping ratio','uniquely mapping ratio'])
    bsmapresult = BsmapResult(resultfilename['bsmap'],param['clip'],param['bamfiles'])#All contents required by the last extend

    sample=0
    datatable=[marker]
    #print(trimresult)
    for orin in originalfilename:
        l = filelabel[sample]
        br = bsmapresult[sample]
        temp=[]
        filenum=0
        if len(orin)==1:
            nn=orin[0]
        else:
            nn=orin[0]+','+orin[1]
        temp.extend([nn,l])
        #if param['qc']:
        temp.append(qcresult[filenum])
        #if param['trim']:
        trim_s = trimresult[sample]
        if trim_s=="/" or len(trim_s)==1:
            temp.append(trimresult[sample][0])
        else:
            temp.append(trimresult[sample][0]+','+trimresult[sample][1])
        filenum=filenum+1
        temp.extend(br)
        sample = sample+1
        datatable.append(temp)
        

    datatableprint = list(map(lambda x: ((''.join(list(map(lambda y:str(y)+'\t',x )))).strip()+'\n'),datatable))
    with open('RESULT/datatable.txt','w') as f:
        f.writelines(datatableprint)
    '''
    Table generated as RESULT/datatable.txt
    '''
    #print(filelabel)
    if param['genome']!=None and len(filelabel)>1:
        bedtools.setparam(param)
        bedtools.makewindow()
        shortnames=list(map(lambda x:x+'.short.bed',meth_cpg_bed_name))
        intersectnames=bedtools.intersect(meth_cpg_bed_name)
        methdic=union(intersectnames)
        columns = ['chrom','start','end']
        #print(methdic)
        methdata=list(map(lambda x:x.split()+methdic[x],methdic))
        #columns.extend(list(map(lambda x:'F'+str(x),list(range(1,len(methdata[0])-2)))))
        columns.extend(filelabel)
        #print(methdata)
        #print(columns)
        df = DataFrame(methdata,columns=columns)
        point_cluster(df,'RESULT/point_cluster.png')
        heatmap(df,'RESULT/heatmap.png')
    abspath = os.path.abspath(__file__)
    pos=abspath.find('computeProcess')
    abspath = abspath[:pos]
    os.system('cp '+abspath+'result.html RESULT/')
    os.system('cp -r '+abspath+'jslib RESULT/')

        
    
if __name__=="__main__":
    param={'name':[['6P1_notrim_val_1.fq.gz','6P2_notrim_val_2.fq.gz'],['6P1_notrim_val_1.fq.gz','6P2_notrim_val_2.fq.gz']], 
           'clip':1, 
           'label':['sample1','sample2'], 
            'window':30, 
            'step':5, 
            'process':1, 
            'ref':'/data/dsun/ref/humanigenome/hg19.fa', 
            'qc':1,
            'trim':1,
            'genome':'hg19',
            'bin':1000000,
            'filter_len':40
        }
    computeProcess(param)   
         
        

        

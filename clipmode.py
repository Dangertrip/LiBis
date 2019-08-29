from mapreduce import *
import multiprocessing
#from Comb_fastq import combine
from utils import *
import sys
import gzip
import os
def cut(step,length_bin,link,i):
    start = step*i
    end = length_bin + start
    if end > len(link):
        end = len(link)
    if start>=end or end-start<=length_bin-step:
        return False,''
    return True,link[start:end]

def writefiles(UnmappedReads,step,length_bin,max_length,outputname):
    # Part_Fastq_Filename = []
    unmapped_file = outputname + '.unmapped.fastq'
    for i in range(max_length):
        filecontent = []
        for readsname in UnmappedReads:
            link = UnmappedReads[readsname]
            mark,cutreads = cut(step,length_bin,link[0],i)
            if not mark: continue
            _,cutquality = cut(step,length_bin,link[1],i)
            filecontent.append('@'+readsname+'&'+str(i+1)+'\n'+cutreads+'\n+\n'+cutquality+'\n')
        if len(filecontent)==0: break
        # name = outputname+'.part'+str(i+1)+'.fastq'

        # Part_Fastq_Filename.append(name)
        with open(unmapped_file,'a') as f:
            f.writelines(filecontent)
    # print(Part_Fastq_Filename)
    return unmapped_file

def clip_process(inputfileinfo,param,given_bam_file):
    # print(inputfileinfo)
    input_file_num = len(inputfileinfo)
    if input_file_num<=0 or input_file_num>2:
        print("Parameter error in "+l)
        sys.exit()
    refpath = param['ref']
    if 'filter_len' in param:
        filter = param['filter_len']
    else:
        filter = 40
    purename = inputfileinfo[0][inputfileinfo[0].rfind('/')+1:]
    
    if not given_bam_file:
        outputname = RemoveFastqExtension(purename)
    else:
        purename_bam = given_bam_file[given_bam_file.rfind('/')+1:]
        outputname = RemoveFastqExtension(purename)+'_'+purename_bam[:-4]


    # Part_Fastq_Filename=[]
    unmapped_file = outputname + '.unmapped.fastq'
    os.system('rm '+unmapped_file)
    # for i in range(24):
    #     name = outputname+'.part'+str(i+1)+'.fastq'
    #     if os.path.exists(name):
    #         os.system('rm '+name)

    if not given_bam_file:
        originallogname = 'BAM_FILE/'+outputname+'_originallog.record'
    else:
        originallogname = ''

    phred=33
    if not given_bam_file:
        if input_file_num==2 :
            command='bsmap -a '+inputfileinfo[0]+' -b '+inputfileinfo[1]+' -z '+str(phred)+' -d '+refpath+' -o '+outputname+'.bam -S 123 -n 1 -r 0 1>>log 2>>'+originallogname
        else:
            command='bsmap -a '+inputfileinfo[0]+' -z '+str(phred)+' -d '+refpath+'  -o '+outputname+'.bam -S 123 -n 1 -r 0 1>>log 2>>'+originallogname
        First_try = Pshell(command)
        First_try.process()

#Test1 done
    if not given_bam_file:
        command = 'samtools view '+outputname+'.bam > '+outputname+'.sam'
    else:
        command = 'samtools view '+given_bam_file+' > '+outputname+'.sam'
    BamFileReader = Pshell(command)

    BamFileReader.process()
    
    
    with open(outputname+".sam") as sam:
    #second column in sam file: 64, mate 1; 128, mate 2;
        samlines = sam.readlines()
    set_sam = {}
    
    for line in samlines:
        temp = line.strip().split()
        m1 = (int(temp[1]) & 64)
        m2 = (int(temp[1]) & 128)
        if m1>0: mate = 1
        elif m2>0: mate = 2
        else: mate = 0
        if temp[0] in set_sam: set_sam[temp[0]]=3
        else: set_sam[temp[0]]=mate 
    del samlines
    
    UnmappedReads = {}
    o=0
    step = param['step']
    length_bin = param['window']
    max_length = 24#50

    
    for filename in inputfileinfo:
        o+=1
        gzmark=False
        if filename.endswith('.gz'):
            f = gzip.open(filename)
            gzmark=True
        else:
            f = open(filename)


        if f:
            while 1:
                if gzmark:
                    line1 = f.readline().decode()
                else:
                    line1 = f.readline()
                if not line1:
                    break
                if gzmark:
                    line2 = f.readline().decode().strip()
                    line3 = f.readline().decode()
                    line4 = f.readline().decode().strip()
                else:
                    line2 = f.readline().strip()
                    line3 = f.readline()
                    line4 = f.readline().strip()
                line1 = line1.strip().split()
                read_name = line1[0][1:]
                pair_file_rank = 0
                if len(line1)>1:
                    if line1[1][0]>='1' and line1[1][0]<='2':
                        pair_file_rank = int(line1[1][0])
                if '/' in read_name or '.' in read_name:
                    split_pos=0
                    if '/' in read_name:
                        split_pos = read_name.find('/')
                    else:
                        split_pos = read_name.find('.')
                    if input_file_num>1:
                        pairnum = read_name[split_pos+1]
                        if pairnum>='1' and pairnum<='2':
                            pair_file_rank = int(pairnum)
                        if pair_file_rank==0:
                            print('Unable to recognize given pairend fastq file format!!')
                            pair_file_rank = -1
                    read_name = read_name[:split_pos]

                string_mark = o
                if (read_name in set_sam):
                    if set_sam[read_name]==0 or set_sam[read_name]==3: continue
                    if set_sam[read_name]==pair_file_rank: continue
                    if pair_file_rank<0: continue

                if input_file_num>1: read_name+='_'+str(pair_file_rank)
                #Maybe the mate search method is buggy. Cuz there are different structures of reads name generated by different sequencing machine.
                # Fixed in issue 1 in github.
                #fastqlines[i] = line1[0]+'_'+line1[1][0]+' '+line1[1]
                UnmappedReads[read_name]=[line2,line4]
                if len(UnmappedReads)>1000000:
                    writefiles(UnmappedReads,step,length_bin,max_length,outputname)
                    UnmappedReads={}
                    # if len(pfn)>len(Part_Fastq_Filename):
                    #     Part_Fastq_Filename=pfn
                    

#We've got a dictionary named UnmappedReads = {readsname:[line1,line2,line3,line4]}

#Change cut funtion into cut(setp,length_bin,string,fileorder), return Available(T/F), reads_fraction
    if len(UnmappedReads)>0:
        writefiles(UnmappedReads,step,length_bin,max_length,outputname)
        # if len(pfn)>len(Part_Fastq_Filename):
        #     Part_Fastq_Filename=pfn
    print('finish')
    f.close()
    del UnmappedReads
    
    #We've got the splited fastq file, filename is stored in Part_Fastq_Filename
    # for i in range(len(Part_Fastq_Filename)):
    command = 'bsmap -a '+unmapped_file+' -z '+str(phred)+' -d '+refpath+'  -o '+unmapped_file+'.bam -S 123 -n 1 -r 0 -R 1>>bsmap_log 2>>bsmap_err'
    Batch_try = Pshell(command)
    Batch_try.process()
    command = 'samtools view '+unmapped_file+'.bam >'+unmapped_file+'.sam'
    Batch_try = Pshell(command)
    Batch_try.process()

   #run bsmap and get bam files named as Part_Fastq_Filename[i].bam
    #import combine to generate the finalfastq

    args={'step':param['step'],
          'binsize':param['window'],
          'filter':filter,
          'outputname':outputname,
          'originalfile':inputfileinfo,
          'mapfilenumber':10
          #'finish':1
         }
    mapreduce_names = reads_map(unmapped_file,args)
    reads_reduce(mapreduce_names,args)
    splitlogname = 'BAM_FILE/'+outputname+'_split_log.record'

    command = 'bsmap -a '+outputname+'_finalfastq.fastq -d '+refpath+' -z '+str(phred)+' -o '+outputname+'_split.bam -S 123 -n 1 -r 0 1>>log 2>> '+splitlogname
    Bam = Pshell(command)
    Bam.process()
    splitfilename = outputname+'_split.bam'

    command='samtools sort -f -@ 4 '+splitfilename+' '+splitfilename+'.sorted.bam'
    filter = Pshell(command)
    filter.process()
    if not given_bam_file:
        command='samtools sort -f -@ 4 '+outputname+'.bam'+' '+outputname+'.sort.bam'
    else:
        command='samtools sort -f -@ 4 '+given_bam_file+' '+outputname+'.sort.bam'
    filter.change(command)
    filter.process()
    command='mv '+outputname+'.sort.bam '+outputname+'.bam'
    filter.change(command)
    filter.process()
    command='mv '+splitfilename+'.sorted.bam '+splitfilename
    filter.change(command)
    filter.process()
    m=Pshell('samtools merge BAM_FILE/'+outputname+'_combine.bam '+outputname+'.bam '+splitfilename)
    m.process()
    command='mv '+outputname+'.bam BAM_FILE/'
    filter.change(command)
    filter.process()
    command='mv '+splitfilename+' BAM_FILE/'
    filter.change(command)
    filter.process()
    return 'BAM_FILE/'+outputname+'_combine.bam',originallogname,splitlogname,[unmapped_file,outputname+'_finalfastq.fastq',outputname+'.sam']
    print("Merge done!\nCreated final bam file called "+outputname+'_combine.bam')

def clipmode(name,param, given_bam_file):
    '''
    When we get the mapping result, we should report 
    mapping ratio, mapped reads number, length distribution,
    original mapping ratio, original mapped reads number,
    new generated splitted reads number, new generated splitted reads length
    '''
    newn,originallog,splitlog,cleanname=clip_process(name,param, given_bam_file)

     if (not 'cleanmode' in param) or param['cleanmode']:
         #Set a clean mode and full mode for clipping mode
         cleanupmess(cleanname)

    return newn,[originallog,splitlog]

def cleanupmess(name):
    unmapped_file,n1,n2 = name
    outputname = n1[:n1.rfind('_')]
    print(outputname)
    os.system('rm '+n1)
    os.system('rm '+n2)
    os.system('rm '+outputname+'.bam.bai')
    os.system('rm '+outputname+'_split.bam.bai')
    os.system('rm '+unmapped_file)
    os.system('rm '+unmapped_file+'.bam')
    os.system('rm '+unmapped_file+'.sam')
    os.system('rm '+unmapped_file+'.bam.bai')
    for i in range(10):
        os.system('rm '+outputname+'_'+str(i)+'.mapreduce')
    


if __name__=="__main__":
    
    cleanupmess(['6P1_notrim_val_1_val_1_test.unmapped.fastq','6P1_notrim_val_1_val_1_test_finalfastq.fastq','6P1_notrim_val_1_val_1_test.sam'])
    #with open("config.txt") as f:
    #    lines = f.readlines()
    #import multiprocessing
    #pool = multiprocessing.Pool(processes=2)
    #for l in lines:
    #    #pool.apply_async(clip_process,(l,))
    #    clip_process(l.strip().split(),{'ref':'/data/dsun/ref/humanigenome/hg19.fa','step':5,'window':30,'cleanmode':True}) #pass file name to clip_process
    #pool.close()
    #pool.join()
    #p = {'ref':'/data/dsun/ref/humanigenome/hg19.fa','step':5,'window':30,'cleanmode':True}
    #clipmode(['Trim/head.fq'],p)

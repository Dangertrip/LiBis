from .mapreduce import *
import multiprocessing
#from Comb_fastq import combine
import pysam
from .utils import *
import sys
import gzip
import os
from .pairend import *

max_file_length = 500 # 24

_index = []
for i in range(max_file_length):
    _index.append(str(i))

reads_cut_buff = {}

def cut(step,length_bin,link,i):
    start = step*i
    end = length_bin + start
    if end > len(link):
        end = len(link)
    if start>=end or end-start<=length_bin-step:
        return False,''
    return True,link[start:end]

def cut_serial(step, length_bin, read_len, max_file_length):
    #read_len = len(link)
    cutted_list = []
    for i in range(max_file_length):
        start = step*i
        end = length_bin + start
        if end > read_len:
            end = read_len
        if start>=end or end-start<=length_bin-step:
            break

        cutted_list.append([start,end])
    return cutted_list


def writefiles(UnmappedReads,step,length_bin,max_file_length,outputname):
    unmapped_file = outputname + '.unmapped.fastq.gz'
    filecontent = []
    for readsname in UnmappedReads:
        link = UnmappedReads[readsname]
        read_length = len(link[0])
        if read_length in reads_cut_buff:
            start_end_list = reads_cut_buff[read_length]
        else:
            start_end_list = cut_serial(step,length_bin, len(link[0]),max_file_length)
            reads_cut_buff[read_length] = start_end_list
        #cutreads = cut_serial(step, length_bin, link[0])
        #cutquality = cut_serial(step, length_bin, link[0])
        #mark,cutreads = cut(step,length_bin,link[0],i)
        #if not mark: continue
        #_,cutquality = cut(step,length_bin,link[1],i)
        for ind, start_end in zip(_index,start_end_list):
            start,end = start_end
            fc = '@'+readsname+'&'+ind+'\n'+link[0][start:end]+'\n+\n'+link[1][start:end]+'\n'
            fc = fc.encode()
            filecontent.append(fc)
        #if len(filecontent)==0: break
        # name = outputname+'.part'+str(i+1)+'.fastq'

        # Part_Fastq_Filename.append(name)
    with gzip.open(unmapped_file,'a') as f:
        f.writelines(filecontent)
    return unmapped_file

def clip_process(inputfileinfo,param,given_bam_file,given_label):
    input_file_num = len(inputfileinfo)
    if input_file_num<=0 or input_file_num>2:
        print("Parameter error in "+l)
        sys.exit()
    refpath = param['ref']
    filter = param['filter_len']
    threads = param['threads']
    purename = inputfileinfo[0][inputfileinfo[0].rfind('/')+1:]
    
    if not given_bam_file:
        outputname = RemoveFastqExtension(purename)
    else:
        purename_bam = given_bam_file[given_bam_file.rfind('/')+1:]
        outputname = RemoveFastqExtension(purename)+'_'+purename_bam[:-4]
        if param['moabs']:
            outputname = purename_bam[:-4]


    # Part_Fastq_Filename=[]
    unmapped_file = outputname + '.unmapped.fastq.gz'
    removeFileIfExist(unmapped_file)
    
    # for i in range(24):
    #     name = outputname+'.part'+str(i+1)+'.fastq'
    #     if os.path.exists(name):
    #         os.system('rm '+name)
    if not given_bam_file:
        output_prefix = outputname
    else:
        output_prefix = given_bam_file[:given_bam_file.rfind('.')]
    
    if param['moabs']:
        store_file_prefix = ''
    else:
        store_file_prefix = 'BAM_FILE/'


    if not given_bam_file:
        originallogname = store_file_prefix+outputname+'_originallog.record'
    else:
        originallogname = ''

    phred=33
    if not given_bam_file:
        if input_file_num==2 :
            command='bsmap -a '+inputfileinfo[0]+' -b '+inputfileinfo[1]+' -z '+str(phred)+' -d '+refpath+' -o '+outputname+'.bam -S 123 -n 1 -r 0 -p ' + threads + ' 1>>LiBis_log 2>>'+originallogname
        else:
            command='bsmap -a '+inputfileinfo[0]+' -z '+str(phred)+' -d '+refpath+'  -o '+outputname+'.bam -S 123 -n 1 -r 0 -p ' + threads + ' 1>>LiBis_log 2>>'+originallogname
        First_try = Pshell(command)
        First_try.process()
    # --------------------------------------------Multiple mapping for multiple matched reads in pairend mapping.----------------------------------
    if 'multiple' in param and param['multiple']:
        ori_bam = outputname+'.bam'
        if given_bam_file:
            ori_bam = given_bam_file
            # pairend_remap(bam_filename,fastq_filename,output_name,refpath, multihits_num=10)
            remove_list = []
            multiple_recovered_bam = pairend_remap(ori_bam, inputfileinfo, outputname, refpath, output_prefix, multihits_num=10)
            remove_list.append(multiple_recovered_bam)
            remove_list.append(multiple_recovered_bam+'.bai')
            multiple_recovered_bam = bam_merge(ori_bam, multiple_recovered_bam)
            for r in remove_list:
                removeFileIfExist(r)
 
    set_sam = {}
    set_sam_pair = set()
    bam_file_name = outputname+'.bam'
    if given_bam_file:
        bam_file_name = given_bam_file
    bam_file = pysam.AlignmentFile(bam_file_name,'rb')
    for line in bam_file:
        qname = line.query_name
        next_name = line.next_reference_name
        next_pos = line.template_length
        if next_name!=None:
            if next_pos>0:
                set_sam_pair.add(qname)
        else:
            if input_file_num==2:
                m1 = (line.flag & 64)
                m2 = (line.flag & 128)
                if m1>0: mate = 1
                elif m2>0: mate = 2
                else: mate = 0
            else:
                mate = 0
            set_sam[qname]=mate
    
    UnmappedReads = {}
    o=0
    step = param['step']
    length_bin = param['window']

    
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
                pair_file_rank = o
                if '/' in read_name:
                    split_pos=0
                    if '/' in read_name:
                        split_pos = read_name.rfind('/')
                    read_name = read_name[:split_pos]

                if read_name in set_sam_pair: continue
                mapped_mate = set_sam.get(read_name,-1)
                if mapped_mate==0 or mapped_mate==pair_file_rank: continue
                #if (read_name in set_sam):
                #    if set_sam[read_name]==0 or set_sam[read_name]==3: continue
                #    if set_sam[read_name]==pair_file_rank: continue
                #    if pair_file_rank<0: continue

                if input_file_num>1: read_name+='_'+str(pair_file_rank)
                # Maybe the mate search method is buggy. Cuz there are different structures of reads name generated by different sequencing machine.
                # Fixed in issue 1 in github.
                UnmappedReads[read_name]=[line2,line4]
                if len(UnmappedReads)>3000000:
                    writefiles(UnmappedReads,step,length_bin,max_file_length,outputname)
                    UnmappedReads={}
                    

    if len(UnmappedReads)>0:
        writefiles(UnmappedReads,step,length_bin,max_file_length,outputname)
    f.close()
    del UnmappedReads
    # sys.exit()    
    # We've got the splited fastq file, filename is stored in Part_Fastq_Filename
    # for i in range(len(Part_Fastq_Filename)):
    command = 'bsmap -a '+unmapped_file+' -d '+refpath+'  -o '+unmapped_file[:-3]+'.bam -S 123 -n 1 -r 0 -R -p ' + threads + ' 1>>LiBis_log 2>>'+store_file_prefix+unmapped_file[:-3]+'_log.txt'
    Batch_try = Pshell(command)
    Batch_try.process()

    # run bsmap and get bam files named as Part_Fastq_Filename[i].bam
    # import combine to generate the finalfastq

    args={'step':param['step'],
          'binsize':param['window'],
          'filter':filter,
          'outputname':outputname,
          'originalfile':inputfileinfo,
          'mapfilenumber':10,
          'report_clip':param['report_clip']
          #'finish':1
         }
    mapreduce_names = reads_map(unmapped_file[:-3],args)
    reads_reduce(mapreduce_names,args)
    splitlogname = store_file_prefix+outputname+'_split_log.record'

    command = 'bsmap -a '+outputname+'_finalfastq.fastq.gz -d '+refpath+' -o '+outputname+'_split.bam -S 123 -n 1 -r 0 -p ' + threads + ' 1>>LiBis_log 2>> '+splitlogname
    Bam = Pshell(command)
    Bam.process()
    splitfilename = outputname+'_split.bam'
    #------------------------------------------------Pairend filter for LiBis rescued reads----------------------------------------------------------
    if 'pairend' in param and param['pairend']:
        ori_bam = outputname+'.bam'
        if given_bam_file:
            ori_bam = given_bam_file
        PEfilter_bam = libis_filter(ori_bam, splitfilename)
        command = 'mv '+PEfilter_bam+' '+splitfilename
        renamer = Pshell(command)
        renamer.process()
        #libis_filter(ori_bam, split_bam)
    #------------------------------------------------------------------------------------------------------------------------------------------------
    
    '''
    command='samtools sort -f -@ '+threads+' '+splitfilename+' '+splitfilename+'.sorted.bam'
    filter = Pshell(command)
    filter.process()
    if not given_bam_file:
        command='samtools sort -f -@ '+threads+' '+outputname+'.bam'+' '+outputname+'.sort.bam'
    else:
        command='samtools sort -f -@ '+threads+' '+given_bam_file+' '+outputname+'.sort.bam'
    filter.change(command)
    filter.process()
    command='mv '+outputname+'.sort.bam '+outputname+'.bam'
    filter.change(command)
    filter.process()
    command='mv '+splitfilename+'.sorted.bam '+splitfilename
    filter.change(command)
    filter.process()
    '''
    if not param["moabs"]:
        filter = Pshell("")
        if given_bam_file:
            command='cp '+given_bam_file+' '+outputname+'.bam'
            filter.change(command)
            filter.process()
        m=Pshell('samtools merge '+store_file_prefix+outputname+'_combine.bam '+outputname+'.bam '+splitfilename)
        m.process()
        command='mv '+outputname+'.bam BAM_FILE/'
        filter.change(command)
        filter.process()
        command='mv '+splitfilename+' BAM_FILE/'
        filter.change(command)
        filter.process()

    return 'BAM_FILE/'+outputname+'_combine.bam',originallogname,splitlogname,[unmapped_file,outputname+'_finalfastq.fastq.gz',outputname+'.sam']
    print("Merge done!")#\nCreated final bam file called "+outputname+'_combine.bam')

def clipmode(name,param, given_bam_file,given_label):
    '''
    When we get the mapping result, we should report 
    mapping ratio, mapped reads number, length distribution,
    original mapping ratio, original mapped reads number,
    new generated splitted reads number, new generated splitted reads length
    '''
    newn,originallog,splitlog,cleanname=clip_process(name,param, given_bam_file,given_label)

    if (not 'cleanmode' in param) or param['cleanmode']:
        #Set a clean mode and full mode for clipping mode
        cleanupmess(cleanname)

    return newn,[originallog,splitlog]


def cleanupmess(name):
    unmapped_file,n1,n2 = name
    outputname = n1[:n1.rfind('_')]
    removeFileIfExist(n1)
    removeFileIfExist(n2)
    removeFileIfExist(outputname+'.bam.bai')
    removeFileIfExist(outputname+'_split.bam.bai')
    removeFileIfExist(unmapped_file)
    removeFileIfExist(unmapped_file[:-3]+'.bam')
    removeFileIfExist(unmapped_file[:-3]+'.sam')
    removeFileIfExist(unmapped_file[:-3]+'.bam.bai')
    for i in range(10):
        removeFileIfExist(outputname+'_'+str(i)+'.mapreduce')
    


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

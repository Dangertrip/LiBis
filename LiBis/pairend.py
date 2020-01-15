import gzip
import pysam
from .utils import Pshell
import os

def pair_valid(chr1,chr2,strand1,strand2,pos1,pos2):
    if chr1!=chr2:
        return False
    if strand1=='++' and strand2!='+-':
        return False
    if strand1=='-+' and strand2!='--':
        return False
    if strand1=='+-' and strand2!='++':
        return False
    if strand1=='--' and strand2!='-+':
        return False
    if abs(int(pos1)-int(pos2))<=1000:
        return True
    return False
 

def writefiles(UnmappedReads, outputname):

    filecontent = []
    for readsname in UnmappedReads:
        link = UnmappedReads[readsname]
        filecontent.append('@'+readsname+'\n'+link[0]+'\n+\n'+link[1]+'\n')
    
    with open(outputname,'a') as f:
        f.writelines(filecontent)
    


def pairend_remap(bam_filename,fastq_filename,output_name,refpath,output_prefix, multihits_num=10):
    '''
        Args:
            bam_filename: string, name of combined bam file
            fastq_filename: string[], list of names of fastq files
        Output:
            1. 
    '''

    #with open(outputname+".sam") as sam:
    #second column in sam file: 64, mate 1; 128, mate 2;
    #    samlines = sam.readlines()
    output_name_fq = output_name + '_multi.fq'
    output_name_bam = output_name + '_multi.bam'
    bamfile = pysam.AlignmentFile(bam_filename, "rb")
    set_sam = {}
    fq_file_num = len(fastq_filename)
    for line in bamfile:
        temp = line.to_dict()
        name = temp['name']
        if temp['next_ref_name']!='*':
            continue
        clip_from = -1
        if '_' in name:
            pos = name.find('_')
            if name[pos+3]=='0':
                clip_from = 1
            else:
                clip_from = 2
            name = name[:pos]
        flag = int(temp['flag'])
        m1 = (flag & 64)
        m2 = (flag & 128)
        if m1>0: mate = 1
        elif m2>0: mate = 2
        else: mate = 0
        if clip_from>0 and fq_file_num>1:
            mate = clip_from
        if name in set_sam: set_sam[name]=3
        else: set_sam[name]=mate 
    bamfile.close()
    
    UnmappedReads = {}
    o=0

    for filename in fastq_filename:
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
                    split_pos = read_name.rfind('/')
                    read_name = read_name[:split_pos]

                if (read_name in set_sam):
                    if set_sam[read_name]==0 or set_sam[read_name]==3: continue
                    if set_sam[read_name]==pair_file_rank: continue
                    if pair_file_rank<0: continue
                else:
                    continue

                #if fq_file_num>1: read_name+='_'+str(pair_file_rank)
                UnmappedReads[read_name]=[line2,line4]
                if len(UnmappedReads)>1000000:
                    writefiles(UnmappedReads,output_name_fq)
                    UnmappedReads={}
        f.close()            

    if len(UnmappedReads)>0:
        writefiles(UnmappedReads,output_name_fq)
    #print('finish')
    del UnmappedReads

    command='bsmap -a '+output_name_fq+' -d '+refpath+'  -o '+output_name_bam+' -w '+str(multihits_num)+' -S 123 -n 1 -r 2 1>>temp_log 2>>'+output_prefix+'_multi_log.txt'
    multiple_mapping = Pshell(command)
    multiple_mapping.process()
    if os.path.exists(output_name_fq):
        os.remove(output_name_fq)
    return output_name_bam

def bam_merge(ori_bam, multi_bam):
    '''
        Args:
            ori_bam: string, name of the original *_combined.bam
            multi_bam: string, name of the new multiple mapping bam
    '''
    
    samfile = pysam.AlignmentFile(multi_bam, "rb") 
    multi_dict = {}
    for line in samfile:
        temp = line.to_dict()
        chr = temp['ref_name']
        pos = temp['ref_pos']
        strand=''
        for tag in temp['tags']:
            if 'ZS' in tag:
                strand = tag[-2:]
        name = temp['name'][:-2]
        if name in multi_dict:
            multi_dict[name].append([chr,pos,strand])
        else:
            multi_dict[name] = [[chr,pos,strand]]

    samfile.close()
    potential_marker = {}
    samfile = pysam.AlignmentFile(ori_bam, "rb")
    for line in samfile:
        temp = line.to_dict()
        if temp['next_ref_name']!='*':
            continue
        name = temp['name']
        if name not in multi_dict:
            continue
        chr = temp['ref_name']
        pos = int(temp['ref_pos'])
        strand=''
        for tag in temp['tags']:
            if 'ZS' in tag:
                strand = tag[-2:]
        
        candidates = multi_dict[name]
        for c in candidates:
            if chr==c[0]:
                if strand=='++' and c[2]!='+-':
                    continue
                if strand=='-+' and c[2]!='--':
                    continue
                if strand=='+-' and c[2]!='++':
                    continue
                if strand=='--' and c[2]!='-+':
                    continue
                if abs(pos-int(c[1]))<=1000:
                    if name in potential_marker:
                        potential_marker[name].append(c)
                    else:
                        potential_marker[name] = [c]
                    
    samfile.close()
        
    multifile = pysam.AlignmentFile(multi_bam, "rb") 
    output_file = pysam.AlignmentFile(multi_bam+'.select.bam','wb',template=multifile)
    for line in multifile:
        temp = line.to_dict()
        name = temp['name']
        #mate = int(name[-1])
        #name = name[:-2]
        if name not in potential_marker:
            continue
        mate = int(name[-1])
        if mate==1:
            flag = 64
        else:
            flag = 128
        #line.flag |= 
        chr = temp['ref_name']
        pos = temp['ref_pos']
        for c in potential_marker:
            if chr==c[0] and pos==c[1]:
                output_file.write(line)
    multifile.close()
    output_file.close()
    return multi_bam+'.select.bam'

def libis_filter(ori_bam, split_bam):
    '''
        Args:
            ori_bam: string, name of the original *_combined.bam
            multi_bam: string, name of the new multiple mapping bam
    '''
    
    samfile = pysam.AlignmentFile(ori_bam, "rb") 
    ori_dict = {}
    for line in samfile:
        temp = line.to_dict()
        if temp['next_ref_name']!='*':
            continue
        chr = temp['ref_name']
        pos = temp['ref_pos']
        strand=''
        for tag in temp['tags']:
            if 'ZS' in tag:
                strand = tag[-2:]
        name = temp['name'][:-2]
        ori_dict[name] = [chr,pos,strand]

    samfile.close()
    split_dict = {}
    candidates = {}
    samfile = pysam.AlignmentFile(split_bam, "rb")
    output_file = pysam.AlignmentFile(split_bam+'.filter.bam','wb',template=samfile)
    number = 0
    for line in samfile:
        temp = line.to_dict()
        name = temp['name']
        #print(name)
        _pos = name.find('_')
        #print(pos)
        #print(name[:pos])
        clip_from = -1
        if name[_pos+3]=='1':
            clip_from = 1
        else:
            clip_from = 2
        chr = temp['ref_name']
        pos = temp['ref_pos']
        strand = ''
        for tag in temp['tags']:
            if 'ZS' in tag:
                strand = tag[-2:]
        name = name[:_pos]
        if name in ori_dict:
            chr1, pos1, strand1 = ori_dict[name]
            if pair_valid(chr,chr1,strand,strand1,pos,pos1):
                output_file.write(line)
                #number += 1
        else:
            if name not in split_dict:
                split_dict[name] = [[clip_from, chr, pos, strand]]
            else:
                split_dict[name].append([clip_from, chr, pos, strand])
    samfile.close()
    for name, c in split_dict.items():
        if len(c)<2: continue
        clip_from = 0
        for cc in c:
            clip_from = clip_from | cc[0]
        if clip_from!=3: continue
        #candidate.append(name)
        mark = False
        cand = set()
        for i in range(len(c)):
            for j in range(i+1,len(c)):
                if c[i][0] | c[j][0]==3 and pair_valid(c[i][1],c[j][1],c[i][3],c[j][3],c[i][2],c[j][2]):
                    cand.add((c[j][1],c[j][2],c[j][3]))
                    cand.add((c[i][1],c[i][2],c[i][3]))
                    mark = True
        if mark:
            candidates[name] = list(cand)
                    #number += 1
    #print(number)
        
    splitfile = pysam.AlignmentFile(split_bam, "rb") 
    for line in splitfile:
        temp = line.to_dict()
        name = temp['name']
        #mate = int(name[-1])
        #name = name[:-2]
        _pos = name.find('_')
        name = name[:_pos]
        if name not in candidates:
            continue
        chr = temp['ref_name']
        pos = temp['ref_pos']
        strand=''
        for tag in temp['tags']:
            if 'ZS' in tag:
                strand = tag[-2:]
        for c in candidates[name]:
            if chr==c[0] and pos==c[1] and strand==c[2]:
                output_file.write(line)
    splitfile.close()
    output_file.close()
    return split_bam+'.filter.bam'
    

if __name__=="__main__":
    import sys
    refpath = '/fdata/scratch/lijin.bio/resource/genome/hg38/hg38.fa'
    #bam = sys.argv[1]#'H1_1A.bam'
    #fq = sys.argv[2:4]#['H1_1A_CKDL190142710-1a_HMFVMDSXX_L2_1_val_1.fq.gz','H1_1A_CKDL190142710-1a_HMFVMDSXX_L2_2_val_2.fq.gz']
    #pairend_remap(bam,fq,bam[:bam.find('.')]+'_multi',refpath)   
    #bam_merge('H1_1D.bam','H1_1D_multi.bam')
    ori_bam, split_bam = sys.argv[1], sys.argv[2]
    libis_filter(ori_bam, split_bam)

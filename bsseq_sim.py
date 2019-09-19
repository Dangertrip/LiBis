import numpy as np
import math
import random
import sys

#==============================Variable setting===================================
readsname = sys.argv[1]
readsnum= int(sys.argv[2])
conversion_ratio=1
methylation_ratio=0.75
read_length=100
chrom_set, chrom_len, chrom_name, genome_len = genome_loader('/data/dsun/ref/humanigenome/hg19.fa')
methy_dict = bed_reader('./hESC.bed')
'''
@Variables:
chrom_set: dictionary, {chromsome_name:chromsome_genome}
chrom_len: int[], chromsome length ordered by chromsome name in fastq file
chrom_name: string[], chromsome names
genome_len: int, length of whole genome
methy_dict: dictionary, {'chr1_10469':0.1(methylation ratio)}
'''
#================================================================================

def bed_reader(filename):
    '''
        read bed file as fixed methylation ratio
    '''
    dict = {}
    with open(filename) as f:
        for line in f:
            line_content = line.strip().split()
            dict[line_content[0]+'_'+line_content[2]] = float(line_content[4])
    return dict

def genome_loader(filename):
    chr = ''
    chrom_set = {}
    seq = ''
    genome_len = 0
    chrom_len = []
    chrom_name = []
    with open(filename) as f:
        for line in f:
            line_content = line.strip()
            if line_content[0]=='>':
                if chr!='':
                    chrom_set[chr] = seq
                    chrom_len.append(len(seq))
                    genome_len += chrom_len[-1]
                    seq=''
                chr = line_content[1:]
                chrom_name.append(chr)
                continue
            seq += line_content
        chrom_set[chr] = seq
        chrom_len.append(len(seq))
        genome_len += chrom_len[-1]
    return chrom_set, chrom_len, chrom_name, genome_len


def fake_read(length):

    base = ['A','T','C','G']
    seq = ''
    for i in range(length):
        pos = random.randint(0,3)
        seq += base[pos]
    return seq


def reverse(read):
    dic={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    r=''
    for rr in read:
        r=dic[rr.upper()]+r
    return r


def bisulfite(chrom, start, read):
    '''
    return bisulfited seq.
    Ignore the reads if get ''
    '''
    r=''
    l = len(read)
    for i in range(1,l-1): #waste 2 bp
        base=read[i].upper()
        if base=='C':
            c = random.random()
            if c<conversion_ratio:
                if read[i+1].upper()=='G':
                    pos = chrom+'_'+str(start+i)
                    if pos not in methy_dict:
                        return ''
                    m = random.random()
                    if m>methy_dict[pos]:
                        base='T'
                else:
                    base='T'
        r=r+base
    return r



'''
@Variables:
chrom_set: dictionary, {chromsome_name:chromsome_genome}
chrom_len: int[], chromsome length ordered by chromsome name in fastq file
chrom_name: string[], chromsome names
genome_len: int, length of whole genome
methy_dict: dictionary, {'chr1_10469':0.1(methylation ratio)}
'''

def random_head_tail():
    finalbed=[]
    fake_length = 20
    real_read_length = read_length - fake_length
    for i in range(readsnum):
        pos = random.randint(0,genome_len)
        chr=0
        real_read_length = read_length - 20
        while pos>chrom_len[chr]:
            pos-=chrom_len[chr]
            chr+=1
        start=pos-1
        end=pos+real_read_length+1
        if pos<1:continue
        if pos>chrom_len[chr]:continue
        read=chrom_set[chrom_name[chr]][start:end]
        r=read
        a1=random.random()
        a2=random.random()
        if a1>0.5:
            r=reverse(read)#Get reads from +/- strand
        r = bisulfite(conversion_ratio,methylation_ratio,r)
        if a2>0.5:
            r=reverse(r)#PCR +/-
        fake_marker=''
        if fake_length>0:
            head_fake_length = random.randint(0,fake_length)
            tail_fake_length = fake_length - head_fake_length
            head = fake_read(head_fake_length)
            tail = fake_read(tail_fake_length)
            r = head + r + tail
            fake_marker = '_'+str(head_fake_length)+'_'+str(tail_fake_length)
        quality = 'E'*read_length
        print('@'+str(i)+'_'+readsname)
        print(r)
        print('+')
        print(quality)
        finalbed.append(chrom_name[chr]+'\t'+str(start+1)+'\t'+str(end-1)+'\t'+str(i)+'_'+readsname+'\t'+str(f)+'\t'+str(t)+'\n')
    with open(readsname+'_simulation.bed','w') as f:
        f.writelines(finalbed)


if __name__=="__main__":
    random_head_tail()
    
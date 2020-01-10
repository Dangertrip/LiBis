from .utils import *

def BsmapOutputExtractor(filename):
    #aligned=[]
    #unique=[]
    #total=[]
    lines = readf(filename)
    countinputfile=sum(list(map(lambda x: 'Input read file' in x,lines)))
    total=0
    unique=0
    aligned=0
    if countinputfile==2:
        for line in lines:
            if "total read pairs:" in line:
                a = line.strip().split()
                total=int(line[line.rfind('total read pairs:')+len('total read pairs:'):].split()[0])                    
                continue
            if "aligned pairs:" in line:
                a = line.strip().split()
                aligned = aligned+int(a[2])*2
                unique = unique + int(a[6])*2
                continue
            if "unpaired read #1:" in line or "unpaired read #2:" in line:
                a = line.strip().split()
                aligned = aligned+int(a[3])
                unique = unique+int(a[7])
    else:
        for line in lines:
            if "total reads:" in line:
                a = line.strip().split()
                # print(a)
                total=int(line[line.rfind('total reads:')+len('total reads:'):].split()[0])
                continue
            if "aligned reads:" in line:
                a = line.strip().split()
                aligned = aligned+int(a[2])
                unique = unique+int(a[6])
    return [total,aligned,unique]


def McallOutputExtractor(filename):
    lines = readf(filename)


def TrimOutputExtractor(filename):
    '''
    Tested
    '''
    ratio=[]
    for nn in filename:
        temp=[]
        for name in nn: 
            lines = readf(name)
            #Total written (filtered):  1,085,317,679 bp (94.2%)
            for line in lines:
                temp = line.strip().split()
                if len(temp)<5: continue
                r=None
                if temp[0]=='Total' and temp[1]=='written' and temp[2]=='(filtered):':
                    r=float(temp[5].strip()[1:-2])
                    break
                temp.append(r)
            ratio.append(temp)
    return ratio

'''
Result from BsmapOutputExtractor:
    total reads, mapped reads, uniquely mapped reads
Result from BsmapResult:
    [[total reads,mapped reads,uniquely mapped reads, clipped reads, unique clipped reads,
    all mapped reads, all uniquely mapped reads, mapping ratio, uniquely mapping ratio],...]
'''
def BsmapResult(filenames,clip,offered_bamfile):
    if clip:
        '''
        We will get two record files in this mode
        '''
        result=[]
        # print(filenames)
        # print(offered_bamfile)
        for (name, processed_bam) in zip(filenames,offered_bamfile):
            ori,spl = name
            ori_info = BsmapOutputExtractor(ori) if processed_bam=='' else [1,1,1]
            spl_info = BsmapOutputExtractor(spl)
            total_reads,mapped_reads,uniquely_mapped_reads=ori_info
            _,clipped_reads,uniquely_clipped_reads=spl_info
            result.append([total_reads,
                mapped_reads,
                uniquely_mapped_reads,
                clipped_reads,
                uniquely_clipped_reads,
                mapped_reads+clipped_reads,
                uniquely_mapped_reads+uniquely_clipped_reads,
                (mapped_reads+clipped_reads)/float(total_reads),
                (uniquely_mapped_reads+uniquely_clipped_reads)/float(total_reads)
            ])
            print(result)
            if processed_bam!='':
                for i in [0,1,4,5,6,7,8]:
                    result[-1][i]='Not Available'
    else:
        '''
        Only one record file need to be processed.
        '''
        if processed_bam!='':
            result.append(['Not Available']*9)
            return result
        result=[]
        for name in filenames:
            ori = name
            ori_info = BsmapOutputExtractor(ori)
            print(ori_info)
            total_reads,mapped_reads,uniquely_mapped_reads=ori_info
            _,clipped_reads,uniquely_clipped_reads=[0,0,0]
            result.append([total_reads,mapped_reads,uniquely_mapped_reads,clipped_reads,
                uniquely_clipped_reads,
                mapped_reads+clipped_reads,uniquely_mapped_reads+uniquely_clipped_reads,
                (mapped_reads+clipped_reads)/float(total_reads),
                (uniquely_mapped_reads+uniquely_clipped_reads)/float(total_reads)])
    return result

if __name__=="__main__":
    r=BsmapResult(['BAM_FILE/temp'],False)
    print(r)

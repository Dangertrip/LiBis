from utils import overlap,Pshell
import utils
import sys
def GetFastqList(joined_reads,Part_Fastq_Filename,step,length_bin):

    nameset={}
    contentset={}
    filenum = len(Part_Fastq_Filename)
    #Generate a dictionary which contains readsname, start file order and extend fraction number
    for name in joined_reads:
        reads_list = joined_reads[name]
        if len(reads_list)==0: continue
        #n = reads_list[0].fqname
        n = name
        nameset[n]=[[read.order,read.sum] for read in reads_list]
        contentset[n]=[['',''] for i in range(len(nameset[n]))]#read_content,read_quality
        #contentset[n]=['','']


    file_order=0
    for name in Part_Fastq_Filename:
        with open(name) as f:
            lines = f.readlines()
        for i in range(0,len(lines),4):
            fqname = lines[i].strip().split()[0][1:]
            if not fqname in nameset: continue
            read = lines[i+1].strip()
            quality = lines[i+3].strip()
            fraction_num=0
            for order,sum in nameset[fqname]:
                #print(order,sum)
                if file_order==order:
                    contentset[fqname][fraction_num][0]=read
                    contentset[fqname][fraction_num][1]=quality
                if file_order>order and file_order<=order+sum:
                    add_length = ((len(read)-1)%step)+1
                    contentset[fqname][fraction_num][0]+=read[-1*add_length:]
                    contentset[fqname][fraction_num][1]+=quality[-1*add_length:]
                fraction_num+=1
        file_order+=1
    return contentset


def combine(outputname,Part_Fastq_Filename,step,length_bin,filter=40):
    Part_Fastq_Filename=Part_Fastq_Filename
    #print(Part_Fastq_Filename)
    cache_length=3
    result={}
    file_order=0
    for name in Part_Fastq_Filename:
        command = 'samtools view '+name+'.bam >'+name+'.sam'
        SamFileMaker = Pshell(command)
        SamFileMaker.process()
        with open(name+'.sam') as f:
        #print(partsamlines[-1])
            for line in f:
                #print(line)
                s = line.strip().split('\t')
                mismatch = int(s[11][s[11].rfind(':')+1:])
                if mismatch>1: continue  #testmismatch
                if not(s[0] in result) or (len(result[s[0]])==0):
                    result[s[0]]=[utils.reads(s[2],int(s[3]),file_order,mismatch)]
                else:
                    temp = utils.reads(s[2],int(s[3]),file_order,mismatch)
                    #reads from cliped mapped bam
                    join_or_not=False
                    read_length = len(s[9])
                    tail_length = ((read_length-1)%step)+1
                    refseq = s[12][-2-tail_length:-2]
                    readsseq = s[9][-tail_length:]
                    strand = s[13][-2:]
                    for reads in result[s[0]]:#Try to join existing seeds
                        if reads.canjoin(temp,step,read_length,length_bin):
                            mis=0
                            for ppp in range(tail_length):
                                if (refseq[ppp]!=readsseq[ppp]):
                                    #Here ++/+-/-+/-- should be considered. C/T or A/G match should be identified.
                                    if strand[0]=='+':
                                        if (refseq[ppp]=='C' and readsseq=='T'): continue
                                    else:
                                        if (refseq[ppp]=='G' and readsseq=='A'): continue
                                    mis+=1
                            reads.join(temp,mis)
                            join_or_not=True
                            break

                    frac_list=result[s[0]]
                    if not join_or_not: #temp reads haven't join any exist reads
                        frac_list.append(temp) #add temp reads to array as new seed
                        #if file_order>2 and len(frac_list)>=s:
                        #    for i in range(len(frac_list)-1,-1,-1):
                        #        read = frac_list[i]
                        #        if file_order-read.order>2:
                        #            if read.getSum()==0 and read.getMismatch()>1:
                        #                frac_list.pop(i)

                        #print(len(result[s[0]]))
        file_order+=1
    #join done
    #filter results: filter1
    for name in result:
        nonjoin_num=0
        reads_list=result[name]
        for i in range(len(reads_list)-1,-1,-1):
            if reads_list[i].getSum()==0 and reads_list[i].getMismatch()>1:
                reads_list.pop(i)#Remove all reads which have more than 1 mistake and never be joined
    #filter results: filter2
    for name in result:
        reads_list = result[name]
        num = len(reads_list)
        del_mark = [0 for i in range(num)]
        for i in range(num):
            for j in range(i+1,num):
                if overlap(result[name][i],result[name][j],step,length_bin):
                    sss = result[name][i].getSum()-result[name][j].getSum()
                    if sss>0: del_mark[j]=1
                    elif sss<0: del_mark[i]=1
                    else:
                        mis = result[name][i].getMismatch()-result[name][j].getMismatch()
                        if mis>0: del_mark[i]=1
                        else: del_mark[j]=1
        #Only keep the best read which has the most extends and the least mismatches.
        for i in range(num-1,-1,-1):
            if del_mark[i]==1:
                reads_list.pop(i)

    fastq_dic = GetFastqList(result,Part_Fastq_Filename,step,length_bin)
    #print(fastq_dic['SRR1248467.53331_1'])
    #print(result['SRR1248467.53331_H'][0])
    #sys.exit()
    with open(outputname+'_finalfastq.fastq','w') as f:
        for name in fastq_dic:
            num=1
            for read,quality in fastq_dic[name]:
                #read,quality = fastq_dic[name][num]
                if (len(read)<filter): continue
                f.write('@'+name+'_'+str(num)+'\n')
                f.write(read+'\n')
                f.write('+\n')
                f.write(quality+'\n')
                num+=1








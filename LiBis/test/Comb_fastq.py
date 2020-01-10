from utils import overlap,Pshell
import utils
import sys
import os

def GetFastqList(joined_reads,Part_Fastq_Filename,step,length_bin,filter,outputname,originalfile):

    nameset={}
    contentset={}
    filenum = len(Part_Fastq_Filename)
    #Generate a dictionary which contains readsname, start file order and extend fraction number
    for name in joined_reads:
        reads_list = joined_reads[name]
        if len(reads_list)==0: continue
        #n = reads_list[0].fqname
        n = name
        nameset[n]=[[read[2],read[4]] for read in reads_list]
        contentset[n]=[['',''] for i in range(len(nameset[n]))]#read_content,read_quality
        #contentset[n]=['','']

    pos_mark=[{},{}]
    for name in nameset:
        readinfo = nameset[name]
        pos=0
        if name[-2:]=='_2':
            pos=1
        if name[-2:]=='_1' or name[-2:]=='_2':
            name = name[:-2]
        for order,sum in readinfo:
            start = (order)*step
            end = start + step*sum + length_bin
            if end-start<filter: continue
            if name in pos_mark[pos]:
                pos_mark[pos][name].append([start,end])
            else:
                pos_mark[pos][name]=[[start,end]]
    #num=0
    #for n in pos_mark[0]:
    #    print(n)
    #    if n in pos_mark[0]:
    #        print(pos_mark[0][n])
    #    num+=1
    #    if num>10: break
    fileorder=0
    result=[]
    for file in originalfile:
        with open(file) as f:
            while True:
                name = f.readline()
                if not name:
                    break
                reads = f.readline()
                _ = f.readline()
                quality = f.readline()
                fqname = name.strip().split()[0][1:]
                if not fqname in pos_mark[fileorder]: continue
                for i in range(len(pos_mark[fileorder][fqname])):
                    start,end = pos_mark[fileorder][fqname][i]
                    s_name = fqname
                    if len(pos_mark[pos])>1:
                        s_name += '_'+str(i)
                    s_read = reads[start:end]
                    s_qua = quality[start:end]
                    s_final = '@'+s_name+'_'+str(i)+'\n'+s_read+'\n'+'+\n'+s_qua+'\n'
                    result.append(s_final)
                if len(result)>5000000:
                    with open(outputname+'_finalfastq.fastq','a') as ff:     
                        ff.writelines(result)
                        result=[]
        fileorder+=1
    if len(result)>0:
        with open(outputname+'_finalfastq.fastq','a') as ff:
            ff.writelines(result)
                    
'''
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
'''

def combine(outputname,Part_Fastq_Filename,step,length_bin,filter,originalfile):
    Part_Fastq_Filename=Part_Fastq_Filename
    #print(Part_Fastq_Filename)
    cache_length=3
    result={}
    file_order=0
    for name in Part_Fastq_Filename:
        command = 'samtools view '+name+'.bam >'+name+'.sam'
        SamFileMaker = Pshell(command)
        #SamFileMaker.process()
        print(name)
        print(len(result))
        with open(name+'.sam') as f:
        #print(partsamlines[-1])
            for line in f:
                #print(line)
                s = line.strip().split('\t')
                #name = s[0].split(':')
                #if len(name)>6:
                #    s[0] = ':'.join(name[3:])
                mismatch = int(s[11][s[11].rfind(':')+1:])
                if mismatch>1: continue  #testmismatch
                if not(s[0] in result) or (len(result[s[0]])==0):
                    #result[s[0]]=[utils.reads(s[2][3:],int(s[3]),file_order,mismatch)]
                    result[s[0]]=[[s[2][3:],int(s[3]),file_order,mismatch,0]]
                else:
                    temp = [s[2][3:],int(s[3]),file_order,mismatch,0]#utils.reads(s[2][3:],int(s[3]),file_order,mismatch)
                    #reads from cliped mapped bam
                    join_or_not=False
                    read_length = len(s[9])
                    tail_length = ((read_length-1)%step)+1
                    refseq = s[12][-2-tail_length:-2]
                    readsseq = s[9][-tail_length:]
                    strand = s[13][-2:]
                    for reads in result[s[0]]:#Try to join existing seeds
                        #if reads.canjoin(temp,step,read_length,length_bin):
                        if readsjoin(reads,temp,step,read_length,length_bin):
                            mis=0
                            for ppp in range(tail_length):
                                if (refseq[ppp]!=readsseq[ppp]):
                                    #Here ++/+-/-+/-- should be considered. C/T or A/G match should be identified.
                                    if strand[0]=='+':
                                        if (refseq[ppp]=='C' and readsseq=='T'): continue
                                    else:
                                        if (refseq[ppp]=='G' and readsseq=='A'): continue
                                    mis+=1
                            #update:
                            #reads.join(temp,mis)
                            reads[3]+=mis#mismatch
                            reads[4] = temp[2]-reads[2]#sum
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
            if reads_list[i][4]==0 and reads_list[i][3]>1:
                reads_list.pop(i)#Remove all reads which have more than 1 mistake and never be joined
    #filter results: filter2
    for name in result:
        reads_list = result[name]
        num = len(reads_list)
        del_mark = [0 for i in range(num)]
        for i in range(num):
            for j in range(i+1,num):
                if overlap(result[name][i],result[name][j],step,length_bin):
                    sss = result[name][i][4]-result[name][j][4]
                    if sss>0: del_mark[j]=1
                    elif sss<0: del_mark[i]=1
                    else:
                        mis = result[name][i][3]-result[name][j][3]
                        if mis>0: del_mark[i]=1
                        else: del_mark[j]=1
        #Only keep the best read which has the most extends and the least mismatches.
        for i in range(num-1,-1,-1):
            if del_mark[i]==1:
                reads_list.pop(i)

    #fastq_dic = 
    GetFastqList(result,Part_Fastq_Filename,step,length_bin,filter,outputname,originalfile)
    #print(fastq_dic['SRR1248467.53331_1'])
    #print(result['SRR1248467.53331_H'][0])
    #sys.exit()
    #with open(outputname+'_finalfastq.fastq','w') as f:
    #    for name in fastq_dic:
    #        num=1
    #        for read,quality in fastq_dic[name]:
    #            #read,quality = fastq_dic[name][num]
    #            if (len(read)<filter): continue
    #            f.write('@'+name+'_'+str(num)+'\n')
    #            f.write(read+'\n')
    #            f.write('+\n')
    #            f.write(quality+'\n')
    #            num+=1

if __name__=="__main__":
    partfilename = ['6P.part1.fastq','6P.part2.fastq','6P.part3.fastq','6P.part4.fastq']
    '''
    ['6P.part10.fastq','6P.part14.fastq','6P.part18.fastq','6P.part21.fastq',
                    '6P.part2.fastq','6P.part6.fastq','6P.part11.fastq','6P.part15.fastq',
                    '6P.part19.fastq','6P.part22.fastq','6P.part3.fastq','6P.part7.fastq',
                    '6P.part12.fastq','6P.part16.fastq','6P.part1.fastq','6P.part23.fastq',
                    '6P.part4.fastq','6P.part8.fastq','6P.part13.fastq','6P.part17.fastq',
                    '6P.part20.fastq','6P.part24.fastq','6P.part5.fastq','6P.part9.fastq']
    '''
    combine('6P',partfilename,5,30,40,['6P_R1_trimmed.fq','6P_R2_trimmed.fq'])






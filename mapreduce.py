from utils import *
import os
import gzip

def reads_map(partfilelist,args):
    mapfilenum = args['mapfilenumber']
    step = args['step']

    file_order=0
    mapreduce_file=[]
    rootfile = partfilelist[0][:partfilelist[0].find('.')]
    dic={}
    for i in range(mapfilenum):
        mapreduce_file.append(rootfile+'_'+str(i)+'.mapreduce')
        if os.path.exists(mapreduce_file[-1]) and (not 'finish' in args):
            os.remove(mapreduce_file[-1])
            print('Delete '+mapreduce_file[-1])
        dic[i]=[]
    if 'finish' in args:
        return mapreduce_file

    for file in partfilelist:
        print(file)
        count=0
        with open(file+'.sam') as f:
            for line in f:
                s = line.strip().split('\t')
                mismatch = int(s[11][s[11].rfind(':')+1:])
                if mismatch>1: continue
                read_length = len(s[9])
                tail_length = ((read_length-1)%step)+1
                refseq = s[12][-2-tail_length:-2]
                readsseq = s[9][-tail_length:]
                strand = s[13][-2:]
                mis=0
                for base in range(tail_length):
                    if (refseq[base]!=readsseq[base]):
                        if strand[0]=='+':
                            if (refseq[base]=='C' and readsseq[base]=='T'): continue
                        else:
                            if (refseq[base]=='G' and readsseq[base]=='A'): continue
                        mis+=1
                tail_mismatch = mis

                hashnum = abs(hash(s[0])) % mapfilenum
                dic[hashnum].append([s[0],s[2][3:],s[3],str(file_order),str(mismatch),str(tail_mismatch),str(read_length)])
                count+=1
                if count>5000000:
                    for i in range(mapfilenum):
                        arr = dic[i]
                        arr = list(map(lambda x:'\t'.join(x)+'\n',arr))
                        with open(mapreduce_file[i],'a') as ff:
                            ff.writelines(arr)
                        dic[i]=[]
                    count=0
        file_order+=1
    if count>0:
        for i in range(mapfilenum):
            arr = dic[i]
            arr = list(map(lambda x:'\t'.join(x)+'\n',arr))
            with open(mapreduce_file[i],'a') as ff:
                ff.writelines(arr)
            dic[i]=[]
        count=0
    return mapreduce_file


def reads_combine(filename,args):
    step = args['step']
    length_bin = args['binsize']
    filter = args['filter']
    outputname = args['outputname']
    originalfile = args['originalfile']
    result={}
    with open(filename) as f:
        for line in f:
            arr = line.strip().split()
            name,chr,startpos,fileorder,mismatch,tail_mismatch,read_length = arr
            startpos = int(startpos)
            fileorder = int(fileorder)
            mismatch = int(mismatch)
            tail_mismatch = int(tail_mismatch)
            read_length = int(read_length)
            if (not name in result) or (len(result[name])==0):
                result[name]=[[chr,startpos,fileorder,mismatch,0]]
            else:
                temp = [chr,startpos,fileorder,mismatch,0]
                #COPIED CODE; NEED TO BE MODIFIED
                #reads from cliped mapped bam
                join_or_not=False
                for reads in result[name]:
                    if reads[3]+tail_mismatch<=1 and readsjoin(reads,temp,step,read_length,length_bin):
                        reads[3]+=tail_mismatch
                        reads[4]=temp[2]-reads[2]
                        join_or_not=True
                        break

                frac_list=result[name]
                if not join_or_not:
                    frac_list.append(temp)
            #print(len(result))
    #Delete short fragments
    print(len(result))
    del_name=[]
    for name in result:
        nonjoin_num=0
        reads_list=result[name]
        for i in range(len(reads_list)-1,-1,-1):
            if reads_list[i][4]<=1:# and reads_list[i][3]>1:
                reads_list.pop(i)
        if len(reads_list)==0:
            del_name.append(name)
    for name in del_name:
        del result[name]
    print(len(result))
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
    return result
    #GetFastqList(result,step,length_bin,filter,outputname,originalfile)

    

def reads_reduce(mapreduce_file,args):
    step = args['step']
    length_bin = args['binsize']
    filter = args['filter']
    outputname = args['outputname']
    originalfile = args['originalfile']
    mapfilenum = args['mapfilenumber']

    totalresult={}
    for i in range(mapfilenum):
        print(str(i)+' start!')
        result=reads_combine(mapreduce_file[i],args)
        totalresult.update(result)
        if len(totalresult)>5000000 or i==mapfilenum-1:
            GetFastqList(totalresult,step,length_bin,filter,outputname,originalfile)
            totalresult={}
        #GetFastqList(result,step,length_bin,filter,outputname,originalfile)
        #totalresult.update(result)
        #print(str(i)+' finished! length='+str(len(totalresult)))
    #GetFastqList(totalresult,step,length_bin,filter,outputname,originalfile)


def GetFastqList(joined_reads,step,length_bin,filter,outputname,originalfile):
    #print(joined_reads)
    nameset={}
    #Generate a dictionary which contains readsname, start file order and extend fraction number
    for name in joined_reads:
        reads_list = joined_reads[name]
        if len(reads_list)==0: continue
        n = name
        nameset[n]=[[read[2],read[4]] for read in reads_list]
        #contentset[n]=[['',''] for i in range(len(nameset[n]))]#read_content,read_quality
    print(len(nameset))
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
    print(len(pos_mark[0]),len(pos_mark[1]))
    del nameset
    #num=0
    #for n in pos_mark[0]:
    #    print(n)
    #    if n in pos_mark[0]:
    #        print(pos_mark[0][n])
    #    num+=1
    #    if num>10: break
    fileorder=0
    #print(pos_mark)
    result=[]
    for file in originalfile:
        gzmark=False
        if file.endswith('gz'):
            gzmark=True
        if not gzmark:
            f = open(file)
        else:
            f = gzip.open(file)

        while True:
            name = f.readline()
            if not name:
                break
            reads = f.readline()
            _ = f.readline()
            quality = f.readline()
            if gzmark:
                name = name.decode()
                reads = reads.decode()
                quality = quality.decode()

            reads = reads.strip()
            quality = quality.strip()
            fqname = name.strip().split()[0][1:]
            if not fqname in pos_mark[fileorder]: continue
            for i in range(len(pos_mark[fileorder][fqname])):
                start,end = pos_mark[fileorder][fqname][i]
                s_name = fqname
                if len(pos_mark[pos])>1:
                    s_name += '_'+str(i)
                s_read = reads[start:end]
                s_qua = quality[start:end]
                s_final = '@'+s_name+'_'+str(fileorder)+'\n'+s_read+'\n'+'+\n'+s_qua+'\n'
                result.append(s_final)
            if len(result)>5000000:
                with open(outputname+'_finalfastq.fastq','a') as ff:     
                    ff.writelines(result)
                    result=[]
        fileorder+=1
        f.close()
    print(len(result))
    if len(result)>0:
        with open(outputname+'_finalfastq.fastq','a') as ff:
            ff.writelines(result)
 

if __name__=='__main__':

    args={'step':5,
          'binsize':30,
          'filter':40,
          'outputname':'2384-EPN-CSF_S1_L001_R1_001_trimmed',
          'originalfile':['2384-EPN-CSF_S1_L001_R1_001_trimmed.fq.gz'],
          'mapfilenumber':10,
          'finish':1
    }

    epn=['2384-EPN-CSF_S1_L001_R1_001_trimmed.part1.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part2.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part3.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part4.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part5.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part6.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part7.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part8.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part9.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part10.fastq',
         '2384-EPN-CSF_S1_L001_R1_001_trimmed.part11.fastq']

    names = reads_map(epn,args)
    print(names)
    reads_reduce(names,args)
#    step = args['step']
##    length_bin = args['binsize']
#    filter = args['filter']
#    outputname = args['outputname']
#    originalfile = args['originalfile']


from utils import *
import os

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
        GetFastqList(result,step,length_bin,filter,outputname,originalfile)
        #GetFastqList(result,step,length_bin,filter,outputname,originalfile)
        #totalresult.update(result)
        print(str(i)+' finished! length='+str(len(totalresult)))
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
    result_start=[]
    result_end=[]
    for file in originalfile:
        with open(file) as f:
            while True:
                name = f.readline()
                if not name:
                    break
                reads = f.readline().strip()
                _ = f.readline()
                quality = f.readline().strip()
                fqname = name.strip().split()[0][1:]
                if not fqname in pos_mark[fileorder]: continue
                for i in range(len(pos_mark[fileorder][fqname])):
                    start,end = pos_mark[fileorder][fqname][i]
                    s_name = fqname
                    if len(pos_mark[pos])>1:
                        s_name += '_'+str(i)
                    s_read = reads[:start]
                    #0-start
                    s_qua = quality[:start]
                    readlen = str(len(reads))
                    s_final = '@'+s_name+'_'+str(fileorder)+':'+readlen+'\n'+s_read+'\n'+'+\n'+s_qua+'\n'
                    result_start.append(s_final)
                    s_read = reads[end:]
                    s_qua = quality[end:]
                    s_final = '@'+s_name+'_'+str(fileorder)+':'+readlen+'\n'+s_read+'\n'+'+\n'+s_qua+'\n'
                    result_end.append(s_final)
                if len(result_start)>5000000:
                    with open(outputname+'_clipped_start.fastq','a') as ff:     
                        ff.writelines(result_start)
                        result_start=[]
                    with open(outputname+'_clipped_end.fastq','a') as ff:
                        ff.writelines(result_end)
                        result_end=[]
        fileorder+=1
    if len(result_start)>0:
        with open(outputname+'_clipped_start.fastq','a') as ff:     
            ff.writelines(result_start)
        with open(outputname+'_clipped_end.fastq','a') as ff:
            ff.writelines(result_end)

 

if __name__=='__main__':

    args={'step':5,
          'binsize':30,
          'filter':40,
          'outputname':'6P',
          'originalfile':['6P_R1_val_1.fq','6P_R2_val_2.fq'],
          'mapfilenumber':10,
          'finish':1
    }
    name_6g=['6G.part1.fastq','6G.part2.fastq','6G.part3.fastq','6G.part4.fastq',
             '6G.part5.fastq','6G.part6.fastq','6G.part7.fastq','6G.part8.fastq',
             '6G.part9.fastq','6G.part10.fastq','6G.part11.fastq','6G.part12.fastq',
             '6G.part13.fastq','6G.part14.fastq','6G.part15.fastq','6G.part16.fastq',
             '6G.part17.fastq','6G.part18.fastq','6G.part19.fastq','6G.part20.fastq',
             '6G.part21.fastq','6G.part22.fastq','6G.part23.fastq','6G.part24.fastq']

    name_6p=['6P.part1.fastq','6P.part2.fastq','6P.part3.fastq','6P.part4.fastq',
             '6P.part5.fastq','6P.part6.fastq','6P.part7.fastq','6P.part8.fastq',
             '6P.part9.fastq','6P.part10.fastq','6P.part11.fastq','6P.part12.fastq',
             '6P.part13.fastq','6P.part14.fastq','6P.part15.fastq','6P.part16.fastq',
             '6P.part17.fastq','6P.part18.fastq','6P.part19.fastq','6P.part20.fastq',
             '6P.part21.fastq','6P.part22.fastq','6P.part23.fastq','6P.part24.fastq']

    name_m6g=['M6G.part1.fastq','M6G.part2.fastq','M6G.part3.fastq','M6G.part4.fastq',
             'M6G.part5.fastq','M6G.part6.fastq','M6G.part7.fastq','M6G.part8.fastq',
             'M6G.part9.fastq','M6G.part10.fastq','M6G.part11.fastq','M6G.part12.fastq',
             'M6G.part13.fastq','M6G.part14.fastq','M6G.part15.fastq','M6G.part16.fastq',
             'M6G.part17.fastq','M6G.part18.fastq','M6G.part19.fastq','M6G.part20.fastq',
             'M6G.part21.fastq','M6G.part22.fastq','M6G.part23.fastq','M6G.part24.fastq']

    names = reads_map(name_6p,args)
    print(names)
    reads_reduce(names,args)
#    step = args['step']
##    length_bin = args['binsize']
#    filter = args['filter']
#    outputname = args['outputname']
#    originalfile = args['originalfile']


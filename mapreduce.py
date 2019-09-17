from utils import *
import os
import gzip

def reads_map(unmapped_file,args):
    mapfilenum = args['mapfilenumber']
    step = args['step']

    file_order=0
    mapreduce_file=[]
    rootfile = unmapped_file[:unmapped_file.find('.')]
    dic={}
    for i in range(mapfilenum):
        mapreduce_file.append(rootfile+'_'+str(i)+'.mapreduce')
        if os.path.exists(mapreduce_file[-1]) and (not 'finish' in args):
            os.remove(mapreduce_file[-1])
            print('Delete '+mapreduce_file[-1])
        dic[i]=[]
    if 'finish' in args:
        return mapreduce_file

    count=0
    with open(unmapped_file+'.sam') as f:
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
            s[0] = s[0].strip()
            read_name = s[0].split('&')[0]
            file_order = int(s[0].split('&')[1]) - 1
            hashnum = abs(hash(read_name)) % mapfilenum
            dic[hashnum].append([read_name,s[2],s[3],str(file_order),str(mismatch),str(tail_mismatch),str(read_length)])
            count+=1
            if count>5000000:
                for i in range(mapfilenum):
                    arr = dic[i]
                    arr = list(map(lambda x:'\t'.join(x)+'\n',arr))
                    with open(mapreduce_file[i],'a') as ff:
                        ff.writelines(arr)
                    dic[i]=[]
                count=0
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
                        # break

                # frac_list=result[name]
                if not join_or_not:
                    result[name].append(temp)
            #print(len(result))
    #Delete short fragments
    # print(len(result))
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
    # print(len(result))
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
    if os.path.exists(outputname+'_finalfastq.fastq'):
        os.system('rm '+outputname+'_finalfastq.fastq')
    totalresult={}
    for i in range(mapfilenum):
        print(str(i)+' start!')
        command = 'sort -k4n,4 -o '+mapreduce_file[i]+' '+mapreduce_file[i]
        p = Pshell(command)
        p.change(command)
        p.process()
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
    if len(joined_reads)==0: return
    nameset={}
    #Generate a dictionary which contains readsname, start file order and extend fraction number
    for name in joined_reads:
        reads_list = joined_reads[name]
        if len(reads_list)==0: continue
        n = name
        # print(n,reads_list)
        nameset[n]=[[read[2],read[4]] for read in reads_list]
        #contentset[n]=[['',''] for i in range(len(nameset[n]))]#read_content,read_quality
    # print(len(nameset))
    pos_mark=[{},{}]
    for name in nameset:
        readinfo = nameset[name]

        # print(name, readinfo)

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
    del nameset
    fileorder=0
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
            if '/' in fqname: # or '.' in fqname:
                # split_pos=0
                # if '/' in fqname:
                split_pos = fqname.rfind('/')
                # else:
                #     split_pos = fqname.find('.')
                fqname = fqname[:split_pos]
            if not fqname in pos_mark[fileorder]: continue
            for i in range(len(pos_mark[fileorder][fqname])):
                start,end = pos_mark[fileorder][fqname][i]
                s_name = fqname
                if len(pos_mark[pos])>1:
                    s_name += '_'+str(i)
                s_read = reads[start:end]
                s_qua = quality[start:end]
                leftdel = start
                rightdel = len(reads)-end
                if rightdel<0:
                    rightdel=0
                s_final = '@'+s_name+'_'+str(fileorder)+'_'+str(leftdel)+'_'+str(rightdel)+'\n'+s_read+'\n'+'+\n'+s_qua+'\n'
                result.append(s_final)
            if len(result)>5000000:
                with open(outputname+'_finalfastq.fastq','a') as ff:     
                    ff.writelines(result)
                    result=[]
        fileorder+=1
        f.close()
    # print(len(result))
    if len(result)>0:
        with open(outputname+'_finalfastq.fastq','a') as ff:
            ff.writelines(result)
 

if __name__=='__main__':

    args={'step':5,
          'binsize':30,
          'filter':40,
          'outputname':'hpv_test',
          'originalfile':['D15_S1_L001_R1_001.fastq.gz'],
          'mapfilenumber':2,
          'finish':1,
    }

    mr_file = ['D15_S1_L001_R1_001_trimmed_0.mapreduce',
               'D15_S1_L001_R1_001_trimmed_1.mapreduce']
    #names = reads_map('6P1_notrim_val_1_val_1_test.unmapped.fastq',args)
    #print(names)
    reads_reduce(mr_file,args)
#    step = args['step']
##    length_bin = args['binsize']
#    filter = args['filter']
#    outputname = args['outputname']
#    originalfile = args['originalfile']


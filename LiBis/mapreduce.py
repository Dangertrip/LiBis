from .utils import *
import os
import gzip
import pysam

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
            #print('Delete '+mapreduce_file[-1])
        dic[i]=[]
    if 'finish' in args:
        return mapreduce_file

    count=0

    with pysam.AlignmentFile(unmapped_file+'.bam','r') as f:
        for line in f:
            #s = line.strip().split('\t')
            mismatch = line.tags[0][1]#int(s[11][s[11].rfind(':')+1:])
            if mismatch>1: continue
            read_length = line.query_length#len(s[9])
            tail_length = ((read_length-1)%step)+1
            refseq = line.tags[1][1][-2-tail_length:-2]#s[12][-2-tail_length:-2]
            readsseq = line.seq[-tail_length:]#s[9][-tail_length:]
            strand = line.tags[2][1]#s[13][-2:]
            mis=0
            for base in range(tail_length):
                if (refseq[base]!=readsseq[base]):
                    if strand[0]=='+':
                        if (refseq[base]=='C' and readsseq[base]=='T'): continue
                    else:
                        if (refseq[base]=='G' and readsseq[base]=='A'): continue
                    mis+=1
            tail_mismatch = mis
            query_name = line.query_name.split('&')
            #s[0] = s[0].strip()
            read_name = query_name[0]
            file_order = int(query_name[1]) 
            hashnum = abs(hash(read_name)) % mapfilenum
            dic[hashnum].append([read_name,line.reference_name,str(line.reference_start),str(file_order),str(mismatch),str(tail_mismatch),str(read_length)])
            count+=1
            if count>10000000:
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
            if name not in result:
                result[name]=[[chr,startpos,fileorder,mismatch,tail_mismatch,read_length]]
            else:
                result[name].append([chr,startpos,fileorder,mismatch,tail_mismatch,read_length])
    print(len(result))
    for name, content_list in result.items():
        content_list = sorted(content_list, key=lambda x:x[2])
        # initialize result_fragment
        chr,startpos,fileorder,mismatch,tail_mismatch,read_length = content_list[0]
        result_fragment = [[chr,startpos,fileorder,mismatch,0]]
        # Fragment combination
        for c in content_list[1:]:
            chr,startpos,fileorder,mismatch,tail_mismatch,read_length = c
            temp = [chr,startpos,fileorder,mismatch,0]
            join_or_not=False
            for reads in result_fragment:
                if reads[3]+tail_mismatch<=1 and readsjoin(reads,temp,step,read_length,length_bin):
                    reads[3]+=tail_mismatch
                    reads[4]=temp[2]-reads[2]
                    join_or_not=True
            if not join_or_not:
                result_fragment.append(temp)
        frag_num = len(result_fragment)
        for i in range(frag_num-1,-1,-1):
            if result_fragment[i][4]<=1 and result_fragment[i][3]>0:
                result_fragment.pop(i)
        frag_num = len(result_fragment)
        del_mark = [0] * frag_num
        for i in range(frag_num):
            for j in range(i+1,frag_num):
                if overlap(result_fragment[i],result_fragment[j],step,length_bin):
                    sss = result_fragment[i][4] - result_fragment[j][4]
                    if sss>0: del_mark[j]=1
                    elif sss<0: del_mark[i]=1
                    else:
                        mis = result_fragment[i][3] - result_fragment[j][3]
                        if mis>0: del_mark[i]=1
                        else: del_mark[j]=1
        for i in range(frag_num-1,-1,-1):
            if del_mark[i]==1:
                result_fragment.pop(i)
        result[name] = result_fragment
    return result


    

def reads_reduce(mapreduce_file,args):
    step = args['step']
    length_bin = args['binsize']
    filter = args['filter']
    outputname = args['outputname']
    originalfile = args['originalfile']
    mapfilenum = args['mapfilenumber']
    report_clip = args['report_clip']
    if report_clip==None:
        report_clip = False
    removeFileIfExist(outputname+'_finalfastq.fastq.gz')
    if report_clip:
        removeFileIfExist(outputname+'_finalfastq_clipped_head.fastq.gz')
        removeFileIfExist(outputname+'_finalfastq_clipped_tail.fastq.gz')
    totalresult=[{},{}]
    for i in range(mapfilenum):
        result=reads_combine(mapreduce_file[i],args)
        
        #print(len(result))

        pos_mark = [{},{}]
        for name, reads_list in result.items():#nameset:
            #readinfo = nameset[name]
            if len(reads_list)==0: continue
            pos = 0
            
            if name[-2:]=='_1' or name[-2:]=='_2':
                if name[-2:]=='_2':
                    pos = 1
                name = name[:-2]
            for read in reads_list:#readinfo:
                order, sum = read[2], read[4]
                start = (order)*step
                end = start + step*sum + length_bin
                if end-start<filter: 
                    continue
                if name in pos_mark[pos]:
                    pos_mark[pos][name].append([start,end])
                else:
                    pos_mark[pos][name]=[[start,end]]
 


        totalresult[0].update(pos_mark[0])
        totalresult[1].update(pos_mark[1])
        content_number = len(totalresult[0])+len(totalresult[1])
        if content_number>5000000 or (i==mapfilenum-1 and content_number>0):
            GetFastqList(totalresult,step,length_bin,filter,outputname,originalfile,report_clip)
            totalresult={}


def GetFastqList(joined_reads,step,length_bin,filter,outputname,originalfile,report_clip):
    pos_mark = joined_reads
    fileorder=0
    result=[]
    result_begin=[]
    result_end=[]
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
            if '/' in fqname[-5:]: # or '.' in fqname:
                split_pos = fqname.rfind('/')
                fqname = fqname[:split_pos]
            if not fqname in pos_mark[fileorder]: continue
            for i in range(len(pos_mark[fileorder][fqname])):
                start,end = pos_mark[fileorder][fqname][i]
                s_name = fqname
                #if len(pos_mark[pos])>1:
                #    s_name += '_'+str(i)
                s_read = reads[start:end]
                if len(s_read)<filter:
                    continue
                s_qua = quality[start:end]
                if report_clip:
                    s_read_begin = reads[:start]
                    s_read_end = reads[end:]
                    s_qua_begin = quality[:start]
                    s_qua_end = quality[:end]
                leftdel = start
                rightdel = len(reads)-end
                if rightdel<0:
                    rightdel=0
                s_final = '@'+s_name+'_'+str(i)+'_'+str(fileorder)+'_'+str(leftdel)+'_'+str(rightdel)+'\n'+s_read+'\n'+'+\n'+s_qua+'\n'
                s_final = s_final.encode()
                if report_clip:
                    s_begin = '@'+s_name+'_'+str(fileorder)+'_'+str(leftdel)+'_'+str(rightdel)+'\n'+s_read_begin+'\n'+'+\n'+s_qua_begin+'\n'
                    s_end = '@'+s_name+'_'+str(fileorder)+'_'+str(leftdel)+'_'+str(rightdel)+'\n'+s_read_end+'\n'+'+\n'+s_qua_end+'\n'
                    s_begin, s_end = s_begin.encode(), s_end.encode()
                    if len(s_read_begin)>0:
                        result_begin.append(s_begin)
                    if len(s_read_end)>0:
                        result_end.append(s_end)
                result.append(s_final)
            if len(result)>2000000:
                with gzip.open(outputname+'_finalfastq.fastq.gz','a') as ff:     
                    ff.writelines(result)
                if report_clip:
                    with gzip.open(outputname+'_finalfastq_clipped_head.fastq.gz','a') as ff:
                        ff.writelines(result_begin)
                    with gzip.open(outputname+'_finalfastq_clipped_tail.fastq.gz','a') as ff:
                        ff.writelines(result_end)
                    result_begin=[]
                    result_end=[]
                result=[]
                
        fileorder+=1
        f.close()
    # print(len(result))
    if len(result)>0:
        with gzip.open(outputname+'_finalfastq.fastq.gz','a') as ff:
            ff.writelines(result)
        if report_clip:
            with gzip.open(outputname+'_finalfastq_clipped_head.fastq.gz','a') as ff:
                ff.writelines(result_begin)
            with gzip.open(outputname+'_finalfastq_clipped_tail.fastq.gz','a') as ff:
                ff.writelines(result_end)


if __name__=='__main__':

    args={'step':5,
          'binsize':30,
          'filter':30,
          'outputname':'mate1_bsmap',
          'originalfile':['mate1.fq.gz','mate2.fq.gz'],
          'mapfilenumber':10,
          #'finish':1,
          'report_clip':1
    }

    mr_file = ['mate1_bsmap_0.mapreduce',
               'mate1_bsmap_1.mapreduce',
               'mate1_bsmap_2.mapreduce',
               'mate1_bsmap_3.mapreduce',
               'mate1_bsmap_4.mapreduce',
               'mate1_bsmap_5.mapreduce',
               'mate1_bsmap_6.mapreduce',
               'mate1_bsmap_7.mapreduce',
               'mate1_bsmap_8.mapreduce',
               'mate1_bsmap_9.mapreduce']
    names = reads_map('mate1_bsmap.unmapped.fastq',args)
    #print(names)
    # reads_reduce(mr_file,args)
    reads_reduce(names, args)
#    step = args['step']
##    length_bin = args['binsize']
#    filter = args['filter']
#    outputname = args['outputname']
#    originalfile = args['originalfile']


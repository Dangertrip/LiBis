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


def combine(read1,read2,strand,step):
    read_length = read2.reference_length
    tail_length = ((read_length-1)%step)+1
    if strand=='++' or strand=='--':
        refseq = read2.get_tag("XR")[-2-tail_length:-2]
        add_seq = read2.seq[-1*tail_length:]
        add_qua = read2.query_qualities[-1*tail_length:]
        read1.query_qualities.extend(add_qua)
        read1_qua = read1.query_qualities
        read1.seq = read1.seq + add_seq
        read1.query_qualities = read1_qua
    else:
        refseq = read2.get_tag("XR")[2:2+tail_length]
        add_seq = read2.seq[:tail_length]
        add_qua = read2.query_qualities[:tail_length]
        add_qua.extend(read1.query_qualities)
        read1.seq = add_seq + read1.seq
        read1.query_qualities = add_qua
        
    read1.reference_start = min(read1.reference_start,read2.reference_start)
    if read2.get_tag("NM")>0:
        new_mismatch = read1.get_tag("NM")
        for i in range(len(refseq)):
            if refseq[i]==add_seq[i]: continue
            if strand[0]=='+':
                if (refseq[i]=='C' and add_seq[i]=='T'): continue
            else:
                if (refseq[i]=='G' and add_seq[i]=='A'): continue
            new_mismatch+=1
        read1.set_tag("NM",new_mismatch)


def fragCombine(param):#,outf):
    frags = param['frags']
    order = param['order']
    args = param['args']
    step = args['step']
    filter = args['filter']
    window = args['binsize']
    order_max_gap = window//step
    #print(len(frags),order)
    strand = []
    mismatch = []
    for fr in frags:
        strand.append(fr.get_tag("ZS"))
        mismatch.append(fr.get_tag("NM"))
        
    head_frag_length_diff = 0
    if order[0]==0:
        head_frag_length_diff = abs(window - frags[0].reference_length)

    result_fragment = [frags[0]]
    result_order = [order[0]]
    result_strand = [strand[0]]
    result_extension = [0]
    for i in range(1,len(frags)):
        join_marker = False
        candidate = frags[i]
        candidate_order = order[i]
        for j, read in enumerate(result_fragment):
            if read.reference_name==candidate.reference_name and strand[i]==result_strand[j]:
                read_order = result_order[j]
                order_diff = candidate_order - read_order
                if strand[i]=='+-' or strand[i]=='-+':
                    order_diff = order_diff - result_extension[j]
                if order_diff>=order_max_gap: continue
                pos_diff = abs(read.reference_start - candidate.reference_start)
                #print(candidate_order,read_order,result_extension[j],pos_diff, order_diff)
                #print(result_extension)
                if order_diff*step==pos_diff or (j==0 and order_diff*step==pos_diff+head_frag_length_diff):
                    join_marker = True
                    combine(read, candidate,result_strand[j],step)
                    update_read(read)
                    result_extension[j] += 1
        if not join_marker:
            result_fragment.append(frags[i])
            result_order.append(order[i])
            result_strand.append(strand[i])
            result_extension.append(0)

    # Finish combining fragments. Now apply filters towards combined fragments:
    # 1. 
    valid_frags = []
    for i, frag in enumerate(result_fragment):
        #print(frag.to_string())
        length = len(frag.seq)
        if length>=filter:# and frag.get_tag("NM")<=length*0.05:
            valid_frags.append([frag,result_order[i],result_extension[i]])

    result_fragment = valid_frags
    valid_frags = []
    del_mark = [0] * len(result_fragment)
    frag_num = len(result_fragment)
    for i in range(frag_num):
        for j in range(i+1,frag_num):
            read1, order1, ext1 = result_fragment[i]
            read2, order2, ext2 = result_fragment[j]
            gap = order2 - order1
            if gap<ext1+order_max_gap:
                ext_diff = ext1-ext2
                if ext_diff>0:
                    del_mark[j]=1
                elif ext_diff<0: del_mark[i]=1
                else:
                    mis_diff = read1.get_tag("NM") - read2.get_tag("NM")
                    if mis_diff>0: del_mark[i]=1
                    else: del_mark[j]=1

    read_length = []
    for i,frag in enumerate(result_fragment):
        if del_mark[i]==0:
            valid_frags.append(result_fragment[i])
            read_length.append(result_fragment[i][0].reference_length)

    frag_rank = 0

    combined_frags = []
    for i in range(len(valid_frags)):
        frag, order, ext = valid_frags[i]
        #frag.query_name = frag.query_name + '_' + str(frag_rank) + '_' + str(order*step) + '_' + str(window+(order+ext)*step)
        name = frag.query_name
        if '_' in name[-5:]:
            name = frag.query_name.split('_')
            frag.query_name = name[0]
            frag.set_tag("MT", name[1])
        else:
            frag.set_tag("MT", '3')
        frag.set_tag("FG", str(frag_rank))
        frag.set_tag("LC", str(order*step))
        frag.set_tag("RC", str(window+(order+ext)*step))
        frag.set_tag('XR', None)
        # print(frag.to_string())
        #outf.write(frag)
        combined_frags.append(frag)
        frag_rank += 1

    return len(combined_frags), read_length, combined_frags 


def pairwise(combined_frags):
    for m1_frag in combined_frags[0]:
        for m2_frag in combined_frags[1]:
            r1 = m1_frag
            r1_mate = r1.get_tag('MT')
            r1_frag = r1.get_tag('FG')
            r1_chr = r1.reference_name
            r1_start = r1.reference_start
            r1_len = r1.reference_length
            r1_end = r1_start + r1_len
            r2 = m2_frag
            r2_mate = r2.get_tag('MT')
            r2_frag = r2.get_tag('FG')
            r2_chr = r2.reference_name
            r2_start = r2.reference_start
            r2_len = r2.reference_length
            r2_end = r2_start + r2_len
            if r1_chr==r2_chr:
                template_length = max(r1_end, r2_end) - min(r1_start, r2_start)
                if template_length<r1.template_length or r1.next_reference_start<0 or r1.next_reference_name!=r1_chr:
                    r1.template_length = template_length
                    r1.next_reference_name = r2_chr
                    r1.next_reference_start = r2_start
                    if r1_mate==1:
                        m_info=64
                    else:
                        m_info=128
                    if r2.flag&16>0:
                        r1.flag = r1.flag|32
                    r1.flag = r1.flag|1
                    r1.flag = r1.flag|m_info
                if template_length<r2.template_length or r2.next_reference_start<0 or r2.next_reference_name!=r2_chr:
                    r2.template_length = template_length
                    r2.next_reference_name = r1_chr
                    r2.next_reference_start = r1_start
                    if r2_mate==1:
                        m_info=64
                    else:
                        m_info=128
                    if r1.flag&16>0:
                        r2.flag = r2.flag|32
                    r2.flag = r2.flag|1
                    r2.flag = r2.flag|m_info
            else:
                if r1.next_reference_start<0:
                    r1.next_reference_name = r2_chr
                    r1.next_reference_start = r2_start
                    if r1_mate==1:
                        m_info=64
                    else:
                        m_info=128
                    if r2.flag&16>0:
                        r1.flag = r1.flag|32
                    r1.flag = r1.flag|1
                    r1.flag = r1.flag|m_info
                if r2.next_reference_start<0:
                    r2.next_reference_name = r1_chr
                    r2.next_reference_start = r1_start
                    if r2_mate==1:
                        m_info=64
                    else:
                        m_info=128
                    if r1.flag&16>0:
                        r2.flag = r2.flag|32
                    r2.flag = r2.flag|1
                    r2.flag = r2.flag|m_info



def unsortedCombine(unmapped_file,args):
    step = args['step']
    outputname = args['outputname']
    window = args['binsize']
    filter = args['filter']
    inputfilenumber = len(args['originalfile'])
    least_frags = (filter-window)//step

    result_fragment = []

    feed_in_parameter = {
        'args' : args,
    }

    reads_num_dist = [0]*10
    reads_len_dist = [0]*200

    with pysam.AlignmentFile(unmapped_file+'.bam','r') as f:
        wf = pysam.AlignmentFile(outputname+"_split.bam", "wb", template=f)
        frags = [[],[]]
        order = [[],[]]
        current_reads_name=''
        for line in f:
            mismatch = line.get_tag('NM')
            if mismatch>line.query_length*0.05: continue
            #print(line.query_name)
            if '&' not in line.query_name[-5:]:
                continue
            read_name, o = line.query_name.split('&')

            line.query_name = read_name
            if inputfilenumber>1:#'_' in line.query_name[-5:]:
                c = line.query_name.split('_')
                pure_name = c[0]
                _mate = int(c[-1])-1
            else:
                pure_name = line.query_name
                _mate = 0

            #if 'E00488:423:HYHFMCCXY:8:1101:14874:1872' not in read_name: continue
            #print(line.to_string())
            #if read_name!=current_reads_name:
            if pure_name!=current_reads_name:
                #print(current_reads_name)
                if current_reads_name!='':
                    combined_frags = [[],[]]
                    
                    # Combined fragments retrieving
                    for m in range(2):
                        if len(frags[m])>least_frags:
                            feed_in_parameter['frags'] = frags[m]
                            feed_in_parameter['order'] = order[m]
                            reads_num, reads_len, combined_frags[m] = fragCombine(feed_in_parameter)#fragCombine(feed_in_parameter,wf)
                            if reads_num>0:
                                reads_num_dist[reads_num] += 1
                                for l in reads_len:
                                    reads_len_dist[l] += 1
                    #Combined fragments pairing
                    if len(combined_frags[0])>0 and len(combined_frags[1])>0:
                        pairwise(combined_frags)
                    
                    #Combined fragments output
                    for m in range(2):
                        for outfrag in combined_frags[m]:
                            wf.write(outfrag)

                    frags = [[],[]]
                    order = [[],[]]
                frags[_mate].append(line)
                order[_mate].append(int(o))
                current_reads_name = pure_name # read_name
            else:
                frags[_mate].append(line)
                order[_mate].append(int(o))
    if len(frags)>least_frags:
        combined_frags = [[],[]]
        for m in range(2):
            if len(frags[m])>least_frags:
                feed_in_parameter['frags'] = frags[m]
                feed_in_parameter['order'] = order[m]
                reads_num, reads_len, combined_frags[m] = fragCombine(feed_in_parameter)
                if reads_num>0:
                    reads_num_dist[reads_num] += 1
                    for l in reads_len:
                        reads_len_dist[l] += 1
        if len(combined_frags[0])>0 and len(combined_frags[1])>0:
            pairwise(combined_frags)
        for m in range(2):
            for outfrag in combined_frags[m]:
                wf.write(outfrag)
        frags = [[],[]]
        order = [[],[]]
    
    return reads_num_dist,reads_len_dist

if __name__=='__main__':

    args={'step':5,
          'binsize':40,
          'filter':46,
          'outputname':'sample1',
          #'originalfile':['mate1.fq.gz','mate2.fq.gz'],
          #'mapfilenumber':10,
          #'finish':1,
          #'report_clip':1
    }
    #unsortedCombine('/data/yyin/LiBis_test/data/simulation/rht_1sample/LiBis_011_w40_unmap/1_1.unmapped',args)
    unsortedCombine('/scratch0/MB4CSF/MB4.R1.fq.ae.unmapped', args)
    #unsortedCombine('test_mr', args)

'''
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
'''

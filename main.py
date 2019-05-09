import argparse
import os
import xml.dom.minidom
from computeProcess import computeProcess
from utils import *
'''
    fastq pairs             k pairs
    sample label            k labels
    clipmode                True/False
    clipsetting             window, step
    ref                     string
    processes               number
'''

def nameprocess(s):
    #split samples using blank. split double end files using comma(,);
    ans=[]
    for f in s:
        if ',' in f:
            ff=f.strip().split(',')
        else:
            ff = [f]
        ans.append(ff)
    return ans


def getxmlcontent(dom,x):
    temp = dom.getElementsByTagName(x)
    if len(temp)==0:
        return None
    else:
        return temp[0].firstChild.data


def text_process(filename):
    '''
    All process for reading config text in XML format.
    This only works for the mode which use file to input parameters
    '''
    dom = xml.dom.minidom.parse(filename)
    #root = dom.documentElement
    dic={}
    tagnames = ['fastq','label','clip','window','step','ref','process','genome','qc','trim','binsize','filter']
    tagcontent = list(map(lambda x:getxmlcontent(dom,x),tagnames))
    for i in range(len(tagnames)):
        name = tagnames[i]
        content = tagcontent[i]
        if name=='fastq':
            dic[name]=nameprocess(content.strip().split())
            '''
            Split fastq file name using comma(,)
            Format will be like S1_1,S1_2 S2 S3 S4_1,S4_2
            feed nameprocess() a name list like above(use split to eliminate the blank)
            it will return [['S1_1','S1_2'],['S2'],['S3'],['S4_1','S4_2']]
            '''
            continue
        if name=='label':
            dic[name]=content.strip().split()
            continue
        if name=='ref':
            dic[name]=content
            continue
        if name=='genome':
            dic[name]=content
        dic[name]=int(content)
        
    return dic

def inputorN(v):
    if v:
        return v
    else:
        return None

def valid(param):
    '''
    check whether parameters are valid.
    '''
    name = param['name']
    for n in name:
        for nn in n:
            if not os.path.exists(nn):
                raise Exception(nn+' not exist!')
    if len(param['label'])!=len(param['name']):
        raise Exception('Number of samples and number of labels should be the same!')
    if not os.path.exists(param['ref']):
        raise Exception(param['ref']+' not exist!')
    if param['clip'] and param['trim']:
        print("Don't need to trim for clipping mode! Set -t to 0")
    if os.path.exists('RESULT'):
        raise "RESULT file exists! Please delete RESULT or change name"
    os.mkdir("RESULT")

def input_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file',type=int, help=r'Required. Enter a number, 0 means using parameter to set up, 1 means using text file to set up',required=True,default=0)
    parser.add_argument('-sf','--settingfile',help='Required. Setting txt file name. Ignore if -f is 0.')
    parser.add_argument('-n','--name',nargs="*",help=r'Required. Fastq file name. ')
    parser.add_argument('-c','--clip',help=r'Clip mode. 0 means close. 1 means open. default: 0',type=int,default=0)
    parser.add_argument('-l','--label',nargs="*",help=r'Required. Labels for samples')
    parser.add_argument('-g','--genome',help='Required. Genome the reference belong to.(Use for plotting) hg18/hg19/mm10/mm9 and so on. Plotting script will not avaliable if leave it blank')
    parser.add_argument('-w','--window',type=int,help=r'Window length for clipping mode, default=1000000',default=30)
    parser.add_argument('-s','--step',type=int,help=r'Step size for clipping mode.',default=5)
    parser.add_argument('-p','--process',type=int,help=r'Process using for one pipeline. Normally bsmap will cost 8 cpu number. So total will be 8p.',default=1)
    parser.add_argument('-r','--ref',help=r'Required. Reference')
    parser.add_argument('-qc','--QualityControl',help=r'Do(1) quality control or not(0)',default=True,required=False)
    parser.add_argument('-t','--trim',help=r"Do(1) trimming or not(0). Don't need to do trimming if you use clip mode.",default=True,required=False)
    parser.add_argument('-b','--binsize',help="Plot setting. Set the bin size for averaging methylation ratio among samples",default=1000000,required=False)
    parser.add_argument('-ft','--filter',help="filter for clipped reads",default=40)

    args = parser.parse_args()
    #print(args.name) 
    if args.file==1:
        param=text_process(args.settingfile)
    else:
        param={'name':(nameprocess(args.name) or None), 
               'clip':args.clip, 
               'label':(args.label or None), 
               'window':int(inputorN(args.window)), 
               'step':int(inputorN(args.step)), 
               'process':int(inputorN(args.step)), 
               'ref':(args.ref or None), 
               'qc':int(args.QualityControl),
               'trim':int(args.trim),
               'genome':args.genome,
               'bin':args.binsize,
               'filter_len':args.filter
              }
    valid(param)
    return param

if __name__=="__main__":
    param=input_args()
    computeProcess(param)



import subprocess

def processline(line,dic,order,cov):
    if line[0]=='#': return
    t = line.strip().split()
    key=(t[0],t[1])
    if len(t)>4 and int(t[4])<cov: return
    if not key in dic:
        if order==0:
            dic[key]=[t[3]]
    else:
        dic[key].append(t[3])

def formdata(files,cov=0,bedfile=''):
    dic={}
    i=0
    for file in files:
        if bedfile=='':
            with open(file) as f:
                for line in f:
                    processline(line,dic,i,cov)
        else:
            p = subprocess.Popen('bedtools intersect -a %s -b %s' %(file,bedfile),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            for line in p.stdout:
                processline(line,dic,i,cov)
        i+=1
    result=[]
    for key in dic:
        if len(dic[key])==len(files):
            result.append(dic[key])
    return result

if __name__=="__main__":
    import sys
    argv = sys.argv[1:]
    final=[argv]
    final.extend(formdata(argv))
    ans=[]
    for f in final:
        s=str(f[0])
        for i in range(1,len(f)):
            s=s+','+str(f[i])
        s=s+'/n'
        ans.append(s)
    with open('data.csv','w') as f:
        f.writelines(ans)

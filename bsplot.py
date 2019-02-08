'''
Tested
'''
import matplotlib
matplotlib.use('Agg')
import numpy as np 
import matplotlib.pyplot as plt 
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import os
from urllib import request
from copy import deepcopy
import matplotlib.cm as cm

def point_cluster(data,outputname,method='PCA'):
    #data: DataFrame
    #contains labelname column, samplename column and data
    d = deepcopy(data)
    #print(di)
    d.sort_values(['chrom','start'])

    windowdata = d.values[:,3:].T
    position = d.values[:,:3]
    label = d.columns[3:]
    #Every sample in a row
    #print(windowdata)
    dim=2
    colors = cm.rainbow(np.linspace(0,1,len(label)))
    if method == 'PCA':
        pca = PCA(n_components=dim)
        x_tr = pca.fit_transform(windowdata)
    if method == 'TSNE':
        tsne = TSNE(n_components=dim)
        x_tr = tsne.fit_transform(windowdata)
    fig,ax = plt.subplots(figsize=(7,7))
    xmin = np.min(x_tr[:,0])
    xmax = np.max(x_tr[:,0])
    ymin = np.min(x_tr[:,1])
    ymax = np.max(x_tr[:,1])
    sample_size = x_tr.shape[0]
    plt.xlim(xmin-0.2*np.abs(xmin),xmax+0.2*np.abs(xmax))
    plt.ylim(ymin-0.2*np.abs(ymin),ymax+0.2*np.abs(ymax))
    markers=['o', '^','v','<','>','1','2', '3','4','8','s','P','p', '*','H','h','x','X','D']
    label_u = np.unique(label)
    for i in range(len(label_u)):
        cc = colors[i]
        l = label_u[i]
        pos = np.where(label==l)
        ma = markers[i]
        plt.scatter(x_tr[pos,0],x_tr[pos,1],c=cc,alpha=0.8,s=50,marker=ma,label=l)

    if method=='PCA':
        plt.xlabel('PC1',fontsize=13)
        plt.ylabel('PC2',fontsize=13)
    if method=='TSNE':
        plt.xlabel('TSNE1',fontsize=13)
        plt.ylabel('TSNE2',fontsize=13)
    plt.legend(loc='best',fontsize=13)
    plt.savefig(outputname)

def heatmap(data,outputname):
    from seaborn import clustermap 
    d = deepcopy(data)
    #print(d)
    d=d.drop(['chrom','start','end'],axis=1)
    sns_plot=clustermap(d)
    sns_plot.ax_row_dendrogram.set_visible(False)
    sns_plot.savefig(outputname)


if __name__=='__main__':
    with open('BED_FILE/head_combine.bam.G.bed.short.bed') as f:
        lines = f.readlines()
    d=[]
    for line in lines:
        temp = line.strip().split()
        temp[-1]=float(temp[-1])
        d.append(temp)
    import random
    for i in range(5):
        for dd in d:
            dd.append(random.random())
    import pandas as pd
    d = pd.DataFrame(d,columns=['chrom','start','end','L1','L2','L3','L4','L4','L4'])
    heatmap(d,'a.pdf')


    


    



    

#!/usr/bin/env python3
import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

VERSION = 3.0

###################################################################################################
# Input
###################################################################################################
parser = argparse.ArgumentParser(description='.')
# IO
parser.add_argument('-i', "--input", help='the input methylation data',)
parser.add_argument('-g', "--groupinfo", help='(optional) the input group information table',)
parser.add_argument('-o', "--output", help='the output directory',)
# system
parser.add_argument('-t', "--threads", type=int, default=1, help='(optional) number of threads',)
parser.add_argument('-s', "--seed", type=int, default=1, help='(optional) set seed for random generator',)
parser.add_argument('-O', "--outputImputed", type=lambda x: (str(x).lower() == 'true'), default=False, help='(optional) True or False, save the CpG methylation matrix with imputed values as imputed.tsv',)
parser.add_argument('-p', "--verbose", type=lambda x: (str(x).lower() == 'true'), default=False, help='(optional) True or False, track the segmentation',)
# DMR
parser.add_argument('-M', "--maxdist", type=int, default=300, help='(optional) maximum distance between two CpG',)
parser.add_argument('-m', "--minCpGs", type=int, default=10, help='(optional) minimum CpGs',)
parser.add_argument('-d', "--minMethDiff", type=float, default=0.1, help='(optional) minimum mean methylation difference',)
parser.add_argument('-r', "--minDMR", type=int, default=5, help='(optional) minimum CpGs with minimum mean methylation difference in a segment',)
parser.add_argument('-X', "--minNonNA", type=int, default=1, help='(optional) minimum samples with non-NA values in each group',)
parser.add_argument('-v', "--valley", type=float, default=0.7, help='(optional) a cutoff for the difference between global and regional methylation differences',)
parser.add_argument('-D', "--minMethDiffHigh", type=float, default=0.5, help='(optional) minimum mean methylation difference for DMTree and GSEA, similar to -d but a higher value will be recommanded to reduce the number of false positive DMRs',)
parser.add_argument('-u', "--clusteringRatio", type=float, default=0.5, help='(optional) maximum ratio of CpGs with minimum difference in a cluster',)
# DMTree
parser.add_argument('-n', "--minNSamples", type=int, default=3, help='(optional) minimum samples in a cluster',)
parser.add_argument('-w', "--minSumDMRs", type=int, default=100, help='(optional) minimum sum of DMR weights to split samples',)
# optional
parser.add_argument('-plot', "--visualization", type=lambda x: (str(x).lower() == 'true'), default=False, help='(optional) plot PCA and heatmap based on DMR methylation',)
parser.add_argument('-anno', "--annotation", help='(optional) hg19 or hg38, use ChIPseeker to annotate the DMRs',)
parser.add_argument('-refs', "--refSeq", help='(optional) reference genome, for sequence annotation',)
parser.add_argument('-gsea', "--genesets", help='(optional) geneset gmt file for GSEA',)
parser.add_argument('-wsup', "--withSupervised", help='(optional) run supervised mode on clusters after unsupervised mode', type=lambda x: (str(x).lower() == 'true'), default=True)
parser.add_argument('--version', action='version', version=f'%(prog)s {VERSION}', help='Get the version of metilene3',)

# hidden
parser.add_argument('-sk', "--skipMetilene", type=lambda x: (str(x).lower() == 'true'), default=False, help=argparse.SUPPRESS)
parser.add_argument('-kt', "--keeptmp", type=lambda x: (str(x).lower() == 'true'), default=False, help=argparse.SUPPRESS)
parser.add_argument('-n0', "--minN0", type=int, default=2, help=argparse.SUPPRESS)



###################################################################################################
# Install
###################################################################################################
def getMetilene():
    os.system("cd "+os.path.realpath(__file__).replace('metilene3.py','')+";make")



###################################################################################################
# Run
###################################################################################################
def preprocess(args, headerfile, ifsup, grpinfo=None):
    # if args.skipMetilene:
    #     return None
    
    if ifsup=='unsup':
        cols = pd.read_table(args.input, nrows=0)
        newcols = list(cols.columns)
        # print(newcols)
        for i in range(len(newcols))[2:]:
            newcols[i] = str(i-2)+'_Sample'+str(i-2)
        cols.columns = newcols
        cols.to_csv(headerfile, sep='\t', index=False)
        
    else:
        cols = pd.read_table(args.input, nrows=0)
        newcols = list(cols.columns)
        
        grp = pd.read_table(grpinfo, index_col='ID')['Group'].astype(str)
        try:
            grp = grp.loc[newcols[2:]]
        except:
            print('ERROR: group information table is not matched!')
            return

        grpid = {}
        j = 0
        for i in sorted(grp.unique()):
            grpid[i] = j
            j += 1
            
        df_grpid = pd.DataFrame(pd.Series(grpid))
        df_grpid.columns = ['Group_ID']
        df_grpid.index.name = 'Group'
        df_grpid.to_csv(args.output+'/group-ID.tsv', sep='\t')
        
        grpdict = grp.map(grpid).to_dict()
        
        # print(newcols)
        for i in range(len(newcols))[2:]:
            newcols[i] = str(grpdict[newcols[i]])+'_Sample'+str(i-2)#+'_'+newcols[i]
        cols.columns = newcols
        cols.to_csv(headerfile, sep='\t', index=False)


def runMetilene(args, headerfile, ifsup):
    if args.skipMetilene:
        return None
    # print(os.path.realpath(__file__))
    if ifsup=='unsup':
        os.system(os.path.realpath(__file__).replace('metilene3.py','metilene')+\
                    " -t "+str(args.threads)+\
                    " -s "+str(args.seed)+\
                    " -p "+str(args.verbose*1)+\
                    
                    " -M "+str(args.maxdist)+\
                    " -m "+str(args.minCpGs)+\
                    " -d "+str(args.minMethDiffHigh)+\
                    " -v "+str(args.valley)+\
                    
                    " -r "+str(args.minDMR)+\
                    " -w "+str(args.minMethDiffHigh)+\
                    " -e "+str(args.clusteringRatio)+\
                    " -q "+str(args.minMethDiffHigh)+\
                    
                    " -H "+headerfile+\
                    " -l 1 "+args.input+" > "+\
                    args.output+'/DMRs-unsupervised.tsv' )

    else:
        os.system(os.path.realpath(__file__).replace('metilene3.py','metilene')+\
                    " -t "+str(args.threads)+\
                    " -s "+str(args.seed)+\
                    " -p "+str(args.verbose*1)+\
                    " -O "+str(args.outputImputed*1)+\
                    
                    " -M "+str(args.maxdist)+\
                    " -m "+str(args.minCpGs)+\
                    " -d "+str(args.minMethDiff)+\
                    " -X "+str(args.minNonNA)+\
                    " -v "+str(args.valley)+\
                    
                    " -r "+str(args.minDMR)+\
                    " -w "+str(args.minMethDiff)+\
                    " -e "+str(args.clusteringRatio)+\
                    " -q "+str(args.minMethDiff)+\
                    
                    " -H "+headerfile+\
                    " -l 1 "+args.input+" > "+\
                    args.output+'/'+args.input.split('/')[-1]+'.aout' )
                    
        if args.outputImputed:
            os.system("grep -v \'//Imputed:\' "+\
            args.output+'/'+args.input.split('/')[-1]+'.aout >' + \
            args.output+'/DMRs.tsv')
            
            os.system("head -n1 "+args.input+" > "+\
                        args.output+'/'+args.input.split('/')[-1]+'.imputed')
                        
            os.system("grep \'//Imputed:\' "+\
            args.output+'/'+args.input.split('/')[-1]+r".aout|sed 's/\/\/Imputed://' >>" + \
            args.output+'/'+args.input.split('/')[-1]+'.imputed')
            
            os.system("rm "+args.output+'/'+args.input.split('/')[-1]+'.aout')
            
        else:
            os.system("mv "+ args.output+'/'+args.input.split('/')[-1]+'.aout ' + \
            args.output+'/DMRs.tsv')


def chipseeker(mout, moutPath, anno):
    if anno in ['hg19','HG19']:
        anno = 'TxDb.Hsapiens.UCSC.hg19.knownGene'
    if anno in ['hg38','HG38']:
        anno = 'TxDb.Hsapiens.UCSC.hg38.knownGene'
    cmd = "require("+anno+");require(ChIPseeker);setwd(\'"+str(os.getcwd())+"\');"+\
    "peakfile=\'"+moutPath+".bed\';"+\
    "txdb<-"+anno+";"+\
    "peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 1000), TxDb=txdb, annoDb=\'org.Hs.eg.db\');"+\
    "write.csv(as.GRanges(peakAnno), \'"+moutPath+".bed.csv\')"
    
    mout.sort_values(['chr','start','stop',])[['chr','start','stop']].to_csv(\
    moutPath+".bed",sep='\t',index=False,header=None)
    
    os.system('Rscript -e \"'+cmd+'\"')
    
    annoed = pd.read_csv(moutPath+".bed.csv", index_col=0)
    
    os.remove(moutPath+".bed")
    os.remove(moutPath+".bed.csv")
    
    annoed.index = annoed['seqnames']+':'+annoed['start'].astype(str)+'-'+annoed['end'].astype(str)
    annoed['anno'] = annoed['annotation'].apply(lambda x:x.split(' (')[0])
    
    for i in ['distanceToTSS','ENSEMBL','SYMBOL','anno']:
        mout[i] = (mout['chr']+':'+(mout['start']+1).astype(str)+'-'+mout['stop'].astype(str)).map(annoed[i])
        
    return mout


def addSeq(mout, refSeq):
    from Bio import SeqIO

    ref = {}
    with open(refSeq) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ref[record.id] = record.seq

    mout['seq'] = mout.apply(lambda x:str(ref[x['chr']][x['start']-1:x['stop']]), axis=1)

    return mout


def processOutput(args, ifsup, anno='F'):
    if ifsup=='unsup':
        moutPath = args.output + '/DMRs-unsupervised.tsv'
    else:
        moutPath = args.output + '/DMRs.tsv'
    mout = pd.read_table(moutPath)
    
    mout = mout.loc[mout['sig.comparison']!='TBC']
    # if args.skipMetilene:
    #     return mout
    if mout.shape[0]<1:
        print("No DMR found!")
        return None
    
    mout['meandiffabs'] = mout['meandiff'].apply(abs)

    def rename_cls_pn(x):
        x = x.replace('0','1').replace('4','3')
        if x[0]=='p':
            x = x.replace('1','x').replace('3','1').replace('x','3')
        x = x[1:]
        return x
    mout['sig.comparison'] = ( (1*(mout['meandiff']>0)).map({1:"p", 0:"n"}) \
                                     +mout['sig.comparison']).apply(rename_cls_pn)
    
    mout['#Hypo'] = mout['sig.comparison'].apply(lambda x:(len(x.split('1'))-1))
    mout['#Int'] = mout['sig.comparison'].apply(lambda x:(len(x.split('2'))-1))
    mout['#Hyper'] = mout['sig.comparison'].apply(lambda x:(len(x.split('3'))-1))
    
    def calmean(a,b,c):
        a = a.split('|')
        b = b.split('|')
        s = 0
        n = 0
        for i in range(len(a)):
            if b[i]==c:
                s += float(a[i])
                n += 1
        try:
            return s/n
        except:
            return None
            
    mout['meanHypo'] = mout.apply(lambda x:calmean(x['mean'],x['sig.comparison'],'1'), axis=1)
    mout['meanInt'] = mout.apply(lambda x:calmean(x['mean'],x['sig.comparison'],'2'), axis=1)
    mout['meanHyper'] = mout.apply(lambda x:calmean(x['mean'],x['sig.comparison'],'3'), axis=1)

    def sigcom2hypo(x):
        y = []
        for i,j in enumerate(x.split('|')):
            if j=='1':
                y.append(i)
        return y

    def sigcom2inter(x):
        y = []
        for i,j in enumerate(x.split('|')):
            if j=='2':
                y.append(i)
        return y

    def sigcom2hyper(x):
        y = []
        for i,j in enumerate(x.split('|')):
            if j=='3':
                y.append(i)
        return y
    
    if ifsup=='unsup':
        sids = [str(i) for i in pd.read_table(args.input, nrows=0).columns[2:]]
        mout['Hypo-samples'] = mout['sig.comparison'].apply(sigcom2hypo).apply(lambda x:','.join(sorted([str(sids[i]) for i in x])))
        mout['Int-samples'] = mout['sig.comparison'].apply(sigcom2inter).apply(lambda x:','.join(sorted([str(sids[i]) for i in x])))\
                                                                                                                            .apply(lambda x:x if x!='' else '-')
        mout['Hyper-samples'] = mout['sig.comparison'].apply(sigcom2hyper).apply(lambda x:','.join(sorted([str(sids[i]) for i in x])))

    else:
        rename_cls = pd.read_table(args.output + '/group-ID.tsv', index_col='Group_ID')['Group'].astype(str).to_dict()

        mout['Hypo-groups'] = mout['sig.comparison'].apply(sigcom2hypo).apply(lambda x:','.join(sorted([rename_cls[i] for i in x])))
        mout['Int-groups'] = mout['sig.comparison'].apply(sigcom2inter).apply(lambda x:','.join(sorted([rename_cls[i] for i in x])))\
                                                                                                                            .apply(lambda x:x if x!='' else '-')
        mout['Hyper-groups'] = mout['sig.comparison'].apply(sigcom2hyper).apply(lambda x:','.join(sorted([rename_cls[i] for i in x])))

    # print('# of processed DMRs:',mout.shape[0])
    if anno == 'T' and args.annotation:
        mout = chipseeker(mout, moutPath, args.annotation)

    if anno == 'T' and args.refSeq:
        mout = addSeq(mout, args.refSeq)

    mout.to_csv(moutPath, index=False, sep='\t')
                    
    return mout


def addDMTree2DMR(args, ifsup, cls, finalCls):
    if ifsup=='unsup':
        moutPath = args.output + '/DMRs-unsupervised.tsv'
    else:
        moutPath = args.output + '/DMRs.tsv'
    mout = pd.read_table(moutPath)
    
    def rev123(x):
        return x.replace('3','x').replace('1','3').replace('x','1')
        
    if ifsup=='unsup':
        def rename_cls_pn2(x):
            if x.count('3')>x.count('1'):
                x = x.replace('1','2')
            elif x.count('3')<=x.count('1'):
                x = x.replace('3','2')
            return x
        mout['sig.comparison.bin'] = mout['sig.comparison'].apply(rename_cls_pn2)

        cls_id = {}
        for i in cls[0]:
            finalCls[i] = i.split('|')
            tmp = pd.crosstab(finalCls['Group'],finalCls[i])
            for j in tmp.columns:
                tmp[j] = (tmp[j]>0).map({True:str(j),False:''})
            cls_id[i] = '|'.join(list(tmp.T.sum().sort_index()))
        for i in cls[0]:
            finalCls[cls_id[i]] = finalCls[i]
        finalCls.drop(columns=cls[0]).to_csv(args.output + '/clusters.tsv', sep='\t')

        def findDMTreeID(a, b, c, d):
            a = a.split('|')
            b = b.split('|')
            num = 0
            all = 0
            for i in range(len(a)):
                if a[i] != '0':
                    if a[i] != b[i]:
                        num += 1
                    all += 1
            if num == 0 and all > 0:
                return d+c
            else:
                return ''
            
        mout['DMTree'] = ''
        for i in cls[0]:
            mout['DMTree'] += mout['sig.comparison.bin'].apply(lambda x:findDMTreeID(i,x,cls_id[i]+',','P'))
            mout['DMTree'] += mout['sig.comparison'].apply(lambda x:findDMTreeID(rev123(i),x,cls_id[i]+',','N'))

    else:
        cls_id = {}
        for i in cls[0]:
            finalCls[i] = i.split('|')
            tmp = pd.crosstab(finalCls['Group'],finalCls[i])
            for j in tmp.columns:
                tmp[j] = (tmp[j]>0).map({True:str(j),False:''})
            cls_id[i] = '|'.join(list(tmp.T.sum().sort_index()))

        def findDMTreeIDsup(a, b, c, d):
            a = a.split('|')
            if ('1' in a) and ('3' in a):
                b = b.split('|')
            else:
                if '3' in a:
                    b = b.replace('1','2').split('|')
                else:
                    b = b.replace('3','2').split('|')
            num = 0
            all = 0
            for i in range(len(a)):
                if a[i] != '0':
                    if a[i] != b[i]:
                        num += 1
                    all += 1
            if num == 0 and all > 0:
                return d+c
            else:
                return ''
            
        mout['DMTree'] = ''
        for i in cls[0]:
            mout['DMTree'] += mout['sig.comparison'].apply(lambda x:findDMTreeIDsup(cls_id[i],x,cls_id[i]+',','P'))
            mout['DMTree'] += mout['sig.comparison'].apply(lambda x:findDMTreeIDsup(rev123(cls_id[i]),x,cls_id[i]+',','N'))

    mout[mout.columns[~mout.columns.str.contains('sig.comparison.bin')]].to_csv(moutPath, index=False, sep='\t')
    return mout



###################################################################################################
# DMR-Freq-based Clustering
###################################################################################################
def recurSplit(arr, ref=0, depth=0, nsep=0, minN=2, minSumDMRs=100, fulltree=False):
    finalList = []
    depthList = []
    weightList = []
    
    def numVS(a):
        return sorted([a.count('1'),a.count('2'),a.count('3')])[1]
    
    if ref == 0:
        arr = arr.groupby('sig.comparison.bin').sum().sort_values(ascending=False)
        ifSig = 0
        for i in range(len(arr)):
            newref = arr.index[0]
            nsep = arr.iloc[0]
            if arr.iloc[i] > minSumDMRs and numVS(arr.index[i]) >= minN:#0.5*nonsep:
                newref = arr.index[i]
                nsep = arr.iloc[i]
                ifSig = 1
                break
        if ifSig:
            finalList.append(newref)
            depthList.append(depth)
            weightList.append(nsep)
        else:
            return None
    else:
        arr = arr.groupby('sig.comparison.bin').sum().sort_values(ascending=False)
        for i in arr.index:
            newref = 0
            if arr[i] < minSumDMRs:
                return ([],[],[])
            if numVS(i)<minN:
                continue
            newref = i
            finalList.append(newref)
            depthList.append(depth)
            nsep = arr[i]
            weightList.append(nsep)
            break
        # print(newref,fulltree)
        if newref==0 and fulltree:
            # print('FT')
            for i in arr.index:
                newref = 0
                if arr[i] < minSumDMRs:
                    return ([],[],[])
                if numVS(i)==0:
                    continue
                newref = i
                finalList.append(newref)
                depthList.append(depth)
                nsep = arr[i]
                weightList.append(nsep)
                break
            # print(newref)

    if newref == 0:
        return (finalList, depthList, weightList)
        
    else:
        def mask(x, y, v):
            x = list(x)
            for i in range(len(y)):
                if y[i] != '|':
                    if y[i] != v:
                        x[i] = '0'
            return ''.join(x)
        masked = {}
        for i in ['1','2','3']:
            masked[i] = arr.copy()
            masked[i].index = pd.Series(arr.index).apply(lambda x:mask(x, newref, i))
            resRS = recurSplit(masked[i], mask(newref, newref, i), depth+1, nsep, minN, minSumDMRs, fulltree)
            finalList += resRS[0]
            depthList += resRS[1]
            weightList+= resRS[2]
            
        return (finalList, depthList, weightList)

    
def plotDMTree(cls, finalCls, reportPath, sids, cmap):
    k = len(cls[0])
    dmrcluster_m = []
    pd.Series(cls[0]).apply(lambda x:dmrcluster_m.append(x.split('|')))
    dmrcluster_m = pd.DataFrame(dmrcluster_m)
    dmrcluster_m = dmrcluster_m.astype(float)
    
    dmrcluster_m.columns = sids
    dmrcluster_m = dmrcluster_m.T
    dmrcluster_m = dmrcluster_m.sort_values(list(range(len(cls[2]))))
    # dmrcluster_m.sort_values(list(range(len(cls[2])))).to_csv(reportPath+'/clusters_detailed.tsv', sep='\t')
    
    clsD = finalCls['Group'].to_dict()
    ids = list(dmrcluster_m.sort_values(list(range(len(cls[2])))).index.map(clsD)+'_'+\
               dmrcluster_m.sort_values(list(range(len(cls[2])))).index+',')+['']
    treestr = ids.copy()
    
    for i in range(len(cls[2])):
        st = {}
        ed = {}
        for j in ['1','2','3']:
            st[j] = dmrcluster_m[i].astype(int).astype(str).sum().find(j)
            ed[j] = dmrcluster_m[i].astype(int).astype(str).sum().rfind(j)
            
            if st[j]!=-1:
                treestr[st[j]] = treestr[st[j]].split(ids[st[j]])[0]+'('+\
                                    ids[st[j]]+treestr[st[j]].split(ids[st[j]])[1]
                treestr[ed[j]] = treestr[ed[j]].split(ids[ed[j]])[0]+ids[ed[j]]+'):'+\
                                    str(cls[2][i])+','+treestr[ed[j]].split(ids[ed[j]])[1]
    
    f = open(reportPath+"DMTree.nwk", "w")
    f.write(''.join(treestr[:-1]).replace(',)',')')[:-1])
    f.close()

    oldkeys = list(cmap.keys())
    for i in oldkeys:
        cmap[clsD[i]+'_'+i] = cmap[i]
    
    from Bio import Phylo
    import matplotlib.pyplot as plt
    
    tree = Phylo.read(reportPath+"DMTree.nwk", "newick")

    f,a = plt.subplots(figsize=[20,len(sids)/5])
    Phylo.draw(tree, axes=a, do_show=False, label_colors=cmap)
    plt.savefig(reportPath+'DMTree.jpg', bbox_inches='tight')
    plt.savefig(reportPath+'DMTree.pdf', bbox_inches='tight')
    
    
def plotClustermap(mout, cls, reportPath, sids, finalCls, cls_full):
    import random
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    import numpy as np
    from matplotlib.patches import Patch
    
    k = len(cls[0])
    dmrcluster_m = []
    pd.Series(cls_full[0]).apply(lambda x:dmrcluster_m.append(x.split('|')))
    dmrcluster_m = pd.DataFrame(dmrcluster_m)
    dmrcluster_m = dmrcluster_m.astype(float)
    
    dmrcluster_m.columns = sids
    dmrcluster_m = dmrcluster_m.T
    dmrcluster_m = dmrcluster_m.sort_values(list(range(len(cls_full[2]))))
    
    dmrmean_m = [] # for pCA
    mout['mean'].apply(lambda x:dmrmean_m.append(x.split('|')))
    dmrmean_m = pd.DataFrame(dmrmean_m)
    dmrmean_m = dmrmean_m.astype(float).T
    dmrmean_m.index = sids
    pca = PCA(n_components=2)
    X = pd.DataFrame(pca.fit_transform(np.array(dmrmean_m)))
    X.index = dmrmean_m.index
    X.to_csv(reportPath+'PCA.tsv', sep='\t')
    pd.DataFrame(pca.explained_variance_ratio_).to_csv(reportPath+'PCA_ratio.tsv', sep='\t')
    # X['k'] = X[0]-X[0].min()
    # X['k'] = X['k']/X['k'].max()
    # dmrcluster_m[k] = dmrcluster_m.index.map(X['k'])
    
    # def td2(a,b):
    #     d = 0.01*abs(a[-1]-b[-1])
    #     for i in range(len(a)-1):
    #         if a[i]!=b[i]:
    #             d += max(cls[1])+1 - cls[1][i]
    #             break
    #     return d
    
    def td2(a,b):
        d = 0
        for i in range(len(a)):
            if a[i]!=b[i]:
                d += max(cls_full[1])+1 - cls_full[1][i]
                break
        return d
    
    # print(dmrcluster_m )
    # cm = sns.clustermap(dmrcluster_m, metric=td2,\
    #                     figsize=[4,4],row_cluster=True,col_cluster=False,\
    #                      dendrogram_ratio=0.2, colors_ratio=0.15, xticklabels=False, yticklabels=False, \
    #                     method='complete', cmap='Spectral_r')
    # lk = cm.dendrogram_row
    # print(lk)
    # print(cls)
    dmrmean_m = []
    
    def editD(a, b):
        a = a.split('|')
        b = b.split('|')
        num = 0
        all = 0
        for i in range(len(a)):
            if a[i] != '0':
                if a[i] != b[i]:
                    num += 1
                all += 1
        if all == 0:
            return 0
        return num/all
        
    denovo_pn2 = mout[['mean','sig.comparison','sig.comparison.bin']]
    denovo_pn2['tmp'] = 0
    c = np.int64(2**61)
    def allrelated(x):
        if x.find('1') > -1:
            return [x,\
                    x.replace('2','3').replace('1','2'),\
                    x.replace('1','3'),\
                    x.replace('1','x').replace('2','1').replace('x','2'),\
                    ]
        else:
            return [x,\
                    x.replace('2','1').replace('3','2'),\
                    x.replace('3','1'),\
                    x.replace('3','x').replace('2','3').replace('x','2'),\
                    ]
    for i in cls[0]:
        for j in allrelated(i):
            denovo_pn2['tmp'] += c*(denovo_pn2['sig.comparison.bin'].apply(lambda x:editD(j,x))==0)
            c = c/2
    denovo_filtered = denovo_pn2.loc[denovo_pn2['tmp']>0].sort_values('tmp', ascending=False)
    denovo_filtered['mean'].apply(lambda x:dmrmean_m.append(x.split('|')))
    # print(denovo_filtered)
    
    dmrmean_m = pd.DataFrame(dmrmean_m)
    dmrmean_m = dmrmean_m.astype(float)
    dmrmean_m.columns = sids
    dmrmean_m = dmrmean_m[dmrcluster_m.index].T
    
    clsD = finalCls['Group']
    random.seed(0)
    clsCD = {}
    for i in clsD.unique():
        clsCD[i] = (random.randint(1,100)/100,random.randint(1,100)/100,random.randint(1,100)/100)
    cmap = {}
    for i in dmrmean_m.index:
        cmap[i] = clsCD[clsD.to_dict()[i]]

    dmrmean_m_rename = dmrmean_m.loc[dmrcluster_m.index].copy()
    dmrmean_m_rename.index = dmrmean_m_rename.index.map(clsD)+' '+dmrmean_m_rename.index
    dmrmean_m_rename.to_csv(reportPath+'heatmap.tsv', sep='\t')
    dmrmean_m_rename_anno = pd.DataFrame(list(denovo_filtered['sig.comparison']),dmrmean_m_rename.columns)
    dmrmean_m_rename_anno.to_csv(reportPath+'heatmap.anno.tsv', sep='\t')
    if dmrmean_m.shape[1]>1:
        cm = sns.clustermap(dmrmean_m_rename,\
            row_colors=[dmrmean_m.index.map(cmap),\
                        (dmrmean_m[0]=='NO').map({False:'white'})],\
            # row_linkage=lk.linkage,\
            col_cluster=False,row_cluster=False,\
            cmap='Spectral_r', figsize=[0.2*len(sids)+0.1*max([len(i) for i in sids]),0.2*len(sids)], dendrogram_ratio=0.000001, xticklabels=False, yticklabels=True, \
            method='ward', cbar_pos=None, vmax=1, vmin=0, center=0.5, colors_ratio=0.02)
        plt.savefig(reportPath+'heatmap.jpg', bbox_inches='tight')
        plt.savefig(reportPath+'heatmap.pdf', bbox_inches='tight')
    else:
        cm = sns.clustermap(dmrmean_m_rename,\
            row_colors=[dmrmean_m.index.map(cmap),\
                        (dmrmean_m[0]=='NO').map({False:'white'})],\
            # row_linkage=lk.linkage,\
            col_cluster=False,row_cluster=False,\
            cmap='Spectral_r', figsize=[0.2*len(sids)+0.1*max([len(i) for i in sids]),0.2*len(sids)], dendrogram_ratio=0.000001, xticklabels=False, yticklabels=True, \
            method='ward', cbar_pos=None, vmax=1, vmin=0, center=0.5, colors_ratio=0.02)
        plt.savefig(reportPath+'heatmap.jpg', bbox_inches='tight')
        plt.savefig(reportPath+'heatmap.pdf', bbox_inches='tight')

    fig, ax0 = plt.subplots(figsize=(3, 3))
    X['grp'] = X.index.map(cmap)
    # print(X)
    sns.scatterplot(x=X[0],y=X[1],c=X['grp'],ax=ax0,s=30)
    
    label_color_dict = clsCD.copy()
    legend_handles = [Patch(color=color, label=label) for label, color in label_color_dict.items()]
    ax0.legend(handles=legend_handles, ncol=1, )
#    ax0.axis('off')
    plt.savefig(reportPath+'PCA.jpg', bbox_inches='tight')
    plt.savefig(reportPath+'PCA.pdf', bbox_inches='tight')

    return cmap


def clustering(mout, args):
    minN0 = args.minN0
    minN = args.minNSamples
    minSumDMRs = args.minSumDMRs
    mindiff_unsup = args.minMethDiffHigh

    def rename_cls_pn2(x):
        if x.count('3')>x.count('1'):
            x = x.replace('1','2')
        elif x.count('3')<=x.count('1'):
            x = x.replace('3','2')
        return x
    mout = mout.loc[(mout['#Hypo']>=minN0)\
                    &(mout['#Hyper']>=minN0)\
                    &(mout['meandiffabs']>mindiff_unsup)]
    mout['sig.comparison.bin'] = mout['sig.comparison'].apply(rename_cls_pn2)
    ranked = mout[['sig.comparison.bin','meandiffabs']].groupby('sig.comparison.bin').\
        sum()['meandiffabs'].sort_values(ascending=False)
    # print(ranked)
    cls = recurSplit(ranked.sort_values(ascending = False), minN=minN, minSumDMRs=minSumDMRs)
    # print(cls, minN, minSumDMRs)
    if cls is None:
        return (None, None)

    reportPath = args.output+'/'
    sids = [str(i) for i in pd.read_table(args.input, nrows=0).columns[2:]]
    
    finalCls = pd.DataFrame([i.split('|') for i in cls[0]]).sum()
    finalCls.index = sids
    rename_cls_id = {}
    j=0
    for i in sorted(finalCls.unique()):
        rename_cls_id[i] = 'G'+str(j)
        j+=1
    finalCls = finalCls.map(rename_cls_id)
    finalCls = pd.DataFrame(finalCls)
    finalCls.columns = ['Group']
    finalCls.index.name = 'ID'
    finalCls.to_csv(reportPath+'clusters.tsv', sep='\t')
    
    if args.visualization:
        # print(args.visualization)
        cls_full = recurSplit(ranked.sort_values(ascending = False), \
                              minN=minN, minSumDMRs=0, fulltree=True)
        cmap = plotClustermap(mout, cls, reportPath, sids, finalCls, cls_full)
        plotDMTree(cls_full, finalCls, reportPath, sids, cmap)
    
    return (finalCls, cls)



###################################################################################################
# GSEA
###################################################################################################
def DMRtable(args, finalCls, mout, unmout=None):
    tables = []
    
    dmrs_list = [mout,]
    if unmout is not None:
        dmrs_list.append(unmout)

    for dmrs in dmrs_list:
        if args.groupinfo:
            table = pd.DataFrame(mout['sig.comparison'].value_counts()[:10])
            table.columns = ['#DMRs']
        else:
            table = pd.DataFrame([dmrs['DMTree'].str.contains(('P'+i+',').replace('|',r'\|')).sum() for i in finalCls.columns[1:]], list(finalCls.columns[1:]))
            table.columns = ['#DMRs_hypo_in_left']
            table['#DMRs_hypo_in_right'] = [dmrs['DMTree'].str.contains(('N'+i+',').replace('|',r'\|')).sum() for i in table.index]
        # print(table)
        
        def decodeSigCmp(x):
            upmlist = {'1':[], '2':[], '3':[], '0':[]}
            x = x.split('|')
            grp_dict = pd.read_table(args.output+'/group-ID.tsv',\
                                    index_col='Group_ID')
            grp_dict = grp_dict['Group'].astype(str).to_dict()
            for i in range(len(x)):
                upmlist[x[i]].append(grp_dict[i])
            return [upmlist, 'Hypo:'+','.join(upmlist['1'])+' - vs - '+'Hyper:'+','.join(upmlist['3'])]
        
        def decodeSigCmpLR(x):
            upmlist = {'1':[], '2':[], '3':[], '0':[]}
            x = x.split('|')
            grp_dict = pd.read_table(args.output+'/group-ID.tsv',\
                                    index_col='Group_ID')
            grp_dict = grp_dict['Group'].astype(str).to_dict()
            for i in range(len(x)):
                upmlist[x[i]].append(grp_dict[i])
            
            if len(upmlist['3'])==0:
                return {'L':upmlist['1'], 'R':upmlist['2']}
            else:
                return {'L':upmlist['2'], 'R':upmlist['3']}
        
        if args.groupinfo:
            table['Hypo'] = [','.join(decodeSigCmp(i)[0]['1']) for i in table.index]
            table['Int'] = [','.join(decodeSigCmp(i)[0]['2']) for i in table.index]
            table['Hyper'] = [','.join(decodeSigCmp(i)[0]['3']) for i in table.index]
        else:
            table['Left_Child'] = [','.join(decodeSigCmpLR(i)['L']) for i in table.index]
            table['Right_Child'] = [','.join(decodeSigCmpLR(i)['R']) for i in table.index]
        
        table.index = range(len(table.index))#[decodeSigCmp(i) for i in table.index]
        tables.append(table)
    
    return tables

def gsea(args, finalCls, mout, unmout=None):
    import gseapy as gp
    
    gseapopup = ''
    tables = []
    
    dmrs_list = [mout,]
    if unmout is not None:
        dmrs_list.append(unmout)

    def decodeSigCmp(x):
        upmlist = {'1':[], '2':[], '3':[], '0':[]}
        x = x.split('|')
        grp_dict = pd.read_table(args.output+'/group-ID.tsv',\
                                index_col='Group_ID')
        grp_dict = grp_dict['Group'].astype(str).to_dict()
        for i in range(len(x)):
            upmlist[x[i]].append(grp_dict[i])
        return [upmlist, 'Hypo:'+','.join(upmlist['1'])+' - vs - '+'Hyper:'+','.join(upmlist['3'])]
    
    def decodeSigCmpLR(x):
        upmlist = {'1':[], '2':[], '3':[], '0':[]}
        x = x.split('|')
        grp_dict = pd.read_table(args.output+'/group-ID.tsv',\
                                index_col='Group_ID')
        grp_dict = grp_dict['Group'].astype(str).to_dict()
        for i in range(len(x)):
            upmlist[x[i]].append(grp_dict[i])
        
        if len(upmlist['3'])==0:
            return {'L':upmlist['1'], 'R':upmlist['2']}
        else:
            return {'L':upmlist['2'], 'R':upmlist['3']}
        
    if args.groupinfo:
        table = pd.DataFrame(mout['sig.comparison'].value_counts()[:10])
        table.columns = ['#DMRs']
        table['Hypo'] = [','.join(decodeSigCmp(i)[0]['1']) for i in table.index]
        table['Int'] = [','.join(decodeSigCmp(i)[0]['2']) for i in table.index]
        table['Hyper'] = [','.join(decodeSigCmp(i)[0]['3']) for i in table.index]

        if args.genesets and args.annotation:
            import gseapy as gp
            j = 0
            for i in table.index:
                gene_sets = args.genesets
                gene_list = list(set(mout.loc[(mout['sig.comparison']==i)&(mout['meandiffabs']>args.minMethDiffHigh)]['SYMBOL'].dropna()))
                
                try:
                    for gs in gene_sets.split(','):
                        enr = gp.enrichr(gene_list=gene_list,
                                    gene_sets=gs,
                                    organism='human',
                                    outdir=args.output+'/GSEA/'+i.replace('|','_'),
                                    cutoff = 1,
                                    format = 'jpg',
                                    )
                except:
                    pass
                    # print("GSEA error:",gene_list)
                            
                fig_path = './GSEA/'+i.replace('|','_')+\
                            "/"+args.genesets.split(',')[0].split('/')[-1]+".human.enrichr.reports.jpg"
                            
                gseapopup += "<div id=\"popupgsea"+str(j)+"\" class=\"popup\"><br>\
                    <button onclick=\"hidePopup('popupgsea"+str(j)+"')\">Close</button>\
                        <p>GSEA"+"</p><img src="+fig_path+" height=\"200\"><br></div>\n"
                j+=1
                
            table['GSEA'] = ["<button onclick=\"showPopup('popupgsea"+str(i)+\
                            "')\">Click to show GSEA results</button>" \
                            for i in range(table.shape[0])]
        
        table.index = range(len(table.index))#[decodeSigCmp(i) for i in table.index]
        tables.append(table)
    else:
        uors = 'sup'
        for dmrs in dmrs_list:
            if dmrs is not None:
                table = pd.DataFrame([dmrs['DMTree'].str.contains(('P'+i+',').replace('|',r'\|')).sum() for i in finalCls.columns[1:]], list(finalCls.columns[1:]))
                table.columns = ['#DMRs_hypo_in_left']
                table['#DMRs_hypo_in_right'] = [dmrs['DMTree'].str.contains(('N'+i+',').replace('|',r'\|')).sum() for i in table.index]
                # print(table)
                
                table['left'] = [','.join(decodeSigCmpLR(i)['L']) for i in table.index]
                table['right'] = [','.join(decodeSigCmpLR(i)['R']) for i in table.index]
    
                if args.genesets and args.annotation:
                    import gseapy as gp
                    j = 0
                    for i in table.index:
                        gene_sets = args.genesets
                        gene_list = list(set(dmrs.loc[(dmrs['DMTree'].str.contains(('P'+i+',').replace('|',r'\|')))&(dmrs['meandiffabs']>args.minMethDiffHigh)]['SYMBOL'].dropna()))
                        
                        try:
                            for gs in gene_sets.split(','):
                                enr = gp.enrichr(gene_list=gene_list,
                                            gene_sets=gs,
                                            organism='human',
                                            outdir=args.output+'/GSEA/'+'P'+uors+i.replace('|','_'),
                                            cutoff = 1,
                                            format = 'jpg',
                                            )
                        except:
                            pass
                            # print("GSEA error:",gene_list)
                                    
                        fig_path = './GSEA/'+'P'+uors+i.replace('|','_')+\
                                    "/"+args.genesets.split(',')[0].split('/')[-1]+".human.enrichr.reports.jpg"
                                    
                        gseapopup += "<div id=\"popupgsea"+'P'+uors+str(j)+"\" class=\"popup\"><br>\
                            <button onclick=\"hidePopup('popupgsea"+'P'+uors+str(j)+"')\">Close</button>\
                                <p>GSEA for hypo in "+','.join(decodeSigCmpLR(i)['L'])+"</p><img src="+fig_path+" height=\"200\"><br></div>\n"
                        j+=1
                        
                    table['GSEA_hypo_in_left'] = ["<button onclick=\"showPopup('popupgsea"+'P'+uors+str(i)+\
                                    "')\">Click to show GSEA results</button>" \
                                    for i in range(table.shape[0])]
                    
                    j = 0
                    for i in table.index:
                        gene_sets = args.genesets
                        gene_list = list(set(dmrs.loc[(dmrs['DMTree'].str.contains(('N'+i+',').replace('|',r'\|')))&(dmrs['meandiffabs']>args.minMethDiffHigh)]['SYMBOL'].dropna()))
                        
                        try:
                            for gs in gene_sets.split(','):
                                enr = gp.enrichr(gene_list=gene_list,
                                            gene_sets=gs,
                                            organism='human',
                                            outdir=args.output+'/GSEA/'+'N'+uors+i.replace('|','_'),
                                            cutoff = 1,
                                            format = 'jpg',
                                            )
                        except:
                            pass
                            # print("GSEA error:",gene_list)
                                    
                        fig_path = './GSEA/'+'N'+uors+i.replace('|','_')+\
                                    "/"+args.genesets.split(',')[0].split('/')[-1]+".human.enrichr.reports.jpg"
                                    
                        gseapopup += "<div id=\"popupgsea"+'N'+uors+str(j)+"\" class=\"popup\"><br>\
                            <button onclick=\"hidePopup('popupgsea"+'N'+uors+str(j)+"')\">Close</button>\
                                <p>GSEA for hypo in "+','.join(decodeSigCmpLR(i)['R'])+"</p><img src="+fig_path+" height=\"200\"><br></div>\n"
                        j+=1
                        
                    table['GSEA_hypo_in_right'] = ["<button onclick=\"showPopup('popupgsea"+'N'+uors+str(i)+\
                                    "')\">Click to show GSEA results</button>" \
                                    for i in range(table.shape[0])]
                
                table.index = range(len(table.index))#[decodeSigCmp(i) for i in table.index]
                tables.append(table)
            else:
                tables.append(None)
            uors = 'unsup'
    
    return (gseapopup, tables)



###################################################################################################
# HTML report
###################################################################################################
def report_unsup(args, start_time, end_time, unmout, finalCls, mout):
    with open(os.path.realpath(__file__).replace('metilene3.py','')+'template_unsup.html', 'r') as template_file:
        template_content = template_file.read()

    final_html = template_content.replace('<h2>Metilene Report for XXX</h2>', '<h2>Metilene Report for '+args.input.split('/')[-1]+'</h2>')
    final_html = final_html.replace('<div>Version: XXX</div><br>', '<div>Version: '+str(VERSION)+'</div><br>')
    final_html = final_html.replace('<div>Command: XXX</div><br>', '<div>Command: '+''.join([i+' ' for i in sys.argv])+'</div><br>')
    final_html = final_html.replace('<div>Parameters: XXX</div><br>', '<div>Parameters: <br>'+str(args).split('Namespace')[-1][1:-1].split(', skipMetilene')[0]+'</div><br>')
    final_html = final_html.replace('<div>Start time: XXX</div><br>', '<div>Start time: '+str(start_time)+'</div>')
    final_html = final_html.replace('<div>End time: XXX</div><br>', 'End time: '+str(end_time)+'</div><br>')

    final_html = final_html.replace('<div>Number of unsupervised DMRs: XXX</div><br>', 'Number of unsupervised DMRs: '+str(unmout.shape[0])+'</div><br>')
    
    final_html = final_html.replace('<div>Number of clusters: XXX</div><br>', 'Number of clusters: '+str(len(finalCls['Group'].unique()))+'</div><br>')
    cls_table_html = finalCls.to_html(escape=False)
    final_html = final_html.replace('<div id="pandas_table_placeholder_cluster"></div>', cls_table_html)

    if not args.visualization:
        final_html = final_html.replace('<button onclick="showPopup(\'popupTree\')">Click to show the figures</button>', '')

    final_html = final_html.replace('<div>Number of supervised DMRs: XXX</div><br>', 'Number of supervised DMRs: '+str(mout.shape[0])+'</div><br>')
    
    if args.genesets and args.annotation:
        gseapopup, tables = gsea(args, finalCls, mout, unmout)
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_sup"></div>', tables[0].to_html(escape=False))
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_unsup"></div>', tables[1].to_html(escape=False))
        final_html = final_html.replace('<div id="gsea_placeholder"></div>', gseapopup)
    if not args.genesets:
        tables = DMRtable(args, finalCls, mout, unmout)
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_sup"></div>', tables[0].to_html(escape=False))
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_unsup"></div>', tables[1].to_html(escape=False))

    with open(args.output+'/report.html', 'w') as final_file:
        final_file.write(final_html)


def report_sup(args, start_time, end_time, mout):
    with open(os.path.realpath(__file__).replace('metilene3.py','')+'template_sup.html', 'r') as template_file:
        template_content = template_file.read()

    final_html = template_content.replace('<h2>Metilene Report for XXX</h2>', '<h2>Metilene Report for '+args.input.split('/')[-1]+'</h2>')
    final_html = final_html.replace('<div>Version: XXX</div><br>', '<div>Version: '+str(VERSION)+'</div><br>')
    final_html = final_html.replace('<div>Command: XXX</div><br>', '<div>Command: '+''.join([i+' ' for i in sys.argv])+'</div><br>')
    final_html = final_html.replace('<div>Parameters: XXX</div><br>', '<div>Parameters: <br>'+str(args).split('Namespace')[-1][1:-1].split(', skipMetilene')[0]+'</div><br>')
    final_html = final_html.replace('<div>Start time: XXX</div><br>', '<div>Start time: '+str(start_time)+'</div>')
    final_html = final_html.replace('<div>End time: XXX</div><br>', 'End time: '+str(end_time)+'</div><br>')

    final_html = final_html.replace('<div>Number of supervised DMRs: XXX</div><br>', 'Number of supervised DMRs: '+str(mout.shape[0])+'</div><br>')
    finalCls = pd.read_table(args.groupinfo, index_col='ID')[['Group']]
    finalCls['Group'] = finalCls['Group'].astype(str)
    if args.genesets and args.annotation:
        gseapopup, tables = gsea(args, finalCls, mout)
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_sup"></div>', tables[0].to_html(escape=False))
        final_html = final_html.replace('<div id="gsea_placeholder"></div>', gseapopup)
    else:
        tables = DMRtable(args, finalCls, mout)
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_sup"></div>', tables[0].to_html(escape=False))
    with open(args.output+'/report.html', 'w') as final_file:
        final_file.write(final_html)

def report_wosup(args, start_time, end_time, unmout, finalCls):
    with open(os.path.realpath(__file__).replace('metilene3.py','')+'template_wosup.html', 'r') as template_file:
        template_content = template_file.read()

    final_html = template_content.replace('<h2>Metilene Report for XXX</h2>', '<h2>Metilene Report for '+args.input.split('/')[-1]+'</h2>')
    final_html = final_html.replace('<div>Version: XXX</div><br>', '<div>Version: '+str(VERSION)+'</div><br>')
    final_html = final_html.replace('<div>Command: XXX</div><br>', '<div>Command: '+''.join([i+' ' for i in sys.argv])+'</div><br>')
    final_html = final_html.replace('<div>Parameters: XXX</div><br>', '<div>Parameters: <br>'+str(args).split('Namespace')[-1][1:-1].split(', skipMetilene')[0]+'</div><br>')
    final_html = final_html.replace('<div>Start time: XXX</div><br>', '<div>Start time: '+str(start_time)+'</div>')
    final_html = final_html.replace('<div>End time: XXX</div><br>', 'End time: '+str(end_time)+'</div><br>')

    final_html = final_html.replace('<div>Number of unsupervised DMRs: XXX</div><br>', 'Number of unsupervised DMRs: '+str(unmout.shape[0])+'</div><br>')
    
    final_html = final_html.replace('<div>Number of clusters: XXX</div><br>', 'Number of clusters: '+str(len(finalCls['Group'].unique()))+'</div><br>')
    cls_table_html = finalCls.to_html(escape=False)
    final_html = final_html.replace('<div id="pandas_table_placeholder_cluster"></div>', cls_table_html)

    if not args.visualization:
        final_html = final_html.replace('<button onclick="showPopup(\'popupTree\')">Click to show the figures</button>', '')

    if args.genesets and args.annotation:
        gseapopup, tables = gsea(args, finalCls, None, unmout)
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_unsup"></div>', tables[1].to_html(escape=False))
        final_html = final_html.replace('<div id="gsea_placeholder"></div>', gseapopup)
    if not args.genesets:
        tables = DMRtable(args, finalCls, unmout)
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_unsup"></div>', tables[0].to_html(escape=False))

    with open(args.output+'/report.html', 'w') as final_file:
        final_file.write(final_html)

def report_nocls(args, start_time, end_time, unmout):
    with open(os.path.realpath(__file__).replace('metilene3.py','')+'template_wosup.html', 'r') as template_file:
        template_content = template_file.read()

    final_html = template_content.replace('<h2>Metilene Report for XXX</h2>', '<h2>Metilene Report for '+args.input.split('/')[-1]+'</h2>')
    final_html = final_html.replace('<div>Version: XXX</div><br>', '<div>Version: '+str(VERSION)+'</div><br>')
    final_html = final_html.replace('<div>Command: XXX</div><br>', '<div>Command: '+''.join([i+' ' for i in sys.argv])+'</div><br>')
    final_html = final_html.replace('<div>Parameters: XXX</div><br>', '<div>Parameters: <br>'+str(args).split('Namespace')[-1][1:-1].split(', skipMetilene')[0]+'</div><br>')
    final_html = final_html.replace('<div>Start time: XXX</div><br>', '<div>Start time: '+str(start_time)+'</div>')
    final_html = final_html.replace('<div>End time: XXX</div><br>', 'End time: '+str(end_time)+'</div><br>')

    final_html = final_html.replace('<div>Number of unsupervised DMRs: XXX</div><br>', 'Number of unsupervised DMRs: '+str(unmout.shape[0])+'</div><br>')
    
    final_html = final_html.replace('<div>Number of clusters: XXX</div><br>', 'No clusters found. </div><br>')
    final_html = final_html.replace('<div id="pandas_table_placeholder_cluster"></div>', '')

    final_html = final_html.replace('<button onclick="showPopup(\'popupTree\')">Click to show the figures</button>', '')
    final_html = final_html.replace('<button onclick="showPopup(\'popupCluster\')">Click to show the Table of clusters</button>', '')

    finalCls = pd.read_table(args.input, nrows=0).T[2:]
    finalCls['Group'] = finalCls.index
    finalCls['Group'] = finalCls['Group'].astype(str)
    if args.genesets and args.annotation:
        gseapopup, tables = gsea(args, finalCls, unmout)
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_unsup"></div>', tables[0].to_html(escape=False))
        final_html = final_html.replace('<div id="gsea_placeholder"></div>', gseapopup)
    if not args.genesets:
        tables = DMRtable(args, finalCls, unmout)
        final_html = final_html.replace('<div id="pandas_table_placeholder_dmr_unsup"></div>', tables[0].to_html(escape=False))

    with open(args.output+'/report.html', 'w') as final_file:
        final_file.write(final_html)

###################################################################################################
# main
###################################################################################################
def checkParams(args):
    msg = None
    if not (args.input and args.output):
        msg = 'ERROR: please provide both input filename and output folder.'
    
    if args.genesets and (not args.annotation):
        msg = 'ERROR: please also provide annotation if you want to run GSEA.'
        
    if args.groupinfo and args.visualization:
        msg = 'ERROR: visualization function is only for unsupervised mode.'

    if args.groupinfo:
        try:
            pd.read_table(args.groupinfo, index_col='ID')['Group'].astype(str)
        except:
            msg = 'ERROR: please check the format of the table of group information and provide a tab-separated tsv file.'

    try:
        pd.read_table(args.input, nrows=1)
    except:
        msg = 'ERROR: please check the format of the input matrix and provide a tab-separated tsv file.'
         
    return msg

def main():
    start_time = time.ctime()
    args = parser.parse_args()
    print(start_time,": Started.")
    # print(args)
    msg = checkParams(args)
    if msg:
        print(msg)
        return None

    if not os.path.isfile(os.path.realpath(__file__).replace('metilene3.py','metilene')):
        getMetilene()
    try:
        os.mkdir(args.output)
    except:
        pass
        
    if args.groupinfo:
        print(time.ctime(),": Running supervised mode...")
        headerfile = args.output+'/'+args.input.split('/')[-1]+'.header'
        preprocess(args, headerfile, 'sup', \
                   args.groupinfo)
        runMetilene(args, headerfile, 'sup')
        mout = processOutput(args, 'sup', anno='T')
        end_time = time.ctime()
        if mout is None:
            print(end_time,": Finished.")
            return
        report_sup(args, start_time, end_time, mout)
        print(end_time,": Finished.")
    else:
        print(time.ctime(),": Running unsupervised mode...")
        headerfile = args.output+'/'+args.input.split('/')[-1]+'.unsup.header'
        preprocess(args, headerfile, 'unsup')
        runMetilene(args, headerfile, 'unsup')
        unmout = processOutput(args, 'unsup', anno='T')
        if unmout is None:
            end_time = time.ctime()
            print(end_time,": Finished.")
            return
        print(time.ctime(),": Clustering...")
        
        finalCls, cls = clustering(unmout, args)
        if finalCls is None:
            print('Warning: No cluster found. Please check the data or use smaller meandiff for clustering.')
            end_time = time.ctime()
            headerfile = args.output+'/'+args.input.split('/')[-1]+'.header'
            finalCls = pd.read_table(args.input, nrows=0).T[2:]
            finalCls['Group'] = finalCls.index
            finalCls['Group_ID'] = range(len(finalCls.index))
            finalCls.to_csv(args.output+'/group-ID.tsv', sep='\t', index=False)
            args.groupinfo = 1
            report_nocls(args, start_time, end_time, unmout)
            os.system("rm "+args.output+"/group-ID.tsv")
            print(end_time,": Finished.")
            return None

        unmout = addDMTree2DMR(args, 'unsup', cls, finalCls)
        
        if args.withSupervised==False:
            end_time = time.ctime()
            headerfile = args.output+'/'+args.input.split('/')[-1]+'.header'
            preprocess(args, headerfile, 'sup', \
                       args.output+'/clusters.tsv')
            report_wosup(args, start_time, end_time, unmout, finalCls.drop(columns=cls[0]))
            print(end_time,": Finished.")
            return None

        print(time.ctime(),": Running supervised mode...")

        headerfile = args.output+'/'+args.input.split('/')[-1]+'.header'
        preprocess(args, headerfile, 'sup', \
                   args.output+'/clusters.tsv')
        runMetilene(args, headerfile, 'sup')
        mout = processOutput(args, 'sup', anno='T')
        if mout is None:
            end_time = time.ctime()
            print(end_time,": Finished.")
            return
        mout = addDMTree2DMR(args, 'sup', cls, finalCls)

        end_time = time.ctime()
        report_unsup(args, start_time, end_time, unmout, finalCls.drop(columns=cls[0]), mout)
        print(end_time,": Finished.")

    if not args.keeptmp:
        if not args.skipMetilene:
            os.system("rm "+args.output+"/*.header")
        
main()

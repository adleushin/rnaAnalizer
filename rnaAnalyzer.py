
import argparse
import xlrd
import regex as re
import pandas as pd
import subprocess
from simplesam import Reader
import gffutils

par=argparse.ArgumentParser(description='This script predicts and draws primer dimers')
par.add_argument('--reads-file','-reads',
                    dest='readsFile',type=str,
                    help='XLSX-file with NGS read pairs. Use this argument only when you want to search primer dimers in this file',
                    required=True)
par.add_argument('--pcrOligs-file','-pcr',
                    dest='pcrOligsFile',type=str,
                    help='input TXT-file with PCR oligs',
                    required=True)
par.add_argument('--output-file','-out',dest='outFile',
                    type=str,help='XLS-file for output',
                    required=True)
par.add_argument('--rtOligs-file','-rt',dest='rtOligsFile',
                    type=str,help='XLS-file for RT oligs',
                    required=True)
par.add_argument('--readSequences','-rs',dest='fastq',
                    type=str,help='Path to the output FASTQ-file with read sequences for mapping',
                    required=True)
par.add_argument('--indexedTranscriptome', '-inT', dest='indexT',
                  type=str,help='Index of transcriptome filename prefix',
                  required=True)
par.add_argument('--indexedGenome', '-inG', dest='indexG',
                  type=str,help='Index of genome filename prefix',
                  required=True)
par.add_argument('--T-sam-file', '-Tsam', dest='Tsam',
                  type=str,help='Path to the output sam file with reads mapped on transcriptome',
                  required=True)
par.add_argument('--G-sam-file', '-Gsam', dest='Gsam',
                  type=str,help='Path to the output sam file with reads mapped on Genome',
                  required=True)
par.add_argument('--transcriptome-GFF', '-tgff', dest='tgff',
                  type=str,help='Transcriptome GFF file',
                  required=True)
par.add_argument('--genome-GFF', '-ggff', dest='ggff',
                  type=str,help='Genome GFF file',
                  required=True)
par.add_argument('--sensetivity', '-sens', dest='sensitivity',
                  type=str,help='To increase sensetivity of sequence mapping: sensetive or very-sensetive.',
                  required=False, default=None)
par.add_argument('--nrows', '-n', dest='nrows',
                  type=str,help='First n rows of xlsx reads file what will be proceeded',
                  required=False, default='all')

args=par.parse_args()

def revComplement(seq):
    nucToRev={'A':'T','G':'C','T':'A','C':'G',
                'a':'t','g':'c','t':'a','c':'g'}
    seq_rc=[]
    for n in seq[::-1]:
        seq_rc.append(nucToRev[n])
    return(''.join(seq_rc))

def getVarForNonNumericalVariants(name,v):
    if name in v:
        if v==name:
            newV=1
        else:
            v=v[v.find(name)+len(name):]
            if '-' in v:
                v.replace('-', '')
            try:
                newV=int(v)
            except ValueError:
                newV=1
    return newV


def getVar(items):
    varDict={'IX':'9','VIII':'8','VII':'7','VI':'6','IV':'5','V':'5','III':'3',\
             'II':'2','I':'1','hBACHd':'4','hBACHc':'3','hBACHb':'2','hBACHa':'1','P5CDhS':'1','P5CDhL':'1',\
                '1a':'1', 'WNT-2B1':'1','WNT-2B2':'2'}
    n=None
    for i in range (len(items)):
        if items[i][0]=='product':
            n=i
    if n!=None:
        item=items[n][1][0]
        if 'transcript variant ' in item:
            v=item[item.find('transcript variant ')+19:]
            if 'X' in v:
                v=10
            try:
                v=int(v)
            except ValueError:
                try:
                    v=ord(v)
                except TypeError:
                    try:
                        v=varDict[v]
                        v=int(v)
                    except KeyError:
                        if 'Tpm' in v:
                            if v=='Tpm3':
                                v=1
                            else:
                                v=v.split('.')[1]
                                v=int(v)
                        elif 'alpha' in v:
                            v=getVarForNonNumericalVariants('alpha',v)
                        elif 'beta' in v:
                            v=getVarForNonNumericalVariants('beta',v)
                        elif 'gamma' in v:
                            v=getVarForNonNumericalVariants('gamma',v)
                        elif ',' in v:
                            v=v.split(',')[0]
                            try:
                                v=int(v)
                            except ValueError:
                                v=1
                        else:
                            try:
                                v=int(v)
                            except ValueError:
                                v=1
        else:
            v=1
    else:
        v=1
    return(v)


def makeResList(res, list):
    if res!=None:
        list.append(res)
    else:
        list.append(None)

def makeRes(name, res):
    if res==None:
        res=str(name)
        return(res)
    else:
       res=res+','+str(name)
       return(res)
    
def makeResGenes(name, res):
    if res==None:
        res=str(name)
        return(res)
    else:
       res=res+'|'+str(name)
       return(res)
    
def searchStringInString(substr,string,errors):
    if substr in string:
        return(True,0,string.find(substr))
    m=re.search('('+substr+'){e<='+str(errors)+'}',
                    string,flags=re.BESTMATCH)
    if m!=None:
        sumErrNum=sum(m.fuzzy_counts)
        return(True,sumErrNum,m.span()[0])
    else:
        return(False,0,-1)
    
def getTGenesDict(samFile):
    tGenesDict={}
    Tsam_file=open(samFile)
    in_Tsam=Reader(Tsam_file)
    for x in in_Tsam:
        if x.rname!='*' :
            if int(x.qname) not in tGenesDict.keys():
                tGenesDict[int(x.qname)]=[]
            tGenesDict[int(x.qname)].append((x.rname, list(x.coords)[:1][0], list(x.coords)[-1:][0]))
    Tsam_file.close()
    return(tGenesDict)

def getGGenesDict(samFile, tGenesDict):
    gGenesDict={}
    Gsam_file=open(samFile)
    in_Gsam=Reader(Gsam_file)
    for x in in_Gsam:
        if int(x.qname) not in tGenesDict.keys() and x.rname!='*' and len(x.rname)<6:
            if int(x.qname) not in gGenesDict.keys():
                gGenesDict[int(x.qname)]=[]
            gGenesDict[int(x.qname)].append((x.rname, list(x.coords)[:1][0], list(x.coords)[-1:][0]))
    Gsam_file.close()
    return(gGenesDict)

def getGenesAndCoords(nRows, tGenesDict, gGenesDict, allExonsTranscriptome, allExonsGenome):
    genes=[]
    cords=[]
    for i in range (1,nRows):
        resGenes=None
        resCords=None
        if i in tGenesDict.keys():
            for mapped in tGenesDict[i]:
                gene=None
                allResGForCoord=None
                resGenDict={}
                try:
                    gene=allExonsTranscriptome[mapped[0]]
                except KeyError:
                    if mapped[0] not in resGenDict:
                        resGenDict[mapped[0]]=[]
                    resGenDict[mapped[0]].append(' unknown') 
                if gene!=None:   
                    for exon in allExonsTranscriptome[mapped[0]]:
                        exonFounded=False
                        for x in range(1,3):
                            if mapped[x] in range(int(allExonsTranscriptome[mapped[0]][exon][0]), int(allExonsTranscriptome[mapped[0]][exon][1])) and exonFounded==False:
                                if mapped[0] not in resGenDict:
                                    resGenDict[mapped[0]]=[]
                                if exon.split('-')[-1] not in resGenDict[mapped[0]]:
                                    resGenDict[mapped[0]].append(exon.split('-')[-1])
                for key in resGenDict.keys():
                    resGForCoord=key+': ex'
                    for n in resGenDict[key]:
                        resGForCoord+=n+','
                    resGForCoord=resGForCoord[:-1]
                    if allResGForCoord==None:
                        allResGForCoord=resGForCoord
                    else:
                        allResGForCoord+=' '+resGForCoord
                resGenes=makeResGenes(allResGForCoord, resGenes)                        
                resCords=makeRes(mapped[0]+' '+str(mapped[1])+'-'+str(mapped[2]), resCords)
        elif i in gGenesDict.keys():
            for mapped in gGenesDict[i]:
                resGenDict={}
                startReadExonFounded=False
                endReadExonFounded=False
                for g in allExonsGenome[mapped[0]]:
                    for exon in allExonsGenome[mapped[0]][g]:
                        for x in range(1,3):
                            if mapped[x] in range(int(allExonsGenome[mapped[0]][g][exon][0]), int(allExonsGenome[mapped[0]][g][exon][1])):
                                if g not in resGenDict:
                                    resGenDict[g]=[]
                                try:
                                    num=int(exon.split('.')[1].split('-')[1])
                                except IndexError:
                                    try:
                                        num=int(exon.split('-')[-1])
                                    except IndexError:
                                        print('cant extract intron number'+exon)
                                if '.' in exon:
                                    num=int(exon.split('.')[1].split('-')[1])
                                else:
                                    num=int(exon.split('-')[-1])
                                if num not in resGenDict[g]:
                                    resGenDict[g].append(num)
                                if x==1:
                                    startReadExonFounded=True
                                else:
                                    endReadExonFounded=True 
                                exonFounded=True
                if startReadExonFounded==False:
                    for g in allExonsGenome[mapped[0]]:
                        if int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[0]][0])<int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[-1]][1]):
                            start=int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[0]][0])
                            end=int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[-1]][1])
                        else:
                            start=int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[-1]][1])
                            end=int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[0]][0])
                        if  mapped[1] in range(start, end):                         
                            if g not in resGenDict:
                                resGenDict[g]=[]
                            if 'intron' not in resGenDict[g]:
                                resGenDict[g].append('intron')
                            startReadExonFounded=True
                if endReadExonFounded==False:
                    for g in allExonsGenome[mapped[0]]:
                        if int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[0]][0])<int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[-1]][1]):
                            start=int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[0]][0])
                            end=int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[-1]][1])
                        else:
                            start=int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[-1]][1])
                            end=int(allExonsGenome[mapped[0]][g][list(allExonsGenome[mapped[0]][g].keys())[0]][0])
                        if  mapped[2] in range(start, end): 
                            if g not in resGenDict:
                                if g not in resGenDict:
                                    resGenDict[g]=[]
                                if 'intron' not in resGenDict[g]:
                                    resGenDict[g].append('intron')
                            endReadExonFounded=True
                if endReadExonFounded==False or startReadExonFounded==False:
                    resGenDict[mapped[0]]='intergenic region'
                allResGForCoord=None
                for key in resGenDict.keys():
                    resGForCoord=None
                    if resGenDict[key]=='intergenic region':
                        resGForCoord=key+': intergenic region'
                    else:
                        if resGenDict[key]==['intron']:
                            resGForCoord=key+': intron'
                        else:
                            resGForCoord='ex'
                            for n in resGenDict[key]:
                                if n!='intron':
                                    resGForCoord+=str(n)+','
                            resGForCoord=resGForCoord[:-1]
                            if resGenDict[key][0]=='intron':
                                resGForCoord='intron'+resGForCoord
                            if resGenDict[key][0]=='intron':
                                resGForCoord+=' intron'
                            resGForCoord=key+': '+resGForCoord
                    if allResGForCoord==None:
                        allResGForCoord=resGForCoord
                    else:
                        allResGForCoord+=' '+resGForCoord
                resGenes=makeResGenes(allResGForCoord, resGenes)
                resCords=makeRes(mapped[0]+' '+str(mapped[1])+'-'+str(mapped[2]), resCords)
        else:
            resCords=None
        makeResList(resGenes, genes)
        makeResList(resCords, cords)
    return(genes, cords)

pcrNucs4s={}
wb=xlrd.open_workbook(args.pcrOligsFile)
ws=wb.sheet_by_index(0)
for i in range(ws.nrows):
    row=ws.row_values(i)
    if row[1][:4] not in pcrNucs4s.keys():
        pcrNucs4s[row[1][:4]]={} 
    pcrNucs4s[row[1][:4]][row[1]]=row[0]

rtOligs={}
wb=xlrd.open_workbook(args.rtOligsFile)
ws=wb.sheet_by_index(0)
for i in range(ws.nrows):
    row=ws.row_values(i)
    rtOligs[row[1]]=[row[1][:6],row[1][(len(row[0])-6):],row[0]]


reads1=[]
reads2=[]
pcr1=[]
pcr2=[]
mism=[]
rt=[]
quantity=[]
dest=[]
pcrLenList=[]
toFastaNum=[]
toFastaSeq1=[]
toFastaSeq2=[]
wb=xlrd.open_workbook(args.readsFile)
ws=wb.sheet_by_index(0)
if args.nrows=='all':
    nRows=ws.nrows
else:
    nRows=int(args.nrows)
for i in range(nRows):
    if i==0:
        continue
    row=ws.row_values(i)
    quantity.append(row[2])
    pcrOlig=None
    pcrOligsNotInTheStart={}
    pcrOligsNotInTheStartMismatches=None
    foundedRts={}
    rtN=None
    resPCR1=None
    resPCR2=None
    resMism=None
    resRt=None
    resDest=None
    pcrLen=[]
    revRtOlig=None
    foundedRt=None
    for n in range(2):
        destPcr=None
        destRt=None
        if n==0:
            reads1.append(row[n])
        else:
            reads2.append(row[n])
        if row[n][:4] in pcrNucs4s.keys():
            for olig in pcrNucs4s[row[n][:4]]:
                findOlig=searchStringInString(olig,row[n],int(len(olig)*0.2))
                if findOlig[0]==True:
                    if findOlig[2]==0:
                        if pcrOlig==None:
                            pcrOlig=olig
                            pcrOligRow=n
                            oligMismatch=findOlig[1]
                        else:
                            if findOlig[1]<oligMismatch:
                                pcrOlig=olig
                                pcrOligRow=n
                                oligMismatch=findOlig[1]
                    else:
                        if pcrOligsNotInTheStart=={}:
                            pcrOligsNotInTheStart[olig]=[findOlig[2],n]
                            pcrOligsNotInTheStartMismatches=findOlig[1]
                        else:
                            if findOlig[1]>pcrOligsNotInTheStartMismatches:
                                pcrOligsNotInTheStart={}
                                pcrOligsNotInTheStart[olig]=[findOlig[2],n]
                                pcrOligsNotInTheStartMismatches=findOlig[1]
                            elif findOlig[1]==pcrOligsNotInTheStartMismatches:
                                pcrOligsNotInTheStart[olig]=[findOlig[2],n]
        for rtOlig in rtOligs.keys():
            for x in range(2):
                if (rtOligs[rtOlig][x]) in row[n]:
                    findRtOlig=searchStringInString(rtOlig,row[n],int(len(rtOlig)*0.2))
                    if findRtOlig[0]==True:
                        foundedRts[findRtOlig[1]]=(rtOligs[rtOlig][2],findRtOlig[2],n,rtOlig)
    if pcrOlig==None:
        for nucs4 in pcrNucs4s:
            for olig in pcrNucs4s[nucs4]:
                findOlig=searchStringInString(olig,row[1],int(len(olig)*0.2))
                if findOlig[0]==True:
                    if findOlig[2]==0:
                        if pcrOlig==None:
                            pcrOlig=olig
                            pcrOligRow=n
                            oligMismatch=findOlig[1]
                        else:
                            print('Founded two pcrOligs')
                        break
                    else:
                        if pcrOligsNotInTheStart=={}:
                            pcrOligsNotInTheStart[olig]=[findOlig[2],n]
                            pcrOligsNotInTheStartMismatches=findOlig[1]
                        else:
                            if findOlig[1]>pcrOligsNotInTheStartMismatches:
                                pcrOligsNotInTheStart={}
                                pcrOligsNotInTheStart[olig]=[findOlig[2],n]
                                pcrOligsNotInTheStartMismatches=findOlig[1]
                            elif findOlig[1]==pcrOligsNotInTheStartMismatches:
                                pcrOligsNotInTheStart[olig]=[findOlig[2],n]
    if foundedRts!={}:
        foundedRts=sorted(foundedRts.items(), key=lambda item:item[0])[0][1]
        foundedRt=foundedRts[0]
        destRt=foundedRts[1]
        revRtOlig=foundedRts[3]
        resRt=makeRes(foundedRt,resRt)
    if pcrOlig!=None:
        destPcr=len(olig)
        if pcrOligRow==0:
            resPCR1=makeRes(pcrNucs4s[pcrOlig[:4]][pcrOlig], resPCR1)
            toFastaNum.append(str(i))
            toFastaSeq1.append(row[0][len(pcrOlig):][:40])
            toFastaSeq2.append(row[1][:40])
        if pcrOligRow==1:
            resPCR2=makeRes(pcrNucs4s[pcrOlig[:4]][pcrOlig], resPCR2)
            toFastaNum.append(str(i))
            toFastaSeq2.append(row[1][len(pcrOlig):][:40])
            toFastaSeq1.append(row[0][:40])
        resMism=makeRes(oligMismatch,resMism)
    if destPcr!=None and destRt!=None:
        resDest=(destRt-destPcr)
    if revRtOlig!=None and pcrOlig!=None:
        if foundedRts[2]==0:
            rowN=1
        else:
            rowN=0
        foundRevRtOlig=searchStringInString(revComplement(revRtOlig),row[rowN],int(len(revRtOlig)*0.2))
        if foundRevRtOlig[0]==True:
            resDest=('contains revcompl rt olig , '+str(foundRevRtOlig[2]-destPcr))
    if n==1 and pcrOlig==None and foundedRt!=None:
        for i in range(2):
            if searchStringInString(row[n][:6],foundedRt,1)[0]==True:
                if i==0:
                    resPCR1=makeRes('RT olig as PCR olig , '+rtOligs[foundedRt][2],resPCR1)
                if i==1:
                    resPCR2=makeRes('RT olig as PCR olig , '+rtOligs[foundedRt][2],resPCR2)
    if pcrOligsNotInTheStart!={}:
        for key in pcrOligsNotInTheStart.keys():
            if pcrOligsNotInTheStart[key][1]==0:
                resPCR1=makeRes(pcrNucs4s[key[:4]][key]+' dest='+str(pcrOligsNotInTheStart[key][0]),resPCR1)
            if pcrOligsNotInTheStart[key][1]==1:
                resPCR2=makeRes(pcrNucs4s[key[:4]][key]+' dest='+str(pcrOligsNotInTheStart[key][0]),resPCR2)
    makeResList(resPCR1, pcr1)
    makeResList(resPCR2, pcr2)
    makeResList(resMism, mism)
    makeResList(resRt, rt)
    makeResList(resDest, dest)
    makeResList(pcrLen,pcrLenList)
faFile1=open('1.'.join(args.fastq.split('.')),'w')
faFile2=open('2.'.join(args.fastq.split('.')),'w')
for i in range(len(toFastaNum)):
    faFile1.write('@'+toFastaNum[i]+'\n'+toFastaSeq1[i]+'\n'+'+'+'\n'+'F'*len(toFastaSeq1[i])+'\n')
    faFile2.write('@'+toFastaNum[i]+'\n'+toFastaSeq2[i]+'\n'+'+'+'\n'+'F'*len(toFastaSeq2[i])+'\n')
faFile1.close()
faFile2.close()

if args.indexT or args.indexG:
    sens='--bowtie2-dp 0 -k 5 --score-min L,0,-0.2'
    if args.sensitivity=='sensetive':
        sens='--bowtie2-dp 1 -k 30 --score-min L,0,-0.5'
    if args.sensitivity=='very-sensetive':
        sens='--bowtie2-dp 2 -a --score-min L,0,-1'
if args.indexT:
    hisat='hisat2 -x '+args.indexT+' -U '+'1.'.join(args.fastq.split('.'))+' -S '+'1.'.join(args.Tsam.split('.'))+' '+sens
    result=subprocess.run(hisat , shell=True)
    hisat='hisat2 -x '+args.indexT+' -U '+'2.'.join(args.fastq.split('.'))+' -S '+'2.'.join(args.Tsam.split('.'))+' '+sens
    result=subprocess.run(hisat , shell=True)
if args.indexG:
    hisat='hisat2 -x '+args.indexG+' -U '+'1.'.join(args.fastq.split('.'))+' -S '+'1.'.join(args.Gsam.split('.'))+' '+sens
    result=subprocess.run(hisat , shell=True)
    hisat='hisat2 -x '+args.indexG+' -U '+'2.'.join(args.fastq.split('.'))+' -S '+'2.'.join(args.Gsam.split('.'))+' '+sens
    result=subprocess.run(hisat , shell=True)

tGenesDict1=getTGenesDict('1.'.join(args.Tsam.split('.')))
tGenesDict2=getTGenesDict('2.'.join(args.Tsam.split('.')))
gGenesDict1=getGGenesDict('1.'.join(args.Gsam.split('.')),tGenesDict1)
gGenesDict2=getGGenesDict('2.'.join(args.Gsam.split('.')),tGenesDict2)

allExonsTranscriptome={}
allExonsGenome={}
genVars={}
errorlist=[]
fn1=gffutils.example_filename(args.tgff)
for feature in gffutils.DataIterator(fn1):
    if feature[2]=="exon":
        if feature[0] not in allExonsTranscriptome.keys():
            allExonsTranscriptome[feature[0]]={}
        allExonsTranscriptome[feature[0]][feature[8].items()[0][1][0]]=(feature[3],feature[4])
        #pseudo=true
fn2=gffutils.example_filename(args.ggff)
for feature in gffutils.DataIterator(fn2):
    if feature[2]=="exon" and feature[0][:2]=='NC' and ('pseudo', ['true']) not in feature[8].items():
        if int(feature[0][6:][:3])==920:
            chr='chrM'
        elif int(feature[0][7])!=0:
            if int((feature[0][7]))==2:
                if int((feature[0][8]))==3:
                   chr='chrX'
                if int((feature[0][8]))==4:
                    chr='chrY'
                if int((feature[0][8]))!=3 and int((feature[0][8]))!=4:
                    chr='chr'+feature[0][7]+feature[0][8]
            else:
                chr='chr'+feature[0][7]+feature[0][8]
        else:
            chr='chr'+feature[0][8]
        if chr not in allExonsGenome.keys():
            allExonsGenome[chr]={}
        for it in feature[8].items():
            if it[0]=='gene':
                var=getVar(feature[8].items())
                if it[1][0]+chr not in genVars.keys():
                    genVars[it[1][0]+chr]=var
                    allExonsGenome[chr][it[1][0]]={}
                    allExonsGenome[chr][it[1][0]][feature[8].items()[0][1][0]]=(feature[3],feature[4])
                else:
                    if var<genVars[it[1][0]+chr]:
                        allExonsGenome[chr][it[1][0]]={}
                        allExonsGenome[chr][it[1][0]][feature[8].items()[0][1][0]]=(feature[3],feature[4])
                        genVars[it[1][0]+chr]=var
                    if var==genVars[it[1][0]+chr]:
                        allExonsGenome[chr][it[1][0]][feature[8].items()[0][1][0]]=(feature[3],feature[4])

genes1,cords1=getGenesAndCoords(nRows,tGenesDict1,gGenesDict1,allExonsTranscriptome,allExonsGenome)
genes2,cords2=getGenesAndCoords(nRows,tGenesDict2,gGenesDict2,allExonsTranscriptome,allExonsGenome)

data = {"R1" : (reads1), "R2" : (reads2), "Quantity" : (quantity), "PCR Olig in R1" : (pcr1), "PCR Olig in R2" : (pcr2), "Missmatches in PCR Olig" : (mism),"RT Olig" : (rt), "Possition of RT Olig" : (dest), "Cordinates R1": (cords1), "Genes R1" : (genes1), "Cordinates R2": (cords2),"Genes R2": (genes2)}
df = pd.DataFrame (data) 
df.to_excel(args.outFile,index=False) 


      
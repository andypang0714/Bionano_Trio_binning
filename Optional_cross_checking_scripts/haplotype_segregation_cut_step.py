__author__ = 'jlee'
import os 
import numpy as np
import pandas as pd
import sys
import subprocess
import inspect as i
from collections import Counter
import linecache
import argparse

parser=argparse.ArgumentParser(
    description='''Arguments to run: ''',
    epilog="""Make sure you have all the arguments!.""")
parser.add_argument('-m', metavar='\b', type=str, help='Map to cut. Binned assembly in trio-binning pipeline')
parser.add_argument('-c', metavar='\b', type=str, help='File with cuts listed')
args=parser.parse_args()
maptobecut=args.m
cutfilepath=args.c

cutfile=pd.read_table(cutfilepath)
freq = Counter(cutfile.cutType)
print(freq)

tobecuts= cutfile[(cutfile.cutType == "cut") | (cutfile.cutType == "alt_homo_same_tobecut")]
tobecuts=tobecuts.sort_values(by=['QryContigID'])
reffreq = Counter(tobecuts.QryContigID)
print(reffreq)
reffreq=pd.DataFrame.from_dict(reffreq, orient='index').reset_index()
reffreq.columns = ["QryContigID",'count']

command0="less {} | grep ^# | wc -l".format(maptobecut)
print(command0)
num=subprocess.check_output(command0, shell=True)
num=int(num)

header=linecache.getline(maptobecut, num-1)

cmap=pd.read_table(maptobecut,sep='\t', header=num-1)
cmap.columns=header.split("\t")
cmap.rename(columns={'#h CMapId': 'CMapId'}, inplace=True)
print(cmap.head())

def eachcontigcut(contignum):     ###For each contig, find all cuts
    print(contignum)
    eachmap=cmap[cmap.CMapId == contignum]
    look=tobecuts[tobecuts.QryContigID ==contignum]
    print(look)
    #print(eachmap)
    cutloc=eachmap[(eachmap['SiteID'].isin(look['QryStartId'])) | (eachmap['SiteID'].isin(look['QryEndID']))]
    cutat = [ '%.3f' % elem for elem in cutloc.Position/1000 ]
    cutat=' '.join(cutat)
    print(cutat)    
    perresult=str(round(contignum,0)) +" "+ str(cutat)
    print("rounding here")
    print(perresult)
    return(perresult)
	
getallcutsdf=reffreq['QryContigID'].astype(int).apply(eachcontigcut)

getallcuts=getallcutsdf.tolist()
print(getallcuts)

printallcuts = ', '.join(getallcuts).replace(",", " ")
print(printallcuts)


if(cmap.CMapId.max()<10000):
    command1="/home/users3/tanantharaman/tools/10195/RefAligner -i {} -break 10000 {} -f -o {}/{}_crosschecked_cut -merge".format(maptobecut, printallcuts, os.path.dirname(maptobecut), os.path.basename(maptobecut).split(".")[0])
    print(command1)
    #os.system(command1)
    subprocess.check_output(command1, shell=True)


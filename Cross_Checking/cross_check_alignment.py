__author__ = 'jlee'
import os
import sys
#import pandas
#import atexit
#import utilities as util
#import traceback
import subprocess
import argparse
#from subprocess import PIPE
#import fnmatch
#print(os.getcwd())
#pip install rpy2==2.3.0

parser=argparse.ArgumentParser(
    description='''Arguments to run: ''',
    epilog="""Make sure you have all the arguments!.""")
parser.add_argument('-p', metavar='\b', type=str, help='Pipeline folder')
parser.add_argument('-r', metavar='\b', type=str, help='Anchor genome map')
parser.add_argument('-b', metavar='\b', type=str, help='Query genome map')
parser.add_argument('-a', metavar='\b', type=str, help='SV calling opt arg')
parser.add_argument('-o', metavar='\b', type=str, help='output')
parser.add_argument('-mask', metavar = '\b', type=str, help='an optional masking bed file')
parser.add_argument('-RefAligner', metavar='\b',type=str, help='RefAligner path')
args=parser.parse_args()

pipeline=args.p
ref=args.r
query=args.b
optarg=args.a
RefAligner=args.RefAligner
outputfolder=args.o

#os.chdir(outputfolder)
#print("changed directory to " + outputfolder)

alignref_final=RefAligner + " -ref "+ ref + " -i " + query + " -o " + outputfolder + "/exp_refineFinal1 -stdout -stderr -maxthreads 112 -M 1 3 -Msave 1 -T 1e-12 -A 8 -S 0.1 -L 60.0 -outlier 3e-4 -endoutlier 3e-2 -outlierBC -res 2.6 -resSD 0.7 -extend 1 -nosplit 2 -biaswt 0.2 -biaswtEnd 0.0 -biaswtOutlier 0.0 -deltaX 6 -deltaY 6 -PVres 2 -PVendoutlier -hashgen 5 7 2.4 1.5 0.05 5.0 1 1 3 -hash -hashdelta 16 10 46 -hashoffset 1 -hashrange 0 -MultiMatchesTotScore 2 12.0 -RefSplit 1e-4 24 1e-5 -RefSplitStitch 1 -mres 1e-3 -hashGrouped 5 7 -hashMultiMatch 30 -hashGC 300 -hashT2 1 -hashkeys 1 -HSDrange 1.0 -insertThreads 4 -ScaleDelta 0.05 2 -ScaleDeltaBPP -hashScaleDelta 2 -xmapchim 14 2000.0 -xmapUnique 14 -xmaplen -AlignRes 2.0 -outlierExtend 6 120 -Kmax 4 -rres 0.9 -resEstimate -MultiMatches 5 -MultiMatchesDelta 50.0 -MultiMatchesFilter 2 -outlierLambda 10.0 -outlierType1 0 -FP 0.2 -FN 0.02 -sf 0.1 -sd 0.0 -sr 0.02 -se 0.2 -MinSF 0.05 -MaxSF 0.1 -MinSD 0.0 -MaxSD 0.0 -MaxSR 0.02 -MaxSE 0.5 -indel -finalsort-sitesdec -RAmem 3 300 -f -maxmem 120 -maxvirtmem 0 -ScaleDelta 0.05 8 -ScaleDeltaBPP -hashScaleDelta 2"

print(alignref_final)
subprocess.call(alignref_final, shell=True)

if args.mask:
    mask=args.mask
    SVcalling="python " + pipeline + "/runSV.py -t " + RefAligner + " -r "+ ref + " -q " + query + " -o " + outputfolder + "/exp_refineFinal1_sv -T 64 -j 64 -E " + outputfolder + "/exp_refineFinal1.errbin -b " + mask + " -a " + optarg +"&>SV_calling.log&"
else:

    SVcalling="python " + pipeline + "/runSV.py -t " + RefAligner + " -r "+ ref + " -q " + query + " -o " + outputfolder + "/exp_refineFinal1_sv -T 64 -j 64 -E " + outputfolder + "/exp_refineFinal1.errbin -a " + optarg +"&>SV_calling.log&"

print(SVcalling)
subprocess.call(SVcalling, shell=True)

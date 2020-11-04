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
parser.add_argument('-r', metavar='\b', type=str, help='Anchor genome map')
parser.add_argument('-b', metavar='\b', type=str, help='Query molecules')
parser.add_argument('-o', metavar='\b', type=str, help='output')
parser.add_argument('-RefAligner', metavar='\b',type=str, help='RefAligner path')
parser.add_argument('-outputmolname', metavar='\b', type=str, help='output name of aligned molecule')
args=parser.parse_args()

ref=args.r
query=args.b
RefAligner=args.RefAligner
outputfolder=args.o
outputmolname=args.outputmolname
#os.chdir(outputfolder)
#print("changed directory to " + outputfolder)
#os.mkdir(outputfolder+"/autonoise")
#os.mkdir(outputfolder+"/alignmolvref")
#print("make "+ outputfolder+"/autonoise")


###Proband molecule aligns to Parental Genome Maps
#1. Run Autonoise0 and Autonoise1
Autonoise0= RefAligner + " -f -i " + query + " -ref " + ref + " -o " + outputfolder + "/autonoise/autoNoise0 -stdout -stderr -MapRate 0.6 1e-11 -usecolor 1 -FP 1.0 -FN 0.10 -sf 0.12 -sd 0.0 -sr 0.02 -res 3.1 -resSD 0.75 -minlen 120 -minsites 9 -minlen 150 -maxInterval 1000 -usecolor 1 -MinSF 0.0 -MaxSF 0.3 -MaxSD 0.12 -MaxSR 0.04 -MaxSE 0.5 -se 0.25 -L 0 -BestRef 1 -BestRefPV 1 -outlier 1e-3 -outlierMax 40. -endoutlier 1e-4 -nosplit 2 -biaswt 0.0 -S -1000.0 -PVres 2 -PVendoutlier -AlignRes 1.5 -rres 0.9 -resEstimate -f -maptype 0 -RAmem 3 300 -maxmem 120 -maxvirtmem 0 -mres 0.9 -T 1e-11 -A 9 -M 1 3 -sort-runindexShuffle 150.0 1 -subset 1 50000 0 4 -maxEnd 90. -minSNRestimate 2.0 2.0 0.5 100 0.9 -bpp 475 -bppScan 15 5 80 -ScaleDelta 0.03 10 -ScaleDeltaBPP -hashScaleDelta 2 -ScanScaling 2 2.0 -MaxSF 0.12 -MaxSD 0.0 -MaxSR 0.02 -MaxSE 0.25 -hashgen 5 4 2.4 1.4 0.05 5.0 1 1 1 -hash -hashdelta 14 10 24 -hashoffset 1 -hashrange 1 -hashGC 300 -hashT2 1 -hashkeys 1 -hashMultiMatch 30 10 -insertThreads 4 -XmapStatWrite " + outputfolder + "/molecule_stats.txt -TotalThreads 96"
print(Autonoise0)
subprocess.call(Autonoise0, shell=True)

Autonoise1=RefAligner + " -f -i " + query + " -ref " + ref + " -o " + outputfolder + "/autonoise/autoNoise1 -stdout -stderr -FP 0.658814 -readparameters " + outputfolder + "/autonoise/autoNoise0.errbin -res 2.544 -bpp 513.83 -sr 0.015502 -sf 0.12 -FN 0.11332 -sd 0.0 -minlen 120 -minsites 9 -minlen 150 -maxInterval 1000 -usecolor 1 -MinSF 0.0 -MaxSF 0.3 -MaxSD 0.12 -MaxSR 0.04 -MaxSE 0.5 -se 0.25 -L 0 -BestRef 1 -BestRefPV 1 -outlier 1e-3 -outlierMax 40. -endoutlier 1e-4 -nosplit 2 -biaswt 0.0 -S -1000.0 -PVres 2 -PVendoutlier -AlignRes 1.5 -rres 0.9 -resEstimate -f -maptype 0 -RAmem 3 300 -maxmem 120 -maxvirtmem 0 -T 1e-11 -A 9 -M 2 4 -resbias 4.0 64 -sort-runindexShuffle 150.0 1 -subset 1 100000 0 10 -ScanScaling 2 -Hash_Bits 23 -hashgen 5 4 2.4 1.4 0.05 5.0 1 1 1 -hash -hashdelta 14 10 24 -hashoffset 1 -hashrange 1 -ScaleDelta 0.03 3 -ScaleDeltaBPP -hashScaleDelta 2 -hashGC 300 -hashT2 1 -hashkeys 1 -hashbest 0 -insertThreads 4 -maxEnd 90. -minsites 9 -maxsites 750 -maxSiteDensity 30 -MinSD 0.0 -biaswt 0.7 -biaswtEnd 0.0 -biaswtOutlier 0.0 -S 0.1 -XmapStatWrite " + outputfolder + "/molecule_stats_1.txt -TotalThreads 96 -maxthreads 96"
print(Autonoise1)
subprocess.call(Autonoise1, shell=True)

alignmolvref=RefAligner + " -ref " + ref + " -i " + outputfolder + "/autonoise/autoNoise1_rescaled.bnx -o " + outputfolder + "/alignmolvref/exp_refineFinal1 -f -stdout -stderr -maxthreads 64 -usecolor 1 -FP 1.0 -FN 0.10 -sf 0.12 -sd 0.0 -sr 0.02 -res 3.1 -resSD 0.75 -readparameters " + outputfolder + "/autonoise/autoNoise1.errbin -T 1e-11 -A 8 -L 80 -S -1000 -nosplit 2 -biaswt 0.0 -biaswtEnd 0.0 -biaswtOutlier 0.0 -res 3.1 -resSD 0.75 -extend 1 -BestRef 1 -BestRefPV 1 -maptype 0 -PVres 2 -PVendoutlier -AlignRes 1.5 -hashrange 1 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 30 5 -insertThreads 4 -finalsort-sitesdec -RAmem 3 30 -f -MinSD 0.0 -minlen 120 -minsites 9 -minlen 150 -maxsites 750 -maxSiteDensity 30 -maxContigSiteDensity 40 -MaxSD 0.12 -MaxSR 0.03 -outlier 0.0001 -endoutlier 0.0001 -sf 0.25 -sd 0.11 -minsites 9 -hashgen 5 3 2.4 1.4 0.05 5.0 1 1 -hash -hashdelta 14 10 24 -mres 0.9 -hashGC 300 -hashT2 1 -insertThreads 4 -rres 0.9 -maxmem 120 -maxvirtmem 0 -mapped " + outputfolder + "/alignmolvref/" + outputmolname

print(alignmolvref)
subprocess.call(alignmolvref, shell=True)

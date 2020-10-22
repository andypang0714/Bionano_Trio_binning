# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:45:43 2020

@author: jlee
"""

#conda activate bionano_minimal
##pip install psutil
##pip install multiprocess
#pip install snakemake
#conda install snakeparse
#pip install pyyaml

import os
import sys
from snakeparse.parser import argparser

configfile:"temp_config.yaml"   ###Snakemake stores it as "config"
print(config)

anchor1=config['f']
anchor2=config['m']
mol=config['b']
RefAligner=config['RefAligner']
outputfolder=config['o']
fname=config['na']
mname=config['nb']

rule done:
    input:
        #bnx1b=outputfolder + "/mol_align_to_" + fname + "/alignmolvref/mol_align_to_" + fname + "_postbinning_molID.txt",
        #bnx2b=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/mol_align_to_" + mname + "_postbinning_molID.txt",
        #malignmentout=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/mol_align_to_" + mname + "_prebinning.bnx",
        #falignmentout= outputfolder + "/mol_align_to_" + fname + "/alignmolvref/mol_align_to_" + fname + "_prebinning.bnx", #{outputfolder}/mol_align_to_{fname}/mol_align_to_{fname}_prebinning.bnx",
        bbnx1=outputfolder + "/mol_align_to_" + fname + "/alignmolvref/mol_align_to_" + fname + "_postbinning.bnx",
        bbnx2=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/mol_align_to_" + mname + "_postbinning.bnx",


rule make_folder:
    input:
        proband_mol=mol,
        father_cmap=anchor1,
        mother_cmap=anchor2,
        name=fname,
        mname=mname
    output:
        output_father=directory(outputfolder + "/mol_align_to_" + fname),
        output_mother=directory(outputfolder + "/mol_align_to_" + mname)
    run:
        os.mkdir(output.output_father)
        os.mkdir(output.output_mother)
        print ("make " + output.output_father)
        print ("make " + output.output_mother)

rule autonoise_alignmolvref_of_proband_mol_to_father_assembly:
    input:
        align_bnx_to_anchor="Scripts/Align_bnx_to_anchor.py",
        proband_mol=mol,
        father_cmap=anchor1,
        outputfolder=outputfolder,

    output:
        falignmentout= outputfolder + "/mol_align_to_" + fname + "/alignmolvref/mol_align_to_" + fname + "_prebinning.bnx", #{outputfolder}/mol_align_to_{fname}/mol_align_to_{fname}_prebinning.bnx",
        all1=outputfolder + "/mol_align_to_" + fname + "/autonoise/autoNoise1_rescaled.bnx",
        xmap1=outputfolder + "/mol_align_to_" + fname + "/alignmolvref/exp_refineFinal1.xmap",
    run:
        shell("echo 'python {input.align_bnx_to_anchor} -r {input.father_cmap} -b {input.proband_mol} -o {outputfolder}/mol_align_to_{fname} -RefAligner {RefAligner} -outputmolname mol_align_to_{fname}_prebinning'")
        shell("python {input.align_bnx_to_anchor} -r {input.father_cmap} -b {input.proband_mol} -o {outputfolder}/mol_align_to_{fname} -RefAligner {RefAligner} -outputmolname mol_align_to_{fname}_prebinning")
        print("alignmolvref for father done")

rule autonoise_alignmolvref_of_proband_mol_to_mother_assembly:
    input:
        align_bnx_to_anchor="Scripts/Align_bnx_to_anchor.py",
        proband_mol=mol,
        mother_cmap=anchor2,
        outputfolder=outputfolder,

    output:
        malignmentout=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/mol_align_to_" + mname + "_prebinning.bnx",
        all2=outputfolder + "/mol_align_to_" + mname + "/autonoise/autoNoise1_rescaled.bnx",
        xmap2=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/exp_refineFinal1.xmap",
    run:
        shell("echo 'python {input.align_bnx_to_anchor} -r {input.mother_cmap} -b {input.proband_mol} -o {outputfolder}/mol_align_to_{mname} -RefAligner {RefAligner} -outputmolname mol_align_to_{mname}_prebinning'")
        shell("python {input.align_bnx_to_anchor} -r {input.mother_cmap} -b {input.proband_mol} -o {outputfolder}/mol_align_to_{mname} -RefAligner {RefAligner} -outputmolname mol_align_to_{mname}_prebinning")
        print("alignmolvref for mother done")

rule trio_binning:
    input:
        trio_binning_rscript="Scripts/trio_binning_20190726_rscript.R",
        outputfolder=outputfolder,
        bnx1=outputfolder + "/mol_align_to_" + fname + "/alignmolvref/mol_align_to_" + fname + "_prebinning.bnx",
        bnx2=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/mol_align_to_" + mname + "_prebinning.bnx",
        xmap1=outputfolder + "/mol_align_to_" + fname + "/alignmolvref/exp_refineFinal1.xmap",
        xmap2=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/exp_refineFinal1.xmap",
        all1=outputfolder + "/mol_align_to_" + fname + "/autonoise/autoNoise1_rescaled.bnx",
        all2=outputfolder + "/mol_align_to_" + mname + "/autonoise/autoNoise1_rescaled.bnx",

    output:
        bnx1b=outputfolder + "/mol_align_to_" + fname + "/alignmolvref/mol_align_to_" + fname + "_postbinning_molID.txt",
        bnx2b=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/mol_align_to_" + mname + "_postbinning_molID.txt"

    run:
        shell("echo 'Rscript {input.trio_binning_rscript} -a {input.bnx1} -b {input.bnx2} -x {input.xmap1} -y {input.xmap2} -m {input.all1} -n {input.all2} > {outputfolder}/trio_binning_summary.txt'"),
        shell("Rscript {input.trio_binning_rscript} -a {input.bnx1} -b {input.bnx2} -x {input.xmap1} -y {input.xmap2} -m {input.all1} -n {input.all2} > {outputfolder}/trio_binning_summary.txt"),
        shell("echo 'done binning_analysis'")

rule extract_father_binned_molecules:
    input:
         bnxID1=outputfolder + "/mol_align_to_" + fname + "/alignmolvref/mol_align_to_" + fname + "_postbinning_molID.txt",

    output:
         bbnx1=outputfolder + "/mol_align_to_" + fname + "/alignmolvref/mol_align_to_" + fname + "_postbinning.bnx",

    run:
        shell("{RefAligner} -selectidf {input.bnxID1} -i {mol} -o {outputfolder}/mol_align_to_{fname}/alignmolvref/mol_align_to_{fname}_postbinning -merge -bnx")
        shell("echo 'done extracting father binned molecules'")

rule extract_mother_binned_molecules:
    input:
         bnxID2=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/mol_align_to_" + mname + "_postbinning_molID.txt",

    output:
         bbnx2=outputfolder + "/mol_align_to_" + mname + "/alignmolvref/mol_align_to_" + mname + "_postbinning.bnx",

    run:
        shell("{RefAligner} -selectidf {input.bnxID2} -i {mol} -o {outputfolder}/mol_align_to_{mname}/alignmolvref/mol_align_to_{mname}_postbinning -merge -bnx")
        shell("echo 'done extracting mother binned molecules'")


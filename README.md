# Trio_binning

## General Information
This project is for trio-binning Bionano molecules. Proband molecules are aligned to parental haplotype aware assemblies to identify paternal and maternal allele molecules. 

## Requirements
RefAligner - downloaded from Bionano Genomics software download page (Bionano Tools in https://bionanogenomics.com/support/software-downloads/).  
Run has to be launched from the directory of make_config.py.\
Currently only 1 run can be done per time. 

Other requirements include:\
python==3.7.3\
pyYAML==5.3.1\
snakemake==5.26.1\
r==3.4.3\
r-optparse=1.6.6\
snakeparse==0.1.0-py_2

## Example Run:
Generate the temp_config.yaml:\
Bionano_Trio_binning$ python make_config.py -na father -nb mother -f father_EXP_REFINEFINAL1.cmap -m mother_EXP_REFINEFINAL1.cmap -b Proband_mol_folder/Proband_all.bnx -o Proband_mol_folder -RefAligner path/to/RefAligner

Launch the run:\
Bionano_Trio_binning$ screen -L snakemake --jobs

## For Help
python make_config.py -h

## Result
The final binned molecules are in `<your ouput folder>`/mol_align_to_\*/alignmolvref/mol_align_to_\*_postbinning.bnx. The metrics for the binning is in the Summary section (at the buttom) in `<your ouput folder>`/trio_binning_summary.txt 

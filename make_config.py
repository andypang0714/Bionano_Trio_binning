###e.g. JLwrittenscripts/trio_binning$ python command.py -na father -nb mother -f father_EXP_REFINEFINAL1.cmap -m mother_EXP_REFINEFINAL1.cmap -b Proband_mol_folder/Proband_all.bnx -o Proband_mol_folder -RefAligner path/to/RefAligner
import argparse
import os.path
import snakemake
import sys
import io
#import configparser
import yaml

#thisdir = os.path.abspath(os.path.dirname(__file__))
#parentdir = os.path.join(thisdir,'..')
#cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):
    parser=argparse.ArgumentParser(
    description='''Arguments to run: ''',
    epilog="""Make sure you have all the arguments!.""")
    parser.add_argument('-f', metavar='\b', type=str, help='anchor1 genome map', required=True)
    parser.add_argument('-na', metavar='\b', type=str, help='anchor1 name', required=True)
    parser.add_argument('-m', metavar='\b', type=str, help='anchor2 genome map', required=True)
    parser.add_argument('-nb', metavar='\b', type=str, help='anchor2 name', required=True)
    parser.add_argument('-b', metavar='\b', type=str, help='query molecules', required=True)
    parser.add_argument('-o', metavar='\b', type=str, help='output folder', required=True)
    parser.add_argument('-RefAligner', metavar='\b',type=str, help='RefAligner path', required=True)
    args=parser.parse_args()

    print(args)

    with io.open('temp_config.yaml', 'w', encoding='utf8') as outfile:
        yaml.dump(vars(args), outfile, default_flow_style=False, allow_unicode=True )

if __name__=='__main__':
    main()




animport glob
import argparse

from genome_vcf_to_maf_functions import *

parser = argparse.ArgumentParser(description='read path and parameters')
parser.add_argument('--samples_file',  metavar='p',help='file containing the list of')
parser.add_argument('--destination',  metavar='d',help='...')
parser.add_argument('--add_to_existing',  metavar='n',help='...',default=0)
#
args = parser.parse_args()
samples_file=args.samples_file 
destination=args.destination
add_to_existing=args.add_to_existing

#python genomes_to_maf.py --samples_file '.....' --destination '...'


samples_metadata=pd.read_table(samples_file,index_col=0)
samples_metadata.index=samples_metadata['Sample_folder']+samples_metadata['Sample_file']
print ("\nTo change the sample list please edit the file: "+samples_file+"\n")
''' MAIN '''

'''alternative code for searching vcf files in fodlers'''
#files=np.hstack([glob.glob(f+"/*vcf*")for f in np.unique(samples_metadata['Sample_folder'])])

chroms=np.hstack([[str(c) for c in np.arange(1,23)],"MT",'X','Y'])
if add_to_existing==0: create_new_MAF_results_file(chroms,destination)
for file in samples_metadata.index:
    
    [skip,samples]=get_sample_first_line(file,samples_metadata)
    
    process_input_chunks(file,skip,samples,chunksize = 10 ** 5,nrows= 10**8,chroms=chroms,destination=destination)
   
    
    
        

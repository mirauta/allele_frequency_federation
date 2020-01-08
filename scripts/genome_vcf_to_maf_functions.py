import numpy as np
import pandas as pd
import os


def funreplace(x): return x.replace('chr','')
def getfrequency(x): return np.array(int(x[0])+int(x[2]))

def get_sample_first_line(file,samples_metadata):
#    print (file)
    temp0=pd.read_table(file, nrows=1000,sep="£££££%%")
    skip=np.where(["#CHROM"==t[0][:6] for t in temp0.values])[0][0]
    
    samples=np.array(samples_metadata.loc[file,'Sample_names'].split(';'))
    print ("\n\n>>>>>>>>>>>>>\n")
    print ("The following samples will be aggregated for "+file)
    print (samples)
    print ("\n>>>>>>>>>>>>>\n\n")
    return[skip,samples]


             
#       
def get_sample_frequency_one_chunk(dat,samples):
 
    dat=dat[(dat['REF']=="A")+(dat['REF']=="C")+(dat['REF']=="G")+(dat['REF']=="T")]
    dat=dat[(dat['ALT']=="A")+(dat['ALT']=="C")+(dat['ALT']=="G")+(dat['ALT']=="T")]
    try:    dat.index=dat['#CHROM'].apply(funreplace).astype("category") 
    except: dat.index=dat['#CHROM'].astype('U').astype("category")   


    dat['freqALT']=0
    for s in samples:
        dat['freqALT']=dat['freqALT']+dat[s].apply(getfrequency)
    for nt in ['A','C','G','T']:
        dat[nt+"sample"]=dat['freqALT']*(dat['ALT']==nt)+(2*samples.shape[0]-dat['freqALT'])*(dat['REF']==nt)
        
    return dat

def process_input_chunks(file,skip,samples,chunksize = 10 ** 1,nrows=10**6,chroms=[],destination=''):
    for ic,dat in enumerate(pd.read_table(file,skiprows=skip+1,\
                                          nrows=nrows,\
                                          chunksize=chunksize,\
                                          sep="\t")):
        if ic==0: dat =dat[1:]
        dat=dat[np.hstack([['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'],samples])]
        print (dat.shape)
        dat=get_sample_frequency_one_chunk(dat,samples)
#        print( "...........1")
            
        for c in np.intersect1d(chroms,dat.index.categories):
            datc=dat.loc[c].set_index('POS')
            datc=datc.iloc[np.unique(datc.index,return_index=1)[1]]
            dfchrom=pd.read_table(destination+"/MAF_output/"+"MAF_chrom_"+str(c)+".tsv",sep='\t',index_col=0)
#            print( "..........."+str(c))
            dfchrom=pd.concat([dfchrom,datc['Asample']],1); 
            dfchrom['A']=dfchrom[['A','Asample']].sum(1); del dfchrom['Asample'];
            for nt in ['C','G','T']:
                dfchrom=pd.concat([dfchrom,datc[nt+'sample']],1); dfchrom[nt]=dfchrom[[nt,nt+'sample']].sum(1); del dfchrom[nt+'sample'];
            dfchrom.to_csv(destination+"/MAF_output/"+"MAF_chrom_"+str(c)+".tsv",sep='\t')       
            print( " wrote ...........chrom "+str(c)+" wrote")
        
        
def create_new_MAF_results_file(chroms,destination):
    
    try:
        if os.path.isdir(destination+"/MAF_output/"):
            print ("Using this output directory: " + destination+"/MAF_output/\n\n")
    except:
        os.mkdir(destination+"/MAF_output/")
        print ("Creation of the directory: " % destination+"/MAF_output/\n\n")
 
    df={}
    for c in chroms:
        df[c]=pd.DataFrame(columns=['A','C','G','T'])
        df[c].to_csv(destination+"/MAF_output/"+"MAF_chrom_"+str(c)+".tsv",sep='\t')

import numpy as np
import pandas as pd
import os
import sys


def funreplace(x): return x.replace('chr','')
def getfrequency(x): return np.array(int(x[0])+int(x[2]))

class MAF():
    
    def __init__(self):
        print(sys.platform)
        if (sys.platform=="win32")|(sys.platform=="win64"):
            self.foldersep=r'\\'[0]
        else:
            self.foldersep="/"
        
    def get_sample_first_line(self,file,samples_metadata):
    #    print (file)
        temp0=pd.read_table(file, nrows=1000,sep="£££££%%",engine="python")
        skip=np.where(["#CHROM"==t[0][:6] for t in temp0.values])[0][0]
        
        try:
            samples=np.array(samples_metadata.loc[file,'Sample_names'].split(';'))
        except:
            samples=np.array(pd.read_table(file,skiprows=skip+1,\
                                              nrows=1000,\
                                              sep="\t").columns[[9]])
            samples_metadata.loc[file,'Sample_names']=samples[0]
        print ("\n\n>>>>>>>>>>>>>\n")
        print ("The following samples will be aggregated for "+file)
        print (samples)
        print ("\n>>>>>>>>>>>>>\n\n")
        return[skip,samples]


                 
    #       
    def get_sample_frequency_one_chunk(self,dat,samples):
     
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

    def process_input_chunks(self,file,skip,samples,chunksize = 10 ** 1,nrows=10**6,chroms=[],destination=''):
        for ic,dat in enumerate(pd.read_table(file,skiprows=skip+1,\
                                              nrows=nrows,\
                                              chunksize=chunksize,\
                                              sep="\t")):
            if ic==0: dat =dat[1:]
            dat=dat[np.hstack([['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'],samples])]
            print (dat.shape)
            dat=self.get_sample_frequency_one_chunk(dat,samples)
    #        print( "...........1")
                
            for c in np.intersect1d(chroms,dat.index.categories):
                datc=dat.loc[c].set_index('POS')
                datc=datc.iloc[np.unique(datc.index,return_index=1)[1]]
                dfchrom=pd.read_table(destination+self.foldersep+"MAF_output"+self.foldersep+"MAF_chrom_"+str(c)+".tsv",sep='\t',index_col=0)
    #            print( "..........."+str(c))
                dfchrom=pd.concat([dfchrom,datc['Asample']],1); 
                dfchrom['A']=dfchrom[['A','Asample']].sum(1); del dfchrom['Asample'];
                for nt in ['C','G','T']:
                    dfchrom=pd.concat([dfchrom,datc[nt+'sample']],1); dfchrom[nt]=dfchrom[[nt,nt+'sample']].sum(1); del dfchrom[nt+'sample'];
                dfchrom.to_csv(destination+self.foldersep+"MAF_output"+self.foldersep+"MAF_chrom_"+str(c)+".tsv",sep='\t')       
                print( "........... chrom "+str(c)+" done")
            


    def create_new_MAF_results_file(self,chroms,destination):
        temp=destination+self.foldersep+"MAF_output"+self.foldersep
        
        if os.path.isdir(temp):
            print ("Using this output directory: " + temp+"\n\n")
            print(1)

        else:
            print (temp)
            os.mkdir(temp)
            print ("Creation of the directory: " +temp+"\n\n")
        print(temp) 
        df={}
        for c in chroms:
            df[c]=pd.DataFrame(columns=['A','C','G','T'])
            df[c].to_csv(destination+self.foldersep+"MAF_output"+self.foldersep+"MAF_chrom_"+str(c)+".tsv",sep='\t')

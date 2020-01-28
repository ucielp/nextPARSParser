# -*- coding: utf-8 -*-
"""
Created on Wed May 31 19:38:02 2017

@author: jwillis
"""

import argparse, sys
import numpy as np #... dont use??? since numpy wont seem to import on the cluster
from datetime import datetime
from termcolor import colored


##################################
def load_tab(tabFile):
    
    loaded_tab = {}
    
    with open(tabFile) as readTab:
        for line in readTab:
            try:
                name, counts = line.split()
            except ValueError:
                sys.stderr.write( '%s\t%s\n' %(len(line.split()), line.split()) )
                continue
            
            if name == ".":
                continue
            
            # will remove ';' at the end and add them back later, in case there are any that do not end in ';'
            counts = counts.strip(';')
            
            if name not in loaded_tab:
                loaded_tab[name] = []
            loaded_tab[name].append( counts )
    
    return loaded_tab
##################################

def filter_tab(loaded_tab, minC, tabFile):
    
    mols_of_interest = ['226_neg_ROX2','ROX2','L1_4_5T1','MinH_4_5T1','2_4_5T1','beg_4_5T1','s2mHotair','mid_4_5T1','H7_4_5T1',
                        'STR_4_5T1','4_5exonT1','PRC2_mHotair','smHotair','HOT2','HOTAIR-ncbiRinn',
                        'VHA_0_mut','VHA_3_mut','VHA_6_mut','hSRA','mSRA','SRA','3UTR_long_SNCA','AK123408','CR610499',
                        'HAR1RDQ860410','B2','U1','TETp4p6','TETp9p9.1','BC041455','AL832444',
                        'ENST00000394991-5utr','ENST00000394991-cds','ENST00000394991-3utr',
                        'ERNST4_ORF_FWDstrand','ERNST4_ORF_REVstrand']

    tab = open('%s_temp' %tabFile, 'w')
#    tab = open('%s' %tabFile, 'w')
    
    all_means = []
    low_means = []
    
    for name in loaded_tab:
        
        counts = ';'.join(loaded_tab[name])
#        counts = [ int(x) for x in counts.split(';') ]
#        print colored(name,'red'), np.mean( [ int(x) for x in counts.split(';') ] )#, counts
        if np.mean( [ int(x) for x in counts.split(';') ] ) < minC and name not in mols_of_interest:
#        if float(sum( [ int(x) for x in counts.split(';') ] )) / len( [ int(x) for x in counts.split(';') ] ) < minC:
            low_means.append( np.mean( [ int(x) for x in counts.split(';') ] ) )
            continue
        
        all_means.append( np.mean( [ int(x) for x in counts.split(';') ] ) )
#        all_means.append( float(sum( [ int(x) for x in counts.split(';') ] ) ) / len( [ int(x) for x in counts.split(';') ] ) )
        
        # write transcripts that pass mincount filter
        tab.write( '%s\t%s;\n' %(name, counts) )
    
    tab.close()
    
    sys.stderr.write('Minimum average count per site threshold = %s\n' %minC)
    sys.stderr.write('Overall average counts per site for retained transcripts  = %s (%s transcripts)\n' %(np.mean( all_means ), len(all_means)))
    sys.stderr.write('Overall average counts per site for discarded transcripts = %s (%s transcripts)\n' %(np.mean( low_means ), len(low_means)))
#    sys.stderr.write('Overall average counts per site = %s\n' %float(sum( all_means )) / len( all_means ))
    
##################################
def main():
    
    parser = argparse.ArgumentParser(description="Reformat tab file output by PARSparser.jar",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", "--minC", dest="minC", default=5, type=int, help="Min average counts for given transcript.")
    parser.add_argument("-t", "--tab", dest="tab", required=True, help="tab file to be reformatted.")
    args = parser.parse_args()

    t0 = datetime.now()
    sys.stderr.write('\n\n## Reformatting tab file %s\n' %args.tab)
    sys.stderr.write('## Starting at %s\n\n' %t0)
    
    loaded_tab = load_tab(args.tab)
    
    filter_tab(loaded_tab, args.minC, args.tab)
    
    t1 = datetime.now() - t0
    sys.stderr.write('\n##Time elapsed during reformatting: %s\n' %t1)
#    for name in new_tab:
#        print colored(name, 'red'), new_tab[name]


##################################

if __name__ == "__main__":
    main()
# -*- coding: UTF-8 -*-

import sys
import os
import seaborn as sns
import numpy as np
import pandas as pd
import h5py
import scipy

# print(sys.argv)

startpos = sys.argv[1]      #zero-based position, in variant genome
endpos = sys.argv[2]        #zero-based position, in variant genome
targetpos = sys.argv[3]     #one-based, in ref genome
substrlen = sys.argv[4]     #upstream/downstream N bases
tombo_dir = sys.argv[5]     #"/data1/yefq/data/fast5/20220703_WGA_twist/processed/20220629_basecalled/workspace/fast5_pass/Omicron_BA.1/fast5_pass_single/all_single_fast5s"
id_dir = sys.argv[6]        #"/data1/yefq/data/fast5/20220703_WGA_twist/processed/20220724_comparison"
outputfile = sys.argv[7]      #"/data1/yefq/data/fast5/20220703_WGA_twist/processed/20220629_basecalled/workspace/fast5_pass/Omicron_BA.1/fast5_pass_single/test"


def tombo_extraction(fast5_path):
    """
    extract signal and substitute true str seq
    
    """
    
    fast5_file = h5py.File(fast5_path, 'r')
    
    if "RawGenomeCorrected_000" not in fast5_file['Analyses']:
        return None, "Error: RawGenomeCorrected_000", None
    
    if 'BaseCalled_template' not in fast5_file['Analyses']['RawGenomeCorrected_000']: 
        return None, "Error: BaseCalled_template", None
    if fast5_file['Analyses']['RawGenomeCorrected_000']['BaseCalled_template'].attrs['status'] != 'success': 
        return None, "Error: status!=success", None
    
    for readname in fast5_file['Raw']['Reads']:
        signal = fast5_file['Raw']['Reads'][readname]['Signal'][:]
        if len(signal) > 0:
            break
    if len(signal) == 0:
        return None, "Error: signal empty", None

    
    #### extract signal&query from fast5
    tombo_events        = fast5_file['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']['Events'][:]
    tombo_event_starts  = [entry[2] for entry in tombo_events]
    tombo_event_lengths = [entry[3] for entry in tombo_events]
    tombo_event_bases   = [entry[4].decode('utf-8') for entry in tombo_events]
    
    # tombo_signal_length = tombo_event_starts[-1] + tombo_event_lengths[-1]
    tombo_signal_length = sum([entry[3] for entry in tombo_events])
    
    # read_start_rel_to_raw (int) – start (in raw signal vector) of assigned bases
    tombo_signal_start  = fast5_file['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']['Events'].attrs['read_start_rel_to_raw'] 
    tombo_singal_end    = tombo_signal_start + tombo_signal_length
    feature_signal      = signal[tombo_signal_start:tombo_singal_end]

    #### amplicon_seq subsititution
    #raw signal length
    


    ## tombo_info
    mapped_chrom  = fast5_file['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']['Alignment'].attrs['mapped_chrom']
    mapped_start  = fast5_file['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']['Alignment'].attrs['mapped_start']
    mapped_end    = fast5_file['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']['Alignment'].attrs['mapped_end']
    mapped_strand = fast5_file['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']['Alignment'].attrs['mapped_strand']
    
    return mapped_chrom, mapped_start, mapped_end, mapped_strand, feature_signal, tombo_event_bases, tombo_event_starts, tombo_event_lengths
    



#==============================20221031 list signals for one region=========================================
ranges1 = [int(startpos),int(endpos)]
ranges = [int(startpos)-int(substrlen),int(endpos)+int(substrlen)]
rangelen = ranges[1]-ranges[0]+1
substrlen = int(substrlen)

plus_del_alllists = [[] for i in range(rangelen)]  #each position corresponding to one list
plus_del_alllabels = ["B" for i in range(rangelen)]

plus_match_alllists = [[] for i in range(rangelen)]  #each position corresponding to one list
plus_match_alllabels = ["B" for i in range(rangelen)]

minus_del_alllists = [[] for i in range(rangelen)]  #each position corresponding to one list
minus_del_alllabels = ["B" for i in range(rangelen)]

minus_match_alllists = [[] for i in range(rangelen)]  #each position corresponding to one list
minus_match_alllabels = ["B" for i in range(rangelen)]


plus_del_file = "target.plus.del.readID.txt"        #"BA.1_guppy_sup.pos9812.plus.del.readID.txt"
plus_match_file = "target.plus.match.readID.txt"    #"BA.1_guppy_sup.pos9812.plus.match.readID.txt"
minus_del_file = "target.minus.del.readID.txt"      #"BA.1_guppy_sup.pos9812.minus.del.readID.txt"
minus_match_file = "target.minus.match.readID.txt"  #"BA.1_guppy_sup.pos9812.minus.match.readID.txt"

print("ReadIDs start to load") 
#load readIDs
f = open(id_dir + '/' + plus_del_file,'r')
lines = f.readlines()
plus_del_readidlists= []
for line in lines:
    line = line.strip()
    plus_del_readidlists.append(line)
f.close()


f = open(id_dir + '/' + plus_match_file,'r')
lines = f.readlines()
plus_match_readidlists= []
for line in lines:
    line = line.strip()
    plus_match_readidlists.append(line)
f.close()


f = open(id_dir + '/' + minus_del_file,'r')
lines = f.readlines()
minus_del_readidlists= []
for line in lines:
    line = line.strip()
    minus_del_readidlists.append(line)
f.close()


f = open(id_dir + '/' + minus_match_file,'r')
lines = f.readlines()
minus_match_readidlists= []
for line in lines:
    line = line.strip()
    minus_match_readidlists.append(line)
f.close()
print("ReadIDs loading finished") 

print("Plus Del starts")
count = 0
for readid in plus_del_readidlists:
    count+=1
    if count%1000 == 0:
        print(count) 
    fast5_filepath = tombo_dir + '/' + readid + '.fast5'
    if len(tombo_extraction(fast5_filepath)) == 8:
        mapped_chrom, mapped_start, mapped_end, mapped_strand, feature_signal, tombo_event_bases, tombo_event_starts, tombo_event_lengths = tombo_extraction(fast5_filepath)
        dat = pd.DataFrame()
        dat['tombo_event_bases'] = tombo_event_bases
        dat['tombo_event_starts'] = tombo_event_starts
        dat['tombo_event_lengths'] = tombo_event_lengths
        dat['signal_list'] = dat.apply(lambda x: feature_signal[x['tombo_event_starts']:(x['tombo_event_starts'] + x['tombo_event_lengths'])],axis=1)
        dat['signal_mean'] = [entry.mean() for entry in dat['signal_list']]
        event_bases_string = "".join(tombo_event_bases)
        if mapped_strand == "+":
            if (mapped_start<=ranges[0]) and (mapped_end >= ranges[1]+1):
                tmp = dat.loc[(ranges[0]-mapped_start):(ranges[1]-mapped_start)]
                for i in range(rangelen):
                    plus_del_alllists[i].append(tmp['signal_list'].iloc[i])
                    plus_del_alllabels[i] = tmp['tombo_event_bases'].iloc[i]

print("Plus Match starts")
count = 0
for readid in plus_match_readidlists:
    count+=1
    if count%1000 == 0:
        print(count) 
    fast5_filepath = tombo_dir + '/' + readid + '.fast5'
    if len(tombo_extraction(fast5_filepath)) == 8:
        mapped_chrom, mapped_start, mapped_end, mapped_strand, feature_signal, tombo_event_bases, tombo_event_starts, tombo_event_lengths = tombo_extraction(fast5_filepath)
        dat = pd.DataFrame()
        dat['tombo_event_bases'] = tombo_event_bases
        dat['tombo_event_starts'] = tombo_event_starts
        dat['tombo_event_lengths'] = tombo_event_lengths
        dat['signal_list'] = dat.apply(lambda x: feature_signal[x['tombo_event_starts']:(x['tombo_event_starts'] + x['tombo_event_lengths'])],axis=1)
        dat['signal_mean'] = [entry.mean() for entry in dat['signal_list']]
        event_bases_string = "".join(tombo_event_bases)
        if mapped_strand == "+":
            if (mapped_start<=ranges[0]) and (mapped_end >= ranges[1]+1):
                tmp = dat.loc[(ranges[0]-mapped_start):(ranges[1]-mapped_start)]
                for i in range(rangelen):
                    plus_match_alllists[i].append(tmp['signal_list'].iloc[i])
                    plus_match_alllabels[i] = tmp['tombo_event_bases'].iloc[i]
                
print("Minus Del starts")
count = 0
for readid in minus_del_readidlists:
    count+=1
    if count%1000 == 0:
        print(count) 
    fast5_filepath = tombo_dir + '/' + readid + '.fast5'
    if len(tombo_extraction(fast5_filepath)) == 8:
        mapped_chrom, mapped_start, mapped_end, mapped_strand, feature_signal, tombo_event_bases, tombo_event_starts, tombo_event_lengths = tombo_extraction(fast5_filepath)
        dat = pd.DataFrame()
        dat['tombo_event_bases'] = tombo_event_bases
        dat['tombo_event_starts'] = tombo_event_starts
        dat['tombo_event_lengths'] = tombo_event_lengths
        dat['signal_list'] = dat.apply(lambda x: feature_signal[x['tombo_event_starts']:(x['tombo_event_starts'] + x['tombo_event_lengths'])],axis=1)
        dat['signal_mean'] = [entry.mean() for entry in dat['signal_list']]
        event_bases_string = "".join(tombo_event_bases)
        if mapped_strand == "-":
            if (mapped_start<=ranges[0]) and (mapped_end >= ranges[1]+1):
                tmp1 = dat.loc[(mapped_end-1-ranges[1]):(mapped_end-1-ranges[0])]
                tmp1 = tmp1.iloc[::-1] #reverse rows
                for i in range(rangelen):
                    minus_del_alllists[i].append(tmp1['signal_list'].iloc[i])
                    minus_del_alllabels[i] = tmp1['tombo_event_bases'].iloc[i]

print("Minus Match starts")
count = 0
for readid in minus_match_readidlists:
    count+=1
    if count%1000 == 0:
        print(count) 
    fast5_filepath = tombo_dir + '/' + readid + '.fast5'
    if len(tombo_extraction(fast5_filepath)) == 8:
        mapped_chrom, mapped_start, mapped_end, mapped_strand, feature_signal, tombo_event_bases, tombo_event_starts, tombo_event_lengths = tombo_extraction(fast5_filepath)
        dat = pd.DataFrame()
        dat['tombo_event_bases'] = tombo_event_bases
        dat['tombo_event_starts'] = tombo_event_starts
        dat['tombo_event_lengths'] = tombo_event_lengths
        dat['signal_list'] = dat.apply(lambda x: feature_signal[x['tombo_event_starts']:(x['tombo_event_starts'] + x['tombo_event_lengths'])],axis=1)
        dat['signal_mean'] = [entry.mean() for entry in dat['signal_list']]
        event_bases_string = "".join(tombo_event_bases)
        if mapped_strand == "-":
            if (mapped_start<=ranges[0]) and (mapped_end >= ranges[1]+1):
                tmp1 = dat.loc[(mapped_end-1-ranges[1]):(mapped_end-1-ranges[0])]
                tmp1 = tmp1.iloc[::-1] #reverse rows
                for i in range(rangelen):
                    minus_match_alllists[i].append(tmp1['signal_list'].iloc[i])
                    minus_match_alllabels[i] = tmp1['tombo_event_bases'].iloc[i]
print("All signals loaded")

plus_match_alllists_len = [[] for i in range(rangelen)]
plus_del_alllists_len = [[] for i in range(rangelen)]
minus_match_alllists_len = [[] for i in range(rangelen)]
minus_del_alllists_len = [[] for i in range(rangelen)]
match_alllists_len = [[] for i in range(rangelen)]
del_alllists_len = [[] for i in range(rangelen)]

for i in range(rangelen):
    for j in range(len(plus_match_alllists[i])):
        plus_match_alllists_len[i].append(len(plus_match_alllists[i][j]))
    for j in range(len(plus_del_alllists[i])):
        plus_del_alllists_len[i].append(len(plus_del_alllists[i][j]))
    for j in range(len(minus_del_alllists[i])):
        minus_del_alllists_len[i].append(len(minus_del_alllists[i][j]))
    for j in range(len(minus_match_alllists[i])):
        minus_match_alllists_len[i].append(len(minus_match_alllists[i][j]))
    match_alllists_len[i] = plus_match_alllists_len[i] + minus_match_alllists_len[i]
    del_alllists_len[i] = plus_del_alllists_len[i] + minus_del_alllists_len[i]


'''per read level'''
#upstream
plus_match_reads_len = [[] for j in range(len(plus_match_alllists[0]))]
plus_del_reads_len = [[] for j in range(len(plus_del_alllists[0]))]
minus_match_reads_len = [[] for j in range(len(minus_match_alllists[0]))]
minus_del_reads_len = [[] for j in range(len(minus_del_alllists[0]))]

plus_match_reads_signals = [[] for j in range(len(plus_match_alllists[0]))]
plus_del_reads_signals = [[] for j in range(len(plus_del_alllists[0]))]
minus_match_reads_signals = [[] for j in range(len(minus_match_alllists[0]))]
minus_del_reads_signals = [[] for j in range(len(minus_del_alllists[0]))]


plus_match_reads_total_len = [[] for j in range(len(plus_match_alllists[0]))]
plus_del_reads_total_len = [[] for j in range(len(plus_del_alllists[0]))]
minus_match_reads_total_len = [[] for j in range(len(minus_match_alllists[0]))]
minus_del_reads_total_len = [[] for j in range(len(minus_del_alllists[0]))]

for j in range(len(plus_match_reads_len)):
    for i in range(rangelen):
        plus_match_reads_len[j].append(len(plus_match_alllists[i][j]))
        plus_match_reads_signals[j].extend(plus_match_alllists[i][j])
    #print(plus_match_reads_signals[j])
    plus_match_reads_total_len[j] = sum(plus_match_reads_len[j])
for j in range(len(plus_del_reads_len)):
    for i in range(rangelen):
        plus_del_reads_len[j].append(len(plus_del_alllists[i][j]))
        plus_del_reads_signals[j].extend(plus_del_alllists[i][j])
    plus_del_reads_total_len[j] = sum(plus_del_reads_len[j])
for j in range(len(minus_match_reads_len)):
    for i in range(rangelen):
        minus_match_reads_len[j].append(len(minus_match_alllists[i][j]))
        minus_match_reads_signals[j].extend(minus_match_alllists[i][j])
    minus_match_reads_total_len[j] = sum(minus_match_reads_len[j])
for j in range(len(minus_del_reads_len)):
    for i in range(rangelen):
        minus_del_reads_len[j].append(len(minus_del_alllists[i][j]))
        minus_del_reads_signals[j].extend(minus_del_alllists[i][j])
    minus_del_reads_total_len[j] = sum(minus_del_reads_len[j]) 

resfile = outputfile
f = open(resfile,'w')
f.write("readid_prefix\tType\tPosition\tBase_Num\tGroup\tGroup_Signal_Length\tGroup_Signals\tGroup_Filtered_Signals\n")
for j in range(len(plus_match_reads_len)):
    tmp = [id for id in plus_match_reads_signals[j] if id >=300 and id <=500]
    res = "target" + "\t" + "Per_Read\t" + str(targetpos) + "\t" + str(rangelen) + "\t" + "plus_match" + "\t" + str(plus_match_reads_total_len[j]) + "\t" + ",".join('%s' %id for id in plus_match_reads_signals[j]) + "\t" + ",".join(['%s' %id for id in tmp]) + "\n"
    f.write(res)
for j in range(len(plus_del_reads_len)):
    tmp = [id for id in plus_del_reads_signals[j] if id >=300 and id <=500]
    res = "target" + "\t" + "Per_Read\t" + str(targetpos) + "\t" + str(rangelen) + "\t" + "plus_del" + "\t" + str(plus_del_reads_total_len[j]) + "\t" + ",".join('%s' %id for id in plus_del_reads_signals[j]) + "\t" + ",".join(['%s' %id for id in tmp]) + "\n"
    f.write(res)
for j in range(len(minus_match_reads_len)):
    tmp = [id for id in minus_match_reads_signals[j] if id >=300 and id <=500]
    res = "target" + "\t" + "Per_Read\t" + str(targetpos) + "\t" + str(rangelen) + "\t" + "minus_match" + "\t" + str(minus_match_reads_total_len[j]) + "\t" + ",".join('%s' %id for id in minus_match_reads_signals[j]) + "\t" + ",".join(['%s' %id for id in tmp]) + "\n"
    f.write(res)
for j in range(len(minus_del_reads_len)):
    tmp = [id for id in minus_del_reads_signals[j] if id >=300 and id <=500]
    res = "target" + "\t" + "Per_Read\t" + str(targetpos) + "\t" + str(rangelen) + "\t" + "minus_del" + "\t" + str(minus_del_reads_total_len[j]) + "\t" + ",".join('%s' %id for id in minus_del_reads_signals[j]) + "\t" + ",".join(['%s' %id for id in tmp]) + "\n"
    f.write(res)
f.close()

print("Job is done")

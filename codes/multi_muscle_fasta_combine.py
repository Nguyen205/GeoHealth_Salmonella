# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 10:34:30 2022

@author: 67535
"""
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i',dest='input',type=str,required=True,help="Define input folder containing separate aligned fasta")
parser.add_argument('-o',dest='output',type=str,required=True,help="Define output file path")
args=parser.parse_args()

dir_path=args.input.strip("'")
aln_file_list=[]
dir_list=os.listdir(dir_path)

for cur_file in dir_list:
    if os.path.splitext(cur_file)[1]=='.fasta':
        aln_file=os.path.join(dir_path,cur_file)
        #print(blast_file)
        aln_file_list.append(aln_file)

marker='>'

for num in range(0,len(aln_file_list)):
    print('Processing '+aln_file_list[num]+' ......')
    marker_loc=[]
    original_file_raw=open(aln_file_list[num])
    original_file=original_file_raw.readlines()
    original_file_raw.close()
    
    for i in range(0,len(original_file)):
        original_file[i]=original_file[i].strip('\n')
        
    for i in range(0,len(original_file)):
        if marker in original_file[i]:
            marker_loc.append(i)
    #print(marker_loc, len(marker_loc), len(original_file))
    
    attach_dict={}
    
    for i in range(0,len(marker_loc)):
        name_with_marker=original_file[marker_loc[i]].split(' ')[0]
        name_wo_marker=name_with_marker.split('>')[1]
        current_loc=marker_loc[i]
        if i==len(marker_loc)-1:
            next_loc=len(original_file)
        else:
            next_loc=marker_loc[i+1]
        title=original_file[current_loc].replace('\n','').replace('\r','')
        #print(title, file=output)
        sequence=''
        for j in range(current_loc+1, next_loc):
            sequence+=original_file[j].replace('\n','').replace('\r','')
        attach_dict[title]=sequence
    
    if num==0:
        start_dict=attach_dict.copy()
    else:
        for key in start_dict:
            start_dict[key]+=attach_dict[key]
        
output_file=open(args.output,'w')
print('Writing output file......')
for key in start_dict:
    output_file.write(key+'\n')
    output_file.write(start_dict[key]+'\n')
output_file.close()
      
        
        
        
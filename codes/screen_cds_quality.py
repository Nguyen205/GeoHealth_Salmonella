# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 17:53:38 2022

@author: 67535
"""

import argparse
import os
import re
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i',dest='input_ID',type=int,required=True,help='Define the ID of the output files obtained from "pick_CDS_hit_all.py"')
parser.add_argument('-f',dest='folder_path',type=str,required=True,help='Define the folder containing all tblastn files (outfmt = 14) (The same folder as the input -b in "pick_CDS_hit_all.py")')
args=parser.parse_args()

xml_path_file=open(args.folder_path.strip("'")+'/tblastn_screening_output/xml_file_paths/xml_file_path_'+str(args.input_ID)+'.txt')
xml_path_list=xml_path_file.readlines()
xml_path_file.close()
for i in range(0,len(xml_path_list)):
    xml_path_list[i]=xml_path_list[i].strip('\n')

srr_list_file=open(args.folder_path.strip("'")+'/tblastn_screening_output/srr_lists_after_first_screening/srr_list_after_first_screening_'+str(args.input_ID)+'.txt')
srr_list=srr_list_file.readlines()
srr_list_file.close()
for i in range(0,len(srr_list)):
    srr_list[i]=srr_list[i].strip('\n')

removed_dict={}

output_path=args.folder_path+'/tblastn_screening_output/fasta_for_align'
if not(os.path.exists(output_path)):
    os.mkdir(output_path)

for xml_path in xml_path_list:
    print('Processing '+xml_path+' ......')
    info_table_raw=pd.DataFrame(columns=['SRA','sequence','identity','coverage'])
    xml_file=open(xml_path)
    xml_file_content=xml_file.read()
    xml_file.close()
    query_name=re.findall(r'<query-title>.+?</query-title>',xml_file_content)[0]
    query_name=query_name.split('<query-title>')[1]
    query_name=query_name.split('</query-title>')[0]
    query_length=re.findall(r'<query-len>.+?</query-len>',xml_file_content)[0]
    query_length=query_length.split('<query-len>')[1]
    query_length=float(query_length.split('</query-len>')[0])
    hit_list=re.findall(r'<Hit>[\s\S]*?</Hit>',xml_file_content)
    for i in range(0,len(hit_list)):
        SRA_raw=re.findall(r'<accession>.+?</accession>',hit_list[i])[0]     #如果以后出bug先检查格式再优先修改这个位置
        SRA=SRA_raw.split('<accession>')[1]
        SRA=SRA.split('_')[0]
        sequence_raw=re.findall(r'<hseq>.+?</hseq>',hit_list[i])[0]
        sequence=sequence_raw.split('<hseq>')[1]
        sequence=sequence.split('</hseq>')[0]
        sequence=sequence.replace('*','.')
        identity_raw=re.findall(r'<identity>.+?</identity>',hit_list[i])[0]
        identity_str=identity_raw.split('<identity>')[1]
        identity=float(identity_str.split('</identity>')[0])
        identity_perc=identity/query_length
        coverage_raw=re.findall(r'<align-len>.+?</align-len>',hit_list[i])[0]
        coverage=coverage_raw.split('<align-len>')[1]
        coverage=float(coverage.split('</align-len>')[0])
        coverage_perc=coverage/query_length
        info_table_raw.loc[i,'SRA']=SRA
        info_table_raw.loc[i,'sequence']=sequence
        info_table_raw.loc[i,'identity']=identity_perc  
        info_table_raw.loc[i,'coverage']=coverage_perc
    info_table=pd.DataFrame(columns=['SRA','sequence','identity','coverage'])
    for i in range(0,len(info_table_raw)):
        for j in range(0,len(srr_list)):
            if info_table_raw.loc[i,'SRA']==srr_list[j]:
                info_sra=info_table_raw.loc[i,'SRA']
                info_sequence=info_table_raw.loc[i,'sequence']
                info_identity=info_table_raw.loc[i,'identity']
                info_coverage=info_table_raw.loc[i,'coverage']
                info_table.loc[i,'SRA']=info_sra
                info_table.loc[i,'sequence']=info_sequence
                info_table.loc[i,'identity']=info_identity
                info_table.loc[i,'coverage']=info_coverage
                break
    print(len(info_table))
    info_table.reset_index(inplace=True)
    info_table.sort_values(by=['identity','coverage'],ascending=False,inplace=True)
    info_table.drop_duplicates(subset=['SRA'],keep="first",inplace=True)
    info_table.reset_index(inplace=True)
    #info_table.to_csv(r'C:\Users\67535\Desktop\info_table_example.csv')
    #info_table.sort_values(by=['SRA'],ascending=True,inplace=True,kind='mergesort')
    flag=0
    for i in range(0,len(info_table)):
        if info_table.loc[i,'identity']<0.75 or info_table.loc[i,'coverage']<0.5:
            removed_dict[query_name]=xml_path
            break
        if i==len(info_table)-1 and info_table.loc[i,'identity']>=0.75 and info_table.loc[i,'coverage']>=0.5:
            flag=1
    if flag==1:
        fasta_output=open(output_path+'/'+query_name+'.fasta','w')
        for i in range(0,len(info_table)):
            fasta_output.write('>'+info_table.loc[i,'SRA']+'\n')
            fasta_output.write(info_table.loc[i,'sequence']+'\n')
        fasta_output.close()

removed_df=pd.DataFrame(removed_dict,index=[0]).T     
removed_df.to_csv(args.folder_path+'/tblastn_screening_output/removed_CDS_after_second_screen.csv')            
        
        
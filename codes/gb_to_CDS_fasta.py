# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 16:56:01 2022

@author: 67535
"""

file=open(r"C:\Users\67535\Desktop\NC_003197.2.gb")
gb_file=file.readlines()
file.close()

CDS_loc_list=[]

#Get all locations starting with "CDS"
for i in range(0,len(gb_file)):
    if 'CDS ' in gb_file[i]:
        #print(gb_file[i])
        CDS_loc_list.append(i)

ID_sequence_dict={}
sequence=''

for i in range(0,len(CDS_loc_list)):
    start_loc=CDS_loc_list[i]
    if i<len(CDS_loc_list)-1:
        end_loc=CDS_loc_list[i+1]
    else:
        end_loc=len(gb_file)
    for j in range(start_loc,end_loc):
        if "/protein_id=" in gb_file[j]:
            protein_ID=gb_file[j].split('="')[1]
            protein_ID=protein_ID.split('"')[0]
            #print(protein_ID)
        if "/translation=" in gb_file[j]:
            #print(gb_file[j].split('="')[1])
            translation_start=j
            for k in range(j,end_loc):
                if gb_file[k].endswith('"\n'):
                    translation_end=k
                    break
    for loc in range(translation_start,translation_end+1):
        partial_sequence=gb_file[loc].strip('\n')
        partial_sequence=partial_sequence.strip(' ')
        sequence+=partial_sequence
    sequence=sequence.split('/translation="')[1]
    sequence=sequence.split('"')[0]
    ID_sequence_dict[protein_ID]=sequence
    sequence=''  
    protein_ID=''
    translation_start=0
    translation_end=0

output_file=open(r"C:\Users\67535\Desktop\NC_003197.2_CDS.fasta",'w')
for key in ID_sequence_dict:
    output_file.writelines('>'+key+'\n')
    output_file.writelines(ID_sequence_dict[key]+'\n')
output_file.close()
    
        
            
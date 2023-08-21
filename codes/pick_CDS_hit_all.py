# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 11:24:17 2022

@author: 67535
"""

import argparse
import os
import re
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-e',dest='exclude_genome',type=int,required=False,default=0,help="Define the maximum number of genomes that can be excluded")
parser.add_argument('-b',dest='blast',type=str,required=True,help="Define the folder containing all tblastn files (outfmt = 14)")
parser.add_argument('-r',dest='reference',type=str,required=True,help="Define the input list containing all SRR IDs included in this analysis")
args=parser.parse_args()

file=open(args.reference.strip("'"))
reference_list=file.readlines()
file.close()
reference_list=list(set(reference_list))
reference_list=sorted(reference_list, key=str.casefold)  
for i in range(0,len(reference_list)):
    reference_list[i]=reference_list[i].strip('\n')

dir_path=args.blast.strip("'")
blast_file_list=[]
dir_list=os.listdir(dir_path)

for cur_file in dir_list:
    if os.path.splitext(cur_file)[1]=='.xml':
        blast_file=os.path.join(dir_path,cur_file)
        #print(blast_file)
        blast_file_list.append(blast_file)

exclude_file_path=''
print('Screening xml files......')
for xml_file_path in blast_file_list:
    xml_file=open(xml_file_path)
    xml_file_content_lines=xml_file.readlines()
    xml_file.close()
    for i in range(0,10):
        if '<xi:include href=' in xml_file_content_lines[i]:
            exclude_file_path=xml_file_path
            break
    break
for i in range(0,len(blast_file_list)):
    if exclude_file_path==blast_file_list[i]:
        blast_file_list.pop(i)
        break

work_xml_file_path_list=[[]]*(args.exclude_genome+1)   
missing_list=[[]]*(args.exclude_genome+1)   
#plot_x=[0]*(args.exclude_genome+1)                      #For plotting useful # of CDS V.S. # of genomes that can be excluded
#plot_y=[0]*(args.exclude_genome+1)
for xml_file_path in blast_file_list:
    print('Processing '+xml_file_path+' ......')
    xml_file=open(xml_file_path)
    xml_file_content=xml_file.read()
    xml_file.close()
    title_list_raw=re.findall(r'<accession>.+?</accession>',xml_file_content)  #之后再出bug就检查格式并和这个位置比对
    title_list=title_list_raw
    for i in range(0,len(title_list)):
        title_list[i]=title_list[i].split('<accession>')[1]
        title_list[i]=title_list[i].split('_')[0]
    title_list=list(set(title_list))
    title_list=sorted(title_list, key=str.casefold)
    for num in range(0,args.exclude_genome+1):
        if len(title_list)==len(reference_list)-num:
            work_xml_file_path_list[num].append(xml_file_path)
            missing_list_temp=list(set(reference_list).difference(set(title_list)))
            print(missing_list_temp)
            for j in range(0,len(missing_list_temp)):
                missing_list[num].append(missing_list_temp[j])
        missing_list[num]=sorted(list(set(missing_list[num])),key=str.casefold)
        work_xml_file_path_list[num]=sorted(list(set(work_xml_file_path_list[num])),key=str.casefold)

output_xml_file_path_list=[[]]*len(work_xml_file_path_list)
for i in range(0,len(work_xml_file_path_list)):
    if i==0:
        output_xml_file_path_list[0]=list(set(work_xml_file_path_list[0]))
    else:
        output_xml_file_path_list[i]=list(set(work_xml_file_path_list[i]).union(output_xml_file_path_list[i-1]))

output_srr_remove_list=[[]]*len(missing_list)
for i in range(0,len(missing_list)):
    if i==0:
        output_srr_remove_list[0]=list(set(missing_list[0]))
    else:
        output_srr_remove_list[i]=list(set(missing_list[i]).union(output_srr_remove_list[i-1]))

output_srr_list=[reference_list]*(args.exclude_genome+1)
for i in range(0,len(output_srr_list)):
    output_srr_list[i]=list(set(output_srr_list[i]).difference(output_srr_remove_list[i]))
                
plot_x=list(range(0,args.exclude_genome+1))
plot_y1=[0]*len(plot_x)
plot_y2=[0]*len(plot_x)

for i in plot_x:
    plot_y1[i]=len(output_xml_file_path_list[i])
    plot_y2[i]=len(output_srr_list[i])

output_path=args.blast+'/tblastn_screening_output'
if not(os.path.exists(output_path)):
    os.mkdir(output_path)

fig,ax1=plt.subplots()
ax2=ax1.twinx()
ax1.plot(plot_x,plot_y1,'s-',color='dodgerblue')
ax2.plot(plot_x,plot_y2,'o-',color='darkorange')
ax1.set_xlabel('Maximum # of genomes to be excluded in each CDS')
ax1.set_ylabel('Total # of CDS can be aligned \n after removing the genomes',color='dodgerblue')
ax2.set_ylabel('Total # of genomes after screening',color='darkorange')
ax2.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig(output_path+'/Statistics.png')
plt.close()

xml_path_out_path=output_path+'/xml_file_paths'
if not(os.path.exists(xml_path_out_path)):
    os.mkdir(xml_path_out_path)

for i in range(0,len(output_xml_file_path_list)):
    xml_out=open(xml_path_out_path+'/xml_file_path_'+str(i)+'.txt','w')
    for j in range(0,len(output_xml_file_path_list[i])):
        xml_out.write(output_xml_file_path_list[i][j]+'\n')
    xml_out.close()
    
srr_list_out_path=output_path+'/srr_lists_after_first_screening'
if not(os.path.exists(srr_list_out_path)):
    os.mkdir(srr_list_out_path)

for i in range(0,len(output_srr_list)):
    srr_out=open(srr_list_out_path+'/srr_list_after_first_screening_'+str(i)+'.txt','w')
    for j in range(0,len(output_srr_list[i])):
        srr_out.write(output_srr_list[i][j]+'\n')
    srr_out.close()


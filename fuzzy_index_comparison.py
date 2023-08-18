import editdistance
import sys
import os

path=sys.argv[1]

x=0
distance_dict={}
first_indexes={}
second_indexes={}

equivalent_indexes={}
all_indexes=set()
for line in open(path+'/indexed_reads.txt'):
    x+=1
    a=line.strip().split('\t')
    index_number1=a[1]
    index_sequence1=a[2]
    index_sequence2=a[3]
    compare_indexes={}
    try:
        first_indexes[index_sequence1][index_sequence2]=index_number1
    except:
        first_indexes[index_sequence1]={}
        first_indexes[index_sequence1][index_sequence2]=index_number1

    try:
        second_indexes[index_sequence2][index_sequence1]=index_number1
    except:
        second_indexes[index_sequence2]={}
        second_indexes[index_sequence2][index_sequence1]=index_number1



    all_indexes.add((index_sequence1,index_sequence2,index_number1))

for index in sorted(list(all_indexes),key=lambda x:int(x[2])):
    index_number=index[2]
    index_sequence1=index[0]
    index_sequence2=index[1]

    try:
        new_number=equivalent_indexes[index_number]
    except:
        new_number=index_number

    for new_index in first_indexes[index_sequence1]:
            compare_index=first_indexes[index_sequence1][new_index]
            if compare_index!=index_number:
                if len(index_sequence2)==len(new_index):
                    dis1=editdistance.eval(index_sequence2,new_index)
                    if dis1<=1:
                        
                        equivalent_indexes[compare_index]=new_number

    for new_index in second_indexes[index_sequence2]:
            compare_index=second_indexes[index_sequence2][new_index]
            if compare_index!=index_number:
                if len(index_sequence1)==len(new_index):
                    dis1=editdistance.eval(index_sequence1,new_index)
                    if dis1<=1:

                        equivalent_indexes[compare_index]=new_number
                


print(len(equivalent_indexes))

out=open(path+'/indexed_reads_mismatch.txt','w')    
for line in open(path+'/indexed_reads.txt'):
    
    a=line.strip().split('\t')            
    try:
        a[1]=equivalent_indexes[a[1]]+'M'
    except:
        pass
    for item in a:
        out.write(item+'\t')
    out.write('\n')

out.close()

os.system('sort -k2,2n -k7,7n '+path+'/indexed_reads_mismatch.txt >'+path+'//indexed_reads_mismatch_sorted.txt')

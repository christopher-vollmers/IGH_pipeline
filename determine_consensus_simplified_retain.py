
import re
import sys
import os
import time



path=sys.argv[1]+'//'
out_all=open(path+'/all','w')
Random_Nucleotides=21

Minimum_SeqPrep_Length=120
Single_Read_Length=300


def reverse_complement(sequence):
  Seq=''
  complement = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq


def determine_consensus_unequal_length(read,path):
    


    fast1=list()
    fast2=list()


    for sequence in read:
        
        seq1=sequence[0].split('~~~~')[0]
        qual1=sequence[1].split('~~~~')[0]

        seq2=sequence[0].split('~~~~')[1]
       	qual2=sequence[1].split('~~~~')[1]


        fast1.append([seq1,qual1,sequence[2],sequence[3],sequence[4],'forward'])
        fast2.append([seq2,qual2,sequence[2],sequence[3],sequence[4],'forward'])


    final_1,qual_final_1,matched_reads_1, error_1,count1,Matrix1,Length1=determine_consensus_equal_length(fast1,Single_Read_Length)
    final_2,qual_final_2,matched_reads_2, error_2,count2,Matrix2,Length2=determine_consensus_equal_length(fast2,Single_Read_Length)

    return final_1,qual_final_1,final_2,qual_final_2


def determine_consensus_known_length(read,SeqPrep):

    fast1=list()
    fast2=list()

    
    for sequence in read:

        seq1=sequence[0].split('~~~~')[0]
        qual1=sequence[1].split('~~~~')[0]

        seq2=sequence[0].split('~~~~')[1]
        qual2=sequence[1].split('~~~~')[1]

        R2=reverse_complement(seq2)

        fast1.append([seq1+' '*(SeqPrep-len(seq1)),qual1+' '*(SeqPrep-len(seq1)),sequence[2],sequence[3],sequence[4],'forward'])
        fast1.append([' '*(SeqPrep-len(seq2))+str(R2),' '*(SeqPrep-len(seq2))+qual2[::-1],sequence[2],sequence[3],sequence[4],'reverse'])
        
        out_all.write(seq1+'\n')
        out_all.write('!'*(SeqPrep-len(seq2))+str(R2)+'\n')


    final_1,qual_final_1,matched_reads_1, error_1,count1,Matrix,Length=determine_consensus_equal_length(fast1,SeqPrep)

    return final_1,qual_final_1,Matrix,Length
            

    
def determine_consensus_equal_length(reads,Length):
    Consensus_Seq=''
    Consensus_Qual=''
    not_complete=0
    total_dis=0

    Matrix={}

    for x in range(0,Length,1):

        Matrix[x]={}
        Base_Count={}
        Base_Count['A']=0
        Base_Count['T']=0
        Base_Count['C']=0
        Base_Count['G']=0
        Base_Count['N']=0
        Base_Count['-']=0
        Base_Count['Total']=0

        Base_Qual={}
        Base_Qual['A']=0
        Base_Qual['T']=0
        Base_Qual['C']=0
        Base_Qual['G']=0
        Base_Qual['N']=0
        Base_Qual['-']=0
        Base_Qual['Total']=0

        Base_Qual_List={}
        Base_Qual_List['A']=[]
        Base_Qual_List['T']=[]
        Base_Qual_List['C']=[]
        Base_Qual_List['G']=[]
        Base_Qual_List['N']=[]
        Base_Qual_List['-']=[]




        for read in reads:
            try:
                if read[0][x]!=' ':
                    Base_Count[read[0][x]]+=1
                    Base_Count['Total']+=1
                    Base_Qual[read[0][x]]+=ord(read[1][x])-33
                    Base_Qual['Total']+=ord(read[1][x])-33
                    Base_Qual_List[read[0][x]].append(ord(read[1][x])-33)
            except:
                pass

        Max=max(Base_Count['A'],Base_Count['T'],Base_Count['C'],Base_Count['G'],Base_Count['N'],Base_Count['-'])
        Total=Base_Count['Total']
        if Max>=1:
            total_dis+=Total-Max
            Winner=[]
            Scores=[]
            for Base in ['A','T','C','G','N','-']:
                Matrix[x][Base]=[]
                for qual in Base_Qual_List[Base]:
                    Matrix[x][Base].append(qual)
                if Base_Count[Base]==Max:
                    Winner.append(Base)
                    Scores.append(Base_Qual[Base])
            if len(Winner)==1:
                Base=Winner[0]
                Consensus_Seq+=Base
                Quality_Score=chr(int(round(((max(0,min(90,Base_Qual[Base]-(Base_Qual['Total']-Base_Qual[Base]))))+33),0)))
                Consensus_Qual+=Quality_Score
            else:
                for Base in Winner:
                    if Base_Qual[Base]==max(Scores):
                        Consensus_Seq+=Base
                        Quality_Score=chr(int(round(((max(0,min(90,Base_Qual[Base]-(Base_Qual['Total']-Base_Qual[Base]))))+33),0)))
                        Consensus_Qual+=Quality_Score
                        break
         
              
    count=len(reads)

    return Consensus_Seq,Consensus_Qual, reads,total_dis, count,Matrix,Length

def prelim_analysis(index,read1,path1,path2,path1_complete):
    SeqPrep=[]
    diversity=set()
    count=0
    for read in read1:

        diversity.add(read[0])
        count+=1
        bar1=read[2]
        bar2=read[3]
        SeqPrep=int(read[4])


    final_1,qual_final_1,final_2,qual_final_2=determine_consensus_unequal_length(read1,path)

    out1.write('@'+str(index)+'\n'+str(final_1)+'\n'+'+'+'\n'+str(qual_final_1)+'\n')
    out2.write('@'+str(index)+'\n'+str(final_2)+'\n'+'+'+'\n'+str(qual_final_2)+'\n')
    return path1, path2, path1_complete

def start_analysis(index,read1,path1,path2,path1_complete):
   
    SeqPrep=[]
    diversity=set()
    count=0
    for read in read1:

        diversity.add(read[0])
        count+=1
        bar1=read[2]
        bar2=read[3]
        SeqPrep=int(read[4])

    final_1,qual_final_1,Matrix,Length=determine_consensus_known_length(read1,SeqPrep)
    out_complete.write('@'+str(index)+'\n'+str(final_1)+'\n'+'+'+'\n'+str(qual_final_1)+'\n')
    qual_complete.write('@'+str(index)+'\t'+str(SeqPrep)+'\t')
#    for x in range(0,Length,1):
        
#        for Base in ['A','T','C','G','N','-']:
#             qual_complete.write(Base+':')            
#             Base_Scores=Matrix[x]
#             for Base_Score in Base_Scores[Base]:
#                 qual_complete.write(str(Base_Score)+',')
#             qual_complete.write('/')
#        qual_complete.write('\t')            
#    qual_complete.write('\n')         




outfile=path+'/consensus_'
     
out1=open(str(outfile)+'R1.txt','w')
out2=open(str(outfile)+'R2.txt','w')
qual_1=open(str(outfile)+'_qual_R1.txt','w')
qual_2=open(str(outfile)+'_qual_R2.txt','w')
out_singletons_1=open(path+'singletons_1','w')
out_singletons_2=open(path+'singletons_2','w')
out_complete=open(str(outfile)+'complete.txt','w')
qual_complete=open(str(outfile)+'qual_complete.txt','w')
path1=0
path1_complete=0
path2=0





molecule=''
counter=0


read1=[]
for line in open(path+'indexed_reads_mismatch_sorted.txt'):
    read=line
    a=read.strip().split('\t')
    index=a[1]
    if molecule=='':
        molecule=index
    bar1=a[2]
    bar2=a[3]
    base1=a[4]
    qual1=a[5]
    position1=-1

    if index==molecule:
        read1.append([base1,qual1,bar1,bar2,position1])
    else:
        prelim_analysis(molecule,read1,path1,path2,path1_complete)
        read1=[]
        read1.append([base1,qual1,bar1,bar2,position1])
        molecule=index

prelim_analysis(index,read1,path1,path2,path1_complete)

out1.close()
out2.close()

os.system('/home/bd1/Downloads/FLASH-1.2.11-Linux-x86_64/flash -r 300 -f 460 -s 50 -d '+path+' '+path+'/consensus_R1.txt '+path+'/consensus_R2.txt')

Length=0
for line in open(path+'out.extendedFrags.fastq'):
    Length+=1

in1=open(path+'out.extendedFrags.fastq','r')
counter=0
length_dict={}
while counter < Length:
    a=in1.readline().strip()    
    b=in1.readline().strip()
    c=in1.readline().strip()
    d=in1.readline().strip()
    index=a.split('@')[1]
    read_length=len(b)
    length_dict[index]=read_length
    counter+=4

for line in open(path+'indexed_reads_mismatch_sorted.txt'):
    read=line
    a=read.strip().split('\t')
    index=a[1]
    if molecule=='':
        molecule=index
    bar1=a[2]
    bar2=a[3]
    base1=a[4]
    qual1=a[5]

    try:
        position1=length_dict[index]
    except:
        position1=-1

    if index==molecule:
        read1.append([base1,qual1,bar1,bar2,position1])
    else:
        if read1[0][4]!=-1:

            start_analysis(molecule,read1,path1,path2,path1_complete)
        read1=[]
        read1.append([base1,qual1,bar1,bar2,position1])
        molecule=index

if read1[0][4]!=-1:
    start_analysis(index,read1,path1,path2,path1_complete)










        




  
















from __future__ import division
import time
import sys
import os


path=sys.argv[1]+'//'


def file_len(fileName):
        i=0
        for line in open(fileName):
            i+=1

        return i

infile1_list=[]
infile2_list=[]


print sorted(os.listdir(path))

os.system('rm '+path+'Combined_R1.fastq')
os.system('rm '+path+'Combined_R2.fastq')

for file1 in sorted(os.listdir(path)):
    if 'fastq' in file1:
            if 'gz' not in file1:
                if 'R1' in file1:
                    infile1_list.append(path+file1)
                elif 'R2' in file1:
                    infile2_list.append(path+file1)

print infile1_list
print infile2_list

if len(infile1_list)==1:
  print '1 file only'
  infile1=infile1_list[0]
  infile2=infile2_list[0]

else:
  out_reads1=open(path+'Combined_R1.fastq','w')
  for element in infile1_list:
    print element
    for line in open(element):
      out_reads1.write(line)
  out_reads1.close()

  out_reads2=open(path+'Combined_R2.fastq','w')
  for element in infile2_list:
    print element
    for line in open(element):
      out_reads2.write(line)
  out_reads2.close()

  infile1=path+'Combined_R1.fastq'
  infile2=path+'Combined_R2.fastq'


print 'processing'
print infile1
print 'and'
print infile2
print 'to'
print path+'indexed_reads.txt'

out=open(path+'indexed_reads.txt','w')

length=file_len(infile1)

compare_indexes={}
coverage_dict={}
in1=open(infile1,'r')
in2=open(infile2,'r')
start=time.time()
x=4
while x<=length+4:

        a=in1.readline()
        b=in1.readline()
        c=in1.readline()
        d=in1.readline()

        e=in2.readline()
        f=in2.readline()
        g=in2.readline()
        h=in2.readline()

        position=a.strip().split(' ')[0]
        index1=b[0:14]
        index2=f[0:14]
        read1=b[14:].strip()+'~~~~'+f[14:].strip()
        qual1=d[14:].strip()+'~~~~'+h[14:].strip()
        index=index1+index2

        try:
            UID_number=compare_indexes[index]
            coverage_dict[index]=coverage_dict[index]+1
        except:
            compare_indexes[index]=int(x/4)
            UID_number=int(x/4)
            coverage_dict[index]=1

        out.write('0\t'+str(UID_number)+'\t'+index1+'\t'+index2+'\t'+read1+'\t'+qual1+'\t'+str(int(x/4))+'\t'+str(position)+'\n')
        x+=4

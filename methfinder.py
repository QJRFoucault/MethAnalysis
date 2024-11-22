import sys
import time
import subprocess
import argparse
import os
parser = argparse.ArgumentParser(prog="Usage :python3 methfinder.py -c <yourcoveragefile> -cx <yourCXfile> -p <1|0> -o <youroutputname> ",description="The script will use the coverage file and CX file produced by bismark and return a bedfile of cytosine position fully methylated (p 1) or fully unmethylated (p 0)", epilog="" )
parser.add_argument("-b", "--bed", dest="bed",  help="bed file (make sure it use the same chromosome identifier than your GFF file)", action='store', type=str)
parser.add_argument("-cx", "--CX", dest="CX",  help="CX file from Bismark", action='store', type=str)
parser.add_argument("-c", "--cov", dest="cov",  help="coverage file from Bismark", action='store', type=str)
parser.add_argument("-o", "--output", dest="output",  help="Output directory ", action='store', type=str)
parser.add_argument("-p", "--percent", dest="percent",  help="percentage of methylation for a base to be detected (1 = 100, 0 = 0) ", action='store', type=str)

args = parser.parse_args()

cov=(("%s")%(args.cov))
cx=(("%s")%(args.CX))
p=(("%s")%(args.percent))
 
percent=0

if p=="1":
    percent==100
elif p=="0":
    percent==0
else :
    print ('Error in p parameter, please select between fully methylated (p 1) or fully unmethylated (p 0)') 
    exit()

file1 =open (cov)



print("preparing")
file2 =open (cx)
j=0
c= sum(1 for i in open (cov,'rb'))
i=0

#create list on only the location of 100% methylated
print("selection of cytosine")
temp=(("%stemp")%(args.output))
print (temp)
with open(temp,'w') as fp:
        
    for line in file1:
            j=j+1
            prog= int((j*100)/c)
            print('[ {0}{1} ] {2}%'.format('#' * prog, ' ' * int(100 - prog), prog),end='\r')
            r=line.split()
            if (int(float(r[3])) == percent  and (int(r[4])+int(r[5]) >= 10))  :
                for lines in file2:
                    t=lines.split()
                    if r[0] ==t[0] and (int(r[1])==int(t[1])):
                        i=i+1
                        fp.write(lines)
                        break
print("selection done, transformation in bedfile")    
k=0

tobed=open(temp)
for line in tobed:
    k=k+1
    prog= int((k*100)/i)
    print('[ {0}{1} ] {2}%'.format('#' * prog, ' ' * int(100 - prog), prog),end='\r')
    t=line.split()
    bed=(("%s.bed")%(temp))
    f = open(bed, 'a')
    out=("%s\t%s\t%s\n")%(t[0],(int(t[1])),(int(t[1])))
    f.write(out)

f.close
os.remove(temp)
print("complete")

import numpy
import pandas
import subprocess
import sys
import argparse
import os

parser = argparse.ArgumentParser(prog="Usage :python3 genefinder.py -b <yourbedfile> -r <referenceGFF> -o <youroutputname> ",description="The script will use the BED file and return a list of genes covering those regions", epilog="" )
parser.add_argument("-b", "--bed", dest="bed",  help="bed file (make sure it use the same chromosome identifier than your GFF file)", action='store', type=str)
parser.add_argument("-r", "--ref", dest="ref",  help="GFF file ", action='store', type=str)
parser.add_argument("-o", "--output", dest="output",  help="Output file ", action='store', type=str)
args = parser.parse_args()
#if any ("-h" in s for s in sys.argv): 
 #   print (parser.print_help())

    
prep=(("bedtools intersect -wb -a \"%s\" -b \"%s\" >./temp")%(args.bed,args.ref))
print(prep)
subprocess.call(prep,shell =True)
o=open(args.output,"w")
temp=open("./temp")
for line in temp:
    t=line.split()
    if len(t)>>3:
        if t[5] =="gene":
            out=("%s\t%s\t%s\t%s\n")%(t[0],t[1],t[2],t[11])
            o.write(out)
o.close
os.remove("./temp")
            
            
			
			

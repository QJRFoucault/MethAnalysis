import os
import sys
import subprocess

#fonctions

def searchingmultsample ( dirlist ) :
    
    samplelist=[]
    print ("ok1", dirlist)  
    analist=input("how many sample do you have ?: ")
    analist=int(analist)
    
    while (len(samplelist)< analist*2):
        if any (".paired.fq" in s for s in dirlist):
            po = [i for i in dirlist if  i.endswith(".paired.fq") ]
            print(po)
            samplelist.extend(po)
            dirlist = [i for i in dirlist if i not in po]
            c=len(samplelist)
            #print(c)

        elif any (".paired.fastq" in s for s in dirlist):
            po = [ i for i in dirlist if i.endswith(".paired.fastq") ]
            print(po)
            samplelist.extend(po)
            dirlist = [i for i in dirlist if i not in po]
            c=len(samplelist)
            #print(c)
            
        else:
            print ("Sorry I was not able able to find your sample, please check the name and/or path to the file and try again" )
            exit()
    print("ok2",dirlist,samplelist)
    input("alright")
    return (samplelist)

def nounpaired (path):
    pa=".paired"
    unp="unpaired"
    newdirlist=[]
    dirlist=[]
    trimdirlist=[]
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.endswith(".fq"):
                dirlist.append(os.path.abspath(os.path.join(root,name)))
            if name.endswith(".fastq"):
                dirlist.append(os.path.abspath(os.path.join(root,name)))
            if (unp in s for s in dirlist):
                newdirlist = [i for i in dirlist if unp not in i ]
            if (pa in s for s in dirlist):
                trimdirlist = [i for i in newdirlist if pa not in i ]
    return (newdirlist,trimdirlist)




# arguments check

if any ("-help" in s for s in sys.argv):
    print ("ok")
    quit()

if any ("/"in s for s in sys.argv):
    for i in range (1, len(sys.argv)):
        if ("/" in  sys.argv[i]):
            path= (sys.argv[i]) 
            print(path)
            save=open("/vol/storage/script/.SavefileBis.txt",'w')
            save.write(('%s\n')%(path))
            save.close()
if any ("-o"in s for s in sys.argv):
    save=open("/vol/storage/script/.SavefileBis.txt",'w')
    save.write(('%s\n')%(path))
    save.close()

samplist,triml=nounpaired(path)
samplist.sort()
triml.sort()
print ("samplist:",samplist)
print("triml:",triml)

#trimming block

trim=open("/vol/storage/script/.samples.txt",'w')
i=0
while (i<(len(triml)-1)): 
    p1=(triml[i])
    p2=(triml[i+1])
   # p=('%s\t%s\n')%(p1,p2)
    trim.write(('%s\t%s\n')%(p1,p2))
    i=i+2
trim.close()
save=open("/vol/storage/script/.SavefileBis.txt",'a')
#save.write('yep\n')
save.close()
# Commands task

didone=open("/vol/storage/script/.SavefileBis.txt",'r').read().splitlines();
path= didone[0]
start= len(didone)
if (start >0  ):
    samp=searchingmultsample(samplist)
    samp.sort()
    samp1=samp [0]
    samp2=samp [1]
    sampname=samp1.split(".",1)[0]
    print(samp,"1")
    print (sampname,"2")
    input("ready?")
task = {
    1:("#perl /vol/bigstorage/tools/autotrim/autotrim.pl -fofn /vol/bigstorage/script/.samples.txt -trim /vol/bigstorage/tools/autotrim/trimming.txt -tt 60 -log %s ")%(path),
    2:("/vol/storage/tool/Bismark/bismark -parallel 5 --genome_folder /vol/storage/data/Reference/ -o %s  -1 %s -2 %s ")% (path,samp1,samp2),
    3:("/vol/storage/tool/Bismark/deduplicate_bismark -p  --bam  %s.paired_bismark_bt2_pe.bam" )%(sampname) ,
    4:("/vol/storage/tool/Bismark/bismark_methylation_extractor -p --parallel 30 --bedGraph --CX --scaffold --cytosine_report --CX  --genome_folder /vol/storage/data/Reference/ %s.paired_bismark_bt2_pe.deduplicated.bam")%(sampname)

        }
print (task)


save=open("/vol/storage/script/.SavefileBis.txt",'a')
i=1
while (i<5):
    cmd=task[i]
    print(cmd)
    save.write(('%s\n')%(task[i]))
    subprocess.call (cmd,shell =True)
    i=i+1
save.close()


# open("/vol/bigstorage/script/.SavefileBis.txt",'r')




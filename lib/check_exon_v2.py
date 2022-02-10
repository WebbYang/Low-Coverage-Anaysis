#argv $low_cov_candidate_file
import sys

#ncbi to chr
ncbi_chr={}
for i in range(1,23):
    ncbi_chr["chr"+str(i)]='NC_0000'+("0"+str(i) if(i<10)  else str(i))
ncbi_chr['chrX']='NC_000023'
ncbi_chr['chrY']='NC_000024'
ncbi_chr['chrM']='NC_012920'
exon_region={}

for i in ncbi_chr.keys():
    exon_region[ncbi_chr[i]]=[]

#set exon flanking bp
flanking=50

#get exon region
gff = sys.argv[3] #"/home/dna/webb/low_coverage_test/samtools_sol/interim_GRCh37.p13_top_level_2017-01-13.gff3"
with open(gff,'r') as r:
    line=r.readline()
    while line:
        if(line.startswith("#")):
            pass
        else:
            info=line.split('\t')
            if(info[2]=="CDS" and info[0][:3]=="NC_"):
                # print([int(info[3])-flanking,int(info[4])+flanking],info[0].split('.')[0])
                exon_region[info[0].split('.')[0]].append([int(info[3]),int(info[4])])

        line=r.readline()


#load low coverage candidate region
candidate=[]
path=sys.argv[1]
with open(path,'r') as r:
    line=r.readline()
    while line:
        if(int(line.split('\t')[1])>10):
            info=line.split('\t')[0].split(':')
            chrom=info[0]
            pos=info[1].split('-')
            candidate.append([chrom,int(pos[0]),int(pos[1])])
        line=r.readline()


# low coverage overlap with exon region
# out_path = path.split('.')[0]+'_exon_region'
out_path = sys.argv[2]
fh = open(out_path,'w')
fh.write("chr\tstart\tend\tpos\n")

for i in candidate:
    #loop over specific chromsome exon
    for index,j in enumerate(exon_region[ncbi_chr[i[0]]]):
        if((j[0]-flanking<i[2]<j[1]+flanking) or (j[0]-flanking<i[1]<j[1]+flanking)):
            #left junction
            if(j[0]>i[2]>j[0]-flanking):
                print(i[0]+":"+str(i[1])+"-"+str(i[2]),"leftjunction"," -"+str(j[0]-i[2]))
                # write exon info
                fh.write(f"{i[0]}\t{j[0]}\t{j[1]}\t-{j[0]-i[2]}\n")
                break
            #right junction
            if(j[1]+flanking>i[1]>j[1]):
                print(i[0]+":"+str(i[1])+"-"+str(i[2]),"rightjunction"," +"+str(i[1]-j[1]))
                # write exon info
                fh.write(f"{i[0]}\t{j[0]}\t{j[1]}\t+{i[1]-j[1]}\n")
                break
            #internal region
            if((j[0]<i[2]<j[1]) or (j[0]<i[1]<j[1])):
                print(i[0]+":"+str(i[1])+"-"+str(i[2]),"internal")
                # write exon info
                fh.write(f"{i[0]}\t{j[0]}\t{j[1]}\tin\n")
                break

fh.close()

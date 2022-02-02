import subprocess
import re
import pickle


REF_ClinVar = '/home/dna/webb/low_coverage_test/samtools_sol/test_vcf/clinvar_20200316.vcf.gz'
REF_gff = ''

def call_cmd(cmd):
    '''
    Call shell command in python
    '''
    res = subprocess.Popen(args=cmd, shell=True, stdout=subprocess.PIPE).communicate()
    return res

def concat_contig(base_coverage, low_coverage=True):
    '''
    concat all low-coverage regions
    '''
    info = [item.split('\t') for item in base_coverage]
    if low_coverage:
        depth_thrd = 10
        info = [item for item in info if int(item[2])<depth_thrd]
    chr_num = info[0][0]
    start = int(info[0][1])
    old = start
    contigs = []
    restart_flag = False
    for item in info[1:]:
        if int(item[1])-old==1:
            old = int(item[1])
            restart_flag = False
        else:
            end = old
            contigs.append(f'{chr_num}:{start}-{end}')
            start = int(item[1])
            old = start
            restart_flag = True
    if restart_flag:
        contigs.append(f'{chr_num}:{start}-{start}')
    else:
        contigs.append(f'{chr_num}:{start}-{old}')
    return contigs

def low_coverage_check(bam_file):
    low_coveage_regions = []
    regions = [f'chr{i}' for i in range(1,23)] + ['chrX']
    for region in regions:
        res = call_cmd(f'samtools depth -Q 1 -r {region} {bam_file}')
        coverage = res[0].decode("utf-8").split('\n')[:-1]
        if len(coverage)>0:
            low_coveage_regions+=concat_contig(coverage)
    return low_coveage_regions

def annotate_clinvar_variant(low_coverage_candidate, annotate_file, REF_VCF=REF_ClinVar):
    with open(low_coverage_candidate, 'r') as rh:
        candidate_rows = rh.readlines()
 
    ClinVar_map = {}

    wh = open(annotate_file, 'w')
    for row in candidate_rows:
        row = row.split('\t')[0] #rstrip()
        res = call_cmd(f'bcftools view -r {row} {REF_VCF} | grep -v ^#')
        lines = res[0].decode("utf-8").split('\n')[:-1]
        for found in lines:
            if 'athogenic;' in found:             
                if row not in ClinVar_map.keys():
                    ClinVar_map[row] = []
                ClinVar_map[row].append(found)

    for k,vs in ClinVar_map.items():
        wh.write(k+'\n')
        for v in vs:
            wh.write(v+'\n')
    wh.close()

    with open(f'{annotate_file}.pickle', 'wb') as handle:
        pickle.dump(ClinVar_map, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return ClinVar_map #count_ClinVar(ClinVar_map)

def count_ClinVar(annotate_file):
    with open(f'{annotate_file}.pickle', 'rb') as rh:
        ClinVar_map = pickle.load(rh)
    single_nucleotide_variant = []
    deletion = []
    others = []

    for k,vs in ClinVar_map.items():
        for val in vs:       
            if 'Deletion' in val:
                deletion.append(val)
            elif 'single_nucleotide_variant' in val:
                single_nucleotide_variant.append(val)
            else:
                others.append(val)
    snv_num = len(set(single_nucleotide_variant))
    deletion_num = len(set(deletion))
    others_num = len(set(others))

    print(snv_num, deletion_num, others_num)

    return snv_num, deletion_num, others_num

def annotate_exon(low_coverage_candidate, annotate_file, REF_VCF=REF_gff):
    cmd = f'python3 lib/check_exon_v2.py {low_coverage_candidate} {annotate_file}.region > {annotate_file}'
    _ = call_cmd(cmd)

def low_coverage_report(low_coverage_candidate, clinvar_annotate_file, exon_annotate_file, low_coverage_file):
    # Load ClinVar low coverage region
    with open(f'{clinvar_annotate_file}.pickle', 'rb') as handle:
        ClinVar_map = pickle.load(handle) 

    # Load Exon low coverage region
    exon_map = {}
    with open(exon_annotate_file, 'r') as fh:
        for line in fh:
            info = line.split()
            exon_map[info[0]] = ' '.join(info[1:])

    to_report = {}
    wh = open(low_coverage_file, 'w')
    with open(low_coverage_candidate, 'r') as fh:
        for line in fh:
            key = line.split('\t')[0]
            if key in ClinVar_map.keys():
                #if 'Deletion' not in ClinVar_map[key]:
                info_list = []
                for row in ClinVar_map[key]:
                    if 'single_nucleotide_variant' in row:
                        rsid = row.split('\t')[2]
                        clnsig = [item[7:] for item in row.split('\t')[7].split(';') if item[:6]=='CLNSIG'][0]
                        info_list.append(f'{rsid} {clnsig}')
                if len(info_list)>0:
                    to_report[key] = [line.rstrip()] + ['clinvar P/LP: '+', '.join(info_list)]
            if key in exon_map.keys():
                if key not in to_report.keys():
                    to_report[key] = [line.rstrip()] + [exon_map[key]]
                # should check deletion region?
                else:
                    to_report[key] += [exon_map[key]]

    for v in to_report.values():
        wh.write('\t'.join(v)+'\n')
    wh.close()
    return

def drop_del(out_dir, sample_id, vcf_del=None):
    '''
    remove 1. exact match rows of candidate low coverage and sample vcf deletion
    remove 2. likely deletion
    remove 3. male chrX with no clinvar snp inside
    '''
    
    if vcf_del is None:
        del_file = f"/home/dna/webb/low_coverage_test/low_coverage_out/{sample_id}/vcf_scan/report_vcf_del.txt"
    else:
        del_file = vcf_del
    del_file2 = f"{out_dir}/{sample_id}_likely_deletion.txt"
    del_file3 = f"{out_dir}/{sample_id}_chrX_snp.txt"

    candidate_file = f"{out_dir}/{sample_id}_coverage_check.txt"
    out_file = f"{out_dir}/{sample_id}_coverage_check_filter.txt"

    # remove 1: exact match vcf deletion region
    with open(candidate_file, 'r') as fh:
        candidate_list = fh.readlines()

    del_k_list, del_k_list_check = [], []
    del_dic = {}
    with open(del_file, 'r') as fh:
        for line in fh:
            if len(line)>5:
                if ":" in line[4:6]:
                    k = line.split('\t')[0]
                    del_dic[k] = []
                else:
                    del_dic[k].append(line.rstrip())

    for candidate in candidate_list:
        candidate_key = candidate.split('\t')[0]
        if candidate_key in del_dic.keys():
            for item in del_dic[candidate_key]:
                item_info = item.split('\t')
                try:
                    len_ref, len_alt = int(len(item_info[3])), int(len(item_info[4]))
                except:
                    print(f"Error: {item_info}")
                del_chr, del_start, del_end = item_info[0],int(item_info[1])+len_alt, int(item_info[1])+len_ref-len_alt
                del_key = del_chr+':'+str(del_start)+'-'+str(del_end)
                #print(f"{key}\t{del_key}")
                if candidate_key==del_key:
                    del_k_list_check.append(candidate)
                    break
                    
    vcf_rm_num = len(del_k_list_check)
    print(f'Exact match sample vcf deletion remove: {vcf_rm_num}')

    # remove 2: likely deletion
    likely_del_rm_num = 0
    with open(del_file2, 'r') as fh:
        for line in fh:
            #if line not in del_k_list:
            if line in del_k_list_check:
                del_k_list.append(line)
            likely_del_rm_num += 1

    print(f'likely deletion: {likely_del_rm_num}')

    # remove 3
    chrX_del_rm_num = 0
    candidate_map = {item.split('\t')[0]:item for item in candidate_list}
    with open(del_file3, 'r') as fh:
        for line in fh:
            k,v = line.rstrip().split('\t')
            if v=='W/O':
                if k not in [item.split('\t')[0] for item in del_k_list]:
                    try:
                        del_k_list.append(candidate_map[k])
                    except:
                        print(f'Not found: {k}, should be benign. Check {sample_id}_clinvar')
                chrX_del_rm_num += 1

        print(f'Remove chrX : {chrX_del_rm_num}')

    for item in del_k_list:
        candidate_list.remove(item)

    with open(out_file, 'w') as fh:
        fh.writelines(candidate_list)

    # update clinvar file
    annotate_file = f"{out_dir}/{sample_id}_clinvar"
    with open(f'{annotate_file}.pickle', 'rb') as rh:
        ClinVar_map = pickle.load(rh)

    to_rm_keys = [item.split('\t')[0] for item in del_k_list]
    for item in to_rm_keys:
        ClinVar_map.pop(item, None)

    with open(f'{annotate_file}.pickle', 'wb') as handle:
        pickle.dump(ClinVar_map, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
def find_deletion(line):
    '''
    return the input line if len(ref) < len(alt)
    '''
    info_list = line.split('\t')
    if len(info_list[3])>len(info_list[4]):
        return line

def parse_vcf(file, func=None):
    '''
    Parse the vcf file with a function to screen the rows 
    '''
    result = []
    h = open(file, 'r')
    line = h.readline()
    start_flag = False
    line_cnt = 0
    while line:
        if start_flag:
            if func is not None:
                row = func(line)
                if row is not None:
                    result.append(row)
            else:
                print('No function input, return content start line #')
                return line_cnt
        else:
            if line[:3]=='#CH':
                #print('Start reading content at: %d'%line_cnt)
                start_flag = True
        line = h.readline()
        line_cnt+=1

    h.close()
    return result

def extend_region(line, length=1):
    chr_, start, end = re.split(':|-|\s', line)[:3]
    start = int(start)-length
    end = int(end)+length
    return f'{chr_}:{start}-{end}'

def likely_deletion_check(bam_dir, out_dir, sample_id):
    '''
    check boundary big drops 
    save to {sample_id}.likely_deletion.txt
    '''
    likely_deletion = []
    dp_threshold = 10 #20
    low_coverage_path = f'{out_dir}/{sample_id}_coverage_check.txt'
    with open(low_coverage_path,'r') as h:
        for item in h:
            info = item.split('\t')
            if 'clinvar P/LP' in info[2]:
                continue
            ext_len = 2
            if int(info[1])>10:
                region = extend_region(info[0], ext_len)
                cmd = f'samtools depth -Q 1 -r {region} {bam_dir}'
                result  = subprocess.Popen(args=cmd, shell=True, stdout=subprocess.PIPE).communicate(0)[0]
                #coverage = !samtools depth -Q 1 -r {region} {bam_dir}
                coverage = result.decode("utf-8").split('\n')[:-1]
                dp_region = [int(base.split('\t')[2]) for base in coverage]
                #if (dp_region[0]-dp_region[1]+dp_region[-1]-dp_region[-2])>dp_threshold:
                shift = ext_len-1
                if (dp_region[shift]-dp_region[ext_len]>dp_threshold) and (dp_region[-ext_len]-dp_region[-ext_len-1]>dp_threshold):
                    likely_deletion.append(item)
                elif (dp_region[0]-dp_region[ext_len]>dp_threshold) and (dp_region[-1]-dp_region[-ext_len-1]>dp_threshold):
                    likely_deletion.append(item)

    with open(f'{out_dir}/{sample_id}_likely_deletion.txt','w') as h:
        for item in likely_deletion:
            h.write(item)

def chrX_check(sex, bam_dir, out_dir, sample_id):
    '''
    if male, check reported chrX snp.
    return {k: Ture} or {k: False}
    True: some clinvar snp inside
    False: no any snp => throw out
    '''
    clinvar_snp_in = {}
    clinvar_chrX = {}
    clinvar_path = f'{out_dir}/{sample_id}_clinvar'
    # gender check first
    male = False if sex =='F' else True
    
    # clivar chrX snp check
    if male:
        with open(clinvar_path,'r') as h:
            for line in h:
                if line[3]=='X':
                    if line[4]==':':
                        k = line.rstrip()
                        clinvar_chrX[k] = []
                    else:
                        clinvar_chrX[k].append(line.split('\t')[1:5])
        for k,v in clinvar_chrX.items():
            #print(item)
            CHECK_FLAG = True
            for vi in v:
                if CHECK_FLAG:
                    pos = int(vi[0])
                    region = f"chrX:{pos}-{pos}"
                    cmd = f'samtools mpileup -r {region} {bam_dir}'
                    result = subprocess.Popen(args=cmd, shell=True, stdout=subprocess.PIPE)
                    bam_region = result.communicate(0)[0].decode("utf-8").split('\n')[:-1]
                    for item in bam_region[1:]:
                        ref_base = vi[2][0]
                        if check_snv(item.split('\t')[4],ref_base):
                            #print(item, ref_base)
                            clinvar_snp_in[k] = 'W snp'#True
                            CHECK_FLAG = False
                            break
            if CHECK_FLAG:
                clinvar_snp_in[k] = 'W/O'#False

    with open(f'{out_dir}/{sample_id}_chrX_snp.txt', 'w') as h:
        for k,v in clinvar_snp_in.items():
            h.write(k+'\t'+str(v)+'\n')
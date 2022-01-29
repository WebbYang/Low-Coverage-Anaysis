import argparse
import os
from lib.utils import *
# call_cmd, low_coverage_check, annotate_clinvar_variant, count_ClinVar, annotate_exon, low_coverage_report
'''
args: bam_dir, subject_id, out_dir
example: python3 low_coverage_analysis.py --id GB0001 --out GB0001

Date:20220127 to do
1. Add ClinVar ID after 'P/LP'
2. male chrX 
3. remove deletion? (already only include snv)
'''

parser = argparse.ArgumentParser()
parser.add_argument('--id',help='input subject id', required=True)
parser.add_argument('--out', help='output directory', required=True)
parser.add_argument('--sex', help='M of F', required=True)

args = parser.parse_args()
subject_id = args.id
out_dir = args.out
sex = args.sex

# input file
bam_dir = '/home/missmi/Missmi_Intelligence/Results/'
bam_file = f'{bam_dir}/{subject_id}/variant_call/{subject_id}.recaled.bam' # realigned
# output files
acmg_bam_file = f'{out_dir}/{subject_id}_acmg.bam'
low_coverage_candidate = f'{out_dir}/{subject_id}.low_coverage_candidate.txt'
annotate_clinvar_file = f'{out_dir}/{subject_id}_clinvar'
annotate_exon_file = f'{out_dir}/{subject_id}_exon'
low_coverage_file = f'{out_dir}/{subject_id}_coverage_check.txt'

# Extract bam
def extract_acmg():   
    cmd_list = [f"mkdir -p {out_dir}",
                f'extract_acmg_bam {acmg_bam_file} {bam_file}']
    for cmd in cmd_list:
        _ = call_cmd(cmd)
    
    print(f'Extract {bam_file} to {out_dir}/{subject_id}_acmg.bam done!')

# Get low-coverage region
def get_low_coverage():
    res = low_coverage_check(acmg_bam_file) #bam?
    with open(low_coverage_candidate, 'w') as wh:
        for line in res:
            start, end = line.split(':')[1].split('-')
            dist = int(end) - int(start) + 1
            wh.write(f'{line}\t{dist}\n')

    print(f'Generate {low_coverage_candidate} done!')

def annotate_region():
    #print('Start annotating ClinVar variants')   
    #_ = annotate_clinvar_variant(low_coverage_candidate, annotate_clinvar_file)
    print('Start annotating exon region') 
    _ = annotate_exon(low_coverage_candidate, annotate_exon_file)

def report():
    low_coverage_report(low_coverage_candidate, annotate_clinvar_file, annotate_exon_file, low_coverage_file)

def filter():
    likely_deletion_check(acmg_bam_file, out_dir, subject_id)
    chrX_check(sex, acmg_bam_file, out_dir, subject_id)
    drop_del(out_dir, subject_id)

def main():
    # extract_acmg()
    # get_low_coverage()
    annotate_region()
    report()
    filter()

if __name__ == "__main__":
    main()
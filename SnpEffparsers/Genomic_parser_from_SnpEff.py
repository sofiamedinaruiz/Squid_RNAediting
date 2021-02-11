
import argparse, os
import vcf
import pandas as pd
import time
import sys
import numpy as np
import pickle

# March 29, 2020 - Sofia Medina 
# Rokhsar lab - UC Berkeley
# Run Test: module load python/2.7-anaconda-4.4 && source activate Cori_new && python ~/SCRIPTS/Python_scripts/RNA_editing/Genomic_parser_from_SnpEff.py -v /global/cscratch1/sd/sofiamr/SQUID_2/VCF_files/Dpe46/Dpe46.mpile.cds.recode.snpEff.vcf -o /global/cscratch1/sd/sofiamr/SQUID_2/VCF_parser_v2/
# It takes 3 min per chromosome.


def parse_args():
    parser = argparse.ArgumentParser(description='VCF parser for RNA transcript variants')
    parser.add_argument("-o", "--output", metavar='STR', help='Directory where output files will be written', type=str, default='VCF_parser/')
    parser.add_argument('-m', '--moduleLoad', metavar='STR', help='Environment requirements', type=str, default="module load python/2.7-anaconda-4.4 && source activate Cori_new && ")
    parser.add_argument("-r", "--ranges", metavar='STR', help='Path to expected Depth ranges', type=str, default='/global/cscratch1/sd/sofiamr/SQUID_2/2sdfit.txt')
    parser.add_argument("-V", "--version", action='version', version='%(prog)s 1.0')
   
    required = parser.add_argument_group('required arguments')
    required.add_argument('-v', "--vcf", help="Path to vcf file, ideally it is the output from snpEff", type=str, default ='')
    args = parser.parse_args()
    
    if (args.vcf.endswith('.vcf') | args.vcf.endswith('.vcf.gz')):
        if os.path.exists(args.vcf):
            print "Input CVF file:", args.vcf
    else:
        print "Please provide a VCF file"
        sys.exit(2) 
    return(args)


def main():
    args = parse_args()
    vcf_file =  args.vcf
    expected_ranges_file = args.ranges
    
    output_dir = ''.join((args.output,'/')).replace('//','/')
    
    if (output_dir == 'VCF_parser/'):
        ensure_dir(output_dir)
        output_file = ''.join((output_dir,vcf_file.split('/')[-1].split('.vcf')[0],'.tab'))
        print "Output directory will be:", output_dir
    else:
        ensure_dir(output_dir)
        output_file = ''.join((output_dir,'/',vcf_file.split('/')[-1].split('.vcf')[0],'.tab'))
        print "Output directory will be:", output_dir
    
    analyze_VCF(vcf_file, output_file, expected_ranges_file) 
    
    print_time('Done running!')
    return()

def convert_to_DF(dictionary_, variable):
        Tst = pd.DataFrame.from_dict(dictionary_, orient='index')
        Tst = Tst.reset_index()
        Tst[variable] = Tst['index'].apply(lambda x: x[1])
        Tst['index'] = Tst['index'].apply(lambda x: x[0])
        Tst = Tst.pivot(columns=variable, index='index').fillna('')
        Tst = Tst.rename(columns={0:variable})
        return(Tst)
        
def parse_snpEff_vcf(vcf_file, Ranges_DP):
    print_time('Start')
    if 'gz' in vcf_file:
        vcf_reader = vcf.Reader(filename=vcf_file)
    else:
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    counter = 0

    
    ####
    variant_R = {}
    variant_N_alt = {}
    variant_Letters = {}
    variant_KIND = {}
    variant_SUBKIND = {}
    variant_Transcript_NAME = {}
    variant_Prot_eff  = {}
    variant_AA_eff  =  {}
    variant_BlosumScore = {}


    variant_GT = {}
    variant_GT_decoded = {}
    variant_AD1 = {}
    variant_AD2 = {}
    variant_DP  = {}
    variant_DP_OK  = {}
    variant_AD1_2 = {}
    variant_FR = {}
    variant_MQ = {}

    VAR_description_ANNOT = {}


    for record in vcf_reader:
        var_id = ''.join((record.CHROM, ":",str(record.POS)))
    
        var_all = {}
        if str(record.ALT) =='[None]': 
            var_all[var_id] = [record.REF]; 
        else: 
            var_all[var_id] = [record.REF] + map(str, record.ALT) 
    
        variant_MQ[var_id,''] = record.INFO['MQ']
    
        for sample_record in record.samples:
            sample_name = sample_record.sample
            GT = sample_record['GT']
            variant_GT[var_id,sample_name] = GT

            if (GT=='./.'):
                variant_GT_decoded[var_id,sample_name] =  GT  
                variant_AD1[var_id,sample_name]=0
                variant_AD2[var_id,sample_name]=0
                if (type(sample_record['AD']) is list )== True:
                    variant_DP[var_id,sample_name] = sum(np.array(sample_record['AD']))
                else:
                    variant_DP[var_id,sample_name] = sample_record['AD']
            else:
                GT_bin = (int(GT[0]),int(GT[-1]))
                if GT == '0/0':
                    variant_GT_decoded[var_id,sample_name] =  '/'.join((str(var_all[var_id][int(GT_bin[0])]),str(var_all[var_id][int(GT_bin[0])])))
                else:
                    variant_GT_decoded[var_id,sample_name] = '/'.join((str(var_all[var_id][int(GT_bin[0])]),str(var_all[var_id][int(GT_bin[1])])))
            
                
                if (type(sample_record['AD']) is list) == False:
                    variant_AD1[var_id,sample_name] =  sample_record['AD']
                    variant_AD2[var_id,sample_name] = 0
                    variant_DP[var_id,sample_name] = sample_record['AD']
                
                    variant_FR[var_id,sample_name] = (sample_record['ADF']>2) &  (sample_record['ADR'] >2)
                if (type(sample_record['AD']) is list )== True:
                    variant_AD1[var_id,sample_name] = sample_record['AD'][int(GT_bin[0])]
                    variant_AD2[var_id,sample_name] = sample_record['AD'][int(GT_bin[1])]
                    variant_DP[var_id,sample_name] = sample_record['AD'][int(GT_bin[0])] + sample_record['AD'][GT_bin[1]]
                    variant_FR[var_id,sample_name] = (sample_record['ADF'][GT_bin[0]]>2) &  (sample_record['ADR'][GT_bin[1]] >2)
                
            variant_DP_OK[var_id,sample_name] = (variant_DP[var_id,sample_name]>=Ranges_DP[sample_name][0]) & ((variant_DP[var_id,sample_name]<=Ranges_DP[sample_name][1]))
                
        if 'ANN' in record.INFO:
            if len(record.ALT)>=0:
                counter= counter+1
            
                prev_variants  = []
                for i in record.INFO['ANN']:
                    Kind, SubKind = obtain_variant_kind(i)
                    nucleotide = i.split('|')[0]
                    variant_Transcript_NAME[var_id,''] = i.split('|')[3]
                    variant_R[var_id,'Ref'] = record.REF
                    variant_N_alt[var_id,'N_Alt']  = len(record.ALT)
                    if (nucleotide in prev_variants) == False :
                        prev_variants.append(nucleotide)
                        variant_Letters[var_id, nucleotide] = nucleotide
                        variant_KIND[var_id,nucleotide] = Kind
                        variant_SUBKIND[var_id,nucleotide] = SubKind
                        if  len(i.split('|')[10])>0:
                            variant_Prot_eff[var_id,nucleotide] = i.split('|')[10].split('.')[1]
                            codon_eff = variant_Prot_eff[var_id,nucleotide].replace("ext*?",'').replace('*','Ter').replace('?','???')
                            variant_AA_eff[var_id,nucleotide] = '>'.join((aa_names[codon_eff[:3]],aa_names[codon_eff[-3:]]))
                            variant_BlosumScore[var_id,nucleotide] =  Blosum62[variant_AA_eff[var_id,nucleotide][0]][variant_AA_eff[var_id,nucleotide][2]]

        
    VAR_description_ANNOT['Nucleotide'] = variant_Letters
    VAR_description_ANNOT['Kind'] = variant_KIND
    VAR_description_ANNOT['SubKind'] = variant_SUBKIND
    VAR_description_ANNOT['Gene_ID'] = variant_Transcript_NAME
    VAR_description_ANNOT['AA_eff'] = variant_AA_eff
    VAR_description_ANNOT['Prot_eff'] = variant_Prot_eff
    VAR_description_ANNOT['Score'] = variant_BlosumScore
    VAR_description_ANNOT['Ref'] = variant_R
    VAR_description_ANNOT['N_Alt'] =variant_N_alt
    ###new:
    VAR_description_ANNOT['GT'] =variant_GT_decoded
    VAR_description_ANNOT['AD1'] =variant_AD1
    VAR_description_ANNOT['AD2'] =variant_AD2
    VAR_description_ANNOT['Dpt_FR'] =variant_FR
    VAR_description_ANNOT['DP'] =variant_DP
    VAR_description_ANNOT['MQ'] =variant_MQ
    VAR_description_ANNOT['DP_OK']  =  variant_DP_OK


    Annotation_DF = pd.DataFrame()
    for i in VAR_description_ANNOT.keys():
        if Annotation_DF.shape[0] == 0:
            Annotation_DF = convert_to_DF(VAR_description_ANNOT[i],i) 
            continue
        else:
            Annotation_DF = pd.merge(Annotation_DF,convert_to_DF(VAR_description_ANNOT[i],i) ,left_index=True, right_index=True, how='outer')
    
    print_time('End')
    return(Annotation_DF)


def obtain_variant_kind(record_annotation):
    kind = record_annotation.split('|')[1].replace('_variant','').replace('_region','').replace('_gene','')
    if kind in kind_rename.keys():
        kind = kind_rename[kind]
    note = ''
    if 'splice' in kind:
        note = kind
        note = note.replace('_',' ').replace('&',' & ')
        kind = "SJ"
    if '5_prime' in kind :
        if kind == '5_prime_UTR':
            note=''
        else:
            note = kind.replace('5_prime_UTR_','').replace('_',' ')
        kind = "5' UTR"
    if '3_prime' in kind :
        if kind == '3_prime_UTR':
            note=''
        else:
            note = kind.replace('3_prime_UTR_','').replace('_',' ')
        kind = "3' UTR"
    if 'missense' in kind:
        note = kind
        kind = "Coding"
    if 'synonymous' in kind:
        note = kind
        kind = "Coding"
    if kind == 'Rec':
        kind = "Coding"
        note = "Rec"
    if kind == 'Syn':
        kind = "Coding"
        note = "Syn"
    if "lost" in kind:
        note = kind.replace('_',' ')
        kind = "Coding"
    if "gained" in kind:
        note = kind.replace('_',' ')
        kind = "Coding"

    return(kind, note)


def analyze_VCF(vcf_file, output_file, expected_ranges_file):
    
    Ranges_DP_cds = pd.read_csv(expected_ranges_file, sep='\t', index_col=0)
    Ranges_DP = {}
    for i in Ranges_DP_cds.index.to_list():
        Ranges_DP[i] = (int(round(Ranges_DP_cds.loc[i]['min'])), int(round(Ranges_DP_cds.loc[i]['max2sd'])))
    
    Edit_Sites_df = parse_snpEff_vcf(vcf_file, Ranges_DP)
    
    ref_specie = 'cliff'
    if  ref_specie in Ranges_DP.keys():
        Edit_Sites_df =  Edit_Sites_df[(Edit_Sites_df.DP.cliff>=Ranges_DP[ref_specie][0]) & (Edit_Sites_df.DP.cliff<=Ranges_DP[ref_specie][1]) & (Edit_Sites_df['MQ']>23) &  (Edit_Sites_df['Dpt_FR'][ref_specie]==True) & (Edit_Sites_df[['DP_OK']]['DP_OK'].T.sum()>=3)]
    
    Edit_Sites_df.to_csv(output_file, sep='\t')
    print "File saved as:", output_file ,'\n\n'
    return(Edit_Sites_df)


def print_time(extra_text=''):
    t = time.localtime()
    print time.strftime("%D - %H:%M:%S", t), extra_text

    
def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
        print "made new dir: ", file_path
    return()

Blosum62 = pd.read_csv('/global/homes/s/sofiamr/IpythonNotebooks/SquidAnalysis/February_2020/blosum62.mtx', sep='\s+', index_col=0).to_dict()

aa_names = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C',
            'Glu':'E','Gln':'Q','Gly':'G','His':'H','Ile':'I',
            'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P',
            'Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
            'Ter':'*','???':'X','Nan':'X'}

kind_rename = {'missense':'Rec',
                   'synonymous':'Syn', 
                   'stop_retained':'Syn',
                   'intron':'Intron',
                   'intergenic':'intergenic', 
                   'downstream':'DS',
                   'upstream':'US'}



if __name__=='__main__':
    main()
    

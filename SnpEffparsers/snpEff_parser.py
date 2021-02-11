import argparse, os
import vcf
import pandas as pd
import time
import sys
import numpy as np
#import pickle


# Feb 7th, 2019 - Sofia Medina 
# Rokhsar lab - UC Berkeley
# To run: module load python/2.7-anaconda-4.4 && source activate Cori_new && python snpEff_parser.py -v /global/cscratch1/sd/sofiamr/SQUID_2/Dpe01-carrieRNA.md.mpile.snpEff.vcf
# It takes 3 min per chromosome.

def parse_args():
    parser = argparse.ArgumentParser(description='VCF parser for RNA transcript variants')
    parser.add_argument("-o", "--output", metavar='STR', help='Directory where output files will be written', type=str, default='VCF_parser/')
    parser.add_argument('-m', '--moduleLoad', metavar='STR', help='Environment requirements', type=str, default="module load python/2.7-anaconda-4.4 && source activate Cori_new && ")
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
    output_dir = ''.join((args.output,'/')).replace('//','/')
    
    if (output_dir == 'VCF_parser/'):
        ensure_dir(output_dir)
        output_file = ''.join((output_dir,vcf_file.split('/')[-1].split('.vcf')[0],'.tab'))
        print "Output directory will be:", output_dir
    else:
        ensure_dir(output_dir)
        output_file = ''.join((output_dir,'/',vcf_file.split('/')[-1].split('.vcf')[0],'.tab'))
        print "Output directory will be:", output_dir
    
    output_multiIndex = ''.join((output_file.split('/')[0], '/MultiIndex_header.txt'))
    analyze_VCF(vcf_file, output_file, output_multiIndex)
    
    print_time('Done running!')
    return()


def select_proper_anot(record):
    for i  in range(0,len(record.INFO['ANN'])):
        if 'A>G' in record.INFO['ANN'][i]:
            return(i)
        else:
            return(0)
    return(0)

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


def parse_record(record_annotation_, Eff_weight):
    flag = 'other'
    impact = Eff_weight[record_annotation_.split('|')[2]]
    transcript_id = record_annotation_.split('|')[3]
    coding_eff = record_annotation_.split('|')[9].replace('c.','')
    prot_eff = record_annotation_.split('|')[10].replace('p.','')
    #ADAR
    if ('T>C' in coding_eff) : flag = 'F'
    if ('A>G' in coding_eff) : flag = 'OK'
    #APOBEC
    if ('G>A' in coding_eff) : flag = 'F'
    if ('C>T' in coding_eff) : flag = 'OK'
        
    (kind, note )=   obtain_variant_kind(record_annotation_)
    return(transcript_id,coding_eff, prot_eff,flag,impact,kind, note)

def parse_snpEff_vcf(vcf_file):
    
    print_time('Start')
    if 'gz' in vcf_file:
        vcf_reader = vcf.Reader(filename=vcf_file)
    else:
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    counter = 0

    
    samples_dp = {}
    VAR_description = {}
    
    var_kinds = {}
    var_kind_info = {}
    var_kind_info_nuc= {}
    var_kind_info_prot= {}
    var_kind_info_note= {}
    var_transcript= {}
    var_impact = {}
    var_flag = {}
    var_ref = {}
    var_alt = {}
    var_depth = {}
    var_G_VAR = {}
    var_Enzyme = {}
    
    for record in vcf_reader:
        counter= counter+1
        if (str(record.ALT) == "[None]") ==False :
#            break
#        else:
            if (sum(record.INFO['DP4'])>0):
                var_id = ''.join((record.CHROM, ":",str(record.POS)))
                
                (transcript_id, coding_eff, prot_eff, flag, impact,kind, note) = ('', '','','','','','')
                if 'ANN' in record.INFO:
                    i = select_proper_anot(record) 
                    record_annotation = record.INFO['ANN'][i]
                    (transcript_id, coding_eff, prot_eff, flag, impact,kind, note) = parse_record(record_annotation, Eff_weight)
                    var_kinds[var_id] = kind
                    var_kind_info_nuc[var_id] = coding_eff
                    var_kind_info_prot[var_id] = prot_eff
                    var_kind_info_note[var_id] = note
                    var_transcript[var_id] = transcript_id
                    var_impact[var_id] = impact
                    var_flag[var_id] = flag
                
                var_ref[var_id] = record.REF
                var_alt[var_id] = record.ALT
                var_depth[var_id] = sum(record.INFO['DP4']) #record.INFO['DP']
                
                index_alt = -1
                var_Enzyme[var_id] = '.'
                var_G_VAR[var_id] = ''.join((str(record.REF[0]),">",str(record.ALT[0])))
                if record.REF == 'A': #A>G = ADAR
                    if 'G' in record.ALT:
                        index_alt =  record.ALT.index('G') + 1
                        var_G_VAR[var_id] = 'A>G'
                        var_Enzyme[var_id] = 'ADAR'
                if record.REF == 'T':  #T>C = ADAR
                    if 'C' in record.ALT:
                        index_alt = record.ALT.index('C') + 1
                        var_G_VAR[var_id] = 'T>C'
                        var_Enzyme[var_id] = 'ADAR'
                if record.REF == 'G': #'G>A' = APOBEC
                    if 'A' in record.ALT:
                        index_alt =  record.ALT.index('A') + 1
                        var_G_VAR[var_id] = 'G>A'
                        var_Enzyme[var_id] = 'APOBEC'
                if record.REF == 'C': #C>T = APOBEC
                    if 'T' in record.ALT:
                        index_alt = record.ALT.index('T') + 1
                        var_G_VAR[var_id] = 'C>T'
                        var_Enzyme[var_id] = 'APOBEC'
                        
                if (index_alt > -1 ):# & (len(record.ALT)==1):
                    for i in range(0,len(record.samples)):
                        if (record.samples[i]['AD'][0]+record.samples[i]['AD'][1])>0:
                            samples_dp[var_id,record.samples[i].sample] = [record.samples[i]['AD'][0],record.samples[i]['AD'][1]] + [record.samples[i]['AD'][index_alt]/float(record.samples[i]['AD'][0]+record.samples[i]['AD'][1])*100] 
                        else:
                            samples_dp[var_id,record.samples[i].sample] = [record.samples[i]['AD'][0],record.samples[i]['AD'][1]]  + [-1] 
                else:
                    for i in range(0,len(record.samples)):
                        if (record.samples[i]['AD'][0]+record.samples[i]['AD'][index_alt])>0:
                            samples_dp[var_id,record.samples[i].sample] = [record.samples[i]['AD'][0],record.samples[i]['AD'][index_alt]]  + [record.samples[i]['AD'][index_alt]/float(record.samples[i]['AD'][0]+record.samples[i]['AD'][index_alt])*100] 
            
    
    if 'ANN' in record.INFO:
        VAR_description['Kind'] = var_kinds
        VAR_description['SubKind'] = var_kind_info_note
        VAR_description['Nuc'] = var_kind_info_nuc
        VAR_description['Prot'] =  var_kind_info_prot
        VAR_description['Gene_ID'] = var_transcript
        VAR_description['I'] = var_impact
        VAR_description['Flag'] = var_flag
        VAR_description['Variant'] = var_G_VAR
        VAR_description['Enzyme'] = var_Enzyme
        
    VAR_description['REF'] = var_ref
    VAR_description['ALT'] = var_alt
    VAR_description['DP'] = var_depth
    print_time('End')
    
    return(samples_dp, VAR_description)

def make_df_from_snpEff_dict(samples_dp, VAR_description):
    print_time('DF start')  
    Tst = pd.DataFrame.from_dict(samples_dp, orient='index')
    Tst = Tst.reset_index()
    Tst = Tst.rename(columns={0:'Ref',1:'Alt',2:'Freq'})
    Tst['sample'] = Tst['index'].apply(lambda x: x[1])
    Tst['index'] = Tst['index'].apply(lambda x: x[0])
    Tst = Tst.set_index('index')
    Tst = Tst.pivot(columns='sample')
    
    
    for i in VAR_description.keys():
        Tst[i] = Tst.index.to_series().apply(lambda x: VAR_description[i][x])
        if 'Kind' in VAR_description.keys():
            Tst['Genic'] = Tst.Kind.apply(lambda x: 'Genic' if ( (x=='Coding') or (x=='Intron') or (x=="5' UTR")  or (x=='SJ') or (x=="3' UTR")) else 'Intergenic')
    if 'Prot' in VAR_description.keys():
        
        Tst.loc[Tst.Kind=='Coding','Prot'] = Tst[Tst.Kind=='Coding'].Prot.apply(lambda x: fix_prot_values(x))
        Tst.loc[Tst.Kind=='Coding','Penalty'] = Tst[Tst.Kind=='Coding'].Prot.apply(lambda x: add_blosum62_value(x))
        Tst.loc[(Tst.Kind=='Coding'),'AAsubs'] = Tst[(Tst.Kind=='Coding')].Prot.apply(lambda x:  '>'.join((aa_names[x[:3]],aa_names[x[-3:]])))
        #Tst['Enzyme'] = Tst.apply(lambda x: edited_by(x))

    Tst['N_Alts'] = Tst.ALT.apply(lambda x: len(x))
    Tst['DP_AgTc'] = Tst['Ref'].T.sum() + Tst['Alt'].T.sum() 
    Tst['Ref'] = Tst.Ref.fillna(0)
    Tst['Alt'] = Tst.Alt.fillna(0)
    Tst['Freq'] =  Tst.Alt.fillna(-1)
    print_time('DF End')  
    
    return(Tst)

def edited_by(x):
    if (x=='A>G') or (x=='T>C'):
        Enzyme = 'ADAR'
    if (x=='T>C') or (x=='A>G'):
        Enzyme = 'APOBEC'
    return(Enzyme)

def fix_prot_values(Prot_pair_aa):
    Prot_pair_aa = Prot_pair_aa.replace('ext*?','')
    Prot_pair_aa = Prot_pair_aa.replace('?','Nan')
    Prot_pair_aa = Prot_pair_aa.replace('*','Ter')
    return(Prot_pair_aa)
    
def add_blosum62_value(Prot_pair_aa):
    Ref_aa = aa_names[Prot_pair_aa[:3]]
    Alt_aa = aa_names[Prot_pair_aa[-3:]]
    return(int(Blosum62[Ref_aa][Alt_aa]))


def analyze_VCF(vcf_file, output_file, output_multiIndex):
    (samples_dp, snpEff_decriptors) = parse_snpEff_vcf(vcf_file)
    Edit_Sites_df = make_df_from_snpEff_dict(samples_dp, snpEff_decriptors)
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
            'Ter':'*','Nan':'X'}
kind_rename = {'missense':'Rec',
                   'synonymous':'Syn', 
                   'stop_retained':'Syn',
                   'intron':'Intron',
                   'intergenic':'intergenic', 
                   'downstream':'DS',
                   'upstream':'US'}
Eff_weight = {'HIGH':3,'MODERATE':2,'LOW':1,"MODIFIER":0}


if __name__=='__main__':
    main()

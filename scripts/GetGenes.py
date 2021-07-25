import argparse
from bioservices import QuickGO
import re

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", dest="vcf", help="Annotated VCF")
parser.add_argument("-g", "--GO", dest="GO", help="Genes and GO terms")
parser.add_argument("-o", "--output", dest="output", help="Output file")
parser.add_argument("-n", "--num_terms", dest="num_terms", help="Number of GO terms to include", type = int)
parser.add_argument("-b", "--BP", dest="BP", help="Return on Biological Process terms", type = bool, default = False)
args = parser.parse_args()

sig_snp_data = args.vcf
GO_data_file = args.GO
outfile = args.output
num_terms = args.num_terms
BP_only = args.BP

print('Retreiving Clunio Gene Model GO IDs')
GO_dict = {}
with open(GO_data_file, 'r') as GO_F:
    for line in GO_F:
        key = line.split('\t')[0]
        value = list(line.split('\t')[1].strip('\n').split(','))
        GO_dict[key] = value
print('Done!')

print('Relating GO IDs to SNPs')
sig_GO = {}
with open(sig_snp_data, 'r') as SNP_F:
    for line in SNP_F:
        if line[0] != '#':
            key = line.split('\t')[0] + '_' + line.split('\t')[1]
            sig_GO[key] = []
            for gene in set(re.findall('CLUMA_[A-Z]+[0-9]+',line)):
                info = {}
                info['gene_id'] = gene
                info['gene_name'] = {*()}
                info['effect'] = {*()}
                info['go_id'] = {*()}
                for ann in line.split(','):
                    if re.match(gene, ann.split('|')[4]):
                        info['gene_name'].add(ann.split('|')[3])
                        info['effect'].add(ann.split('|')[1])
                for id in info['gene_id']:
                    try:
                        for term in GO_dict[gene]:
                            info['go_id'].add(term)
                    except KeyError:
                        info['go_id'].add('NA')
                sig_GO[key].append(info)
print('Done!')


if BP_only:
    print('Building Dictionary of GO terms')
    g = QuickGO(verbose=False)

    all_terms_list = []
    for key in sig_GO:
        for gene in sig_GO[key]:
            for ID in gene['go_id']:
                all_terms_list.append(ID)

    all_terms = set(all_terms_list)

    all_terms_dic = {}
    for i, term in enumerate(all_terms):
        all_terms_dic[term] = g.get_go_terms(term)
        if int((i/len(all_terms))*100)%10 == 0:
            print(int((i/len(all_terms))*100),'% complete')
    print('Done!')
    print('Writing output')
    with open(outfile, 'w') as f:
        for snp in sig_GO:
            f.write('###' + snp + '###\n')
            for gene in sig_GO[snp]:
                if len(gene['gene_name']) > 0:
                    f.write('name:' + '\t'.join(gene['gene_name']) + '\n')
                    f.write('effect:' + '\t'.join(gene['effect'])+ '\n')
                    try:
                        try:
                            BPs = [term for term in gene['go_id'] if all_terms_dic[term][0]['aspect'] == 'biological_process']
                        except (KeyError,TypeError):
                            BPs = ['NA']
                        top_ids = [term[1] for term in sorted(zip([item.split(':')[1] for item in BPs],BPs), reverse=True)[:5]]
                        top_terms = [all_terms_dic[term][0]['name'] for term in top_ids]
                        f.write('terms:' + '\t'.join(top_terms)+ '\n')
                    except IndexError:
                        f.write('terms: NA'+ '\n')
                    f.write('\n')
    print('Done!')

else:
    g = QuickGO(verbose=False)
    print('Writing output')
    with open(outfile, 'w') as f:
        for snp in sig_GO:
            f.write('###' + snp + '###\n')
            for gene in sig_GO[snp]:
                if len(gene['gene_name']) > 0:
                    f.write('name:' + '\t'.join(gene['gene_name']) + '\n')
                    f.write('effect:' + '\t'.join(gene['effect'])+ '\n')
                    try:
                        top_ids = [term[1] for term in sorted(zip([item.split(':')[1] for item in gene['go_id']],gene['go_id']), reverse=True)[:5]]
                        top_terms = [g.get_go_terms(term)[0]['name'] for term in top_ids]
                        f.write('terms:' + '\t'.join(top_terms)+ '\n')
                    except IndexError:
                        f.write('terms: NA'+ '\n')
                    f.write('\n')
    print('Done!')

import argparse
import requests
from io import StringIO

parser = argparse.ArgumentParser(
   prog='parse_data',
   description='script that takes input hotspots csv file, and filters for '
   'genes that are represented above a given threshold in file for patients '
   'with pancreatic cancer'
   )
parser.add_argument('-i', '--input', metavar='i',
   help='input file for analysis from '
   'https://www.cancerhotspots.org/#/download')
parser.add_argument('-t', '--thresh', metavar='t', default=0.05,
   required=False, help='threshold for mutated gene frequency in patient '
   'population, prior to adjustment to account for gene size')

args = vars(parser.parse_args())
input_file = args['input']
thresh = float(args['thresh'])

with open(input_file, 'r') as f:
   cols = f.readline().split(',')
   # initialize dictionaries
   gene_dict = {}
   gene_count = {}
   protein_length = {}
   for line in f.readlines():
      line = line.split(',')
      gene = line[0]
      pos = line[1]
      # find number of pancreas samples from field with all organs present,
      # separated by '|'
      organs = line[12].split('|')
      for o in organs:
         organ_split = o.split(':')
         if organ_split[0] == 'pancreas':
            pan_cases = int(organ_split[2])
            pan_total = int(organ_split[1])
      # if gene already has a parsed mutation, but another position is mutated
      # in gene, add to gene_dict and increase case count in gene_count
      if gene in gene_dict and pos not in gene_dict[gene]:
         gene_dict[gene].append(pos)
         gene_count[gene] += pan_cases
      # if gene has no parsed mutations yet, add to gene_dict and initialize
      # key in gene_count with number of mutated cases associated
      elif gene not in gene_dict:
         gene_dict[gene] = [pos]
         gene_count[gene] = pan_cases
   # after parsing all genes, filter out cases accounting for protein size
   for g in gene_count:
      # only consider genes that already pass threshold
      if gene_count[g] > pan_total * thresh:
         # get amino acid sequence using kegg
         # query kegg with gene symbol to get corresponding kegg id, then
         # query kegg with specific gene id to get protein information
         api_url = 'https://rest.kegg.jp/find/hsa/{}'.format(g)
         response = requests.get(api_url)
         data = StringIO(response.text)
         for i in data.readlines():
            first_gene = i.split('\t')[1].split(',')[0]
            if(first_gene == g):
               kegg_gene = i.split('\t')[0]
         get_url = 'https://rest.kegg.jp/get/{}/aaseq'.format(kegg_gene)
         response = requests.get(get_url)
         data = StringIO(response.text)
         aa_length = 0
         # read protein sequence to calculate length
         for line in data:
            if(line[0] != '>' and len(line) == 61):
               aa_length += 60
            elif(line[0] != '>'):
               aa_length += len(line)-1
         protein_length[g] = aa_length
   # get average length of all proteins found
   average_length = sum(list(protein_length.values()))/len(list(protein_length.values()))
   output_list = []
   for g in gene_count:
      if gene_count[g] > pan_total * thresh:
         # divide threshold number of samples by average length of gene protein product
         adjusted_n = (pan_total * thresh)/average_length
         # divide number of samples with gene mutated by length of protein product
         adjusted_g = gene_count[g] / protein_length[g]
         if adjusted_g > adjusted_n:
            output_list.append(g)

# get all pathways associated with significant genes
path_dict = {}
path_names = {}
for gene in output_list:
   # query kegg with gene name to get kegg id
   api_url = 'https://rest.kegg.jp/find/hsa/{}'.format(gene)
   response = requests.get(api_url)
   data = StringIO(response.text)
   for i in data.readlines():
      first_gene = i.split('\t')[1].split(',')[0]
      if(first_gene == gene):
         kegg_gene = i.split('\t')[0]
   # query kegg with kegg id to get gene data
   get_url = 'https://rest.kegg.jp/get/{}'.format(kegg_gene)
   response = requests.get(get_url)
   data = StringIO(response.text)
   # parse gene data for pathway information, adding to dictionary
   pathways = False
   for line in data:
      if line[0:7] == 'PATHWAY':
         pathways = True
      elif line[0:7] != '       ':
         pathways = False
      if pathways:
         try:
            trimmed_line = line.split('hsa')[1]
            path = 'hsa{}'.format(trimmed_line.split('  ')[0])
            name = trimmed_line.split('  ')[1][:-1]
            if path in path_dict:
               path_dict[path] += 1
            else:
               path_dict[path] = 1
               path_names[path] = name
         except:
            print('poor data: ',line)
# sort pathway dictionary in descending order of included genes
sorted_dict = sorted(path_dict.items(), key=lambda item: item[1], reverse=True)
# print number of genes along with each pathway connected to at least one gene
print('{} total genes'.format(len(output_list)))
for path in sorted_dict:
   print('{}\t{}\t{}'.format(path[1],path[0],path_names[path[0]]))

print('\nGenes:')
# output final list of genes, formatted for pasting into query_pathway argument
for g in output_list:
   print(g,end=' ')
import argparse
import requests
from io import StringIO
from scipy.stats import hypergeom

parser = argparse.ArgumentParser(
   prog='query_pathway',
   description='script that takes input genes and pathway, and queries '
   'KEGG to '
   )
parser.add_argument('-g', '--genes', metavar='g', nargs='+',
   help='input list of genes for analysis')
parser.add_argument('-p', '--path', metavar='t',
   help='pathway id for analysis (ex. "hsa05200")')

args = vars(parser.parse_args())
genes = args['genes']
pathway_of_interest = args['path']

# total number of human protein-coding genes stored on kegg
# https://www.genome.jp/kegg-bin/show_organism?org=hsa
kegg_total = 20537

path_genes = []
for gene in genes:
   # query gene name to find kegg id, using first gene to identify
   api_url = 'https://rest.kegg.jp/find/hsa/{}'.format(gene)
   response = requests.get(api_url)
   data = StringIO(response.text)
   for i in data.readlines():
      first_gene = i.split('\t')[1].split(',')[0]
      if(first_gene == gene):
         kegg_gene = i.split('\t')[0]
   # query kegg with gene id to get gene data
   get_url = 'https://rest.kegg.jp/get/{}'.format(kegg_gene)
   response = requests.get(get_url)
   data = StringIO(response.text)
   pathways = False
   # iterate through data to get pathway data
   for line in data:
      if line[0:7] == 'PATHWAY':
         pathways = True
      elif line[0:7] != '       ':
         pathways = False
      if pathways:
         try:
            trimmed_line = line.split('hsa')[1]
            path = 'hsa{}'.format(trimmed_line.split('  ')[0])
            if path == pathway_of_interest:
               path_genes.append(kegg_gene)
         except:
            print('poor data: ',line)

# retrieve pathway data from kegg
api_url = 'https://rest.kegg.jp/get/{}'.format(pathway_of_interest)
response = requests.get(api_url)
data = StringIO(response.text)

# boolean to determine if current pathway data contains gene information
genes_section = False
# number of genes involved in pathway
total_genes = 0

# iterate through data, counting number of total genes in pathway
for line in data:
   if line[0:4] == 'GENE':
      genes_section = True
   if line[0:4] != '    ' and line[0:4] != 'GENE':
      genes_section = False
   if genes_section:
      total_genes+=1

# perform hypergeometric test
# number of sampled successes = length of path_genes dictionary
# number of total genes in population = kegg_total
# number of total successes = length of genes list
# number of total genes in sample = total_genes
print('{} genes of {} input genes in pathway'.format(len(path_genes),len(genes)))
print('{} genes in pathway of {} represented in KEGG'.format(total_genes, kegg_total))
p = hypergeom.sf(len(path_genes)-1,kegg_total,len(genes),total_genes)
print('p-value: ',p)

# generate text file for use in coloring pathway maps
# on webpage for given pathway (ex. https://www.kegg.jp/pathway/hsa05200),
# click 'plus' symbol next to color on left-hand sidebar, and paste this
# list to highlight genes from the dataset in orange
with open('pathway_{}.txt'.format(pathway_of_interest), 'w+') as f:
   for gene in path_genes:
      print('{} orange'.format(gene),file=f)

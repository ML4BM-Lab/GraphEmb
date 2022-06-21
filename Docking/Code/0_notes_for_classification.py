
#df_uni2pdb = pd.DataFrame(columns=['uniprotid', 'gene_id', 'pdbid', 'method', 'resolution', 'chains', 'keywords', 'EC_Code'])
# for uniprotid in tqdm(all_prot, desc='Retrieving info for protein'):
#     uniprotid
#     url = f'https://www.uniprot.org/uniprot/{uniprotid}.xml'
#     try:   
#         r = requests.get(url)
#         data = xmltodict.parse(r.content)
#     except:
#         print(f'error trying to parse  {uniprotid}')    
#     data_entry = data['uniprot']['entry']
#     if 'PDB' in [i.get('@type',0) for i in data_entry['dbReference']]:
#         EC_Code = '-'
#         if data_entry['protein'].get('recommendedName') and data_entry['protein'].get('recommendedName', {}).get('ecNumber'):
#             ec_entry = data_entry['protein']['recommendedName']['ecNumber']
#             if type(ec_entry) == str:
#                 EC_Code = data_entry['protein']['recommendedName']['ecNumber']
#             elif type(ec_entry) == list:
#                 if ec_entry[0] != dict:
#                     EC_Code = ec_entry[0]
#                 else:
#                     EC_Code = [element['#text'] for element in data_entry['protein']['recommendedName']['ecNumber']][0]
#             else:
#                 EC_Code = data_entry['protein']['recommendedName']['ecNumber'].get('#text')
#         if data_entry.get('gene') and data_entry.get('gene', {}).get('name'):
#             gene_name = data_entry['gene']['name']
#             if type(gene_name) is list:
#                 if type(gene_name[0]) == dict:
#                     gene_id = [text.get('#text') for text in data_entry['gene']['name']]
#                 else:
#                     gene_id = gene_name
#                     print(gene_id)
#             else:
#                 gene_id = gene_name.get('#text')
#         #keywords
#         keywords = None
#         if type(data_entry['keyword']) == list:
#             keywords = [key['#text'] for key in data_entry['keyword']]
#         # pdb information
#         for reference in data_entry['dbReference']:
#             if reference['@type'] == 'PDB':
#                 pdbid = reference['@id']
#                 props = reference['property']
#                 for prop in props:
#                     if prop['@type'] == 'method': 
#                         method = prop['@value']
#                     elif prop['@type'] == 'resolution':
#                         resolution = prop['@value']
#                     elif prop['@type'] == 'chains':
#                         chains = prop['@value']
#                 df_uni2pdb.loc[len(df_uni2pdb.index)] = [uniprotid, gene_id, pdbid, method, resolution, chains, keywords, EC_Code]

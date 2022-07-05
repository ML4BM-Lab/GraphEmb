from ast import keyword
import requests
from requests.adapters import HTTPAdapter, Retry

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


url ='https://rest.uniprot.org/uniprotkb/search?fields=accession%2Clength%2Cec%2Ckeyword&format=tsv&query=ion%20channel%20AND%20%28model_organism%3A9606%29%20AND%20%28proteins_with%3A2%29%20AND%20%28proteins_with%3A3%29%20AND%20%28proteins_with%3A6%29%20AND%20%28existence%3A1%29%20AND%20%28length%3A%5B401%20TO%20600%5D%29&size=500'

url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Clength%2Cec&format=tsv&query=GPCR%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29%20AND%20%28annotation_score%3A2%29&size=500'
interactions = {}

url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Clength%2Cec%2Ckeyword&format=tsv&query=GPCR%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29%20AND%20%28annotation_score%3A2%29&size=500'

url_all_human = ''

dict_unis = {}

for batch, total in get_batch(url):
    for line in batch.text.splitlines()[1:]:
        print('----'*3)
        uniprotid, length, ecode, keys = line.split('\t')
        uniprotid
        length
        ecode
        keys





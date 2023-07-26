import collections
collections.Callable = collections.abc.Callable

from bs4 import BeautifulSoup
import requests
import pandas as pd

def get_page(url):
    response = requests.get(url, headers = {'User-Agent': 'Custom'})
    if not response.ok:
        print('Status code:', response.status_code)
        raise Exception('Failed to load page {}'.format(url))
    page_content = response.text
    doc = BeautifulSoup(page_content, 'lxml')
    return doc

def getPDBIDs(url):
    table = str(get_page(url).find_all('table'))
    return pd.read_html(table)

def createMetadata():
    NBmeta = getPDBIDs('https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/nanobodies/?all=true')[0]
    NBmeta = NBmeta[NBmeta['Antigens'] == 'protein']
    NBmeta.to_csv('Data/Nanobody/Nanobody_metadata.csv')

    Abmeta = getPDBIDs('https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/search/?all=true')[0]
    Abmeta = Abmeta[Abmeta['Antigens'] == 'protein']
    Abmeta.to_csv('Data/Antibody/Antibody_metadata.csv')

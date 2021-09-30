import requests
from friendlylog import colored_logger as log

API = 'https://rest.ensembl.org'
OLD_API = 'https://grch37.rest.ensembl.org'


def vep(input, species='human', input_type='hgvs', GRCh37=False, refseq=False):
    req = f"{API if not GRCh37 else OLD_API}/vep/{species}/{input_type}/{input}"
    params = {
        'canonical': True,
        'hgvs': True,
        'SpliceRegion': True,
        'ccds': True,
        'tls': True,
        'xref_refseq': True
    }
    if refseq:
        params['refseq'] = True

    log.debug(f"Request: {req} {params}")
    resp = requests.get(req, headers={
        'Content-Type': 'application/json'
    }, params=params)
    if resp.ok:
        log.debug(resp.json())
        return resp.json()
    else:
        log.error(f"Response error {resp.status_code}: {resp.content}")
        resp.raise_for_status()

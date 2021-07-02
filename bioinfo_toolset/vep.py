import requests
from friendlylog import colored_logger as log

API = 'https://rest.ensembl.org'
OLD_API = 'https://grch37.rest.ensembl.org/'


def vep(input, species='human', input_type='hgvs', GRCh37=False):
    req = f"{API if not GRCh37 else OLD_API}/vep/{species}/{input_type}/{input}"
    log.debug(f"Request: {req}")
    resp = requests.get(req, headers={
        'Content-Type': 'application/json'
    }, params={
        'canonical': True
    })
    if resp.ok:
        log.debug(resp.json())
        return resp.json()
    else:
        log.error(f"Response error {resp.status_code}: {resp.content}")
        resp.raise_for_status()

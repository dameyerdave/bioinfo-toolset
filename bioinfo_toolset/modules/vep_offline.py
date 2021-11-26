from pathlib import Path
from os import makedirs, remove
from os.path import join, splitext
from tempfile import mkstemp
import json
import re
from friendlylog import colored_logger as log
from datetime import datetime as dt
from threading import Thread
from sys import stdout
from urllib.request import urlretrieve
import tarfile
import gzip
from tqdm import tqdm
from subprocess import Popen, PIPE, STDOUT
from shutil import move

import docker

IMAGE = 'ensemblorg/ensembl-vep'
VEP_DATA = join(Path.home(), 'vep_data')
DATA = '/opt/vep/.vep'
RELEASE = '104'
INSTALLED = Path(join(VEP_DATA, '.installed'))

CACHE = {
    'GRCh37': f"http://ftp.ensembl.org/pub/release-{RELEASE}/variation/indexed_vep_cache/homo_sapiens_merged_vep_{RELEASE}_GRCh37.tar.gz",
    'GRCh38': f"http://ftp.ensembl.org/pub/release-{RELEASE}/variation/indexed_vep_cache/homo_sapiens_merged_vep_{RELEASE}_GRCh38.tar.gz",
}

# The last GRCh37 fasta file is in release 75:
FASTA = {
    'GRCh37': 'http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz',
    'GRCh38': f"http://ftp.ensembl.org/pub/release-{RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
}

INDICATORS = ['\\', '|', '/', '|']
indicator_idx = 0


def next_indicator():
    global indicator_idx
    indicator_idx = (indicator_idx + 1) % len(INDICATORS)
    return INDICATORS[indicator_idx]


class OfflineVep():
    def __init__(self):
        self.client = docker.client.from_env()
        makedirs(VEP_DATA, exist_ok=True)

    def __rsync_and_extract(self, url, force=False):
        # replace http with rsync
        url = re.sub(r'http(s)?://', r'rsync://', url)
        # change the location of the file
        url = re.sub(r'rsync://ftp.ensembl.org/',
                     r'rsync://ftp.ensembl.org/ensembl/', url)
        local_file = join(VEP_DATA, re.sub(r'.*/', '', url))
        if force or not Path(local_file).is_file():
            command = ['rsync', '-ahP', url, VEP_DATA]
            log.info(f"Executing '{' '.join(command)}'...")
            with Popen(command, shell=False,
                       stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as rsync:
                nl = r'\n'
                for line in rsync.stdout:
                    print(
                        f"{re.sub(r'.*/', '', url)}: {re.sub(nl, '', line)}", end='\r')
                    # stdout.flush()
            if rsync.returncode == 0:
                self.__unzip(local_file, VEP_DATA)
            else:
                raise Exception(
                    f"Error rsyncing file '{url}'': {rsync.returncode}")
        else:
            log.warning(
                f"File {local_file} already exists, skipping.")

    def __unzip(self, file, dir):
        try:
            if splitext(file)[1] == '.gz':
                log.info(f"Gunzip {file}...")
                with gzip.open(file, 'r') as gz:
                    file_content = gz.read()
                    file = splitext(file)[0]
                    new_file = join(dir, re.sub(r'.*/', '', file))
                    if not Path(new_file).is_file():
                        with open(file, 'wb') as decompressed:
                            decompressed.write(file_content)
                        move(file, new_file)
                    file = new_file
            if splitext(file)[1] == '.tar':
                log.info(f"Untar {file} to {dir}...")
                with tarfile.open(file) as tar:
                    tar.extractall(dir)
        except Exception as ex:
            raise Exception(f"Cannot unzip file {file}: {ex}")

    def __download_and_extract(self, url, force=False):
        local_file = join(VEP_DATA, re.sub(r'.*/', '', url))

        if force or not Path(local_file).is_file():
            try:
                with tqdm(desc=re.sub(r'.*/', '', url), unit='B', unit_scale=True) as pbar:
                    def report_hook(_, block_size, total_size):
                        pbar.total = total_size
                        pbar.update(block_size)
                    urlretrieve(url, local_file, report_hook)
                self.__unzip(local_file, VEP_DATA)
            except Exception as ex:
                raise Exception(
                    f"Error downloading or extracting file '{url}': {ex}")
        else:
            log.warning(
                f"File {local_file} already exists, skipping.")

    def populate_cache(self):
        if not INSTALLED.is_file():
            log.info('Populating VEP cache...')
            try:
                self.__rsync_and_extract(CACHE['GRCh37'], force=True)
                self.__rsync_and_extract(CACHE['GRCh38'], force=True)
                self.__rsync_and_extract(FASTA['GRCh37'], force=True)
                self.__rsync_and_extract(FASTA['GRCh38'], force=True)
                INSTALLED.touch()
            except Exception as ex:
                raise Exception(f"Error populating cache: {ex}")
        else:
            log.info('VEP cache is already polpulated, doing nothing.')

    def __convert_path(self, local_file):
        return local_file.replace(VEP_DATA, DATA)

    def evaluate(self, input: list, GRCh37=False, **kwargs):
        assembly = 'GRCh38' if not GRCh37 else 'GRCh37'
        (_, input_file) = mkstemp(suffix='.vcf', dir=VEP_DATA)
        with open(input_file, 'w') as inf:
            tab = '\t'
            for line in input:
                # print('--line--', f"{re.sub(r'[ ]+', tab, line)}\n")
                inf.write(f"{re.sub(r'[ ]+', tab, line)}\n")
        (_, output_file) = mkstemp(suffix='.json', dir=VEP_DATA)
        vep_args = {
            'cache': True,
            'offline': True,
            'merged': True,
            'use_given_ref': True,
            'assembly': assembly,
            'fasta': join(DATA, splitext(re.sub(r'.*/', '', FASTA[assembly]))[0]),
            'canonical': True,
            'hgvs': True,
            'hgvsg': True,
            'symbol': True,
            'xref_refseq': True,
            'af_gnomad': True,
            'input_file': self.__convert_path(input_file),
            'format': 'vcf',
            'output_file': self.__convert_path(output_file),
            'json': True,
            'force_overwrite': True,
            'dir_cache': f"{DATA}/",
            'dir_plugins': f"{DATA}/Plugins/",
            'sift': 'b',
            'uniprot': True,
            'max_af': True,
            'variant_class': True
        }
        vep_args.update(kwargs)
        command = './vep '
        command += ' '.join([f"--{key} {value if value != True else ''}" for key,
                             value in vep_args.items() if value])

        self.run(command)

        with open(output_file, 'r') as outf:
            ret = []
            for line in outf.readlines():
                # print(line)
                ret.append(json.loads(line))

        remove(input_file)
        remove(output_file)
        # if not isinstance(ret, list):
        #     ret = [ret]
        return ret

    def run(self, command):
        log.debug(f"Running command: {command}...")

        def run_container():
            container = self.client.containers.run(
                IMAGE,
                name=f"vep_{dt.now().strftime('%Y-%m-%d_%H_%M_%S')}",
                auto_remove=True,
                remove=True,
                tty=True,
                stdout=True,
                stderr=True,
                stream=True,
                detach=True,
                volumes={VEP_DATA: {'bind': DATA, 'mode': 'rw'}},
                command=command
            )
            for output in container.logs(stdout=True, stderr=True,
                                         stream=True, timestamps=True, follow=True):
                stdout.write(output.decode('utf-8'))

            container.wait()

        runner = Thread(target=run_container)
        runner.start()
        runner.join()

# pipenv install "jsonschema<4.0" is required
# from bioinfo_toolset.modules.helper import dj
import docker
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy

SEQREPO_PORT = 5051
SEQREPO_IMAGE = 'biocommons/seqrepo-rest-service'
SEQREPO_DATA = '/usr/local/share/seqrepo/'
SEQREPO_REST_SERVICE_URL = f"http://localhost:{SEQREPO_PORT}/seqrepo"


class VRS():
    instance = None

    @classmethod
    def get_instance(cls):
        if cls.instance is None:
            cls.instance = VRS()
        return cls.instance

    def __init__(self):
        self.client = docker.client.from_env()
        self.start_seqrepo_rest()
        self.seqrepo = SeqRepoRESTDataProxy(
            base_url=SEQREPO_REST_SERVICE_URL)

    def start_seqrepo_rest(self):
        CONTAINER_NAME = 'seqrepo_rest'
        container = None
        try:
            container = self.client.containers.get(CONTAINER_NAME)
            if container.attrs['State']['Status'] != 'running':
                container.start()
        except docker.errors.NotFound:
            container = self.client.containers.run(
                SEQREPO_IMAGE,
                name=CONTAINER_NAME,
                auto_remove=False,
                remove=False,
                tty=False,
                stdout=False,
                stderr=False,
                stream=False,
                detach=True,
                volumes={SEQREPO_DATA: {
                    'bind': '/usr/local/share/seqrepo/', 'mode': 'ro'}},
                ports={'5000/tcp': ('127.0.0.1', SEQREPO_PORT)}
            )
        if container is None:
            raise Exception(
                "Problem running docker container ${CONTAINER_NAME}")

    def identify(self, sequence_id: str, allele: str, start: int, end: int = None):
        # if only start is provided it is basically a one xxx change and we
        # recalculate the start/end values
        if end is None or start == end:
            end = start
            start = start - 1
        # sequence_id = 'GRCh38:19'
        # print('sequence_id', sequence_id)
        # print(self.seqrepo.translate_sequence_identifier(sequence_id, "ga4gh"))
        try:
            _sequence_id = self.seqrepo.translate_sequence_identifier(sequence_id, "ga4gh")[
                0]
            _interval = models.SimpleInterval(start=start, end=end)
            _location = models.SequenceLocation(
                sequence_id=_sequence_id, interval=_interval)
            # _location['_id'] = ga4gh_identify(_location)
            _state = models.SequenceState(sequence=allele)
            _allele = models.Allele(location=_location, state=_state)
            # dj(_allele)
        except Exception as ex:
            raise Exception(
                f"Cannot query local seqrepo server. Check if the container 'seqrepo_local' is running: {ex}")

        return ga4gh_identify(_allele)

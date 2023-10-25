from enum import Enum,auto

class MetricsGroup(Enum):
    SURF = auto()
    RSCU = auto()
    CPGX = auto()
    CPG = auto()
    DRSCU = auto()
    GERP = auto()
    SYNVEP = auto()

class GenomeReference(Enum):
    hg19 = auto()
    hg38 = auto()
    
class VariantObject():

    _chr = ''
    _pos = ''
    _alt = ''
    _ref = ''
    _genome = 'hg38'

    def __init__(self,chr:str = '',pos:str = '', ref:str = '',alt:str = '',genome:GenomeReference = GenomeReference.hg38) -> None:
        self._chr = chr
        self._pos = pos
        self._ref = ref
        self._alt = alt
        self._genome = genome
        
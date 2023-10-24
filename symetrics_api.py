import numpy as np
import pandas as pd
import sqlite3
from sqlite3 import Error
import abc
import logging
from enum import Enum,auto
from sklearn.preprocessing import StandardScaler

class MetricsGroup(Enum):
    SURF = auto()
    RSCU = auto()
    CPGX = auto()
    CPG = auto()
    DRSCU = auto()
    GERP = auto()
    SYNVEP = auto()

class ISymetrics(abc.ABC):

    @abc.abstractclassmethod
    def connect_to_database():
       pass

    @abc.abstractclassmethod
    def get_silva_score():
        pass

    @abc.abstractclassmethod
    def get_surf_score():
        pass

    @abc.abstractclassmethod
    def get_synvep_score():
        pass

    @abc.abstractclassmethod
    def get_prop_score():
        pass

class Symetrics(ISymetrics):


    _db = None
    _conn = None
    
    def __init__(self, db) -> None:
        self._db = db
        self._conn = self.connect_to_database()    
        if self._conn != None:
            logging.info(f"Connection to f{self._db} is successful")

    def connect_to_database(self):

        conn = None
        try:
            conn = sqlite3.connect(self._db)
            return conn
        except Error as e:
            logging.error(f"Connection to f{self._db} failed")

        return conn
    
    def get_silva_score(self,chr = '',pos = '', ref = '',alt = ''):
        
        silva_scores = None
        try:
            # dont forget silva is hg19
            silva_cursor = self._conn.cursor()
            silva_query = f'SELECT "#chrom" AS CHR,pos AS POS,ref AS REF,alt AS ALT,gene AS GENE,"#RSCU" AS RSCU,dRSCU,"#GERP++" AS GERP,"#CpG?" AS CPG,CpG_exon AS CPGX FROM SILVA WHERE "#chrom" = {chr} AND pos = {pos} AND ref = "{ref}" AND alt = "{alt}"'
            silva_cursor.execute(silva_query)
            silva_rows = silva_cursor.fetchall()
            silva_scores = pd.DataFrame(silva_rows)
            silva_scores = silva_scores.to_dict(orient='records')

        except Error as e:
            logging.error(f"Connection to {self._db} failed")
        
        return silva_scores

    def get_surf_score(self,chr = '',pos = '', ref = '',alt = ''):
        

        surf_scores = None
        try:
            # SURF is hg38
            surf_cursor = self._conn.cursor()
            surf_query = f'SELECT CHR,POS,REF,ALT,GENE,SURF FROM SURF WHERE CHR = {chr} AND POS = {pos} AND REF = "{ref}" AND ALT = "{alt}"'
            surf_cursor.execute(surf_query)
            surf_rows = surf_cursor.fetchall()
            surf_scores = pd.DataFrame(surf_rows)
            surf_scores = surf_scores.to_dict(orient='records')
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return surf_scores
    
    def get_synvep_score(self,chr = '',pos = '', ref = '',alt = ''):
        

        synvep_scores = None
        try:
            # SURF is hg38
            synvep_cursor = self._conn.cursor()
            synvep_query = f'SELECT chr as CHR,pos_GRCh38 as POS,ref as REF,alt as ALT, HGNC_gene_symbol as GENE,synVep as SYNVEP FROM SYNVEP WHERE chr = {chr} AND pos_GRCh38 = {pos} AND ref = "{ref}" AND alt = "{alt}"'
            synvep_cursor.execute(synvep_query)
            synvep_rows = synvep_cursor.fetchall()
            synvep_scores = pd.DataFrame(synvep_rows)
            synvep_scores = synvep_scores.to_dict(orient='records')
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return synvep_scores


    def get_prop_score(self,group = MetricsGroup.SYNVEP.name,gene = ''):
        
        scores = None
        scaler = StandardScaler()


        if group in MetricsGroup.__members__:
            match group:
                case MetricsGroup.SYNVEP.name:
                    scores = pd.read_csv(f"data/{group}_DATA.csv")
                    scaled_scores = scaler.fit_transform(scores[['z_ne']])
                    scores['scaled_z'] = scaled_scores
                    scores = scores[scores.GENE ==  gene]
                    scores = scores[['GENE','pval_ne','fdr_ne','z_ne','scaled_z']]
                    scores.columns = ['GENE','PVAL','FDR','SYMETRIC_SCORE','NORM_SYMETRIC_SCORE']
                    scores['GROUP'] = group
                    scores = scores.to_dict(orient='records')
                case MetricsGroup.SURF.name:
                    scores = pd.read_csv(f"data/{group}_DATA.csv")
                    scaled_scores = scaler.fit_transform(scores[['z']])
                    scores['scaled_z'] = scaled_scores
                    scores = scores[scores.GENES ==  gene]
                    scores = scores[['GENES','pval','fdr','z','scaled_z']]
                    scores.columns = ['GENE','PVAL','FDR','SYMETRIC_SCORE','NORM_SYMETRIC_SCORE']
                    scores['GROUP'] = group
                    scores = scores.to_dict(orient='records')
                case _:
                    scores = pd.read_csv(f"data/{group}_DATA.csv")
                    scaled_scores = scaler.fit_transform(scores[['z']])
                    scores['scaled_z'] = scaled_scores
                    scores = scores[scores.GENE ==  gene]
                    scores = scores[['GENE','pval','fdr','z','scaled_z']]
                    scores.columns = ['GENE','PVAL','FDR','SYMETRIC_SCORE','NORM_SYMETRIC_SCORE']
                    scores['GROUP'] = group
                    scores = scores.to_dict(orient='records')
        else:
            logging.error(f'Group: {group} is not valid')       
    
        return scores

if __name__ == "__main__":

    symetrics_db = Symetrics("symetrics.db")
    metrics = []
    for member_name in MetricsGroup.__members__:
        result = symetrics_db.get_prop_score(group = member_name,gene = 'A1BG')
        metrics.append(result)
    print(metrics)
    a = symetrics_db.get_silva_score('7','91763673','C','A')
    b = symetrics_db.get_surf_score('2','232536581','A','T')
    c = symetrics_db.get_surf_score('2','232536581','A','T')
    print(a)
    print(b)

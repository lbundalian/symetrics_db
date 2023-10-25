import numpy as np
import pandas as pd
import sqlite3
from sqlite3 import Error
import abc
import logging
from enum import Enum,auto
from sklearn.preprocessing import StandardScaler
from .src.datastruct import *


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
    def get_spliceai_score():
        pass

    @abc.abstractclassmethod
    def get_prop_score():
        pass

    @abc.abstractclassmethod
    def get_gnomad_data():
        pass

    @abc.abstractclassmethod
    def get_gnomad_constraints():
        pass


    @abc.abstractclassmethod
    def liftover():
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
    
    def get_silva_score(self,variant: VariantObject):

        silva_scores = None
        try:
            # dont forget silva is hg19
            silva_cursor = self._conn.cursor()
            silva_query = f'SELECT "#chrom" AS CHR,pos AS POS,ref AS REF,alt AS ALT,gene AS GENE,"#RSCU" AS RSCU,dRSCU,"#GERP++" AS GERP,"#CpG?" AS CPG,CpG_exon AS CPGX FROM SILVA WHERE "#chrom" = {variant._chr} AND pos = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
            silva_cursor.execute(silva_query)
            silva_rows = silva_cursor.fetchall()
            silva_scores = silva_rows[0]
            silva_scores = {
                "CHR": silva_scores[0],
                "POS": silva_scores[1],
                "REF": silva_scores[2],
                "ALT": silva_scores[3],
                "GENE": silva_scores[4],
                "RSCU": silva_scores[5],
                "dRSCU": silva_scores[6],
                "GERP": silva_scores[7],
                "CPG": silva_scores[8],
                "CPGX": silva_scores[9]
            }

        except Error as e:
            logging.error(f"Connection to {self._db} failed")
        
        return silva_scores

    def get_surf_score(self,variant: VariantObject):
        

        surf_scores = None
        try:
            # SURF is hg38
            surf_cursor = self._conn.cursor()
            surf_query = f'SELECT CHR,POS,REF,ALT,GENE,SURF FROM SURF WHERE CHR = {variant._chr} AND POS = {variant._pos} AND REF = "{variant._ref}" AND ALT = "{variant._alt}"'
            surf_cursor.execute(surf_query)
            surf_rows = surf_cursor.fetchall()
            surf_scores = surf_rows[0]
            surf_scores = {
                "CHR": surf_scores[0],
                "POS": surf_scores[1],
                "REF": surf_scores[2],
                "ALT": surf_scores[3],
                "SURF": surf_scores[4]

            }
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return surf_scores
    
    def get_synvep_score(self,variant: VariantObject):
        
        synvep_scores = None
        try:
            # synvep is hg38 (pos_GRCh38) abd hg19 (pos)
            synvep_cursor = self._conn.cursor()
            synvep_query = ''
            if variant._genome == GenomeReference.hg38.name:
                synvep_query = f'SELECT chr as CHR,pos_GRCh38 as POS,ref as REF,alt as ALT, HGNC_gene_symbol as GENE,synVep as SYNVEP FROM SYNVEP WHERE chr = {variant._chr} AND pos_GRCh38 = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
            elif variant._genome == GenomeReference.hg19.name:
                synvep_query = f'SELECT chr as CHR,pos as POS,ref as REF,alt as ALT, HGNC_gene_symbol as GENE,synVep as SYNVEP FROM SYNVEP WHERE chr = {variant._chr} AND pos_GRCh38 = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
            synvep_cursor.execute(synvep_query)
            synvep_rows = synvep_cursor.fetchall()
            synvep_scores = synvep_rows[0]
            synvep_scores = {
                "CHR": synvep_scores[0],
                "POS": synvep_scores[1],
                "REF": synvep_scores[2],
                "ALT": synvep_scores[3],
                "GENE": synvep_scores[4],
                "SYNVEP": synvep_scores[5]

            }
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return synvep_scores

    def get_spliceai_score(self, variant: VariantObject):
                
        spliceai_score = None
        try:
            # synvep is hg38 (pos_GRCh38) abd hg19 (pos)
            spliceai_cursor = self._conn.cursor()
            spliceai_query = f'SELECT chr as CHR,pos as POS,ref as REF,alt as ALT, INFO FROM SPLICEAI WHERE chr = {variant._chr} AND pos = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
            spliceai_cursor.execute(spliceai_query)
            spliceai_rows = spliceai_cursor.fetchall()
            spliceai_score = pd.DataFrame(spliceai_rows)
            spliceai_score.columns = ['CHR','POS','REF','ALT','INFO']
            if not spliceai_score.empty:
                vcf_header = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"
                vcf_header = vcf_header.split('|')
                spliceai_score[vcf_header] = spliceai_score['INFO'].str.split('|', expand=True)
                spliceai_score['MAX_DS'] = spliceai_score.apply(lambda row: max(row['DS_AG'],row['DS_AL'],row['DS_DG'],row['DS_DL']), axis=1)
                spliceai_score = spliceai_score[['CHR','POS','REF','ALT','MAX_DS']]
                spliceai_score = spliceai_score.to_dict(orient='records')

        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return spliceai_score


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

    def get_gnomad_data(self, variant: VariantObject):

        gnomad_conn = None
        gnomad_data = None

        if variant._genome == GenomeReference.hg19.name:
            gnomad_conn = self.connect_to_database('data/gnomad2/gnomad_db.sqlite3')
        elif variant._genome == GenomeReference.hg38.name:
            gnomad_conn = self.connect_to_database('data/gnomad3/gnomad_db.sqlite3')

        try:
            gnomad_cursor = gnomad_conn.cursor()
            gnomad_query = f'SELECT chr as CHR,pos as POS,ref as REF,alt as ALT, AC, AN, AF FROM gnomad_db WHERE chr = {variant._chr} AND pos = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
            gnomad_cursor.execute(gnomad_query)
            gnomad_rows = gnomad_cursor.fetchall()
            gnomad_data = gnomad_rows[0]
            gnomad_data = {
                "CHR": gnomad_data[0],
                "POS": gnomad_data[1],
                "REF": gnomad_data[2],
                "ALT": gnomad_data[3],
                "AC": gnomad_data[4],
                "AN": gnomad_data[5],
                "AF": gnomad_data[6]
            }

        except Error as e:
            logging.error(f"Connection to Gnomad failed")
    
        return gnomad_data

    def get_gnomad_constraints(self,data='',gene=''):
        
        gnomad_data =  None
        gnomad_data = pd.read_csv(data,sep="\t")
        gnomad_data = gnomad_data[['gene','transcript','syn_z','mis_z','lof_z','pLI']]
        gnomad_data =  gnomad_data[gnomad_data.gene == gene]
        gnomad_data = gnomad_data.to_dict(orient='records')
        return gnomad_data

    def liftover(self,variant: VariantObject):

        liftover_variant = None
        try:
            # synvep is hg38 (pos_GRCh38) abd hg19 (pos)
            
            synvep_cursor = self._conn.cursor()
            synvep_query = ''
            if variant._genome == GenomeReference.hg38:
                new_reference = GenomeReference.hg19
                synvep_query = f'SELECT chr as CHR,pos as POS,ref as REF,alt as ALT, HGNC_gene_symbol as GENE,synVep as SYNVEP FROM SYNVEP WHERE chr = {variant._chr} AND pos_GRCh38 = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
            elif variant._genome == GenomeReference.hg19:
                new_reference = GenomeReference.hg38
                synvep_query = f'SELECT chr as CHR,pos_GRCh38 as POS,ref as REF,alt as ALT, HGNC_gene_symbol as GENE,synVep as SYNVEP FROM SYNVEP WHERE chr = {variant._chr} AND pos = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
            synvep_cursor.execute(synvep_query)
            synvep_rows = synvep_cursor.fetchall()
            variant_info = synvep_rows[0]
            liftover_variant = VariantObject(
                chr=variant_info[0],
                pos=variant_info[1],
                ref=variant_info[2],
                alt=variant_info[3],
                genome=new_reference
            )
            
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return liftover_variant
        

# if __name__ == "__main__":

#     symetrics_db = Symetrics("data/symetrics.db")
#     symetrics_db.get_gnomad_constraints("data/gnomad.tsv",'A1BG')
#     # metrics = []
#     # for member_name in MetricsGroup.__members__:
#     #     result = symetrics_db.get_prop_score(group = member_name,gene = 'A1BG')
#     #     metrics.append(result)
    
#     test_variant = VariantObject(chr='7',pos='91763673',ref='C',alt='A',genome=GenomeReference.hg19)
#     test_variant_liftover = symetrics_db.liftover(test_variant)
#     #a = symetrics_db.get_silva_score('7','91763673','C','A')
#     #b = symetrics_db.get_surf_score('2','232536581','A','T')
#     c = symetrics_db.get_spliceai_score(test_variant_liftover)
#     print(c)

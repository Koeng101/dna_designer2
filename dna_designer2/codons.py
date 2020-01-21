import sqlite3
from Bio import SeqIO
import random
import re
from jsonschema import validate


class CodonDatabase:
    def __init__(self,database=None,default_organism=None,rare_cutoff=None):
        if database != None:
            self.conn = sqlite3.connect(database)
        else:
            self.conn = sqlite3.connect(':memory:')
        self.default_organism=default_organism
        if rare_cutoff is None:
            self.rare_cutoff = .1
        else:
            self.rare_cutoff=rare_cutoff
        self.build_database()
            
    codon_table_validator = {
          "$schema": "http://json-schema.org/draft-07/schema#",
          #"$id": "http://example.com/product.schema.json",
          "title": "Codon Table Schema",
          "description": "A JSON validator schema for Codon Tables",
          "type": "object",
            "properties": {
                "id": {"description": "Unique identifier for the table. Usually a GenBank ID.",
                      "type": "string"},
                "organism": {"description": "Name of organism for the table. Human readable.", "type": "string"},
                "description": {"description": "Description of the table or of the organism used to build the table. Human readable."},
                "transl_table": {"description": "Translation table that the codon table uses", "type": "integer"},
                "codons": {"description": "Each individual codon in the table", "type": "array", "items": {
                    "type": "object",
                    "description": "A single codon in a table",
                    "properties": {"codon": {"description": "Uppercase 3 letter non-degenerative DNA code", "pattern": "^[ATGC]*$", "type": "string",
                                            "maxLength":3,"minLength":3},
                                  "amino_acid": {"description": "The amino acid coded by the codon. Uppercase one latter code.", "pattern":"^[ARNDBCEQZGHILKMFPSTWYV*X]$",
                                                "type":"string", "maxLength":1,"minLength":1},
                                  "codon_count": {"description": "Count of codon occurrence in all genes in given organism or GenBank file.", "type": "integer"}
                                  }
                }}
            }
        }

    def build_database(self):
        # Database schema of the Codon Database
        CREATE_DB ="""
            CREATE TABLE IF NOT EXISTS 'amino_acids' (
                'amino_acid' TEXT PRIMARY KEY
            );

            CREATE TABLE IF NOT EXISTS 'organisms' (
                'id' TEXT PRIMARY KEY,
                'description' TEXT NOT NULL,
                'organism' TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS 'transl_tables' (
                'transl_table' INT PRIMARY KEY
            );

            CREATE TABLE IF NOT EXISTS 'codons' (
                'codon' TEXT NOT NULL,
                'codon_count' INTEGER NOT NULL DEFAULT 0,
                'amino_acid' TEXT NOT NULL,
                'transl_table' INTEGER NOT NULL,
                'organism' TEXT NOT NULL,
                FOREIGN KEY('amino_acid') REFERENCES 'amino_acids'('amino_acid'),
                FOREIGN KEY('transl_table') REFERENCES 'transl_tables'('transl_table'),
                FOREIGN KEY('organism') REFERENCES 'organisms'('id')
                PRIMARY KEY('codon', 'amino_acid', 'organism')
            );

            """
        for x in CREATE_DB.split(';'):
            self.conn.cursor().execute(x)

        amino_acids = 'ARNDBCEQZGHILKMFPSTWYV*X'
        self.conn.cursor().executemany('INSERT OR IGNORE INTO amino_acids(amino_acid) VALUES(?)', [(x,) for x in amino_acids])
        self.conn.commit()

        # transl_tables are for 3NF standardization

    def choose_codon(self,aa,organism_id=None,banned=[],gc_bias=None,rare_cutoff=None):
        if organism_id is None:
            organism_id=self.default_organism
        if self.conn.cursor().execute('SELECT * FROM organisms WHERE id = ?',(organism_id,)).fetchone() == None:
            raise ValueError('Organism_id {} not available'.format(organism_id))
        
        if rare_cutoff is None:
            rare_cutoff=self.rare_cutoff

        # Check that the amino acid is valid
        valid_aa = 'ARNDBCEQZGHILKMFPSTWYV*'
        if aa not in valid_aa:
            raise ValueError('{} not valid amino acid ({})'.format(aa,valid_aa))

        # Check that all codons aren't banned
        codons_per_aa = self.conn.cursor().execute('SELECT COUNT(*) FROM codons WHERE amino_acid = ? AND organism = ?', (aa,organism_id)).fetchone()[0]
        if len(banned) >= codons_per_aa:
            raise ValueError('Too many banned codons')

        # 1. Build a table with a certain amino acid and a certain codon without any of the banned codons
        # 2. Take cumulative probability and, using a random number input from python, select a codon
        random_select = """
            WITH codons_selective AS ( WITH codon_gc AS
            (
                   SELECT codon, -- 2. Selects a table with desired amino acid, organism, and no codons from banned list. Also gives gc_content of each codon
                          amino_acid,
                          organism,
                          codon_count,
                          ( ( Length(Replace(( Replace(codon, 'T', '') ), 'A', '')) ) * 1.0 ) / Length( codon) AS gc_content
                   FROM   (
                                 SELECT * -- 1. Removes rare codons that occur fewer than 10% of the time, though this is configurable.
                                 FROM   codons
                                 WHERE  ((
                                                      codon_count /
                                                      (
                                                               SELECT   Sum(codon_count)*1.0
                                                               FROM     codons
                                                               GROUP BY amino_acid)) > ?)) AS common_codons
                   WHERE  amino_acid = ?
                   AND    organism = ?
                   AND    codon NOT IN {})
            SELECT * -- 3. Selects for either high GC or low GC codons
            FROM   codon_gc {})
            SELECT   codon, -- Calculates cumulative probability of each codon based on occurrence, then chooses one, given a random number 0-1
                     Cast(
                            (
                            SELECT Sum(codon_count)
                            FROM   codons_selective AS c
                            WHERE  c.codon <= codons_selective.codon) AS FLOAT) /
                     (
                            SELECT Sum(codon_count)
                            FROM   codons_selective) AS probability
            FROM     codons_selective
            ORDER BY Abs(probability - Abs(Random() / 18446744073709551616 + 0.5)) -- Number comes from https://www.sqlite.org/lang_corefunc.html#random
                     ASC limit 1
            """ # https://stackoverflow.com/questions/50534961/sqlite-how-to-select-rows-based-on-probability-via-integer-value

        if gc_bias == 'GC':
            gc_bias_string = 'WHERE  gc_content = (SELECT Max(gc_content) FROM codon_gc)'
        if gc_bias == 'AT':
            gc_bias_string = 'WHERE  gc_content = (SELECT Min(gc_content) FROM codon_gc)'
        else:
            gc_bias_string = ''
        try:
            # Select a codon. First, add in a variable number of ? for banned codons, then run above statement. Fetch first codon of first result
            r = self.conn.cursor().execute(random_select.format('(' + ','.join(["?" for _ in banned]) + ')',gc_bias_string), (rare_cutoff,aa,organism_id,)+tuple(banned)).fetchone()[0]
        except Exception as e:
            # If there is a value error, it is likely that the organism is not in the database
            raise ValueError('Organism or amino acid not found in database')
        return r

    def optimize_sequence(self,protein_seq,organism_id):
        if organism_id is None:
            organism_id = self.default_organism
        return ''.join([self.choose_codon(aa,organism_id=organism_id) for aa in protein_seq])

    def build_from_record(self,record):
        c = self.conn.cursor()
        c.execute('INSERT OR IGNORE INTO organisms(id,description,organism) VALUES(?,?,?)', (record.id, record.description,record.annotations['organism']))
        self.default_organism = record.id
        for feature in record.features:
            if 'translation' in feature.qualifiers:
                translation = feature.qualifiers['translation'][0] + '*'
                seq = feature.extract(record).seq
                if 'transl_table' in feature.qualifiers:
                    transl_table = int(feature.qualifiers['transl_table'][0])
                else:
                    transl_table = 0
                c.execute('INSERT OR IGNORE INTO transl_tables(transl_table) VALUES(?)', (transl_table,))
                c.executemany('INSERT OR IGNORE INTO codons(amino_acid,codon,organism,transl_table) VALUES(?,?,?,?)', [(aa,str(seq[i*3:i*3+3]),record.id,transl_table) for i,aa in enumerate(translation)])
                c.executemany('UPDATE codons SET codon_count = codon_count + 1 WHERE amino_acid = ? AND codon = ? AND organism = ?', [(aa,str(seq[i*3:i*3+3]),record.id) for i,aa in enumerate(translation)])
        self.conn.commit()

    def build_from_genbank(self,genbank_file):
        for record in SeqIO.parse(genbank_file, "genbank"):
            self.build_from_record(record)
        return self

    def codon_to_aa(self,codon,organism_id=None):
        if organism_id is None:
            organism_id = self.default_organism
        return self.conn.cursor().execute('SELECT amino_acid FROM codons WHERE codon = ? AND organism = ? ORDER BY codon_count DESC', (codon,organism_id)).fetchone()[0]

    def dump_codon_database(self):
        return '\n'.join([line for line in self.conn.iterdump()])

    def available_codons(self,aa,banned,organism_id=None,rare_cutoff=None):
        if organism_id is None:
            organism_id = self.default_organism
        if rare_cutoff is None:
            rare_cutoff = self.rare_cutoff
        return [{'codon':x[0],'weight':x[1]} for x in self.conn.cursor().execute("""SELECT codon, 1.0*codon_count/(SELECT Sum(codon_count) FROM codons WHERE amino_acid = ?) as percent FROM codons WHERE amino_acid = ? AND organism = ? AND percent > ? AND codon NOT IN ({}) """.format('(' + ','.join(["?" for _ in banned]) + ')'), (aa,aa,organism_id,rare_cutoff)+tuple(banned)).fetchall()]

    def export_table(self,organism_id=None):
        if organism_id == None:
            organism_id = self.default_organism
            if organism_id == None:
                raise ValueError('No organism_id given')
        rows = self.conn.cursor().execute("SELECT o.id,o.description,o.organism,c.codon,c.codon_count,c.amino_acid,c.transl_table FROM organisms as o JOIN codons as c on c.organism=o.id WHERE id = ?",(organism_id,)).fetchall()
        organism_json = {'id':rows[0][0],'description':rows[0][1],'organism':rows[0][2],'transl_table':rows[0][6]}
        organism_json['codons'] = [{'codon':row[3],'codon_count':row[4],'amino_acid':row[5]} for row in rows]
        return organism_json
    
    def import_table(self,json_table):
        validate(instance=json_table,schema=self.codon_table_validator)
        transl_table = json_table['transl_table']
        organism = json_table['id']
        c = self.conn.cursor()
        c.execute('INSERT OR IGNORE INTO organisms(id,description,organism) VALUES(?,?,?)',(organism,json_table['description'],json_table['organism']))
        c.execute('INSERT OR IGNORE INTO transl_tables(transl_table) VALUES(?)', (transl_table,))
        c.executemany('INSERT OR IGNORE INTO codons(amino_acid,codon,codon_count,organism,transl_table) VALUES(?,?,?,?,?)', [(codon['amino_acid'],codon['codon'],codon['codon_count'],organism,transl_table) for codon in json_table['codons']])
        self.conn.commit()
    
class Triplet:
    def __init__(self,table,codon,organism_id=None,rare_cutoff=None):
        self.table = table
        self.codon = codon.upper()
        self.last_codon = self.codon
        self.aa = self.table.codon_to_aa(codon,organism_id)
        self.organism_id = organism_id
        self.banned = [codon]
        self.rare_cutoff = rare_cutoff
    
    def change_codon(self, gc_bias=None):
        new_codon = self.table.choose_codon(self.aa,self.organism_id,self.banned,gc_bias=gc_bias)
        self.last_codon = self.codon
        self.codon = new_codon
        self.banned.append(new_codon)
        return new_codon
    
    def available_codons(self):
        return self.table.available_codons(self.aa,self.banned,self.organism_id,self.rare_cutoff)
        
    def __str__(self):
        return self.codon


class CodingSequence:
    def __init__(self,table,sequence,organism_id=None):
        # Add in verification
        self.triplets = [Triplet(table,x,organism_id) for x in re.findall('...',sequence)]
        
    def __str__(self):
        return ''.join([str(x) for x in self.triplets])

    def aa_seq(self):
        return ''.join([x.aa for x in self.triplets])

import sqlite3
import re
import pandas as pd
import requests
from io import StringIO

def build_enzyme_database(db=':memory:',references=False,type_iis_database_url='http://rebase.neb.com/rebase/link_itype2',type_iis_reference_database_url='http://rebase.neb.com/rebase/link_type2ref'):

    type_iis_database = StringIO(requests.get(type_iis_database_url).content.decode('utf-8'))
    if references == True:
        type_iis_reference_database = requests.get(type_iis_reference_database_url).content.decode('utf-8')

    conn = sqlite3.connect(db)
    c = conn.cursor()

    CREATE_DB = """
    CREATE TABLE IF NOT EXISTS 'enzymes' (
        'name' TEXT PRIMARY KEY,
        'prototype' INTEGER NOT NULL,
        'commercial_source' TEXT,
        'recognition_site' TEXT  NOT NULL,
        'methylation_site' TEXT NOT NULL,
        FOREIGN KEY('recognition_site') REFERENCES 'recognition_sites'('site'),
        FOREIGN KEY('methylation_site') REFERENCES 'methylation_sites'('site')
    );

    CREATE TABLE IF NOT EXISTS 'recognition_sites' (
        'site' TEXT PRIMARY KEY
    );

    CREATE TABLE IF NOT EXISTS 'methylation_sites' (
        'site' TEXT PRIMARY KEY
    )
    """
    for x in CREATE_DB.split(';'):
        c.execute(x)
        
    if references == True:
        CREATE_DB_REFERENCES= """
            CREATE TABLE IF NOT EXISTS 'paper_references' (
                'reference_number' INTEGER PRIMARY KEY,
                'date' INTEGER,
                'reference_string' TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS 'authors' (
                'name' TEXT PRIMARY KEY
            );

            CREATE TABLE IF NOT EXISTS 'authors_references' (
                'author' TEXT NOT NULL,
                'paper_reference' INT NOT NULL,
                FOREIGN KEY('author') REFERENCES 'authors'('name'),
                FOREIGN KEY('paper_reference') REFERENCES 'paper_references'('reference_number'),
                PRIMARY KEY('author','paper_reference')
            );

            CREATE TABLE IF NOT EXISTS 'enzymes_references' (
                'enzyme' TEXT NOT NULL,
                'paper_reference' INT NOT NULL,
                FOREIGN KEY('enzyme') REFERENCES 'enzymes'('id'),
                FOREIGN KEY('paper_reference') REFERENCES 'paper_references'('reference_number'),
                PRIMARY KEY('enzyme','paper_reference')
            )"""
        for x in CREATE_DB_REFERENCES.split(';'):
            c.execute(x)


    if references == True:
        # Insert Reference data
        data = type_iis_reference_database.split('\n \n')
        for ref in data[3:]:

            # Parse reference number
            num = int(ref.split('.')[0].split('\n')[-1])


            s = ref.split('.',1)[1].replace('\n','').replace('\t','')
            # Parse authors
            words = s.split('(',1)[0].split(',')
            authors = [x for x in list(map(str.strip, [",".join(words[i:i+2]) for i in range(0, len(words), 2)])) if x != 'Unpublished observations.' and x!= '' and 'Patent' not in x]

            # Parse date
            if 'Patent' in s:
                date = s.strip()[-5:-1]
            elif '(' in s:
                date = s[s.find('(')+len('('):s.find('(')+5]
            else:
                date = None

            # Insert reference
            c.execute('INSERT OR IGNORE INTO paper_references(reference_number,reference_string,date) VALUES(?,?,?)', (num,s,date))
            c.executemany('INSERT OR IGNORE INTO authors(name) VALUES(?)', [(x,) for x in authors])
            c.executemany('INSERT OR IGNORE INTO authors_references(author,paper_reference) VALUES(?,?)', [(x,num) for x in authors])

            if num == 2261: # This is hardcoded, will have to be updated as new versions come out
                break

    # Insert Enzyme data
    for i,row in pd.read_csv(type_iis_database,sep='\t',skiprows=9,names=['enzyme_name','prototype','recognition_sequence','methylation_site','commercial_source','references']).iterrows():
        # Insert recognition sites or methylation sites if they aren't already in database
        c.execute('INSERT OR IGNORE INTO recognition_sites(site) VALUES(?)', (row['recognition_sequence'].strip(),))
        c.execute('INSERT OR IGNORE INTO methylation_sites(site) VALUES(?)', (str(row['methylation_site']).strip(),))

        # Insert enzyme       
        c.execute('INSERT OR IGNORE INTO enzymes(name,prototype,commercial_source,recognition_site,methylation_site) VALUES(?,?,?,?,?)', (row['enzyme_name'], 0 if type(row['prototype']) == str else 1, row['commercial_source'],row['recognition_sequence'].strip(),str(row['methylation_site']).strip()))
        if references == True:
            c.executemany('INSERT OR IGNORE INTO enzymes_references(enzyme,paper_reference) VALUES(?,?)', [(row['enzyme_name'], x) for x in row['references'].split(',')])
    conn.commit()
    
    return conn

def dna_to_regex_str(recognition_site:str):
    regex_nucleotide = { #https://opisthokonta.net/?p=549
            "B": "[CGTBSKY]",
            "D": "[AGTDRWK]",
            "H": "[ACTHMYW]",
            "K": "[GTK]",
            "M": "[ACM]",
            "N": "[ACGTBDHKMNRSVWY]",
            "R": "[AGR]",
            "S": "[CGS]",
            "V": "[ACGVMSR]",
            "W": "[ATW]",
            "Y": "[CTY]"}
    plain_recognition_site = re.sub("[^a-zA-Z]+", "", recognition_site.upper())
    regex_recognition_site = ''.join([regex_nucleotide[nucleotide] if nucleotide not in 'ATGC' else nucleotide for nucleotide in plain_recognition_site])
    return regex_recognition_site

def enzymes_list_to_regex(enzyme_database_conn,enzyme_names:list):
    return [{"enzyme":row[0],"regex_site":dna_to_regex_str(row[1])} for row in enzyme_database_conn.cursor().execute('SELECT name,recognition_site FROM enzymes WHERE name IN {}'.format('(' + ','.join(["?" for _ in enzyme_names]) + ')'),tuple(enzyme_names)).fetchall()]

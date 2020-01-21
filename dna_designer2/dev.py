import codons
c = codons.CodonDatabase('example.db')
c.optimize_sequence('MA','NC_000964.3')
c.optimize_sequence('MA','NC_000')

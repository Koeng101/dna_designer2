# dna_designer2

DNA Designer 2 is a library made to codon optimize genes for synthesis. 

Internally, codon tables are saved in a SQLite database, which is then queried to optimize genes. After optimization, genes are checked for elements that affect sequencing, such as GC content and homopolymers, which are then intelligently fixed to produce the final gene product.  

Developed by Keoni Gandall.

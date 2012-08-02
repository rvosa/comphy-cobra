# various directories within project
SCRIPT=script
LIB=lib
DATA=data
RAWDATA=$(DATA)/pythonsequences/

# raw fasta files
RAWFILES := $(wildcard $(RAWDATA)/*.raw)

# cleaned fasta files
FASTAFILES = $(patsubst %.raw,%.fas,$(RAWFILES))

# this so that perl includes $(LIB) in search path
PERL=perl -I$(LIB)

# mapping between different labels
TAXAMAP=$(DATA)/excel/taxa.csv

fasta : $(FASTAFILES)

# clean up nick's raw fasta files
$(FASTAFILES) : %.fas : %.raw
	cat $< | $(PERL) $(SCRIPT)/filter_frameshifts.pl -c $(TAXAMAP) \
        | $(PERL) $(SCRIPT)/filter_sparse_codons.pl \
        | $(PERL) $(SCRIPT)/filter_short_seqs.pl -c $(TAXAMAP) -i OPHIHANN -i ANOLCARO \
        | $(PERL) $(SCRIPT)/filter_duplicate_seqs.pl > $@
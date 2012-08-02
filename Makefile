# executables that should be on the PATH or specified here
MB=mb
PHYML=phyml

# various directories within project
SCRIPT=script
LIB=lib
DATA=data
RAWDATA=$(DATA)/pythonsequences/

# raw fasta files
RAWFILES := $(wildcard $(RAWDATA)/*.raw)

# cleaned fasta files
FASTAFILES = $(patsubst %.raw,%.fas,$(RAWFILES))

# mrbayes input files
MRBAYESFILES = $(patsubst %.fas,%.mb,$(FASTAFILES))

# mrbayes output consensus trees
MCMCTREES = $(patsubst %.mb,%.mb.con.tre,$(MRBAYESFILES))

# phyml output ml trees
PHYMLTREES = $(patsubst %.mb.con.tre,%.phylip_phyml_tree.txt,$(MCMCTREES))

# phyml trees with collapsed monophyletic taxa
COLLAPSEDTREES = $(patsubst %.phylip_phyml_tree.txt,%.ctree,$(PHYMLTREES))

# files for the paml analysis
PAMLSEQS  = $(patsubst %.ctree,%.pamlseq,$(COLLAPSEDTREES))
PAMLTREES = $(patsubst %.ctree,%.pamltree,$(COLLAPSEDTREES))
PAMLCTLS  = $(patsubst %.pamltree,%.pamlctl,$(PAMLTREES))
PAMLOUTS  = $(patsubst %.pamlctl,%.pamlout,$(PAMLCTLS))
PAMLRESULT=$(DATA)/pamlresult.txt

# this so that perl includes $(LIB) in search path
PERL=perl -I$(LIB)

# mapping between different labels
TAXAMAP=$(DATA)/excel/taxa.csv

# commands for MrBayes analysis
MRBAYESBLOCK=$(DATA)/bayesian/commandblock.mb

fasta : $(FASTAFILES)

genetrees : $(COLLAPSEDTREES)

paml : $(PAMLRESULT)

# clean up nick's raw fasta files
$(FASTAFILES) : %.fas : %.raw
	cat $< | $(PERL) $(SCRIPT)/filter_frameshifts.pl -c $(TAXAMAP) \
        | $(PERL) $(SCRIPT)/filter_sparse_codons.pl \
        | $(PERL) $(SCRIPT)/filter_short_seqs.pl -c $(TAXAMAP) -i OPHIHANN -i ANOLCARO \
        | $(PERL) $(SCRIPT)/filter_duplicate_seqs.pl > $@

# converts fasta files to nexus files for mrbayes
$(MRBAYESFILES) : %.mb : %.fas
	$(PERL) $(SCRIPT)/make_mrbayes_data.pl -i $< -c $(TAXAMAP) -s fasta -d dna -m $(MRBAYESBLOCK) > $@

# run mrbayes
$(MCMCTREES) : %.mb.con.tre : %.mb
	$(MB) $<

# run phyml
$(PHYMLTREES) : %.phylip_phyml_tree.txt : %.mb.con.tre
	$(PERL) $(SCRIPT)/fas2phylip.pl -i $*.fas -c $(TAXAMAP) > $*.phylip
	$(PERL) $(SCRIPT)/nexus2newick.pl -i $< > $*.dnd
	$(PHYML) -i $*.phylip -u $*.dnd -s BEST

# collapse monophyletic clades, which are either in-paralogs or multiple
# genbank records for the same locus and the same species
$(COLLAPSEDTREES) : %.ctree : %.phylip_phyml_tree.txt
	$(PERL) $(SCRIPT)/collapse_monophyletic.pl -i $< -f newick -l phylip \
		-c $(TAXAMAP) -s ANOLCARO -s OPHIHANN > $@

# generate trees with labeled internal nodes for PAML codeml
$(PAMLTREES) : %.pamltree : %.ctree
	$(PERL) $(SCRIPT)/make_paml_tree.pl -c $(TAXAMAP) -i $< > $@

# generate PAML's version of phylip files
$(PAMLSEQS) : %.pamlseq : %.ctree
	$(PERL) $(SCRIPT)/make_paml_seqs.pl -i $*.phylip -t $< > $@

# generate paml control files
$(PAMLCTLS) : %.pamlctl : %.pamltree
	$(PERL) $(SCRIPT)/make_paml_ctl.pl -t $< -s $*.pamlseq -o $*.pamlout > $@

# run codonml
$(PAMLOUTS) : %.pamlout : %.pamlctl
	$(CODONML) $<

# parse codonml output
$(PAMLRESULT) : $(PAMLSEQS) $(PAMLOUTS)
	$(PERL) $(SCRIPT)/parse_paml_out.pl $(PAMLOUTS) > $@
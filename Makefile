# executables that should be on the PATH or specified here
MB=mb
PHYML=phyml
CODONML=codeml
PAUP=paup

# various directories within project
SCRIPT=script
LIB=lib
DATA=data
RAWDATA=$(DATA)/pythonsequences/
SOURCETREES=$(DATA)/sourcetrees

# this so that perl includes $(LIB) in search path
PERL=perl -I$(LIB)
JAVA=java

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

# phyloxml translations of ml trees, for sdi
PHYLOXMLGENETREES = $(patsubst %.ctree,%.phyloxml,$(COLLAPSEDTREES))

# trees produced by archeopteryx
SDITREES = $(patsubst %.phyloxml,%.sdi,$(PHYLOXMLGENETREES))

# files for the paml analysis
PAMLSEQS  = $(patsubst %.ctree,%.pamlseq,$(COLLAPSEDTREES))
PAMLTREES = $(patsubst %.ctree,%.pamltree,$(COLLAPSEDTREES))
PAMLCTLS  = $(patsubst %.pamltree,%.pamlctl,$(PAMLTREES))
PAMLOUTS  = $(patsubst %.pamlctl,%.pamlout,$(PAMLCTLS))
PAMLRESULT=$(DATA)/pamlresult.txt

# nexus files for upload to datamonkey
NEXUSFILES = $(patsubst %.ctree,%.nex,$(COLLAPSEDTREES))

# variables for NCBI taxonomy tree
NCBISTEM=phyliptree
NCBITREE=$(DATA)/$(NCBISTEM).phy
NCBIMRP=$(SOURCETREES)/$(NCBISTEM).dat

# variables for supertree
SUPERTREE=$(DATA)/supertree
SUPERMRPSTEM=MRP_matrix
SUPERMRP=$(SUPERTREE)/$(SUPERMRPSTEM).nex
MRPOUTGROUP=mrp_outgroup
RATCHETSETUP=$(SUPERTREE)/setup.nex
RATCHETCOMMANDS=ratchet.nex
RATCHETCOMMANDSABS=$(SUPERTREE)/$(RATCHETCOMMANDS)
RATCHETRESULT=$(SUPERTREE)/mydata.tre
RATCHETFILES=$(SUPERTREE)/mydata.tre $(SUPERTREE)/mydata.tmp $(RATCHETCOMMANDSABS) $(SUPERMRP) $(SUPERTREE)/paupratchet.log
SPECIESPHYLOXML=$(DATA)/speciestree.xml

# file patterns for supertree
NEXMLFILES := $(wildcard $(SOURCETREES)/*.xml)
MRPMATRICES = $(patsubst %.xml,%.dat,$(NEXMLFILES))

# mapping between different labels
TAXAMAP=$(DATA)/excel/taxa.csv

# commands for MrBayes analysis
MRBAYESBLOCK=$(DATA)/bayesian/commandblock.mb

fasta : $(FASTAFILES)

genetrees : $(PHYLOXMLGENETREES)

paml : $(PAMLRESULT)

nexus : $(NEXUSFILES)

supertree : $(SPECIESPHYLOXML)

treebase :
	$(PERL) $(SCRIPT)/fetch_trees.pl -d $(SOURCETREES) -c $(TAXAMAP) -v -v -v -v

sdi : $(SPECIESPHYLOXML) $(PHYLOXMLGENETREES) $(SDITREES)

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

# generates phyloxml trees
$(PHYLOXMLGENETREES) : %.phyloxml : %.ctree
	$(PERL) $(SCRIPT)/phyloxml.pl -s $* -f newick -c $(TAXAMAP) -e ctree > $@

# generates nexus files
$(NEXUSFILES) : %.nex : %.ctree
	$(PERL) $(SCRIPT)/make_nexus.pl \
        --treefile=$< \
        --treeformat=newick \
        --labels \
        --datafile=$*.phylip \
        --dataformat=phylip \
        --datatype=dna > $@

# converts downloaded source trees to mrp matrices
$(MRPMATRICES) : %.dat : %.xml
	$(PERL) $(SCRIPT)/treebase2mrp.pl -i $< -f nexml -c $(TAXAMAP) > $@

# converts the NCBI common tree to mrp matrix
$(NCBIMRP) :
	$(PERL) $(SCRIPT)/ncbi2mrp.pl -i $(NCBITREE) -f newick -c $(TAXAMAP) -v -v -v -v > $@

# concatenates mrp matrices from treebase and ncbi common tree to nexus file
$(SUPERMRP) : $(NCBIMRP) $(MRPMATRICES)
	$(PERL) $(SCRIPT)/concat_tables.pl -d $(SOURCETREES) -c $(TAXAMAP) \
        -o $(MRPOUTGROUP) > $@

# appends command blocks to mrp nexus file
$(RATCHETCOMMANDSABS) : $(SUPERMRP)
	$(PERL) $(SCRIPT)/make_ratchet_commands.pl -s $(RATCHETSETUP) -r 200 -p 15 -m $(SUPERMRP) > $(RATCHETCOMMANDSABS)
	$(PERL) $(SCRIPT)/make_ratchet_footer.pl --constraint $(NCBITREE) \
        -f newick -o $(MRPOUTGROUP) --csv $(TAXAMAP) -r $(RATCHETCOMMANDS) \
        >> $(SUPERMRP)

# runs the parsimony ratchet
$(RATCHETRESULT) : $(RATCHETCOMMANDSABS)
	cd $(SUPERTREE) && $(PAUP) $(SUPERMRPSTEM).nex && cd -

# makes consensus species tree in phyloxml format
$(SPECIESPHYLOXML) : $(RATCHETRESULT)
	$(PERL) $(SCRIPT)/make_consensus.pl -i $(RATCHETRESULT) -c $(TAXAMAP) \
        -o $(MRPOUTGROUP) > $@

# runs the sdi analysis, this requires that forester.jar is in the $CLASSPATH
$(SDITREES) : %.sdi : %.phyloxml
	$(JAVA) org.forester.application.sdi $< $(SPECIESPHYLOXML) $@


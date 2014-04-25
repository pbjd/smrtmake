SHELL = /bin/bash
ROOTDIR = /mnt/data3/vol54/2590675/0013/Analysis_Results
QUEUE ?= secondary
GENOME_SIZE ?= 5000000
CHUNK_SIZE ?= 3
NPROC ?= 15
LOCALTMP ?= /scratch
SPLITBESTN = $(shell echo $$(( 30/$(CHUNK_SIZE)+2 )))

TASKDIRS = filter correct assemble polish

SMRTETC = $(SEYMOUR_HOME)/analysis/etc
PWD = $(shell pwd)

BAXFILES := $(wildcard $(ROOTDIR)/*.bax.h5)

CHUNKS := $(shell seq 1 $(CHUNK_SIZE))

# Chunks are maintained throughout workflow. Need to add some elasticity to handle
# larger datasets
BAXFOFNS := $(foreach c,$(CHUNKS),$(shell printf "input.chunk%03dof%03d.fofn" $(c) $(CHUNK_SIZE)))
REGFOFNS := $(BAXFOFNS:input.%=filtered_regions.%)
SUBFASTA := $(BAXFOFNS:input.%.fofn=filtered_subreads.%.fasta)
SUBLENGTHS := $(BAXFOFNS:input.%.fofn=filtered_subreads.%.lengths)
LONGFASTA := $(BAXFOFNS:input.%.fofn=filtered_longreads.%.fasta)
MAPPEDM4 := $(BAXFOFNS:input.%.fofn=seeds.%.m4)
CORRECTED := $(BAXFOFNS:input.%.fofn=corrected.%.fasta)
CMPH5 := $(BAXFOFNS:input.%.fofn=aligned_reads.%.cmp.h5)

REFERENCE = reference/sequence/reference.fasta
CUTOFF = filtered_longreads.cutoff
QUERYFOFN = filtered_subreads.fofn
MAPPEDFOFN = seeds.m4.fofn

QSUB = qsub -S /bin/bash -cwd -b y -j y -sync y -V -q $(QUEUE)

all : assembly

assembly : polished_assembly.fasta

## Assembly polishing ##

polished_assembly.fasta : aligned_reads.cmp.h5
	$(QSUB) -N polish -pe smp $(NPROC) variantCaller.py -P $(SMRTETC)/algorithm_parameters/2014-03 \
	-v -j $(NPROC) --algorithm=quiver $< -r reference/sequence/reference.fasta -o corrections.gff \
	-o $@ -o polished_assembly.fastq.gz

aligned_reads.cmp.h5 : $(CMPH5)
	$(QSUB) -N merge assertCmpH5NonEmpty.py --debug $^
	cmph5tools.py  -vv merge --outFile=$@ $^
	cmph5tools.py -vv sort --deep --inPlace $@
	#h5repack -f GZIP=1 $@ $@_TMP && mv $@_TMP $@

$(CMPH5) : aligned_reads.%.cmp.h5 : input.%.fofn filtered_regions.%.fofn $(REFERENCE)
	$(QSUB) -N res.$* -pe smp $(NPROC) pbalign.py $< reference $@ --seed=1 --minAccuracy=0.75 \
	--minLength=50 --algorithmOptions=\"-useQuality -minMatch 12 -bestn 10 -minPctIdentity 70.0\" \
	--hitPolicy=randombest --tmpDir=$(LOCALTMP) -vv --nproc=$(NPROC) --regionTable=$(word 2,$^)

$(REFERENCE) : draft_assembly.fasta
	referenceUploader --skipIndexUpdate -c -n "reference" -p . -f $<  --saw="sawriter -blt 8 -welter" --samIdx="samtools faidx"

##

## Assembly draft ##

draft_assembly.fasta : asm.finished
	@echo '#!/bin/bash' > make_draft.sh
	@echo 'tmp=$$(mktemp -d -p $(LOCALTMP))' >> make_draft.sh
	@echo -n 'tigStore -g celera-assembler.gkpStore -t celera-assembler.tigStore 1 -d properties -U ' >> make_draft.sh
	@echo '| awk '\''BEGIN{t=0}$$1=="numFrags"{if($$2>1){print t, $$2}t++}'\'' | sort -nrk2,2 > unitig.lst' >> make_draft.sh
	@echo 'tmp=$$tmp cap=$(PWD)/celera-assembler utg=$(PWD)/unitig.lst cns=$(PWD)/$@ nproc=$(NPROC) pbutgcns_wf.sh' >> make_draft.sh
	@chmod +x make_draft.sh
	$(QSUB) -pe smp $(NPROC) ./make_draft.sh

asm.finished : asm.spec asm.frg
	$(QSUB) -pe smp $(NPROC) runCA -d . -p celera-assembler -s $^
	touch asm.finished

asm.spec : corrected.fasta
	runCASpecWriter.py  -vv --bitTable=$(SMRTETC)/celeraAssembler/bitTable \
	--interactiveTmpl=$(SMRTETC)/cluster/SGE/interactive.tmpl \
	--smrtpipeRc=$(SMRTETC)/smrtpipe.rc --genomeSize=$(GENOME_SIZE) --defaultFrgMinLen=500 \
	--xCoverage=20 --ovlErrorRate=0.06 --ovlMinLen=40 --merSize=14 --corrReadsFasta=corrected.fasta \
	--specOut=$@ --sgeName=asm --gridParams="useGrid:0, scriptOnGrid:0, frgCorrOnGrid:0, ovlCorrOnGrid:0" \
	--maxSlotPerc=1 $(SMRTETC)/celeraAssembler/unitig.spec

asm.frg : $(CORRECTED)
	fastqToCA -technology sanger -type sanger -libraryname corr $(patsubst %,-reads %,$(^:.fasta=.fastq)) > $@

corrected.fasta : $(CORRECTED)
	cat $^ > $@

##

## Correction (optimizations available here) ##
corrected : $(CORRECTED) ;

# make get away with escaping in the right places
$(CORRECTED) : corrected.%.fasta : seeds.%.m4 $(MAPPEDFOFN) filtered_subreads.fasta
	@echo '#!/bin/bash' > correct.$*.sh
	@echo 'tmp=$$(mktemp -d -p $(LOCALTMP))' >> correct.$*.sh
	@echo -n 'mym4=$(PWD)/$< allm4=$(PWD)/$(word 2,$^) subreads=$(PWD)/$(word 3, $^) ' >> correct.$*.sh
	@echo -n 'bestn=24 nproc=$(NPROC) fasta=$(PWD)/$@ fastq=$(PWD)/$(@:.fasta=.fastq) ' >> correct.$*.sh
	@echo 'tmp=$$tmp pbdagcon_wf.sh' >> correct.$*.sh
	@echo 'rm -rf $$tmp' >> correct.$*.sh
	@chmod +x correct.$*.sh
	$(QSUB) -pe smp $(NPROC) ./correct.$*.sh

$(MAPPEDFOFN) : $(MAPPEDM4)
	echo $^ | sed 's/ /\n/g' > $@

filtered_subreads.fasta : $(SUBFASTA)
	cat $^ > $@
##

## Read overlap using BLASR ##
$(MAPPEDM4) : seeds.%.m4 : filtered_longreads.%.fasta $(QUERYFOFN)
	$(QSUB) -pe smp $(NPROC) blasr $(QUERYFOFN) $< -out $@ -m 4 -nproc $(NPROC) -bestn 12 -noSplitSubreads \
	-maxScore -1000 -maxLCPLength 16 -minMatch 14

$(QUERYFOFN) : $(SUBFASTA)
	echo $^ | sed 's/ /\n/g' > $@
##

## Generating the long seed reads for mapping ##
$(LONGFASTA) : filtered_longreads.%.fasta : filtered_subreads.%.fasta filtered_subreads.%.lengths $(CUTOFF)
	awk -v len=$$(cat $(CUTOFF)) '($$1 < len ){ print $$2 } ' $(word 2,$^) | fastaremove $< stdin > $@

$(CUTOFF) : $(SUBLENGTHS)
	sort -nrmk1,1 $^ | awk '{t+=$$1;if(t>=$(GENOME_SIZE)*30){print $$1;exit;}}' > $@
	
$(SUBLENGTHS) : filtered_subreads.%.lengths : filtered_subreads.%.fasta
	fastalength $< | sort -nrk1,1 > $@ 
##

## Extracting subreads (avoidable with some work) ##
subreads : $(SUBFASTA) ;

$(SUBFASTA) : filtered_subreads.%.fasta : filtered_regions.%.fofn input.%.fofn
	$(QSUB) pls2fasta -trimByRegion -regionTable $^ $@
##

## Filtering ##
filter : $(REGFOFNS) ;

$(REGFOFNS) : filtered_regions.%.fofn : input.%.fofn
	$(QSUB) filter_plsh5.py --filter='MinReadScore=0.80,MinSRL=500,MinRL=100' --trim='True' --outputFofn=$@ $<
##

## Initial chunking ##
$(BAXFOFNS) : input.fofn
	awk 'BEGIN{c=1}{print $$0 > sprintf("input.chunk%03dof%03d.fofn", c++, $(CHUNK_SIZE))}' $<
##

input.fofn : $(BAXFILES)
	echo $^ | sed 's/ /\n/g' > $@

$(BAXFILES) : ;

clean :
	rm -f *.fofn
	rm -f *.fast[a,q]
	rm -f *.lengths
	rm -f *.cutoff
	rm -f *.rgn.h5
	rm -f *.m4
	rm -f asm.*
	rm -rf celera-assembler.*
	rm -rf [0-4]-*
	rm -rf runCA-logs
	rm -rf reference*
	rm -f unitig.lst
	rm -f aligned*
	rm -f polished*
	rm -f corrections.gff
	rm -f *.sh

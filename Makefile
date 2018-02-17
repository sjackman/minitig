# Parameters of minitig
k=32
w=28
f=3

# Number of threads
t=16

# Compress a file
gzip=pigz -p$t

# Report run time and memory usage using zsh
export SHELL=zsh -e -o pipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J

all: lint test

.DELETE_ON_ERROR:
.SECONDARY:

################################################################################
# Test

# Check the source code for errors with Pylint.
lint:
	pylint --rcfile=.pylintrc minitig

# Test Minitig.
test: \
	mt.pe.bfc.minitig.fa \
	mt.pe.bfc.minitig.mt.sort.bam.bai \
	mt.pe.bfc.minitig.gfa.png \
	mt.pe.minitig.fa \
	mt.pe.minitig.mt.sort.bam.bai \
	mt.pe.minitig.gfa.png

# Download the human mitochondrial genome.
mt.fa:
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz \
	| gunzip -c | tr -d N | seqtk seq >$@

# Convert Postscript to PDF.
%.pdf: %.ps
	pstopdf -o $@ $<

################################################################################
# Samtools

# Simulate paired-end reads using wgsim.
%.pe.fq.gz: %.fa
	wgsim -e 0.002 -r 0 -d 400 -s 100 -1 150 -2 150 -S 99 -N 3000 $< $*.pe.1.fq $*.pe.2.fq
	seqtk mergepe $*.pe.1.fq $*.pe.2.fq | $(gzip) >$@
	rm -f $*.pe.1.fq $*.pe.2.fq

# Index a BAM file with samtools.
%.sort.bam.bai: %.sort.bam
	samtools index $<

################################################################################
# BFC

# Correct reads.
%.bfc.fq.gz: %.fq.gz
	bfc -t$t $< | $(gzip) >$@

################################################################################
# Minimap2

# Map sequences to the reference.
%.$(ref).sort.bam: $(ref).fa %.fa
	minimap2 -a -x asm10 $^ | samtools sort -o $@

# Map reads to the reference.
%.$(ref).sort.bam: $(ref).fa %.fq.gz
	minimap2 -a -x sr $^ | samtools sort -o $@

################################################################################
# Miniasm

# Create a dot plot.
%.minidot.ps: %.paf.gz
	minidot $< >$@

################################################################################
# Minitig

# Index a FASTA file.
%.minitig.json: %.fa
	./minitig index -k$k -w$w $< >$@

# Index a FASTQ file.
%.minitig.json: %.fq.gz
	gunzip -c $< | ./minitig index -k$k -w$w - >$@

# Create a graph of minimizers of a FASTA file.
%.minitig.graph.gv: %.fa
	./minitig graph -k$k -w$w $< >$@

# Create a graph of minimizers of a FASTQ file.
%.minitig.graph.gv: %.fq.gz
	gunzip -c $< | ./minitig graph -k$k -w$w - >$@

# Assemble a FASTQ file.
%.minitig.fa %.minitig.gfa: %.fq.gz
	gunzip -c $< | ./minitig assemble -k$k -w$w -f$f -g $*.minitig.gfa - >$*.minitig.fa

################################################################################
# Bandage

# Render an assembly graph to PNG.
%.gfa.png: %.gfa
	bandage image $< $@

# Render an assembly graph to SVG.
%.gfa.svg: %.gfa
	bandage image $< $@

################################################################################
# Graphviz

# Render a graph using dot.
%.gv.pdf: %.gv
	dot -Tpdf -o $@ $<

all: lint test

# Parameters of minitig
k=32
w=48

# Number of threads
t=16

# Compress a file
gzip=pigz -p$t

# Report run time and memory usage using zsh
export SHELL=zsh -e -o pipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J

.DELETE_ON_ERROR:
.SECONDARY:

################################################################################
# Test

# Check the source code for errors with Pylint.
lint:
	pylint --rcfile=.pylintrc minitig

# Test Minitig.
test: mt.pe.minitig.gv.pdf

# Download the human mitochondrial genome.
mt.fa:
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz \
	| gunzip -c | tr -d N | seqtk seq >$@

# Convert Postscript to PDF.
%.pdf: %.ps
	pstopdf -o $@ $<

# Simulate paired-end reads using wgsim.
%.pe.fq.gz: %.fa
	wgsim -e 0 -r 0 -d 400 -s 100 -1 150 -2 150 -S 99 -N 3000 $< $*.pe.1.fq $*.pe.2.fq
	seqtk mergepe $*.pe.1.fq $*.pe.2.fq | $(gzip) >$@
	rm -f $*.pe.1.fq $*.pe.2.fq

# Correct reads using BFC.
%.bfc.fq.gz: %.fq.gz
	bfc -t$t $< | $(gzip) >$@

# Map reads to the reference using minimap2.
%.$(ref).sort.bam: $(ref).fa %.fq.gz
	minimap2 -a -x sr $^ | samtools sort -o $@

# Index a BAM file with samtools.
%.sort.bam.bai: %.sort.bam
	samtools index $<

# Create a dot plot using minidot from miniasm.
%.minidot.ps: %.paf.gz
	minidot $< >$@

# Index a FASTA file with Minitig.
%.minitig.json: %.fa
	./minitig index -k$k -w$w $< >$@

# Index a FASTQ file with Minitig.
%.minitig.json: %.fq.gz
	gunzip -c $< | ./minitig index -k$k -w$w - >$@

# Map sequences with Minitig.
%.minitig.gv: %.fa
	./minitig map -k$k -w$w $< >$@

# Map sequences with Minitig.
%.minitig.gv: %.fq.gz
	gunzip -c $< | ./minitig map -k$k -w$w - >$@

# Convert Minitig JSON to GraphViz.
%.minitig.tsv: %.minitig.json
	jq . $< | sed 's/^  /{ "u": /;s/,$$/ },/;s/^]$$/} ]/' | mlr --ijson --otsvlite cat >$@

# Render a graph with GraphViz.
%.gv.pdf: %.gv
	dot -Tpdf -o $@ $<

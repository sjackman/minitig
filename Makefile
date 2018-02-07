all: lint test

# Parameters of minitig
w=64
k=32

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
test: mt.minitig

# Download the human mitochondrial genome.
mt.fa:
	curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz \
	| gunzip -c | seqtk seq >$@

# Simulate paired-end reads using wgsim.
%.pe.fq.gz: %.fa
	wgsim -r 0 -d 400 -s 100 -1 150 -2 150 -S 1 -N 3000 $< $*.pe.1.fq $*.pe.2.fq
	seqtk mergepe $*.pe.1.fq $*.pe.2.fq | $(gzip) >$@
	rm -f $*.pe.1.fq $*.pe.2.fq

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

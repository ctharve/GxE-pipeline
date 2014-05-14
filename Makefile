## 3.13.14
## rpique [2014-03-13  8:07 ] I changed the resource requirments for sorting. 
## Make .bams for system plate 4
genomeindex=/wsu/home/groups/piquelab/data/RefGenome/hg19.fa.gz
fastqs=../fastqs
sortedBams=$(shell find ./ -name '*sorted.bam' | sed 's/_[SL].*\.bam//g' | sort | uniq)
Qsub.ppn=16
Qsub.q=masxq
Qsub.N=bams


all: 
	echo $(sortedBams)

.SECONDARY:

%.Qsub:
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q $(Qsub.q) -l nodes=1:ppn=$(Qsub.ppn) -N $(Qsub.N) -o $@ -e $@.e
	sleep 1;

all_bams: $(patsubst $(fastqs)/%_R1.fastq.gz,%.bam.Qsub,$(wildcard $(fastqs)/*_R1.fastq.gz))

%.bam: $(fastqs)/%_R1.fastq.gz $(fastqs)/%_R2.fastq.gz
	bwa mem -t $(Qsub.ppn) $(genomeindex) $^  | samtools view -Sb1 - > $@
	samtools sort $@ $*_sorted
#	samtools sort -@$(Qsub.ppn) -m 1000000000 $@ $*_sorted
	samtools index $*_sorted.bam

all_clean: $(patsubst %, %_clean.bam.Qsub,$(sortedBams))

%_merged.bam: %_[SL]*_sorted.bam 
	samtools merge $@ $^
	samtools index $@

%_quality.bam: %_merged.bam		
	samtools view -b1 -q10 $^ > $@
	samtools index $@

%_clean.bam: %_quality.bam
	samtools view -b -q10 $^ | samtools rmdup - $@
	samtools index $@

clean:
	rm *merged*
	rm *quality*

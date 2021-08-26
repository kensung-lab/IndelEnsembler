# IndelEnsembler
IndelEnsembler is an ensemble method for identifying deletions (DELs), tandem duplications (DUPs) and insertions (INSs) (either novel or due to transposition) from next generation sequencing data. It merges calls from different callers: Lumpy, Manta, SurVIndel and TranSurVeyor.

### Preparing the reference
The reference fasta file should be indexed by both bwa and samtools. For example, assuming the file is hg19.fa, you should run

```
$ bwa index hg19.fa
$ samtools faidx hg19.fa
```
Although not mandatory, SurVIndel will generally give higher quality results if a simple repeats file is provided. This can normally be downloaded from the simpleRepeats table in UCSC. The header must be removed and only the chromosome, the start, the end and the period columns must be retained, i.e.:

```
$ cat downloaded-file | grep -v "#" | cut -f2,3,4,6 > file-for-survindel.bed
```
Alternatively, you can run [TRF](https://tandem.bu.edu/trf/trf.html) and use the provided trf-to-bed.sh, i.e.:

```
## Generated the repeats file
$ cd /path/to/SurVIndel
$ cat trf-output.dat | ./trf-to-bed.sh > /path/to/simple/repeats/file
```
### Preparing the BAM file
The BAM files should be coordinate sorted, indexed, and should contain the MC and MQ tags. MC and MQ tags can be added using Picard FixMateInformation (http://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation).

Supposing file.bam is the file resulting from the alignment:

```
$ java -jar picard.jar FixMateInformation I=file.bam
$ samtools sort file.bam > sorted.bam
$ samtools index sorted.bam
```

### Dependencies
1. Manta  
Manta's source code was cloned from the repository at [HERE](https://github.com/Illumina/manta/releases/download/v1.6.0/). We used manta-1.6.0-0. 

```
$ wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
$ tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2
$ cd manta
# Add Manta to your PATH
$ export PATH=/path/to/manta/bin:$PATH
# run Manta
$ python2 /path/to/manta/bin/configManta.py --bam /path/to/bamfile --referenceFasta /path/to/reference/fasta --runDir /an/empty/working/directory
$ /an/empty/working/directory/runWorkflow.py -m local -j 40
$ gunzip /an/empty/working/directory/results/variants/candidateSV.vcf.gz
```
2. Lumpy  
Lumpy's source code was cloned from the repository at [HERE](https://github.com/arq5x/lumpy-sv). We used LUMPY-0.2.13.

```
$ git clone --recursive https://github.com/arq5x/lumpy-sv.git
$ cd lumpy-sv
# Add Lumpy to your PATH
$ export PATH=/path/to/lumpy-sv/bin:$PATH
# run Lumpy
$ samtools view -b -F 1294 /path/to/bamfile > /path/to/discordants/bamfile
$ samtools view -h /path/to/bamfile | /path/to/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | /path/to/samtools view -Sb -  > /path/to/splitters/bamfile
$ samtools sort /path/to/discordants/bamfile -o /path/to/sorted/discordants/bamfile
$ samtools sort /path/to/splitters/bamfile -o /path/to/sorted/splitters/bamfile
$ /path/to/lumpy-sv/bin/lumpyexpress -B /path/to/bamfile -S /path/to/sorted/splitters/bamfile -D /path/to/sorted/discordants/bamfile -o /an/empty/working/directory
```
3. SurVIndel  
SurVIndel's source code was cloned from the repository at [HERE](https://github.com/Mesh89/SurVIndel).

```
$ git clone https://github.com/Mesh89/SurVIndel.git
$ cd SurVIndel
$ ./build_htslib.sh
$ cmake -DCMAKE_BUILD_TYPE=Release . && make
# run SurVIndel
$ python /path/to/SurVIndel/surveyor.py /path/to/bamfile /an/empty/working/directory	/path/to/reference/fasta --threads 40 --samtools /path/to/samtools --bwa /path/to/bwa --simple-rep /path/to/simple/repeats/file
$ ./filter /path/to/working/directory alpha-value score-cutoff min-size simple-repeats
# The value of alpha-value, score-cutoff and min-size, you can refer (https://github.com/Mesh89/SurVIndel)
```
4. TranSurVeyor  
TranSurVeyor's source code was cloned from the repository at [HERE](https://github.com/Mesh89/TranSurVeyor).

```
$ git clone https://github.com/Mesh89/TranSurVeyor.git
$ cd TranSurVeyor
$ ./build_htslib.sh
$ cmake -DCMAKE_BUILD_TYPE=Release . && make
# run TranSurVeyor
$ python surveyor.py /path/to/bamfile /an/empty/working/directory	/path/to/reference/fasta --threads 40 --samtools /path/to/samtools --bwa /path/to/bwa --maxTRAsize 10000
$ ./filter /path/to/working/directory
```


### Download and Usage
Installing IndelEnsembler is easy. You can download and uncompress the IndelEnsembler package using wget or through git.



```
# download the IndelEnsembler
$ wget https://github.com/kensung-lab/IndelEnsembler/archive/refs/heads/main.zip
or
$ git clone https://github.com/kensung-lab/IndelEnsembler.git
$ ./build.sh
# Usage
$ cd IndelEnsembler
$ vi pipeline.sh
# Change the path of ref_genome and repeats
$ bash pipeline.sh /path/to/manta/result /path/to/lumpy/result /path/to/survindel/result /path/to/transurveyor/result /an/empty/working/directory 200
```

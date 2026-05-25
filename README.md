# F64
## License

This software is dual-licensed:

Academic version:
- Open-source under the GNU GPL v3 for academic and non-commercial use

Industrial version:
- Commercial license required for proprietary, industrial, or SaaS usage

For commercial licensing, contact: Alar.Aints@gmail.com, Alar.Aints@ut.ee


F64 © Alar Aints 2024, 2025, 2026
a program for filtering sequence reads


Options: 

    -F 	--filterFile 		Filter binary file(s) ( none to many )
    -f 	--flatFile		flat filter file (optional) ( none or one )
    -O 	--outPath 		outPath where the filtered outfiles go (optional)
    -o 	--outFile		output file for build merge and refine
    -i 	--input 		input files 
    -p 	--paired		paired or single files [2]
    -L 	--readLimit		readlLimit, 0 reads the whole files [0]
    -m 	--minLength		min length of reads to keep [40]
    -s 	--strictness		strictness: that many hits will trigger skip [1]
    -d 	--disableReset		disable automatic reset of parameters -m and -s (Reset function is On by default)
    -D 	--density		scanning density [3]
    -K 	--K			length of kmer (for build) [25]
    -R 	--tableReport		report index table specifics to stderr
    -j 	--json			print json report of filter data to stdout
    -t  --time			print various time data to stderr
       	--printPositive		print positivley identified reads instead of unrecognized reads
       	--stdin1
       	--stdin2
	      --stdout1
	      --stdout2

Commands:	Options:

    stats	-F -f 
		reads filterfiles ( -F and -f ) (no import) and reports the K values and number of indices. 

    merge	-F -f -o -R -j
		imports binfiles ( -F ) and merges their data. All files must have the same K value. If -f name is provided, output is written in flat hash format. If the -f file exists, the program exits with message: "Flat Files can not be updated or merged". If -o is specified, and exists, a warning is printed : "Outfile is replaced, not updated !". If both -f and -o are provided, -f is written as flat file, and -o is written as binfile. By default, the data are saved in bin format as the first -F file. If the file exists, a new version is created. If -F name is provided, but not found (and -o is not provided), the filename is treated as outfile, else the program exits with message "binfile -F not found. Check the filenames." 

    build	-F -f -o -K -i -R -j 
		imports binfiles ( -F ) and merges their data. All files must have the same K value. If -f name is provided, output is written in flat hash format. If the -f file exists, the program exits with message: "Flat Files can not be updated or merged". If -o is specified, and exists, a warning is printed : "Outfile is replaced, not updated !". If both -f and -o are provided, -f is written as flat file, and -o is written as binfile. By default, the data are saved in bin format as the first -F file. If the file exists, a new version is created. If -F name is provided, but not found (and -o is not provided), the filenameis treated as outfile, else the program exits with message "binfile -F not found. Check the filenames." Then, imports -i files and adds their data. Duplicate data are skipped automatically. fastq and fasta format are automatically detected. All files are analyzed independently. If "R1" exists in fastq file name, the reads are reverse-transcribed before analysis.

    refine	-F -f -o -i -j
		imports binfiles ( -F ) and merges their data. All files must have the same K value. If -f name is provided, output is written in flat hash format. If the -f file exists, the program exits with message: "Flat Files can not be updated or merged". If -o is specified, and exists, a warning is printed : "Outfile is replaced, not updated !". If both -f and -o are provided, -f is written as flat file, and -o is written as binfile. By default, the data are saved in bin format as the first -F file. If the file exists, a new version is created. If -F name is provided, but not found (and -o is not provided), the filenameis treated as outfile, else the program exits with message "binfile -F not found. Check the filenames." Then, imports -i files and REMOVES their data from index table. fastq and fasta format are automatically detected. All files are analyzed independently. All reads are also reverse-transcribed before analysis. 

    filter 	-F -f -L -p -m -s -D -d -R -j -O --stdin1 --stdin2 --stdout1 --stdout2 --printPositive
		default behaviour (no command needed). Imports filter files ( -F or -f ). All files must have the same K value. Then, fastq files (-i) are analyzed incrementally. Readlimit (-L) can be used to analyze only a certain number of reads. By defaut, all reads are analyzed. Fastq files can be paired reads or single reads, this can be specified with parameter -p [1/2]. By default, paired reads are assumed. Reads shorter or equal to cutoff length -m are skipped (default 40). Strictness (-s) means how many hash matches triggers skipping (default 1). Setting -d disables automatic resetting of parameters -m and -s based on average read length of the first 100k reads. Output files' destination can be specified with -O. By default, output goes to the same folder as input. Output file names are automatically derived: ".Filtered" is added to the filename before the last "." . 

Usage examples:

    F64 build -K 25 -i org1.fasta STRAINS/*.fasta org1extra.fasta org1Mito.fasta org1_strainX_R1.fastq org1_strainX_R2.fastq -o Org1.25.bin 
    F64 build -F Org1.25.bin -i more.fasta more2.fasta   // -> Org1.25.bin.v1
    F64 build -F NegCont.25.bin -i NegCont_R1.fastq NegCont_R2.fastq -K 25
    F64 merge -F Org1.25.bin.v1 NegCont.25.bin -o Exp1.25.bin -f Exp1.25.flat
    F64 refine -F Exp1.25.bin -o Exp1.25.refined.bin -f Exp1.25.refined.flat -i bact1.fasta bact2.fasta metag_R1.fastq metag_R2.fastq
    F64 -f Exp1.25.refined.flat -p 1 -i test1.fastq  // -> test1.filtered.fastq
    F64 filter -f Exp1.25.refined.flat --stdin1 <(zcat R1.fq.gz) --stdin2 <(zcat R2.fq.gz) --stdout1 >(gzip > R1.filtered.fq.gz) --stdout2 >(gzip > R2.filtered.fq.gz) --json
    F64 -f Exp1.25.refined.flat -p 2 -i Exp1/*.fastq
    F64 -F Exp1.25.refined.bin SILVA138.25.bin -d -s 2 -m 55 -D 1 -i Exp1/*.fastq -O Result2/


System requirements

Disk space: 
the index files are stored as binary int collections. Full human index is 27 GB, SILVA is 2.9 GB. For repeated routine use, the index files can be saved in flat hash format. These load quicker, but can not be merged or updated. 

Memory:
up to 220GB for full human db in bin format. (120GB of RAM are required for binfile import and merging. For conversion to cache-efficient flat hash format, another 100 GB is required. The 120GB is then freed. For flat file 100 GB is needed). Minimally, more than 16 GB is required. 

Dependencies: 
F64 has no dependencies. F64 is written in plain C. It needs to be compiled using gcc -std=gnu99 -pthread options. Low complexity regions are handled like everything else (no masking). 

File formats:
Filter files can be stored in two binary formats: bin and flat. Bin files can be merged and updated incrementally and converted to flat format. Flat files can not be updated or merged. Bin files must be entered as -F arguments (none to many), flat file must be entered as -f argument (none or one). Either -F or -f must be given. The formats are verified internally.
For build, fasta and fastq files can be used, can be mixed. File type is detected automatically.
For filtration, fasta and fastq files can be entered. Paired and single files can not be mixed. 

MacOS notice:
F64 runs on MacOS and Linux. 

Network connectivity:
F64 has no network connectivity functions. 

Multihreading: 
Optimally, F64 runs on 6-8 threads using mpmc architecture. 

Behaviour
The program first loads the filter file(s) containing 64-bit hash values of the sequences to be removed. All these files must have been prepared using the same kmer length. The files are merged in memory. For large files containing billions of indices 120GB of RAM are required. For conversion to cache-efficient flat hash format, another 100 GB is required. The 120GB is then freed. Fastq files are analyzed incrementally. If the input folder contains no other files, and they are named sequentially, -i FOLDER/*.* can be used. Paired reads are analyzed in the order read2(sense) -> read2(antisense) -> read1(sense)-> read1(antisense). When skip is triggered, the cycle is broken. Either read can trigger the skipping of the pair. "_Filtered" is added to the file names, before the last "." character. If the output file already exists, the process moves to the next pair (or single file). Paired files are assumed. Reads containing near-terminal "N" are truncated in the beginnig or end respectively, if the remaining sequence is not shorter than cutoff (-m). If the remaining read contains more than 1/50 N, it is skipped. Reads not recognized by the filter are printed to output files. Progress report and error messages are printed to stderr. Json reports are printed to stdout.

Run optimization

F64 is designed for fastq read filtration, for example removing human, or ribosomal, (or both) reads from metatranscriptome samples. Strictness 1 allows efficient removal of even slightly mutated target reads, however, some non-targeted reads are also removed. (less than 200reads / 10M,  in human - mouse mixture, k=25). Increasing kmer length or strictness reduces false positive rate, but increases sensitivity to mutations. When many files with different length reads are processed, it is advisable to allow ( no -d ) automatic read length detection, and strictness and cutoff optimization. Cutoff is set to 70% of the average read length of the first 100k reads, (min 40 max 70); strictness is calculated as s = (int)Avg.length / k / 2 (min 1). Processing files in batches is preferrable, because filterfile loading and conversion takes about 5 minutes (bin files) or 15 sec (flat files), while processing itself proceeds at 50 M reads/min (mixed orientation, 10% retainable reads) including file loading and N detection, or up to 70 M reads/min (sense orientation, 0.5% retainable reads). 
Many index binfiles can be loaded for one run, they are merged. (For example Human genome + SILVA + experiment-specific negative controls).
For repeated routine use, the index files can be merged and saved in flat hash format. These load quicker, but can not be further updated or merged. 

Data management flexibility:
Unlike existing tools that require monolithic database construction, our method supports incremental index file updates and merging, enabling flexible and efficient integration of diverse reference datasets.

Modularity
We support modular, incrementally extensible reference index libraries. In practice, this means that new data can be added to existing index files without having to recalculate the whole library. Small, custom libraries can be created with ease, containing experiment-specific, or lab-specific negative control data, from fasta, or directly from fastq files. These small libraries can be merged seamlessly with larger libraries on the go, to provide customizable flexibility for each experiment. Updates create new versions of index library files, allowing strict validation and versioning.

Scalability
Our indexing algorithm supports true 64-bit index storage in 32-bit table, delivering robust discrimination power required for metagenomic analysis in a memory-optimized, yet scalable format. Merging large libraries containing billions of indices is not limited to 4.2 B (UINT_MAX), while keeping RAM requirement in the 120 GB range. Essentially, the filter is only limited by hardware.

Batch processing
The program supports batch input, making processing entire folder contents hassle-free. The output can be directed to a user-specified destination, or left in the source folder. 

Auto-calibration
By default, the program sets read length cutoff and strictness parameters based on the average length of the first 100k reads. This is particularly useful when analyzing large number of files with different read lengths. Auto-calibration can be overriden by the user, and is turned off when using --stdin.


Performance metrics:

Speed
Processing speed is 50 - 60 M reads/min, including file loading, reverse transcription, detection, N removal and output file writing. Index loading and conversion proceeds at 7 GB/min (bin files) or 2 GB/sec (flat files). For filtering  simulated datasets containing 10M human genomic reads and 1M E.coli reads, speed gains of 16x and 43x were observed, compared to Kraken2 and BWA-mem, respectively, not accounting for hands-on time. 

Accuracy
For unmutated test reads, the program displays Sensitivity = 1, which means every single targeted read is removed. Mutation burden prevents complete removal. If a read can not be identified, it is left untouched. Specificity, or retention of non-targeted reads, depends on the kmer length and strictness parameters. For distantly related organisms, such as H.sapiens - E.coli mixture, sensitivity and specificity are 1. (Kraken2: Sensitivity = 0.984011, Specificity = 0.999977; BWA-mem : Sensitivity = 0.9935316, Specificity = 0.99996) For Human -  Mouse genomic mixture, sensitivity is between 0.997-0.999, spcificity is between 0.99998 and 1 .  Increasing k and s will improve specificity, but increase retention of mutated targeted reads. The high accuracy is achieved by our unique library refinement function.



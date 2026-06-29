#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
// #include <sys/sysinfo.h>
// #include <sys/sysctl.h>
#include <getopt.h>
#include <ctype.h>
#include <pthread.h>
#include <sched.h>
#include <stdatomic.h>

    int FilesCount=0; 
    char** fnameHS = NULL;  // -F Filter binary file(s)               
    char* outfile1 = NULL;  //    Automatically derived from input files
    char* outfile2 = NULL;  //    Automatically derived from input files
    char* outfile = NULL;   // -o outfile for build, merge and refine
    char* OutPath = NULL;   // -O outfiles destination (Optional. By default the outfiles go to the same folder as input files.)
    int readlimit = 0;      // -L readlimit, 0 reads the whole files
    int cutoff = 40;        // -m min length of reads to keep
    int strictness=1;       // -s that many hits will trigger skip
    int DD = 3;             // -D density of scanning. typically -D 5 is enough
    int disableReset=0;     // -d disable automatic reset of parameters
    char** fnamesInput=NULL;// -i input files 
    int pairs = 2;          // -p single or paired files
    int FR = 0;             // -K length of kmer
    int TableReport = 0;    // -R report index table specifics
    char* FlatFile = NULL;  // -f read/save Filter file in flat hash format
    int json = 0;           // -j print json report of filter data
    int takeTime = 0;       // -t measure time
    int detectN = 0;        // -N detect Ns
    int InputCount=0;   
    int filter = 0;
    int stats = 0;
    int build = 0;
    int merge = 0;
    int printPositive = 0;
    int refine = 0;
    int help = 0;
    int use_stdin = 0;
    char* stdin1_path = NULL;
    char* stdin2_path = NULL;
    char* stdout1_path = NULL;
    char* stdout2_path = NULL;

    int N_WORKERS = 1;       // determined automatically
    int N_BUILDERS = 1;     //  determined automatically


    // int version = 2;        // added input file buffering
    // int version = 3;        // added FlatHash lookup; 
    // int version = 4;        // handling of single or paired read files, removed Forward reads option. All reads are analyzed in both orientations.
    //int version = 5;         // change of order: Scanning, N detection, scan again
    //int version = 6;         // CLI API
    //int version = 7;         // streaming and multithreading
    //int version = 8;         // multithreading and worker batch
    //int version = 10;        // write buffering
    //int version = 11;        // Added back stats, build, OutPath, readlimit, -m and -s reset, -d, -p and folder/*.* input
    // int version = 12.2;     // Added multithreading to build
    //int version = 13;        // Added flat hash format saving option
    // int version = 14;       // Added positive selection option
    //int version = 15;        // Suport for fasta file analysis
    int version = 16;          // added refine command
    
    char* msg = "\n© Alar Aints 2024-2026\n\n Usage:\n\n\t -F <FilterFile(s).bin> <F2.bin> <etc.bin> \n\t -f <flatFile.flat>\n\t -o <Filter.out.bin>\n\t -d disableReset\n\t -m minLength [int] \n\t -s strictness [int] \n\t -p single[1]/paired[2]\n\t -K kmer length [int] \n\t -R report specifics\n\t -j json report\n\t -v version\n\t -h help\n\t -L readLimit [int] \n\t -D scanDensity [int] \n\t -O OutPath/\n\t --printPositive\n\t --stdin1 <path>\n\t --stdin2 <path>\n\t --stdout1 <path>\n\t --stdout2 <path>\n\t -i <File1.fastq> <File2.fastq> <etc...>\n\n";
    char* msgHelp = "type \"F64 help\" for usage examples.\n";
    char* msg2 = "\n\
    // CLI Usage Examples:\n\
    // \n\
    // F64 build -K 25 -i org1.fasta STRAINS/*.fasta org1extra.fasta org1Mito.fasta org1_strainX_R1.fastq org1_strainX_R2.fastq -o Org1.25.bin \n\
    // F64 build -F Org1.25.bin -i more.fasta more2.fasta   // -> Org1.25.bin.v1\n\
    // F64 build -F NegCont.25.bin -i NegCont_R1.fastq NegCont_R2.fastq -K 25\n\
    // F64 merge -F Org1.25.bin.v1 NegCont.25.bin -o Exp1.25.bin -f Exp1.25.flat\n\
    // F64 refine -F Exp1.25.bin -o Exp1.25.refined.bin -f Exp1.25.refined.flat -i bact1.fasta bact2.fasta metag_R1.fastq metag_R2.fastq\n\
    // F64 -f Exp1.25.refined.flat -p 1 -i test1.fastq  // -> test1.filtered.fastq\n\
    // F64 filter -f Exp1.25.refined.flat --stdin1 <(zcat R1.fq.gz) --stdin2 <(zcat R2.fq.gz) --stdout1 >(gzip > R1.filtered.fq.gz) --stdout2 >(gzip > R2.filtered.fq.gz) --json\n\
    // F64 -f Exp1.25.refined.flat -p 2 -d -i Exp1/*.fastq\n\
    // F64 -F Exp1.25.refined.bin SILVA138.25.bin -i Exp1/*.fastq -O Result2/\n\
    \n";


  /*
 * F64 - High-performance read filtering engine
 *
 * Copyright (C) 2026 Alar Aints / Tartu University
 * - Open-source under the GNU GPL v3 for academic and non-commercial use
 * - Commercial license required for proprietary, industrial, or SaaS usage

 * For commercial licensing, contact: Alar.Aints@gmail.com, Alar.Aints@ut.ee

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
 static struct option long_opts[] = {
        {"filter", no_argument, &filter, 1},
        {"stats", no_argument, &stats, 1},
        {"build", no_argument, &build, 1},
        {"merge", no_argument, &merge, 1},
        {"printPositive", no_argument, &printPositive, 1},
        {"refine", no_argument, &refine, 1},
        {"input", required_argument, 0, 'i'},
        {"filterFile", required_argument, 0, 'F'},
        {"outPath", required_argument, 0, 'O'},
        {"outFile", required_argument, 0, 'o'},
        {"paired", required_argument, 0, 'p'},
        {"readLimit", required_argument, 0, 'L'},
        {"minLength", required_argument, 0, 'm'},
        {"strictness", required_argument, 0, 's'},
        {"density", required_argument, 0, 'D'},
        {"K", required_argument, 0, 'K'},
        //{"Workers", required_argument, 0, 'T'},
        //{"Builders", required_argument, 0, 'B'},
        {"tableReport", no_argument, 0, 'R'},
        {"detectNs", no_argument, 0, 'N'},
        {"disableReset", no_argument, 0, 'd'},
        {"flatFile", required_argument, 0, 'f'},
        //{"index", no_argument, 0, 'x'},
        {"json", no_argument, 0, 'j'},
        {"time", no_argument, 0, 't'},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 1},
        {"stdin1", required_argument, 0, 1001},
        {"stdin2", required_argument, 0, 1002},
        {"stdout1", required_argument, 0, 1003},
        {"stdout2", required_argument, 0, 1004},
        {0,0,0,0}
    };

void parse_args(int argc, char** argv) {
    int opt;
    int long_index;

    while ((opt = getopt_long(argc, argv, "hF:L:D:m:s:o:O:tdjvf:RNi:p:K:", long_opts, &long_index)) != -1) {
        switch (opt) {
        case 0:
            // long option with flag was processed
            break;
        case 'h':
            printf("F64.%d\n%s\n%s\n",version,msg,msgHelp);
            exit(EXIT_SUCCESS);  
        case 'F': {
            fnameHS = calloc(argc, sizeof(char*));
            fnameHS[0]= optarg;
            FilesCount++;
            // Collect following arguments that are not options
            while (optind < argc && argv[optind][0] != '-') {
                fnameHS[FilesCount] = argv[optind];
                FilesCount++;
                optind++;
            }
            break;
        }
        case 'i': {
            fnamesInput = calloc(argc, sizeof(char*));
            fnamesInput[0]= optarg;
            InputCount++;
            // Collect following arguments that are not options
            while (optind < argc && argv[optind][0] != '-') {
                fnamesInput[InputCount] = argv[optind];
                InputCount++;
                optind++;
            }
            break;
        }
        case 'L':
            readlimit = atoi(optarg);
            break;  
        case 'D':
            DD = atoi(optarg);
            break;   
        case 'm':
            cutoff = atoi(optarg);
            break;
        case 'p':
            pairs = atoi(optarg);
            break;  
        case 's':
            strictness = atoi(optarg);
            break;
        case 'K':
            FR = atoi(optarg);
            break;
        case 'R':
            TableReport = 1;
            break;
        case 'N':
            detectN = 1;
            break;    
        case 'd':
            disableReset = 1;
            break;
        case 't':
            takeTime = 1;
            break;
        case 'f':
            FlatFile = optarg;
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'O':
            OutPath = optarg;
            break;
        case 'j':
            json = 1;
            break;
        case 'v':
            printf("F64 version %d",version);;
            break;
        case 1001:
            use_stdin = 1;
            stdin1_path = optarg;
            break;
        case 1002:
            use_stdin = 1;
            stdin2_path = optarg;
            break;
        case 1003:
            stdout1_path = optarg;
            break;
        case 1004:
            stdout2_path = optarg;
            break;
        default:
            fprintf(stderr, "F64.%d\n%s\n",version,msg);
            exit(EXIT_FAILURE);
        }
    }
}

// MurmurHash2-64A by Austin Appleby
#if defined(_MSC_VER)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else   // defined(_MSC_VER)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

uint64_t MurmurHash64A ( const void * key, int len, uint64_t seed )
{
  const uint64_t m = BIG_CONSTANT(0xc6a4a7935bd1e995);
  const int r = 47;

  uint64_t h = seed ^ (len * m);

  const uint64_t * data = (const uint64_t *)key;
  const uint64_t * end = data + (len/8);

  while(data != end)
  {
    uint64_t k = *data++;

    k *= m; 
    k ^= k >> r; 
    k *= m; 
    
    h ^= k;
    h *= m; 
  }

  const unsigned char * data2 = (const unsigned char*)data;

  switch(len & 7)
  {
  case 7: h ^= (uint64_t)(data2[6]) << 48;
  case 6: h ^= (uint64_t)(data2[5]) << 40;
  case 5: h ^= (uint64_t)(data2[4]) << 32;
  case 4: h ^= (uint64_t)(data2[3]) << 24;
  case 3: h ^= (uint64_t)(data2[2]) << 16;
  case 2: h ^= (uint64_t)(data2[1]) << 8;
  case 1: h ^= (uint64_t)(data2[0]);
          h *= m;
  };
 
  h ^= h >> r;
  h *= m;
  h ^= h >> r;

  return h;
}
// end of MurmurHash2-64A 

const char PT[] = "ACGTN";
char* revTrans(char* seq) {
    int seq_len = strlen(seq);
    char* rt = malloc((seq_len + 1) * sizeof(char));
    int i = seq_len - 1;
    for (; i >= 0; i--) {
        char r = seq[i];
        int j = 0;
        for (; j < 5; j++) {
            if (PT[j] == r) {
                break;
            }
        }
        rt[seq_len - 1 - i] = (j >= 4) ? r : PT[3 - j]; // j==5: base not in ACGTN (N/IUPAC/other) -> keep unchanged (was OOB PT[3-5])
    }
    rt[seq_len] = '\0';
    return rt;
}

void printAeg(int aeg) { // igaks juhuks...
    int sekund = aeg % 60;
    int jk = aeg / 60;
    int minut = jk % 60;
    jk = jk / 60;
    int tund = jk % 24;
    int pv = jk / 24;
    fprintf(stderr,"%d_%02d:%02d:%02d\n", pv, tund, minut, sekund);
}

uint32_t htSize = 0xffffffff; // UINT_MAX = 4 294 967 295 = 2^32 - 1
uint64_t osSize = 4294967296; //(uint64_t)htSize + 1;
_Atomic size_t counter=0; // table allocations
_Atomic size_t entryCount=0; // all indices
uint64_t n_values = 0;
uint32_t* offsets = NULL;
uint32_t* values = NULL;
uint32_t** htHS = NULL;

#define HT_LOCK_COUNT 4096
static pthread_mutex_t ht_locks[HT_LOCK_COUNT];

static inline pthread_mutex_t *ht_lock_for(uint32_t index) {
    return &ht_locks[index & (HT_LOCK_COUNT - 1)];
}

void init_ht_locks(void) {
    for (int i = 0; i < HT_LOCK_COUNT; i++) {
        pthread_mutex_init(&ht_locks[i], NULL);
    }
}

static inline int bucket_insert_unique_locked(uint32_t index, uint32_t hsValue, uint32_t **htHS) {
    if (htHS[index] == NULL) {
        htHS[index] = malloc(2 * sizeof(uint32_t));
        if (htHS[index] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        htHS[index][0] = 2;
        htHS[index][1] = hsValue;
        return 2;
    }

    uint32_t isize = htHS[index][0];
    for (uint32_t ii = 1; ii < isize; ii++) {
        if (hsValue == htHS[index][ii]) {
            return 0;
        }
    }

    uint32_t *tmp = realloc(htHS[index], (isize + 1) * sizeof(uint32_t));
    if (tmp == NULL) {
        fprintf(stderr, "Realloc failed\n");
        exit(EXIT_FAILURE);
    }
    htHS[index] = tmp;
    htHS[index][0] = isize + 1;
    htHS[index][isize] = hsValue;
    return 1;
}

static inline int bucket_remove_value_locked(uint32_t index, uint32_t hsValue, uint32_t **htHS) {
    if (htHS[index] == NULL) {
        return 0;
    }

    uint32_t num = htHS[index][0];
    for (uint32_t ii = 1; ii < num; ii++) {
        if (htHS[index][ii] != hsValue) {
            continue;
        }

        if (num == 2) {
            free(htHS[index]);
            htHS[index] = NULL;
            return 2;
        }

        uint32_t mov = num - ii - 1;
        memmove(&htHS[index][ii], &htHS[index][ii + 1], mov * sizeof(uint32_t));
        htHS[index][0]--;
        return 1;
    }

    return 0;
}

void insertHash(int FR, size_t klen, char* target, uint32_t** htHS) {
    const size_t kmer_len = (size_t)FR;
    if (FR <= 0 || klen < kmer_len) {
        return;
    }

    size_t new_bucket_count = 0;
    size_t new_entry_count = 0;
    size_t segment_start = 0;

    while (segment_start < klen) {
        // Hash contiguous non-N segments directly instead of rescanning each
        // candidate k-mer window for ambiguous bases.
        while (segment_start < klen &&
               (target[segment_start] == 'N' || target[segment_start] == 'n')) {
            segment_start++;
        }

        size_t segment_end = segment_start;
        while (segment_end < klen &&
               target[segment_end] != 'N' && target[segment_end] != 'n') {
            segment_end++;
        }

        if (segment_end - segment_start >= kmer_len) {
            size_t limit = segment_end - kmer_len + 1;
            for (size_t i = segment_start; i < limit; i++) {
                uint64_t frN = MurmurHash64A(target + i, FR, 0);
                uint32_t index = (uint32_t)(frN >> 32);
                uint32_t hsValue = (uint32_t)(frN & 0xFFFFFFFF);

                pthread_mutex_t *lock = ht_lock_for(index);
                pthread_mutex_lock(lock);
                int insert_state = bucket_insert_unique_locked(index, hsValue, htHS);
                pthread_mutex_unlock(lock);

                if (insert_state != 0) {
                    new_entry_count++;
                    if (insert_state == 2) {
                        new_bucket_count++;
                    }
                }
            }
        }

        if (segment_end == klen) {
            break;
        }
        segment_start = segment_end + 1;
    }

    if (new_bucket_count) {
        atomic_fetch_add(&counter, new_bucket_count);
    }
    if (new_entry_count) {
        atomic_fetch_add(&entryCount, new_entry_count);
    }
}

    void write_outfile (char* outfile, int FR, uint32_t counter, uint32_t** htHS){
        FILE *Ffile = fopen(outfile, "wb");
        if (Ffile == NULL) {
            fprintf(stderr, "Error opening file for writing.\n");
            exit(EXIT_FAILURE);
        }
        
        fwrite("F64B", sizeof(uint32_t), 1, Ffile);                 // magic int
        fwrite(&(FR), sizeof(uint32_t), 1, Ffile);                  // K-mer length
        fwrite(&(counter), sizeof(uint32_t), 1, Ffile);             // number of table entries
        for (uint32_t i = 0; i < htSize; i++) {
            if(htHS[i] != NULL){
                uint32_t Entrysize = htHS[i][0];
                fwrite(&htHS[i][0],sizeof(uint32_t),1,Ffile);             // the number of ints to follow; 'Entrysize' in reader
                fwrite(&i,sizeof(uint32_t),1,Ffile);                      // the ht index
                fwrite(&htHS[i][1],sizeof(uint32_t),Entrysize-1,Ffile);   // the entries
            }
        }
        fclose(Ffile);
        fprintf(stderr, "outfile written: %s\n",outfile);           // FYI: Human Genome: fileSize 23 GB; memory size 98 GB
        return;
    }

void writeFlatOutfile (char* outfile){
    FILE *Ffile = fopen(outfile, "wb");
    if (Ffile == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
        exit(EXIT_FAILURE);
    }
    fwrite("F64F", sizeof(uint32_t), 1, Ffile);                      // magic int
    fwrite(&(FR), sizeof(uint32_t), 1, Ffile);                        // K-mer value
    fwrite(&osSize, sizeof(uint64_t), 1, Ffile);                       // size of *offsets
    fwrite(offsets,sizeof(uint32_t),osSize,Ffile);                   // the offsets
    fwrite(&n_values, sizeof(uint64_t), 1, Ffile);                     // the number of values
    fwrite(values,sizeof(uint32_t),n_values,Ffile);                   // the values
    fclose(Ffile);
    fprintf(stderr, "flat outfile written: %s\n",outfile);          
    return;
}
void readFlatfile (char* file){
    FILE *Ffile = fopen(file, "rb");
    if (Ffile == NULL) {
        fprintf(stderr, "Error opening flat file.\n");
        exit(EXIT_FAILURE);
    }
    uint64_t fTest = 0; // 4294967296;
    char HeaderF [4];
    fread(&HeaderF, sizeof(uint32_t), 1, Ffile);                       // F64Flat header
    if(memcmp(HeaderF, "F64F", 4) != 0){
        fprintf(stderr, "Wrong header ! %s is not F64 file.\n",file);
        exit(EXIT_FAILURE);
    }
    fread(&FR, sizeof(uint32_t), 1, Ffile);                          // K-mer value
    if(FR > 5000){
        fprintf(stderr, "Wrong K value ! %s is not F64 file.\n",file);
        exit(EXIT_FAILURE);
    }
    fread(&fTest,sizeof(uint64_t), 1, Ffile);                      // test file type
    if(fTest != osSize){ 
        fprintf(stderr, "Wrong osSize ! %s is not F64 flat file.\n",file);
        exit(EXIT_FAILURE);
    }
    offsets = calloc( osSize , sizeof(uint32_t));                  // malloc *offsets
    size_t nread = fread(offsets, sizeof(uint32_t), osSize, Ffile);     // read *offsets
    if (nread != osSize) {
        if (feof(Ffile)) {
            fprintf(stderr, "Unexpected end of file\n");
        } else if (ferror(Ffile)) {
            perror("fread failed");
        }
    }
    fread(&n_values, sizeof(uint64_t), 1, Ffile);                     // the number of values
    values = calloc( (n_values+1), sizeof(uint32_t));                // malloc *values
    fread(values,sizeof(uint32_t),n_values,Ffile);                   // the values
    fclose(Ffile);
    //fprintf(stderr, "outfile written: %s\n",file);          
    return;
}

    void convert_to_FlatHash(){
        time_t now = time(NULL);
        offsets = malloc( (osSize + 1) * sizeof(uint32_t)); // +1: offsets[osSize] is written below
        //if(offsets) {fprintf(stderr,"Offsets allocated: %llu \n",osSize); }
        //else {fprintf(stderr,"Offsets allocation error: \n"); exit(1);}
        // First pass: compute total size
        n_values = 0;
        for(uint32_t i = 0; i < htSize; i++) {
            offsets[i] = n_values;
            if(htHS[i]) {
                n_values += (htHS[i][0] - 1); // skip size field
            }
        }
        offsets[osSize] = n_values;
        values = malloc((n_values+1) * sizeof(uint32_t));
        //fprintf(stderr,"Values allocated: %llu \n",n_values);
        // Second pass: fill values
        size_t pos = 0;
        for(uint32_t i = 0; i < htSize; i++) {
            if( ! htHS[i] ) {continue;}
            uint32_t size = htHS[i][0];
            memcpy(&values[pos], &htHS[i][1], (size - 1) * sizeof(uint32_t));
            pos += (size - 1);    
        }
        time_t fin = time(NULL);
        if(takeTime){
            int Aeg = difftime(fin,now);
            fprintf(stderr,"FlatHash converison time: ");
            printAeg(Aeg);
            printf("\n");
        }
        return;
    }

static void add_bucket_report_count(uint64_t **counts_ptr, size_t *counts_size, uint32_t bucket_len) {
    if (bucket_len >= *counts_size) {
        size_t new_size = *counts_size ? *counts_size : 16;
        while (bucket_len >= new_size) {
            new_size *= 2;
        }
        uint64_t *tmp = realloc(*counts_ptr, new_size * sizeof(uint64_t));
        if (tmp == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        memset(tmp + *counts_size, 0, (new_size - *counts_size) * sizeof(uint64_t));
        *counts_ptr = tmp;
        *counts_size = new_size;
    }
    (*counts_ptr)[bucket_len]++;
}

static void print_bucket_report(uint64_t *counts, size_t counts_size) {
    fprintf(stderr,"\nhtHS report: %d-mer\n\n",FR);
    fprintf(stderr,"Length\tOccurrence\n");
    for(size_t i = 0; i < counts_size; i++){
        if(counts[i] > 0 ){
            fprintf(stderr,"%zu\t%12llu\n",i,(unsigned long long)counts[i]);
        }
    }
}

void report_hashtable (uint32_t*** htHS_ptr){ 
    if (htHS_ptr == NULL || *htHS_ptr == NULL) {
        fprintf(stderr, "\nhtHS report unavailable: table is not loaded\n");
        return;
    }
    uint32_t** htHS = *htHS_ptr;
    size_t counterTsize = atomic_load(&counter) + 1;
    if (counterTsize < 16) {
        counterTsize = 16;
    }
    uint64_t* counterT = calloc(counterTsize,sizeof(uint64_t));
    if (counterT == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for(uint32_t i = 0;i<htSize;i++){
        if(htHS[i] == NULL){
            add_bucket_report_count(&counterT, &counterTsize, 0);
            continue;
        }
        uint32_t es = htHS[i][0]-1;
        add_bucket_report_count(&counterT, &counterTsize, es);
    }
    print_bucket_report(counterT, counterTsize);
    free(counterT);
    counterT = NULL;
    return;
}
uint64_t filestatRes[3];

uint64_t* fileStat(char* file) {
    
    uint64_t count = 0;
    //char c;
    int cc = 0; // character count per line
    //int ll = 0; // line length
    //int maxll = 0;
    uint64_t fl = 0; // fasta line count
    uint64_t tc = 0; // total characters
    size_t line_length = 0;
    ssize_t read;
    char* line = NULL;
    // Open the file
    FILE *fp = fopen(file, "r");

    // Check if file is open
    if (fp == NULL) {
        fprintf(stderr,"Could not open file %s\n", file);
        return 0;
        exit(1);
    }

    // Count the number of lines
    while ((read = getline(&line, &line_length, fp)) != -1 ) {
        if (line[0] == '>') {
            fl++;
            count++;
            } // fasta count
        else{
            cc = strlen(line)-1;
            tc+=cc;
            count++;
        }
    }
    filestatRes[0]=tc;          // total characters
    filestatRes[1]=count;       // total lines
    filestatRes[2]=fl;          // fasta lines
    // Close the file
    fclose(fp);
    // fprintf(stderr,"Number of lines: %d", count);
    // fprintf(stderr,"Longest line: %d bp\n",maxll);
    return filestatRes;
}

char* openFasta ( char* fileIn){

    FILE* file = fopen(fileIn, "r");
    if (file == NULL) {
        fprintf(stderr,"Error opening file\n");
        return 0;
    }
         // --- 1. Find file size ---
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    rewind(file);  // Go back to start

        // --- 2. Allocate buffer ---
    char *buffer = malloc(size + 3);
    if (!buffer) {
        perror("Memory allocation failed");
        fclose(file);
        return 0;
    }
         // --- 3. Read entire file into memory ---
    size_t line_length;
    ssize_t read;
    char* line = NULL;
    size_t pos = 0;
    while ((read = getline(&line, &line_length, file)) != -1) {
        if (line[0] == '>'){buffer[pos++] ='\n';
            for (int i = 0; line[i] != '\0'; i++) {
                buffer[pos++] = line[i];
            }
        }
        else { 
            for (int i = 0; line[i] != '\0'; i++) {
                if (isalpha(line[i])){
                    if(line[i] >= 'a' && line[i] <= 'z'){line[i] -= 32;} // ensure upper case
                    buffer[pos++] = line[i];
                }
            }
        }        
    }
    fclose(file);
    free(line);
    buffer[pos++] = '\n'; 
    buffer[pos++] = '\0';  // Null-terminate to be safe
    return buffer;
}

char* build_stream_output_path(const char *fallback_name) {
    if (OutPath == NULL || OutPath[0] == '\0') {
        return strdup(fallback_name);
    }

    size_t dir_len = strlen(OutPath);
    int needs_slash = (OutPath[dir_len - 1] != '/');
    size_t total = dir_len + needs_slash + strlen(fallback_name) + 1;
    char *path = malloc(total);
    if (path == NULL) {
        fprintf(stderr, "Memory allocation failed for output path\n");
        exit(EXIT_FAILURE);
    }

    strcpy(path, OutPath);
    if (needs_slash) {
        strcat(path, "/");
    }
    strcat(path, fallback_name);
    return path;
}

char* derive_output_path(const char *input_name, const char *fallback_name) {
    if (input_name != NULL &&
        strcmp(input_name, "-") != 0 &&
        strncmp(input_name, "/dev/fd/", 8) != 0) {
        size_t in_len = strlen(input_name);
        char *out_name = malloc(in_len + 2);
        if (out_name == NULL) {
            fprintf(stderr, "Memory allocation failed for output filename\n");
            exit(EXIT_FAILURE);
        }
        strcpy(out_name, input_name);
        strcat(out_name, "F");
        return out_name;
    }

    return build_stream_output_path(fallback_name);
}
int dir_exists_and_writable(const char *path) {
    struct stat st;

    // Check if it exists
    if (stat(path, &st) != 0) {
        return 0;  // does not exist
    }

    // Check if it's a directory
    if (!S_ISDIR(st.st_mode)) {
        return 0;  // not a directory
    }

    // Check if writable
    if (access(path, W_OK) != 0) {
        return 0;  // not writable
    }

    return 1;  // exists, is directory, and writable
}
FILE* open_input_stream(const char *path, const char *label) {
    FILE *fp = NULL;

    if (path == NULL) {
        fprintf(stderr, "Missing input for %s\n", label);
        exit(EXIT_FAILURE);
    }

    if (strcmp(path, "-") == 0) {
        fp = stdin;
    } else {
        fp = fopen(path, "r");
    }

    if (fp == NULL) {
        fprintf(stderr, "Error opening %s: %s\n", label, path);
        exit(EXIT_FAILURE);
    }

    return fp;
}

FILE* open_output_stream(const char *path, const char *label) {
    FILE *fp = NULL;

    if (path == NULL) {
        fprintf(stderr, "Missing output for %s\n", label);
        exit(EXIT_FAILURE);
    }

    if (strcmp(path, "-") == 0) {
        fp = stdout;
    } else {
        fp = fopen(path, "w");
    }

    if (fp == NULL) {
        fprintf(stderr, "Error opening %s: %s\n", label, path);
        exit(EXIT_FAILURE);
    }

    return fp;
}

char* derive_filtered_output_path(const char *input_name, const char *fallback_name) {
    const char *suffix = printPositive ? ".positive.fastq" : ".filtered.fastq";
    size_t suffix_len = strlen(suffix);

    if (input_name == NULL ||
        strcmp(input_name, "-") == 0 ||
        strncmp(input_name, "/dev/fd/", 8) == 0) {
        return build_stream_output_path(fallback_name);
    }

    const char *basename = strrchr(input_name, '/');
    basename = (basename != NULL) ? basename + 1 : input_name;

    const char *dir = NULL;
    size_t dir_len = 0;
    if (OutPath != NULL && OutPath[0] != '\0') {
        dir = OutPath;
        dir_len = strlen(OutPath);
        if (! dir_exists_and_writable(OutPath)) {
            mkdir(OutPath, 0777);
    }
    } else if (basename != input_name) {
        dir = input_name;
        dir_len = (size_t)(basename - input_name);
    }

    const char *dot = strrchr(basename, '.');
    size_t stem_len = (dot != NULL) ? (size_t)(dot - basename) : strlen(basename);
    int needs_slash = (dir_len > 0 && dir[dir_len - 1] != '/');

    char *path = malloc(dir_len + needs_slash + stem_len + suffix_len + 1);
    if (path == NULL) {
        fprintf(stderr, "Memory allocation failed for output filename\n");
        exit(EXIT_FAILURE);
    }

    size_t pos = 0;
    if (dir_len > 0) {
        memcpy(path + pos, dir, dir_len);
        pos += dir_len;
    }
    if (needs_slash) {
        path[pos++] = '/';
    }
    memcpy(path + pos, basename, stem_len);
    pos += stem_len;
    memcpy(path + pos, suffix, suffix_len + 1);

    return path;
}

int lookup_flat(char* line, int lLen) {
    //if(lLen < cutoff){return 1;}
    int vver = 0;   // verdict !
    char* revLine=NULL;  
    char* target=NULL;
    int hFound=0;
    size_t scanLen=lLen-FR+1;
    for(int d = 1; d<=DD; d+=2){              // scanning density  
        size_t step = FR/d;                     
        for(int t = 0; t<2;t++){            // orientation 
            //int hFound=0;                            
            if(t==0){target=line;}
            else{target=revLine; if(target==NULL){revLine=revTrans(line);target=revLine;}}
            for(size_t i = d-1;i<scanLen;i+=step){
                uint64_t frN     = MurmurHash64A(target+i, FR, 0);
                uint32_t index   = (uint32_t)(frN >> 32);
                uint32_t hsValue = (uint32_t)(frN & 0xFFFFFFFF);                 
                // FlatHash specific lookup:
                size_t start = offsets[index];
                size_t end = offsets[index+1];
                //
                for(size_t ii = start; ii < end; ii++) {
                    if(values[ii] == hsValue) {
                        hFound+=1;
                        break; // for ii
                    }
                }
                if(hFound >= strictness){vver=1;break;} // for i 
            }
            if(vver){break;} // for t
        }
        if(vver){break;} // for d
    }
    if(revLine){free(revLine);revLine=NULL;}
    return vver;
}

    void refineFilter(char* line, size_t lLen){
        if(lLen < (size_t)FR){return ;}
        char* revLine=NULL;  
        char* target=NULL;
        size_t scanLen = lLen-FR+1;
        for(int t = 0; t<2;t++){  
            if(t==0){target=line;}
            else{revLine=revTrans(line);target=revLine;}
            for(size_t i = 0; i<scanLen;i++ ){
                uint64_t frN     = MurmurHash64A(target+i, FR, 0);
                uint32_t index   = (uint32_t)(frN >> 32);
                uint32_t hsValue = (uint32_t)(frN & 0xFFFFFFFF);
                pthread_mutex_t *lock = ht_lock_for(index);
                pthread_mutex_lock(lock);
                int remove_state = bucket_remove_value_locked(index, hsValue, htHS);
                pthread_mutex_unlock(lock);

                if (remove_state != 0) {
                    atomic_fetch_sub(&entryCount, 1);
                    if (remove_state == 2) {
                        atomic_fetch_sub(&counter, 1);
                    }
                }                
            }            
        }
        free(revLine);
        return;
    }
/*
int fastqRefine( char* Fbuffer ){
    int n = 1;
    int m = -1;

    char* line = Fbuffer;
    size_t lLen ;
    for (char *p = Fbuffer; *p; ++p) {
      if (*p == '\n') {
        *p = '\0';  // Replace newline with terminator // Now 'line' points to one line     
        if(readlimit && m == readlimit){break;}
        int nL = n % 4;
        switch(nL){
            case 1:
                if( line[0] != '@'){
                    fprintf(stderr,"File1 corrupt: line %d : %s\n",n,line); EXIT_FAILURE;
                }  
                m++;
                break;              
            case 2:
                lLen = strlen(line);
                refineFilter(line, lLen);        
                break;
            case 3:
                //RIGHTt2[m] = line ;   
                break;     
            case 0:
                //RIGHTq[m] = line ; 
                break;
        }
        n++;
        line = p + 1;  // Move to start of next line
      }
    }
    return m+1 ;
}
*/
int Ndetect (char* input, char* input2, size_t lLen) {
    
    int skip =0;
                int pos = -1;
                for (int i = 0; i < FR+1; i++) { // scan the beginning
                    if (input[i] == 'N') {
                        pos = i;
                    }
                }
                if (pos >= 0) {//mxmx++;
                    int remain = lLen - (pos + 2);
                    if (remain < cutoff) {
                        skip=1; //Nex++;
                        //continue;
                    }
                    else { 
                        input=&input[pos+2];
                        input2=&input2[pos+2];
                        lLen = remain;
                    }
                }
                pos = -1;
                for (int i = lLen - FR; i < (int)lLen; i++) { // scan the end
                    if (input[i] == 'N') {
                        pos = i;
                        break;
                    }
                }
                if (pos >= 0) {//mxmx++;
                    if (pos - 1  < cutoff) {
                        skip=1;//Nex++;
                        //continue;
                    }
                    else {
                        input[pos - 1] = '\0';
                        input2[pos - 1] = '\0';
                    }
                    //lLen=pos-1;
                } 
                // once more total scan
                int N = 0;
                for (int i = 0; i < (int)lLen; i++) { // scan all
                    if (input[i] == 'N') {
                        N++;
                    }
                } 
                if (N > 0) {
                    int r = lLen / N;
                    if( r < 50){skip=1;}
                }
    return skip;
}

    int import_binfiles (uint32_t*** htHS_ptr, char** fnameHS, int outfileSpecified){ 
        uint32_t** htHS = *htHS_ptr;
        if (htHS == NULL) {
            htHS = calloc(htSize , sizeof(uint32_t*));
            *htHS_ptr = htHS;
        }
        int fret=0;
        // read the Ffilter files to htHS 
        for( int fff = 0; fff<FilesCount;fff++){  // many filter files option
            
            if(fnameHS[fff]==NULL){continue;}
            time_t LoadStart = time(NULL);
            char* infile = fnameHS[fff];
            
            uint32_t* buffer = NULL;
            size_t buffer_capacity = 0;
            int Kmer = 0;
            uint32_t tableSize=0;
            uint32_t entrySize=0;
            FILE* file = fopen(infile, "rb");
            if (file == NULL) {
                if( (!outfileSpecified) && ( !FlatFile ) && (build || merge || refine)){
                    fprintf(stderr, "binfile not found. File treated as outfile\n");
                    outfile = infile;
                    return fret;
                }
                else{
                    fprintf(stderr, "binfile -F not found. Check the filenames.\n");
                    exit(EXIT_FAILURE);
                }
                //
            }
            // fprintf(stderr,"Importing: %s :",infile); 
            char binHeader [4];
            fread(&binHeader, sizeof(uint32_t), 1, file);
                if(memcmp(binHeader, "F64B", 4) != 0){
                    fprintf(stderr, "Wrong header ! %s is not F64 binfile.\n",infile);
                    exit(EXIT_FAILURE);
                }
            // read the FR first
            if(FR==0){
                fread(&(FR), sizeof(uint32_t), 1, file);
                if(FR > 5000){
                    fprintf(stderr, "Wrong K value ! %s is not F64 binfile.\n",infile);
                    exit(EXIT_FAILURE);
                }
            }
            else {
                fread(&(Kmer), sizeof(uint32_t), 1, file); 
                if( Kmer != FR){
                    fprintf(stderr, "Error: Files have different Kmer lenght: FR1 = %d ; FR%d = %d \n",FR,fff+1, Kmer);
                    exit(EXIT_FAILURE);
                } 
            }           
            // Read the size of the table 
            //  fread(&tableSize, sizeof(uint32_t), 1, file);
            if (fread(&tableSize, sizeof(uint32_t), 1, file) != 1) {
                fprintf(stderr, "Error reading table size\n");
                exit(EXIT_FAILURE);
            }  
            //fprintf(stderr, "Table size: %u\n",tableSize);    
            for (uint32_t i = 0; i < tableSize; i++) {
                fread(&(entrySize), sizeof(uint32_t), 1, file);
                if (entrySize > buffer_capacity) {
                    uint32_t *tmp = realloc(buffer, entrySize * sizeof(uint32_t));
                    if (!tmp) {
                        fprintf(stderr, "Memory allocation failed\n");
                        exit(EXIT_FAILURE);
                    }
                    buffer = tmp;
                    buffer_capacity = entrySize;
                }
                // buffer = realloc(buffer,(entrySize)*sizeof(uint32_t));
                fread((buffer),sizeof(uint32_t),entrySize,file);
                uint32_t index = buffer[0];
                if(htHS[index] == NULL){
                    htHS[index] = malloc(entrySize*sizeof(uint32_t));
                    htHS[index][0]=entrySize;
                    memcpy(&htHS[index][1], &buffer[1], (entrySize - 1) * sizeof(uint32_t) );
                    atomic_fetch_add(&counter, 1);
                    //for(uint32_t ii=1;ii<entrySize;ii++){ htHS[index][ii]=buffer[ii];}  
                } 
                else{
                    uint32_t oldSize = htHS[index][0];
                    uint32_t newSize = oldSize + entrySize - 1;
                    uint32_t *tmp = realloc(htHS[index], newSize * sizeof(uint32_t));
                    if (!tmp) {
                        fprintf(stderr, "Realloc failed\n");
                        exit(EXIT_FAILURE);
                    }
                    htHS[index] = tmp;
                    //htHS[index] = realloc(htHS[index],( oldSize + entrySize )*sizeof(uint32_t));
                    htHS[index][0] = newSize;
                    memcpy(&htHS[index][oldSize], &buffer[1], (entrySize - 1) * sizeof(uint32_t) );
                    // for(uint32_t ii=oldSize;ii < (oldSize + entrySize);ii++){  htHS[index][ii]=buffer[ii]; }  
                }
                atomic_fetch_add(&entryCount, entrySize-1);      
            }
            free(buffer);buffer=NULL;
            fclose(file);
            time_t now = time(NULL);
            if(takeTime){
                time_t Aeg = difftime(now, LoadStart);
                fprintf(stderr,"Loaded: %s; %zu indices\n",infile,atomic_load(&entryCount) );
                fprintf(stderr,"FilterFile %d read time = ",fff); // FYI: FilterFile read time = 0_00:02:49
                printAeg(Aeg);
                fprintf(stderr,"\n");
            }
            fret++;
        } // end of Filters import
        return fret;
    }

#define QSIZE 65536  // must be power of 2
#define MAX_LINE 256

typedef struct {
    char r1[4][MAX_LINE];
    char r2[4][MAX_LINE];
    int keep;
    uint64_t id;
    char pad[64];
} paired_read_t;

typedef struct {
    atomic_uint seq;
    void *data;
    char pad[64];
} mpmc_cell_t;

typedef struct {
    mpmc_cell_t *buffer;
    unsigned int mask;
    atomic_uint enqueue_pos;
    atomic_uint dequeue_pos;
} mpmc_queue_t;

int read_fastq_record(FILE *fp, char out[4][MAX_LINE]) {
    if (fp == NULL) {
        for (int i = 0; i < 4; i++) {
            out[i][0] = '\0';
        }
        return 1;
    }

    for (int i = 0; i < 4; i++) {
        if (!fgets(out[i], MAX_LINE, fp)){
            for (int j = i; j < 4; j++) {
                out[j][0] = '\0';
            }
            return 0;
        }
    }
    return 1;
}
/*
static inline void cpu_relax() {
    __asm__ volatile("pause" ::: "memory");  // x86
}
*/
static inline void cpu_relax() {
#if defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();

#elif defined(__aarch64__) || defined(__arm__)
    __asm__ volatile("yield");

#else
    ;
#endif
}

mpmc_queue_t read_queue;
mpmc_queue_t write_queue;

void mpmc_init(mpmc_queue_t *q) {
    q->buffer = malloc(sizeof(mpmc_cell_t) * QSIZE);
    q->mask = QSIZE - 1;

    for (unsigned int i = 0; i < QSIZE; i++) {
        atomic_store(&q->buffer[i].seq, i);
    }

    atomic_store(&q->enqueue_pos, 0);
    atomic_store(&q->dequeue_pos, 0);
}

void mpmc_destroy(mpmc_queue_t *q) {
    free(q->buffer);
    q->buffer = NULL;
    q->mask = 0;
}


int mpmc_enqueue(mpmc_queue_t *q, void *data) {
    mpmc_cell_t *cell;
    unsigned int pos = atomic_load(&q->enqueue_pos);

    while (1) {
        cell = &q->buffer[pos & q->mask];
        unsigned int seq = atomic_load(&cell->seq);
        int diff = (int)seq - (int)pos;

        if (diff == 0) {
            if (atomic_compare_exchange_weak(&q->enqueue_pos, &pos, pos + 1))
                break;
        } else if (diff < 0) {
            return 0; // full
        } else {
            pos = atomic_load(&q->enqueue_pos);
        }
    }

    cell->data = data;
    atomic_store(&cell->seq, pos + 1);

    return 1;
}

void* mpmc_dequeue(mpmc_queue_t *q) {
    mpmc_cell_t *cell;
    unsigned int pos = atomic_load(&q->dequeue_pos);

    while (1) {
        cell = &q->buffer[pos & q->mask];
        unsigned int seq = atomic_load(&cell->seq);
        int diff = (int)seq - (int)(pos + 1);

        if (diff == 0) {
            if (atomic_compare_exchange_weak(&q->dequeue_pos, &pos, pos + 1))
                break;
        } else if (diff < 0) {
            return NULL; // empty
        } else {
            pos = atomic_load(&q->dequeue_pos);
        }
    }

    void *data = cell->data;
    atomic_store(&cell->seq, pos + q->mask + 1);

    return data;
}

// Wait wrapper
paired_read_t* dequeue_block(mpmc_queue_t *q) {
    paired_read_t *item;
    int spins = 0;
    while ((item = mpmc_dequeue(q)) == NULL) {
        
        //fprintf(stderr, "Dequeue spinning to enqueue poison pill\n");
        if (spins < 100) {
            cpu_relax();          // fast spin
        } else if (spins < 1000) {
            sched_yield();        // light yield
        } else {
            struct timespec ts = {.tv_sec = 0, .tv_nsec = 1000};
            nanosleep(&ts, NULL); // backoff
        }
        spins++;
        
        //sched_yield();
    }
    return item;
}

// Memory pool
#define POOL_SIZE 100000

typedef struct {
    paired_read_t *items;
    paired_read_t **recycled;
    unsigned int recycled_capacity;
    atomic_uint top;
    atomic_uint recycled_head;
    atomic_uint recycled_tail;
    int multi_producer_recycle;
    pthread_mutex_t recycled_lock;
} pool_t;

pool_t pool;

// Init pool
void pool_init(int multi_producer_recycle) {
    pool.items = malloc(sizeof(paired_read_t) * POOL_SIZE);
    pool.recycled_capacity = POOL_SIZE + 1;
    pool.recycled = malloc(sizeof(paired_read_t*) * pool.recycled_capacity);
    atomic_store(&pool.top, 0);
    atomic_store(&pool.recycled_head, 0);
    atomic_store(&pool.recycled_tail, 0);
    pool.multi_producer_recycle = multi_producer_recycle;
    if (pool.multi_producer_recycle) {
        pthread_mutex_init(&pool.recycled_lock, NULL);
    }
}

static inline paired_read_t* pool_try_reuse() {
    if (pool.multi_producer_recycle) {
        paired_read_t *item = NULL;

        pthread_mutex_lock(&pool.recycled_lock);
        unsigned int tail = atomic_load_explicit(&pool.recycled_tail, memory_order_relaxed);
        unsigned int head = atomic_load_explicit(&pool.recycled_head, memory_order_acquire);

        if (tail != head) {
            item = pool.recycled[tail];
            tail++;
            if (tail == pool.recycled_capacity) {
                tail = 0;
            }
            atomic_store_explicit(&pool.recycled_tail, tail, memory_order_release);
        }
        pthread_mutex_unlock(&pool.recycled_lock);
        return item;
    }

    unsigned int tail = atomic_load_explicit(&pool.recycled_tail, memory_order_relaxed);
    unsigned int head = atomic_load_explicit(&pool.recycled_head, memory_order_acquire);

    if (tail == head) {
        return NULL;
    }

    paired_read_t *item = pool.recycled[tail];
    tail++;
    if (tail == pool.recycled_capacity) {
        tail = 0;
    }

    atomic_store_explicit(&pool.recycled_tail, tail, memory_order_release);
    return item;
}

static inline int pool_owns(paired_read_t *p) {
    uintptr_t addr = (uintptr_t)p;
    uintptr_t start = (uintptr_t)pool.items;
    uintptr_t end = start + sizeof(paired_read_t) * POOL_SIZE;
    return addr >= start && addr < end;
}

// Allocate pool
paired_read_t* pool_alloc() {
    paired_read_t *reused = pool_try_reuse();
    if (reused != NULL) {
        return reused;
    }

    unsigned int idx = atomic_fetch_add(&pool.top, 1);
    if (idx >= POOL_SIZE) return malloc(sizeof(paired_read_t));
    return &pool.items[idx];
}

// Return records to the pool after the writer has finished with them.
void pool_free(paired_read_t *p) {
    if (p == NULL) {
        return;
    }

    if (!pool_owns(p)) {
        free(p);
        return;
    }

    if (pool.multi_producer_recycle) {
        pthread_mutex_lock(&pool.recycled_lock);

        unsigned int head = atomic_load_explicit(&pool.recycled_head, memory_order_relaxed);
        unsigned int next = head + 1;
        if (next == pool.recycled_capacity) {
            next = 0;
        }

        unsigned int tail = atomic_load_explicit(&pool.recycled_tail, memory_order_acquire);
        if (next == tail) {
            pthread_mutex_unlock(&pool.recycled_lock);
            fprintf(stderr, "Pool recycle ring overflow\n");
            exit(EXIT_FAILURE);
        }

        pool.recycled[head] = p;
        atomic_store_explicit(&pool.recycled_head, next, memory_order_release);
        pthread_mutex_unlock(&pool.recycled_lock);
        return;
    }

    unsigned int head = atomic_load_explicit(&pool.recycled_head, memory_order_relaxed);
    unsigned int next = head + 1;
    if (next == pool.recycled_capacity) {
        next = 0;
    }

    unsigned int tail = atomic_load_explicit(&pool.recycled_tail, memory_order_acquire);
    if (next == tail) {
        fprintf(stderr, "Pool recycle ring overflow\n");
        exit(EXIT_FAILURE);
    }

    pool.recycled[head] = p;
    atomic_store_explicit(&pool.recycled_head, next, memory_order_release);
}

void pool_destroy() {
    if (pool.multi_producer_recycle) {
        pthread_mutex_destroy(&pool.recycled_lock);
    }
    free(pool.recycled);
    free(pool.items);
    pool.recycled = NULL;
    pool.items = NULL;
    pool.recycled_capacity = 0;
    pool.multi_producer_recycle = 0;
}

// Paired Streaming  Reader args
typedef struct {
    FILE *fp1;
    FILE *fp2;
    int n_consumers;
} reader_args_t;

// writer args:
typedef struct {
    FILE *out1;
    FILE *out2;
    int n_workers;
} writer_args_t;


static paired_read_t poison_pill;
_Atomic uint64_t total_reads;
_Atomic uint64_t reads_written_total;
_Atomic uint64_t reads_filtered_total;
int reverseR = 0; 
int compRead1(char *infile){
    int rev = 0;
    fprintf(stderr,"Analyzing: %s \n",infile);
    if(strstr(infile,"R1") != 0){rev=1;}
    return rev;
}

// Reader thread (zero-copy)
void* reader_thread(void *arg) {
    reader_args_t *a = (reader_args_t*)arg;
    uint64_t id = 0;
    
    while (1) {
        paired_read_t *pr = pool_alloc();

        if (!read_fastq_record(a->fp1, pr->r1)) {
            pool_free(pr);
            break;
        }
        if (!read_fastq_record(a->fp2, pr->r2)) {
            pool_free(pr);
            break;
        }
        pr->keep = 1;
        pr->id = id;
        id++;
        int reached_limit = (readlimit > 0 && id >= (uint64_t)readlimit);
        int spins = 0;
        while (!mpmc_enqueue(&read_queue, pr)) {
            //fprintf(stderr, "Reader waiting for free read_queue slot\n");

            if (spins < 100) {
                cpu_relax();          // fast spin
            } else if (spins < 1000) {
                sched_yield();        // light yield
            } else {
                struct timespec ts = {.tv_sec = 0, .tv_nsec = 1000};
                nanosleep(&ts, NULL); // backoff
            }
            spins++;
        }
        if (reached_limit) {
            break;
        }
        
    }

    atomic_store(&total_reads, id);
    //fprintf(stderr, "Reader exiting after %llu pairs\n", (unsigned long long)id);

    // poison pills
    for (int i = 0; i < a->n_consumers; i++){
        
        int spins = 0;
        while (!mpmc_enqueue(&read_queue, &poison_pill)) {
            //fprintf(stderr, "Reader waiting (spinning) to enqueue poison pill\n");
            if (spins < 100) {
                cpu_relax();          // fast spin
            } else if (spins < 1000) {
                sched_yield();        // light yield
            } else {
                struct timespec ts = {.tv_sec = 0, .tv_nsec = 1000};
                nanosleep(&ts, NULL); // backoff
            }
            spins++;
        } 
    }
    return NULL;
}

#define BATCH_SIZE 32768
#define FASTA_CHUNK_SIZE (1 << 23)  // 8 MB

static inline void enqueue_blocking(mpmc_queue_t *q, void *item) {
    int spins = 0;
    while (!mpmc_enqueue(q, item)) {
        if (spins < 100) {
            cpu_relax();
        } else if (spins < 1000) {
            sched_yield();
        } else {
            struct timespec ts = {.tv_sec = 0, .tv_nsec = 1000};
            nanosleep(&ts, NULL);
        }
        spins++;
    }
}

static inline void enqueue_fasta_chunks(char *seq, size_t len) {
    if (len == 0) {
        return;
    }

    size_t overlap = 0;
    if (FR > 1 && (size_t)(FR - 1) < FASTA_CHUNK_SIZE) {
        overlap = (size_t)(FR - 1);
    }
    size_t step = FASTA_CHUNK_SIZE - overlap;

    for (size_t start = 0; start < len; start += step) {
        enqueue_blocking(&read_queue, seq + start);
        if (len - start <= FASTA_CHUNK_SIZE) {
            break;
        }
    }
}

static inline size_t fasta_chunk_len(const char *chunk) {
    const char *stop = memchr(chunk, '\0', FASTA_CHUNK_SIZE);
    if (stop != NULL) {
        return (size_t)(stop - chunk);
    }
    return FASTA_CHUNK_SIZE;
}

void* worker_thread(void *arg) {
    (void)arg;
    paired_read_t *batch[BATCH_SIZE];
    int has_pair = (pairs == 2);

    while (1) {
        int count = 0;
        int saw_poison = 0;

        // dequeue batch
        for (int i = 0; i < BATCH_SIZE; i++) {
            paired_read_t *pr = mpmc_dequeue(&read_queue);
            if (!pr) break;

            if (pr == &poison_pill) {  // poison
                saw_poison = 1;
                break;
            }

            batch[count++] = pr;
        }

        if (count == 0) {
            if (saw_poison) {
                while (!mpmc_enqueue(&write_queue, &poison_pill)) {
                    cpu_relax();
                }
                return NULL;
            }
            //sched_yield();
            cpu_relax();
            continue;
        }

        // process batch
        for (int i = 0; i < count; i++) {
            paired_read_t *pr = batch[i];  
            if(has_pair) {
                size_t lLen2 = strcspn(pr->r2[1], "\n");
                if(lLen2 < (size_t)cutoff){pr->keep = 0;continue;}
                int skip2 = lookup_flat(pr->r2[1], lLen2);
                if(skip2){pr->keep = 0;continue;}
                if(detectN){skip2 = Ndetect(pr->r2[1],pr->r2[3],lLen2);
                    if(skip2){pr->keep = 0;continue;}
                }
            }
            size_t lLen1 = strcspn(pr->r1[1], "\n");
            if(lLen1 < (size_t)cutoff){pr->keep = 0;continue;}       
            int skip1 = lookup_flat(pr->r1[1], lLen1);
            if(skip1){pr->keep = 0;continue;}
            if(detectN){skip1 = Ndetect(pr->r1[1],pr->r1[3],lLen1);
                if(skip1){pr->keep = 0;}
            }     
        }
        if(printPositive){
            for (int i = 0; i < count; i++) {
                paired_read_t *pr = batch[i];
                pr->keep = !pr->keep;
            }
        }
        // enqueue batch
        for (int i = 0; i < count; i++) {
            while (!mpmc_enqueue(&write_queue, batch[i])) {
                //fprintf(stderr, "Worker waiting for free write_queue slot\n");
                cpu_relax();
            }
        }

        if (saw_poison) {
            while (!mpmc_enqueue(&write_queue, &poison_pill)) {
                cpu_relax();
            }
            return NULL;
        }
    }
}

void* fq_builder_thread(void *arg) {
    (void)arg;
    paired_read_t *batch[BATCH_SIZE];

    while (1) {
        int count = 0;
        int saw_poison = 0;

        // dequeue batch
        for (int i = 0; i < BATCH_SIZE; i++) {
            paired_read_t *pr = mpmc_dequeue(&read_queue);
            if (!pr) break;

            if (pr == &poison_pill) {  // poison
                saw_poison = 1;
                break;
            }

            batch[count++] = pr;
        }

        if (count == 0) {
            if (saw_poison) {
                return NULL;
            }
            //sched_yield();
            cpu_relax();
            continue;
        }

        // process batch
        for (int i = 0; i < count; i++) {
            paired_read_t *pr = batch[i];
            size_t lLen1 = strcspn(pr->r1[1], "\n");  
            if(build) {                   
                if(reverseR==0){
                    insertHash(FR,lLen1,pr->r1[1],htHS);
                }
                else {
                    pr->r1[1][lLen1]='\0';
                    char* revLine = revTrans(pr->r1[1]);
                    insertHash(FR,lLen1,revLine,htHS); 
                    free(revLine);revLine=NULL;
                }
            }
            if(refine){
                refineFilter(pr->r1[1],lLen1);
            }
            pool_free(pr);
        }

        if (saw_poison) {
            return NULL;
        }
    }
}

void* fa_builder_thread(void *arg) {
    (void)arg;

    while (1) {
        // FASTA chunks are large enough that one-at-a-time dequeueing gives
        // better cross-thread load sharing than a large local prefetch batch.
        char *chunk = (char*)mpmc_dequeue(&read_queue);
        if (!chunk) {
            cpu_relax();
            continue;
        }

        if (chunk == (char*)&poison_pill) {
            return NULL;
        }

        size_t llen = fasta_chunk_len(chunk);
        if(build){insertHash(FR, llen, chunk, htHS);}
        if(refine){refineFilter(chunk, llen); }
    }
}

#define MAX_PENDING (1 << 17)  // 131072
#define MASK (MAX_PENDING - 1)
#define WRITE_BUF_SIZE (1 << 20) // 1MB

typedef struct {
    char *buf;
    size_t pos;
} write_buffer_t;

typedef struct {
    paired_read_t **items;
    size_t capacity;
} pending_buffer_t;

static inline void wb_flush(FILE *f, write_buffer_t *wb) {
    fwrite(wb->buf, 1, wb->pos, f);
    wb->pos = 0;
}

static inline void wb_puts(FILE *f, write_buffer_t *wb, const char *s) {
    size_t len = strlen(s);

    if (wb->pos + len >= WRITE_BUF_SIZE) {
        wb_flush(f, wb);
    }

    memcpy(wb->buf + wb->pos, s, len);
    wb->pos += len;
}

static void pending_init(pending_buffer_t *pending, size_t capacity) {
    pending->items = calloc(capacity, sizeof(*pending->items));
    if (pending->items == NULL) {
        fprintf(stderr, "Memory allocation failed for writer pending buffer\n");
        exit(EXIT_FAILURE);
    }
    pending->capacity = capacity;
}

static void pending_destroy(pending_buffer_t *pending) {
    free(pending->items);
    pending->items = NULL;
    pending->capacity = 0;
}

static void pending_grow(pending_buffer_t *pending, uint64_t next_id, uint64_t incoming_id) {
    size_t new_capacity = pending->capacity;
    uint64_t distance = incoming_id - next_id;

    // Workers can finish far enough out of order to lap the initial ring.
    // Grow before inserting so read ids never alias the same pending slot.
    while (distance >= new_capacity) {
        if (new_capacity > ((size_t)-1) / 2) {
            fprintf(stderr, "ERROR: writer pending buffer cannot grow for read id %llu\n",
                    (unsigned long long)incoming_id);
            exit(EXIT_FAILURE);
        }
        new_capacity *= 2;
    }

    paired_read_t **new_items = calloc(new_capacity, sizeof(*new_items));
    if (new_items == NULL) {
        fprintf(stderr, "Memory allocation failed while growing writer pending buffer\n");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < pending->capacity; i++) {
        paired_read_t *cur = pending->items[i];
        if (cur == NULL) {
            continue;
        }
        if (cur->id < next_id || cur->id - next_id >= new_capacity) {
            fprintf(stderr, "ERROR: stale pending read id %llu while growing writer buffer\n",
                    (unsigned long long)cur->id);
            exit(EXIT_FAILURE);
        }
        size_t slot = cur->id & (new_capacity - 1);
        if (new_items[slot] != NULL) {
            fprintf(stderr, "ERROR: duplicate pending read ids %llu and %llu\n",
                    (unsigned long long)new_items[slot]->id,
                    (unsigned long long)cur->id);
            exit(EXIT_FAILURE);
        }
        new_items[slot] = cur;
    }

    free(pending->items);
    pending->items = new_items;
    pending->capacity = new_capacity;
}

static void pending_ensure_capacity(pending_buffer_t *pending, uint64_t next_id, uint64_t incoming_id) {
    if (incoming_id < next_id) {
        fprintf(stderr, "ERROR: stale read id %llu, writer already emitted through %llu\n",
                (unsigned long long)incoming_id,
                (unsigned long long)next_id);
        exit(EXIT_FAILURE);
    }

    if (incoming_id - next_id >= pending->capacity) {
        pending_grow(pending, next_id, incoming_id);
    }
}

void* writer_thread(void *arg) {
    writer_args_t *a = (writer_args_t*)arg;
    write_buffer_t wb1 = {malloc(WRITE_BUF_SIZE), 0};
    write_buffer_t wb2 = {NULL, 0};
    pending_buffer_t pending;
    uint64_t next_id = 0;
    int finished = 0;
    uint64_t reads_written = 0;
    uint64_t reads_filtered = 0;
    paired_read_t *pr = NULL;
    int has_pair = (a->out2 != NULL);

    pending_init(&pending, MAX_PENDING);

    if (has_pair) {
        wb2.buf = malloc(WRITE_BUF_SIZE);
    }
    while (1) {
        pr = dequeue_block(&write_queue);
        
        if (pr == &poison_pill) {
            finished++;
            // fprintf(stderr, "Writer got poison pill (%d/%d)\n", finished, a->n_workers);
            if (finished == a->n_workers) {
                uint64_t final = atomic_load(&total_reads);
                //fprintf(stderr, "Writer draining %llu pairs\n", (unsigned long long)final);
                while (next_id < final) {
                    size_t pending_slot = next_id & (pending.capacity - 1);
                    paired_read_t *cur = pending.items[pending_slot];
                    if (cur != NULL && cur->id != next_id) {
                        fprintf(stderr, "ERROR: writer pending mismatch: expected id %llu, found id %llu\n",
                                (unsigned long long)next_id,
                                (unsigned long long)cur->id);
                        exit(EXIT_FAILURE);
                    }

                    if (cur) {
                        if (cur->keep) {
                            for (int i = 0; i < 4; i++) {
                                wb_puts(a->out1, &wb1, cur->r1[i]);
                                if (has_pair) {
                                    wb_puts(a->out2, &wb2, cur->r2[i]);
                                }
                            }
                            reads_written++;
                        }else {
                            reads_filtered++;
                        }

                    pending.items[pending_slot] = NULL;
                    pool_free(cur);
                    next_id++;
                    } else {
                    // wait for missing earlier work
                        
                        int spins = 0;
            
                        //fprintf(stderr, "Writer waiting (spinning)\n");
                        if (spins < 100) {
                            cpu_relax();          // fast spin
                        } else if (spins < 1000) {
                            sched_yield();        // light yield
                        } else {
                            struct timespec ts = {.tv_sec = 0, .tv_nsec = 1000};
                            nanosleep(&ts, NULL); // backoff
                        }
                        spins++;
                    }
                }
                //fprintf(stderr, "Writer exiting\n");
                break;
            }
            continue;
        }
        pending_ensure_capacity(&pending, next_id, pr->id);
        uint64_t slot = pr->id & (pending.capacity - 1);

        if (pending.items[slot] != NULL) {
            fprintf(stderr, "ERROR: duplicate writer slot %llu for read ids %llu and %llu\n",
                    (unsigned long long)slot,
                    (unsigned long long)pending.items[slot]->id,
                    (unsigned long long)pr->id);
            exit(1);
        }
        pending.items[slot] = pr;

        // flush in order
        while (1) {
            size_t pending_slot = next_id & (pending.capacity - 1);
            paired_read_t *cur = pending.items[pending_slot];
            if (!cur) break;
            if (cur->id != next_id) {
                fprintf(stderr, "ERROR: writer pending mismatch: expected id %llu, found id %llu\n",
                        (unsigned long long)next_id,
                        (unsigned long long)cur->id);
                exit(EXIT_FAILURE);
            }

            pending.items[pending_slot] = NULL;

            if (cur->keep) {
                for (int i = 0; i < 4; i++) {
                    wb_puts(a->out1, &wb1, cur->r1[i]);
                    if (has_pair) {
                        wb_puts(a->out2, &wb2, cur->r2[i]);
                    }
                }
                reads_written++;
            } else {
                reads_filtered++;
            }

            pool_free(cur);
            next_id++;
        }
    }

    wb_flush(a->out1, &wb1);
    if (has_pair) {
        wb_flush(a->out2, &wb2);
    }
    atomic_store(&reads_written_total, reads_written);
    atomic_store(&reads_filtered_total, reads_filtered);
    free(wb1.buf);
    free(wb2.buf);
    pending_destroy(&pending);
  
    return NULL;
    
}


int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "%s: Usage: F64 <command> [options]\n",argv[0]);
        return 1;
    }

    if(argv[1][0] != '-'){                                  // if command is given
        if(strcmp(argv[1], "stats") == 0){stats = 1;}
        else if(strcmp(argv[1], "filter") == 0){ filter = 1;}
        else if(strcmp(argv[1], "build") == 0){ build = 1;}
        else if(strcmp(argv[1], "merge") == 0){ merge = 1;}
        else if(strcmp(argv[1], "refine") == 0){ refine = 1;}
        else if(strcmp(argv[1], "help") == 0){ help = 1;}
        else {fprintf(stderr, "Unknown command: %s\n",argv[1]); exit(1);}
        parse_args(argc - 1, argv + 1);
        fprintf(stderr, "\nRunning: %s: %s \n",argv[0],argv[1]);
    }
    else {parse_args(argc, argv);
        fprintf(stderr, "\nRunning: %s filter\n",argv[0]);
    }

    if(help){fprintf(stderr,"%s",msg2);exit(0);}

    time_t algusAeg = time(NULL);
    int Aeg = 0;
    time_t now = time(NULL);

  // licencing time test
    time_t nyyd;
    struct tm expiry = {0};
    // Get current time
    time(&nyyd);
    // Set expiry date (YYYY, MM-1, DD)
    expiry.tm_year = 2027 - 1900;  // years since 1900
    expiry.tm_mon  = 4;            // May (0 = Jan)
    expiry.tm_mday = 1;            // 1st
    // Convert expiry date to time_t
    time_t expiry_time = mktime(&expiry);
    // Compare
    if (difftime(nyyd, expiry_time) > 0) {
        printf("Program expired.\n");
        exit(1);
    }
  // end time test

    //n_values = 0;
    int cpus = 1;

    char* env = getenv("SLURM_CPUS_PER_TASK");
    if (env) {
        cpus = atoi(env);
    } else {
        cpus = sysconf(_SC_NPROCESSORS_ONLN);
    }

    N_WORKERS = cpus - 2;
    if (N_WORKERS < 1) {N_WORKERS = 1;}  // safety
    N_BUILDERS =cpus-1;
    if (N_BUILDERS < 1) {N_BUILDERS = 1;}

  // memory
    /*
    // for macOSX
    int64_t mem;
    size_t len = sizeof(mem);
    sysctlbyname("hw.memsize", &mem, &len, NULL, 0);

    uint64_t Mmem = (uint64_t) mem / (1024 * 1024);
    fprintf(stderr,"Total RAM: %lld MB\n", mem / (1024 * 1024));
    if (Mmem < 16777216){
        fprintf(stderr,"Not enough RAM to run F64 !\n");
        exit(0);
    }
    */



    if (stats) {  // only read filterfiles, no import
        n_values = 0;
        if(FlatFile){
            FILE* file = fopen(FlatFile, "rb");
            if (file == NULL) {
                fprintf(stderr, "Error opening FlatFile for reading.\n");
                exit(EXIT_FAILURE);
            }
            char flatHeader [4];
            fread(&flatHeader, sizeof(uint32_t), 1, file);            // read FR 
            if(memcmp(flatHeader, "F64F", 4) != 0){
                fprintf(stderr, "Wrong header ! %s is not f64flat file.\n",FlatFile);
                exit(EXIT_FAILURE);
            }  
            fread(&(FR), sizeof(uint32_t), 1, file);            // read FR 
            if(FR > 5000){
                fprintf(stderr, "Wrong K value ! %s is not F64 file.\n",FlatFile);
                exit(EXIT_FAILURE);
            }
            uint64_t fTest = 0;
            fread(&fTest,sizeof(uint64_t), 1, file);                      // test file type
            if(fTest != osSize){ 
                fprintf(stderr, "Wrong file type ! %s is not F64 flat file.\n",FlatFile);
                exit(EXIT_FAILURE);
            }
            fseek(file, osSize*sizeof(uint32_t),SEEK_CUR );     // skip offsets    
            fread(&n_values, sizeof(uint64_t), 1, file);        // the number of values
            fclose(file);   
            printf("{\n");
            printf("  \"FlatFile\": %s,\n", FlatFile);
            printf("  \"K\": %d,\n", FR);
            printf("  \"Indices\": %llu,\n", n_values);
            printf("}\n");
        }
        
        for( int fff = 0; fff<FilesCount;fff++){  // many filter files option
            if(fnameHS[fff]==NULL){continue;}
            //time_t LoadStart = time(NULL);
            char* infile = fnameHS[fff];
            
            entryCount = 0;
            //fprintf(stderr,"Importing: %s \n",infile); 
            uint32_t tableSize=0;
            uint32_t entrySize=0;  
            FILE* file = fopen(infile, "rb");
            if (file == NULL) {
                fprintf(stderr, "Error opening binfile for reading.\n");
                exit(EXIT_FAILURE);
            }
            char binHeader[4];
            fread(&binHeader, sizeof(uint32_t), 1, file);            // read FR 
            if(memcmp(binHeader, "F64B", 4) != 0){
                fprintf(stderr, "Wrong header ! %s is not F64bin file.\n",infile);
                exit(EXIT_FAILURE);
            }         
            fread(&(FR), sizeof(uint32_t), 1, file);            // read FR 
            if(FR > 5000){
                fprintf(stderr, "Wrong K value ! %s is not F64 file.\n",infile);
                exit(EXIT_FAILURE);
            }
            // fprintf(stderr,"  \"K\": %d,\n", FR);
            fread(&tableSize, sizeof(uint32_t), 1, file);       // read Table size
            uint64_t* counterT = NULL;
            size_t counterTsize = 0;
            if(TableReport){
                counterTsize = 16;
                counterT = calloc(counterTsize,sizeof(uint64_t));
                if (counterT == NULL) {
                    fprintf(stderr, "Memory allocation failed\n");
                    exit(EXIT_FAILURE);
                }
                counterT[0] = (uint64_t)htSize - tableSize;
            }
            for (uint32_t i = 0; i < tableSize; i++) {
                fread(&(entrySize), sizeof(uint32_t), 1, file); // read Entry size
                     // count indices   
                entryCount+=entrySize-1; 
                if(TableReport){
                    add_bucket_report_count(&counterT, &counterTsize, entrySize-1);
                }
                fseek(file,entrySize*sizeof(uint32_t),SEEK_CUR);
            }
            fclose(file);  
             
            printf("{\n");
            printf("  \"FilterFile %d\": %s,\n",fff+1, infile);
            printf("  \"K\": %d,\n", FR);
            printf("  \"Indices\": %zu,\n", entryCount);
            printf("}\n");
            n_values+=entryCount; 
            if(TableReport){
                print_bucket_report(counterT, counterTsize);
                free(counterT);
            }
        }
        if(FilesCount > 0){
            printf("{\n");
            printf("  \"TotalIndices\": %llu,\n", n_values);
            printf("}\n");
        }
        if(takeTime){
            now = time(NULL);
            int Aeg = difftime(now, algusAeg);
            fprintf(stderr,"Read time = ");    
            printAeg(Aeg);
            fprintf(stderr,"\n");
        }
        exit(0);
    } 
    int outfileSpecified = 0;
    if(build || merge || refine){
        
        if(outfile){ outfileSpecified = 1;}
        if(! outfile){
            outfile=fnameHS[0];
        }
        size_t outFileLen =strlen(outfile);
        outfile=strdup(outfile); 

        atomic_store(&counter, 0);
        atomic_store(&entryCount, 0);
        init_ht_locks();
        htHS = calloc(htSize, sizeof(uint32_t*)); fprintf(stderr,"htHS initialized\n");
           
        if(fnameHS){
            int ret = import_binfiles(&htHS,fnameHS,outfileSpecified);
            if(ret == 0){fprintf(stderr,"Updating...\n");}
        }
        if(FlatFile){
            if (access(FlatFile, F_OK) == 0){ // the output file exists
                fprintf(stderr,"Flat Files can not be updated or merged\n");
                exit(0);
            }
            if(strstr(FlatFile,".f64flat")==0){
                size_t flatFileLen =strlen(FlatFile);
                FlatFile=strdup(FlatFile);
                FlatFile=realloc(FlatFile, (flatFileLen + 11)*sizeof(char)); 
                strncat(FlatFile,".f64flat",8);
                FlatFile[flatFileLen + 9]='\0';
            }
        }

        if(! FlatFile && access(outfile, F_OK) != 0){// the output file does not exist
            if(strstr(outfile,".f64bin")==0){
                outfile=realloc(outfile, (outFileLen + 11)*sizeof(char)); 
                strncat(outfile,".f64bin",6);
                outfile[outFileLen + 7]='\0';
            }
        }
        // check if the file exists
        if (access(outfile, F_OK) == 0){ // the output file exists
            if(outfileSpecified == 0){
                       // if outfile exists, change the name
                int outFileLen = strlen(outfile);
                char* version = ".v";
                char* vFound=strstr(outfile,version);
                char* extension = strstr(outfile,"f64");
                if(vFound==0){
                    size_t nameLen = 0;
                    if(extension){nameLen = extension - outfile;}
                    outfile=realloc(outfile,(outFileLen+4)*sizeof(char) );
                    outfile[nameLen]='\0';
                    strncat(outfile,".v1",3);
                    outfile[nameLen+3]='\0';
                    strncat(outfile,".f64bin",6);
                    outfile[nameLen+9]='\0';
                    outFileLen=strlen(outfile);
                }
                else {
                    size_t scanStart = vFound+2 - outfile;
                    int versionLen = 1;
                    while (outfile[scanStart++] >= '0' && outfile[scanStart++]<='9'){versionLen++;}
                    char* versionString = calloc(versionLen,sizeof(char));
                    strncpy(versionString,outfile+scanStart,versionLen);
                    int versionNr=atoi(versionString);
                    versionNr++;
                    char* newNr=calloc(5,sizeof(char));
                    snprintf(newNr,4,"%d",versionNr);
                    int newNrLen=strlen(newNr);
                    outfile=realloc(outfile,(outFileLen+newNrLen+1)*sizeof(char));
                    outfile[scanStart-2]='\0';
                    strncat(outfile,newNr,newNrLen);
                    outfile[scanStart-2+newNrLen]='\0';
                    strncat(outfile,"f64bin",6);
                    outfile[scanStart-2+newNrLen+6]='\0';
                    outFileLen =strlen(outfile);
                }
            }
            else{
                fprintf(stderr,"Outfile is replaced, not updated !\n");
            }
        }
        
        if (merge){
            if(TableReport){
                report_hashtable (&htHS);
            }
            if(json){
                printf("{\n");
                printf("  \"FileOut\": %s,\n", outfile);
                printf("  \"K\": %u,\n", FR);
                printf("  \"indices\": %zu,\n", atomic_load(&entryCount));
                printf("  \"time_sec\": %u,\n", Aeg); 
                printf("}\n");

            }
            if(FlatFile){    
                convert_to_FlatHash();
                writeFlatOutfile(FlatFile);
                fprintf(stderr,"Outfile: %s ; %llu indices\n",FlatFile,n_values);
            }

            if(outfileSpecified || (! FlatFile )){    
                write_outfile (outfile,FR,(uint32_t)atomic_load(&counter),htHS);
            }
            
            exit(0);
        }
       
        // build and refine:
        // fprintf(stderr,"Builder Threads: %d\n",N_BUILDERS);
        if( FR == 0){FR = 25;}
        for( int fff = 0; fff<InputCount;fff++){  // start of batch cycle, reading fasta and fastq
            // time_t Cstart = time(NULL);
            char* infile = fnamesInput[fff];
            if( !infile ){fprintf(stderr,"No files to convert\n");exit(0); }
            size_t line_length = 0;
            ssize_t read;
            char* line = NULL;
            int fileType = 0;

            FILE* file1 = fopen(infile, "r");
            if (file1 == NULL) {fprintf(stderr,"Error opening file1\n");return 1;}
            read = getline(&line, &line_length, file1);
            fclose(file1);
            if (read == -1 || line == NULL) {
                fprintf(stderr,"Empty or unreadable input: %s\n",infile);
                free(line);
                return 1;
            }
            if(line[0]=='@'){fileType=1;}
            else if(line[0]=='>'){fileType=0;}
            else {
                fprintf(stderr,"Unknown file type: %s\n",infile);
                free(line);
                return 1;
            }
            free(line);
            line = NULL;

            // fastq files
            if(fileType == 1) { // fastq files, analyzed one by one as they come
                
                FILE* file1 = fopen(infile, "r");
                if (file1 == NULL) {fprintf(stderr,"Error opening file1\n");return 1;}
                FILE* file2 = NULL;
                reverseR = compRead1(infile);
                atomic_store(&total_reads, 0);                
                mpmc_init(&read_queue);                
                pool_init(1);
                //fprintf(stderr,"Queues initiated\n");
                
                pthread_t reader, builders[N_BUILDERS];
                reader_args_t* rargs = malloc(sizeof(reader_args_t));
                rargs->fp1 = file1;   
                rargs->fp2 = NULL;
                rargs->n_consumers = N_BUILDERS;
                pthread_create(&reader, NULL, reader_thread, rargs);
                for (int i = 0; i < N_BUILDERS; i++){
                    pthread_create(&builders[i], NULL, fq_builder_thread, NULL);
                }
                               
                pthread_join(reader, NULL);
                for (int i = 0; i < N_BUILDERS; i++){
                    pthread_join(builders[i], NULL);
                }

                if (file1 != NULL && file1 != stdin) { fclose(file1); }
                if (file2 != NULL && file2 != stdin) { fclose(file2); }
                
                free(rargs);
                mpmc_destroy(&read_queue);
                pool_destroy();
                
            }
            // end of fastq 

            if(fileType==0){ // fasta files
                fprintf(stderr,"Analyzing: %s \n:",infile);
                char* FBuff = openFasta(infile);
                if (FBuff == NULL) {
                    fprintf(stderr,"Error reading fasta file\n");
                    exit(EXIT_FAILURE);
                }
                int m = -1;
                char* line = FBuff;
                mpmc_init(&read_queue);
                // fprintf(stderr," fa Builder Threads: %d\n",N_BUILDERS);
                pthread_t builders[N_BUILDERS];
                for (int i = 0; i < N_BUILDERS; i++) {
                    pthread_create(&builders[i], NULL, fa_builder_thread, NULL);
                }
                for (char *p = FBuff; *p; ++p) {
                    if (*p == '\n') {
                        *p = '\0';  // Replace newline with terminator
                        size_t lLen = strlen(line);
                        if(lLen < 2){line = p + 1;continue;}
                        if(readlimit && m==readlimit){break;}
                        if (line[0] == '>') {m++;line = p + 1;continue;}
                        else{
                            enqueue_fasta_chunks(line, lLen);
                        }
                        line = p + 1;
                    }
                }
                for (int i = 0; i < N_BUILDERS; i++) {
                    enqueue_blocking(&read_queue, &poison_pill);
                }
                for (int i = 0; i < N_BUILDERS; i++) {
                    pthread_join(builders[i], NULL);
                }
                mpmc_destroy(&read_queue);
                free(FBuff);FBuff=NULL;
                //fileStat(infile);
                //int numFCont = m; //filestatRes[2] % 0xffffffff;
                fprintf(stderr,"Fasta reads: %d\n",m+1);
            }
            // end of fasta
            //readlimit=0;
        }

        // fprintf(stderr,"htHS filled\n");
        time_t now = time(NULL);
        if(takeTime){
            int Aeg = difftime(now, algusAeg);
            fprintf(stderr,"htHS fill time = ");                // FYI: htHS fill time = 0_00:17:20 for human genome
            printAeg(Aeg);
            fprintf(stderr,"\n");
        }
        fprintf(stderr,"Indices: %zu\n",atomic_load(&entryCount));      // FYI: Fragments: 2 004 547 527
    
        // write htHS to a file as binary data
        if(FlatFile){
            convert_to_FlatHash();
            writeFlatOutfile(FlatFile);
            fprintf(stderr,"Outfile: %s ; %llu indices\n",FlatFile,n_values);
        }
        if(outfileSpecified || (! FlatFile )) {
            write_outfile (outfile,FR,(uint32_t)atomic_load(&counter),htHS);
        }
        // table report
        if(TableReport){
            report_hashtable (&htHS);
        }

        if(json){
            time_t now = time(NULL);
            Aeg = difftime(now, algusAeg);
            printf("{\n");
            printf("  \"FileOut\": %s,\n", outfile);
            printf("  \"K\": %u,\n", FR);
            printf("  \"indices\": %zu,\n", atomic_load(&entryCount));
            printf("  \"time_sec\": %u,\n", Aeg); 
            printf("}\n");

        }

        exit(0);
    }
    // end of build || refine || merge

    // filter == Default behaviour
    
    if(! FlatFile){
      //initialize htHS
        int ret = 0;
        ret = import_binfiles (&htHS, fnameHS, outfileSpecified);
        now = time(NULL);
        Aeg = difftime(now, algusAeg);
        if(ret){
            fprintf(stderr,"Loaded %d files; %zu indices\n",ret, atomic_load(&entryCount));
        }
        if(ret != FilesCount){ fprintf(stderr," %d files missing\n",FilesCount - ret); EXIT_FAILURE;}
        if(takeTime){
            fprintf(stderr,"Total loading time = "); 
            printAeg(Aeg);
            fprintf(stderr,"\n");
        }
        if(TableReport){
            report_hashtable (&htHS);
        }
      //  build_flat hash from **htHS for cache-efficient lookup
        convert_to_FlatHash();
        //fprintf(stderr, "Values: %llu\n",n_values);

      // erase old htHS
        for(uint32_t i=0;i<htSize;i++){if(htHS[i]){free(htHS[i]);htHS[i]=NULL;}}
        free(htHS);htHS=NULL;   
        /*
        time_t conversionTime = time(NULL);
        Aeg = difftime(conversionTime,now);
        fprintf(stderr,"FlatHash conversion time = "); 
        printAeg(Aeg);
        fprintf(stderr,"\n");
        */
    // end of FlatHash building
    }
    else if(FlatFile && FilesCount){ 
        fprintf(stderr,"Flat files and bin files can not be merged");
        exit(0);
    }
    else{
        readFlatfile(FlatFile);
        fprintf(stderr,"Imported: %s\n",FlatFile);
        if(takeTime){
            time_t flatReadTime = time(NULL);
            Aeg = difftime(flatReadTime,algusAeg);
            fprintf(stderr,"FlatHash Load time = "); 
            printAeg(Aeg);
            fprintf(stderr,"\n");
        }
    }

    time_t process_time = time(NULL);
    char* fname1 = NULL;
    char* fname2 = NULL;

    if(InputCount == 1){pairs = 1;}
    if (pairs != 1 && pairs != 2) {
        fprintf(stderr, "Error: -p must be 1 for single-end input or 2 for paired input.\n");
        exit(EXIT_FAILURE);
    }
    if (stdout2_path != NULL && pairs != 2) {
        fprintf(stderr, "Error: --stdout2 can only be used with paired input.\n");
        exit(EXIT_FAILURE);
    }
    if (pairs == 2 &&
        stdout1_path != NULL &&
        stdout2_path != NULL &&
        strcmp(stdout1_path, "-") == 0 &&
        strcmp(stdout2_path, "-") == 0) {
        fprintf(stderr, "Error: --stdout1 and --stdout2 cannot both write to stdout.\n");
        exit(EXIT_FAILURE);
    }
    if (use_stdin) {
        if (stdin1_path == NULL) {
            fprintf(stderr, "Error: --stdin1 must be provided for streaming input.\n");
            exit(EXIT_FAILURE);
        }
        if (pairs == 2 && stdin2_path == NULL) {
            fprintf(stderr, "Error: both --stdin1 and --stdin2 must be provided for paired streaming input.\n");
            exit(EXIT_FAILURE);
        }
        disableReset = 1;
    } else {
        if (InputCount == 0) {
            fprintf(stderr, "Error: provide input with -i <R1.fastq> [R2.fastq ...] or use --stdin1/--stdin2.\n");
            exit(EXIT_FAILURE);
        }
        if (pairs == 2 && (InputCount % 2) != 0) {
            fprintf(stderr, "Error: paired mode (-p 2) requires an even number of input files.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    int cycle_limit = use_stdin ? 1 : InputCount;
    int cycle_step = use_stdin ? 1 : pairs;

    /// input files cycle
  for( int fff = 0; fff<cycle_limit;fff+=cycle_step){  // start of batch cycle
    time_t Cstart = time(NULL);
    outfile1 = NULL;
    outfile2 = NULL;
    fname1 = use_stdin ? stdin1_path : fnamesInput[fff];
    fname2 = NULL;
    if (pairs == 2) {
        fname2 = use_stdin ? stdin2_path : fnamesInput[fff+1];
    }

    if (fname2) {
        outfile1 = stdout1_path ? strdup(stdout1_path) : derive_filtered_output_path(fname1, printPositive ? "stdin1.Positive.fastq" : "stdin1.Filtered.fastq");
        outfile2 = stdout2_path ? strdup(stdout2_path) : derive_filtered_output_path(fname2, printPositive ? "stdin2.Positive.fastq" : "stdin2.Filtered.fastq");
    } else {
        outfile1 = stdout1_path ? strdup(stdout1_path) : derive_filtered_output_path(fname1, printPositive ? "stdin.Positive.fastq" : "stdin.Filtered.fastq");
    }

    if (stdout1_path == NULL && access(outfile1, F_OK) == 0){ // the output file exists
        fprintf(stderr,"Already processed: %s \n",outfile1);
        free(outfile1);
        free(outfile2);
        outfile1 = NULL;
        outfile2 = NULL;
        continue;
    }
    FILE *fp1 = NULL;
    FILE *filep = NULL;
    if (use_stdin) {
        fp1 = open_input_stream(fname1, "input 1");
        filep = fp1;
    } else {
        filep = fopen(fname1, "r");
        if (filep == NULL) {
            fprintf(stderr, "Error opening input file: %s\n", fname1);
            exit(EXIT_FAILURE);
        }
    }
    int firstchar;
    int fileType = 1;
    do {
        firstchar = fgetc(filep);
    } while (firstchar == '\n' || firstchar == '\r' || firstchar == ' ' || firstchar == '\t');
    if (firstchar == EOF) {
        printf("Empty file or error\n");
    } else if (firstchar == '>') {
        fileType=0;
        disableReset = 1;
        cycle_step = 1;
        //printf("FASTA format (starts with '>')\n");
    } else if (firstchar == '@') {
        fileType=1;
        // printf("FASTQ format (starts with '@')\n");
    } else {
        printf("Unknown format\n");
    }
    if (firstchar != EOF && ungetc(firstchar, filep) == EOF) {
        fprintf(stderr, "Error resetting input stream after format check: %s\n", fname1);
        exit(EXIT_FAILURE);
    }
    

    if(disableReset == 0){          // check the readlenghts and adapt parameters
       
        char* line = NULL;
        size_t line_length = 0;
        ssize_t read;
        int totalL=0;
        int m = 0;
        int n = 0;
        while (m < 100000 && (read = getline(&line, &line_length, filep)) != -1) {
            if(n % 4 == 1){
                int lLen = (int)read;
                while (lLen > 0 && (line[lLen - 1] == '\n' || line[lLen - 1] == '\r')) {
                    lLen--;
                }
                totalL+=lLen;
                m++;
            }
            n++;
        }
        free(line);
        if (m == 0) {
            fprintf(stderr,"Error: no FASTQ reads found while calculating read length.\n");
            exit(EXIT_FAILURE);
        }
        int averageReadLength = totalL/m;
        fprintf(stderr,"Analyzed %d reads\n",m);
        fprintf(stderr,"Average Read length: %d\n",averageReadLength);
        cutoff = (int)(averageReadLength*0.7);
        if(cutoff < 40){cutoff=40;}
        if(cutoff > 99){cutoff=70;}
        strictness=FR > 0 ? averageReadLength/FR/2 : 1;
        if(strictness < 1){strictness=1;}
        fprintf(stderr,"Read length cutoff set to %d\n",cutoff);
        fprintf(stderr,"Strictness set to %d\n",strictness);

    }
    if (!use_stdin) {
        fclose(filep);
    } else if (fileType == 0) {
        fprintf(stderr, "Error: streaming FASTA input is not supported; use FASTQ with --stdin1/--stdin2 or process FASTA with -i.\n");
        if (fp1 != NULL && fp1 != stdin) { fclose(fp1); }
        exit(EXIT_FAILURE);
    }
    // single thread process for fasta files.
        if(fileType == 0){
            if(readlimit == 0){ // get the inputfile linecount via wc -l    
                fileStat(fname1);
                readlimit = (uint32_t)filestatRes[2];
            }
            char** TITLE = calloc(readlimit+1, sizeof(char*));
            char** seqs = calloc(readlimit+1, sizeof(char*));
            char* FBuff=openFasta(fname1);
            if (FBuff == NULL) {
                fprintf(stderr,"Error reading fasta file\n");
                exit(EXIT_FAILURE);
            }
            int m = -1;
            char* line = FBuff;
                
            for (char *p = FBuff; *p; ++p) {
                if (*p == '\n') {
                    *p = '\0';  // Replace newline with terminator
                    size_t lLen = strlen(line);
                    if(lLen < 2){line = p + 1;continue;}
                    
                    if (line[0] == '>') {
                        m++;
                        if(readlimit && m==readlimit){break;}
                        TITLE[m]=line;line = p + 1;continue;
                    }
                    else{
                        if(lLen < (size_t)cutoff){line = p + 1;continue;}
                        int skip = lookup_flat(line, lLen);
                        if(printPositive){
                            if(skip){
                                seqs[m]=line;

                            }
                        }
                        else{
                            if(!skip){
                                seqs[m]=line;
                            }
                        }
                    }
                    line = p + 1;
                }
            }
            FILE *out1 = open_output_stream(outfile1, "output 1");
            // print out
            for(int i = 0; i<readlimit;i++){
                if(seqs[i]){
                    fprintf(out1,"%s\n%s\n",TITLE[i],seqs[i]);
                }
            }
            if (out1 != stdout) { fclose(out1); }
            // free
            free(FBuff);FBuff=NULL;
            free(TITLE);TITLE=NULL;
            free(seqs);seqs=NULL;
            
            readlimit = 0;
        }
    else {
    if (fp1 == NULL) {
        fp1 = open_input_stream(fname1, "input 1");
    }
    FILE *fp2 = NULL;
    if( fname2) {fp2 = open_input_stream(fname2, "input 2");}

    FILE *out1 = open_output_stream(outfile1, "output 1");
    FILE *out2 = NULL;
    if( fname2) {out2 = open_output_stream(outfile2, "output 2");}
    if (out1 == NULL || (fname2 && out2 == NULL)) {
        fprintf(stderr, "Error opening output files.\n");
        exit(EXIT_FAILURE);
    }

    atomic_store(&total_reads, 0);
    atomic_store(&reads_written_total, 0);
    atomic_store(&reads_filtered_total, 0);
    mpmc_init(&read_queue);
    mpmc_init(&write_queue);
    pool_init(0);
    //fprintf(stderr,"Queues initiated\n");
    //fprintf(stderr,"Worker Threads: %d\n",N_WORKERS);
    pthread_t reader, workers[N_WORKERS], writer;
    reader_args_t* rargs = malloc(sizeof(reader_args_t));
    rargs->fp1 = fp1;   
    rargs->fp2 = fp2;
    rargs->n_consumers = N_WORKERS;
    pthread_create(&reader, NULL, reader_thread, rargs);
    //fprintf(stderr,"Reader created\n");
    writer_args_t* wargs= malloc(sizeof(writer_args_t));
    wargs->n_workers=N_WORKERS;
    for (int i = 0; i < N_WORKERS; i++){
        pthread_create(&workers[i], NULL, worker_thread, NULL);
    }
    wargs-> out1 = out1;  
    wargs-> out2 = out2;
    //fprintf(stderr,"Workers created\n");
    pthread_create(&writer, NULL, writer_thread, wargs);
    //fprintf(stderr,"Writer created\n");
    
    pthread_join(reader, NULL);
    for (int i = 0; i < N_WORKERS; i++){
        pthread_join(workers[i], NULL);
    }
    pthread_join(writer, NULL);


    //fprintf(stderr,"Reads processed:  %llu \n", (unsigned long long)atomic_load(&total_reads));
    //fprintf(stderr,"Reads skipped:  %llu \n", (unsigned long long)atomic_load(&reads_filtered_total));
    //fprintf(stderr,"Reads kept:  %llu \n", (unsigned long long)atomic_load(&reads_written_total));

    if (fp1 != NULL && fp1 != stdin) { fclose(fp1); }
    if (fp2 != NULL && fp2 != stdin) { fclose(fp2); }
    if (out1 != NULL && out1 != stdout) { fclose(out1); }
    if (out2 != NULL && out2 != stdout) { fclose(out2); }
    free(rargs);
    free(wargs);
    mpmc_destroy(&read_queue);
    mpmc_destroy(&write_queue);
    pool_destroy();
} // end of filetype else
    now = time(NULL);
    Aeg = difftime(now, Cstart);
    //fprintf(stderr,"Processed: %s\n",fname1);
    if(takeTime){
        fprintf(stderr,"File Processing time = ");
        printAeg(Aeg);
        fprintf(stderr,"\n");
    }
    if(json){
        printf("{\n");
        printf("  \"K\": %d,\n", FR);
        printf("  \"FileOut1\": \"%s\",\n", outfile1 ? outfile1 : "");
        printf("  \"FileOut2\": \"%s\",\n", outfile2 ? outfile2 : "");
        printf("  \"Reads processed\": %llu,\n",(unsigned long long)atomic_load(&total_reads));
        printf("  \"Reads skipped\": %llu,\n",(unsigned long long)atomic_load(&reads_filtered_total));
        printf("  \"Reads kept\": %llu,\n",(unsigned long long)atomic_load(&reads_written_total));
        printf("  \"time_sec\": %u,\n", Aeg); 
        printf("}\n");
    }
    atomic_store(&reads_written_total, 0);
    atomic_store(&reads_filtered_total, 0);
    atomic_store(&total_reads, 0);
    free(outfile1);
    free(outfile2);
    outfile1 = NULL;
    outfile2 = NULL;
  } // end of infile cycle
    now = time(NULL);
    if(takeTime){
        Aeg = difftime(now, process_time);
        fprintf(stderr,"Processing time = ");
        printAeg(Aeg);
        fprintf(stderr,"\n");
        Aeg = difftime(now, algusAeg);
        fprintf(stderr,"Total time = ");
        printAeg(Aeg);
        fprintf(stderr,"\n");
    }
    return 0;
}

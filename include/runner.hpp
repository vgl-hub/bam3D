#ifndef runner_hpp
#define runner_hpp

#include <htslib/sam.h>
#include <map>

struct UserInputBam3D : UserInput { // additional input
	bool hist_none        =false;
	bool hist_global      = false;
	bool hist_by_chrom    = false;

	uint8_t decompression_threads = 4;
	uint8_t histogram_mode = hist_none; // histogram mode

};

struct ReadStats {
    uint64_t readN   = 0;
    uint64_t qc_fail = 0;
    uint64_t unmapped = 0;
    uint64_t secondary = 0;
    uint64_t supplementary = 0;
    uint64_t primary = 0;
    uint64_t mapQ0 = 0;

	long double mean_insert=0;
	long double quadratic_mean=0;

	double error_rate=0;
};

struct PairStats {
    uint64_t pairN        = 0;
    uint64_t proper_pairs = 0;
    uint64_t good_pairs   = 0;
	bool good_read1=false; //dovrebbe essere locale
	bool good_read2=false; //dovrebbe essere locale

    uint64_t UMone_sided  = 0;
    uint64_t UMtwo_sided  = 0;
    uint64_t duplicated   = 0;
    uint64_t sameCr       = 0;  // cis

    uint64_t read1 = 0;
    uint64_t read2 = 0;
};

enum class Maptype :uint8_t {N=0, U=1, M=2};

struct QnameStats { 
	Maptype type1;
	Maptype type2;

	uint64_t UU=0;
	uint64_t RU=0;
	uint64_t UR=0;
	uint64_t WW=0;
	uint64_t DD=0;
	uint64_t MU=0;
	uint64_t MM=0;
	uint64_t MR=0;
	uint64_t NM=0;
	uint64_t NU=0;
	uint64_t NR=0;
	uint64_t NN=0;
};

class Runner {
    
    UserInputBam3D userInput;
	ReadStats readStats;
	PairStats pairStats;
	QnameStats qnameStats;
    
public:
    
    void loadInput(UserInputBam3D userInput);
	void bam_sorting();
	long double update_mean_tlen(long double,uint64_t, bam1_t*);
	long double update_quadratic_mean_tlen(long double,uint64_t, bam1_t*);
	double error_rate(uint64_t,uint64_t);
	void qname_group(bam1_t*,std::string&,std::vector<bam1_t*> &);
	void qname_stats(std::vector<bam1_t*> &);
	void flag_inspector(bam1_t*);
	void histo_global_distance(std::unordered_map<uint64_t,uint64_t>&);
	void histo_chrom_distance(std::map<uint32_t,std::unordered_map<uint64_t,uint64_t>>&); 
	void processReads(samFile* , bam_hdr_t* , bam1_t*);
	void output();
	void run();
    
};

#endif /* runner_hpp */


/*  
DOMANDE:
-se il mate di una coppia per un errore non è stato segnato la coppia è invalida ma non sapremo mai quale e non ha senso fare pairN/2(forse con grandi numeri è inutile) 
-ma il pos e il pos del mate non hanno di mezzo la sequenza del frammento stesso? non dovrebbe essere pos dell'ultima base mappata del primo frammento e inizio del mate?

DOMANDE Y:
-i duplicati. sono read doppie? se il file viene prima passato da markdup/picard la flag è segnata giusta? mi fido? se no dovrei fare una passata iniziale (fare il lavoro di markdup)
-WW e R ,DD,parametri di calcolo

MIGLIORAMENTI:
-eliminare variabili globali (metterle come membri pivati della classe runner), sistemare le variabili!!
inserire da termianle le statistiche che si vogliono fare e inizializzare solo le variabili che servono
-file di output
-file di imput(dove si spiega quali imput sono leciti e cosa fanno)
-i pairN se non sono pari, controllo sulle statistiche per non rompere il tool !WARNING


TO DO:
-chiedere info
-statistiche pair type più complesse
-pushare le modifiche nel fork e sulla repo principale
-leggere info sul cluster
-rendere parallelizzabile il codice
-studiare python
*/
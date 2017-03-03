// Wave A-CGH Correction Algorithm implementation (Lepretre et al, Nucleic Acids Res. 2010 Apr;38(7):e94)
// Author : Sylvain Mareschal <mareschal@ovsa.fr>

#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include "R_WACA.h"


// Routine registration (executed at package loading)
void attribute_visible R_init_cghRA(DllInfo *info) {
	// R_WACA argument definition
	R_NativePrimitiveArgType R_WACA_types[] = {
		STRSXP, STRSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
	};
	
	// C() calls
	R_CMethodDef cMethods[] = {
		{"R_WACA", (DL_FUNC) &R_WACA, 11, R_WACA_types},
		{NULL, NULL, 0}
	};
	
	// Register all entry points
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
	
	// Prevent access to all non-registered entry points on Linux, except "attribute_visible" ones
	// Remember to update "src/cghRA-win.def" file for Windows !
	R_useDynamicSymbols(info, FALSE);
}


// Preprocessing of an alphabet for FASTA sequence reading
typedef struct Alphabet {
	int* codeIsAllowed;    // for all unsigned char, 1 if it is allowed, 0 otherwise
	char* codeReplaceBy;   // 'char' to use as code for all unsigned char (may be itself)
} Alphabet;


// Produces an Alphabet structure from an array of allowed characters
// Assumes lowercase and uppercase to be similar
Alphabet* newAlphabet(char* codes, int codeCount) {
	int c;
	Alphabet* object = (Alphabet*) malloc(sizeof(struct Alphabet));
	object->codeIsAllowed = (int*) calloc(UCHAR_MAX+1, sizeof(int));
	object->codeReplaceBy = (char*) calloc(UCHAR_MAX+1, sizeof(char));
	for(c=0; c<codeCount; c++) {
		object->codeIsAllowed[ (unsigned char) tolower(codes[c]) ] = 1;
		object->codeIsAllowed[ (unsigned char) toupper(codes[c]) ] = 1;
		object->codeReplaceBy[ (unsigned char) tolower(codes[c]) ] = toupper(codes[c]);
		object->codeReplaceBy[ (unsigned char) toupper(codes[c]) ] = toupper(codes[c]);
	}
	return object;
}


// Free an Alphabet pointer
void freeAlphabet(Alphabet** alpha) {
	free((*alpha)->codeIsAllowed);
	(*alpha)->codeIsAllowed = NULL;
	free((*alpha)->codeReplaceBy);
	(*alpha)->codeReplaceBy = NULL;
	free((*alpha));
	(*alpha) = NULL;
}


// Counts DNA letters in a single FASTA sequence
// Assumes first line is header, and ignores not allowed characters
int singleFASTA_length(char* fileName, Alphabet* alpha) {
	int i;
	int sequenceLength = 0;
	char* buffer = (char*) malloc(1024*sizeof(char));
	
	FILE* file = fopen(fileName, "r");
	if(fgets(buffer, 1024, file) != NULL) {
		while(fgets(buffer, 1024, file) != NULL) {
			for(i=0; buffer[i] != '\n'; i++) {
				sequenceLength += alpha->codeIsAllowed[ (unsigned char) buffer[i] ];
			}
		
			// To keep R-GUI awoken
			if(sequenceLength % 1000000 == 0) {
				R_ProcessEvents();
				R_CheckUserInterrupt();
			}
		}
	}
	fclose(file);
	free(buffer);
	
	return sequenceLength;
}


// Reads DNA letters from a single FASTA file (needs the sequence length)
// Assumes first line is header, and ignores not allowed characters
// Returned string begins with a '-' and ends with a '\0' (no strlen() !)
char* singleFASTA_sequence(char* fileName, int length, Alphabet* alpha) {
	int i;
	char* sequence = (char*) malloc((length+2)*sizeof(char));
	char* buffer = (char*) malloc(1024*sizeof(char));
	
	// Coordinates begin at 1
	int s = 1;
	sequence[0] = '-';
	
	FILE* file = fopen(fileName, "r");
	if(fgets(buffer, 1024, file) != NULL) {
		while(fgets(buffer, 1024, file) != NULL) {
			for(i=0; buffer[i] != '\n'; i++) {
				sequence[s] = alpha->codeReplaceBy[ (unsigned char) buffer[i] ];
				s++;
			}
		
			// To keep R-GUI awoken
			if(s % 1000000 == 0) {
				R_ProcessEvents();
				R_CheckUserInterrupt();
			}
		}
	}
	sequence[s] = '\0';
	fclose(file);
	free(buffer);
	
	return sequence;
}


// GC frequency in a window of a large sequence
// "start" and "end" begin at 1, both boundaries are comprised
double freqGC(char* sequence, int start, int end) {
	int i, c;
	char* ptr = &(sequence[start]);                          // 'ptr' starts at 'start' in 'sequence'
	int* counts = (int*) calloc(UCHAR_MAX+1, sizeof(int));   // to store all (unsigned) char counts
	for(i=start; i<=end; i++, ptr++) {                       // 'i' only for ending loop, 'ptr' always points to the considered character
		c = (unsigned char) *ptr;                            // converts the letter code into an index (needs to be positive)
		counts[c]++;
	}
	double GC = counts[(unsigned char)'G'] + counts[(unsigned char)'C'] + counts[(unsigned char)'S'];
	double AT = counts[(unsigned char)'A'] + counts[(unsigned char)'T'] + counts[(unsigned char)'W'];
	free(counts);
	return GC / (GC+AT);
}


// Structure for restriction sites
typedef struct Site {
	int length;       // Total base count in the pattern
	char* sequence;   // Bases of the pattern
	int left;         // Width of the left fragment
} Site;


// Returns the last position in the fragment, crawling forward or backward
// 'to' is never reached, stops before
int crawl(Site* sites, int siteCount, char* chrom, int from, int to, int way) {
	int c, site, s;
	for(c=from; c != to; c += way) {                                            // "c" loops between "from" and "to" by "way"
		for(site=0; site < siteCount; site++) {                                     // "site" loops over restriction sites
			for(s=0; chrom[c+s] == sites[site].sequence[s]; s++) {                      // while "chrom" and "site" match, s loops forward (always)
				if(s == sites[site].length-1) {                                             // reached the end of "site", complete match ...
					if(way == 1) {
						return c + sites[site].left - 1;                                        // ... while crawling forward
					} else {
						return c + sites[site].left;                                            // ... while crawling backward
					}
				}
			}
		}
	}
	return to;                                                                  // No match, 'to' reached
}


// wGCprobe, wGC150 and wGC500 computer
// "probeStarts" and "probeEnds" begin at 1, both boundaries are comprised
void wGCxxx(double* target, char* chrom, int chromLength, int window, int probeCount, int* probeStarts, int* probeEnds) {
	int probe, start, end;
	for(probe=0; probe<probeCount; probe++) {
		// Boundaries
		start = probeStarts[probe] - window;
		end = probeEnds[probe] + window;
		if(start < 1) {
			start = 1;
		}
		if(end > chromLength) {
			end = chromLength;
		}
		
		// Computation
		target[probe] = freqGC(chrom, start, end);
		
		// To keep R-GUI awoken
		R_ProcessEvents();
		R_CheckUserInterrupt();
	}
}


// wGCfrag and wFragSize computer
// "probeStarts" and "probeEnds" begin at 1, both boundaries are comprised
// Fragmentation site searching begins at the middle of the probe
void frag(double* wGCfrag, int* wFragSize, char** rawSites, int siteCount, char* chrom, int chromLength, int probeCount, int* probeStarts, int* probeEnds) {
	int site, i, s, probe, left, right, start;
	
	// Fragmentation sites
	Site* sites = (Site*) malloc(siteCount*sizeof(struct Site));
	for(site=0; site<siteCount; site++) {
		// Sequence length
		sites[site].length = 0;
		for(i=0; rawSites[site][i] != '\0'; i++) {
			switch(rawSites[site][i]) {
				case 'A': case 'a': case 'C': case 'c': case 'G': case 'g': case 'T': case 't': sites[site].length++;
			}
		}
		
		// Sequence copy
		sites[site].sequence = (char*) malloc(sites[site].length * sizeof(char));
		sites[site].left = 0;
		s = 0;
		for(i=0; rawSites[site][i] != '\0'; i++) {
			switch(rawSites[site][i]) {
				case 'A': case 'a': sites[site].sequence[s] = 'A'; s++; break;
				case 'C': case 'c': sites[site].sequence[s] = 'C'; s++; break;
				case 'G': case 'g': sites[site].sequence[s] = 'G'; s++; break;
				case 'T': case 't': sites[site].sequence[s] = 'T'; s++; break;
				case '|':           sites[site].left = s; break;
			}
		}
	}
	
	for(probe=0; probe<probeCount; probe++) {
		// Fragment boundaries
		start = (probeEnds[probe] + probeStarts[probe]) / 2;
		right = crawl(sites, 2, chrom, start+1, chromLength+1, 1);
		left = crawl(sites, 2, chrom, start, 0, -1);
		
		// Computations
		wFragSize[probe] = right - left + 1;
		wGCfrag[probe] = freqGC(chrom, left, right);
		
		// To keep R-GUI awoken
		R_ProcessEvents();
		R_CheckUserInterrupt();
	}
	
	// Cleaning
	for(site=0; site<siteCount; site++) {
		free(sites[site].sequence);
	}
	free(sites);
}


// WACA bias computation, R interface
// IUPAC DNA codes are handled properly
void R_WACA(
	char** chromFile,   // Single chromosome at a time
	char** rawSites,    // Multiple fragmentation sites, e.g. "AC|GT"
	int* siteCount,     // Length of "sites"
	int* probeCount,
	int* probeStarts,   // Begins at 1, both boundaries are comprised
	int* probeEnds,     // Begins at 1, both boundaries are comprised
	double* wGC150,
	double* wGC500,
	double* wGCprobe,
	double* wGCfrag,
	int* wFragSize
	) {
	// Chromosome importation
	Alphabet* DNA = newAlphabet("ACGTNRYMSKWBDHV", 15);
	int chromLength = singleFASTA_length(chromFile[0], DNA);
	char* chrom = singleFASTA_sequence(chromFile[0], chromLength, DNA);
	freeAlphabet(&DNA);

	// Bias computations
	wGCxxx(wGCprobe, chrom, chromLength, 0, probeCount[0], probeStarts, probeEnds);
	wGCxxx(wGC150, chrom, chromLength, 150000, probeCount[0], probeStarts, probeEnds);
	wGCxxx(wGC500, chrom, chromLength, 500000, probeCount[0], probeStarts, probeEnds);
	frag(wGCfrag, wFragSize, rawSites, siteCount[0], chrom, chromLength, probeCount[0], probeStarts, probeEnds);
	
	// Cleaning
	free(chrom);
}

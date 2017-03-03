// Wave A-CGH Correction Algorithm implementation (Lepretre et al, Nucleic Acids Res. 2010 Apr;38(7):e94)
// Author : Sylvain Mareschal <mareschal@ovsa.fr>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/*
	Compute the 5 bias values required by WACA for all provided CGH probes
*/

void R_WACA(
	char** chromFile,   // Single chromosome at a time
	char** sites,       // Multiple fragmentation sites, e.g. "AC|GT"
	int* siteCount,     // Length of "sites"
	int* probeCount,
	int* probeStarts,   // 0 based, half open, only for "chromFile"
	int* probeEnds,     // 0 based, half open, only for "chromFile"
	double* wGC150,
	double* wGC500,
	double* wGCprobe,
	double* wGCfrag,
	int* wFragSize
);


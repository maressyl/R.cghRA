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

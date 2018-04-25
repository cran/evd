#include "header.h"
#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[] = {
    {"nlgpd",  (DL_FUNC) &nlgpd,  6},
    {"nlpp", (DL_FUNC) &nlpp, 8},
    {"clusters",  (DL_FUNC) &clusters,  6},
	
	{"ccop",  (DL_FUNC) &ccop,  11},
	
	{"rbvlog_shi",  (DL_FUNC) &rbvlog_shi,  3},
    {"rbvalog_shi", (DL_FUNC) &rbvalog_shi, 4},
    {"rmvlog_tawn",  (DL_FUNC) &rmvlog_tawn,  4},
	{"rmvalog_tawn",  (DL_FUNC) &rmvalog_tawn,  6},
    {"rbvlog", (DL_FUNC) &rbvlog, 3},
    {"rbvalog",  (DL_FUNC) &rbvalog,  4},
	{"rbvhr",  (DL_FUNC) &rbvhr,  3},
    {"rbvneglog", (DL_FUNC) &rbvneglog, 3},
    {"rbvaneglog",  (DL_FUNC) &rbvaneglog,  4},
	{"rbvbilog",  (DL_FUNC) &rbvbilog,  4},
    {"rbvnegbilog", (DL_FUNC) &rbvnegbilog, 4},
    {"rbvct",  (DL_FUNC) &rbvct,  4},
	{"rbvamix",  (DL_FUNC) &rbvamix,  4},
	
	{"nlgev", (DL_FUNC) &nlgev, 6},
	{"nslmvalog",  (DL_FUNC) &nslmvalog,  13},
	
	{"nlbvlog", (DL_FUNC) &nlbvlog, 13},
    {"nlbvalog",  (DL_FUNC) &nlbvalog,  15},
	{"nlbvhr",  (DL_FUNC) &nlbvhr,  13},
    {"nlbvneglog", (DL_FUNC) &nlbvneglog, 13},
    {"nlbvaneglog",  (DL_FUNC) &nlbvaneglog,  15},
	{"nlbvbilog",  (DL_FUNC) &nlbvbilog,  14},
    {"nlbvnegbilog", (DL_FUNC) &nlbvnegbilog, 14},
    {"nlbvct",  (DL_FUNC) &nlbvct,  14},
	{"nlbvamix",  (DL_FUNC) &nlbvamix,  14},
	
	{"nllbvclog", (DL_FUNC) &nllbvclog, 12},
    {"nllbvcalog",  (DL_FUNC) &nllbvcalog,  14},
	{"nllbvchr",  (DL_FUNC) &nllbvchr,  12},
    {"nllbvcneglog", (DL_FUNC) &nllbvcneglog, 12},
    {"nllbvcaneglog",  (DL_FUNC) &nllbvcaneglog,  14},
	{"nllbvcbilog",  (DL_FUNC) &nllbvcbilog,  13},
    {"nllbvcnegbilog", (DL_FUNC) &nllbvcnegbilog, 13},
    {"nllbvcct",  (DL_FUNC) &nllbvcct,  13},
	{"nllbvcamix",  (DL_FUNC) &nllbvcamix,  13},
	
	{"nllbvplog", (DL_FUNC) &nllbvplog, 13},
	{"nllbvphr",  (DL_FUNC) &nllbvphr,  13},
    {"nllbvpneglog", (DL_FUNC) &nllbvpneglog, 13},
	{"nllbvpbilog",  (DL_FUNC) &nllbvpbilog,  14},
    {"nllbvpnegbilog", (DL_FUNC) &nllbvpnegbilog, 14},
    {"nllbvpct",  (DL_FUNC) &nllbvpct,  14},
	
    {NULL, NULL, 0}
};	   

void R_init_evd(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}



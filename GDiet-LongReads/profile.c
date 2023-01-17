#include "profile.h"
#include <stdint.h>
#include <stdio.h>

#ifdef PROFILE
atomic_uint_least64_t pf_indexing           = 0;
atomic_uint_least64_t pf_pattern_alignment  = 0;
atomic_uint_least64_t pf_seeding            = 0;
atomic_uint_least64_t pf_voting             = 0;
atomic_uint_least64_t pf_sequence_alignment = 0;

void print_profile() {
	fprintf(stderr, "[PROFILING] indexing time: %lu ns\n", pf_indexing);
	fprintf(stderr, "[PROFILING] pattern alignment time: %lu ns\n", pf_pattern_alignment);
	fprintf(stderr, "[PROFILING] seeding time: %lu ns\n", pf_seeding);
	fprintf(stderr, "[PROFILING] voting time: %lu ns\n", pf_voting);
	fprintf(stderr, "[PROFILING] sequence alignment time: %lu ns\n", pf_sequence_alignment);
}

#else
void print_profile() {}
#endif

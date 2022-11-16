#ifndef PROFILE_H
#define PROFILE_H

void print_profile();

#ifdef PROFILE
#include <stdatomic.h>
#include <time.h>

extern atomic_uint_least64_t pf_indexing;
extern atomic_uint_least64_t pf_pattern_alignment;
extern atomic_uint_least64_t pf_seeding;
extern atomic_uint_least64_t pf_voting;
extern atomic_uint_least64_t pf_sequence_alignment;

#define PROF_INIT struct timespec start, end
#define PROF_START clock_gettime(CLOCK_MONOTONIC, &start)
#define PROF_END(prof)                                                                                                 \
	{                                                                                                              \
		clock_gettime(CLOCK_MONOTONIC, &end);                                                                  \
		atomic_fetch_add(&prof, (end.tv_sec - start.tv_sec) * 1000000000 + (end.tv_nsec - start.tv_nsec));     \
	}

#else
#define PROF_INIT
#define PROF_START
#define PROF_END(prof)
#endif

#endif

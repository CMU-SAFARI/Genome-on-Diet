#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"
#include <math.h>

unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

static inline unsigned get_real_location(const unsigned location, const unsigned pattern_len, const int ones,
                                         const int *const ones_loc, const unsigned shift) {
	return (location / ones) * pattern_len + ones_loc[location % ones] + shift;
}

static inline uint64_t hash64(uint64_t key, uint64_t mask) {
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x) { q->a[((q->count++) + q->front) & 0x1f] = x; }

static inline int tq_shift(tiny_queue_t *q) {
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */

#ifdef __AVX512DQ__
#include <immintrin.h>

static uint8_t first_pos_array[256] = {
    0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2,
    0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0,
    1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1,
    0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0,
    2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3,
    0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
    1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0};

static uint8_t last_pos_array[256] = {
    0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};

static inline uint64_t extract_u64(__m512i vec, uint8_t pos) {
	switch (pos) {
		case 0:
			return _mm_extract_epi64(_mm512_extracti64x2_epi64(vec, 0), 0);
		case 1:
			return _mm_extract_epi64(_mm512_extracti64x2_epi64(vec, 0), 1);
		case 2:
			return _mm_extract_epi64(_mm512_extracti64x2_epi64(vec, 1), 0);
		case 3:
			return _mm_extract_epi64(_mm512_extracti64x2_epi64(vec, 1), 1);
		case 4:
			return _mm_extract_epi64(_mm512_extracti64x2_epi64(vec, 2), 0);
		case 5:
			return _mm_extract_epi64(_mm512_extracti64x2_epi64(vec, 2), 1);
		case 6:
			return _mm_extract_epi64(_mm512_extracti64x2_epi64(vec, 3), 0);
		case 7:
			return _mm_extract_epi64(_mm512_extracti64x2_epi64(vec, 3), 1);
	}
	return 0;
}

static inline void pack(const char *str, uint64_t *c, uint64_t *c_rev, uint8_t *n_mask, const unsigned location,
                        const unsigned shift, const unsigned pattern_len, const int k, const int ones,
                        const int *const ones_loc, const unsigned len) {
	*c      = 0;
	*c_rev  = 0;
	*n_mask = 0;

	for (unsigned j = 0; j < 8; j++) {
		unsigned real_location = get_real_location(j + location, pattern_len, ones, ones_loc, shift);
		if (real_location >= len) {
			return;
		}
		int _c = seq_nt4_table[(uint8_t)str[real_location]];
		if (_c == 4) {
			*n_mask |= (1 << j);
		} else {
			*c |= _c << (14 - 2 * j);
			*c_rev |= (3ULL ^ _c) << (2 * (k - 8 + j));
		}
	}
}

static inline __m512i hash64_avx(__m512i key, uint64_t MASK) {
	const __m512i mask = _mm512_set1_epi64(MASK);

	key = _mm512_and_epi64(_mm512_add_epi64(_mm512_andnot_epi64(key, mask), _mm512_slli_epi64(key, 21)), mask);
	key = _mm512_xor_epi64(key, _mm512_srli_epi64(key, 24));
	key = _mm512_and_epi64(
	    _mm512_add_epi64(_mm512_add_epi64(key, _mm512_slli_epi64(key, 3)), _mm512_slli_epi64(key, 8)), mask);
	key = _mm512_xor_epi64(key, _mm512_srli_epi64(key, 14));
	key = _mm512_and_epi64(
	    _mm512_add_epi64(_mm512_add_epi64(key, _mm512_slli_epi64(key, 2)), _mm512_slli_epi64(key, 4)), mask);
	key = _mm512_xor_epi64(key, _mm512_srli_epi64(key, 28));
	key = _mm512_and_epi64(_mm512_add_epi64(key, _mm512_slli_epi64(key, 31)), mask);
	return key;
}

typedef struct {
	uint64_t hash; // size: 2*K
	uint32_t loc;
	uint8_t str;
} seed_t;

void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p, const char *Z,
               int W) {
	uint64_t mask = (1ULL << 2 * k) - 1;

	assert(len > 0 && (w > 0 && w < 256) &&
	       (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice

	// char pattern[10] = "111000111"; // determine the pattern
	// int pattern_len = 9;

	const char *pattern = Z; // determine the pattern
	int pattern_len     = W; // determine length of pattern

	int ones = 0;     // amount of ones in the pattern
	int ones_loc[40]; // the maximum number of ones in the pattern is 40 for now

	for (int g = 0; g < pattern_len; ++g) {
		if (pattern[g] == '1') {
			ones_loc[ones] = g;
			++ones; // count the ones in the pattern
		}
	}

	// Len after appliying the pattern
	unsigned diet_len        = (len / pattern_len) * ones;
	unsigned diet_len_offset = len % pattern_len;
	for (unsigned i = 0; i < diet_len_offset; i++) {
		if (pattern[i] == '1') {
			diet_len += 1;
		}
	}

	/*
	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
	        char *new_genome = 0;
	        new_genome  = (char *)malloc(sizeof(char) * diet_len); // allocate memory for the shortened
	sequence int current = 0; // gives the position of the base in the new_genome for (int f = 0; f < len;
	++f) {
	                // printf("%d_%c_%d_%d\n",pattern[((f-shift)%pattern_len)], str[f],f,
	                // ((f-shift)%pattern_len));
	                if (pattern[((f - shift) % pattern_len)] == '1') {
	                        new_genome[current] = str[f];
	                        ++current;
	                }
	        }
	        kv_resize(mm128_t, km, *p, p->n + diet_len / w);
	        // Print the new Genome to check
	        printf("\nNew Read: ");
	        for (int pp = 0; pp < diet_len; ++pp) {
	                printf("%c", new_genome[pp]);
	                // if ((pp+1) % k == 0) printf("\t");
	        }
	        printf("\n");
	        free(new_genome);
	}
	*/

	const __m512i kmer_shift = _mm512_set_epi64(16, 14, 12, 10, 8, 6, 4, 2);
	const __m512i c_shift    = _mm512_set_epi64(0, 2, 4, 6, 8, 10, 12, 14);
	const __m512i mask_vec   = _mm512_set1_epi64(mask);
	const unsigned buf_size  = (w / 8) + ((w % 8) != 0);

	for (uint32_t extracted_len = 0; extracted_len < diet_len;) {
		__m512i kmer[2]  = {_mm512_setzero_si512()};
		seed_t minimizer = {.hash = UINT64_MAX, .loc = k - 2, .str = 0};
		__m512i hash_buf[32];
		__mmask8 mask_buf[32];
		__mmask8 strand_buf[32];
		unsigned buf_pos = 0;

		uint32_t len_early_exit = diet_len - extracted_len;
		// const uint8_t *seq_off     = seq + extracted_len;
		unsigned early_exit_offset = 0;

		int first_window    = 1;
		__m512i minimum_vec = _mm512_set1_epi64(minimizer.hash);
		for (uint32_t i = 0; i < len_early_exit; i += 8) {

			uint64_t seq_packed, seq_packed_rev;
			uint8_t n_mask;
			pack(str, &seq_packed, &seq_packed_rev, &n_mask, extracted_len + i, 0, pattern_len, k, ones,
			     ones_loc, len);

			// Compute the + and - kmers
			__m512i c = _mm512_set1_epi64(seq_packed);
			c         = _mm512_srlv_epi64(c, c_shift);
			kmer[0]   = _mm512_set1_epi64(_mm_extract_epi64(_mm512_extracti64x2_epi64(kmer[0], 3), 1));
			kmer[0]   = _mm512_sllv_epi64(kmer[0], kmer_shift);
			kmer[0]   = _mm512_or_epi64(kmer[0], c);
			kmer[0]   = _mm512_and_epi64(kmer[0], mask_vec);

			__m512i c_rev = _mm512_set1_epi64(seq_packed_rev);
			c_rev         = _mm512_sllv_epi64(c_rev, c_shift);
			kmer[1]       = _mm512_set1_epi64(_mm_extract_epi64(_mm512_extracti64x2_epi64(kmer[1], 3), 1));
			kmer[1]       = _mm512_srlv_epi64(kmer[1], kmer_shift);
			kmer[1]       = _mm512_or_epi64(kmer[1], c_rev);
			kmer[1]       = _mm512_and_epi64(kmer[1], mask_vec);

			// If there are ambiguous bases
			if (n_mask != 0) {
				// If there are multiple `N`s
				unsigned first_n_pos = first_pos_array[n_mask];
				unsigned last_n_pos  = last_pos_array[n_mask];
				if (i + first_n_pos <= k - 1) {
					len_early_exit = i + last_n_pos + 1;
					break;
				}
				len_early_exit = (i + first_n_pos > len_early_exit) ? len_early_exit : i + first_n_pos;
				early_exit_offset = last_n_pos - first_n_pos + 1;
			}

			// Skip before the first k-mers
			if (i + 8 < k) {
				continue;
			}

			// Skip symetric k-mers
			__mmask8 skip_mask = _mm512_cmpneq_epu64_mask(kmer[0], kmer[1]);

			// Skip before the first k-mers
			if (i < k - 1) {
				skip_mask &= UINT8_MAX << (k - 1 - i);
			}

			// Skip after the last k-mers
			if (i + 8 > len_early_exit) {
				skip_mask &= UINT8_MAX >> (i + 8 - len_early_exit);
			}

			// Compute the hashes
			__mmask8 strand_mask    = _mm512_cmpgt_epu64_mask(kmer[0], kmer[1]);
			__m512i min_kmer_strand = _mm512_mask_blend_epi64(strand_mask, kmer[0], kmer[1]);
			__m512i hash            = hash64_avx(min_kmer_strand, mask);

			/*
			uint64_t debug_vec[8];
			_mm512_store_epi64(debug_vec, hash);
			for (unsigned k = 0; k < 8; k++) {
			        printf("0x%010lx\t%lu\n", debug_vec[k], i + k + offset);
			}
			printf("n_mask: %x\n", n_mask);
			printf("minimum: %lx, loc: %lu\n", minimizer.hash, minimizer.loc + offset);
			puts("");
			*/

			__mmask8 in_window_mask = skip_mask;

			if (i + 7 - minimizer.loc >= w) {
				in_window_mask &= UINT8_MAX >> (i + 7 - minimizer.loc - (w - 1));
			}
			__mmask8 cmp_mask = _mm512_mask_cmple_epu64_mask(in_window_mask, hash, minimum_vec);
			while (cmp_mask != 0) {
				unsigned minimizer_pos = first_pos_array[cmp_mask];
				// Only store the minimizer after the first window
				if (i + minimizer_pos >= k + w - 1) {
					// If previous minimizer in the first window, search for all the seeds
					// with the same hash in the first window
					if (first_window) {
						first_window = 0;
						if (minimizer.hash != UINT64_MAX) {

							// Offset to compute the location of the k-mers
							uint32_t base_pos = (k - 1) & (UINT32_MAX << 3);
							// Buffer index of the last minimizer of the window
							unsigned minimizer_buf = (minimizer.loc - base_pos) >> 3;
							for (int j = 0; j < minimizer_buf; j++) {
								__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
								    mask_buf[j], minimum_vec, hash_buf[j]);
								while (cmp_mask != 0) {
									unsigned minimizer_pos =
									    first_pos_array[cmp_mask];
									__mmask8 select_mask = 1 << minimizer_pos;
									int strand =
									    ((select_mask & strand_buf[j]) != 0);
									uint32_t real_loc = get_real_location(
									    base_pos + 8 * j + minimizer_pos +
									        extracted_len,
									    pattern_len, ones, ones_loc, 0);
									mm128_t min_tmp =
									    (mm128_t){.x = minimizer.hash << 8 | k,
									              .y = (uint64_t)rid << 32 |
									                   real_loc << 1 | strand};
									kv_push(mm128_t, km, *p, min_tmp);
									mask_buf[j] &= ~select_mask;
									cmp_mask = _mm512_mask_cmpeq_epu64_mask(
									    mask_buf[j], minimum_vec, hash_buf[j]);
								}
							}

							__m512i last_hash;
							uint8_t *last_mask;
							uint8_t last_strand_mask;
							if (minimizer.loc < i) {
								last_hash        = hash_buf[minimizer_buf];
								last_mask        = &mask_buf[minimizer_buf];
								last_strand_mask = strand_buf[minimizer_buf];
							} else {
								last_hash        = hash;
								last_mask        = &skip_mask;
								last_strand_mask = strand_mask;
							}
							// Mask the k-mers after the current minimizer
							unsigned minimizer_pos = (minimizer.loc - base_pos) & 0x7;
							__mmask8 in_window_mask =
							    (UINT8_MAX >> (8 - minimizer_pos)) & *last_mask;
							__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    in_window_mask, minimum_vec, last_hash);
							while (cmp_mask != 0) {
								unsigned minimizer_pos = first_pos_array[cmp_mask];
								__mmask8 select_mask   = 1 << minimizer_pos;
								int strand = ((select_mask & last_strand_mask) != 0);
								uint32_t real_loc =
								    get_real_location(base_pos + 8 * minimizer_buf +
								                          minimizer_pos + extracted_len,
								                      pattern_len, ones, ones_loc, 0);
								mm128_t min_tmp = (mm128_t){
								    .x = minimizer.hash << 8 | k,
								    .y = (uint64_t)rid << 32 | real_loc << 1 | strand};
								kv_push(mm128_t, km, *p, min_tmp);
								in_window_mask &= ~select_mask;
								*last_mask &= ~select_mask;
								cmp_mask = _mm512_mask_cmpeq_epu64_mask(
								    in_window_mask, minimum_vec, last_hash);
							}

							// Mask the current minimizer in skip_mask
							if (minimizer.loc >= i) {
								skip_mask &= ~(1 << minimizer_pos);
							}
						}
					}
					if (minimizer.hash != UINT64_MAX) {
						uint32_t real_loc = get_real_location(minimizer.loc + extracted_len,
						                                      pattern_len, ones, ones_loc, 0);

						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};

						kv_push(mm128_t, km, *p, min_tmp);
					}
					__mmask8 select_mask = 1 << minimizer_pos;
					minimizer            = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                                .loc  = i + minimizer_pos,
					                                .str  = ((select_mask & strand_mask) != 0)};
					skip_mask &= ~select_mask;
					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(skip_mask, hash, minimum_vec);
				} else {
					__mmask8 select_mask = 1 << minimizer_pos;
					minimizer            = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                                .loc  = i + minimizer_pos,
					                                .str  = ((select_mask & strand_mask) != 0)};

					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(
					       skip_mask & (UINT8_MAX << (minimizer_pos + 1)), hash, minimum_vec);
				}
			}

			// If the last fetched k-mer is not in the first window
			if (first_window && i + 7 >= k + w - 1) {
				first_window = 0;
				if (minimizer.hash != UINT64_MAX) {
					uint32_t base_pos = ((k - 1) >> 3) << 3;
					// Buffer index of the last minimizer of the window
					unsigned minimizer_buf = (minimizer.loc - base_pos) >> 3;
					for (int j = 0; j < minimizer_buf; j++) {
						__mmask8 cmp_mask =
						    _mm512_mask_cmpeq_epu64_mask(mask_buf[j], minimum_vec, hash_buf[j]);
						while (cmp_mask != 0) {
							unsigned minimizer_pos = first_pos_array[cmp_mask];
							__mmask8 select_mask   = 1 << minimizer_pos;
							int strand             = ((select_mask & strand_buf[j]) != 0);
							uint32_t real_loc      = get_real_location(
							         base_pos + 8 * j + minimizer_pos + extracted_len,
							         pattern_len, ones, ones_loc, 0);
							mm128_t min_tmp = (mm128_t){.x = minimizer.hash << 8 | k,
							                            .y = (uint64_t)rid << 32 |
							                                 real_loc << 1 | strand};
							kv_push(mm128_t, km, *p, min_tmp);
							mask_buf[j] &= ~select_mask;
							cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    mask_buf[j], minimum_vec, hash_buf[j]);
						}
					}

					__m512i last_hash;
					uint8_t last_mask;
					uint8_t last_strand_mask;
					if (minimizer.loc < i) {
						last_hash        = hash_buf[minimizer_buf];
						last_mask        = mask_buf[minimizer_buf];
						last_strand_mask = strand_buf[minimizer_buf];
					} else {
						last_hash        = hash;
						last_mask        = skip_mask;
						last_strand_mask = strand_mask;
					}

					// Mask the k-mers after the current minimizer
					unsigned minimizer_pos  = (minimizer.loc - base_pos) & 0x7;
					__mmask8 in_window_mask = (UINT8_MAX >> (8 - minimizer_pos)) & last_mask;
					__mmask8 cmp_mask =
					    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, last_hash);
					while (cmp_mask != 0) {
						unsigned minimizer_pos = first_pos_array[cmp_mask];
						__mmask8 select_mask   = 1 << minimizer_pos;
						int strand             = ((select_mask & last_strand_mask) != 0);

						uint32_t real_loc = get_real_location(base_pos + 8 * minimizer_buf +
						                                          minimizer_pos + extracted_len,
						                                      pattern_len, ones, ones_loc, 0);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | strand};
						kv_push(mm128_t, km, *p, min_tmp);
						in_window_mask &= ~select_mask;
						cmp_mask = _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec,
						                                        last_hash);
					}
				}
			}

			// If the current window doens't cover the last fetched k-mer
			while (minimizer.loc + w < len_early_exit && i + 7 - minimizer.loc >= w) {
				// Find the Minimimum in the next window

				// Mask the k-mers after the next window
				in_window_mask = skip_mask & (UINT8_MAX >> (i + 7 - minimizer.loc - w));

				// Mask the k-mers before the current minimizer
				unsigned buf_index = buf_pos;
				mask_buf[buf_index] &= UINT8_MAX << ((buf_size << 3) + 1 + minimizer.loc - i);
				if (minimizer.loc >= i - ((buf_size - 1) << 3)) {
					mask_buf[(buf_index + 1) % buf_size] &=
					    UINT8_MAX << (((buf_size - 1) << 3) + 1 + minimizer.loc - i);
				}

				uint64_t minimum = UINT64_MAX;
				for (unsigned j = 0; j < buf_size; j++) {
					uint64_t min_tmp =
					    _mm512_mask_reduce_min_epu64(mask_buf[buf_index], hash_buf[buf_index]);
					if (min_tmp < minimum) {
						minimum = min_tmp;
					}
					buf_index = (buf_index + 1) % buf_size;
				}
				uint64_t min_tmp = _mm512_mask_reduce_min_epu64(in_window_mask, hash);
				if (min_tmp < minimum) {
					minimum = min_tmp;
				}

				minimum_vec = _mm512_set1_epi64(minimum);

				// Get its location(s)
				if (minimum == UINT64_MAX) {
					if (minimizer.hash != UINT64_MAX) {
						uint32_t real_loc = get_real_location(minimizer.loc + extracted_len,
						                                      pattern_len, ones, ones_loc, 0);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};

						kv_push(mm128_t, km, *p, min_tmp);
					}
					minimizer.hash = UINT64_MAX;
					minimizer.loc  = minimizer.loc + w;
				} else {
					unsigned loc = i - 8 * buf_size;
					for (unsigned j = 0; j < buf_size; j++) {
						__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
						    mask_buf[buf_index], minimum_vec, hash_buf[buf_index]);
						while (cmp_mask != 0) {
							unsigned minimizer_pos = first_pos_array[cmp_mask];
							__mmask8 select_mask   = 1 << minimizer_pos;
							if (minimizer.hash != UINT64_MAX) {
								uint32_t real_loc =
								    get_real_location(minimizer.loc + extracted_len,
								                      pattern_len, ones, ones_loc, 0);
								mm128_t min_tmp =
								    (mm128_t){.x = minimizer.hash << 8 | k,
								              .y = (uint64_t)rid << 32 | real_loc << 1 |
								                   minimizer.str};
								kv_push(mm128_t, km, *p, min_tmp);
							}
							minimizer = (seed_t){
							    .hash = minimum,
							    .loc  = loc + minimizer_pos,
							    .str  = ((strand_buf[buf_index] & select_mask) != 0)};
							mask_buf[buf_index] &= ~select_mask;
							cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    mask_buf[buf_index], minimum_vec, hash_buf[buf_index]);
						}
						loc += 8;
						buf_index = (buf_index + 1) % buf_size;
					}
					__mmask8 cmp_mask =
					    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, hash);
					while (cmp_mask != 0) {
						unsigned minimizer_pos = first_pos_array[cmp_mask];
						__mmask8 select_mask   = 1 << minimizer_pos;
						if (minimizer.hash != UINT64_MAX) {
							uint32_t real_loc =
							    get_real_location(minimizer.loc + extracted_len,
							                      pattern_len, ones, ones_loc, 0);
							mm128_t min_tmp = (mm128_t){.x = minimizer.hash << 8 | k,
							                            .y = (uint64_t)rid << 32 |
							                                 real_loc << 1 | minimizer.str};
							kv_push(mm128_t, km, *p, min_tmp);
						}
						minimizer = (seed_t){.hash = minimum,
						                     .loc  = i + minimizer_pos,
						                     .str  = ((strand_mask & select_mask) != 0)};
						skip_mask &= ~select_mask;
						in_window_mask &= ~select_mask;
						cmp_mask =
						    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, hash);
					}
				}

				// Check if there are minimizers inside of the new window
				in_window_mask = skip_mask;
				if (i + 7 - minimizer.loc >= w) {
					in_window_mask &= UINT8_MAX >> (i + 7 - minimizer.loc - (w - 1));
				}
				__mmask8 cmp_mask = _mm512_mask_cmple_epu64_mask(in_window_mask, hash, minimum_vec);
				while (cmp_mask != 0) {
					unsigned minimizer_pos = first_pos_array[cmp_mask];
					__mmask8 select_mask   = 1 << minimizer_pos;
					if (minimizer.hash != UINT64_MAX) {
						uint32_t real_loc = get_real_location(minimizer.loc + extracted_len,
						                                      pattern_len, ones, ones_loc, 0);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};
						kv_push(mm128_t, km, *p, min_tmp);
					}
					minimizer = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                     .loc  = i + minimizer_pos,
					                     .str  = ((strand_mask & select_mask) != 0)};
					skip_mask &= ~select_mask;
					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(skip_mask, hash, minimum_vec);
				}
			}
			hash_buf[buf_pos]   = hash;
			mask_buf[buf_pos]   = skip_mask;
			strand_buf[buf_pos] = strand_mask;
			buf_pos             = (buf_pos + 1) % buf_size;
		}
		if (len_early_exit >= w + k - 1 && minimizer.hash != UINT64_MAX) {
			uint32_t real_loc =
			    get_real_location(minimizer.loc + extracted_len, pattern_len, ones, ones_loc, 0);
			mm128_t min_tmp = (mm128_t){.x = minimizer.hash << 8 | k,
			                            .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};
			kv_push(mm128_t, km, *p, min_tmp);
		}
		extracted_len += len_early_exit + early_exit_offset;
	}
}

static inline unsigned mm_sketch2_sub(void *km, const char *str, const unsigned diet_len, const int w, const int k,
                                      const uint32_t rid, mm128_v *p, const unsigned shift, const unsigned max_seeds,
                                      const int ones, const int *const ones_loc, const unsigned pattern_len,
                                      const uint64_t shift1, const uint64_t mask, const unsigned len) {
	unsigned min_counter     = 0;
	const __m512i kmer_shift = _mm512_set_epi64(16, 14, 12, 10, 8, 6, 4, 2);
	const __m512i c_shift    = _mm512_set_epi64(0, 2, 4, 6, 8, 10, 12, 14);
	const __m512i mask_vec   = _mm512_set1_epi64(mask);
	const unsigned buf_size  = (w / 8) + ((w % 8) != 0);

	for (uint32_t extracted_len = 0; extracted_len < diet_len;) {
		__m512i kmer[2]  = {_mm512_setzero_si512()};
		seed_t minimizer = {.hash = UINT64_MAX, .loc = k - 2, .str = 0};
		__m512i hash_buf[32];
		__mmask8 mask_buf[32];
		__mmask8 strand_buf[32];
		unsigned buf_pos = 0;

		uint32_t len_early_exit = diet_len - extracted_len;
		// printf("len_early_exit: %u\n", len_early_exit);
		// const uint8_t *seq_off     = seq + extracted_len;
		unsigned early_exit_offset = 0;

		int first_window    = 1;
		__m512i minimum_vec = _mm512_set1_epi64(minimizer.hash);
		for (uint32_t i = 0; i < len_early_exit; i += 8) {

			uint64_t seq_packed, seq_packed_rev;
			uint8_t n_mask;
			pack(str, &seq_packed, &seq_packed_rev, &n_mask, extracted_len + i, shift, pattern_len, k, ones,
			     ones_loc, len);

			// Compute the + and - kmers
			__m512i c = _mm512_set1_epi64(seq_packed);
			c         = _mm512_srlv_epi64(c, c_shift);
			kmer[0]   = _mm512_set1_epi64(_mm_extract_epi64(_mm512_extracti64x2_epi64(kmer[0], 3), 1));
			kmer[0]   = _mm512_sllv_epi64(kmer[0], kmer_shift);
			kmer[0]   = _mm512_or_epi64(kmer[0], c);
			kmer[0]   = _mm512_and_epi64(kmer[0], mask_vec);

			__m512i c_rev = _mm512_set1_epi64(seq_packed_rev);
			c_rev         = _mm512_sllv_epi64(c_rev, c_shift);
			kmer[1]       = _mm512_set1_epi64(_mm_extract_epi64(_mm512_extracti64x2_epi64(kmer[1], 3), 1));
			kmer[1]       = _mm512_srlv_epi64(kmer[1], kmer_shift);
			kmer[1]       = _mm512_or_epi64(kmer[1], c_rev);
			kmer[1]       = _mm512_and_epi64(kmer[1], mask_vec);

			// If there are ambiguous bases
			if (n_mask != 0) {
				// If there are multiple `N`s
				unsigned first_n_pos = first_pos_array[n_mask];
				unsigned last_n_pos  = last_pos_array[n_mask];
				if (i + first_n_pos <= k - 1) {
					len_early_exit = i + last_n_pos + 1;
					break;
				}
				len_early_exit = (i + first_n_pos > len_early_exit) ? len_early_exit : i + first_n_pos;
				early_exit_offset = last_n_pos - first_n_pos + 1;
			}

			// Skip before the first k-mers
			if (i + 8 < k) {
				continue;
			}

			// Skip symetric k-mers
			__mmask8 skip_mask = _mm512_cmpneq_epu64_mask(kmer[0], kmer[1]);

			// Skip before the first k-mers
			if (i < k - 1) {
				skip_mask &= UINT8_MAX << (k - 1 - i);
			}

			// Skip after the last k-mers
			if (i + 8 > len_early_exit) {
				skip_mask &= UINT8_MAX >> (i + 8 - len_early_exit);
				// printf("len_early_exit: %u, i: %u, mask: %x\n", len_early_exit, i,
				// skip_mask);
			}

			// Compute the hashes
			__mmask8 strand_mask    = _mm512_cmpgt_epu64_mask(kmer[0], kmer[1]);
			__m512i min_kmer_strand = _mm512_mask_blend_epi64(strand_mask, kmer[0], kmer[1]);
			__m512i hash            = hash64_avx(min_kmer_strand, mask);

			/*
			uint64_t debug_vec[8];
			_mm512_store_epi64(debug_vec, hash);
			for (unsigned k = 0; k < 8; k++) {
			        printf("0x%010lx\t%u\n", debug_vec[k], i + k);
			}
			printf("n_mask: %x\n", n_mask);
			printf("minimum: %lx, loc: %u\n", minimizer.hash, minimizer.loc);
			puts("");
			*/

			__mmask8 in_window_mask = skip_mask;

			if (i + 7 - minimizer.loc >= w) {
				in_window_mask &= UINT8_MAX >> (i + 7 - minimizer.loc - (w - 1));
			}
			__mmask8 cmp_mask = _mm512_mask_cmple_epu64_mask(in_window_mask, hash, minimum_vec);
			while (cmp_mask != 0) {
				unsigned minimizer_pos = first_pos_array[cmp_mask];
				// Only store the minimizer after the first window
				// printf("m_pos: %u\n", minimizer_pos);
				if (i + minimizer_pos >= k + w - 1) {
					// If previous minimizer in the first window, search for all the seeds
					// with the same hash in the first window
					if (first_window) {
						// printf("hash: %lx\n", minimizer.hash << 8 | k);
						first_window = 0;
						if (minimizer.hash != UINT64_MAX) {

							// Offset to compute the location of the k-mers
							uint32_t base_pos = (k - 1) & (UINT32_MAX << 3);
							// Buffer index of the last minimizer of the window
							unsigned minimizer_buf = (minimizer.loc - base_pos) >> 3;
							for (int j = 0; j < minimizer_buf; j++) {
								__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
								    mask_buf[j], minimum_vec, hash_buf[j]);
								while (cmp_mask != 0) {
									unsigned minimizer_pos =
									    first_pos_array[cmp_mask];
									__mmask8 select_mask = 1 << minimizer_pos;
									int strand =
									    ((select_mask & strand_buf[j]) != 0);
									uint32_t loc = base_pos + 8 * j +
									               minimizer_pos + extracted_len;
									uint32_t real_loc = get_real_location(
									    loc, pattern_len, ones, ones_loc, shift);
									mm128_t min_tmp =
									    (mm128_t){.x = minimizer.hash << 8 | k,
									              .y = (uint64_t)rid << 32 |
									                   real_loc << 1 | strand};
									kv_push(mm128_t, km, *p, min_tmp);
									min_counter++;
									if (min_counter == max_seeds) {
										return min_counter;
									}
									mask_buf[j] &= ~select_mask;
									cmp_mask = _mm512_mask_cmpeq_epu64_mask(
									    mask_buf[j], minimum_vec, hash_buf[j]);
								}
							}

							__m512i last_hash;
							uint8_t *last_mask;
							uint8_t last_strand_mask;
							if (minimizer.loc < i) {
								last_hash        = hash_buf[minimizer_buf];
								last_mask        = &mask_buf[minimizer_buf];
								last_strand_mask = strand_buf[minimizer_buf];
							} else {
								last_hash        = hash;
								last_mask        = &skip_mask;
								last_strand_mask = strand_mask;
							}
							// Mask the k-mers after the current minimizer
							unsigned minimizer_pos = (minimizer.loc - base_pos) & 0x7;
							__mmask8 in_window_mask =
							    (UINT8_MAX >> (8 - minimizer_pos)) & *last_mask;
							__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    in_window_mask, minimum_vec, last_hash);
							while (cmp_mask != 0) {
								unsigned minimizer_pos = first_pos_array[cmp_mask];
								__mmask8 select_mask   = 1 << minimizer_pos;
								int strand   = ((select_mask & last_strand_mask) != 0);
								uint32_t loc = base_pos + 8 * minimizer_buf +
								               minimizer_pos + extracted_len;
								uint32_t real_loc = get_real_location(
								    loc, pattern_len, ones, ones_loc, shift);
								mm128_t min_tmp = (mm128_t){
								    .x = minimizer.hash << 8 | k,
								    .y = (uint64_t)rid << 32 | real_loc << 1 | strand};
								kv_push(mm128_t, km, *p, min_tmp);
								min_counter++;
								if (min_counter == max_seeds) {
									return min_counter;
								}
								in_window_mask &= ~select_mask;
								*last_mask &= ~select_mask;
								cmp_mask = _mm512_mask_cmpeq_epu64_mask(
								    in_window_mask, minimum_vec, last_hash);
							}

							// Mask the current minimizer in skip_mask
							if (minimizer.loc >= i) {
								skip_mask &= ~(1 << minimizer_pos);
							}
						}
					}
					if (minimizer.hash != UINT64_MAX) {
						uint32_t loc = minimizer.loc + extracted_len;
						uint32_t real_loc =
						    get_real_location(loc, pattern_len, ones, ones_loc, shift);

						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};

						kv_push(mm128_t, km, *p, min_tmp);
						min_counter++;
						if (min_counter == max_seeds) {
							return min_counter;
						}
					}
					__mmask8 select_mask = 1 << minimizer_pos;
					minimizer            = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                                .loc  = i + minimizer_pos,
					                                .str  = ((select_mask & strand_mask) != 0)};
					skip_mask &= ~select_mask;
					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(skip_mask, hash, minimum_vec);
				} else {
					__mmask8 select_mask = 1 << minimizer_pos;
					minimizer            = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                                .loc  = i + minimizer_pos,
					                                .str  = ((select_mask & strand_mask) != 0)};

					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(
					       skip_mask & (UINT8_MAX << (minimizer_pos + 1)), hash, minimum_vec);
				}
			}

			// If the last fetched k-mer is not in the first window
			if (first_window && i + 7 >= k + w - 1) {
				first_window = 0;
				if (minimizer.hash != UINT64_MAX) {
					uint32_t base_pos = ((k - 1) >> 3) << 3;
					// Buffer index of the last minimizer of the window
					unsigned minimizer_buf = (minimizer.loc - base_pos) >> 3;
					for (int j = 0; j < minimizer_buf; j++) {
						__mmask8 cmp_mask =
						    _mm512_mask_cmpeq_epu64_mask(mask_buf[j], minimum_vec, hash_buf[j]);
						while (cmp_mask != 0) {
							unsigned minimizer_pos = first_pos_array[cmp_mask];
							__mmask8 select_mask   = 1 << minimizer_pos;
							int strand             = ((select_mask & strand_buf[j]) != 0);
							uint32_t loc = base_pos + 8 * j + minimizer_pos + extracted_len;
							uint32_t real_loc =
							    get_real_location(loc, pattern_len, ones, ones_loc, shift);
							mm128_t min_tmp = (mm128_t){.x = minimizer.hash << 8 | k,
							                            .y = (uint64_t)rid << 32 |
							                                 real_loc << 1 | strand};
							kv_push(mm128_t, km, *p, min_tmp);
							min_counter++;
							if (min_counter == max_seeds) {
								return min_counter;
							}
							mask_buf[j] &= ~select_mask;
							cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    mask_buf[j], minimum_vec, hash_buf[j]);
						}
					}

					__m512i last_hash;
					uint8_t last_mask;
					uint8_t last_strand_mask;
					if (minimizer.loc < i) {
						last_hash        = hash_buf[minimizer_buf];
						last_mask        = mask_buf[minimizer_buf];
						last_strand_mask = strand_buf[minimizer_buf];
					} else {
						last_hash        = hash;
						last_mask        = skip_mask;
						last_strand_mask = strand_mask;
					}

					// Mask the k-mers after the current minimizer
					unsigned minimizer_pos  = (minimizer.loc - base_pos) & 0x7;
					__mmask8 in_window_mask = (UINT8_MAX >> (8 - minimizer_pos)) & last_mask;
					__mmask8 cmp_mask =
					    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, last_hash);
					while (cmp_mask != 0) {
						unsigned minimizer_pos = first_pos_array[cmp_mask];
						__mmask8 select_mask   = 1 << minimizer_pos;
						int strand             = ((select_mask & last_strand_mask) != 0);

						uint32_t loc =
						    base_pos + 8 * minimizer_buf + minimizer_pos + extracted_len;
						uint32_t real_loc =
						    get_real_location(loc, pattern_len, ones, ones_loc, shift);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | strand};
						kv_push(mm128_t, km, *p, min_tmp);
						min_counter++;
						if (min_counter == max_seeds) {
							return max_seeds;
						}

						in_window_mask &= ~select_mask;
						cmp_mask = _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec,
						                                        last_hash);
					}
				}
			}

			// If the current window doens't cover the last fetched k-mer
			while (minimizer.loc + w < len_early_exit && i + 7 - minimizer.loc >= w) {
				// Find the Minimimum in the next window

				// Mask the k-mers after the next window
				in_window_mask = skip_mask & (UINT8_MAX >> (i + 7 - minimizer.loc - w));

				// Mask the k-mers before the current minimizer
				unsigned buf_index = buf_pos;
				mask_buf[buf_index] &= UINT8_MAX << ((buf_size << 3) + 1 + minimizer.loc - i);
				if (minimizer.loc >= i - ((buf_size - 1) << 3)) {
					mask_buf[(buf_index + 1) % buf_size] &=
					    UINT8_MAX << (((buf_size - 1) << 3) + 1 + minimizer.loc - i);
				}

				uint64_t minimum = UINT64_MAX;
				for (unsigned j = 0; j < buf_size; j++) {
					uint64_t min_tmp =
					    _mm512_mask_reduce_min_epu64(mask_buf[buf_index], hash_buf[buf_index]);
					if (min_tmp < minimum) {
						minimum = min_tmp;
					}
					buf_index = (buf_index + 1) % buf_size;
				}
				uint64_t min_tmp = _mm512_mask_reduce_min_epu64(in_window_mask, hash);
				if (min_tmp < minimum) {
					minimum = min_tmp;
				}

				minimum_vec = _mm512_set1_epi64(minimum);

				// Get its location(s)
				if (minimum == UINT64_MAX) {
					if (minimizer.hash != UINT64_MAX) {
						uint32_t loc = minimizer.loc + extracted_len;
						uint32_t real_loc =
						    get_real_location(loc, pattern_len, ones, ones_loc, shift);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};

						kv_push(mm128_t, km, *p, min_tmp);
						min_counter++;
						if (min_counter == max_seeds) {
							return min_counter;
						}
					}
					minimizer.hash = UINT64_MAX;
					minimizer.loc  = minimizer.loc + w;
				} else {
					unsigned loc = i - 8 * buf_size;
					for (unsigned j = 0; j < buf_size; j++) {
						__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
						    mask_buf[buf_index], minimum_vec, hash_buf[buf_index]);
						while (cmp_mask != 0) {
							unsigned minimizer_pos = first_pos_array[cmp_mask];
							__mmask8 select_mask   = 1 << minimizer_pos;
							if (minimizer.hash != UINT64_MAX) {
								uint32_t real_loc = get_real_location(
								    minimizer.loc + extracted_len, pattern_len, ones,
								    ones_loc, shift);
								mm128_t min_tmp =
								    (mm128_t){.x = minimizer.hash << 8 | k,
								              .y = (uint64_t)rid << 32 | real_loc << 1 |
								                   minimizer.str};
								kv_push(mm128_t, km, *p, min_tmp);
								min_counter++;
								if (min_counter == max_seeds) {
									return min_counter;
								}
							}
							minimizer = (seed_t){
							    .hash = minimum,
							    .loc  = loc + minimizer_pos,
							    .str  = ((strand_buf[buf_index] & select_mask) != 0)};
							mask_buf[buf_index] &= ~select_mask;
							cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    mask_buf[buf_index], minimum_vec, hash_buf[buf_index]);
						}
						loc += 8;
						buf_index = (buf_index + 1) % buf_size;
					}
					__mmask8 cmp_mask =
					    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, hash);
					while (cmp_mask != 0) {
						unsigned minimizer_pos = first_pos_array[cmp_mask];
						__mmask8 select_mask   = 1 << minimizer_pos;
						if (minimizer.hash != UINT64_MAX) {
							uint32_t real_loc =
							    get_real_location(minimizer.loc + extracted_len,
							                      pattern_len, ones, ones_loc, shift);
							mm128_t min_tmp = (mm128_t){.x = minimizer.hash << 8 | k,
							                            .y = (uint64_t)rid << 32 |
							                                 real_loc << 1 | minimizer.str};
							kv_push(mm128_t, km, *p, min_tmp);
							min_counter++;
							if (min_counter == max_seeds) {
								return min_counter;
							}
						}
						minimizer = (seed_t){.hash = minimum,
						                     .loc  = i + minimizer_pos,
						                     .str  = ((strand_mask & select_mask) != 0)};
						skip_mask &= ~select_mask;
						in_window_mask &= ~select_mask;
						cmp_mask =
						    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, hash);
					}
				}

				// Check if there are minimizers inside of the new window
				in_window_mask = skip_mask;
				if (i + 7 - minimizer.loc >= w) {
					in_window_mask &= UINT8_MAX >> (i + 7 - minimizer.loc - (w - 1));
				}
				__mmask8 cmp_mask = _mm512_mask_cmple_epu64_mask(in_window_mask, hash, minimum_vec);
				while (cmp_mask != 0) {
					unsigned minimizer_pos = first_pos_array[cmp_mask];
					__mmask8 select_mask   = 1 << minimizer_pos;
					if (minimizer.hash != UINT64_MAX) {
						uint32_t loc = minimizer.loc + extracted_len;
						uint32_t real_loc =
						    get_real_location(loc, pattern_len, ones, ones_loc, shift);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};
						kv_push(mm128_t, km, *p, min_tmp);
						min_counter++;
						if (min_counter == max_seeds) {
							return min_counter;
						}
					}
					minimizer = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                     .loc  = i + minimizer_pos,
					                     .str  = ((strand_mask & select_mask) != 0)};
					skip_mask &= ~select_mask;
					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(skip_mask, hash, minimum_vec);
				}
			}
			hash_buf[buf_pos]   = hash;
			mask_buf[buf_pos]   = skip_mask;
			strand_buf[buf_pos] = strand_mask;
			buf_pos             = (buf_pos + 1) % buf_size;
		}
		if (len_early_exit >= w + k - 1 && minimizer.hash != UINT64_MAX) {
			uint32_t loc      = minimizer.loc + extracted_len;
			uint32_t real_loc = get_real_location(loc, pattern_len, ones, ones_loc, shift);
			mm128_t min_tmp   = (mm128_t){.x = minimizer.hash << 8 | k,
			                              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};
			kv_push(mm128_t, km, *p, min_tmp);
			min_counter++;
			if (min_counter == max_seeds) {
				return min_counter;
			}
		}
		extracted_len += len_early_exit + early_exit_offset;
	}
	return min_counter;
}
unsigned mm_sketch3(void *km, const char *str, const unsigned len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p,
                    const char *Z, int W, int shift2, uint32_t MAX_NB_SEEDS) {
	uint64_t mask = (1ULL << 2 * k) - 1;

	assert(len > 0 && (w > 0 && w < 256) &&
	       (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice

	// char pattern[10] = "111000111"; // determine the pattern
	// int pattern_len = 9;

	const char *pattern = Z; // determine the pattern
	int pattern_len     = W; // determine length of pattern

	int ones = 0;     // amount of ones in the pattern
	int ones_loc[40]; // the maximum number of ones in the pattern is 40 for now

	for (int g = 0; g < pattern_len; ++g) {
		if (pattern[g] == '1') {
			ones_loc[ones] = g;
			++ones; // count the ones in the pattern
		}
	}

	int shift = 0;
	if (shift2 < 0)
		shift = 0;
	else
		shift = shift2;

	// Len after appliying the pattern
	unsigned diet_len        = ((len - shift) / pattern_len) * ones;
	unsigned diet_len_offset = (len - shift) % pattern_len;
	for (unsigned i = 0; i < diet_len_offset; i++) {
		if (pattern[i] == '1') {
			diet_len += 1;
		}
	}

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		char *new_genome;
		new_genome = (char *)malloc(sizeof(char) * diet_len); // allocate memory for the shortened sequence
		// Print the new Genome to check
		fprintf(stderr, "New Read:\n");
		for (unsigned i = 0; i < diet_len; i++) {
			unsigned real_location = get_real_location(i, pattern_len, ones, ones_loc, shift);
			fprintf(stderr, "%c", str[real_location]);
		}
		fprintf(stderr, "\n");
		free(new_genome);
	}

	const __m512i kmer_shift = _mm512_set_epi64(16, 14, 12, 10, 8, 6, 4, 2);
	const __m512i c_shift    = _mm512_set_epi64(0, 2, 4, 6, 8, 10, 12, 14);
	const __m512i mask_vec   = _mm512_set1_epi64(mask);
	const unsigned buf_size  = (w / 8) + ((w % 8) != 0);

	for (uint32_t extracted_len = 0; extracted_len < diet_len;) {
		__m512i kmer[2]  = {_mm512_setzero_si512()};
		seed_t minimizer = {.hash = UINT64_MAX, .loc = k - 2, .str = 0};
		__m512i hash_buf[32];
		__mmask8 mask_buf[32];
		__mmask8 strand_buf[32];
		unsigned buf_pos = 0;

		uint32_t len_early_exit = diet_len - extracted_len;
		// printf("len_early_exit: %u\n", len_early_exit);
		// const uint8_t *seq_off     = seq + extracted_len;
		unsigned early_exit_offset = 0;

		int first_window    = 1;
		__m512i minimum_vec = _mm512_set1_epi64(minimizer.hash);
		for (uint32_t i = 0; i < len_early_exit; i += 8) {

			uint64_t seq_packed, seq_packed_rev;
			uint8_t n_mask;
			pack(str, &seq_packed, &seq_packed_rev, &n_mask, extracted_len + i, shift, pattern_len, k, ones,
			     ones_loc, len);

			// Compute the + and - kmers
			__m512i c = _mm512_set1_epi64(seq_packed);
			c         = _mm512_srlv_epi64(c, c_shift);
			kmer[0]   = _mm512_set1_epi64(_mm_extract_epi64(_mm512_extracti64x2_epi64(kmer[0], 3), 1));
			kmer[0]   = _mm512_sllv_epi64(kmer[0], kmer_shift);
			kmer[0]   = _mm512_or_epi64(kmer[0], c);
			kmer[0]   = _mm512_and_epi64(kmer[0], mask_vec);

			__m512i c_rev = _mm512_set1_epi64(seq_packed_rev);
			c_rev         = _mm512_sllv_epi64(c_rev, c_shift);
			kmer[1]       = _mm512_set1_epi64(_mm_extract_epi64(_mm512_extracti64x2_epi64(kmer[1], 3), 1));
			kmer[1]       = _mm512_srlv_epi64(kmer[1], kmer_shift);
			kmer[1]       = _mm512_or_epi64(kmer[1], c_rev);
			kmer[1]       = _mm512_and_epi64(kmer[1], mask_vec);

			// If there are ambiguous bases
			if (n_mask != 0) {
				// If there are multiple `N`s
				unsigned first_n_pos = first_pos_array[n_mask];
				unsigned last_n_pos  = last_pos_array[n_mask];
				if (i + first_n_pos <= k - 1) {
					len_early_exit = i + last_n_pos + 1;
					break;
				}
				len_early_exit = (i + first_n_pos > len_early_exit) ? len_early_exit : i + first_n_pos;
				early_exit_offset = last_n_pos - first_n_pos + 1;
			}

			// Skip before the first k-mers
			if (i + 8 < k) {
				continue;
			}

			// Skip symetric k-mers
			__mmask8 skip_mask = _mm512_cmpneq_epu64_mask(kmer[0], kmer[1]);

			// Skip before the first k-mers
			if (i < k - 1) {
				skip_mask &= UINT8_MAX << (k - 1 - i);
			}

			// Skip after the last k-mers
			if (i + 8 > len_early_exit) {
				skip_mask &= UINT8_MAX >> (i + 8 - len_early_exit);
				// printf("len_early_exit: %u, i: %u, mask: %x\n", len_early_exit, i,
				// skip_mask);
			}

			// Compute the hashes
			__mmask8 strand_mask    = _mm512_cmpgt_epu64_mask(kmer[0], kmer[1]);
			__m512i min_kmer_strand = _mm512_mask_blend_epi64(strand_mask, kmer[0], kmer[1]);
			__m512i hash            = hash64_avx(min_kmer_strand, mask);

			/*
			uint64_t debug_vec[8];
			_mm512_store_epi64(debug_vec, hash);
			for (unsigned k = 0; k < 8; k++) {
			        printf("0x%010lx\t%u\n", debug_vec[k], i + k);
			}
			printf("n_mask: %x\n", n_mask);
			printf("minimum: %lx, loc: %u\n", minimizer.hash, minimizer.loc);
			puts("");
			*/

			__mmask8 in_window_mask = skip_mask;

			if (i + 7 - minimizer.loc >= w) {
				in_window_mask &= UINT8_MAX >> (i + 7 - minimizer.loc - (w - 1));
			}
			__mmask8 cmp_mask = _mm512_mask_cmple_epu64_mask(in_window_mask, hash, minimum_vec);
			while (cmp_mask != 0) {
				unsigned minimizer_pos = first_pos_array[cmp_mask];
				// Only store the minimizer after the first window
				// printf("m_pos: %u\n", minimizer_pos);
				if (i + minimizer_pos >= k + w - 1) {
					// If previous minimizer in the first window, search for all the seeds
					// with the same hash in the first window
					if (first_window) {
						// printf("hash: %lx\n", minimizer.hash << 8 | k);
						first_window = 0;
						if (minimizer.hash != UINT64_MAX) {

							// Offset to compute the location of the k-mers
							uint32_t base_pos = (k - 1) & (UINT32_MAX << 3);
							// Buffer index of the last minimizer of the window
							unsigned minimizer_buf = (minimizer.loc - base_pos) >> 3;
							for (int j = 0; j < minimizer_buf; j++) {
								__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
								    mask_buf[j], minimum_vec, hash_buf[j]);
								while (cmp_mask != 0) {
									unsigned minimizer_pos =
									    first_pos_array[cmp_mask];
									__mmask8 select_mask = 1 << minimizer_pos;
									int strand =
									    ((select_mask & strand_buf[j]) != 0);
									uint32_t loc = base_pos + 8 * j +
									               minimizer_pos + extracted_len;
									uint32_t real_loc = get_real_location(
									    loc, pattern_len, ones, ones_loc, shift);
									mm128_t min_tmp =
									    (mm128_t){.x = minimizer.hash << 8 | k,
									              .y = (uint64_t)rid << 32 |
									                   real_loc << 1 | strand};
									kv_push(mm128_t, km, *p, min_tmp);
									if (p->n == MAX_NB_SEEDS) {
										return real_loc;
									}
									mask_buf[j] &= ~select_mask;
									cmp_mask = _mm512_mask_cmpeq_epu64_mask(
									    mask_buf[j], minimum_vec, hash_buf[j]);
								}
							}

							__m512i last_hash;
							uint8_t *last_mask;
							uint8_t last_strand_mask;
							if (minimizer.loc < i) {
								last_hash        = hash_buf[minimizer_buf];
								last_mask        = &mask_buf[minimizer_buf];
								last_strand_mask = strand_buf[minimizer_buf];
							} else {
								last_hash        = hash;
								last_mask        = &skip_mask;
								last_strand_mask = strand_mask;
							}
							// Mask the k-mers after the current minimizer
							unsigned minimizer_pos = (minimizer.loc - base_pos) & 0x7;
							__mmask8 in_window_mask =
							    (UINT8_MAX >> (8 - minimizer_pos)) & *last_mask;
							__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    in_window_mask, minimum_vec, last_hash);
							while (cmp_mask != 0) {
								unsigned minimizer_pos = first_pos_array[cmp_mask];
								__mmask8 select_mask   = 1 << minimizer_pos;
								int strand   = ((select_mask & last_strand_mask) != 0);
								uint32_t loc = base_pos + 8 * minimizer_buf +
								               minimizer_pos + extracted_len;
								uint32_t real_loc = get_real_location(
								    loc, pattern_len, ones, ones_loc, shift);
								mm128_t min_tmp = (mm128_t){
								    .x = minimizer.hash << 8 | k,
								    .y = (uint64_t)rid << 32 | real_loc << 1 | strand};
								kv_push(mm128_t, km, *p, min_tmp);
								if (p->n == MAX_NB_SEEDS) {
									return real_loc;
								}
								in_window_mask &= ~select_mask;
								*last_mask &= ~select_mask;
								cmp_mask = _mm512_mask_cmpeq_epu64_mask(
								    in_window_mask, minimum_vec, last_hash);
							}

							// Mask the current minimizer in skip_mask
							if (minimizer.loc >= i) {
								skip_mask &= ~(1 << minimizer_pos);
							}
						}
					}
					if (minimizer.hash != UINT64_MAX) {
						uint32_t loc = minimizer.loc + extracted_len;
						uint32_t real_loc =
						    get_real_location(loc, pattern_len, ones, ones_loc, shift);

						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};

						kv_push(mm128_t, km, *p, min_tmp);
						if (p->n == MAX_NB_SEEDS) {
							return real_loc;
						}
					}
					__mmask8 select_mask = 1 << minimizer_pos;
					minimizer            = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                                .loc  = i + minimizer_pos,
					                                .str  = ((select_mask & strand_mask) != 0)};
					skip_mask &= ~select_mask;
					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(skip_mask, hash, minimum_vec);
				} else {
					__mmask8 select_mask = 1 << minimizer_pos;
					minimizer            = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                                .loc  = i + minimizer_pos,
					                                .str  = ((select_mask & strand_mask) != 0)};

					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(
					       skip_mask & (UINT8_MAX << (minimizer_pos + 1)), hash, minimum_vec);
				}
			}

			// If the last fetched k-mer is not in the first window
			if (first_window && i + 7 >= k + w - 1) {
				first_window = 0;
				if (minimizer.hash != UINT64_MAX) {
					uint32_t base_pos = ((k - 1) >> 3) << 3;
					// Buffer index of the last minimizer of the window
					unsigned minimizer_buf = (minimizer.loc - base_pos) >> 3;
					for (int j = 0; j < minimizer_buf; j++) {
						__mmask8 cmp_mask =
						    _mm512_mask_cmpeq_epu64_mask(mask_buf[j], minimum_vec, hash_buf[j]);
						while (cmp_mask != 0) {
							unsigned minimizer_pos = first_pos_array[cmp_mask];
							__mmask8 select_mask   = 1 << minimizer_pos;
							int strand             = ((select_mask & strand_buf[j]) != 0);
							uint32_t loc = base_pos + 8 * j + minimizer_pos + extracted_len;
							uint32_t real_loc =
							    get_real_location(loc, pattern_len, ones, ones_loc, shift);
							mm128_t min_tmp = (mm128_t){.x = minimizer.hash << 8 | k,
							                            .y = (uint64_t)rid << 32 |
							                                 real_loc << 1 | strand};
							kv_push(mm128_t, km, *p, min_tmp);
							if (p->n == MAX_NB_SEEDS) {
								return real_loc;
							}
							mask_buf[j] &= ~select_mask;
							cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    mask_buf[j], minimum_vec, hash_buf[j]);
						}
					}

					__m512i last_hash;
					uint8_t last_mask;
					uint8_t last_strand_mask;
					if (minimizer.loc < i) {
						last_hash        = hash_buf[minimizer_buf];
						last_mask        = mask_buf[minimizer_buf];
						last_strand_mask = strand_buf[minimizer_buf];
					} else {
						last_hash        = hash;
						last_mask        = skip_mask;
						last_strand_mask = strand_mask;
					}

					// Mask the k-mers after the current minimizer
					unsigned minimizer_pos  = (minimizer.loc - base_pos) & 0x7;
					__mmask8 in_window_mask = (UINT8_MAX >> (8 - minimizer_pos)) & last_mask;
					__mmask8 cmp_mask =
					    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, last_hash);
					while (cmp_mask != 0) {
						unsigned minimizer_pos = first_pos_array[cmp_mask];
						__mmask8 select_mask   = 1 << minimizer_pos;
						int strand             = ((select_mask & last_strand_mask) != 0);

						uint32_t loc =
						    base_pos + 8 * minimizer_buf + minimizer_pos + extracted_len;
						uint32_t real_loc =
						    get_real_location(loc, pattern_len, ones, ones_loc, shift);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | strand};
						kv_push(mm128_t, km, *p, min_tmp);
						if (p->n == MAX_NB_SEEDS) {
							return real_loc;
						}

						in_window_mask &= ~select_mask;
						cmp_mask = _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec,
						                                        last_hash);
					}
				}
			}

			// If the current window doens't cover the last fetched k-mer
			while (minimizer.loc + w < len_early_exit && i + 7 - minimizer.loc >= w) {
				// Find the Minimimum in the next window

				// Mask the k-mers after the next window
				in_window_mask = skip_mask & (UINT8_MAX >> (i + 7 - minimizer.loc - w));

				// Mask the k-mers before the current minimizer
				unsigned buf_index = buf_pos;
				mask_buf[buf_index] &= UINT8_MAX << ((buf_size << 3) + 1 + minimizer.loc - i);
				if (minimizer.loc >= i - ((buf_size - 1) << 3)) {
					mask_buf[(buf_index + 1) % buf_size] &=
					    UINT8_MAX << (((buf_size - 1) << 3) + 1 + minimizer.loc - i);
				}

				uint64_t minimum = UINT64_MAX;
				for (unsigned j = 0; j < buf_size; j++) {
					uint64_t min_tmp =
					    _mm512_mask_reduce_min_epu64(mask_buf[buf_index], hash_buf[buf_index]);
					if (min_tmp < minimum) {
						minimum = min_tmp;
					}
					buf_index = (buf_index + 1) % buf_size;
				}
				uint64_t min_tmp = _mm512_mask_reduce_min_epu64(in_window_mask, hash);
				if (min_tmp < minimum) {
					minimum = min_tmp;
				}

				minimum_vec = _mm512_set1_epi64(minimum);

				// Get its location(s)
				if (minimum == UINT64_MAX) {
					if (minimizer.hash != UINT64_MAX) {
						uint32_t loc = minimizer.loc + extracted_len;
						uint32_t real_loc =
						    get_real_location(loc, pattern_len, ones, ones_loc, shift);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};

						kv_push(mm128_t, km, *p, min_tmp);
						if (p->n == MAX_NB_SEEDS) {
							return real_loc;
						}
					}
					minimizer.hash = UINT64_MAX;
					minimizer.loc  = minimizer.loc + w;
				} else {
					unsigned loc = i - 8 * buf_size;
					for (unsigned j = 0; j < buf_size; j++) {
						__mmask8 cmp_mask = _mm512_mask_cmpeq_epu64_mask(
						    mask_buf[buf_index], minimum_vec, hash_buf[buf_index]);
						while (cmp_mask != 0) {
							unsigned minimizer_pos = first_pos_array[cmp_mask];
							__mmask8 select_mask   = 1 << minimizer_pos;
							if (minimizer.hash != UINT64_MAX) {
								uint32_t real_loc = get_real_location(
								    minimizer.loc + extracted_len, pattern_len, ones,
								    ones_loc, shift);
								mm128_t min_tmp =
								    (mm128_t){.x = minimizer.hash << 8 | k,
								              .y = (uint64_t)rid << 32 | real_loc << 1 |
								                   minimizer.str};
								kv_push(mm128_t, km, *p, min_tmp);
								if (p->n == MAX_NB_SEEDS) {
									return real_loc;
								}
							}
							minimizer = (seed_t){
							    .hash = minimum,
							    .loc  = loc + minimizer_pos,
							    .str  = ((strand_buf[buf_index] & select_mask) != 0)};
							mask_buf[buf_index] &= ~select_mask;
							cmp_mask = _mm512_mask_cmpeq_epu64_mask(
							    mask_buf[buf_index], minimum_vec, hash_buf[buf_index]);
						}
						loc += 8;
						buf_index = (buf_index + 1) % buf_size;
					}
					__mmask8 cmp_mask =
					    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, hash);
					while (cmp_mask != 0) {
						unsigned minimizer_pos = first_pos_array[cmp_mask];
						__mmask8 select_mask   = 1 << minimizer_pos;
						if (minimizer.hash != UINT64_MAX) {
							uint32_t real_loc =
							    get_real_location(minimizer.loc + extracted_len,
							                      pattern_len, ones, ones_loc, shift);
							mm128_t min_tmp = (mm128_t){.x = minimizer.hash << 8 | k,
							                            .y = (uint64_t)rid << 32 |
							                                 real_loc << 1 | minimizer.str};
							kv_push(mm128_t, km, *p, min_tmp);
							if (p->n == MAX_NB_SEEDS) {
								return real_loc;
							}
						}
						minimizer = (seed_t){.hash = minimum,
						                     .loc  = i + minimizer_pos,
						                     .str  = ((strand_mask & select_mask) != 0)};
						skip_mask &= ~select_mask;
						in_window_mask &= ~select_mask;
						cmp_mask =
						    _mm512_mask_cmpeq_epu64_mask(in_window_mask, minimum_vec, hash);
					}
				}

				// Check if there are minimizers inside of the new window
				in_window_mask = skip_mask;
				if (i + 7 - minimizer.loc >= w) {
					in_window_mask &= UINT8_MAX >> (i + 7 - minimizer.loc - (w - 1));
				}
				__mmask8 cmp_mask = _mm512_mask_cmple_epu64_mask(in_window_mask, hash, minimum_vec);
				while (cmp_mask != 0) {
					unsigned minimizer_pos = first_pos_array[cmp_mask];
					__mmask8 select_mask   = 1 << minimizer_pos;
					if (minimizer.hash != UINT64_MAX) {
						uint32_t loc = minimizer.loc + extracted_len;
						uint32_t real_loc =
						    get_real_location(loc, pattern_len, ones, ones_loc, shift);
						mm128_t min_tmp =
						    (mm128_t){.x = minimizer.hash << 8 | k,
						              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};
						kv_push(mm128_t, km, *p, min_tmp);
						if (p->n == MAX_NB_SEEDS) {
							return real_loc;
						}
					}
					minimizer = (seed_t){.hash = extract_u64(hash, minimizer_pos),
					                     .loc  = i + minimizer_pos,
					                     .str  = ((strand_mask & select_mask) != 0)};
					skip_mask &= ~select_mask;
					minimum_vec = _mm512_set1_epi64(minimizer.hash);
					cmp_mask    = _mm512_mask_cmple_epu64_mask(skip_mask, hash, minimum_vec);
				}
			}
			hash_buf[buf_pos]   = hash;
			mask_buf[buf_pos]   = skip_mask;
			strand_buf[buf_pos] = strand_mask;
			buf_pos             = (buf_pos + 1) % buf_size;
		}
		if (len_early_exit >= w + k - 1 && minimizer.hash != UINT64_MAX) {
			uint32_t loc      = minimizer.loc + extracted_len;
			uint32_t real_loc = get_real_location(loc, pattern_len, ones, ones_loc, shift);
			mm128_t min_tmp   = (mm128_t){.x = minimizer.hash << 8 | k,
			                              .y = (uint64_t)rid << 32 | real_loc << 1 | minimizer.str};
			kv_push(mm128_t, km, *p, min_tmp);
			if (p->n == MAX_NB_SEEDS) {
				return real_loc;
			}
		}
		extracted_len += len_early_exit + early_exit_offset;
	}
	return len;
}

#else

void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p, const char *Z,
               int W) {
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = {UINT64_MAX, UINT64_MAX};
	// tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) &&
	       (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	// memset(&tq, 0, sizeof(tiny_queue_t));

	// int pattern[9] = {1,1,1, 0,0,0, 1,1,1}; // determine the pattern
	// int pattern_len = 9; // determine length of pattern

	const char *pattern = Z; // determine the pattern
	int pattern_len     = W; // determine length of pattern

	// calculate the new_size
	// int new_size;
	int ones = 0;     // amount of ones in the pattern
	int ones_loc[40]; // the maximum number of ones in the pattern is 40 for now
	// char *new_genome   = 0;
	for (int g = 0; g < pattern_len; ++g) {
		if (pattern[g] == '1') {
			ones_loc[ones] = g;
			++ones; // count the ones in the pattern
		}
	}

	// Len after appliying the pattern
	unsigned diet_len        = (len / pattern_len) * ones;
	unsigned diet_len_offset = len % pattern_len;
	for (unsigned i = 0; i < diet_len_offset; i++) {
		if (pattern[i] == '1') {
			diet_len += 1;
		}
	}

	/*if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
	        new_size = (len / pattern_len)*ones; // determine the new size
	        if(len % pattern_len != 0){ // loop if the pattern is interrupted in the middle
	                for(int m = 0; m < (len % pattern_len); ++m){
	                        if(pattern[m] == '1') ++new_size;
	                }
	        }
	        new_genome = (char *)malloc(sizeof(char) * new_size); // allocate memory for the shortened
	sequence int current = 0; // gives the position of the base in the new_genome for(int f = 0; f < len;
	++f){ if(pattern[(f%pattern_len)] == '1'){ new_genome[current] = str[f];
	                        ++current;
	                }
	        }

	        kv_resize(mm128_t, km, *p, p->n + new_size/w);
	        printf("============Sketch1: \n");
	        printf("Old Ref Genome: %s\n",str);
	        printf("New Ref Genome: ");
	        for(int pp = 0; pp < new_size; ++pp){
	                printf("%c", new_genome[pp]);
	        }
	        printf("\n");
	}*/

	for (i = l = buf_pos = min_pos = 0; i < diet_len; ++i) {
		unsigned real_location = get_real_location(i, pattern_len, ones, ones_loc, 0);
		int c                  = seq_nt4_table[(uint8_t)str[real_location]];
		mm128_t info           = {UINT64_MAX, UINT64_MAX};
		if (c < 4) { // not an ambiguous base
			int z;
			// TODO::: enable HPC later on
			/*if (is_hpc) {
			        int skip_len = 1;
			        if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
			                for (skip_len = 2; i + skip_len < len; ++skip_len)
			                        if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
			                                break;
			                i += skip_len - 1; // put $i at the end of the current homopolymer run
			        }
			        tq_push(&tq, skip_len);
			        kmer_span += skip_len;
			        if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else */
			kmer_span = l + 1 < k ? l + 1 : k;
			kmer[0]   = (kmer[0] << 2 | c) & mask;             // forward k-mer
			kmer[1]   = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
			++l;
			if (kmer[0] == kmer[1]) {
				info.x = UINT64_MAX;
			} else {                               // skip "symmetric k-mers" as we don't know it strand
				z = kmer[0] < kmer[1] ? 0 : 1; // strand
				if (l >= k && kmer_span < 256) {
					info.x = hash64(kmer[z], mask) << 8 | kmer_span;
					/*if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
					        printf("kmer[0] and kmer[1]\t");
					        printbits(kmer[0],64);
					        printf(", ");
					        printbits(kmer[1],64);
					        printf(", ");
					        printf("%d ", z);
					        printbits(info.x, 64);
					        printf("\n");
					}*/

					// printf(" new i: %d, old i: %d\n",i-k+1, old_location-k+1);
					// printf(" new i: %d, old i: %d\n",i, old_location);

					info.y = (uint64_t)rid << 32 | (uint32_t)real_location << 1 | z;
				}
			}
		} else {
			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				kv_push(mm128_t, km, *p, min);
			}
			l = 0, /*tq.count = tq.front = 0,*/ kmer_span = 0;
		}
		buf[buf_pos] = info;   // need to do this here as appropriate buf_pos and buf[buf_pos] are needed
		                       // below printf("k-mers: "); printf("%016llx\n", buf[buf_pos].x);
		                       // printbits(buf[buf_pos].x,24); // print w k-mers
		                       // printf("\n");
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) {
				// printbits(min.x, 64);
				// printf("\tLocation: ");
				// printf("%d\n",(min.y&0xff)>>1);
				kv_push(mm128_t, km, *p, min);
			}
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				// printbits(min.x, 64);
				// printf("\tLocation: ");
				// printf("\n%d",(min.y&0xff)>>1);
				kv_push(mm128_t, km, *p, min);
			}
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w;
			     ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x)
					min     = buf[j],
					min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) {  // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) {
						// printbits(buf[j].x, 64);
						// printf("\tLocation: ");
						// printf("%d\n",(buf[j].y&0xff)>>1);
						kv_push(mm128_t, km, *p, buf[j]);
					}
				}
				for (j = 0; j <= buf_pos; ++j) {
					if (min.x == buf[j].x && min.y != buf[j].y) {
						// printbits(buf[j].x, 64);
						// printf("\tLocation: ");
						// printf("%d\n",(buf[j].y&0xff)>>1);
						kv_push(mm128_t, km, *p, buf[j]);
					}
				}
			}
		}

		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because
			                                     // identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					// printbits(buf[j].x, 64);
					// printf("\tLocation: ");
					// printf("%d\n",(buf[j].y&0xff)>>1);
					kv_push(mm128_t, km, *p, buf[j]);
				}
			}
			for (j = 0; j < buf_pos; ++j) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					// printbits(buf[j].x, 64);
					// printf("\tLocation: ");
					// printf("%d\n",(buf[j].y&0xff)>>1);
					kv_push(mm128_t, km, *p, buf[j]);
				}
			}
		}

		if (++buf_pos == w) buf_pos = 0;
	}
	if (l > w + k - 1 && min.x != UINT64_MAX) {
		// printbits(min.x, 64);
		// printf("\tLocation: ");
		// printf("%d\n",(min.y&0xff)>>1);
		kv_push(mm128_t, km, *p, min);
	}
	// if (mm_dbg_flag & MM_DBG_PRINT_SEED) free(new_genome);
}

static inline unsigned mm_sketch2_sub(void *km, const char *str, const unsigned diet_len, const int w, const int k,
                                      const uint32_t rid, mm128_v *p, const unsigned shift, const unsigned max_seeds,
                                      const int ones, const int *const ones_loc, const unsigned pattern_len,
                                      const uint64_t shift1, const uint64_t mask, unsigned len) {
	uint64_t kmer[2] = {0, 0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = {UINT64_MAX, UINT64_MAX};
	memset(buf, 0xff, w * 16);

	unsigned min_counter = 0;

	for (i = l = buf_pos = min_pos = 0; i < diet_len; ++i) {

		unsigned real_location = get_real_location(i, pattern_len, ones, ones_loc, shift);
		int c                  = seq_nt4_table[(uint8_t)str[real_location]];
		mm128_t info           = {UINT64_MAX, UINT64_MAX};
		if (c < 4) { // not an ambiguous base
			int z;
			// TODO:: Enable HPC
			/*
		if (is_hpc) {
			int skip_len = 1;
			if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
			        for (skip_len = 2; i + skip_len < len; ++skip_len)
			                if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
			                        break;
			        i += skip_len - 1; // put $i at the end of the current homopolymer run
			}
			tq_push(&tq, skip_len);
			kmer_span += skip_len;
			if (tq.count > k) kmer_span -= tq_shift(&tq);
		} else*/
			kmer_span = l + 1 < k ? l + 1 : k;
			kmer[0]   = (kmer[0] << 2 | c) & mask;             // forward k-mer
			kmer[1]   = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
			++l;
			if (kmer[0] == kmer[1]) {
				info.x = UINT64_MAX;
			} else {                               // skip "symmetric k-mers" as we don't know it strand
				z = kmer[0] < kmer[1] ? 0 : 1; // strand
				if (l >= k && kmer_span < 256) {
					info.x = hash64(kmer[z], mask) << 8 | kmer_span;
					info.y = (uint64_t)rid << 32 | (uint32_t)real_location << 1 | z;
				}
			}
		} else {
			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				kv_push(mm128_t, km, *p, min);
				min_counter++;
				if (min_counter == max_seeds) {
					return min_counter;
				}
			}
			l = 0, /*tq.count = tq.front = 0,*/ kmer_span = 0;
		}
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are
		                     // needed below

		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) {
				kv_push(mm128_t, km, *p, min);
				min_counter++;
				if (min_counter == max_seeds) {
					return min_counter;
				}
			}
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window

			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				kv_push(mm128_t, km, *p, min);
				min_counter++;
				if (min_counter == max_seeds) {
					return min_counter;
				}
			}
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w;
			     ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) {
					min     = buf[j];
					min_pos = j; // >= is important s.t. min is always the closest k-mer
				}
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) {  // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) {
						kv_push(mm128_t, km, *p, buf[j]);
						min_counter++;
						if (min_counter == max_seeds) {
							return min_counter;
						}
					}
				}
				for (j = 0; j <= buf_pos; ++j) {
					if (min.x == buf[j].x && min.y != buf[j].y) {
						kv_push(mm128_t, km, *p, buf[j]);
						min_counter++;
						if (min_counter == max_seeds) {
							return min_counter;
						}
					}
				}
			}
		}

		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because
			                                     // identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					kv_push(mm128_t, km, *p, buf[j]);
					min_counter++;
					if (min_counter == max_seeds) {
						return min_counter;
					}
				}
			}
			for (j = 0; j < buf_pos; ++j) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					kv_push(mm128_t, km, *p, buf[j]);
					min_counter++;
					if (min_counter == max_seeds) {
						return min_counter;
					}
				}
			}
		}
		buf_pos = (buf_pos == w - 1) ? 0 : buf_pos + 1;
	}
	if (l >= w + k - 1 && min.x != UINT64_MAX) {
		kv_push(mm128_t, km, *p, min);
		min_counter++;
		if (min_counter == max_seeds) {
			return min_counter;
		}
	}
	return min_counter;
}

unsigned mm_sketch3(void *km, const char *str, const unsigned len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p,
                    const char *Z, int W, int shift2, uint32_t MAX_NB_SEEDS) {
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = {UINT64_MAX, UINT64_MAX};
	// tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) &&
	       (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	// memset(&tq, 0, sizeof(tiny_queue_t));

	// char pattern[10] = "111000111"; // determine the pattern
	// int pattern_len = 9;

	const char *pattern = Z; // determine the pattern
	int pattern_len     = W; // determine length of pattern

	int ones = 0;     // amount of ones in the pattern
	int ones_loc[40]; // the maximum number of ones in the pattern is 40 for now

	for (int g = 0; g < pattern_len; ++g) {
		if (pattern[g] == '1') {
			ones_loc[ones] = g;
			++ones; // count the ones in the pattern
		}
	}
	int shift = 0;
	if (shift2 < 0)
		shift = 0;
	else
		shift = shift2;

	// Len after appliying the pattern
	unsigned diet_len        = ((len - shift) / pattern_len) * ones;
	unsigned diet_len_offset = (len - shift) % pattern_len;
	for (unsigned i = 0; i < diet_len_offset; i++) {
		if (pattern[i] == '1') {
			diet_len += 1;
		}
	}

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		char *new_genome;
		new_genome = (char *)malloc(sizeof(char) * diet_len); // allocate memory for the shortened sequence
		// Print the new Genome to check
		fprintf(stderr, "New Read:\n");
		for (unsigned i = 0; i < diet_len; i++) {
			unsigned real_location = get_real_location(i, pattern_len, ones, ones_loc, shift);
			fprintf(stderr, "%c", str[real_location]);
		}
		fprintf(stderr, "\n");
		free(new_genome);
	}

	l       = 0;
	buf_pos = 0;
	min_pos = 0;

	// printf("k: %u, w: %u\n", k, w);
	for (i = 0; i < diet_len; ++i) {
		unsigned real_location = get_real_location(i, pattern_len, ones, ones_loc, shift);
		int c                  = seq_nt4_table[(uint8_t)str[real_location]];
		mm128_t info           = {UINT64_MAX, UINT64_MAX};
		if (c < 4) { // not an ambiguous base
			int z;
			// TODO:: Enable HPC
			/*
		if (is_hpc) {
			int skip_len = 1;
			if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
			        for (skip_len = 2; i + skip_len < len; ++skip_len)
			                if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
			                        break;
			        i += skip_len - 1; // put $i at the end of the current homopolymer run
			}
			tq_push(&tq, skip_len);
			kmer_span += skip_len;
			if (tq.count > k) kmer_span -= tq_shift(&tq);
		} else*/
			kmer_span = l + 1 < k ? l + 1 : k;
			kmer[0]   = (kmer[0] << 2 | c) & mask;             // forward k-mer
			kmer[1]   = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
			++l;
			if (kmer[0] == kmer[1]) {
				info.x = UINT64_MAX;
			} else {                               // skip "symmetric k-mers" as we don't know it strand
				z = kmer[0] < kmer[1] ? 0 : 1; // strand
				if (l >= k && kmer_span < 256) {
					info.x = hash64(kmer[z], mask) << 8 | kmer_span;
					// printf("infox: %lu, buf_pos: %u\n", hash64(kmer[z], mask), buf_pos);

					// printf(" new i: %d, old i: %d\n",i-k+1, old_location-k+1);
					// printf(" new i: %d, old i: %d\n",i, old_location);
					// info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
					info.y = (uint64_t)rid << 32 | (uint32_t)real_location << 1 | z;
					// printf("i: %u\t old_location: %u\n", i, old_location);
				}
			}
		} else {
			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				kv_push(mm128_t, km, *p, min);
				if (p->n == MAX_NB_SEEDS) {
					return (uint32_t)(min.y >> 1);
				}
			}
			l = 0, /*tq.count = tq.front = 0,*/ kmer_span = 0;
		}
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are
		                     // needed below

		/*
		if (min.x == 0x21be060c && (min.y >> 1) == 4467) {
		        printf("new_min: %lx, buf_pos: %u, min_pos:%u\n", info.x >> 8, buf_pos, min_pos);
		}
		*/
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) {
				//					printbits(min.x, 64);
				//					printf("\tLocation: ");
				//					printf("%d",(min.y&0xff)>>1);
				//					printf("\tStrand:
				//%c","+-"[min.y&0x00000001]); printf("\n");
				kv_push(mm128_t, km, *p, min);
				if (p->n == MAX_NB_SEEDS) {
					return (uint32_t)(min.y >> 1);
				}
			}
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window

			if (l >= w + k - 1 && min.x != UINT64_MAX) {
				//					printbits(min.x, 64);
				//					printf("\tLocation: ");
				//					printf("%d",(min.y&0xff)>>1);
				//					printf("\tStrand:
				//%c","+-"[min.y&0x00000001]); printf("\n");
				kv_push(mm128_t, km, *p, min);
				if (p->n == MAX_NB_SEEDS) {
					return (uint32_t)(min.y >> 1);
				}
			}
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w;
			     ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) {
					min     = buf[j];
					min_pos = j; // >= is important s.t. min is always the closest k-mer
				}
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) {  // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) {
						//							printbits(buf[j].x,
						// 64);
						// printf("\tLocation: ");
						// printf("%d",(buf[j].y&0xff)>>1);
						// printf("\tStrand:
						//%c","+-"[buf[j].y&0x00000001]);
						// printf("\n");
						kv_push(mm128_t, km, *p, buf[j]);
						if (p->n == MAX_NB_SEEDS) {
							return (uint32_t)(buf[j].y >> 1);
						}
					}
				}
				for (j = 0; j <= buf_pos; ++j) {
					if (min.x == buf[j].x && min.y != buf[j].y) {
						//							printbits(buf[j].x,
						// 64);
						// printf("\tLocation: ");
						// printf("%d",(buf[j].y&0xff)>>1);
						// printf("\tStrand:
						//%c","+-"[buf[j].y&0x00000001]);
						// printf("\n");
						kv_push(mm128_t, km, *p, buf[j]);
						if (p->n == MAX_NB_SEEDS) {
							return (uint32_t)(buf[j].y >> 1);
						}
					}
				}
			}
		}

		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because
			                                     // identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					//						printbits(buf[j].x,
					// 64);
					// printf("\tLocation:
					// "); printf("%d",(buf[j].y&0xff)>>1);
					// printf("\tStrand:
					//%c","+-"[buf[j].y&0x00000001]);
					// printf("\n");
					kv_push(mm128_t, km, *p, buf[j]);
					if (p->n == MAX_NB_SEEDS) {
						return (uint32_t)(buf[j].y >> 1);
					}
				}
			}
			for (j = 0; j < buf_pos; ++j) {
				if (min.x == buf[j].x && buf[j].y != min.y) {
					//						printbits(buf[j].x,
					// 64);
					// printf("\tLocation:
					// "); printf("%d",(buf[j].y&0xff)>>1);
					// printf("\tStrand:
					//%c","+-"[buf[j].y&0x00000001]);
					// printf("\n");
					kv_push(mm128_t, km, *p, buf[j]);
					if (p->n == MAX_NB_SEEDS) {
						return (uint32_t)(buf[j].y >> 1);
					}
				}
			}
		}
		buf_pos = (buf_pos == w - 1) ? 0 : buf_pos + 1;
	}
	if (l >= w + k - 1 && min.x != UINT64_MAX) {
		//		printbits(min.x, 64);
		//		printf("\tLocation: ");
		//		printf("%d",(min.y&0xff)>>1);
		//		printf("\tStrand: %c","+-"[min.y&0x00000001]);
		//		printf("\n");
		kv_push(mm128_t, km, *p, min);
		if (p->n == MAX_NB_SEEDS) {
			return (uint32_t)(min.y >> 1);
		}
	}
	return len;
}

#endif

mm_pattern_t mm_sketch2(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p,
                        const char *Z, int W, const float max_seeds) {
	const char *pattern = Z; // determine the pattern
	int pattern_len     = W; // determine length of pattern
	// char pattern[10] = "111000111"; // determine the pattern
	// int pattern_len = 9;

	// int f = 0;

	uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1;

	assert(len > 0 && (w > 0 && w < 256) &&
	       (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice

	mm_pattern_t mm_pattern       = {.n = pattern_len};
	mm_pattern.shift_seeds_number = (uint32_t *)kmalloc(km, pattern_len * sizeof(uint32_t));

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "Read before pattern  : %s\n", str);
	}

	int ones = 0;     // amount of ones in the pattern
	int ones_loc[40]; // the maximum number of ones in the pattern is 40 for now

	for (int g = 0; g < pattern_len; ++g) {
		if (pattern[g] == '1') {
			ones_loc[ones] = g;
			++ones; // count the ones in the pattern
		}
	}

	unsigned len_crop;
	uint32_t max_nb_seeds;

	if (max_seeds < 1) {
		len_crop     = (unsigned)((float)max_seeds * len);
		max_nb_seeds = UINT32_MAX;
	} else {
		len_crop     = len;
		max_nb_seeds = max_seeds;
	}

	for (int shift = 0; shift < pattern_len; ++shift) {
		// Len after appliying the pattern
		unsigned diet_len        = ((len_crop - shift) / pattern_len) * ones;
		unsigned diet_len_offset = (len_crop - shift) % pattern_len;
		for (unsigned i = 0; i < diet_len_offset; i++) {
			if (pattern[i] == '1') {
				diet_len += 1;
			}
		}

		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			unsigned diet_len        = ((len - shift) / pattern_len) * ones;
			unsigned diet_len_offset = (len - shift) % pattern_len;
			for (unsigned i = 0; i < diet_len_offset; i++) {
				if (pattern[i] == '1') {
					diet_len += 1;
				}
			}
			char *new_genome;
			new_genome =
			    (char *)malloc(sizeof(char) * diet_len); // allocate memory for the shortened sequence
			// Print the new Genome to check
			fprintf(stderr, "Read after pattern(%u) (len: %u):\n", shift, diet_len);
			for (unsigned i = 0; i < diet_len; i++) {
				unsigned real_location = get_real_location(i, pattern_len, ones, ones_loc, shift);
				fprintf(stderr, "%c", str[real_location]);
			}
			fprintf(stderr, "\n");
			free(new_genome);
		}
		mm_pattern.shift_seeds_number[shift] =
		    mm_sketch2_sub(km, str, diet_len, w, k, rid, p, shift, max_nb_seeds, ones, ones_loc, pattern_len,
		                   shift1, mask, len);

		if (max_nb_seeds == UINT32_MAX) {
			len_crop     = len;
			max_nb_seeds = mm_pattern.shift_seeds_number[shift];
		}
	}
	return mm_pattern;
}

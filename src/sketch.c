#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <inttypes.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"
#include <math.h>

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
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

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
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

void printBinary (uint64_t n)
{
	for (int i = 63; 0 <= i; i--) {
		printf("%c", (n & (1 << i)) ? '1' : '0');
		if (i%4==0) printf("_");
	}
}

void printbits(uint64_t number, unsigned int num_bits_to_print)
{
	if (number || num_bits_to_print > 0) {
		printbits(number >> 1, num_bits_to_print - 1);
		printf("%d", number & 1);
	}
}

void bits2let(int number){
	int *a[24] = {};
	int bin = 0;
	int i;

	int bits = (log(number)/log(2)) + 1;

	for(i = 0; 0 < number; ++i){
		a[i] = number%2;
		number=number/2;
	}
	for(i = bits; i < 24; ++i)
		a[i] = 0;
	printf("bits2let: ");
	/*for(i -= 1; i >= 0; --i)
		printf("%d", a[i]);
	printf("\n");*/

	//printf("letters: ");
	for(i = 23; i >= 0; i -= 2){
		bin = ((int)a[i])*2 + ((int)a[i-1]);
		switch(bin){
			case 0: printf("A"); break;
			case 1: printf("C"); break;
			case 2: printf("G"); break;
			case 3: printf("T"); break;
			default: printf("E"); break;
		}
	}
	printf("\n");
}

void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p, const char *Z, int W)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));

	//int pattern[9] = {1,1,1, 0,0,0, 1,1,1}; // determine the pattern
	//int pattern_len = 9; // determine length of pattern

	char *pattern = Z; // determine the pattern
	int pattern_len = W; // determine length of pattern

	// calculate the new_size
	int new_size;
	int ones = 0; // amount of ones in the pattern
	int *ones_loc = (int *)malloc( 40*sizeof(int)); // the maximum number of ones in the pattern is 40 for now
	int ones_loc_index=1;
	int old_location=0;
	int ones_index=0;
	int offset=0;
	char *new_genome = 0;
	for(int g = 0; g < pattern_len; ++g){
		if(pattern[g] == '1') {
			ones_loc[ones]=g;
			++ones; // count the ones in the pattern
		}
	}

	/*if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		new_size = (len / pattern_len)*ones; // determine the new size
		if(len % pattern_len != 0){ // loop if the pattern is interrupted in the middle
			for(int m = 0; m < (len % pattern_len); ++m){
				if(pattern[m] == '1') ++new_size;
			}
		}
		new_genome = (char *)malloc(sizeof(char) * new_size); // allocate memory for the shortened sequence
		int current = 0; // gives the position of the base in the new_genome
		for(int f = 0; f < len; ++f){
			if(pattern[(f%pattern_len)] == '1'){
				new_genome[current] = str[f];
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

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {

		if(pattern[(i%pattern_len)] == '1'){

			int c = seq_nt4_table[(uint8_t)str[i]];
			mm128_t info = { UINT64_MAX, UINT64_MAX };
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
				kmer_span = l + 1 < k? l + 1 : k;
				kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
				kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
				if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
				z = kmer[0] < kmer[1]? 0 : 1; // strand
				++l;
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
					ones_index=(ones_loc_index) % ones;
					old_location= offset + ((ones_index==0)?ones_loc[ones-1] : ones_loc[ones_index-1]);

					//printf(" new i: %d, old i: %d\n",i-k+1, old_location-k+1);
					//printf(" new i: %d, old i: %d\n",i, old_location);

					info.y = (uint64_t)rid<<32 | (uint32_t)old_location<<1 | z;
					ones_loc_index++;
					offset = offset + ((ones_index==0)?pattern_len:0);
				}
			} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
			buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
			//printf("k-mers: ");
			//printf("%016llx\n", buf[buf_pos].x);
			//printbits(buf[buf_pos].x,24); // print w k-mers
			//printf("\n");
			if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
				for (j = buf_pos + 1; j < w; ++j){
					if (min.x == buf[j].x && buf[j].y != min.y){
						//printbits(buf[j].x, 64);
						//printf("\tLocation: ");
						//printf("%d\n",(buf[j].y&0xff)>>1);
						kv_push(mm128_t, km, *p, buf[j]);
					}
				}
				for (j = 0; j < buf_pos; ++j){
					if (min.x == buf[j].x && buf[j].y != min.y){
						//printbits(buf[j].x, 64);
						//printf("\tLocation: ");
						//printf("%d\n",(buf[j].y&0xff)>>1);
						kv_push(mm128_t, km, *p, buf[j]);
					}
				}
			}
			if (info.x <= min.x) { // a new minimum; then write the old min
				if (l >= w + k && min.x != UINT64_MAX){
					//printbits(min.x, 64);
					//printf("\tLocation: ");
					//printf("%d\n",(min.y&0xff)>>1);
					kv_push(mm128_t, km, *p, min);
				}
				min = info, min_pos = buf_pos;
			} else if (buf_pos == min_pos) { // old min has moved outside the window
				if (l >= w + k - 1 && min.x != UINT64_MAX){
					//printbits(min.x, 64);
					//printf("\tLocation: ");
					//printf("\n%d",(min.y&0xff)>>1);
					kv_push(mm128_t, km, *p, min);
				}
				for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
					if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
				for (j = 0; j <= buf_pos; ++j)
					if (min.x >= buf[j].x) min = buf[j], min_pos = j;
				if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
					for (j = buf_pos + 1; j < w; ++j){ // these two loops make sure the output is sorted
						if (min.x == buf[j].x && min.y != buf[j].y){
							//printbits(buf[j].x, 64);
							//printf("\tLocation: ");
							//printf("%d\n",(buf[j].y&0xff)>>1);
							kv_push(mm128_t, km, *p, buf[j]);
						}
					}
					for (j = 0; j <= buf_pos; ++j){
						if (min.x == buf[j].x && min.y != buf[j].y){
							//printbits(buf[j].x, 64);
							//printf("\tLocation: ");
							//printf("%d\n",(buf[j].y&0xff)>>1);
							kv_push(mm128_t, km, *p, buf[j]);
						}
					}
				}
			}
			if (++buf_pos == w) buf_pos = 0;
		}
	}
	if (min.x != UINT64_MAX){
		//printbits(min.x, 64);
		//printf("\tLocation: ");
		//printf("%d\n",(min.y&0xff)>>1);
		kv_push(mm128_t, km, *p, min);
	}
	//if (mm_dbg_flag & MM_DBG_PRINT_SEED) free(new_genome);
	free(ones_loc);
}


mm_pattern_versions *mm_sketch2(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p, const char *Z, int W, int max_seeds)
{
	char *pattern = Z; // determine the pattern
	int pattern_len=W; // determine length of pattern
	//char pattern[10] = "111000111"; // determine the pattern
	//int pattern_len = 9; 

	int new_size=0;
	int current = 0; // gives the position of the base in the new_genome
	//int f = 0;

	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1;
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));

	for(int f = 0; f < len;++f){
		if(pattern[(f % pattern_len)] == '1'){
			++new_size;
		}
	}

	kv_resize(mm128_t, km, *p, pattern_len*(p->n + new_size/w));

	mm_pattern_versions *mm_pattern_versions1=0;
	mm_pattern_versions1 = (mm_pattern_versions*)calloc(pattern_len*sizeof(uint32_t)+1, sizeof(mm_pattern_versions));
	mm_pattern_versions1->n_shift_seeds_number  = pattern_len;

	char *new_genome=0;

	for(int shift = 0; shift < pattern_len; ++shift){

		uint64_t kmer[2] = {0,0};
		kmer_span = 0;

		if(shift>0){
			new_size=0;
			for(int f = shift; f < len;++f){
				if(pattern[((f-shift)%pattern_len)] == '1'){
					++new_size;
				}
			}
		}


		if (shift ==0)
			new_genome = (char *)malloc(sizeof(char) * new_size); // allocate memory for the shortened sequence shift 0

		current = 0;

		//printf("Sketch2\t shift: %d\tpattern length: %d\t new_size: %d\tlen: %d\n", shift, pattern_len, new_size, len);
		//printf("Pattern: ");
		for(int f = shift; f < len;++f){
			//printf("%c",pattern[((f-shift)%pattern_len)]);
			//printf("%d_%c_%d_%d\n",pattern[((f-shift)%pattern_len)], str[f],f, ((f-shift)%pattern_len));
			if(pattern[((f-shift)%pattern_len)] == '1'){
				new_genome[current] = str[f];
				++current;
			}
		}
		//printf("\n");

		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			//Print the new Genome to check
			printf("Read before pattern  : %s\n",str);
			printf("Read after pattern(%d): ", shift);
			for(int pp = 0; pp < new_size; ++pp){
				printf("%c", new_genome[pp]);
			}
			printf("\n");
		}

		for (i = l = buf_pos = min_pos = 0; i < new_size; ++i) { // k + k/pattern_len = k*pattern_len + k = k*(pattern_len + 1)
			int c = seq_nt4_table[(uint8_t)new_genome[i]];
			//printf("%c\t%d\n", new_genome[i],c); // print the k-mers
			//if ((i+1) % k == 0) printf(" %d\n", k); // give the length of a k-mer
			mm128_t info = { UINT64_MAX, UINT64_MAX };
			if (c < 4) { // not an ambiguous base
				int z;
				if (is_hpc) {
					int skip_len = 1;
					if (i + 1 < new_size && seq_nt4_table[(uint8_t)new_genome[i + 1]] == c) {
						for (skip_len = 2; i + skip_len < new_size; ++skip_len)
							if (seq_nt4_table[(uint8_t)new_genome[i + skip_len]] != c)
								break;
						i += skip_len - 1; // put $i at the end of the current homopolymer run
					}
					tq_push(&tq, skip_len);
					kmer_span += skip_len;
					if (tq.count > k) kmer_span -= tq_shift(&tq);
				} else kmer_span = l + 1 < k? l + 1 : k;
				kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
				kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
				if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
				z = kmer[0] < kmer[1]? 0 : 1; // strand
				++l;
				if (l >= k && kmer_span < 256) {
					info.x = hash64(kmer[z], mask) << 8 | kmer_span;
					//printbits(kmer[z],64);
					//printf("\n");
					//printbits(kmer[z],24);
					//printf("\n");
					//printf("kmer[0]\t");
					//printbits(kmer[0],64);
					//printf("\n");
					//printf("kmer[1]\t");
					//printbits(kmer[1],24);
					//printf("\n");
					//printbits(info.x,64);
					//printf("\n");
					//bits2let(kmer[0]);
					//printf("\n");
					info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
				}
			} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
			buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below

			//printf("%016llx\n", buf[buf_pos].x);
			//printbits(buf[buf_pos].x,24); // print w k-mers
			//printf("\n");
			if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
				for (j = buf_pos + 1; j < w; ++j){
					if (min.x == buf[j].x && buf[j].y != min.y){
						//printbits(buf[j].x, 64);
						//printf("\n");
						kv_push(mm128_t, km, *p, buf[j]);
						mm_pattern_versions1->shift_seeds_number[shift]++;
					}
				}
				for (j = 0; j < buf_pos; ++j){
					if (min.x == buf[j].x && buf[j].y != min.y){
						//printbits(buf[j].x, 64);
						//printf("\n");
						kv_push(mm128_t, km, *p, buf[j]);
						mm_pattern_versions1->shift_seeds_number[shift]++;
					}
				}
			}
			if (info.x <= min.x) { // a new minimum; then write the old min
				if (l >= w + k && min.x != UINT64_MAX){
					//printbits(min.x, 64);
					//printf("\n");
					kv_push(mm128_t, km, *p, min);
					mm_pattern_versions1->shift_seeds_number[shift]++;
				}
				min = info, min_pos = buf_pos;
			} else if (buf_pos == min_pos) { // old min has moved outside the window
				if (l >= w + k - 1 && min.x != UINT64_MAX){
					//printbits(min.x, 64);
					//printf("\n");
					kv_push(mm128_t, km, *p, min);
					mm_pattern_versions1->shift_seeds_number[shift]++;
				}
				for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
					if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
				for (j = 0; j <= buf_pos; ++j)
					if (min.x >= buf[j].x) min = buf[j], min_pos = j;
				if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
					for (j = buf_pos + 1; j < w; ++j){ // these two loops make sure the output is sorted
						if (min.x == buf[j].x && min.y != buf[j].y){
							//printbits(buf[j].x, 64);
							//printf("\n");
							kv_push(mm128_t, km, *p, buf[j]);
							mm_pattern_versions1->shift_seeds_number[shift]++;
						}
					}
					for (j = 0; j <= buf_pos; ++j){
						if (min.x == buf[j].x && min.y != buf[j].y){
							//printbits(buf[j].x, 64);
							//printf("\n");
							kv_push(mm128_t, km, *p, buf[j]);
							mm_pattern_versions1->shift_seeds_number[shift]++;
						}
					}
				}
			}
			if (++buf_pos == w) buf_pos = 0;
			if (mm_pattern_versions1->shift_seeds_number[shift]>= max_seeds) goto SKIPSEEDS;
		}
		SKIPSEEDS:
		if (min.x != UINT64_MAX){
			//printbits(min.x, 64);
			//printf("\n");
			kv_push(mm128_t, km, *p, min);
			mm_pattern_versions1->shift_seeds_number[shift]++;
		}

		//printf("shift_seeds_number[shift]: %d\n",mm_pattern_versions1->shift_seeds_number[shift]);
	}
	free(new_genome);
	return mm_pattern_versions1;
}



void mm_sketch3(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p, const char *Z, int W, int shift2)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));

	//char pattern[10] = "111000111"; // determine the pattern
	//int pattern_len = 9; 

	char *pattern = Z; // determine the pattern
	int pattern_len=W; // determine length of pattern

	// calculate the new_size
	int new_size = 0;
	int ones = 0; // amount of ones in the pattern
	int *ones_loc = (int *)malloc(40*sizeof(int)); // the maximum number of ones in the pattern is 40 for now
	int ones_loc_index=1;
	int old_location=0;
	int ones_index=0;
	int offset=0;

	for(int g = 0; g < pattern_len; ++g){
		if(pattern[g] == '1') {
			ones_loc[ones]=g;
			++ones; // count the ones in the pattern
		}
	}
	int shift =0;
	if (shift2<0) shift=0;
	else shift=shift2;
	char *new_genome = 0;

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {

		for(int f = shift; f < len;++f){
			if(pattern[((f-shift)%pattern_len)] == '1'){
				++new_size;
			}
		}
		new_genome = (char *)malloc(sizeof(char) * new_size); // allocate memory for the shortened sequence
		int current = 0; // gives the position of the base in the new_genome
		for(int f = shift; f < len;++f){
			//printf("%d_%c_%d_%d\n",pattern[((f-shift)%pattern_len)], str[f],f, ((f-shift)%pattern_len));
			if(pattern[((f-shift)%pattern_len)] == '1'){
				new_genome[current] = str[f];
				++current;
			}
		}
		kv_resize(mm128_t, km, *p, p->n + new_size/w);
		// Print the new Genome to check
		printf("\nNew Ref Genome: ");
		for(int pp = 0; pp < new_size; ++pp){
			printf("%c", new_genome[pp]);
			//if ((pp+1) % k == 0) printf("\t");
		}
		printf("\n");
	}

	l = 0;
	buf_pos = 0;
	min_pos = 0;
	for (i = shift ; i < len; ++i) {
		if(pattern[((i-shift)%pattern_len)] == '1'){

			int c = seq_nt4_table[(uint8_t)str[i]];
			mm128_t info = { UINT64_MAX, UINT64_MAX };
			if (c < 4) { // not an ambiguous base
				int z;
				//TODO:: Enable HPC
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
				kmer_span = l + 1 < k? l + 1 : k;
				kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
				kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
				if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
				z = kmer[0] < kmer[1]? 0 : 1; // strand
				++l;
				if (l >= k && kmer_span < 256) {
					info.x = hash64(kmer[z], mask) << 8 | kmer_span;

					ones_index=(ones_loc_index) % ones;
					old_location= offset + ((ones_index==0)?ones_loc[ones-1] : ones_loc[ones_index-1]);

					//printf(" new i: %d, old i: %d\n",i-k+1, old_location-k+1);
					//printf(" new i: %d, old i: %d\n",i, old_location);
					//info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
					info.y = (uint64_t)rid<<32 | (uint32_t)old_location<<1 | z;
					ones_loc_index++;
					offset = offset + ((ones_index==0)?pattern_len:0);
				}
			} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
			buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below

			if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
				for (j = buf_pos + 1; j < w; ++j){
					if (min.x == buf[j].x && buf[j].y != min.y){
//						printbits(buf[j].x, 64);
//						printf("\tLocation: ");
//						printf("%d",(buf[j].y&0xff)>>1);
//						printf("\tStrand: %c","+-"[buf[j].y&0x00000001]);
//						printf("\n");
						kv_push(mm128_t, km, *p, buf[j]);
					}
				}
				for (j = 0; j < buf_pos; ++j){
					if (min.x == buf[j].x && buf[j].y != min.y){
//						printbits(buf[j].x, 64);
//						printf("\tLocation: ");
//						printf("%d",(buf[j].y&0xff)>>1);
//						printf("\tStrand: %c","+-"[buf[j].y&0x00000001]);
//						printf("\n");
						kv_push(mm128_t, km, *p, buf[j]);
					}
				}
			}
			if (info.x <= min.x) { // a new minimum; then write the old min
				if (l >= w + k && min.x != UINT64_MAX){
//					printbits(min.x, 64);
//					printf("\tLocation: ");
//					printf("%d",(min.y&0xff)>>1);
//					printf("\tStrand: %c","+-"[min.y&0x00000001]);
//					printf("\n");
					kv_push(mm128_t, km, *p, min);
				}
				min = info, min_pos = buf_pos;
			} else if (buf_pos == min_pos) { // old min has moved outside the window
				if (l >= w + k - 1 && min.x != UINT64_MAX){
//					printbits(min.x, 64);
//					printf("\tLocation: ");
//					printf("%d",(min.y&0xff)>>1);
//					printf("\tStrand: %c","+-"[min.y&0x00000001]);
//					printf("\n");
					kv_push(mm128_t, km, *p, min);
				}
				for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
					if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
				for (j = 0; j <= buf_pos; ++j)
					if (min.x >= buf[j].x) min = buf[j], min_pos = j;
				if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
					for (j = buf_pos + 1; j < w; ++j){ // these two loops make sure the output is sorted
						if (min.x == buf[j].x && min.y != buf[j].y){
//							printbits(buf[j].x, 64);
//							printf("\tLocation: ");
//							printf("%d",(buf[j].y&0xff)>>1);
//							printf("\tStrand: %c","+-"[buf[j].y&0x00000001]);
//							printf("\n");
							kv_push(mm128_t, km, *p, buf[j]);
						}
					}
					for (j = 0; j <= buf_pos; ++j){
						if (min.x == buf[j].x && min.y != buf[j].y){
//							printbits(buf[j].x, 64);
//							printf("\tLocation: ");
//							printf("%d",(buf[j].y&0xff)>>1);
//							printf("\tStrand: %c","+-"[buf[j].y&0x00000001]);
//							printf("\n");
							kv_push(mm128_t, km, *p, buf[j]);
						}
					}
				}
			}
			if (++buf_pos == w) buf_pos = 0;
		}
	}
	if (min.x != UINT64_MAX){
//		printbits(min.x, 64);
//		printf("\tLocation: ");
//		printf("%d",(min.y&0xff)>>1);
//		printf("\tStrand: %c","+-"[min.y&0x00000001]);
//		printf("\n");
		kv_push(mm128_t, km, *p, min);
	}
	if (mm_dbg_flag & MM_DBG_PRINT_SEED) free(new_genome);
	free(ones_loc);
}
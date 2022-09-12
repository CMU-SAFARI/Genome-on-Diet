#include <string.h>
#include <assert.h>
#include "ksw2.h"

#ifdef __SSE2__
#include <emmintrin.h>

#ifdef KSW_SSE2_ONLY
#undef __SSE4_1__
#endif

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

#ifdef KSW_CPU_DISPATCH
#ifdef __SSE4_1__
void exact_match_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez, bool *exact_match, int *mismatch_cnt)
#else
void exact_match_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez, bool *exact_match, int *mismatch_cnt)
#endif
#else
void exact_match_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez, bool *exact_match, int *mismatch_cnt)
#endif // ~KSW_CPU_DISPATCH
{
									 
    if (m <= 0 || qlen <= 0 || tlen <= 0) return;


    *exact_match = true;

		__m128i sq_new, st_new, tmp, mask;	
		__m128i neq;
        __m128i vcmp;
        uint16_t vmask;

		for (int t = 0; t + 16 <= qlen; t += 16) {
			sq_new = _mm_loadu_si128((__m128i*)&target[t]);
			st_new = _mm_loadu_si128((__m128i*)&query[t]);


            vcmp = _mm_cmpeq_epi8(sq_new, st_new);       // PCMPEQB
            vmask = _mm_movemask_epi8(vcmp);    // PMOVMSKB
            if (vmask == 0xffff)
                *exact_match &= true;
            else
                *exact_match &= false;

           // *mismatch_cnt += (16 - (__builtin_popcount(vmask)));

           // printf("mismatch: %d \n", *mismatch_cnt);

         //  printf("stage exatch match: %s", *exact_match ? "true \n": "false \n");
          if(! *exact_match)
               break;
		}

        if(*exact_match){

            int mul = qlen / 16;
            int res = qlen % 16;

            uint8_t target_tail[16], query_tail[16];
            for(int i = 0; i < 16; i++){

                if(i < res){
                    target_tail[i] = target[16 * mul + i];
                    query_tail[i] = query[16 * mul + i];
                }
                else {
                    target_tail[i] = 0;
                    query_tail[i] = 0;
                }
            }

            sq_new = _mm_loadu_si128((__m128i*)&target_tail[0]);
            st_new = _mm_loadu_si128((__m128i*)&query_tail[0]);
            vcmp = _mm_cmpeq_epi8(sq_new, st_new);       // PCMPEQB
            vmask = _mm_movemask_epi8(vcmp);    // PMOVMSKB
            if (vmask == 0xffff)
                *exact_match &= true;
            else
                *exact_match &= false;

        }

           // *mismatch_cnt += (16 - (__builtin_popcount(vmask)));
            *mismatch_cnt = 0;

     //  printf("%s", *exact_match ? "exact match!! \n" : "not exact \n");

}
#endif // __SSE2__

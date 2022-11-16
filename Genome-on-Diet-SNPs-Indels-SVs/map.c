#include "kalloc.h"
#include "khash.h"
#include "kthread.h"
#include "kvec.h"
#include "mmpriv.h"
#include "profile.h"
#include "sdust.h"
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ksw2.h"
#include "ksw2_extd2_avx.h"
#include <inttypes.h>

#define mm_seq4_set(s, i, c) ((s)[(i) >> 3] |= (uint32_t)(c) << (((i)&7) << 2))
#define mm_seq4_get(s, i) ((s)[(i) >> 3] >> (((i)&7) << 2) & 0xf)

struct mm_tbuf_s {
	void *km;
	int rep_len, frag_gap;
};

mm_tbuf_t *mm_tbuf_init(void) {
	mm_tbuf_t *b;
	b = (mm_tbuf_t *)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b) {
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

void *mm_tbuf_get_km(mm_tbuf_t *b) { return b->km; }

int concatenate_cigars(mm_reg1_t *const rstart, const mm_reg1_t rend, const uint8_t *const qseq, const uint8_t str,
                       const uint32_t read_len, const mm_idx_t *const mi, void *km, const uint32_t sc_mch,
                       const uint32_t sc_mis, const uint32_t gapo1, const uint32_t gape1, const uint32_t gapo2,
                       const uint32_t gape2) {
	const uint32_t tstart = rstart->rs;
	const uint32_t tend   = rend.re;

	const uint32_t tstart_junc = rend.rs;
	const uint32_t tend_junc   = rstart->re;

	const uint32_t qstart = str ? read_len - rstart->qe : rstart->qs;
	const uint32_t qend   = str ? read_len - rend.qs : rend.qe;

	const uint32_t qstart_junc = str ? read_len - rend.qe : rend.qs;
	const uint32_t qend_junc   = str ? read_len - rstart->qs : rstart->qe;

	// Check if they cross either on query or target
	if (tend_junc <= tstart_junc && qend_junc <= qstart_junc) {
		return 1;
	}

	// Check if one is not within the other
	if (tend_junc >= tend || tstart >= tstart_junc) {
		return 1;
	}

	if (qend_junc >= qend || qstart >= qstart_junc) {
		return 1;
	}

	const uint32_t size_start = rstart->re - rstart->rs;
	const uint32_t size_end   = rend.re - rend.rs;
	const uint32_t tseq_size  = size_start > size_end ? size_start : size_end;
	uint8_t *tseq             = (uint8_t *)kmalloc(km, tseq_size * sizeof(uint8_t));
	unsigned juncq;
	unsigned junct;
	unsigned cigar_pos;
	int score;
	/*
	fprintf(stderr, "First part:\n");
	for (unsigned i = 0; i < rstart->p->n_cigar; i++) {
	        fprintf(stderr, "%d%c", rstart->p->cigar[i] >> 4, "MIDNSH"[rstart->p->cigar[i] & 0xf]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "Second part:\n");
	for (unsigned i = 0; i < rend.p->n_cigar; i++) {
	        fprintf(stderr, "%d%c", rend.p->cigar[i] >> 4, "MIDNSH"[rend.p->cigar[i] & 0xf]);
	}
	fprintf(stderr, "\n");
	*/

	if (qend_junc > qstart_junc) {
		mm_idx_getseq2(mi, 0, rstart->rid, tstart, tend_junc, tseq);
		const uint32_t juncture_len = qend_junc - qstart_junc;
		int *const al_start_a       = (int *)kmalloc(km, juncture_len * sizeof(int));
		int *const al_end_a         = (int *)kmalloc(km, juncture_len * sizeof(int));

		// Populate al_start_a array
		int al_score       = 0;
		uint32_t toff      = 0;
		uint32_t qoffstart = qstart;
		for (uint32_t i = 0; i < rstart->p->n_cigar; i++) {
			uint32_t op = rstart->p->cigar[i] & 0xf, len = rstart->p->cigar[i] >> 4;
			// fprintf(stderr, "op %u\n", op);
			if (op == MM_CIGAR_MATCH) {
				/*
				fprintf(stderr, "Read:\n");
				fprintf(stderr, "%u= %u\n", len, qoffstart);
				for (unsigned b = 0; b < len; b++) {
				        fprintf(stderr, "%c", "ACGTN"[qseq[qoffstart + b]]);
				}
				fprintf(stderr, "\nReff:\n");
				for (unsigned b = 0; b < len; b++) {
				        fprintf(stderr, "%c", "ACGTN"[tseq[toff + b]]);
				}
				fputs("\n", stderr);
				*/
				for (unsigned j = 0; j < len; j++) {
					if (qoffstart + j >= qstart_junc) {
						al_start_a[qoffstart + j - qstart_junc] = al_score;
					}
					if (qseq[qoffstart + j] == tseq[toff + j]) {
						al_score += sc_mch;
					} else {
						al_score -= sc_mis;
					}
				}
				qoffstart += len;
				toff += len;
			} else if (op == MM_CIGAR_INS) {
				const uint32_t p1 = gapo1 + len * gape1;
				const uint32_t p2 = gapo2 + len * gape2;
				if (qoffstart + len <= qstart_junc) {
					al_score -= (p1 < p2) ? p1 : p2;
				} else if (qoffstart < qstart_junc) {
					unsigned o;
					unsigned e;
					if (p1 < p2) {
						o = gapo1;
						e = gape1;
					} else {
						o = gapo2;
						e = gape2;
					}
					al_score -= o + e * (qstart_junc - qoffstart);
					// fprintf(stderr, "%u\n", len - soff);
					for (unsigned j = 0; j < qoffstart + len - qstart_junc; j++) {
						al_start_a[j] = al_score;
						al_score -= e;
					}
				} else {
					unsigned o;
					unsigned e;
					if (p1 < p2) {
						o = gapo1;
						e = gape1;
					} else {
						o = gapo2;
						e = gape2;
					}
					al_start_a[qoffstart - qstart_junc] = al_score;
					al_score -= o + e;
					for (unsigned j = 1; j < len; j++) {
						al_start_a[qoffstart + j - qstart_junc] = al_score;
						al_score -= e;
					}
				}
				qoffstart += len;
			} else if (op == MM_CIGAR_DEL) {
				const uint32_t p1 = gapo1 + len * gape1;
				const uint32_t p2 = gapo2 + len * gape2;
				al_score -= (p1 < p2) ? p1 : p2;
				toff += len;
			} else if (op == MM_CIGAR_N_SKIP) {
				toff += len;
			}
		}

		/*
		fputs("\n", stderr);
		for (unsigned i = 0; i < juncture_len; i++) {
		        fprintf(stderr, "<%d>", al_start_a[i]);
		}
		fputs("\n", stderr);
		*/

		mm_idx_getseq2(mi, 0, rend.rid, tstart_junc, tend, tseq);
		// Populate al_end array
		toff             = 0;
		uint32_t qoffend = qstart_junc;
		al_score         = rend.score;
		for (uint32_t i = 0; i < rend.p->n_cigar && qoffend <= qend_junc; i++) {
			uint32_t op = rend.p->cigar[i] & 0xf, len = rend.p->cigar[i] >> 4;
			// fprintf(stderr, "op %u\n", op);
			if (op == MM_CIGAR_MATCH) {
				/*
				fprintf(stderr, "Read:\n");
				fprintf(stderr, "%u= %u\n", len, qoffend);
				for (unsigned b = 0; b < len; b++) {
				        fprintf(stderr, "%c", "ACGTN"[qseq[qoffend + b]]);
				}
				fprintf(stderr, "\nReff:\n");
				for (unsigned b = 0; b < len; b++) {
				        fprintf(stderr, "%c", "ACGTN"[tseqend[toff + b]]);
				}
				fputs("\n", stderr);
				*/
				for (unsigned j = 0; j < len; j++) {
					if (qoffend + j < qend_junc) {
						if (qseq[qoffend + j] == tseq[toff + j]) {
							al_score -= sc_mch;
						} else {
							al_score += sc_mis;
						}
						al_end_a[qoffend + j - qstart_junc] = al_score;
					} else {
						break;
					}
				}
				qoffend += len;
				toff += len;
			} else if (op == MM_CIGAR_INS) {
				const uint32_t p1 = gapo1 + len * gape1;
				const uint32_t p2 = gapo2 + len * gape2;
				unsigned o;
				unsigned e;
				if (p1 < p2) {
					o = gapo1;
					e = gape1;
				} else {
					o = gapo2;
					e = gape2;
				}
				al_score += o;
				for (unsigned j = 0; j < len; j++) {
					if (qoffend + j < qend_junc) {
						al_score += e;
						al_end_a[qoffend + j - qstart_junc] = al_score;
					} else {
						break;
					}
				}
				qoffend += len;
			} else if (op == MM_CIGAR_DEL) {
				const uint32_t p1 = gapo1 + len * gape1;
				const uint32_t p2 = gapo2 + len * gape2;
				al_score += (p1 < p2) ? p1 : p2;
				toff += len;
			} else if (op == MM_CIGAR_N_SKIP) {
				toff += len;
			}
		}
		kfree(km, tseq);

		/*
		fputs("\n", stderr);
		for (unsigned i = 0; i < juncture_len; i++) {
		        fprintf(stderr, "<%d>", al_end_a[i]);
		}
		fputs("\n", stderr);
		*/

		// Find the best junction point to maximise the alignment score
		int max_score = al_start_a[0] + al_end_a[0];
		juncq         = 0;
		for (unsigned start = 1; start < juncture_len; start++) {
			int total_score = al_start_a[start] + al_start_a[start];
			if (total_score > max_score) {
				max_score = total_score;
				juncq     = start;
			}
		}
		score = max_score;
		// fprintf(stderr, "juncture: %u, max_nb_match %u\n", juncture, max_score);
		juncq += qstart_junc;
		kfree(km, al_start_a);
		kfree(km, al_end_a);

		// Alocate new CIGAR
		if (rstart->p->n_cigar + rend.p->n_cigar + 1 + sizeof(mm_extra_t) / 4 > rstart->p->capacity) {
			rstart->p->capacity = rstart->p->n_cigar + rend.p->n_cigar + 1 + sizeof(mm_extra_t) / 4;
			kroundup32(rstart->p->capacity);
			rstart->p = (mm_extra_t *)realloc(rstart->p, rstart->p->capacity * 4);
		}

		qoffstart = qstart;
		uint32_t i;
		uint32_t toffs = rstart->rs;
		for (i = 0; i < rstart->p->n_cigar; i++) {
			uint32_t op = rstart->p->cigar[i] & 0xf, len = rstart->p->cigar[i] >> 4;
			if (op == MM_CIGAR_MATCH) {
				if (qoffstart + len >= juncq) {
					const uint32_t new_len = juncq - qoffstart;
					// fprintf(stderr, "len: %u, new_len %u\n", len, new_len);
					rstart->p->cigar[i] = MM_CIGAR_MATCH | (new_len << 4);
					qoffstart += new_len;
					toffs += new_len;
					i++;
					break;
				}
				qoffstart += len;
				toffs += len;
			} else if (op == MM_CIGAR_INS) {
				if (qoffstart + len >= juncq) {
					// move the juncture
					juncq = qoffstart;
					break;
				}
				qoffstart += len;
			} else if (op == MM_CIGAR_DEL) {
				toffs += len;
			} else if (op == MM_CIGAR_N_SKIP) {
				toffs += len;
			}
		}
		junct     = toffs;
		cigar_pos = i;
	} else {
		const uint32_t juncture_len = tend_junc - tstart_junc;
		int *const al_start_a       = (int *)kmalloc(km, juncture_len * sizeof(int));
		int *const al_end_a         = (int *)kmalloc(km, juncture_len * sizeof(int));
		mm_idx_getseq2(mi, 0, rstart->rid, tstart, tend_junc, tseq);

		uint32_t toff          = 0;
		uint32_t qoffstart     = qstart;
		int al_score           = 0;
		const uint32_t sofft_s = tstart_junc - tstart;
		for (uint32_t i = 0; i < rstart->p->n_cigar; i++) {
			uint32_t op = rstart->p->cigar[i] & 0xf, len = rstart->p->cigar[i] >> 4;
			// fprintf(stderr, "op %u\n", op);
			if (op == MM_CIGAR_MATCH) {
				/*
				fprintf(stderr, "Read:\n");
				fprintf(stderr, "%u= %u\n", len, qoffstart);
				for (unsigned b = 0; b < len; b++) {
				        fprintf(stderr, "%c", "ACGTN"[qseq[qoffstart + b]]);
				}
				fprintf(stderr, "\nReff:\n");
				for (unsigned b = 0; b < len; b++) {
				        fprintf(stderr, "%c", "ACGTN"[tseq[toff + b]]);
				}
				fputs("\n", stderr);
				*/
				for (unsigned j = 0; j < len; j++) {
					if (toff + j >= sofft_s) {
						al_start_a[toff + j - sofft_s] = al_score;
					}
					if (qseq[qoffstart + j] == tseq[toff + j]) {
						al_score += sc_mch;
					} else {
						al_score -= sc_mis;
					}
				}
				qoffstart += len;
				toff += len;
			} else if (op == MM_CIGAR_DEL) {
				const uint32_t p1 = gapo1 + len * gape1;
				const uint32_t p2 = gapo2 + len * gape2;
				if (toff + len <= sofft_s) {
					al_score -= (p1 < p2) ? p1 : p2;
				} else if (toff < sofft_s) {
					unsigned o;
					unsigned e;
					if (p1 < p2) {
						o = gapo1;
						e = gape1;
					} else {
						o = gapo2;
						e = gape2;
					}
					al_score -= o + e * (sofft_s - toff);
					// fprintf(stderr, "%u\n", len - soff);
					for (unsigned j = 0; j < toff + len - sofft_s; j++) {
						al_start_a[j] = al_score;
						al_score -= e;
					}
				} else {
					unsigned o;
					unsigned e;
					if (p1 < p2) {
						o = gapo1;
						e = gape1;
					} else {
						o = gapo2;
						e = gape2;
					}
					al_start_a[toff - sofft_s] = al_score;
					al_score -= o + e;
					for (unsigned j = 1; j < len; j++) {
						al_start_a[toff + j - sofft_s] = al_score;
						al_score -= e;
					}
				}
				toff += len;
			} else if (op == MM_CIGAR_INS) {
				const uint32_t p1 = gapo1 + len * gape1;
				const uint32_t p2 = gapo2 + len * gape2;
				al_score -= (p1 < p2) ? p1 : p2;
				qoffstart += len;
			} else if (op == MM_CIGAR_N_SKIP) {
				toff += len;
			}
		}

		/*
		fputs("\n", stderr);
		for (unsigned i = 0; i < juncture_len; i++) {
		        fprintf(stderr, "<%d>", al_start_a[i]);
		}
		fputs("\n", stderr);
		*/

		mm_idx_getseq2(mi, 0, rend.rid, rend.rs, rend.re, tseq);
		toff                   = 0;
		uint32_t qoffend       = qstart_junc;
		al_score               = 0;
		const uint32_t eofft_s = tend_junc - tstart_junc;

		for (uint32_t i = 0; i < rend.p->n_cigar && toff <= eofft_s; i++) {
			uint32_t op = rend.p->cigar[i] & 0xf, len = rend.p->cigar[i] >> 4;
			// fprintf(stderr, "op %u\n", op);
			if (op == MM_CIGAR_MATCH) {
				/*
				fprintf(stderr, "Read:\n");
				fprintf(stderr, "%u= %u\n", len, qoffend);
				for (unsigned b = 0; b < len; b++) {
				        fprintf(stderr, "%c", "ACGTN"[qseq[qoffend + b]]);
				}
				fprintf(stderr, "\nReff:\n");
				for (unsigned b = 0; b < len; b++) {
				        fprintf(stderr, "%c", "ACGTN"[tseqend[toff + b]]);
				}
				fputs("\n", stderr);
				*/
				for (unsigned j = 0; j < len; j++) {
					if (toff + j < eofft_s) {
						if (qseq[qoffend + j] == tseq[toff + j]) {
							al_score -= sc_mch;
						} else {
							al_score += sc_mis;
						}
						al_end_a[toff + j] = al_score;
					} else {
						break;
					}
				}
				qoffend += len;
				toff += len;
			} else if (op == MM_CIGAR_DEL) {
				const uint32_t p1 = gapo1 + len * gape1;
				const uint32_t p2 = gapo2 + len * gape2;
				unsigned o;
				unsigned e;
				if (p1 < p2) {
					o = gapo1;
					e = gape1;
				} else {
					o = gapo2;
					e = gape2;
				}
				al_score += o;
				for (unsigned j = 0; j < len; j++) {
					if (toff + j < eofft_s) {
						al_score += e;
						al_end_a[toff + j] = al_score;
					} else {
						break;
					}
				}
				toff += len;
			} else if (op == MM_CIGAR_INS) {
				const uint32_t p1 = gapo1 + len * gape1;
				const uint32_t p2 = gapo2 + len * gape2;
				al_score += (p1 < p2) ? p1 : p2;
				qoffend += len;
			} else if (op == MM_CIGAR_N_SKIP) {
				toff += len;
			}
		}
		/*
		fputs("\n", stderr);
		for (unsigned i = 0; i < juncture_len; i++) {
		        fprintf(stderr, "<%d>", al_end_a[i]);
		}
		fputs("\n", stderr);
		*/
		kfree(km, tseq);

		// Find the best junction point to maximise the alignment score
		int max_score = al_start_a[0] + al_end_a[0];
		junct         = 0;
		for (unsigned start = 1; start < juncture_len; start++) {
			int total_score = al_start_a[start] + al_start_a[start];
			if (total_score > max_score) {
				max_score = total_score;
				junct     = start;
			}
		}
		// fprintf(stderr, "juncture_t: %u, max_nb_match %u\n", junct, max_score);
		score = max_score;

		junct += tstart_junc;
		kfree(km, al_start_a);
		kfree(km, al_end_a);

		// Alocate new CIGAR
		if (rstart->p->n_cigar + rend.p->n_cigar + 1 + sizeof(mm_extra_t) / 4 > rstart->p->capacity) {
			rstart->p->capacity = rstart->p->n_cigar + rend.p->n_cigar + 1 + sizeof(mm_extra_t) / 4;
			kroundup32(rstart->p->capacity);
			rstart->p = (mm_extra_t *)realloc(rstart->p, rstart->p->capacity * 4);
		}

		qoffstart = qstart;
		uint32_t i;
		uint32_t toffs = rstart->rs;
		for (i = 0; i < rstart->p->n_cigar; i++) {
			uint32_t op = rstart->p->cigar[i] & 0xf, len = rstart->p->cigar[i] >> 4;
			if (op == MM_CIGAR_MATCH) {
				if (toffs + len >= junct) {
					const uint32_t new_len = junct - toffs;
					// fprintf(stderr, "len: %u, new_len %u\n", len, new_len);
					rstart->p->cigar[i] = MM_CIGAR_MATCH | (new_len << 4);
					qoffstart += new_len;
					toffs += new_len;
					i++;
					break;
				}
				qoffstart += len;
				toffs += len;
			} else if (op == MM_CIGAR_DEL) {
				if (toffs + len >= junct) {
					// move the juncture
					junct = toffs;
					break;
				}
				toffs += len;
			} else if (op == MM_CIGAR_INS) {
				qoffstart += len;
			} else if (op == MM_CIGAR_N_SKIP) {
				toffs += len;
			}
		}
		juncq     = qoffstart;
		cigar_pos = i;
	}

	// fprintf(stderr, "juncq: %u junct %u\n", juncq, junct);

	uint32_t toffe   = rend.rs;
	uint32_t qoffend = qstart_junc;
	unsigned i       = cigar_pos;
	int crossed      = 0;
	for (uint32_t j = 0; j < rend.p->n_cigar; j++) {
		uint32_t op = rend.p->cigar[j] & 0xf, len = rend.p->cigar[j] >> 4;
		// fprintf(stderr, "c: %u%c\n", len, "MIDNSH"[op]);
		if (op == MM_CIGAR_MATCH) {
			if (crossed) {
				rstart->p->cigar[i] = rend.p->cigar[j];
				i++;
			}
			qoffend += len;
			toffe += len;
		} else if (op == MM_CIGAR_INS) {
			if (crossed) {
				rstart->p->cigar[i] = rend.p->cigar[j];
				i++;
			}
			qoffend += len;
		} else if (op == MM_CIGAR_DEL) {
			if (crossed) {
				rstart->p->cigar[i] = rend.p->cigar[j];
				i++;
			}
			toffe += len;
		} else if (op == MM_CIGAR_N_SKIP) {
			if (crossed) {
				rstart->p->cigar[i] = rend.p->cigar[j];
				i++;
			}
			toffe += len;
		}

		if (crossed == 0 && qoffend >= juncq && toffe >= junct) {
			const uint32_t tar_len = toffe - junct;
			const uint32_t que_len = qoffend - juncq;
			if (que_len > tar_len) {
				// need to insert
				const uint32_t len = que_len - tar_len;
				const uint32_t p1  = gapo1 + len * gape1;
				const uint32_t p2  = gapo2 + len * gape2;
				score -= (p1 < p2) ? p1 : p2;
				rstart->p->cigar[i] = MM_CIGAR_INS | (len << 4);
				i++;
				if (tar_len != 0) {
					rstart->p->cigar[i] = MM_CIGAR_MATCH | (tar_len << 4);
					i++;
				}
			} else if (que_len < tar_len) {
				// need to delete
				const uint32_t len = tar_len - que_len;
				const uint32_t p1  = gapo1 + len * gape1;
				const uint32_t p2  = gapo2 + len * gape2;
				score -= (p1 < p2) ? p1 : p2;
				rstart->p->cigar[i] = MM_CIGAR_DEL | (len << 4);
				i++;
				if (que_len != 0) {
					rstart->p->cigar[i] = MM_CIGAR_MATCH | (que_len << 4);
					i++;
				}
			} else {
				// only match
				rstart->p->cigar[i] = MM_CIGAR_MATCH | ((tar_len) << 4);
				i++;
			}
			crossed = 1;
		}
	}
	rstart->p->n_cigar  = i;
	rstart->p->dp_score = score;
	rstart->score       = score;

	/*
	for (unsigned i = 0; i < rstart->p->n_cigar; i++) {
	        fprintf(stderr, "%d%c", rstart->p->cigar[i] >> 4, "MIDNSH"[rstart->p->cigar[i] & 0xf]);
	        if (rstart->p->cigar[i] >> 4 > 50000) {
	                exit(0);
	        }
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "start.qs: %u start.qe: %u, end.qs: %u, end.qe: %u\n", rstart->qs, rstart->qe, rend.qs,
	        rend.qe);
	        */
	if (str) {
		rstart->qs = rend.qs;
	} else {
		rstart->qe = rend.qe;
	}
	rstart->re = rend.re;
	return 0;
}

void mm_update_extra(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e,
                     int is_eqx, int log_gap);

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres) {
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb  = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t *)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y >> 1, span = a[j].x & 0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (int32_t)dreg[u] <= s)
			++u;
		if (u < n_dreg && (int32_t)(dreg[u] >> 32) < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && (int32_t)(dreg[v] >> 32) < e;
			     ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > (int32_t)(dreg[v] >> 32) ? s : dreg[v] >> 32;
				int ee = e < (int32_t)dreg[v] ? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span >> 1)
				a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else
			a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static inline unsigned collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, const int qlen,
                                          const char *seq, mm128_v *mv, const char *Z, int W, int shift,
                                          uint32_t MAX_NB_SEEDS) {
	int n = 0;
	mv->n = 0;
	// TODO: HPC should be disabled
	unsigned extracted_len =
	    mm_sketch3(km, seq, qlen, mi->w, mi->k, 0, mi->flag & MM_I_HPC, mv, Z, W, shift, MAX_NB_SEEDS);
	/*
	printf("nb_seeds %lu\n", mv->n);
	printf("extracted_len %u\n", extracted_len);
	printf("len %u\n", qlen);
	for (unsigned j = 0; j < mv->n; j++) {
	        printf(".x: %lx\t .y: %lu\n", mv->a[j].x, mv->a[j].y >> 1 & UINT32_MAX);
	}
	*/
	if (opt->sdust_thres > 0) // mask low-complexity minimizers
		mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlen, seq, opt->sdust_thres);
	return extracted_len;
}

static inline mm_pattern_t collect_minimizers2(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, const int qlen,
                                               const char *seq, mm128_v *mv, const char *Z, int W,
                                               const float max_seeds) {
	mv->n = 0;
	return mm_sketch2(km, seq, qlen, mi->w, mi->k, 0, mi->flag & MM_I_HPC, mv, Z, W, max_seeds);
	// fprintf(stderr, "Collect_minimizers2:\n");
	// printf("shift_seeds_number: %d, size: %d\n",mm_pattern_versions1->shift_seeds_number[0],
	// mm_pattern_versions1->n_shift_seeds_number);
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

static inline int skip_seed(int flag, uint64_t r, const mm_seed_t *q, const char *qname, int qlen, const mm_idx_t *mi,
                            int *is_self) {
	*is_self = 0;
	if (qname && (flag & (MM_F_NO_DIAG | MM_F_NO_DUAL))) {
		const mm_idx_seq_t *s = &mi->seq[r >> 32];
		int cmp;
		cmp = strcmp(qname, s->name);
		if ((flag & MM_F_NO_DIAG) && cmp == 0 && (int)s->len == qlen) {
			if ((uint32_t)r >> 1 == (q->q_pos >> 1)) return 1; // avoid the diagnonal anchors
			if ((r & 1) == (q->q_pos & 1))
				*is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag & MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY | MM_F_REV_ONLY)) {
		if ((r & 1) == (q->q_pos & 1)) { // forward strand
			if (flag & MM_F_REV_ONLY) return 1;
		} else {
			if (flag & MM_F_FOR_ONLY) return 1;
		}
	}
	return 0;
}

// loc_t:
// loc: |- 0 -|--- CHROM_ID ---|--- LOC ---|
//          63               32           0
// target can be < 0
typedef struct {
	uint64_t target;
	uint32_t query;
} loc_t;

void heap_sort(void *km, loc_t **const queries, loc_t **const buf, unsigned *const pos, unsigned nb_units) {
	if (nb_units <= 1) {
		return;
	}

	loc_t *src = *queries;
	loc_t *tgt = *buf;

	// Create the heap
	unsigned heap_size = nb_units;
	mm128_t *heap      = (mm128_t *)kmalloc(km, heap_size * sizeof(mm128_t));

	heap[0].x = src[0].target;
	heap[0].y = 0;
	for (unsigned i = 1; i < heap_size; i++) {
		heap[i].x = src[pos[i - 1]].target;
		heap[i].y = (uint64_t)i << 32;
	}

	ks_heapmake_heap(heap_size, heap);
	for (unsigned i = 0; heap_size > 0; i++) {
		unsigned unit     = heap[0].y >> 32;
		unsigned base_pos = (unit == 0) ? 0 : pos[unit - 1];
		tgt[i]            = (loc_t){.target = heap[0].x, .query = src[base_pos + (uint32_t)heap[0].y].query};
		if ((uint32_t)heap[0].y + base_pos < pos[unit] - 1) {
			heap[0].y++;
			heap[0].x = src[base_pos + (uint32_t)heap[0].y].target;
		} else {
			heap[0] = heap[heap_size - 1];
			heap_size--;
		}
		ks_heapdown_heap(0, heap_size, heap);
	}
	kfree(km, heap);

	*buf     = src;
	*queries = tgt;
}

// Branchless merge
static inline void merge_locations(const loc_t *const src, loc_t *const tgt, const unsigned l_1, const unsigned l_2) {
	unsigned i_1 = 0, i_2 = 0, i_tgt = 0;

	const loc_t *const src_1 = src;
	const loc_t *const src_2 = src + l_1;

	while (i_1 < l_1 && i_2 < l_2) {
		loc_t in_1 = src_1[i_1];
		loc_t in_2 = src_2[i_2];

		int flag          = (in_1.target - in_2.target) >> 63;
		tgt[i_tgt].target = flag * in_1.target + (1 - flag) * in_2.target;
		tgt[i_tgt].query  = flag * in_1.query + (1 - flag) * in_2.query;

		i_tgt++;
		i_1 = i_1 + flag;
		i_2 = i_2 + (1 - flag);
	}

	if (i_1 < l_1) {
		unsigned len = l_1 - i_1;
		memcpy(&tgt[i_tgt], &src_1[i_1], len * sizeof(loc_t));
	} else {
		unsigned len = l_2 - i_2;
		memcpy(&tgt[i_tgt], &src_2[i_2], len * sizeof(loc_t));
	}
}

static inline void merge_sort(loc_t **const queries, loc_t **const buf, unsigned *const pos, unsigned nb_units) {
	if (nb_units <= 1) {
		if (nb_units == 0) {
			pos[0] = 0;
		}
		return;
	}

	loc_t *src = *queries;
	loc_t *tgt = *buf;
	for (;;) {
		for (unsigned i = 0; i < nb_units / 2; i++) {
			const unsigned offset   = (i == 0) ? 0 : pos[2 * i - 1];
			const loc_t *const _src = src + offset;
			loc_t *const _tgt       = tgt + offset;
			const unsigned l1       = pos[2 * i] - offset;
			const unsigned l2       = pos[2 * i + 1] - pos[2 * i];
			merge_locations(_src, _tgt, l1, l2);

			pos[i] = pos[2 * i + 1];
		}

		// If there is an odd number of units
		if (nb_units % 2) {
			const unsigned offset = pos[nb_units - 2];
			const unsigned len    = pos[nb_units - 1] - pos[nb_units - 2];

			const loc_t *const _src = src + offset;
			loc_t *const _tgt       = tgt + offset;

			memcpy(_tgt, _src, len * sizeof(loc_t));

			pos[nb_units / 2] = pos[nb_units - 1];
		}

		nb_units = nb_units / 2 + nb_units % 2;

		if (nb_units == 1) {
			*queries = tgt;
			*buf     = src;
			return;
		}

		// Swap the buffers
		loc_t *tmp = src;
		src        = tgt;
		tgt        = tmp;
	}
}

static void collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname,
                              const mm128_v *mv, int qlen, unsigned tmp_extracted_len, loc_t **a_for, loc_t **a_rev,
                              unsigned *n_a_for, unsigned *n_a_rev, int heap_sort_en) {
	int i, n_m;
	mm_seed_t *m;
	int64_t n_a;
	m = mm_collect_matches2(km, &n_m, qlen, max_occ, opt->max_max_occ, opt->occ_dist, mi, mv, &n_a);

	*a_for = (loc_t *)kmalloc(km, n_a * sizeof(loc_t));
	*a_rev = (loc_t *)kmalloc(km, n_a * sizeof(loc_t));

	loc_t *loc_for = *a_for;
	loc_t *loc_rev = *a_rev;
	loc_t *buf     = (loc_t *)kmalloc(km, n_a * sizeof(loc_t));

	unsigned len_for      = 0;
	unsigned len_rev      = 0;
	unsigned nb_seeds_for = 0;
	unsigned nb_seeds_rev = 0;

	unsigned *pos_for = (unsigned *)kmalloc(km, n_m * sizeof(unsigned));
	unsigned *pos_rev = (unsigned *)kmalloc(km, n_m * sizeof(unsigned));

	for (i = 0; i < n_m; ++i) {
		mm_seed_t *q           = &m[i];
		uint32_t nb_locs       = q->n;
		const uint64_t *r      = q->cr;
		uint32_t query_len_for = 0;
		uint32_t query_len_rev = 0;

		loc_t *buf_for = &loc_for[len_for];
		loc_t *buf_rev = &loc_rev[len_rev];

		for (uint32_t k = 0; k < nb_locs; ++k) {
			int32_t is_self;
			if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self)) continue;
			uint32_t qpos     = q->q_pos >> 1;
			unsigned str      = (r[k] & 1) ^ (q->q_pos & 1);
			uint32_t loc      = (uint32_t)r[k] >> 1;
			uint64_t chrom_id = r[k] >> 32;
			if (str) { // reverse strand
				loc                           = loc + qpos - q->q_span + 1;
				buf_rev[query_len_rev].target = (chrom_id << 32) | loc;
				buf_rev[query_len_rev].query  = qpos;
				query_len_rev++;
			} else { // forward strand
				loc                           = loc + tmp_extracted_len - qpos;
				buf_for[query_len_for].target = (chrom_id << 32) | loc;
				buf_for[query_len_for].query  = qpos;
				query_len_for++;
			}
		}

		if (query_len_for > 0) {
			len_for += query_len_for;
			pos_for[nb_seeds_for] = len_for;
			nb_seeds_for++;
		}

		if (query_len_rev > 0) {
			len_rev += query_len_rev;
			pos_rev[nb_seeds_rev] = len_rev;
			nb_seeds_rev++;
		}
	}

	kfree(km, m);
	if (nb_seeds_for == 0) {
		*n_a_for = 0;
	} else {
		*n_a_for = pos_for[nb_seeds_for - 1];
		// fprintf(stderr, "for, nb_seeds: %u, nb_locs: %u\n", nb_seeds_for, *n_a_for);
		if (heap_sort_en) {
			heap_sort(km, a_for, &buf, pos_for, nb_seeds_for);
		} else {
			merge_sort(a_for, &buf, pos_for, nb_seeds_for);
		}
	}

	if (nb_seeds_rev == 0) {
		*n_a_rev = 0;
	} else {
		*n_a_rev = pos_rev[nb_seeds_rev - 1];
		// fprintf(stderr, "rev, nb_seeds: %u, nb_locs: %u\n", nb_seeds_for, *n_a_for);
		if (heap_sort_en) {
			heap_sort(km, a_rev, &buf, pos_rev, nb_seeds_rev);
		} else {
			merge_sort(a_rev, &buf, pos_rev, nb_seeds_rev);
		}
	}
	kfree(km, pos_for);
	kfree(km, pos_rev);
	kfree(km, buf);
	return;
}
#define sort_key_loc_t(a) ((a).target)
KRADIX_SORT_INIT(loc_t, loc_t, sort_key_loc_t, 8)

static void collect_seed_hits_radix(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi,
                                    const char *qname, const mm128_v *mv, int qlen, unsigned tmp_extracted_len,
                                    loc_t **a_for, loc_t **a_rev, unsigned *n_a_for, unsigned *n_a_rev) {

	int i, n_m;
	mm_seed_t *m;
	int64_t n_a;
	m = mm_collect_matches2(km, &n_m, qlen, max_occ, opt->max_max_occ, opt->occ_dist, mi, mv, &n_a);

	*a_for = (loc_t *)kmalloc(km, n_a * sizeof(loc_t));
	*a_rev = (loc_t *)kmalloc(km, n_a * sizeof(loc_t));

	loc_t *queries_for = *a_for;
	loc_t *queries_rev = *a_rev;

	unsigned len_for      = 0;
	unsigned len_rev      = 0;
	unsigned nb_seeds_for = 0;
	unsigned nb_seeds_rev = 0;

	for (i = 0; i < n_m; ++i) {
		mm_seed_t *q           = &m[i];
		uint32_t nb_locs       = q->n;
		const uint64_t *r      = q->cr;
		uint32_t query_len_for = 0;
		uint32_t query_len_rev = 0;

		loc_t *buf_for = &queries_for[len_for];
		loc_t *buf_rev = &queries_rev[len_rev];

		for (uint32_t k = 0; k < nb_locs; ++k) {
			int32_t is_self;
			if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self)) continue;
			uint32_t qpos     = q->q_pos >> 1;
			unsigned str      = (r[k] & 1) ^ (q->q_pos & 1);
			uint32_t loc      = (uint32_t)r[k] >> 1;
			uint64_t chrom_id = r[k] >> 32;
			if (str) { // reverse strand
				loc                           = loc + qpos;
				buf_rev[query_len_rev].target = (chrom_id << 32) | loc;
				buf_rev[query_len_rev].query  = qpos;
				query_len_rev++;
			} else { // forward strand
				loc                           = loc + tmp_extracted_len - qpos;
				buf_for[query_len_for].target = (chrom_id << 32) | loc;
				buf_for[query_len_for].query  = qpos;
				query_len_for++;
			}
		}

		if (query_len_for > 0) {
			len_for += query_len_for;
			nb_seeds_for++;
		}

		if (query_len_rev > 0) {
			len_rev += query_len_rev;
			nb_seeds_rev++;
		}
	}

	kfree(km, m);
	*n_a_for = len_for;
	if (nb_seeds_for > 1) {
		radix_sort_loc_t(*a_for, *a_for + len_for);
	}

	*n_a_rev = len_rev;
	if (nb_seeds_rev > 1) {
		radix_sort_loc_t(*a_rev, *a_rev + len_rev);
	}
	return;
}

typedef struct vt_t {
	uint32_t chrom_id;
	int32_t first_target_loc;
	int32_t last_target_loc;
	uint32_t first_query_loc;
	uint32_t last_query_loc;
	unsigned score;
	struct vt_t *next;
	mm_reg1_t r;
	unsigned str : 1;
	unsigned concat : 1;
	unsigned valid : 1;
} vt_t;

typedef struct {
	unsigned nb_seqs;
	vt_t *seqs;
} vt_v;

static inline void vote(const loc_t *const loc, const unsigned len, const int str, vt_v *const vt,
                        const uint32_t vt_distance, const int32_t tmp_extracted_len, const unsigned vt_max_nb_locations,
                        const uint32_t coverage_threshold) {
	if (len == 0) {
		return;
	}

	vt_t *const seqs = vt->seqs;
	unsigned out_len = vt->nb_seqs;

	unsigned counter = 1;
	const uint64_t loc_tmp =
	    str ? (loc[0].target - loc[0].query) : loc[0].target - (tmp_extracted_len - loc[0].query);
	uint64_t first_target_loc = loc_tmp;
	uint64_t last_target_loc  = loc_tmp;
	uint32_t first_query_loc  = loc[0].query;
	uint32_t last_query_loc   = loc[0].query;
	uint64_t ref_loc          = loc[0].target;
	loc_t cur;

	for (unsigned i = 1; i < len; i++) {
		cur = loc[i];

		// If curent location in the range, increase the counter and update the values
		if (cur.target - ref_loc <= vt_distance) {
			counter++;
			if (cur.query < first_query_loc) {
				first_query_loc = cur.query;
				ref_loc         = cur.target;
			}
			if (cur.query > last_query_loc) {
				last_query_loc = cur.query;
			}
			uint64_t loc = str ? (cur.target - cur.query) : cur.target - (tmp_extracted_len - cur.query);
			if (loc > last_target_loc) {
				last_target_loc = loc;
			}
			if (loc < first_target_loc) {
				first_target_loc = loc;
			}
		}
		// Else check the number of votes we have and if we are above the coverage threshold
		else {
			if (last_query_loc - first_query_loc > coverage_threshold) {
				if (out_len == vt_max_nb_locations) {
					// If not enough votes, we just continue
					if (seqs[out_len - 1].score >= counter) {
						uint64_t loc_tmp = str ? (cur.target - cur.query)
						                       : cur.target - (tmp_extracted_len - cur.query);
						first_target_loc = loc_tmp;
						last_target_loc  = loc_tmp;
						first_query_loc  = cur.query;
						last_query_loc   = cur.query;
						ref_loc          = cur.target;
						counter          = 1;
						continue;
					}
				} else {
					out_len++;
				}
				if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
					fprintf(stderr, "counter: %d, target: [%d, %d], query: [%d, %d], vt_dist: %u\n",
					        counter, (int32_t)(first_target_loc & UINT32_MAX),
					        (int32_t)(last_target_loc & UINT32_MAX), first_query_loc,
					        last_query_loc, vt_distance);
				}
				uint32_t chrom_id = first_target_loc >> 32;
				seqs[out_len - 1] = (vt_t){.chrom_id         = chrom_id,
				                           .first_target_loc = (uint32_t)first_target_loc,
				                           .last_target_loc  = (uint32_t)last_target_loc,
				                           .first_query_loc  = (uint32_t)first_query_loc,
				                           .last_query_loc   = (uint32_t)last_query_loc,
				                           .str              = str,
				                           .score            = counter};
				for (unsigned k = out_len - 1; k > 0; k--) {
					if (seqs[k].score > seqs[k - 1].score) {
						vt_t tmp    = seqs[k];
						seqs[k]     = seqs[k - 1];
						seqs[k - 1] = tmp;
					} else {
						break;
					}
				}
			}
			const uint64_t loc_tmp =
			    str ? (cur.target - cur.query) : cur.target - (tmp_extracted_len - cur.query);
			first_target_loc = loc_tmp;
			last_target_loc  = loc_tmp;
			first_query_loc  = cur.query;
			last_query_loc   = cur.query;
			ref_loc          = cur.target;
			counter          = 1;
		}
	}

	if (last_query_loc - first_query_loc > coverage_threshold) {
		if (out_len == vt_max_nb_locations) {
			if (seqs[out_len - 1].score >= counter) {
				vt->nb_seqs = out_len;
				return;
			}
		} else {
			out_len++;
		}
		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			fprintf(stderr, "counter: %d, target: [%d, %d], query: [%d, %d], vt_dist: %u\n", counter,
			        (int32_t)(first_target_loc & UINT32_MAX), (int32_t)(last_target_loc & UINT32_MAX),
			        first_query_loc, last_query_loc, vt_distance);
		}
		uint32_t chrom_id = first_target_loc >> 32;
		seqs[out_len - 1] = (vt_t){.chrom_id         = chrom_id,
		                           .first_target_loc = (uint32_t)first_target_loc,
		                           .last_target_loc  = (uint32_t)last_target_loc,
		                           .first_query_loc  = (uint32_t)first_query_loc,
		                           .last_query_loc   = (uint32_t)last_query_loc,
		                           .str              = str,
		                           .score            = counter};
		for (unsigned k = out_len - 1; k > 0; k--) {
			if (seqs[k].score > seqs[k - 1].score) {
				vt_t tmp    = seqs[k];
				seqs[k]     = seqs[k - 1];
				seqs[k - 1] = tmp;
			} else {
				break;
			}
		}
	}
	vt->nb_seqs = out_len;
}

void vote_2(const loc_t *const loc, const unsigned len, const int str, vt_t *const vt, const unsigned vt_distance,
            const int32_t tmp_extracted_len, const uint32_t min, const uint32_t max) {
	if (len == 0) {
		return;
	}

	vt_t best_vt = *vt;

	unsigned counter = 1;
	const uint64_t loc_tmp =
	    str ? (loc[0].target - loc[0].query) : loc[0].target - (tmp_extracted_len - loc[0].query);
	uint64_t first_target_loc = loc_tmp;
	uint64_t last_target_loc  = loc_tmp;
	uint32_t first_query_loc  = loc[0].query;
	uint32_t last_query_loc   = loc[0].query;
	uint64_t ref_loc          = loc[0].target;
	loc_t cur;

	for (unsigned i = 1; i < len; i++) {
		cur = loc[i];

		// If curent location in the range, increase the counter and update the values
		if (cur.target - ref_loc <= vt_distance) {
			if (cur.query < max && cur.query > min) {
				counter++;
				if (cur.query < first_query_loc) {
					first_query_loc = cur.query;
					ref_loc         = cur.target;
				}
				if (cur.query > last_query_loc) {
					last_query_loc = cur.query;
				}
				uint64_t loc =
				    str ? (cur.target - cur.query) : cur.target - (tmp_extracted_len - cur.query);
				if (loc > last_target_loc) {
					last_target_loc = loc;
				}
				if (loc < first_target_loc) {
					first_target_loc = loc;
				}
			}
		}
		// Else check the number of votes we have
		else {
			if (counter > best_vt.score) {
				if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
					fprintf(stderr,
					        "VT2 counter: %d, target: [%d, %d], query: [%d, %d], vt_dist: %u\n",
					        counter, (int32_t)(first_target_loc & UINT32_MAX),
					        (int32_t)(last_target_loc & UINT32_MAX), first_query_loc,
					        last_query_loc, vt_distance);
				}

				uint32_t chrom_id = first_target_loc >> 32;
				best_vt           = (vt_t){.chrom_id         = chrom_id,
                                                 .first_target_loc = (uint32_t)first_target_loc,
                                                 .last_target_loc  = (uint32_t)last_target_loc,
                                                 .first_query_loc  = (uint32_t)first_query_loc,
                                                 .last_query_loc   = (uint32_t)last_query_loc,
                                                 .str              = str,
                                                 .score            = counter};
			}
			const uint64_t loc_tmp =
			    str ? (cur.target - cur.query) : cur.target - (tmp_extracted_len - cur.query);
			first_target_loc = loc_tmp;
			last_target_loc  = loc_tmp;
			first_query_loc  = cur.query;
			last_query_loc   = cur.query;
			ref_loc          = cur.target;
			counter          = 1;
		}
	}
	if (counter > best_vt.score) {
		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			fprintf(stderr, "VT2 counter: %d, target: [%d, %d], query: [%d, %d], vt_dist: %u\n", counter,
			        (int32_t)(first_target_loc & UINT32_MAX), (int32_t)(last_target_loc & UINT32_MAX),
			        first_query_loc, last_query_loc, vt_distance);
		}

		uint32_t chrom_id = first_target_loc >> 32;
		best_vt           = (vt_t){.chrom_id         = chrom_id,
                                 .first_target_loc = (uint32_t)first_target_loc,
                                 .last_target_loc  = (uint32_t)last_target_loc,
                                 .first_query_loc  = (uint32_t)first_query_loc,
                                 .last_query_loc   = (uint32_t)last_query_loc,
                                 .str              = str,
                                 .score            = counter};
	}
	*vt = best_vt;
}

void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs,
                 mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname) {
	unsigned qlen_sum = 0;
	mm128_v mv        = {0, 0, 0};

	int Pattern_Num_Ones = 0;
	for (int f = 0; f < opt->pattern_len; ++f) {
		if (opt->pattern[f] == '1') {
			++Pattern_Num_Ones;
		}
	}

	for (unsigned i = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	PROF_INIT;
	PROF_START;

	// sketch2 or collect_minimizer2
	int shift = 0;
	mm_pattern_t mm_pattern =
	    collect_minimizers2(b->km, opt, mi, qlens[0], seqs[0], &mv, opt->pattern, opt->pattern_len, opt->max_seeds);
	shift = mm_get_shift(b->km, mi, mv, mm_pattern); // & or * or nothing
	// printf("Shift: %d\n", shift);
	// printf("%d\n\n\n\n", sizeof(shift_seeds_number)/sizeof(shift_seeds_number[0]));

	kfree(b->km, mm_pattern.shift_seeds_number);
	PROF_END(pf_pattern_alignment);

	vt_v vt = {nb_seqs : 0};
	vt.seqs = (vt_t *)kmalloc(b->km, (opt->vt_nb_loc + 2) * sizeof(vt_t));

	uint32_t max_nb_seeds =
	    (opt->flag & MM_F_FRAG_MODE) ? ((opt->max_frag_len == 0) ? 800 : opt->max_frag_len) : UINT32_MAX;

	PROF_START;
	unsigned tmp_extracted_len = collect_minimizers(b->km, opt, mi, qlens[0], seqs[0], &mv, opt->pattern,
	                                                opt->pattern_len, shift, max_nb_seeds);
	if (opt->q_occ_frac > 0.0f) mm_seed_mz_flt(b->km, &mv, opt->mid_occ, opt->q_occ_frac); // checks freq
	loc_t *a_for;
	loc_t *a_rev;
	unsigned n_a_for, n_a_rev;

	if (opt->flag & MM_F_RADIX_SORT) {
		collect_seed_hits_radix(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, tmp_extracted_len, &a_for,
		                        &a_rev, &n_a_for, &n_a_rev);
	} else {
		collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, tmp_extracted_len, &a_for, &a_rev,
		                  &n_a_for, &n_a_rev, opt->flag & MM_F_HEAP_SORT);
	}
	kfree(b->km, mv.a);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "RS n_a_for: %u, n_a_rev: %u\n", n_a_for, n_a_rev);
		for (unsigned i = 0; i < n_a_for; ++i) {
			fprintf(stderr, "SD\t%s\t%d\t+\t%u\n", mi->seq[a_for[i].target >> 32].name,
			        (int32_t)a_for[i].target + 1 - tmp_extracted_len, a_for[i].query);
		}
		for (unsigned i = 0; i < n_a_rev; ++i) {
			fprintf(stderr, "SD\t%s\t%d\t-\t%u\n", mi->seq[a_rev[i].target >> 32].name,
			        (uint32_t)a_rev[i].target + 1, a_rev[i].query);
		}
	}
	PROF_END(pf_seeding);
	PROF_START;

	const uint32_t coverage_threshold = (float)qlen_sum * opt->vt_cov;
	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "VT: cov_threshold: %u, df1: %f, f: %f\n", coverage_threshold, opt->vt_df1, opt->vt_f);
	}
	vote(a_for, n_a_for, 0, &vt, opt->vt_dis, tmp_extracted_len, opt->vt_nb_loc, coverage_threshold);
	vote(a_rev, n_a_rev, 1, &vt, opt->vt_dis, tmp_extracted_len, opt->vt_nb_loc, coverage_threshold);
	PROF_END(pf_voting);

	if (vt.nb_seqs == 0) {
		kfree(b->km, vt.seqs);
		return;
	}

	// density filter
	unsigned nb_seqs_df = 0;
	for (unsigned i = 0; i < vt.nb_seqs; i++) {
		if ((float)vt.seqs[i].score >
		    opt->vt_df1 * (float)(vt.seqs[i].last_target_loc - vt.seqs[i].first_target_loc)) {
			vt.seqs[i] = vt.seqs[nb_seqs_df];
			nb_seqs_df++;
		}
	}

	vt.nb_seqs = nb_seqs_df;
	if (vt.nb_seqs == 0) {
		kfree(b->km, vt.seqs);
		return;
	}

	// filter out the bad seqs
	// & adjust the first location for each potential
	// & check if we cover the entire read
	const unsigned bw                  = opt->bw;
	uint32_t start                     = qlen_sum;
	uint32_t end                       = 0;
	const unsigned filtering_threshold = (float)vt.seqs[0].score * opt->vt_f;
	for (unsigned i = 0; i < vt.nb_seqs; i++) {
		if (vt.seqs[i].score < filtering_threshold) {
			vt.nb_seqs = i;
			// filtered   = 1;
			break;
		}
		vt.seqs[i].first_query_loc -= (mi->k - 1);
		vt.seqs[i].first_target_loc -= (mi->k - 1);
		vt.seqs[i].next   = NULL;
		vt.seqs[i].concat = 0;
		// To avoid overflow during alignment
		if (vt.seqs[i].last_query_loc - vt.seqs[i].first_query_loc + 0.5 * bw <
		    vt.seqs[i].last_target_loc - vt.seqs[i].first_target_loc) {
			vt.seqs[i].last_target_loc = vt.seqs[i].first_target_loc + vt.seqs[i].last_query_loc -
			                             vt.seqs[i].first_query_loc + 0.5 * bw;
		}
		if (vt.seqs[i].first_query_loc < start) {
			start = vt.seqs[i].first_query_loc;
		}
		if (vt.seqs[i].last_query_loc > end) {
			end = vt.seqs[i].last_query_loc;
		}
	}

	// Perform a second round of voting if we don't cover the entire read
	if (start > coverage_threshold) {
		vt_t vt2 = {.score = 0};
		vote_2(a_for, n_a_for, 0, &vt2, opt->vt_dis, tmp_extracted_len, 0, start);
		vote_2(a_rev, n_a_rev, 1, &vt2, opt->vt_dis, tmp_extracted_len, 0, start);
		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			fprintf(stderr, "start: %u\n", start);
			const char *chrom_name = mi->seq[vt2.chrom_id].name;
			fprintf(stderr, "VT2: score: %u, chrom: %s, loc: %u, df2: %f\n", vt2.score, chrom_name,
			        vt2.first_target_loc, opt->vt_df2);
		}
		vt2.first_query_loc -= (mi->k - 1);
		vt2.first_target_loc -= (mi->k - 1);
		if ((float)vt2.score > opt->vt_df2 * (float)(vt2.last_target_loc - vt2.first_target_loc)) {
			if (vt2.last_query_loc - vt2.first_query_loc + 0.5 * bw <
			    vt2.last_target_loc - vt2.first_target_loc) {
				vt2.last_target_loc =
				    vt2.first_target_loc + vt2.last_query_loc - vt2.first_query_loc + 0.5 * bw;
			}
			vt.seqs[vt.nb_seqs] = vt2;
			vt.nb_seqs++;
		}
	}
	if (qlen_sum - end > coverage_threshold) {
		vt_t vt2 = {.score = 0};
		vote_2(a_for, n_a_for, 0, &vt2, opt->vt_dis, tmp_extracted_len, end, qlen_sum);
		vote_2(a_rev, n_a_rev, 1, &vt2, opt->vt_dis, tmp_extracted_len, 0, start);
		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			const char *chrom_name = mi->seq[vt2.chrom_id].name;
			fprintf(stderr, "VT2: score: %u, chrom: %s, loc: %u, df2: %f\n", vt2.score, chrom_name,
			        vt2.first_target_loc, opt->vt_df2);
		}
		vt2.first_query_loc -= (mi->k - 1);
		vt2.first_target_loc -= (mi->k - 1);
		if ((float)vt2.score > opt->vt_df2 * (float)(vt2.last_target_loc - vt2.first_target_loc)) {
			if (vt2.last_query_loc - vt2.first_query_loc + 0.5 * bw <
			    vt2.last_target_loc - vt2.first_target_loc) {
				vt2.last_target_loc =
				    vt2.first_target_loc + vt2.last_query_loc - vt2.first_query_loc + 0.5 * bw;
			}
			vt.seqs[vt.nb_seqs] = vt2;
			vt.nb_seqs++;
		}
	}
	kfree(b->km, a_for);
	kfree(b->km, a_rev);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "VT n: %u, len: %u\n", vt.nb_seqs, qlen_sum);
		for (unsigned i = 0; i < vt.nb_seqs; ++i) {
			const vt_t potential    = vt.seqs[i];
			const char *chrom_name  = mi->seq[potential.chrom_id].name;
			const int32_t chrom_len = mi->seq[potential.chrom_id].len;
			fprintf(stderr, "VT\t%s (len: %u)\t[%u, %u]\t%c\t[%u, %u]\t%u\n", chrom_name, chrom_len,
			        potential.first_target_loc, potential.last_target_loc, "+-"[potential.str],
			        potential.first_query_loc, potential.last_query_loc, potential.score);
		}
	}

	/*
	// DEBUGING
	kfree(b->km, vt.seqs);
	return;
	*/

	const unsigned max_max_gap = opt->max_max_gap;
	const unsigned max_min_gap = opt->max_min_gap;
	// Check which potential sequences could be concatenated
	for (unsigned i = 0; i < vt.nb_seqs; i++) {
		// First choose the best potential sequence
		vt_t *const s1 = &vt.seqs[i];
		for (unsigned j = 0; j < vt.nb_seqs; j++) {
			if (j != i) {
				vt_t *const s2 = &vt.seqs[j];
				if (s2->concat == 0 && s1->str == s2->str && s1->chrom_id == s2->chrom_id) {
					if (s1->str) {
						if (s2->last_query_loc < s1->first_query_loc &&
						    s1->last_target_loc > s2->first_target_loc &&
						    s1->first_target_loc < s2->first_target_loc) {
							if (s2->last_query_loc + max_max_gap > s1->first_query_loc) {
								if (s1->next == NULL) {
									s1->next = s2;
								} else if (s2->last_query_loc >
								           s1->next->last_query_loc) {
									s1->next = s2;
								}
							}
						} else if (s2->last_query_loc < s1->first_query_loc &&
						           s1->last_target_loc < s2->first_target_loc) {
							if ((s2->last_query_loc + max_min_gap > s1->first_query_loc ||
							     s1->last_target_loc + max_min_gap >
							         s2->first_target_loc) &&
							    s2->last_query_loc + max_max_gap > s1->first_query_loc &&
							    s1->last_target_loc + max_max_gap > s2->first_target_loc) {
								if (s1->next == NULL) {
									s1->next = s2;
								} else if (s2->last_query_loc >
								           s1->next->last_query_loc) {
									s1->next = s2;
								}
							}
						} else if (s2->last_query_loc > s1->first_query_loc &&
						           s1->last_target_loc < s2->first_target_loc &&
						           s2->last_query_loc < s1->last_query_loc &&
						           s2->first_query_loc < s1->first_query_loc) {
							if (s1->last_target_loc + max_max_gap > s2->first_target_loc) {
								if (s1->next == NULL) {
									s1->next = s2;
								} else if (s2->last_query_loc <
								           s1->next->last_query_loc) {
									s1->next = s2;
								}
							}
						}
					} else {
						if (s1->last_query_loc < s2->first_query_loc &&
						    s1->last_target_loc > s2->first_target_loc &&
						    s1->first_target_loc < s2->first_target_loc) {
							if (s1->last_query_loc + max_max_gap > s2->first_query_loc) {
								if (s1->next == NULL) {
									s1->next = s2;
								} else if (s2->first_query_loc <
								           s1->next->first_query_loc) {
									s1->next = s2;
								}
							}
						} else if (s1->last_query_loc < s2->first_query_loc &&
						           s1->last_target_loc < s2->first_target_loc) {
							if ((s1->last_query_loc + max_min_gap > s2->first_query_loc ||
							     s1->last_target_loc + max_min_gap >
							         s2->first_target_loc) &&
							    s1->last_target_loc + max_max_gap > s2->first_target_loc &&
							    s1->last_query_loc + max_max_gap > s2->first_query_loc) {
								if (s1->next == NULL) {
									s1->next = s2;
								} else if (s2->first_query_loc <
								           s1->next->first_query_loc) {
									s1->next = s2;
								}
							}
						} else if (s1->last_query_loc > s2->first_query_loc &&
						           s1->last_target_loc < s2->first_target_loc &&
						           s1->first_query_loc < s2->first_query_loc &&
						           s1->last_query_loc < s2->last_query_loc) {
							if (s1->last_target_loc + max_max_gap > s2->first_target_loc) {
								if (s1->next == NULL) {
									s1->next = s2;
								} else if (s2->first_query_loc <
								           s1->next->first_query_loc) {
									s1->next = s2;
								}
							}
						}
					}
				}
			}
		}
		// Second, adjust the boundaries of the sequences
		if (s1->next != NULL) {
			vt_t *const s2 = s1->next;
			s2->concat     = 1;
			if (s1->str) {
				if (s2->last_query_loc < s1->first_query_loc &&
				    s1->last_target_loc < s2->first_target_loc) {
					const uint32_t diffq = s1->first_query_loc - s2->last_query_loc;
					const uint32_t difft = s2->first_target_loc - s1->last_target_loc;
					const uint32_t min   = difft > diffq ? diffq : difft;
					s2->last_query_loc += min;
					s1->last_target_loc += min;
					s1->first_query_loc -= min;
					s2->first_target_loc -= min;
				}
			} else {
				if (s1->last_query_loc < s2->first_query_loc &&
				    s1->last_target_loc < s2->first_target_loc) {
					const uint32_t diffq = s2->first_query_loc - s1->last_query_loc;
					const uint32_t difft = s2->first_target_loc - s1->last_target_loc;
					const uint32_t min   = difft > diffq ? diffq : difft;
					s1->last_query_loc += min;
					s1->last_target_loc += min;
					s2->first_query_loc -= min;
					s2->first_target_loc -= min;
				}
			}
			if (s2->last_target_loc < s1->last_target_loc) {
				s1->last_target_loc = s2->last_target_loc - 1;
			}
		}
	}

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "AVT n: %u, len: %u\n", vt.nb_seqs, qlen_sum);
		for (unsigned i = 0; i < vt.nb_seqs; ++i) {
			const vt_t potential    = vt.seqs[i];
			const char *chrom_name  = mi->seq[potential.chrom_id].name;
			const int32_t chrom_len = mi->seq[potential.chrom_id].len;
			fprintf(stderr, "AVT\t%s (len: %u)\t[%u, %u]\t%c\t[%u, %u]\t%u\tc:%u\n", chrom_name, chrom_len,
			        potential.first_target_loc, potential.last_target_loc, "+-"[potential.str],
			        potential.first_query_loc, potential.last_query_loc, potential.score, potential.concat);
		}
	}

	PROF_START;

	// Get the read sequence (either only one strand or both if required)
	int str_for      = 0;
	int str_rev      = 0;
	uint32_t max_len = 0;
	for (unsigned i = 0; i < vt.nb_seqs; i++) {
		uint32_t len = vt.seqs[i].last_target_loc - vt.seqs[i].first_target_loc + 1;
		if (len > max_len) {
			max_len = len;
		}
		if (vt.seqs[i].str) {
			str_rev = 1;
		} else {
			str_for = 1;
		}
	}

	uint8_t *qs_for = NULL;
	uint8_t *qs_rev = NULL;

	extern unsigned char seq_nt4_table[256];
	if (!str_rev) {
		qs_for = (uint8_t *)kmalloc(b->km, qlen_sum * sizeof(uint8_t));
		for (size_t j = 0; j < qlen_sum; j++) {
			qs_for[j] = seq_nt4_table[(uint8_t)(seqs[0][j])];
		}
	} else if (!str_for) {
		qs_rev = (uint8_t *)kmalloc(b->km, qlen_sum * sizeof(uint8_t));
		for (size_t j = 0; j < qlen_sum; j++) {
			qs_rev[qlen_sum - j - 1] = seq_nt4_table[(uint8_t)(seqs[0][j])] ^ 0x3;
		}
	} else {
		qs_for = (uint8_t *)kmalloc(b->km, qlen_sum * sizeof(uint8_t));
		qs_rev = (uint8_t *)kmalloc(b->km, qlen_sum * sizeof(uint8_t));
		for (size_t j = 0; j < qlen_sum; j++) {
			qs_for[j]                = seq_nt4_table[(uint8_t)(seqs[0][j])];
			qs_rev[qlen_sum - j - 1] = qs_for[j] ^ 0x3;
		}
	}

	uint8_t *qseq;
	uint8_t *tseq;

	if (qlen_sum > 300) {
		tseq = (uint8_t *)kmalloc(b->km, max_len * sizeof(uint8_t));
	} else {
		tseq = (uint8_t *)kmalloc(b->km, (max_len + qlen_sum) * sizeof(uint8_t));
	}

	for (unsigned i = 0; i < vt.nb_seqs; i++) {
		vt.seqs[i].valid   = 1;
		int str            = vt.seqs[i].str;
		unsigned target_id = vt.seqs[i].chrom_id;
		;
		uint32_t target_start = vt.seqs[i].first_target_loc;
		uint32_t target_end   = vt.seqs[i].last_target_loc;
		uint32_t query_start;
		uint32_t query_end;
		if (str) {
			query_end   = qlen_sum - 1 - vt.seqs[i].first_query_loc;
			query_start = qlen_sum - 1 - vt.seqs[i].last_query_loc;
		} else {
			query_start = vt.seqs[i].first_query_loc;
			query_end   = vt.seqs[i].last_query_loc;
		}
		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			const char *chrom_name = mi->seq[vt.seqs[i].chrom_id].name;
			fprintf(stderr, "BE\t%s, [%u, %u[ (chrom_len: %u) -> '%c' [%u, %u[ (read_len: %u)\n",
			        chrom_name, target_start, target_end, mi->seq[target_id].len, "+-"[str], query_start,
			        query_end, qlen_sum);
		}

		// Get the ref sequence
		if (qlen_sum > 300) {
			/*
			   //TODO: skiping
			if (vt.seqs[i].first_query_loc == vt.seqs[i].last_query_loc) {
			        if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			                fprintf(stderr, "SKIPPED\n");
			        }
			        continue;
			}
			*/
			qseq = str ? &qs_rev[query_start] : &qs_for[query_start];
		} else {
			int32_t chrom_len = mi->seq[target_id].len;
			if (target_start < query_start) {
				query_start -= target_start;
				target_start = 0;
			} else {
				target_start -= query_start;
				query_start = 0;
			}
			if (chrom_len + query_end < qlen_sum + target_end) {
				query_end += chrom_len - target_end - 1;
				target_end = chrom_len - 1;
			} else {
				target_end += qlen_sum - query_end - 1;
				query_end = qlen_sum - 1;
			}
			qseq = str ? &qs_rev[query_start] : &qs_for[query_start];
		}
		uint32_t qlen = query_end - query_start + 1;
		uint32_t tlen = target_end - target_start + 1;
		if (str) {
			uint32_t tmp = qlen_sum - 1 - query_start;
			query_start  = qlen_sum - 1 - query_end;
			query_end    = tmp;
		}
		mm_idx_getseq2(mi, 0, target_id, target_start, target_end + 1, tseq);

		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			if (str) {
				fprintf(stderr, "Read: str: -, [%u, %u], len: %u\n", qlen_sum - 1 - query_end,
				        qlen_sum - 1 - query_start, qlen);
			} else {
				fprintf(stderr, "Read: str: +, [%u, %u], len: %u\n", query_start, query_end, qlen);
			}
			for (unsigned b = 0; b < qlen; b++) {
				fprintf(stderr, "%c", "ACGTN"[qseq[b]]);
			}
			fprintf(stderr, "\nRef: %s [%d, %d], len: %u\n", mi->seq[target_id].name, target_start,
			        target_end, tlen);
			for (unsigned b = 0; b < tlen; b++) {
				fprintf(stderr, "%c", "ACGTN"[tseq[b]]);
			}
			fputs("\n", stderr);
		}

		// Alignment
		ksw_extz_t ez;
		int sc_mch = opt->a;
		int sc_mis = opt->b;
		int g = sc_mch, bb = sc_mis < 0 ? sc_mis : -sc_mis; // g>0 and b<0
		int8_t mat[25] = {g,  bb, bb, bb, 0,  bb, g, bb, bb, 0, bb, bb, g,
		                  bb, 0,  bb, bb, bb, g,  0, 0,  0,  0, 0,  0};
		memset(&ez, 0, sizeof(ksw_extz_t));
		int flag = KSW_EZ_APPROX_MAX;

		// Exact match filtering (GenStore 2022)
		bool exact_match = false;
		int mismatch_cnt = 0;

		if (qlen_sum < 300 && qlen == tlen) { // TODO: check if short read preset is used

			exact_match_sse(b->km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->bw, opt->zdrop,
			                opt->end_bonus, flag, &ez, &exact_match, &mismatch_cnt);

			/*
			exact_match=true;
			for (size_t j = 0; j <qlen_sum; ++j) {
			        if (qs[j] != ts[j]){
			                exact_match=false;
			                break;
			        }
			}
			*/
			if (exact_match) {
				// *exact_cnt = *exact_cnt + 1;
				if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
					fprintf(stderr, "Query: ");
					for (size_t j = 0; j < qlen; j++)
						fprintf(stderr, "%c", "ACGTN"[qseq[j]]);
					fprintf(stderr, "\nTarget: ");
					for (size_t j = 0; j < tlen; j++)
						fprintf(stderr, "%c", "ACGTN"[tseq[j]]);
				}

				ksw_reset_extz(&ez);
				/*
				ez.score = 0;
				for (size_t j = 0; j < qlen_sum; ++j) {
				        if (qseq[j] >= 4 || tseq[j] >= 4) ez.score += opt->e;
				        else
				                ez.score += qseq[j] == tseq[j] ? opt->a : -opt->b;
				}
				*/
				ez.score = qlen_sum * sc_mch;
				ez.cigar =
				    ksw_push_cigar(b->km, &ez.n_cigar, &ez.m_cigar, ez.cigar, MM_CIGAR_MATCH, qlen);

				if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
					fprintf(stderr, "\nCigar: ");
					for (unsigned j = 0; j < ez.n_cigar; ++j)
						fprintf(stderr, "%d%c", ez.cigar[j] >> 4, "MIDNSH"[ez.cigar[j] & 0xf]);
					fprintf(stderr, "\n");
				}
			}
		}

		if (!exact_match) {
			/*printf("bw: %u, bw_frac: %f, bw_min: %u, bw_max: %u\n", bw, opt->bw_frac,
			   opt->bw_min, opt->bw_max);*/
			// fprintf(stderr, "bw: %u\n", bw);
#ifdef __AVX512BW__
			ksw_extd2_avx512(b->km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->e2, bw,
			                 opt->zdrop, opt->end_bonus, flag, &ez);
#else
			ksw_extd2_sse(b->km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->e2, bw,
			              opt->zdrop, opt->end_bonus, flag, &ez);
#endif
		}

		if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
			fprintf(stderr, "AL_SCORE: %d\n", ez.score);
		}
		// Discard the location if score is too low
		if (ez.score == KSW_NEG_INF) {
			kfree(b->km, ez.cigar);
			vt.seqs[i].valid = 0;
			continue;
		}

		mm_reg1_t r_tmp = {.rid   = target_id,
		                   .score = ez.score,
		                   .qs    = query_start,
		                   .qe    = query_end + 1,
		                   .rs    = target_start,
		                   .re    = target_end + 1,
		                   .rev   = str};

		if (ez.n_cigar >= 0) {
			uint32_t capacity = ez.n_cigar + sizeof(mm_extra_t) / 4;
			kroundup32(capacity);
			r_tmp.p           = (mm_extra_t *)calloc(capacity, 4);
			r_tmp.p->capacity = capacity;
			r_tmp.p->n_cigar  = ez.n_cigar;
			memcpy(r_tmp.p->cigar, ez.cigar, (r_tmp.p->n_cigar) * 4);
		} else {
			r_tmp.p = (mm_extra_t *)calloc(sizeof(mm_extra_t), 1);
		}
		r_tmp.p->dp_score = ez.score;
		kfree(b->km, ez.cigar);

		mm_update_extra(&r_tmp, qseq, tseq, mat, opt->q, opt->e, opt->flag & MM_F_EQX, !(opt->flag & MM_F_SR));

		uint32_t clip_len[2];
		clip_len[0] = r_tmp.rev ? qlen_sum - r_tmp.qe : r_tmp.qs;
		clip_len[1] = r_tmp.rev ? r_tmp.qs : qlen_sum - r_tmp.qe;

		// Discard the location if these conditions are not met
		if (!(clip_len[0] < qlen_sum && clip_len[1] < qlen_sum)) {
			free(r_tmp.p);
			vt.seqs[i].valid = 0;
			continue;
		}

		vt.seqs[i].r = r_tmp;
	}

	// Concatenate the records
	for (unsigned i = 0; i < vt.nb_seqs; i++) {
		while (vt.seqs[i].valid && vt.seqs[i].next != NULL && vt.seqs[i].next->valid) {
			if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
				fprintf(stderr, "CONQ[%u, %u] || [%u, %u]\n", vt.seqs[i].r.qs, vt.seqs[i].r.qe,
				        vt.seqs[i].next->r.qs, vt.seqs[i].next->r.qe);
				fprintf(stderr, "CONT [%u, %u] || [%u, %u]\n", vt.seqs[i].r.rs, vt.seqs[i].r.re,
				        vt.seqs[i].next->r.rs, vt.seqs[i].next->r.re);
			}
			if (concatenate_cigars(&vt.seqs[i].r, vt.seqs[i].next->r, vt.seqs[i].str ? qs_rev : qs_for,
			                       vt.seqs[i].str, qlen_sum, mi, b->km, opt->a, opt->b, opt->q, opt->e,
			                       opt->q2, opt->e2) == 0) {
				free(vt.seqs[i].next->r.p);
				vt.seqs[i].next->valid = 0;
				vt.seqs[i].next        = vt.seqs[i].next->next;
			} else {
				vt.seqs[i].next = NULL;
			}
		}
	}

	// Count the number of remaining sequences
	unsigned nb_out = 0;
	for (unsigned i = 0; i < vt.nb_seqs; i++) {
		if (vt.seqs[i].valid) {
			nb_out++;
		}
	}

	if (nb_out > 0) {
		// Alocate the output buffer and sort the sequences according to the alignment score
		mm_reg1_t *r;
		r = (mm_reg1_t *)malloc(nb_out * sizeof(mm_reg1_t));

		unsigned out_pos = 0;
		for (unsigned i = 0; i < vt.nb_seqs; i++) {
			if (vt.seqs[i].valid) {
				r[out_pos] = vt.seqs[i].r;
				unsigned j = out_pos;
				out_pos++;
				while (j > 0) {
					if (r[j].score > r[j - 1].score) {
						mm_reg1_t r_tmp = r[j];
						r[j]            = r[j - 1];
						r[j - 1]        = r_tmp;
					} else {
						break;
					}
				}
			}
		}
		const unsigned max_nb_sec = (opt->flag & MM_F_NO_PRINT_2ND) ? 0 : opt->best_n;
		// Set Mapq
		mm_set_sam_params(b->km, nb_out, r, qlen_sum, opt->a, max_nb_sec);
		*regs = r;
	} else {
		*regs = NULL;
	}

	kfree(b->km, vt.seqs);
	kfree(b->km, tseq);
	kfree(b->km, qs_for);
	kfree(b->km, qs_rev);
	*n_regs = nb_out;
	PROF_END(pf_sequence_alignment);

	/*
	// Limit the memory
	if (b->km) {
	        km_stat_t kmst;
	        km_stat(b->km, &kmst);
	        if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
	                fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum,
	kmst.capacity, kmst.n_cores, kmst.largest);
	        // assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
	        if (kmst.largest > 1U << 28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc))
	{ if (mm_dbg_flag & MM_DBG_PRINT_QNAME) fprintf(stderr, "[W::%s] reset thread-local memory after
	read %s\n", __func__, qname); km_destroy(b->km); b->km = km_init();
	        }
	}
	*/
}

/*
mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const
mm_mapopt_t *opt, const char *qname) { mm_reg1_t *regs; mm_map_frag(mi, 1, &qlen, &seq, n_regs, &regs,
b, opt, qname); return regs;
}
*/

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int n_processed, n_threads, n_fp;
	int64_t mini_batch_size;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;

	int n_parts;
	uint32_t *rid_shift;
	FILE *fp_split, **fp_parts;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
	int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
	step_t *s = (step_t *)_data;
	int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MM_MAX_SEG];
	double t     = 0.0;
	mm_tbuf_t *b = s->buf[tid];
	assert(s->n_seg[i] <= MM_MAX_SEG);
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME) {
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
		t = realtime();
	}
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori >> 1 & 1)) || (j == 1 && (pe_ori & 1))))
			mm_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
	}
	if (s->p->opt->flag & MM_F_INDEPEND_SEG) {
		for (j = 0; j < s->n_seg[i]; ++j) {
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off + j], &s->reg[off + j], b,
			            s->p->opt, s->seq[off + j].name);
			s->rep_len[off + j]  = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt,
		            s->seq[off].name);
		for (j = 0; j < s->n_seg[i]; ++j) {
			s->rep_len[off + j]  = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	}
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori >> 1 & 1)) || (j == 1 && (pe_ori & 1)))) {
			int k, t;
			mm_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mm_reg1_t *r = &s->reg[off + j][k];
				t            = r->qs;
				r->qs        = qlens[j] - r->qe;
				r->qe        = qlens[j] - t;
				r->rev       = !r->rev;
			}
		}
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QT\t%s\t%d\t%.6f\n", s->seq[off].name, tid, realtime() - t);
}

static void merge_hits(step_t *s) {
	int f, i, k0, k, max_seg = 0, *n_reg_part, *rep_len_part, *frag_gap_part, *qlens;
	void *km;
	FILE **fp              = s->p->fp_parts;
	const mm_mapopt_t *opt = s->p->opt;

	km = km_init();
	for (f = 0; f < s->n_frag; ++f)
		max_seg = max_seg > s->n_seg[f] ? max_seg : s->n_seg[f];
	qlens         = CALLOC(int, max_seg + s->p->n_parts * 3);
	n_reg_part    = qlens + max_seg;
	rep_len_part  = n_reg_part + s->p->n_parts;
	frag_gap_part = rep_len_part + s->p->n_parts;
	for (f = 0, k = k0 = 0; f < s->n_frag; ++f) {
		k0 = k;
		for (i = 0; i < s->n_seg[f]; ++i, ++k) {
			int j, l, t, rep_len = 0;
			qlens[i] = s->seq[k].l_seq;
			for (j = 0, s->n_reg[k] = 0; j < s->p->n_parts; ++j) {
				mm_err_fread(&n_reg_part[j], sizeof(int), 1, fp[j]);
				mm_err_fread(&rep_len_part[j], sizeof(int), 1, fp[j]);
				mm_err_fread(&frag_gap_part[j], sizeof(int), 1, fp[j]);
				s->n_reg[k] += n_reg_part[j];
				if (rep_len < rep_len_part[j]) rep_len = rep_len_part[j];
			}
			s->reg[k] = CALLOC(mm_reg1_t, s->n_reg[k]);
			for (j = 0, l = 0; j < s->p->n_parts; ++j) {
				for (t = 0; t < n_reg_part[j]; ++t, ++l) {
					mm_reg1_t *r = &s->reg[k][l];
					uint32_t capacity;
					mm_err_fread(r, sizeof(mm_reg1_t), 1, fp[j]);
					r->rid += s->p->rid_shift[j];
					if (opt->flag & MM_F_CIGAR) {
						mm_err_fread(&capacity, 4, 1, fp[j]);
						r->p           = (mm_extra_t *)calloc(capacity, 4);
						r->p->capacity = capacity;
						mm_err_fread(r->p, r->p->capacity, 4, fp[j]);
					}
				}
			}
			if (!(opt->flag & MM_F_SR) && s->seq[k].l_seq >= opt->rank_min_len)
				mm_update_dp_max(s->seq[k].l_seq, s->n_reg[k], s->reg[k], opt->rank_frac, opt->a,
				                 opt->b);
			for (j = 0; j < s->n_reg[k]; ++j) {
				mm_reg1_t *r = &s->reg[k][j];
				if (r->p)
					r->p->dp_max2 = 0; // reset ->dp_max2 as mm_set_parent() doesn't clear
					                   // it; necessary with mm_update_dp_max()
				r->subsc = 0;              // this may not be necessary
				r->n_sub = 0;              // n_sub will be an underestimate as we don't see all the
				                           // chains now, but it can't be accurate anyway
			}
			mm_hit_sort(km, &s->n_reg[k], s->reg[k], opt->alt_drop);
			mm_set_parent(km, opt->mask_level, opt->mask_len, s->n_reg[k], s->reg[k], opt->a * 2 + opt->b,
			              opt->flag & MM_F_HARD_MLEVEL, opt->alt_drop);
			if (!(opt->flag & MM_F_ALL_CHAINS)) {
				mm_select_sub(km, opt->pri_ratio, s->p->mi->k * 2, opt->best_n, 0, opt->max_gap * 0.8,
				              &s->n_reg[k], s->reg[k]);
				mm_set_sam_pri(s->n_reg[k], s->reg[k]);
			}
			mm_set_mapq(km, s->n_reg[k], s->reg[k], opt->min_chain_score, opt->a, rep_len,
			            !!(opt->flag & MM_F_SR));
		}
		if (s->n_seg[f] == 2 && opt->pe_ori >= 0 && (opt->flag & MM_F_CIGAR))
			mm_pair(km, frag_gap_part[0], opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, &s->n_reg[k0],
			        &s->reg[k0]);
	}
	free(qlens);
	km_destroy(km);
}

static void *worker_pipeline(void *shared, int step, void *in) {
	int i, j, k;
	pipeline_t *p = (pipeline_t *)shared;
	if (step == 0) { // step 0: read sequences
		int with_qual    = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
		int with_comment = !!(p->opt->flag & MM_F_COPY_COMMENT);
		// int frag_mode    = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
		step_t *s;
		s = (step_t *)calloc(1, sizeof(step_t));
		if (p->n_fp > 1)
			s->seq =
			    mm_bseq_read_frag2(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
		else
			s->seq = mm_bseq_read3(p->fp[0], p->mini_batch_size, with_qual, with_comment, 0, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t **)calloc(p->n_threads, sizeof(mm_tbuf_t *));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mm_tbuf_init();
			s->n_reg   = (int *)calloc(5 * s->n_seq, sizeof(int));
			s->seg_off = s->n_reg + s->n_seq; // seg_off, n_seg, rep_len and frag_gap are
			                                  // allocated together with n_reg
			s->n_seg    = s->seg_off + s->n_seq;
			s->rep_len  = s->n_seg + s->n_seq;
			s->frag_gap = s->rep_len + s->n_seq;
			s->reg      = (mm_reg1_t **)calloc(s->n_seq, sizeof(mm_reg1_t *));
			for (i = 1, j = 0; i <= s->n_seq; ++i)
				if (i == s->n_seq || !mm_qname_same(s->seq[i - 1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag]     = i - j;
					s->seg_off[s->n_frag++] = j;
					j                       = i;
				}
			return s;
		} else
			free(s);
	} else if (step == 1) { // step 1: map
		if (p->n_parts > 0)
			merge_hits((step_t *)in);
		else
			kt_for(p->n_threads, worker_for, in, ((step_t *)in)->n_frag);
		return in;
	} else if (step == 2) { // step 2: output
		void *km           = 0;
		step_t *s          = (step_t *)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i)
			mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				if (p->opt->split_prefix && p->n_parts == 0) { // then write to temporary files
					mm_err_fwrite(&s->n_reg[i], sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->rep_len[i], sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->frag_gap[i], sizeof(int), 1, p->fp_split);
					for (j = 0; j < s->n_reg[i]; ++j) {
						mm_reg1_t *r = &s->reg[i][j];
						mm_err_fwrite(r, sizeof(mm_reg1_t), 1, p->fp_split);
						if (p->opt->flag & MM_F_CIGAR) {
							mm_err_fwrite(&r->p->capacity, 4, 1, p->fp_split);
							mm_err_fwrite(r->p, r->p->capacity, 4, p->fp_split);
						}
					}
				} else if (s->n_reg[i] > 0) { // the query has at least one hit
					for (j = 0; j < s->n_reg[i]; ++j) {
						mm_reg1_t *r = &s->reg[i][j];
						assert(!r->sam_pri || r->id == r->parent);
						if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent) continue;
						if (p->opt->flag & MM_F_OUT_SAM)
							mm_write_sam3(&p->str, mi, t, i - seg_st, j, s->n_seg[k],
							              &s->n_reg[seg_st],
							              (const mm_reg1_t *const *)&s->reg[seg_st], km,
							              p->opt->flag, s->rep_len[i]);
						else
							mm_write_paf3(&p->str, mi, t, r, km, p->opt->flag,
							              s->rep_len[i]);
						mm_err_puts(p->str.s);
					}
				} else if ((p->opt->flag & MM_F_PAF_NO_HIT) ||
				           ((p->opt->flag & MM_F_OUT_SAM) &&
				            !(p->opt->flag & MM_F_SAM_HIT_ONLY))) { // output an empty
					                                            // hit, if requested
					if (p->opt->flag & MM_F_OUT_SAM)
						mm_write_sam3(&p->str, mi, t, i - seg_st, -1, s->n_seg[k],
						              &s->n_reg[seg_st],
						              (const mm_reg1_t *const *)&s->reg[seg_st], km,
						              p->opt->flag, s->rep_len[i]);
					else
						mm_write_paf3(&p->str, mi, t, 0, 0, p->opt->flag, s->rep_len[i]);
					mm_err_puts(p->str.s);
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_reg[i]; ++j)
					free(s->reg[i][j].p);
				free(s->reg[i]);
				free(s->seq[i].seq);
				free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
				if (s->seq[i].comment) free(s->seq[i].comment);
			}
		}
		free(s->reg);
		free(s->n_reg);
		free(s->seq); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory
		              // leak here
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0,
			        cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
	return 0;
}

static mm_bseq_file_t **open_bseqs(int n, const char **fn) {
	mm_bseq_file_t **fp;
	int i, j;
	fp = (mm_bseq_file_t **)calloc(n, sizeof(mm_bseq_file_t *));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = mm_bseq_open(fn[i])) == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j)
				mm_bseq_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads) {
	int i, pl_threads;
	pipeline_t pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp   = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.n_threads       = n_threads > 1 ? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	if (opt->split_prefix) pl.fp_split = mm_split_init(opt->split_prefix, idx);
	pl_threads = n_threads == 1 ? 1 : (opt->flag & MM_F_2_IO_THREADS) ? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);

	free(pl.str.s);
	if (pl.fp_split) fclose(pl.fp_split);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads) {
	return mm_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int mm_split_merge(int n_segs, const char **fn, const mm_mapopt_t *opt, int n_split_idx) {
	int i;
	pipeline_t pl;
	mm_idx_t *mi;
	if (n_segs < 1 || n_split_idx < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp   = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt             = opt;
	pl.mini_batch_size = opt->mini_batch_size;

	pl.n_parts   = n_split_idx;
	pl.fp_parts  = CALLOC(FILE *, pl.n_parts);
	pl.rid_shift = CALLOC(uint32_t, pl.n_parts);
	pl.mi = mi = mm_split_merge_prep(opt->split_prefix, n_split_idx, pl.fp_parts, pl.rid_shift);
	if (pl.mi == 0) {
		free(pl.fp_parts);
		free(pl.rid_shift);
		return -1;
	}
	for (i = n_split_idx - 1; i > 0; --i)
		pl.rid_shift[i] = pl.rid_shift[i - 1];
	for (pl.rid_shift[0] = 0, i = 1; i < n_split_idx; ++i)
		pl.rid_shift[i] += pl.rid_shift[i - 1];
	if (opt->flag & MM_F_OUT_SAM)
		for (i = 0; i < (int32_t)pl.mi->n_seq; ++i)
			printf("@SQ\tSN:%s\tLN:%d\n", pl.mi->seq[i].name, pl.mi->seq[i].len);

	kt_pipeline(2, worker_pipeline, &pl, 3);

	free(pl.str.s);
	mm_idx_destroy(mi);
	free(pl.rid_shift);
	for (i = 0; i < n_split_idx; ++i)
		fclose(pl.fp_parts[i]);
	free(pl.fp_parts);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	mm_split_rm_tmp(opt->split_prefix, n_split_idx);
	return 0;
}

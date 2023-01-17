CFLAGS=		-g -Wall -O2 -Wc++-compat #-fsanitize=address -static-libasan #-Wextra
CPPFLAGS=	-DHAVE_KALLOC -DPROFILE
INCLUDES=
OBJS=		kthread.o kalloc.o misc.o bseq.o sdust.o options.o index.o \
			lchain.o align.o hit.o seed.o format.o pe.o esterr.o splitidx.o \
			ksw2_ll_sse.o profile.o
PROG=		GDiet
PROG_EXTRA=	sdust minimap2-lite
LIBS=		-lm -lz -lpthread

VPATH= src

ifeq ($(arm_neon),) # if arm_neon is not defined
ifeq ($(sse2only),) # if sse2only is not defined
	OBJS+=ksw2_extz2_sse41.o exact_match_sse41.o exact_match_sse2.o ksw2_extd2_sse41.o ksw2_exts2_sse41.o ksw2_extz2_sse2.o ksw2_extd2_sse2.o ksw2_exts2_sse2.o ksw2_dispatch.o
else                # if sse2only is defined
	OBJS+=ksw2_extz2_sse.o exact_match_sse.o ksw2_extd2_sse.o ksw2_exts2_sse.o
endif
else				# if arm_neon is defined
	OBJS+=ksw2_extz2_neon.o exact_match_neon.o ksw2_extd2_neon.o ksw2_exts2_neon.o
    INCLUDES+=-Isse2neon
ifeq ($(aarch64),)	#if aarch64 is not defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -mfpu=neon -fsigned-char
else				#if aarch64 is defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -fsigned-char
endif
endif

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

ifneq ($(tsan),)
	CFLAGS+=-fsanitize=thread
	LIBS+=-fsanitize=thread
endif

.PHONY:all extra clean depend
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@


all: $(PROG) $(PROG)_avx $(PROG)_sketch_avx $(PROG)_ksw2_avx

extra:all $(PROG_EXTRA)


$(PROG):main.o libminimap2.a
		$(CC) $(CFLAGS) main.o -o $@ -L. -lminimap2 $(LIBS)

$(PROG)_avx:main.o libminimap2_avx.a
		$(CC) $(CFLAGS) main.o -o $@ -L. -lminimap2_avx $(LIBS)

$(PROG)_sketch_avx:main.o libminimap2_sketch_avx.a
		$(CC) $(CFLAGS) main.o -o $@ -L. -lminimap2_sketch_avx $(LIBS)

$(PROG)_ksw2_avx:main.o libminimap2_ksw2_avx.a
		$(CC) $(CFLAGS) main.o -o $@ -L. -lminimap2_ksw2_avx $(LIBS)

libminimap2_ksw2_avx.a:$(OBJS) map_avx.o ksw2_extd2_avx.o sketch.o
	$(AR) -csru $@ $^

libminimap2_avx.a:$(OBJS) sketch_avx.o map_avx.o ksw2_extd2_avx.o
	$(AR) -csru $@ $^

libminimap2_sketch_avx.a:$(OBJS) sketch_avx.o map.o
	$(AR) -csru $@ $^

libminimap2.a:$(OBJS) map.o sketch.o
		$(AR) -csru $@ $^

sdust:sdust.c kalloc.o kalloc.h kdq.h kvec.h kseq.h ketopt.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< kalloc.o -o $@ -lz

sketch_avx.o: sketch.c kvec.h kalloc.h mmpriv.h minimap.h bseq.h kseq.h
	$(CC) -c $(CFLAGS) -mavx512dq $(CPPFLAGS) $(INCLUDES) $< -o $@

map_avx.o: map.c kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h kseq.h ksw2_extd2_avx.h khash.h ksort.h
	$(CC) -c $(CFLAGS) -mavx512bw $(CPPFLAGS) $(INCLUDES) $< -o $@

ksw2_extd2_avx.o: ksw2_extd2_avx.c ksw2.h ksw2_extd2_avx.h
	$(CC) -c $(CFLAGS) -mavx512bw $(CPPFLAGS) $(INCLUDES) $< -o $@

# SSE-specific targets on x86/x86_64

ifeq ($(arm_neon),)   # if arm_neon is defined, compile this target with the default setting (i.e. no -msse2)
ksw2_ll_sse.o:ksw2_ll_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 $(CPPFLAGS) $(INCLUDES) $< -o $@
endif

ksw2_extz2_sse41.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extz2_sse2.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

exact_match_sse41.o:exact_match_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

exact_match_sse2.o:exact_match_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_extd2_sse41.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extd2_sse2.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_exts2_sse41.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_exts2_sse2.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_dispatch.o:ksw2_dispatch.c ksw2.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

# NEON-specific targets on ARM

ksw2_extz2_neon.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

ksw2_extd2_neon.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

ksw2_exts2_neon.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

exact_match_neon.o:exact_match_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

# other non-file targets

clean:
		rm -fr gmon.out *.o  $(PROG) $(PROG)_avx $(PROG)_sketch_avx $(PROG)_ksw2_avx $(PROG_EXTRA) *~ *.a *.dSYM build dist mappy*.so mappy.c python/mappy.c mappy.egg*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

align.o: minimap.h mmpriv.h bseq.h kseq.h ksw2.h kalloc.h
bseq.o: bseq.h kvec.h kalloc.h kseq.h
esterr.o: mmpriv.h minimap.h bseq.h kseq.h
example.o: minimap.h kseq.h
format.o: kalloc.h mmpriv.h minimap.h bseq.h kseq.h
hit.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h khash.h
index.o: kthread.h bseq.h minimap.h mmpriv.h kseq.h kvec.h kalloc.h khash.h
index.o: ksort.h
kalloc.o: kalloc.h
ksw2_extd2_sse.o: ksw2.h kalloc.h
ksw2_exts2_sse.o: ksw2.h kalloc.h
ksw2_extz2_sse.o: ksw2.h kalloc.h
exact_match_sse.o: ksw2.h kalloc.h
ksw2_ll_sse.o: ksw2.h kalloc.h
kthread.o: kthread.h
lchain.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h krmq.h
main.o: bseq.h minimap.h mmpriv.h kseq.h ketopt.h
map.o: kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h kseq.h ksw2_extd2_avx.h
map.o: khash.h ksort.h
misc.o: mmpriv.h minimap.h bseq.h kseq.h ksort.h
options.o: mmpriv.h minimap.h bseq.h kseq.h
pe.o: mmpriv.h minimap.h bseq.h kseq.h kvec.h kalloc.h ksort.h
sdust.o: kalloc.h kdq.h kvec.h sdust.h
seed.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h ksort.h
sketch.o: kvec.h kalloc.h mmpriv.h minimap.h bseq.h kseq.h
splitidx.o: mmpriv.h minimap.h bseq.h kseq.h
profile.o: profile.h

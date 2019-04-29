# Makefile
# for compiling fdk_mex

# common compile flags
aflags = -std=c99 \
	-O3 -ffast-math -fexpensive-optimizations \
        -Wall -Wshadow -Wpointer-arith -W -Werror -Wmissing-prototypes \
        -DUse_thread \
	-DUse_fdk_mex

# OS specific compile flags
ifeq ($(Arch), mi64) # mac 64-bit
	oflags = -arch x86_64 -DUse_ncore_sysctl # -DUse_aff_mac1
	libs =
	mextop = ~fessler/l/matlab64/extern
	mexinc = -I$(mextop)/include
	mexsuf = maci64
	mexgnu = -DMATLAB_HOST -DMmex
	mexcom = $(mextop)/../bin/mex CFLAGS='$(aflags) $(oflags) $(mexgnu)'

else # linux 64 bit
	oflags = -DNeed_uint
	libs = -lpthread -lm
	# user may need to edit the following line to set up path to matlab:
	mextop = /usr/local/matlab/extern
	mexinc = -I$(mextop)/include
	mexsuf = a64
	mexgnu = -D_GNU_SOURCE -UNeed_uint -fPIC -DMATLAB_HOST -DMmex
	mexcom = $(mextop)/../bin/mex CFLAGS='$(aflags) $(oflags) $(mexgnu)'

endif

CFLAGS = $(aflags) $(oflags)

goal = fdk_mex.mex$(mexsuf)


all:
	@echo 'choose "get" (jf only) or "fdk" to make $(goal)'

#
# jf: copy files from master source directories to these local copies
#
home = $(HOME)
jf_cbc = $(home)/l/src/util/cbct/
jf_umx = $(home)/l/src/util/mex/
get:
	cp -pf \
	\
		$(jf_cbc)/cbct,def.h \
		$(jf_cbc)/cbct,mask2.c \
	\
		$(jf_umx)/def,fdk.h \
		$(jf_umx)/fdk,mex.c \
		$(jf_umx)/fdk,st.c \
		$(jf_umx)/fdk,ts.c \
		$(jf_umx)/fdk,ts,t.c \
	.

# compile it all together

def = ../def
incs = -I. -I$(def)

allc = ./*.c \
		$(def)/mexarg.c \
		$(def)/jf,thread1.c
allh = ./*.h \
		$(def)/defs-env.h \
		$(def)/def,mexarg.h \
		$(def)/jf,mex,def.h \
		$(def)/jf,thread1.h

fdk:	$(goal)

$(goal):	$(allh) $(allc)
	$(mexcom) -o $@ $(allc) $(incs) $(libs)

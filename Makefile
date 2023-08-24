# Examples:
# make jumpmat
# make NBITS=512

ifndef NBITS
   $(info WARNING: NBITS not defined. Using default value: 128)
   NBITS=128
endif
$(info NBITS: $(NBITS))

ifndef TESTU01_DIR
   $(info WARNING: TESTU01_DIR not defined. Using default value: ../testu01/install)
   TESTU01_DIR=../testu01/install
endif

ifdef MKLROOT
    $(info MKLROOT: $(MKLROOT))
endif

ifneq ("$(wildcard $(TESTU01_DIR)/include/TestU01.h)","")
    TESTU01_AVAIL = 1
    $(info TESTU01_DIR: $(TESTU01_DIR))
else
    TESTU01_AVAIL = 0
    $(info TestU01.h header file NOT found)
endif

PLATFORM := $(shell uname -s)
$(info PLATFORM: $(PLATFORM))

ifeq ($(NBITS), 512)
   SIMD=-mavx512f -mavx512bw -mavx512dq
else ifeq ($(NBITS), 256)
   SIMD=-mavx2
else ifeq ($(NBITS), 128)
   SIMD=-msse4.2
endif

BINDIR=bin-$(NBITS)
LOGDIR=logs/testu01

COMMONFLAGS = -c -O3 $(SIMD)

SFMT_FLAGS = -DSFMT_MEXP=19937 -DHAVE_SSE2

CFLAGS += $(COMMONFLAGS)
CPPFLAGS += $(COMMONFLAGS) -O3 -std=c++17 -Iinclude $(SIMD)

HEADERS := $(wildcard include/*.h)
$(info HEADERS: $(HEADERS))


CPP_SRC=$(wildcard src/*.cpp)
$(info C++ files: $(CPP_SRC))

CPP_OBJ=$(patsubst src/%.cpp,$(BINDIR)/%.cpp.obj,$(CPP_SRC))
$(info C++ obj: $(CPP_OBJ))

DATDIR=dat
POWERS=00009 00100 19933 19934 19935 19936
MT_JUMP_7Z=$(DATDIR)/mt/F19937.7z $(patsubst %,$(DATDIR)/mt/F%.7z,$(POWERS))
SFMT_JUMP_7Z=$(patsubst %,$(DATDIR)/sfmt/F%.7z,$(POWERS))
ALL_JUMP_7Z=$(MT_JUMP_7Z) $(SFMT_JUMP_7Z)
JUMP_TARGETS=$(patsubst %.7z,%.bits,$(ALL_JUMP_7Z))
$(info JUMP MATRIX FILES: $(JUMP_TARGETS))

MT_OBJ = $(BINDIR)/mt19937ar.c.obj
SFMT_OBJ = $(BINDIR)/SFMT.c.obj

TARGETS := $(patsubst src/%.cpp,$(BINDIR)/%.exe,$(CPP_SRC))
ifeq ($(TESTU01_AVAIL), 0)
    TARGETS := $(filter-out $(BINDIR)/testu01.exe, $(TARGETS))
endif
$(info TARGETS: $(TARGETS))

all: $(TARGETS)

matrix : $(JUMP_TARGETS)

%.bits : %.7z
	7za e -o$(@D) -y $< > /dev/null
	touch $@

#dat/%.hmat : dat/%.bits $(BINDIR)/encoder.exe
#	$(BINDIR)/encoder.exe -i $< -o $@

# extra compilation flags specific files
$(BINDIR)/perf.cpp.obj $(BINDIR)/test.cpp.obj : CPPFLAGS += $(SFMT_FLAGS)
ifdef MKLROOT
    $(BINDIR)/perf.cpp.obj : CPPFLAGS += -DTEST_MKL -I$(MKLROOT)/include
endif
$(BINDIR)/testu01.cpp.obj : CPPFLAGS += -I$(TESTU01_DIR)/include

$(BINDIR)/%.cpp.obj : src/%.cpp $(HEADERS) Makefile | $(BINDIR)
	g++ $(CPPFLAGS) -o $@ $<

$(MT_OBJ) : mt19937-original/mt19937ar.c Makefile | $(BINDIR)
	gcc $(CFLAGS) -o $@ $<

$(SFMT_OBJ) : CFLAGS += -DHAVE_SSE2
$(SFMT_OBJ) : SFMT-src-1.5.1/SFMT.c Makefile | $(BINDIR)
	gcc $(CFLAGS) $(SFMT_FLAGS) -o $@ $<

# extra dependencies and flags for specific executable
$(BINDIR)/test.exe $(BINDIR)/perf.exe : $(MT_OBJ) $(SFMT_OBJ)
$(BINDIR)/testu01.exe :	LFLAGS += -L$(TESTU01_DIR)/lib -ltestu01 -lprobdist -lmylib -lm
ifdef MKLROOT
    $(BINDIR)/perf.exe : LFLAGS += -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif

$(BINDIR)/%.exe : $(BINDIR)/%.cpp.obj
	g++ -o $@ $^ $(LFLAGS)

BITS=32 128 256 512
TESTLOGS=$(patsubst %,$(LOGDIR)/SmallCrush-%.log,$(BITS)) $(patsubst %,$(LOGDIR)/Crush-%.log,$(BITS)) $(patsubst %,$(LOGDIR)/BigCrush-%.log,$(BITS))
$(info TESTLOGS: $(TESTLOGS))

testu01logs: $(TESTLOGS)

$(LOGDIR)/SmallCrush-%.log : | $(BINDIR)/testu01.exe
	$(BINDIR)/testu01.exe -b $* -m 0 > $@

$(LOGDIR)/Crush-%.log : | $(BINDIR)/testu01.exe
	$(BINDIR)/testu01.exe -b $* -m 1 > $@

$(LOGDIR)/BigCrush-%.log : | $(BINDIR)/testu01.exe
	$(BINDIR)/testu01.exe -b $* -m 2 > $@

.PHONY: clean
clean:
	rm -rf bin-*

# use -p for multithreading
$(BINDIR):
	mkdir -p $(BINDIR)

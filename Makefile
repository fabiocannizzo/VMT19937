# Examples:
# - make NBITS=512  (default to 128)
# - make TESTU01_DIR=/path/to/testu01/install (default to ../testu01/install)
# - make MKLROOT=/path/to/mkl (default to /opt/intel/oneapi/mkl/latest/)

ifndef NBITS
   $(info WARNING: NBITS not defined. Using default value: 128)
   NBITS=128
endif
$(info NBITS: $(NBITS))

ifndef MKLROOT
   $(info WARNING: MKLROOT not defined: using default path.)
   MKLROOT=/opt/intel/oneapi/mkl/latest/
endif

# Check if the directory exists
ifeq ($(wildcard $(MKLROOT)),)
   # Code to run if the directory does NOT exist
   $(info The MKL include directory was not found at $(MKLROOT). Disabling MKL.)
   MKLROOT :=
else
   $(info MKLROOT: $(MKLROOT))
   $(info NOTE: reemmber to export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(MKLROOT))
endif

ifndef TESTU01_DIR
   $(info WARNING: TESTU01_DIR not defined: using default path.)
   TESTU01_DIR=../testu01/install
endif

ifneq ("$(wildcard $(TESTU01_DIR)/include/TestU01.h)","")
    TESTU01_AVAIL = 1
    $(info TESTU01_DIR: $(TESTU01_DIR))
else
    TESTU01_AVAIL = 0
    $(info TestU01.h header file NOT found at $(TESTU01_DIR)/include/)
endif

PLATFORM := $(shell uname -s)
$(info PLATFORM: $(PLATFORM))

CYGWIN := $(findstring CYGWIN,$(PLATFORM))

CC:=gcc
CXX:=g++

$(info CXX: $(CXX))
$(info CC: $(CC))

ifeq ($(NBITS), 512)
   SIMD=-mavx512f -mavx512bw -mavx512dq
else ifeq ($(NBITS), 256)
   SIMD=-mavx2
else ifeq ($(NBITS), 128)
   SIMD=-msse4.2
endif

BINDIR=bin-$(NBITS)-$(CC)
$(info BINDIR: $(BINDIR))

LOGDIR=logs/testu01

COMMONFLAGS = -c -O3 $(SIMD)

SFMT_FLAGS = -DSFMT_MEXP=19937 -DHAVE_SSE2

# clear flags
CFLAGS :=
CPPFLAGS :=

CFLAGS += $(COMMONFLAGS)
CPPFLAGS += $(COMMONFLAGS) -O3 -std=c++20 -Iinclude

HEADERS := $(wildcard include/*.h)
$(info HEADERS: $(HEADERS))


CPP_SRC=$(wildcard src/*.cpp)
$(info C++ files: $(CPP_SRC))

CPP_WITH_MAIN=$(shell grep -l "int main" $(CPP_SRC))
CPP_WITHOUT_MAIN=$(filter-out $(CPP_WITH_MAIN), $(CPP_SRC))

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

TARGETS := $(patsubst src/%.cpp,$(BINDIR)/%.exe,$(CPP_WITH_MAIN))
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
ifneq ($(MKLROOT),)
     $(BINDIR)/perf.cpp.obj : CPPFLAGS += -I$(MKLROOT)/include/
endif
$(BINDIR)/testu01.cpp.obj : CPPFLAGS += -I$(TESTU01_DIR)/include

$(BINDIR)/%.cpp.obj : src/%.cpp $(HEADERS) Makefile | $(BINDIR)
	$(CXX) $(CPPFLAGS) -o $@ $<

$(MT_OBJ) : mt19937-original/mt19937ar.c Makefile | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $<

$(SFMT_OBJ) : SFMT-src-1.5.1/SFMT.c Makefile | $(BINDIR)
	$(CC) $(CFLAGS) $(SFMT_FLAGS) -o $@ $<

# extra dependencies and flags for specific executable
$(BINDIR)/test.exe $(BINDIR)/perf.exe : $(MT_OBJ) $(SFMT_OBJ)
$(BINDIR)/perf.exe : $(BINDIR)/cpu.cpp.obj
$(BINDIR)/testu01.exe :	LFLAGS += -L$(TESTU01_DIR)/lib -ltestu01 -lprobdist -lmylib -lm
ifneq ($(MKLROOT),)
#    $(BINDIR)/perf.exe : LFLAGS += -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
    $(BINDIR)/perf.exe : LFLAGS += -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
endif

$(BINDIR)/%.exe : $(BINDIR)/%.cpp.obj
	$(CXX) -o $@ $^ $(LFLAGS)

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

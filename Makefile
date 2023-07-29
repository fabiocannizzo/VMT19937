ifndef $(NBITS)
   $(info WARNING: NBITS not defined. Using default value: 128)
   NBITS=128
endif

ifndef $(TESTU01_DIR)
   $(info WARNING: TESTU01_DIR not defined. Using default value: ../testu01/install)
   TESTU01_DIR=../testu01/install
endif

ifneq ("$(wildcard $(TESTU01_DIR)/include/TestU01.h)","")
    TESTU01_AVAIL = 1
    $(info TestU01.h header file found)
else
    TESTU01_AVAIL = 0
    $(info TestU01.h header file NOT found)
endif

$(info NBITS: $(NBITS))
$(info TESTU01_DIR: $(TESTU01_DIR))


ifeq ($(NBITS), 512)
   SIMD=-mavx512f -mavx512dq
else ifeq ($(NBITS), 256)
   SIMD=-mavx2
else ifeq ($(NBITS), 128)
   SIMD=-msse4.2
endif

BINDIR=bin-$(NBITS)

COMMONFLAGS = -c -O3 $(SIMD)

SFMT_FLAGS = -DSFMT_MEXP=19937

CFLAGS += $(COMMONFLAGS)
CPPFLAGS += $(COMMONFLAGS) -O3 -std=c++17 -Iinclude $(SIMD)

HEADERS := $(wildcard include/*.h)
$(info HEADERS: $(HEADERS))


CPP_SRC=$(wildcard src/*.cpp)
$(info C++ files: $(CPP_SRC))

CPP_OBJ=$(patsubst src/%.cpp,$(BINDIR)/%.cpp.obj,$(CPP_SRC))
$(info C++ obj: $(CPP_OBJ))

MT_OBJ = $(BINDIR)/mt19937ar.c.obj
SFMT_OBJ = $(BINDIR)/sfmt.c.obj

TARGETS := $(patsubst src/%.cpp,$(BINDIR)/%.exe,$(CPP_SRC))
ifeq ($(TESTU01_AVAIL), 0)
    TARGETS := $(filter-out src/testu01.cpp, $(TARGETS))
endif
$(info TARGETS: $(TARGETS))

all: $(TARGETS)

dat/%.bits : dat/%.7z
	7za e -odat -y $< > /dev/null
	touch $@

dat/%.hmat : dat/%.bits $(BINDIR)/encoder.exe
	$(BINDIR)/encoder.exe -i $< -o $@

# extra compilation flags for perf.cpp
$(BINDIR)/perf.cpp.obj : CPPFLAGS += $(SFMT_FLAGS)
$(BINDIR)/testu01.cpp.obj : CPPFLAGS += -I$(TESTU01_DIR)/include

$(BINDIR)/%.cpp.obj : src/%.cpp $(HEADERS) Makefile | $(BINDIR)
	g++ $(CPPFLAGS) -o $@ $<

$(MT_OBJ) : mt19937-original/mt19937ar.c Makefile | $(BINDIR)
	gcc $(CFLAGS) -o $@ $<

$(SFMT_OBJ) : SFMT-src-1.5.1/sfmt.c Makefile | $(BINDIR)
	gcc $(CFLAGS) $(SFMT_FLAGS) -o $@ $<

# extra dependencies for test.exe
$(BINDIR)/test.exe : | dat/F00010.bits dat/F19937.bits
$(BINDIR)/test.exe : $(MT_OBJ) $(SFMT_OBJ)
$(BINDIR)/perf.exe : $(MT_OBJ) $(SFMT_OBJ)
$(BINDIR)/testu01.exe : | dat/F19933.bits dat/F19934.bits dat/F19935.bits
$(BINDIR)/testu01.exe :	LFLAGS += -L$(TESTU01_DIR)/lib -ltestu01 -lprobdist -lmylib -lm


$(BINDIR)/%.exe : $(BINDIR)/%.cpp.obj
	g++ -o $@ $^ $(LFLAGS)


.PHONY: clean
clean:
	rm -rf bin-*

# use -p for multithreading
$(BINDIR):
	mkdir -p $(BINDIR)

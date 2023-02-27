INCLUDE_DIRS=include

NBITS ?= 128
$(info NBITS: $(NBITS))

ifeq ($(NBITS), 512)
   SIMD=-mavx512f -mavx512dq
else ifeq ($(NBITS), 256)
   SIMD=-mavx2
else ifeq ($(NBITS), 128)
   SIMD=-msse4.2
endif

BINDIR=bin-$(NBITS)

COMMONFLAGS = -c -O3 $(SIMD) -I$(INCLUDE_DIRS)

CFLAGS += $(COMMONFLAGS)
CPPFLAGS += $(COMMONFLAGS) -O3 $(SIMD)


HEADERS := $(shell find . -name "*.h")
$(info HEADERS: $(HEADERS))


CPP_SRC=$(wildcard src/*.cpp)
$(info C++ files: $(CPP_SRC))

CPP_OBJ=$(patsubst src/%.cpp,$(BINDIR)/%.cpp.obj,$(CPP_SRC))
$(info C++ obj: $(CPP_OBJ))

MT_OBJ = $(BINDIR)/mt19937ar.c.obj

TARGETS := $(patsubst src/%.cpp,$(BINDIR)/%.exe,$(CPP_SRC))
$(info TARGETS: $(TARGETS))

all: $(TARGETS)

$(BINDIR)/%.cpp.obj : src/%.cpp $(HEADERS) Makefile $(BINDIR)
	g++ $(CPPFLAGS) -o $@ $<

$(MT_OBJ) : mt19937-original/mt19937ar.c $(HEADERS) Makefile $(BINDIR)
	gcc $(CFLAGS) -o $@ $<

$(BINDIR)/%.exe : $(BINDIR)/%.cpp.obj $(HEADERS) Makefile $(MT_OBJ) $(BINDIR)
	g++ $(LFLAGS) -o $@ $(MT_OBJ) $<


.PHONY: clean
clean:
	rm -rf bin-*

# use -p for multithreading
$(BINDIR):
	mkdir -p $(BINDIR)
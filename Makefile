# Examples:
# - make NBITS=512  (default to 128)
# - make TESTU01_DIR=/path/to/testu01/install (default to ../testu01/install)
# - make MKLROOT=/path/to/mkl (default to automatic discovery)

ifndef NBITS
   NBITS := native
endif
$(info NBITS: $(NBITS))

# Detect compiler
ifeq ($(CXX),)
    CXX := g++
endif
ifeq ($(CC),)
    CC := gcc
endif
export CXX
export CC

# Check if we are using MSVC (cl.exe)
IS_MSVC := $(findstring cl,$(CXX))
ifeq ($(IS_MSVC),)
    # Check if cl is available in path if CXX is just 'cl' or empty
    ifeq ($(shell where cl.exe 2>NUL),)
        IS_MSVC :=
    else
        IS_MSVC := 1
        CXX := cl
        CC := cl
    endif
endif

# --- Platform & Compiler Specific Flags ---
ifeq ($(IS_MSVC),1)
    $(info Compiler: MSVC)
    COMMON_FLAGS := /O2 /MD /EHsc /Zi
    CXX_ONLY_FLAGS := /std:c++20
    C_FLAG := /c
    I_FLAG := /I
    D_FLAG := /D
    OBJ_EXT := .obj
    EXE_EXT := .exe
    OUT_OBJ := /Fo:
    OUT_EXE := /Fe:

    ifeq ($(NBITS), native)
        $(info WARNING: NBITS=native not supported with MSVC. Defaulting to 128)
        NBITS := 128
    endif
    ifeq ($(NBITS), 512)
        SIMD := /arch:AVX512
    else ifeq ($(NBITS), 256)
        SIMD := /arch:AVX2
    else ifeq ($(NBITS), 128)
        SIMD :=
    endif

    # MKL Discovery for MSVC
    ifndef MKLROOT
        # Try to find mkl.h in INCLUDE path or common locations
        MKL_H_FOUND := $(shell where mkl.h 2>NUL)
        ifneq ($(MKL_H_FOUND),)
            $(info MKL found in PATH/INCLUDE)
            MKL_AVAIL := 1
        else
            # Try vcpkg default path if MKLROOT not set
            VCPKG_MKL := <vcpkg-root>/installed/x64-windows
            ifneq ("$(wildcard $(VCPKG_MKL)/include/mkl.h)","")
                MKLROOT := $(VCPKG_MKL)
                MKL_AVAIL := 1
            endif
        endif
    endif

    ifneq ($(MKLROOT),)
        MKL_AVAIL := 1
        MKL_INC := /I"$(MKLROOT)/include"
        MKL_LIB_DIR := /LIBPATH:"$(MKLROOT)/lib"
        MKL_LIBS := mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib Advapi32.lib
    endif

    LFLAGS := /link Advapi32.lib
    SFMT_FLAGS := /D SFMT_MEXP=19937 /D HAVE_SSE2
else
    $(info Compiler: GCC/Clang)
    # Detect Architecture
    ARCH ?= $(shell uname -m)
    $(info Architecture: $(ARCH))

    COMMON_FLAGS := -O3 -pthread
    CXX_ONLY_FLAGS := -std=c++20
    C_FLAG := -c
    I_FLAG := -I
    D_FLAG := -D
    OBJ_EXT := .o
    EXE_EXT := .exe
    OUT_OBJ := -o
    OUT_EXE := -o

    LFLAGS := -pthread
    ifeq ($(ARCH), aarch64)
        # Check if userland is 32-bit
        USERLAND_BITS ?= $(shell getconf LONG_BIT)
        ifeq ($(USERLAND_BITS), 32)
            SIMD := -mfpu=neon -mfloat-abi=hard
        else
            SIMD := -march=armv8-a+simd
        endif
        SFMT_FLAGS := -DSFMT_MEXP=19937 -DHAVE_NEON
        EXTRA_CXXFLAGS += -Wno-psabi
    else ifeq ($(ARCH), armv7l)
        SIMD := -mfpu=neon -mfloat-abi=hard
        SFMT_FLAGS := -DSFMT_MEXP=19937 -DHAVE_NEON
        EXTRA_CXXFLAGS += -Wno-psabi
    else
        ifeq ($(NBITS), native)
            SIMD := -march=native
        else ifeq ($(NBITS), 512)
            SIMD := -mavx512f -mavx512bw -mavx512dq
        else ifeq ($(NBITS), 256)
            SIMD := -mavx2
        else ifeq ($(NBITS), 128)
            SIMD := -msse4.2
        endif
        SFMT_FLAGS := -DSFMT_MEXP=19937 -DHAVE_SSE2
    endif

    # MKL Discovery for GCC
    ifndef MKLROOT
        MKLROOT := /opt/intel/oneapi/mkl/latest
    endif

    ifneq ("$(wildcard $(MKLROOT)/include/mkl.h)","")
        MKL_AVAIL := 1
        MKL_INC := -I$(MKLROOT)/include
        ifneq ("$(wildcard $(MKLROOT)/lib/intel64)","")
            MKL_LIB_DIR := -L$(MKLROOT)/lib/intel64
        else
            MKL_LIB_DIR := -L$(MKLROOT)/lib
        endif
        MKL_LIBS := -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
    endif
endif
# --- TestU01 Discovery ---
ifndef TESTU01_DIR
   TESTU01_DIR := ../testu01/install
endif

ifneq ("$(wildcard $(TESTU01_DIR)/include/TestU01.h)","")
    TESTU01_AVAIL := 1
    ifeq ($(IS_MSVC),1)
        TESTU01_INC := /I"$(TESTU01_DIR)/include"
        TESTU01_LIB_DIR := /LIBPATH:"$(TESTU01_DIR)/lib"
        TESTU01_LIBS := testu01.lib probdist.lib mylib.lib
    else
        TESTU01_INC := -I$(TESTU01_DIR)/include
        TESTU01_LIB_DIR := -L$(TESTU01_DIR)/lib
        TESTU01_LIBS := -ltestu01 -lprobdist -lmylib -lm
    endif
endif

# --- Paths & Files ---
BINDIR := bin-$(NBITS)-$(notdir $(CXX))
$(info BINDIR: $(BINDIR))

CPPFLAGS := $(I_FLAG)include
CXXFLAGS := $(COMMON_FLAGS) $(CXX_ONLY_FLAGS) $(SIMD) $(EXTRA_CXXFLAGS)
CFLAGS := $(COMMON_FLAGS) $(SIMD)
ifneq ($(IS_MSVC),1)
    CXXFLAGS += -MMD -MP
    CFLAGS   += -MMD -MP
endif

MAKEFILE_DEPS := Makefile

CPP_SRC := $(wildcard src/*.cpp)
# We identify files with main() to determine targets
# On Windows/MSVC, grep might not be available, so we assume these files have main
CPP_WITH_MAIN := src/perf.cpp src/test.cpp src/demo.cpp src/encoder.cpp src/jump.cpp src/testu01.cpp
CPP_WITHOUT_MAIN := src/cpu.cpp

MT_OBJ := $(BINDIR)/mt19937ar$(OBJ_EXT)
SFMT_OBJ := $(BINDIR)/SFMT$(OBJ_EXT)
CPU_OBJ := $(BINDIR)/cpu$(OBJ_EXT)

TARGETS := $(patsubst src/%.cpp,$(BINDIR)/%$(EXE_EXT),$(CPP_WITH_MAIN))
ifneq ($(TESTU01_AVAIL), 1)
    TARGETS := $(filter-out $(BINDIR)/testu01$(EXE_EXT), $(TARGETS))
endif

# --- Rules ---

all: $(TARGETS)

$(BINDIR):
	mkdir -p $(BINDIR)

$(MT_OBJ): mt19937-original/mt19937ar.c $(MAKEFILE_DEPS) | $(BINDIR)
	$(CC) $(CFLAGS) $(C_FLAG) $(OUT_OBJ)$@ $<

$(SFMT_OBJ): SFMT-src-1.5.1/SFMT.c $(MAKEFILE_DEPS) | $(BINDIR)
	$(CC) $(CFLAGS) $(SFMT_FLAGS) $(C_FLAG) $(OUT_OBJ)$@ $<

$(BINDIR)/%$(OBJ_EXT): src/%.cpp $(MAKEFILE_DEPS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(C_FLAG) $(OUT_OBJ)$@ $<

# Specific flags for objects
$(BINDIR)/perf$(OBJ_EXT) $(BINDIR)/test$(OBJ_EXT): CPPFLAGS += $(SFMT_FLAGS)
ifeq ($(MKL_AVAIL),1)
    $(BINDIR)/perf$(OBJ_EXT): CPPFLAGS += $(MKL_INC)
endif
ifeq ($(TESTU01_AVAIL),1)
    $(BINDIR)/testu01$(OBJ_EXT): CPPFLAGS += $(TESTU01_INC)
endif

# Executables
$(BINDIR)/test$(EXE_EXT): $(BINDIR)/test$(OBJ_EXT) $(MT_OBJ) $(SFMT_OBJ)
	$(CXX) $(OUT_EXE)$@ $^ $(LFLAGS)

$(BINDIR)/perf$(EXE_EXT): $(BINDIR)/perf$(OBJ_EXT) $(MT_OBJ) $(SFMT_OBJ) $(CPU_OBJ)
	$(CXX) $(OUT_EXE)$@ $^ $(LFLAGS) $(MKL_LIB_DIR) $(MKL_LIBS)

$(BINDIR)/testu01$(EXE_EXT): $(BINDIR)/testu01$(OBJ_EXT)
	$(CXX) $(OUT_EXE)$@ $^ $(LFLAGS) $(TESTU01_LIB_DIR) $(TESTU01_LIBS)

$(BINDIR)/%$(EXE_EXT): $(BINDIR)/%$(OBJ_EXT)
	$(CXX) $(OUT_EXE)$@ $^ $(LFLAGS)

.PHONY: clean
clean:
	rm -rf bin-*

-include $(wildcard $(BINDIR)/*.d)

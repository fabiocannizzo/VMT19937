# Build examples:
#   make                              # NBITS=native, ISA auto-detected
#   make NBITS=128                    # SSE4.2 (x86) or NEON (ARM)
#   make NBITS=256                    # AVX2 (x86) or SVE-256 (ARM)
#   make NBITS=512                    # AVX-512+VL (x86)
#   make NBITS=256 ISA=avx512vl       # 256-bit virtual regs with AVX-512VL instructions
#   make CXX=cl NBITS=128             # MSVC SSE4.2
#   make TESTU01_DIR=/path/to/testu01/install
#   make MKLROOT=/path/to/mkl

ifndef NBITS
   NBITS := native
endif
DEBUG ?= 0
$(info DEBUG: $(DEBUG))

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
ifneq ($(findstring cl,$(CXX)),)
    IS_MSVC := 1
    CC := cl
else ifeq ($(shell where cl.exe 2>/dev/null),)
    IS_MSVC :=
else
    IS_MSVC := 1
    CXX := cl
    CC := cl
endif

# ============================================================
# Platform & Compiler Specific Flags
# ============================================================

ifeq ($(IS_MSVC),1)
    $(info Compiler: $(shell cl 2>&1 | head -1))

    ifeq ($(DEBUG),1)
        COMMON_FLAGS := /Od /MDd /EHsc /Zi /FS
    else
        COMMON_FLAGS := /O2 /MD /EHsc /Zi /FS
    endif
    CXX_ONLY_FLAGS := /std:c++20
    C_FLAG    := /c
    I_FLAG    := /I
    D_FLAG    := /D
    OBJ_EXT   := .obj
    EXE_EXT   := .exe
    OUT_OBJ   := /Fo:
    OUT_EXE   := /Fe:
    PDB_FLAGS  = /Fd:$(basename $@).pdb

    # Native ISA detection: compile and run cpuid_probe.c using cl.exe
    ifeq ($(NBITS), native)
        _ISA_PROBE := $(strip $(shell \
            cl /nologo /W0 /O2 /Fo:cpuid_probe_tmp.obj /Fe:cpuid_probe_tmp.exe \
               src/cpuid_probe.c >/dev/null 2>&1 \
            && ./cpuid_probe_tmp.exe 2>/dev/null || echo sse42; \
            rm -f cpuid_probe_tmp.obj cpuid_probe_tmp.exe >/dev/null 2>&1))
        ISA ?= $(_ISA_PROBE)
    else
        ifeq ($(NBITS), 512)
            ISA ?= avx512vl
        else ifeq ($(NBITS), 256)
            ISA ?= avx2
        else
            ISA ?= sse42
        endif
    endif

    # ISA → MSVC /arch flag
    ifeq ($(ISA), avx512vl)
        SIMD := /arch:AVX512
    else ifeq ($(ISA), avx2)
        SIMD := /arch:AVX2
    else
        SIMD :=
    endif

    # SIMD_N_BITS_VAL: for native derive from ISA; for explicit NBITS use NBITS
    ifeq ($(NBITS), native)
        ifeq ($(ISA), avx512vl)
            SIMD_N_BITS_VAL := 512
        else ifeq ($(ISA), avx2)
            SIMD_N_BITS_VAL := 256
        else
            SIMD_N_BITS_VAL := 128
        endif
    else
        SIMD_N_BITS_VAL := $(NBITS)
    endif

    LFLAGS := /link Advapi32.lib
    ifeq ($(DEBUG),1)
        LFLAGS += /DEBUG
    endif
    SFMT_FLAGS := /D SFMT_MEXP=19937 /D HAVE_SSE2

    # MKL Discovery for MSVC
    ifndef MKLROOT
        MKL_H_FOUND := $(shell where mkl.h 2>/dev/null)
        ifneq ($(MKL_H_FOUND),)
            $(info MKL found in PATH/INCLUDE)
            MKL_AVAIL := 1
        else
            VCPKG_MKL := <vcpkg-root>/installed/x64-windows
            ifneq ("$(wildcard $(VCPKG_MKL)/include/mkl.h)","")
                MKLROOT := $(VCPKG_MKL)
                MKL_AVAIL := 1
            endif
        endif
    endif
    ifneq ($(MKLROOT),)
        MKL_AVAIL := 1
        MKL_INC     := /I"$(MKLROOT)/include"
        MKL_LIB_DIR := /LIBPATH:"$(MKLROOT)/lib"
        MKL_LIBS    := mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib Advapi32.lib
    endif

else
    # -------- GCC / Clang --------
    $(info Compiler: $(shell $(CXX) --version | head -1))
    ARCH ?= $(shell uname -m)
    $(info Architecture: $(ARCH))

    ifeq ($(DEBUG),1)
        COMMON_FLAGS := -O0 -g
    else
        COMMON_FLAGS := -O3
    endif
    CXX_ONLY_FLAGS := -std=c++20
    C_FLAG    := -c
    I_FLAG    := -I
    D_FLAG    := -D
    OBJ_EXT   := .o
    EXE_EXT   := .exe
    OUT_OBJ   := -o
    OUT_EXE   := -o
    LFLAGS    :=

    ifeq ($(NBITS), native)
        # Probe the compiler to identify the native ISA
        _DETECT := $(shell echo | $(CXX) -march=native -dM -E -x c++ - 2>/dev/null)
        ifneq (,$(findstring __AVX512VL__,$(_DETECT)))
            ISA            ?= avx512vl
            SIMD_N_BITS_VAL := 512
        else ifneq (,$(findstring __AVX2__,$(_DETECT)))
            ISA            ?= avx2
            SIMD_N_BITS_VAL := 256
        else ifneq (,$(findstring __SSE4_2__,$(_DETECT)))
            ISA            ?= sse42
            SIMD_N_BITS_VAL := 128
        else ifneq (,$(findstring __ARM_FEATURE_SVE,$(_DETECT)))
            ISA            ?= sve256
            SIMD_N_BITS_VAL := 256
        else ifneq (,$(findstring __ARM_NEON,$(_DETECT)))
            ISA            ?= neon
            SIMD_N_BITS_VAL := 128
        else
            ISA            ?= scalar
            SIMD_N_BITS_VAL := 32
        endif
        SIMD := -march=native
    else
        SIMD_N_BITS_VAL := $(NBITS)
        # Derive default ISA from ARCH + NBITS
        ifeq ($(ARCH), aarch64)
            ifeq ($(NBITS), 256)
                ISA ?= sve256
            else
                ISA ?= neon
            endif
        else ifeq ($(ARCH), armv7l)
            ISA ?= neon
        else
            # x86 / x86-64
            ifeq ($(NBITS), 512)
                ISA ?= avx512vl
            else ifeq ($(NBITS), 256)
                ISA ?= avx2
            else
                ISA ?= sse42
            endif
        endif
        # ISA → GCC flags (non-native)
        ifeq ($(ISA), avx512vl)
            SIMD := -mavx512f -mavx512vl -mavx512bw -mavx512dq
        else ifeq ($(ISA), avx2)
            SIMD := -mavx2
        else ifeq ($(ISA), sse42)
            SIMD := -msse4.2
        else ifeq ($(ISA), neon)
            ifeq ($(ARCH), armv7l)
                SIMD := -mfpu=neon -mfloat-abi=hard
            else
                SIMD := -march=armv8-a+simd
            endif
        else ifeq ($(ISA), sve256)
            SIMD := -march=armv8-a+sve
        endif
    endif

    # SFMT flags and ARM-specific extras
    ifeq ($(ISA), neon)
        SFMT_FLAGS := -DSFMT_MEXP=19937 -DHAVE_NEON
        EXTRA_CXXFLAGS += -Wno-psabi
    else ifeq ($(ISA), sve256)
        SFMT_FLAGS := -DSFMT_MEXP=19937 -DHAVE_NEON
        EXTRA_CXXFLAGS += -Wno-psabi
    else
        SFMT_FLAGS := -DSFMT_MEXP=19937 -DHAVE_SSE2
    endif

    # MKL Discovery for GCC
    ifndef MKLROOT
        MKLROOT := /opt/intel/oneapi/mkl/latest
    endif
    ifneq ("$(wildcard $(MKLROOT)/include/mkl.h)","")
        MKL_AVAIL := 1
        MKL_INC   := -I$(MKLROOT)/include
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
        TESTU01_INC     := /I"$(TESTU01_DIR)/include"
        TESTU01_LIB_DIR := /LIBPATH:"$(TESTU01_DIR)/lib"
        TESTU01_LIBS    := testu01.lib probdist.lib mylib.lib
    else
        TESTU01_INC     := -I$(TESTU01_DIR)/include
        TESTU01_LIB_DIR := -L$(TESTU01_DIR)/lib
        TESTU01_LIBS    := -ltestu01 -lprobdist -lmylib -lm
    endif
endif

# --- Paths & Files ---
$(info NBITS: $(if $(filter native,$(NBITS)),native ($(SIMD_N_BITS_VAL)),$(SIMD_N_BITS_VAL)))
$(info ISA: $(ISA))
BINDIR := bin-$(SIMD_N_BITS_VAL)-$(ISA)-$(notdir $(CXX))$(if $(filter 1,$(DEBUG)),-dbg,)
$(info BINDIR: $(BINDIR))

# Always pass SIMD_N_BITS explicitly so simd_config.h #ifndef guard works
CPPFLAGS  := $(I_FLAG)include $(D_FLAG)SIMD_N_BITS=$(SIMD_N_BITS_VAL)
CXXFLAGS  := $(COMMON_FLAGS) $(CXX_ONLY_FLAGS) $(SIMD) $(EXTRA_CXXFLAGS)
CFLAGS    := $(COMMON_FLAGS) $(SIMD)
ifneq ($(IS_MSVC),1)
    CXXFLAGS += -MMD -MP
    CFLAGS   += -MMD -MP
endif

MAKEFILE_DEPS := Makefile

CPP_SRC := $(wildcard src/*.cpp)
# Files with main() — cpuid_probe.c is C only, excluded from C++ targets
CPP_WITH_MAIN    := src/perf.cpp src/test.cpp src/demo.cpp src/encoder.cpp src/jump.cpp src/testu01.cpp
CPP_WITHOUT_MAIN := src/cpu.cpp

MT_OBJ   := $(BINDIR)/mt19937ar$(OBJ_EXT)
MT64_OBJ := $(BINDIR)/mt19937-64$(OBJ_EXT)
SFMT_OBJ := $(BINDIR)/SFMT$(OBJ_EXT)
CPU_OBJ  := $(BINDIR)/cpu$(OBJ_EXT)

TARGETS := $(patsubst src/%.cpp,$(BINDIR)/%$(EXE_EXT),$(CPP_WITH_MAIN))
ifneq ($(TESTU01_AVAIL), 1)
    TARGETS := $(filter-out $(BINDIR)/testu01$(EXE_EXT), $(TARGETS))
endif
ifeq ($(IS_MSVC),1)
    TARGETS := $(filter-out $(BINDIR)/testu01$(EXE_EXT), $(TARGETS))
endif

# --- Rules ---

all: $(TARGETS)

$(BINDIR):
	mkdir -p $(BINDIR)

$(MT_OBJ): mt19937-original/mt19937ar.c $(MAKEFILE_DEPS) | $(BINDIR)
	$(CC) $(CFLAGS) $(if $(IS_MSVC),/FS) $(C_FLAG) $(OUT_OBJ)$@ $(if $(IS_MSVC),$(PDB_FLAGS)) $<

$(MT64_OBJ): mt19937-original/mt19937-64.c $(MAKEFILE_DEPS) | $(BINDIR)
	$(CC) $(CFLAGS) $(if $(IS_MSVC),/FS) $(D_FLAG)MT19937_64_NO_MAIN $(C_FLAG) $(OUT_OBJ)$@ $(if $(IS_MSVC),$(PDB_FLAGS)) $<

$(SFMT_OBJ): SFMT-src-1.5.1/SFMT.c $(MAKEFILE_DEPS) | $(BINDIR)
	$(CC) $(CFLAGS) $(if $(IS_MSVC),/FS) $(SFMT_FLAGS) $(C_FLAG) $(OUT_OBJ)$@ $(if $(IS_MSVC),$(PDB_FLAGS)) $<

$(BINDIR)/%$(OBJ_EXT): src/%.cpp $(MAKEFILE_DEPS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(if $(IS_MSVC),/FS) $(CPPFLAGS) $(C_FLAG) $(OUT_OBJ)$@ $(if $(IS_MSVC),$(PDB_FLAGS)) $<

# Specific flags for objects
ifneq ($(IS_MSVC),1)
$(BINDIR)/jump$(OBJ_EXT): CXXFLAGS += -pthread
$(BINDIR)/jump$(EXE_EXT): LFLAGS   += -pthread
endif
$(BINDIR)/perf$(OBJ_EXT) $(BINDIR)/test$(OBJ_EXT): CPPFLAGS += $(SFMT_FLAGS)
ifeq ($(MKL_AVAIL),1)
    $(BINDIR)/perf$(OBJ_EXT): CPPFLAGS += $(MKL_INC)
endif
ifeq ($(TESTU01_AVAIL),1)
    $(BINDIR)/testu01$(OBJ_EXT): CPPFLAGS += $(TESTU01_INC)
endif

# Executables
$(BINDIR)/test$(EXE_EXT): $(BINDIR)/test$(OBJ_EXT) $(MT_OBJ) $(MT64_OBJ) $(SFMT_OBJ)
	$(CXX) $(OUT_EXE)$@ $^ $(LFLAGS)

$(BINDIR)/perf$(EXE_EXT): $(BINDIR)/perf$(OBJ_EXT) $(MT_OBJ) $(MT64_OBJ) $(SFMT_OBJ) $(CPU_OBJ)
	$(CXX) $(OUT_EXE)$@ $^ $(LFLAGS) $(MKL_LIB_DIR) $(MKL_LIBS)

$(BINDIR)/testu01$(EXE_EXT): $(BINDIR)/testu01$(OBJ_EXT)
	$(CXX) $(OUT_EXE)$@ $^ $(LFLAGS) $(TESTU01_LIB_DIR) $(TESTU01_LIBS)

$(BINDIR)/%$(EXE_EXT): $(BINDIR)/%$(OBJ_EXT)
	$(CXX) $(OUT_EXE)$@ $^ $(LFLAGS)

.PHONY: clean
clean:
	rm -rf bin-*

-include $(wildcard $(BINDIR)/*.d)

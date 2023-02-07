INCLUDE_DIR=./

COMMONFLAGS = -c -O3 -msse4.1 -I$(INCLUDE_DIR)

CFLAGS += $(COMMONFLAGS)
CPPFLAGS += $(COMMONFLAGS) -O3 -msse4.1 -I$(INCLUDE_DIR)


HEADERS := $(shell find $(INCLUDE_DIR) -name "*.h")
$(info HEADERS: $(HEADERS))


C_SRC=$(wildcard *.c)
CPP_SRC=$(wildcard *.cpp)
$(info C files: $(C_SRC))
$(info C++ files: $(CPP_SRC))

C_OBJ=$(patsubst %.c,%.c.obj,$(C_SRC))
CPP_OBJ=$(patsubst %.cpp,%.cpp.obj,$(CPP_SRC))
$(info C obj: $(C_OBJ))
$(info C++ obj: $(CPP_OBJ))


TARGETS := $(patsubst %.cpp,%.exe,$(CPP_SRC))
$(info TARGETS: $(TARGETS))

%.cpp.obj : %.cpp $(HEADERS) Makefile
	g++ $(CPPFLAGS) -o $@ $<

%.c.obj : %.c $(HEADERS) Makefile
	gcc $(CFLAGS) -o $@ $<

%.exe : %.cpp.obj $(C_OBJ) $(HEADERS) Makefile
	g++ $(LFLAGS) -o $@ $(C_OBJ) $<

all: $(TARGETS)

.PHONY: clean
clean:
	rm *.exe *.obj

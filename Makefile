# Written by Ivan Pogrebnyak

# TODO: fix problem with making in multiple threads

CXX := g++
STD := -std=c++1z
DF := $(STD)
CF := $(STD) -Wall -fmax-errors=3 #-Og #-flto
LF := $(STD) #-flto

ROOT_INCDIR := $(shell root-config --incdir)
ROOT_CFLAGS := -pthread -Wno-deprecated-declarations -m64 -I$(ROOT_INCDIR)
ROOT_LIBS   := -L$(shell root-config --libdir) -lCore -lRIO -lHist -lTree -lTreePlayer -pthread -lm -ldl -rdynamic

DF += -I$(ROOT_INCDIR)
CF += $(ROOT_CFLAGS)
LF += $(ROOT_LIBS)

L_select += -lboost_program_options

SRC := src
BIN := bin
BLD := .build

SRCS := $(shell find $(SRC) -type f -name '*.cc')
DEPS := $(patsubst $(SRC)%.cc,$(BLD)%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' $(SRCS)
EXES := $(patsubst $(SRC)%.cc,$(BIN)%,$(shell $(GREP_EXES)))

NODEPS := clean
.PHONY: all clean

all: $(EXES)

#Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DEPS)
endif

$(DEPS): $(BLD)/%.d: $(SRC)/%.cc | $(BLD)
	$(CXX) $(DF) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o: | $(BLD)
	$(CXX) $(CF) $(C_$*) -c $(filter %.cc,$^) -o $@

$(BIN)/%: $(BLD)/%.o | $(BIN)
	$(CXX) $(filter %.o,$^) -o $@ $(LF) $(L_$*)

$(BLD) $(BIN):
	mkdir $@

clean:
	@rm -rfv $(BLD) $(EXES)

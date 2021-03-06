SILENCE_WARNINGS := 0
EXECUTION_PREFIX = 
PROCESS_NUM = 2
ifeq ($(MODE),debug)
	CC = clang++
	CFLAGS := -O0 -g -pg -Weverything -Werror -Wno-error=padded -ferror-limit=1 -Wno-c++98-compat
	LFLAGS := -pg
else ifeq ($(MODE),gcc)
	CC = gcc
	CFLAGS = -O3 -march=native -g -ffast-math -Werror -Wall -Wextra -pedantic -Wno-long-long -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations -Wconversion
else ifeq ($(MODE),parallel)
	EXECUTION_PREFIX = mpirun -n $(PROCESS_NUM)
	CC = mpiCC
	CFLAGS = -DPARALLEL -std=c++11 -Wfatal-errors -O3 -march=native -ffast-math
else ifeq ($(MODE),debug_parallel)
	EXECUTION_PREFIX = mpirun -n $(PROCESS_NUM)
	CC = mpiCC
	CFLAGS = -DPARALLEL -std=c++11 -Wfatal-errors -ggdb -march=native
else
	CC = clang++
	MODE = production
	CFLAGS := -O2 -Weverything -Werror -Wno-error=padded -ferror-limit=1 -ffast-math -Wno-c++98-compat
	ifeq ($(SILENCE_WARNINGS),1)
		CFLAGS := -O2
	endif
	LFLAGS :=
endif
PYTHON_PATH = python
BINARY = $(TEMP_FOLDER)/test.exe
SOURCE_FILES = $(shell find $(FUJIN_ROOT)/source -name '*.*pp')
EMPTY := 
TEST_SOURCE_LIST := $(shell find . -name 'test.cpp')
TEMP1 := $(subst /test.cpp,$(EMPTY),$(TEST_SOURCE_LIST))
TEMP2 := $(subst ./,/,$(TEMP1))
TEMP_FOLDER_LIST := $(addprefix $(MODE),$(TEMP2))
BINARY_LIST := $(addsuffix /test.exe,$(TEMP_FOLDER_LIST))
OBJECT_LIST := $(addsuffix /test.o,$(TEMP_FOLDER_LIST))
RESULT_LIST := $(addsuffix /test_passed.res,$(TEMP_FOLDER_LIST))
NORMAL_TERMINATION_LIST := $(addsuffix /test_terminated_normally.res,$(TEMP_FOLDER_LIST))
TEST_ROOT := $(shell pwd)

.PRECIOUS: $(NORMAL_TERMINATION_LIST) $(RESULT_LIST)

all: $(RESULT_LIST)
	@echo "Done"

$(RESULT_LIST): $(MODE)/%/test_passed.res: $(MODE)/%/test_terminated_normally.res %/test.py
	cd $(MODE)/$* && $(PYTHON_PATH) $(TEST_ROOT)/$*/test.py && tail ./test_passed.res

$(NORMAL_TERMINATION_LIST): $(MODE)/%/test_terminated_normally.res: $(MODE)/%/test.exe
	cd $(MODE)/$* && $(EXECUTION_PREFIX) ./test.exe

$(BINARY_LIST): $(MODE)/%/test.exe: $(MODE)/%/test.o $(FUJIN_ROOT)/library_$(MODE)/libfujin.a
	echo $(LD_LIBRARY_PATH)
	$(CC) $(MODE)/$*/test.o $(FUJIN_ROOT)/library_$(MODE)/libfujin.a -o $(MODE)/$*/test.exe $(LFLAGS) -L $(HDF5_LIB_PATH) -lhdf5 -lhdf5_cpp -lz

$(OBJECT_LIST): $(MODE)/%/test.o: %/test.cpp $(FUJIN_ROOT)/library_$(MODE)/libfujin.a | $(MODE)/%
	$(CC) -c $(CFLAGS) $*/test.cpp -o $(MODE)/$*/test.o -I $(FUJIN_ROOT)/source
	$(CC) -MM $(CFLAGS) $*/test.cpp -o $(MODE)/$*/test.d -I $(FUJIN_ROOT)/source

$(TEMP_FOLDER_LIST): %:
	mkdir -p $*

$(FUJIN_ROOT)/library_$(MODE)/libfujin.a: $(SOURCE_FILES)
	$(MAKE) -C $(FUJIN_ROOT) MODE=$(MODE)

.PHONY: clean all

clean:
	rm -rf $(MODE)

SOURCE_DIR := source
SOURCES := $(shell find $(SOURCE_DIR) -name '*.cpp')
STRIPPED_SOURCES := $(notdir $(SOURCES))
LIB_FILE := libfujin.a
CC := clang++
ifeq ($(MODE),debug)
	CFLAGS := -O0 -g -pg -Weverything -Werror -ferror-limit=1 -Wno-error=padded
	LFLAGS := -pg
else
	MODE = production
	CFLAGS := -O2 -Weverything -Werror -ferror-limit=1 -Wno-error=padded
	LFLAGS :=
endif

LIBRARY_FOLDER := library_$(MODE)
OBJECTS := $(patsubst %.cpp,$(LIBRARY_FOLDER)/%.o,$(STRIPPED_SOURCES))

$(LIBRARY_FOLDER)/$(LIB_FILE): $(OBJECTS)
	ar cr $@ $^

$(OBJECTS): $(LIBRARY_FOLDER)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p $(LIBRARY_FOLDER)
	$(CC) -c $(CFLAGS) $< -o $@
	$(CC) -MM $(CFLAGS) $< -o $(LIBRARY_FOLDER)/$*.d
	@sed 's,\(\w*\)\.o,$(LIBRARY_FOLDER)/\1.o,g' -i $(LIBRARY_FOLDER)/$*.d

-include $(OBJECTS:.o=.d)

.PHONY: clean

clean:
	rm -rf $(LIBRARY_FOLDER)

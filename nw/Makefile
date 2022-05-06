#Commands
CC=gcc
RM=rm -f

#flags
CFLAGS=-Wall -g3 -I include/ -Wall -Wextra
OFLAGS=-O2 -march=native -mtune=native
# OFLAGS=-Ofast -march=native -mtune=native -funroll-loops -finline-functions -ftree-loop-vectorize -ftree-vectorize
LDFLAGS=
# LDFLAGS=-lm

# Config
SRC=src
DEPS=build
TARGET=nw bin

#Files
NW_SOURCES=$(SRC)/alignment.c 
NW_OBJECTS=$(NW_SOURCES:$(SRC)%.c=$(DEPS)%.o)

NWBIN_SOURCES=$(SRC)/bin_alignment.c 
NWBIN_OBJECTS=$(NWBIN_SOURCES:$(SRC)%.c=$(DEPS)%.o)


all: $(TARGET)

$(DEPS)/%.o: $(SRC)/%.c
	@mkdir -p $(DEPS)
	$(CC) $(CFLAGS) $(OFLAGS) -c $< -o $@

nw: $(NW_OBJECTS)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(LDFLAGS)

bin: $(NWBIN_OBJECTS)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(LDFLAGS)

display: $(SRC)/display.c
	$(CC) $(CFLAGS) $? -o $@

clean:
	$(RM) $(NW_OBJECTS) $(NWBIN_OBJECTS)
	$(RM) $(TARGET) $(DEPS)/*

.PHONY: clean all

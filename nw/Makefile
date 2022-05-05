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
TARGET=nw

#Files
LBM_SOURCES=$(SRC)/alignment.c 
LBM_OBJECTS=$(LBM_SOURCES:$(SRC)%.c=$(DEPS)%.o)


all: $(TARGET)

$(DEPS)/%.o: $(SRC)/%.c
	@mkdir -p $(DEPS)
	$(CC) $(CFLAGS) $(OFLAGS) -c $< -o $@

nw: $(LBM_OBJECTS)
	$(CC) $(CFLAGS) $(OFLAGS) $^ -o $@ $(LDFLAGS)

display: $(SRC)/display.c
	$(CC) $(CFLAGS) $? -o $@

clean:
	$(RM) $(LBM_OBJECTS)
	$(RM) $(TARGET) $(DEPS)/*

.PHONY: clean all

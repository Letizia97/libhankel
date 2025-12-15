
# Compiler & flags
CC      = gcc
CFLAGS  = -Wall -Iinclude -fPIC
LDFLAGS =
# Use pkg-config for GSL if present; fallback to explicit libs
LDLIBS  = $(shell pkg-config --libs gsl 2>/dev/null)
ifeq ($(LDLIBS),)
LDLIBS  = -lgsl -lgslcblas -lm
endif

# Directories
SRC_DIR   = src
INC_DIR   = include
TEST_DIR  = test

# Library
LIB_SHARED = liblibhankel.so

# Sources/objects
SRCS  = $(wildcard $(SRC_DIR)/*.c)
OBJS  = $(SRCS:.c=.o)

# Test executable (built from all .c in test/)
TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)
TEST_OBJS = $(TEST_SRCS:.c=.o)
TEST_BIN  = test_libhankel

.PHONY: all lib test clean

# Default: build lib and test app
all: lib test

# --- Library ---
lib: $(LIB_SHARED)

$(LIB_SHARED): $(OBJS)
	$(CC) -shared $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(INC_DIR)/libhankel.h
	$(CC) $(CFLAGS) -c $< -o $@

# --- Test executable ---
test: $(TEST_BIN)

$(TEST_BIN): $(TEST_OBJS) $(LIB_SHARED)
	$(CC) $(CFLAGS) -L. -o $@ $(TEST_OBJS) -llibhankel $(LDLIBS)

$(TEST_DIR)/%.o: $(TEST_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# --- Clean ---
clean:
	rm -f $(OBJS) $(LIB_SHARED) $(TEST_OBJS) $(TEST_BIN)

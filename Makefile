# compiler and compiler flags
CC = gcc
CFLAGS = -O2 -Wall -Wextra -march=native

# output binary
#  EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  EXEDIR = bin
  EXE = $(EXEDIR)/bicubic-interpolation

# source files
  SRC += main.c
  SRC += InputFunction.c

#  MKLINC += `pkg-config --cflags mkl`
  LIB += -lm
#  LIB += `pkg-config --libs mkl`

all: $(EXE) Makefile

$(SRC:.c=.o): $(SRC)
	$(CC) $(CFLAGS) $(MKLINC) $? -c

$(EXE): $(SRC:.c=.o)
	$(CC) $(CFLAGS) $(LIB) $^ -o $(EXE)

clean:
	rm -f $(SRC:.c=.o) $(EXE)

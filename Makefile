  CC=gcc
  LD=${CC} 
  CFLAGS= -Wall -pedantic -std=c99 -O3 -D_GNU_SOURCE_ -g  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DDBGLEVEL=0 -DPROGNFO -I -Lsrc -DPCGRNG
  LDFLAGS=-lm -lpthread 
  CTAGS=ctags > tags
  LIBS=

  METSEGOBJ = src/manopt.o\
	src/memory.o\
	src/debug.o\
	src/info.o\
	src/stringutils.o\
	src/fileio.o\
	src/mathematics.o\
	src/sort.o\
	src/vstack.o\
	src/segmentstack.o\
	src/mtc.o\
	src/metseg.o

all: metilene 

metilene: ${METSEGOBJ}
	gcc $(CFLAGS) ${METSEGOBJ} -o metilene $(LDFLAGS)

clean:
	rm -rf src/*.o metilene

.PHONY: all clean

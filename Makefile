
# change this TARGET to compile your own programs
TARGET = main
SHELL  = /bin/sh

DBGFLAGS   = -g -Wall
CFLAGS     = -O3 -Wall -ffp-contract=off -fno-strict-aliasing
EXTRACFLAGS=
#CFLAGS    = $(DBGFLAGS)

CC         = /usr/bin/gcc
CPPFLAGS   = -I. -I/NOBACKUP/ericahl/pari-lib/include
LD         = /usr/bin/gcc
LDFLAGS    = -O3 -Wall -ffp-contract=off -fno-strict-aliasing    -Wl,--export-dynamic
MODLD      = /usr/bin/gcc
MODLDFLAGS = -shared  $(CFLAGS) $(DLCFLAGS) -Wl,-shared
EXTRAMODLDFLAGS = -lc -lm -L/NOBACKUP/ericahl/pari-lib/lib -lpari
EXTRALIBS  =

RUNPTH     = -Wl,-rpath "/NOBACKUP/ericahl/pari-lib/lib"
DLCFLAGS   = -fPIC
LIBS       = -lm -L/NOBACKUP/ericahl/pari-lib/lib -lpari

RM = rm -f

OBJS = $(TARGET).o
DYN = lib$(TARGET).so
ALL = $(TARGET) $(DYN)

all: $(ALL)

dynlib: $(DYN)

$(DYN): $(OBJS)
		$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) $(OBJS) $(EXTRAMODLDFLAGS)

$(TARGET): $(OBJS)
		$(LD) -o $@ $(LDFLAGS) $(EXTRACFLAGS) $< $(EXTRALIBS) $(RUNPTH) $(LIBS)

%.o: %.c
		$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
clean:
		-$(RM) *.o $(ALL)



# #eight

# # change this TARGET to compile your own programs
# TARGET = main
# SHELL  = /bin/sh

# DBGFLAGS   = -g -Wall
# CFLAGS     = -O3 -Wall -ffp-contract=off -fno-strict-aliasing
# EXTRACFLAGS=
# #CFLAGS    = $(DBGFLAGS)

# CC         = /usr/bin/gcc
# CPPFLAGS   = -I. -I/NOBACKUP/ericahl/pari-lib/include
# LD         = /usr/bin/gcc
# LDFLAGS    = -O3 -Wall -ffp-contract=off -fno-strict-aliasing    -Wl,--export-dynamic
# MODLD      = /usr/bin/gcc
# MODLDFLAGS = -shared  $(CFLAGS) $(DLCFLAGS) -Wl,-shared
# EXTRAMODLDFLAGS = -lc -lm -L/NOBACKUP/ericahl/pari-lib/lib -lpari
# EXTRALIBS  =

# RUNPTH     = -Wl,-rpath "/NOBACKUP/ericahl/pari-lib/lib"
# DLCFLAGS   = -fPIC
# LIBS       = -lm -L/NOBACKUP/ericahl/pari-lib/lib -lpari

# RM = rm -f

# OBJS = $(TARGET).o
# DYN = lib$(TARGET).so
# ALL = $(TARGET) $(DYN)

# all: $(ALL)

# dynlib: $(DYN)

# $(DYN): $(OBJS)
# 		$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) $(OBJS) $(EXTRAMODLDFLAGS)

# $(TARGET): $(OBJS)
# 		$(LD) -o $@ $(LDFLAGS) $(EXTRACFLAGS) $< $(EXTRALIBS) $(RUNPTH) $(LIBS)

# %.o: %.c
# 		$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
# clean:
# 		-$(RM) *.o $(ALL)

nine

# change this TARGET to compile your own programs
TARGET = main
SHELL  = /bin/sh

DBGFLAGS   = -g -Wall
CFLAGS     = -O3 -Wall -ffp-contract=off -fno-strict-aliasing
EXTRACFLAGS=
#CFLAGS    = $(DBGFLAGS)

CC         = /usr/bin/gcc
CPPFLAGS   = -I. -I/NOBACKUP/ericahl/pari/include
LD         = /usr/bin/gcc
LDFLAGS    = -O3 -Wall -ffp-contract=off -fno-strict-aliasing    -Wl,--export-dynamic
MODLD      = /usr/bin/gcc
MODLDFLAGS = -shared  $(CFLAGS) $(DLCFLAGS) -Wl,-shared
EXTRAMODLDFLAGS = -lc -lm -L/NOBACKUP/ericahl/pari/lib -lpari
EXTRALIBS  =

RUNPTH     = -Wl,-rpath "/NOBACKUP/ericahl/pari/lib"
DLCFLAGS   = -fPIC
LIBS       = -lm -L/NOBACKUP/ericahl/pari/lib -lpari

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

# # Mac version

# # change this TARGET to compile your own programs
# TARGET = main
# SHELL  = /bin/sh

# DBGFLAGS   = -g -Wall
# CFLAGS     = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer
# EXTRACFLAGS=
# #CFLAGS    = $(DBGFLAGS)

# # Various linkers use different flags to force static compilation. Choose
# # the one which is relevant for your installation.
# #
# # Solaris ld (global)
# #STATIC    = -dn

# # Solaris ld (toggle: no shared object accepted until -B dynamic is seen
# #STATIC    = -B static

# # gcc
# #STATIC    = -static

# CC         = /usr/bin/gcc
# CPPFLAGS   = -I. -I/usr/local/include -I/opt/homebrew/Cellar/pari/2.13.4/include
# LD         = /usr/bin/gcc
# LDFLAGS    = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer    -Wl,-search_paths_first 
# MODLD      = /usr/bin/gcc
# MODLDFLAGS = -bundle -undefined dynamic_lookup $(CFLAGS) $(DLCFLAGS)
# EXTRAMODLDFLAGS = 
# EXTRALIBS  =

# RUNPTH     = 
# DLCFLAGS   = -fPIC
# LIBS       = -L/usr/local/lib -lpari -L/opt/homebrew/Cellar/pari/2.13.4/lib -lpari

# RM = rm -f

# OBJS = $(TARGET).o
# # DYN = lib$(TARGET).dylib
# ALL = $(TARGET)-sta # $(TARGET)-dyn $(DYN)

# # dft: $(TARGET)-dyn

# all: $(ALL)

# sta: $(TARGET)-sta

# # dyn: $(TARGET)-dyn

# # dynlib: $(DYN)

# # $(DYN): $(OBJS)
# # 	$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) $(OBJS) $(EXTRAMODLDFLAGS)

# $(TARGET)-sta: $(OBJS)
# 	$(LD) -o $@ $(LDFLAGS) $(EXTRACFLAGS) $< $(EXTRALIBS) $(STATIC) $(LIBS)

# # $(TARGET)-dyn: $(OBJS)
# # 	$(LD) -o $@ $(LDFLAGS) $(EXTRACFLAGS) $< $(EXTRALIBS) $(RUNPTH) $(LIBS)

# %.o: %.c
# 	$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
# clean:
# 	-$(RM) *.o $(ALL) main-sta
# # massey-dyn libmassey.dylib


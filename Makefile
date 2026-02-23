CC      = gcc
CFLAGS  = -O2
BINDIR  = $(HOME)/bin

TARGETS = \
	$(BINDIR)/PCdrift \
	$(BINDIR)/polcap \
	$(BINDIR)/tracklos \
	$(BINDIR)/elipsDrift \
	$(BINDIR)/edgEllipse \
	$(BINDIR)/twoDdrift \
	$(BINDIR)/twoDout \
	$(BINDIR)/twoDvar \
	$(BINDIR)/oneDspark \
	$(BINDIR)/edgDist \
	$(BINDIR)/plotedge \
	$(BINDIR)/polcapfit \
	$(BINDIR)/polcaprot \
	$(BINDIR)/sprkplt \
	$(BINDIR)/sprkfig \
	$(BINDIR)/twoDelsp \
	$(BINDIR)/twoDgrid

.PHONY: all clean

all: $(BINDIR) $(TARGETS)

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/PCdrift:    PCdrift.c    PCdrift.h plotfunc.h sparkfunc.h usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lfftw3 -lcpgplot -lm

$(BINDIR)/polcap:     polcap.c     usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lm

$(BINDIR)/tracklos:   tracklos.c   usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lm

$(BINDIR)/elipsDrift: elipsDrift.c sparkfunc.h usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lcpgplot -lm

$(BINDIR)/edgEllipse: edgEllipse.c sparkfunc.h usefunc.h ellipse_fit.h
	$(CC) $(CFLAGS) $< -o $@ -lcpgplot -lm

$(BINDIR)/twoDdrift:  twoDdrift.c  usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lcpgplot -lm

$(BINDIR)/twoDout:    twoDout.c    usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lcpgplot -lm

$(BINDIR)/twoDvar:    twoDvar.c    usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lcpgplot -lm

$(BINDIR)/oneDspark:  oneDspark.c  usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lm

$(BINDIR)/edgDist:    edgDist.c    usefunc.h ellipse_fit.h
	$(CC) $(CFLAGS) $< -o $@ -lm

$(BINDIR)/plotedge:   plotedge.c   usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lm

$(BINDIR)/polcapfit:  polcapfit.c  usefunc.h ellipse_fit.h
	$(CC) $(CFLAGS) $< -o $@ -lm

$(BINDIR)/polcaprot:  polcaprot.c  usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lm

$(BINDIR)/sprkplt:    sprkplt.c    sparkfunc.h usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lcpgplot -lm

$(BINDIR)/sprkfig:    sprkfig.c    sparkfunc.h usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lcpgplot -lm

$(BINDIR)/twoDelsp:   twoDelsp.c   usefunc.h ellipse_fit.h
	$(CC) $(CFLAGS) $< -o $@ -lm

$(BINDIR)/twoDgrid:   twoDgrid.c   usefunc.h
	$(CC) $(CFLAGS) $< -o $@ -lm

clean:
	rm -f $(TARGETS)

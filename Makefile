PREFIX              = ~/.local/bin
HC                  = ghc
HFLAGS              = -W -O2

all: hoptics

hoptics:
	cd src && $(HC) $(HFLAGS) --make Hoptics.hs

clean:
	cd src && rm -f *.hi *.o Hoptics

install:
	cd src && cp Hoptics $(PREFIX)/hoptics

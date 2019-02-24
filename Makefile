TARGET=rhydrogen
TARGETMAIN=main
TARGETMAINSRC=$(TARGETMAIN).cxx
TARGETSRC=$(TARGET).cxx
TARGETHDR=$(TARGET).h
CMDOPT=cmdline
CMDOPTSRC=$(CMDOPT).c
CMDOPTHDR=$(CMDOPT).h
CC=g++
CFLAGS=`root-config --cflags`
RLINKFLAGS=`root-config --libs`
GSLCFLAGS=`gsl-config --cflags`
GSLLINKFLAGS=`gsl-config --libs`
OBJS=$(TARGET).o $(CMDOPT).o $(TARGETMAIN).o

all: $(OBJS)
	g++ $(RLINKFLAGS) $(OBJS) -o $(TARGET) 

$(OBJS): $(TARGETSRC) $(CMDOPTSRC) $(TARGETMAINSRC) 
	g++ $(CFLAGS) -c $(TARGETSRC)
	g++ $(CFLAGS) -c $(TARGETMAINSRC)
	g++ -I. -c $(CMDOPTSRC)


clean:
	rm *.o $(TARGET)

tar: 
#	rm Hydrogen.tar.gz
	tar -czf RHydrogen.tar.gz rhydrogen.cxx rhydrogen.h main.cxx cmdline.c cmdline.h README Makefile

# Use also with gsl
#
#all: $(OBJS)
#	g++ $(RLINKFLAGS) $(GSLLINKFLAGS) $(OBJS) -o $(TARGET) 
#
#$(OBJS): $(DICTFILE) $(TARGETSRC) $(CMDOPTSRC) $(TARGETMAINSRC) 
#	g++ $(CFLAGS) $(GSLCFLAGS) -c $(DICTFILE)
#	g++ $(CFLAGS) $(GSLCFLAGS) -c $(TARGETSRC)
#	g++ $(CFLAGS) $(GSLCFLAGS) -c $(TARGETMAINSRC)
#	g++ -I. -c $(CMDOPTSRC)

OBJS = main.o 2dMemory.o calc_eh.o calcH.o calcWLPcoef.o composeMatA.o Fourier.o Jr_test.o Jr.o Laguerre.o PECboundary.o toTimeDomain.o update_beta.o update_sum.o

main: $(OBJS)
	g++ -o $@ $(OBJS)

%.o: %.cpp const.h const_region.h const_sq.h const_time.h head.h PMLparameter.h r_theta.h
	g++ -c -Wall -O2 $<
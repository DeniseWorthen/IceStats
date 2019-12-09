#CDF=/apps/netcdf/4.6.1/intel/16.1.150/
CDF=/apps/netcdf/4.7.0/intel/18.0.5.274
#####################################################################
# compiler options
#####################################################################
#FOPT = -O
FOPT = -C
#FOPT = -convert big_endian
#FOPT = -p

F90 = ifort
#F90 = ifort -warn

#####################################################################
# 
#####################################################################

opt1 = -Duse_m6c5
#opt2 = -Duse_cfsv2
#opt2 = -Duse_sis2
#opt2 = -Duse_cpc

optall = $(opt1) $(opt2) $(opt3) $(opt4)

OBJS = param.o charstrings.o cdf.o variablelist.o regmask_regrid_north.o grdvar.o stats.o caldata.o runparams.o icestats.o tm_secs_from_bc.o write_cdf.o

makeit: $(OBJS) 
	$(F90) $(FOPT) -o makeit $(OBJS) -L$(CDF)/lib -lnetcdff -lnetcdf 

%.o: %.F90
	$(F90) $(FOPT) $(optall) -c -I$(CDF)/include $<
	cpp $(optall) -I$(CDF)/include $*.F90>$*.i

clean:
	/bin/rm -f makeit *.o *.i *.mod

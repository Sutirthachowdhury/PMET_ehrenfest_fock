FC=gfortran
#FLAGS = -ofast -cpp -ffree-line-length-none -fbacktrace -g -DDEBUG -D"ASSERT(x)=call assert(x, 'x', __LINE__, __FILE__)"
FLAGS = -lfftw3 -llapack -lblas -fbounds-check #-DDEBUG -D"ASSERT(x)=call assert(x, 'x', __LINE__, __FILE__)"     
SRC= global.f90 pimc.f90 gran.f90  parameter_input.f90 init_mapping_sampling.f90 mapping_pot.f90 histogram.f90 mapverlet.f90 plot_potential.f90 gethel.f90 getdhel.f90 run_traj.f90  init_nuclear_sampling.f90 

OBJS  =  ${SRC:.f90=.o}

all:  $(OBJS)
	$(FC) $(OBJS) mapping-rpmd-RT_old_code.f90 -o rpmd.exe $(FLAGS) 


%.o : %.f90
	$(FC) -c $(FLAGS) $< -o $@

clean:
	rm -rf *.o 


# this project sucks !!!
# no, try to work harder!!!!!!
#yes I agree.

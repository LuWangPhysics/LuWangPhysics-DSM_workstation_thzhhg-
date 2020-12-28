.SUFFIXES: .c .cpp .o .ex

# the compiler: gcc for C program, define as g++ for C++
  # compiler flags:
  #  -g    adds debugging information to the executable file
  #  -Wall turns on most, but not all, compiler warnings
  CC =g++ -fopenmp   -ggdb  -Wall  -std=c++14  -Wfatal-errors 
  #-std=c++11
  # -framework Python#for python
  #-fopenmp # for openmp
  #for compile without mpi, change mpic++ to g++, mpirun
  # g++ -g -Wall 


INCLUDE= -I/DSM_THzHHg_LuWang/include.hpp
INCLUDE+= -I /usr/local/include/eigen3/
INCLUDE+= -I /build/thirdparty/eigen/src/eigen/ 



LIBS = -L   /build/thirdparty/fftw/src/fftw/  -lfftw3 -lm  


 

.cpp.o:
	$(CC) $(INCLUDE) -c $< 

.c.o:
	gcc $(INCLUDE) -c $< 

.o.ex:
	@echo g++ ... -o $@ $< ... $(OBJS) ... $(LIBS) 
	@$(CC) -o $@ $< $(OBJS) $(LIBS)

######################################################################
# targets
######################################################################
OBJS  = 

objs: $(OBJS)

clean:
	rm -f *.o *.ex *~


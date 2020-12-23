.SUFFIXES: .c .cpp .o .ex

# the compiler: gcc for C program, define as g++ for C++
    # compiler flags:
  #  -g    adds debugging information to the executable file
  #  -Wall turns on most, but not all, compiler warnings
  CC =g++ -fopenmp -ggdb -framework Python -Wall  -std=c++14  -Wfatal-errors 
  #-std=c++11
  # -framework Python#for python
  #-fopenmp # for openmp
  #for compile without mpi, change mpic++ to g++, mpirun
  # g++ -g -Wall 


  # to include external libraries to your code


INCLUDE= -I /Users/luwang1/Documents/posdoc_related/NTU/third_order_thz_cpp/
INCLUDE+= -I /Users/luwang1/Documents/posdoc_related/NTU/third_order_thz_cpp/external_lib/eigen/eigen-master/
INCLUDE+= -I /Applications/Xcode.app/Contents/Developer/Library/Frameworks/Python3.framework/Versions/3.7/Headers/
INCLUDE+= -I /Applications/Xcode.app/Contents/Developer/Library/Frameworks/Python3.framework/Versions/3.7/lib/python3.7/importlib/
INCLUDE+=-I /usr/local/Cellar/python/3.7.7/Frameworks/Python.framework/Versions/3.7/include/python3.7m/
INCLUDE+= -I  /Users/luwang1/Library/Python/3.7/lib/python/site-packages/


LIBS = -L  /Users/luwang1/Documents/posdoc_related/NTU/third_order_thz_cpp/external_lib/fftw/fftw-3.3.8/ -lfftw3_threads -lfftw3 -lm  



 



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


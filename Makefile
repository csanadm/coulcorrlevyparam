all: calc_func.exe
clean:
	rm calc_func.o coulcorr_param.o calc_func.exe
calc_func.exe: calc_func.o coulcorr_param.o
	g++ $^ -o $@ 
calc_func.o: calc_func.cc coulcorr_param.h
	g++ -c calc_func.cc
coulcorr_param.o: coulcorr_param.cc coulcorr_param.h
	g++ -c coulcorr_param.cc

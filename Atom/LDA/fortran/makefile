main: main.f radial.f LDA.f
	ifort -c main.f
	ifort -c radial.f
	ifort -c LDA.f
	ifort -g -traceback -check LDA.o main.o radial.o -o main -llapack
	#ifort main.o fzero.o fun.o radial.o LDA.o -o main -llapack

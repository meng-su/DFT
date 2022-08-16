main: main.f fzero.f radial.f LDA.f ps.f besseljzero.f
	ifort -c main.f
	ifort -c radial.f
	ifort -c LDA.f
	ifort -c fzero.f
	ifort -c besseljzero.f
	ifort -c ps.f
	ifort -g -traceback -check LDA.o main.o fzero.o radial.o ps.o besseljzero.o -o main -llapack
	#ifort main.o fzero.o fun.o radial.o LDA.o -o main -llapack

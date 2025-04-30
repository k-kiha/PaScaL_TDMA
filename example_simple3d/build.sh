mpiifort -c main.f90 -I../include 
mpiifort *.o -I../include -L../lib -lpascal_tdma
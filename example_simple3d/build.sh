mpiifort -c -I../include main.f90 
mpiifort -I../include -L../lib -lpascal_tdma *.o
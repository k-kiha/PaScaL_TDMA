mpiifort -O3 -fPIC -fp-model=precise -xMIC-AVX512 -c main.f90 -I../include 
mpiifort -O3 -fPIC -fp-model=precise -xMIC-AVX512 *.o -I../include -L../lib -lpascal_tdma
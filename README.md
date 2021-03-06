# PaScaL_TDMA

Parallel and Scalable Library for TriDiagonal Matrix Algorithm

PaScal_TDMA provides an efficient and scalable computational procedure to solve many tridiagonal systems in multi-dimensional partial differential equations. The modified Thomas algorithm proposed by Laszlo et al.(2016) and the newly designed communication scheme have been used to reduce the communication overhead in solving many tridiagonal systems.

This library is for both single and many tridiagonal systems of equations. The main algorithm for a tridiagonal matrix consists of the following five steps: 

- (1) Transform the partitioned submatrices in the tridiagonal systems into modified submatrices:
        Each computing core transforms the partitioned submatrices in the tridiagonal systems of equations into the modified forms by applying modified Thomas algorithm.
- (2) Construct reduced tridiagonal systems from the modified submatrices:
        The reduced tridiagonal systems are constructed by collecting the first and last rows of the modified submatrices from each core using MPI_Ialltoallw.
- (3) Solve the reduced tridiagonal systems:
        The reduced tridiagonal systems constructed in Step 2 are solved by applying the Thomas algorithm.
- (4) Distribute the solutions of the reduced tridiagonal systems:
        The solutions of the reduced tridiagonal systems in Step 3 are distributed to each core using MPI_Ialltoallw.
        This communication is an exact inverse of the communication in Step 2.
- (5) Update the other unknowns in the modified tridiagonal systems:
        The remaining unknowns in the modified submatrices in Step 1 are solved in each computing core with the solutions obtained in Step 3 and Step 4.
    
Step 1 and Step 5 are similar to the method proposed by Laszlo et al.(2016) which uses parallel cyclic reduction (PCR) algorithm to build and solve the reduced tridiagonal systems. Instead of using the PCR, we develop an all-to-all communication scheme using the MPI_Ialltoall function after the modified Thomas algorithm is executed. The number of coefficients for the reduced tridiagonal systems are greatly reduced, so we can avoid the communication bandwidth problem, which is a main bottleneck for all-to-all communications. Our algorithm is also distinguished from the work of Mattor et al. (1995) which assembles the undetermined coefficients of the temporary solutions in a single processor using MPI_Gather, where load imbalances are serious.

# Authors
- Kiha Kim (k-kiha@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University
- Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
- Jung-Il Choi (jic@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University

# References
For more information, please the reference paper (in preparation) and [Multi-Physics Modeling and Computation Lab.](https://www.mpmc.yonsei.ac.kr/)
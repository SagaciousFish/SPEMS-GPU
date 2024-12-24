## SPEMS-GPU: GPU-Accelerated Parallel solver for EMS.

The current repo mainly contains the implementation of:

1. recursEMS:  a parallel EMS solver build upon a recursion tree model;
2. parEMS: EMS Solver with Introducing Concurrent Data Structure.
3. recursEMS-GPU: GPU-Accelerated recursion tree model for EMS.
4. parEMS-GPU: GPU-Accelerated concurrent data structure for EMS.

__Note that, this repo is still under construction, and it may not work well for all other real-world datasets.__

__Please feel free to contact us if you find there are any bugs.__

#### Compilation

There are some requirements for building the completed SPEMS-GPU: `cmake 3.30+` and `nvcc 12.6+`.

To compile SPEMS-GPU, you may use the following commands:

```
TODO: add the compilation commands
``` 

#### Running applications in SPEMS-GPU

The applications take the DNA datasets as input as well as some parallel execution arguments. For example:

```
TODO: add the usage
``` 

#### Input datasets

Input DNA datasets should be in form of plain text files, containing a list of DNA sequences.

## References

- J. Qiu and A. Ebnenasir, “Exploring scalable parallelization for edit distance-based motif search,” IEEE/ACM Trans.
  Comput. Biol. Bioinformatics, vol. 20, no. 2, pp. 1587–1593, Mar. 2023, doi:
  10.1109/TCBB.2022.3208867. (https://github.com/AutoPalSys/SPEMS)

- S. Pal, P. Xiao, and S. Rajasekaran, “Efficient sequential and parallel algorithms for finding edit distance based
  motifs,” BMC genomics, vol. 17, no. 4, pp. 315–326, 2016.
  (https://github.com/soumitrakp/ems2)

- P. Xiao, X. Cai, and S. Rajasekaran, “Ems3: An improved algorithm for finding edit-distance based motifs,” IEEE/ACM
  Transactions on Computational Biology and Bioinformatics, 2020
  (https://github.com/xiaopeng9371/EMS3)

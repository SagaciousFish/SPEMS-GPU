## SPEMS: Scalability-sensitive Parallel solver for EMS.

The current repo mainly contains the implementation of:
1. recursEMS:  a parallel EMS solver build upon a recursion tree model;
2. parEMS: EMS Solver with Introducing Concurrent Data Structure.

__Note that, this repo is still under construction and it may not work well for all other real-world datasets.__

__Please feel free to contact us if you find there are any bugs.__

#### Compilation

There are some requirements for building the completed SPEMS: `cmake3.2+` and `gcc4.9+`.
And you also need to install TBB (https://github.com/oneapi-src/oneTBB) before you complie the codes.

To compile recursEMS, just run `make` in the `./recursEMS` directory;

To compile parEMS, you may use the following commands:
```
$ cd ./parEMS
$ mkdir build && cd build && cmake ..
$ make
``` 

#### Running applications in SPEMS

The applications take the DNA datasets as input as well as some parallel execution arguments. For example:

```
$ ./ems -s recusEMS -l vaule-of-l -d value-of-d -t #cores path-to-input-graph
$ ./parEMS  -l vaule-of-l -d value-of-d -t #cores path-to-input-graph
``` 

#### Input datasets

Input DNA datasets should be in form of plain text files, containing a list of DNA sequences.

## References

- S. Pal, P. Xiao, and S. Rajasekaran, “Efficient sequential and parallel algorithms for finding edit distance based motifs,” BMC genomics, vol. 17, no. 4, pp. 315–326, 2016. 
(https://github.com/soumitrakp/ems2)

- P. Xiao, X. Cai, and S. Rajasekaran, “Ems3: An improved algorithm for finding edit-distance based motifs,” IEEE/ACM Transactions on Computational Biology and Bioinformatics, 2020 
(https://github.com/xiaopeng9371/EMS3)

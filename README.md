# Implementation of the GWSAT algorithm

C++ implementation of GWSAT, a stochastic local search SAT algorithm proposed by McAllester, Selman and Kautz.
For information see their paper: _Evidence for Invariants in Local Search, AAAI'97_

## How the algorithm works 

The algorithm can be used to test satisfiability of propositional formulas in CNF.
In the beginning the algorithm creates a random assignment of the variables.  
More TODO.

## Compilation and usage

A Makefile is included to compile the source code with g++. Apart from the C++11 standary library there are no further dependencies.

The program can be run on the command line, where input formulas must be passed as files in the DIMACS format. If the formula is satisfiable the solver will output the number of iterations it took to find a satisfiying assignment. Otherwise the solver will run until a maximal number of iterations is reached.

    $./GWSAT formula.cnf

Furthermore the maximum number of iterations (option -mi),  maximum number of restarts (-mr), noise parameter (-n) and initial RNG seed (-s) can be specified as command line options. For example:

    $./GWSAT formula.cnf -mi 100000007 -mr 42 -n 60 -s 137984113 

## Notes and copyright
This program was created as part of a university project. Many implementation choices and code simplifications were made to fit the specific requirements of the further project. Depending on your own requirements it may be advisable to modify the program.  

&copy; Daniel Kinzler

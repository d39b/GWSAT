# Implementation of the GWSAT algorithm

C++ implementation of GWSAT, a stochastic local search SAT algorithm proposed by McAllester, Selman and Kautz.
For more information see their paper: _Evidence for Invariants in Local Search, AAAI'97_

## How the algorithm works

The algorithm can be used to test satisfiability of propositional formulas in CNF. Starting with a random initial assignment, in each step the algorithm randomly selects an unsatisfied clause and flips the value of a variable within that clause. There are two ways to choose the variable to be flipped. Either the variable is selected randomly or the variable that results in the highest change in the number of satisfied clauses is chosen. Which selection method is used in a particular step is determined randomly according to a fixed noise probability. If no satisfying assignment could be found after a certain number of steps, the search can be restarted with a new random assignment.

## Compilation and usage

A makefile is included to compile the source code with g++. Apart from the C++11 standard library there are no further dependencies.

The program can be run on the command line, where input formulas must be passed as files in the DIMACS format. If the formula is satisfiable the solver will output the number of iterations it took to find a satisfying assignment. Otherwise the solver will run until a maximum number of iterations is reached.

> ./GWSAT formula.cnf

Furthermore the maximum number of iterations (option -mi),  maximum number of restarts (-mr), noise parameter (-n, integer between 0 and 100) and initial RNG seed (-s) can be specified as command line options. For example:

> ./GWSAT formula.cnf -mi 100000007 -mr 42 -n 60 -s 137984113

## License  

[MIT](LICENSE)

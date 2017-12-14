#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <random>

/*
Author: Daniel Kinzler
*/


/*
Implementation of GWSAT, a stochastic local search SAT algorithm proposed by McAllester/Selman/Kautz.
For more information see their paper:
Evidence for Invariants in Local Search
appearing in the Proceedings of the Fourteenth National Conference on Artificial Intelligence (AAAI-97)
*/


//type used for all necessary values such as number of iterations, rng seeds, variables etc.
//for larger formulas (e.g. >1000 variables) 32bit integers may be too small to hold some values
typedef long long int long_t;

/*
Wrapper class for generating random numbers that internally uses a Mersenne Twister 19937 generator.

The constructor can take a seed parameter to initialize the RNG.
*/
class Random {
    private: 
        std::mt19937* random_generator;
        
    public:
        Random() {
            random_generator = new std::mt19937();
        }

	Random(long_t seed) {
		random_generator = new std::mt19937(seed);
	}
        
        ~Random() {
            delete random_generator;
        }
	
        //returns a uniformly distributed random integer between a and b (inclusive)
        long_t uniform_int(long_t a, long_t b) {
            std::uniform_int_distribution<long_t> distribution(a,b);
            return distribution(*random_generator);
        }
};


/*
Solver superclass a specific SAT Solver should extend.
Contains fields and methods for the representation of a propositional formula. 

Variables and clauses of a formula are numbered starting with 0.
A formula is represented as an array of clauses, where a clause is in turn represented as an array of literals.
A positive literal for variable i is represented by the number 2*i + 1 and a negative literal by the number 2*i.
This allows for the efficient manipulation of literals using bitwise operations, e.g. to extract the variable from a literal l a simple bitshift (l >> 1) is sufficient.

A SatSolver instance should be initialized by calling the constructor with the number of variables and clauses of the input formula as argmuents.  
After that the actual clauses and literals can be added using the methods add_clause() and add_lit_to_clause().
To solve another formula, a new SatSolver instance should be created.
*/
class SatSolver {
	
	protected:			
		//array of clauses; each clause is an array of literal numbers
		long_t ** clauses;
		
		//array containing the length (number of literals) of each clause
		long_t * clause_size;
		
		//number of variables of the formula
		long_t num_var;
		
		//number of clauses of the formula
		long_t num_clause;
		
		// array of 0/1 values representing the truth value of each variable
		long_t * assignment;	
	
	public:
		SatSolver (long_t num_var, long_t num_clause) : num_var(num_var), num_clause(num_clause) {clauses = new long_t*[num_clause]; clause_size = new long_t[num_clause]; assignment = new long_t[num_var];} 
       
               
		virtual ~SatSolver () {}	
	 
		//allocates memory for a new clause, the clause can later on be referenced by the given index
	        //arguments are not checked for validity since these methods are only used by the DIMACS parser method 
		void add_clause(long_t length, long_t index) {
		    clauses[index] = new long_t[length];
		    clause_size[index] = length;
	        }
               
		//adds a given literal to a clause
		void add_lit_to_clause(long_t clause_index, long_t lit_index, long_t lit) {
                        clauses[clause_index][lit_index] = lit;
                }
		
		void print_assignment() {
			for(int i = 0; i < num_var; i++) {
				std::cout << "var " << i << ": " << assignment[i] << "\n";
			}
		}
		
		//subclasses should override this method 
                virtual bool solve() = 0;
        	
};

/*
Subclass implementing the GWSAT solver.

After calling the constructor, the solver should be initialized with the input formula (see description of superclass SatSolver).
When all clauses have been added, the init() function must be called.
Then solve() or solve(max_restarts) can be called to try to find a solution (these methods leave the solver in a consitent state, i.e. solve() functions can be called repeatedly to (maybe) get different satisfying assignments )
*/
class GWSAT : public SatSolver {
	private:
		//array containing the indices of all unsatisified clauses under the current assignment
		long_t * unsatisfied_clauses;

		//idx[j] is the index at which clause j is stored in the array unsatisfied_clauses 
		//i.e if j is the index of an unsatisifed clause then unsatisfied_clauses[idx[j]] = j;
		long_t * idx
;
		//array containing the number of satisfied literals in each clause
		long_t * true_lit_count;
	
		//an array for each literal containing the indices of the clauses that contain this literal 	
		long_t ** clauses_of_lit;
		
		//array containing the lengths of the arrays of clauses_of_lit
		long_t * col_lengths;
	
		//maximal number of iterations before a restart with a new random assignment	
		long_t max_iterations;
		
		//number between 0 and 100 corresponding to the probability that the variable to be flipped is picked at random
		long_t rand_prob;

		//number of possible literals, should be 2 times the number of variables
		long_t num_lit;

		//array storing the variables with minimal make-break count 
		long_t * min_indices;

		//counter for the total number of iterations of this solver instance	
                long_t num_iterations = 0;
           
		//random number generator
                Random * rand;
	
	
	public:
		GWSAT (long_t num_var, long_t num_clause, long_t max_iterations, long_t rand_prob) : SatSolver(num_var,num_clause), max_iterations(max_iterations), rand_prob(rand_prob) {
			num_lit = num_var << 1;
			min_indices = new long_t[num_lit];
			rand = new Random();
		}
		
		GWSAT (long_t num_var, long_t num_clause, long_t max_iterations, long_t rand_prob, long_t rng_seed) : SatSolver(num_var,num_clause), max_iterations(max_iterations), rand_prob(rand_prob) {
			num_lit = num_var << 1;
			min_indices = new long_t[num_lit];
			rand = new Random(rng_seed);
		}		

		~GWSAT() {
			delete[] assignment;
			delete[] unsatisfied_clauses;
			delete[] idx;
			delete[] true_lit_count;
			delete[] clauses_of_lit;
			delete[] col_lengths;
			delete[] min_indices;
			for(long_t i = 0; i < num_clause; ++i) {
				delete clauses[i];
			}
			
			delete[] clauses;
			delete[] clause_size;
                        delete rand;
		}
	
                void generate_random_assignment() {
	            for(long_t i = 0; i < num_var; i++) {
                        assignment[i] = rand -> uniform_int(0,1);
                    }    
	        }
    
		//returns the number of unsatisfied clauses for the current assignment
		//initializes the arrays true_lit_count, unsatisfied_clauses and idx with the correct values
                long_t find_unsatisfied_clauses() {
	        	long_t num_unsatisfied = 0;

		        for(long_t i = 0; i < num_clause; ++i) {
			    true_lit_count[i] = 0;
			    bool clause_satisfied = false;
			
			    for(long_t j = 0; j < clause_size[i]; ++j) {
			        long_t lit = clauses[i][j];
				long_t var = lit >> 1;
				long_t sign = lit & 1;

				if(sign == assignment[var]) {
					clause_satisfied = true;
					true_lit_count[i]++;
				}

			    }			
			    if(!clause_satisfied) {
				num_unsatisfied++;
				unsatisfied_clauses[num_unsatisfied] = i;
				idx[i] = num_unsatisfied;
			    }

		        }
		
		        return num_unsatisfied;

	        }


	//must be called after all clauses have been added, to initialize all necessary variables  
        void init() {
		unsatisfied_clauses = new long_t[num_clause+1];
		idx = new long_t[num_clause];
		
		true_lit_count = new long_t[num_clause];

		col_lengths = new long_t[num_lit]();
		for(long_t i = 0; i < num_clause; ++i) {
			for(long_t j = 0; j < clause_size[i]; ++j) {
				long_t lit = clauses[i][j];
				col_lengths[lit]++;
			}
		}
		clauses_of_lit = new long_t*[num_lit];
		for(long_t i = 0; i < num_lit; ++i) {
			clauses_of_lit[i] = new long_t[col_lengths[i]];
			col_lengths[i] = 0;
		}		
		for(long_t i = 0; i < num_clause; ++i) {
			for(long_t j = 0; j < clause_size[i]; ++j) {
				long_t lit = clauses[i][j];
				clauses_of_lit[lit][col_lengths[lit]] = i;
				col_lengths[lit]++;
			}
		}
	}	
	
        long_t get_num_iterations() {
            return num_iterations;
        }
	
	//tries to find a satisfying assignment using at most max_iterations flips
	//returns true iff a satisfying assignment has been found
	bool solve() {
                generate_random_assignment();
		
		//num_u will always contain the number of unsatisified clauses for the current assignment
		long_t num_u = find_unsatisfied_clauses();
		
                for(long_t i = 0; i < max_iterations; ++i) {
                        num_iterations++;
			
			//if there are no unsatisfied clauses, we have found a satisfying assignment	
                        if(num_u == 0) {
				return true;
			}
			
			//choose an unsatisfied clause at random
			long_t pr_index = rand -> uniform_int(1,num_u);
			long_t* c = clauses[unsatisfied_clauses[pr_index]];
			
			
			long_t min = num_clause;
			long_t min_var = -1;
		        long_t min_count = 0;

                        long_t p = rand -> uniform_int(1,100);
                        
			//choose the variable to be flipped in the clause at random
			if(p <= rand_prob) {
                            long_t lit = c[rand -> uniform_int(0,clause_size[pr_index]-1)];
                            min_var = lit >> 1;                             
                        } else {
			//choose the variable to be flipped as the variable with the lowest break-make count
			        
				//compute the break-make count for each literal in the clause
				for(long_t j = 0; j < clause_size[pr_index]; ++j) {
			            long_t lit = c[j];
				    long_t neg_lit = lit ^ 1;		
				    long_t var = lit >> 1;		
				    long_t score = 0;
				    
				    //if the number of true literals in a clause, that contains lit, is 0, the clause would become satisfied if we flip lit, thus adding 1 to the make count of lit  
				    for(long_t k = 0; k < col_lengths[lit]; ++k) {
				        long_t clause_index = clauses_of_lit[lit][k];
					if(true_lit_count[clause_index] == 0) {
							score--;
					}			
				    }
				    //if (not lit) is the only true literal in a clause, that clause would become unsatisfied if we flip the truth value of lit, thus adding 1 to the break count of lit
				    for(long_t k = 0; k < col_lengths[neg_lit]; ++k) {
					long_t clause_index = clauses_of_lit[neg_lit][k];
					if(true_lit_count[clause_index] == 1) {
						score++;
					
						if(score > min) {
							break;
						}
							
					}
			            }
	            		    if(score < min) {
					    min = score;
                                            min_indices[0] = var;
                                            min_count = 1;
				    } else if(score == min) {
                                        min_indices[min_count] = var;
                                        min_count++;
                                    }								
			        }
				
				//select a variable with minimal break-make count
                                if(min_count == 1) {
                                    min_var = min_indices[0];
                                } else {
				    //if there is more than one literal in the clause that has minimal break-make count, we choose one at random
                                    pr_index = rand -> uniform_int(0,min_count-1);
                                    min_var = min_indices[pr_index];
                                }
			}

			//flip the truth value of the chosen variable	
			assignment[min_var] = assignment[min_var] ^ 1;
			long_t lit = (min_var << 1) + assignment[min_var];
			long_t neg_lit = lit ^ 1;
		
			//update true_lit_count and unsatisifed_clauses for unsatisfied clauses that become satisfied
			for(long_t k = 0; k < col_lengths[lit]; ++k) {
				long_t clause_index = clauses_of_lit[lit][k];
				true_lit_count[clause_index]++;
				if(true_lit_count[clause_index] == 1) {
					long_t tmp_idx = unsatisfied_clauses[num_u];
					unsatisfied_clauses[idx[clause_index]] = tmp_idx;
					idx[tmp_idx] = idx[clause_index];
					idx[clause_index] = 0;
					num_u--;
				}			
			}
			//update true_lit_count and unsatisfied_clauses for satisfied clauses that become unsatisfied
			for(long_t k = 0; k < col_lengths[neg_lit]; ++k) {
				long_t clause_index = clauses_of_lit[neg_lit][k];
				true_lit_count[clause_index]--;
				if(true_lit_count[clause_index] == 0) {
					num_u++;
					unsatisfied_clauses[num_u] = clause_index;
					idx[clause_index] = num_u;
				}
			}
			

		}
	
		//check if formula is satisfied after last iteration	
                if(num_u == 0) {
                    return true;
                }
                return false;		

	}

	//tries to find a satisfying assignment with at most max_restarts restarts
        bool solve(long_t max_restarts) {
            long_t sat = false;
            for(long_t i = 0; i < max_restarts; i++) {
                sat = solve();
                if(sat) {
                    break;
                }
            }
            return sat;
        }
       
	//method that tests if the values in the various array are correctly maintained 
        long_t static_correctness_check() {
            long_t num_unsatisfied = 0;
            for(long_t i = 0; i < num_clause; i++) {
                bool sat = false;
                long_t tlc = 0;
                long_t* c = clauses[i];
                for(long_t j = 0; j < clause_size[i]; j++) {
                    long_t lit = c[j];
                    long_t var = lit >> 1;
                    long_t sign = lit & 1;
                    if(assignment[var] == sign) {
                        sat = true;
                        tlc++;
                    }
                }
                if(tlc != true_lit_count[i]) {
                    std::cout << "error in true lit count";
                }
                if(!sat){
                    if(unsatisfied_clauses[idx[i]] != i) {
                        std::cout << "error in unsatisfied clauses";
                    }
                    num_unsatisfied++;
                }
            }             
            return num_unsatisfied;
        }
        	
		
};

//parse a DIMACS file, create a corresponding solver instance with the given parameters, then try to find a solution
//returns the number of iterations it took to find a satisfying assignment or -1 if no such assignment could be found
long_t run_on_file(std::string filename, long_t max_iterations, long_t max_restarts, long_t rand_prob, long_t rng_seed) {
	std::ifstream infile(filename);
	
	if(infile.fail()) {
		std::cout << "Error, could not open file: " << filename << "\n";
		return -1;
	}
	
	std::string line;
	
	
	bool header_found = false;
	long_t num_var;
	long_t num_clause = 0;
	long_t * tmp = NULL;

        GWSAT * solver = NULL;		
	long_t clause_count = 0;	
	
	long_t * literals = NULL;

	//parse the DIMACS file
	while(std::getline(infile,line)) {
		if(line[0] != 'c' && line[0] != '0' && line[0] != '%' && (!header_found || clause_count < num_clause)) {
			std::istringstream iss(line);
			if(line[0] == 'p' && !header_found) {
				std::string t;
				iss >> t;
				iss >> t;
				iss >> num_var;
				iss >> num_clause;
				tmp = new long_t[2*num_var];			
				if(rng_seed == -1) {
					solver = new GWSAT(num_var,num_clause,max_iterations,rand_prob);
				} else {
					solver = new GWSAT(num_var,num_clause,max_iterations,rand_prob, rng_seed);
				}
				literals = new long_t[2*num_var];
				header_found = true;
			} else {
				if(!header_found) {
					std::cout << "Error while parsing: header not found\n";
					return -1;	
				} else {
					
					long_t index = 0;
					while(index < num_var && (iss >> tmp[index])) {
						if(tmp[index] == 0) {
							break;
						}
						index++;
					}		
					solver -> add_clause(index,clause_count);	
					
					for(long_t i = 0; i < index; ++i) {
						long_t lit = tmp[i];
						if(lit > 0) {
							lit = (lit-1)*2 + 1;
						} else {
							lit = ((lit*(-1))-1)*2;
						}
						solver -> add_lit_to_clause(clause_count, i, lit);
					}	
					clause_count++;
							
				}	



			}
		}
	}
	solver -> init();

	bool solution_found = solver -> solve(max_restarts);
        long_t result = -2;
        if(solution_found) {
            result = solver -> get_num_iterations();
        }	

	if(tmp != NULL) {
		delete[] tmp;
	}
	if(literals != NULL) {
		delete[] literals;
	}
	if(solver != NULL) {
		delete solver;
        }
        return result;
}

void print_help_message() {
    std::cout << "Usage:\n\n";
    std::cout << "gwsat filename [OPTIONS]\n";
    std::cout << "Additional options:\n\n";
    std::cout << "-mi NUM       set NUM as the maximum number of iterations\n";
    std::cout << "-mr NUM       set NUM as the maximum number of restarts\n";
    std::cout << "-n NUM        set NUM as the noise parameter, NUM should be a number between 0 and 100\n";
    std::cout << "-s NUM	set NUM as the initial seed of the random number generator\n\n";
    std::cout << "Default number of iterations is 2^63-1\n";
    std::cout << "Default number of restarts is 1\n";
    std::cout << "Default noise parameter is 60\n";
}


//see print_help_message() function for usage
//checks the command line arguments for validity and calls the appropriate methods to solve the given SAT instances
//output is a single line containing the number of iterations it took to find a satisfying assignment
int main(int argc, char* argv[]) {
	
	long_t rand_prob = 60;
        long_t max_iterations = 1;
	max_iterations = ~(max_iterations << 63);
	long_t max_restarts = 1;
	long_t rng_seed = -1;	
        
        if(argc < 2) {
	    std::cout << "Error: missing arguments\n\n";
            print_help_message();
        } else {
            std::string filename(argv[1]);
            long_t arg_index = 2;
            while(arg_index+1 < argc) {
            	std::string option(argv[arg_index]);
                std::string option_arg(argv[arg_index+1]);
                    
                long_t num;
		try {
			num = std::stoll(option_arg);                    
    		} catch(std::exception& e) {
                    	std::cout << "Invalid argument " << option_arg << ", positive integer expected\n\n";
			print_help_message();
			return 0;	
		}

                if(option == "-mi") {
                	max_iterations = num;            
                } else if(option == "-mr") {
                        max_restarts = num;
                } else if(option == "-n") {
                        rand_prob = num;
                } else if(option == "-s") {
			rng_seed = num;	
		} else {
                        print_help_message();
                        return 0;
                }

                arg_index += 2;
            }                     
            
	    /*  	       	
	    std::cout << "running on input: " << filename << "\n";
            std::cout << "max iterations: " << max_iterations << "\n";
            std::cout << "max restarts: " << max_restarts << "\n";
            std::cout << "noise: " << rand_prob << "\n";                
	    std::cout << "seed: " << rng_seed << "\n";	    	             
	    */
		
            bool solution_found = false;
                    
	    long_t num_iterations = run_on_file(filename,max_iterations,max_restarts,rand_prob, rng_seed);
	
            if(num_iterations == -1) {
            	return 0;	
	    } else if(num_iterations > 0) {
           	solution_found = true; 	 	       
            } else {
                num_iterations = max_iterations*max_restarts;
	    }

	    if(solution_found) {
		std::cout << num_iterations << "\n";
	    } else {
		std::cout << -1 << "\n";
	    }
   	    
	    std::cout << "seed: " << rng_seed << "\n"; 
			 
        }     
	
	return 0;
}

#ifndef __SOLUTION_H__
#define __SOLUTION_H__

class CSolution {
private:
public:
    CDiscretization disc;  // spatial discretization class
    CPreconditioner prec;  // precondtioner class 
    CSolver solv;          // linear and nonlinear solvers
    //ofstream filestr;     // storing residual norms
    
    // constructor 
    CSolution(string filein, string fileout, Int mpiprocs, Int mpirank, Int ompthreads, Int omprank, Int backend)   
       : disc(filein, fileout, mpiprocs, mpirank, ompthreads, omprank, backend),
         prec(disc, backend),
         solv(disc, backend) { };        
    
    // destructor        
    ~CSolution(){  }; 

    // solve steady-state problems
    void SteadyProblem(ofstream &out, Int backend);    
            
    // advance the solution to the next time step using DIRK and BDF schemes
    void TimeStepping(ofstream &out, dstype time, Int istep, Int backend);
        
    // solve time-dependent problems
    void UnsteadyProblem(ofstream &out, Int backend);    
        
    // solve time-dependent problems
    void DIRK(ofstream &out, Int backend);    
    
    // precompute some quantities
    void InitSolution(Int backend);    
    
    // solve both steady-state and time-dependent problems
    void SolveProblem(ofstream &out, Int backend);    
    
    // save solutions in binary files
    void SaveSolutions(Int backend);    

    // save solutions in binary files
    void SaveSolutionsOnBoundary(Int backend);    
    
    // save solutions in binary files
    void SaveNodesOnBoundary(Int backend);    
    
    // read solutions from binary files
    void ReadSolutions(Int backend);        
    
    // save output in binary files
    void SaveOutputDG(Int backend);        
    
    // save output in binary files
    void SaveOutputCG(Int backend);        
};

#endif        



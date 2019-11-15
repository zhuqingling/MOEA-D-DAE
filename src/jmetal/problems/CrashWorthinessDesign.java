//  CrashWorthinessDesign.java
//
//  Author:
//       Yi Xiang <antonio@lcc.uma.es>


/**
 * This is a real world problem presented in the following paper:
 * Multiobjective optimization for crash safety design of vehicles using stepwise regression model
   Xingtao Liao & Qing Li & Xujing Yang &
    Weigang Zhang & Wei Li
 */

package jmetal.problems;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

/**
 * Class representing problem Water
 */
public class CrashWorthinessDesign extends Problem {                     

 /**
  * Constructor.
  * Creates a default instance of the CrashWorthinessDesign problem.
  * @param solutionType The solution type must "Real" or "BinaryReal".
  */
  public CrashWorthinessDesign(String solutionType) {
    numberOfVariables_   = 5 ;
    numberOfObjectives_  = 3 ;
    numberOfConstraints_ = 0 ;
    problemName_         = "CrashWorthinessDesign";
	        
    upperLimit_ = new double[numberOfVariables_];
    lowerLimit_ = new double[numberOfVariables_];

    for (int var = 0; var < numberOfVariables_; var++){
      lowerLimit_[var] = 1.0;
      upperLimit_[var] = 3.0;
    } // for
	        
    if (solutionType.compareTo("BinaryReal") == 0)
      solutionType_ = new BinaryRealSolutionType(this) ;
    else if (solutionType.compareTo("Real") == 0)
    	solutionType_ = new RealSolutionType(this) ;
    else {
    	System.out.println("Error: solution type " + solutionType + " invalid") ;
    	System.exit(-1) ;
    }  
 } // CrashWorthinessDesign
	
  /**
   * Evaluates a solution
   * @param solution The solution to evaluate
   * @throws JMException 
   */
  public void evaluate(Solution solution) throws JMException {         
    double [] x = new double[5] ; // 5 decision variables
    double [] f = new double[3] ; // 3 functions
    
    x[0] = solution.getDecisionVariables()[0].getValue();
    x[1] = solution.getDecisionVariables()[1].getValue();
    x[2] = solution.getDecisionVariables()[2].getValue();
    x[3] = solution.getDecisionVariables()[3].getValue();
    x[4] = solution.getDecisionVariables()[4].getValue();
    
    // First function
    f[0] = 1640.2823 + 2.3573285 * x[0] + 2.3220035 * x[1] + 4.5688768* x[2] 
    		+ 7.7213633 * x[3] + 4.4559504 * x[4];
    // Second function
    f[1] = 6.5856 + 1.15 * x[0] - 1.0427 * x[1] + 0.9738 * x[2] + 0.8364 * x[3]
    		- 0.3695 *  x[0] * x[3] + 0.0861*  x[0] * x[4] + 0.3628*  x[1] * x[3]
    		- 0.1106 *  x[0] * x[0] - 0.3437 *  x[2] * x[2] 
    		+ 0.1764 *  x[3] * x[3];
    // Third function
    f[2] = - 0.0551 + 0.0181 * x[0] + 0.1024 * x[1]
    		+ 0.0421 * x[2] - 0.0073 * x[0]* x[1]
    		+ 0.024 * x[1]* x[2] - 0.0118* x[1]* x[3]
    		-0.0204 * x[2]* x[3] - 0.008 * x[2]* x[4]
    		- 0.0241 * x[1]* x[1] + 0.0109 * x[3]* x[3] ;
  
             
    solution.setObjective(0,f[0]);    
    solution.setObjective(1,f[1]);
    solution.setObjective(2,f[2]);

  } // evaluate

  /** 
   * Evaluates the constraint overhead of a solution 
   * @param solution The solution
   * @throws JMException 
   */  
  public void evaluateConstraints(Solution solution) throws JMException {
  }
} // Water

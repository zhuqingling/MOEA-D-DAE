//  CEC2009_UF3.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.problems.cec2009Competition;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

/**
 * Class representing problem CEC2009_UF3
 */
public class CF5 extends Problem {
    
 /** 
  * Constructor.
  * Creates a default instance of problem CEC2009_UF3 (30 decision variables)
  * @param solutionType The solution type must "Real" or "BinaryReal".
  */
  public CF5(String solutionType) throws ClassNotFoundException {
    this(solutionType, 10); // 30 variables by default
  } // CEC2009_CF4
  
 /**
  * Creates a new instance of problem CEC2009_UF3.
  * @param numberOfVariables Number of variables.
  * @param solutionType The solution type must "Real" or "BinaryReal".
  */
  public CF5(String solutionType, Integer numberOfVariables) {
    numberOfVariables_  = numberOfVariables;
    numberOfObjectives_ =  2;
    numberOfConstraints_=  1;
    problemName_        = "CEC2009_CF5";

    upperLimit_ = new double[numberOfVariables_];
    lowerLimit_ = new double[numberOfVariables_];

    // Establishes upper and lower limits for the variables
    lowerLimit_[0] = 0.0 ;
    upperLimit_[0] = 1.0 ;
    for (int var = 1; var < numberOfVariables_; var++){
    	lowerLimit_[var] = -2.0;
    	upperLimit_[var] = 2.0;
    } //for

    if (solutionType.compareTo("BinaryReal") == 0)
    	solutionType_ = new BinaryRealSolutionType(this) ;
    else if (solutionType.compareTo("Real") == 0)
    	solutionType_ = new RealSolutionType(this) ;
    else {
    	System.out.println("Error: solution type " + solutionType + " invalid") ;
    	System.exit(-1) ;
    }  
  } // CEC2009_UF3
  
  private static double MYSIGN(double x) {
	  return(x>0.0?1.0:-1.0);
  }
  
  /** 
   * Evaluates a solution.
   * @param solution The solution to evaluate.
   * @throws JMException 
   */
  public void evaluate(Solution solution) throws JMException {
    Variable[] decisionVariables  = solution.getDecisionVariables();
    
    double [] x = new double[numberOfVariables_];
    double [] f = new double[numberOfObjectives_];
    for (int i = 0; i < numberOfVariables_; i++)
    	x[i] = decisionVariables[i].getValue() ;

    double sum1, sum2, yj;
	sum1   = sum2   = 0.0;
    
    for (int j = 2 ; j <= numberOfVariables_; j++) {
    	if (j % 2 == 1) {
			yj    = x[j-1] - 0.8*x[0]*Math.cos(6.0*Math.PI*x[0] + j*Math.PI/numberOfVariables_);
			sum1 += 2.0*yj*yj - Math.cos(4.0*Math.PI*yj) + 1.0;
		} else {
			yj = x[j-1] - 0.8*x[0]*Math.sin(6.0*Math.PI*x[0] + j*Math.PI/numberOfVariables_);
			if (j==2) {
				sum2 += yj < 1.5-0.75*Math.sqrt(2.0) ? Math.abs(yj) : (0.125+(yj-1)*(yj-1));
			} else {
				sum2 += 2.0*yj*yj - Math.cos(4.0*Math.PI*yj) + 1.0;
			}
		}
    }
    
    f[0] = x[0]		  + sum1;
	f[1] = 1.0 - x[0] + sum2;
    // constraint evaluation
    double constraint;
    constraint = x[1] - 0.8*x[0]*Math.sin(6.0*x[0]*Math.PI+2.0*Math.PI/numberOfVariables_) - 0.5*x[0] + 0.25;
    double total = 0.0;
    int number = 0;
    if (constraint < 0.0) {
        total += -constraint;
        number++;
        solution.setConstraint(0, -constraint);
    }else {
    	solution.setConstraint(0,0.0);
    }
    solution.setOverallConstraintViolation(total);
//  solution.setOverallConstraintViolation(total/(double)numberOfConstraints_);
    solution.setNumberOfViolatedConstraint(number);

    for (int i = 0; i < numberOfObjectives_; i++)
    	solution.setObjective(i, f[i]);
  } // evaluate
} // CEC2009_UF3

//  CEC2009_UF2.java
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
 * Class representing problem CEC2009_UF2
 */
public class CF2 extends Problem {
	
 /**
  * Constructor.
  * Creates a default instance of problem CEC2009_UF2 (30 decision variables)
  * @param solutionType The solution type must "Real" or "BinaryReal".
  */
  public CF2(String solutionType) throws ClassNotFoundException {
    this(solutionType, 10); // 30 variables by default
  } // CEC2009_UF2

 /**
  * Creates a new instance of problem CEC2009_UF2.
  * @param numberOfVariables Number of variables.
  * @param solutionType The solution type must "Real" or "BinaryReal".
  */
  public CF2(String solutionType, Integer numberOfVariables) {
    numberOfVariables_  = numberOfVariables;
    numberOfObjectives_ =  2;
    numberOfConstraints_=  1;
    problemName_        = "CEC2009_CF2";

    upperLimit_ = new double[numberOfVariables_];
    lowerLimit_ = new double[numberOfVariables_];

    lowerLimit_[0] = 0.0 ;
    upperLimit_[0] = 1.0 ;
    for (int var = 1; var < numberOfVariables_; var++){
    	lowerLimit_[var] = -1.0;
    	upperLimit_[var] = 1.0;
    } //for

    if (solutionType.compareTo("BinaryReal") == 0)
    	solutionType_ = new BinaryRealSolutionType(this) ;
    else if (solutionType.compareTo("Real") == 0)
    	solutionType_ = new RealSolutionType(this) ;
    else {
    	System.out.println("Error: solution type " + solutionType + " invalid") ;
    	System.exit(-1) ;
    }
  } // CEC2009_UF2
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

  	int count1, count2;
	double sum1, sum2, yj, N, a, t;
	N = 2.0; a = 1.0;
	sum1   = sum2   = 0.0;
	count1 = count2 = 0;

    for (int j = 2 ; j <= numberOfVariables_; j++) {
    	yj = x[j-1] - Math.sin(6.0*Math.PI*x[0] + j*Math.PI/numberOfVariables_);
		if (j % 2 == 1) {
			yj = x[j-1] - Math.sin(6.0*Math.PI*x[0] + j*Math.PI/numberOfVariables_);
			sum1  += yj*yj;
			count1++;
		} else {
			yj = x[j-1] - Math.cos(6.0*Math.PI*x[0] + j*Math.PI/numberOfVariables_);
			sum2  += yj*yj;
			count2++;
		}
    }
    f[0] = x[0]		        + 2.0*sum1 / (double)count1;
	f[1] = 1.0 - Math.sqrt(x[0]) + 2.0*sum2 / (double)count2;
	t	 = f[1] + Math.sqrt(f[0]) - a*Math.sin(N*Math.PI*(Math.sqrt(f[0])-f[1]+1.0)) - 1.0;
    // constraint evaluation
    double constraint;
    constraint = MYSIGN(t)*Math.abs(t)/(1+Math.exp(4.0*Math.abs(t)));
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
} // CEC2009_UF2

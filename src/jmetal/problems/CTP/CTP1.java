//  CTP2.java
//
//  Author:
//  wenji li Email: wenji_li@126.com

package jmetal.problems.CTP;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;


/**
 * Class representing problem CTP1
 */
public class CTP1 extends Problem {
    /**
     * Creates an  CTP1 problem (10 variables and 2 objectives)
     *
     * @param solutionType The solution type must "Real" or "BinaryReal".
     */
    public CTP1(String solutionType) throws ClassNotFoundException {
        this(solutionType, 4, 2);
    } // CTP1

    /**
     * Creates an CTP1 problem instance
     *
     * @param numberOfVariables  Number of variables
     * @param numberOfObjectives Number of objective functions
     * @param solutionType       The solution type must "Real" or "BinaryReal".
     */
    public CTP1(String solutionType,
                Integer numberOfVariables,
                Integer numberOfObjectives) throws ClassNotFoundException {
        numberOfVariables_ = numberOfVariables.intValue();
        numberOfObjectives_ = numberOfObjectives.intValue();
        numberOfConstraints_ = 2;
        problemName_ = "CTP1";

        lowerLimit_ = new double[numberOfVariables_];
        upperLimit_ = new double[numberOfVariables_];
        lowerLimit_[0] = 0.0;
        upperLimit_[0] = 1.0;

        for (int var = 1; var < numberOfVariables; var++) {
            lowerLimit_[var] = -5.0;
            upperLimit_[var] = 5.0;
        } //for
        if (solutionType.compareTo("BinaryReal") == 0)
            solutionType_ = new BinaryRealSolutionType(this);
        else if (solutionType.compareTo("Real") == 0)
            solutionType_ = new RealSolutionType(this);
        else {
            System.out.println("Error: solution type " + solutionType + " invalid");
            System.exit(-1);
        }
    }

    /**
     * Evaluates a solution
     *
     * @param solution The solution to evaluate
     * @throws JMException
     */


    public void evaluate(Solution solution) throws JMException {
        Variable[] gen = solution.getDecisionVariables();

        double[] x = new double[numberOfVariables_];
        double[] f = new double[numberOfObjectives_];

        for (int i = 0; i < numberOfVariables_; i++)
            x[i] = gen[i].getValue();

        double g = 1.0;
        for (int i = 1; i < numberOfVariables_; i++)
            g += x[i] * x[i] - 10.0 * Math.cos(4.0 * Math.PI * x[i]) + 10.0;

        f[0] = x[0];
        f[1] = g * Math.exp(-f[0] / g);

        for (int i = 0; i < numberOfObjectives_; i++)
            solution.setObjective(i, f[i]);
    } // evaluate

    /**
     * Evaluates the constraint overhead of a solution
     *
     * @param solution The solution
     * @throws JMException
     */
    public void evaluateConstraints(Solution solution) throws JMException {
        double f1 = solution.getObjective(0);
        double f2 = solution.getObjective(1);
        double constraint[] = new double[numberOfConstraints_];
        constraint[0] = f2 - 0.858 * Math.exp(-0.541 * f1);
        constraint[1] = f2 - 0.728 * Math.exp(-0.295 * f1);
        double total = 0.0;
        int number = 0;
        for (int i = 0; i < numberOfConstraints_; i++) {
            if (constraint[i] < 0.0) {
            	total = total - constraint[i];
                number++;
                solution.setConstraint(i,-constraint[i]);
            }else {
            	solution.setConstraint(i,0.0);
            }
        }
        solution.setOverallConstraintViolation(total);
//      solution.setOverallConstraintViolation(total/(double)numberOfConstraints_);
        solution.setNumberOfViolatedConstraint(number);
    } // evaluateConstraints

}


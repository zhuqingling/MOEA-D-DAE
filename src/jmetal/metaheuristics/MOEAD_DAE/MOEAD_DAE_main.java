//  Author: Qingling Zhu
//
package jmetal.metaheuristics.MOEAD_DAE;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.ProblemFactory;
import jmetal.problems.CTP.*;
import jmetal.problems.LIRCMOP.*;
import jmetal.problems.DTLZ.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.JMException;

import java.util.HashMap;

public class MOEAD_DAE_main {
  /**
   * @param args Command line arguments. The first (optional) argument specifies 
   *      the problem to solve.
   * @throws JMException 
   * @throws IOException 
   * @throws SecurityException 
   * Usage: three options
   *      - jmetal.metaheuristics.moead.MOEAD_main
   *      - jmetal.metaheuristics.moead.MOEAD_main problemName
   *      - jmetal.metaheuristics.moead.MOEAD_main problemName ParetoFrontFile
   * @throws ClassNotFoundException 
 
   */
  public static void main(String [] args) throws JMException, SecurityException, ClassNotFoundException {
		int runtimes=30;
		long[] estimatedTime=new long[runtimes];

		Problem problem=null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover ; // Crossover operator
		Operator mutation; // Mutation operator

		HashMap parameters; // Operator parameters

		QualityIndicator indicators; // Object to get quality indicators

		// Logger object and file to store log messages

		indicators = null;
		
		for(int i=0;i<runtimes;i++){
		
			if (args.length == 1) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
				indicators = new QualityIndicator(problem, args[1]);
			} // if
			else { // Default problem
				problem = new LIRCMOP4("Real");
			}

    algorithm = new MOEAD_DAE(problem);
    double F=0.5, CR=1.0;
    algorithm.setInputParameter("populationSize",100);
    algorithm.setInputParameter("maxEvaluations",1500*100);
    int SBX_flag=0; // 0: use DE operator; 1: use SBX operator
    F=0.5;CR=1.0;
    
    // Crossover operator 
    parameters = new HashMap() ;
    if(SBX_flag==1) {
    	parameters.put("probability", 0.2) ;
        parameters.put("distributionIndex", 20.0) ;
        crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
    }else {
    	parameters.put("CR", CR) ;
        parameters.put("F", F) ;
        crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters); 
    }
    
    // Mutation operator
    parameters = new HashMap();
    parameters.put("probability", 1.0/problem.getNumberOfVariables()) ;
    parameters.put("distributionIndex", 20.0) ;
    mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);                    
    
    algorithm.addOperator("crossover",crossover);
    algorithm.addOperator("mutation",mutation);
    
    // Execute the Algorithm
    long initTime = System.currentTimeMillis();
    SolutionSet population = algorithm.execute();
    estimatedTime[i] = System.currentTimeMillis() - initTime;
    population.printObjectivesToFile(problem.getName()+"_"+problem.getNumberOfObjectives()+"_"+"T"+(i+1));
    System.out.println(problem.getName()+"_"+problem.getNumberOfObjectives()+"_"+"T"+(i+1));
	}
  } //main
} // MOEAD_main

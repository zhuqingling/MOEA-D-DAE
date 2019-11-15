//  Original Version
// Better than MOEAD_TrialProcedureI

package jmetal.metaheuristics.MOEAD_DAE;
import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.qualityIndicator.util.MetricsUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;

public class MOEAD_DAE extends Algorithm {

  private int populationSize_;
  private double [] utility_;
  /**
   * Stores the population
   */
  private SolutionSet population_;
  private SolutionSet population_around;
  private SolutionSet archive;
  /**
   * Z vector (ideal point)
   */
  double[] z_;
  /**
   * Lambda vectors
   */
  //Vector<Vector<Double>> lambda_ ; 
  double[][] lambda_;
  /**
   * T: neighbour size
   */
  int T_;
  int H_ = 13;
  /**
   * Neighborhood
   */
  int[][] neighborhood_;
  /**
   * delta: probability that parent solutions are selected from neighbourhood
   */
  double delta_;
  double whole_cv;
  /**
   * nr: maximal number of solutions replaced by each child solution
   */
  int nr_;
  int TrialFlag=0;
  double sigma;
  double sigma_min;
  double avg_fit;
  String functionType_;
  int evaluations_;
  /**
   * Operators
   */
  Operator crossover_;
  Operator mutation_;

  MetricsUtil utils_ = new MetricsUtil();
  double change=Double.MAX_VALUE;

  /** 
   * Constructor
   * @param problem Problem to solve
   */
  public MOEAD_DAE(Problem problem) {
    super (problem) ;
    functionType_ = "_TCHE1";
  } // MOEA/D-DAE

  public SolutionSet execute() throws JMException, ClassNotFoundException {
    int maxEvaluations;
    
    evaluations_ = 0;
    maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
    populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();

    population_ = new SolutionSet(populationSize_);

    population_around=new SolutionSet(populationSize_);
    archive = new SolutionSet(populationSize_);

    utility_     = new double[populationSize_];

    T_ = 20;
    delta_ = 0.9;
    nr_ = 2;
    neighborhood_ = new int[populationSize_][T_];

    z_ = new double[problem_.getNumberOfObjectives()];
    lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];

    crossover_ = operators_.get("crossover");
    mutation_ = operators_.get("mutation");  // default: polynomial mutation

    // STEP 1. Initialization
    // STEP 1.1. Compute euclidean distances between weight vectors and find T
    initUniformWeight();
    
    initNeighborhood();

    // STEP 1.2. Initialize population
    initPopulation();

    // STEP 1.3. Initialize z_
    initIdealPoint(population_);
    whole_cv=0.0;
    double max_cv=Double.NEGATIVE_INFINITY;
    int feasible=0;
    double whole_fit=0.0;
    for (int i = 0; i < populationSize_; i++) {
    	whole_fit+=fitnessFunction(population_.get(i),lambda_[i]);
    	double cv = population_.get(i).getOverallConstraintViolation();
        if(cv==0)
    	    feasible +=1;
        if(cv>max_cv)
    	    max_cv=cv;
        whole_cv=whole_cv+cv;
    } // for
    avg_fit=(whole_fit/(double)populationSize_);
    int gen = 0;
    double f_r=(double)feasible/(double)populationSize_;
    
    double epsilon=max_cv*problem_.getNumberOfConstraints();
    double interval=1.0/3.0*(double)maxEvaluations/(double)populationSize_;
    sigma_min = 1.0-Math.pow(1.0/(epsilon+1.0e-10), 1.0/interval);
    SolutionSet offspring_pop=new SolutionSet(populationSize_);
    // STEP 2. Update
    do {
    	List<Integer> order = tour_selection(10);
        for (int i = 0; i < order.size(); i++) {
    	int n = order.get(i);
        int type;
        double rnd = PseudoRandom.randDouble();

        // STEP 2.1. Mating selection based on probability
        if (rnd < delta_) // if (rnd < realb)    
        {
          type = 1;   // neighborhood
        } else {
          type = 2;   // whole population
        }
        Vector<Integer> p = new Vector<Integer>();

        // STEP 2.2. Reproduction
        Solution child=new Solution();
        
        if(crossover_.getClass().getSimpleName().equalsIgnoreCase("SBXCrossover")) {
        	Solution[] parents = new Solution[2];
        	matingSelection(population_, p, n, 1, type);
        	parents[0] = population_.get(n);//CURRENT
        	parents[1] = population_.get(p.get(0));
            
	        Solution[] offSpring = (Solution[]) crossover_.execute(parents);
	        child=offSpring[0];
        }else if(crossover_.getClass().getSimpleName().equalsIgnoreCase("DifferentialEvolutionCrossover")){
        	Solution[] parents = new Solution[3];
        	matingSelection(population_, p, n, 2, type);
        	parents[0] = population_.get(p.get(0));
        	parents[1] = population_.get(p.get(1));
        	parents[2] = population_.get(n);//CURRENT
            
            // Apply DE crossover 
	        child = (Solution) crossover_.execute(new Object[]{population_.get(n), parents});
        }
        mutation_.execute(child);
        // Evaluation
        problem_.evaluate(child);
        problem_.evaluateConstraints(child); 
        evaluations_++;
        // STEP 2.3. Repair. Not necessary
        // STEP 2.4. Update z_
        updateReference(child);
        // STEP 2.5. Update of solutions
        if(child.getOverallConstraintViolation()==0) {
        	offspring_pop.add(child);
        }
        updateProblem(child, n, type, epsilon);
        if(TrialFlag==1) {
        	updateProblem(child, n, type);
        }
      }// for
      gen++;
      if(offspring_pop.size()>10 || (offspring_pop.size()>0 && archive.size()<1)) {
    	  ArchiveUpdate(offspring_pop);
    	  offspring_pop.clear();
      }
      feasible=0;
  	  for (int i = 0; i < populationSize_; i++) {
  		  if(population_.get(i).getOverallConstraintViolation()==0)
  			  feasible=feasible+1;
  	  }
  	  f_r=(double)feasible/(double)populationSize_;
  	  sigma = Math.max(sigma_min, f_r);
      if((TrialFlag==0 || TrialFlag==2)) {
	    	if(f_r<0.95){
	    		epsilon=(1.0-sigma)*epsilon;//improve diversity for the population
	    	}else {
	    		if(TrialFlag==0){
	    			//TrialProcedure
	    			TrialFlag=1;epsilon=1.0e+30;gen=1;
	    			//copy current population
	    			population_around.clear();
			    	for (int i = 0; i < populationSize_; i++) {
			    		population_around.add(new Solution(population_.get(i)));
			    	}
	    	    }else if(TrialFlag==2) {
	    	    	epsilon= max_cv;
	    	    	sigma_min = 1.0-Math.pow(1.0/(epsilon+1.0e-10), 1.0/interval);
	    	    }else {
		    		System.out.println("error and need to stop");
		    		return null;
	    	    }
	    	}
      }
      
      if(gen%10==0)
      {
    	  double whole_old=whole_cv;
          whole_cv=0.0;
          whole_fit=0.0;
          for (int i = 0; i < populationSize_; i++) {
        	whole_fit+=fitnessFunction(population_.get(i),lambda_[i]);
            whole_cv=whole_cv+Math.abs(population_.get(i).getOverallConstraintViolation());
          } // for
          avg_fit=whole_fit/(double)populationSize_;
          change=(double)Math.abs(whole_cv-whole_old)/((double)whole_cv+1.0e-10);
    	  if(TrialFlag==0) {//When trapped into local optimal
	          if((change>0&&change<1.0e-5)&&whole_cv>0.1*max_cv) {
	        	  //TrialProcedure
		    	  TrialFlag=1;epsilon=1.0e+30;gen=1;
		    	  //copy current population
		    	  population_around.clear();
			    	for (int i = 0; i < populationSize_; i++) {
			    		population_around.add(new Solution(population_.get(i)));
			    	}
	          }
    	  }else if(TrialFlag==1){
    		  if(change>0&&change<1.0e-5) {//increase f_r by change
    			TrialFlag=2;
		    	max_cv=0.0;
		    	for (int i = 0; i < populationSize_; i++) {
		    		if(Math.abs(population_around.get(i).getOverallConstraintViolation())>max_cv)
		    			max_cv=Math.abs(population_around.get(i).getOverallConstraintViolation());
		    	}
		    	population_.clear();
		    	for (int i = 0; i < populationSize_; i++) {
		    		population_.add(new Solution(population_around.get(i)));
		    	}
		    	epsilon=max_cv;
		    	sigma_min = 1.0-Math.pow(1.0/(epsilon+1.0e-10), 1.0/interval);
		    	//re-set the utility
		    	
		    	//re-initialize the ideal point
		    	initIdealPoint(population_);
		    }
          }
      }
    } while (evaluations_ < maxEvaluations);
    return archive;
  }
  public void ArchiveUpdate(SolutionSet offspring) {
	  SolutionSet union=archive.union(offspring);
	  archive.clear();
	//1. find nonDomianted solutions (Corner sort)
	  ParetoDominanceRanking ranking = new ParetoDominanceRanking(union);
	  SolutionSet archive_temp=ranking.getSubfront(0);
	  
	  int j=0;
	  int m=problem_.getNumberOfObjectives();
	  archive_temp.Suppress();
	  //
	  if(archive_temp.size()>=populationSize_) {
		  int l=archive_temp.size();
		  if(m==2) {
			  double[] f_max=new double[m];
			  double[] f_min=new double[m];
			  for(j=0;j<m;j++) {
				  f_max[j]=-1.0e+30;
				  f_min[j]=1.0e+30;
			  }
			  for(int i=0;i<l;i++) {
				  for(j=0;j<m;j++) {
					  if(archive_temp.get(i).getObjective(j)>f_max[j]) {
						  f_max[j]=archive_temp.get(i).getObjective(j);
					  }
					  if(archive_temp.get(i).getObjective(j)<f_min[j]) {
						  f_min[j]=archive_temp.get(i).getObjective(j);
					  }
				  }
			  }
			  while(archive_temp.size()>populationSize_) {
				  //crowding distance assignment
				  l=archive_temp.size();
				  //sort archive according to the first objective value
				  archive_temp.sort(new jmetal.util.comparators.ObjectiveComparator(0, false));//decending order ==false
				  double d_min=1.0e+30;
				  int min_idx=-1;
				  archive_temp.get(0).setFitness(1.0e+30);
				  archive_temp.get(l-1).setFitness(1.0e+30);
				  for(int i=1;i<l-1;i++) {
					  double d=0.0;
					  for(j=0;j<m;j++) {
						  d=d+Math.abs(archive_temp.get(i+1).getObjective(j)-archive_temp.get(i-1).getObjective(j))/(f_max[j]-f_min[j] + Double.MIN_VALUE);
					  }
					  archive_temp.get(i).setFitness(d);
					  if(d<d_min) {
						  d_min=d;
						  min_idx=i;
					  }
				  }
				  //remove the most crowded solution
				  archive_temp.remove(min_idx);
			  }
			  for(int i=0;i<archive_temp.size();i++) {
				  archive.add(archive_temp.get(i));
			  }
		  }else {//calculate the crowding distance with niching method in BCE
			  double [][] normalizedFront ;
			  double [][] PF = archive_temp.writeObjectivesToMatrix();//N*m
			  // STEP 1. Obtain the maximum and minimum values of the Pareto front
			  double [] maximumValue = utils_.getMaximumValues(PF, m);
			  double [] minimumValue = utils_.getMinimumValues(PF, m);
			  
			  // STEP 2. Get the normalized front and true Pareto fronts
			  normalizedFront = utils_.getNormalizedFront(PF,maximumValue,minimumValue);//N*m
			  
			  double [] solutionI, solutionJ;

			  //The matrix of distances
			  double [][] d_m = new double [l][l];
			  double [][] d_t = new double [l][l];
			  //-> Calculate the distances
			  for (int i = 0; i < l; i++){
				  d_m[i][i] = 0.0;
				  d_t[i][i] = 0.0;
			      solutionI = normalizedFront[i];
			      for (j = i + 1; j < normalizedFront.length; j++){
			    	  solutionJ = normalizedFront[j];
			    	  d_m[i][j] = utils_.distance(solutionI,solutionJ);
			    	  d_m[j][i] = d_m[i][j];
			    	  d_t[i][j] = d_m[i][j];
			    	  d_t[j][i] = d_t[i][j];
			      } // for
			  } // for  
			  int k = 2;
			  //int k = (int) Math.sqrt(solutionSet_.size()/2.0);
			  double sum_d3=0.0;
			  for (int i = 0; i < d_t.length; i++) {
			    Arrays.sort(d_t[i]);
			    sum_d3 = sum_d3 + d_t[i][k]; // Calcule de D(i) distance         
			  }// for
			  double r = sum_d3/(double)d_t.length;
			  
			  ArrayList<Integer> [] neighbor = new ArrayList [l];
			  double[] R = new double[l];
			  for(int i=0;i<l;i++) {
				  neighbor[i] = new ArrayList<Integer>();
			  }
			  for(int i=0;i<l;i++){
				  R[i]=1.0;
			  	  for(j=0; j<l;j++){
			  		  if(d_m[i][j] <= r && d_m[i][j]>0.0){
			  			  R[i] = R[i] * d_m[i][j]/r;
			  			  neighbor[i].add(j);//solution j is a neighbor of solution i
			  		  }
			  	  }
			  	  archive_temp.get(i).setFitness(R[i]);
			  }
			  ArrayList<Integer> removed_idx= new ArrayList<Integer>();
			  while(l-removed_idx.size()>populationSize_) {
				  //find the most crowded solution index in R
				  double crowded=1.0e+30;
				  int crowded_idx=-1;
				  for(int i=0;i<l;i++) {
					  if(!removed_idx.contains(i)&&R[i]<crowded) {
						  crowded=R[i];
						  crowded_idx=i;
					  }
				  }
				  //update the neighbors of the removed solution
				  for(int i=0;i<l;i++) {
					  if(neighbor[i].contains(crowded_idx)) {
						  R[i]=R[i]*r/d_m[i][crowded_idx];
						  int idx=neighbor[i].indexOf(crowded_idx);
//						  System.out.println(crowded_idx+"="+neighbor[i].get(idx));
						  neighbor[i].remove(idx);
						  archive_temp.get(i).setFitness(R[i]);
					  }
				  }
				  //add the most crowded index to removed_idx
				  removed_idx.add(crowded_idx);
			  }
			  for(int i=0;i<l;i++) {
				  if(!removed_idx.contains(i)) {
					  archive.add(archive_temp.get(i));
				  }
			  }
		  }
	  }else {//archive_temp.size is less than N
		  for(int i=0;i<archive_temp.size();i++) {
			  archive.add(archive_temp.get(i));
		  }
		  //assign fitness for the archive???
	  }
  }
  public List<Integer> tour_selection(int depth)
  {
	// selection based on utility
	List<Integer> selected = new ArrayList<Integer>();
	List<Integer> candidate = new ArrayList<Integer>();

//	for(int k=0; k<problem_.getNumberOfObjectives(); k++)
//     selected.add(k);   // WARNING! HERE YOU HAVE TO USE THE WEIGHT PROVIDED BY QINGFU (NOT SORTED!!!!)
 
	for(int n=0; n<populationSize_; n++)
		candidate.add(n);  // set of unselected weights

	while(selected.size()<(int)(populationSize_/5.0))
	{
	    //int best_idd = (int) (rnd_uni(&rnd_uni_init)*candidate.size()), i2;
           int best_idd = (int) (PseudoRandom.randDouble()*candidate.size());
           //System.out.println(best_idd);
           int i2;
	    int best_sub = candidate.get(best_idd);
           int s2;
           for(int i=1; i<depth; i++)
           {
               i2  = (int) (PseudoRandom.randDouble()*candidate.size());
               s2  = candidate.get(i2);
               //System.out.println("Candidate: "+i2);
               if(utility_[s2]>utility_[best_sub])
               {
                   best_idd = i2;
                   best_sub = s2;
               }
           }
	    selected.add(best_sub);
	    candidate.remove(best_idd);
	}
       return selected;
   }
  
  public void initUniformWeight() { // init lambda vectors
		int nw = 0;
		if (problem_.getNumberOfObjectives() == 2) {
			lambda_[nw][0] = 0.0;
			lambda_[nw][1] = 1.0;
			nw++;
			lambda_[nw][0] = 1.0;
			lambda_[nw][1] = 0.0;
			nw++;
			for (int n = 1; n < populationSize_-1; n++) {
				double a = 1.0 * n / (populationSize_ - 1);
				lambda_[nw][0] = a;
				lambda_[nw][1] = 1 - a;
				nw++;
			} // for
		} // if
		else {
			int i, j;
			lambda_[nw][0] = 0.0;
			lambda_[nw][1] = 0.0;
			lambda_[nw][2] = 1.0;
			nw++;
			lambda_[nw][0] = 0.0;
			lambda_[nw][1] = 1.0;
			lambda_[nw][2] = 0.0;
			nw++;
			lambda_[nw][0] = 1.0;
			lambda_[nw][1] = 0.0;
			lambda_[nw][2] = 0.0;
			nw++;
			for (i = 0; i <= H_; i++) {
				for (j = 0; j <= H_; j++) {
					if (i + j <= H_) {
						if((i==0 && j==0) || (i==0&&j==H_) || (i==H_&&j==0))
							continue;
						lambda_[nw][0] = (double) (1.0 * i) / H_;
						lambda_[nw][1] = (double) (1.0 * j) / H_;
						lambda_[nw][2] = (double) (1.0 * (H_ - i - j) / H_);
						nw++;
					} // if
				} // for
			} // for
		} // else

		if (nw != populationSize_) {
			System.out.println(nw + "---" + (populationSize_));
			System.out.println("ERROR: population size <> #weights");
			System.exit(0);
		}
		//Apply the WS-transformation on the generated weight vectors
//		if (problem_.getNumberOfObjectives()>2)
		for (int i=0;i<populationSize_;i++){
			double prod = 1.0, sum = 0.0;
			for (int j=0;j<problem_.getNumberOfObjectives();j++){
				prod = prod * lambda_[i][j];
			}
			if(prod != 0.0){
				for (int j=0;j<problem_.getNumberOfObjectives();j++){
					sum = sum + 1.0/lambda_[i][j];
				}
				for (int j=0;j<problem_.getNumberOfObjectives();j++){
					lambda_[i][j] = 1.0/lambda_[i][j]/sum;
				}
			}else{
				for (int j=0;j<problem_.getNumberOfObjectives();j++){
					sum = sum + 1.0/(lambda_[i][j]+0.0000001);
				}
				for (int j=0;j<problem_.getNumberOfObjectives();j++){
					lambda_[i][j] = 1.0/(lambda_[i][j]+0.0000001)/sum;
				}
			}
		}
	} // initUniformWeight
  
  /**
   * 
   */
  public void initNeighborhood() {
    double[] x = new double[populationSize_];
    int[] idx = new int[populationSize_];

    for (int i = 0; i < populationSize_; i++) {
      // calculate the distances based on weight vectors
      for (int j = 0; j < populationSize_; j++) {
        x[j] = Utils.distVector(lambda_[i], lambda_[j]);
        //x[j] = dist_vector(population[i].namda,population[j].namda);
        idx[j] = j;
      //System.out.println("x["+j+"]: "+x[j]+ ". idx["+j+"]: "+idx[j]) ;
      } // for

      // find 'niche' nearest neighboring subproblems
      Utils.minFastSort(x, idx, populationSize_, T_);
      //minfastsort(x,idx,population.size(),niche);

        System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
    } // for
  } // initNeighborhood

  /**
   * 
   */
  public void initPopulation() throws JMException, ClassNotFoundException {
    for (int i = 0; i < populationSize_; i++) {
      Solution newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);
      problem_.evaluateConstraints(newSolution);
      evaluations_++;
      population_.add(newSolution) ;
      utility_[i]=1;
    } // for
  } // initPopulation

  /**
   * 
   */
  void initIdealPoint(SolutionSet pop) throws JMException, ClassNotFoundException {
    for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
      z_[i] = 1.0e+30;
    } // for

    for (int i = 0; i < pop.size(); i++) {
      updateReference(pop.get(i));
//      utility_[i]=1;
    } // for
  } // initIdealPoint

  /**
   * 
   */
  public void matingSelection(SolutionSet pop, Vector<Integer> list, int cid, int size, int type) {
    // list : the set of the indexes of selected mating parents
    // cid  : the id of current subproblem
    // size : the number of selected mating parents
    // type : 1 - neighborhood; otherwise - whole population
    int ss;
    int r;
    int p;

    ss = neighborhood_[cid].length;
    while (list.size() < size) {
      if (type == 1) {
        r = PseudoRandom.randInt(0, ss - 1);
        p = neighborhood_[cid][r];
      //p = population[cid].table[r];
      } else {
        p = PseudoRandom.randInt(0, populationSize_ - 1);
      }
      boolean flag = true;
      for (int i = 0; i < list.size(); i++) {
//        if (list.get(i) == p) // p is in the list
    	if (pop.get(list.get(i)).getDecisionVariables().equals(pop.get(p).getDecisionVariables()))
        {
          flag = false;
          break;
        }
      }

      //if (flag) list.push_back(p);
      if (flag) {
        list.addElement(p);
      }
    }
  } // matingSelection

  /**
   * 
   * @param individual
   */
  void updateReference(Solution individual) {
    for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
      if (individual.getObjective(n) < z_[n]) {
        z_[n] = individual.getObjective(n);
      }
    }
  } // updateReference
  /**
   * @param individual
   * @param id
   * @param type
   */
  void updateProblem(Solution indiv, int id, int type, double epsilon) {
    // indiv: child solution
    // id:   the id of current subproblem
    // type: update solutions in - neighborhood (1) or whole population (otherwise)
    int size;
    int time;
    time = 0;

    if (type == 1) {
      size = neighborhood_[id].length;
    } else {
      size = population_.size();
    }
    int[] perm = new int[size];

    Utils.randomPermutation(perm, size);

    for (int i = 0; i < size; i++) {
      int k;
      if (type == 1) {
        k = neighborhood_[id][perm[i]];
      } else {
        k = perm[i];      // calculate the values of objective function regarding the current subproblem
      }
      double offspring_fit=fitnessFunction(indiv, lambda_[k]);
      double child_constraint=Math.abs(indiv.getOverallConstraintViolation());
      double parent_constraint=Math.abs(population_.get(k).getOverallConstraintViolation());
      if((child_constraint<=epsilon && parent_constraint<=epsilon) || child_constraint==parent_constraint) {
	      double f1, f2;
	      f1 = fitnessFunction(population_.get(k), lambda_[k]);
	      f2 = offspring_fit;
	      if (f2 < f1) {
	        population_.replace(k, new Solution(indiv));
	        time++;
	        double delta = (f1 - f2)/(f1+1.0e-10);
			if (delta > 0.001) {
				utility_[k] = 1.0;
			}else {
				double uti = 0.95 + (0.05 * delta / 0.001) * utility_[k];
//				double uti = 0.95 * (1.0 + delta / threshold) * utility_[k];//0.95 + (0.05 * delta / 0.001) * utility_[k];
				utility_[k] = uti < 1.0 ? uti : 1.0;
			}
	      }
      }else {
    	  double f1=fitnessFunction(population_.get(k), lambda_[k]);
    	  double f2=fitnessFunction(indiv, lambda_[k]);
	      f1 = sigma*f1+(1.0-sigma)*avg_fit*parent_constraint;
	      f2 = sigma*f2+(1.0-sigma)*avg_fit*child_constraint;
	
	      if (f2 < f1) {
	        population_.replace(k, new Solution(indiv));
	        time++;
	        double delta = (f1 - f2)/(f1+1.0e-10);
			if (delta > 0.001) {
				utility_[k] = 1.0;
			}else {
				double uti = 0.95 + (0.05 * delta / 0.001) * utility_[k];
//				double uti = 0.95 * (1.0 + delta / threshold) * utility_[k];//0.95 + (0.05 * delta / 0.001) * utility_[k];
				utility_[k] = uti < 1.0 ? uti : 1.0;
			}
	      }
      }
      // the maximal number of solutions updated is not allowed to exceed 'limit'
      if (time >= nr_) {
        return;
      }
    }
  } // updateProblem
void updateProblem(Solution indiv, int id, int type) {
  int size;
  int time;
  time = 0;
  if (type == 1) {
    size = neighborhood_[id].length;
  } else {
    size = population_around.size();
  }
  int[] perm = new int[size];

  Utils.randomPermutation(perm, size);

  for (int i = 0; i < size; i++) {
    int k;
    if (type == 1) {
      k = neighborhood_[id][perm[i]];
    } else {
      k = perm[i];      // calculate the values of objective function regarding the current subproblem
    }
    double f1, f2;
    double c1, c2;
    f1 = fitnessFunction(population_around.get(k), lambda_[k]);
    f2 = fitnessFunction(indiv, lambda_[k]);
    c1 = Math.abs(population_around.get(k).getOverallConstraintViolation());
    c2 = Math.abs(indiv.getOverallConstraintViolation());
    double e=1.0/populationSize_;
    if (e*f2+(1.0-e)*avg_fit*c2<e*f1+(1.0-e)*avg_fit*c1) {
  	  population_around.replace(k, new Solution(indiv));
  	  time++;
    }
    if (time >= nr_) {
        return;
      }
  }
} // updateProblem
	double fitnessFunction(Solution individual, double[] lamda) {
		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				// double diff = Math.abs(individual.getObjective(n)
				// - this.idealPoint[n]);

				double diff = Math.abs(individual.getObjective(n) - z_[n]);
				double feval;
				if (lamda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lamda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			return maxFun;
		} else if (functionType_.equals("_TCHE2")) {
            double maxFun = -1.0e+30;

            for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
                double diff = Math.abs(individual.getObjective(i) - z_[i]);

                double feval;
                if (lamda[i] == 0) {
                    feval = diff / 0.000001;
                } else {
                    feval = diff / lamda[i];
                }
                if (feval > maxFun) {
                    maxFun = feval;
                }
            } // for
            return maxFun;
        } else if (functionType_.equals("_WSUM")) {

			double sum = 0;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				sum += (lamda[n]) * individual.getObjective(n);
			}
			return sum;

		} // if
		else if (functionType_.equals("_NBI")) {
			int i;
			double d1, d2, nl;
			double theta = 5.0;
			double fin;

			d1 = d2 = nl = 0.0;
			for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
				d1 += (individual.getObjective(i) - z_[i]) * lamda[i];
				nl += Math.pow(lamda[i], 2.0);
			}
			d1 = Math.abs(d1) / Math.sqrt(nl);
			if (nl == 0.0) {
				System.out
						.println("ERROR: dived by zero(bad weihgted vector)\n");
				System.exit(0);
			}
			for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
				d2 += Math.pow((individual.getObjective(i) - z_[i]) - (d1 * lamda[i]), 2.0);
			}
			d2 = Math.sqrt(d2);
			fin = (d1 + theta * d2);
			return fin;
		}
		else {
			System.out.println("SDMOPSO.fitnessFunction: unknown type "
					+ functionType_);
			return 0;
		}
	} // fitnessEvaluation
} // MOEAD
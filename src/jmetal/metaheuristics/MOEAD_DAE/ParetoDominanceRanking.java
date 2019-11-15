//  Ranking.java
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

package jmetal.metaheuristics.MOEAD_DAE;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.OverallConstraintViolationComparator;

import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * This class implements some facilities for ranking solutions.
 * Given a <code>SolutionSet</code> object, their solutions are ranked 
 * according to scheme proposed in NSGA-II; as a result, a set of subsets 
 * are obtained. The subsets are numbered starting from 0 (in NSGA-II, the 
 * numbering starts from 1); thus, subset 0 contains the non-dominated 
 * solutions, subset 1 contains the non-dominated solutions after removing those
 * belonging to subset 0, and so on.
 */
public class ParetoDominanceRanking {
  
  /**
   * The <code>SolutionSet</code> to rank
   */
  private SolutionSet   solutionSet_ ;
  
  /**
   * An array containing all the fronts found during the search
   */
  private SolutionSet[] ranking_  ;
  
  /**
   * stores a <code>Comparator</code> for dominance checking
   */
//  private static final Comparator dominance_ = new DominanceComparator();
  
  /**
   * stores a <code>Comparator</code> for Overal Constraint Violation Comparator
   * checking
   */
//  private static final Comparator constraint_ = new OverallConstraintViolationComparator();
    
  /** 
   * Constructor.
   * @param solutionSet The <code>SolutionSet</code> to be ranked.
   * TrailFlag==1
   */       
  public ParetoDominanceRanking(SolutionSet solutionSet) {        
    solutionSet_ = solutionSet ;

    // dominateMe[i] contains the number of solutions dominating i        
    int [] dominateMe = new int[solutionSet_.size()];

    // iDominate[k] contains the list of solutions dominated by k
    List<Integer> [] iDominate = new List[solutionSet_.size()];

    // front[i] contains the list of individuals belonging to the front i
    List<Integer> [] front = new List[solutionSet_.size()+1];
        
    // flagDominate is an auxiliar encodings.variable
    int flagDominate;    

    // Initialize the fronts 
    for (int i = 0; i < front.length; i++)
      front[i] = new LinkedList<Integer>();
        
    /*
    //-> Fast non dominated sorting algorithm
    for (int p = 0; p < solutionSet_.size(); p++) {
    // Initialice the list of individuals that i dominate and the number
    // of individuals that dominate me
      iDominate[p] = new LinkedList<Integer>();
      dominateMe[p] = 0;            
      // For all q individuals , calculate if p dominates q or vice versa
      for (int q = 0; q < solutionSet_.size(); q++) {
        flagDominate =constraint_.compare(solutionSet.get(p),solutionSet.get(q));
        if (flagDominate == 0) {
          flagDominate =dominance_.compare(solutionSet.get(p),solutionSet.get(q));                                
        }
        
        if (flagDominate == -1) {
          iDominate[p].add(new Integer(q));
        } else if (flagDominate == 1) {
                    dominateMe[p]++;   
        }
      }
            
      // If nobody dominates p, p belongs to the first front
      if (dominateMe[p] == 0) {
        front[0].add(new Integer(p));
        solutionSet.get(p).setRank(0);
      }            
    }
    */
    
    //-> Fast non dominated sorting algorithm
    // Contribution of Guillaume Jacquenot
    for (int p = 0; p < solutionSet_.size(); p++) {
    // Initialize the list of individuals that i dominate and the number
    // of individuals that dominate me
      iDominate[p] = new LinkedList<Integer>();
      dominateMe[p] = 0;
    }
    for (int p = 0; p < (solutionSet_.size()-1); p++) {
      // For all q individuals , calculate if p dominates q or vice versa
      for (int q = p+1; q < solutionSet_.size(); q++) {
    	flagDominate =Dominancecompare(solutionSet.get(p),solutionSet.get(q));
        if (flagDominate == -1)
        {
          iDominate[p].add(q);
          dominateMe[q]++;
        }
        else if (flagDominate == 1)
        {
          iDominate[q].add(p);
          dominateMe[p]++;
        }
      }
      // If nobody dominates p, p belongs to the first front
    }
    for (int p = 0; p < solutionSet_.size(); p++) {
      if (dominateMe[p] == 0) {
        front[0].add(p);
        solutionSet.get(p).setRank(0);
      }
    }    
    
    //Obtain the rest of fronts
    int i = 0;
    Iterator<Integer> it1, it2 ; // Iterators
    while (front[i].size()!= 0) {
      i++;
      it1 = front[i-1].iterator();
      while (it1.hasNext()) {
        it2 = iDominate[it1.next()].iterator();
        while (it2.hasNext()) {
          int index = it2.next();
          dominateMe[index]--;
          if (dominateMe[index]==0) {
            front[i].add(index);
            solutionSet_.get(index).setRank(i);
          }
        }
      }
    }
    //<-
        
    ranking_ = new SolutionSet[i];
    //0,1,2,....,i-1 are front, then i fronts
    for (int j = 0; j < i; j++) {
      ranking_[j] = new SolutionSet(front[j].size());
      it1 = front[j].iterator();
      while (it1.hasNext()) {
                ranking_[j].add(solutionSet.get(it1.next()));
      }
    }
    
  } // Ranking

  /**
   * Returns a <code>SolutionSet</code> containing the solutions of a given rank. 
   * @param rank The rank
   * @return Object representing the <code>SolutionSet</code>.
   */
  public SolutionSet getSubfront(int rank) {
    return ranking_[rank];
  } // getSubFront

  /** 
  * Returns the total number of subFronts founds.
  */
  public int getNumberOfSubfronts() {
    return ranking_.length;
  } // getNumberOfSubfronts
  /**
   * Compares two solutions.
   * @param object1 Object representing the first <code>Solution</code>.
   * @param object2 Object representing the second <code>Solution</code>.
   * @return -1, or 0, or 1 if solution1 dominates solution2, both are 
   * non-dominated, or solution1  is dominated by solution22, respectively.
   */
   public int Dominancecompare(Object object1, Object object2) {
     if (object1==null)
       return 1;
     else if (object2 == null)
       return -1;
     
     Solution solution1 = (Solution)object1;
     Solution solution2 = (Solution)object2;

     int dominate1 ; // dominate1 indicates if some objective of solution1 
                     // dominates the same objective in solution2. dominate2
     int dominate2 ; // is the complementary of dominate1.

     dominate1 = 0 ; 
     dominate2 = 0 ;
     
     int flag; //stores the result of the comparison

     /*
     if (solution1.getOverallConstraintViolation()!= 
         solution2.getOverallConstraintViolation() &&
        (solution1.getOverallConstraintViolation() < 0) ||         
        (solution2.getOverallConstraintViolation() < 0)){            
       return (overallConstraintViolationComparator_.compare(solution1,solution2));
     }
    */
     
     // Equal number of violated constraints. Applying a dominance Test then
     double value1, value2;
     for (int i = 0; i < solution1.getNumberOfObjectives(); i++) {
       value1 = solution1.getObjective(i);
       value2 = solution2.getObjective(i);
       if (value1 < value2) {
         flag = -1;
       } else if (value1 > value2) {
         flag = 1;
       } else {
         flag = 0;
       }
       
       if (flag == -1) {
         dominate1 = 1;
       }
       
       if (flag == 1) {
         dominate2 = 1;           
       }
     }
             
     if (dominate1 == dominate2) {            
       return 0; //No one dominate the other
     }
     if (dominate1 == 1) {
       return -1; // solution1 dominate
     }
     return 1;    // solution2 dominate   
   } // compare
} // Ranking

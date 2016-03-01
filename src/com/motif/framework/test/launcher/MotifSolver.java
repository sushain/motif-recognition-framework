/**
* <code>MotifSolver</code> is a test Motif-solver solver
*
* @author pandit [sushain.pandit@gmail.com]
*
* @version 1.0
*
* Revision Log
*
* Date Version Description
**********************************************************************
*
* 15th Nov 1.0 First cut at MotifSolver
*
**********************************************************************
*/

package com.motif.framework.test.launcher;

import com.motif.framework.contract.IMotifSolvingEngine;
import com.motif.framework.factory.MotifSolvingEngineFactory;


public class MotifSolver {

	   public static void main(final String [] args) {
		   
		   String params [] = {"resources/testcase_aaa.txt", "35", "-1", "1", "20", "1", "true"};
		   
		   try {
			   
			   IMotifSolvingEngine instance = MotifSolvingEngineFactory.newPatternBranchingEngine(params);
			   
			   instance.processSequenceData();
		   
		   } catch (final Exception fileExp) {
			   
			   fileExp.printStackTrace();
		   }
		}
}

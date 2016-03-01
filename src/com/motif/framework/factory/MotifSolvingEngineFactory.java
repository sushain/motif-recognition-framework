/**
* <code>MotifSolvingEngineFactory</code> is a factory class that is 
* responsible for creating different Motif-solving engines.
* 
* The default Solving engine is <code>PatternBranchingEngine</code>.
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
* 15th Nov 1.0 First cut at MotifSolvingEngineFactory
*
**********************************************************************
*/

package com.motif.framework.factory;

import java.io.IOException;
import com.motif.framework.impl.patternbranching.engines.PatternBranchingEngine;


public final class MotifSolvingEngineFactory {


	public static PatternBranchingEngine newPatternBranchingEngine (final String [] params) throws IOException {
		
		return PatternBranchingEngine.getInstance(params);
	}
}
/**
* <code>IMotifSolvingEngine</code> defines the contract for any Motif
* Solving Engine.
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
* 15th Nov 1.0 First cut at IMotifSolvingEngine
*
**********************************************************************
*/

package com.motif.framework.contract;

import java.io.IOException;


public interface IMotifSolvingEngine {
	
	public void populateSequenceData(final String inputFile) throws IOException;
	
	public void processSequenceData();
}
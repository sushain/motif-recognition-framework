/**
* <code>TestSequenceGenerator</code>
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
* 15th Nov 1.0 First cut at TestSequenceGenerator
*
**********************************************************************
*/

package com.motif.framework.test.launcher;

import java.io.IOException;

import com.motif.framework.impl.patternbranching.engines.SequenceGeneratingEngine;


public class TestSequenceGenerator {

	   public static void main(final String [] args) {
		   
		   String params [] = {"resources/testcase_aaa.txt", "20", "2000", "35", "10", "true"};
		   
		   
		   SequenceGeneratingEngine instance = new SequenceGeneratingEngine (params);
		   
		   try {
			   
			   instance.generateSequence();
		   
		   } catch (final IOException fileExp) {
			   
			   fileExp.printStackTrace();
		   }
		}
}

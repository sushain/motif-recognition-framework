/**
* <code>Util</code>
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
* 15th Nov 1.0 First cut at the Utilities class
*
**********************************************************************
*/

package com.motif.framework.impl.patternbranching.util;

public class Util {

	public static final int WORST_SCORE = 10000;
	
	public static final int MAXIMUM_SEQUENCE_LENGTH = 5000; 

	public static final int MAXIMUM_NUMBER_OF_SEQUENCES = 21;

	public static final int MAXIMUM_MOTIF_LENGTH = 40; 

	public static final int MAXIMUM_NUMBER_OF_MUTATIONS = 10;

	public static final int MAXIMUM_NUMBER_OF_BEST_NEIGHBORS = 5;

	public static final int MAXIMUM_NUMBER_OF_TOP_MOTIFS = 5;
	
	public static char getNucleotide(final int representativeNumber) {

		char nucleotide = '\0';

		switch (representativeNumber) {

		case 0:
			nucleotide = 'A';
			break;

		case 1:
			nucleotide = 'C';
			break;

		case 2:
			nucleotide = 'G';
			break;

		case 3:
			nucleotide = 'T';
			break;
		}

		return (nucleotide);
	}

	public static int getNumericValue(final char nucleotide) {

		if ((new String("" + nucleotide).toUpperCase()).equals("A")) {

			return (0);

		} else if ((new String("" + nucleotide).toUpperCase()).equals("C")) {

			return (1);

		} else if ((new String("" + nucleotide).toUpperCase()).equals("G")) {

			return (2);

		} else if ((new String("" + nucleotide).toUpperCase()).equals("T")) {

			return (3);

		} else {

			return (-1);
		}
	}

	public static int randomize(final int multiplier) {

		float randomNumber = multiplier * ((float) Math.random());

		return (int) randomNumber;
	}
}

/**
* <code>SequenceGeneratingEngine</code>
*
* Reference - DOI: 10.1093/bioinformatics/btg1072, Price A., Ramabhadran S., Pevzner P.A. 2003.
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
* 15th Nov 1.0 First cut at SequenceGeneratingEngine
*
* 18th Nov 1.1 Finished coding and testing
*
**********************************************************************
*/

package com.motif.framework.impl.patternbranching.engines;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import com.motif.framework.impl.patternbranching.util.Util;

public class SequenceGeneratingEngine {

	private int numberOfSequences;

	private int lengthOfSequence;

	private int motifLength;

	private int maxDistance;

	private int [ ] sequence;

	private int [ ] motif;

	private int motifSite;

	private int unmutatedPosition;

	private String outputFile;

	private boolean debugFlag;

	public SequenceGeneratingEngine(final String [ ] params) {

		this.outputFile = params[0];

		this.numberOfSequences = Integer.parseInt(params[1]);

		this.lengthOfSequence = Integer.parseInt(params[2]);

		this.motifLength = Integer.parseInt(params[3]);

		this.maxDistance = Integer.parseInt(params[4]);

		this.sequence = new int [lengthOfSequence];

		this.motif = new int [motifLength];

		this.debugFlag = new Boolean(params[5]);
	}

	public void generateSequence () throws IOException {

		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));

		// Generate random motif
		for (int motifNucleotideIndex = 0; motifNucleotideIndex < motifLength; motifNucleotideIndex++) {

			motif[motifNucleotideIndex] = Util.randomize(4);

			if (debugFlag) {

				System.out.print(Util
						.getNucleotide(motif[motifNucleotideIndex]));
			}
		}
		
		System.out.println();

		// Generate sequences
		for (int sequenceIndex = 0; sequenceIndex < numberOfSequences; sequenceIndex++) {

			for (int sequenceNucleotideIndex = 0; sequenceNucleotideIndex < lengthOfSequence; sequenceNucleotideIndex++) {

				sequence[sequenceNucleotideIndex] = Util.randomize(4);
			}

			// Insert the motif
			motifSite = Util.randomize(lengthOfSequence - motifLength + 1);

			if (debugFlag) {

				System.out.println((sequenceIndex + 1) + " - " + motifSite);
			}

			for (int motifNucleotideIndex = 0; motifNucleotideIndex < motifLength; motifNucleotideIndex++) {

				sequence[motifSite + motifNucleotideIndex] = motif[motifNucleotideIndex];
			}

			// Mutate the motif
			for (int hammingDistanceCounter = 0; hammingDistanceCounter < maxDistance; hammingDistanceCounter++) {

				unmutatedPosition = Util.randomize(motifLength);

				while (sequence[motifSite + unmutatedPosition] != motif[unmutatedPosition]) {

					unmutatedPosition = Util.randomize(motifLength);
				}

				sequence[motifSite + unmutatedPosition] = (sequence[motifSite
						+ unmutatedPosition]
						+ Util.randomize(3) + 1) % 4;
			}

			for (int sequenceNucleotideIndex = 0; sequenceNucleotideIndex < lengthOfSequence; sequenceNucleotideIndex++) {

				writer.append(Util
						.getNucleotide(sequence[sequenceNucleotideIndex]));
			}

			writer.append("\n");
		}

		writer.flush();

		writer.close();
	}

	public int getNumberOfSequences () {

		return numberOfSequences;
	}

	public void setNumberOfSequences (final int numberOfSequences) {

		this.numberOfSequences = numberOfSequences;
	}

	public int getLengthOfSequence () {

		return lengthOfSequence;
	}

	public void setLengthOfSequence (final int lengthOfSequence) {

		this.lengthOfSequence = lengthOfSequence;
	}

	public int getMotifLength () {

		return motifLength;
	}

	public void setMotifLength (final int motifLength) {

		this.motifLength = motifLength;
	}

	public int getMaxDistance () {

		return maxDistance;
	}

	public void setMaxDistance (final int maxDistance) {

		this.maxDistance = maxDistance;
	}

	public boolean getDebugFlag () {

		return debugFlag;
	}

	public void setDebugFlag (final boolean debugFlag) {

		this.debugFlag = debugFlag;
	}
}

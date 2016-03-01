/**
* <code>PatternBranchingEngine</code>
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
* 15th Nov 1.0 First cut at PatternBranchingEngine
*
* 16th Nov 1.1 Modified code to add worker code
*
* 19th Nov 1.2 Modified code to add heuristic speed-up
* 
**********************************************************************
*/

package com.motif.framework.impl.patternbranching.engines;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import com.motif.framework.contract.IMotifSolvingEngine;
import com.motif.framework.impl.patternbranching.util.Util;


public class PatternBranchingEngine implements IMotifSolvingEngine {

	private int [] sequenceLength; 

	private int numberOfSequences; 

	private int motifLength; 

	private int numberOfMutations; 

	private int numberOfCandidateNeighbors; 

	private int numberOfTopMotifs; 

	private int nref;

	private int distanceThreshold;

	private int[][] sequenceMatrix; 

	private int[][][][] distance;

	private int[][][][] totalDistance;

	private int[][][] setOfCorrectedMotifs; 

	private int[][] correctionMatrix;

	private int[] rfrom;

	private int ref;

	private int[][] goodpos; 

	private int[] positionsBetterThanDistanceThreshold;

	private int[] M;

	private int motifScore;

	private int[][] setOfTopMotifs;

	private int[] setOfTopMotifScores;

	private int ncorr;

	private int posA;

	private int worstMotifScore;

	private int positionOfWorstMotifScore;

	private boolean debugFlag;
	
	private static PatternBranchingEngine instance;
	
	private PatternBranchingEngine(final String[] params) throws IOException {
		
		this.sequenceLength = new int[Util.MAXIMUM_NUMBER_OF_SEQUENCES];
		
        this.sequenceMatrix = new int[Util.MAXIMUM_NUMBER_OF_SEQUENCES][Util.MAXIMUM_SEQUENCE_LENGTH]; 

        this.distance = new int[Util.MAXIMUM_NUMBER_OF_MUTATIONS + 1][Util.MAXIMUM_NUMBER_OF_BEST_NEIGHBORS][Util.MAXIMUM_NUMBER_OF_SEQUENCES][Util.MAXIMUM_SEQUENCE_LENGTH + 1];

		this.totalDistance = new int[Util.MAXIMUM_NUMBER_OF_MUTATIONS][Util.MAXIMUM_NUMBER_OF_BEST_NEIGHBORS][Util.MAXIMUM_MOTIF_LENGTH][4];

		this.setOfCorrectedMotifs = new int[Util.MAXIMUM_NUMBER_OF_MUTATIONS + 1][Util.MAXIMUM_NUMBER_OF_BEST_NEIGHBORS][Util.MAXIMUM_MOTIF_LENGTH]; 

		this.correctionMatrix = new int[Util.MAXIMUM_NUMBER_OF_MUTATIONS][Util.MAXIMUM_NUMBER_OF_BEST_NEIGHBORS];

		this.rfrom = new int[Util.MAXIMUM_NUMBER_OF_BEST_NEIGHBORS];

		this.goodpos = new int[Util.MAXIMUM_NUMBER_OF_SEQUENCES][Util.MAXIMUM_SEQUENCE_LENGTH + 1];

		this.positionsBetterThanDistanceThreshold = new int[Util.MAXIMUM_NUMBER_OF_SEQUENCES];

		this.M = new int[Util.MAXIMUM_MOTIF_LENGTH];

		this.setOfTopMotifs = new int[Util.MAXIMUM_NUMBER_OF_TOP_MOTIFS][Util.MAXIMUM_MOTIF_LENGTH];

		this.setOfTopMotifScores = new int[Util.MAXIMUM_NUMBER_OF_TOP_MOTIFS];
		
		populateSequenceData(params[0]); 

		this.nref = this.numberOfSequences;

		this.motifLength = Integer.parseInt(params[1]);

		this.numberOfMutations = Integer.parseInt(params[2]);

		this.numberOfCandidateNeighbors = Integer.parseInt(params[3]);

		this.nref = Integer.parseInt(params[4]);

		this.numberOfTopMotifs = Integer.parseInt(params[5]);

		this.debugFlag = new Boolean(params[6]);
	}

	public static PatternBranchingEngine getInstance (final String [] params) throws IOException {
	
		if (instance == null) {
			
			instance = new PatternBranchingEngine (params);
		}
		
		return instance;
	}
	
	public void processSequenceData() {

		int bestocc, dist, bestdist, pos;

		generateNumberOfMutationsAndDistanceThreshold (); // Borrowed from Price A., Ramabhadran S., Pevzner P.A. 2003

		initializeBestMotifs (); 
		
		long startTime = System.currentTimeMillis();
		
		for (ref = 0; ref < nref; ref++) {
	
			for (posA = 0; posA < sequenceLength[ref] - motifLength + 1; posA++) {
				
				this.motifScore = Util.WORST_SCORE;

				buildDistance(); 
				
				for (int motifIndex = 0; motifIndex < motifLength; motifIndex++) {
					
					setOfCorrectedMotifs[0][0][motifIndex] = sequenceMatrix[ref][posA + motifIndex];
				}

				for (ncorr = 0; ncorr <= numberOfMutations; ncorr++) {
					
					patternBranchingWorker(); 
					
					if (ncorr < numberOfMutations) {
						
						calculateSetOfBestNeighbors();  // Borrowed from Price A., Ramabhadran S., Pevzner P.A. 2003
					}
				}
				
				updateBestMotifs(); 
			}
		}

		long endTime = System.currentTimeMillis();

		System.out.println("Runtime = " + (endTime - startTime) + " milliseconds");
		
		motifScore = Util.WORST_SCORE;

		int motifIndex;
		
		System.out.println("Recognized Motif - ");
		
		for (int topMotifIndex = 0; topMotifIndex < numberOfTopMotifs; topMotifIndex++) {
	
			for (motifIndex = 0; motifIndex < motifLength; motifIndex++) {

				System.out.print(Util.getNucleotide(setOfTopMotifs[topMotifIndex][motifIndex]));
			}

			System.out.println("\nMotif Score - " + setOfTopMotifScores[topMotifIndex]);
			
			System.out.println("Motif positions in sequences are as follows: -");

			if (debugFlag) {

				for (int sequenceIndex = 0; sequenceIndex < numberOfSequences; sequenceIndex++) {
					
					bestocc = -1;
					
					bestdist = 1000;
					
					for (pos = 0; pos < sequenceLength[sequenceIndex] - this.motifLength + 1; pos++) {
						
						dist = 0;
						
						for (motifIndex = 0; motifIndex < this.motifLength; motifIndex++) {
							
							if (setOfTopMotifs[topMotifIndex][motifIndex] != sequenceMatrix[sequenceIndex][pos + motifIndex]) {
								
								dist++;
							}
						}
						
						if (dist >= bestdist) {
							
							continue;	
						}
						
						bestdist = dist;
						
						bestocc = pos;
					}
					
					System.out.println((sequenceIndex + 1) + " - " + bestocc);
				}
				
				System.out.println();
			}
		}
	}

	private void initializeBestMotifs () {

		for (int motifIndex = 0; motifIndex < numberOfTopMotifs; motifIndex++) {

			setOfTopMotifScores[motifIndex] = Util.WORST_SCORE;
		}

		worstMotifScore = Util.WORST_SCORE;

		positionOfWorstMotifScore = 0;
	}
	
	private void patternBranchingWorker () {

		int[][] alpha = new int[this.motifLength][this.numberOfMutations];

		int i1, ci1, pos, i, x, g, rr, thisr, a;

		int totaldist0, totaldist0i;

		a = distanceThreshold - ncorr;
		
		thisr = this.numberOfCandidateNeighbors;
		
		if (ncorr == 0) {
		
			thisr = 1;
		}

		if (ncorr > numberOfMutations) {
			
			System.exit(1);
		}

		for (rr = 0; rr < thisr; rr++) {
			
			totaldist0 = 0;
			
			if (ncorr < numberOfMutations) {
			
				for (i1 = 0; i1 < motifLength; i1++) {
				
					for (ci1 = 0; ci1 < 4; ci1++) {
					
						totalDistance[ncorr][rr][i1][ci1] = (numberOfSequences - 1) * (a + 1);
					}
				}
			}

			for (i = 0; i < numberOfSequences; i++) {
				
				totaldist0i = 1000;
			
				if (i == ref) {
	
					continue;
				}
				
				for (i1 = 0; i1 < motifLength; i1++) {
					
					for (ci1 = 0; ci1 < 4; ci1++) {
						
						alpha[i1][ci1] = a + 1;
					}
				}
					
				if (ncorr > 0)  {
	
					for (g = 0; g < positionsBetterThanDistanceThreshold[i]; g++) {
						
						pos = goodpos[i][g];
						
						distance[ncorr][rr][i][pos] = distance[ncorr - 1][rfrom[rr]][i][pos];
						
						if (setOfCorrectedMotifs[ncorr - 1][rfrom[rr]][correctionMatrix[ncorr - 1][rr]] != sequenceMatrix[i][correctionMatrix[ncorr - 1][rr] + pos]) {
						
							distance[ncorr][rr][i][pos] -= 1;
						}
						
						if (setOfCorrectedMotifs[ncorr][rr][correctionMatrix[ncorr - 1][rr]] != sequenceMatrix[i][correctionMatrix[ncorr - 1][rr] + pos]) {
						
							distance[ncorr][rr][i][pos] += 1;
						}
						
						if (distance[ncorr][rr][i][pos] < totaldist0i) {
						
							totaldist0i = distance[ncorr][rr][i][pos];
						}
	
						if (ncorr < numberOfMutations) {
	
							if (distance[ncorr][rr][i][pos] > distanceThreshold - ncorr) {
								
								continue;
							}
								
							for (i1 = 0; i1 < motifLength; i1++) {
								
								ci1 = sequenceMatrix[i][i1 + pos];
								
								if (setOfCorrectedMotifs[ncorr][rr][i1] == ci1) {
								
									continue;
								}
								
								
								if (alpha[i1][ci1] > distance[ncorr][rr][i][pos]) {
									
									totalDistance[ncorr][rr][i1][ci1] += distance[ncorr][rr][i][pos] - alpha[i1][ci1];
									
									alpha[i1][ci1] = distance[ncorr][rr][i][pos];
								}
							}
						}
					}
					
				} else {
					
					for (pos = 0; pos < sequenceLength[i] - motifLength + 1; pos++) {
						
						if (distance[ncorr][rr][i][pos] < totaldist0i) {
						
							totaldist0i = distance[ncorr][rr][i][pos];
						}
						
						if (distance[ncorr][rr][i][pos] > distanceThreshold - ncorr) {
							
							continue;
						}
						
						for (i1 = 0; i1 < motifLength; i1++) {
						
							ci1 = sequenceMatrix[i][i1 + pos];
							
							if (setOfCorrectedMotifs[ncorr][rr][i1] == ci1) {
								
								continue;
							}
							
							if (alpha[i1][ci1] > distance[ncorr][rr][i][pos]) {
								
								totalDistance[ncorr][rr][i1][ci1] += distance[ncorr][rr][i][pos] - alpha[i1][ci1];
								
								alpha[i1][ci1] = distance[ncorr][rr][i][pos]; 
							}
						}
					}
				}
	
				if (totaldist0i > numberOfMutations) {
					
					totaldist0++;
				}
			} 
			
			if (totaldist0 < motifScore) {
				
				motifScore = totaldist0;
				
				for (x = 0; x < motifLength; x++) {
				
					M[x] = setOfCorrectedMotifs[ncorr][rr][x];
				}
			}	
		}
			
		return;
	}

	private void calculateSetOfBestNeighbors () {  // Borrowed from Price A., Ramabhadran S., Pevzner P.A. 2003

		int thisr;
		
		int[][] topr = new int[this.numberOfCandidateNeighbors][3]; 
		
		int[] toprscore = new int[this.numberOfCandidateNeighbors];
		
		int oldrr, rr, worstofrscore, worstofrpos, i1, ci1, x;

		thisr = this.numberOfCandidateNeighbors;
		
		if (ncorr == 0) {
			
			thisr = 1;
		}
		
		for (rr = 0; rr < this.numberOfCandidateNeighbors; rr++) {
			
			toprscore[rr] = Util.WORST_SCORE;	
		}
		
		worstofrscore = Util.WORST_SCORE;
		
		worstofrpos = 0;
		
		for (oldrr = 0; oldrr < thisr; oldrr++) {
			
			for (i1 = 0; i1 < motifLength; i1++) {
				
				for (ci1 = 0; ci1 < 4; ci1++) {
					
					if (totalDistance[ncorr][oldrr][i1][ci1] < worstofrscore) {
						
						toprscore[worstofrpos] = totalDistance[ncorr][oldrr][i1][ci1];
						
						topr[worstofrpos][0] = oldrr;
						
						topr[worstofrpos][1] = i1;
						
						topr[worstofrpos][2] = ci1;
						
						worstofrpos = 0;
						
						worstofrscore = toprscore[worstofrpos];
						
						for (rr = 1; rr < this.numberOfCandidateNeighbors; rr++) {
						
							if (toprscore[rr] > worstofrscore) {
							
								worstofrscore = toprscore[rr];
								
								worstofrpos = rr;
							}
						}
					}
				}
			}
		}
		
		for (rr = 0; rr < this.numberOfCandidateNeighbors; rr++) {
			
			rfrom[rr] = topr[rr][0];
			
			for (x = 0; x < motifLength; x++) {
				
				setOfCorrectedMotifs[ncorr + 1][rr][x] = setOfCorrectedMotifs[ncorr][rfrom[rr]][x];
			}
			
			i1 = topr[rr][1];
			
			ci1 = topr[rr][2];
			
			setOfCorrectedMotifs[ncorr + 1][rr][i1] = ci1;
			
			correctionMatrix[ncorr][rr] = i1;
		}
	}
	
	private void updateBestMotifs() {

		if (motifScore >= worstMotifScore) {

			return;
		}

		setOfTopMotifScores[positionOfWorstMotifScore] = motifScore;

		for (int motifIndex = 0; motifIndex < motifLength; motifIndex++) {

			setOfTopMotifs[positionOfWorstMotifScore][motifIndex] = M[motifIndex];
		}

		worstMotifScore = setOfTopMotifScores[0];

		positionOfWorstMotifScore = 0;

		for (int index = 1; index < numberOfTopMotifs; index++) {

			if (setOfTopMotifScores[index] <= worstMotifScore) {

				continue;
			}

			worstMotifScore = setOfTopMotifScores[index];

			positionOfWorstMotifScore = index;
		}
	}

	public void populateSequenceData(final String inputFile) throws IOException {

		int i, pos = 0;

		String readSequence;

		char[] sequenceChars;

		BufferedReader reader = new BufferedReader(new FileReader(inputFile));

		i = 0;

		while ((readSequence = reader.readLine()) != null) {

			sequenceChars = readSequence.toCharArray();

			pos = 0;

			for (int count = 0; count < sequenceChars.length; count++) {

				sequenceMatrix[i][pos++] = Util
						.getNumericValue(sequenceChars[count]);
			}

			sequenceLength[i] = pos;

			i++;
		}

		numberOfSequences = i;
	}

	private void buildDistance() {

		int pos;

		for (int index = 0; index < numberOfSequences; index++) {
			
			positionsBetterThanDistanceThreshold[index] = 0;
		}

		if (posA == 0) {
			
			for (int index = 0; index < numberOfSequences; index++) {
				
				if (index == ref) {
				
					continue;
				}
				
				for (pos = 0; pos < sequenceLength[index] - motifLength + 1; pos++) {
					
					distance[0][0][index][pos] = 0;
					
					for (int x = 0; x < motifLength; x++) {
					
						if (sequenceMatrix[index][x + pos] != sequenceMatrix[ref][posA+ x]) {
						
							distance[0][0][index][pos] += 1;
						}
					}
					
					if (distance[0][0][index][pos] > distanceThreshold) {
					
						continue;
					}
					
					goodpos[index][positionsBetterThanDistanceThreshold[index]] = pos;
					
					positionsBetterThanDistanceThreshold[index] += 1;
				}
			}
			
		} else {
			
			for (int index = 0; index < numberOfSequences; index++) {
			
				if (index == ref) {
					
					continue;
				}
				
				for (pos = sequenceLength[index] - motifLength; pos > 0; pos--) {
					
					distance[0][0][index][pos] = distance[0][0][index][pos - 1];
					
					if (sequenceMatrix[index][motifLength - 1 + pos] != sequenceMatrix[ref][posA+ motifLength - 1]) {
						
						distance[0][0][index][pos] += 1;
					}
					
					if (sequenceMatrix[index][pos - 1] != sequenceMatrix[ref][posA - 1]) {
					
						distance[0][0][index][pos] -= 1;
					}
					
					if (distance[0][0][index][pos] > distanceThreshold) {
						
						continue; 
					}
					
					goodpos[index][positionsBetterThanDistanceThreshold[index]] = pos;
					
					positionsBetterThanDistanceThreshold[index] += 1;
				}
				
				pos = 0;
				
				distance[0][0][index][pos] = 0;
				
				for (int x = 0; x < motifLength; x++) {
					
					if (sequenceMatrix[index][x + pos] != sequenceMatrix[ref][posA + x]) {
					
						distance[0][0][index][pos] += 1;
					}
				}
				
				if (distance[0][0][index][pos] <= distanceThreshold) {
					
					goodpos[index][positionsBetterThanDistanceThreshold[index]] = pos;
					
					positionsBetterThanDistanceThreshold[index] += 1;
				}
			}
		}
	}
	
	private void generateNumberOfMutationsAndDistanceThreshold() { // Borrowed from Price A., Ramabhadran S., Pevzner P.A. 2003

		this.numberOfMutations = (3 * this.motifLength) / 10;

		if (motifLength <= 17) {

			distanceThreshold = (motifLength + 1) / 2;

		} else if (motifLength <= 20) {

			distanceThreshold = (motifLength + 2) / 2;

		} else if (motifLength <= 23) {

			distanceThreshold = (motifLength + 3) / 2;

		} else if (motifLength <= 26) {

			distanceThreshold = (motifLength + 4) / 2;

		} else if (motifLength <= 29) {

			distanceThreshold = (motifLength + 5) / 2;

		} else {

			distanceThreshold = (motifLength + 6) / 2;
		}
	}
}

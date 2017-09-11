/*
# GAPE - Genetic Algorithm for Parameter Estimation of SIFTER tool
#
# Created by Eng. (C) Julian Camilo Castañeda Alonso and Msc. Tania Andrea Rodriguez Quiñones on August 2017.
# Copyright (c) 2017. Eng. (C) Julian Camilo Castañeda Alonso and Msc. Tania Andrea Rodriguez Quiñones. Universidad Antonio Narino. All rights reserved.
#
# GAPE is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License v3.0 found in the LICENSE file in the root directory of this project.
*/

package gape.markov_model; 

/**
 * This class represents an Individual for Genetic Algorithm with the parameters for SIFTER
 * @author Eng. (C) Julian Camilo Castañeda Alonso - Msc. Tania Andrea Rodriguez Quiñones
 *
 */

public class MarkovModel {
	
	/**
	 * This method build a Transtion Matrix
	 * @param numStates
	 * @return matrix in an vector
	 */
	public static double[][] createTransitionMatrix(int numStates) {
		double [][] transitionMatrix = new double[numStates][numStates];
		double [] sumRowTransitionMatrix = new double[numStates];
		
		double sum = 0; 
		
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++) {
				transitionMatrix[i][j] = Math.random();
				sum += transitionMatrix[i][j];
			}
			sumRowTransitionMatrix[i] = sum;
			sum = 0;
		}
		
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++) {
				transitionMatrix[i][j] = transitionMatrix[i][j] / sumRowTransitionMatrix[i];
			}
		}
		
		return transitionMatrix;
	}
	
	/**
	 * This method build a Alpha Matrix
	 * @param numStates
	 * @return matrix in an vector
	 */
	public static double[] createAlpha(int numStates) {
		double [] alpha = new double[numStates];
		
		for (int i = 0; i < alpha.length; i++) {
			alpha[i] = Math.random();
		}
		
		return alpha;
	}
	
	/**
	 * This method build a Sigma Matrix with speciation and duplication events
	 * @return matrix in an vector
	 */

	public static double[] createSigma() {
		double [] phi = new double[2];
		
		double random_1 = Math.random();
		double random_2 = Math.random();
		
		while (random_2 <= random_1) {
			random_2 = Math.random();
		}
				
		phi[0] = random_1;
		phi[1] = random_2;
		
		return phi;
	}
	
}

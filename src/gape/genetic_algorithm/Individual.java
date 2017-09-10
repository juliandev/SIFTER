/*
# GAPE - Genetic Algorithm for Parameter Estimation of SIFTER tool
#
# Created by Eng. (C) Julian Camilo Castañeda Alonso and Msc. Tania Andrea Rodriguez Quiñones on August 2017.
# Copyright (c) 2017. Eng. (C) Julian Camilo Castañeda Alonso and Msc. Tania Andrea Rodriguez Quiñones. Universidad Antonio Narino. All rights reserved.
#
# GAPE is free software: you can redistribute it and/or modify it under the terms of the 
# Apache License 2.0 found in the LICENSE file in the root directory of this project.
*/

package gape.genetic_algorithm;

/**
 * 
 * @author Eng. (C) Julian Camilo Castañeda Alonso - Msc. Tania Andrea Rodriguez Quiñones
 *
 */

public class Individual {
	
	// Attributes
	private double[][] transitionMatrix;
	private double[] alpha;
	private double[] sigma;
	private double fitness;
	
	/**
	 * Constructor of the class that requires the genetic representation of new individual
	 * @param transitionMatrix
	 * @param alpha
	 * @param sigma
	 * @param fitness
	 */
	
	public Individual(double[][] transitionMatrix, double[] alpha, double[] sigma, double fitness) {
		this.transitionMatrix = transitionMatrix; // Gene 1
		this.alpha = alpha; // Gene 2
		this.sigma = sigma; // Gene 3
		this.fitness = fitness;
	}
		
	/**
	 * @return the transitionMatrix
	 */
	public double[][] getTransitionMatrix() {
		return transitionMatrix;
	}

	/**
	 * @param transitionMatrix the transitionMatrix to set
	 */
	public void setTransitionMatrix(double[][] transitionMatrix) {
		this.transitionMatrix = transitionMatrix;
	}

	/**
	 * @return the alpha
	 */
	public double[] getAlpha() {
		return alpha;
	}

	/**
	 * @param alpha the alpha to set
	 */
	public void setAlpha(double[] alpha) {
		this.alpha = alpha;
	}

	/**
	 * @return the sigma
	 */
	public double[] getSigma() {
		return sigma;
	}

	/**
	 * @param sigma the sigma to set
	 */
	public void setSigma(double[] sigma) {
		this.sigma = sigma;
	}

	/**
	 * @return the fitness
	 */
	public double getFitness() {
		return fitness;
	}

	/**
	 * @param fitness the fitness to set
	 */
	public void setFitness(double fitness) {
		this.fitness = fitness;
	}

	/**
	 * @return the String of Tansition Matrix
	 */
	public String getStringTransitionMatrix() {
		
		String matrix = "";
		
		for (int i = 0; i < transitionMatrix.length; i++) {
			
			for (int j = 0; j <transitionMatrix[i].length; j++) {
				
				matrix += transitionMatrix[i][j] + "\t";
				
			}
			
			matrix += "\n";
			
		}
		
		return matrix;	
		
	}
	
	/**
	 * @return the String of Alpha Matrix
	 */
	public String getStringAlphaMatrix() {
		
		String matrix = "";
		
		for (int i = 0; i < alpha.length; i++) {
			matrix += alpha[i] + " ";
		}
		
		return matrix;	
		
	}
	
	/**
	 * @return the String of Sigma Matrix
	 */
	public String getStringSigmaMatrix() {
		
		String matrix = "";
		
		for (int i = 0; i < sigma.length; i++) {
			matrix += sigma[i] + " ";
		}
		
		return matrix;	
		
	}
	
}
/*
# GAPE - Genetic Algorithm for Parameter Estimation of SIFTER tool
#
# Created by Eng. (C) Julian Camilo Castañeda Alonso, Msc. Carlos Andres Sierra and Msc. Tania Andrea Rodriguez Quiñones on August 2017.
# Copyright (c) 2017. Eng. (C) Julian Camilo Castañeda Alonso, Msc. Carlos Andres Sierra and Msc. Tania Andrea Rodriguez Quiñones. Universidad Antonio Narino. All rights reserved.
#
# GAPE is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License v3.0 found in the LICENSE file in the root directory of this project.
*/

package gape.genetic_algorithm;

/**
 * This class represents fitness function to genetic algorithm
 * 
 * @author Eng. (C) Julian Camilo Castañeda Alonso, Msc. Carlos Andres Sierra and Msc. Tania Andrea Rodriguez Quiñones
 *
 */

public class FitnessFunction {

	/**
	 * Default constructor
	 */
	public FitnessFunction() {

	}

	/**
	 * This method calculate fitness of individual
	 * 
	 * @param transitionMatrix
	 * @param alphaArray
	 * @param phiArray
	 * @return
	 */
	public double fitness(double[][] transitionMatrix, double[] alphaArray, double[] phiArray) {

		double fitness = 0;
		double alphaM = alphaArray(alphaArray);
		double phiM = sigmaArray(phiArray);
		double transitionM = transitionMatrix(transitionMatrix);

		fitness = transitionM + alphaM + phiM;

		return fitness;

	}

	/**
	 * This method calculate the Log-Likelihood function to the Transition
	 * Matrix
	 * 
	 * @param transitionMatrix
	 * @return the Log-Likelihood of Transtion Matrix
	 */
	private double transitionMatrix(double[][] transitionMatrix) {

		double sum = 0;

		int j = 1;

		for (int i = 1; i < transitionMatrix.length; i++) {
			sum += Math.log(transitionMatrix[i - 1][j]);
			j++;
		}

		sum += Math.log(transitionMatrix[0][0]);

		return sum;
	}

	/**
	 * This method calculate the Log-Likelihood function to the Alpha Matrix in
	 * one dimension
	 * 
	 * @param alphaArray
	 * @return the Log-Likelihood of Alpha Array
	 */
	private double alphaArray(double[] alphaArray) {

		double sum = 0;

		for (int i = 0; i < alphaArray.length; i++) {
			sum += Math.log(alphaArray[i]);
		}

		return sum;
	}

	/**
	 * This method calculate the difference between speciation and duplication
	 * events
	 * 
	 * @param sigmaArray
	 * @return the difference between the events
	 */
	private double sigmaArray(double[] sigmaArray) {

		return sigmaArray[1] - sigmaArray[0];

	}

}

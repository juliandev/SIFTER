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

import java.util.Random;
import java.util.Vector;

/**
 * 
 * @author Eng. (C) Julian Camilo Castañeda Alonso, Msc. Carlos Andres Sierra and Msc. Tania Andrea Rodriguez Quiñones
 *
 */

public class GeneticOperator {

	/**
	 * Default constructor
	 */
	public GeneticOperator() {

	}

	/**
	 * This method implements Mutation strategy for Transtion Matrix of Markov
	 * Models
	 * 
	 * @param parent_1
	 * @param parent_2
	 * @return sons in an vector
	 */
	public Vector<double[][]> mutationTransitionMatrix(double[][] parent_1, double[][] parent_2) {

		Random rnd = new Random();

		// Vector of sons's chromosomes
		Vector<double[][]> sons = new Vector<>();

		double[][] firstSon = parent_1.clone();
		double[][] secondSon = parent_2.clone();

		// The first child is produced

		// Selected row
		int x = rnd.nextInt(parent_1.length);

		// Position of the row to change
		int y = rnd.nextInt(parent_1.length);

		// Point to change number from first father
		double random = rnd.nextDouble();

		// Verifies that the random number is less than the number selected in
		// the row
		while (random > firstSon[x][y]) {
			random = rnd.nextDouble();
		}

		// Subtract the random number from the number selected in the row
		firstSon[x][y] -= random;

		// Divide the random number by the number of states - 1
		double difference = random / (parent_1.length - 1);

		// Recalculates the odds of the other states by adding up the difference
		for (int j = 0; j < firstSon.length; j++) {

			if (j != y) {

				firstSon[x][j] += difference;

			}

		}

		// The second child is produced

		// Selected row
		x = rnd.nextInt(parent_2.length);

		// Position of the row to change
		y = rnd.nextInt(parent_2.length);

		// Point to change number from second father
		random = rnd.nextDouble();

		// Verifies that the random number is less than the number selected in
		// the row
		while (random > secondSon[x][y]) {
			random = rnd.nextDouble();
		}

		// Subtract the random number from the number selected in the row
		secondSon[x][y] -= random;

		// Divide the random number by the number of states - 1
		difference = random / (parent_2.length - 1);

		// Recalculates the odds of the other states by adding up the difference
		for (int j = 0; j < secondSon.length; j++) {

			if (j != y) {

				secondSon[x][j] += difference;

			}

		}

		sons.add(firstSon);
		sons.add(secondSon);

		return sons;
	}

	/**
	 * This method implements traditional Mutation strategy for Alpha Matrix
	 * 
	 * @param parent_1
	 * @param parent_2
	 * @return sons in an vector
	 */
	public Vector<double[][]> mutationAlpha(double[] parent_1, double[] parent_2) {

		Random rnd = new Random();

		// Vector of sons's chromosomes
		Vector<double[][]> sons = new Vector<>();

		double[][] childs = new double[parent_1.length][parent_1.length];

		childs[0] = parent_1.clone();
		childs[1] = parent_2.clone();

		// Points to exchange two positions
		int random_1 = rnd.nextInt(parent_1.length);
		int random_2 = rnd.nextInt(parent_2.length);
		double temp = 0;

		// Verifies that the two random numbers are different
		while (random_1 == random_2) {
			random_2 = rnd.nextInt(parent_2.length);
		}

		// Swaps the two selected positions in child 1
		temp = childs[0][random_1];
		childs[0][random_1] = childs[0][random_2];
		childs[0][random_2] = temp;

		// Swaps the two selected positions in child 2
		temp = childs[1][random_1];
		childs[1][random_1] = childs[1][random_2];
		childs[1][random_2] = temp;

		sons.add(childs);

		return sons;

	}

	/**
	 * This method implements traditional Mutation strategy for Sigma Matrix
	 * 
	 * @param parent_1
	 * @param parent_2
	 * @return sons in an vector
	 */
	public Vector<double[][]> mutationSigma(double[] parent_1, double[] parent_2) {

		Random rnd = new Random();

		// Vector of sons's chromosomes
		Vector<double[][]> sons = new Vector<>();

		double[][] childs = new double[parent_1.length][parent_1.length];

		childs[0] = parent_1.clone();
		childs[1] = parent_2.clone();

		// Position to change the value of first child
		int random = rnd.nextInt(2);

		// Generate a random number in the selected position
		childs[0][random] = rnd.nextDouble();

		// For the first child if the event of the first position is greater
		// than the value of the second position it exchanges them
		if (childs[0][0] > childs[0][1]) {
			double temp = childs[0][0];
			childs[0][0] = childs[0][1];
			childs[0][1] = temp;
		}

		// Position to change the value of second child
		random = rnd.nextInt(2);

		// Generate a random number in the selected position
		childs[1][random] = rnd.nextDouble();

		// For the first child if the event of the first position is greater
		// than the value of the second position it exchanges them
		if (childs[1][0] > childs[1][1]) {
			double temp = childs[1][0];
			childs[1][0] = childs[1][1];
			childs[1][1] = temp;
		}

		sons.add(childs);

		return sons;
	}

	/**
	 * This method implements Crossover strategy for Transtion Matrix of Markov
	 * Models
	 * 
	 * @param parent_1
	 * @param parent_2
	 * @return sons in an vector
	 */
	public Vector<double[][]> crossoverTransitionMatrix(double[][] parent_1, double[][] parent_2) {

		Random rnd = new Random();

		// Vector of sons's chromosomes
		Vector<double[][]> sons = new Vector<>();

		double[][] firstSon = parent_1.clone();
		double[][] secondSon = parent_2.clone();

		// First row selected to swap
		int point_1 = rnd.nextInt(parent_1.length);

		// Second row selected to swap
		int point_2 = rnd.nextInt(parent_1.length);

		// Exchange the selected row between the first and the second parent
		double[] temp = firstSon[point_1];
		firstSon[point_1] = secondSon[point_2];
		secondSon[point_2] = temp;

		sons.addElement(firstSon);
		sons.addElement(secondSon);

		return sons;

	}

	/**
	 * This method implements traditional Crossover strategy
	 * 
	 * @param parent_1
	 * @param parent_2
	 * @return sons in an vector
	 */
	public Vector<double[][]> crossoverAlpha(double[] parent_1, double[] parent_2) {

		Random rnd = new Random();

		// Vector of sons's chromosomes
		Vector<double[][]> sons = new Vector<>();

		double[][] childs = new double[parent_1.length][parent_1.length];

		childs[0] = parent_1.clone();
		childs[1] = parent_2.clone();

		// Segment of the matrix to cross
		int point = rnd.nextInt(parent_1.length);

		double temp = 0;

		// Cross segments between first and second parent
		for (int i = 0; i < point; i++) {
			temp = childs[0][i];
			childs[0][i] = childs[1][i];
			childs[1][i] = temp;
		}

		sons.add(childs);

		return sons;

	}

	/**
	 * This method implements traditional Crossover strategy
	 * 
	 * @param parent_1
	 * @param parent_2
	 * @return sons in an vector
	 */
	public Vector<double[][]> crossoverSigma(double[] parent_1, double[] parent_2) {

		Random rnd = new Random();

		// Vector of sons's chromosomes
		Vector<double[][]> sons = new Vector<>();

		double[][] childs = new double[parent_1.length][parent_1.length];

		childs[0] = parent_1.clone();
		childs[1] = parent_2.clone();

		// Segment of the matrix to cross
		int point = rnd.nextInt(parent_1.length);

		// Cross segments between first and second parent
		double temp = childs[0][point];
		childs[0][point] = childs[1][point];
		childs[1][point] = temp;

		// For the first child if the event of the first position is greater
		// than the value of the second position it exchanges them
		if (childs[0][0] > childs[0][1]) {

			temp = childs[0][0];
			childs[0][0] = childs[0][1];
			childs[0][1] = temp;

		}

		// For the second child if the event of the first position is greater
		// than the value of the second position it exchanges them
		if (childs[1][0] > childs[1][1]) {

			temp = childs[1][0];
			childs[1][0] = childs[1][1];
			childs[1][1] = temp;

		}

		sons.add(childs);

		return sons;

	}

}

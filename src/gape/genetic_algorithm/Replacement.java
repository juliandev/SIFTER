/*
# GAPE - Genetic Algorithm for Parameter Estimation of SIFTER tool
#
# Created by Eng. (C) Julian Camilo Castañeda Alonso and Msc. Tania Andrea Rodriguez Quiñones on August 2017.
# Copyright (c) 2017. Eng. (C) Julian Camilo Castañeda Alonso and Msc. Tania Andrea Rodriguez Quiñones. Universidad Antonio Narino. All rights reserved.
#
# GAPE is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License v3.0 found in the LICENSE file in the root directory of this project.
*/

package gape.genetic_algorithm;

import java.util.Random;

/**
 * 
 * @author Eng. (C) Julian Camilo Castañeda Alonso - Msc. Tania Andrea Rodriguez Quiñones
 *
 */

public class Replacement {

	Random rnd = new Random();

	/**
	 * Default constructor
	 */
	public Replacement() {

	}

	/**
	 * This method implements traditional Steady-State strategy
	 * 
	 * @param father_1
	 * @param father_2
	 * @param son_1
	 * @param son_2
	 * @return invididuals an a vector
	 */
	public Individual[] steadyState(Individual father_1, Individual father_2, Individual son_1, Individual son_2) {

		Individual[] winners = new Individual[2];

		if (father_1.getFitness() < father_2.getFitness()) {

			if (son_1.getFitness() < son_2.getFitness()) {

				winners[0] = this.roulette(father_2, son_2); // Best father -
																// Best son
				winners[1] = this.roulette(father_1, son_1); // Worst father -
																// Worst son

			} else {

				winners[0] = this.roulette(father_2, son_1); // Best father -
																// Best son
				winners[1] = this.roulette(father_1, son_2); // Worst father -
																// Worst son

			}
		} else {

			if (son_1.getFitness() < son_2.getFitness()) {

				winners[0] = this.roulette(father_1, son_2); // Best father -
																// Best son
				winners[1] = this.roulette(father_2, son_1); // Worst father -
																// Worst son

			} else {

				winners[0] = this.roulette(father_1, son_1); // Best father -
																// Best son
				winners[1] = this.roulette(father_2, son_2); // Worst father -
																// Worst son

			}
		}

		return winners;

	}

	/**
	 * This method implements traditional Roullete strategy
	 * 
	 * @param player_1
	 * @param player_2
	 * @return an individual
	 */
	private Individual roulette(Individual player_1, Individual player_2) {

		// Fitness of the two individuals
		double totalFitness = player_1.getFitness() + player_2.getFitness();

		// Normalizes the fitness of the first individual
		double point = 1 - (player_1.getFitness() / totalFitness);

		// Generate a random number to select one of the two individuals
		return rnd.nextDouble() < point ? player_1 : player_2;
	
	}

}

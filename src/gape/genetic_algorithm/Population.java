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

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Random;
import java.util.Vector;
import gape.markov_model.MarkovModel;
import gape.output.PrintFiles;

/**
 * 
 * @author Eng. (C) Julian Camilo Castañeda Alonso - Msc. Tania Andrea Rodriguez Quiñones
 *
 */

public class Population {

	// Attributes
	private int numStates = 0;
	private int populationSize = 0;
	private int iterations = 0;
	private double rateMutation = 0;
	private String idFamily = null;
	private Random rnd = null;
	private Vector<Individual> population = null;
	private int[] rowNames = null;

	FitnessFunction fitnessFunction = null;
	GeneticOperator geneticOperator = null;
	Replacement replacement = null;
	BufferedWriter bw = null;

	/**
	 * Constructor of the class
	 * 
	 * @param populationSize
	 * @param iterations
	 * @param rateMutation
	 * @param rowNames
	 * @param idFamily
	 */
	public Population(int populationSize, int iterations, double rateMutation, int[] rowNames, String idFamily) {

		bw = new BufferedWriter(new OutputStreamWriter(System.out));

		try {

			bw.write("\nInitializing Genetic Algorithm...\n");
			bw.flush();

		} catch (IOException e) {

			e.printStackTrace();

		}

		this.numStates = rowNames.length;
		this.populationSize = populationSize;
		this.iterations = iterations;
		this.rateMutation = rateMutation;
		this.fitnessFunction = new FitnessFunction();
		this.geneticOperator = new GeneticOperator();
		this.replacement = new Replacement();
		this.replacement = new Replacement();
		this.rnd = new Random();
		this.population = new Vector<>();
		this.rowNames = rowNames;
		this.idFamily = idFamily;

		try {

			bw.write("Building initial population...\n");
			bw.flush();

			this.initialPopulation();

			bw.write("Start Genetic Algorithm...\n");
			bw.flush();

			this.generations();

		} catch (IOException e) {

			e.printStackTrace();

		}

	}

	/**
	 * This method initializes the population of algorithm
	 */
	private void initialPopulation() {

		double fitness = 0;

		for (int i = 0; i < populationSize; i++) {
			double[][] individualTransitionMatrix = MarkovModel.createTransitionMatrix(numStates);
			double[] individualAlpha = MarkovModel.createAlpha(numStates);
			double[] individualPhi = MarkovModel.createSigma();
			fitness = fitnessFunction.fitness(individualTransitionMatrix, individualAlpha, individualPhi);
			this.population.add(new Individual(individualTransitionMatrix, individualAlpha, individualPhi, fitness));
		}

	}

	/**
	 * This method contains the different generations of algorithm
	 */
	private void generations() {

		Individual father_1 = null;
		Individual father_2 = null;
		Vector<double[][]> sonsGenotipes = null;
		Individual[] sons = new Individual[2];
		Individual[] newGeneration = new Individual[2];
		Vector<Individual> nextGeneration = null;
		double tempFitness;

		// Iterations or generations
		for (int i = 0; i < this.iterations; i++) {

			System.gc();

			nextGeneration = new Vector<Individual>();

			// Repeat operation to build new generation
			for (int j = 0; j < (populationSize / 2); j++) {

				father_1 = tournament();  
				father_2 = tournament();

				// Apply genetic operator to obtain sons
				if (rnd.nextDouble() < rateMutation) {

					int random = this.rnd.nextInt(3) + 1;

					if (random == 1) {

						sonsGenotipes = this.geneticOperator.mutationTransitionMatrix(father_1.getTransitionMatrix(),
								father_2.getTransitionMatrix());

						tempFitness = this.fitnessFunction.fitness(sonsGenotipes.get(0), father_1.getAlpha(),
								father_1.getSigma());
						sons[0] = new Individual(sonsGenotipes.get(0), father_1.getAlpha(), father_1.getSigma(),
								tempFitness);

						tempFitness = this.fitnessFunction.fitness(sonsGenotipes.get(1), father_2.getAlpha(),
								father_2.getSigma());
						sons[1] = new Individual(sonsGenotipes.get(1), father_2.getAlpha(), father_2.getSigma(),
								tempFitness);

					} else if (random == 2) {

						sonsGenotipes = this.geneticOperator.mutationAlpha(father_1.getAlpha(), father_2.getAlpha());

						tempFitness = this.fitnessFunction.fitness(father_1.getTransitionMatrix(),
								sonsGenotipes.get(0)[0], father_1.getSigma());
						sons[0] = new Individual(father_1.getTransitionMatrix(), sonsGenotipes.get(0)[0],
								father_1.getSigma(), tempFitness);

						tempFitness = this.fitnessFunction.fitness(father_2.getTransitionMatrix(),
								sonsGenotipes.get(0)[1], father_2.getSigma());
						sons[1] = new Individual(father_1.getTransitionMatrix(), sonsGenotipes.get(0)[1],
								father_1.getSigma(), tempFitness);

					} else if (random == 3) {

						sonsGenotipes = this.geneticOperator.mutationSigma(father_1.getSigma(), father_2.getSigma());

						tempFitness = this.fitnessFunction.fitness(father_1.getTransitionMatrix(), father_1.getAlpha(),
								sonsGenotipes.get(0)[0]);
						sons[0] = new Individual(father_1.getTransitionMatrix(), father_1.getAlpha(),
								sonsGenotipes.get(0)[0], tempFitness);

						tempFitness = this.fitnessFunction.fitness(father_2.getTransitionMatrix(), father_2.getAlpha(),
								sonsGenotipes.get(0)[1]);
						sons[1] = new Individual(father_1.getTransitionMatrix(), father_1.getAlpha(),
								sonsGenotipes.get(0)[1], tempFitness);

					}

				} else {

					int random = this.rnd.nextInt(3) + 1;

					if (random == 1) {

						sonsGenotipes = this.geneticOperator.crossoverTransitionMatrix(father_1.getTransitionMatrix(),
								father_2.getTransitionMatrix());

						tempFitness = this.fitnessFunction.fitness(sonsGenotipes.get(0), father_1.getAlpha(),
								father_1.getSigma());
						sons[0] = new Individual(sonsGenotipes.get(0), father_1.getAlpha(), father_1.getSigma(),
								tempFitness);

						tempFitness = this.fitnessFunction.fitness(sonsGenotipes.get(1), father_2.getAlpha(),
								father_2.getSigma());
						sons[1] = new Individual(sonsGenotipes.get(1), father_2.getAlpha(), father_2.getSigma(),
								tempFitness);

					} else if (random == 2) {

						sonsGenotipes = this.geneticOperator.crossoverAlpha(father_1.getAlpha(), father_2.getAlpha());

						tempFitness = this.fitnessFunction.fitness(father_1.getTransitionMatrix(),
								sonsGenotipes.get(0)[0], father_1.getSigma());
						sons[0] = new Individual(father_1.getTransitionMatrix(), sonsGenotipes.get(0)[0],
								father_1.getSigma(), tempFitness);

						tempFitness = this.fitnessFunction.fitness(father_2.getTransitionMatrix(),
								sonsGenotipes.get(0)[1], father_2.getSigma());
						sons[1] = new Individual(father_1.getTransitionMatrix(), sonsGenotipes.get(0)[1],
								father_1.getSigma(), tempFitness);

					} else if (random == 3) {

						sonsGenotipes = this.geneticOperator.crossoverSigma(father_1.getSigma(), father_2.getSigma());

						tempFitness = this.fitnessFunction.fitness(father_1.getTransitionMatrix(), father_1.getAlpha(),
								sonsGenotipes.get(0)[0]);
						sons[0] = new Individual(father_1.getTransitionMatrix(), father_1.getAlpha(),
								sonsGenotipes.get(0)[0], tempFitness);

						tempFitness = this.fitnessFunction.fitness(father_2.getTransitionMatrix(), father_2.getAlpha(),
								sonsGenotipes.get(0)[1]);
						sons[1] = new Individual(father_1.getTransitionMatrix(), father_1.getAlpha(),
								sonsGenotipes.get(0)[1], tempFitness);

					}

				}

				/*newGeneration = this.replacement.steadyState(father_1, father_2, sons[0], sons[1]);

				nextGeneration.add(newGeneration[0]);
				nextGeneration.add(newGeneration[1]);*/
				nextGeneration.add(sons[0]);
				nextGeneration.add(sons[1]);
			}

			this.population = nextGeneration;

		}

	}

	public Individual uniform() {

		Individual father = null;

		int candidate = rnd.nextInt(this.populationSize);

		father = this.population.get(candidate);

		return father;
	}

	/**
	 * This method implements Tournament strategy based in traditional Roulette
	 * strategy
	 * 
	 * @return an individual
	 */
	private Individual tournament() {

		Individual father = null;

		int candidate_1, candidate_2, candidate_3, candidate_4;

		candidate_1 = rnd.nextInt(populationSize);
		candidate_2 = rnd.nextInt(populationSize);
		candidate_3 = rnd.nextInt(populationSize);
		candidate_4 = rnd.nextInt(populationSize);

		father = tournament(this.population.get(candidate_1), this.population.get(candidate_2),
				this.population.get(candidate_3), this.population.get(candidate_4));

		return father;

	}

	/**
	 * @param player_1
	 * @param player_2
	 * @return an individual
	 */
	private Individual tournament(Individual player_1, Individual player_2, Individual player_3, Individual player_4) {
		System.out.println(player_1.getFitness() + "\t" + player_2.getFitness() + "\t" + player_3.getFitness() + "\t" + player_4.getFitness());
		Individual winner = roulette(roulette(player_1, player_2), roulette(player_3, player_4));
		System.out.println("Ganador: " + winner.getFitness());
		return winner;
	}

	/**
	 * This method implements traditional Roullete strategy
	 * 
	 * @param player_1
	 * @param player_2
	 * @return an individual
	 */
	private Individual roulette(Individual player_1, Individual player_2) {

		double totalFitness = player_1.getFitness() + player_2.getFitness();

		double point = 1 - (player_1.getFitness() / totalFitness);
		
		return rnd.nextDouble() < point ? player_1 : player_2;
	}

	// TODO
	/**
	 * This method selects the best individual from the population
	 * 
	 * @return the best individual of population
	 */
	public String getBest() {

		int index = 0;

		double fitness = this.population.get(0).getFitness();

		for (int i = 0; i < this.populationSize; i++) {

			if (this.population.get(i).getFitness() > fitness) {

				index = i;

				fitness = this.population.get(i).getFitness();

			}

		}

		String best = "Transition Matrix: \n" + this.population.get(index).getStringTransitionMatrix();
		best += "\nAlpha:" + this.population.get(index).getStringAlphaMatrix();
		best += "\n\nPhi:" + this.population.get(index).getStringSigmaMatrix();
		best += "\n\nFitness:" + this.population.get(index).getFitness();

		PrintFiles files = new PrintFiles();
		files.printTransitionMatrix("output/infer-" + this.idFamily + ".fx",
				this.population.get(index).getTransitionMatrix(), this.rowNames);
		files.printAlpha("output/alpha-" + this.idFamily + ".fx", this.population.get(index).getAlpha());
		files.printSigma("output/scale-" + this.idFamily + ".fx", this.population.get(index).getSigma());
		files.printBestIndividual("output/BestIndividual.txt", best);

		return best;

	}

}

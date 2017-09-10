/**
 * Copyright 2003-2005 Barbara Engelhardt (bee@cs.berkeley.edu)
 *
 * This module handles probabilities in each node of the DAG.
 * Used in ProbabilisticDAG.
 *
 * $Id: ProbabilisticDAGNode.java,v 1.2 2005/05/27 06:56:01 sprite Exp $
 * @author Barbara Engelhardt.
 */

package sifter.components;

import java.util.*;
//import sifter_components.*;

public class ProbabilisticDAGNode extends Object {

	private double prior;
	private double evidence;
	private double likelihood;
	private double additionalPosterior;
	private Vector<Integer> descendantLeaves;
	private double subsetProbability;

	public ProbabilisticDAGNode(double pri) {
		prior = pri;
		//System.out.println("Setting initial prior to: "+pri);
		evidence = 0;
		likelihood = 0.0;
		additionalPosterior = 0.0;
		descendantLeaves = null;
		subsetProbability = -1.0;
	}

	public ProbabilisticDAGNode() {
		prior = 0.0;
		evidence = 0;
		likelihood = 0.0;
		additionalPosterior = 0.0;
		descendantLeaves = null;
		subsetProbability = -1.0;
	}

	public Object clone() {
		return this.clone();
	}

	public void setPrior(double p) {
		prior = p;
		//System.out.println("Setting prior to: "+p);
	}

	public double getInitialPrior(int l, double rValue) {
		if (prior > 0.0)
			return prior;

		double pri = 0.0;

		for (int i = 0; i < l; i++) {
			pri = pri + (GenericMathFunctions.choose(l, i + 1) / Math.pow(rValue, i + 1));
		}

		setPrior(pri);
		//System.out.println("Calculating prior: "+pri+", for l = "+l);
		return pri;
	}
	public double getPrior() {
		return prior;
	}

	// Equivalent function to evidence,
	// except non-direct evidence included.
	public void addLikelihood(double l) {
		likelihood = 1 - ((1 - likelihood) * (1 - l));
	}

	public void setPosterior(double l) {
		if (prior > 0)
			likelihood = l / prior;
		else {
			prior = 1;
			likelihood = l;
		}

		additionalPosterior = 0;
	}

	public double getLikelihood() {
		return likelihood;
	}

	public void setLikelihood(double l) {
		likelihood = l;
	}

	public void addPosterior(double p) {
		additionalPosterior += p;
	}

	public double getPosterior() {
		return ((likelihood*prior) + additionalPosterior);
	}

	// If evidence is 0, as it is to start,
	// then this is just e for the first
	// piece of evidence added.
	public void addEvidence(double e) {
		evidence = 1 - ((1 - evidence) * (1 - e));
		addLikelihood(e);
	}

	public double getEvidence() {
		return evidence;
	}

	public void resetProbs() {
		evidence = 0.0;
		likelihood = 0.0;
		additionalPosterior = 0.0;
	}

	public void setLeaves(Vector<Integer> l, double rValue) {
		if (l.size() == 0) {
			getInitialPrior(1, rValue);
			return;
		}
		else
			getInitialPrior(l.size(), rValue);

		//System.out.println("Setting leaves, prior: "+nprob.getPrior()+", size: "+l.size());
		addDescendantLeaves(l);
	}

	public void addDescendantLeaf(Integer d) {
		if (descendantLeaves == null)
			descendantLeaves = new Vector<Integer>();

		descendantLeaves.add(d);
	}

	public void addDescendantLeaves(Vector<Integer> l) {
		descendantLeaves = l;

	}

	public Vector<Integer> getDescendantLeaves() {
		return descendantLeaves;
	}

	public int getNumDescendantLeaves() {
		if (descendantLeaves == null)
			return 0;
		else
			return descendantLeaves.size();
	}

	public double getSubsetProbability() {
		return subsetProbability;
	}

	public void setSubsetProbability(double sp) {
		subsetProbability = sp;
	}

}

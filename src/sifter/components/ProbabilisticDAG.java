/**
 * Copyright 2003-2005 Barbara Engelhardt (bee@cs.berkeley.edu)
 *
 * DAG data structure with probability ideas in it.
 * It really should be generic
 *
 * @author Barbara Engelhardt
 */

package sifter.components;

import java.util.*;


public class ProbabilisticDAG extends Object {
	// indexed by id number.
	public Hashtable<Integer, Node> nodeList;
	public Hashtable<Integer, ProbabilisticDAGNode> probList;
	public double rValue;
	public Vector<Integer> leaves;

	public ProbabilisticDAG() {
		nodeList = new Hashtable<Integer, Node>();
		probList = new Hashtable<Integer, ProbabilisticDAGNode>();
		rValue = 0.0;
		leaves = null;
	}

	public boolean contains(int id) {
		return(nodeList.containsKey(new Integer(id)));
	}


	public void setRValue(double r) {
		rValue = r;
	}

	public void addEvidence(int id, double evidence) {
		Integer node = new Integer(id);

		if (!probList.containsKey(node))
			return;

		ProbabilisticDAGNode n = probList.get(node);
		n.addEvidence(evidence);
	}

	public void addNode(String name, int id) {
		Node n = new Node(name, id);
		ProbabilisticDAGNode np = new ProbabilisticDAGNode();
		Integer intID = new Integer(id);
		nodeList.put(intID, n);
		probList.put(intID, np);
	}

	public void addParent(int childId, int parentId) {
		Integer child = new Integer(childId);
		Integer parent = new Integer(parentId);

		if (!nodeList.containsKey(child) ||
		    !nodeList.containsKey(parent))
			return;

		Node ch = nodeList.get(child);
		Node p = nodeList.get(parent);
		ch.addInEdge(p);
		p.addOutEdge(ch);
	}

	// Remove node: removes a node from the DAG
	// and connects the children with the parents
	// appropriately:
	// for all parents and all children of node,
	// remove child/parent (respectively)
	// then addParent with the two.
	public void removeNode(int rNode) {
		Integer node = new Integer(rNode);

		if (!nodeList.containsKey(node))
			return;

		Node n = nodeList.get(node);
		Vector<Node> parents = n.getParents();
		Vector<Node> children = n.getChildren();

		if (parents != null) {
			for (int i = 0; i < parents.size(); i++) {
				Node p = parents.elementAt(i);

				if (p != null) {
					p.removeChild(n);
				}
			}
		}

		if (children != null) {
			for (int i = 0; i < children.size(); i++) {
				Node c = children.elementAt(i);

				if (c != null) {
					c.removeParent(n);
				}
			}
		}

		if (parents != null && children != null) {
			for (int i = 0; i < parents.size(); i++) {
				Node p = nodeList.get(parents.elementAt(i));

				if (p != null) {
					for (int j = 0; j < children.size(); j++) {
						Node c = nodeList.get(children.elementAt(j));

						if (c != null) {
							addParent(c.getId(), p.getId());
						}
					}
				}
			}
		}

		nodeList.remove(node);

		if (!probList.containsKey(node))
			return;

		probList.remove(node);
		leaves = null;
	}

	// for all parents remove child
	// then remove node & descendants from nodeList
	public void removeNodeNDescendants(int rNode) {
		Integer node = new Integer(rNode);

		if (!nodeList.containsKey(node))
			return;

		Node n = nodeList.get(node);
		Vector<Node> parents = n.getParents();

		if (parents != null) {
			for (int i = 0; i < parents.size(); i++) {
				Node p = nodeList.get(parents.elementAt(i));

				if (p != null) {
					p.removeChild(n);
				}
			}
		}

		Vector<Integer> dg = getDescendantGraph(rNode, false);

		if (dg != null) {
			for (int i = 0; i < dg.size(); i++) {
				nodeList.remove(dg.elementAt(i));
			}
		}

		nodeList.remove(node);
	}

	// could stop if already in master list to save time.
	// and automatically remove multiples.
	public Vector<Integer> getInducedSubgraph(int nodeId, Vector<Integer> v) {
		Integer node = new Integer(nodeId);

		if (!nodeList.containsKey(node))
			return v;

		if (v.contains(node))
			return v;

		Node n = nodeList.get(node);
		Vector<Node> parents = n.getParents();
		v.add(new Integer(nodeId));

		if (parents == null)
			return v;

		for (int i = 0; i < parents.size(); i++) {
			v = getInducedSubgraph(parents.elementAt(i), v);
		}

		//v = removeMultiples(v);
		return v;
	}

	public Vector<Integer> getInducedSubgraph(Node node, Vector<Integer> v) {
		Vector<Node> parents = node.getParents();
		Integer n = new Integer(node.getId());

		if (v.contains(n))
			return v;

		v.add(n);

		if (parents == null)
			return v;

		for (int i = 0; i < parents.size(); i++) {
			v = getInducedSubgraph(parents.elementAt(i), v);
		}

		return v;
	}

	// Runs in breadth-first search-order over parents
	// Vector returned includes the node which seeded the function.
	// Currently prunes duplicates!
	public Vector<Integer> getInducedSubgraph(int fn) {
		Vector<Integer> isg = new Vector<Integer>();
		Integer node = new Integer(fn);

		if (nodeList.containsKey(node)) {
			isg.add(node);
		}
		else
			return isg;

		int index = 0;

		while (index < isg.size()) {
			Integer dNode = (Integer)isg.elementAt(index);

			// Check for duplicates; remove this one if found
			if (isg.indexOf(dNode) != isg.lastIndexOf(dNode)) {
				isg.removeElementAt(index);
				index--;
			}
			else
				if (nodeList.containsKey(dNode)) {
					Vector<Node> parents = nodeList.get(dNode).getParents();

					if (parents != null) {
						for (int i = 0; i < parents.size(); i++) {
							isg.add(new
							        Integer(parents.elementAt(i).getId()));
						}
					}
				}

			index++;
		}

		return isg;
	}

	// Runs in breadth-first search-order over children
	// Currently prunes duplicates!
	// Does not include self.
	public Vector<Integer> getDescendantGraph(int fn, boolean setLeaves) {
		Vector<Integer> dg = new Vector<Integer>();
		Node n = null;
		Integer node = new Integer(fn);

		if (nodeList.containsKey(node)) {
			n = nodeList.get(node);

			if (n.getChildren() != null) {
				Vector<Node> c = n.getChildren();

				for (int i = 0; i < c.size(); i++) {
					dg.add(new Integer(c.elementAt(i).getId()));
				}
			}
			else {
				return dg;
			}
		}
		else
			return dg;

		int index = 0;

		while (index < dg.size()) {
			Integer dNode = dg.elementAt(index);

			// Check for duplicates; remove this one if found
			if (dg.indexOf(dNode) != dg.lastIndexOf(dNode)) {
				dg.removeElementAt(index);
				index--;
			}
			else
				if (nodeList.containsKey(dNode)) {
					Vector<Node> children = nodeList.get(dNode).getChildren();

					if (children != null) {
						for (int i = 0; i < children.size(); i++) {
							dg.add(dg.size(),
							       new Integer(children.elementAt(i)
							                   .getId()));
						}
					}
				}

			index++;
		}

		return dg;
	}

	public void propagateEvidenceDownward(int fn, double e) {
		Vector<Integer> dg = getDescendantGraph(fn, true);
		double oldPrior = getPrior(fn);

		for (int i = 0; i < dg.size(); i++) {
			Integer prob = dg.elementAt(i);
			double newPrior = getPrior(prob.intValue());

			if (probList.containsKey(prob)) {
				probList.get(prob)
				.addLikelihood(e*newPrior / oldPrior);
				//System.out.println("Adding "+e*newPrior/oldPrior);
			}
		}
	}

	// Given a node which is NOT a leaf, propagates
	// all the leaf probabilities upwards.
	// Should be combined with propagating the
	// internal nodes upwards too...
	public void propagateLeavesUpwards(Integer node) {
		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);
			Vector<Integer> dl = getDescendantLeaves(node);

			if (dl == null || dl.size() < 1)
				return;

			double totalLikelihood = 0.0;
			double totalPrior = 0.0;
			double currentPrior = np.getPrior();

			for (int i = 0; i < dl.size(); i++) {
				if (probList.containsKey(dl.elementAt(i))) {
					ProbabilisticDAGNode ch =
					  probList.get(dl.elementAt(i));
					double singlePrior = ch.getPrior();
					double chLikelihood = ch.getLikelihood();
					//totalLikelihood += chLikelihood;
					totalLikelihood += (chLikelihood * singlePrior);
					totalPrior += singlePrior;
				}
			}

			if (totalPrior > 0)
				np.setPosterior(totalLikelihood / totalPrior);
			else
				np.setPosterior(currentPrior * totalLikelihood);
		}
	}

	// Accessor Functions for external classes
	// for the probabilities in a node
	public double getPrior(int pfName) {
		Integer node = new Integer(pfName);

		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);

			if (np.getNumDescendantLeaves() == 0)
				getDescendantLeaves(node);

			return np.getPrior();
		}

		return 0;
	}

	public double getPosterior(int pfName) {
		Integer node = new Integer(pfName);

		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);
			return np.getPosterior();
		}

		return 0;
	}

	public double getLikelihood(int pfName) {
		Integer node = new Integer(pfName);

		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);
			return np.getLikelihood();
		}

		return 0;

	}

	public void setLikelihood(int pfName, double l) {
		Integer node = new Integer(pfName);

		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);
			np.setLikelihood(l);
		}
	}

	public double getEvidence(int pfName) {
		Integer node = new Integer(pfName);

		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);
			return np.getEvidence();
		}

		return 0;
	}

	public String getName(int pfName) {
		Integer node = new Integer(pfName);

		if (nodeList.containsKey(node)) {
			Node n = nodeList.get(node);
			return n.getName();
		}

		return null;
	}

	public int getNumberOfNodes() {
		return nodeList.size();
	}

	public int getNumDescendantLeaves(int pfName) {
		Integer node = new Integer(pfName);

		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);
			return np.getNumDescendantLeaves();
		}

		return 0;
	}
	// Always updates the node with the
	// descendant leaves
	public Vector<Integer> getDescendantLeaves(Integer node) {
		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);
			Vector<Integer> dl = np.getDescendantLeaves();

			while (dl == null) {
				Vector<Integer> d = getDescendantGraph(node.intValue(), false);
				dl = getAllLeaves(d);
				np.setLeaves(dl, rValue);
			}

			return dl;
		}

		return null;
	}

	public void propagateInternalUpward(int pfName) {
		ProbabilisticDAGNode np = null;
		Node n = null;
		Integer node = new Integer(pfName);

		if (probList.containsKey(node) && nodeList.containsKey(node)) {
			np = probList.get(node);
			n = nodeList.get(node);

			// don't propagate for leaves -- already done above.
			if (n.getChildren() == null || n.getChildren().size() < 1)
				return;

			Vector<Integer> is = getDescendantGraph(node.intValue(), true);

			if (is == null)
				return;

			if (!containsInternalEvidence(is))
				return;

			double currentLikelihood = np.getLikelihood();

			for (int i = 0; i < is.size(); i++) {
				if (probList.containsKey(is.elementAt(i))) {
					ProbabilisticDAGNode pnp =
					  probList.get(is.elementAt(i));

					if (pnp.getNumDescendantLeaves() == 0)
						getDescendantLeaves(is.elementAt(i));

					// Does this line make any sense?
					np.addPosterior(currentLikelihood * pnp.getPrior());
				}
			}
		}
	}

	// Returns true if some internal node in this set
	// has evidence associated with it.
	private boolean containsInternalEvidence(Vector<Integer> is) {
		for (int i = 0; i < is.size(); i++) {
			Integer node = is.elementAt(i);

			if (getEvidence(node.intValue()) > 0
			    && getNumChildren(node.intValue()) != 0)
				return true;
		}

		return false;
	}

	// This function is never used?!?
	public double getIntroducedPrior(Node n, ProbabilisticDAGNode np) {
		double introducedL = 0;
		Vector<Node> ch = n.getChildren();
		Vector<Integer> leaves = np.getDescendantLeaves();
		int leafSize = 0;

		if (leaves != null)
			leafSize = leaves.size();

		if (ch == null)
			return (1.0 / rValue);

		Vector<Vector<Integer>> leavesLeft = new Vector<Vector<Integer>>();

		for (int i = 0; i < ch.size(); i++) {
			ProbabilisticDAGNode np2 = probList.get(ch.elementAt(i));
			Vector<Integer> l2 = np2.getDescendantLeaves();
			int l2Size = 0;

			if (l2 != null)
				l2Size = l2.size();

			if (l2Size == leafSize)
				return 0.0;

			leavesLeft.add(i, l2);
		}

		// For each child, if this child shares any
		// i leaves in common with any other child (left)
		// then don't add that set's prior.
		Vector<Integer> poly = getPolynomialVector(leafSize);

		for (int i = 0; i < ch.size() ; i++) {
			EliminateNotInCommon(poly, leavesLeft.elementAt(i),
			                     (Vector<Vector<Integer>>)leavesLeft.subList(i,
			                         leavesLeft.size()));
			// might be size-1
		}

		return introducedL;
	}

	public void EliminateNotInCommon(Vector<Integer> poly, Vector<Integer> leafSet,
	                                 Vector<Vector<Integer>> setLeafSet) {
		// First, make a list of items this leaf set has
		// in common with the remaining leaf sets
		Vector<Vector<Integer>> commonSets = new Vector<Vector<Integer>>();

		for (int i = 0; i < setLeafSet.size(); i++) {
			Vector<Integer> v = getSingleCommonSubset(leafSet,
			                    setLeafSet.elementAt(i));

			if (v != null)
				commonSets.add(v);
		}

		// Next, go through common sets. For each
		// common set, check to see whether it has
		// any elements in common with the priors
		// that have already been eliminated. Otherwise
		// use them to eliminate priors.
		for (int i = 0; i < commonSets.size(); i++) {
			Vector<Integer> cset = commonSets.elementAt(i);
			Vector<Integer> pThis = getPolynomialVector(cset.size());

			for (int j = 0; j < i; j++) {
				cset = eliminateFromPolySet(cset, pThis,
				                            commonSets.elementAt(j));
			}

			poly = eliminateFromPolySet(poly, cset);
		}
	}

	// this function takes two vectors, the first is the
	// polynomial coefficients, the second is the list of
	// leaf nodes that have already been accounted for
	// up to their numbers. it returns the polynomial
	// vector without the smaller sets made into a polynomial
	// Note: does not deal with x^0 coefficients, or sluff vars.
	public Vector<Integer> eliminateFromPolySet(Vector<Integer> poly, Vector<Integer> cset) {
		Vector<Integer> cpoly = getPolynomialVector(cset.size());

		for (int i = 0; i < cpoly.size(); i++) {
			Integer coe = poly.elementAt(i + 1);
			poly.removeElementAt(i + 1);
			poly.add(i + 1, new Integer(coe.intValue()
			                            - cpoly.elementAt(i + 1)
			                            .intValue()));
		}

		return poly;
	}

	// This function takes three vectors. The first is the
	// set common leaves, the second is the polynomial for
	// this common set and the third is one of the common
	// sets that have already been eliminated.
	// Note: not quite sure this is all right.
	public Vector<Integer> eliminateFromPolySet(Vector<Integer> cset, Vector<Integer> poly, Vector<Integer> nextSet) {
		Vector<Integer> intersect = GenericMathFunctions.intersection(cset, nextSet);

		if (intersect != null)
			cset.removeAll(intersect);

		return cset;
	}

	// returns a vector with elements at i corresponding to
	// integer parameters for x^i, including sluff variables.
	public Vector<Integer> getPolynomialVector(int level) {
		Vector<Integer> poly = new Vector<Integer>();
		poly.add(new Integer(0));

		// to deal with sluff variable
		if (level > 1)
			poly.add(new Integer(level + 1));
		else
			poly.add(new Integer(level));

		for (int i = 1; i < level; i++) {
			poly.add(new Integer((int)GenericMathFunctions.choose(level, i + 1)));
		}

		return poly;
	}

	//Returns an array with c elements which are the polynomial
	// coefficients for the polynomial describing a single
	// leaf's set containment possibility.
	public int[] getLeafPolynomialVector(int leaves) {
		int[] leafv = new int[leaves];

		for (int i = 0; i < leaves; i++) {
			leafv[i] = (int)GenericMathFunctions.choose(leaves - 1, i);
		}

		return leafv;
	}

	// will return null when there are no elements in common
	// between the two vectors
	public Vector<Integer> getSingleCommonSubset(Vector<Integer> l1, Vector<Integer> l2) {
		Vector<Integer> common = null;

		for (int i = 0; i < l1.size(); i++) {
			if (l2.contains(l1.elementAt(i))) {
				if (common == null)
					common = new Vector<Integer>();

				common.add(l1.elementAt(i));
			}
		}

		return common;
	}


	public int getTreeDistance(int leafA, int leafB) {
		int dist = 0;
		int maxLevel = 0;
		Integer nodeA = new Integer(leafA);
		Integer nodeB = new Integer(leafB);

		if (nodeList.containsKey(nodeA) && nodeList.containsKey(nodeB)) {
			Vector<Integer> intersect =
			  GenericMathFunctions.intersection(getInducedSubgraph(leafA),
			                                    getInducedSubgraph(leafB));

			for (int i = 0; i < intersect.size(); i++) {
				Integer n = (Integer)intersect.elementAt(i);
				int l = getLevel(n.intValue());

				if (l > maxLevel)
					maxLevel = l;
			}

			int lA = getLevel(leafA);
			int lB = getLevel(leafB);
			//System.out.println("A level: "+lA+", B level: "
			//		       +lB+", intersect level "+maxLevel);
			dist = (Math.max(lA - maxLevel, 1) + Math.max(lB - maxLevel, 1));
		}

		return dist;
	}

	/*
	  private Vector removeMultiples(Vector v)
	  {
	Vector vu = new Vector();
	for(int i = 0; i < v.size(); i++) {
	    Object o = v.elementAt(i);
	    vu.addElement(o);
	    while(v.contains(o)) {
		v.removeElement(o);
	    }
	}
	return vu;
	}*/

	public void addAdditional(int nodeId, Object o) {
		Integer node = new Integer(nodeId);

		if (!nodeList.containsKey(node))
			return;

		Node n = nodeList.get(node);
		n.addAdditional(o);
	}

	public Object getAdditional(int nodeId) {
		Integer node = new Integer(nodeId);

		if (!nodeList.containsKey(node))
			return null;

		Node n = nodeList.get(node);
		return n.getAdditional();
	}

	public void addLevel(int nodeId, int l) {
		Integer node = new Integer(nodeId);

		if (!nodeList.containsKey(node))
			return;

		Node n = nodeList.get(node);
		n.addLevel(l);
	}

	public int getLevel(int nodeId) {
		Integer node = new Integer(nodeId);

		if (!nodeList.containsKey(node))
			return -1;

		Node n = nodeList.get(node);
		return n.getLevel();
	}

	public int getNumParents(int nodeId) {
		Integer node = new Integer(nodeId);

		if (!nodeList.containsKey(node))
			return -1;

		Node n = nodeList.get(node);
		Vector<Node> v = n.getParents();

		if (v == null)
			return 0;

		return v.size();
	}

	public Vector<Node> getParents(int nodeId) {
		Integer node = new Integer(nodeId);

		if (!nodeList.containsKey(node))
			return null;

		Node n = nodeList.get(node);
		Vector<Node> v = n.getParents();
		return v;
	}

	public int getNumChildren(int nodeId) {
		Integer node = new Integer(nodeId);

		if (!nodeList.containsKey(node))
			return -1;

		Node n = nodeList.get(node);
		Vector<Node> v = n.getChildren();

		if (v == null)
			return 0;

		return v.size();
	}

	public int getNumLeaves() {
		int leaves = 0;
		Enumeration<Integer> keys = nodeList.keys();

		while (keys.hasMoreElements()) {
			if (getNumChildren(keys.nextElement().intValue()) == 0)
				leaves++;
		}

		return leaves;
	}

	public Vector<Integer> getAllLeaves() {
		if (leaves != null)
			return leaves;

		Vector<Integer> leavesl = new Vector<Integer>();
		Enumeration<Integer> keys = nodeList.keys();

		while (keys.hasMoreElements()) {
			Integer n = keys.nextElement();

			if (getNumChildren((n).intValue()) == 0) {
				leavesl.add(n);

			}
		}

		leaves = leavesl;
		return leavesl;
	}

	public Vector<Integer> getAllLeaves(Vector<Integer> descendants) {
		Vector<Integer> leaves = new Vector<Integer>();

		for (int i = 0; i < descendants.size(); i++) {
			Integer n = descendants.elementAt(i);

			if (getNumChildren((n).intValue()) == 0) {
				leaves.add(n);

			}
		}

		return leaves;
	}

	public Hashtable<Integer, ProbabilisticDAGNode> removeProbabilityInstance() {
		Hashtable<Integer, ProbabilisticDAGNode> plist = new Hashtable<Integer, ProbabilisticDAGNode>();
		Enumeration<Integer> probs = probList.keys();

		while (probs.hasMoreElements()) {
			Integer node = probs.nextElement();

			Object npObject = probList.get(node).clone();
			ProbabilisticDAGNode np = null;

			if (npObject instanceof ProbabilisticDAGNode) {
				np = (ProbabilisticDAGNode) npObject;
			}

			plist.put(node, np);
		}

		return plist;
	}

	// goes through nodes and makes sure that their evidence
	// likelihoods, and posteriors are all 0.
	public void zeroProbabilities() {
		Enumeration<Integer> probs = probList.keys();

		while (probs.hasMoreElements()) {
			Integer node = probs.nextElement();
			probList.get(node).resetProbs();
		}
	}

	// I am not sure we need this function, but here
	// is what we would use either way.
	public void reInsertProbabilityInstance(Hashtable<Integer, ProbabilisticDAGNode> probs) {
		probList = probs;
	}

	public void setPosterior(int nodeId, double likelihood) {
		Integer node = new Integer(nodeId);

		if (probList.containsKey(node)) {
			ProbabilisticDAGNode np = probList.get(node);

			if (np == null)
				System.out.println("Can't set posterior in ProbabilisticDAG.setPosterior");
			else
				np.setPosterior(likelihood);
		}
	}

	// This is an _approximation_ to two leaf node subsets
	// It checks to see if evidence associated with the
	// leaf and each (for three, pair of) other leaves
	// have evidence associated with them. If so, the
	// evidence is combined; otherwise it is added to
	// the appropriate leaf (pair) probability.

	// This function maps the integer onto the vector
	// of leaves that it has stored.
	public double getLeafSubsetProbability(Integer leaf, Vector<Integer> leaves) {
		// First check to see if leaf has its own
		// subset probability there already...

		ProbabilisticDAGNode leafNP = null;

		if (probList.containsKey(leaf)) {
			leafNP = probList.get(leaf);
			double sp = leafNP.getSubsetProbability();

			if (sp >= 0.0)
				return sp;
		}
		else {
			System.out.println("In getLeafSubsetProbability (ProbabilisticDAG): Cannot find leaf in probNode set: " + leaf);
		}


		// Leaf's own index is its single element subset
		double[] doubleElement = new double[leaves.size()];
		int leafIndex = 0;

		if (leaves.contains(leaf)) {
			leafIndex = leaves.indexOf(leaf);
		}
		else {
			System.out.println("In getLeafSubsetProbability (ProbabilisticDAG): Can't find leaf in set of leaves: " + leaf);
			return 0.0;
		}

		// Initialize
		for (int i = 0; i < leaves.size(); i++) {
			doubleElement[i] = 0.0;
		}

		// Check to see if there is evidence in leaf node too;
		double evidence = leafNP.getEvidence();

		if (evidence > 0.0) {
			doubleElement[leafIndex] = evidence;
		}

		// Find elements of induced subtree with evidence
		// If evidence exists, add it to the appropriate
		// leaves element.
		Vector<Integer> subgraph = getInducedSubgraph(leaf.intValue());

		for (int i = 0; i < subgraph.size(); i++) {
			Integer node = (Integer)subgraph.elementAt(i);

			if (probList.containsKey(node)) {
				ProbabilisticDAGNode n = probList.get(node);
				evidence = n.getEvidence();

				if (evidence > 0.0) {
					Vector<Integer> internalLeaves = n.getDescendantLeaves();

					if (internalLeaves == null) {
						System.out.println("In getLeafSubsetProbability (ProbabilisticDAG): no descendant leaves for node: " + node);
					}
					else {
						for (int j = 0; j < internalLeaves.size(); j++) {
							if (leaves.contains(internalLeaves.elementAt(j))) {
								int index =
								  leaves.indexOf(internalLeaves.elementAt(j));
								doubleElement[index] =
								  1.0
								  - ((1.0 - doubleElement[index])
								     * (1.0 - evidence));
							}
							else {
								System.out.println("In getLeafSubsetProbability (ProbabilisticDAG): internal leaves and descendant leaves do not match for node:" + node);
							}
						}
					}
				}
			}
		}

		double totalEvidence = 0.0;
		double singlePrior = getPolynomialMultiplier(1);
		double doublePrior = getPolynomialMultiplier(2);

		for (int i = 0; i < doubleElement.length; i++) {
			if (i != leafIndex) {
				totalEvidence += doubleElement[i] * doublePrior;
			}
			else {
				totalEvidence += doubleElement[i] * singlePrior;
			}
		}

		// Set it here so we don't have to calculate it
		// in the future.
		leafNP.setSubsetProbability(totalEvidence);
		return totalEvidence;
	}

	private double getPolynomialMultiplier(int desLeaves) {
		double total = 0;

		for (int i = 1; i <= desLeaves; i++) {
			double mult = GenericMathFunctions.choose(desLeaves, i);

			// again accounting for sluff variables.
			if (mult == desLeaves - 1)
				mult = mult + 1;

			total += mult * (1 / Math.pow(rValue, (double)i));
		}

		return total;
	}

	public class Node {
		// This is a basic DAG, taking on an object if needed.
		private String name;
		private int id;
		private Vector<Node> inEdges;
		private Vector<Node> outEdges;
		private Object additional;
		private int level; // minimum level. not really required.
		//private ProbabilisticDAGNode nprob;

		public Node(String nameS, int idS) {
			name = nameS;
			id = idS;
			inEdges = null;
			outEdges = null;
			additional = null;
			level = -1;
		}

		public int getId() {
			return id;
		}

		public int getLevel() {
			return level;
		}

		//if the level is being added, and there is already
		// one there, then take the smaller of the two.
		public void addLevel(int l) {

			if (level > -1) {
				if (level < l)
					return;
				else
					level = l;
			}

			level = l;
			return;
		}

		// I hope this is passed in Java
		// as a pointer. Otherwise this
		// type of edge is wrong.
		public void addInEdge(Node n) {
			if (inEdges == null)
				inEdges = new Vector<Node>();

			if (!inEdges.contains(n))
				inEdges.add(n);

			inEdges.trimToSize();
		}

		public void addOutEdge(Node n) {
			if (outEdges == null)
				outEdges = new Vector<Node>();

			if (!outEdges.contains(n))
				outEdges.add(n);

			outEdges.trimToSize();
		}

		public Vector<Node> getParents() {
			return inEdges;
		}

		public Vector<Node> getChildren() {
			return outEdges;
		}

		public void addAdditional(Object o) {
			additional = o;
		}

		public Object getAdditional() {
			return additional;
		}

		public String getName() {
			return name;
		}

		public void removeParent(Node n) {
			while (inEdges != null && inEdges.contains(n)) {
				inEdges.remove(n);
				inEdges.trimToSize();
			}

			if (inEdges.size() < 1)
				inEdges = null;

		}

		public void removeChild(Node n) {
			while (outEdges != null && outEdges.contains(n)) {
				outEdges.remove(n);
				outEdges.trimToSize();
			}

			if (outEdges.size() < 1)
				outEdges = null;
		}

		public void printOutNode() {
			System.out.println(name + ", " + id + ", " + inEdges + ", "
			                   + outEdges + ", " + additional);
		}

	}

	/**
	 * Method keys.
	 * @return Enumeration
	 */
	public Enumeration<Integer> keys() {
		if (nodeList != null)
			return nodeList.keys();

		return null;
	}

	public void printOutDAG() {
		Enumeration<Integer> nodes = nodeList.keys();

		while (nodes.hasMoreElements()) {
			Node n = nodeList.get(nodes.nextElement());
			n.printOutNode();
		}
	}

	public String toString() {
		String retval = "";

		Enumeration<Integer> nodes = nodeList.keys();

		while (nodes.hasMoreElements()) {
			Node n = nodeList.get(nodes.nextElement());
			retval += n.name + "\n";
		}

		return retval;
	}
	/**
	 * Method removeNodeParent.
	 * @param name
	 * @param node
	 */
	public void removeNodeParent(Integer name, Node node) {
		//if(!nodeList.containsKey(node)) return;
		Node n = nodeList.get(name);

		if (n != null)
			n.removeParent(node);
	}

	/**
	 *
	 */
	public void aPrioriEvidence() {
		// TODO Auto-generated method stub
		double total = 1.0;
		int noEvCount = 0;
		leaves = getAllLeaves();

		for (int i = 0; i < leaves.size(); i++) {
			double l = getLikelihood(leaves.elementAt(i).intValue());

			if (l > 0)
				total *= l;
			else
				noEvCount++;
		}

		double a = 1 / rValue;

		if (noEvCount > 0) {
			double rest = Math.pow(1.0 / rValue, leaves.size()) / total;
			a = Math.pow(rest, 1.0 / (double)noEvCount);
		}

		for (int i = 0; i < leaves.size(); i++) {
			double l = getLikelihood(leaves.elementAt(i).intValue());

			if (l <= 0) {
				ProbabilisticDAGNode n =
				  probList.get(leaves.elementAt(i));
				n.setLikelihood(a);
			}
		}
	}
}

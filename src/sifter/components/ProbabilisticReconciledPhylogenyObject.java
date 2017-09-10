/**
 * Data structure for reconciled phylogeny, but also does
 * whole propagation up and down tree.
 * Manipulates probabilities.
 *
 * Copyright 2003-2005 Barbara Engelhardt (bee@cs.berkeley.edu)
 * @author Barbara Engelhardt (bee@cs.berkeley.edu)
 */

package sifter.components;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;
import java.util.Vector;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;


public class ProbabilisticReconciledPhylogenyObject {
	private int root;
	private Vector<Node> tree;
	private Random rand;
	//private double leafPrior;
	private double singleLeafPrior;
	private double[] priorPolynomial;
	private int alignmentLength;
	private Hashtable<Node, double[]> gammasFinal;
	private Hashtable<Node, double[]> deltasFinal;
	//private Hashtable facListsRepository;
	//private int maxCheckSum;

	private static String SPECIATION_PARAMS = "species";
	private static String DUPLICATION_PARAMS = "duplication";

	public ProbabilisticReconciledPhylogenyObject(int alignmentLen) {
		root = -1;
		tree = null;
		rand = new Random();
		//leafPrior = 0;
		singleLeafPrior = 0;
		priorPolynomial = null;
		//facListsRepository = new Hashtable();
		//maxCheckSum = 11;
		alignmentLength = alignmentLen;
	}

	/* iterator with breadth first traversal of tree */
	public Vector<Node> getBFVector() {
		if (root < 0)
			return null;

		Vector<Node> bfVector = new Vector<Node>();
		bfVector.add(tree.elementAt(root));
		Vector<Node> nextLayer = tree.elementAt(root).getChildren();

		while (!nextLayer.isEmpty()) {
			Node n = nextLayer.elementAt(0);
			bfVector.add(n);
			nextLayer.removeElementAt(0);

			if (n.isLeaf() == false)
				nextLayer.addAll(n.getChildren());
		}

		return bfVector;
	}

	// iterator with depth first
	// traversal of tree
	public Vector<Node> getDFVector() {
		if (root < 0)
			return null;

		Vector<Node> dfVector = new Vector<Node>();
		dfVector.add(tree.elementAt(root));
		Vector<Node> nextLayer = tree.elementAt(root).getChildren();

		while (!nextLayer.isEmpty()) {
			Node n = nextLayer.elementAt(0);
			dfVector.add(n);
			nextLayer.removeElementAt(0);

			if (n.isLeaf() == false)
				nextLayer.addAll(0, n.getChildren());
		}

		return dfVector;
	}


	///////////////////////////////////////////////////////////
	//   Compiling statistics of tree
	///////////////////////////////////////////////////////////

	public int getNumNodes() {
		if (tree != null)
			return tree.size();

		return 0;
	}

	///////////////////////////////////////////////////////////
	//   Reading in Reconciled Tree
	///////////////////////////////////////////////////////////


	public void createReconciled(String filename) {
		tree = new Vector<Node>();

//		String Dataset = null;
		if (filename != null) {
//			try{
//				filename="./data/test_nhx_xml/reconciled-pf00685.nhx";
//				Dataset = fileAsString(filename);
//			} catch (IOException e) {
//				System.err.println("Couldn't get data file "
//						+filename+" from the dir");
//				System.exit(1);
//			}
		}
		else {
			System.out.println("Error in createReconciled: "
			                   + "reconciled tree filename is null");
			System.exit(1);
		}

//		Dataset = Dataset.toLowerCase();

		// Modify dataset a bit

		//System.out.println("The input (modified) dataset is:\n" + Dataset);
		try {
			readInSpecies(filename);
		}
		catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}

		for (int i = 0; i < tree.size(); i++) {
			if (tree.elementAt(i).isRoot()) {
				root = i;
				return;
			}
		}

	}

	public Vector<Node> getTree() {
		return tree;
	}

	public void readInSpecies(String filename) throws Exception {
		//maptoTree(input);

//		// alternative way
		tree = new Vector<Node>();

		PhylogenyParser parser = new PhyloXmlParser();
		File file = new File(filename);
		final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
		final Phylogeny[] phylogenies = factory.create(file, parser);

		if ((phylogenies == null) || (phylogenies.length == 0)) {
			throw new Exception("Unable to parse phylogeny from file: " + file);
		}

		Phylogeny phylogeny = phylogenies[0];

		for (final PhylogenyNodeIterator it = phylogeny.iteratorPreorder(); it.hasNext();) {
			final PhylogenyNode node = it.next();

			Node p = node.getParent() != null ? this.getNode(node.getParent().getNodeId()) : null;
			String nodeName = node.getNodeName() != "" ? node.getNodeName().toLowerCase() : "Node" + node.getNodeId();
			Node n = new Node(nodeName, p);
			n.setParentDistance(node.getDistanceToParent(), alignmentLength);
			boolean isDuplication = node.getNodeData().getEvent() != null ? node.getNodeData().getEvent().isDuplication() : false;
			n.setDuplicate(isDuplication);
			tree.addElement(n);

			if (p != null)
				p.addChildren(n);
		}

		for (int i = 0; i < tree.size(); i++) {
			if (tree.elementAt(i).isRoot()) {
				root = i;
				break;
			}
		}

	}

	public void maptoTree(String doc) {
		if (!tree.isEmpty())
			tree.removeAllElements();

		//int Begin = doc.indexOf("(");
		int Begin = -1;
		int current = 0;
		Node currentNode = null;
		Node parent = null;
		int next = 0;

		while (doc.length() > 1) {
			doc = doc.substring(Begin + 1, doc.length());
			// Begin can either be the next space, ) or (
			int nextParen = doc.indexOf("(");
			int nextSpace = doc.indexOf(",");
			int nextClose = doc.indexOf(")");
			int lastFinal = doc.indexOf(";");

			if (nextParen == -1)
				nextParen = lastFinal;

			if (nextSpace == -1)
				nextSpace = lastFinal;

			if (nextClose == -1)
				nextClose = lastFinal;

			next = Math.min(nextParen, nextSpace);
			next = Math.min(next, nextClose);

			if (next < 0 || next == lastFinal) {
				if (doc.indexOf(']') > 0)
					parseClose(doc.substring(0, doc.indexOf(']')),
					           currentNode);

				return;
			}

			Begin = 0;
			//System.out.println(next+": "+doc);

			// Pull the next first node if it exists
			//if(next - Begin > 1) {
			String name = doc.substring(0, next);
			Node n = null;

			if (doc.indexOf(':') != 0 && doc.indexOf('[') != 0) {
				if (name.indexOf('[') >= 0) {
					name = name.substring(0, name.indexOf('['));
				}

				if (name.equals("")) {
					name = "Node" + current;
					current++;
				}

				n = new Node(name, parent);

				if (parent != null)
					parent.addChildren(n);

				//System.out.println(doc.substring(0, next));
				tree.addElement(n);
			}

			if (next == nextParen) {
				parent = n;
			}

			// might have to add non-zero length cond.
			if (next == nextSpace) {
				if (currentNode != null &&
				    (doc.indexOf(':') == 0 || doc.indexOf('[') == 0))
					parseSpace(doc.substring(0, doc.indexOf(']')),
					           currentNode);
				else
					parseSpace(doc.substring(0, doc.indexOf(']')), n);
			}

			if (next == nextClose) {
				if (currentNode != null && currentNode.getParent() == null)
					System.out.println("Mismatched parens in species file");
				else {
					if (doc.indexOf(':') == 0) {
						parseClose(doc.substring(0, doc.indexOf(']')),
						           currentNode);
					}
					else {
						parseClose(doc.substring(0, doc.indexOf(']')), parent);
					}

					currentNode = parent;
					parent = parent.getParent();

					if (doc.indexOf(',') == 0) {
						doc = doc.substring(1);
					}
				}
			}

			Begin = next;
		}
	}

	public void parseClose(String close, Node c) {
		if (((String)c.getNodeID()).startsWith("Node") != true) { //"OST6_MOUSE"))
			System.out.println("[DEBUG] looking at ost6_mouse");
		}

		if (close.indexOf(':') == 0 || close.indexOf('[') == 0) {
			if (close.indexOf(':') < close.indexOf('[')) {
				String dist = close.substring(close.indexOf(':') + 1,
				                              close.indexOf('['));
				double parentDist = new Double(dist).doubleValue();
				c.setParentDistance(parentDist, alignmentLength);
			}
			else {
				c.setParentDistance(-1, alignmentLength);
			}

			if (close.indexOf("d=") > 0) {
				String dup = close.substring(close.indexOf("d=") + 2,
				                             close.indexOf("d=") + 3);

				if (dup.equals("y")) {
					c.setDuplicate(true);
				}
			}

			//System.out.println("Closed final for node "+c.getNodeID());
		}
		else
			if ((close.indexOf(':') > 0 && close.indexOf('[') != 0)
			    && close.indexOf(':') < close.indexOf('[')) {
				c.setParentDistance(Double.parseDouble(
				                      close.substring(close.indexOf(':') + 1,
				                                      close.indexOf('['))),
				                    alignmentLength);

				if (close.indexOf("d=") > 0) {
					String dup = close.substring(close.indexOf("d=") + 2,
					                             close.indexOf("d=") + 3);

					if (dup.equals("y")) {
						c.setDuplicate(true);
					}
				}
			}
			else
				if ((close.indexOf(':') > 0 && close.indexOf('[') != 0)) {
					c.setParentDistance(-1, alignmentLength);

					if (close.indexOf("d=") > 0) {
						String dup = close.substring(close.indexOf("d=") + 2,
						                             close.indexOf("d=") + 3);

						if (dup.equals("y")) {
							c.setDuplicate(true);
						}
					}
				}
				else {
					c.setParentDistance(1.0, alignmentLength);

					if (close.indexOf("d=") > 0) {
						String dup = close.substring(close.indexOf("d=") + 2,
						                             close.indexOf("d=") + 3);

						if (dup.equals("y")) {
							c.setDuplicate(true);
						}
					}

				}

	}

	public void parseSpace(String close, Node c) {
//		if (((String)c.getNodeID()).startsWith("ost6_mouse"))
//		{
//			System.out.println("[DEBUG] looking at ost6_mouse");
//		}

		if (close.indexOf(':') >= 0) {
			if (close.indexOf(':') < close.indexOf('[')) {
				String dist = close.substring(close.indexOf(':') + 1,
				                              close.indexOf('['));
				double parentDist = new Double(dist).doubleValue();
				c.setParentDistance(parentDist, alignmentLength);
			}
			else {
				c.setParentDistance(-1, alignmentLength);
			}

			if (close.indexOf("d=") > 0) {
				String dup = close.substring(close.indexOf("d=") + 2,
				                             close.indexOf("d=") + 3);

				if (dup.equals("y")) {
					c.setDuplicate(true);
				}
			}
		}
	}

	public void printTree() {
		Vector<Node> t = tree;
		System.out.println("The species doc tree looks like:\n");
		int i;
		int size = t.size();

		for (i = 0; i < size; i++) {
			String value = (String)t.elementAt(i).getNodeID();
			String nodenum = (new Integer(i)).toString();
			double dist = t.elementAt(i).getParentDistance();
			boolean dup = t.elementAt(i).hasDuplication();

			if (t.elementAt(i).getParent() != null) {
				System.out.println("(Node " + nodenum + "): "
				                   + value + ", parent: "
				                   + ((String)t.elementAt(i)
				                      .getParent().getNodeID())
				                   + " dist: " + dist + " dup?: " + dup);
			}
			else {
				System.out.println("(Node " + nodenum + "): "
				                   + value + " dup?: " + dup);
			}
		}
	}
	
	// TODO
	public void printTreeAsFile(String filename) {
		
		File o = new File(filename);
    	
    	//Create output file
    	FileWriter fw;
        if (!o.exists()) 
        {
            try 
            {
                o.createNewFile();
            } 
            catch (IOException ex) 
            {
                System.out.println("\n\nFile " + o);
                System.exit(0);
            }
        }
        BufferedWriter bw = null;
        
        try 
    	{
    		fw = new FileWriter(o.getAbsoluteFile(), true);
    		bw = new BufferedWriter(fw);
    		Vector<Node> t = tree;
    		int i;
    		int size = t.size();
    		
    		for (i = 0; i < size; i++) {
    			String value = (String)t.elementAt(i).getNodeID();
    			String nodenum = (new Integer(i)).toString();
    			double dist = t.elementAt(i).getParentDistance();
    			boolean dup = t.elementAt(i).hasDuplication();

    			if (t.elementAt(i).getParent() != null) {
    				bw.write(value + ", parent: " + ((String)t.elementAt(i).getParent().getNodeID()) + 
    						", probabilities: " + Arrays.toString(t.elementAt(i).getLocalProbabilities()) + " dist: " + dist + 
    						" dup?: " + dup + " " + " leaf?: " + t.elementAt(i).isLeaf() + "\n");
    				bw.flush();
    			}
    			else {
    				System.out.println("(Node " + nodenum + "): "
    				                   + value + " dup?: " + dup);
    			}
    		}
    	}
        catch (IOException e) 
    	{
			e.printStackTrace();
		}
	}

	@SuppressWarnings("resource")
	public String fileAsString(String fname) throws IOException {
		File pn = new File(fname);
		FileReader frdr = new FileReader(pn);
		StringBuffer buf = new StringBuffer();
		int ch;

		do {
			ch = frdr.read();

			if (ch != -1 && ch != '\n')
				buf.append((char)ch);
		}
		while (ch != -1);

		return buf.toString();
	}

	public void stringToFile(String input, String fname) throws IOException {
		File outputFile = new File(fname);
		FileWriter out = new FileWriter(outputFile);
		int c;
		int i = 0;
		System.out.println("begin write file\n");
		c = input.charAt(i);

		//   while (c != -1 & c != '\r') //for linux
		while (c != -1 & i < input.length() - 1) { //for windows and linux
			out.write(c);
			i = i + 1;
			c = input.charAt(i);
		}

		System.out.println("finish write: " + fname + "\n");
		out.close();
	}

	/* converts the tree in memory to a NXH file */
	public void treeToNHXFile(String fname) throws IOException {
		File outputFile = new File(fname);
		FileWriter out = new FileWriter(outputFile);

		System.out.println("begin write file\n");

		Node n = tree.elementAt(0);
		Integer prevLevel = new Integer(0);
		Integer currentLevel = null;

		if (n.getParent() == null) {
			// root, don't print out name
			//out.write("(");
			Vector<Node> nextList = n.getChildren();
			Vector<Integer> nextLevel = new Vector<Integer>();
			int levelIndex = 0;

			for (int i = 0; i < nextList.size(); i++) {
				nextLevel.add(new Integer(1));
			}

			while (!nextList.isEmpty()) {
				Node current = nextList.elementAt(0);
				currentLevel = nextLevel.elementAt(levelIndex);
				nextList.removeElementAt(0);
				levelIndex++;

				if (currentLevel.compareTo(prevLevel) == 0) {
					out.write(",");
					prevLevel = currentLevel;
				}
				else
					if (currentLevel.compareTo(prevLevel) < 0) {
						int backtrack = prevLevel.intValue() -
						                currentLevel.intValue();

						for (int j = 0; j < backtrack; j++) {
							out.write("):5");
						}

						prevLevel = new Integer(currentLevel.intValue()
						                        - backtrack);

						if (currentLevel.compareTo(prevLevel) >= 0) {
							out.write(",");
							prevLevel = currentLevel;
						}
					}
					else {
						out.write("(");
						prevLevel = currentLevel;
					}

				if (current.getChildren() != null) {
					nextList.addAll(0, current.getChildren());

					for (int i = 0; i < current.getChildren().size(); i++) {
						nextLevel.add(levelIndex,
						              new Integer(currentLevel.intValue() + 1));
					}
				}
				else {
					out.write(((String)current.getNodeID()).toUpperCase());
				}
			}
		}

		for (int i = 0; i < currentLevel.intValue() - 1; i++) {
			// TODO: make correct branch length
			out.write("):5");
		}

		out.write(");");
		System.out.println("finish write: " + fname + "\n");
		out.close();
	}


	/////////////////////////////////////////////////////////
	//  Accessor functions                                 //
	/////////////////////////////////////////////////////////

	public boolean hasNode(String nID) {
		String n2 = nID.toLowerCase();

		for (int i = 0; i < tree.size(); i++) {
			Node n = tree.elementAt(i);

			if (((String)n.getNodeID()).startsWith(n2)) {
				return true;
			}
		}

		System.out.println("In hasNode (ProbabilisticReconciledPhylogenyObject): "
		                   + "did not find node in tree with ID: " + n2);
		return false;
	}

	public Node getNode(String nID) {
		String n2 = nID.toLowerCase();

		for (int i = 0; i < tree.size(); i++) {
			Node n = tree.elementAt(i);

			if (((String)n.getNodeID()).startsWith(n2)) {
				return n;
			}
		}

		System.out.println("In getNode (ProbabilisticReconciledPhylogenyObject): "
		                   + "did not find node in tree with ID: " + n2);
		return null;
	}
	
	public Vector<Node> getNodes(String nID) {
		String n2 = nID.toLowerCase();

		Vector<Node> ns = new Vector<Node>();
		
		
		for (int i = 0; i < tree.size(); i++) {
			Node n = tree.elementAt(i);

			if (((String)n.getNodeID()).startsWith(n2)) {
				ns.add(n);
			}
		}
		
		if (ns.isEmpty()){
			System.out.println("In getNode (ProbabilisticReconciledPhylogenyObject): "
		                   + "did not find node in tree with ID: " + n2);
			return null;
		}else{
			return ns;			
		}
	}
	
	

	public Node getNode(int nID) {
		String n2 = "Node" + nID;

		for (int i = 0; i < tree.size(); i++) {
			Node n = tree.elementAt(i);

			if (((String)n.getNodeID()).startsWith(n2)) {
				return n;
			}
		}

		System.out.println("In getNode (ProbabilisticReconciledPhylogenyObject): "
		                   + "did not find node in tree with ID: " + n2);
		return null;
	}


	/**
	 * Method getNode.
	 * @param root
	 * @return Node
	 */
	/*
	  private Node getNode(int num) {
	Node n = (Node)tree.elementAt(num);
	if(n.getParent() == null) return null;
	else return n;
	}*/

	
	
	public void setNodeEvidenceProbabilities(String nID, double[] ep) {
		if (ep == null)
			return;

		Vector<Node> ns = getNodes(nID);

		if (ns == null)
			return;
		
		for (int i = 0; i < ns.size(); i++) {
			Node n = ns.elementAt(i);
			//System.out.println(n+nID+","+(double[])ep.clone()+","+((double[])ep.clone()!= null));
			n.setLocalProbabilities(ep);
		}
	}
	
	/*
	  public void setNodeEvidenceProbabilities(String nID, double[] ep) {
		if (ep == null)
			return;

		Node n = getNode(nID);

		if (n == null)
			return;
		System.out.println(n+nID+","+(double[])ep.clone()+","+((double[])ep.clone()!= null));
		n.setLocalProbabilities(ep);
	}*/

	//////////////////////////////////////////////////////////
	//  Node Class declaration                              //
	//////////////////////////////////////////////////////////

	public class Node {
		public Object obj;
		public boolean fixed;
		public boolean duplication;
		public Node parent;
		public double parentDistance;
		public Vector<Node> children;
		private double[] localProbabilities;

		public Node(Object o) {
			obj = o;
			fixed = false;
			parent = null;
			parentDistance = 1;
			children = null;
			duplication = false;
			localProbabilities = null;
		}

		/**
		 * Constructor Node.
		 * @param string
		 * @param p
		 */
		public Node(Object value, Node p) {
			parent = p;
			parentDistance = 1;
			duplication = false;
			obj = parseNHXName((String)value);
			children = null;
			localProbabilities = null;
		}

		private String parseNHXName(String v) {
			String n = v;
			int colonIndex = v.indexOf(':');

			if (colonIndex > 0) {

				n = v.substring(0, colonIndex);

				if (colonIndex + 1 < v.indexOf('[')) {
					String dist = v.substring(colonIndex + 1, v.indexOf('['));
					double distTemp = (new Double(dist)).doubleValue();
					setParentDistance(distTemp, alignmentLength);
				}
				else {
					setParentDistance(-1, alignmentLength);
				}
			}

			return n;
		}

		public void setNodeID(Object id, boolean f) {
			if (fixed == true)
				return;

			obj = id;
			fixed = f;
		}

		public Object getNodeID() {
			return obj;
		}

		public void addParent(Node p, double pd) {
			parent = p;
			p.addChildren(this);
			setParentDistance(pd, alignmentLength);
		}

		public Node getParent() {
			return parent;
		}

		public boolean isRoot() {
			return(parent == null);
		}

		public boolean isLeaf() {
			return(children == null);
		}

		public void setDuplicate(boolean dup) {
			duplication = dup;
		}

		public boolean hasDuplication() {
			return duplication;
		}

		public void addChildren(Node c) {
			if (children == null)
				children = new Vector<Node>();

			children.add(c);
			children.trimToSize();
		}

		public double getParentDistance() {
			return parentDistance;
		}

		// This is the unfortunate place where I am
		// putting the distance function
		// as the d_i parameter goes down, uncertainty (p(x=1))
		// goes down.
		public void setParentDistance(double pd, int alignLength) {
			if ((pd >= 1.0 || pd == 0.0) && alignLength > 0) {
				parentDistance = (double)(pd) / (double)alignLength;
			}
			else
				if (pd < 1.0 && pd > 0.0) {
					parentDistance = pd;
				}

			//highest/lowest it should get
			//if(parentDistance > 1.0) parentDistance = 1.0;
			if (parentDistance <= 0.0 && alignLength > 0)
				parentDistance = 1.0 / (double)alignLength;

			//parentDistance = Math.sqrt(parentDistance);
			if (parentDistance > 1.0)
				parentDistance = 1.0;
		}

		public Vector<Node> getChildren() {
			return children;
		}

		public void setLocalProbabilities(double[] lp) {
			localProbabilities = (double[])lp.clone();
		}

		public double[] getLocalProbabilities() {
			return localProbabilities;
		}

		public void removeLocalProbabilities() {
			localProbabilities = null;
		}

		public boolean hasLocalProbabilities() {
			return(localProbabilities != null);
		}

		/**
		 *
		 */
		public void printParentsChildren() {
			// TODO Auto-generated method stub
			System.out.print(obj + " has parent " + parent.getNodeID()
			                 + " and children ");

			if (children != null) {

				for (int i = 0; i < children.size(); i++) {
					System.out.print(children.elementAt(i) + " ");
				}
			}

			System.out.println();
		}


		public String toString() {
			String retval = "";

			retval += obj.toString();

			return retval;
		}

	}

	public int throwDart(double total, double[] board) {
		double eye = rand.nextDouble() * total;
		double area = 0.0;

		for (int i = 0; i < board.length; i++) {
			area += board[i];

			if (area > eye)
				return i + 1;
		}

		return (board.length - 1);
	}


	public int throwDartNegExp(double total, double[] board) {
		double eye = rand.nextDouble() * total;
		double area = 0.0;

		for (int i = 0; i < board.length; i++) {
			area += (1 - Math.exp(board[i]));

			if (area > eye)
				return i;
		}

		return (board.length - 1);
	}

	public void initializePriorPolynomial(int len) {
		if (priorPolynomial != null)
			return;

		priorPolynomial = new double[len];

		for (int i = 0; i < len; i++) {
			priorPolynomial[i] = Math.pow(singleLeafPrior, i + 1)
			                     * choose(len, i + 1);
		}
	}

	public double getPrior(ProteinFunctionMarkovState ps) {
		int numWithout = 0;

		for (int i = 0; i < ps.length(); i++) {
			if (ps.elementAt(i) == 1) {
				numWithout++;
			}
		}

		if (numWithout == 0)
			return 0.0;

		return (Math.pow(singleLeafPrior, numWithout));
	}

	/**
	 * Method setSingleNodeSampleLikelihoods.
	 */
	public void setSingleNodeSampleLikelihoods(double lprior,
	    double singlelprior) {
		//leafPrior = lprior;
		singleLeafPrior = singlelprior;
	}

	/**
	 * performs exact inference in the tree, given the speciation
	 * transition matrix and the duplication transition matrix.
	 */
	@SuppressWarnings("unused")
	public Hashtable<Node, double[]> propagateExactThroughoutTree(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx,
	    Hashtable<String, Double> scaleParams) {
		//if(true) return null;
		//Propagate up tree
		double[] evidenceProbs = null;
		double[] gamma = null;
		double[] deltaRoot = null;
		double[] prior = new double[1];
		prior[0] = 1;
		Node rootNode = null;
		Vector<Node> messages = new Vector<Node>();
		Hashtable<Node, double[]> gammaList = new Hashtable<Node, double[]>();
		Hashtable<Node, double[]> deltaList = new Hashtable<Node, double[]>();

		// bee bmc added line
		// precompute the matrix exponential
		pfx.buildMarkovTransitionRateMatrix();

		//Get children with evidence first.
		//since the tree is constructed depth first,
		//this should be the correct way to go through the tree.
		for (int node = tree.size() - 1; node >= 0; node--) {
			// Get the name of the node
			Node n = tree.elementAt(node);
			Node parent = n.getParent();
			//if (parent !=null){
			//	System.out.println("good");
			//}
			if (n == null) {
				System.out.println("Can't find node in getNewNodeSampleFinal");
			}

			//Case 1: Leaf with evidence
			//sample directly from those local probabilities

			//System.out.println(node+","+n.getNodeID()+","+n.hasLocalProbabilities()+","+n.isLeaf());
			
			if (n.hasLocalProbabilities() && n.isLeaf()) {
				evidenceProbs = n.getLocalProbabilities();

				for (int i = 0; i < evidenceProbs.length; i++) {
					if (evidenceProbs[i] > 0) {
						evidenceProbs[i] = logSafe(evidenceProbs[i]);
					}
				}

				gammaList.put(n, evidenceProbs);
				messages.add(n);
				
				
				//System.out.println("Added gamma for: "+n.getNodeID());
			} // end of leaf with evidence
			else
				if (n.isLeaf()) {
					gammaList.put(n, prior);
					messages.add(n);
					//if(verbose)
					// System.out.println("Added gamma for; "+n.getNodeID());
				} // end of leaf without evidence
		}// end of initial leaf search

		while (!messages.isEmpty()) {
			boolean done = false;
			int count = 0;

			while (!done && messages.size() > 0) {
				count = (count + 1) % messages.size();
				//bee good debug line
				//System.out.print("Messages left gamma ("+count+") ");
				//for(int i = 0; i < messages.size(); i++){
				//	System.out.print(((Node)messages.elementAt(i)).getNodeID()+" ");
				//}
				//System.out.println();
				Node parent = messages.elementAt(count).getParent();

				if (parent != null) {
					Vector<Node> children = (Vector<Node>)parent.getChildren();

					if (messagesContainsAllChildren(messages, children)) {
						double inUse = ((Double)scaleParams.get(SPECIATION_PARAMS)).doubleValue();

						if (parent.hasDuplication())
							inUse = ((Double)scaleParams.get(DUPLICATION_PARAMS)).doubleValue();

						//System.out.println("Working on "+parent.getNodeID());
						gamma = gammaPropBinaryExp(gammaList, children,
						                           pfx, inUse);

						if (gamma.length == 1)
							gamma = prior;

						messages.removeAll(children);

						if (!parent.isRoot()) {
							gammaList.put(parent, gamma);
							messages.add(parent);
							//if(verbose)
							//System.out.println("Added gamma for "+parent.getNodeID());
						}
						else { // For the root
							gammaList.put(parent, gamma);
							deltaRoot = new double[gamma.length];

							for (int r = 0; r < gamma.length; r++)
								deltaRoot[r] = 0.0;

							deltaList.put(parent, deltaRoot);
							//if(verbose)
							//System.out.println("Added gamma for root"+gamma);
							messages.clear();
							rootNode = parent;
							done = true;
						}
					} // finished that message
				}
			}
		}

		// Propagating down from root.
		messages.add(rootNode);

		while (!messages.isEmpty()) {
			// bee helpful debug line
			//System.out.print("Messages left delta ("+messages.size()+") ");
			//for(int i = 0; i < messages.size(); i++) {
			//	System.out.print(((Node)messages.elementAt(i)).getNodeID()+" ");
			//}
			//System.out.println();
			// bee end helpful debug line

			Node parent = messages.elementAt(0);
			Vector<Node> children = parent.getChildren();
			double inUse = ((Double)scaleParams.get(SPECIATION_PARAMS)).doubleValue();

			if (parent.hasDuplication())
				inUse = ((Double)scaleParams.get(DUPLICATION_PARAMS)).doubleValue();

			for (int i = 0; i < children.size(); i++) {
				double[] delta0;
				delta0 = deltaPropBinaryExp(gammaList, children,
				                            children.elementAt(i),
				                            pfx, inUse,
				                            deltaList.get(parent));
				deltaList.put(children.elementAt(i), delta0);
				//if(verbose)
				//	System.out.println("Added delta for "
				//	                   + children.elementAt(i).getNodeID());

				if (!children.elementAt(i).isLeaf())
					messages.add(children.elementAt(i));
			}

			messages.remove(parent);
		}

		//Have gamma/delta. Multiply the two vectors and normalize
		Hashtable<Node, double[]> posteriors = new Hashtable<Node, double[]>();

		for (int node = tree.size() - 1; node >= 0; node--) {
			Node n = tree.elementAt(node);
			System.out.println("Pringing out "+n.getNodeID()); // TODO
			posteriors.put(n, gammaDeltaBinaryLog(gammaList.get(n),
			                                      deltaList.get(n)));
		}

		gammasFinal = gammaList;
		deltasFinal = deltaList;
		return posteriors;
	}

	// delta is defined as exp{-\sum{ \theta_{m,n}^{d_i} x_{\pi_i}^m}
	private double[] deltaPropBinaryExp(Hashtable<Node, double[]> gammaList,
	                                    Vector<Node> children, Node node,
	                                    GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx,
	                                    double scale_param,
	                                    double[] delta1) {
		int len2 = delta1.length;
		double[] delta = new double[len2];
		double[] notDelta = new double[len2];
		double distance = node.getParentDistance();
		ProteinFunctionMarkovState psParent = new ProteinFunctionMarkovState(len2, pfx.maxFunctions());
		double[] gammaParent = new double[len2];
		// product over sums of gamma functions
		// to incorporate siblings annotations
		int index = children.indexOf(node);
		children.remove(index);
		//System.out.println("Distance: "+distance+", rate: "
		//		   +scale_param+", yn: "+ pfx.getDelta(1,1));
		gammaParent = gammaPropBinaryExp(gammaList, children,
		                                 pfx, scale_param);
		children.add(index, node);

		if (gammaParent.length == 1) {
			gammaParent = new double[len2];

			for (int i = 0; i < gammaParent.length; i++) {
				gammaParent[0] = 0.0; // log of 1

			}

		}

		// incorporate deltas of parent in position-appropriate way.
		while (psParent.hasNext()) {
			psParent.getNextNonZeros();
			double[] deltaScale = new double[delta.length];
			double[] deltaNotScale = new double[delta.length];

			for (int j = 0; j < deltaScale.length; j++) {
				deltaScale[j] = 1.0;
				deltaNotScale[j] = 1.0;
			}

			double deltaProd = 1.0;

			for (int i = 0; i < psParent.length(); i++) {
				// incorporate delta of parent
				if (psParent.elementAt(i) == 1) {
					deltaProd *= Math.exp(delta1[i]) * Math.exp(gammaParent[i]);
				}
				else {
					double negGammaParent = 1 - Math.exp(gammaParent[i]);

					if (negGammaParent <= 0.0)
						negGammaParent = 1.0;

					double negDelta1 = 1 - Math.exp(delta1[i]);

					if (negDelta1 <= 0.0)
						negDelta1 = 1.0;

					deltaProd *= negDelta1 * negGammaParent;
				}

				//System.out.println("Deltaprod: "+deltaProd+", "
				//		   +Math.exp(delta1[i])+", "
				//		   +Math.exp(gammaParent[i]));
			}

			// incorporate deltas of parent in position-appropriate way.
			//psParent.printPowerSet();
			for (int j = 0; j < delta.length; j++) {
				double scale = probChildGivenParentME(1, j,
				                                      psParent, pfx,
				                                      distance,
				                                      scale_param);
				deltaScale[j] *= scale * deltaProd;
				deltaNotScale[j] *= (1.0 - scale) * deltaProd;
				//System.out.println("Scale ("+j+"): "+scale);
			}

			// Sum over all possible parents
			for (int j = 0; j < delta.length; j++) {
				delta[j] += deltaScale[j];
				notDelta[j] += deltaNotScale[j];
				deltaScale[j] = 1.0;
				deltaNotScale[j] = 1.0;
			}
		} // done with all possible parents

		// normalize, log
		for (int j = 0; j < delta.length; j++) {
			//System.out.println("Delta (from parent): "+delta[j]
			// +", not delta: "+notDelta[j]);
			delta[j] = logSafe(delta[j] / (delta[j] + notDelta[j]));
			//System.out.println("total (from parent):"+Math.exp(delta[j]));
		}

		return delta;
	}


	// The conditional probabilities in the model are
	// computed here.
	// yunesj: unused
	//    private double probChildGivenParent(int child, int functionIndex,
	//					ProteinFunctionMarkovState parentsSet,
	//					GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx, double distance,
	//					double rate)
	//    {
	//	double parentProd = 1.0;
	//	// is immediate parent 1?
	//	if(parentsSet.elementAt(functionIndex) == 1) {
	//	    parentProd = (Math.exp(-pfx.getDelta(functionIndex,
	//						functionIndex))
	//			  * (1.0 - Math.exp(-distance*rate)));
	//	    if(child == 1)
	//		parentProd = 1.0 - parentProd;
	//	} else { // immediate parent 0
	//	    parentProd = (1.0 - Math.exp(-distance * rate));
	//	    for(int k = 0; k < parentsSet.length(); k++) {
	//		if(parentsSet.elementAt(k) == 1)
	//		    parentProd *= Math.exp(-pfx.getDelta(k,functionIndex));
	//	    }
	//	    if(child == 0) {
	//		parentProd = 1.0 - parentProd;
	//	    }
	//	}
	//	return parentProd;
	//    }

	// The conditional probabilities in the model are
	// computed here.
	// Takes the exponentiated matrix (with rate and distance)
	// and returns the appropriate index.
	private double probChildGivenParentME(int child, int functionIndex,
	                                      ProteinFunctionMarkovState parentsSet,
	                                      GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx,
	                                      double rate, double distance) {
		if (child == 1) {
			return pfx.getExpProb(parentsSet.setIndex(), functionIndex,
			                      rate, distance);
		}
		else
			return (1.0 - pfx.getExpProb(parentsSet.setIndex(), functionIndex,
			                             rate, distance));
	}

	/* Computes gamma for a single node
	 */
	private double[] gammaPropBinaryExp(Hashtable<Node, double[]> gammaList,
	                                    Vector<Node> children, GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx,
	                                    double scale) {
		Vector<double[]> gammaSet = new Vector<double[]>();
		ProteinFunctionMarkovState ps = null;
		double distance = 0.0;
		Vector<Node> prunedChildren = new Vector<Node>();

		// pull out children with interesting gammas
		for (int i = 0; i < children.size(); i++) {
			Node child = children.elementAt(i);
			double[] gammaChild = gammaList.get(child);

			if (gammaChild.length > 1 && gammaChild[0] <= 0) {
				gammaSet.add(gammaChild);
				prunedChildren.add(child);
				//System.out.println("  Adding:	"+child.getNodeID());
			}
		}

		// If interesting gamma children list is size 0,
		// exit, setting gamma to uninteresting.
		if (gammaSet.size() == 0) {
			double[] gamma = new double[1];
			gamma[0] = 1.0;
			return gamma;
		}

		// gammaSet.size() >= 1 now, get length of gamma
		int len0 = gammaSet.elementAt(0).length;
		double[] gamma = new double[len0];
		double[] gammaNot = new double[len0];
		double[] gammaSum = new double[len0];
		double[] gammaNotSum = new double[len0];
		ps = new ProteinFunctionMarkovState(len0, pfx.maxFunctions());

		for (int i = 0; i < gamma.length; i++) {
			gamma[i] = 1.0;
			gammaSum[i] = 0.0;
		}

		// for each possible values over all of the children i
		// sum over p(X_{\pi i}|x_i)p(x_i|D_i)
		// where the first expression is a function computed
		// elsewhere, and the second is gamma.
		//ProteinFunctionMarkovState psChildren = new ProteinFunctionMarkovState(gammaSet.size(), len0);
		ProteinFunctionMarkovState psChildren = new ProteinFunctionMarkovState(len0, pfx.maxFunctions());
		double[] gammaTemp = new double[len0];
		double[] gammaTempNot = new double[len0];

		for (int i = 0; i < gamma.length; i++) {
			gammaNot[i] = gamma[i];
			gammaNotSum[i] = gammaSum[i];
		}

		Vector<Node> child = new Vector<Node>();

		for (int k = 0; k < prunedChildren.size(); k++) {
			psChildren.Reset();
			child.clear();
			child.add(prunedChildren.elementAt(k));
			distance =
			  prunedChildren.elementAt(k).getParentDistance();

			while (psChildren.hasNext()) {
				psChildren.getNextNonZeros();
				//System.out.print("Child: ");
				//psChildren.printPowerSet();

				probParentsGivenKidsME(gammaTemp, child,
				                       psChildren, pfx, scale, distance);

				for (int i = 0; i < gammaTemp.length; i++) {
					gammaTempNot[i] = 1.0 - gammaTemp[i];
					//System.out.println("Gamma temp"+i+", "+gammaTemp[i]);
				}

				// incorporate the gamma probabilities at the
				// children nodes now
				//System.out.println("gammaTemp: "+gammaTemp[0]+","
				//		   +gammaTemp[1]);
				for (int i = 0; i < gamma.length; i++) {
					//System.out.println("\ni= "+i+" gamma ="
					//+gammaTemp[i]+" notGamma = "+gammaTempNot[i]);
					//psChildren.printPowerSet();
					for (int j = 0; j < psChildren.length(); j++) {
						double[] gammaj = gammaSet.elementAt(k);
						double gammajTerm = 0;
						gammajTerm =
						  Math.exp(gammaj[psChildren.functionIndex(j)]);

						if (psChildren.elementAt(j) == 1) {
							gammaTemp[i] *= (gammajTerm);
							gammaTempNot[i] *= gammajTerm;
						}
						else {
							gammaTemp[i] *= (1.0 - gammajTerm);
							gammaTempNot[i] *= (1.0 - gammajTerm);
						}

						//System.out.println("gammajterm: "+gammajTerm);
					}

					//System.out.println("i="+i+", GammaTemp: "+gammaTemp[i]+", gammaTempNot: "+gammaTempNot[i]);
					gammaSum[i] += gammaTemp[i];
					gammaNotSum[i] += gammaTempNot[i];
				}

				//psChildren.printPowerSet();
				ps.getNext();
			} // done with the while loop over all children settings

			for (int i = 0; i < gamma.length; i++) {
				gamma[i] *= gammaSum[i];
				gammaNot[i] *= gammaNotSum[i];
				//System.out.println("For child "+k+",gammaSum "+gammaSum[i]+" gammaNotSum "+gammaNotSum[i]+" distance "+distance);
				//System.out.println("   and "+i+", gamma "+gamma[i]+" gammaNot "+gammaNot[i]);
				gammaSum[i] = 0.0;
				gammaNotSum[i] = 0.0;
			}
		}// done with product over children

		for (int i = 0; i < gamma.length; i++) {
			//System.out.println(" "+gamma[i]+" "+gammaNot[i]+"  "+(gamma[i]/(gamma[i]+gammaNot[i])));
			gamma[i] = logSafe(gamma[i] / (gamma[i] + gammaNot[i]));
		}

		//if(gamma.length > 0) System.out.println();
		return gamma;
	}


	// Compute the probability of the parent(s) given the kids.
	// Denominator is p(K+,K-), numerator is p(P_i,K+,K-)
	// This function is derived from heckerman90, equation 11,
	// but assumes that the D variables are *not* mutually independent
	// faster, debugged.
	/*
	  private void probParentsGivenKids(double[] parents, Vector children,
				      ProteinFunctionMarkovState psChildren,
				      GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx, double rate,
				      double distance)
	  {
	double[] negParents = new double[parents.length];
	for(int i = 0; i < parents.length; i++) {
	    parents[i] = 0.0;
	    negParents[i] = 0.0;
	}
	// Iterate through this setting of kids
	ProteinFunctionMarkovState parentsSet = new ProteinFunctionMarkovState(parents.length, pfx.maxFunctions());
	while(parentsSet.hasNext()) {
	    parentsSet.getNext();
	    //parentsSet.printPowerSet();
	    double parentProd = 1.0; // product for a particular i
	    // Product over all positive power set and all negative kids
	    for(int j = 0; j < psChildren.length(); j++) {
		parentProd *= probChildGivenParent(psChildren.elementAt(j),
						  j, parentsSet, pfx,
						  distance, rate);
		//System.out.println("delta: "+pfx.getDelta(k, j)+", scale_param: "+scale_param+", distances: "+distances[childIndex]);
		//System.out.println("     d: "+pfx.getDelta(i,j)+", dist: "+distances[childIndex]);

	    } // summation over all children
	    // add onto the appropriate parent
	    parentProd *= getPrior(parentsSet);
	    for(int j = 0; j < parents.length; j++) {
		//System.out.println("   adding "+(negParentsProd[j])+", "+(parentsProd[j]));
		if(parentsSet.elementAt(j) == 1) {
		    parents[j] += parentProd;
		} else {
		    negParents[j] += parentProd;
		}
	    }
	}
	// Normalize the parent messages
	for(int i = 0; i < parents.length; i++) {
	    parents[i] = parents[i]/ (parents[i]+negParents[i]);
	}
	}*/


	// Compute the probability of the parent(s) given the kids.
	// Denominator is p(K+,K-), numerator is p(P_i,K+,K-)
	// This function is derived from heckerman90, equation 11,
	// but assumes that the D variables are *not* mutually independent
	// faster, not debugged.
	private void probParentsGivenKidsME(double[] parents, Vector<Node> children,
	                                    ProteinFunctionMarkovState psChildren,
	                                    GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx, double rate,
	                                    double distance) {
		double[] negParents = new double[parents.length];

		for (int i = 0; i < parents.length; i++) {
			parents[i] = 0.0;
			negParents[i] = 0.0;
		}

		// Iterate through this setting of kids
		ProteinFunctionMarkovState parentsSet = new ProteinFunctionMarkovState(parents.length, pfx.maxFunctions());
		double parentProd = 1.0;

		while (parentsSet.hasNext()) {
			parentsSet.getNext();
			//parentsSet.printPowerSet();
			// Product over all positive power set and all negative kids
			parentProd = pfx.getExpProb(parentsSet.setIndex(),
			                            psChildren.setIndex(),
			                            rate * distance);
			//System.out.println("parents set index: "+parentsSet.setIndex());
			//System.out.println("Parent prod: "+parentProd);
			// add onto the appropriate parent
			parentProd *= getPrior(parentsSet);

			for (int j = 0; j < parents.length; j++) {
				if (parentsSet.elementAt(j) == 1) {
					parents[j] += parentProd;
					//System.out.println(" adding "+(parentProd)+" to "+j);
				}
				else {
					negParents[j] += parentProd;
					//System.out.println(" adding neg "+(parentProd)+" to "+j);
				}
			}
		}

		// Normalize the parent messages
		for (int i = 0; i < parents.length; i++) {
			//System.out.println("Parents: "+parents[i]);
			parents[i] = parents[i] / (parents[i] + negParents[i]);
		}
	}

	/**
	 * @param messages
	 * @param children
	 * @return
	 */
	private boolean messagesContainsAllChildren(Vector<Node> messages,
	    Vector<Node> children) {
		for (int i = 0; i < children.size(); i++) {
			if (!messages.contains(children.elementAt(i)))
				return false;
		}

		return true;
	}

	/**
	 * @param gamma
	 * @param delta
	 */
	private double[] gammaDeltaBinaryLog(double[] gamma, double[] delta) {
		double[] posterior = new double[delta.length];

		if (gamma.length != delta.length) { // gamma was only a prior
			for (int i = 0; i < delta.length; i++) {
				posterior[i] = (Math.exp(delta[i]))
				               / ((Math.exp(delta[i])) + (1 - Math.exp(delta[i])));
			}
		}
		else
			if (delta[0] >= 0) { // in the case of the root node
				for (int i = 0; i < gamma.length; i++) {
					posterior[i] = (Math.exp(gamma[i]))
					               / ((Math.exp(gamma[i])) + (1 - Math.exp(gamma[i])));
				}
			}
			else {
				for (int i = 0; i < gamma.length && i < delta.length; i++) {
					posterior[i] = Math.exp(gamma[i] + delta[i]) /
					               (Math.exp(delta[i] + gamma[i]) +
					                ((1 - Math.exp(gamma[i])) * (1 - Math.exp(delta[i]))));
				}
			}

		return posterior;
	}

	public Hashtable<Node, double[]> getGammas() {
		return gammasFinal;
	}

	public Hashtable<Node, double[]> getDeltas() {
		return deltasFinal;
	}

	////////////////////////////////////////////////////
	// Log Likelihood
	////////////////////////////////////////////////////

	/* calculates the log likelihood of the model given the
	 * sets of parameters
	 */

	// yunesj: unused and subfunctions have a bug
	//    public double getLogLikelihood(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx, Hashtable<String, Double> scaleParams,
	//				   double rate)
	//    {
	//	//Propagate likelihood messages up tree
	//	double[] prior = new double[1];
	//	prior[0] = 1;
	//	Vector<Node> messages = new Vector<Node>();
	//	Hashtable<Node, double[]> messageList = new Hashtable<Node, double[]>();
	//	int messageLength = 1;
	//	double logLikelihood = 0.0;
	//	//Get children with evidence first.
	//	//since the tree is constructed depth first,
	//	//this should be the correct way to go through the tree.
	//	for(int node = tree.size()-1; node >= 0; node--) {
	//	    // Get the name of the node
	//	    Node n = tree.elementAt(node);
	//	    if(n == null) {
	//		System.out.println("Can't find node in getLogLikelihood");
	//	    }
	//	    //Case 1: Leaf with evidence
	//	    //sample directly from those local probabilities
	//	    if(n.hasLocalProbabilities() && n.isLeaf()) {
	//		double[] evidenceProbs = n.getLocalProbabilities();
	//		for(int i = 0; i < evidenceProbs.length; i++){
	//		    //System.out.println(evidenceProbs[i]);
	//		    if(evidenceProbs[i] > 0) {
	//			evidenceProbs[i] = logSafe(evidenceProbs[i]);
	//		    }
	//		}
	//		messageList.put(n, evidenceProbs);
	//		messages.add(n);
	//		if(evidenceProbs.length > messageLength)
	//		    messageLength = evidenceProbs.length;
	//		//System.out.println("Added message for: "+n.getNodeID());
	//	    } // end of leaf with evidence
	//	    else if(n.isLeaf()){
	//		messageList.put(n, prior);
	//		messages.add(n);
	//		//System.out.println("Added message for; "+n.getNodeID());
	//	    } // end of leaf without evidence
	//	    else {
	//	    } // end of internal node
	//	}// end of initial leaf search
	//
	//	boolean done = false;
	//	while(!messages.isEmpty()) {
	//	    done = false;
	//	    int count = 0;
	//	    while(!done && messages.size()>0) {
	//		count = (count+1)%messages.size();
	//		Node parent = messages.elementAt(count).getParent();
	//		if(parent != null) {
	//		    Vector<Node> children = (Vector<Node>)parent.getChildren();
	//		    if(messagesContainsAllChildren(messages, children)) {
	//			double inUse = ((Double)scaleParams.get(SPECIATION_PARAMS)).doubleValue();
	//			if(parent.hasDuplication())
	//			    inUse = ((Double)scaleParams.get(DUPLICATION_PARAMS)).doubleValue();
	//			//System.out.println("Working on "+parent.getNodeID());
	//			if(!parent.isRoot()) {
	//			    // Something not right about two message lists
	//			    double[] message = singleLikelihood(messageList,
	//								messageList,
	//								children,
	//								pfx, inUse,
	//								messageLength,
	//								rate);
	//			    if(message.length == 1) message = prior;
	//			    messageList.put(parent, message);
	//			    messages.add(parent);
	//			}
	//			else { // For the root
	//			    logLikelihood = rootLikelihood(messageList,
	//							   children, pfx,
	//							   inUse,
	//							   messageLength,
	//							   rate);
	//			    messages.clear();
	//			}
	//			messages.removeAll(children);
	//		    } // finished that message
	//		}
	//	    }
	//	}
	//	System.out.println("Log likelihood of model: "+logLikelihood);
	//	return logLikelihood;
	//    }

	/* propagates log likelihood calculation up from children messages
	 * to parents, in the usual way
	 */

	// yunesj: these two functions have a bug (psChildren is always null)
	// and are unused
	//    private double[] singleLikelihood(Hashtable<Node, double[]> message1List,
	//				      Hashtable<Node, double[]> message0List, Vector<Node> children,
	//				      GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx, double scale_param,
	//				      int messageLength, double rate)
	//    {
	//	int psLength = messageLength;
	//	messageLength = (int)Math.round(Math.pow(2, messageLength)) - 1;
	//	double[] message = new double[messageLength];
	//	//Iterate through this setting of kids
	//	for(int i = 0; i < messageLength; i++) {
	//	    message[i] = 0.0;
	//	}
	//	boolean success = false;
	//	ProteinFunctionMarkovState psChildren = null;
	//	ProteinFunctionMarkovState psParent = new ProteinFunctionMarkovState(psLength, pfx.maxFunctions());
	//	for(int j = 0; j < children.size(); j++) {
	//	    double[] childMessage =
	//		message1List.get(children.elementAt(j));
	//	    if(childMessage != null
	//	       && childMessage.length >= 1 && childMessage[0] <= 0.0) {
	//		//if(childMessage == null)
	//		//    childMessage = new ProteinFunctionMarkovState(childMessage.length, maxFunctions);
	//		double distance =
	//		    children.elementAt(j).getParentDistance();
	//		int parentPSCount = 0;
	//		psParent.Reset();
	//		while(psParent.hasNext()) {
	//		    psParent.getNextNonZeros();
	//		    int psCount = 0;
	//		    double parentSum = 1.0;
	//		    psChildren.Reset();
	//		    while(psChildren.hasNext()) {
	//			psChildren.getNextNonZeros();
	//			double parentsProd = childMessage[psCount];
	//			for(int k = 0; k < psChildren.length(); k++) {
	//			  parentsProd +=
	//			  logSafe(probChildGivenParent(psChildren.elementAt(k),
	//						       k,
	//						       psParent,
	//						       pfx,
	//						       distance,
	//						       rate));
	//			}
	//			if(parentSum > 0) parentSum = parentsProd;
	//			else parentSum = GenericMathFunctions.getLogSum(parentSum,
	//							    parentsProd);
	//			psCount++;
	//		    }
	//		    message[parentPSCount] += parentSum;
	//		    parentPSCount++;
	//		}
	//		success = true;
	//	    }
	//	}
	//	//System.out.println("message at new node ");
	//	//for(int i = 0; i < message.length; i++)
	//	//    System.out.println("  "+message[i]);
	//	if(success) return message;
	//	// otherwise, none of children had evidence; parent also meaningless
	//	double[] prior = new double[1];
	//	prior[0] = 1.0;
	//	return prior;
	//    }
	//
	//    /* propagates messages below root up to root, then marginalizes
	//     * out the root node too, according to non-marginally independent
	//     * probabilites
	//     */
	//    private double rootLikelihood(Hashtable<Node, double[]> messageList, Vector<Node> children,
	//				  GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx, double scale_param,
	//				  int messageLength, double rate)
	//    {
	//	int psLength = messageLength;
	//	messageLength = (int)Math.round(Math.pow(2, messageLength)) - 1;
	//	double message = 0.0;
	//	//double[] message = new double[messageLength];
	//	//Iterate through this setting of kids
	//	//for(int i = 0; i < messageLength; i++) {
	//	//    message[i] = 0.0;
	//	//}
	//	boolean success = false;
	//	ProteinFunctionMarkovState psChildren = null;
	//	ProteinFunctionMarkovState psParent = new ProteinFunctionMarkovState(psLength, pfx.maxFunctions());
	//	for(int j = 0; j < children.size(); j++) {
	//	    double[] childMessage =
	//		messageList.get(children.elementAt(j));
	//	    if(childMessage != null
	//	       && childMessage.length >= 1 && childMessage[0] <= 0.0) {
	//		//if(childMessage == null)
	//		//    childMessage = new ProteinFunctionMarkovState(childMessage.length, maxFunctions);
	//		double distance =
	//		    children.elementAt(j).getParentDistance();
	//		int parentPSCount = 0;
	//		psParent.Reset();
	//		while(psParent.hasNext()) {
	//		    psParent.getNextNonZeros();
	//		    int psCount = 0;
	//		    double parentSum = 1.0;
	//		    psChildren.Reset();
	//		    while(psChildren.hasNext()) {
	//			psChildren.getNextNonZeros();
	//			double parentsProd = childMessage[psCount];
	//			for(int k = 0; k < psChildren.length(); k++) {
	//			  parentsProd +=
	//			  logSafe(probChildGivenParent(psChildren.elementAt(k),
	//						       k,
	//						       psParent,
	//						       pfx,
	//						       distance,
	//						       rate));
	//			}
	//			if(parentSum > 0) parentSum = parentsProd;
	//			else parentSum = GenericMathFunctions.getLogSum(parentSum,
	//							    parentsProd);
	//			psCount++;
	//		    }
	//		    message += parentSum;
	//		    parentPSCount++;
	//		}
	//		success = true;
	//	    }
	//	}
	//	//System.out.println("message at new node ");
	//	//for(int i = 0; i < message.length; i++)
	//	//    System.out.println("  "+message[i]);
	//	if(success) return message;
	//	// otherwise, none of children had evidence; parent also meaningless
	//	System.out.println("In LogLikelihood: didn't find any evidence");
	//	return(0.0);
	//    }


	/////////////////////////////////////////////////////////
	// Expectation Maximization (at least the M step)
	/////////////////////////////////////////////////////////

	public boolean maximizationStep(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx0,
	                                Hashtable<String, Double> scaleParams,
	                                Hashtable<Node, double[]> posteriors, double rho,
	                                double cutoff,
	                                int iteration) {
		//double logLikelihood = getLogLikelihood(pfx0, scaleParams, rate);

		// Perform gradent ascent for M step of EM
		int i = 0;
		double delta = 1.0;

		while (delta > cutoff && i < 1) {
			// do update
			delta = singleCompleteMaxUpdateME(pfx0, scaleParams,
			                                  posteriors, rho, iteration);
			System.out.println("M-step gradient iteration: " + iteration
			                   + ", delta: " + delta);
			//+", likelihood: " +logLikelihood);
			i++;
		}

		if (delta < cutoff)
			return true;

		return false;
		//System.out.println("Gradient iteration: "+i
		//		+", delta: "+delta+", likelihood: "+logLikelihood);
	}

	// yunesj: unused
	//    public double singleCompleteMaxUpdate(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfxT,
	//					  Hashtable scaleParams,
	//					  Hashtable posteriors, double rho,
	//					  int iteration)
	//    {
	//	double[][] sumGradients =
	//	    new double[pfxT.getRows()][pfxT.getColumns()];
	//	Vector<Node> evidence = findTreeWithEvidence();
	//	Hashtable<String, Double> scaleGradients = new Hashtable<String, Double>();
	//	Enumeration scaleNames = scaleParams.keys();
	//	while(scaleNames.hasMoreElements()) {
	//	    String scaleName = (String)scaleNames.nextElement();
	//	    scaleGradients.put(scaleName, new Double(0.0));
	//	}
	//
	//	double deltaScale = 0.0;
	//	for(int i = 0; i < evidence.size(); i++) {
	//	    Node n = evidence.elementAt(i);
	//	    if(!n.isRoot()) {
	//		//System.out.println("working on node: "+n.getNodeID());
	//		double inUse =
	//		    ((Double)scaleParams.get(SPECIATION_PARAMS)).doubleValue();
	//		if(n.getParent() != null && n.getParent().hasDuplication())
	//		   inUse =
	//		   ((Double)scaleParams.get(DUPLICATION_PARAMS)).doubleValue();
	//		double scaleGradient;
	//		scaleGradient =
	//		    maximizeSingleParentChild(pfxT, inUse,
	//					      (double[])posteriors.get(n),
	//					      (double[])posteriors
	//					      .get(n.getParent()),
	//					      n.getParentDistance(),
	//					      sumGradients,
	//					      iteration);
	//		scaleGradient = scaleGradient*rho;
	//		if(n.getParent() != null && n.getParent().hasDuplication()) {
	//		    double dupGradient =
	//        	   ((Double)scaleParams.get(DUPLICATION_PARAMS)).doubleValue();
	//		    if(dupGradient + scaleGradient > 0) {
	//			scaleGradients.put(DUPLICATION_PARAMS,
	//					   new Double(dupGradient
	//						      + scaleGradient));
	//			deltaScale += scaleGradient;
	//		    }
	//		} else {
	//		    double speGradient =
	//        	   ((Double)scaleParams.get(SPECIATION_PARAMS)).doubleValue();
	//		    if(speGradient + scaleGradient > 0) {
	//			scaleGradients.put(SPECIATION_PARAMS,
	//					   new Double(speGradient
	//						      + scaleGradient));
	//			deltaScale += scaleGradient;
	//		    }
	//		}
	//	    }
	//	}
	//	System.out.println("Sum gradients ");
	//	// Update pfx matrix
	//	double delta = 0.0;
	//	for(int m = 0; m < sumGradients.length; m++) {
	//	    for(int n = 0; n < sumGradients[m].length; n++) {
	//		double gradient = (rho * sumGradients[m][n]);
	//		pfxT.addDeltaPositive(m, n, gradient);
	//		delta += Math.abs(gradient);
	//		//System.out.println("m,n "+m+","+n+": "+gradient);
	//	    }
	//	}
	//	// Update scale params
	//	scaleNames = scaleParams.keys();
	//	while(scaleNames.hasMoreElements()) {
	//	    String scaleName = (String)scaleNames.nextElement();
	//	    Double scaleGrad = scaleGradients.get(scaleName);
	//	    pfxT.setScale(scaleName, scaleGrad.doubleValue());
	//            scaleParams = pfxT.getScale();
	//	}
	//	//System.out.println("Delta: "+delta+", deltascale: "+deltaScale);
	//	return (Math.abs(deltaScale) + delta);
	//    }

	public double singleCompleteMaxUpdateME(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfxT,
	                                        Hashtable<String, Double> scaleParams,
	                                        Hashtable<Node, double[]> posteriors, double rho,
	                                        int iteration) {
		double[][] sumGradients =
		  new double[pfxT.getRows()][pfxT.getColumns()];
		double[] alphaGradients =
		  new double[pfxT.getColumns()];
		double[][] phiGradients =
		  new double[pfxT.getRows()][pfxT.getColumns()];
		double[] phi0Gradient =
		  new double[pfxT.getColumns()];
		LAPACKMatrixOperationsWrapper phiGradientsLAPACKMatrixOperationsWrapper = new LAPACKMatrixOperationsWrapper(pfxT.getSizeExpMatrix(),
		    pfxT.getSizeExpMatrix());
		precomputePhiGradient(pfxT, phiGradients, phiGradientsLAPACKMatrixOperationsWrapper,
		                      phi0Gradient);
		//printMatrix(phiGradients);

		Vector<Node> evidence = findTreeWithEvidence();
		Hashtable<String, Double> scaleGradients = new Hashtable<String, Double>();
		Hashtable<String, Integer> scaleCounts = new Hashtable<String, Integer>();
		Enumeration<String> scaleNames = scaleParams.keys();

		while (scaleNames.hasMoreElements()) {
			String scaleName = scaleNames.nextElement();
			scaleGradients.put(scaleName, new Double(0.0));
			scaleCounts.put(scaleName, new Integer(0));
		}

		double deltaScale = 0.0;
		//double epsilon = 1.0E-18;
		double epsilon = 0.01;

		for (int m = 0; m < sumGradients.length; m++) {
			for (int n = 0; n < sumGradients[m].length; n++) {
				sumGradients[m][n] = 0.0;
			}
		}

		for (int i = 0; i < evidence.size(); i++) {
			Node n = evidence.elementAt(i);

			if (!n.isRoot()) {
				//System.out.println("working on node: "+n.getNodeID());
				double inUse =
				  scaleParams.get(SPECIATION_PARAMS).doubleValue();

				if (n.getParent() != null && n.getParent().hasDuplication())
					inUse =
					  scaleParams.get(DUPLICATION_PARAMS).doubleValue();

				double scaleGradient;
				scaleGradient =
				  maximizeSingleParentChildME(pfxT, inUse,
				                              (double[])posteriors
				                              .get(n.getParent()),
				                              (double[])posteriors.get(n),
				                              n.getParentDistance(),
				                              sumGradients,
				                              alphaGradients,
				                              phiGradients,
				                              phiGradientsLAPACKMatrixOperationsWrapper,
				                              phi0Gradient,
				                              iteration);
				// scale times step size
				scaleGradient = scaleGradient * rho;

				if (n.getParent() != null && n.getParent().hasDuplication()) {
					double dupGradient =
					  scaleGradients.get(DUPLICATION_PARAMS).doubleValue();
					scaleGradients.put(DUPLICATION_PARAMS,
					                   new Double(dupGradient
					                              + scaleGradient));
					scaleCounts.put(DUPLICATION_PARAMS,
					                new Integer(1 + scaleCounts.get(DUPLICATION_PARAMS).intValue()));
				}
				else {
					double speGradient =
					  scaleGradients.get(SPECIATION_PARAMS).doubleValue();
					scaleGradients.put(SPECIATION_PARAMS,
					                   new Double(speGradient
					                              + scaleGradient));
					scaleCounts.put(SPECIATION_PARAMS,
					                new Integer(1 + scaleCounts.get(SPECIATION_PARAMS).intValue()));
				}
			}
		}

		// Update pfx matrix
		double delta = 0.0;
		elementwiseDivide(sumGradients, phiGradients);
		//printMatrix(sumGradients);
		// do the diagonal elements separately (project on simplex)
		/*
		double total = 0.0;
		for(int m = 0; m < sumGradients.length; m++) {
		  double gradient = (rho * sumGradients[m][m]);
		  double newVal = pfxT.getDelta(m,m)+gradient;
		  newVal = (newVal>epsilon) ? newVal : epsilon;
		  total += newVal;
		}
		// Projected part of gradient ascent
		for(int m = 0; m < sumGradients.length; m++) {
		  if(total <= (double)pfxT.getRows())
		total = 1;
		  if(total >= (double)pfxT.getRows())
		total /= (double)pfxT.getRows();
		  // diagonal times step size
		  double gradient = (rho * sumGradients[m][m]);
		  double oldVal = pfxT.getDelta(m,m);
		  double newVal = oldVal+gradient;
		  newVal = (newVal>epsilon) ? (newVal/total) : epsilon;
		  pfxT.setDelta(m,m, newVal);
		  delta += Math.abs(newVal-oldVal);
		}
		    double deltadiag = delta;
		 */
		double deltadiag = 0;
		double[] total = new double[sumGradients.length];

		for (int m = 0; m < sumGradients.length; m++) {
			total[m] = 0.0;

			for (int n = 0; n < sumGradients[m].length; n++) {
				double gradient = (rho * sumGradients[m][n]);
				double newVal = pfxT.getDelta(m, n) + gradient;
				newVal = (newVal > epsilon) ? newVal : epsilon;
				total[m] += newVal;
			}
		}

		// Projected part of gradient ascent
		for (int m = 0; m < sumGradients.length; m++) {
			if (total[m] <= (double)pfxT.getRows())
				total[m] = 1;

			if (total[m] > (double)pfxT.getRows())
				total[m] /= (double)pfxT.getRows();

			for (int n = 0; n < sumGradients[m].length; n++) {
				double gradient = (rho * sumGradients[m][n]);
				double oldVal = pfxT.getDelta(m, n);
				double newVal = oldVal + gradient;
				newVal = (newVal > epsilon) ? (newVal / total[m]) : epsilon;
				pfxT.setDelta(m, n, newVal);
				delta += Math.abs(newVal - oldVal);

				if (n == m)
					deltadiag += Math.abs(newVal - oldVal);
			}
		}

		/*for(int m = 0; m < sumGradients.length; m++) {
		  for(int n = 0; n < sumGradients[m].length; n++) {
		if(m != n) {
		    // Off-diagona times step size
		    double gradient = (rho * sumGradients[m][n]);
		    if(pfxT.addDeltaPositive(m, n, gradient))
			delta += Math.abs(gradient);
		    else {
			delta += pfxT.getDelta(m,n)-epsilon;
			pfxT.setDelta(m,n,epsilon);
		    }
		}
		//System.out.println("m,n "+m+","+n+": "+gradient);
		  }
		  }*/
		// Update alpha parameters
		double deltaAlpha = 0.0;
		//elementwiseDivide(alphaGradients, phi0Gradient);
		//printVector(alphaGradients);
		/*for(int m = 0; m < alphaGradients.length; m++) {
		  // alpha times step size
		  double gradient = (rho * alphaGradients[m]);
		  pfxT.addAlphaPositive(m, gradient);
		  deltaAlpha += Math.abs(gradient);
		  //System.out.println("m "+m+": "+gradient);
		  }*/
		double alphaTotal = 0.0;

		for (int m = 0; m < alphaGradients.length; m++) {
			double gradient = (rho * alphaGradients[m]);
			double newVal = pfxT.getAlpha(m) + gradient;
			newVal = (newVal > epsilon) ? newVal : epsilon;
			alphaTotal += newVal;
		}

		// Projected part of gradient ascent
		for (int m = 0; m < alphaGradients.length; m++) {
			if (alphaTotal <= (double)pfxT.getRows())
				alphaTotal = 1;

			if (alphaTotal > (double)pfxT.getRows())
				alphaTotal /= (double)pfxT.getRows();

			// diagonal times step size
			double gradient = (rho * alphaGradients[m]);
			double oldVal = pfxT.getAlpha(m);
			double newVal = oldVal + gradient;
			newVal = (newVal > epsilon) ? (newVal / alphaTotal) : epsilon;
			pfxT.setAlpha(m, newVal);
			deltaAlpha += Math.abs(newVal - oldVal);
		}

		// Update scale params
		scaleNames = scaleParams.keys();

		while (scaleNames.hasMoreElements()) {
			String scaleName = scaleNames.nextElement();
			int scaleCount = scaleCounts.get(scaleName).intValue();
			double scaleGrad = (scaleGradients.get(scaleName).doubleValue() / (double)scaleCount);

			if (scaleGrad + pfxT.getScale(scaleName) > epsilon) {
				pfxT.setScale(scaleName, scaleGrad + pfxT.getScale(scaleName));
				deltaScale += Math.abs(scaleGrad);
			}
			else {
				pfxT.setScale(scaleName, epsilon);
			}

			scaleParams = pfxT.getScale();
			//System.out.println("scaleGrad: "+scaleGrad);
			//System.out.println("scaleCount: "+scaleCount);
		}

		//printVector(alphaGradients);
		//printMatrix(sumGradients);
		System.out.println("Delta: " + delta + ", deltaAlpha: "
		                   + deltaAlpha);
		System.out.println("deltascale: " + deltaScale + ", Deltadiag: " + deltadiag);
		return (Math.abs(deltaScale) + deltaAlpha + delta);
	}

	public double maximizeSingleParentChild(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx,
	                                        double scaleParam,
	                                        double[] posteriorsParent,
	                                        double[] posteriorsChild,
	                                        double distance,
	                                        double[][] sumGradients,
	                                        int iteration) {
		// Length terms
		double scaleGradient = 0.0;

		for (int n = 0; n < sumGradients.length; n++) {
			for (int m = 0; m < sumGradients[n].length; m++) {
				if (m != n) {
					double gradientStep =
					  singleYMNGradientStep(m, n, distance,
					                        posteriorsParent,
					                        posteriorsChild,
					                        pfx,
					                        scaleParam);
					sumGradients[m][n] += gradientStep;
					sumGradients[n][m] += gradientStep; // for symmetry
				}
				else {
					double gradientStep =
					  singleYNNGradientStep(n, distance,
					                        posteriorsParent,
					                        posteriorsChild,
					                        pfx,
					                        scaleParam);
					sumGradients[n][n] += gradientStep;
				}
			}

			scaleGradient += singleScaleGradientStep(n, distance,
			                 posteriorsParent,
			                 posteriorsChild,
			                 pfx, scaleParam);
		}

		return scaleGradient;
	}

	@SuppressWarnings("unused")
	public double maximizeSingleParentChildME(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx,
	    double scaleParam,
	    double[] posteriorsParent,
	    double[] posteriorsChild,
	    double distance,
	    double[][] sumGradients,
	    double[] alphaGradients,
	    double[][] phiGradient,
	    LAPACKMatrixOperationsWrapper phiGradientLAPACKMatrixOperationsWrapper,
	    double[] phi0Gradient,
	    int iteration) {
		// assumes initialized sumGradients/alphaGradients
		if (sumGradients == null || alphaGradients == null) {
			System.out.print("Error: in maximizeSingleParentChildME, must ");
			System.out.print("initialize sumGradients and alphaGradients ");
			System.out.println("before call");
			System.exit(0);
		}

		// Iterate over power set of parents/children
		ProteinFunctionMarkovState parent = new ProteinFunctionMarkovState(posteriorsParent.length,
		                               pfx.maxFunctions());
		ProteinFunctionMarkovState child = new ProteinFunctionMarkovState(posteriorsParent.length,
		                              pfx.maxFunctions());

		double scaleGradient = 0.0;
		int diffIndex = 0;

		//System.out.println("One node");

		//printMatrix(sumGradients);
		double prod = 1;
		double noParent = 1.0;

		for (int i = 0; i < posteriorsParent.length; i++)
			noParent *= (1.0 - posteriorsParent[i]);

		double noChild = 1.0;

		for (int j = 0; j < posteriorsChild.length; j++)
			noChild *= (1.0 - posteriorsChild[j]);

		double checksum = 0.0;

		while (child.hasNext()) {
			child.getNextNonZeros();
			prod = noParent;
			double currentChild = 1.0;

			if (child.getCurrentSum() == 1) {
				for (int j = 0; j < child.length(); j++) {
					if (child.elementAt(j) == 1)
						currentChild *= posteriorsChild[j];
					else
						currentChild *= (1.0 - posteriorsChild[j]);
				}

				checksum += currentChild;

				for (int j = 0; j < child.length(); j++) {
					if (child.elementAt(j) == 1) {
						alphaGradients[j] += (currentChild * noParent
						                      * scaleParam * distance);
						alphaGradients[j] -= (noChild * noParent
						                      * scaleParam * distance);
					}
				}
			}
		}

		//System.out.println("Checksum: "+(checksum+noChild));

		int childCount = 0;
		int parentCount = 1;
		child.Reset();

		while (parent.hasNext()) {
			parent.getNextNonZeros();
			double currentParent = 1.0;

			for (int i = 0; i < parent.length(); i++) {
				if (parent.elementAt(i) == 1)
					currentParent *= posteriorsParent[i];
				else
					currentParent *= (1.0 - posteriorsParent[i]);
			}

			double sameChild = 1.0;

			for (int j = 0; j < child.length(); j++) {
				if (child.elementAt(j) == 1)
					sameChild *= posteriorsChild[j];
				else
					sameChild *= (1.0 - posteriorsChild[j]);
			}

			childCount = 0;
			child.Reset();
			double identityProd = currentParent * sameChild * distance * scaleParam;

			while (child.hasNext()) {
				child.getNext();
				//parent.printPowerSet();
				//child.printPowerSet();
				double currentChild = 1.0;
				diffIndex = parent.offByOne(child);

				if (diffIndex != -1) {
					for (int j = 0; j < child.length(); j++) {
						if (child.elementAt(j) == 1)
							currentChild *= posteriorsChild[j];
						else
							currentChild *= (1.0 - posteriorsChild[j]);
					}

					prod = currentChild * currentParent
					       * distance * scaleParam;

					for (int i = 0; i < parent.length(); i++) {
						if (parent.elementAt(i) == 1) {
							for (int j = 0; j < child.length(); j++) {
								if (i == j && child.elementAt(j) == 0) {
									sumGradients[i][i] += prod;
									sumGradients[i][i] -= identityProd;
								}
								else
									if (i != j && child.elementAt(j) == 1
									    && parent.elementAt(j) == 0) {
										sumGradients[i][j] += prod;
										sumGradients[j][i] += prod;
										sumGradients[i][j] -= identityProd;
										sumGradients[j][i] -= identityProd;
										// Alpha updates
										alphaGradients[j] += prod;
										alphaGradients[j] -= identityProd;

									}
							}
						}
					}

					// Update scaleGradient
					scaleGradient += (currentChild * currentParent
					                  * distance
					                  * pfx.getQElement(parentCount, childCount));
					//System.out.println("Adding: "+(currentChild*currentParent*distance*pfx.getQElement(parentCount, childCount)));
				}

				childCount++;
			}

			parentCount++;
		}

		//System.out.println("Scale gradient before:"+scaleGradient);
		return scaleGradient;
	}

	public double maximizeSingleParentChildScale(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx,
	    double scaleParam,
	    double[] posteriorsParent,
	    double[] posteriorsChild,
	    double distance,
	    double[][] sumGradients) {
		// Length terms
		double scaleGradient = 0.0;

		for (int n = 0; n < sumGradients.length; n++) {
			scaleGradient += singleScaleGradientStep(n, distance,
			                 posteriorsParent,
			                 posteriorsChild,
			                 pfx, scaleParam);
		}

		return scaleGradient;
	}

	/* computes a gradient computation, putting the
	 * result into the gradient parameter,
	 * for a specific setting of parents/children.
	 */
	public double singleGradientStep(int m, int n, double distance,
	                                 double[] posteriorsParent,
	                                 double[] posteriorsChild,
	                                 GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject inUse,
	                                 double scaleParam) {
		double denomProd = 1.0;

		for (int k = 0; k < posteriorsParent.length; k++) {
			denomProd *= Math.pow((1.0 - (1.0 / inUse.getLength(k, n))),
			                      (distance * scaleParam
			                       * posteriorsParent[k]));
		}

		double lengthMN = inUse.getLength(m, n);
		// First term
		double gradient = distance * scaleParam
		                  * (1.0 / (lengthMN * lengthMN))
		                  * posteriorsParent[n];
		// Second term
		gradient -= scaleParam * distance * posteriorsChild[m]
		            * posteriorsParent[n] * (1.0 / lengthMN * lengthMN)
		            / (1.0 - denomProd);
		//System.out.println(", "+gradient);
		return gradient;
	}

	/* computes denominators of gradient computation,
	 * putting the result into the gradient parameter,
	 * for a specific setting of parents/children.
	 */
	public double singleYMNGradientStep(int m, int n, double distance,
	                                    double[] posteriorsParent,
	                                    double[] posteriorsChild,
	                                    GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject inUse,
	                                    double scaleParam) {
		double prodTerm = singleProdTerm(n, posteriorsParent, inUse);
		double denomAProd = singleDenomA(n, distance, posteriorsParent,
		                                 inUse, scaleParam, prodTerm);
		double denomBProd = singleDenomB(n, distance, posteriorsParent,
		                                 inUse, scaleParam, prodTerm);
		double gradient = posteriorsChild[n]
		                  * ((1.0 - posteriorsParent[n])
		                     * posteriorsParent[m]
		                     * (prodTerm * (Math.exp(-distance * scaleParam) - 1.0)))
		                  / denomAProd;
		//System.out.println("posterior child: "+posteriorsChild[n]+", posterior parent: "+posteriorsParent[n]);
		//System.out.println("Gradient part "+gradient);
		gradient = gradient + ((1.0 - posteriorsChild[n])
		                       * (1.0 - posteriorsParent[n])
		                       * posteriorsParent[m]
		                       * (1.0 - Math.exp(-distance * scaleParam))
		                       * prodTerm / denomBProd);
		//System.out.println("Gradient full "+gradient);
		return gradient;
	}

	/* computes denominators of gradient computation,
	 * putting the result into the gradient parameter,
	 * for a specific setting of parents/children.
	 */
	public double singleScaleGradientStep(int n, double distance,
	                                      double[] posteriorsParent,
	                                      double[] posteriorsChild,
	                                      GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject inUse,
	                                      double scaleParam) {
		double prodTerm = singleProdTerm(n, posteriorsParent, inUse);
		double denomAProd = singleDenomA(n, distance, posteriorsParent,
		                                 inUse, scaleParam, prodTerm);
		double denomBProd = singleDenomB(n, distance, posteriorsParent,
		                                 inUse, scaleParam, prodTerm);
		double gradient = posteriorsChild[n]
		                  * (((1.0 - posteriorsParent[n]) * distance
		                      * prodTerm * Math.exp(-distance * scaleParam))
		                     - (posteriorsParent[n] * Math.exp(-inUse.getDelta(n, n))
		                        * distance * Math.exp(-distance * scaleParam)))
		                  / denomAProd;
		gradient = gradient
		           + ((1.0 - posteriorsChild[n])
		              * (((1.0 - posteriorsParent[n]) * (-distance)
		                  * prodTerm * Math.exp(-distance * scaleParam))
		                 + (posteriorsParent[n] * Math.exp(-inUse.getDelta(n, n))
		                    * distance * Math.exp(-distance * scaleParam)))
		              / denomBProd);
		//System.out.println(", "+gradient);
		return gradient;
	}

	/* computes denominators of gradient computation,
	 * putting the result into the gradient parameter,
	 * for a specific setting of parents/children.
	 */
	/* computes denominators of gradient computation,
	 * putting the result into the gradient parameter,
	 * for a specific setting of parents/children.
	 */
	public double singleYNNGradientStep(int n, double distance,
	                                    double[] posteriorsParent,
	                                    double[] posteriorsChild,
	                                    GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject inUse,
	                                    double scaleParam) {
		double prodTerm = singleProdTerm(n, posteriorsParent, inUse);
		double denomAProd = singleDenomA(n, distance, posteriorsParent,
		                                 inUse, scaleParam, prodTerm);
		double denomBProd = singleDenomB(n, distance, posteriorsParent,
		                                 inUse, scaleParam, prodTerm);
		double gradient = posteriorsChild[n]
		                  * (posteriorsParent[n]
		                     * (Math.exp(-inUse.getDelta(n, n))
		                        * (1.0 - Math.exp(-distance * scaleParam))))
		                  / denomAProd;
		gradient = gradient + ((1.0 - posteriorsChild[n])
		                       * posteriorsParent[n]
		                       * (Math.exp(-distance * scaleParam) - 1.0)
		                       * Math.exp(-inUse.getDelta(n, n)) / denomBProd);
		//System.out.println(", "+gradient);
		return gradient;
	}

	private double singleProdTerm(int n, double[] posteriorsParent,
	                              GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject inUse) {
		double prodTerm = 1.0;

		for (int k = 0; k < posteriorsParent.length; k++) {
			if (n != k)
				prodTerm *= Math.exp(-inUse.getDelta(n, k)
				                     * posteriorsParent[k]);
		}

		return prodTerm;
	}

	private double singleDenomA(int n, double distance,
	                            double[] posteriorsParent,
	                            GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject inUse, double scaleParam,
	                            double prodTerm) {
		double denomAProd = posteriorsParent[n] *
		                    (1.0 - (Math.exp(-inUse.getDelta(n, n))
		                            * (1.0 - Math.exp(-distance * scaleParam))));
		denomAProd = ((1.0 - posteriorsParent[n])
		              * (1.0 - Math.exp(-distance * scaleParam))
		              * prodTerm) + denomAProd;
		return denomAProd;
	}

	private double singleDenomB(int n, double distance,
	                            double[] posteriorsParent,
	                            GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject inUse, double scaleParam,
	                            double prodTerm) {
		double denomBProd = posteriorsParent[n] *
		                    Math.exp(-inUse.getDelta(n, n))
		                    * (1.0 - Math.exp(-distance * scaleParam));
		denomBProd = ((1.0 - posteriorsParent[n])
		              * (1.0 - ((1.0 - Math.exp(-distance * scaleParam))
		                        * prodTerm))) + denomBProd;
		return denomBProd;
	}

	/* returns the vector of tree nodes that have evidence or
	 * have a descendent leaf with evidence */
	public Vector<Node> findTreeWithEvidence() {
		Vector<Node> evidence = new Vector<Node>();
		Vector<Node> messages = new Vector<Node>();

		for (int node = tree.size() - 1; node >= 0; node--) {
			// Get the name of the node
			Node n = tree.elementAt(node);

			//Case 1: Leaf with evidence
			if (n.hasLocalProbabilities() && n.isLeaf()) {
				messages.add(n);
				evidence.add(n);
			} // end of leaf with evidence
		}// end of initial leaf search

		while (!messages.isEmpty()) {
			Node parent = messages.elementAt(0).getParent();

			if (parent != null
			    && !parent.isRoot()
			    && !evidence.contains(parent)) {
				evidence.add(parent);
				messages.add(parent);
			}

			messages.remove(0);
		}

		return evidence;
	}


	@SuppressWarnings("unused")
	public void precomputePhiGradient(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx,
	                                  double[][] phiGradients,
	                                  LAPACKMatrixOperationsWrapper qGradients,
	                                  double[] phi0Gradient) {
		ProteinFunctionMarkovState parent = new ProteinFunctionMarkovState(pfx.getRows(), pfx.maxFunctions());
		ProteinFunctionMarkovState child = new ProteinFunctionMarkovState(pfx.getRows(), pfx.maxFunctions());

		if (phiGradients == null || phi0Gradient == null || qGradients == null) {
			System.out.print("Error: in precomputeQGradient, must ");
			System.out.print("initialize phiGradients and phi0Gradient ");
			System.out.println("before call");
			System.exit(0);
		}

		for (int i = 0; i < parent.length(); i++) {
			phi0Gradient[i] = 0.0;

			for (int j = 0; j < child.length(); j++) {
				phiGradients[i][j] = 0.0;
			}
		}

		//qGradients.zero();
		int count = 0;
		int diffIndex = 0;

		while (child.hasNext()) {
			child.getNextNonZeros();

			// cycle through alphas to find gradients
			if (child.getCurrentSum() == 1) {
				for (int j = 0; j < child.length(); j++) {
					if (child.elementAt(j) == 1 && pfx.getAlpha(j) != 0) {
						phi0Gradient[j] += 1;
						//qGradients.add(0, count, 1);
					}
				}
			}

			count++;
		}

		child.Reset();
		int parentCount = 1;
		int childCount = 0;

		while (parent.hasNext()) {
			parent.getNextNonZeros();
			child.Reset();
			childCount = 0;

			while (child.hasNext()) {
				// important: since the left column is not included
				// in the counts, don't include it in the scaling gradient.
				child.getNext();
				//child.printPowerSet();
				diffIndex = parent.offByOne(child);

				if (diffIndex != -1 && diffIndex != parent.length()) {
					for (int i = 0; i < parent.length(); i++) {
						if (parent.elementAt(i) == 1) {
							// Cycle through, add gradient
							for (int j = 0; j < child.length(); j++) {
								if (i == j && child.elementAt(j) == 0) {
									phiGradients[i][i] += 1;
								}
								else
									if (i != j && child.elementAt(j) == 1
									    && parent.elementAt(j) == 0) {
										phiGradients[i][j] += 1;
										phiGradients[j][i] += 1;
										phi0Gradient[j] += 1;
										//qGradients.add(parentCount, childCount, 1);
									}
							}
						}
					}
				}

				childCount++;
			}

			parentCount++;
		}

		//System.out.println("PhiGradients:");
		//printMatrix(phiGradients);
		//qGradients.print();
	}


	/////////////////////////////////////////////////////////
	// Helper functions; move to Util
	/////////////////////////////////////////////////////////

	public double logSafe(double x) {
		if (x <= 0) {
			System.out.println("Warning: Underflow in log: log(" + x + ")");
			return -10000000;
		}
		else
			return Math.log(x);
	}

	public double choose(int n, int k) {
		double total = 1;
		double divisor = 1;
		k = Math.max(k, n - k);

		for (int i = n; i > k; i--) {
			total *= (double)i;
		}

		for (int i = 1; i <= (n - k); i++) {
			divisor *= (double)i;
		}

		total /= divisor;
		return total;
	}

	public double oneMinusLogSafe(double a) {
		return(logSafe(1 - Math.exp(a)));
	}


	public void elementwiseDivide(double[][] a, double[][] b) {
		if (a.length != b.length) {
			System.out.println("Error in elementwiseDivide: "
			                   + "rows do not match");
			return;
		}

		for (int i = 0; i < a.length; i++) {
			if (a[0].length != b[0].length) {
				System.out.println("Error in elementwiseDivide: "
				                   + "columns do not match");
				return;
			}

			for (int j = 0; j < a[0].length; j++) {
				a[i][j] = a[i][j] / b[i][j];
			}
		}
	}

	public void elementwiseDivide(double[] a, double[] b) {
		if (a.length != b.length) {
			System.out.println("Error in elementwiseDivide: "
			                   + "lengths do not match");
			return;
		}

		for (int i = 0; i < a.length; i++) {
			a[i] = a[i] / b[i];
		}
	}

	public void printVector(double[] v) {
		for (int i = 0; i < v.length; i++) {
			System.out.print(v[i] + " ");
		}

		System.out.println();
	}

	public void printMatrix(double[][] m) {
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[i].length; j++) {
				System.out.print(m[i][j] + " ");
			}

			System.out.println();
		}

		System.out.println();
	}

}

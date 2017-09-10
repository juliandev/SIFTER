/**
 * Data structure for reading and holding the parameters .fx file.
 * Mostly retrieving and setting values --- does some error checking.
 *
 * Copyright 2003-2005 Barbara Engelhardt (bee@cs.berkeley.edu)
 * @author Barbara Engelhardt
 */

package sifter.components;

import java.util.*;

import sifter.components.LAPACKMatrixOperationsWrapper;

import java.io.*;

public class GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject extends Object {

	public int rows;
	public int columns;
	
	public Vector<Vector<Parameter>> matrixRows;
	private Vector<Integer> goTerms;
	private double s;
	private Hashtable<String, Double> scaleParams;
	private Vector<Double> alphaParams;
	private LAPACKMatrixOperationsWrapper transRateMatrix;
	private LAPACKMatrixOperationsWrapper exlAPACKMatrixOperationsWrapperTemp;
	private LAPACKMatrixOperationsWrapper exlAPACKMatrixOperationsWrapperSummary;
	private ProteinFunctionMarkovState childPS;
	private int maxFunctions;
	private double currentScaling;
	private boolean summaryNotCurrent;

	private static double probEpsilon = 1.0E-18;

	/* initializer; sets rows and columns to all 0 */
	public GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject(int r, int c) {
		rows = r;
		columns = c;
		matrixRows = new Vector<Vector<Parameter>>();

		for (int i = 0; i < rows; i++) {
			Vector<Parameter> col = new Vector<Parameter>();

			for (int j = 0; j < columns; j++) {
				col.add(new Parameter(0.0));
			}

			matrixRows.add(i, new Vector<Parameter>(columns));
		}

		matrixRows.trimToSize();
		maxFunctions = matrixRows.size();
		goTerms = null;
		s = -1.0;
		scaleParams = null;
		transRateMatrix = null;
		exlAPACKMatrixOperationsWrapperTemp = null;
		exlAPACKMatrixOperationsWrapperSummary = null;
		childPS = null;
	}

	@SuppressWarnings("unused")
	public GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject(String infile) {
		//System.out.println("GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject(" + infile + ")");
		String str = null;
		int lineno = 0;
		BufferedReader fin = null;
		Double d = null;
		s = -1.0;
		transRateMatrix = null;
		exlAPACKMatrixOperationsWrapperTemp = null;
		exlAPACKMatrixOperationsWrapperSummary = null;
		childPS = null;
		matrixRows = new Vector<Vector<Parameter>>();

		if (goTerms == null)
			goTerms = new Vector<Integer>();

		try {
			fin = new BufferedReader(new FileReader(infile));

			try {
				while ((str = fin.readLine()) != null) {
					int pos = 0;

					if (str.indexOf("#") == 0) {
						str = "";
					}
					else
						if (str.indexOf("#") > 0) {
							str = str.substring(0, str.indexOf("#"));
						}

					if (str.length() > 0 && str.indexOf('\t') != 0) {
						if (str.indexOf('\t') > -1) {
							Vector<Parameter> singleRow = new Vector<Parameter>();
							str = str.substring(str.indexOf('\t') + 1,
							                    str.length());

							while (str != null && str.length() > 0) {
								if (str.indexOf('\t') > -1) {
									d = new Double(str.substring(0,
									                             str.indexOf('\t')));
									str = str.substring(str.indexOf('\t') + 1,
									                    str.length());
								}
								else {
									d = new Double(str);
									str = null;
								}

								singleRow.add(new Parameter(d.doubleValue()));
								pos++;
							}
							
							matrixRows.add(singleRow);
							columns = singleRow.size();
						}
					}
					else
						if (str.length() > 0) {
							str = str.substring(1, str.length());
							Integer index = null;

							while (str != null && str.length() > 0) {
								if (str.indexOf('\t') >= 0) {
									index = (new Integer(str.substring(0,
									                                   str.indexOf('\t'))));
									str = str.substring(str.indexOf('\t') + 1,
									                    str.length());
								}
								else {
									index = new Integer(str);
									str = null;
								}

								goTerms.add(index);
								pos++;
							}
						}

					lineno++;
				}
			}
			catch (Exception e) {
				System.out.println("Error in pFunTransMatrix: " + str
				                   + " " + lineno + " " + e.getMessage());
				System.exit(1);
			}

			fin.close();
		}
		catch (Exception ioe) {
			System.err.println("pFunTransMatrix: " + infile + " " +
			                   ioe.getMessage());
			System.exit(1);
		}

		rows = matrixRows.size();
		matrixRows.trimToSize();
		goTerms.trimToSize();
		maxFunctions = matrixRows.size();
	}

	public int getRows() {
		return rows;
	}

	public int getColumns() {
		return columns;
	}

	public int maxFunctions() {
		return maxFunctions;
	}

	public void setMaxFunctions(int mf) {
		if (matrixRows != null && mf <= matrixRows.size())
			maxFunctions = mf;
		else
			if (matrixRows != null)
				maxFunctions = matrixRows.size();
	}

	// Delta: p(x = 1|pi(x) = 1);
	public double getDelta(int r, int c) {
		if (r < 0 || c < 0 || r >= rows || c >= columns) {
			System.out.println("Indices exceed matrix dimensions: " + r + ", " + c);
			System.exit(1);
			//return -1;
		}
		
		return ((Parameter)matrixRows.elementAt(r).elementAt(c)).getDelta();
	}

	// Delta: p(x = 0|pi(x) = 1) = 1/length
	public double getLength(int r, int c) {
		if (r < 0 || c < 0 || r >= rows || c >= columns) {
			System.out.println("Indices exceed matrix dimensions: " + r + ", " + c);
			System.exit(1);
			//return -1;
		}

		return(((Parameter)
		        matrixRows.elementAt(r).elementAt(c)).getLength());
	}

	// Gets rows to make copying easier, since clone is being weird.
	@SuppressWarnings("unchecked")
	public Vector<Parameter> getRow(int r) {
		Vector<Parameter> retval = null;

		Object vectorOfParameterObject = matrixRows.elementAt(r).clone();

		if (vectorOfParameterObject instanceof Vector<?>) {

			retval = (Vector<Parameter>) vectorOfParameterObject;
		}

		return retval;
	}

	public void setRow(Vector<Parameter> v, int r) {
		matrixRows.remove(r);
		matrixRows.add(r, v);
	}

	public void setRowConsistent(Vector<Parameter> v, int r) {
		Vector<Parameter> ro = matrixRows.elementAt(r);
		double tot = 0;

		// Find ratio of old row and new row
		for (int j = 0; j < ro.size() - 1; j++) {
			tot += 1 - ro.elementAt(j).getDelta();
		}

		double tot2 = 0;

		for (int j = 0; j < v.size() - 1; j++) {
			tot2 += 1 - v.elementAt(j).getDelta();
		}

		// Re"normalize" new row to sum to the value of old row.
		for (int j = 0; j < v.size() - 1; j++) {
			Parameter p = v.elementAt(j);
			p.setDelta(1 - ((1 - p.getDelta())*tot / tot2));
		}

		matrixRows.remove(r);
		matrixRows.add(r, v);
	}

	/* sets the parameter at row r column c to v */
	public void setDelta(int r, int c, double v) {
		if (r >= 0 && r < rows && c >= 0 && c < columns) {
			((Parameter)matrixRows.elementAt(r)
			 .elementAt(c)).setDelta(v);
		}
		else
			System.out.println("Error: Can't set transition matrix at row "
			                   + r + ", column " + c);
	}

	public void addDelta(double[][] v) {
		if (v.length == rows && v[0].length == columns) {
			for (int r = 0; r < rows; r++) {
				for (int c = 0; c < columns; c++) {
					((Parameter)matrixRows.elementAt(r)
					 .elementAt(c)).addDelta(v[r][c]);
				}
			}
		}
		else
			System.out.println("Error: Number of rows, columns do not match");
	}

	public void addDeltaPositive(double[][] v) {
		if (v.length == rows && v[0].length == columns) {
			for (int r = 0; r < rows; r++) {
				for (int c = 0; c < columns; c++) {
					((Parameter)matrixRows.elementAt(r)
					 .elementAt(c)).addDeltaPositive(v[r][c]);
				}
			}
		}
		else
			System.out.println("Error: Number of rows, columns do not match");
	}

	public void addLengthPositive(double[][] v) {
		if (v.length == rows && v[0].length == columns) {
			for (int r = 0; r < rows; r++) {
				for (int c = 0; c < columns; c++) {
					((Parameter)matrixRows.elementAt(r)
					 .elementAt(c)).addLengthPositive(v[r][c]);
				}
			}
		}
		else
			System.out.println("Error: Number of rows, columns do not match");
	}

	/* converts the transition matrix in 1-q form to theta form */
	public void convertToTheta() {
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				((Parameter)matrixRows.elementAt(r)
				 .elementAt(c)).convertToTheta();
			}
		}
	}

	/* converts the transition matrix from theta form to 1-q */
	public void convertToOneMinusQ() {
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				((Parameter)matrixRows.elementAt(r)
				 .elementAt(c)).convertToOneMinusQ();
			}
		}
	}

	public void addDelta(int r, int c, double v) {
		if (r >= 0 && r < rows && c >= 0 && c < columns) {
			((Parameter)matrixRows.elementAt(r)
			 .elementAt(c)).addDelta(v);
		}
		else
			System.out.println("Error: Can't set numerator at row "
			                   + r + ", column " + c);
	}

	public boolean addDeltaPositive(int r, int c, double v) {
		if (r >= 0 && r < rows && c >= 0 && c < columns) {
			return ((Parameter)matrixRows.elementAt(r)
			        .elementAt(c)).addDeltaPositive(v);
		}
		else
			System.out.println("Error: Can't set numerator at row "
			                   + r + ", column " + c);

		return false;
	}

	public void addAlphaPositive(int index, double v) {
		if (index >= 0 && index < alphaParams.size()) {
			double newAlpha =
			  alphaParams.elementAt(index).doubleValue() + v;

			if (newAlpha > 0)
				alphaParams.set(index, new Double(newAlpha));
		}
		else
			System.out.println("Error: Can't set alphaParam at index "
			                   + index);
	}

	public void addLengthPositive(int r, int c, double v) {
		if (r >= 0 && r < rows && c >= 0 && c < columns) {
			((Parameter)matrixRows.elementAt(r)
			 .elementAt(c)).addLengthPositive(v);
		}
		else
			System.out.println("Error: Can't add length positive at "
			                   + r + ", column " + c);
	}

	public void initializeToUniform(double epsilon) {
		double selfV = 1.0 - epsilon * (columns - 1);

		if (selfV < 0.0) {
			epsilon = epsilon / (1.0 - selfV);
		}

		for (int i = 0; i < rows; i++) {
			Vector<Parameter> v = matrixRows.elementAt(i);

			for (int j = 0; j < columns; j++) {
				if (i != j)
					v.add(j, new Parameter(epsilon));
				else
					v.add(j, new Parameter(selfV));
			}
		}
	}


	public void initializeToOnes() {
		for (int i = 0; i < rows; i++) {
			Vector<Parameter> v = matrixRows.elementAt(i);

			for (int j = 0; j < columns; j++) {
				v.add(j, new Parameter(1));
			}
		}
	}

	public void initializeToZero() {
		for (int i = 0; i < rows; i++) {
			Vector<Parameter> v = matrixRows.elementAt(i);

			for (int j = 0; j < columns; j++) {
				v.add(j, new Parameter(0));
			}
		}
	}

	// This is the initialization that is used
	@SuppressWarnings("unused")
	public void initializeToDistanceSimple(int row,
	                                       int[] distances,
	                                       double scale) {
		double total = 0;
		s = scale;

		if (row < 0 || row > rows) {
			System.out.println("Indices exceed matrix dimensions: " + row);
			System.exit(1);
			//return;
		}

		Vector<Parameter> v = matrixRows.elementAt(row);

		for (int i = 0; i < columns; i++) {
			if (i < distances.length) {
				int dist = distances[i];

				if (dist < 1)
					dist = 1;

				if (dist > 20)
					dist = 20;

				if (i != row) {
					double par = scale + ((double)dist / (double)scale);
					v.add(i, new Parameter(par));
					total += par;
				}
				else {
					double par = scale;
					v.add(i, new Parameter(par));
					total += par;
				}
			}
		}
	}

	// This is the initialization that is used for rate matrix
	@SuppressWarnings("unused")
	public void initializeToDistanceSimpleRate(int row,
	    int[] distances,
	    double scale) {
		double total = 0;
		s = scale;

		if (row < 0 || row > rows) {
			System.out.println("Indices exceed matrix dimensions: " + row);
			System.exit(1);
			//return;
		}

		Vector<Parameter> v = matrixRows.elementAt(row);

		for (int i = 0; i < columns; i++) {
			if (i < distances.length) {
				int dist = distances[i];

				if (dist < 1)
					dist = 1;

				if (dist > 20)
					dist = 20;

				if (i != row) {
					double par = 1.0;
					v.add(i, new Parameter(par));
					total += par;
				}
				else {
					double par = 0.5;
					v.add(i, new Parameter(par));
					total += par;
				}
			}
		}
	}


	// This just sets the names associated with each column/row
	// mainly used for printing purposes (now).
	public void setGOTerms(Vector<Integer> names) {
		goTerms = names;
	}

	public Vector<Integer> getGOTerms() {
		return goTerms;
	}

	public void addDelta(int i, int j, int s) {
		((Parameter)matrixRows.elementAt(i).elementAt(j)).addDelta(s);
	}

	public Object clone() {
		return this.clone();
	}

	/* prints out the current matrix, given a filename */
	public void printOutGOTermConversionRateMatrix(String filename) {
		PrintStream fout = null;
		
		try {
			fout =  new PrintStream(new FileOutputStream(new File(filename)));
		} catch (Exception ioe) {
			System.err.println("Transition Matrix Output Error: " + filename + " " +ioe.toString());
			System.exit(1);
		}	
		
		fout.println("# scale parameter: " + s);
		
		if (goTerms != null) {
			for (int i = 0; i < goTerms.size(); i++) {
				fout.print("\t" + goTerms.elementAt(i));
			}

			fout.println();
		}

		for (int i = 0; i < rows && i < goTerms.size(); i++) {
			// need to print out name of position!
			if (goTerms != null)
				fout.print(goTerms.elementAt(i));

			for (int j = 0; j < columns; j++) {
				fout.print("\t" + getDelta(i, j));
			}

			fout.println();
		}

		fout.close();
	}

	/* prints out the current matrix, given a filename */
	public void printOutGOTermConversionRateMatrix() {
		System.out.println("# scale parameter: " + s);

		if (goTerms != null) {
			for (int i = 0; i < goTerms.size(); i++) {
				System.out.print("\t" + goTerms.elementAt(i));
			}

			System.out.println();
		}

		for (int i = 0; i < rows && i < goTerms.size(); i++) {
			// need to print out name of position!
			if (goTerms != null)
				System.out.print(goTerms.elementAt(i));

			for (int j = 0; j < columns; j++) {
				System.out.print("\t" + getDelta(i, j));
			}

			System.out.println();
		}
	}

	/* reads in scale parameters from file */
	public boolean readInScale(String filename) {
		scaleParams = GenericMathFunctions.readInHashtable(filename);
		return(scaleParams != null);
	}

	/* prints out the scale parameters, given a filename */
	public void printOutScale(String filename) {
		if (scaleParams == null)
			return;

		PrintStream fout;

		try {
			fout =  new PrintStream(new FileOutputStream(new File(filename)));
			Enumeration<String> scaleNames = scaleParams.keys();

			while (scaleNames.hasMoreElements()) {
				String scaleName = scaleNames.nextElement();
				fout.println(scaleName + "\t" + scaleParams.get(scaleName));
			}

			fout.close();
		}
		catch (Exception ioe) {
			System.err.println("PrintOutScaleParams: " + filename + " " +
			                   ioe.getMessage());
			System.exit(1);
		}

	}

	/* gets hashtable with scale parameters */
	public Hashtable<String, Double> getScale() {
		return scaleParams;
	}

	public double getScale(String name) {
		if (scaleParams.containsKey(name)) {
			return scaleParams.get(name).doubleValue();
		}

		return 0.0;
	}

	/* sets values of the scale parameter */
	public void setScale(String name, double value) {
		if (scaleParams == null)
			scaleParams = new Hashtable<String, Double>();

		scaleParams.put(name, new Double(value));
	}

	/* reads in alpha parameters from file */
	public boolean readInAlpha(String filename) {
		alphaParams = GenericMathFunctions.readInVector(filename);
		return(alphaParams != null);
	}

	/* prints out the alpha parameters, given a filename */
	public void printOutAlpha(String filename) {
		if (alphaParams == null)
			return;

		PrintStream fout;

		try {
			fout =  new PrintStream(new FileOutputStream(new File(filename)));

			for (int i = 0; i < alphaParams.size(); i++) {
				fout.println(alphaParams.elementAt(i) + "\t");
			}

			fout.close();
		}
		catch (Exception ioe) {
			System.err.println("PrintOutAlphaParams: " + filename + " " +
			                   ioe.getMessage());
			System.exit(1);
		}

	}

	/* prints out the alpha parameters */
	public void printOutAlpha() {
		if (alphaParams == null)
			return;

		for (int i = 0; i < alphaParams.size(); i++) {
			System.out.println(alphaParams.elementAt(i) + "\t");
		}
	}

	/* gets vector with alpha parameters */
	public Vector<Double> getAlpha() {
		return alphaParams;
	}

	/* gets double from alpha parameters */
	public double getAlpha(int i) {
		return ((alphaParams.elementAt(i)).doubleValue());
	}

	/* sets values of the alpha parameter */
	public void setAlpha(int index, double value) {
		if (alphaParams == null)
			alphaParams = new Vector<Double>();

		if (alphaParams.size() <= index) {
			for (int i = alphaParams.size(); i < index; i++)
				alphaParams.add(i, new Double(0.0));

			alphaParams.add(index, new Double(value));
		}
		else
			alphaParams.set(index, new Double(value));
	}

	private class Parameter extends Object {
		public double delta; // P(x=1 | pi(x) = 1);

		public Parameter(double d) {
			delta = d;
		}

		public double getDelta() {
			return delta;
		}

		public double getLength() {
			return ((1.0) / (double)(1 - delta));
		}

		public void setDelta(double d) {
			delta = d;
		}

		public void setLength(double l) {
			delta = (1.0 - ((double)1.0 / (double)l));
		}

		public void addDelta(double x) {
			delta += x;
		}

		public boolean addDeltaPositive(double x) {
			delta += x;

			if (delta < 0.0) {
				delta -= x;
				return false;
			}

			return true;
		}

		public void addLengthPositive(double l) {
			double len = getLength();
			len += l;
			System.out.println("length: " + len + ", change: " + l);

			if (len < 1.0) {
				len -= l;
			}
			else {
				setLength(len);
			}
		}

		public void convertToTheta() {
			delta = -GenericMathFunctions.logSafe(delta);
		}

		public void convertToOneMinusQ() {
			delta = Math.exp(-delta);
		}
	}

	/////////////////////////////////////////
	// Functions for interfacing with LAPACKMatrixOperationsWrapper
	/////////////////////////////////////////

	// Full version of power set
	// This function builds the rate transition matrix for
	// the Markov chain.
	// If this is changed, also change precomputePhiGradient in ProbabilisticReconciledPhylogenyObject
	public void buildMarkovTransitionRateMatrix() {
		ProteinFunctionMarkovState parent = new ProteinFunctionMarkovState(rows, maxFunctions);
		ProteinFunctionMarkovState child  = new ProteinFunctionMarkovState(rows, maxFunctions);
		System.out.println("Parent power set size: " + parent.powerSetSize());
		System.out.println("Child power set size: " + child.powerSetSize());

		if (transRateMatrix == null)
			transRateMatrix = new LAPACKMatrixOperationsWrapper(parent.powerSetSize(), child.powerSetSize());

		double sum = 0;
		double prod = 1;
		int count = 1;

		while (child.hasNext()) {
			child.getNextNonZeros();
			prod = 1;

			if (child.getCurrentSum() == 1) {
				for (int i = 0; i < child.length(); i++) {
					if (child.elementAt(i) == 1)
						prod *= getAlpha(i);
				}
			}
			else {
				prod = 0;
			}

			transRateMatrix.set(0, count++, prod);
			sum += prod;
		}

		transRateMatrix.set(0, 0, -sum); // so that row sums to 0

		int childCount = 0;
		int parentCount = 1;
		child.Reset();
		double oneSum = 0.0;

		while (parent.hasNext()) {
			parent.getNextNonZeros();
			sum = 0;
			childCount = 0;
			child.Reset();

			while (child.hasNext()) {
				child.getNext();
				//child.printPowerSet();
				oneSum = 0;
				int diffIndex = parent.offByOne(child);

				if (diffIndex != -1 && diffIndex != parent.length()) {
					for (int i = 0; i < parent.length(); i++) {
						if (parent.elementAt(i) == 1) {
							for (int j = 0; j < child.length(); j++) {
								if (i == j && child.elementAt(j) == 0)
									oneSum += getDelta(i, i);
								else
									if (i != j && child.elementAt(j) == 1
									    && parent.elementAt(j) == 0) {
										oneSum += getDelta(i, j);
										oneSum += getAlpha(j);
									}
							}
						}
					}
				}

				if (parentCount != childCount)
					sum += oneSum;

				transRateMatrix.set(parentCount, childCount++, oneSum);
			}

			transRateMatrix.set(parentCount, parentCount, -sum); //ensure rows sum to 0
			//System.out.println("Setting: "+parentCount+","+parentCount+" to "+(-sum));
			parentCount++;
		}

		//printOutMatrix();
		//printOutAlpha();
		transRateMatrix.print();
		//System.out.println("Probability Matrix: (exp to 1)");
		//transRateMatrix.rowsSumToZero();
		//transRateMatrix.print();
		//transRateMatrix.sumRows();
		//System.out.println("Printing matrix exp (1.0):");
		computeMatrixExp(1.0, false);
		//exlAPACKMatrixOperationsWrapperTemp.print();
		//computeMatrixExp(0.05, false);
		//System.out.println("Printing matrix exp (0.05):");
		//exlAPACKMatrixOperationsWrapperTemp.print();
		//exlAPACKMatrixOperationsWrapperSummary.print();
		exlAPACKMatrixOperationsWrapperTemp = null;
		currentScaling = 0.0;
		summaryNotCurrent = true;
	}

	// Takes power set indices i and j
	// and returns element of Q
	public double getQElement(int i, int j) {
		if (transRateMatrix == null)
			buildMarkovTransitionRateMatrix();

		//System.out.println("In getQElement");
		return transRateMatrix.get(i, j);
	}

	// put into matrix that holds parent ps index, child index
	// and normalize (for getting rid of 00 state)
	public void computeMatrixExp(double scale, boolean summary) {
		if (currentScaling != scale) {
			//Date d = new Date();
			//long time1 = d.getTime();
			exlAPACKMatrixOperationsWrapperTemp = transRateMatrix.matrixExponential(scale);
			//d = new Date();
			//long time2 = d.getTime();
			//System.out.println("Time for matrix exponential: "+(time2-time1));
			//System.out.println("Scaling: "+scale+", "+currentScaling);
			summaryNotCurrent = true;
			currentScaling = scale;
		}

		if (!summary ||
		    (currentScaling == scale && !summaryNotCurrent))
			return;

		//System.out.println("recomputing summary");
		int functions = matrixRows.size();

		if (childPS == null)
			childPS = new ProteinFunctionMarkovState(functions, maxFunctions);
		else
			childPS.Reset();

		if (exlAPACKMatrixOperationsWrapperSummary == null)
			exlAPACKMatrixOperationsWrapperSummary = new LAPACKMatrixOperationsWrapper(transRateMatrix.rows, functions);
		else
			exlAPACKMatrixOperationsWrapperSummary.zero();

		//exlAPACKMatrixOperationsWrapperTemp.print();
		int setIndex = 0;

		while (childPS.hasNext()) {
			childPS.getNextNonZeros();
			setIndex++;

			//childPS.getNext();
			for (int j = 0; j < childPS.length(); j++) {
				if (childPS.elementAt(j) == 1) {
					for (int i = 0; i < exlAPACKMatrixOperationsWrapperTemp.rows; i++) {
						//System.out.println("in computeMatrixExp");
						exlAPACKMatrixOperationsWrapperSummary.add(i, j,
						                     exlAPACKMatrixOperationsWrapperTemp.get(i, setIndex));
						//System.out.println("Summarizing "+i+", "+setIndex
						//		   +" to "+i+","+j
						//		   +" and prob "
						//		   +exlAPACKMatrixOperationsWrapperTemp.get(i,setIndex));
					}

					//exlAPACKMatrixOperationsWrapperSummary.print();
				}
			}
		}

		//exlAPACKMatrixOperationsWrapperSummary.print();
		summaryNotCurrent = false;
	}

	public LAPACKMatrixOperationsWrapper getMatrixExp(double rate, double distance) {
		computeMatrixExp(rate*distance, false);
		return exlAPACKMatrixOperationsWrapperTemp;
	}

	public LAPACKMatrixOperationsWrapper getQMatrix() {
		if (transRateMatrix == null)
			buildMarkovTransitionRateMatrix();

		return transRateMatrix;
	}

	// grab appropriate
	public double getExpProb(int parentPSIndex, int childIndex,
	                         double rate, double distance) {
		if (exlAPACKMatrixOperationsWrapperSummary == null || rate*distance != currentScaling)
			computeMatrixExp(rate*distance, true);

		//System.out.println("in getExpProb");
		double p = exlAPACKMatrixOperationsWrapperSummary.get(parentPSIndex, childIndex);

		if (p > 0)
			return p;
		else
			return probEpsilon;
	}

	// grab appropriate
	public double getExpProb(int parentPSIndex, int childPSIndex, double t) {
		if (transRateMatrix == null)
			buildMarkovTransitionRateMatrix();

		//Date d = new Date();
		//long time = d.getTime();
		if (exlAPACKMatrixOperationsWrapperTemp == null || t != currentScaling)
			computeMatrixExp(t, false);

		//d = new Date();
		//long time2 = d.getTime();
		//System.out.println("in getExpProb3: "+(time2-time)+", t = "+t+", currentscaling: "+currentScaling);
		double p = exlAPACKMatrixOperationsWrapperTemp.get(parentPSIndex, childPSIndex);

		if (p > 0)
			return p;
		else
			return probEpsilon;
	}

	public int getSizeExpMatrix() {
		if (transRateMatrix != null)
			return transRateMatrix.rows;
		else
			return 0;
	}

}

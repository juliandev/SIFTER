/**
 * SIFTER project. This implements the high level core
 * functionality for the inference and learning methods.
 * $Id: Sifter.java,v 1.2 2005/05/27 06:56:01 sprite Exp $
 *
 * Copyright (C) 2010, Barbara Engelhardt (bee@compbio.berkeley.edu)
 *
 * @author Barbara Engelhardt (primary investigator and author)
 * @author Steven R. Chan (later code hacking, documentation, GUI)
 * @version 1.0
 */

package sifter.main;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import gape.genetic_algorithm.Population;
import gape.output.PrintFiles;
import sifter.components.ExpectationMaximizationObject;
import sifter.components.GOOntologyWrapper;
import sifter.components.GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject;
import sifter.components.GenericMathFunctions;
import sifter.components.PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser;
import sifter.components.ProbabilisticReconciledPhylogenyObject;
import sifter.components.ProteinAnnotationObjectWithProbabilityMethods;
import sifter.components.SifterExecutionSettingsContainer;
import sifter.components.ProbabilisticReconciledPhylogenyObject.Node;

import java.io.File;

public class SifterPipelineObject {
	/**
	 * Holds user's preferences.
	 */
	private SifterExecutionSettingsContainer settings;
	
	/**
	 * Arbitrary constants not provided by user
	 */
	private static final class CONSTANTS {
		public static final int ITERATIONS_DEFAULT = 4000;
		public static final double STEPSIZE_DEFAULT = 0.01;
		public static final double CUTOFF_DEFAULT = 0.000115;
		
		public static final int FOLDS_DEFAULT = 0;
		public static final int TRUNCATION_DEFAULT = 10;

		public static final double SPECIATION_RATE = 0.03;
		public static final double DUPLICATION_RATE = 0.05;
	}
	
	// Global to all modes:
	private boolean em;
	private String familyFilename;
	private String scaleParamsFilename;
	private String alphaParamsFilename;
	private int setfolds;
	private int truncation;
	private String outputDirectory;

	private Hashtable<String, ProteinAnnotationObjectWithProbabilityMethods> proteinList; // Indexed by swissprot numbers. KEY: swissprotID (not name); VALUE: PFunprotein
	private PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser family; // Indexed by pfam ID's

	/** Once we get a new settings object, initialize any variables
	 * that we care about.
	 * For example, I know SifterPipelineObject cares about verbose, family name, etc.
	 */
	private void initDefaultSettings() {
		if (settings.getSetting("truncation") == null)
			settings.setSetting("truncation", CONSTANTS.TRUNCATION_DEFAULT);
		if (settings.getSetting("folds") == null)
			settings.setSetting("folds", new Integer(CONSTANTS.FOLDS_DEFAULT));
		if (settings.getSetting("stepsize") == null)
			settings.setSetting("stepsize", new Double(CONSTANTS.STEPSIZE_DEFAULT));			
		if (settings.getSetting("cutoff") == null)
			settings.setSetting("cutoff", new Double(CONSTANTS.CUTOFF_DEFAULT));
		if (settings.getSetting("iterations") == null)
			settings.setSetting("iterations", new Integer(CONSTANTS.ITERATIONS_DEFAULT));
		
		this.em = ((Boolean)this.settings.getSetting("em")).booleanValue();
		this.familyFilename = (String) this.settings.getSetting("familyFilename");
		this.scaleParamsFilename = (String) this.settings.getSetting("scaleParamsFilename");
		this.alphaParamsFilename = (String) this.settings.getSetting("alphaParamsFilename");
		this.setfolds = ((Integer)this.settings.getSetting("folds")).intValue();
		this.truncation = ((Integer)this.settings.getSetting("truncation")).intValue();
		
		String outputFile = (String)settings.getSetting("output");
		File f = new File(outputFile);
		f = new File(f.getParent());
		
		this.outputDirectory = f.getAbsolutePath();
	}
	
	public SifterPipelineObject(SifterExecutionSettingsContainer settings) throws ClassNotFoundException {
		this.settings = settings;
		initDefaultSettings();
		
		proteinList = null;
		family = null;
		
		String proteinFilename	 = (String)settings.getSetting("pli");
		String phylogenyFilename = (String)settings.getSetting("phylo");
		System.out.println(phylogenyFilename);
		String fxFileName		 = (String)settings.getSetting("fx");
		String scaleFilename	 = (String)settings.getSetting("scale");
		String alphaFilename	 = (String)settings.getSetting("alpha");
		
		String ontology = ((String) this.settings.getSetting("ontology"));
		boolean noIEA = !((Boolean)this.settings.getSetting("iea")).booleanValue();
		
		String runmode = (String)settings.getSetting("runmode");
		
		GOOntologyWrapper pfgodag = buildDatasetGODAG(proteinFilename, ontology, noIEA);
		
		// The "Generate" step creates parameter files that can be loaded for subsequent
		// inference calls. 
		if (runmode.equals("generate")) {
			/**
			 * Generates a set of input parameters for the inference problem. All
			 * parameters must be >= 0. The lower the theta_m,n parameters, the
			 * more likely function m will mutate to function n. The s parameter
			 * is a rate of mutation.
			 */
			if ((Boolean)settings.getSetting("verbose")) {
				System.out.println("Mode: Generating Network Parameters.");
				System.out.println("**********************************************");
			}
			
			double scale = 20.0;
			int numLeaves = pfgodag.getNumLeaves();
			System.out.println("There are " + numLeaves + " candidate functions.");
			
			// Build the transition matrix
			GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx = new GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject(numLeaves, numLeaves);
			
			pfx.setScale("species", CONSTANTS.SPECIATION_RATE);
			pfx.setScale("duplication", CONSTANTS.DUPLICATION_RATE);
			for (int i = 0; i < numLeaves; i++)
				pfx.setAlpha(i, 1.0);
			
			// Uses the GO DAG to name
			Vector<Integer> leaves = pfgodag.probabilisticDAG.getAllLeaves();
			int[] dist = new int[leaves.size() + 1];
			
			// Note that unless EM is done, there's no reason to have
			// the names stored in the transition matrix file. But
			// we leave the option of EM open, so we provide the GO terms of 
			Vector<Integer> rowNames = new Vector<Integer>();
			for (int i = 0; i < leaves.size(); i++) {
				Integer nodei = leaves.elementAt(i);
	
				for (int j = 0; j < leaves.size(); j++) {
					if (i == j)
						dist[i] = 0;
					else
						dist[j] = pfgodag.probabilisticDAG.getTreeDistance(nodei.intValue(), leaves.elementAt(j).intValue());
				}
	
				pfx.initializeToDistanceSimpleRate(i, dist, scale);
				rowNames.add(leaves.elementAt(i));
			}
			pfx.setGOTerms(rowNames);
			
			String ConvMatrixOutFile = this.outputDirectory + "/" + fxFileName;
			pfx.printOutGOTermConversionRateMatrix(ConvMatrixOutFile);
			
			String ScaleOutFile = this.outputDirectory + "/" + scaleFilename;
			pfx.printOutScale(ScaleOutFile);
			
			String AlphaOutFile = this.outputDirectory + "/" + alphaFilename;
			pfx.printOutAlpha(AlphaOutFile);
			
			//System.out.println("Output written to " + ConvMatrixOutFile + ", " + ScaleOutFile + ", " + AlphaOutFile);
		} else {
			//System.out.println(phylogenyFilename);
			// Build the graphical model
			buildGraphicalModel(pfgodag, phylogenyFilename, noIEA);
			
			if (runmode.equals("default")) {
				if ((Boolean)this.settings.getSetting("verbose")) {
					System.out.println("Mode: Use Exact Network Parameters.");
					System.out.println("**************************************************");
				}
				
				inferWithTransitionMatrixExact(pfgodag, getFamily());
			} else if (runmode.equals("xvalidate")) {
				if ((Boolean)this.settings.getSetting("verbose")) {
					System.out.println("Mode: Cross Validation, Estimating Network Parameters.");
					System.out.println("**********************************************");
				}
				
				crossValidation(pfgodag, getFamily());
			} else if (runmode.equals("expectationmaximization")) {
				if ((Boolean)settings.getSetting("verbose")) {
					System.out.println("Mode: Estimating Network Parameters.");
					System.out.println("**********************************************");
				}
				long startTime = System.currentTimeMillis();
				
				estimateParameters(pfgodag, getFamily());
				
				long endTime = System.currentTimeMillis() - startTime;
				
				System.out.println(endTime);
				
			} else if (runmode.equals("gaparameterestimation")) {
				System.out.println("GAPE");
				System.out.println(settings.getSetting("population"));
				System.out.println(settings.getSetting("generations"));
				System.out.println(settings.getSetting("geneticoperators"));
				System.out.println(this.getFamily().getID());
				
				int populationSize = (int)settings.getSetting("population");
				int iterations = (int)settings.getSetting("generations");
				int[] rowNames = pfgodag.probabilisticDAG.getAllLeaves().stream().mapToInt((Integer i) -> i.intValue()).toArray(); 
				double geneticoperators = (double)settings.getSetting("geneticoperators");
				String idFamily = (String) this.settings.getSetting("family");
				
				long startTime = System.nanoTime(); //Start time counter
				
				Population population = new Population(populationSize, iterations, geneticoperators, rowNames, idFamily);
				population.getBest();
				
				double elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
				
				PrintFiles files = new PrintFiles();
				
				files.printTimeExecution("output/BestIndividual.txt", String.valueOf(elapsedTimeInSec));
				
				System.out.println(elapsedTimeInSec);
			}
		}
	}

	
	/*****************************************************************
	      ONLY USED BY -lib
	 ******************************************************************/

	// Only called from testLibrary() and testLibrarySmall().
	// Otherwise also seems useless.
//	public void printOutProteinList(String outFileName, Hashtable pList)
//	{
//		PrintStream fout;
//		try {
//			fout =  new PrintStream(new FileOutputStream(new File(outFileName)));
//
//			Enumeration pnames = pList.keys();
//			while(pnames.hasMoreElements()) {
//				ProteinAnnotationObjectWithProbabilityMethods pfp = (ProteinAnnotationObjectWithProbabilityMethods)pList.get(pnames.nextElement());
//				pfp.printProteinToFile(fout);
//			}
//			fout.close();
//		}
//		catch (Exception ioe) {
//			System.err.println("PrintOutProteinList: " + outFileName + " " +
//					ioe.getMessage());
//			System.exit(1);
//		}
//	}

	//////////////////////////////////////////////////
	// Ontology building section                    //
	//////////////////////////////////////////////////

	/**
	 * Builds DAG: Most of the time she uses this (string-calling).
	 * All this does is take an input file "inFiles", stick it in a
	 * vector, and pass it on like a hot potato to another overloaded
	 * buildDatasetGODAG function (in this case, completely
	 * passing all of its own parameters to the overloaded f(x)
	 * except ignoring the "outfileName" argument.)
	 *
	 * STEP 1.
	 * TODO: Check to see if this calls other overloaded
	 * buildDatasetGODAG() functions.
	 *
	 * @param inFamilyFile
	 * @param goFile
	 * @param noIEA
	 * @return
	 */

	/**
	 * STEP 2 in the buildDatasetGODAG potato-passing process.
	 * does the actual work of building a model from the files.
	 * In general, she uses the single file ... calls ImportGODataFromFile
	 *
	 * @param files
	 * @param goFile
	 * @param noIEA
	 * @return
	 * @throws ClassNotFoundException 
	 */
	public GOOntologyWrapper buildDatasetGODAG(String file, String goFile, boolean noIEA) throws ClassNotFoundException {
		GOOntologyWrapper pfgodag = new GOOntologyWrapper(goFile, true, true, true);
		pfgodag.setSettingsObject(this.settings);
		
		pfgodag.readGOFile();
		
		Vector<Integer> fns = null; // evidence functions;
		Vector<String> methods = null; // methods list
		Vector<Double> confidences = null;
		proteinList = new Hashtable<String, ProteinAnnotationObjectWithProbabilityMethods>();
		
		// counting the xml file with GO dag stats
		// print out at the end.
		PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser proteinFam = new PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser(file);
		Vector<ProteinAnnotationObjectWithProbabilityMethods> p = proteinFam.readInFromXMLFile(file);
		
		if (p != null) {
			//if ((Boolean)this.settings.getSetting("verbose"))
			//	System.out.println(file);

			fns = new Vector<Integer>();
			methods = new Vector<String>();
			confidences = new Vector<Double>();

			for (int i = 0; i < p.size(); i++) {
				ProteinAnnotationObjectWithProbabilityMethods pro = p.elementAt(i);

				//if(pro.getGONumber() != null) { // bee added for BMC experiments
				if (pro.getGONumbers() != null) {
					fns.addAll(pro.getGONumbers());
					methods.addAll(pro.getMocs());

					if (pro.getProbabilities() == null) {
						Vector<Double> ones = new Vector<Double>();

						for (int one_index = 0; one_index < methods.size(); ++one_index) {
							ones.add(Double.valueOf(1));
						}

						pro.setProbabilities(ones);
					}

					confidences.addAll(pro.getProbabilities());

				}

				if (!proteinList.contains(pro.getID())) {
					proteinList.put(pro.getID(), pro);
				}
			}

			// should be noIEA, but tally is only IDA
			// don't need confidences here?
			pfgodag.tallyFunctions(fns, methods);
		} else {
			System.out.println("No data in pli?");
		}

		family = proteinFam;
		pfgodag.pruneZeroHitsAndLeaves();
		// functions is now a GO DAG with original leaves removed if they did not appear directly in PLI
		pfgodag.padSingletonLeaves();

		// REFACTORME: gets R values based on leaves, but not all leaves have evidence, and not all evidence is in leaves...
		pfgodag.findRValue();

		Enumeration<String> proteinIDs = proteinList.keys();

		while (proteinIDs.hasMoreElements()) {
			// iterate through protein list,
			// incorporating evidence
			String pID = proteinIDs.nextElement();
			ProteinAnnotationObjectWithProbabilityMethods pro = proteinList.get(pID);

			if (pro.getFunctionRatio() != null) {
				pro.addRatioLikelihoods();
			}
			else
				if (pro.getGONumbers() != null) {
					methods = pro.getMocs();
					Vector<Integer> functions = pro.getGONumbers();
					pfgodag.incorporateEvidenceShort(functions,
					                                 methods, confidences, noIEA);

					if ((Boolean)this.settings.getSetting("verbose")) {
						for (int i = 0; i < methods.size(); i++) {
							System.out.println("Evidence "
							                   + proteinFam.getProteinNameFromID(pID).toUpperCase() + " "
							                   + methods.elementAt(i) + " "
							                   + functions.elementAt(i));
						}
					}

					double[] leaflike = pfgodag.pullOutLeafLikelihoods();
					pfgodag.clearDAGProbabilities();
					pro.addLikelihoods(leaflike);
				}
		}

		if ((Boolean)this.settings.getSetting("verbose"))
			System.out.println("Finished building GO graphs");

		return pfgodag;
	}

	/**
	 * Called by parseTestNetworkInputExact() and a WHOLE bunch of
	 * other functions!
	 *
	 * Tries comparing the reconciled tree with the one in GODAG.
	 * They're never merged. With the set of likelihoods from the GODAG,
	 * associate each protein with the leaf in the tree.
	 *
	 * TODO: Could have error checking ...
	 **/
	public void buildGraphicalModel(GOOntologyWrapper pfgodag, String nexfile, boolean noIEA) {
		// Read in their reconciled tree (not the GODAG)
		//
		ProbabilisticReconciledPhylogenyObject t = new ProbabilisticReconciledPhylogenyObject(family.getMaxAlignment());
		String reconciledFilename = nexfile;
		t.createReconciled(reconciledFilename);

		if ((Boolean)this.settings.getSetting("verbose"))
			t.printTree();

		// Add evidence at leaves with evidence
		// PFam-type names (not swiss prot ids)
		Vector<String> proteins = family.getProteinNames();

		for (int i = 0; i < proteins.size(); i++) {
			// src: Does comparison. Does not do error-checking.
			//match protein with reconciled tree leaf
			String currentProtein = proteins.elementAt(i);
			ProteinAnnotationObjectWithProbabilityMethods pfp = proteinList.get(family.getProteinName(currentProtein));

			if (t.hasNode(currentProtein) && pfp.getLeafLikelihoods() != null && (pfp.hasProperEvidence(settings))) {
				t.setNodeEvidenceProbabilities(currentProtein, pfp.getLeafLikelihoods());
				
				/*System.out.println(currentProtein);
				for (int ii = 0; ii <  pfp.getLeafLikelihoods().length; ii++) {
					System.out.println(pfp.getLeafLikelihoods()[ii]);
				}*/
				
				
			}
		}

		// and attach them to the family.
		// src: makes sure that each protein in this reconciled tree
		// has the appropriate likelihoods from the GODAG.
		t.setSingleNodeSampleLikelihoods(pfgodag.getLeafSubsetPrior(), pfgodag.getSingleLeafPrior());

		// SETs the tree (
		family.addTree(t);
		
		// TODO
		if ((Boolean)this.settings.getSetting("verbose")) {
			System.out.println("Escribiendo árbol en archivo");
		}
		
		t.printTreeAsFile("phylogeneticTree.txt");
	}

	/** Accessor for familyList.
	 * @return Hashtable: a list of families
	 */
	public PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser getFamily() {
		return family;
	}

	/**
	 * This is the main function for inferences that she uses;
	 * not the other inferWithTransitionMatrix() functions.
	 * learnTransitionMatrix().
	 *
	 * Actually exponential in terms of the number of functions in the tree
	 * (even though it's linear in the number of proteins).
	 * TODO: To run it on larger datasets: figure out how to speed this up,
	 *
	 * @param pfgodag
	 * @param hashtable
	 */

	public Hashtable<Node, double[]> inferWithTransitionMatrixExact(GOOntologyWrapper pfgodag,
	    PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser fam) {
		// iterate through protein list,
		// incorporating evidence

		GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx = new GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject(familyFilename);
		pfx.readInAlpha(alphaParamsFilename);
		pfx.readInScale(scaleParamsFilename);
		pfx.setMaxFunctions(truncation); // truncation here
		System.out.println("Setting truncation level to " + truncation);
		ExpectationMaximizationObject pfl = new ExpectationMaximizationObject(pfx, proteinList, pfgodag, settings);
		pfl.setSettingsObject(settings);
		pfl.setFamily(fam);
		return(pfl.inferPosteriorsExact(fam));
		// Note: this prints to the file provided by "--output"; using "prettyPrintResults" function in ExpectationMaximizationObject
		//fam.getTree().getLogLikelihood(pfx, scaleParams);
	}

	/**
	 * @param pfgodag
	 * @param hashtable
	 */
	public double estimateParameters(GOOntologyWrapper pfgodag, PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser fam) {
		// iterate through protein list,
		// incorporating evidence
		GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx = new GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject(familyFilename);
		pfx.readInScale(scaleParamsFilename);
		pfx.readInAlpha(alphaParamsFilename);

		ExpectationMaximizationObject pfl = new ExpectationMaximizationObject(pfx, proteinList,
		    pfgodag, settings);
		pfl.setSettingsObject(settings);
		pfl.setFamily(fam);
		pfl.ExpectationMaximization(pfx, fam);
		return 0.0;
	}

	/**
	 * crossValidation
	 *
	 */

	public double crossValidation(GOOntologyWrapper pfgodag, PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser fam) {
		return (crossValidation(pfgodag, fam, setfolds));
	}

	public double crossValidation(GOOntologyWrapper pfgodag, PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser fam, int folds) {
		// iterate through protein list,
		// incorporating evidence
		int correctGO = 0;
		int totalGO = 0;
		Vector<Node> toRemove = new Vector<Node>();
		int lpsize = 0;

		Vector<String> proteinNames = fam.getProteinNames();

		if (folds > 0)
			GenericMathFunctions.randomizeVector(proteinNames);

		// Build the list of proteins to remove
		ProbabilisticReconciledPhylogenyObject tree = fam.getTree();

		for (int i = 0; i < proteinNames.size(); i++) {
			//match protein with reconciled tree leaf
			String currentProtein = proteinNames.elementAt(i);
			ProteinAnnotationObjectWithProbabilityMethods p = proteinList.get(fam.getProteinName(currentProtein));

			if (p != null && tree.hasNode(currentProtein)) {
				ProbabilisticReconciledPhylogenyObject.Node n = tree.getNode(currentProtein);

				if (n.hasLocalProbabilities() && n.isLeaf()) {
					toRemove.addElement(n);
					lpsize = n.getLocalProbabilities().length;
				}
			}
		}

		if (folds == 0)
			folds = toRemove.size();

		int foldSize = (int)Math.floor((double)toRemove.size() / (double)folds);
		int foldSizeRemainder = toRemove.size() % folds;
		System.out.println("Fold size: " + foldSize);
		System.out.println("Number of holdouts: " + toRemove.size());
		int startIndex = 0;
		int endIndex = 0;
		double[][] lps = new double[toRemove.size()][lpsize];

		for (int i = 0; i < folds; i++) {
			startIndex = endIndex;
			endIndex = (endIndex + foldSize > toRemove.size()) ?
			           (toRemove.size()) : (endIndex + foldSize);

			if (i < foldSizeRemainder)
				++endIndex;

			for (int j = startIndex; j < endIndex; j++) {
				ProbabilisticReconciledPhylogenyObject.Node n = toRemove.elementAt(j);
				lps[j] = n.getLocalProbabilities();
				n.removeLocalProbabilities();

				if ((Boolean)this.settings.getSetting("verbose"))
					System.out.println("X Val: removing evidence for "
					                   + n.getNodeID());
			}

			// Run EM like usual
			GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx = new GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject(familyFilename);
			pfx.setMaxFunctions(truncation); // truncation here
			System.out.println("Setting truncation level to " + truncation);

			pfx.readInScale(scaleParamsFilename);
			pfx.readInAlpha(alphaParamsFilename);

			ExpectationMaximizationObject pfl =
			  new ExpectationMaximizationObject(pfx, proteinList,
			                                    pfgodag, settings);
			pfl.setSettingsObject(settings);
			pfl.setFamily(fam);
			Hashtable<Node, double[]> posteriors;

			if (em) {
				posteriors = pfl.ExpectationMaximization(pfx, fam);
			}
			else {
				posteriors = inferWithTransitionMatrixExact(pfgodag, fam);
			}

			// Check to see if we got held-out ones correct
			for (int k = startIndex; k < endIndex; k++) {
				ProbabilisticReconciledPhylogenyObject.Node n = toRemove.elementAt(k);
				double[] lp = lps[k];

				if (lp == null)
					System.out.println("Error: don't have local probabilities for " + n.getNodeID());

				if (posteriors.containsKey(n)) {
					double[] d = (double[])posteriors.get(n);
					System.out.print("x-val (" + n.getNodeID() + ") ");

					for (int di = 0; di < d.length; di++) {
						System.out.print(d[di] + " ");
					}

					System.out.println();
					double maxPrior = 0.0;
					double maxPosterior = 0.0;
					int[] maxPriorIndex = new int[1];
					int maxPosteriorIndex = 0;

					for (int j = 0; j < d.length; j++) {
						//System.out.println("LP["+j+"] = "+lp[j]);
						double lpj = lp[j];

						try {
							if (d[j] > maxPosterior) {
								maxPosterior = d[j];
								maxPosteriorIndex = j;
							}
							else
								if (GenericMathFunctions.areEqual(d[j], maxPosterior)) {
									maxPosteriorIndex = -1;
								}

							if (lpj < 0)
								lpj = Math.exp(lpj);

							if (lpj > maxPrior) {
								maxPrior = lpj;
								maxPriorIndex[0] = j;
							}
							else
								if (GenericMathFunctions.areEqual(lpj, maxPrior)) {
									int[] maxPriorIndexT =
									  new int[maxPriorIndex.length+1];

									for (int pi = 0; pi < maxPriorIndex.length; pi++)
										maxPriorIndexT[pi] = maxPriorIndex[pi];

									maxPriorIndexT[maxPriorIndex.length] = j;
								}
						}
						catch (ArrayIndexOutOfBoundsException e) {
							continue;
						}
					}

					System.out.println("In x-validation: "
					                   + n.getNodeID() + " node name, "
					                   + d.length + " = d length, "
					                   + maxPosteriorIndex
					                   + " = max posterior index, "
					                   + maxPriorIndex.length + " = max prior index, "
					                   + maxPosterior + " = max posterior, "
					                   + maxPrior + " = max prior");

					boolean thisCorrect = false;

					for (int pi = 0; pi < maxPriorIndex.length; pi++) {
						if (maxPosteriorIndex == maxPriorIndex[pi] &&
						    (maxPosteriorIndex != -1)) {
							correctGO++;
							thisCorrect = true;
						}
					}

					if (!thisCorrect) {
						System.out.print("Missed " + (n.getNodeID())
						                 + ", predicted "
						                 + maxPosteriorIndex
						                 + ", real");

						for (int pi = 0; pi < maxPriorIndex.length; pi++) {
							System.out.print(" " + maxPriorIndex[pi]);
						}

						System.out.println();
					}

					totalGO++;
				}

				tree.setNodeEvidenceProbabilities((String)n.getNodeID(),
				                                  lp);
			}
		}

		System.out.println("Cross-validation results: "
		                   + ((double)correctGO / (double)totalGO)
		                   + " (" + correctGO + " out of " + totalGO + ")");
		
		System.out.println("correctGO: " + correctGO + "\ttotalGO: " + totalGO);
		return ((double)correctGO / (double)totalGO);
	}

	/** Should we ever want to use an already existing SifterExecutionSettingsContainer object,
	 * we can direct our pointer to that.
	 * @param settings SifterExecutionSettingsContainer object to set to.
	 * @see Sifter.main()
	 */
	public void setSettingsObject(SifterExecutionSettingsContainer settings) {
		this.settings = settings;
		initDefaultSettings();
	}
}

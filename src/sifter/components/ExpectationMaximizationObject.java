/**
 *  ExpectationMaximizationObject coordinates inference and EM.
 *
 * Copyright 2003-2005 Barbara Engelhardt (bee@cs.berkeley.edu)
 * @author Barbara Engelhardt
 */

package sifter.components;

import java.util.Enumeration;
import java.util.Hashtable;
//import java.util.Random;
import java.util.Vector;

import sifter.components.ProbabilisticReconciledPhylogenyObject;
import sifter.components.ProbabilisticReconciledPhylogenyObject.Node;

import java.io.*;
import java.util.Date;

public class ExpectationMaximizationObject {
	private GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx;
	//private GOOntologyWrapper pfgodag;
	private PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser family;
	private Hashtable<Node, double[]> posteriors;
	private Hashtable<String, ProteinAnnotationObjectWithProbabilityMethods> proteins;
	//private Random rand;
	private SifterExecutionSettingsContainer settings;
	private boolean verbose;
	private Date ts;

	public ExpectationMaximizationObject(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pf,
	                                     Hashtable<String, ProteinAnnotationObjectWithProbabilityMethods> prots, GOOntologyWrapper pfgo,
	                                     SifterExecutionSettingsContainer s) {
		pfx = pf;
		family = null;
		//pfgodag = pfgo;
		proteins = prots;
		posteriors = null;
		//rand = new Random();
		setSettingsObject(s);
		checkFunctionsMatch(pf, pfgo);
		//checkFunctionsMatch(pfD, pfgo);
		ts = new Date();
	}

	public ExpectationMaximizationObject(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pf,
	                                     GOOntologyWrapper pfgo, SifterExecutionSettingsContainer s) {
		pfx = pf;
		family = null;
		posteriors = null;
		//pfgodag = pfgo;
		proteins = null;
		//rand = new Random();
		setSettingsObject(s);
		checkFunctionsMatch(pf, pfgo);
		//checkFunctionsMatch(pfD, pfgo);
		ts = new Date();
	}

	public void checkFunctionsMatch(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx, GOOntologyWrapper pfgo) {
		Vector<Integer> namesPFX = pfx.getGOTerms();
		Vector<Integer> namesDAG = pfgo.getLeafNames();
		boolean error = false;

		if (namesPFX.size() != namesDAG.size())
			error = true;
		else {
			for (int i = 0; i < namesPFX.size(); i++) {
				if (!namesDAG.contains(namesPFX.elementAt(i)))
					error = true;
			}
		}

		if (error) {
			System.out.println("*********** Error **********");
			System.out.print("Transition matrix functions do not match");
			System.out.println(" functions in pruned GO DAG");
			System.out.print("Transition file: [");

			for (int i = 0; i < namesPFX.size(); i++) {
				if (i + 1 < namesPFX.size())
					System.out.print(namesPFX.elementAt(i) + ",");
				else
					System.out.print(namesPFX.elementAt(i));
			}

			System.out.println("]");
			System.out.print("Pruned GO DAG: [");

			for (int i = 0; i < namesDAG.size(); i++) {
				if (i + 1 < namesDAG.size())
					System.out.print(namesDAG.elementAt(i) + ",");
				else
					System.out.print(namesDAG.elementAt(i));
			}

			System.out.println("]");
			System.exit(1);
		}
	}

	public void setFamily(PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser fam) {
		family = fam;
	}

	// Performs EM, used to estimate parameters.
	// Output currently put into output/
	// and returns Hashtable with posteriors from latest iteration
	// of parameter estimation.
	public Hashtable<Node, double[]> ExpectationMaximization(GOTermConversionMatrixAndMarkovTransitionMatrixInputAndMathObject pfx0,
	    PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser pfFam) {
		//System.out.println("Starting to learn transition matrix");
		//rootSamples = null;
		int iterations = ((Integer)this.settings.getSetting("iterations")).intValue();
		//int iterations = Integer.parseInt((String)this.settings.getSetting("iterations"));
		double rho = ((Double)this.settings.getSetting("stepsize")).doubleValue();
		double cutoff = ((Double)this.settings.getSetting("cutoff")).doubleValue();
		printOutTransitionMatrices(0);
		posteriors = null;
		boolean done = false;
		//Hashtable removed;
		int i = 0;

		for (i = 0; i < iterations && !done; i++) {
			inferPosteriorsExact(pfFam);
			ts = new Date();
			long starttime = ts.getTime();
			done = pfFam.getTree().maximizationStep(pfx0, pfx0.getScale(),
			                                        posteriors,
			                                        rho, cutoff, i);
			ts = new Date();
			long maxtime = ts.getTime();
			System.out.println("Max time: " + (maxtime - starttime));
			printOutTransitionMatrices(i + 1);
		}

		inferPosteriorsExact(pfFam);
		System.out.println("Printing transition matrix/scale parameters to "
		                   + "output/pfxIteration" + i + ".fx, output/scale"
		                   + i + ".fx");
		return posteriors;
	}

	/**
	 * @param removed
	 */
	/*
	  private void returnOneTenth(Hashtable removed) {
	Enumeration nodes = removed.keys();
	while(nodes.hasMoreElements()) {
	    String node = (String)nodes.nextElement();
	    ProbabilisticReconciledPhylogenyObject tree = family.getTree();
	    tree.setNodeEvidenceProbabilities(node,
					      (double[])removed.get(node));
	}
	}*/

	/**
	 * @return
	 */
	/*
	  private Hashtable removeOneTenth() {
	Hashtable v = new Hashtable();
	Vector tree = family.getTree().getTree();
	for(int node = tree.size()-1; node >= 0; node--) {
	    // Get the name of the node
	    ProbabilisticReconciledPhylogenyObject.Node n = (ProbabilisticReconciledPhylogenyObject.Node)tree.elementAt(node);
	    //Case 1: Leaf with evidence
	    if(n.hasLocalProbabilities() && n.isLeaf()) {
		double tenth = rand.nextDouble();
		if(tenth < 0.1) {
		    v.put(n.getNodeID(), n.getLocalProbabilities());
		    n.removeLocalProbabilities();
		}
	    }
	}
	return v;
	}*/

	/**
	 * Method printOutTransitionMatrices.
	 */
	private void printOutTransitionMatrices(int iter) {
		pfx.printOutGOTermConversionRateMatrix("output/pfxIteration" + iter + ".fx");
		pfx.printOutScale("output/scale" + iter + ".fx");
		pfx.printOutAlpha("output/alpha" + iter + ".fx");
	}

	/**
	 * @param pfx
	 * @param fams
	 */
	@SuppressWarnings("unused")
	public Hashtable<Node, double[]> inferPosteriorsExact(PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser pfFam) {
		ts = new Date();
		long starttime = ts.getTime();
		posteriors = pfFam.getTree().propagateExactThroughoutTree(pfx, pfx.getScale());
		ts = new Date();
		long endtime = ts.getTime();
		//System.out.println("Inference time: " + (endtime - starttime));
		Hashtable<Node, double[]> gammas = pfFam.getTree().getGammas();
		Hashtable <Node, double[]> deltas = pfFam.getTree().getDeltas();
		Enumeration<ProbabilisticReconciledPhylogenyObject.Node> posts = posteriors.keys();

		while (posts.hasMoreElements()) {
			ProbabilisticReconciledPhylogenyObject.Node n = (posts.nextElement());

			if (posteriors.containsKey(n)) {
				double[] d = posteriors.get(n);
				double[] gamma = gammas.get(n);
				double[] delta = deltas.get(n);

				// Print out the gamma/deltas for each protein
				if (verbose) {
					System.out.println(n.getNodeID() + ": " + d);
					System.out.print("   gamma ");

					for (int i = 0; i < gamma.length; i++) {
						System.out.print(Math.exp(gamma[i]) + " ");
					}

					System.out.println();
					System.out.print("   delta ");

					for (int i = 0; i < delta.length; i++) {
						System.out.print(Math.exp(delta[i]) + " ");
					}

					System.out.println();
				}
			}
			else
				System.out.println("Doesn't contain key " + n.getNodeID());
		}

		prettyPrintResults(pfx.getGOTerms(), posteriors);

		if (verbose)
			prettyPrintProteinsValidation(pfx.getGOTerms(),
			                              posteriors);
		
		prettyPrintProteinsValidation(pfx.getGOTerms(), posteriors);

		//prettyPrintValidation(pfx.getNames(), posteriors);
		return(posteriors);
	}

	/**
	 *
	 */
	private void prettyPrintResults(Vector<Integer> names, Hashtable<Node, double[]> posteriors) {
		int[] order = new int[names.size()];
		String outFileName = (String)this.settings.getSetting("output");

		PrintStream fout;

		try {
			fout =  new PrintStream(new FileOutputStream(new File(outFileName)));
			fout.print("#Names");
			for (Integer goId : names) {
				fout.print("\t" + "GO:"+String.format("%07d", goId));
			}
			fout.print("\n");

			//System.out.println(names);
			Enumeration<Node> posts = posteriors.keys();

			while (posts.hasMoreElements()) {
				ProbabilisticReconciledPhylogenyObject.Node n = (posts.nextElement());

				if (posteriors.containsKey(n) && n.isLeaf()) { // lose n.isLeaf() to print internal nodes too.
					//System.out.print(n.getNodeID());
					fout.print(n.getNodeID().toString().toUpperCase());
					double[] d = posteriors.get(n);
					for (int i = 0; i < order.length; i++)
						order[i] = i;
					
					//double max = 0;
					//int maxIndex = 0;
					for (int i = 0; i < order.length; i++) {
						try {
							//System.out.print("\t" + d[order[i]]);
							fout.print("\t" + d[order[i]]);
							//if (d[order[i]] > max) {
							//	max = d[order[i]];
							//	maxIndex = order[i];
							//}
						}
						catch (ArrayIndexOutOfBoundsException e) {
							continue;
						}
					}
					//System.out.println("\t" + names.elementAt(maxIndex));
					//fout.println("\t" + names.elementAt(maxIndex));
					fout.print("\n");
				}
			}

			fout.close();
		}
		catch (Exception ioe) {
			System.err.println("PrintOutProteinList: " + outFileName + " " + ioe.toString());
			System.exit(1);
		}
	}

	//	 Quick check to see how many we got right.
	/*
	  private double prettyPrintValidation(Vector names, Hashtable posteriors) {
	// TODO Auto-generated method stub

	int[] order = new int[names.size()];
	int totalGO = 0;
	int correctGO = 0;
	int numIDA = 0;
	int numIDAp = 0;
	int numIDACorrect = 0;
	int numIDACorrectp = 0;
	int numIMP = 0;
	int numIMPp = 0;
	int numIMPCorrect = 0;
	int numIMPCorrectp = 0;
	int numIEA = 0;
	int numIEAp = 0;
	int numIEACorrect = 0;
	int numIEACorrectp = 0;
	int numNAS = 0;
	int numNASp = 0;
	int numNASCorrect = 0;
	int numNASCorrectp = 0;
	int numTAS = 0;
	int numTASp = 0;
	int numTASCorrect = 0;
	int numTASCorrectp = 0;
	int numISS = 0;
	int numISSp = 0;
	int numISSCorrect = 0;
	int numISSCorrectp = 0;
	int numIPI = 0;
	int numIPIp = 0;
	int numIPICorrect = 0;
	int numIPICorrectp = 0;
	int numIEP = 0;
	int numIEPp = 0;
	int numIEPCorrect = 0;
	int numIEPCorrectp = 0;
	int numND = 0;
	int numNDp = 0;
	int numNDCorrect = 0;
	int numNDCorrectp = 0;
	int numCorrect = 0;
	int totalP = 0;
	int totalPCorrect = 0;

	// PFam-type names (not swiss prot ids)
	Vector proteinNames = family.getProteinNames();
	for(int i = 0; i < proteinNames.size(); i++) {
	    //match protein with reconciled tree leaf
	    String currentProtein = (String)proteinNames.elementAt(i);
	    ProteinAnnotationObjectWithProbabilityMethods p = (ProteinAnnotationObjectWithProbabilityMethods)
		proteins.get(family.getProteinName(currentProtein));
	    ProbabilisticReconciledPhylogenyObject tree = family.getTree();
	    if(p != null && tree.hasNode(currentProtein)) {
		ProbabilisticReconciledPhylogenyObject.Node n = tree.getNode(currentProtein);
		if(posteriors.containsKey(n)) {
		    double[] d = (double[])posteriors.get(n);
		    for(int j = 0; j < order.length; j++) order[j] = j;
		    double max = 0;
		    int maxIndex = 0;
		    for(int j = 0; j < order.length; j++) {
	                      try {
	                          if(d[order[j]] > max) {
	                              max = d[order[j]];
	                              maxIndex = order[j];
	                          }
	                      } catch (ArrayIndexOutOfBoundsException e) {
	                          continue;
	                      }
		    }
		    if(p.getGOReal() != null) {
			if(p.getGOReal().contains(names.elementAt(maxIndex))){
			    correctGO++;
			} else if (verbose){
			    System.out.println("Missed "+(n.getNodeID())
					       +", predicted "
					       +names.elementAt(maxIndex)
					       +", real "+p.getGOReal());
			}
			totalGO++;
		    }
		    // All of the existing GO annotations
		    Vector goAnnot = p.getGONumber();
		    Vector MOC = p.getMOC();
		    if(goAnnot != null) {
			if(MOC.contains("IDA")) {
			    numIDAp++;
			    numIDA += countElementsInVector(MOC, "IDA");
			    numCorrect =
				countCorrectElementsInVector(MOC, "IDA",
							     goAnnot,
					     names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numIDACorrect += numCorrect;
				numIDACorrectp += 1;
			    }
			}
			if(MOC.contains("IMP")) {
			    numIMPp++;
			    numIMP += countElementsInVector(MOC, "IMP");
			    numCorrect =
				countCorrectElementsInVector(MOC, "IMP",
							     goAnnot,
						    names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numIMPCorrect += numCorrect;
				numIMPCorrectp += 1;
			    }
			}
			if(MOC.contains("IEA")) {
			    numIEAp++;
			    numIEA += countElementsInVector(MOC, "IEA");
			    numCorrect =
				countCorrectElementsInVector(MOC, "IEA",
							     goAnnot,
						    names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numIEACorrect += numCorrect;
				numIEACorrectp += 1;
			    }
			}
			if(MOC.contains("TAS")) {
			    numTASp++;
			    numTAS += countElementsInVector(MOC, "TAS");
			    numCorrect
				= countCorrectElementsInVector(MOC, "TAS",
							       goAnnot,
						    names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numTASCorrect += numCorrect;
				numTASCorrectp += 1;
			    }
			}
			if(MOC.contains("NAS")) {
			    numNASp++;
			    numNAS += countElementsInVector(MOC, "NAS");
			    numCorrect = countCorrectElementsInVector(MOC,
								      "NAS",
								      goAnnot,
						   names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numNASCorrect += numCorrect;
				numNASCorrectp += 1;
			    }
			}
			if(MOC.contains("ISS")) {
			    numISSp++;
			    numISS += countElementsInVector(MOC, "ISS");
			    numCorrect
				= countCorrectElementsInVector(MOC, "ISS",
							       goAnnot,
						    names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numISSCorrect += numCorrect;
				numISSCorrectp += 1;
			    }
			}
			if(MOC.contains("IPI")) {
			    numIPIp++;
			    numIPI += countElementsInVector(MOC, "IPI");
			    numCorrect =
				countCorrectElementsInVector(MOC, "IPI",
							     goAnnot,
						    names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numIPICorrect += numCorrect;
				numIPICorrectp += 1;
			    }
			}
			if(MOC.contains("IEP")) {
			    numIEPp++;
			    numIEP += countElementsInVector(MOC, "IEP");
			    numCorrect =
				countCorrectElementsInVector(MOC, "IEP",
							     goAnnot,
						    names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numIEPCorrect += numCorrect;
				numIEPCorrectp += 1;
			    }
			}
			if(MOC.contains("ND")) {
			    numNDp++;
			    numND += countElementsInVector(MOC, "ND");
			    numCorrect =
				countCorrectElementsInVector(MOC, "ND",
							     goAnnot,
						    names.elementAt(maxIndex));
			    if(numCorrect >= 1) {
				numNDCorrect += numCorrect;
				numNDCorrectp += 1;
			    }
			}
		    }
		}
	    } else {
		System.out.println("Couldnt get protein "
				   +family.getProteinName(currentProtein));
	    }
	}
	totalP = numIEAp+numTASp+numNASp+numISSp+numIPIp+numIEPp;
	totalPCorrect = numIEACorrectp + numTASCorrectp
	    + numNASCorrectp + numISSCorrectp+ numIPICorrectp+numIEPCorrectp;
	//System.out.println(names);
	if(verbose)
	    System.out.println("Total (GOReal): "
			       +((double)correctGO/(double)totalGO)
			       +" out of "+totalGO+" ("+(totalGO-correctGO)
			       +" wrong)");
	if(verbose)
	    System.out.println(proteins.size()+"\t"+numIDA+"\t"
			       +numIDAp+"\t"
			       +numIDACorrect+"\t"+numIDACorrectp+"\t"
			       +numIMP+"\t"+numIMPp+"\t"
			       +numIMPCorrect+"\t"+numIMPCorrectp+"\t"
			       +numIEA+"\t"+numIEAp+"\t"+numIEACorrect
			       +"\t"+numIEACorrectp+"\t"
			       +numTAS+"\t"+numTASp+"\t"
			       +numTASCorrect+"\t"+numTASCorrectp+"\t"
			       +numNAS+"\t"+numNASp+"\t"+numNASCorrect
			       +"\t"+numNASCorrectp+"\t"
			       +numISS+"\t"+numISSp
			       +"\t"+numISSCorrect+"\t"+numISSCorrectp+"\t"
			       +numIPI+"\t"+numIPIp+"\t"+numIPICorrect+"\t"
			       +numIPICorrectp+"\t"
			       +numIEP+"\t"+numIEPp+"\t"
			       +numIEPCorrect+"\t"+numIEPCorrectp+"\t"
			       +numND+"\t"+numNDp+"\t"
			       +numNDCorrect+"\t"+numNDCorrectp);
	if(totalP > 0)
	    if(verbose)
		System.out.println("Total correct: "
				   +((double)totalPCorrect/(double)totalP));
	if(totalP == 0) return 1.0;
	if(totalGO == 0) {
	    return((double)totalPCorrect/(double)totalP);
	}
	return ((double)correctGO/(double)totalGO);
	  }
	 */

	// Quick check to print out the details of the validation
	// proteins.
	private double prettyPrintProteinsValidation(Vector<Integer> names,
	    Hashtable<Node, double[]> posteriors) {
		int[] order = new int[names.size()];
		int totalGO = 0;
		int correctGO = 0;
		int numIDA = 0;
		int numIDAp = 0;
		int numIMP = 0;
		int numIMPp = 0;
		int numIEA = 0;
		int numIEAp = 0;
		int numIEACorrect = 0;
		int numIEACorrectp = 0;
		int numNAS = 0;
		int numNASp = 0;
		int numNASCorrect = 0;
		int numNASCorrectp = 0;
		int numTAS = 0;
		int numTASp = 0;
		int numTASCorrect = 0;
		int numTASCorrectp = 0;
		int numISS = 0;
		int numISSp = 0;
		int numISSCorrect = 0;
		int numISSCorrectp = 0;
		int numIPI = 0;
		int numIPIp = 0;
		int numIPICorrect = 0;
		int numIPICorrectp = 0;
		int numIEP = 0;
		int numIEPp = 0;
		int numIEPCorrect = 0;
		int numIEPCorrectp = 0;
		int numND = 0;
		int numNDp = 0;
		int numNDCorrect = 0;
		int numNDCorrectp = 0;
		int numCorrect = 0;
		int totalP = 0;
		int totalPCorrect = 0;

		// PFam-type names (not swiss prot ids)
		Vector<String> proteinNames = family.getProteinNames();
		System.out.println("***************************");

		for (int i = 0; i < proteinNames.size(); i++) {
			//match protein with reconciled tree leaf
			String currentProtein = (String)proteinNames.elementAt(i);
			ProteinAnnotationObjectWithProbabilityMethods p = (ProteinAnnotationObjectWithProbabilityMethods) proteins.get(family.getProteinName(currentProtein));
			ProbabilisticReconciledPhylogenyObject tree = family.getTree();

			if (p != null && tree.hasNode(currentProtein)) {
				ProbabilisticReconciledPhylogenyObject.Node n = tree.getNode(currentProtein);

				if (posteriors.containsKey(n)) {
					double[] d = posteriors.get(n);

					for (int j = 0; j < order.length; j++)
						order[j] = j;

					double max = 0;
					int maxIndex = 0;

					for (int j = 0; j < order.length; j++) {
						try {
							if (d[order[j]] > max) {
								max = d[order[j]];
								maxIndex = order[j];
							}
						}
						catch (ArrayIndexOutOfBoundsException e) {
							continue;
						}
					}

					if (p.getGOReal() != null) {
						if (p.getGOReal().contains(names.elementAt(maxIndex))) {
							correctGO++;
						}
						else
							if (verbose) {
								System.out.println("Missed " + (n.getNodeID())
								                   + ", predicted "
								                   + names.elementAt(maxIndex)
								                   + ", real " + p.getGOReal());
							}

						totalGO++;
					}

					// All of the existing GO annotations
					Vector<Integer> goAnnot = p.getGONumbers();
					Vector<String> MOC = p.getMocs();

					if (goAnnot != null) {
						if (MOC.contains("IDA")) {
							numIDAp++;
							numIDA += countElementsInVector(MOC, "IDA");
						}

						if (MOC.contains("IMP")) {
							numIMPp++;
							numIMP += countElementsInVector(MOC, "IMP");
						}

						if (goAnnot.size() > 0) {
							System.out.print(p.getName() + "\t");
							listPredictedElementsInVector(MOC, goAnnot,
							                              names.elementAt(maxIndex));
							System.out.print("\t" + max + "\t");
							listOtherElementsInVector(MOC, goAnnot,
							                          names.elementAt(maxIndex));
							System.out.println();
						}

						if (MOC.contains("IEA")) {
							numIEAp++;
							numIEA += countElementsInVector(MOC, "IEA");
							numCorrect =
							  countCorrectElementsInVector(MOC, "IEA",
							                               goAnnot,
							                               names.elementAt(maxIndex));

							if (numCorrect >= 1) {
								numIEACorrect += numCorrect;
								numIEACorrectp += 1;
							}
						}

						if (MOC.contains("TAS")) {
							numTASp++;
							numTAS += countElementsInVector(MOC, "TAS");
							numCorrect =
							  countCorrectElementsInVector(MOC, "TAS",
							                               goAnnot,
							                               names.elementAt(maxIndex));

							if (numCorrect >= 1) {
								numTASCorrect += numCorrect;
								numTASCorrectp += 1;
							}
						}

						if (MOC.contains("NAS")) {
							numNASp++;
							numNAS += countElementsInVector(MOC, "NAS");
							numCorrect =
							  countCorrectElementsInVector(MOC, "NAS",
							                               goAnnot,
							                               names.elementAt(maxIndex));

							if (numCorrect >= 1) {
								numNASCorrect += numCorrect;
								numNASCorrectp += 1;
							}
						}

						if (MOC.contains("ISS")) {
							numISSp++;
							numISS += countElementsInVector(MOC, "ISS");
							numCorrect = countCorrectElementsInVector(MOC,
							             "ISS",
							             goAnnot,
							             names.elementAt(maxIndex));

							if (numCorrect >= 1) {
								numISSCorrect += numCorrect;
								numISSCorrectp += 1;
							}
						}

						if (MOC.contains("IPI")) {
							numIPIp++;
							numIPI += countElementsInVector(MOC, "IPI");
							numCorrect =
							  countCorrectElementsInVector(MOC, "IPI",
							                               goAnnot,
							                               names.elementAt(maxIndex));

							if (numCorrect >= 1) {
								numIPICorrect += numCorrect;
								numIPICorrectp += 1;
							}
						}

						if (MOC.contains("IEP")) {
							numIEPp++;
							numIEP += countElementsInVector(MOC, "IEP");
							numCorrect =
							  countCorrectElementsInVector(MOC, "IEP",
							                               goAnnot,
							                               names.elementAt(maxIndex));

							if (numCorrect >= 1) {
								numIEPCorrect += numCorrect;
								numIEPCorrectp += 1;
							}
						}

						if (MOC.contains("ND")) {
							numNDp++;
							numND += countElementsInVector(MOC, "ND");
							numCorrect =
							  countCorrectElementsInVector(MOC, "ND",
							                               goAnnot,
							                               names.elementAt(maxIndex));

							if (numCorrect >= 1) {
								numNDCorrect += numCorrect;
								numNDCorrectp += 1;
							}
						}
					}
				}
			}
			else {
				System.out.println("Couldnt get protein "
				                   + family.getProteinName(currentProtein));
			}
		}

		totalP = numIEAp + numTASp + numNASp + numISSp + numIPIp + numIEPp;
		totalPCorrect = numIEACorrectp + numTASCorrectp + numNASCorrectp
		                + numISSCorrectp + numIPICorrectp + numIEPCorrectp;

		//System.out.println(names);
		if (verbose)
			System.out.println("Total (GOReal): " + ((double)correctGO
			                   / (double)totalGO)
			                   + " out of " + totalGO + " ("
			                   + (totalGO - correctGO) + " wrong)");
		
		System.out.println("Total (GOReal): " + ((double)correctGO
                / (double)totalGO)
                + " out of " + totalGO + " ("
                + (totalGO - correctGO) + " wrong)");

		if (verbose)
			System.out.println(proteins.size() + "\t" + numIDA + "\t"
			                   + numIDAp + "\t" + numIMP + "\t" + numIMPp + "\t"
			                   + numIEA + "\t" + numIEAp + "\t"
			                   + numIEACorrect + "\t" + numIEACorrectp + "\t"
			                   + numTAS + "\t" + numTASp + "\t"
			                   + numTASCorrect + "\t" + numTASCorrectp + "\t"
			                   + numNAS + "\t" + numNASp
			                   + "\t" + numNASCorrect + "\t" + numNASCorrectp + "\t"
			                   + numISS + "\t" + numISSp
			                   + "\t" + numISSCorrect + "\t" + numISSCorrectp + "\t"
			                   + numIPI + "\t" + numIPIp
			                   + "\t" + numIPICorrect + "\t" + numIPICorrectp + "\t"
			                   + numIEP + "\t" + numIEPp
			                   + "\t" + numIEPCorrect + "\t" + numIEPCorrectp + "\t"
			                   + numND + "\t" + numNDp
			                   + "\t" + numNDCorrect + "\t" + numNDCorrectp);

		if (totalP > 0)
			if (verbose)
				System.out.println("Total correct: "
				                   + ((double)totalPCorrect / (double)totalP));
				
				System.out.println("---------------------------------------------------");
				System.out.println("totalPCorrect: " + totalPCorrect + "\ttotalP: " + totalP);
				System.out.println("---------------------------------------------------");
				
				System.out.println("Total correct: "
		                   + ((double)totalPCorrect / (double)totalP));
				
		if (totalP == 0)
			return 1.0;

		return((double)totalPCorrect / (double)totalP);
	}

	/* Prints out leaf nodes in respective clusters of
	 * function prediction.
	 */
	/*
	  private void printClustersOfFunctions(Vector names,
					    Hashtable posteriors)
	  {
	int[] order = new int[names.size()];

	// PFam-type names (not swiss prot ids)
	Vector proteinNames = family.getProteinNames();
	ProbabilisticReconciledPhylogenyObject tree = family.getTree();
	System.out.println("***************************");
	for(int k = 0; k < names.size(); k++) {
	    String currentFn = (String)names.elementAt(k);
	    System.out.println("Function: "+currentFn);
	    for(int i = 0; i < proteinNames.size(); i++) {
		//match protein with reconciled tree leaf
		String currentProtein = (String)proteinNames.elementAt(i);
		ProteinAnnotationObjectWithProbabilityMethods p = (ProteinAnnotationObjectWithProbabilityMethods)
		    proteins.get(family.getProteinName(currentProtein));
		if(p != null && tree.hasNode(currentProtein)) {
		    ProbabilisticReconciledPhylogenyObject.Node n = tree.getNode(currentProtein);
		    if(n.isLeaf() && posteriors.containsKey(n)) {
			double[] d = (double[])posteriors.get(n);
			for(int j = 0; j < order.length; j++) order[j] = j;
			double max = 0;
			int maxIndex = 0;
			for(int j = 0; j < order.length; j++) {
			    try {
				if(d[order[j]] > max) {
				    max = d[order[j]];
				    maxIndex = order[j];
				}
			    } catch (ArrayIndexOutOfBoundsException e) {
				continue;
			    }
			}
			if(maxIndex == k) {
			    System.out.println(n.getNodeID()+"\t"+max);
			}
		    }
		}
	    }
	}
	}*/


	/**
	 * @param moc
	 * @param string
	 * @param goAnnot
	 * @param object
	 * @return
	 */
	private int countCorrectElementsInVector(Vector<String> moc, String method,
	    Vector<Integer> goAnnot,
	    Object annotation) {
		int total = 0;

		for (int i = 0; i < moc.size(); i++) {
			if (annotation.equals(goAnnot.elementAt(i))) {
				if (method.equals(moc.elementAt(i))) {
					total++;
				}
			}
		}

		return total;
	}

	private void listPredictedElementsInVector(Vector<String> moc, Vector<Integer> goAnnot,
	    Integer annotation) {
		for (int i = 0; i < moc.size(); i++) {
			if (((String)moc.elementAt(i)).equals("IDA")) {
				System.out.print(goAnnot.elementAt(i) + "\tIDA\t");
				return;
			}

			if (((String)moc.elementAt(i)).equals("IMP")) {
				System.out.print(goAnnot.elementAt(i) + "\tIMP\t");
				return;
			}
		}

		for (int i = 0; i < moc.size(); i++) {
			if (annotation.equals(goAnnot.elementAt(i))) {
				System.out.print(goAnnot.elementAt(i) + "\t"
				                 + moc.elementAt(i) + "\t");
				return;
			}
		}

		return;
	}

	private void listOtherElementsInVector(Vector<String> moc, Vector<Integer> goAnnot,
	                                       Integer annotation) {
		for (int i = 0; i < moc.size(); i++) {
			if (((String)moc.elementAt(i)).equals("IDA")) {
			}
			else
				if (((String)moc.elementAt(i)).equals("IMP")) {
				}
				else
					if (!annotation.equals(goAnnot.elementAt(i))) {
						System.out.print(goAnnot.elementAt(i) + "\t"
						                 + moc.elementAt(i) + "\t");
					}
		}

		return;
	}

	/*
	  private boolean correctElementsInVector(Vector moc, Vector goAnnot,
					    Object annotation) {
	if(moc.contains("IDA") || moc.contains("IMP")) return true;
	for(int i = 0; i < moc.size(); i++){
	    if(annotation.equals(goAnnot.elementAt(i))) {
		return true;
	    }
	}
	return false;
	}*/

	/**
	 * @param moc
	 * @param string
	 * @return
	 */
	private int countElementsInVector(Vector<String> moc, String method) {
		int total = 0;

		for (int i = 0; i < moc.size(); i++) {
			if (method.equals(moc.elementAt(i))) {
				total++;
			}
		}

		return total;
	}

	/**
	 * @param goAnnot
	 * @param moc
	 * @param names
	 */
	/*
	  private void pruneGOAnnotations(Vector goAnnot, Vector moc, Vector names)
	  {
	Vector toRemove = new Vector();
	for(int i = 0; i < goAnnot.size(); i++) {
	    if(!names.contains(goAnnot.elementAt(i))) {
		toRemove.add(0, new Integer(i));
	    }
	    else if((moc.contains("IDA") || moc.contains("IMP"))
		    && !(moc.elementAt(i).equals("IDA")
			 || moc.elementAt(i).equals("IMP"))) {
		toRemove.add(0, new Integer(i));
	    }
	}
	for(int i = 0; i < toRemove.size(); i++) {
	    goAnnot.remove(((Integer)toRemove.elementAt(i)).intValue());
	    moc.remove(((Integer)toRemove.elementAt(i)).intValue());
	}
	}*/


	/** Should we ever want to use an already existing SifterExecutionSettingsContainer object,
	 * we can direct our pointer to that.
	 * @param settings SifterExecutionSettingsContainer object to set to.
	 * @see Sifter.main()
	 */
	public void setSettingsObject(SifterExecutionSettingsContainer settings) {
		this.settings = settings;
		initSettings();
	}

	/** Once we get a new settings object, initialize
	 * any variables that we care about.
	 * For example, I know ExpectationMaximizationObject cares about verbose, family name, etc.
	 */
	private void initSettings() {
		this.verbose =
		  ((Boolean)this.settings.getSetting("verbose")).booleanValue();
	}


}

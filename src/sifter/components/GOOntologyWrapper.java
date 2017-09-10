/**
 * Data structure for the GO (gene ontology) DAG.
 *
 * Copyright 2003-2005 Barbara Engelhardt (bee@cs.berkeley.edu)
 * @author Barbara Engelhardt
 */

package sifter.components;

import java.util.*;

import sifter.components.ProbabilisticDAG.Node;

//import org.apache.regexp.*;
import java.io.*;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public class GOOntologyWrapper {

	// Denotes whether these members of the proteins
	// should be pulled out of the SwissProt DB.
	//private static boolean nameP;
	private String goFileName;

	private Vector<Entry> functionList; // list of most recent pnums
	//parents at each level.
	public ProbabilisticDAG probabilisticDAG;   // DAG of functions
	private int lineno;
	private Entry currentEntry;
	//private static ProteinAnnotationObjectWithProbabilityMethods currentNode;
	private double rValue;
	private Vector<Integer> leafList;

	private SifterExecutionSettingsContainer settings;
	private boolean verbose;

	public GOOntologyWrapper(String fileName, boolean pName, boolean pGo, boolean pMoc) {
		functionList = new Vector<Entry>();
		probabilisticDAG = new ProbabilisticDAG();
		lineno = 0;
		goFileName = fileName;
		rValue = (double) - 1.0;
		leafList = null;
	}

	public void parseOntologyLine(String s) {
		// ProbabilisticDAG functions: seems to add it stack-like to
		// "functionList" and then
		// really really add it to ProbabilisticDAG.
		int level = -1;
		int level2;
		level = s.indexOf('%') - 1;
		level2 = s.indexOf('<') - 1;

		if ((level2 > level && (level2 >= 0 && level < 0)) || (level2 < level && level2 >= 0 && level >= 0))
			level = level2;

		if (level < 0) {
			if (verbose) System.out.println("Error: no percents or" + " less thans in GO ontology line: "+lineno);
			System.out.println(s);
			return;
		}

		s = s.substring(level + 2, s.length());
		s = pullOneEntry(s, level);
		
		Entry e = currentEntry;
		
		// add parent as previous level
		Entry p = null;
		if (functionList.size() >= level && level - 1 >= 0)
			p = functionList.get(level - 1);
		if (p != null)
			probabilisticDAG.addParent(e.getNumber(), p.getNumber());

		// add as new level
		functionList.add(level, e);

		// if any more %, create new, add as parent.
		while (s.indexOf('%') >= 0 || s.indexOf('<') >= 0) {
			if (s.indexOf('%') >= 0) {
				s = pullOneEntry(s.substring(s.indexOf('%') + 2), level - 1);
				p = currentEntry;
			} else {
				s = pullOneEntry(s.substring(s.indexOf('<') + 2), level - 1);
				p = currentEntry;
			}

			if (p != null)
				probabilisticDAG.addParent(e.getNumber(), p.getNumber());
		}
	}

	private String pullOneEntry(String s, int level) {
		Entry e = null;
		int semicolon = s.indexOf(';');
		String name = s.substring(0, semicolon - 1);
		s = s.substring(semicolon + 1, s.length());
		semicolon = s.indexOf(':');
		int id = new Integer(s.substring(semicolon + 1, semicolon + 8)).intValue();
		s = s.substring(semicolon + 8);
		// not parsing EC numbers here, but could!
		// check if its already been seen.
		e = new Entry(id, name, level);

		if (!probabilisticDAG.contains(e.getNumber())) {
			probabilisticDAG.addNode(e.getName(), e.getNumber());
		}

		probabilisticDAG.addLevel(e.getNumber(), level);

		if (e != null)
			currentEntry = e;

		return s;
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked", "unused" })
	public void readGOFile() throws ClassNotFoundException {
		// http://code.google.com/p/variationtoolkit/wiki/GeneOntologyDbManager
		// https://bitbucket.org/xerial/sqlite-jdbc
		// http://sqlitebrowser.sourceforge.net/development.html
		Class.forName("org.sqlite.JDBC");
		
		try {
			Connection connection = null;
			connection = DriverManager.getConnection("jdbc:sqlite:"+goFileName);
			Statement statement = connection.createStatement();
			statement.setQueryTimeout(30);  // set timeout to 30 sec.
			
			Stack term_stack = new Stack();
			Stack level_stack = new Stack();
			
			String topElem = "unused name ; GO:0003674"; // MF ontology top node is GO:0003674.
			int topLevel = 0;
			
			// Add initial element			
			level_stack.push(new Integer(topLevel));
			pullOneEntry(topElem, topLevel);
			term_stack.push(currentEntry);
			
			int c = 0;
			
			while (!term_stack.empty()) {
				// Get "level" of this element
				Entry p = (Entry)term_stack.pop();
				int pLevel = ((Integer)level_stack.pop()).intValue();
				functionList.add(pLevel, p);
				
				// add as new level
				
				ResultSet rs = statement.executeQuery("select acn,rel,target from TERM2REL where target=\"GO:"+String.format("%07d", p.goNum)+"\" and rel=\"is_a\"");
				
				while (rs.next()) {
					c += 1;
					// Each rs.getString(1) is child to parent rs.getString(3).
					// System.out.println(rs.getString(1)+" is child to parent: "+rs.getString(3));
					pullOneEntry("unused name ; "+rs.getString(1), pLevel + 1);
					term_stack.push(currentEntry); 
					level_stack.push(pLevel + 1);
					
					probabilisticDAG.addParent(currentEntry.getNumber(), p.getNumber());
					
					//if (c % 100 == 0 )
					//	System.out.println("Slow");
				}
			}
			
		} catch (SQLException e) {
			System.err.println(e.getMessage());
		}
	}
	
	// Often called by SifterPipelineObject.java's buildDatasetGODAG().
	// Basically, it takes the .ont file and puts it in its own data structure.
	public void readGOFileOld() {
		String str = null;
		BufferedReader fin;
		
		try {
			fin = new BufferedReader(new FileReader(goFileName));

			try {
				// Read each line of the .ont Gene Ontology filename
				while ((str = fin.readLine()) != null) {
					lineno++;

					// If the line isn't blank and it's not a comment,
					if (str.length() > 0 && !str.startsWith("!")) {
						parseOntologyLine(str);
					}
				}
			}
			catch (Exception e) {
				System.out.println("Error in readGOFile: " + str + " " + lineno + " " + e.getMessage());
				System.exit(1);
			}

			fin.close();
		}
		catch (Exception ioe) {
			// TODO: Why is this exception generalized?
			// I think it's for FileNotFoundError.
			System.err.println("readGOFile: " + goFileName + " " +
			                   ioe.getMessage());
			// TODO: Make error message friendlier.
			System.exit(1); // TODO: Have this propagate upwards.
			//In a GUI, I can't have this completely bring down the program.
		}

		//functions.printOutDAG();
	}


	//////////////////////////////////////////////////////////
	//  Statistics-Gathering section                        //
	//////////////////////////////////////////////////////////

	// Sent a vector with functions, a vector with experimental methods,
	// this function tallies all the counts from the induced subgraph
	// of each of the functions, and associated counts for each
	// experimental method.
//    public void tallyObservedFunctions(Vector<Integer> fns, Vector<String> methods)
//    {
//        String m;
//        int fn;
//        for(int i = 0; i < fns.size(); i++) {
//            fn = ((Integer)fns.elementAt(i)).intValue();
//            m = null;
//            if(methods.size() > i) m = (String)methods.get(i);
//            if(m == null) {
//                if (verbose) System.out.println("Function without a method: "
//						+ fn);
//            }
//            tallyOneFunction(fn, m);
//        }
//    }

//    private void tallyOneFunction(int fn, String m)
//    {
//        Vector<Integer> isg = new Vector<Integer>();
//        isg = functions.getInducedSubgraph(fn, isg);
//        for(int i = 0; i < isg.size(); i++) {
//            fn = isg.elementAt(i).intValue();
//            if(!functions.contains(fn)) { // should this be containskey?
//                if (verbose)
//		    System.out.println("FunctionNode does not exist: " + fn);
//            }
//            Vector<Integer> ex = ((Vector<Integer>)functions.getAdditional(fn));
//            //System.out.println(fn);
//            if(ex == null) {
//                ex = new Vector<Integer>(15);
//                for(int j = 0; j < 15; j++) ex.add(new Integer(0));
//            }
//            if(m.equals("IEA"))
//		ex.set(0,new Integer(ex.elementAt(0).intValue()+1));
//            else if(m.equals("IMP"))
//		ex.set(1,new Integer(ex.elementAt(1).intValue()+1));
//            else if(m.equals("IGI"))
//		ex.set(2,new Integer(ex.elementAt(2).intValue()+1));
//            else if(m.equals("IPI"))
//		ex.set(3,new Integer(ex.elementAt(3).intValue()+1));
//            else if(m.equals("ISS"))
//		ex.set(4,new Integer(ex.elementAt(4).intValue()+1));
//            else if(m.equals("IDA"))
//		ex.set(5,new Integer(ex.elementAt(5).intValue()+1));
//            else if(m.equals("IEP"))
//		ex.set(6,new Integer(ex.elementAt(6).intValue()+1));
//            else if(m.equals("IEA"))
//		ex.set(7,new Integer(ex.elementAt(7).intValue()+1));
//            else if(m.equals("TAS"))
//		ex.set(8,new Integer(ex.elementAt(8).intValue()+1));
//            else if(m.equals("NAS"))
//		ex.set(9,new Integer(ex.elementAt(9).intValue()+1));
//            else if(m.equals("RCA"))
//		ex.set(10,new Integer(ex.elementAt(10).intValue()+1));
//            else if(m.equals("IGC"))
//		ex.set(11,new Integer(ex.elementAt(11).intValue()+1));
//            else if(m.equals("ND"))
//		ex.set(12,
//		       new Integer(ex.elementAt(12).intValue()+1));
//            else if(m.equals("IC"))
//		ex.set(13,
//		       new Integer(ex.elementAt(13).intValue()+1));
//            else if(m.equals("P"))
//		ex.set(14,
//		       new Integer(ex.elementAt(14).intValue()+1));
//            else if(m.equals("E"))
//		ex.set(15,
//		       new Integer(ex.elementAt(15).intValue()+1));
//            else if(m.equals("NR"))
//		ex.set(16,
//		       new Integer(ex.elementAt(16).intValue()+1));
//            else if(m.equals("GOR"))
//		ex.set(17,
//		       new Integer(ex.elementAt(17).intValue()+1));
//            else
//                if (verbose)
//		    System.out.println("Methods not recognized: "+m);
//            functions.addAdditional(fn, ex.clone());
//        }
//    }

	// prints this out into a nice, matlab-readable file
	@SuppressWarnings("rawtypes")
	public void printOutGOTally(String fileName) {
		PrintStream fout;
		fileName = fileName + ".pdt";

		try {
			fout =  new PrintStream(new FileOutputStream(new File(fileName)));

			Enumeration<Integer> pnames = probabilisticDAG.keys();

			while (pnames.hasMoreElements()) {
				int pfName = pnames.nextElement().intValue();
				int total = 0;

				if (!probabilisticDAG.contains(pfName)) {
					if (verbose)
						System.out.println("FunctionNode does not exist: "
						                   + pfName);
				}

				int level = probabilisticDAG.getLevel(pfName);
				int pCount = probabilisticDAG.getNumParents(pfName);
				int cCount = probabilisticDAG.getNumChildren(pfName);
				Vector ex = (Vector)probabilisticDAG.getAdditional(pfName);

				if (ex == null) {
					fout.println(pfName + " " + level + " " + pCount + " "
					             + cCount + " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
				}
				else {
					fout.print(pfName + " " + level + " " + pCount + " " + cCount + " ");

					for (int i = 0; i < 15; i++) {
						fout.print(ex.elementAt(i) + " ");
						total += ((Integer)ex.elementAt(i)).intValue();
					}

					fout.println(total);
				}
			}

			fout.close();
		}
		catch (Exception ioe) {
			System.err.println("PrintOutGOTally: " + fileName + " " +
			                   ioe.getMessage());
			System.exit(1);
		}
	}


	//////////////////////////////////////////////////
	// Ontology building section                    //
	//////////////////////////////////////////////////

	// REFACTORME: unused
//    public void tallyInducedFunctions(Vector fns, Vector ms)
//    {
//        String m;
//        int fn;
//        for(int i = 0; i < fns.size(); i++) {
//            fn = ((Integer)fns.elementAt(i)).intValue();
//            m = null;
//            if(ms.size() > i) m = (String)ms.get(i);
//            if(m == null) {
//                if (verbose)
//		    System.out.println("Function without a method: "+fn);
//            }
//            Vector<Integer> isg = functions.getInducedSubgraph(fn);
//            Vector<Integer> dg = functions.getDescendantGraph(fn, false);
//            isg.addAll(dg);
//            for(int j = 0; j < isg.size(); j++) {
//                fn = ((Integer)isg.elementAt(j)).intValue();
//                if(!functions.contains(fn)) { // should this be containskey?
//                    if (verbose)
//			System.out.println("FunctionNode does not exist: "+fn);
//                }
//                Integer ex = (Integer)functions.getAdditional(fn);
//                //System.out.println(fn);
//                if(ex == null) {
//                    ex = new Integer(1);
//                } else {
//                    ex = new Integer(ex.intValue()+1);
//                }
//                functions.addAdditional(fn, ex);
//            }
//        }
//    }

	// Only ever called from SifterPipelineObject.buildDatasetGODAG().
	// fns - vectorlist of all of the GO function codes from the XML file.
	// ms - vectorlist of all of the methods from the XML file.
	public void tallyFunctions(Vector<Integer> fns, Vector<String> ms) {
		String m;
		int fn;

		for (int i = 0; i < fns.size(); i++) {
			// fn - function
			fn = ((Integer)fns.elementAt(i)).intValue();
			// m  - method
			m = null;

			if (ms.size() > i) {
				m = (String)ms.get(i);
			}

			if (m == null) {
				if (verbose)
					System.out.println("Function without a method: " + fn);
			}
			else
				if (hasProperBackgroundSettings(m)) {
					Integer ex = (Integer)probabilisticDAG.getAdditional(fn);

					//if (verbose) System.out.println("Adding function: "+fn);
					if (ex == null) {
						// add ...
						ex = new Integer(1);
					}
					else {
						ex = new Integer(ex.intValue() + 1);
					}

					// add this additional function, this is a count
					probabilisticDAG.addAdditional(fn, ex);
					//System.out.print(ex + " ");
				}

			// so it goes through and puts a count on the GO DAG
			// which proteins are used.
		}
	}

	private boolean hasProperBackgroundSettings(String m) {
		if (((Boolean)settings.getSetting("bg")).booleanValue()) {
			System.err.println("bg is not updated / implemented for new experimental evidence codes in hasProperBackgroundSettings");
			System.exit(1);
		}

		// jeff: it used to do this. (see Sifter::buildOptions for details on the bg option)
//    	setp = ((Boolean)settings.getSetting("nas")).booleanValue();
//    	if(setp && m.equalsIgnoreCase("NAS") && !setb) return true;
//    	setp = ((Boolean)settings.getSetting("igi")).booleanValue();
//    	if(setp && m.equalsIgnoreCase("IGI")) return true;


		return hasProperSettings(m);
	}

	private boolean hasProperSettings(String m) {
		
		if (isNumeric(m))
			return ((Boolean)settings.getSetting("floats")).booleanValue();

		
		String moc_lowercase = m.toLowerCase();
		return ((Boolean)settings.getSetting(moc_lowercase)).booleanValue();
	}

	public void tallyFunctions(Vector<Integer> fns) {
		int fn;

		for (int i = 0; i < fns.size(); i++) {
			fn = ((Integer)fns.elementAt(i)).intValue();
			Integer ex = (Integer)probabilisticDAG.getAdditional(fn);

			//System.out.println(fn);
			if (ex == null) {
				ex = new Integer(1);
			}
			else {
				ex = new Integer(ex.intValue() + 1);
			}

			probabilisticDAG.addAdditional(fn, ex);
		}
	}

	public void tallyInducedFunctions(Vector<Integer> fns) {
		int fn;

		for (int i = 0; i < fns.size(); i++) {
			fn = ((Integer)fns.elementAt(i)).intValue();
			Vector<Integer> isg = probabilisticDAG.getInducedSubgraph(fn);
			Vector<Integer> dg = probabilisticDAG.getDescendantGraph(fn, false);
			isg.addAll(dg);

			for (int j = 0; j < isg.size(); j++) {
				fn = ((Integer)isg.elementAt(j)).intValue();

				if (!probabilisticDAG.contains(fn)) { // should this be containskey?
					if (verbose)
						System.out.println("FunctionNode does not exist: " + fn);
				}

				Integer ex = (Integer)probabilisticDAG.getAdditional(fn);

				//System.out.println(fn);
				if (ex == null) {
					ex = new Integer(1);
				}
				else {
					ex = new Integer(ex.intValue() + 1);
				}

				probabilisticDAG.addAdditional(fn, ex);
			}
		}
	}

	public void pruneZeroHits() {
		Enumeration<Integer> pnames = probabilisticDAG.keys();

		while (pnames.hasMoreElements()) {
			int pfName = pnames.nextElement().intValue();

			if (!probabilisticDAG.contains(pfName)) {
				if (verbose)
					System.out.println("FunctionNode does not exist: " + pfName);
			}

			Integer ex = (Integer)probabilisticDAG.getAdditional(pfName);

			if (ex == null) {
				probabilisticDAG.removeNode(pfName);
			}
		}
	}

	// jeff: remove leaves from the GO dag if there it does not appear directly in the PLI
	// jeff shouldn't we propagate downward first?
	// jeff current leaves might not have evidence, but have smaller depth
	public void pruneZeroHitsAndLeaves() {
		Vector<Integer> leaves = probabilisticDAG.getAllLeaves();

		for (int i = 0; i < leaves.size(); i++) {
			int pfName = leaves.elementAt(i).intValue();

			if (!probabilisticDAG.contains(pfName)) {
				if (verbose)
					System.out.println("FunctionNode does not exist: " + pfName);
			}

			Integer ex = (Integer)probabilisticDAG.getAdditional(pfName);

			if (ex == null) {
				Vector<Node> p = probabilisticDAG
				                 .getParents(leaves.elementAt(i).intValue());
				probabilisticDAG.removeNode(pfName);

				if (p != null) {
					for (int j = 0; j < p.size(); j++) {
						int node = ((ProbabilisticDAG.Node)p.elementAt(j)).getId();

						if (node > -1 && probabilisticDAG.getNumChildren(node) < 1) {
							leaves.add(new Integer(node));
						}
					}
				}
			}
		}

		// If there's only one function, add another
		// so that inference proceeds normally.
		if (leaves.size() == 1) {
			probabilisticDAG.addNode("NOFUNCTION", 0);
			// REFACTORME: this was found to be broken after templatizing
//	    Vector p = ((ProbabilisticDAG.Node)leaves.elementAt(0)).getParents();
//	    if(p != null && p.size() > 0) {
//		functions.addParent(0, ((Integer)p.elementAt(0)).intValue());
//	    }
		}
	}

	/* If the family has a singleton leaf in the
	 * GO DAG, add a new "other" leaf as a sibling
	 * of the singleton leaf, run normally */
	public void padSingletonLeaves() {
		if (probabilisticDAG.getNumLeaves() == 0) {
			System.out.println("************** Error **************");
			System.out.println("No ontology terms matching those found in the ontology file. Make sure that (a) ontology type matches and (b) .pli file contains at least one ontology term");
			System.exit(1);
		}

		if (probabilisticDAG.getNumLeaves() > 1)
			return;
		
	}

	// prints this out into a nice, matlab-readable file
	public void printOutInducedTally(String fileName) {
		PrintStream fout;
		fileName = fileName + ".pdt";
		
		try {
			fout =  new PrintStream(new FileOutputStream(new File(fileName)));
			Enumeration<Integer> pnames = probabilisticDAG.keys();

			while (pnames.hasMoreElements()) {
				int pfName = pnames.nextElement().intValue();

				if (!probabilisticDAG.contains(pfName)) {
					if (verbose)
						System.out.println("FunctionNode does not exist: "
						                   + pfName);
				}

				int level = probabilisticDAG.getLevel(pfName);
				int pCount = probabilisticDAG.getNumParents(pfName);
				int cCount = probabilisticDAG.getNumChildren(pfName);
				Integer ex = (Integer)probabilisticDAG.getAdditional(pfName);

				if (ex == null) {
					fout.println(pfName + " " + level + " " + pCount + " " + cCount + " 0");
				}
				else {
					fout.println(pfName + " " + level + " "
					             + pCount + " " + cCount + " " + ex);
				}
			}

			fout.close();
		}
		catch (Exception ioe) {
			System.err.println("PrintOutGOTally: " + fileName + " " +
			                   ioe.getMessage());
			System.exit(1);
		}
	}

	// prints this out into a nice, matlab-readable file
	public void printOutInducedLeafTally(String fileName) {
		PrintStream fout;
		fileName = fileName + ".pdt";

		try {
			fout =  new PrintStream(new FileOutputStream(new File(fileName)));
			Vector<Integer> leaves = probabilisticDAG.getAllLeaves();

			for (int i = 0; i < leaves.size(); i++) {
				int pfName = leaves.elementAt(i).intValue();

				if (!probabilisticDAG.contains(pfName)) {
					if (verbose)
						System.out.println("FunctionNode does not exist: "
						                   + pfName);
				}

				int level = probabilisticDAG.getLevel(pfName);
				int pCount = probabilisticDAG.getNumParents(pfName);
				double prior = probabilisticDAG.getPrior(pfName);
				double posterior = probabilisticDAG.getPosterior(pfName);
				double likelihood = probabilisticDAG.getLikelihood(pfName);
				fout.println(pfName + " " + i + " " + level + " "
				             + pCount + " " + prior + " " + posterior + " " + likelihood);
			}

			fout.close();
		}
		catch (Exception ioe) {
			System.err.println("PrintOutGOTally: " + fileName + " " +
			                   ioe.getMessage());
			System.exit(1);
		}
	}


	// prints this out into a nice, matlab-readable file
	public void printOutInducedProbabilities(String fileName) {
		PrintStream fout;
		fileName = fileName + ".pdt";

		try {
			fout =  new PrintStream(new FileOutputStream(new File(fileName)));
			Enumeration<Integer> pnames = probabilisticDAG.keys();

			while (pnames.hasMoreElements()) {
				int pfName = pnames.nextElement().intValue();

				if (!probabilisticDAG.contains(pfName)) {
					if (verbose)
						System.out.println("FunctionNode does not exist: "
						                   + pfName);
				}

				int level = probabilisticDAG.getLevel(pfName);
				int pCount = probabilisticDAG.getNumParents(pfName);
				int cCount = probabilisticDAG.getNumChildren(pfName);
				double prior = probabilisticDAG.getPrior(pfName);
				double posterior = probabilisticDAG.getPosterior(pfName);
				double likelihood = probabilisticDAG.getLikelihood(pfName);
				double evidence = probabilisticDAG.getEvidence(pfName);
				int leaves = probabilisticDAG.getNumDescendantLeaves(pfName);
				fout.println(pfName + " " + level + " "
				             + pCount + " " + cCount + " " + " "
				             + leaves + " " + prior + " "
				             + likelihood + " " + posterior + " " + evidence);
			}

			fout.close();
		}
		catch (Exception ioe) {
			System.err.println("PrintOutGOTally: " + fileName + " " +
			                   ioe.getMessage());
			System.exit(1);
		}
	}

	// calculates r value from paper.
	// According to Matlab, this calculation gets flakey
	// if leaves < 20. Anything above that and it is
	// accurate up to approx 5 digits.
	public void findRValue() {
		// The additional one is the "sluff" variable
		// connected to all nodes (counts as single extra leaf).
		int leaves = probabilisticDAG.getNumLeaves();
		double r = 1;

		for (int i = 2; i <= leaves; i++) {
			double t = (double) 2.2605 * (double)((i * i) - i);
			r = r + ((double)1.0 / t);
		}

		rValue = r * (double)leaves;

		if (verbose)
			System.out.println("Multiplier: " + r
			                   + ", Leaves: " + leaves
			                   + ", R-value: " + rValue);

		Vector<Integer> tempLeaves = probabilisticDAG.getAllLeaves();

		for (int i = 0; i < tempLeaves.size(); i++) {
			if (verbose)
				System.out.println("Leaf: " + tempLeaves.elementAt(i));
		}

		if (probabilisticDAG != null)
			probabilisticDAG.setRValue(rValue);
	}

	// Accessor function
	public double getRValue() {
		if (rValue < 0)
			findRValue();

		return rValue;
	}

	public double getPriorProbability(int fn) {
		return probabilisticDAG.getPrior(fn);
	}

	// Given a DAG, take all evidence from pli
	// with annotations for specific proteins. put all annotations into leaves.
	public void incorporateEvidenceShort(Vector<Integer> fs, Vector<String> ms, Vector<Double> confidences, boolean noIEA) {
		// First put actual evidence in the appropriate positions.
		for (int i = 0; i < fs.size(); i++) {
			int pfName = ((Integer)fs.elementAt(i)).intValue();
			String method = (String)ms.elementAt(i);

			if (probabilisticDAG.contains(pfName)
			    && hasProperSettings(method)) {
				double evidence = convertEvidence(method);
				evidence *= confidences.elementAt(i);

				if (verbose)
					System.out.println("Adding evidence: "
					                   + evidence + " to " + pfName);

				probabilisticDAG.addEvidence(pfName, evidence);
				probabilisticDAG.propagateEvidenceDownward(pfName, evidence);
			}
		}

		synchronizeLikelihoods();
		probabilisticDAG.aPrioriEvidence();
	}

	public Hashtable<Integer, ProbabilisticDAGNode> switchOutProtein() {
		Hashtable<Integer, ProbabilisticDAGNode> pi = probabilisticDAG.removeProbabilityInstance();
		probabilisticDAG.zeroProbabilities();
		return pi;
	}

	public void clearDAGProbabilities() {
		probabilisticDAG.zeroProbabilities();
	}

	public double[] pullOutLeafLikelihoods() {
		Vector<Integer> leaves = getAllLeaves();
		double[] probs = new double[leaves.size()];

		for (int i = 0; i < leaves.size(); i++) {
			probs[i] = probabilisticDAG
			           .getLikelihood(leaves.elementAt(i).intValue());
		}

		return probs;
	}

	/* sets the likelihood for leaves with evidence */
	public void synchronizeLikelihoods() {
		Vector<Integer> leaves = getAllLeaves();
		double[] probs = new double[leaves.size()];
		double prior = getLeafSubsetPrior();
		double ltemp = 0.0;
		double totall = 0.0;

		//double remainder = 0.0;
		for (int i = 0; i < leaves.size(); i++) {
			ltemp = probabilisticDAG
			        .getLikelihood(leaves.elementAt(i).intValue());

			if (ltemp > 0)
				totall = (1.0 - (1.0 - totall) * (1.0 - ltemp));
		}

		//remainder = getNumLeaves()*prior - totall;
		//if(remainder > 0) {
		//   remainder = remainder/(double)leaves.size();

		prior = (1 - totall) * prior;

		for (int i = 0; i < leaves.size(); i++) {
			ltemp = probabilisticDAG
			        .getLikelihood(leaves.elementAt(i).intValue());
			probs[i] = ltemp;
			probs[i] = 1.0 - ((1.0 - probs[i]) * (1.0 - prior));
			probabilisticDAG.setLikelihood(leaves.elementAt(i).intValue(),
			                        probs[i]);
		}
	}

	public void insertLeafLikelihoods(double[] ll) {
		Vector<Integer> leaves = getAllLeaves();

		for (int i = 0; i < ll.length - 1; i++) {
			probabilisticDAG.setLikelihood(leaves.elementAt(i).intValue(),
			                        ll[i]);
		}
	}
	
	@SuppressWarnings("unused")
	public static boolean isNumeric(String str)  
	{  
	  try  
	  {  
	    double d = Double.parseDouble(str);  
	  }  
	  catch(NumberFormatException nfe)  
	  {  
	    return false;  
	  }  
	  return true;  
	}

	private double convertEvidence(String m) {
		// This stuff shouldn't be hard-coded this way
		// GOA:
		
		double ec;
		if (isNumeric(m))
		{
			ec=Double.parseDouble(m);
			System.out.format("The evidence code is %f \n",ec);	
			return ec;
			
		}
		
		if (m.equalsIgnoreCase("exp")) return 0.9;
		if (m.equalsIgnoreCase("ida")) return 0.9;
		if (m.equalsIgnoreCase("ipi")) return 0.8;
		if (m.equalsIgnoreCase("imp")) return 0.8;
		if (m.equalsIgnoreCase("igi")) return 0.8;
		if (m.equalsIgnoreCase("iep")) return 0.4;
		
		if (m.equalsIgnoreCase("tas")) return 0.9;
		if (m.equalsIgnoreCase("nas")) return 0.3;
		
		if (m.equalsIgnoreCase("iss")) return 0.4;
		if (m.equalsIgnoreCase("iso")) return 0.4;
		if (m.equalsIgnoreCase("isa")) return 0.4;
		if (m.equalsIgnoreCase("ism")) return 0.4;
		if (m.equalsIgnoreCase("igc")) return 0.4;
		if (m.equalsIgnoreCase("iba")) return 0.4;
		if (m.equalsIgnoreCase("ibd")) return 0.4;
		if (m.equalsIgnoreCase("ikr")) return 0.4;
		if (m.equalsIgnoreCase("ird")) return 0.4;
		if (m.equalsIgnoreCase("rca")) return 0.4;
		
		if (m.equalsIgnoreCase("ic")) return 0.4;
		if (m.equalsIgnoreCase("nd")) return 0.3;
		
		if (m.equalsIgnoreCase("iea")) return 0.04;
		
		if (m.equalsIgnoreCase("nr")) return 0.3;
		
		// Barbara Inventions
		if (m.equals("P")) return 0.2; // no idea what p is
		if (m.equals("E")) return 0.6; // no idea what e is (evidence)
		if (m.equals("GOR")) return 0.9; // from reading papers, my entry
		
		// Jeff's crap:
		if (m.equalsIgnoreCase("gene3d")) return 0.002;
		if (m.equalsIgnoreCase("hamap")) return 0.02;
		if (m.equalsIgnoreCase("panther")) return 0.001;
		if (m.equalsIgnoreCase("pfam")) return 0.0002;
		if (m.equalsIgnoreCase("pirsf")) return 0.00002;
		if (m.equalsIgnoreCase("prints")) return 0.2;
		if (m.equalsIgnoreCase("prodom")) return 0.0;
		if (m.equalsIgnoreCase("prosite_patterns")) return 0.2;
		if (m.equalsIgnoreCase("prosite_profiles")) return 0.002;
		if (m.equalsIgnoreCase("smart")) return 0.0006;
		if (m.equalsIgnoreCase("superfamily")) return 0.0;
		if (m.equalsIgnoreCase("tigrfams")) return 0.0;
		if (m.equalsIgnoreCase("genemania")) return 0.2;
		if (m.equalsIgnoreCase("pagosub")) return 0.01;
		if (m.equalsIgnoreCase("protfun")) return 0.01;
		
		if (verbose)
			System.out.println("Methods not recognized: " + m);
		
		return 0;
	}

	public double getPosteriorProbability(int fn) {
		return probabilisticDAG.getPosterior(fn);
	}

	public Vector<Integer> getLeafNames() {
		Vector<Integer> leaves = probabilisticDAG.getAllLeaves();
		Vector<Integer> rowNames = new Vector<Integer>();

		for (int i = 0; i < leaves.size(); i++) {
			rowNames.add(leaves.elementAt(i));
		}

		return rowNames;
	}

	public int getNumLeaves() {
		return probabilisticDAG.getNumLeaves();
	}

	////////////////////////////////////////////////////////////
	// Accessing DAG information for inference/learning       //
	////////////////////////////////////////////////////////////

	// This is an approximation to two leaf node subsets
	// It checks to see if evidence associated with the
	// leaf and each (for three, pair of) other leaves
	// have evidence associated with them. If so, the
	// evidence is combined; otherwise it is added to
	// the appropriate leaf (pair) probability.

	// This function maps the integer onto the vector
	// of leaves that it has stored.
	public double getLeafSubsetProbability(int leaf) {
		Vector<Integer> leaves = getAllLeaves();

		if (leaf >= 0 && leaf < leaves.size()) {
			Integer leafID = leaves.elementAt(leaf);
			return probabilisticDAG.getLeafSubsetProbability(leafID, leaves);
		}

		if (verbose)
			System.out.println("In GOOntologyWrapper: entered wrong id for leaf: " + leaf);
		
		return 0;
	}

	// Gets the prior probability of any leaf being present
	// in the function vector; this should be the same for
	// all leaf functions.
	public double getLeafSubsetPrior() {
		double lprior = 0.0;
		int[] leafv = probabilisticDAG.getLeafPolynomialVector(getNumLeaves());
		
		for (int i = 0; i < leafv.length; i++) {
			lprior += ((double)leafv[i] / (double)Math.pow(rValue, i + 1));
		}
		
		return lprior;
	}

	public double getSingleLeafPrior() {
		return 1 / rValue;
	}

	////////////////////////////////////////////////////////////
	// Local class definition for an entry into the GO DAG    //
	////////////////////////////////////////////////////////////

	private class Entry extends Object {
		public int goNum;
		public String goName;

		public Entry(int num, String name, int level) {
			goNum =	num;

			while (name.indexOf('\\') >= 0) {
				name = name.substring(0, name.indexOf('\\'))
				       + name.substring(name.indexOf('\\') + 1,
				                        name.length());
			}

			goName = name;
		}

		public String getName() {
			return goName;
		}

		public int getNumber() {
			return goNum;
		}

	}

	/**
	 * Method getAllLeaves.
	 * @return Vector
	 */
	public Vector<Integer> getAllLeaves() {
		if (leafList == null)
			leafList = probabilisticDAG.getAllLeaves();

		return leafList;
	}

	/**
	 * Method propagateLeavesUpwards.
	 */
	/*
	private void propagateLeavesUpwards()
	{
	    Enumeration pnames = functions.keys();
	    while(pnames.hasMoreElements()) {
	        Integer name = (Integer)pnames.nextElement();
	        functions.propagateLeavesUpwards(name);
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


	/** Once we get a new settings object, initialize any variables that we care about.
	 * For example, I know SifterPipelineObject cares about verbose, family name, etc.
	 */
	private void initSettings() {
		this.verbose = ((Boolean)this.settings.getSetting("verbose")).booleanValue();
	}

	public String toString() {
		String retval = "";

		retval += settings.toString();

		return retval;
	}

}

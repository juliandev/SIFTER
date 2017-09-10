/**
 * Copyright 2003-2005 Barbara Engelhardt (bee@cs.berkeley.edu)
 *
 * A single protein.
 * Used in each node of the reconciled tree.
 *
 * @author Barbara Engelhardt
 */

package sifter.components;

import java.io.PrintStream;
import java.util.Vector;

public class ProteinAnnotationObjectWithProbabilityMethods {

	private String id;
	private String name;
	private String function;
	private String origin;
	private String ecNum;
	private String pFam;
	private Vector<Integer> goNumbers;
	private Vector<String> mocs;
	private Vector<Integer> goReal;
	private String species;
	private String sequence;
	//private Hashtable fnProbabilities;
	private double[] fnLikelihoods;
	private Vector<Double> functionRatio;
	private Vector<Double> probabilities;

	public ProteinAnnotationObjectWithProbabilityMethods(String idN) {
		id = idN;
		name = null;
		function = null;
		origin = null;
		ecNum = null;
		pFam = null;
		goNumbers = null;
		mocs = null;
		goReal = null;
		species = null;
		sequence = null;
		//fnProbabilities = null;
		fnLikelihoods = null;
		functionRatio = null;
	}

	public String getID() {
		return id;
	}

	public String getName() {
		return name;
	}

	public void setName(String nameN) {
		if (nameN != null || name == null) {
			name = nameN;
		}
	}

	public String getFunction() {
		return function;
	}

	public void addFunction(String moreFunction) {
		if (function == null) {
			function = new String(moreFunction);
		}
		else {
			function += (" " + moreFunction);
		}
	}

	public String getOrigin() {
		return origin;
	}

	public void setOrigin(String o) {
		origin = o;
	}

	public String getSequence() {
		return sequence;
	}

	public void addSequence(String moreSequence) {
		if (sequence == null) {
			sequence = new String(moreSequence);
		}
		else {
			sequence += (moreSequence);
		}
	}

	public void addLikelihoods(double[] likelihoods) {
		fnLikelihoods = likelihoods;
	}

	public void addRatioLikelihoods() {
		if (functionRatio == null)
			return;

		fnLikelihoods = new double[functionRatio.size()];
		double total = 0;

		for (int i = 0; i < fnLikelihoods.length; i++) {
			fnLikelihoods[i] =
			  functionRatio.elementAt(i).doubleValue();
			total += fnLikelihoods[i];
		}

		for (int i = 0; i < fnLikelihoods.length; i++) {
			fnLikelihoods[i] = fnLikelihoods[i] / total;
		}

	}

	public double[] getLeafLikelihoods() {
		return fnLikelihoods;
	}

	public String getECNumber() {
		return ecNum;
	}

	public void setECNumber(String ecNumN) {
		if (ecNumN != null || ecNum == null) {
			ecNum = ecNumN;
		}
	}

	public String getSpecies() {
		return species;
	}

	public void setSpecies(String speciesN) {
		if (speciesN != null || species == null) {
			species = speciesN;
		}
	}

	public void setGOReal(Vector<Integer> gor) {
		if (goReal == null) {
			goReal = new Vector<Integer>();
		}

		goReal.addAll(gor);
	}

	public Vector<Integer> getGOReal() {
		return goReal;
	}

	public void setFunctionRatio(Vector<Double> fnRatio) {
		if (functionRatio == null)
			functionRatio = new Vector<Double>();

		functionRatio.addAll(fnRatio);
	}

	public Vector<Double> getFunctionRatio() {
		return functionRatio;
	}

	public String getPFam() {
		return pFam;
	}

	public void setPFam(String pFamN) {
		if (pFamN != null || pFam == null) {
			pFam = pFamN;
		}
	}

	public Vector<Integer> getGONumbers() {
		return goNumbers;
	}

	public void setGONumber(Integer goN) {
		if (goNumbers == null) {
			goNumbers = new Vector<Integer>();
		}

		if (goN != null) {
			goNumbers.add(goN);
		}
	}

	public void setGONumbers(Vector<Integer> goNs) {
		//if(name.indexOf("SCHPO") == -1) // Added for FUNGI experiments
		goNumbers = goNs;
	}

	// MOC = method of categorization
	public Vector<String> getMocs() {
		return mocs;
	}

	public void setMOC(String mocN) {
		if (mocs == null) {
			mocs = new Vector<String>();
		}

		if (mocN != null) {
			mocs.add(mocN);
		}
	}

	public void setMOCs(Vector<String> mocNs) {
		//if(name.indexOf("SCHPO") == -1) // Added for FUNGI experiments
		mocs = mocNs;
	}

	public void setProbabilities(Vector<Double> probabilities) {
		this.probabilities = probabilities;
	}

	public Vector<Double> getProbabilities() {
		return probabilities;
	}

	public boolean hasGORealAnnotation(Integer function) {
		if (getGOReal() == null)
			return false;

		Vector<Integer> gra = getGOReal();

		if (gra.contains(function))
			return true;

		return false;
	}

	public int numGORealAnnotation(int function) {
		if (getGOReal() == null)
			return 0;

		return (getGOReal().size());
	}

	public boolean hasAnnotation(Integer function, String mocLevel) {
		if (getGONumbers() == null)
			return false;

		Vector<Integer> gon = getGONumbers();
		Vector<String> moc = getMocs();

		if (gon.contains(function) && moc.contains(mocLevel)) {
			for (int i = 0; i < gon.size(); i++) {
				if (gon.elementAt(i) == function
				    && moc.elementAt(i).equals(mocLevel)) {
					return true;
				}
			}
		}

		return false;
	}

	public static boolean isNumeric(String str)  
	{  
	  try  
	  {  
	    @SuppressWarnings("unused")
		double d = Double.parseDouble(str);  
	  }  
	  catch(NumberFormatException nfe)  
	  {  
	    return false;  
	  }  
	  return true;  
	}
	
	public boolean hasProperEvidence(SifterExecutionSettingsContainer settings) {
		Vector<String> mocs = getMocs();
		for (String moc : mocs) {

			if (isNumeric(moc)){
				boolean isProper = ((Boolean)settings.getSetting("floats")).booleanValue();
				if (isProper)  { 
					
					
					return true;
				}
			}else{
				String moc_lowercase = moc.toLowerCase();
				boolean isProper = ((Boolean)settings.getSetting(moc_lowercase)).booleanValue();
	
				if (isProper) {
					return true;
				}
			}
		}

		return false;
	}


	// bee-- add in list of possible functions
	public boolean existsAnnotation(Integer function, String mocLevel) {
		if (getGONumbers() == null)
			return false;

		Vector<String> moc = getMocs();

		if (moc.contains(mocLevel)) {
			return true;
		}

		return false;
	}

	public void printProteinToFile(PrintStream fout) {
		printProteinToFile(fout, null);
	}

	public void printProteinToFile(PrintStream fout, String alignment) {

		fout.println("\t<Protein>");

		if (getName() != null)
			fout.println("\t\t<ProteinName>" + getName() + "</ProteinName>");

		if (getID() != null)
			fout.println("\t\t<ProteinNumber>" + getID() + "</ProteinNumber>");

		if (getOrigin() != null)
			fout.println("\t\t<ProteinLocation>" + getOrigin() + "</ProteinLocation>");

		if (getSpecies() != null)
			fout.println("\t\t<SpeciesName>" + getSpecies() + "</SpeciesName>");

		if (getECNumber() != null)
			fout.println("\t\t<ECNumber>" + getECNumber() + "</ECNumber>");

		if (getPFam() != null)
			fout.println("\t\t<PFamNumber>" + getPFam() + "</PFamNumber>");

		if (getGONumbers() != null)
			fout.println("\t\t<GONumber>" + getGONumbers() + "</GONumber>");

		if (getMocs() != null)
			fout.println("\t\t<MOC>" + getMocs() + "</MOC>");

		if (getSequence() != null)
			fout.println("\t\t<Sequence>" + getSequence() + "</Sequence>");

		if (alignment != null)
			fout.println("\t\t<Alignment>" + alignment + "</Alignment>");

		if (getFunction() != null)
			fout.println("\t\t<Function>" + getFunction() + "</Function>");

		fout.println("\t</Protein>");

	}

	public String toString() {
		String retval = "";

		retval += name + ", ";
		retval += "goNumbers: " + goNumbers.toString();

		return retval;
	}

}

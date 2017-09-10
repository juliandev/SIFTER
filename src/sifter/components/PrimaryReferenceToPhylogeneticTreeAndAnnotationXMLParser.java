/**
 * Reads in XML files:
 * Protein family with go annotations ".xml":
 * proteins/final/proteinfamily*.pli"
 * Returns set of proteins.
 * TODO: should be made to work with phyloXML files.
 *
 * Copyright 2003-2005 Barbara Engelhardt (bee@cs.berkeley.edu)
 * @author Barbara Engelhardt
 */

package sifter.components;

import java.util.*;
import java.io.*;
import java.util.regex.*;
//import org.apache.regexp.RE;


public class PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser {
	private Hashtable<String, String> proteinAlignments;
	private Hashtable<String, String> proteinNames; // key: pfam name value: swiss prot name
	private Hashtable<String, String> proteinIDs; // key: SPID: swiss prot name
	private Vector<String> proteinList; // key: swissprot names
	private String pFamNum;
	private String pFamName;
	private String pFamFunction;
	private String GTree;
	private String STree;
	private ProbabilisticReconciledPhylogenyObject phyTree;
	private int maxAlignment;

	// For parsing from XML file
	private String proteinName;
	private ProteinAnnotationObjectWithProbabilityMethods currentProtein;
	private Vector<ProteinAnnotationObjectWithProbabilityMethods> tempProteins;
	public static Pattern re = null;
	public static Pattern reStartLoner = null;
	public static Pattern reEndLoner = null;


	public PrimaryReferenceToPhylogeneticTreeAndAnnotationXMLParser(String id) {
		pFamNum = id.trim();
		pFamName = null;
		pFamFunction = null;
		GTree = null;
		STree = null;
		proteinNames = new Hashtable<String, String>();
		proteinIDs = new Hashtable<String, String>();
		//proteinIDs = null;
		proteinAlignments = new Hashtable<String, String>();
		proteinList = null;
		maxAlignment = 0;

		proteinName = null;
		currentProtein = null;
		tempProteins = null;
	}

	public void setName(String name) {
		if (name != null) {
			pFamName = name.trim();
		}
	}

	public String getName() {
		return pFamName;
	}

	public String getID() {
		return pFamNum.trim();
	}

	public void setProteinName(String spName, String pfName) {
		if (!proteinNames.containsKey(pfName)) {
			proteinNames.put(pfName, spName);
			proteinIDs.put(spName, pfName);
		}
		else {
			proteinNames.remove(pfName);
			proteinIDs.remove(spName);
			proteinNames.put(pfName, spName);
			proteinIDs.put(spName, pfName);
		}
	}

	public String getProteinName(String pfName) {
		if (proteinNames.containsKey(pfName)) {
			return proteinNames.get(pfName);
		}
		else
			return null;
	}

	public String getProteinNameFromID(String pID) {
		//Enumeration pnames = proteinIDs.keys();
		//while(pnames.hasMoreElements()) {
		//  String pfName = (String)pnames.nextElement();
		//System.out.println("Key: "+pfName+", element: "+proteinIDs.get(pfName));
		//}
		if (proteinIDs.containsKey(pID)) {
			return proteinIDs.get(pID).toLowerCase();
		}
		else
			return null;
	}

	// Returns a vector of the Names (i.e., swissprot IDs) of
	// the proteins associated with this family.
	public Vector<String> getProteins() {
		if (proteinList != null)
			return proteinList;

		proteinList = new Vector<String>();
		Enumeration<String> pnames = proteinNames.keys();

		while (pnames.hasMoreElements()) {
			String pfName = pnames.nextElement();
			String spName = proteinNames.get(pfName);

			if (!proteinList.contains(spName))
				proteinList.add(spName);
		}

		return proteinList;
	}

	// Same as above, but returns a Vector of pfam names
	public Vector<String> getProteinNames() {
		Vector<String> proteinNList = new Vector<String>();
		Enumeration<String> pnames = proteinNames.keys();

		while (pnames.hasMoreElements()) {
			String pfName = pnames.nextElement();

			if (!proteinNList.contains(pfName))
				proteinNList.add(pfName);
		}

		return proteinNList;
	}

	public void addPFamFunction(String pfFunction) {
		if (pFamFunction == null) {
			pFamFunction = pfFunction;
		}
		else {
			pFamFunction += pfFunction;
		}
	}

	public String getPFamFunction() {
		return pFamFunction;
	}

	// pfName is the name
	public void setProteinAlignment(String pfName, String alignment) {
		if (proteinNames.containsKey(pfName)) {
			String spName = proteinNames.get(pfName);

			if (!proteinAlignments.contains(spName)) {
				proteinAlignments.put(spName, alignment);
			}
			else {
				proteinAlignments.remove(spName);
				proteinAlignments.put(spName, alignment);
			}
		}
		else {
			System.out.println("Problem setting an alignment: protein " + pfName + " does not have a mapping");
			System.exit(1);
		}
	}

	public String getProteinAlignment(String spName) {
		if (proteinAlignments.contains(spName)) {
			return proteinAlignments.get(spName);
		}
		else
			return null;
	}

	public String getSpeciesTree() {
		return STree;
	}

	public String getGeneTree() {
		return GTree;
	}

	public void setMaxAlignment(int len) {
		if (maxAlignment < len)
			maxAlignment = len;
	}

	public int getMaxAlignment() {
		return maxAlignment;
	}

	// REFACTORME: unused
//    private int tooFewProteins(Hashtable proteins)
//    {
//	int pcount = 0;
//	Enumeration<String> pnames = getProteins().elements();
//	while(pnames.hasMoreElements()) {
//	    if(proteins.containsKey(pnames.nextElement())) {
//		pcount++;
//	    }
//	}
//	return pcount;
//    }


	public void addTree(ProbabilisticReconciledPhylogenyObject t) {
		phyTree = t;
	}

	public ProbabilisticReconciledPhylogenyObject getTree() {
		return phyTree;
	}


	////////////////////////////////////////////////////////////////
	// Input from XML file                                        //
	////////////////////////////////////////////////////////////////

	/** Imports XML file information into a Vector data structure.
	 * not stored locally.
	 * NYI: describe data structure.
	 * @return Returns a vector of PFunProteins that were found in the family.
	 */
	public Vector<ProteinAnnotationObjectWithProbabilityMethods> readInFromXMLFile(String fileName) {
		BufferedReader fin;
		tempProteins = new Vector<ProteinAnnotationObjectWithProbabilityMethods>();
		proteinList = new Vector<String>();
		setupRegex();
		
		try {
			File fr = new File(fileName);

			if (!fr.exists())
				return null;

			fin = new BufferedReader(new FileReader(fileName));
			String str = "Empty";
			int lnum = 1;

			try {
				while ((str = fin.readLine()) != null) {
					lnum++;
					processRegexResults(str);   // I'm confused... where does tempProteins come in, if at all?
				}
			}
			catch (Exception e) {
				System.out.println("Error in readInFromXMLFile1: " + str + " " + lnum + " " + e.getMessage());
			}

			fin.close();
		}
		catch (Exception ioe) {
			System.err.println("Error in readInFromXMLFile2: " + fileName + " " +
			                   ioe.getMessage());
			System.exit(1);
		}

		return tempProteins;
	}

	public static void setupRegex() {
		String notdelim = "([^<]+)";
		String notdelimxml = "([^<^>^/]+)";
		String regex1 = "<" + notdelimxml + ">" + notdelim + "</" + notdelimxml + ">";
		String regex2 = "<" + notdelimxml + ">";
		String regex3 = "</" + notdelimxml + ">";

		re = Pattern.compile(regex1);
		reStartLoner = Pattern.compile(regex2);
		reEndLoner = Pattern.compile(regex3);
	}

	// The following looks like it just parses XML, and that's it...
	// this could be better served by an XML-parsing library.
	// Just calls processContents(), processLonerEndTag(),
	// and processLonerStartTag().
	//
	public void processRegexResults(String s) {
		// processContents() stores stuff in instance variable currentProtein.
		// processendLonerTag() calls addProtein().
		String paren0, paren1, paren2, paren3;

		if (s.indexOf('>') >= 0) {
			paren0 = s.substring(s.indexOf('<') + 1, s.indexOf('>'));

			if (paren0.equals("SpeciesTree"))
				return;

			if (paren0.equals("SpeciesName"))
				return;

			if (paren0.equals("PFamNumber"))
				return;

			if (paren0.equals("Sequence"))
				return;

			if (s.lastIndexOf('<') > s.indexOf('<') &&
			    s.lastIndexOf('>') > s.indexOf('>')) {
				processContents(paren0, s.substring(s.lastIndexOf('<') + 2, s.lastIndexOf('>')), s.substring(s.indexOf('>') + 1, s.lastIndexOf('<')));
				return;
			}
			else
				if (paren0.indexOf('/') >= 0) {
					processLonerEndTag(paren0.substring(1, paren0.length()));
					return;
				}
				else {
					processLonerStartTag(paren0);
					return;
				}
		}

		Matcher matcher = reStartLoner.matcher(s);

		if (matcher.matches()) {
			paren0 = matcher.group(0);
			paren1 = matcher.group(1);

			if (paren1.equals("SpeciesTree")) {
				processContents("SpeciesTree", "SpeciesTree",
				                s.substring(s.indexOf('>') + 1,
				                            s.lastIndexOf('<')));
				return;
			}
		}

		matcher = re.matcher(s);
		Matcher matcherEndLoner = reEndLoner.matcher(s);
		Matcher matcherStartLoner = reStartLoner.matcher(s);

		if (matcher.matches()) {
			paren0 = matcher.group(0);
			paren1 = matcher.group(1);
			paren2 = matcher.group(2);
			paren3 = matcher.group(3);

			if (paren3 != null) {
				processContents(paren1, paren3, paren2);
			}
		}
		else
			if (matcherEndLoner.matches()) {
				paren0 = matcherEndLoner.group(0);
				paren1 = matcherEndLoner.group(1);
				processLonerEndTag(paren1);
			}
			else
				if (matcherStartLoner.matches()) {
					paren0 = matcherStartLoner.group(0);
					paren1 = matcherStartLoner.group(1);
					processLonerStartTag(paren1);
				}
	}

	// ONLY called from processRegexResults()
	public void processContents(String startTag, String endTag,
	                            String content) {

		if (!startTag.equals(endTag)) {
			System.out.println("Error in processContents: Start and End tags do not match: " + startTag + ", " + endTag);
			return;
		}

		if (startTag.equals("FamilyID")) {
			pFamNum = content;
		}
		else if (startTag.equals("FamilyName")) {
			setName(content);
		} else if (startTag.equals("ProteinName")) {
			proteinName = content;
		} else if (startTag.equals("ProteinNumber")) {
			currentProtein = new ProteinAnnotationObjectWithProbabilityMethods(content);

			if (proteinName != null) {
				currentProtein.setName(proteinName);
				proteinNames.put(proteinName, content);
				proteinIDs.put(content, proteinName);
				proteinList.add(content);
			}
		} else if (startTag.equals("ProteinLocation")) {
			currentProtein.setOrigin(content);
		} else if (startTag.equals("PFamNumber")) {
			currentProtein.setPFam(content);
		} else if (startTag.equals("SpeciesName")) {
			currentProtein.setSpecies(content);
		} else if (startTag.equals("ECNumber")) {
			currentProtein.setECNumber(content);
		} else if (startTag.equals("GOReal")) {
			currentProtein.setGOReal(parseVectorIntegers(content));
		} else if (startTag.equals("FunctionRatio")) {
			currentProtein.setFunctionRatio(parseVectorDoubles(content));
		} else if (startTag.equals("GOAssumed")) {
			currentProtein.setGOReal(parseVectorIntegers(content));
		} else if (startTag.equals("GOSP")) {
			currentProtein.setGOReal(parseVectorIntegers(content));
		} else if (startTag.equals("GOT")) {
			currentProtein.setGOReal(parseVectorIntegers(content));
		} else if (startTag.equals("GONumber")) {
			currentProtein.setGONumbers(parseVectorIntegers(content));
		} else if (startTag.equals("MOC")) {
			currentProtein.setMOCs(parseVector(content));
		//System.out.println("Got: "+currentProtein.getMOC());
		} else if (startTag.equals("Probability")) {
			currentProtein.setProbabilities(parseVectorDoubles(content));
		//System.out.println("Got: "+currentProtein.getMOC());
		} else if (startTag.equals("Sequence")) {
			//currentProtein.addSequence(content);
		} else if (startTag.equals("Alignment")) {
			setMaxAlignment(content.length());
			//if(proteinName != null)
			//	setProteinAlignment(proteinName, content);
			//else
			//    setProteinAlignment(currentProtein.getID(), content);
			//proteinName = null;
			//} else if (startTag.equals("SpeciesTree")) {
			//  STree = content;
		} else if (startTag.equals("GeneTree")) {
				GTree = content;
		} else if (startTag.equals("Structure")) {
				//GTree = content;
		} else if (startTag.equals("PubMed")) {
				//GTree = content;
		} else if (startTag.equals("Information")) {
				//GTree = content;
		} else if (startTag.equals("GOAssumed")) {
				//GTree = content;
		}
		else if (startTag.equals("GOSP")) {
				//GTree = content;
		}
		else if (startTag.equals("GOT")) {
				//GTree = content;
		}
		else if (startTag.equals("GOReal")) {
				//GTree = content;
		}
		else if (startTag.equals("SwissProtHits")) {
				// dont want to do this here.
		} else if (startTag.equals("GeneName")) {
			// dont want to do this here.
		} else if (startTag.equals("FunctionDescription")) {
				// dont want to do this here.

				//} else if (startTag.equals("Language")) {
				//  lang = content;
				//} else if (startTag.equals("Language")) {
				//  lang = content;
		} else {
			System.out.println("Unknown tag: " + startTag);
		}
	}

	// ONLY called from processRegexResults()
	public static void processLonerStartTag(String tag) {

	}

	// ONLY called from processRegexResults()
	public void processLonerEndTag(String tag) {
		if (tag.equals("Protein") && currentProtein != null) {
			addProtein(currentProtein);
			currentProtein = null;
			proteinName = null;
		}
	}

	// Called only from ~ processRegExp()
	public static Vector<String> parseVector(String s) {
		if (s.indexOf('[') < 0)
			return null;

		Vector<String> v = new Vector<String>();
		s = s.substring(s.indexOf('[') + 1, s.length() - 1);
		StringTokenizer st = new StringTokenizer(s, " ,");

		while (st.hasMoreTokens()) {
			v.add(st.nextToken());
		}

		return v;
	}

	// Called only from ~ processRegExp()
	public static Vector<Integer> parseVectorIntegers(String s) {
		if (s.indexOf('[') < 0)
			return null;

		Vector<Integer> v = new Vector<Integer>();
		s = s.substring(s.indexOf('[') + 1, s.length() - 1);
		StringTokenizer st = new StringTokenizer(s, " ,");

		while (st.hasMoreTokens()) {
			v.add(new Integer(st.nextToken()));
		}

		return v;
	}

	// Called only from ~ processRegExp()
	public static Vector<Double> parseVectorDoubles(String s) {
		if (s.indexOf('[') < 0)
			return null;

		Vector<Double> v = new Vector<Double>();
		s = s.substring(s.indexOf('[') + 1, s.length() - 1);
		StringTokenizer st = new StringTokenizer(s, " ,");

		while (st.hasMoreTokens()) {
			v.add(new Double(st.nextToken()));
		}

		return v;
	}


	// ONLY CALLED BY LonerEndTag()
	// and adds it to tempProtein().
	public void addProtein(ProteinAnnotationObjectWithProbabilityMethods p) {
		tempProteins.add(p);
	}

	// REFACTORME: unused
//    public void printOutFamily(String fileName, Hashtable proteins)
//    {
//	PrintStream fout;
//	fileName = fileName+"_"+pFamNum+".pfi";
//	int foundProteins = tooFewProteins(proteins);
//	if(foundProteins < 2) return;
//	try {
//	    fout =  new PrintStream(new FileOutputStream(new File(fileName)));
//
//	    fout.println("<Family>");
//	    fout.println("\t<FamilyID> "+pFamNum+" </FamilyID>\n");
//	    if(pFamName != null)
//		fout.println("\t<FamilyName> "+pFamName+" </FamilyName>\n");
//
//	    Enumeration<String> pnames = proteinAlignments.keys();
//	    while(pnames.hasMoreElements()) {
//		String spName = pnames.nextElement();
//		ProteinAnnotationObjectWithProbabilityMethods pfp = (ProteinAnnotationObjectWithProbabilityMethods)proteins.get(spName);
//		if(pfp != null)
//		    pfp.printProteinToFile(fout, proteinAlignments.get(spName));
//	    }
//	    if(GTree != null) {
//		fout.println("\t<GeneTree>"+GTree+"</GeneTree>");
//		fout.println("\t<GTMethod>"+GTMethod+"<GTMethod>");
//	    }
//	    if(STree != null) {
//		fout.println("\t<SpeciesTree>"+STree+"</SpeciesTree>");
//		fout.println("\t<STSource>"+STSource+"</STSource>");
//	    }
//	    if(pFamFunction != null)
//		fout.println("\t<FamilyFunction>"+pFamFunction+"</FamilyFunction>");
//	    fout.println("\t<SwissProtHits>"+foundProteins+", "+proteinList.size()+"</SwissProtHits>");
//	    fout.println("</Family>");
//
//	    fout.close();
//        }
//        catch (Exception ioe) {
//            System.err.println("PrintOutProteinFamily: " + fileName + " " +
//			       ioe.getMessage());
//            System.exit(1);
//        }
//    }

}

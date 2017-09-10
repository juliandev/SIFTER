/**
 * SIFTER project.
 * This file handles the command-line interface *and nothing else*.
 * 
 * Copyright (c) 2010, Barbara Engelhardt (bee@compbio.berkeley.edu)
 *
 * @author Barbara Engelhardt (primary investigator and author)
 * @author Steven R. Chan (later code hacking, documentation, GUI)
 * @version 1.0
 */

package sifter.main;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.cli.UnrecognizedOptionException;

import sifter.components.SifterExecutionSettingsContainer;

// Insert constants here?


public class Sifter {
	public static void main(String[] args) throws ParseException, IllegalArgumentException, FileNotFoundException, ClassNotFoundException {
		// Collect and parse settings, failing if necessary.
		SifterExecutionSettingsContainer settings = parseCommandLineOptions(args);
		
		new SifterPipelineObject(settings);
	}
	
	/**
	 * A repository of settings.
	 */
	static SifterExecutionSettingsContainer settings;
	
	public static final boolean FILE_EXISTENCE_REQUIRED = true;
	public static final boolean FILE_EXISTENCE_OPTIONAL = false;
	
	private static SifterExecutionSettingsContainer parseCommandLineOptions(String[] args) {
		settings = new SifterExecutionSettingsContainer();
		
		CommandLineParser parser = new PosixParser();
		Options options = buildOptions();
		CommandLine line = null; 
		
		try {
			line = parser.parse(options, args);
			args = line.getArgs();
		} catch (UnrecognizedOptionException e) {
			System.err.println("ERROR: Unrecognized options. ("+ e.getMessage() + ")");
			printHelp();
		} catch (MissingArgumentException e) {
			System.err.println("ERROR: Incomplete arguments. Some arguments require values.");
			printHelp();
		} catch (ParseException e) {
			e.printStackTrace();
		}
		
		// Make sure we have exactly one argument for FAMILY setting.
		String family = "";
		
		if (args.length != 1) {
			System.err.println("ERROR: Please specify run name: ");
			
			for(int i = 0; i < args.length; i++)
				System.out.println(args[i]);
			
			printHelp();
			System.exit(1);
		} else {
			family = args[0];
			settings.setSetting("family", family);
			//if ((Boolean)settings.getSetting("verbose"))
				System.out.println("Family set to " + family);
		}
		
		/**
		 * FILENAME INFORMATION
		 * [option, setting, description, help]
		 */
		// These are all completely arbitrary and shouldn't be hard-coded this way.
		
		String[] FILE_OUTPUT_DEFAULTS	= {"output/default.rdata"};
		String[] FILE_OUTPUT			= {"output", "output", "Resulting output file", "File in which Sifter outputs its final results."};
		String[] FILE_PLI_DEFAULTS		= {"proteinfamily_" + family + ".pli", "proteinfamily_" + family.toUpperCase() + ".pli", family + ".pli", family.toUpperCase() + ".pli"};
		String[] FILE_PLI				= {"pli", "proteinFilename", "Protein family file (.pli)", "To obtain a .pli file, use 'pfam2pli' script. See README.txt."};
		String[] FILE_PHYLOXML_DEFAULTS = {"reconciled_" + family +".xml", "reconciled_" + family +".xml", "reconciled_" + family.toUpperCase()+".xml", "reconciled_" + family.toUpperCase() +".xml", family +".xml", family.toUpperCase() +".xml"};
		String[] FILE_PHYLOXML			= {"phylo", "reconciledFilename (.xml)", "Reconciled tree in PhyloXML", "To obtain a reconciled tree file, use 'pli2tree' script. See README.txt."};
		
		String[] FILE_FX_DEFAULTS		= {"infer-" + family + ".fx"};
		String[] FILE_FX				= {"fx", "familyFilename", "Parameter file (.fx)", "To obtain a .fx file, run 'sifter --generate'. See README.txt."};
		String[] FILE_SCALE_DEFAULTS	= {"scale-" + family + ".fx"};
		String[] FILE_SCALE				= {"sfx", "scaleParamsFilename", "Scale parameter file (.fx)", "To obtain a .fx scale file, run with --generate. See README.txt."};
		String[] FILE_ALPHA_DEFAULTS	= {"alpha-" + family + ".fx"};
		String[] FILE_ALPHA				= {"afx", "alphaParamsFilename", "Alpha parameter file (.fx)", "To obtain a .fx scale file, run with --generate. See README.txt."};
		
		String[] FILE_ONTOLOGY_DEFAULTS = {"function.ontology"};
		String[] FILE_ONTOLOGY			= {"ontology", "ontology", "Functions gene ontology database", "Get function.ontology from www.geneontology.org."};
		
		// We insist on having these be specified/exist:
		try {
			prepareFilename(FILE_ONTOLOGY, FILE_ONTOLOGY_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED);
			prepareFilename(FILE_OUTPUT, FILE_OUTPUT_DEFAULTS, family, line);
		} catch (FileNotFoundException e) {
				System.err.println("ERROR: File not found.\n   " 
						+ e.getMessage() + ".\n" 
						+ "   To specify files, please use the options "
						+ "listed under Sifter's help by entering:\n" 
						+ "    java -jar sifter.jar --help");
		}
		
		// Now parsing command lines:
		try {
			// Generic options
			if (line.hasOption("help")) {
				printHelp();
				System.exit(0);
			}
	
			// Verbosity
			if (line.hasOption("v")) {
				settings.setSetting("verbose", new Boolean(true));
				System.out.println("**************************************************");
				System.out.println("JOB SUMMARY");
			} else {
				settings.setSetting("verbose", new Boolean(false));
			}
			
			// Determine what run-mode the program should use.
			// "g" is "generate"
			if (line.hasOption("g")) {
				settings.setSetting("runmode", new String("generate"));
				
				if (line.hasOption("folds"))
					settings.setSetting("folds", Integer.valueOf(line.getOptionValue("folds")));
				
				if (line.hasOption("truncation"))
					settings.setSetting("truncation", Integer.valueOf(line.getOptionValue("truncation")));
				
				
				settings.setSetting("pli",   prepareFilename(FILE_PLI, FILE_PLI_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("fx",    prepareFilename(FILE_FX, FILE_FX_DEFAULTS, family, line));
				settings.setSetting("scale", prepareFilename(FILE_SCALE, FILE_SCALE_DEFAULTS, family, line));
				settings.setSetting("alpha", prepareFilename(FILE_ALPHA, FILE_ALPHA_DEFAULTS, family, line));
			}
			// "X-Validate" run mode
			else if (line.hasOption("xval")) {
				settings.setSetting("runmode", new String("xvalidate"));
				 
				if (line.hasOption("iter"))
					settings.setSetting("iterations",Integer.getInteger(line.getOptionValue("iter")));
				
				if (line.hasOption("step"))
					settings.setSetting("stepsize", Double.valueOf(line.getOptionValue("step")));
				
				if (line.hasOption("cutoff"))
					settings.setSetting("cutoff", Double.valueOf(line.getOptionValue("cutoff")));
					
				if (line.hasOption("folds"))
					settings.setSetting("folds", Integer.valueOf(line.getOptionValue("folds")));
				
				if (line.hasOption("truncation"))
					settings.setSetting("truncation", Integer.valueOf(line.getOptionValue("truncation")));
				
				if (line.hasOption("em"))
					settings.setSetting("em", new Boolean(true));
				
				
				settings.setSetting("pli",   prepareFilename(FILE_PLI, FILE_PLI_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("phylo", prepareFilename(FILE_PHYLOXML, FILE_PHYLOXML_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("fx",    prepareFilename(FILE_FX, FILE_FX_DEFAULTS, family, line));
				settings.setSetting("scale", prepareFilename(FILE_SCALE, FILE_SCALE_DEFAULTS, family, line));
				settings.setSetting("alpha", prepareFilename(FILE_ALPHA, FILE_ALPHA_DEFAULTS, family, line));
			}
			// Else if run mode is expectation-maximization
			else if (line.hasOption("em")) {
				settings.setSetting("runmode", new String("expectationmaximization"));
				
				if (line.hasOption("iter"))
					settings.setSetting("iterations", Integer.valueOf(line.getOptionValue("iter")));
				
				if (line.hasOption("step"))
					settings.setSetting("stepsize", Double.valueOf(line.getOptionValue("step")));
				
				if (line.hasOption("folds"))
					settings.setSetting("folds", Integer.valueOf(line.getOptionValue("folds")));
				
				if (line.hasOption("truncation"))
					settings.setSetting("truncation", Integer.valueOf(line.getOptionValue("truncation")));
				
				if (line.hasOption("cutoff"))
					settings.setSetting("cutoff", Double.valueOf(line.getOptionValue("cutoff")));
				
				settings.setSetting("pli",   prepareFilename(FILE_PLI, FILE_PLI_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("phylo", prepareFilename(FILE_PHYLOXML, FILE_PHYLOXML_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("fx",    prepareFilename(FILE_FX, FILE_FX_DEFAULTS, family, line));
				settings.setSetting("scale", prepareFilename(FILE_SCALE, FILE_SCALE_DEFAULTS, family, line));
				settings.setSetting("alpha", prepareFilename(FILE_ALPHA, FILE_ALPHA_DEFAULTS, family, line));
			}
			// Else if run mode is GAPE
			else if (line.hasOption("gape")) {
				settings.setSetting("runmode", new String("gaparameterestimation"));
				
				if (line.hasOption("population"))
					settings.setSetting("population", Integer.valueOf(line.getOptionValue("population")));
				
				if (line.hasOption("generations"))
					settings.setSetting("generations", Integer.valueOf(line.getOptionValue("generations")));
				
				if (line.hasOption("geneticoperators"))
					settings.setSetting("geneticoperators", Double.valueOf(line.getOptionValue("geneticoperators")));
				
				settings.setSetting("pli",   prepareFilename(FILE_PLI, FILE_PLI_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("phylo", prepareFilename(FILE_PHYLOXML, FILE_PHYLOXML_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("fx",    prepareFilename(FILE_FX, FILE_FX_DEFAULTS, family, line));
				settings.setSetting("scale", prepareFilename(FILE_SCALE, FILE_SCALE_DEFAULTS, family, line));
				settings.setSetting("alpha", prepareFilename(FILE_ALPHA, FILE_ALPHA_DEFAULTS, family, line));
			}
			// Otherwise go to the default run mode.
			else {
				settings.setSetting("runmode", new String("default"));
				
				if (line.hasOption("truncation"))
					settings.setSetting("truncation", Integer.valueOf(line.getOptionValue("truncation")));
				
				
				if (line.hasOption("folds")) 
					settings.setSetting("folds", Integer.valueOf(line.getOptionValue("folds")));
				
				settings.setSetting("pli",   prepareFilename(FILE_PLI, FILE_PLI_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("reconciledFilename", prepareFilename(FILE_PHYLOXML, FILE_PHYLOXML_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("phylo", prepareFilename(FILE_PHYLOXML, FILE_PHYLOXML_DEFAULTS, family, line, FILE_EXISTENCE_REQUIRED));
				settings.setSetting("fx",    prepareFilename(FILE_FX, FILE_FX_DEFAULTS, family, line));
				settings.setSetting("scale", prepareFilename(FILE_SCALE, FILE_SCALE_DEFAULTS, family, line));
				settings.setSetting("alpha", prepareFilename(FILE_ALPHA, FILE_ALPHA_DEFAULTS, family, line));
			}
			
			// The following is a completely retarded way to specify which
			// evidence should be included.
			
			settings.setSetting("bg", new Boolean(line.hasOption("bg")));
	
			//settings.setSetting("imp", new Boolean(true));
			//settings.setSetting("ida", new Boolean(true));  
	
			settings.setSetting("exp", new Boolean(line.hasOption("exp")));
			settings.setSetting("ida", new Boolean(line.hasOption("ida")));
			settings.setSetting("ipi", new Boolean(line.hasOption("ipi")));
			settings.setSetting("imp", new Boolean(line.hasOption("imp")));
			settings.setSetting("igi", new Boolean(line.hasOption("igi")));
			settings.setSetting("iep", new Boolean(line.hasOption("iep")));
			
			settings.setSetting("tas", new Boolean(line.hasOption("tas")));
			settings.setSetting("nas", new Boolean(line.hasOption("nas")));
			
			settings.setSetting("iss", new Boolean(line.hasOption("iss")));
			settings.setSetting("iso", new Boolean(line.hasOption("iso")));
			settings.setSetting("isa", new Boolean(line.hasOption("isa")));
			settings.setSetting("ism", new Boolean(line.hasOption("ism")));
			settings.setSetting("igc", new Boolean(line.hasOption("igc")));
			settings.setSetting("iba", new Boolean(line.hasOption("iba")));
			settings.setSetting("ibd", new Boolean(line.hasOption("ibd")));
			settings.setSetting("ikr", new Boolean(line.hasOption("ikr")));
			settings.setSetting("ird", new Boolean(line.hasOption("ird")));
			settings.setSetting("rca", new Boolean(line.hasOption("rca")));

			settings.setSetting("ic", new Boolean(line.hasOption("ic")));
			settings.setSetting("nd", new Boolean(line.hasOption("nd")));
			
			settings.setSetting("iea", new Boolean(line.hasOption("iea")));
			settings.setSetting("nr", new Boolean(line.hasOption("nr")));
			settings.setSetting("floats", new Boolean(line.hasOption("floats")));
			// Rest is jeff:
			settings.setSetting("gene3d", new Boolean(line.hasOption("gene3d")));
			settings.setSetting("hamap", new Boolean(line.hasOption("hamap")));
			settings.setSetting("panther", new Boolean(line.hasOption("panther")));
			settings.setSetting("pfam", new Boolean(line.hasOption("pfam")));
			settings.setSetting("pirsf", new Boolean(line.hasOption("pirsf")));
			settings.setSetting("prints", new Boolean(line.hasOption("prodom")));
			settings.setSetting("prodom", new Boolean(line.hasOption("prodom")));
			settings.setSetting("prosite_patterns", new Boolean(line.hasOption("prosite_patterns")));
			settings.setSetting("prosite_profiles", new Boolean(line.hasOption("prosite_profiles")));
			settings.setSetting("smart", new Boolean(line.hasOption("smart")));
			settings.setSetting("superfamily", new Boolean(line.hasOption("superfamily")));
			settings.setSetting("tigrfams", new Boolean(line.hasOption("tigrfams")));
			
			settings.setSetting("genemania", new Boolean(line.hasOption("genemania")));
			settings.setSetting("pagosub", new Boolean(line.hasOption("pagosub")));
			settings.setSetting("protfun", new Boolean(line.hasOption("protfun")));
			
 			/*
			if (line.hasOption("goe")) {
				settings.setSetting("exp", true);
				settings.setSetting("ida", true);
				settings.setSetting("ipi", true);
				settings.setSetting("imp", true);
				settings.setSetting("igi", true);
				settings.setSetting("iep", true);
			}
			
			if (line.hasOption("goc")) {
				settings.setSetting("iss", true);
				settings.setSetting("iso", true);
				settings.setSetting("isa", true);
				settings.setSetting("ism", true);
				settings.setSetting("igc", true);
				settings.setSetting("rca", true);
			}
			
			if (line.hasOption("goa")) {
				settings.setSetting("tas", true);
				settings.setSetting("nas", true);
			}
			
			if (line.hasOption("gocu")) {
				settings.setSetting("ic", true);
				settings.setSetting("nd", true);
			}
			
			if (line.hasOption("goal")) {
				settings.setSetting("imp", true);
				settings.setSetting("ida", true);
				settings.setSetting("exp", true);
				settings.setSetting("ida", true);
				settings.setSetting("ipi", true);
				settings.setSetting("imp", true);
				settings.setSetting("igi", true);
				settings.setSetting("iep", true);
				settings.setSetting("iss", true);
				settings.setSetting("iso", true);
				settings.setSetting("isa", true);
				settings.setSetting("ism", true);
				settings.setSetting("igc", true);
				settings.setSetting("rca", true);
				settings.setSetting("tas", true);
				settings.setSetting("nas", true);
				settings.setSetting("ic", true);
				settings.setSetting("nd", true);
				settings.setSetting("iea", true);
				settings.setSetting("nr", true);
			}
			
			if (line.hasOption("ip")) {
				settings.setSetting("gene3d", true);
				settings.setSetting("hamap", true);
				settings.setSetting("panther", true);
				settings.setSetting("pfam",true);
				settings.setSetting("pirsf", true);
				settings.setSetting("prints", true);
				settings.setSetting("prodom", true);
				settings.setSetting("prosite_patterns", true);
				settings.setSetting("prosite_profiles", true);
				settings.setSetting("smart", true);
				settings.setSetting("superfamily", true);
				settings.setSetting("tigrfams", true);
			}
			
			if (line.hasOption("sink")) {
				settings.setSetting("imp", true);
				settings.setSetting("ida", true);
				settings.setSetting("exp", true);
				settings.setSetting("ida", true);
				settings.setSetting("ipi", true);
				settings.setSetting("imp", true);
				settings.setSetting("igi", true);
				settings.setSetting("iep", true);
				settings.setSetting("iss", true);
				settings.setSetting("iso", true);
				settings.setSetting("isa", true);
				settings.setSetting("ism", true);
				settings.setSetting("igc", true);
				settings.setSetting("rca", true);
				settings.setSetting("tas", true);
				settings.setSetting("nas", true);
				settings.setSetting("ic", true);
				settings.setSetting("nd", true);
				settings.setSetting("iea", true);
				settings.setSetting("nr", true);
				
				settings.setSetting("gene3d", true);
				settings.setSetting("hamap", true);
				settings.setSetting("panther", true);
				settings.setSetting("pfam",true);
				settings.setSetting("pirsf", true);
				settings.setSetting("prints", true);
				settings.setSetting("prodom", true);
				settings.setSetting("prosite_patterns", true);
				settings.setSetting("prosite_profiles", true);
				settings.setSetting("smart", true);
				settings.setSetting("superfamily", true);
				settings.setSetting("tigrfams", true);
				
				settings.setSetting("genemania", true);
				settings.setSetting("pagosub", true);
				settings.setSetting("protfun", true);
			}
			*/
		} catch (FileNotFoundException e) {
			System.err.println("FILE ERROR: \n   " 
					+ e.getMessage() + ".\n" 
					+ "   To specify files, please use the options "
					+ "listed under Sifter's help by entering:\n" 
					+ "    java -jar sifter.jar --help");
		}
		
		return settings;
	}

	private static String prepareFilename(String[] fileinfo, String[] defaultFilenames, String family, CommandLine line)
	throws FileNotFoundException {
		return prepareFilename(fileinfo, defaultFilenames, family, line, FILE_EXISTENCE_OPTIONAL);
	}

	/** Prepares filenames for use by the rest of the program 
	 * (like initGenerate).
	 * Checks for its existence as well.
	 * 
	 * @param fileinfo Fileinfo arrays, such as FILE_FX coded in main().
	 * @param defaultFilenames Default file name arrays, such as that coded in main().
	 * @param family A protein family name, such as 'pf20005'
	 * @param line Command-line options
	 * @param existence_required Does the file need to exist?
	 * @return String final filename to use
	 * @throws FileNotFoundException If user specifies an invalid file, 
	 *         or the default file does not exist, throw an exception.
	 */
	private static String prepareFilename(String[] fileinfo, String[] defaultFilenames, String family, CommandLine line, boolean existence_required) throws FileNotFoundException {
		//for (int i = 0; i < defaultFilenames.length; i++) System.out.println("File default: '"+defaultFilenames[i]+"'");
		String optionName = fileinfo[0];
		String settingName = fileinfo[1];
		String desc = fileinfo[2];
		String help = fileinfo[3];    
		String userFilename = null;
		String finalFilename = null;
		
		if (line.hasOption(optionName)) { 
			userFilename = line.getOptionValue(optionName);
		}
		
		Filefinder f = new Filefinder();
		
		try {
			finalFilename = f.seek(desc, help, userFilename, defaultFilenames);
		} catch (FileNotFoundException e) {
			if (existence_required) {
				throw e;
			} else {
				// What if the user specifies an invalid file
				//and we don't need it to exist?
				if (userFilename != null) {
					finalFilename = userFilename;
				} else {
					finalFilename = defaultFilenames[0];
					
				}
			}
		}
		
		settings.setSetting(settingName, finalFilename);
		if ((Boolean)settings.getSetting("verbose")) {
			System.out.println(desc + " ("+settingName+"): " + finalFilename);
		}
		return finalFilename;
	}

	/**
	 * Prints command-line help.
	 * @see sifter.main
	 */
	private static void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("java -jar sifter.jar [OPTIONS] FAMILYNAME", buildOptions());
	}

	/** Builds the set of possible options for this program.
	 * @return An Options file that can be used by the CommandLine parser.
	 */
	@SuppressWarnings("static-access")
	private static Options buildOptions() {
		Options res = new Options();
		
		res.addOption("v", "verbose", false, "Verbose operation.");
		res.addOption("g", "generate", false, "(Run mode) Generates a set of input parameters for the inference problem.");
		res.addOption("em", "em", false, "(Run mode) Perform EM to estimate parameters");
		res.addOption("xval", "xvalidation", false, "(Run mode) Use cross-validation with EM.");
		
		// GAPE
		// first param = runmode
		// second param = arg console 
		res.addOption("gape", "GAPE", false, "(Run mode) Use GAPE to estimate parameters.");
		res.addOption(OptionBuilder.withLongOpt("population").withDescription("Number of individuals for genetic algorithm").withArgName("number").hasArg().create("population"));
		res.addOption(OptionBuilder.withLongOpt("generations").withDescription("Number of generations for genetic algorithm").withArgName("number").hasArg().create("generations"));
		res.addOption(OptionBuilder.withLongOpt("geneticoperators").withDescription("Rate of genetic operators for genetic algorithm").withArgName("number").hasArg().create("geneticoperators"));
		// --- 
		
		res.addOption(OptionBuilder.withLongOpt("output").withDescription("Set output file (default: output_directory/default.rdata)").withArgName("filename").hasArg().create("output"));
		res.addOption(OptionBuilder.withLongOpt("protein").withDescription("Set protein file (default: proteins/proteinfamily_<FAMILY>.pli)").withArgName("filename").hasArg().create("pli"));
		res.addOption(OptionBuilder.withLongOpt("reconciled").withDescription("Set reconciled .xml tree (default: reconciled/reconciled_<FAMILY>.xml)").withArgName("filename").hasArg().create("phylo"));
		res.addOption(OptionBuilder.withLongOpt("familyfile").withDescription("Set family .fx parameter filename (default: data/infer-<FAMILY>.fx)").withArgName("filename").hasArg().create("fx"));
		res.addOption(OptionBuilder.withLongOpt("scale").withDescription("Set family .fx scale filename (default: data/scale-<FAMILY>.fx)").withArgName("filename").hasArg().create("sfx"));
		res.addOption(OptionBuilder.withLongOpt("alpha").withDescription("Set family .fx alpha filename (default: data/alpha-<FAMILY>.fx)").withArgName("filename").hasArg().create("afx"));
		res.addOption(OptionBuilder.withLongOpt("ontology").withDescription("Specify which ontology file you want (default: \"data/function.ontology\")").withArgName("filename").hasArg().create("ontology"));
		res.addOption(OptionBuilder.withLongOpt("help").withDescription("Show help for arguments. (More help is available via README.txt)").create());
		res.addOption(OptionBuilder.withLongOpt("iter").withDescription("Number of iterations. At the moment, this applies only to EM. ").withArgName("number").hasArg().create("iter"));
		res.addOption(OptionBuilder.withLongOpt("step").withDescription("Step size for gradient ascent in EM (M-step)").withArgName("number").hasArg().create("step"));
		res.addOption(OptionBuilder.withLongOpt("cutoff").withDescription("Cutoff delta for gradient ascent in EM (M-step)").withArgName("number").hasArg().create("cutoff"));
		res.addOption(OptionBuilder.withLongOpt("folds").withDescription("Number of folds in cross validation, leave-one-out is 0").withArgName("number").hasArg().create("folds"));
		res.addOption(OptionBuilder.withLongOpt("truncation").withDescription("Number of functions to truncate to in approximation").withArgName("number").hasArg().create("truncation"));
		
		
		res.addOption("exp", "with-exp", false, "(Experiment) Use GOA protein annotations inferred from experiment.");
		res.addOption("ida", "with-ida", false, "(Experiment) Use GOA protein annotations inferred from direct assay.");
		res.addOption("ipi", "with-ipi", false, "(Experiment) Use GOA protein annotations from those inferred from physical interaction.");
		res.addOption("imp", "with-imp", false, "(Experiment) Use GOA protein annotations from those inferred from mutant phenotype.");
		res.addOption("igi", "with-igi", false, "(Experiment) Use GOA protein annotations from those inferred from genetic interaction.");
		res.addOption("iep", "with-iep", false, "(Experiment) Use GOA protein annotations from those inferred from expression profiles.");
		
		res.addOption("tas", "with-tas", false, "(Author Statement) Use GOA protein annotations from traceable author statements.");
		res.addOption("nas", "with-nas", false, "(Author Statement) Use GOA protein annotations from non-traceable author statements.");
		
		res.addOption("iss", "with-iss", false, "(Computational) Use GOA protein annotations from those inferred from sequence similarity.");
		res.addOption("iso", "with-iso", false, "(Computational) Use GOA protein annotations from those inferred from sequence orthology.");
		res.addOption("isa", "with-isa", false, "(Computational) Use GOA protein annotations from those inferred from sequence alignment.");
		res.addOption("ism", "with-ism", false, "(Computational) Use GOA protein annotations from those inferred from sequence model.");
		res.addOption("igc", "with-igc", false, "(Computational) Use GOA protein annotations from those inferred from genomic context.");
		res.addOption("iba", "with-iba", false, "(Computational) Use GOA protein annotations from those inferred from biological aspect of ancestor.");
		res.addOption("ibd", "with-ibd", false, "(Computational) Use GOA protein annotations from those inferred from biological aspect of descendant.");
		res.addOption("ikr", "with-ikr", false, "(Computational) Use GOA protein annotations from those inferred from key residues.");
		res.addOption("ird", "with-ird", false, "(Computational) Use GOA protein annotations from those inferred from rapid divergence.");
		res.addOption("rca", "with-rca", false, "(Computational) Use GOA protein annotations from those reconstructed from computational analyses.");
		
		res.addOption("ic", "with-ic", false, "(Curator Statement) Use GOA protein annotations from those inferred from curator");
		res.addOption("nd", "with-nd", false, "(Curator Statement) Use GOA protein annotations from those inferred from curator with no biological data available");
		
		res.addOption("iea", "with-iea", false, "(Automatically Assigned) Use GOA protein annotations inferred by electronic annotation.");
		res.addOption("nr", "with-nr", false, "(Obsolete) Use protein annotations with source not recorded.");
		
		res.addOption("floats", "with-floats", false, "Use floating numbers as confidence values.");
		
		
		
		// REFACTORME: interpro annotations
		res.addOption("gene3d", "with-genethreed", false, "Use gene3d predictions in InterPro.");
		res.addOption("hamap", "with-hamap", false, "Use hamap predictions in InterPro.");
		res.addOption("panther", "with-panther", false, "Use panther predictions in InterPro.");
		res.addOption("pfam", "with-pfam", false, "Use Pfam predictions in InterPro.");
		res.addOption("pirsf", "with-pirsf", false, "Use PIRSF predictions in InterPro.");
		res.addOption("prints", "with-prints", false, "Use prints predictions in InterPro.");
		res.addOption("prodom", "with-prodom", false, "Use ProDom predictions in InterPro.");
		res.addOption("prosite_patterns", "with-prosite-patterns", false, "Use PROSITE pattern predictions in InterPro.");
		res.addOption("prosite_patterns", "with-prosite-profiles", false, "Use PROSITE profile predictionss in InterPro.");
		res.addOption("smart", "with-smart", false, "Use SMART predicitons in InterPro.");
		res.addOption("superfamily", "with-superfamily", false, "Use SUPERFAMILY predictions in InterPro.");
		res.addOption("tigrfams", "with-tigr", false, "Use TIGRFAMs predictions in InterPro.");
		
		// REFACTORME: other "atomic" annotations (enables one MOC)
		res.addOption("genemania", "with-genemania", false, "Use protein predictions from genemania.");
		res.addOption("pagosub", "with-pagosub", false, "Use protein predictions from pagosub.");
		res.addOption("protfun", "with-protfun", false, "Use protein predictions from ProtFun.");
		
		/*
		// REFACTORME: psuedo annotations (enables more than one MOC without 1000 CL arguments)
		res.addOption("goe", "with-go-exp", false, "Use experimental go annotations.");
		res.addOption("goc", "with-go-comp", false, "Use computational go annotations.");
		res.addOption("goa", "with-go-auth", false, "Use author go annotations.");
		res.addOption("gocu", "with-go-curator", false, "Use curator go annotations.");
		res.addOption("goal", "with-go-all", false, "Use all go annotations.");
		res.addOption("ip", "with-interpro", false, "Use protein predictions from interpro.");
		res.addOption("sink", "with-all", false, "Use all evidence.");
		*/
		
		return res;
	}

	
	//  filefinder
	//  for the Sifter project
	//
	//  Created by Steven Chan < www.stevenchan.us > on 10/16/05.
	//
	// Usage:
	//  Filefinder finder = new Filefinder(".");
	//  String[] default_values = {"molecular.ont", "molecular.ontology"};  // not case-sensitive
	//  finder.seek("Molecular function ontology file",
	//              "Downloadable from http://geneontology.org",
	//              "user-molecular.ont",
	//              default_values);
	private static class Filefinder {
		private HashMap<String, File> directoryMap;

		/**
		* Constructor. Prepares internal variables, such as the list of files in the specified directory.
		 * @param startDirectory The directory to search for files.
		 **/
		public Filefinder(String startDirectory) {
			this.directoryMap = this.hashDirectory(startDirectory);
			//this.startDirectory = startDirectory;
			// System.out.println( this.directoryMap.toString() ); // TODO - REMOVE LATER
		}
		/**
		 * Constructor. If no parameters, assume that the starting directory is the current working directory.
		 **/
		public Filefinder() {
			this(".");
		}

		/**
		 * Test code. This is NOT necessary since Filefinder's purpose is not to be run.
		 * If executed with "java Filefinder", however, it will run a few tests.
		 * Feel free to delete this method if it causes problems.
		 **/
		/*public static void main(String args[]) {
			// MAIN-RUNNING CODE
			//   This code is not necessary since it's not the main program
			//   being run, but useful for running tests.

			Filefinder f = new Filefinder("/Users/steven/sifter");

			try {
				String[] defaults = {"mol-function.ont"};
				System.out.println(f.seek("molecular function ontology file", "Downloadable from ONTOLOGY", null, defaults));
			}
			catch (FileNotFoundException e) {
				System.out.println("File not found: " + e.getMessage());
			}

			try {
				String[] defaults = {"mol-function.ont"};
				System.out.println(f.seek("molecular function ontology file", "Downloadable from ONTOLOGY", "nonexistent-file", defaults));
			}
			catch (FileNotFoundException e) {
				System.out.println("File not found: " + e.getMessage());
			}

			try {
				String[] defaults = {"function.ontology"};
				System.out.println(f.seek("molecular function ontology file", "Downloadable from ONTOLOGY", null, defaults));
			}
			catch (FileNotFoundException e) {
				System.out.println("File not found: " + e.getMessage());
			}

		}*/


		public String seek(String description,
		                   String help,
		                   String uservalue,
		                   String[] default_filenames) throws FileNotFoundException {
			String res;

			// If the user specifies a filename, then do something with it.
			if (uservalue != null) {
				if ((new File(uservalue).exists())) {
					return uservalue;
				}
				else {
					// Throw a new exception saying that the following files could work.
					try {
						res = fileIn(default_filenames, this.directoryMap);
					}
					catch (FileNotFoundException e) {
						// If the default files don't work, then give them the help message.
						throw new FileNotFoundException(uservalue + " does not exist. " + help);
					}

					throw new FileNotFoundException(uservalue + " does not exist. This file may work, however: " + res + ".");

				}
			}
			else {
				// If the user does not specify anything, then try to find files from the default filenames array.
				try {
					res = fileIn(default_filenames, this.directoryMap);
				}
				catch (FileNotFoundException e) {
					// If the default files don't work, then give them the help message.
					throw new FileNotFoundException(description + " does not exist. " + help);
				}
				
				return res;
			}
		}

		/**
		 * Creates a HashMap representing all the files in the specified directory.
		 * Searches recursively through the directory and adds key "filename", actual file object, into
		 * the resulting HashMap.
		 * I'm using File objects b/c that's the most effective and direct way to determine
		 * whether a given File object represents an actual file or a directory,
		 * since we're recursing.
		 *
		 * @param startDirectory The directory to search for files.
		 * @return HashMap object with String keys and File objects for values
		 **/
		private HashMap<String, File> hashDirectory(String startDirectory) {
			HashMap<String, File> res = new HashMap<String, File>();

			File[] list = new File(startDirectory).listFiles();

			for (int i = 0; i < list.length; i++) {
				File item = list[i];
				// System.out.println(item);    // Print the item.

				// If it's a directory,
				//    recursively call hashDirectory() on it. Save it in hash.
				//    Add all those hashes to res.
				if (item.isDirectory()) {
					HashMap<String, File> recursedMap = hashDirectory(item.getPath());
					res.putAll(recursedMap);
				}
				else {
					// Add the lower-cased filename to hash, plus the file object.
					String lowerCaseName = item.getName().toLowerCase();
					res.put(lowerCaseName, item);
				}
			}

			return res;
		}

		/**
		 * Given an array of possible filenames (not pathnames!) and a List of all the files,
		 * tries to see if there are matches in the List.
		 * (Does not use case sensitivity.)
		 *
		 * @param possibilities
		 * @param fileMap
		 * @return String containing a perfectly valid filename
		 * @throws FileNotFoundException
		 **/
		private String fileIn(String[] possibilities, HashMap<String, File> fileMap) throws FileNotFoundException {
			// At this point, we've exhausted our list of possibilities,
			// so complain by throwing a filenotfoundexception.
			// TODO: pass possibilities array to the new Exception. Can we do this in Java?
			throw new FileNotFoundException();
		}

	}
}

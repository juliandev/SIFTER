# SIFTER
### Parameter estimation for phylogeny-based protein function prediction using SIFTER (Statistical Inference of Function Through Evolutionary Relationships)

SIFTER is a statistical approach to predicting protein function that uses a protein family's phylogenetic tree, as the natural structure for representing protein relationships, overlaid with all known protein functions in the family.

This package provides a genetic algorithm for parameter estimation of SIFTER algorithm. 

The parameter estimation of SIFTER through Genetic Algorithms called GAPE (Genetic Algorithm for Parameter Estimation) is developed by Julian Castañeda and Tania Rodriguez at Department of Systems Engineering and Computing, Antonio Narino University.

The original SIFTER algorithm is developed by Barbara E Engelhardt. 
Original paper:
- Engelhardt BE, Jordan MI, Srouji JR, Brenner SE. 2011. 
Genome-scale phylogenetic function annotation of large and diverse protein families. Genome Research 21:1969-1980. doi:10.1101/gr.104687.109 

Please cite the last paper:
- Sahraeian SME, Luo KR, Brenner SE (2015)
SIFTER-T: A scalable and optimized framework for the SIFTER phylogenomic method of probabilistic protein domain annotation. BioTechniques 58:140-142. 
doi: 10.2144/000114266

Other previous developers: Philip Johnson, Steven R. Chan, Micheal Souza.

You can also use the SIFTER webserver at [Berkeley](http://sifter.berkeley.edu) to access online the predictions on 16,863,537 proteins across 232,403 species.

## Usage
usage: `java -jar sifter.jar [OPTIONS] FAMILYNAME`

    -sfx,--scale <filename>             Set family .fx scale filename
                                        (default: data/scale-<FAMILY>.fx)
    -afx,--alpha <filename>             Set family .fx alpha filename
                                        (default: data/alpha-<FAMILY>.fx)
    -fx,--familyfile <filename>         Set family .fx parameter
                                        filename (default: data/infer-<FAMILY>.fx)
    -exp,--with-exp                     (Experiment) Use GOA protein
                                        annotations inferred from experiment.
    -floats,--with-floats               Use floating numbers as confidence values.
    -genemania,--with-genemania         Use protein predictions from genemania.
    -gene3d,--with-genethreed           Use gene3d predictions in InterPro.
    -hamap,--with-hamap                 Use hamap predictions in InterPro.
    -iba,--with-iba                     (Computational) Use GOA
                                        protein annotations from those inferred
                                        from biological aspect of ancestor.
                                        
    -ibd,--with-ibd                     (Computational) Use GOA
                                        protein annotations from those inferred
                                        from biological aspect of
                                        descendant.
    -ic,--with-ic                       (Curator Statement) Use GOA
                                        protein annotations from those
                                        inferred from curator
    -ida,--with-ida                     (Experiment) Use GOA protein
                                        annotations inferred from direct assay.
    -iea,--with-iea                     (Automatically Assigned) Use
                                        GOA protein annotations inferred by
                                        electronic annotation.
    -iep,--with-iep                     (Experiment) Use GOA protein
                                        annotations from those inferred from
                                        expression profiles.
    -igc,--with-igc                     (Computational) Use GOA protein
                                        annotations from those inferred from
                                        genomic context.
    -igi,--with-igi                     (Experiment) Use GOA protein annotations
                                        from those inferred from genetic
                                        interaction.
    -ikr,--with-ikr                     (Computational) Use GOA protein
                                        annotations from those inferred 
                                        from key residues.
    -imp,--with-imp                     (Experiment) Use GOA protein annotations
                                        from those inferred from mutant phenotype.
    -ipi,--with-ipi                     (Experiment) Use GOA protein
                                        annotations from those inferred from
                                        physical interaction.
    -ird,--with-ird                     (Computational) Use GOA protein
                                        annotations from those inferred from
                                        rapid divergence.
    -isa,--with-isa                     (Computational) Use GOA protein 
                                        annotations from those inferred
                                        from sequence alignment.
    -ism,--with-ism                     (Computational) Use GOA protein
                                        annotations from those inferred from
                                        sequence model.
    -iso,--with-iso                     (Computational) Use GOA
                                        protein annotations from those inferred
                                        from sequence orthology.
    -iss,--with-iss                     (Computational) Use GOA protein
                                        annotations from those inferred from
                                        sequence similarity.
    -nd,--with-nd                       (Curator Statement) Use GOA
                                        protein annotations from those inferred
                                        from curator with no biological data
                                        available
    -nr,--with-nr                       (Obsolete) Use protein annotations with
                                        source not recorded.
    -folds,--folds <number>             Number of folds in cross validation,
                                        leave-one-out is 0
    -pagosub,--with-pagosub             Use protein predictions from pagosub.
    -panther,--with-panther             Use panther predictions in InterPro.
    -pfam,--with-pfam                   Use Pfam predictions in InterPro.
    -pirsf,--with-pirsf                 Use PIRSF predictions in InterPro.
    -prints,--with-prints               Use prints predictions in InterPro.
    -prodom,--with-prodom               Use ProDom predictions in InterPro.
    -prosite_patterns,                  Use PROSITE profile predictions in
    --with-prosite-profiles             InterPro.
    -protfun,--with-protfun             Use protein predictions from ProtFun.
    -rca,--with-rca                     (Computational) Use GOA
                                        protein annotations from those
                                        reconstructed from computational analyses.
    -smart,--with-smart                 Use SMART predicitons in InterPro.
    -step,--step <number>               Step size for gradient ascent
                                        in EM (M-step)
    -cutoff,--cutoff <number>           Cutoff delta for gradient ascent in EM
                                        (M-step)
    -superfamily,--with-superfamily     Use SUPERFAMILY predictions in InterPro.
    -tas,--with-tas                     (Author Statement) Use GOA
                                        protein annotations from traceable 
                                        author statements.
    -tigrfams,--with-tigr               Use TIGRFAMs predictions in InterPro.
    -nas,--with-nas                     (Author Statement) Use GOA
                                        protein annotations from non-traceable
                                        author statements.
    -xval,--xvalidation                 (Run mode) Use cross-validation with EM.
    --help                              Show help for arguments. (More help is
                                        available via README.txt)
    -em,--em                            (Run mode) Perform EM to
                                        estimate parameters
    -g,--generate                       (Run mode) Generates a set of
                                        input parameters for the inference problem.
    -iter,--iter <number>               Number of iterations. At the
                                        moment, this applies only to EM.
    -ontology,--ontology <filename>     Specify which ontology file
                                        you want (default:"data/function.ontology")
    -output,--output <filename>         Set output file (default:
                                        output_directory/default.rdata)
    -phylo,--reconciled <filename>      Set reconciled .xml tree (default:
                                        reconciled/reconciled_<FAMILY>.xml)
    -pli,--protein <filename>           Set protein file (default:
                                        proteins/proteinfamily_<FAMILY>.pli)
    -truncation,--truncation <number>   Number of functions to
                                        truncate to in approximation
    -gape,--GAPE                        (Run mode) Use GAPE to estimate 
                                        parameters.
    -generations,--generations <number> Number of generations for genetic 
                                        algorithm.
    -geneticoperators,                  Rate of genetic operators
    --geneticoperators <number>         for genetic algorithm.
    -population,--population <number>   Number of individuals for
                                        genetic algorithm.
    -v,--verbose                        Verbose operation.
**Note about command line options:** use option name with `--`, except in the case of verbose (use `-v`).

## Setup
To configure the files, databases, ontologies, etc. you can use the scripts developed by Sahraeian et. al, see the **README** in the scripts/ directory.

After you have successfully:
1. generated a phylogeny (and put it in <SIFTER>/reconciled/reconciled-<FAMILY>.nex)
2. generated a .pli file (and put it in <SIFTER>/proteins/proteinfamily-<FAMILY>.pli)
3. download the appropriate ontology (from Gene Ontology, renamed if necessary and placed in <SIFTER>/data/function.ontology)

Made the java code: `java -jar sifter.jar <FAMILY> --generate` which will generate the parameter files with __GAPE__.

`java -jar sifter.jar <FAMILY> -v`

Then play around with learning parameters, different datasets/phylogenies, cross validation, etc.

We have included files for a family called "test" to run here.

## Debbuging
Please send any problems/comments/questions about **SIFTER** to bee@compbio.berkeley.edu.

Please send any problems/comments/questions about **GAPE** to julicastaneda@uan.edu.co.

## SIFTER output
The output for SIFTER (running inference) is a tab-delimited file (default: output/default.rdata) with the following columns:

`<NODE NAME> <POSTERIOR FN1> ... <POSTERIOR FNm> <MAX POSTERIOR PREDICTION>`

The order of the functions is identical to the order in the transition matrix parameter file (infer-<FAMILY>.fx): see initial row for specific order.

The output for SIFTER (running GAPE) is the scale and transition matrix parameter files. Where they are exactly is output at the end of GAPE.

The output for SIFTER (running leave-one-out cross-validation) is in the command line: the percentage of left-out elements that are correct according to their annotations. However, these are often wrong, so I'd would advise you to double check them manually from the output (I was lazy enough not to consider ties in posterior probabilities).
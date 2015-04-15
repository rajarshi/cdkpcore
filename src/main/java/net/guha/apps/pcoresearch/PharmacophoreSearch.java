package net.guha.apps.pcoresearch;

import org.apache.commons.cli.*;
import org.openscience.cdk.CDK;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.ConformerContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.iterator.IteratingMDLConformerReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.pharmacophore.*;
import org.xml.sax.SAXException;

import java.io.*;
import java.text.DecimalFormat;
import java.util.List;

/**
 * @author Rajarshi Guha
 */
public class PharmacophoreSearch {
    private boolean verbose = false;
    private boolean details = false;
    private boolean annotate = false;

    private String ofilename = null;
    private String ifilename = null;
    private String qfilename = null;
    private String qname = null;

    private BufferedWriter report = null;
    private MDLV2000Writer writer;
    private PharmacophoreMatcher matcher;
    private static final String PCORE_VERSION = "1.3.1";

    DecimalFormat formatter = new DecimalFormat("0.00");

    public String getQname() {
        return qname;
    }

    public void setQname(String qname) {
        this.qname = qname;
    }

    public String getOfilename() {
        return ofilename;
    }

    public void setOfilename(String ofilename) {
        this.ofilename = ofilename;
    }

    public String getIfilename() {
        return ifilename;
    }

    public void setIfilename(String ifilename) {
        this.ifilename = ifilename;
    }

    public String getQfilename() {
        return qfilename;
    }

    public void setQfilename(String qfilename) {
        this.qfilename = qfilename;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void setDetails(boolean details) {
        this.details = details;
    }

    public void setAnnotate(boolean annotate) {
        this.annotate = annotate;
    }

    public boolean getAnnotate() {
        return this.annotate;
    }

    private String getHitFileName(String qfile, String ifile) {
        String[] toks1 = qfile.split("\\.");
        String[] toks2 = ifile.split("\\.");
        return toks1[0] + "." + toks2[0] + ".sdf";
    }

    public void initialize() throws IOException, CDKException {
        List<PharmacophoreQuery> queries = PharmacophoreUtils.readPharmacophoreDefinitions(qfilename);
        if (queries.size() == 0) throw new CDKException("No queries found in " + getQfilename());

        PharmacophoreQuery query = null;
        if (qname == null) {
            query = queries.get(0);
            qname = (String) query.getProperty(CDKConstants.TITLE);
        } else {
            for (PharmacophoreQuery q : queries) {
                String title = (String) q.getProperty(CDKConstants.TITLE);
                if (title != null && title.equals(qname)) {
                    query = q;
                    break;
                }
            }
        }
        if (query == null && qname != null)
            throw new CDKException("Query named '" + qname + "' was not found in " + qfilename);

        matcher = new PharmacophoreMatcher(query);
        writer = new MDLV2000Writer(new FileWriter(ofilename));
        report = new BufferedWriter(new FileWriter("report.txt"));
        report.write("Serial\tTitle\tNconf\tNhit\n");
    }

    public void doSingleSearch() throws IOException, CDKException {
        IteratingSDFReader reader = new IteratingSDFReader(
                new FileReader(new File(ifilename)), DefaultChemObjectBuilder.getInstance()
        );

        int nmol = 0;
        int nhit = 0;
        int nskip = 0;

        long timeStart = System.currentTimeMillis();

        while (reader.hasNext()) {
            IAtomContainer container = reader.next();

            if (!GeometryTools.has3DCoordinates(container)) {
                nskip++;
                continue;
            }

            try {
                Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(),
                        Cycles.vertexShort());
                aromaticity.apply(container);
            } catch (CDKException e) {
                nskip++;
                continue;
//                throw new CDKException("Error in aromaticity detection");
            }

            boolean matched;
            try {
                matched = matcher.matches(container);
            } catch (CDKException e) {
                nskip++;
                continue;
            }

            if (matched) {
                nhit++;
                //TODO: This is a hack, since the matcher does not have the matching
                // patoms before this is called. So the USA matches gives a NPE. Need
                // to fix in the CDK
                matcher.getMatchingPharmacophoreAtoms();
                List<List<PharmacophoreAtom>> matches = matcher.getUniqueMatchingPharmacophoreAtoms();

                // loop over each of the matched
                for (List<PharmacophoreAtom> match : matches) {
                    for (PharmacophoreAtom patom : match) {
                        IAtom pseudoAtom = DefaultChemObjectBuilder.getInstance().newInstance(IAtom.class, "Xe");
                        pseudoAtom.setPoint3d(patom.getPoint3d());
                        container.addAtom(pseudoAtom);
                    }
                }
                try {
                    writer.writeMolecule(container);
                } catch (Exception e) {
                    throw new CDKException("ERROR: Problem writing a hit to disk [title = " + container.getProperty(CDKConstants.TITLE) + "]");
                }
            }

            report.write(nmol + "\t" + container.getProperty(CDKConstants.TITLE) + "\tNA\t" + matched + "\n");
            if (matched && details) {
                List<List<IBond>> matchingPBonds = matcher.getMatchingPharmacophoreBonds();
                int mcount = 0;
                for (List<IBond> bondMatch : matchingPBonds) {
                    report.write("MATCH " + (++mcount) + ": ");
                    for (IBond constraint : bondMatch) {
                        if (constraint instanceof PharmacophoreBond) {
                            PharmacophoreBond pbond = (PharmacophoreBond) constraint;
                            double dist = pbond.getBondLength();
                            IAtom group1 = pbond.getAtom(0);
                            IAtom group2 = pbond.getAtom(1);
                            report.write("(" + group1.getSymbol() + "," +
                                    group2.getSymbol() + "," + formatter.format(dist) + ") ");
                        } else if (constraint instanceof PharmacophoreAngleBond) {
                            PharmacophoreAngleBond pbond = (PharmacophoreAngleBond) constraint;
                            double angle = pbond.getBondLength();
                            IAtom group1 = pbond.getAtom(0);
                            IAtom group2 = pbond.getAtom(1);
                            IAtom group3 = pbond.getAtom(2);
                            report.write("(" + group1.getSymbol() + "," +
                                    group2.getSymbol() + "," + group3.getSymbol() +
                                    "," + formatter.format(angle) + ") ");

                        }

                    }
                    report.write("\n");
                }
                report.write("\n");
            }

            nmol++;
            if (verbose && nmol % 100 == 0)
                System.out.print("\rINFO: Processed " + nmol + " [hits = " + nhit + " skip = " + nskip + "]");
        }
        writer.close();
        report.close();
        long timeEnd = System.currentTimeMillis();
        double elapsed = ((timeEnd - timeStart) / 1000.0);
        double avg = elapsed / (double) (nmol + nskip);
        if (verbose) {
            System.out.println("\nINFO: Processed " + (nmol + nskip) + " molecules in " + formatter.format(elapsed) + "s " + "[" +
                    formatter.format(avg) + " s/mol] " +
                    "and got " + nhit + " hits");
        }
    }

    public void doConfSearch() throws IOException, CDKException {
        IteratingMDLConformerReader reader = new IteratingMDLConformerReader(
                new FileReader(new File(ifilename)), DefaultChemObjectBuilder.getInstance()
        );


        int nmol = 0;
        int nhit = 0;
        int nskip = 0;

        long timeStart = System.currentTimeMillis();

        ConformerContainer confContainer;
        while (reader.hasNext()) {
            confContainer = (ConformerContainer) reader.next();

            // we need to do aromaticity detection, since the
            // atom container we get back is just a ref, doing
            // the arom detection on it should set it for all
            // other conformers in the container
            IAtomContainer tmp = confContainer.get(0);
            try {
                Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(),
                        Cycles.vertexShort());
                aromaticity.apply(tmp);
            } catch (CDKException e) {
                nskip++;
                continue;
//                throw new CDKException("Error in aromaticity detection");
            }

            boolean firstTime = true;
            int confhits = 0;
            for (IAtomContainer conf : confContainer) {
                boolean status;
                if (firstTime) {
                    status = matcher.matches(conf, true);
                    firstTime = false;
                } else {
                    status = matcher.matches(conf, false);
                }
                if (status) try {
                    nhit++;
                    confhits++;

                    matcher.getMatchingPharmacophoreAtoms();

                    // loop over each of the matched and add dummy atoms (using Xe)
                    // Note that the IAtomContainer object obtained from confContainer
                    // should not be modified - hence we have to make a clone
                    IAtomContainer confClone = (IAtomContainer) conf.clone();
                    List<List<PharmacophoreAtom>> matches = matcher.getUniqueMatchingPharmacophoreAtoms();
                    for (List<PharmacophoreAtom> match : matches) {
                        for (PharmacophoreAtom patom : match) {
                            IAtom pseudoAtom = DefaultChemObjectBuilder.getInstance().newInstance(IAtom.class, "Xe");
                            pseudoAtom.setPoint3d(patom.getPoint3d());
                            confClone.addAtom(pseudoAtom);
                        }
                    }
                    writer.writeMolecule(confClone);
                } catch (Exception e) {
                    throw new CDKException("ERROR: problem writing a hit to disk [title = " + confContainer.getTitle() + "]");
                }
            }

            report.write(nmol + "\t" + confContainer.getTitle() + "\t" + confContainer.size() + "\t" + confhits + "\n");

            nmol++;
            if (verbose && nmol % 100 == 0)
                System.out.print("\rINFO: Processed " + nmol + " [hits = " + nhit + " skip = " + nskip + "]");
        }
        writer.close();
        report.close();
        long timeEnd = System.currentTimeMillis();
        double elapsed = ((timeEnd - timeStart) / 1000.0);
        double avg = elapsed / (double) (nmol + nskip);
        if (verbose) {
            System.out.println("\nINFO: Processed " + (nmol + nskip) + " molecules in " + formatter.format(elapsed) + "s " + "[" +
                    formatter.format(avg) + " s/mol] " +
                    "and got " + nhit + " hits");
        }
    }

    public boolean isVerbose() {
        return verbose;
    }

    public static void main(String[] args) throws IOException, CDKException {

        Options options = new Options();
        options.addOption("h", "help", false, "Print this message");
        options.addOption("v", "verbose", false, "Verbose output");
        options.addOption("a", "annotate", false, "Annotates the output SD file with pharmacophore groups." +
                " The result is that the each molecule in the hit file will have pseudo atoms" +
                " representing the pharmacophore groups that match the query. " +
                " The pseudo atoms are indicated by using the Xe symbol in the SD file." +
                " Currently only works for the non-conformer mode.");
        options.addOption("d", "details", false, "Detailed output. This implies verbose output and in" +
                " addition prints the details (such as exact distances) for each hit");
        options.addOption("V", "version", false, "Version");
        options.addOption("c", "conf", false, "Input file is conformer data. If this is the" +
                " case then the conformers for a given molecule should be contiguous and have the same" +
                " title. Currently conformer detection by isomorphism is not supported");
        options.addOption(OptionBuilder.withLongOpt("sdfile").withArgName("file").hasArg().withDescription("Input file. Can be a set of unique molecules or a " +
                "collection of conformers for a set of molecules")
                .create("sdfile"));
        options.addOption(OptionBuilder.withLongOpt("ofile").withArgName("file")
                .hasArg()
                .withDescription("Output file. Default is hits.sdf")
                .create("ofile"));
        options.addOption(OptionBuilder.withLongOpt("query").withArgName("file")
                .hasArg()
                .withDescription("The pharmacophore query in XML format")
                .create("query"));
        options.addOption(OptionBuilder.withLongOpt("qname").withArgName("name")
                .hasArg()
                .withDescription("If the query file contains multiple queries you can specify" +
                        " a query by its name field. If this argument" +
                        " is not provided the first query in the file is used")
                .create("qname"));
        options.addOption(OptionBuilder.withLongOpt("validate").withArgName("file")
                .hasArg()
                .withDescription("The query XML file to validate. Currently does not check for incorrect SMARTS definitions")
                .create("validate"));
        options.addOption(OptionBuilder.withLongOpt("align").withArgName("method")
                .hasArg()
                .withDescription("If specified, hits are aligned using one of the following algorithms: first, jvd")
                .create("align"));


        CommandLine line = null;
        try {
            CommandLineParser parser = new PosixParser();
            line = parser.parse(options, args);
        } catch (ParseException exception) {
            System.err.println("Unexpected exception: " + exception.toString());
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("\njava -jar CDKPsearch.jar",
                    "Pharmacophore searching based on the CDK. The code currently handles distances and angle constraints. Pharmacophore queries can be provided in the form of XML definition files\n\n" +
                            "Rajarshi Guha <rajarshi.guha@gmail.com>\n",
                    options, "");
            System.exit(-1);
        }

        PharmacophoreSearch ps = new PharmacophoreSearch();

        boolean annotate = false;
        boolean useConfs = false;

        if (line.hasOption("validate")) {
            String qfilename = line.getOptionValue("validate");
            try {
                ValidateQuery vq = new ValidateQuery();
                if (vq.validate(qfilename)) {
                    System.out.println("INFO: " + qfilename + " is a valid query file");
                }
            } catch (SAXException e) {
                System.out.println("ERROR: " + qfilename + " is an invalid query file");
                System.out.println("  " + e.getMessage());
            }
            System.exit(0);
        }

        if (line.hasOption("version") || line.hasOption("V")) {
            System.out.println("CDKPSearch version = " + PCORE_VERSION + " JRE " + System.getProperty("java.version") + " CDK " + CDK.getVersion());
            System.exit(-1);
        }

        if (line.hasOption("details") || line.hasOption("d")) {
            ps.setDetails(true);
            ps.setVerbose(true);
        }

        if (line.hasOption("verbose") || line.hasOption("v")) ps.setVerbose(true);
        if (line.hasOption("conf") || line.hasOption("c")) useConfs = true;
        if (line.hasOption("sdfile")) ps.setIfilename(line.getOptionValue("sdfile"));
        if (line.hasOption("ofile")) ps.setOfilename(line.getOptionValue("ofile"));
        if (line.hasOption("qname")) ps.setQname(line.getOptionValue("qname"));
        if (line.hasOption("query") || line.hasOption("q")) ps.setQfilename(line.getOptionValue("query"));
        if (line.hasOption("annotate") || line.hasOption("a")) ps.setAnnotate(true);

        if (ps.getIfilename() == null || ps.getQfilename() == null ||
                line.hasOption("h") || line.hasOption("help")) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PharmacophoreSearch", options);
            System.exit(-1);
        }

        // some simple checking
        File checker = new File(ps.getIfilename());
        if (!checker.exists()) {
            System.out.println(ps.getIfilename() + " does not exist!");
            System.exit(-1);
        }
        checker = new File(ps.getQfilename());
        if (!checker.exists()) {
            System.out.println(ps.getQfilename() + " does not exist!");
            System.exit(-1);
        }
        if (ps.getOfilename() == null) ps.setOfilename(ps.getHitFileName(ps.getQfilename(), ps.getIfilename()));
        ps.initialize();

        if (ps.isVerbose()) {
            System.out.println("INFO: Hits will go to " + ps.getOfilename());
            System.out.println("INFO: Using " + ps.getQname() + " from " + ps.getQfilename());
        }

        if (useConfs) {
            if (ps.isVerbose()) System.out.println("INFO: Will process as conformers");
            ps.doConfSearch();
        } else {
            if (ps.isVerbose()) System.out.println("INFO: Will not process conformers");
            ps.doSingleSearch();
        }
    }


}

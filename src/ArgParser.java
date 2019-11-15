package eu.ernstthuer;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

/**
 * Created by ethur on 7/26/16.
 * Argparser parses arguments
 */


class ArgParser {
    //private static ArgumentParser parser = ArgumentParsers.newArgumentParser("Checksum").defaultHelp(true).description("PARSparser");
    private static ArgumentParser parser = ArgumentParsers.newFor("Checksum").build().defaultHelp(true).description("PARSParser");

    //static List<FileHandler> fileList = new ArrayList<>();
    private double minqual;
    private int mincount;
    private int offset;
    private int qualStringency;
    private boolean samflag;
    private String outfileLocation;

    private eu.ernstthuer.BEDHandler bedHandler;
    private eu.ernstthuer.BamHandler bamHandler;

    eu.ernstthuer.BEDHandler getBedHandler() {
        return bedHandler;
    }

    eu.ernstthuer.BamHandler getBamHandler() {
        return bamHandler;
    }

    double getMinqual() {
        return minqual;
    }

    public int getMincount() {
        return mincount;
    }

    public boolean getSamflag() {
        return samflag;
    }

    public int getOffset() {
        return offset;
    }

    public String getOutfileLocation() {
        return outfileLocation;
    }

    ArgParser(String[] args) {


        // ToDo  put bam as mandatory
        parser.addArgument("-a", "--bam")
                .help("bam file (sorted)").required(false).dest("inBAM");
        parser.addArgument("-b", "--bed")
                .help("input file in BED file format").required(true).dest("inBED");
        parser.addArgument("-o", "--offset")
                .help("reads are counted for base upstream, default -1 (since version 0.67b) process before cutting site").required(false).setDefault(0).dest("offset");
        parser.addArgument("-out", "--outfile")
                .help("output in tsv file format containing Transcript identities and coverage per position").required(false).setDefault("Output_Glory_To_Ernst.tsv").dest("out");
        parser.addArgument("-q", "--minqual")
                .help("min mapping quality,  default 0").required(false).setDefault(0.0).dest("minqual");
        parser.addArgument("-m", "--mincount")
                .help("min TOTAL counts for given transcript, default 5").required(false).setDefault(5).dest("mincount");
        parser.addArgument("-i", "--ignore").help("analyze transcripts without name,  default true").required(false).setDefault(true).dest("ignore");
        parser.addArgument("-S", "--SamFlags")
                .help("Consider sam-flags for multi mapping exclusion,  default true checks SAM flags for multi mapping indication").required(false).dest("samflag");
        parser.addArgument("-Q", "--QualityStringency")
                .help("Set Stringency for output quality.  1 ignores CIGAR string problems. 2 ignores mapping issues and CIGAR, default 0 cuts >1 mutations in CIGAR and uses quality thresolds").required(false).setDefault(0).dest("stringency");

        Namespace ns = null;

        try {
            ns = parser.parseArgs(args);

            try {
                minqual = Double.parseDouble(ns.get("minqual"));
            }catch (ClassCastException e){
                minqual = ns.get("minqual");
            }
            try {
                mincount = Integer.parseInt(ns.get("mincount"));
            }catch (ClassCastException e){
                mincount = ns.get("mincount");
            }

            try {
                offset = Integer.parseInt(ns.get("offset"));
            }catch (ClassCastException e){
                offset = ns.get("offset");
            }

            try {
                samflag = Boolean.parseBoolean(ns.get("samflag"));
            }catch (ClassCastException e){
                samflag = ns.get("samflag");
            }

            try {
                qualStringency = Integer.parseInt(ns.get("stringency"));
            }catch (ClassCastException e){
                qualStringency = ns.get("stringency");
            }


            switch (qualStringency) {
                case 0:
                    System.out.println("[STATUS] Default,  running with normal stringency, cutting CIGAR based Stringency and quality threshold");
                    break;
                case 1:
                    System.out.println("[STATUS] Loosening Stringency, ignoring CIGAR strings");
                    break;
                case 2:
                    System.out.println("[STATUS] No quality control on BAM file, allows multimapped reads");
                    break;
                default:
                    System.out.println("[STATUS]  could not decide on Quality, running Default,  running with normal stringency, cutting CIGAR based Stringency and quality threshold");
                    qualStringency = 0;
                    break;
            }

            //System.out.println(offset + "  " + mincount + "  " + minqual);

            outfileLocation = ns.get("out");

        } catch (NumberFormatException num) {
            System.out.println("[ERROR] Please add numerical values to the parameters -o (offset) and -m (minCount) ");
            num.printStackTrace();
        } catch (ArgumentParserException e) {
            parser.handleError(e);
            System.exit(1);
        }

        try {
            bedHandler = new eu.ernstthuer.BEDHandler(ns.get("inBED").toString(), "BED", ns.get("ignore"));
            bamHandler = new eu.ernstthuer.BamHandler(ns.get("inBAM").toString(), "BAM", Boolean.parseBoolean(ns.get("samflag")), offset, qualStringency);
        } catch (NullPointerException e) {
            System.out.println(" File not found exception ");
        }


        //System.out.println(ns.get("inGFF").toString());

/*        // String locale, String type, String feature, String direction))

        FastaHandler inFasta = new FastaHandler(ns.get("inFasta").toString(), "FASTA", "Input");
        FastaHandler outFasta = new FastaHandler(ns.get("mOut").toString(), "FASTA", "Output");
        CSVHandler finalOut = new CSVHandler(ns.get("outFinal").toString(), "VCF", "Output");

        if (existingSNPinfo) {
            CSVHandler vcfInput = new CSVHandler(ns.get("VCFIN").toString(), "VCF", "INPUT");
            fileList.add(vcfInput);
        }


        fileList.add(gffreader);
        fileList.add(inFasta);
        fileList.add(outFasta);
        fileList.add(finalOut);*/
    }

}



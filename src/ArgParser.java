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
    private static ArgumentParser parser = ArgumentParsers.newArgumentParser("Checksum").defaultHelp(true).description("PARSparser");
    //static List<FileHandler> fileList = new ArrayList<>();
    private double minqual;
    private int mincount;
    private int offset;
    private int qualStringency;
    private String outfileLocation;

    private BEDHandler bedHandler;
    private BamHandler bamHandler;

    BEDHandler getBedHandler() {
        return bedHandler;
    }

    BamHandler getBamHandler() {
        return bamHandler;
    }

    double getMinqual() {
        return minqual;
    }

    public int getMincount() {
        return mincount;
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
                .help("reads are counted for base upstream, default -1 (since version 0.67b) process before cutting site").required(false).setDefault(-1).dest("offset");
        parser.addArgument("-out", "--outfile")
                .help("output in tsv file format containing Transcript identities and coverage per position").required(false).setDefault("Output_Glory_To_Ernst.tsv").dest("out");
        parser.addArgument("-q", "--mapq")
                .help("min mapping quality,  default 0").required(false).setDefault(0).dest("mapq");
        parser.addArgument("-m", "--mincount")
                .help("min TOTAL counts for given transcript, default 5").required(false).setDefault(5).dest("mincount");
        parser.addArgument("-i", "--ignore").help("analyze transcripts without name,  default true").required(false).setDefault(true).dest("ignore");
        parser.addArgument("-S", "--SamFlags")
                .help("COnsider sam-flags for multi mapping exclusion,  default true checks SAM flags for multi mapping indication").required(false).setDefault(true).dest("samflag");
        parser.addArgument("-Q", "--QualityStringency")
                .help("Set Stringency for output quality.  1 ignores CIGAR string problems. 2 ignores mapping issues and CIGAR, default 0 cuts >1 mutations in CIGAR and uses quality thresolds").required(false).setDefault(0).dest("stringency");

        Namespace ns = null;

        try {
            ns = parser.parseArgs(args);

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

            //System.out.println(offset + "  " + mincount);

            outfileLocation = ns.get("out");

        } catch (NumberFormatException num) {
            System.out.println("[ERROR] Please add numerical values to the parameters -o (offset) and -m (minCount) ");
            num.printStackTrace();
        } catch (ArgumentParserException e) {
            parser.handleError(e);
            System.exit(1);
        }

        try {
            bedHandler = new BEDHandler(ns.get("inBED").toString(), "BED", ns.get("ignore"));
            bamHandler = new BamHandler(ns.get("inBAM").toString(), "BAM", ns.get("samflag"), offset, qualStringency);
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



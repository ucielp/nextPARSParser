package eu.ernstthuer;

import java.util.ArrayList;
import java.util.List;

public class Main {


    public static void main(String[] args) {
        System.out.println("[STATUS] Running nextPARS parser v0.69");

        // get arguments
        ArgParser parser = new ArgParser(args);
        //System.out.println("[STATUS] Running with parameters  offset : " + parser.getOffset() + "   minCount for Coverage : " + parser.getMincount() + "  minQual for reads : " + parser.getMinqual());
        System.out.println("[STATUS] Running with parameters  offset : " + parser.getOffset() + "   minCount for Coverage : " + parser.getMincount() + "   SAMflag : " + parser.getSamflag() + "  minQual for reads : " + parser.getMinqual());

        // Handlers take care of the input files Bed and BAMoverage
        BEDHandler bedHandler = parser.getBedHandler();
        BamHandler bamHandler = parser.getBamHandler();

        // transfer the gene list to main
        ArrayList<Gene> geneArrayList = bedHandler.getGeneList();

        // load the chromosomeList.  This speeds up the BAM handling by limiting positions to chromosomes.
        bamHandler.setChromosomeArrayList(bedHandler.getChromosomeList());

        // configure the bam file handler
        bamHandler.setMinqual(parser.getMinqual());

        // read the gene List for coverage
        bamHandler.readLocus(geneArrayList);

        // write to file
        OutfileWriter outfileWriter = new OutfileWriter(parser.getOutfileLocation(), geneArrayList, parser.getMincount());

    }
}

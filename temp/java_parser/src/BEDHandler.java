package eu.ernstthuer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import static java.lang.Integer.parseInt;

/**
 * Created by ethur on 5/11/17.
 * Bedfile generator   parses BED text input
 *
 */
class BEDHandler extends FileHandler {


    private String locale;
    private String type;
    private static String origin;
    private String direction;
    private String feature;
    private String[] lineList;
    private ArrayList<Gene> geneList = new ArrayList<>();
    private ArrayList<Chromosome> chromosomeList = new ArrayList<>();



    BEDHandler(String locale, String type, boolean ignoreDots) {

        super(locale, type, "INPUT");
        this.locale = locale;
        this.origin = null;
        this.direction = direction; // BED is always input anyways
        try {
            this.lineList = openGFF(locale, ignoreDots);
        } catch (IOException e) {
            System.out.println("[ERROR] BED file not found");
            System.out.println(e);
        }
        //System.out.println("Total features " + lineList.length);
        // transfer text input to Arraylist of objects

        geneList = populateGeneList(this.lineList);

    }

     ArrayList<Gene> getGeneList() {
        return geneList;
    }

    String[] openGFF(String locale, boolean ignoreDots) throws IOException {
        ArrayList<String> outList = new ArrayList<String>();
        //System.out.println("This is where the file is " + direction);
        try (BufferedReader br = new BufferedReader(new FileReader(locale))) {

            String sCurrentLine;






            while ((sCurrentLine = br.readLine()) != null) {


                if(ignoreDots &&  sCurrentLine.split("\t")[3].charAt(0) != '.' )
                    outList.add(sCurrentLine);

            }



        } catch (IOException e) {
            e.printStackTrace();
        }

        String[] stockArr = new String[outList.size()];
        String[] linesOfFeatures = outList.toArray(stockArr);

        return linesOfFeatures;

    }

    ArrayList<Chromosome> getChromosomeList() {
        return chromosomeList;
    }

    ArrayList<String> chromosome_Temp = new ArrayList<>();


    ArrayList<Gene> populateGeneList(String[] featureList) {

        ArrayList<Gene> outList = new ArrayList<>();


        for (int i = 0; i < featureList.length; i++) {
            if (!featureList[i].startsWith("#")) {
                String[] row = featureList[i].split("\t");

                //ref,start,end,name,score,strand = l.split("\t")[:6]


                String chromosome = row[0];

                if(! chromosome_Temp.contains(chromosome)){
                    chromosomeList.add(new Chromosome(chromosome));
                    chromosome_Temp.add(chromosome);
                }

                String name = row[3];


                char orientation = row[5].charAt(0);
                //String description = (row[8]);
                int start = parseInt(row[1]);  //start and stop position are read as String
                int stop = parseInt(row[2]);

                double score;
                try {
                    score = Double.parseDouble(row[4]);
                } catch (NumberFormatException e) {
                    score = 0.0;
                  // value is missing in fields,  irrelevant, set to zero
                }

                try {


                    Gene gene = new Gene(chromosome, start, stop, orientation, name, score);

                    outList.add(gene);

                } catch (Exception e) {
                    System.out.println(e);
                }
                //System.out.println(outList.size());//geneList.add(newGene);

            }
            //System.out.println(outList.size());
        }
        System.out.println("[STATUS] " + outList.size() + " Transcripts added");
        return outList;

    }


}

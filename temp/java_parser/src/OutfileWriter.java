package eu.ernstthuer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by ethur on 5/12/17.
 * simply write to file
 * uses the Genes .toString method
 */
class OutfileWriter {

    private String locale;
    private ArrayList<Gene> resultList;
    private double minimumQuality;

    BufferedWriter bw = null;

    OutfileWriter(String locale, ArrayList<Gene> resultList, double minimumQuality) {
        this.locale = locale;
        this.resultList = resultList;
        this.minimumQuality = minimumQuality;
        writeToFile();
    }

    public void setMinimumQuality(int minimumQuality) {

        this.minimumQuality = minimumQuality;

    }

    private void writeToFile() {


        File file = new File(locale);

        BufferedWriter bw = null;

        try {
            if (!file.exists()) {
                file.createNewFile();
            }

            FileWriter fw = new FileWriter(file);
            bw = new BufferedWriter(fw);


            int geneCount = 0;
            int writeCount = 0;
            for (Gene gene : resultList) {

                geneCount++;
                if (gene.qualityScoring() >= minimumQuality) {
                    bw.write(gene.toString() + '\n');
                    writeCount++;
                }

            }
            System.out.println("[FINISHED] Output " + writeCount + " transcripts written to file in total");

        } catch (IOException e) {
            System.out.println("[ERROR] Error in writing to file  " + e.getCause());

        } finally {
            try {
                if (bw != null)
                    bw.close();
            } catch (Exception ex) {
                System.out.println("[ERROR] Error in closing the BufferedWriter" + ex);
            }
        }
    }
}




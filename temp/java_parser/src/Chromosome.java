package eu.ernstthuer;

import java.util.ArrayList;

/**
 * Chromosomes to store partial position Lists
 * helps keep the comparison low
 */
 class Chromosome {

    private String name;
    private ArrayList<Position> localPositionList = new ArrayList<>();

     Chromosome(String name) {
        this.name = name;
        ArrayList<Position> localPositionList = new ArrayList<>();
    }

     String getName() {
        return name;
    }

     ArrayList<Position> getLocalPositionList() {
        return localPositionList;
    }
}

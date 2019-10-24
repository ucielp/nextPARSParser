package eu.ernstthuer;

/**
 * Created by ethur on 5/11/17.
 */
class Position implements Comparable<Position> {


    /**
     * Store observed positions as comparable objects
     */

    private int position;
    private int coverage = 0;
    private char orientation;

    private String chromosome;

    Position(int position, String chromosome) {
        this.position = position;
        this.chromosome = chromosome;
        this.orientation = '.';
    }

    Position(int position, String chromosome, char Orientation) {
        this.position = position;
        this.chromosome = chromosome;
        this.orientation = Orientation;
    }

    int getCoverage() {
        return coverage;
    }

    int getPosition() {
        return position;
    }

    void setCoverage(int coverage) {
        this.coverage = coverage;
    }

    String getChromosome() {
        return chromosome;
    }

    void incrementCoverage() {
        coverage++;
    }

    char getOrientation() {
        return orientation;
    }

    @Override
    public int compareTo(Position o) {
        if (this.position == o.position && this.chromosome == o.chromosome) {
            return 1;
        } else {
            return 0;
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Position position1 = (Position) o;

        if (position != position1.position) return false;
        if (orientation != position1.orientation) return false;
        return chromosome != null ? chromosome.equals(position1.chromosome) : position1.chromosome == null;
    }

    @Override
    public int hashCode() {
        int result = position;
        result = 31 * result + (int) orientation;
        result = 31 * result + (chromosome != null ? chromosome.hashCode() : 0);
        return result;
    }

    void factorOffset(int offset){
        this.position = (this.position + offset);
    }

}

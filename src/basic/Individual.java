package basic;

/**
 * Each value of the chromosome is normalized into the range [0,1]
 *
 * @author cuonglv.hust@gmail.com
 * @date 24/02/2021
 *
 */
public class Individual implements Comparable<Individual> {

    private static int counter = 0; //Đếm id của cá thể
    private int individual_id;

    private int dim;
    private double[] chromosome; //Gen
    private double fitness;

    public Individual(int dim) {
        this.individual_id = Individual.counter++;
        this.dim = dim;
        this.chromosome = new double[dim];
        this.fitness = Double.MAX_VALUE;
    }

    public Individual(double[] chromosome) {
        this.individual_id = Individual.counter++;
        this.dim = chromosome.length;
        this.chromosome = chromosome.clone();
        this.fitness = Double.MAX_VALUE;
    }

    public void randomInit() {
        for (int i = 0; i < dim; i++) {
            this.chromosome[i] = Params.rand.nextDouble();
        }
    }

    public void setFitness(double value) {
        this.fitness = value;
    }

    public double getFitness() {
        return this.fitness;
    }

    public int getID() {
        return this.individual_id;
    }
    public void setID(int ID){
        this.individual_id = ID;
    }

    public int getDim() {
        return dim;
    }

    public void setDim(int dim) {
        this.dim = dim;
    }

    public double[] getChromosome() {
        return chromosome;
    }

    public void setChromosome(double[] chromosome) {
        this.chromosome = chromosome;
    }

    public double getGene(int index) {
        return this.chromosome[index];
    }

    public void setGene(int pos, double val) {
        this.chromosome[pos] = Math.max(0, Math.min(1, val));
    }

    @Override
    public int compareTo(Individual o) {
        return Double.valueOf(this.fitness).compareTo(o.fitness);
    }
}

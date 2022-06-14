package benchmark;

/**
 * @author cuonglv.hust@gmail.com
 * @date 24/02/2021
 *
 */
public abstract class Function {

    public int dim;	// dimension
    public double[] LB;
    public double[] UB;

    public double[] bias;
    public double[][] matrix;

    public Function() {

    }

    public Function(int dim, double[] shift, double[][] rotation, double[] lb, double[] ub) {
        this.dim = dim;
        this.bias = shift;
        this.matrix = rotation;
        this.UB = ub;
        this.LB = lb;
    }

    public abstract double getValue(double[] solution);

}

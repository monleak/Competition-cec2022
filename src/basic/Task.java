package basic;

import benchmark.Function;

/**
 * @author cuonglv.hust@gmail.com
 * @date 24/02/2021
 *
 */
public class Task {

    public int task_id;
    public Function function;

    public Task(int task_id, Function function) {
        this.task_id = task_id;
        this.function = function;
    }

    public Task(int task_id) {
        this.task_id = task_id;
    }

    public void setFunction(Function function) {
        this.function = function;
    }

    public double calculateFitnessValue(double[] x) {

        if (Params.countEvals > Params.maxEvals) {
            return Double.MAX_VALUE;
        }
        Params.countEvals++;

        double[] de_normalized = new double[this.function.dim]; //Dãy số thu được sau khi mã hóa NST từ không gian chung ra không gian riêng
        for (int i = 0; i < this.function.dim; i++) {
            de_normalized[i] = x[i] * (this.function.UB[i] - this.function.LB[i]) + this.function.LB[i];
        }

        double f = function.getValue(de_normalized);
        return f;
    }

}

package basic;

import java.util.Random;

/**
 * @author cuonglv.hust@gmail.com
 * @date 24/02/2021
 *
 */
public class Params {
    //Lưu các tham số
    public static Random rand; //Được khởi tạo theo các seed
    public static int countEvals; //Đếm số lần đánh giá
    public static int maxEvals; //Số lần đánh giá tối đa

    public static int recordsNum; //Số lần ghi lại kết quả
    public static final int REPT = 30; //Số lần chạy
    public static final int MAX_EVALS_PER_TASK = 100000;
    public static final double EPSILON = 5e-7; //điểm tối ưu (coi như bằng 0)

    public static final int H = 30; //Số phần tử của mem_cr, mem_f
    public static final double BEST_RATE = 0.11; //top.size() gấp bao nhiêu lần pop_size 0.11
    public static final double ARC_RATE = 5; //archive.size() gấp bao nhiêu lần pop_size 5
    public static final int MIN_POP_SIZE = 4; //Số cá thể tối thiểu mỗi tác vụ
    public static final int MAX_POP_SIZE = 100; //Số cá thể tối đa của mỗi tác vụ

    public static final double C = 0.02; //Sử dụng để update Rmp
    public static final double INIT_RMP = 0.5;

    public static final double NC = 2;


    public static String MANY_TASKS_BENCHMARKS_PATH = "SO-Manytask-Benchmarks/Tasks/Benchmark_";
    public static String COMPLEX_TASKS_BENCHMARKS_PATH = "SO-Complex-Benchmarks/Tasks/Benchmark_";
}

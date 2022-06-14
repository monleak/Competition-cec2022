package basic;

import java.util.ArrayList;

/**
 * @author cuonglv.hust@gmail.com
 * @date 24/02/2021
 *
 */
public class Problem {
    //Đại diện cho 1 ID benchmark
    public int DIM; //Số chiều trong không gian chung
    public int TASKS_NUM; //Số lượng tác vụ có trong ID
    public ArrayList<Task> tasks; //Danh sách các tác vụ

    public Problem() {
        tasks = new ArrayList<Task>();
        DIM = 0;
        TASKS_NUM = 0;
    }

    public void addTask(Task task) {
        tasks.add(task);
        TASKS_NUM++;
        DIM = Math.max(DIM, task.function.dim);
    }

    public Task getTask(int task_id) {
        return tasks.get(task_id);
    }
}

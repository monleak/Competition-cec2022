package basic;

import java.util.ArrayList;

public class ProblemManager {

    public ArrayList<Problem> problems;

    public ProblemManager(ArrayList<Problem> problems) {
        this.problems = problems;
    }

    public Problem getProblem(int index) {
        return problems.get(index);
    }
}

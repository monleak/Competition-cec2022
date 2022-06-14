package ls;

import basic.Individual;
import basic.Task;

public interface LocalSearch {

    Individual search(Individual start_point, int fes, Task task);

}

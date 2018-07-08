from Optimization.Tasks.UnconstrainedOptimization import UnconstrainedOptimization
from Optimization.Terminators.MaxTimeTerminator import MaxTimeTerminator
from Numerical_Objects.Vector import Vector

from math import *


task_1 = {
    'f': '1 * x ** 2',
    'variables': ['x']
}
f_1 = UnconstrainedOptimization.from_dict(task_1)
v_1 = Vector({'x': 0.0})
area_1 = {
    'x': (-10.0, 10.0)
}

task_2 = {
    'f': '1 * x ** 2 + 2 * y ** 2',
    'variables': ['x', 'y']
}
f_2 = UnconstrainedOptimization.from_dict(task_2)
v_2 = Vector({'x': 0.0, 'y': 0.0})
area_2 = {
    'x': (-10.0, 10.0),
    'y': (-10.0, 10.0)
}

task_3 = {
    'f': '1 * x ** 2 + 2 * y ** 2 + 3 * z ** 2',
    'variables': ['x', 'y', 'z']
}
f_3 = UnconstrainedOptimization.from_dict(task_3)
v_3 = Vector({'x': 0.0, 'y': 0.0, 'z': 0.0})
area_3 = {
    'x': (-10.0, 10.0),
    'y': (-10.0, 10.0),
    'z': (-10.0, 10.0)
}

mt = MaxTimeTerminator('m:1,s:30')


def solution_delta(v1, v2):
    s1 = v1.reduce_to_solution()
    s2 = v2.reduce_to_solution()

    keys = s1.keys() | s2.keys()
    delta = 0.0
    for k in keys:
        delta += (s1[k] - s2[k]) ** 2
    delta = sqrt(delta)

    return delta


class Verifier:
    def __init__(self, tolerance=1e-3, attempts=5):
        self.test_functions = [f_1, f_2, f_3]
        self.search_area = [area_1, area_2, area_3]
        self.solutions = [v_1, v_2, v_3]
        self.tolerance = tolerance
        self.attempts = attempts

    def verify(self, algorithms, logger):
        for i_f, f in enumerate(self.test_functions):
            logger.info('>>> Testing function {0}/{1}'.format(i_f + 1, len(self.test_functions)))
            success = False
            for i_a, a in enumerate(algorithms):
                logger.info('>>> >>> Using config {0}/{1}'.format(i_a + 1, len(algorithms)))
                for j in range(self.attempts):
                    logger.info('>>> >>> >>> Attempt {0}/{1}'.format(j + 1, self.attempts))
                    x_opt = a.work(f, self.search_area[i_f], mt)
                    if solution_delta(x_opt, self.solutions[i_f]) < self.tolerance:
                        success = True
                        break
                if success:
                    break
            if not success:
                return False
        return True
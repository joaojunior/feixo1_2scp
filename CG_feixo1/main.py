import sys
import cplex

import my_io
import colect_time

class SetCoveringProblem(object):
    def __init__(self, costs, incidence_matrix, columns):
        self.costs = costs
        self.incidence_matrix = incidence_matrix
        self.columns = columns
        self.num_variables = len(self.costs)
        self.names_variables = []

    def create_cplex_problem(self):
        self.c = cplex.Cplex()
        self.c.set_results_stream(None)
        self.c.set_log_stream(None)
        self.__add_variables()
        self.__create_initial_constraints()
        self.__add_function_objective()
        self.c.set_problem_type(self.c.problem_type.LP)
        self.c.parameters.timelimit.set(2*3600)

    def __add_variables(self):
        self.names_variables = []
        for i in xrange(self.num_variables):
            self.names_variables.append('x%d' %(i+1))
        self.c.variables.add(names = self.names_variables, types = [self.c.variables.type.integer] * self.num_variables,
                            ub=[1]*self.num_variables)

    def __create_initial_constraints(self):
        linear_constraints = []
        for linear_constraint in self.incidence_matrix:
            linear_constraints.append(cplex.SparsePair(ind = self.names_variables, val = linear_constraint))
        self.c.linear_constraints.add(lin_expr = linear_constraints,
                                senses = ["G"] * len(linear_constraints),
                                rhs = [1] *len(linear_constraints))

    def __add_function_objective(self):
        objective = []
        for i in xrange(self.num_variables):
            objective.append((self.names_variables[i], self.costs[i]))
        self.c.objective.set_linear(objective)
        self.c.objective.set_sense(self.c.objective.sense.minimize)

    def create_cg_separate(self, ro=0.01):
        self.cg_separate = cplex.Cplex()
        self.cg_separate.set_results_stream(None)
        self.cg_separate.set_log_stream(None)
        self.variables_names_alpha = ['alpha0']
        self.variables_names_u = []
        for i in xrange(self.num_variables):
            self.variables_names_alpha.append(('alpha%d' %(i+1)))
        for i in xrange(len(self.incidence_matrix)):
            self.variables_names_u.append(('u%d' %(i+1)))
        self.cg_separate.variables.add(names = self.variables_names_alpha, types = [self.cg_separate.variables.type.integer]*len(self.variables_names_alpha))
        self.cg_separate.variables.add(names = self.variables_names_u, types = [self.cg_separate.variables.type.continuous]*len(self.variables_names_u),
                                        lb=[0]*len(self.variables_names_u),ub=[1-ro]*len(self.variables_names_u))
        linear_constraints = []
        for i in xrange(len(self.columns)):
            linear_constraints.append(cplex.SparsePair(ind = ['alpha%d' %(i+1)] + self.variables_names_u, val = [1] + self.columns[i]))
            linear_constraints.append(cplex.SparsePair(ind = ['alpha%d' %(i+1)] + self.variables_names_u, val = [1] + self.columns[i]))
        linear_constraints.append(cplex.SparsePair(ind = ['alpha0'] + self.variables_names_u, val = [1] + [-1] * len(self.variables_names_u)))
        linear_constraints.append(cplex.SparsePair(ind = ['alpha0'] + self.variables_names_u, val = [1] + [-1] * len(self.variables_names_u)))
        self.cg_separate.linear_constraints.add(lin_expr = linear_constraints,
                            senses = ["G", "L"] * (len(linear_constraints)/2),
                            rhs = [0, 1-ro] * (len(linear_constraints)/2))
                            
        self.cg_separate.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = self.variables_names_alpha, val = [1] + [-1]*(len(self.variables_names_alpha)-1))],
                            senses = ["L"],
                            rhs = [0])
        self.cg_separate.objective.set_sense(self.cg_separate.objective.sense.minimize)

    def add_constraint_in_problem(self, alphas):
        self.c.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = self.names_variables, val = alphas[1])],
                                senses = ["G"],
                                rhs = [alphas[0]])

    def resolve_cg_separate(self, x, w=0.0001):
        objective = [('alpha0', -1)]
        for idx in xrange(len(x)):
            objective.append((self.variables_names_alpha[idx+1],x[idx]))
        for name_variable_u in self.variables_names_u:
            objective.append((name_variable_u,-w))
        self.cg_separate.objective.set_linear(objective)
        self.cg_separate.parameters.mip.tolerances.uppercutoff.set(-0.05)
        self.cg_separate.parameters.mip.limits.solutions.set(5)
        self.cg_separate.parameters.mip.tolerances.mipgap.set(0.5)
        self.cg_separate.parameters.timelimit.set(2*60)
        self.cg_separate.parameters.mip.limits.nodes.set(50000)
        self.cg_separate.parameters.emphasis.mip.set(4)
        self.cg_separate.solve()
        solution = self.cg_separate.solution
        result = []
        if solution.get_status() == 101 or solution.get_status() == 102 or solution.get_status() == 105 or solution.get_status() == 107 or solution.get_status() == 104:
            for name_variable_alpha in self.variables_names_alpha[1:]:
                result.append(solution.get_values(name_variable_alpha))
            return solution.get_values('alpha0'), result
        else:
            return (0, [0] * self.num_variables)

def verify_list_parameter_is_integer(x):
    result = True
    for i in x:
        if i != int(i):
            return False
    return result

def transform_incidence2column(incidence_matrix):
    numbers_column = len(incidence_matrix[0])
    columns = []
    for i in xrange(numbers_column):
        column = []
        for item in incidence_matrix:
            column.append(-item[i])
        columns.append(column)
    return columns
    
def execute(name_file):
    costs, incidence_matrix = my_io.read_file_format_or_library(name_file)
    columns = transform_incidence2column(incidence_matrix)
    set_covering = SetCoveringProblem(costs=costs, incidence_matrix=incidence_matrix,columns=columns)
    set_covering.create_cplex_problem()
    set_covering.c.solve()
    solution = set_covering.c.solution
    alphas = (1, [0]*set_covering.num_variables)
    constraint_add = []
    relax_linear = solution.get_objective_value()
    start = colect_time.cpu_time()
    while (not verify_list_parameter_is_integer(solution.get_values()) and alphas != (0, [0]*set_covering.num_variables) and colect_time.cpu_time() - start <= 7200):
        set_covering.create_cg_separate()
        alphas = set_covering.resolve_cg_separate(solution.get_values())
        if alphas != (0, [0]*set_covering.num_variables):
            constraint_add.append(alphas)
            set_covering.add_constraint_in_problem(alphas)
            set_covering.c.solve()
    end = colect_time.cpu_time()
    print name_file, relax_linear, solution.get_objective_value(), len(constraint_add), end - start
    
if __name__ == '__main__':
    if len(sys.argv) == 2:
      name_file = sys.argv[1]
      execute(name_file)

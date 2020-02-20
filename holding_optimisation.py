import sympy as sy

def print_expression(exp_name, exp):
    print(exp_name + ": ")
    print(exp)

# Inputs
h_0 = sy.MatrixSymbol("h_0", 2, 1)
alpha = sy.MatrixSymbol("alpha", 2, 1)
v_1_1, v_1_2, v_2_2 = sy.symbols("v_1_1, v_1_2, v_2_2")
# Covarience
V =  sy.Matrix([[v_1_1, v_1_2], [v_1_2, v_2_2]])
c = sy.MatrixSymbol("c", 2, 1)
C = sy.Matrix([[c[0,0], 0], [0, c[1,0]]])
sigma = sy.symbols("sigma")

# Variables
t_0, t_1 = sy.symbols("t_0 t_1")
t = sy.Matrix([[t_0], [t_1]])
lam = sy.symbols("lambda")
h = h_0 + t

objective = alpha.T * h - t.T * C * t
constraint = h.T * V * h
utility = sy.simplify((objective - lam * constraint)[0,0])
print_expression("utility", utility)
diff_system = [sy.simplify(sy.diff(utility, t_0)), sy.simplify(sy.diff(utility, t_1))]
print_expression("system", diff_system)
t_solutions = sy.nonlinsolve(diff_system, [t_0, t_1])

# Prepare a 4 degree polynome solution
a_0, a_1, a_2, a_3, a_4 = sy.symbols("a_0 a_1 a_2 a_3 a_4")
eq = sy.Eq(a_4 * lam ** 4 + a_3 * lam ** 3 + a_2 * lam ** 2 + a_1 * lam ** 1 + a_0)
polynome_solutions = sy.solveset(eq, lam)
print_expression("polynome_solutions", polynome_solutions)

candidates_count = 0
max_objective_result = None

for sol in t_solutions:
    sol_vector = sy.Matrix([[sy.simplify(sol[0])], [sy.simplify(sol[1])]])
    print_expression("t solution", sol_vector)
    variance = sy.simplify((sol_vector.T * V * sol_vector)[0,0])
    print_expression("variance", variance)
    print(variance)
    denominator = 1 / variance.args[1]
    denominator = sy.collect(sy.expand(denominator), lam)
    print_expression("variance denominator", denominator)
    variance_normalised = variance * denominator
    variance_normalised = sy.collect(sy.expand(variance_normalised), lam)
    print_expression("variance_normalised", variance_normalised)
    polynome = variance_normalised - sigma * denominator
    for polynome_sol in polynome_solutions:
        lambda_solution = polynome_sol.subs([
            (a_0, polynome.coeff(lam, 0)),
            (a_1, polynome.coeff(lam, 1)),
            (a_2, polynome.coeff(lam, 2)),
            (a_3, polynome.coeff(lam, 3)),
            (a_4, polynome.coeff(lam, 4))
        ])
        objective_result = objective.subs(t, sol_vector.subs(lam, lambda_solution))
        candidates_count = candidates_count + 1
        print_expression("objective_result (" + str(candidates_count) + ")", objective_result)
        if max_objective_result is not None:
            max_objective_result = sy.Max(max_objective_result, objective_result)
        else:
            max_objective_result = objective_result

print_expression("max_objective_result", max_objective_result)

# Tests
print_expression("test_1", max_objective_result.subs(
    [
        (h_0, sy.Matrix([[0], [0]])),
        (alpha, sy.Matrix([[0], [0]])),
        (v_1_1, 1),
        (v_1_2, 0),
        (v_2_2, 0),
        (c, sy.Matrix([[1], [0]])),
        (sigma, 1)
    ]
))


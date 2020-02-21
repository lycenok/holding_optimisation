import sympy as sy


def print_expression(exp_name, exp):
    print(exp_name + ": ")
    print(exp)


# Inputs
h_0 = sy.MatrixSymbol("h_0", 2, 1)
# TODO: delete
h_0 = sy.Matrix([[0], [0]])
alpha = sy.MatrixSymbol("alpha", 2, 1)
# TODO: delete
alpha = sy.Matrix([[1], [1]])
v_1_1, v_1_2, v_2_2 = sy.symbols("v_1_1, v_1_2, v_2_2")
# Covarience
V = sy.Matrix([[v_1_1, v_1_2], [v_1_2, v_2_2]])
# TODO: delete
V = V.subs([
    (v_1_1, 1),
    (v_1_2, 0.5),
    (v_2_2, 1)
])
c = sy.MatrixSymbol("c", 2, 1)
C = sy.Matrix([[c[0,0], 0], [0, c[1,0]]])
# TODO: delete
C = sy.Matrix([[1,0], [0,1]])
sigma = sy.symbols("sigma")
# TODO: delete
sigma = 1

# Variables
t_0, t_1 = sy.symbols("t_0 t_1")
t = sy.Matrix([[t_0], [t_1]])
lam = sy.symbols("lambda")
h = h_0 + t

objective = (alpha.T * h - t.T * C * t)[0,0]
constraint = sy.simplify(sy.expand((h.T * V * h)[0,0]))
utility = sy.simplify(objective - lam * constraint)
print_expression("utility", utility)
diff_system = [sy.simplify(sy.diff(utility, t_0)), sy.simplify(sy.diff(utility, t_1))]
print_expression("system", diff_system)
t_solutions = sy.linsolve(diff_system, [t_0, t_1])

# Prepare a 4 degree polynome solution
a_0, a_1, a_2, a_3, a_4 = sy.symbols("a_0 a_1 a_2 a_3 a_4")
eq = sy.Eq(a_4 * (lam ** 4) + a_3 * (lam ** 3) + a_2 * (lam ** 2) + a_1 * (lam ** 1) + a_0, 0)
# TODO: do only real numbers sy.Reals (raises an excception now when evaluating objective)
polynome_solutions = sy.solveset(eq, lam)
eq = eq.subs([(a_4, 1), (a_3, 0), (a_2, 0), (a_1, 0), (a_0, -1)])
print_expression("eq", eq.subs([(a_4, 1), (a_3, 0), (a_2, 0), (a_1, 0), (a_0, -1)]))
print_expression("polynome_solutions", polynome_solutions)

candidates = []
max_objective_result = None

for sol in t_solutions:
    sol_vector = sy.Matrix([[sy.simplify(sol[0])], [sy.simplify(sol[1])]])
    print_expression("t solution", sol_vector)
    variance = sy.simplify((sol_vector.T * V * sol_vector)[0,0])
    print_expression("variance", variance)
    if len(variance.args) > 1:
        denominator = 1 / variance.args[1]
    else:
        denominator = 1
    denominator = sy.collect(sy.expand(denominator), lam)
    print_expression("variance denominator", denominator)
    variance_normalised = variance * denominator
    variance_normalised = sy.collect(sy.expand(variance_normalised), lam)
    print_expression("variance_normalised", variance_normalised)
    polynome = sy.collect(variance_normalised - sigma * denominator, lam)
    print_expression("polynome", polynome)
    print_expression("polynome.coeff(lam, 0)", polynome.coeff(lam, 0))
    print_expression("polynome.coeff(lam, 2)", polynome.coeff(lam, 2))
    for polynome_sol in polynome_solutions:
        lambda_solution = polynome_sol.subs([
            (a_0, polynome.coeff(lam, 0)),
            (a_1, polynome.coeff(lam, 1)),
            (a_2, polynome.coeff(lam, 2)),
            (a_3, polynome.coeff(lam, 3)),
            (a_4, polynome.coeff(lam, 4))
        ])
        sol_vector_pure = sol_vector.subs(lam, lambda_solution)
        objective_result = objective.subs([(t_0, sol_vector_pure[0,0]), (t_1, sol_vector_pure[1,0])])
        candidates.append(objective_result)
        # TODO: apply maximum and check if it is real
        print_expression("lambda_solution (" + str(len(candidates)) + ")", lambda_solution)
        print_expression("sol_vector_pure (" + str(len(candidates)) + ")", sol_vector_pure)
        print_expression("objective_result (" + str(len(candidates)) + ")", objective_result)

# Tests
for i in range(len(candidates)):
    print_expression("test_1 (" + str(i) + ")", candidates[i].subs(
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


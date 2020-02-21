import sympy as sy

lam = sy.symbols("lambda")
a_0, a_1, a_2, a_3 = sy.symbols("a_0 a_1 a_2 a_3")
eq = sy.Eq(a_3 * (lam ** 3) + a_2 * (lam ** 2) + a_1 * (lam ** 1) + a_0, 0)
polynome_solutions = sy.solve(eq, lam)
print(polynome_solutions)
for solution in  polynome_solutions:
    print(solution.subs([(a_3, 0), (a_2, 1), (a_1, 0), (a_0, -1)]))
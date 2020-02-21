import sympy as sy


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
c = sy.Matrix([[1], [3]])
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

# y = sqrt(C) * t
# t = sqrt(C)**-1 * y
y = sy.MatrixSymbol("y", 2, 1)
Cs = sy.Matrix([[sy.sqrt(c[0,0]), 0], [0, sy.sqrt(c[1,0])]])

t = Cs.inv() * y
h = Cs.inv() * y + h_0
# h_0 = Cs.inv() * h_0_2
h_0_2 = Cs * h_0
objective = (alpha.T * (Cs.inv() * y + h_0) - (Cs.inv() * y).T * C * (Cs.inv() * y))[0,0]
# objective = alpha.T * Cs.Inv() * y + y.T * y

constraint = sy.simplify(sy.expand(((Cs.inv() * y + Cs.inv() * h_0_2).T * V * (Cs.inv() * y + Cs.inv() * h_0_2))[0,0]))
constraint = sy.simplify(sy.expand(((Cs.inv() * (y + h_0_2)).T * V * (Cs.inv() * (y + h_0_2)))[0,0]))
V_2 = Cs * V * Cs.inv()
constraint = sy.simplify(sy.expand( (y + h_0_2).T * V_2 * (y + h_0_2) ))[0,0]

print (objective)
print (constraint)

print (V_2)
(P, D) = V_2.diagonalize()
print (P)
print (D)

# y = P.inv() * x
x = sy.MatrixSymbol("x", 2, 1)
objective = alpha.T * Cs.inv() * P.inv() * x - x.T * x
D_2 = sy.Matrix([[sy.sqrt(D[0,0]), 0], [0, sy.sqrt(D[1,1])]])
y_for_constraint = D_2 * P * (P.inv() * x + h_0_2)
y_for_constraint = D_2 * ( x + P * h_0_2 )
constraint = (y_for_constraint.T * y_for_constraint)[0,0]

# z = D_2 * ( x + P * h_0_2 )
# x = D_2**(-1) * ( z - D_2 * h_0_2 )
z = sy.MatrixSymbol("z", 2, 1)
objective = alpha.T * Cs.inv() * P.inv() * D_2.inv() * z - (D_2.inv() * ( z - D_2 * h_0_2 )).T * (D_2.inv() * ( z - D_2 * h_0_2 ))
objective = ((alpha.T * Cs.inv() * P.inv() * D_2.inv() - 2 * h_0_2.T * D_2.inv()) * z - (D_2.inv() * z).T * (D_2.inv() * z))[0,0]
constraint = (z.T * z)[0,0] - sigma

print (objective)
print (constraint)


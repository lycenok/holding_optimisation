# Holding Optimisation

Optimise portfolio holdings for given predicted price for two assets and quadratic transactions cost

Problem:

```c1,c2,sigma,h_a(0),h_b(0)

r_a(1),r_b(1)

...

r_a(N),r_b(N)

PnL(n) = h_a(n) * r_a(n) + h_b(n) * r_b(n) - c1 * |t_a(n)| - c2 * t_a(n)^2 - c2 & |t_b(n)| - c2 * t_b(n)^2
t_a(n) = h_a(n) - h_a(n-1)
t_b(n) = h_b(n) - h_b(n-1)

1 <= n <= N

M(PnL) = (PnL(1) + ... + PnL(N)) / N
Var(PnL) = ((PnL(1) - M(PnL))^2 + (PnL(2)) - M(PnL))^2 + ... (PnL(N) - M(PnL))^2    

Find maximum of PnL(n) with the constraint Var(PnL) <= sigma^2,
t_a(1), ... t_a(n),t_b(1), ... t_b(n) in R (real numbers)

```


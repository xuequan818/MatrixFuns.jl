# Limitations

- No optimal strategy for selecting `sep`, the initial separation distance.

- No optimal strategy for choosing `scale`, the characteristic scale of variations of $f$. In practice, we select it so that the coefficient of `x` becomes unity. For instance, for `f(x) = exp(100x)`, we typically set `scale = 1/100`.

- Evaluating highly singular functions with large matrices may become numerically unstable. In order to enforce very tight Taylor‐series error tolerances, the eigenvalues can become excessively separated, which in turn can produce large errors when solving the Sylvester equation. 

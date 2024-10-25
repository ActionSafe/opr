An R-package for Fitting ordinally observed panel count regression models
(under development...)
# Usage
```R
fit = ordinalReg(formula = OrdinalSurv(id, time, count_obs, lower, upper) ~ x1 + x2,
           data = dd, method = "ML", se = "Fisher", control = list(max_iter = 100))
```

# BSW 0.1.2

## New
- Added option `conswitch` in `bsw()`:
  - `1` (default): all possible min/max combinations of predictors (bounds predictions).
  - `0`: raw design matrix (suitable for risk factor identification only).
- New function `obj_value()` to compute the log-likelihood value used in the Armijo line search.
- Implemented a robust bootstrap procedure `bootBSW()` for `bsw` objects. The function completes execution even if individual bootstrap samples fail to converge.
- Implemented variable selection methods (`variable_selection_bsw()`), supporting both backward elimination and forward selection. User can specify significance level (`alpha`) and maximum iterations (`maxit`).

## Improved
- Optimized Hessian matrix calculation by replacing explicit loops with vectorized operations.
- Armijo line search implemented to ensure stable convergence.    

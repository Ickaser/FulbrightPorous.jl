@quickactivate :FulbrightPores

const OP = OptimizationProblem
const NLSP = NonlinearLeastSquaresProblem
const BFGS = Optim.BFGS
using LineSearches
alg = BFGS(linesearch=LineSearches.BackTracking())
const LP = LyoPronto
pythonplot()


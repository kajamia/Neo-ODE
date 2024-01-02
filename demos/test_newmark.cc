#include <iostream>

#include <nonlinfunc.h>
#include <ode.h>

using namespace Neo_ODE;
using namespace Neo_CLA;
using namespace std;

class RHS : public NonlinearFunction
{
  size_t DimX() const override { return 1; }
  size_t DimF() const override { return 1; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = -x(0);
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df(0, 0) = -1;
  }
};


int main()
{
  double tend = 2*M_PI;
  int steps = 100;
  Vector<> x { 1, };
  Vector<> dx { 0., };
  auto rhs = make_shared<RHS>();
  auto mass = make_shared<IdentityFunction>(1);
  /* SolveODE_Newmark(tend, steps, x, dx, rhs, mass,
                   [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x(0) << endl; }
                   ); */

  SolveODE_Newmark(tend, steps, x, dx, rhs, mass,
                   [](double t, VectorView<double> x) { cout << x(0) << ","; }
                   );
}

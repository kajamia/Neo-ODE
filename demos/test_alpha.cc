#include <iostream>
#include <fstream>
#include <memory>

#include <nonlinfunc.h>
#include <ode.h>

using namespace Neo_ODE;
using namespace Neo_CLA;
using namespace std;

// the pendulum with a length constraint

// Lagrange = -f*y + lam*(x*x+y*y-1)
// dLagrange
class dLagrange : public NonlinearFunction
{
  size_t DimX() const override { return 3; }
  size_t DimF() const override { return 3; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = 2*x(0)*x(2);
    f(1) = 2*x(1)*x(2) - 1;
    f(2) = x(0)*x(0)+x(1)*x(1)-1;
    
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df(0,0) = 2*x(2);
    df(0,1) = 0;
    df(0,2) = 2*x(0);

    df(1,0) = 0;
    df(1,1) = 2*x(2);
    df(1,2) = 2*x(1);

    df(2,0) = 2*x(0);
    df(2,1) = 2*x(1);
    df(2,2) = 0;
  }
};

// langrange = -(f1*y1 + f2*y2) + l1*(x1^2 + y1^2 - 1) + l2*((x1-x2)^2 + (y1 - y2)^2 - 1)
// with f1=f2=1
// rhs = Lagrange derived
// Norbert has notes on the maths behind it
class dLagrangeDoublePendulum : public NonlinearFunction
{
  size_t DimX() const override { return 6; }
  size_t DimF() const override { return 6; }

  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    // x1…x(0)
    // y1…x(1)
    // x2…x(2)
    // y2…x(3)
    // l1…x(4)
    // l2…x(5)

    f = {2*(x(4)*x(0) + x(5)*x(0) - x(5)*x(2)),
         -1 + 2*(x(4)*x(1) + x(5)*x(1) - x(3)*x(5)),
         2*(-x(5)*x(0) + x(2)*x(5)),
         -1 + 2*(-x(1)*x(5) + x(3)*x(5)),
         x(0)*x(0) + x(1)*x(1) - 1,
         (x(0) - x(2))*(x(0) - x(2)) + (x(1) - x(3))*(x(1) - x(3)) - 1
         };
  }
  void EvaluateDeriv (VectorView<double> x_, MatrixView<double> df) const override
  {
    VectorView x = x_.Slice(0, 2);
    VectorView y = x_.Slice(1, 2);
    VectorView l = x_.Range(4, 6);

    df = {2*(l(0) + l(1)), 0, -2*l(1), 0, 2*x(0), 2*(x(0) - x(1)),
          0, 2*(l(0) + l(1)), 0, -2*l(1), 2*y(0), 2*(y(0) - y(1)),
          -2*l(1), 0, 2*l(1), 0, 0, -2*(x(0) - x(1)),
          0, -2*l(1), 0, 2*l(1), 0, -2*(y(0) - y(1)),
          2*x(0), 2*y(0), 0, 0, 0, 0,
          2*(x(0) - x(1)), 2*(y(0) - y(1)), -2*(x(0) - x(1)), -2*(y(0) - y(1)), 0, 0
          };
  }
};

int main()
{
  double tend = 2*2*M_PI;
  double steps = 1000;
  /* Vector<double> x { 1, 0, 2, 0, 0, 0 };
  Vector<double> dx { 0, 0, 0, 0, 0, 0 };
  Vector<double> ddx { 0, 0, 0, 0, 0, 0 };
  auto rhs = make_shared<dLagrangeDoublePendulum>();
  auto mass = make_shared<Projector>(6, 0, 4); */
  Vector<double> x { 1, 0, 0 };
  Vector<double> dx { 0, 0, 0 };
  Vector<double> ddx { 0, 0, 0 };
  auto rhs = make_shared<dLagrange>();
  auto mass = make_shared<Projector>(3, 0, 2);

  
  
  /* SolveODE_Alpha (tend, steps, 0.8, x, dx, ddx, rhs, mass, 
                   // [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x(0) << " " << x(1) << " " << x(2) << endl; }
                   [](double t, VectorView<double> x) { cout << t << " " << x(0) << " " << x(1) << " " << x(2) << endl; }                   
                   ); */
  
  SolveODE_Alpha (tend, steps, 0.8, x, dx, ddx, rhs, mass,
                   [](double t, VectorView<double> x){ cout << x(0) << ","; });
}

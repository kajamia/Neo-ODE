#include <iostream>
#include <cmath>

#include <nonlinfunc.h>
#include <ode.h>
#include <vector.h>

using namespace Neo_ODE;
using namespace Neo_CLA;
using namespace std;


class MassSpring : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};

void test_mass_spring()
{
  double tend = 4*M_PI;
  int steps = 100;
  Vector<> y { 1, 0 };
  auto rhs = make_shared<MassSpring>();
  
  SolveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { cout << "IE " << t << "  " << y(0) << " " << y(1) << endl; });
  
  SolveODE_EE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { cout << "EE " << t << "  " << y(0) << " " << y(1) << endl; });
  
  SolveODE_CN(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { cout << "CN " << t << "  " << y(0) << " " << y(1) << endl; });
}

void test_exponential()
{
  Vector<> y{1, 0};
  SolveODE_IE(0.5, 100, y, make_shared<ConstantFunction>(Vector<> {1, 2}),
              [](double t, VectorView<double> y) { cout << "IE " << t << "  " << y(0) << " " << y(1) << endl; });
  
  y = {0, 1};

  SolveODE_EE(0.5, 100, y, make_shared<ConstantFunction>(Vector<> {1, 2}),
              [](double t, VectorView<double> y) { cout << "EE " << t << "  " << y(0) << " " << y(1) << endl; });
  
  y = {0, 1};

  SolveODE_CN(0.5, 100, y, make_shared<ConstantFunction>(Vector<> {1, 2}),
              [](double t, VectorView<double> y) { cout << "CN " << t << "  " << y(0) << " " << y(1) << endl; });
}


int main()
{
  //test_mass_spring();
  test_exponential();
}

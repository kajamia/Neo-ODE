#include <iostream>
#include <cmath>
#include <pybind11/pybind11.h>

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

Matrix<> test_mass_spring()
{ 
  // 3 methods, 100 values each -> 300 values
  Matrix<> all_y (3*100, 2);
  all_y.Row(0) = {1, 0};
  all_y.Row(100) = {1, 0};
  all_y.Row(200) = {1, 0};

  double tend = 4*M_PI;
  int steps = 100;
  auto rhs = make_shared<MassSpring>();

  SolveODE_IE(tend, steps, all_y.Rows(0, 100), rhs);
  SolveODE_EE(tend, steps, all_y.Rows(100, 100), rhs);
  SolveODE_CN(tend, steps, all_y.Rows(200, 100), rhs);

  return all_y;  
}

PYBIND11_MODULE(ode, m)
{
  m.doc() = "just a test of ode methods";

  m.def("test_mass_spring", &test_mass_spring);
}


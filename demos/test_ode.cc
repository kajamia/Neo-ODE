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

Matrix<double, ColMajor> test_exponential_py()
{ 
  // 3 methods, 100 values each -> 300 values
  Matrix<double, ColMajor> all_y (2, 3*100);

  // every y argument needs to be {1, 0}
  all_y.Row(0) = 1;
  all_y.Row(1) = 0;

  for (int i = 0; i < 100; i += 3)
  {
    SolveODE_IE(0.5, 100, all_y.Col(i), make_shared<ConstantFunction>(Vector<> {1, 2}));

    SolveODE_EE(0.5, 100, all_y.Col(i + 1), make_shared<ConstantFunction>(Vector<> {1, 2}));

    SolveODE_CN(0.5, 100, all_y.Col(i + 2), make_shared<ConstantFunction>(Vector<> {1, 2}));
  }

  return all_y;  
}

PYBIND11_MODULE(ode, m)
{
  m.doc() = "just a test of ode methods";

  m.def("test_exponential", &test_exponential_py);
}


void test_exponential()
{
  Vector<> y{0, 1};
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

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

  // cout.precision(10);
  
  SolveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { cout << "IE " << t << " \t " << y(0) << " \t " << y(1) << endl; });
  
  SolveODE_EE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { cout << "EE " << t << " \t " << y(0) << " \t " << y(1) << endl; });
  
  SolveODE_CN(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { cout << "CN " << t << " \t " << y(0) << " \t " << y(1) << endl; });

  Matrix<> all_y_ie(100, 2);
  all_y_ie.Row(0) = {1, 0};
  SolveODE_IE(tend, steps, all_y_ie.View(), rhs);
  std::cout << all_y_ie << std::endl;

  Matrix<> all_y_ee(100, 2);
  all_y_ee.Row(0) = {1, 0};
  SolveODE_EE(tend, steps, all_y_ee.View(), rhs);
  std::cout << all_y_ee << std::endl;

  Matrix<> all_y_cn(100, 2);
  all_y_cn.Row(0) = {1, 0};
  SolveODE_CN(tend, steps, all_y_cn.View(), rhs);
  std::cout << all_y_cn << std::endl;
}

int main()
{
  test_mass_spring();
}

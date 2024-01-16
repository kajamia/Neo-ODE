#include <cmath>

#include <nonlinfunc.h>
#include <ode.h>

using namespace Neo_ODE;

class Electric: public NonlinearFunction
{
  double R_; // resistivity of resistor
  double C_; // capacity of capacitor

 public:
  Electric(double R, double C) : R_(R), C_(C) {};

  size_t DimX() const override {return 2;}
  size_t DimF() const override {return 2;}

  void Evaluate(VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = (std::cos(100*M_PI*x(1)) - x(0))/(R_*C_); // voltage at capacitor
    f(1) = 1; // time
  }

  void EvaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
  {
    df(0, 0) = -1/(R_*C_);
    df(1, 0) = 0;
    df(0, 1) = (-std::sin(100*M_PI*x(1)*100*M_PI))/(R_*C_);
    df(1, 1) = 1;
  }  
};

int main()
{
  double tend = 0.05;
  int steps = 1000;
  Vector<> y { 0, 0 };
  std::cout << "0,";

  auto rhs = make_shared<Electric>(100, 1e-6);
  
  SolveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << y(0) << ","; }); //{ std::cout << "IE " << t << " \t " << std::cos(100 * M_PI * t) << " \t " << y(1) << std::endl; });
}
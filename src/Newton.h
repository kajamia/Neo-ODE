#ifndef Newton_h
#define Newton_h

#include "nonlinfunc.h"
#include "matrix.h"

namespace Neo_ODE
{

  void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<> res(func->DimF());
    Matrix<> fprime(func->DimF(), func->DimX());

    // std::cout << "x = " << x << std::endl;
    for (int i = 0; i < maxsteps; i++)
      {
        func->Evaluate(x, res);
        // std::cout << "res = " << res << std::endl;
        // cout << "|res| = " << L2Norm(res) << endl;
        func->EvaluateDeriv(x, fprime);
        // std::cout << "fprime = " << fprime << std::endl;
        fprime = Inverse(fprime);
        // std::cout << "inv fprime = " << fprime << std::endl;
        x -= fprime*res;
        // std::cout << "new x = " << x << std::endl;

        double err= L2Norm(res);
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif

#ifndef ODE_h
#define ODE_h

#include <functional>
#include <exception>
#include <memory>

#include "Newton.h"


namespace Neo_ODE
{
  
  // implicit Euler method for dy/dt = rhs(y)
  void SolveODE_IE(double tend, int steps,
                   VectorView<double> y, shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;

    auto yold = make_shared<ConstantFunction>(y);
    auto ynew = make_shared<IdentityFunction>(y.Size());
    auto equ = ynew-yold - dt * rhs;

    double t = 0;

    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, y);
        yold->Set(y);
        t += dt;
        if (callback) callback(t, y);
      }
  }

  
  // explicit Euler method for dy/dt = rhs(y)
  void SolveODE_EE(double tend, int steps,
                   VectorView<double> y, shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double, VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    if (rhs->DimX() != y.Size() || rhs->DimX() != y.Size()){throw invalid_argument("rhs does not have the right dimensions"); }

    double t = 0;
    Vector<double> tmp(y.Size());

    for (int i = 0; i < steps; i++)
    {
      y(y.Size() - 1) = t;
      rhs->Evaluate(y, tmp);
      y = y + dt*tmp;

      t += dt;
      if (callback) callback(t, y);
    }
  }

  // Crank-Nicholson method
  void SolveODE_CN(double tend, int steps,
                   VectorView<double> y, shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double, VectorView<double>)> callback = nullptr)
  {
    // h
    double dt = tend/steps;
    double t = 0;

    // f(t_i, y_i)
    Vector<double> tmp(y.Size());
    auto tmpfunc = make_shared<ConstantFunction>(tmp);

    // set up equation
    auto yold = make_shared<ConstantFunction>(y); // y_i
    auto ynew = make_shared<IdentityFunction>(y.Size()); // y_{i+1}
    auto equ = ynew-yold - (dt/2) * (tmpfunc + rhs);

    for (int i = 0; i < steps; i++)
    { 
      // calculate f(t_i, y_i)
      rhs->Evaluate(y, tmp);
      tmp *= dt;
      tmp += y;

      // set up equation
      equ = ynew-yold - (dt/2) * (tmpfunc + rhs);
      // solve equation
      NewtonSolver (equ, y);
      yold->Set(y);

      t += dt;
      if (callback) callback(t, y);
    }
  }
  
  
  
  // Newmark and generalized alpha:
  // https://miaodi.github.io/finite%20element%20method/newmark-generalized/
  
  // Newmark method for  mass*d^2x/dt^2 = rhs
  void SolveODE_Newmark(double tend, int steps,
                        VectorView<double> x, VectorView<double> dx,
                        shared_ptr<NonlinearFunction> rhs,   
                        shared_ptr<NonlinearFunction> mass,  
                        std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    double gamma = 0.5;
    double beta = 0.25;

    Vector<> a(x.Size());
    Vector<> v(x.Size());

    auto xold = make_shared<ConstantFunction>(x);
    auto vold = make_shared<ConstantFunction>(dx);
    auto aold = make_shared<ConstantFunction>(x);

    rhs->Evaluate (xold->Get(), aold->Get());
    
    auto anew = make_shared<IdentityFunction>(a.Size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    auto equ = Compose(mass, anew) - Compose(rhs, xnew);

    double t = 0;
    for (int i = 0; i < steps; i++)            
      {
        NewtonSolver (equ, a);
        xnew -> Evaluate (a, x);
        vnew -> Evaluate (a, v);

        xold->Set(x);
        vold->Set(v);
        aold->Set(a);
        t += dt;
        if (callback) callback(t, x);
      }
    dx = v;
  }




  // Generalized alpha method for M d^2x/dt^2 = rhs
  void SolveODE_Alpha (double tend, int steps, double rhoinf,
                       VectorView<double> x, VectorView<double> dx, VectorView<double> ddx,
                       shared_ptr<NonlinearFunction> rhs,   
                       shared_ptr<NonlinearFunction> mass,  
                       std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    double alpham = (2*rhoinf-1)/(rhoinf+1);
    double alphaf = rhoinf/(rhoinf+1);
    double gamma = 0.5-alpham+alphaf;
    double beta = 0.25 * (1-alpham+alphaf)*(1-alpham+alphaf);

    Vector<> a(x.Size());
    Vector<> v(x.Size());

    auto xold = make_shared<ConstantFunction>(x);
    auto vold = make_shared<ConstantFunction>(dx);
    auto aold = make_shared<ConstantFunction>(ddx);
    // rhs->Evaluate (xold->Get(), aold->Get()); // solve with M ???
    
    auto anew = make_shared<IdentityFunction>(a.Size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    // auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - Compose(rhs, (1-alphaf)*xnew+alphaf*xold);
    auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - (1-alphaf)*Compose(rhs,xnew) - alphaf*Compose(rhs, xold);

    double t = 0;
    a = ddx;

    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, a);
        xnew -> Evaluate (a, x);
        vnew -> Evaluate (a, v);

        xold->Set(x);
        vold->Set(v);
        aold->Set(a);
        t += dt;
        if (callback) callback(t, x);
      }
    dx = v;
    ddx = a;
  }

  

}


#endif

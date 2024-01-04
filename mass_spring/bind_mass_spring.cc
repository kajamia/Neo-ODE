// #include <sstream>
// #include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "mass_spring.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<Mass<3>>);
PYBIND11_MAKE_OPAQUE(std::vector<Fix<3>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spring>);

PYBIND11_MODULE(mass_spring, m) {
    m.doc() = "mass-spring-system simulator";

    // not that elegant, but necessary
    py::class_<Vector<double>> (m, "Vector", py::buffer_protocol())
      .def(py::init<size_t>(),
           py::arg("size"), "create vector of given size")

      .def("__len__", &Vector<double>::Size,
           "return size of vector")
      
      .def("__setitem__", [](Vector<double> & self, int i, double v) {
        if (i < 0) i += self.Size();
        if (i < 0 || i >= self.Size()) throw py::index_error("vector index out of range");
        self(i) = v;
      })
      .def("__getitem__", [](Vector<double> & self, int i) { return self(i); })
      
      .def("__setitem__", [](Vector<double> & self, py::slice inds, double val)
      {
        size_t start, stop, step, n;
        if (!inds.compute(self.Size(), &start, &stop, &step, &n))
          throw py::error_already_set();
        self.Range(start, stop).Slice(0,step) = val;
      })
      
      .def("__add__", [](Vector<double> & self, Vector<double> & other)
      { return Vector<double> (self+other); })

      .def("__rmul__", [](Vector<double> & self, double scal)
      { return Vector<double> (scal*self); })
      
      .def("__str__", [](const Vector<double> & self)
      {
        std::stringstream str;
        str << self;
        return str.str();
      })

      .def_buffer([](Vector<double> & self) -> py::buffer_info {
        return py::buffer_info(
          self.Data(),
          sizeof(double),
          py::format_descriptor<double>::format(),
          1,
          {self.Size()},
          {sizeof(double) * self.Dist()}
        );})

      .def(py::pickle(
        [](Vector<double> & self) { // __getstate__
            /* return a tuple that fully encodes the state of the object */
          return py::make_tuple(self.Size(),
                                py::bytes((char*)(void*)&self(0), self.Size()*sizeof(double)));
        },
        [](py::tuple t) { // __setstate__
          if (t.size() != 2)
            throw std::runtime_error("should be a 2-tuple!");

          Vector<double> v(t[0].cast<size_t>());
          py::bytes mem = t[1].cast<py::bytes>();
          std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()), v.Size()*sizeof(double));
          return v;
        }))
    ;

    py::class_<Vec<3, double>> (m, "Vec3", py::buffer_protocol())
      .def("__setitem__", [](Vec<3, double> & self, int i, double v) {
        if (i < 0) i += self.Size();
        if (i < 0 || i >= self.Size()) throw py::index_error("vector index out of range");
        self(i) = v;
      })
      .def("__getitem__", [](Vec<3, double> & self, int i) { return self(i); })
      .def("__str__", [](const Vec<3, double> & self)
      {
        std::stringstream str;
        str << self;
        return str.str();
      });




    py::class_<Mass<2>> (m, "Mass2d")
      ;
      
    m.def("Mass", [](double m, std::array<double,2> p)
    {
      return Mass<2>{m, { p[0], p[1] }};
    });

    
    py::class_<Mass<3>> (m, "Mass3d")
      .def_property("mass",
                    [](Mass<3> & m) { return m.mass; },
                    [](Mass<3> & m, double mass) { m.mass = mass; })
      .def_property_readonly("pos",
                             [](Mass<3> & m) { return m.pos; });
    ;

    
    m.def("Mass", [](double m, std::array<double,3> p)
    {
      return Mass<3>{m, { p[0], p[1], p[2] }};
    });

    
    py::class_<Fix<3>> (m, "Fix3d")
      .def_property_readonly("pos",
                             [](Fix<3> & f) { return f.pos; });
    

    m.def("Fix", [](std::array<double,3> p)
    {
      return Fix<3>{ { p[0], p[1], p[2] } };
    });

    py::class_<Connector> (m, "Connector");

    py::class_<Spring> (m, "Spring")
      .def(py::init<double, double, std::array<Connector,2>>())
      .def_property_readonly("connections",
                             [](Spring & s) { return s.connections; })
      ;

    
    py::bind_vector<std::vector<Mass<3>>>(m, "Masses3d");
    py::bind_vector<std::vector<Fix<3>>>(m, "Fixes3d");
    py::bind_vector<std::vector<Spring>>(m, "Springs");        
    
    
    py::class_<MassSpringSystem<2>> (m, "MassSpringSystem2d")
      .def(py::init<>())
      .def("Add", [](MassSpringSystem<2> & mss, Mass<2> m) { return mss.AddMass(m); })
      ;
      
        
    py::class_<MassSpringSystem<3>> (m, "MassSpringSystem3d")
      .def(py::init<>())
      .def("__str__", [](MassSpringSystem<3> & mss) {
        stringstream sstr;
        sstr << mss;
        return sstr.str();
      })
      .def_property("gravity", [](MassSpringSystem<3> & mss) { return mss.Gravity(); },
                    [](MassSpringSystem<3> & mss, std::array<double,3> g) { mss.SetGravity(Vec<3>{g[0],g[1],g[2]}); })
      .def("Add", [](MassSpringSystem<3> & mss, Mass<3> m) { return mss.AddMass(m); })
      .def("Add", [](MassSpringSystem<3> & mss, Fix<3> f) { return mss.AddFix(f); })
      .def("Add", [](MassSpringSystem<3> & mss, Spring s) { return mss.AddSpring(s); })            
      .def_property_readonly("masses", [](MassSpringSystem<3> & mss) -> auto& { return mss.Masses(); })
      .def_property_readonly("fixes", [](MassSpringSystem<3> & mss) -> auto& { return mss.Fixes(); })
      .def_property_readonly("springs", [](MassSpringSystem<3> & mss) -> auto& { return mss.Springs(); })            
      .def("__getitem__", [](MassSpringSystem<3> mss, Connector & c) {
        if (c.type==Connector::FIX) return py::cast(mss.Fixes()[c.nr]);
        else return py::cast(mss.Masses()[c.nr]);
      })
      
      .def("GetState", [] (MassSpringSystem<3> & mss) {
        Vector<> x(3*mss.Masses().size());
        Vector<> dx(3*mss.Masses().size());
        Vector<> ddx(3*mss.Masses().size());
        mss.GetState (x, dx, ddx);
        return x;
      })
      ;
    

    m.def("Simulate", [](MassSpringSystem<3> & mss, double tend, size_t steps) {
      Vector<> x(3*mss.Masses().size());
      Vector<> dx(3*mss.Masses().size());
      Vector<> ddx(3*mss.Masses().size());
      mss.GetState (x, dx, ddx);
      
      auto mss_func = make_shared<MSS_Function<3>> (mss);
      auto mass = make_shared<IdentityFunction> (x.Size());      
      
      SolveODE_Alpha(tend, steps, 0.8, x, dx, ddx, mss_func, mass);

      mss.SetState (x, dx, ddx);
    });
      
    
}

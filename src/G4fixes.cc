#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/python.hpp>

#include <iostream>
#include <boost/python.hpp>

#include <G4VisExecutive.hh>

void G4fixes_init(int verbose=0)
{
   if (verbose)
      G4cout << "G4fixes_init - successfully loaded" << G4endl;
}

BOOST_PYTHON_MODULE(libG4fixes)
{
   using boost::python::class_;

   boost::python::def("G4fixes_init", &G4fixes_init);

   class_<G4VisExecutive, G4VisExecutive*, boost::noncopyable>
         ("G4VisExecutive",
          "wrapper around G4VisManager that automatically detects and "
          "registers visualization drivers supported by the G4 library build")
      .def("Initialize", &G4VisExecutive::Initialize)
   ;
}

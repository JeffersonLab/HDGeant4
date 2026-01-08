//
// g4py.cc -- python bindings for HDGeant4 user classes
//            using the Boost.Python C++ interface.
//

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/python.hpp>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "DANA/DEvent.h"
#include <DANA/DApplication.h>
#include <GlueXUserOptions.hh>
#include <GlueXDetectorConstruction.hh>
#include <GlueXPrimaryGeneratorAction.hh>
#include <GlueXRunAction.hh>
#include <GlueXEventAction.hh>
#include <GlueXSteppingAction.hh>
#include <GlueXSteppingVerbose.hh>
#include <GlueXPhysicsList.hh>

#include <G4SystemOfUnits.hh>
#include <G4OpenGLViewer.hh>
#include <G4TransportationManager.hh>

// We must wrap an abstract C++ type so python 
// knows how to override pure virtual methods.

struct CB_G4VUserParallelWorld : 
       G4VUserParallelWorld, boost::python::wrapper<G4VUserParallelWorld>
{
   CB_G4VUserParallelWorld(const G4String name)
    : G4VUserParallelWorld(name) {}

   void Construct() {
     this->get_override("Construct")(); 
   }
};

// We must wrap the DApplication class because its only constructor
// has a signature that is not supported by Boost.Python.

struct CB_DApplication : public DApplication {
   CB_DApplication() : DApplication(0,0) {}
};

// Wrap overloaded GlueXUserOptions::Find methods using functions
// with unique names. This shows how to pass stl containers back
// and forth between C++ and python using Boost.Python converters.

typedef std::map<int, std::string> GlueXUserOptions_string_map;
int (GlueXUserOptions::*GlueXUserOptions_Find_string)
    (const char *name, GlueXUserOptions_string_map &value) const
    = &GlueXUserOptions::Find;
typedef std::map<int, double> GlueXUserOptions_double_map;
int (GlueXUserOptions::*GlueXUserOptions_Find_double)
    (const char *name, GlueXUserOptions_double_map &value) const
    = &GlueXUserOptions::Find;
typedef std::map<int, int> GlueXUserOptions_int_map;
int (GlueXUserOptions::*GlueXUserOptions_Find_int)
    (const char *name, GlueXUserOptions_int_map &value) const
    = &GlueXUserOptions::Find;

// Add a python utility function to fake a user "picking" click
// on a particular 3D location in the geometry, to enable probing
// the contents of the geometry and magnetic fields from the
// python command prompt. Of course, the same code can be used
// from a C++ user function if so desired.

void pickPoint3D(G4double x_cm, G4double y_cm, G4double z_cm)
{
   // Act as if the user had "picked" on point x,y,z (specified in cm)
   // from within a graphical view in the OGLSX visualization system.

   G4OpenGLViewerPickMap map;
   G4ThreeVector v(x_cm * cm, y_cm * cm, z_cm * cm);
   map.setPickCoordinates3D(v);
   G4cout << map.print();
}

// More python utility functions to access properties not reachable
// through the standard g4py exposed interfaces.

G4Navigator *GetNavigator(int world)
{
   // Look up the navigator by numerical index, where world=0 is the
   // navigator on the default layer, and parallel layers are world>0
   // returning zero if the parallel world does not exist.
 
   G4TransportationManager *tmanager = G4TransportationManager::
                                       GetTransportationManager();
   if (world < (int)tmanager->GetNoWorlds())
      return tmanager->GetNavigator(tmanager->GetWorldsIterator()[world]);
   else
      return 0;
}

bool GetGlobalExitNormal(int world, G4ThreeVector *point,
                                    G4ThreeVector *norm)
{
   // Forwards the request to the GetGlobalExitNormal method of the
   // navigator corresponding to world, and returns the unit normal
   // vector if the method returns valid=true, otherwise returns null.
 
   G4Navigator *nav = GetNavigator(world);
   if (nav) {
      bool valid;
      *norm = nav->GetGlobalExitNormal(*point, &valid);
      return valid;
   }
   return 0;
}

// Create a python module containing all of the G4 user classes
// that are needed to run the HDGeant simulation from python.
// Here it is named libhdgeant4 (happens to also be the name of
// the shared library) but when it is loaded from python, the
// name of the module is HDGeant4.

BOOST_PYTHON_MODULE(libhdgeant4)
{
   using boost::python::class_;
   using boost::python::enum_;
   using boost::python::def;

   class_<CB_DApplication, CB_DApplication*>
         ("DApplication", 
          "singleton class holding configuration data for current run, "
          "part of the standard jana framework",
          boost::python::init<>())
      .def("Init", &JApplication::Initialize)
   ;

   enum_<jerror_t>("jerror_t")
      .value("NOERROR", NOERROR)
      .value("UNKNOWN_ERROR", UNKNOWN_ERROR)
      .value("MAX_EVENT_PROCESSORS_EXCEEDED", MAX_EVENT_PROCESSORS_EXCEEDED)
      .value("ERROR_OPENING_EVENT_SOURCE", ERROR_OPENING_EVENT_SOURCE)
      .value("ERROR_CLOSING_EVENT_SOURCE", ERROR_CLOSING_EVENT_SOURCE)
      .value("NO_MORE_EVENTS_IN_SOURCE", NO_MORE_EVENTS_IN_SOURCE)
      .value("NO_MORE_EVENT_SOURCES", NO_MORE_EVENT_SOURCES)
      .value("EVENT_NOT_IN_MEMORY", EVENT_NOT_IN_MEMORY)
      .value("EVENT_SOURCE_NOT_OPEN", EVENT_SOURCE_NOT_OPEN)
      .value("OBJECT_NOT_AVAILABLE", OBJECT_NOT_AVAILABLE)
      .value("DEVENT_OBJECT_DOES_NOT_EXIST", DEVENT_OBJECT_DOES_NOT_EXIST)
      .value("MEMORY_ALLOCATION_ERROR", MEMORY_ALLOCATION_ERROR)
      .value("RESOURCE_UNAVAILABLE", RESOURCE_UNAVAILABLE)
      .value("VALUE_OUT_OF_RANGE", VALUE_OUT_OF_RANGE)
      .value("INFINITE_RECURSION", INFINITE_RECURSION)
      .value("UNRECOVERABLE_ERROR",  UNRECOVERABLE_ERROR)
      .value("FILTER_EVENT_OUT", FILTER_EVENT_OUT)
   ;

   class_<GlueXUserOptions_string_map, GlueXUserOptions_string_map*>
         ("GlueXUserOptions_string_map",
          "C++ std::map<int, std::string> as a python indexable object")
      .def(boost::python::map_indexing_suite<GlueXUserOptions_string_map, true>())
   ;
   class_<GlueXUserOptions_double_map, GlueXUserOptions_double_map*>
         ("GlueXUserOptions_double_map",
          "C++ std::map<int, double> as a python indexable object")
      .def(boost::python::map_indexing_suite<GlueXUserOptions_double_map, true>())
   ;
   class_<GlueXUserOptions_int_map, GlueXUserOptions_int_map*>
         ("GlueXUserOptions_int_map",
          "C++ std::map<int, int> as a python indexable object")
      .def(boost::python::map_indexing_suite<GlueXUserOptions_int_map, true>())
   ;

   class_<GlueXUserOptions, GlueXUserOptions*>
         ("GlueXUserOptions",
          "general class for fetching and holding user-defined run parameters")
      .def("ReadControl_in", &GlueXUserOptions::ReadControl_in)
      .def("GetInstance", &GlueXUserOptions::GetInstance,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .staticmethod("GetInstance")
      .def("Find", GlueXUserOptions_Find_string)
      .def("Find", GlueXUserOptions_Find_double)
      .def("Find", GlueXUserOptions_Find_int)
   ;

   class_<GlueXDetectorConstruction, GlueXDetectorConstruction*,
          boost::python::bases<G4VUserDetectorConstruction> >
          ("GlueXDetectorConstruction",
           "sets up the GlueX detector geometry and fields for simulation")
      .def(boost::python::init<G4String>())
      .def("Construct", &GlueXDetectorConstruction::Construct,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetUniformField", &GlueXDetectorConstruction::SetUniformField)
      .def("SetMaxStep", &GlueXDetectorConstruction::SetMaxStep)
      .def("GetParallelWorldCount", &GlueXDetectorConstruction::GetParallelWorldCount)
      .def("GetParallelWorldName", &GlueXDetectorConstruction::GetParallelWorldName)
      .def("GetParallelWorldVolume", &GlueXDetectorConstruction::GetParallelWorldVolume,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("RegisterParallelWorld", &G4VUserDetectorConstruction::RegisterParallelWorld)
      .def("ConstructParallelGeometries", &G4VUserDetectorConstruction::ConstructParallelGeometries)
   ;

   class_<CB_G4VUserParallelWorld, CB_G4VUserParallelWorld*, boost::noncopyable>
         ("G4VUserParallelWorld",
          "base class for user implementations of parallel geometries",
          boost::python::init<G4String>())
      .def("GetName", &G4VUserParallelWorld::GetName)
      .def("Construct", boost::python::pure_virtual(&G4VUserParallelWorld::Construct))
   ;

   class_<GlueXParallelWorld, GlueXParallelWorld*,
          boost::python::bases<G4VUserParallelWorld> >
          ("GlueXParallelWorld", 
           "GlueX implementation of parallel mass geometry",
           boost::python::init<G4String, G4LogicalVolume*>())
      .def("Construct", &GlueXParallelWorld::Construct)
   ;

   class_<GlueXPhysicsList, GlueXPhysicsList*,
          boost::python::bases<G4VUserPhysicsList> >
         ("GlueXPhysicsList",
          "a custom physics list for GlueX simulation",
           boost::python::init<GlueXDetectorConstruction*>())
   ;

   class_<GlueXPrimaryGeneratorAction, GlueXPrimaryGeneratorAction*,
          boost::python::bases<G4VUserPrimaryGeneratorAction> >
         ("GlueXPrimaryGeneratorAction",
          "the generator of interaction events to be simulated",
           boost::python::init<>())
   ;

   class_<GlueXEventAction, GlueXEventAction*,
          boost::python::bases<G4UserEventAction> >
         ("GlueXEventAction",
          "encapsulates actions to take at start and end of each event")
   ;

   class_<GlueXRunAction, GlueXRunAction*,
          boost::python::bases<G4UserRunAction> >
         ("GlueXRunAction",
          "encapsulates actions to take at start and end of each run",
           boost::python::init<GlueXPhysicsList*>())
   ;

   class_<GlueXSteppingAction, GlueXSteppingAction*,
          boost::python::bases<G4UserSteppingAction>, boost::noncopyable>
         ("GlueXSteppingAction",
          "encapsulates actions to take at start and end of each step",
           boost::python::init<>())
   ;

   class_<G4SteppingVerbose, G4SteppingVerbose*>
         ("G4SteppingVerbose",
          "Geant4 standard stepping action class with verbose output")
   ;

   class_<GlueXSteppingVerbose, GlueXSteppingVerbose*,
          boost::python::bases<G4SteppingVerbose> >
         ("GlueXSteppingVerbose",
          "encapsulates actions to take at start and end of each step (verbose)")
   ;

   def("pickPoint3D", pickPoint3D);
   def("GetNavigator", GetNavigator,
       boost::python::return_value_policy<boost::python::reference_existing_object>());
   def("GetGlobalExitNormal", GetGlobalExitNormal);
}

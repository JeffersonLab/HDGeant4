//
// class implementation for
//   * GlueXUniformMagField
//   * GlueXMappedMagField
//   * GlueXComputedMagField
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include <iostream>
#include <fstream>

#include "GlueXMagneticField.hh"
#include "GlueXUserOptions.hh"
#include "GlueXDetectorConstruction.hh"
#include "HddmOutput.hh"

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>

#include <HDGEOMETRY/DMagneticFieldMapFineMesh.h>
#include <HDGEOMETRY/DMagneticFieldMapNoField.h>
#include <HDGEOMETRY/DMagneticFieldMapConst.h>
#include <HDGEOMETRY/DMagneticFieldMapPS2DMap.h>
#include <HDGEOMETRY/DMagneticFieldMapPSConst.h>

// Implementation code for class GlueXUniformMagField

GlueXUniformMagField::GlueXUniformMagField(G4ThreeVector B, G4double unit,
                                           const G4AffineTransform &xform)
  : G4UniformMagField(G4ThreeVector()),
    fXform(xform)
{
   // uniform magnetic field constructor, constant field value B required
   // as input. Factor unit converts B to standard G4 magnetic field units.
   // Transform xform performs immediate active rotation on B, and is not
   // saved. Translation in xform is ignored.
 
   xform.ApplyAxisTransform(B);
   SetMagField(B, unit);
}

GlueXUniformMagField::~GlueXUniformMagField()
{ }

GlueXUniformMagField::GlueXUniformMagField(const GlueXUniformMagField &src)
  : G4UniformMagField(src.GetConstantFieldValue()),
    fXform(src.fXform)
{ }

GlueXUniformMagField& 
GlueXUniformMagField::operator=(const GlueXUniformMagField &src)
{
   // assignment operator
 
   fXform = src.fXform;
   SetMagField(src.GetConstantFieldValue(), 1);
   return *this;
}
      
void GlueXUniformMagField::SetMagFieldTransform(const G4AffineTransform &xform)
{
   // change the rotation applied to the field vector

   G4ThreeVector B = GetConstantFieldValue();
   fXform.Inverse().ApplyAxisTransform(B);
   xform.ApplyAxisTransform(B);
   SetFieldValue(B);
   fXform = xform;
}

const G4AffineTransform& 
GlueXUniformMagField::GetMagFieldTransform() const
{
   // return the rotation applied to the field vector

   return fXform;
}

void GlueXUniformMagField::SetMagField(G4ThreeVector B, G4double unit)
{
   // supplementary magnetic field setter, allows user to specify factor
   // unit that converts input field B into G4 standard field units.
   // Standard magnetic field setter method SetFieldValue works as well
   // for cases when B is already in G4 units.

   fXform.ApplyAxisTransform(B);
   SetFieldValue(B *= unit);
}

G4ThreeVector GlueXUniformMagField::GetMagField(G4double unit) const
{
   // supplementary magnetic field getter, allows user to specify factor
   // unit that is used to convert internal G4 field values into any
   // desired units, eg. tesla.

   G4ThreeVector B = GetConstantFieldValue();
   if (B[2] == 0)
      B[2] = 1e-99; // avoid divide-by-zero by excluding |B|=0
   return B *= unit;
}


// Implementation code for class GlueXMappedMagField

GlueXMappedMagField::GlueXMappedMagField(G4double Bmax, G4double unit,
                                         const G4AffineTransform &xform)
 : fUnit(unit),
   fBmax(Bmax),
   fXform(xform),
   fGridtype(0)
{ 
   // mapped magnetic field constructor, maximum field value Bmax required
   // as input. Factor unit converts B (Bmax and field components that will
   // be read from the map) to standard G4 magnetic field units. Transform
   // xform provides 2 functions: its inverse is used to transform the user
   // coordinates used to access the map into local field map coordinates,
   // and then once the map components are extracted from the map, an active
   // rotation from xform is performed on B before it is returned to the user.
   // Note that one or more calls to either method AddCartesianGrid() or
   // AddCylindricalGrid() followed by ReadMapFile() must be issued before
   // the object is ready for calls to GetFieldValue() or GetMagField().
 
   fXfinv = xform.Inverse();
}

GlueXMappedMagField::~GlueXMappedMagField()
{ }

GlueXMappedMagField::GlueXMappedMagField(const GlueXMappedMagField &src)
 : G4MagneticField(src)
{
   // copy constructor

   *this = src;
}

GlueXMappedMagField& 
GlueXMappedMagField::operator=(const GlueXMappedMagField &src)
{
   // assignment operator, map copies are deep so use sparingly
 
   fUnit = src.fUnit;
   fBmax = src.fBmax;
   fXform = src.fXform;
   fXfinv = src.fXfinv;
   fGridtype = src.fGridtype;
   for (int dim=0; dim < 3; ++dim)
   {
      fDim[dim] = src.fDim[dim];
      fOrder[dim] = src.fOrder[dim];
   }
   fGrid = src.fGrid;
   fFieldMap = src.fFieldMap;
   return *this;
}
      
void GlueXMappedMagField::SetMagFieldTransform(const G4AffineTransform &xform)
{
   // provides a means to modify the placement of the mapped field in
   // the geometry without having to reinitialize the map.
 
   fXform = xform;
   fXfinv = xform.Inverse();
}

const G4AffineTransform& GlueXMappedMagField::GetMagFieldTransform() const
{
   // returns a copy of the current map placement transform
 
   return fXform;
}
 
int GlueXMappedMagField::AddCartesianGrid(const int axsamples[4],
                                          const int axorder[4],
                                          const int axsense[4],
                                          const G4double axlower[4],
                                          const G4double axupper[4])
{
   // specifies that the field map is specified on the nodes of a Cartesian
   // grid, what the dimensions of the grid are (axsamples), how the 3D
   // array of grid points is ordered into a linear list (axorder), whether
   // or not to reverse the sign of the field coordinates when they are read
   // from the map (axsense), and the lower (axlower) and upper (axupper)
   // coordinates of the bounding box of the grid. Multiple calls are allowed
   // to AddCartesianGrid, to enable the same map to be replicated several
   // times (with different axsense inversions in each case, as needed) so
   // that symmetries in the field geometry can be exploited to minimize the
   // total size of the map and eliminate redundancies.
   //     Input arrays all have the same component structure: x=1, y=2, z=3
   // and the 0'th component is unused. Units for axlower and axupper are cm.
   // Values in axorder[1..3] should be in the range 1..3. Values in axsense
   // should be either +1 or -1. Repeated calls to AddCartesianGrid must all
   // agree in the contents of axsamples and axorder, but can differ in the
   // remaining arguments. This is assumed but not checked.
 
   fGridtype = 1;
   struct field_map_grid_t grid;
   for (int dim=0; dim < 3; ++dim)
   {
      fDim[dim] = axsamples[dim + 1];
      fOrder[dim] = axorder[dim + 1];
      grid.sense.push_back(axsense[dim + 1]);
      grid.lower.push_back(axlower[dim + 1]);
      grid.upper.push_back(axupper[dim + 1]);
   }
   fGrid.push_back(grid);
   return 1;
}

int GlueXMappedMagField::AddCylindricalGrid(const int axsamples[4],
                                            const int axorder[4],
                                            const int axsense[4],
                                            const G4double axlower[4],
                                            const G4double axupper[4])
{
   // specifies that the field map is specified on the nodes of a cylindrical
   // grid, what the dimensions of the grid are (axsamples), how the 3D
   // array of grid points is ordered into a linear list (axorder), whether
   // or not to reverse the sign of the field coordinates when they are read
   // from the map (axsense), and the lower (axlower) and upper (axupper)
   // coordinates of the bounding box of the grid. Multiple calls are allowed
   // to AddCylindricalGrid, to enable the same map to be replicated several
   // times (with different axsense inversions in each case, as needed) so
   // that symmetries in the field geometry can be exploited to minimize the
   // total size of the map and eliminate redundancies.
   //     Input arrays all have the same component structure: rho=1, phi=2, z=3
   // and the 0'th component is unused. Units for axlower and axupper are cm
   // for rho and z, and radians for phi. Values in axorder[1..3] should be in
   // the range 1..3. Values in axsense should be either +1 or -1. Repeated
   // calls to AddCylindricalGrid must all agree in the contents of axsamples
   // and axorder, but can differ in the remaining arguments. This is assumed
   // but not checked.
 
   fGridtype = 2;
   struct field_map_grid_t grid;
   for (int dim=0; dim < 3; ++dim)
   {
      fDim[dim] = axsamples[dim + 1];
      fOrder[dim] = axorder[dim + 1];
      grid.sense.push_back(axsense[dim + 1]);
      grid.lower.push_back(axlower[dim + 1]);
      grid.upper.push_back(axupper[dim + 1]);
   }
   fGrid.push_back(grid);
   return 1;
}

int GlueXMappedMagField::ReadMapFile(const char *mapS)
{
   // reads field values from an input text file and stores them in an 
   // internal table for interpolation according to a preconceived spatial
   // grid. The grid coordinates are not contained in the file, and are
   // specified separately through calls to methods AddCartesianGrid() or
   // AddCylindricalGrid(). A map is assumed to be either Cartesian or
   // cylindrical, and invoking both methods on the same map will produce
   // nonsense. The input file should contain one triplet of floating point
   // numbers on each line separated by spaces, each representing the
   // magnetic field components at a given grid site ordered as Bx By Bz
   // for a Cartesian grid or Brho Bphi Bz for a cylindrical grid. It is
   // assumed (but not checked) that the number of lines in the file equals
   // the product of the three dimensions of the grid. Note that at least
   // one call to either AddCartesianGrid() or AddCylindricalGrid() must
   // have occurred on this object prior to the invocation of ReadMapFile.

   std::ifstream mapfile(mapS);
   if (! mapfile.good())
   {
      G4cerr << "GlueXMappedMagField::ReadMapFile error - "
             << "open failed on input map file " << mapS
             << G4endl;
      return 0;
   }

   while (mapfile.good())
   {
      union field_map_entry_t entry;
      mapfile >> entry.cart.x >> entry.cart.y >> entry.cart.z;
      if (mapfile.good())
      {
         fFieldMap.push_back(entry);
      }
      else {
         break;
      }
   }
   return 1;
}

G4ThreeVector GlueXMappedMagField::GetMagField(const G4double point[4],
                                               G4double unit) const
{
   // interpolates the magnetic field map at point, and returns the result
   // scaled to the requested unit (eg. tesla). Position point[0..2] gives
   // the Cartesian coordinates at which the map is to be evaluated. The
   // time stored in point[3] is presently ignored. The position is transformed
   // into local map coordinates, and a 3D linear interpolation based on
   // central estimates for the local field gradient is used. Finally, the
   // field is rotated back into user coordinates and returned in the form
   // of a Cartesian 3-vector.
 
   G4ThreeVector p(point[0], point[1], point[2]);
   fXfinv.ApplyPointTransform(p);
   double u[3];
   if (fGridtype == 1)
   {
      u[0] = p[0] / cm;
      u[1] = p[1] / cm;
      u[2] = p[2] / cm;
   }
   else if (fGridtype == 2)
   {
      u[0] = sqrt(p[0] * p[0] + p[1] * p[1]) / cm;
      u[1] = atan2(p[1], p[0]);
      u[2] = p[2] / cm;
   }
   else
   {
      return G4ThreeVector(0, 0, 1e-99);
   }

   double B[3] = {0,0,0};
   std::vector<struct field_map_grid_t>::const_iterator iter;
   for (iter = fGrid.begin(); iter != fGrid.end(); ++iter)
   {
      double unorm[3];
      unorm[0] = (u[0] - iter->lower[0]) / (iter->upper[0] - iter->lower[0]);
      unorm[1] = (u[1] - iter->lower[1]) / (iter->upper[1] - iter->lower[1]);
      unorm[2] = (u[2] - iter->lower[2]) / (iter->upper[2] - iter->lower[2]);
      if (unorm[0] < 0 || unorm[0] > 1 ||
          unorm[1] < 0 || unorm[1] > 1 ||
          unorm[2] < 0 || unorm[2] > 1)
      {
         break;
      }
      double ur[3], dur[3];
      int ui[3], ui0[3], ui1[3];
      for (int dim = 0; dim < 3; ++dim)
      {
         ur[dim] = unorm[dim] * (fDim[dim] - 1);
         ui[dim] = ceil(ur[dim] - 0.5);
         ui0[dim] = (ui[dim] > 0)? ui[dim] - 1 : ui[dim];
         ui1[dim] = (ui[dim] < fDim[dim] - 1)? ui[dim] + 1 : ui[dim];
         dur[dim] = (ur[dim] - ui[dim]) / (ui1[dim] - ui0[dim] + 1e-20);
      }
      double grad[3][3];
      double mapfield0[3], mapfield1[3];
      lookup_field(ui0[0], ui[1], ui[2], mapfield0);
      lookup_field(ui1[0], ui[1], ui[2], mapfield1);
      grad[0][0] = mapfield1[0] - mapfield0[0];
      grad[1][0] = mapfield1[1] - mapfield0[1];
      grad[2][0] = mapfield1[2] - mapfield0[2];
      lookup_field(ui[0], ui0[1], ui[2], mapfield0);
      lookup_field(ui[0], ui1[1], ui[2], mapfield1);
      grad[0][1] = mapfield1[0] - mapfield0[0];
      grad[1][1] = mapfield1[1] - mapfield0[1];
      grad[2][1] = mapfield1[2] - mapfield0[2];
      lookup_field(ui[0], ui[1], ui0[2], mapfield0);
      lookup_field(ui[0], ui[1], ui1[2], mapfield1);
      grad[0][2] = mapfield1[0] - mapfield0[0];
      grad[1][2] = mapfield1[1] - mapfield0[1];
      grad[2][2] = mapfield1[2] - mapfield0[2];

      lookup_field(ui[0], ui[1], ui[2], B);
      B[0] += grad[0][0] * dur[0] + grad[0][1] * dur[1] + grad[0][2] * dur[2];
      B[1] += grad[1][0] * dur[0] + grad[1][1] * dur[1] + grad[1][2] * dur[2];
      B[2] += grad[2][0] * dur[0] + grad[2][1] * dur[1] + grad[2][2] * dur[2];
      B[0] *= iter->sense[0];
      B[1] *= iter->sense[1];
      B[2] *= iter->sense[2];
   }

   G4ThreeVector Bvec;
   if (fGridtype == 1)
   {
      Bvec.set(B[0], B[1], B[2]);
   }
   else if (fGridtype == 2)
   {
      Bvec.setRhoPhiZ(B[0], B[1], B[2]);
   }
   fXform.ApplyAxisTransform(Bvec);
   if (Bvec[2] == 0)
      Bvec[2] = 1e-99; // avoid divide-by-zero by excluding |B|=0
   return Bvec *= fUnit / unit;
}

void GlueXMappedMagField::GetFieldValue(const G4double point[4],
                                        G4double Bfield[3] ) const
{
   // implements the standard G4MagneticField::GetFieldValue() interface

   G4ThreeVector Bvec = GetMagField(point, 1);
   Bfield[0] = Bvec[0];
   Bfield[1] = Bvec[1];
   Bfield[2] = Bvec[2];
}

int GlueXMappedMagField::lookup_field(int i1, int i2, int i3,
                                      G4double mapvalue[3]) const
{
   // private helper method to perform lookup of the field at a single
   // site in the field map. Many such accesses are performed in a 
   // single interpolation of the field, so it should be reasonably
   // efficient.

   int ivec[3] = {i1, i2, i3};
   int iord[3] = {ivec[fOrder[0]], ivec[fOrder[1]], ivec[fOrder[2]]};
   int nord[3] = {fDim[fOrder[0]], fDim[fOrder[1]], fDim[fOrder[2]]};
   int index = nord[1] * (nord[0] * iord[0] + iord[1]) + iord[2];
   if (index < (int)fFieldMap.size())
   {
      mapvalue[0] = fFieldMap[index].cart.x;
      mapvalue[1] = fFieldMap[index].cart.y;
      mapvalue[2] = fFieldMap[index].cart.z;
      return 1;
   }
   return 0;
}


// Implementation code for class GlueXComputedMagField
// Here a "computed" field covers any case where the map from position to
// field value is returned by an external function, whether that function
// implements the field as a formula or looks it up in a database. This is
// the choice to use for any case where the field is not uniform, and the
// map is not stored in the specific format specified for HDDS field maps.

GlueXComputedMagField::GlueXComputedMagField(G4double Bmax, G4double unit,
                                             const G4AffineTransform &xform)
 : fUnit(unit),
   fBmax(Bmax),
   fXform(xform),
   fJanaFieldMap(0),
   fJanaFieldMapPS(0)
{
   // computed magnetic field constructor, maximum field value Bmax required
   // as input. Factor unit converts B (Bmax and field components that will
   // be read from the map) to standard G4 magnetic field units. Transform
   // xform provides 2 functions: its inverse is used to transform the user
   // coordinates used to access the map into local field map coordinates,
   // and then once the map components are extracted from the map, an active
   // rotation from xform is performed on B before it is returned to the user.
   // Note that method SetFunction() must be issued before the object is ready
   // for calls to GetFieldValue() or GetMagField().
}

GlueXComputedMagField::~GlueXComputedMagField()
{
   if (fJanaFieldMap)
      delete fJanaFieldMap;
   if (fJanaFieldMapPS)
      delete fJanaFieldMapPS;
}
      
GlueXComputedMagField::GlueXComputedMagField(const GlueXComputedMagField &src)
 : G4MagneticField(src)
{
   // copy constructor

   *this = src;
}

GlueXComputedMagField&
GlueXComputedMagField::operator=(const GlueXComputedMagField &src)
{
   // assignment operator

   fBmax = src.fBmax;
   fUnit = src.fUnit;
   fXform = src.fXform;
   fFunction = src.fFunction;
   fJanaFieldMap = 0;
   DMagneticFieldMapFineMesh *mapFineMesh = 0;
   DMagneticFieldMapNoField *mapNoField = 0;
   DMagneticFieldMapConst *mapConstField = 0;
   if (src.fJanaFieldMap) {
      mapFineMesh = dynamic_cast<DMagneticFieldMapFineMesh*>(src.fJanaFieldMap);
      mapNoField = dynamic_cast<DMagneticFieldMapNoField*>(src.fJanaFieldMap);
      mapConstField = dynamic_cast<DMagneticFieldMapConst*>(src.fJanaFieldMap);
      if (mapFineMesh)
         fJanaFieldMap = new DMagneticFieldMapFineMesh(*mapFineMesh);
      else if (mapNoField)
         fJanaFieldMap = new DMagneticFieldMapNoField(*mapNoField);
      else if (mapConstField)
         fJanaFieldMap = new DMagneticFieldMapConst(*mapConstField);
   }
   fJanaFieldMapPS = 0;
   DMagneticFieldMapPS2DMap *mapPS2D = 0;
   DMagneticFieldMapPSConst *mapPSConst = 0;
   if (src.fJanaFieldMapPS) {
      mapPS2D = dynamic_cast<DMagneticFieldMapPS2DMap*>(src.fJanaFieldMapPS);
      mapPSConst = dynamic_cast<DMagneticFieldMapPSConst*>(src.fJanaFieldMapPS);
      if (mapPS2D)
	     fJanaFieldMapPS = new DMagneticFieldMapPS2DMap(*mapPS2D);
      else if (mapPSConst)
         fJanaFieldMapPS = new DMagneticFieldMapPSConst(*mapPSConst);
   }
   return *this;
}

void GlueXComputedMagField::SetFunction(std::string function)
{
   // Pull a DMagneticField object from the JANA framework
   // and store a pointer to it in this object. Note that
   // GlueXComputedMagField never actually owns the object.
   // The function argument supports a few special values:
   //    1) gufld_db(r,B) - refers to the GlueX solenoid field
   //    2) gufld_ps(r,B) - refers to the GlueX pair spectrometer field
   // Lookup of the appropriate magnetic field map in the ccdb
   // is done using a database key whose value is fetched from the
   // global GlueXUserOptions object. Hard-coded defaults are also
   // provided below, but normally the user should specify the 
   // correct field map key as an input option, eg. through the
   // control.in file.

   extern jana::JApplication *japp;
   if (japp == 0) {
      G4cerr << "Error in GlueXComputedMagField::SetFunction - "
             << "jana global DApplication object not set, "
             << "cannot continue." << G4endl;
      exit(-1);
   }
   jana::JParameterManager *jpars = japp->GetJParameterManager();
   if (jpars == 0) {
      G4cerr << "Error in GlueXComputedMagField::SetFunction - "
             << "jana global DApplication object returns null "
             << "JParameterManager object, cannot continue." << G4endl;
      exit(-1);
   }
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXComputedMagField::SetFunction - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   int run_number = HddmOutput::getRunNo();
   jana::JCalibration *jcalib = japp->GetJCalibration(run_number);

   if (function == "gufld_db(r,B)") {

      // load the solenoid field map
  
      std::map<int, std::string> type_opts;
      std::map<int, std::string> map_opts;
      GlueXUserOptions::GetInstance()->Find("BFIELDTYPE", type_opts);
      GlueXUserOptions::GetInstance()->Find("BFIELDMAP", map_opts);
      if (type_opts.find(1) == type_opts.end() || type_opts[1] == "CalibDB") {

         // if the magnetic field is specified in control.in
         // then use that value instead of the CCDB values
 
         if (map_opts.find(1) != map_opts.end()) {
            fJanaFieldMap = new DMagneticFieldMapFineMesh(japp, run_number,
                                                          map_opts[1]);
         }

         // otherwise see if we can load the name of the 
         // magnetic field map to use from the ccdb itself

         else {
            std::map<std::string, std::string> map_name;
            const char *map_key = "/Magnets/Solenoid/solenoid_map";
            if (jcalib->GetCalib(map_key, map_name)) {
               G4cerr << "Error in GlueXComputedMagField::SetFunction - "
                      << "failed to figure out which magnetic field map "
                      << "to use for the solenoid in this simulation, "
                      << "cannot continue." << G4endl;
               exit(-1);
            }
            else if (map_name.find("map_name") != map_name.end()) {
	       if (map_name["map_name"] == "NoField")
	          fJanaFieldMap = new DMagneticFieldMapNoField(japp);
	       else
	          fJanaFieldMap = new DMagneticFieldMapFineMesh(japp,
                                      run_number, map_name["map_name"]);
	    }
            else {
               G4cerr << "Error in GlueXComputedMagField::SetFunction - "
                      << "no solenoid magnetic field map specified for "
                      << "this run in either the simulation options "
                      << "or in ccdb, cannot continue." << G4endl;
	       exit(-1);
	    }
         }
      }
      else if (type_opts[1] =="NoField") {
         fJanaFieldMap = new DMagneticFieldMapNoField(japp);
      }
      else if (type_opts[1] == "Const") {
         const GlueXDetectorConstruction *geom =
               GlueXDetectorConstruction::GetInstance();
         double field_T = (geom)? geom->GetUniformField(tesla) : 1.9;
         if (type_opts.find(2) != type_opts.end()) {
            field_T = std::atof(type_opts[2].c_str());
         }
         fJanaFieldMap = new DMagneticFieldMapConst(0.0, 0.0, field_T);
      }
      else {
         G4cerr << "Error in GlueXComputedMagField::SetFunction - "
                << "unknown DMagneticFieldMap type " << type_opts[1]
                << ", cannot continue." << G4endl;
         exit(-1);
      }
   }

   else if (function == "gufld_ps(r,B)") {

      // load the pair spectrometer field map

      std::map<int, std::string> type_opts;
      std::map<int, std::string> map_opts;
      GlueXUserOptions::GetInstance()->Find("PSBFIELDTYPE", type_opts);
      GlueXUserOptions::GetInstance()->Find("PSBFIELDMAP", map_opts);
      if (type_opts.find(1) == type_opts.end() || type_opts[1] == "CalibDB") {
 
         // if the magnetic field is specified in control.in
         // then use that value instead of the CCDB values
 
         if (map_opts.find(1) != map_opts.end()) {
            fJanaFieldMapPS = new DMagneticFieldMapPS2DMap(japp, run_number,
                                                           map_opts[1]);
         }
         else {

           // see if we can load the name of the magnetic 
           // field map to use from the calib DB

            map<string,string> map_name;
            const char *map_key = "/Magnets/PairSpectrometer/ps_magnet_map";
            if (jcalib->GetCalib(map_key, map_name)) {
               G4cerr << "Error in GlueXComputedMagField::SetFunction - "
                      << "failed to figure out which magnetic field map "
                      << "to use for the pair spectrometer, "
                      << "cannot continue." << G4endl;
               exit(-1);
            }
            else if (map_name.find("map_name") != map_name.end()) {
	       fJanaFieldMapPS = new DMagneticFieldMapPS2DMap(japp,
                                     run_number, map_name["map_name"]);
	    }
            else {
               G4cerr << "Error in GlueXComputedMagField::SetFunction - "
                      << "no pair spectrometer magnetic field map "
                      << "specified for this run either in the options "
                      << "or in ccdb, cannot continue." << G4endl;
	       exit(-1);
	    }
         }
      }
      else if (type_opts[1] == "Const") {
         double field_T = 1.64;
         if (type_opts.find(2) != type_opts.end()) {
            field_T = std::atof(type_opts[2].c_str());
         }
         fJanaFieldMapPS = new DMagneticFieldMapPSConst(0.0, field_T, 0.0);
      }
      else {
         G4cerr << "Error in GlueXComputedMagField::SetFunction - "
                << "unknown DMagneticFieldMapPS type " << type_opts[1]
                << ", cannot continue." << G4endl;
         exit(-1);
      }
   }

   fFunction = function;
}

std::string GlueXComputedMagField::GetFunction() const
{
   // returns the string specified in a prior call to SetFunction,
   // or the null string

   return fFunction;
}

void GlueXComputedMagField::SetMagFieldTransform(const G4AffineTransform &xform)
{
   // provides a means to modify the placement of the mapped field in
   // the geometry without having to reinitialize the map.
 
   fXform = xform;
   fXfinv = xform.Inverse();
}

const G4AffineTransform& GlueXComputedMagField::GetMagFieldTransform() const
{
   // returns a copy of the current map placement transform
 
   return fXform;
}
 
G4ThreeVector GlueXComputedMagField::GetMagField(const G4double point[4],
                                                 G4double unit) const
{
   // looks up the magnetic field at the specified point, and returns
   // the result scaled to the requested unit (eg. gauss). The first 3
   // elements of point are the Cartesian coordinates in the G4 world
   // coordinate system at which the map is to be evaluated. The time
   // stored in point[3] is presently ignored. The position is transformed
   // into local map coordinates for lookup, and then the field is rotated
   // back into user coordinates and returned in the requested units.
 
   G4ThreeVector p(point[0] / cm, point[1] / cm, point[2] / cm);
   fXfinv.ApplyPointTransform(p);
   double B[3] = {0,0,0};
   DMagneticFieldMap *mapso = fJanaFieldMap;
                     //dynamic_cast <DMagneticFieldMap* const> (fJanaFieldMap);
   DMagneticFieldMapPS *mapps = fJanaFieldMapPS;
                     //dynamic_cast <DMagneticFieldMapPS* const> (fJanaFieldMapPS);
   if (mapso)
      mapso->GetField(p[0], p[1], p[2], B[0], B[1], B[2]);
   else if (mapps)
      mapps->GetField(p[0], p[1], p[2], B[0], B[1], B[2]);
   G4ThreeVector Bvec(B[0], B[1], B[2]);
   fXform.ApplyAxisTransform(Bvec);
   if (Bvec[2] == 0)
      Bvec[2] = 1e-99; // avoid divide-by-zero by excluding |B|=0
   return Bvec *= fUnit / unit;
}

void GlueXComputedMagField::GetFieldValue(const G4double point[4],
                                          G4double *Bfield ) const
{
   // implements the standard G4MagneticField::GetFieldValue() interface

   G4ThreeVector Bvec = GetMagField(point, 1);
   Bfield[0] = Bvec[0];
   Bfield[1] = Bvec[1];
   Bfield[2] = Bvec[2];
}

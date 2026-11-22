//
// HddsG4Builder - class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//

#define APP_NAME "HddsG4Builder"

#include <HddsG4Builder.hh>
#include <G4Box.hh>
#include <G4Trd.hh>
#include <G4Trap.hh>
#include <G4Para.hh>
#include <G4Tubs.hh>
#include <G4Cons.hh>
#include <G4Hype.hh>
#include <G4Torus.hh>
#include <G4Orb.hh>
#include <G4Sphere.hh>
#include <G4Polycone.hh>
#include <G4Polyhedra.hh>
#include <G4EllipticalTube.hh>
#include <G4SystemOfUnits.hh>
#include <G4Mag_UsualEqRhs.hh>
#include <G4MagIntegratorStepper.hh>
#include <G4ExactHelixStepper.hh>
#include <G4HelixMixedStepper.hh>
#include <G4ClassicalRK4.hh>
#include <G4ChordFinder.hh>
#include <G4PVPlacement.hh>
#include <G4PVDivision.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4VisAttributes.hh>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLStringTokenizer.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>

using namespace xercesc;

#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <list>
#include <map>

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

#ifdef LINUX_CPUTIME_PROFILING
extern CPUtimer timer;
#endif

HddsG4Builder::HddsG4Builder() : fWorldVolume(0) { }

HddsG4Builder::HddsG4Builder(const HddsG4Builder &src)
{
   fWorldVolume = src.fWorldVolume;
   fElements = src.fElements;
   fMaterials = src.fMaterials;
   fMagneticRegions = src.fMagneticRegions;
   fLogicalVolumes = src.fLogicalVolumes;
   fPhysicalVolumes = src.fPhysicalVolumes;
   fSensitiveVolumes = src.fSensitiveVolumes;
   fRotations = src.fRotations;
   fCurrentMother = src.fCurrentMother;
   fCurrentPlacement = src.fCurrentPlacement;
   fCurrentDivision = src.fCurrentDivision;
   fCurrentPhiCenter = src.fCurrentPhiCenter;
}

HddsG4Builder::~HddsG4Builder() { }

int HddsG4Builder::createMaterial(DOMElement* el)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createMaterial: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int imate = CodeWriter::createMaterial(el);

   if (fSubst.fBrewList.size() == 0)
   {
      XString matS = fSubst.getName();
      XString symS = fSubst.getSymbol();
      double A = fSubst.getAtomicWeight();
      double Z = fSubst.getAtomicNumber();
      double dens = fSubst.getDensity();
      fElements[imate] = new G4Element(matS,symS,Z,A*g/mole);
      fMaterials[imate] = new G4Material(matS,Z,A*g/mole,dens*g/cm3);
   }
   else
   {
      XString matS = fSubst.getName();
      double dens = fSubst.getDensity();
      int ncomp = fSubst.fBrewList.size();
      fMaterials[imate] = new G4Material(matS,dens*g/cm3,ncomp);
      std::list<Substance::Brew>::iterator iter;
      for (iter = fSubst.fBrewList.begin();
           iter != fSubst.fBrewList.end(); iter++)
      {
         int subimate = iter->sub->fUniqueID;
         std::map<int,G4Element*>::iterator elemit = fElements.find(subimate);
         if (iter->natoms)
         {
            fMaterials[imate]->AddElement(elemit->second,iter->natoms);
         }
         else if (elemit != fElements.end())
         {
            fMaterials[imate]->AddElement(elemit->second,iter->wfact);
         }
         else
         {
            fMaterials[imate]->AddMaterial(fMaterials[subimate],iter->wfact);
         }
      }
   }

   DOMNodeList* propList = el->getElementsByTagName(X("optical_properties"));
   if (propList->getLength() > 0)
   {
      DOMElement* propEl = (DOMElement*)propList->item(0);
      DOMNodeList* specL = propEl->getElementsByTagName(X("specify"));
      int len = specL->getLength();
      std::vector<double> Ephot;
      std::vector<double> abslen;
      std::vector<double> effic;
      std::vector<double> rindex;
      std::vector<double> smooth;
      std::vector<double> reflect;
      for (int i=0; i < len; ++i)
      {
         DOMElement* specEl = (DOMElement*)specL->item(i);
         XString valS;
         valS = specEl->getAttribute(X("E"));
         Ephot.push_back(atof(S(valS))*eV);
         valS = specEl->getAttribute(X("refindex"));
         rindex.push_back(atof(S(valS)));
         valS = specEl->getAttribute(X("abslen"));
         abslen.push_back(atof(S(valS))*cm);
         valS = specEl->getAttribute(X("smooth"));
         smooth.push_back(atof(S(valS)));
         valS = specEl->getAttribute(X("reflect"));
         reflect.push_back(atof(S(valS)));
         valS = specEl->getAttribute(X("effic"));
         effic.push_back(atof(S(valS)));
      }
      G4MaterialPropertiesTable *mpt = new G4MaterialPropertiesTable();
      if (rindex[0] > 0)
      {
         mpt->AddProperty("RINDEX",&Ephot[0],&rindex[0],len);
      }
      if (abslen[0] > 0)
      {
         mpt->AddProperty("ABSLENGTH",&Ephot[0],&abslen[0],len);
      }
      if (smooth[0] > 0)
      {
         mpt->AddProperty("POLISH",&Ephot[0],&smooth[0],len);
      }
      if (reflect[0] > 0)
      {
         mpt->AddProperty("REFLECTIVITY",&Ephot[0],&reflect[0],len);
      }
      if (effic[0] > 0)
      {
         mpt->AddProperty("EFFICIENCY",&Ephot[0],&effic[0],len);
      }
      fMaterials[imate]->SetMaterialPropertiesTable(mpt);
   }
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
   return imate;
}

int HddsG4Builder::createSolid(DOMElement* el, Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createSolid: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif

   int ivolu = CodeWriter::createSolid(el,ref);
   XString nameS(el->getAttribute(X("name")));
   int imate = fSubst.fUniqueID;

   Units unit;
   unit.getConversions(el);

   G4VSolid *solid;
   XString shapeS(el->getTagName());
   if (shapeS == "box")
   {
      double xl, yl, zl;
      XString xyzS(el->getAttribute(X("X_Y_Z")));
      std::stringstream listr(xyzS);
      listr >> xl >> yl >> zl;
      solid = new G4Box(S(nameS),
                        xl/2 *cm/unit.cm,
                        yl/2 *cm/unit.cm,
                        zl/2 *cm/unit.cm
                       );
   }
   else if (shapeS == "tubs")
   {
      double ri, ro, zl, phi0, dphi;
      XString riozS(el->getAttribute(X("Rio_Z")));
      std::stringstream listr(riozS);
      listr >> ri >> ro >> zl;
      XString profS(el->getAttribute(X("profile")));
      listr.clear(), listr.str(profS);
      listr >> phi0 >> dphi;
      solid = new G4Tubs(S(nameS),
                         ri * cm/unit.cm,
                         ro * cm/unit.cm,
                         zl/2 * cm/unit.cm,
                         phi0 * deg/unit.deg,
                         dphi * deg/unit.deg
                        );
   }
   else if (shapeS == "eltu")
   {

      double rx, ry, zl;
      XString rxyzS(el->getAttribute(X("Rxy_Z")));
      std::stringstream listr(rxyzS);
      listr >> rx >> ry >> zl;
      solid = new G4EllipticalTube(S(nameS),
                                   rx * cm/unit.cm,
                                   ry * cm/unit.cm,
                                   zl/2 * cm/unit.cm
                                  );
   }
   else if (shapeS == "trd")
   {
      double xm, ym, xp, yp, zl;
      XString xyzS(el->getAttribute(X("Xmp_Ymp_Z")));
      std::stringstream listr(xyzS);
      listr >> xm >> xp >> ym >> yp >> zl;
      xm += 1e-100;
      xp += 1e-100;
      ym += 1e-100;
      yp += 1e-100;
      zl += 1e-100;
      double alph_xz, alph_yz;
      XString incS(el->getAttribute(X("inclination")));
      listr.clear(), listr.str(incS);
      listr >> alph_xz >> alph_yz;
      double x = tan(alph_xz/unit.rad);
      double y = tan(alph_yz/unit.rad);
      double r = sqrt(x*x + y*y);
      solid = new G4Trap(S(nameS),
                         zl/2 * cm/unit.cm,
                         atan2(r,1) * rad,
                         atan2(y,x) * rad,
                         ym/2 * cm/unit.cm,
                         xm/2 * cm/unit.cm,
                         xm/2 * cm/unit.cm,
                         0,
                         yp/2 * cm/unit.cm,
                         xp/2 * cm/unit.cm,
                         xp/2 * cm/unit.cm,
                         0
                        );
   }
   else if (shapeS == "pcon")
   {
      double phi0, dphi;
      XString profS(el->getAttribute(X("profile")));
      std::stringstream listr(profS);
      listr >> phi0 >> dphi;
      DOMNodeList* planeList = el->getElementsByTagName(X("polyplane"));
      int nplanes = planeList->getLength();
      std::vector<double> zPlane;
      std::vector<double> rInner;
      std::vector<double> rOuter;
      double zlast = -1e30;
      double zeps = 0;
      for (unsigned int p = 0; p < planeList->getLength(); p++)
      {
         double ri, ro, zl;
         DOMNode* node = planeList->item(p);
         DOMElement* elem = (DOMElement*) node;
         XString riozS(elem->getAttribute(X("Rio_Z")));
         std::stringstream listr1(riozS);
         listr1 >> ri >> ro >> zl;
         if (zl < zlast)
         {
            G4cerr << APP_NAME << " error: Please re-order the polyplanes"
                 << " of volume " << S(nameS) << " so that the z-values"
                 << " are non-decreasing"
                 << G4endl;
            exit(1);
         }
         zl = (zl < zlast + zeps)? zlast + zeps : zl;
         zlast = zl;
         zPlane.push_back(zl * cm/unit.cm);
         rInner.push_back(ri * cm/unit.cm);
         rOuter.push_back(ro * cm/unit.cm);
      }
      solid = new G4Polycone(S(nameS),
                             phi0 * deg/unit.deg,
                             dphi * deg/unit.deg,
                             nplanes,
                             &zPlane[0],
                             &rInner[0],
                             &rOuter[0]
                            );
   }
   else if (shapeS == "pgon")
   {
      int segments;
      XString segS(el->getAttribute(X("segments")));
      segments = atoi(S(segS));
      double phi0, dphi;
      XString profS(el->getAttribute(X("profile")));
      std::stringstream listr(profS);
      listr >> phi0 >> dphi;
      DOMNodeList* planeList = el->getElementsByTagName(X("polyplane"));
      int nplanes = planeList->getLength();
      std::vector<double> zPlane;
      std::vector<double> rInner;
      std::vector<double> rOuter;
      double zlast = -1e30;
      double zeps = 0;
      for (unsigned int p = 0; p < planeList->getLength(); p++)
      {
         double ri, ro, zl;
         DOMNode* node = planeList->item(p);
         DOMElement* elem = (DOMElement*) node;
         XString riozS(elem->getAttribute(X("Rio_Z")));
         std::stringstream listr1(riozS);
         listr1 >> ri >> ro >> zl;
         if (zl < zlast)
         {
            G4cerr << APP_NAME << " error: Please re-order the polyplanes"
                   << " of volume " << S(nameS) << " so that the z-values"
                   << " are non-decreasing."
                   << G4endl;
            exit(1);
         }
         zl = (zl < zlast + zeps)? zlast + zeps : zl;
         zlast = zl;
         zPlane.push_back(zl * cm/unit.cm);
         rInner.push_back(ri * cm/unit.cm);
         rOuter.push_back(ro * cm/unit.cm);
      }
      solid = new G4Polyhedra(S(nameS),
                              phi0 * deg/unit.deg,
                              dphi * deg/unit.deg,
                              segments,
                              nplanes,
                              &zPlane[0],
                              &rInner[0],
                              &rOuter[0]
                             );
   }
   else if (shapeS == "cons")
   {
      double rim, rip, rom, rop, zl;
      XString riozS(el->getAttribute(X("Rio1_Rio2_Z")));
      std::stringstream listr(riozS);
      listr >> rim >> rom >> rip >> rop >> zl;
      double phi0, dphi;
      XString profS(el->getAttribute(X("profile")));
      listr.clear(), listr.str(profS);
      listr >> phi0 >> dphi;
      solid = new G4Cons(S(nameS),
                         zl/2 * cm/unit.cm,
                         rim * cm/unit.cm,
                         rom * cm/unit.cm,
                         rip * cm/unit.cm,
                         rop * cm/unit.cm,
                         phi0 * deg/unit.deg,
                         dphi * deg/unit.deg
                        );
   }
   else if (shapeS == "sphere")
   {
      double ri, ro;
      XString rioS(el->getAttribute(X("Rio")));
      std::stringstream listr(rioS);
      listr >> ri >> ro;
      double theta0, theta1;
      XString polarS(el->getAttribute(X("polar_bounds")));
      listr.clear(), listr.str(polarS);
      listr >> theta0 >> theta1;
      double phi0, dphi;
      XString profS(el->getAttribute(X("profile")));
      listr.clear(), listr.str(profS);
      listr >> phi0 >> dphi;
      solid = new G4Sphere(S(nameS),
                           ri * cm/unit.cm,
                           ro * cm/unit.cm,
                           phi0 * deg/unit.deg,
                           dphi * deg/unit.deg,
                           theta0 * deg/unit.deg,
                           theta1 * deg/unit.deg
                          );
   }
   else
   {
      G4cerr << APP_NAME << " error: volume " << S(nameS)
             << " should be one of the valid shapes, not " << S(shapeS)
             << G4endl;
      exit(1);
   }

   vpair_t newvol(ivolu,0);
   G4FieldManager* fieldmgr = fFieldManagers[ref.fRegionID];
   fLogicalVolumes[newvol] = new G4LogicalVolume(solid,fMaterials[imate],
                                                 S(nameS),fieldmgr);
   G4Material* mate = fMaterials[imate];
   double dens = mate->GetDensity();
   double dmod = (int)(dens*97.345/(g/cm3)) % 20;
   double red = 1-3/(3+dens /(g/cm3));
   double green = 1-5/(5+dmod);
   double blue = 1-green;
   G4VisAttributes* attr = new G4VisAttributes(G4Colour(red,green,blue));
   fLogicalVolumes[newvol]->SetVisAttributes(attr);

   XString sensiS(el->getAttribute(X("sensitive")));
   if (sensiS == "true")
   {
      fSensitiveVolumes[ivolu] = fLogicalVolumes[newvol];
   }

   G4MaterialPropertiesTable *mpt = fMaterials[imate]->GetMaterialPropertiesTable();
   if (mpt)
   {
      G4MaterialPropertyVector *poli_vector = mpt->GetProperty("POLISH");
      G4MaterialPropertyVector *refl_vector = mpt->GetProperty("REFLECTIVITY");
      G4MaterialPropertyVector *refi_vector = mpt->GetProperty("RINDEX");
      G4MaterialPropertyVector *absl_vector = mpt->GetProperty("ABSLENGTH");
      G4MaterialPropertyVector *effi_vector = mpt->GetProperty("EFFICIENCY");
      if (refi_vector == 0) {
         if (refl_vector == 0) {
            G4cerr << APP_NAME << " error: material " << mate->GetName()
                   << " needs to have either refindex defined (dielectrics)"
                   << " or else reflect (metals) but it has neither!"
                   << G4endl;
            exit(1);
         }
         if (absl_vector != 0) {
            G4cerr << APP_NAME << " error: material " << mate->GetName()
                   << " has abslen defined but no refindex; defining abslen"
                   << " requires that refindex also be specified!"
                   << G4endl;
            exit(1);
         }
         if (poli_vector != 0 || effi_vector != 0) {
            G4OpticalSurface *surface = new G4OpticalSurface(S(nameS));
            surface->SetType(dielectric_metal);
            surface->SetModel(glisur);
            surface->SetMaterialPropertiesTable(mpt);
            if (poli_vector != 0) {
               double polish = poli_vector->GetMaxValue();
               surface->SetPolish(polish);
               surface->SetFinish(ground);
            }
            else {
               surface->SetFinish(polished);
            }
            new G4LogicalSkinSurface(S(nameS),fLogicalVolumes[newvol],surface);
         }
      }
      else if (refl_vector == 0) {
         if (poli_vector != 0) {
            G4OpticalSurface *surface = new G4OpticalSurface(S(nameS));
            surface->SetType(dielectric_dielectric);
            surface->SetModel(glisur);
            surface->SetMaterialPropertiesTable(mpt);
            if (poli_vector != 0) {
               double polish = poli_vector->GetMaxValue();
               surface->SetPolish(polish);
               surface->SetFinish(ground);
            }
            else {
               surface->SetFinish(polished);
            }
            new G4LogicalSkinSurface(S(nameS),fLogicalVolumes[newvol],surface);
         }
      }
      else {
         G4cerr << APP_NAME << " error: material " << mate->GetName()
                << " has optical properties, but neither refindex nor reflect"
                << " data are provided, you need either one or the other!"
                << G4endl;
         exit(1);
      }
   }

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
   return ivolu;
}

int HddsG4Builder::createRotation(Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createRotation: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int irot = CodeWriter::createRotation(ref);

   if (irot > 0)
   {
      std::vector<double> omega = ref.getRotation();
      fRotations[irot] = new G4RotationMatrix();
      fRotations[irot]->rotateX(omega[0]);
      fRotations[irot]->rotateY(omega[1]);
      fRotations[irot]->rotateZ(omega[2]);
      fRotations[irot]->invert();
   }
   else if (fRotations.count(0) == 0)
   {
      fRotations[0] = new G4RotationMatrix();
   }

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
   return irot;
}

int HddsG4Builder::createVolume(DOMElement* el, Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createVolume: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int icopy = CodeWriter::createVolume(el,ref);

   if (fPending)
   {
      XString myvoluS(el->getAttribute(X("HDDSvolu")));
      XString motherS(fRef.fMother->getAttribute(X("HDDSvolu")));
      int myvoluI = atoi(S(myvoluS));
      int motherI = atoi(S(motherS));
      G4ThreeVector origin(fRef.fOrigin[0]*cm,
                           fRef.fOrigin[1]*cm,
                           fRef.fOrigin[2]*cm);
      int irot = fRef.fRotation;

      // Apply fix-up in case we are placing the volume inside a phi division
      // because of a flaw in the way geant4 handles this special case.

      G4RotationMatrix *rot = fRotations[irot];
      if (fCurrentPhiCenter[motherI] != 0)
      {
         rot = new G4RotationMatrix(*rot);
         rot->rotateZ(-fCurrentPhiCenter[motherI]);
         origin.rotateZ(fCurrentPhiCenter[motherI]);
      }

      std::map<vpair_t,G4LogicalVolume*>::iterator mine;
      for (mine = fLogicalVolumes.find(vpair_t(myvoluI,0));
           mine != fLogicalVolumes.end() && mine->first.first == myvoluI;
           ++mine)
      {
         int ilayer = mine->first.second;
         std::map<vpair_t,G4LogicalVolume*>::iterator moms;
         moms = addNewLayer(motherI,fRef.fRelativeLayer + ilayer);
         G4PVPlacement *pvol = new G4PVPlacement(rot, origin,
                                                 mine->second,
                                                 mine->second->GetName(),
                                                 moms->second,0,icopy);
#ifdef CHECK_OVERLAPS_MM
         pvol->CheckOverlaps(1000,CHECK_OVERLAPS_MM);
#endif
#ifdef DEBUG_PLACEMENT
         XString nameS(el->getAttribute(X("name")));
         G4cout << "volume " << nameS 
                << "->" << mine->second->GetName()
                << "->" << mine->second->GetSolid()->GetName()
                << " being placed in mother " << moms->second->GetName()
                << "->" << moms->second->GetSolid()->GetName() << " at "
                << fRef.fOrigin[0]*cm << ","
                << fRef.fOrigin[1]*cm << "," 
                << fRef.fOrigin[2]*cm
                << G4endl;
#endif

         if (ilayer == 0) {
            vpair_t mycopy(myvoluI,icopy);
            fPhysicalVolumes[mycopy] = pvol;
            fCurrentPlacement[myvoluI] = fPhysicalVolumes.find(mycopy);
            fCurrentMother[myvoluI] = moms;
         }
      }

      fPending = false;
   }

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
   return icopy;
}

std::map<HddsG4Builder::vpair_t,G4LogicalVolume*>::iterator
HddsG4Builder::addNewLayer(int volume_id, int layer)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "addNewLayer: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif

   // Check if the volume already exists on the specified layer,
   // and if so, return that, otherwise obtain the maximum layer
   // index that exists already for this volume.
 
   int max_layer = -1;
   std::map<vpair_t,G4LogicalVolume*>::iterator volex;
   for (volex = fLogicalVolumes.find(vpair_t(volume_id,0));
        volex != fLogicalVolumes.end() &&
        volex->first.first == volume_id;
        ++volex)
   {
      max_layer = volex->first.second;
      if (max_layer == layer)
      {
         return volex;
      }
   }
   if (max_layer == -1)
   {
      G4cerr << "HddsG4Builder::addNewLayer called for volume "
             << volume_id << " which does not (yet) exist!"
             << G4endl
             << "Unable to continue, giving up."
             << G4endl;
      exit(1);
   }

   // The volume does not exist yet on the specified layer, so we
   // want to create a reflection of the existing layer 0 volume
   // and place it on the specified layer.  Before that can happen,
   // the volume's mother must exist on the specified layer.
   // Use recursion to ensure that this is always true.
 
   std::map<vpair_t,G4LogicalVolume*>::iterator moms;
   if (fCurrentMother.find(volume_id) != fCurrentMother.end())
   {
      moms = fCurrentMother[volume_id];
      moms = addNewLayer(moms->first.first, moms->first.second + layer);
   }
   else
   {
      moms = fLogicalVolumes.end();
   }

   // Create a new copy of the existing max_layer volume.  If the
   // new layer index is higher than max_layer then set the interior
   // filling material to null in the max_layer copy, otherwise 
   // the new copy should have a null material.

   vpair_t newvol(volume_id,layer);
   vpair_t oldvol(volume_id,max_layer);
   G4VSolid *solid = fLogicalVolumes[oldvol]->GetSolid()->Clone();
   G4String name(fLogicalVolumes[vpair_t(volume_id,0)]->GetName());
   std::stringstream str;
   str << name << "::" << layer;
   solid->SetName(str.str());
   G4Material *material = fLogicalVolumes[oldvol]->GetMaterial();
   G4FieldManager *fieldmgr = fLogicalVolumes[oldvol]->GetFieldManager();
   if (layer > max_layer)
   {
      fLogicalVolumes[newvol] = new G4LogicalVolume(solid,material,
                                                    str.str(),fieldmgr);
      fLogicalVolumes[oldvol]->SetMaterial(0);
   }
   else
   {
      fLogicalVolumes[newvol] = new G4LogicalVolume(solid,0,str.str(),fieldmgr);
   }

   // If this is the top-level (world) volume then it cannot be
   // placed here, so simply return here.

   if (moms == fLogicalVolumes.end())
   {
      return fLogicalVolumes.find(newvol);
   }

   // Now we have a copy of the mother volume on the desired layer,
   // and we have a new instance of our volume to place inside her.
   // Two cases are supported in the present code, although more
   // could be easily added.
 
   //   case 1: G4PVPlacement 
 
   if (fCurrentPlacement.find(volume_id) != fCurrentPlacement.end())
   {
      std::map<vpair_t,G4VPhysicalVolume*>::iterator placement;
      placement = fCurrentPlacement[volume_id];
      G4PVPlacement *player0 = (G4PVPlacement*)placement->second;
      vpair_t mycopy(volume_id,placement->first.second);
      G4PVPlacement *playerN;
      playerN = new G4PVPlacement(player0->GetRotation(),
                                  player0->GetTranslation(),
                                  fLogicalVolumes[newvol],
                                  str.str(),
                                  moms->second,0,
                                  placement->first.second);
#ifdef CHECK_OVERLAPS_MM
      playerN->CheckOverlaps(1000,CHECK_OVERLAPS_MM);
#else
      (void)playerN;
#endif
#ifdef DEBUG_PLACEMENT 
      G4cout << "volume " << str.str() 
             << "->" << fLogicalVolumes[newvol]->GetName()
             << "->" << fLogicalVolumes[newvol]->GetSolid()->GetName()
             << " being placed in mother " << moms->second->GetName()
             << "->" << moms->second->GetSolid()->GetName()
             << G4endl;
#endif
   }

   //   case 2: G4PVDivision
 
   if (fCurrentDivision.find(volume_id) != fCurrentDivision.end())
   {
      std::map<vpair_t,G4VPhysicalVolume*>::iterator division;
      division = fCurrentDivision[volume_id];
      G4PVDivision *division0 = (G4PVDivision*)division->second;
      EAxis axis;
      G4int ndiv;
      G4double width;
      G4double offset;
      G4bool consuming;
      division0->GetReplicationData(axis,ndiv,width,offset,consuming);
      // axis returned by GetReplicationData is sometimes wrong!
      axis = division0->GetDivisionAxis();
      vpair_t mydiv(volume_id,division->first.second);
      G4PVDivision *divisionN;
      divisionN = new G4PVDivision(str.str(),
                                   fLogicalVolumes[newvol],
                                   moms->second,
                                   axis, ndiv, width, offset);
      (void)divisionN;
      fLogicalVolumes[newvol]->SetVisAttributes(new G4VisAttributes(false));
#ifdef DEBUG_PLACEMENT 
      G4cout << ndiv << " copies of division " << str.str() 
             << "->" << fLogicalVolumes[newvol]->GetName()
             << "->" << fLogicalVolumes[newvol]->GetSolid()->GetName()
             << " being placed in mother " << moms->second->GetName()
             << "->" << moms->second->GetSolid()->GetName()
             << " offset by " << offset << " on axis " << axis
             << " and repeating every " << width
             << G4endl;
#endif
   }

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
   return fLogicalVolumes.find(newvol);
}

int HddsG4Builder::createDivision(XString& divStr, Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createDivision: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int ndiv = CodeWriter::createDivision(divStr,ref);

   XString motherS(ref.fMother->getAttribute(X("HDDSvolu")));
   XString myvoluS(ref.fPartition.divEl->getAttribute(X("HDDSvolu")));
   int motherI = atoi(S(motherS));
   int myvoluI = atoi(S(myvoluS));
   std::map<vpair_t,G4LogicalVolume*>::iterator moms;
   moms = fLogicalVolumes.find(vpair_t(motherI,0));
   G4VSolid* solid = moms->second->GetSolid()->Clone();
   solid->SetName(divStr);
   EAxis axis;
   double width,offset;

   if (ref.fPartition.axis == "x")
   {
      axis = kXAxis;
      width = ref.fPartition.step * cm;
      offset = ref.fPartition.offset * cm;
   }
   else if (ref.fPartition.axis == "y")
   {
      axis = kYAxis;
      width = ref.fPartition.step * cm;
      offset = ref.fPartition.offset * cm;
   }
   else if (ref.fPartition.axis == "z")
   {
      axis = kZAxis;
      width = ref.fPartition.step * cm;
      offset = ref.fPartition.offset * cm;
   }
   else if (ref.fPartition.axis == "rho")
   {
      axis = kRho;
      width = ref.fPartition.step * cm;
      offset = ref.fPartition.offset * cm;
   }
   else if (ref.fPartition.axis == "phi")
   {
      axis = kPhi;
      width = ref.fPartition.step * deg;
      offset = ref.fPartition.offset * deg;
   }
   else
   {
      XString nameS(ref.fMother->getAttribute(X("name")));
      G4cerr << APP_NAME << " error: Volume division is requested "
             << " for volume " << S(nameS) << " but division axis "
             << "\"" << ref.fPartition.axis << "\""
             << " is not currently supported."
             << G4endl;
      exit(1);
   }

   vpair_t mothvol(moms->first);
   vpair_t myvol(myvoluI,moms->first.second);
   G4Material* material = moms->second->GetMaterial();
   G4FieldManager* fieldmgr = fFieldManagers[ref.fRegionID];
   fLogicalVolumes[myvol] = new G4LogicalVolume(solid,material,
                                                S(divStr),fieldmgr);
   fLogicalVolumes[myvol]->SetVisAttributes(new G4VisAttributes(false));
   vpair_t mydiv(myvoluI,0);
   fPhysicalVolumes[mydiv] = new G4PVDivision(S(divStr),
                                              fLogicalVolumes[myvol],
                                              fLogicalVolumes[mothvol],
                                              axis,ndiv,width,offset);
   fCurrentDivision[myvoluI] = fPhysicalVolumes.find(mydiv);
   fCurrentMother[myvoluI] = moms;

   // recompute the limits of the solid representing this division
 
   G4VPVParameterisation *param = fPhysicalVolumes[mydiv]
                                     ->GetParameterisation();
   if (dynamic_cast<G4Box*>(solid))
   {
      param->ComputeDimensions(*(G4Box*)solid,0,fPhysicalVolumes[mydiv]);
   }
   else if (dynamic_cast<G4Tubs*>(solid))
   {
      param->ComputeDimensions(*(G4Tubs*)solid,0,fPhysicalVolumes[mydiv]);
      double start = ((G4Tubs*)solid)->GetStartPhiAngle();
      double delta = ((G4Tubs*)solid)->GetDeltaPhiAngle();
      fCurrentPhiCenter[myvoluI] = start + delta/2;
   }
   else if (dynamic_cast<G4Trd*>(solid))
   {
      param->ComputeDimensions(*(G4Trd*)solid,0,fPhysicalVolumes[mydiv]);
   }
   else if (dynamic_cast<G4Trap*>(solid))
   {
      param->ComputeDimensions(*(G4Trap*)solid,0,fPhysicalVolumes[mydiv]);
   }
   else if (dynamic_cast<G4Cons*>(solid))
   {
      param->ComputeDimensions(*(G4Cons*)solid,0,fPhysicalVolumes[mydiv]);
      double start = ((G4Cons*)solid)->GetStartPhiAngle();
      double delta = ((G4Cons*)solid)->GetDeltaPhiAngle();
      fCurrentPhiCenter[myvoluI] = start + delta/2;
   }
   else if (dynamic_cast<G4Sphere*>(solid))
   {
      param->ComputeDimensions(*(G4Sphere*)solid,0,fPhysicalVolumes[mydiv]);
      double start = ((G4Sphere*)solid)->GetStartPhiAngle();
      double delta = ((G4Sphere*)solid)->GetDeltaPhiAngle();
      fCurrentPhiCenter[myvoluI] = start + delta/2;
   }
   else if (dynamic_cast<G4Orb*>(solid))
   {
      param->ComputeDimensions(*(G4Orb*)solid,0,fPhysicalVolumes[mydiv]);
   }
   else if (dynamic_cast<G4Torus*>(solid))
   {
      param->ComputeDimensions(*(G4Torus*)solid,0,fPhysicalVolumes[mydiv]);
      double start = ((G4Torus*)solid)->GetSPhi();
      double delta = ((G4Torus*)solid)->GetDPhi();
      fCurrentPhiCenter[myvoluI] = start + delta/2;
   }
   else if (dynamic_cast<G4Para*>(solid))
   {
      param->ComputeDimensions(*(G4Para*)solid,0,fPhysicalVolumes[mydiv]);
   }
   else if (dynamic_cast<G4Polycone*>(solid))
   {
      param->ComputeDimensions(*(G4Polycone*)solid,0,fPhysicalVolumes[mydiv]);
      G4PolyconeHistorical *params;
      params = ((G4Polycone*)solid)->GetOriginalParameters();
      double start = params->Start_angle;
      double delta = params->Opening_angle;
      fCurrentPhiCenter[myvoluI] = start + delta/2;
   }
   else if (dynamic_cast<G4Polyhedra*>(solid))
   {
      param->ComputeDimensions(*(G4Polyhedra*)solid,0,fPhysicalVolumes[mydiv]);
      G4PolyhedraHistorical *params;
      params = ((G4Polyhedra*)solid)->GetOriginalParameters();
      double start = params->Start_angle;
      double delta = params->Opening_angle;
      fCurrentPhiCenter[myvoluI] = start + delta/2;
   }
   else if (dynamic_cast<G4Hype*>(solid))
   {
      param->ComputeDimensions(*(G4Hype*)solid,0,fPhysicalVolumes[mydiv]);
   }
   else
   {
      G4cerr << "HddsG4Builder::createDivision called for unsupported "
             << "solid type!" << G4endl
             << "Unable to continue, giving up."
             << G4endl;
      exit(1);
   }
#ifdef DEBUG_PLACEMENT
   G4cout << ndiv << " copies of division " << divStr 
          << "->" << fLogicalVolumes[myvol]->GetName()
          << "->" << fLogicalVolumes[myvol]->GetSolid()->GetName()
          << " being placed in mother " << fLogicalVolumes[mothvol]->GetName()
          << "->" << fLogicalVolumes[mothvol]->GetSolid()->GetName()
          << " offset by " << offset << " on axis " << axis
          << " and repeating every " << width
          << G4endl;
#endif

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
   return ndiv;
}

int HddsG4Builder::createRegion(DOMElement* el, Refsys& ref)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createRegion: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   int iregion = CodeWriter::createRegion(el,ref);

   XString motherS(ref.fMother->getAttribute(X("HDDSvolu")));
   int motherI = atoi(S(motherS));
   std::map<vpair_t,G4LogicalVolume*>::iterator moms;
   moms = fLogicalVolumes.find(vpair_t(motherI,0));

   G4double *R = (G4double*)ref.fMRmatrix;
   G4double *O = ref.fMOrigin;
   G4AffineTransform xform(CLHEP::HepRotation(CLHEP::HepRep3x3(R)),
                           CLHEP::Hep3Vector(O[0],O[1],O[2]));

   if (ref.fRegion)
   {
      DOMNodeList* noBfieldL;
      DOMNodeList* uniBfieldL;
      DOMNodeList* mapBfieldL;
      DOMNodeList* compBfieldL;
      DOMNodeList* swimL;
      noBfieldL = ref.fRegion->getElementsByTagName(X("noBfield"));
      uniBfieldL = ref.fRegion->getElementsByTagName(X("uniformBfield"));
      mapBfieldL = ref.fRegion->getElementsByTagName(X("mappedBfield"));
      compBfieldL = ref.fRegion->getElementsByTagName(X("computedBfield"));
      swimL = ref.fRegion->getElementsByTagName(X("swim"));
      double maxArcStep = 0;
      XString methodS("helix");
      if (swimL->getLength() > 0)
      {
         DOMElement* swimEl = (DOMElement*)swimL->item(0);
         methodS = swimEl->getAttribute(X("method"));
         Units unit;
         unit.getConversions(swimEl);
         XString stepS(swimEl->getAttribute(X("maxArcStep")));
         maxArcStep = atof(S(stepS))/unit.rad;
      }
      if (noBfieldL->getLength() > 0)
      {
         fMagneticRegions[iregion] = 0;
         G4ThreeVector Bvec(0,0,1e-99);
         G4UniformMagField *fld = new G4UniformMagField(Bvec);
         G4Mag_EqRhs *eqn = new G4Mag_UsualEqRhs(fld);
         G4MagIntegratorStepper *stepper = new G4ExactHelixStepper(eqn);
         G4ChordFinder *cfinder = new G4ChordFinder(fld, 1e+99, stepper);
         fFieldManagers[iregion] = new G4FieldManager(fld, cfinder);
      }
      else if (uniBfieldL->getLength() > 0)
      {
         Units funit;
         DOMElement* uniBfieldEl = (DOMElement*)uniBfieldL->item(0);
         funit.getConversions(uniBfieldEl);
         XString bvecS(uniBfieldEl->getAttribute(X("Bx_By_Bz")));
         std::stringstream str(S(bvecS));
         double B[3];
         str >> B[0] >> B[1] >> B[2];
         double u = kilogauss/funit.kG;
         G4ThreeVector Bvec(B[0],B[1],B[2]);
         GlueXUniformMagField *fld = new GlueXUniformMagField(Bvec,u,xform);
         fMagneticRegions[iregion] = fld;
         G4Mag_EqRhs *eqn = new G4Mag_UsualEqRhs(fld);
         G4MagIntegratorStepper *stepper = new G4ExactHelixStepper(eqn);
         G4ChordFinder *cfinder = new G4ChordFinder(fld, 0.01, stepper);
         fFieldManagers[iregion] = new G4FieldManager(fld, cfinder);
      }
      else if (mapBfieldL->getLength() > 0)
      {
         Units funit;
         DOMElement* mapBfieldEl = (DOMElement*)mapBfieldL->item(0);
         funit.getConversions(mapBfieldEl);
         XString bvecS(mapBfieldEl->getAttribute(X("maxBfield")));
         std::stringstream str(S(bvecS));
         double Bmax;
         str >> Bmax;
         double u = kilogauss/funit.kG;
         GlueXMappedMagField *fld = new GlueXMappedMagField(Bmax,u,xform);
         fMagneticRegions[iregion] = fld;
         G4Mag_EqRhs *eqn = new G4Mag_UsualEqRhs(fld);
         G4MagIntegratorStepper *stepper;
         if (methodS == "RungeKutta") {
            stepper = new G4ClassicalRK4(eqn);
         }
         else {
            stepper = new G4HelixMixedStepper(eqn);
         }
         G4ChordFinder *cfinder = new G4ChordFinder(fld, 0.01, stepper);
         if (maxArcStep > 0) {
            double rmin = 0.1 / (0.03 * Bmax * u) * meter;
            double max_miss = rmin * (1 - cos(maxArcStep / 2));
            cfinder->SetDeltaChord(max_miss);
         }
         fFieldManagers[iregion] = new G4FieldManager(fld, cfinder);
      }
      else if (compBfieldL->getLength() > 0)
      {
         Units funit;
         DOMElement* compBfieldEl = (DOMElement*)compBfieldL->item(0);
         funit.getConversions(compBfieldEl);
         XString bvecS(compBfieldEl->getAttribute(X("maxBfield")));
         std::stringstream str(S(bvecS));
         double Bmax;
         str >> Bmax;
         double u = kilogauss/funit.kG;
         GlueXComputedMagField *fld = new GlueXComputedMagField(Bmax,u,xform);
         fld->SetFunction(XString(compBfieldEl->getAttribute(X("function"))));
         fMagneticRegions[iregion] = fld;
         G4Mag_EqRhs *eqn = new G4Mag_UsualEqRhs(fld);
         G4MagIntegratorStepper *stepper;
         if (methodS == "RungeKutta") {
            stepper = new G4ClassicalRK4(eqn);
         }
         else {
            stepper = new G4HelixMixedStepper(eqn);
         }
         G4ChordFinder *cfinder = new G4ChordFinder(fld, 0.01, stepper);
         if (maxArcStep > 0) {
            double rmin = 0.1 / (0.03 * Bmax * u) * meter;
            double max_miss = rmin * (1 - cos(maxArcStep / 2));
            cfinder->SetDeltaChord(max_miss);
         }
         fFieldManagers[iregion] = new G4FieldManager(fld, cfinder);
      }
   }

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
   return iregion;
}

void HddsG4Builder::createSetFunctions(DOMElement* el, const XString& ident)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createSetFunctions: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createSetFunctions(el,ident);

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
}

void HddsG4Builder::createGetFunctions(DOMElement* el, const XString& ident)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createGetFunctions: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createGetFunctions(el,ident);

   std::vector<int> table;
   std::vector<int> start;
   start.push_back(0);

   XString funcNameStr;
   XString identCaps(ident);
   identCaps[0] = toupper(identCaps[0]);
   funcNameStr = "get" + identCaps;
   for (int ivolu = 1; ivolu <= Refsys::fVolumes; ivolu++)
   {
      start.push_back(0);
      int ncopy = Refsys::fIdentifierTable[ivolu]["copy counter"].back();
      std::map<std::string,std::vector<int> >::iterator idlist = 
                  Refsys::fIdentifierTable[ivolu].find(ident);
      if (idlist != Refsys::fIdentifierTable[ivolu].end())
      {
         if (ncopy != (int)idlist->second.size())
         {
            G4cerr << APP_NAME << " warning: volume " << ivolu
                   << " has " << ncopy << " copies, but "
                   << idlist->second.size() << " " 
                   << ident << " identifiers!" << G4endl;
            for (int idx = 0; idx < (int)idlist->second.size(); idx++)
            {
               G4cerr << idlist->second[idx]  << " ";
               if (idx/20*20 == idx)
                  G4cerr << G4endl;
            }
            G4cerr << G4endl;
         }
         start[ivolu] = table.size() + 1;
         for (int icopy = 0; icopy < ncopy; icopy++)
         {
            table.push_back(idlist->second[icopy]);
         }
      }
      else
      {
         start[ivolu] = 0;
      }
   }

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
}

void HddsG4Builder::createMapFunctions(DOMElement* el, const XString& ident)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "createMapFunctions: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   CodeWriter::createMapFunctions(el,ident);

   if (el == 0)
   {
      return;
   }
   DOMNodeList* regionsL = el->getOwnerDocument()->getDocumentElement()
                             ->getElementsByTagName(X("regions"));
   if (regionsL->getLength() == 0)
   {
      return;
   }
   DOMElement* regionsEl = (DOMElement*)regionsL->item(0);
   DOMNodeList* regionL = regionsEl->getElementsByTagName(X("region"));
   for (int ireg=0; ireg < (int)regionL->getLength(); ++ireg)
   {
      DOMElement* regionEl = (DOMElement*)regionL->item(ireg);
      DOMNodeList* mapfTagL = regionEl->getElementsByTagName(X("mappedBfield"));
      if (mapfTagL->getLength() == 0)
         continue;
      DOMElement* mapfEl = (DOMElement*)mapfTagL->item(0);
      XString nameS(regionEl->getAttribute(X("name")));
      XString iregionS(regionEl->getAttribute(X("HDDSregion")));
      if (iregionS.size() == 0)
         continue;
      int iregion = atoi(S(iregionS));
      GlueXMappedMagField *magfield = (GlueXMappedMagField*)
                                      fMagneticRegions[iregion];
      int axorder[] = {0,0,0,0};
      int axsamples[] = {0,0,0,0};

      XString gridtype;
      DOMNodeList* gridL = mapfEl->getElementsByTagName(X("grid"));
      int ngrid;
      for (ngrid = 0; ngrid < (int)gridL->getLength(); ++ngrid)
      {
         int axsense[] = {1,1,1,1};
         double axlower[4], axupper[4];
         DOMElement* gridEl = (DOMElement*)gridL->item(ngrid);
         XString typeS(gridEl->getAttribute(X("type")));
         if (gridtype.size() > 0 && typeS != gridtype)
         {
            G4cerr << APP_NAME << " error: mappedBfield in region " << S(nameS)
                   << " superimposes incompatible grid types." << G4endl;
            exit(1);
         }
         gridtype = typeS;

         DOMNodeList* samplesL = gridEl->getElementsByTagName(X("samples"));
         if (samplesL->getLength() != 3)
         {
            G4cerr << APP_NAME << " error: mappedBfield in region " << S(nameS)
                   << " does not have samples for three axes." << G4endl;
            exit(1);
         }

         for (int iax = 1; iax <= 3; ++iax)
         {
            DOMElement* sampleEl = (DOMElement*)samplesL->item(iax-1);
            XString nS(sampleEl->getAttribute(X("n")));
            XString axisS(sampleEl->getAttribute(X("axis")));
            XString boundsS(sampleEl->getAttribute(X("bounds")));
            XString senseS(sampleEl->getAttribute(X("sense")));
            Units sunit;
            double bound[2];
            sunit.getConversions(sampleEl);
            std::stringstream listr(boundsS);
            listr >> bound[0] >> bound[1];
            int iaxis=0;
            if (gridtype == "cartesian")
            {
               if (axisS == "x" && 
                  (axorder[0] == 0 || axorder[0] == iax))
               {
                  iaxis = 1;
                  axorder[0] = iax;
                  bound[0] /= sunit.cm;
                  bound[1] /= sunit.cm;
               }
               else if (axisS == "y" && 
                       (axorder[1] == 0 || axorder[1] == iax))
               {
                  iaxis = 2;
                  axorder[1] = iax;
                  bound[0] /= sunit.cm;
                  bound[1] /= sunit.cm;
               }
               else if (axisS == "z" &&
                       (axorder[2] == 0 || axorder[2] == iax))
               {
                  iaxis = 3;
                  axorder[2] = iax;
                  bound[0] /= sunit.cm;
                  bound[1] /= sunit.cm;
               }
               else
               {
                  G4cerr << APP_NAME << " error: grid in region " << S(nameS)
                         << " contains an incompatible set of samples." 
                         << G4endl;
                  exit(1);
               }
            }
            else if (gridtype == "cylindrical")
            {
               if (axisS == "r" && axorder[0] == 0)
               {
                  iaxis = 1;
                  axorder[0] = iax;
                  bound[0] /= sunit.cm;
                  bound[1] /= sunit.cm;
               }
               else if (axisS == "phi" && axorder[1] == 0)
               {
                  iaxis = 2;
                  axorder[1] = iax;
                  bound[0] /= sunit.rad;
                  bound[1] /= sunit.rad;
               }
               else if (axisS == "z" && axorder[2] == 0)
               {
                  iaxis = 3;
                  axorder[2] = iax;
                  bound[0] /= sunit.cm;
                  bound[1] /= sunit.cm;
               }
               else
               {
                  G4cerr << APP_NAME << " error: grid in region " << S(nameS)
                         << " contains an incompatible set of samples." 
                         << G4endl;
                  exit(1);
               }
            }

            int n = atoi(S(nS));
            if (axsamples[iaxis] == 0 || axsamples[iaxis] == n)
            {
               axsamples[iaxis] = n;
            }
            else
            {
               G4cerr << APP_NAME << " error: mappedBfield in region " << S(nameS)
                      << " combines incompatible grid elements." 
                      << G4endl;
               exit(1);
            }

            axlower[iaxis] = bound[0];
            axupper[iaxis] = bound[1];
            if (senseS == "reverse")
            {
               axsense[iaxis] = -1;
            }
         }

         if (gridtype == "cartesian")
         {
            magfield->AddCartesianGrid(axsamples, axorder, axsense,
                                       axlower, axupper);
         }
         else if (gridtype == "cylindrical")
         {
            magfield->AddCylindricalGrid(axsamples, axorder, axsense,
                                         axlower, axupper);
         }
         else
         {
            G4cerr << APP_NAME << " error: unrecognized grid type " 
                   << S(gridtype) << G4endl;
            exit(1);
         }
      }

      XString mapS(mapfEl->getAttribute(X("map")));
      XString encS(mapfEl->getAttribute(X("encoding")));
      if (encS != "utf-8")
      {
         G4cerr << APP_NAME << " error: mappedBfield in region " << S(nameS)
                << " uses unsupported encoding " << encS << G4endl;
         exit(1);
      }
      else if (mapS.substr(0,7) != "file://")
      {
         G4cerr << APP_NAME << " error: mappedBfield in region " << S(nameS)
                << " uses unsupported map URL " << mapS << G4endl;
         exit(1);
      }
      mapS.erase(0,7);
      magfield->ReadMapFile(S(mapS));
   }

   // Apply a post-build fix to ensure that every point in the geometry
   // has a consistent magnetic field on all layers, using recursion.

   addReflections(1);

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
}

void HddsG4Builder::translate(DOMElement* topel)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "translate: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif

   CodeWriter::translate(topel);
   DOMDocument* document = topel->getOwnerDocument();
   DOMElement* worldEl = document->getElementById(X("everything"));
   XString topvolS(worldEl->getAttribute(X("envelope")));
   DOMElement* topvolEl = document->getElementById(X(topvolS));
   XString ivoluS(topvolEl->getAttribute(X("HDDSvolu")));
   if (ivoluS.size() == 0)
   {
      G4cerr << APP_NAME << " error: top-level composition named "
             << "\"everything\" was not found in input document!"
             << G4endl;
      exit(1);
   }
   fWorldVolume = atoi(S(ivoluS));

#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
}

G4LogicalVolume* HddsG4Builder::getWorldVolume(int parallel) const
{
   int worlds = 0;
   std::map<vpair_t,G4LogicalVolume*>::const_iterator paraworld;
   for (paraworld = fLogicalVolumes.find(vpair_t(fWorldVolume,0));
        paraworld != fLogicalVolumes.end() &&
        paraworld->first.first == fWorldVolume;
        ++paraworld, ++worlds)
   {}
   if (parallel < worlds)
   {
      for (int p=0; p <= parallel; ++p, --paraworld) {}
      return paraworld->second;
   }
   return 0;
}

void HddsG4Builder::addReflections(int volume_id)
{
#ifdef LINUX_CPUTIME_PROFILING
   std::stringstream timestr;
   timestr << "addReflections: " << timer.getUserTime() << ", "
           << timer.getSystemTime() << ", " << timer.getRealTime();
   timer.resetClocks();
#endif
   std::vector<int> layer;
   std::map<vpair_t,G4LogicalVolume*>::iterator witer;
   for (witer = fLogicalVolumes.find(vpair_t(1,0));
        witer != fLogicalVolumes.end() && witer->first.first == 1;
        ++witer)
   {
      layer.push_back(witer->first.second);
   }

   std::map<int,G4FieldManager*> lastManager;
   std::map<int,G4FieldManager*>::iterator miter;
   std::map<vpair_t,G4LogicalVolume*>::iterator mine;
   for (mine = fLogicalVolumes.find(vpair_t(volume_id,0));
        mine != fLogicalVolumes.end() && mine->first.first == volume_id;
        ++mine)
   {
      G4FieldManager* fieldmgr = mine->second->GetFieldManager();
      for (int child = 0; child < (int)mine->second->GetNoDaughters(); ++child) {
         G4VPhysicalVolume* pvol = mine->second->GetDaughter(child);
         G4LogicalVolume* lvol = pvol->GetLogicalVolume();
         G4FieldManager* mgr = lvol->GetFieldManager();
         int childId = getVolumeId(lvol);
         miter = lastManager.find(childId);
         if ((miter != lastManager.end() && mgr != miter->second) ||
             (miter == lastManager.end() && mine->first.second != 0 &&
              mgr != 0 && mgr != fieldmgr))
         {
            G4cerr << "HddsG4Builder::addReflections warning - "
                   << "local magnetic field assigned to volume "
                   << lvol->GetName()
                   << " placed on a buried layer in the geometry."
                   << G4endl
                   << "This is a violation of the assumption made "
                   << "in the Geant4 tracking algorithm that fields "
                   << "at a given point are the same on all layers."
                   << G4endl
                   << "Expect inconsistent results."
                   << G4endl;
         }
         else {
            lastManager[childId] = mgr;
         }

         // Use recursion to visit every child logical volume once

         if (miter == lastManager.end()) {
            addReflections(childId);
         }

         // Any volume in the hierarchy that has a different magnetic field
         // manager than its mother volume needs to be reflected onto all
         // layers to maintain consistency of the field value at all points.

         if (mgr != 0 && mgr != fieldmgr) {
            vpair_t cpy(childId, pvol->GetCopyNo());
            std::map<vpair_t,G4VPhysicalVolume*>::iterator piter;
            piter = fPhysicalVolumes.find(cpy);
            if (piter == fPhysicalVolumes.end()) {
               G4cerr << "HddsG4Builder::addReflections error - "
                      << "physical volume " << pvol->GetName()
                      << " copy " <<  cpy.second
                      << " not found in physical volume lookup table, "
                      << "cannot continue." << G4endl;
               exit(1);
            }
            else if (piter->second->IsReplicated()) {
               fCurrentDivision[childId] = piter;
            }
            else {
               fCurrentPlacement[childId] = piter;
            }
            fCurrentMother[childId] = mine;
            std::vector<int>::iterator literi;
            for (literi = layer.begin(); literi != layer.end(); ++literi) {
               if (*literi > mine->first.second) {
#if DEBUG_REFLECTION
                  G4cout << "HddsG4Builder::addReflections info - "
                         << "reflected volume " << lvol->GetName()
                         << " copy " << cpy.second
                         << " in mother " << mine->second->GetName()
                         << " onto layer " << *literi
                         << " because child field manager " << mgr
                         << " is different from mother field manager " << fieldmgr
                         << G4endl;
#endif
                  addNewLayer(childId, *literi);
               }
            }
         }
      }
   }
#ifdef LINUX_CPUTIME_PROFILING
   timestr << " ( " << timer.getUserDelta() << " ) ";
   G4cerr << timestr.str() << G4endl;
#endif
}

int HddsG4Builder::getVolumeId(G4LogicalVolume* vol) const
{
   std::map<vpair_t,G4LogicalVolume*>::const_iterator iter;
   for (iter = fLogicalVolumes.begin(); iter != fLogicalVolumes.end(); ++iter)
   {
      if (iter->second == vol) {
         return iter->first.first;
      }
   }
   G4cerr << "HddsG4Builder::getVolumeId error - "
          << "logical volume " << vol->GetName()
          << " not found in volume lookup table, "
          << "cannot continue." << G4endl;
   exit(1);
}

G4LogicalVolume* HddsG4Builder::getVolume(const G4String volname) const
{
   std::map<vpair_t,G4LogicalVolume*>::const_iterator iter;
   for (iter = fLogicalVolumes.begin(); iter != fLogicalVolumes.end(); ++iter)
   {
      if (iter->second->GetName() == volname) {
         return iter->second;
      }
   }
   return 0;
}

const std::map<int, G4LogicalVolume*> HddsG4Builder::getSensitiveVolumes() const
{
   return fSensitiveVolumes;
}

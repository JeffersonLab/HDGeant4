//
// GlueXDetectorConstruction class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
 
#include "GlueXSensitiveDetectorDIRC.hh"
#include "GlueXSensitiveDetectorCDC.hh"
#include "GlueXSensitiveDetectorFDC.hh"
#include "GlueXSensitiveDetectorSTC.hh"
#include "GlueXSensitiveDetectorBCAL.hh"
#include "GlueXSensitiveDetectorFCAL.hh"
#include "GlueXSensitiveDetectorFCALinsert.hh"
#include "GlueXSensitiveDetectorGCAL.hh"
#include "GlueXSensitiveDetectorCCAL.hh"
#include "GlueXSensitiveDetectorECAL.hh"
#include "GlueXSensitiveDetectorFTOF.hh"
#include "GlueXSensitiveDetectorCERE.hh"
#include "GlueXSensitiveDetectorFMWPC.hh"
#include "GlueXSensitiveDetectorUPV.hh"
#include "GlueXSensitiveDetectorPSC.hh"
#include "GlueXSensitiveDetectorPS.hh"
#include "GlueXSensitiveDetectorTPOL.hh"
#include "GlueXSensitiveDetectorCTOF.hh"

#include "GlueXDetectorConstruction.hh"
#include "GlueXDetectorMessenger.hh"
#include "GlueXMagneticField.hh"
#include "HddmOutput.hh"

#include "G4Version.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4HelixMixedStepper.hh"
#include "G4ClassicalRK4.hh"

#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"

#include "G4ios.hh"

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

#include "XString.hpp"
#include "XParsers.hpp"

#include <string.h>
#include <libgen.h>
#include <errno.h>

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

#define APP_NAME "HDGeant4"

G4Mutex GlueXDetectorConstruction::fMutex = G4MUTEX_INITIALIZER;
std::list<GlueXDetectorConstruction*> GlueXDetectorConstruction::fInstance;

GlueXDetectorConstruction::GlueXDetectorConstruction(G4String hddsFile)
: fMaxStep(0),
  fUniformField(0),
  fpMagneticField(0),
  fGeometryXML(0)
{
   G4AutoLock barrier(&fMutex);
   fInstance.push_back(this);

   // Initialize the class that implements custom interactive
   // commands under the command prefix /hdgeant4/.
 
   fpDetectorMessenger = new GlueXDetectorMessenger(this);

   // Read the geometry description from a HDDS document
   // (Hall D Detector Specification) and build a Geant4
   // model of it in memory.

   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException& toCatch)
   {
      XString message(toCatch.getMessage());
      G4cerr << APP_NAME << " - error during xercesc initialization!"
             << G4endl << S(message) << G4endl;
      exit(1);
   }

   DOMDocument* document;
   if (hddsFile.size() > 0)
   {
      int size=100;
      char *saved_cwd = new char[size];
      while (getcwd(saved_cwd, size) == 0)
      {
         delete [] saved_cwd;
         saved_cwd = new char[size *= 2];
      }
      XString xmlFile = hddsFile.c_str();
      char *dirpath = new char[hddsFile.size() + 2];
      if (chdir(dirname(strcpy(dirpath, hddsFile.c_str()))) != 0) {
          G4cerr << "Error: Could not change directory to " << dirpath << G4endl;
          exit(1);
      }
      document = buildDOMDocument(xmlFile,false);
      if (chdir(saved_cwd) != 0) {
          G4cerr << "Error: Could not return to original directory!" << G4endl;
          exit(1);
      }
      delete [] saved_cwd;
      delete [] dirpath;
   }
   else if (getenv("JANA_GEOMETRY_URL"))
   {
#ifndef FORCE_HDDS_FILES_PARSING
      std::string url = getenv("JANA_GEOMETRY_URL");
      int run = HddmOutput::getRunNo();
      fGeometryXML = new HddsGeometryXML(url, run);
      last_md5_checksum = fGeometryXML->GetChecksum();
      document = fGeometryXML->getDocument();
#else
      int size=100;
      char *saved_cwd = new char[size];
      while (getcwd(saved_cwd, size) == 0)
      {
         delete [] saved_cwd;
         saved_cwd = new char[size *= 2];
      }
      XString hddsdir = getenv("HDDS_HOME");
      chdir(hddsdir.c_str());
      XString xmlFile = "main_HDDS.xml";
      document = buildDOMDocument(xmlFile,false);
      chdir(saved_cwd);
      delete [] saved_cwd;
#endif
   }
   else {
      G4cerr << APP_NAME << " - no hdds geometry file specified!"
             << " Cannot continue" << G4endl;
      exit(9);
   }
   if (document == 0)
   {
      G4cerr << APP_NAME << " - error parsing HDDS document, "
             << "cannot continue" << G4endl;
      exit(9);
   }

   try {
      (void)document->getDocumentElement();
   }
   catch (DOMException& e) {
      G4cerr << APP_NAME << " - Woops " << e.msg << G4endl;
      exit(9);
   }

   DOMElement* rootEl = document->getElementById(X("everything"));
   if (rootEl == 0)
   {
      G4cerr << APP_NAME << " - error scanning HDDS document, " << G4endl
             << "  no element named \"everything\" found, "
             << "cannot continue" << G4endl;
      exit(9);
   }

   try {
      fHddsBuilder.translate(rootEl);
   }
   catch (const DOMException &e) {
      G4cerr << APP_NAME << " - error scanning HDDS document, "
             << XString(e.getMessage()) << G4endl;
      exit(1);
   }

   XMLPlatformUtils::Terminate();
}

GlueXDetectorConstruction::
GlueXDetectorConstruction(const GlueXDetectorConstruction &src)
 : fHddsBuilder(src.fHddsBuilder)
{
   G4AutoLock barrier(&fMutex);
   fInstance.push_back(this);
   fMaxStep = src.fMaxStep;
   fUniformField = src.fUniformField;
   fpMagneticField = src.fpMagneticField; // shallow copy, sharing magfield
   fpDetectorMessenger = src.fpDetectorMessenger;  // shallow copy, sharing
}

GlueXDetectorConstruction::~GlueXDetectorConstruction()
{
   G4AutoLock barrier(&fMutex);
   fInstance.remove(this);
   delete fpDetectorMessenger;             
}

const GlueXDetectorConstruction *GlueXDetectorConstruction::GetInstance()
{
   // Generally one only needs a single instance of this object
   // per process, and this static method lets any component in the
   // application obtain the primary instance, if any. If none has
   // yet been constructed, it returns zero.

   G4AutoLock barrier(&fMutex);
   if (fInstance.size() > 0)
      return *fInstance.begin();
   return 0;
}

const HddsG4Builder *GlueXDetectorConstruction::GetBuilder()
{
   // Return a const pointer to the internal HddsG4Builder object.

   G4AutoLock barrier(&fMutex);
   if (fInstance.size() > 0)
      return &(*fInstance.begin())->fHddsBuilder;
   return 0;
}

void GlueXDetectorConstruction::SetUniformField(G4double field_T)
{
  // This method embeds the entire Hall D (spectrometer + beam line)
  // in an external uniform magnetic field. Used only for testing,
  // this is not how to simulate the GlueX spectrometer field!

  fUniformField = field_T;
  G4cout << "setting up a new field manager with a uniform field "
         << " of " << field_T << " Tesla." << G4endl;
  G4ThreeVector B(0,0,fUniformField);
  G4AffineTransform xform;
  fpMagneticField = new GlueXUniformMagField(B,tesla,xform);
}

void GlueXDetectorConstruction::SetMaxStep(G4double step_mm)
{
  fMaxStep = (step_mm > 0.)? step_mm*mm : 0;
}

G4VPhysicalVolume* GlueXDetectorConstruction::Construct()
{
   G4LogicalVolume *worldvol = fHddsBuilder.getWorldVolume();
   G4cout << "Root geometry volume " << worldvol->GetName();
   worldvol->SetName("World");
   G4cout << " configured as " << worldvol->GetName() << G4endl;
   worldvol->SetVisAttributes(new G4VisAttributes(false));
   return new G4PVPlacement(0, G4ThreeVector(), worldvol, "World", 0, 0, 0);
}

void GlueXDetectorConstruction::ConstructSDandField()
{
   G4RunManager::RMType rmtype = G4RunManager::GetRunManager()->
                                               GetRunManagerType();

   // Magnetic field objects that were created during geometry building
   // will need to be reconstructed here in order to be properly
   // configured to run in multithreaded mode.
   if (rmtype != G4RunManager::sequentialRM)
      CloneF();

   G4SDManager* SDman = G4SDManager::GetSDMpointer();
   GlueXSensitiveDetectorCDC* cdcHandler = 0;
   GlueXSensitiveDetectorFDC* fdcHandler = 0;
   GlueXSensitiveDetectorSTC* stcHandler = 0;
   GlueXSensitiveDetectorBCAL* bcalHandler = 0;
   GlueXSensitiveDetectorFCAL* fcalHandler = 0;
   GlueXSensitiveDetectorFCALinsert* fcalInsertHandler = 0;
   GlueXSensitiveDetectorGCAL* gcalHandler = 0;
   GlueXSensitiveDetectorCCAL* ccalHandler = 0;
   GlueXSensitiveDetectorECAL* ecalHandler = 0;
   GlueXSensitiveDetectorFTOF* ftofHandler = 0;
   GlueXSensitiveDetectorDIRC* dircHandler = 0;
   GlueXSensitiveDetectorCERE* cereHandler = 0;
   GlueXSensitiveDetectorFMWPC* fmwpcHandler = 0;
   GlueXSensitiveDetectorUPV* upvHandler = 0;
   GlueXSensitiveDetectorPSC* pscHandler = 0;
   GlueXSensitiveDetectorPS* psHandler = 0;
   GlueXSensitiveDetectorTPOL* tpolHandler = 0;
   GlueXSensitiveDetectorCTOF* ctofHandler = 0;

   // During geometry building, certain logical volumes were marked as
   // sensitive by adding them to a list. Now we need to go down that
   // list and pick out an explicit sensitive volume class to handle
   // each one. If any are missing, report error.
   const std::map<int, G4LogicalVolume*>
   svolMap = fHddsBuilder.getSensitiveVolumes();
   std::map<int, G4LogicalVolume*>::const_iterator iter;
   for (iter = svolMap.begin(); iter != svolMap.end(); ++iter) {
      G4String volname = iter->second->GetName();
      if (volname == "STLA" || volname == "STRA") {
         if (cdcHandler == 0) {
            cdcHandler = new GlueXSensitiveDetectorCDC("cdc");
            SDman->AddNewDetector(cdcHandler);
         }
         iter->second->SetSensitiveDetector(cdcHandler);
      }
      else if (volname == "FDA1" || volname == "FDA2" ||
               volname == "FDA3" || volname == "FDA4")
      {
         if (fdcHandler == 0) {
            fdcHandler = new GlueXSensitiveDetectorFDC("fdc");
            SDman->AddNewDetector(fdcHandler);
         }
         iter->second->SetSensitiveDetector(fdcHandler);
      }
      else if (volname == "STRC") {
         if (stcHandler == 0) {
            stcHandler = new GlueXSensitiveDetectorSTC("stc");
            SDman->AddNewDetector(stcHandler);
         }
         iter->second->SetSensitiveDetector(stcHandler);
      }
      else if (volname == "BM01" || volname == "BM02" ||
               volname == "BM03" || volname == "BM04" ||
               volname == "BM05" || volname == "BM06" ||
               volname == "BMF7" || volname == "BMF8" ||
               volname == "BMF9" || volname == "BMFA")
      {
         if (bcalHandler == 0) {
            bcalHandler = new GlueXSensitiveDetectorBCAL("bcal");
            SDman->AddNewDetector(bcalHandler);
            // Also add support beam BCL0 as a "sensitive" element
            // to catch incident particles as they enter the BCAL
            G4LogicalVolume *bcl0 = GlueXDetectorConstruction::GetBuilder()->
                                    getVolume("BCL0");
            if (bcl0 != 0) {
               bcl0->SetSensitiveDetector(bcalHandler);
            }
            else {
               G4cerr << "Warning from GlueXDetectorConstruction"
                      << "::ConstructSDandField - "
                      << "special BCal volume BCL0 not found"
                      << " in geometry definition, bcalTruthShower"
                      << " information will not be generated."
                      << G4endl;
            }
         }
         iter->second->SetSensitiveDetector(bcalHandler);
      }
      else if (volname == "LGBL" || volname == "LGLG") {
         if (fcalHandler == 0) {
            fcalHandler = new GlueXSensitiveDetectorFCAL("fcal");
            SDman->AddNewDetector(fcalHandler);
         }
         iter->second->SetSensitiveDetector(fcalHandler);
      }
     else if (volname == "LTB1") {
         if (fcalHandler == 0) {
            fcalInsertHandler = new GlueXSensitiveDetectorFCALinsert("fcalinsert");
            SDman->AddNewDetector(fcalInsertHandler);
         }
         iter->second->SetSensitiveDetector(fcalInsertHandler);
      }
      else if (volname == "GCAL") {
         if (gcalHandler == 0) {
            gcalHandler = new GlueXSensitiveDetectorGCAL("gcal");
            SDman->AddNewDetector(gcalHandler);
         }
         iter->second->SetSensitiveDetector(gcalHandler);
      }
      else if (volname == "LTBL") {
         if (ccalHandler == 0) {
            ccalHandler = new GlueXSensitiveDetectorCCAL("ccal");
            SDman->AddNewDetector(ccalHandler);
         }
         iter->second->SetSensitiveDetector(ccalHandler);
      }
      else if (volname == "XTBL") {
         if (ecalHandler == 0) {
            ecalHandler = new GlueXSensitiveDetectorECAL("ecal");
            SDman->AddNewDetector(ecalHandler);
         }
         iter->second->SetSensitiveDetector(ecalHandler);
      }
      else if (volname == "FTOC" || volname == "FTOX" || 
               volname == "FTOH" || volname == "FTOL")
      {
         if (ftofHandler == 0) {
            ftofHandler = new GlueXSensitiveDetectorFTOF("ftof");
            SDman->AddNewDetector(ftofHandler);
         }
         iter->second->SetSensitiveDetector(ftofHandler);
      }
      else if (volname == "CTOF")
      {
         if (ctofHandler == 0) {
            ctofHandler = new GlueXSensitiveDetectorCTOF("ctof");
            SDman->AddNewDetector(ctofHandler);
         }
         iter->second->SetSensitiveDetector(ctofHandler);
      }
      // radiator volume: BNNM (NN = bar number 0-47 and M is sub-bar character A-D)
      else if (volname == "PIXV" || 
               (volname[0] == 'B' && 
                10*((int)volname[1]-48)+(int)volname[2]-48 >= 0 &&
                10*((int)volname[1]-48)+(int)volname[2]-48 < 48))
      {  // this is nasty, but it works
         if (dircHandler == 0) {
            dircHandler = new GlueXSensitiveDetectorDIRC("dirc");
            SDman->AddNewDetector(dircHandler);
         }
         iter->second->SetSensitiveDetector(dircHandler);
      }
      else if (volname == "WM1N" || volname == "WM2N" || volname == "WM1S" || volname == "WM2S" ||
               volname == "FTMN" || volname == "FTMS" ||
               volname == "TM1N" || volname == "TM2N" || volname == "TM3N" ||
               volname == "TM1S" || volname == "TM2S" || volname == "TM3S" ||
               volname == "SM1N" || volname == "SM2N" || volname == "SM1S" || volname == "SM2S" ||
               volname == "OWDG" || volname(0,2) == "AG" )
      {
         if (dircHandler == 0) {
           dircHandler = new GlueXSensitiveDetectorDIRC("dirc");
           SDman->AddNewDetector(dircHandler);
         }
         iter->second->SetSensitiveDetector(dircHandler);
      }
      else if (volname == "CERW" || volname == "CPPC") {
         if (cereHandler == 0) {
            cereHandler = new GlueXSensitiveDetectorCERE("cere");
            SDman->AddNewDetector(cereHandler);
         }
         iter->second->SetSensitiveDetector(cereHandler);
      }
      else if (volname == "CPPG") {
         if (fmwpcHandler == 0) {
            fmwpcHandler = new GlueXSensitiveDetectorFMWPC("fmwpc");
            SDman->AddNewDetector(fmwpcHandler);
         }
         iter->second->SetSensitiveDetector(fmwpcHandler);
      }
      else if (volname == "UPVP" || volname == "UPVC") {
         if (upvHandler == 0) {
            upvHandler = new GlueXSensitiveDetectorUPV("upv");
            SDman->AddNewDetector(upvHandler);
         }
         iter->second->SetSensitiveDetector(upvHandler);
      }
      else if (volname == "PSCO") {
         if (pscHandler == 0) {
            pscHandler = new GlueXSensitiveDetectorPSC("psc");
            SDman->AddNewDetector(pscHandler);
         }
         iter->second->SetSensitiveDetector(pscHandler);
      }
      else if (volname(0,3) == "PTS") {
         if (tpolHandler == 0) {
            tpolHandler = new GlueXSensitiveDetectorTPOL("tpol");
            SDman->AddNewDetector(tpolHandler);
         }
         iter->second->SetSensitiveDetector(tpolHandler);
      }
      else if (volname == "PSF1" || volname == "PSF2") {
         if (psHandler == 0) {
            psHandler = new GlueXSensitiveDetectorPS("ps");
            SDman->AddNewDetector(psHandler);
         }
         iter->second->SetSensitiveDetector(psHandler);
      }
      else {
         G4cerr << "Warning from GlueXDetectorConstruction"
                << "::ConstructSDandField - "
                << "unsupported sensitive volume " << volname
                << " found in geometry definition."
                << G4endl;
      }
   }
}

void GlueXDetectorConstruction::CloneF()
{
   typedef std::map<G4FieldManager*,G4FieldManager*> FMtoFMmap;
   FMtoFMmap masterToWorker;
   G4LogicalVolumeStore* const logVolStore = G4LogicalVolumeStore::GetInstance();
   assert(logVolStore != NULL);
   for (G4LogicalVolumeStore::const_iterator iter = logVolStore->begin();
        iter != logVolStore->end();
        ++iter)
   {
      G4LogicalVolume *lvol = *iter;
      G4FieldManager* masterFM = lvol->GetFieldManager();
      G4FieldManager* clonedFM = 0;
      if (masterFM) {
         FMtoFMmap::iterator fmFound = masterToWorker.find(masterFM);
         if (fmFound == masterToWorker.end()) {

            // First time we see this FM, let's clone and remember...

            G4ChordFinder *cfinder = masterFM->GetChordFinder();
#if G4VERSION_NUMBER >= 1040
            G4VIntegrationDriver *midriver = cfinder->GetIntegrationDriver();
            double stepMinimum = 1e-2;
#else
            G4MagInt_Driver *midriver = cfinder->GetIntegrationDriver();
            double stepMinimum = midriver->GetHmin();
#endif
            G4MagIntegratorStepper *stepper = midriver->GetStepper();
            const G4Field *field = masterFM->GetDetectorField();

            G4MagneticField *field_copy;
            if (dynamic_cast<const GlueXUniformMagField*>(field)) {
               GlueXUniformMagField &orig = *(GlueXUniformMagField*)field;
               field_copy = new GlueXUniformMagField(orig);
            }
            else if (dynamic_cast<const GlueXMappedMagField*>(field)) {
               GlueXMappedMagField &orig = *(GlueXMappedMagField*)field;
               field_copy = new GlueXMappedMagField(orig);
            }
            else if (dynamic_cast<const GlueXComputedMagField*>(field)) {
               GlueXComputedMagField &orig = *(GlueXComputedMagField*)field;
               field_copy = new GlueXComputedMagField(orig);
            }
            else if (dynamic_cast<const G4UniformMagField*>(field)) {
               G4UniformMagField &orig = *(G4UniformMagField*)field;
               field_copy = new G4UniformMagField(orig);
            }
            else {
               G4cerr << "GlueXDetectorConstruction::CloneF error - "
                      << "unknown G4MagneticField class found "
                      << "attached to detector volume " << lvol->GetName()
                      << ", cannot continue!" << G4endl;
               exit(1);
            }
            G4Mag_UsualEqRhs *eqn_copy = new G4Mag_UsualEqRhs(field_copy);
            G4MagIntegratorStepper *stepper_copy;
            if (dynamic_cast<const G4ExactHelixStepper*>(stepper)) {
               stepper_copy = new G4ExactHelixStepper(eqn_copy);
            }
            else if (dynamic_cast<const G4ClassicalRK4*>(stepper)) {
               stepper_copy = new G4ClassicalRK4(eqn_copy);
            }
            else if (dynamic_cast<const G4HelixMixedStepper*>(stepper)) {
               stepper_copy = new G4HelixMixedStepper(eqn_copy);
            }
            else {
               G4cerr << "GlueXDetectorConstruction::CloneF error - "
                      << "unknown G4MagIntegratorStepper class found "
                      << "attached to detector volume " << lvol->GetName()
                      << ", cannot continue!" << G4endl;
               exit(1);
            }
            G4ChordFinder *cfinder_copy = new G4ChordFinder(field_copy,
                                                            stepMinimum,
                                                            stepper_copy);
            clonedFM = new G4FieldManager(field_copy, cfinder_copy);
            masterToWorker[masterFM] = clonedFM;
         }
         else {

            // We have already seen this FM attached to a different 
            // LogicalVolume, let's re-use the previous clone
 
            clonedFM = (*fmFound).second;
         }
         lvol->SetFieldManager(clonedFM, false);
      }
   }
}

int GlueXDetectorConstruction::GetParallelWorldCount() const
{
   for (int worlds=0; ; ++worlds) {
      if (fHddsBuilder.getWorldVolume(worlds) == 0) {
         return --worlds;
      }
   }
   return 0;
}

G4String GlueXDetectorConstruction::GetParallelWorldName(int paraIndex) const
{
   std::stringstream worldStr;
   worldStr << "ParallelWorld" << paraIndex;
   return worldStr.str();
}

G4LogicalVolume*
 GlueXDetectorConstruction::GetParallelWorldVolume(int paraIndex) const
{
   return fHddsBuilder.getWorldVolume(paraIndex);
}

G4ThreeVector
GlueXDetectorConstruction::GetMagneticField(G4ThreeVector pos, double unit)
const
{
   // Utility function for use by other simulation components,
   // returns the magnetic field at an arbitrary location in
   // the geometry. If geometry has not yet been constructed
   // the value returned is always zero.

   G4TransportationManager *tmanager;
   G4Navigator *navigator;
   G4VPhysicalVolume *pvol;
   G4LogicalVolume *lvol;
   G4FieldManager *fieldmgr;
   const G4Field *field;
   if ((tmanager = G4TransportationManager::GetTransportationManager()) &&
       (navigator = tmanager->GetNavigatorForTracking()) &&
       (pvol = navigator->LocateGlobalPointAndSetup(pos)) &&
       (lvol = pvol->GetLogicalVolume()) &&
       (fieldmgr = lvol->GetFieldManager()) &&
       (field = fieldmgr->GetDetectorField()) )
   {
      double Bfield[3];
      double xglob[4] = {pos[0], pos[1], pos[2], 0};
      field->GetFieldValue(xglob, Bfield);
      G4ThreeVector B(Bfield[0], Bfield[1], Bfield[2]);
      return B / unit;
   }
   return G4ThreeVector();
}



GlueXParallelWorld::GlueXParallelWorld(const GlueXParallelWorld &src)
 : G4VUserParallelWorld(src.fWorldName), fTopVolume(src.fTopVolume)
{ }

GlueXParallelWorld::~GlueXParallelWorld() { }

void GlueXParallelWorld::Construct()
{
   G4VPhysicalVolume* ghostWorld = GetWorld();
   G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
   for (int child = fTopVolume->GetNoDaughters() - 1; child >= 0; --child) {
      G4VPhysicalVolume* pvol = fTopVolume->GetDaughter(child);
      fTopVolume->RemoveDaughter(pvol);
      worldLogical->AddDaughter(pvol);
   }
   if (worldLogical->GetNoDaughters()) {
      G4cout << "Additional geometry layer configured as "
             << ghostWorld->GetName() << G4endl;
   }
   else {
      G4cout << "Additional geometry layer configured as "
             << ghostWorld->GetName() << " is EMPTY!" << G4endl
             << "This is probably due to a geometry bug -- "
             << "PLEASE REPORT THIS TO THE AUTHORS" << G4endl;
   }
}

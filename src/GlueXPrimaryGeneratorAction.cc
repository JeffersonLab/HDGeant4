//
// class implementation for GlueXPrimaryGeneratorAction
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXPrimaryGenerator.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserOptions.hh"
#include "G4OpticalPhoton.hh"
#include "HddmOutput.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <JANA/Geometry/JGeometryManager.h>
#include <JANA/Calibrations/JCalibrationManager.h>



typedef GlueXPrimaryGeneratorAction::source_type_t source_type_t;
typedef GlueXPrimaryGeneratorAction::single_particle_gun_t particle_gun_t;

source_type_t GlueXPrimaryGeneratorAction::fSourceType = SOURCE_TYPE_NONE;

std::ifstream *GlueXPrimaryGeneratorAction::fHDDMinfile = 0;
hddm_s::istream *GlueXPrimaryGeneratorAction::fHDDMistream = 0;

G4ParticleTable *GlueXPrimaryGeneratorAction::fParticleTable = 0;

particle_gun_t GlueXPrimaryGeneratorAction::fGunParticle;

double GlueXPrimaryGeneratorAction::fBeamBackgroundRate = 0;
double GlueXPrimaryGeneratorAction::fBeamBackgroundGateStart = 0;
double GlueXPrimaryGeneratorAction::fBeamBackgroundGateStop = 0;
double GlueXPrimaryGeneratorAction::fL1triggerTimeSigma = 10 * ns;
double GlueXPrimaryGeneratorAction::fRFreferencePlaneZ = 65 * cm;
double GlueXPrimaryGeneratorAction::fTargetCenterZ = 65 * cm;
double GlueXPrimaryGeneratorAction::fTargetLength = 29.9746 * cm;
double GlueXPrimaryGeneratorAction::fBeamEndpointEnergy = 12 * GeV;
double GlueXPrimaryGeneratorAction::fBeamPeakEnergy = 9 * GeV;

G4Mutex GlueXPrimaryGeneratorAction::fMutex = G4MUTEX_INITIALIZER;
std::list<GlueXPrimaryGeneratorAction*> GlueXPrimaryGeneratorAction::fInstance;

double GlueXPrimaryGeneratorAction::DIRC_LUT_X[48] = {0};
double GlueXPrimaryGeneratorAction::DIRC_BAR_Y[48] = {0};
double GlueXPrimaryGeneratorAction::DIRC_LUT_Z = 0;
double GlueXPrimaryGeneratorAction::DIRC_QZBL_DY = 0;
double GlueXPrimaryGeneratorAction::DIRC_QZBL_DZ = 0;
double GlueXPrimaryGeneratorAction::DIRC_OWDG_DZ = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCS_FDTH_X = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCS_FDTH_Z = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCN_FDTH_X = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCN_FDTH_Z = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCN_FDTH1_Y = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCN_FDTH2_Y = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCN_FDTH3_Y = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCS_FDTH4_Y = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCS_FDTH5_Y = 0;
double GlueXPrimaryGeneratorAction::DIRC_LED_OBCS_FDTH6_Y = 0;

//--------------------------------------------
// GlueXPrimaryGeneratorAction (constructor)
//--------------------------------------------

GlueXPrimaryGeneratorAction::GlueXPrimaryGeneratorAction()
 : fParticleGun(0),
   fCobremsGeneration(0),
   fPhotonBeamGenerator(0),
   fPrimaryGenerator(0),
   fBeamvertex_activated(0)
{
   G4AutoLock barrier(&fMutex);
   fInstance.push_back(this);

   // Initializaton is driven by the control.in file, which
   // gets read and parsed only once, by the first constructor.

   fParticleGun = new GlueXParticleGun();

   if (fSourceType == SOURCE_TYPE_HDDM) {
      clone_photon_beam_generator();
      fPrimaryGenerator = new GlueXPrimaryGenerator(fHDDMistream);
      return;
   }
   else if (fSourceType == SOURCE_TYPE_COBREMS_GEN) {
      clone_photon_beam_generator();
      return;
   }
   else if (fSourceType == SOURCE_TYPE_PARTICLE_GUN) {
      fParticleGun->SetParticleDefinition(fGunParticle.partDef);
      return;
   }

   fParticleTable = G4ParticleTable::GetParticleTable();
   
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXPrimaryGeneratorAction constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   // get positions for LUT from XML geometry
   std::map<int, int> dirclutpars;
   std::map<int, int> dircledpars;

   if (user_opts->Find("DIRCLUT", dirclutpars)) {
      
      int runno = HddmOutput::getRunNo();
      extern JApplication *japp;
      if (japp == 0) {
         G4cerr << "Error in GlueXPrimaryGeneratorAction constructor - "
           << "jana global DApplication object not set, "
           << "cannot continue." << G4endl;
         exit(-1);
      }
      JGeometry *jgeom = japp->GetService<JGeometryManager>()->GetJGeometry(runno);
      if (japp == 0) {   // dummy
         jgeom = 0;
         G4cout << "DIRC: ALL parameters loaded from ccdb" << G4endl;
      }
      
      std::vector<double> DIRC;
      std::vector<double> DRCC;
      std::vector<double> DCML00_XYZ;
      std::vector<double> DCML01_XYZ;
      std::vector<double> DCML10_XYZ;
      std::vector<double> DCML11_XYZ;
      std::vector<double> WNGL00_XYZ;
      std::vector<double> WNGL01_XYZ;
      std::vector<double> WNGL10_XYZ;
      std::vector<double> WNGL11_XYZ;
      std::vector<double> OWDG_XYZ;
      jgeom->Get("//section/composition/posXYZ[@volume='DIRC']/@X_Y_Z", DIRC);
      jgeom->Get("//composition[@name='DIRC']/posXYZ[@volume='DRCC']/@X_Y_Z", DRCC);
      jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='DCML00']/@X_Y_Z", DCML00_XYZ);
      jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='DCML01']/@X_Y_Z", DCML01_XYZ);
      jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='DCML10']/@X_Y_Z", DCML10_XYZ);
      jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='DCML11']/@X_Y_Z", DCML11_XYZ);
      jgeom->Get("//composition[@name='DCML00']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL00_XYZ);
      jgeom->Get("//composition[@name='DCML01']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL01_XYZ);
      jgeom->Get("//composition[@name='DCML10']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL10_XYZ);
      jgeom->Get("//composition[@name='DCML11']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL11_XYZ);
      jgeom->Get("//composition[@name='DCML11']/posXYZ[@volume='WNGL']/@X_Y_Z", WNGL11_XYZ);
      jgeom->Get("//trd[@name='OWDG']/@Xmp_Ymp_Z", OWDG_XYZ);
      DIRC_LUT_Z = (DIRC[2] + DRCC[2] + DCML01_XYZ[2] + 0.8625) * cm;
      DIRC_QZBL_DY = 3.5 * cm;   // nominal width to generate LUT
      DIRC_QZBL_DZ = 1.725 * cm; // nominal thickness to generate LUT
      DIRC_OWDG_DZ = OWDG_XYZ[4];

      // set array of bar positions
      for (int i=0; i<48; i++) {
         std::vector<double> DCBR_XYZ;
         if (i<12) {
            std::stringstream geomDCML10;
            geomDCML10 << "//composition[@name='DCML10']/posXYZ[@volume='DCBR" 
                  << std::setfill('0') << std::setw(2) << i << "']/@X_Y_Z"; 
            jgeom->Get(geomDCML10.str(), DCBR_XYZ);
            DIRC_BAR_Y[i] = (DCML10_XYZ[1] - DCBR_XYZ[1]) * cm;
            DIRC_LUT_X[i] = (DIRC[0] + DRCC[0] + DCML10_XYZ[0] - WNGL10_XYZ[0] + DIRC_OWDG_DZ) * cm;
         }
         else if (i<24) {
            std::stringstream geomDCML11;
            geomDCML11 << "//composition[@name='DCML11']/posXYZ[@volume='DCBR" 
                  << std::setfill('0') << std::setw(2) << i << "']/@X_Y_Z"; 
            jgeom->Get(geomDCML11.str(), DCBR_XYZ);
            DIRC_BAR_Y[i] = (DCML11_XYZ[1] - DCBR_XYZ[1]) * cm;
            DIRC_LUT_X[i] = (DIRC[0] + DRCC[0] + DCML11_XYZ[0] - WNGL11_XYZ[0] + DIRC_OWDG_DZ) * cm;
         }
         else if (i<36) {
            std::stringstream geomDCML01;
            geomDCML01 << "//composition[@name='DCML01']/posXYZ[@volume='DCBR" 
                  << std::setfill('0') << std::setw(2) << i << "']/@X_Y_Z"; 
            jgeom->Get(geomDCML01.str(), DCBR_XYZ);
            DIRC_BAR_Y[i] = (DCML01_XYZ[1] + DCBR_XYZ[1]) * cm;
            DIRC_LUT_X[i] = (DIRC[0] + DRCC[0] + DCML01_XYZ[0] + WNGL01_XYZ[0] - DIRC_OWDG_DZ) * cm;
         }
         else if (i<48) {
            std::stringstream geomDCML00;
            geomDCML00 << "//composition[@name='DCML00']/posXYZ[@volume='DCBR" 
                  << std::setfill('0') << std::setw(2) << i << "']/@X_Y_Z"; 
            jgeom->Get(geomDCML00.str(), DCBR_XYZ);
            DIRC_BAR_Y[i] = (DCML00_XYZ[1] + DCBR_XYZ[1]) * cm;
            DIRC_LUT_X[i] = (DIRC[0] + DRCC[0] + DCML00_XYZ[0] + WNGL00_XYZ[0] - DIRC_OWDG_DZ) * cm;
         }            
      }
   }

   if (user_opts->Find("DIRCLED", dircledpars)) {
      int runno = HddmOutput::getRunNo();
      extern JApplication *japp;
      if (japp == 0) {
         G4cerr << "Error in GlueXPrimaryGeneratorAction constructor - "
           << "jana global DApplication object not set, "
           << "cannot continue." << G4endl;
         exit(-1);
      }
      JGeometry *jgeom = japp->GetService<JGeometryManager>()->GetJGeometry(runno);
      if (japp == 0) {   // dummy
         jgeom = 0;
         G4cout << "DIRC: ALL parameters loaded from ccdb" << G4endl;
      }
      std::vector<double> DIRC;
      std::vector<double> DRCC;
      std::vector<double> OBCS_XYZ;
      std::vector<double> OBCN_XYZ;
      std::vector<double> MRAN_XYZ;
      std::vector<double> MRAS_XYZ;
      std::vector<double> WM1N_XYZ;
      std::vector<double> WM1S_XYZ;
      std::vector<double> WM1N_BOX_XYZ;
      std::vector<double> WM1S_BOX_XYZ;

      jgeom->Get("//section/composition/posXYZ[@volume='DIRC']/@X_Y_Z", DIRC);
      jgeom->Get("//composition[@name='DIRC']/posXYZ[@volume='DRCC']/@X_Y_Z", DRCC);
      jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='OBCN']/@X_Y_Z", OBCN_XYZ);
      jgeom->Get("//composition[@name='DRCC']/posXYZ[@volume='OBCS']/@X_Y_Z", OBCS_XYZ);
      jgeom->Get("//composition[@name='OBCN']/posXYZ[@volume='MRAN']/@X_Y_Z", MRAN_XYZ);
      jgeom->Get("//composition[@name='OBCS']/posXYZ[@volume='MRAS']/@X_Y_Z", MRAS_XYZ);
      jgeom->Get("//composition[@name='MRAN']/posXYZ[@volume='WM1N']/@X_Y_Z", WM1N_XYZ);
      jgeom->Get("//composition[@name='MRAS']/posXYZ[@volume='WM1S']/@X_Y_Z", WM1S_XYZ);
      jgeom->Get("//box[@name='WM1N']/@X_Y_Z", WM1N_BOX_XYZ);
      jgeom->Get("//box[@name='WM1S']/@X_Y_Z", WM1S_BOX_XYZ);

      DIRC_LED_OBCN_FDTH_X  = (DIRC[0] + DRCC[0] + OBCN_XYZ[0] + MRAN_XYZ[0] +
                               WM1N_XYZ[0] + WM1N_BOX_XYZ[0]/2. + 1.27) * cm;
      DIRC_LED_OBCS_FDTH_X  = (DIRC[0] + DRCC[0] + OBCS_XYZ[0] + MRAS_XYZ[0] +
                               WM1S_XYZ[0] - WM1S_BOX_XYZ[0]/2. - 1.27) * cm;

      DIRC_LED_OBCN_FDTH_Z  = (DIRC[2] + DRCC[2] + OBCN_XYZ[2] + MRAN_XYZ[2] +
                               WM1N_XYZ[2] - WM1N_BOX_XYZ[2]/2. ) * cm;
      DIRC_LED_OBCS_FDTH_Z  = (DIRC[2] + DRCC[2] + OBCS_XYZ[2] + MRAS_XYZ[2] +
                               WM1S_XYZ[2] - WM1S_BOX_XYZ[2]/2. ) * cm;

      DIRC_LED_OBCN_FDTH1_Y = (DIRC[1] + DRCC[1] + OBCN_XYZ[1] + MRAN_XYZ[1] +
                               WM1N_XYZ[1] - WM1N_BOX_XYZ[1]/2. + 17.235932) * cm;
      DIRC_LED_OBCN_FDTH2_Y = (DIRC[1] + DRCC[1] + OBCN_XYZ[1] + MRAN_XYZ[1] +
                               WM1N_XYZ[1] - WM1N_BOX_XYZ[1]/2. + 17.235932 + 31.800038 ) * cm;
      DIRC_LED_OBCN_FDTH3_Y = (DIRC[1] + DRCC[1] + OBCN_XYZ[1] + MRAN_XYZ[1] +
                               WM1N_XYZ[1] - WM1N_BOX_XYZ[1]/2. + 17.235932 + 31.800038 * 2. ) * cm;
      DIRC_LED_OBCS_FDTH4_Y = (DIRC[1] + DRCC[1] + OBCS_XYZ[1] + MRAS_XYZ[1] +
                               WM1S_XYZ[1] + WM1S_BOX_XYZ[1]/2. - 17.235932) * cm;
      DIRC_LED_OBCS_FDTH5_Y = (DIRC[1] + DRCC[1] + OBCS_XYZ[1] + MRAS_XYZ[1] +
                               WM1S_XYZ[1] + WM1S_BOX_XYZ[1]/2. - 17.235932 - 31.800038 ) * cm;
      DIRC_LED_OBCS_FDTH6_Y = (DIRC[1] + DRCC[1] + OBCS_XYZ[1] + MRAS_XYZ[1] +
                               WM1S_XYZ[1] + WM1S_BOX_XYZ[1]/2. - 17.235932 - 31.800038 * 2. ) * cm;

   }

   std::map<int,std::string> infile;
   std::map<int,double> beampars;
   std::map<int,double> genbeampars;
   std::map<int,double> kinepars;

   // Three event source options are supported:
   // 1) external generator, hddm input stream source
   // 2) internal coherent bremsstrahlung beam generator
   // 3) internal particle gun generator
 
   if (user_opts->Find("INFILE", infile) ||
       user_opts->Find("INFI", infile))
   {
      fHDDMinfile = new std::ifstream(infile[1].c_str());
      if (!fHDDMinfile->is_open()) {
         G4cerr << "GlueXPrimaryGeneratorAction error: "
                << "Unable to open HDDM input file: " << infile[1]
                << G4endl;
         exit(-1);
      }
      fHDDMistream = new hddm_s::istream(*fHDDMinfile);
      G4cout << "Opened input file: " << infile[1] << G4endl;
      std::map<int,int> skippars;
      if (user_opts->Find("SKIP", skippars))
      {
         if (skippars[1] > 0) 
         {
            fHDDMistream->skip(skippars[1]);
            G4cout << "skipped first " << skippars[1] << " input events." << G4endl;
         }
      }
      fPrimaryGenerator = new GlueXPrimaryGenerator(fHDDMistream);
      fSourceType = SOURCE_TYPE_HDDM;
   }

   else if (user_opts->Find("BEAM", beampars) ||
            user_opts->Find("GENBEAM", genbeampars))
   {
      fSourceType = SOURCE_TYPE_COBREMS_GEN;
   }

   else if (user_opts->Find("DIRCLUT", dirclutpars))
   {
      fGunParticle.geantType = 0;
      fGunParticle.pdgType = 999999;
      fGunParticle.partDef = fParticleTable->FindParticle("opticalphoton");
      fGunParticle.deltaR = 0;
      fGunParticle.deltaZ = 0;
      fGunParticle.mom = 3.18 * eV;

      fGunParticle.deltaMom = 0;
      fGunParticle.deltaTheta = 0;
      fGunParticle.deltaPhi = 0;
      fParticleGun->SetParticleDefinition(fGunParticle.partDef);
 
      fSourceType = SOURCE_TYPE_PARTICLE_GUN;
   }

   else if (user_opts->Find("DIRCLED", dircledpars))
   {
      fGunParticle.geantType = 0;
      fGunParticle.pdgType = 999999;
      fGunParticle.partDef = fParticleTable->FindParticle("opticalphoton");
      fGunParticle.deltaR = 0;
      fGunParticle.deltaZ = 0;
      fGunParticle.mom = 3.0613 * eV;

      fGunParticle.deltaMom = 0;
      fGunParticle.deltaTheta = 0;
      fGunParticle.deltaPhi = 0;
      fParticleGun->SetParticleDefinition(fGunParticle.partDef);

      fSourceType = SOURCE_TYPE_PARTICLE_GUN;
   }

   else if (user_opts->Find("KINE", kinepars))
   {
      if (kinepars[1] == 1000 || kinepars[1] == 148 || kinepars[1] == 48) {
         fGunParticle.geantType = 0;
         fGunParticle.pdgType = 999999;
         fGunParticle.partDef = fParticleTable->FindParticle("geantino");
      }
      else if (kinepars[1] == 1001) {
         fGunParticle.geantType = 0;
         fGunParticle.pdgType = 999999;
         fGunParticle.partDef = fParticleTable->FindParticle("chargedgeantino");
      }
      else if (kinepars[1] == 1050) {
         fGunParticle.geantType = 0;
         fGunParticle.pdgType = 999999;
         fGunParticle.partDef = fParticleTable->FindParticle("opticalphoton");
      }            
      else {
         if (kinepars[1] > 100)
            fGunParticle.geantType = kinepars[1] - 100;
         else
            fGunParticle.geantType = kinepars[1];
         fGunParticle.pdgType = ConvertGeant3ToPdg(fGunParticle.geantType);
         fGunParticle.partDef = fParticleTable->FindParticle(fGunParticle.pdgType);
      }
      if (fGunParticle.partDef == 0) {   
         G4cerr << "GlueXPrimaryGeneratorAction constructor error - "
                << "Unknown GEANT particle type: " << kinepars[1] 
                << " was specified in the control.in file." << G4endl;
         exit(-1);
      }
      fParticleGun->SetParticleDefinition(fGunParticle.partDef);

      double x(0), y(0), z(65 * cm);
      std::map<int,double> scappars;
      if (user_opts->Find("SCAP", scappars)) {
         x = scappars[1] * cm;
         y = scappars[2] * cm;
         z = scappars[3] * cm;
      }
      fGunParticle.pos.set(x,y,z);
      std::map<int,double> tgtwidthpars;
      if (user_opts->Find("tgtwidth", tgtwidthpars)) {
         fGunParticle.deltaR = tgtwidthpars[1] * cm;
         fGunParticle.deltaZ = tgtwidthpars[2] * cm;
      }
      else {
         fGunParticle.deltaR = 0;
         fGunParticle.deltaZ = 0;
      }
      fGunParticle.plogOption = 0;
      std::map<int,int> plogpars;
      if (user_opts->Find("PLOG", plogpars)) {
         fGunParticle.plogOption = plogpars[1];
      }
      fGunParticle.tlogOption = 0;
      std::map<int,int> tlogpars;
      if (user_opts->Find("TLOG", tlogpars)) {
         fGunParticle.tlogOption = tlogpars[1];
      }

      fGunParticle.mom = kinepars[2] * GeV;
      if (kinepars[1] > 100) {
         fGunParticle.theta = kinepars[3] * degree;
         fGunParticle.phi = kinepars[4] * degree;
         fGunParticle.deltaMom = kinepars[5] * GeV;
         fGunParticle.deltaTheta = kinepars[6] * degree;
         fGunParticle.deltaPhi = kinepars[7] * degree;
      }
      else {
         fGunParticle.deltaMom = 0;
         fGunParticle.theta = 90 * degree;
         fGunParticle.deltaTheta = 180 * degree;
         fGunParticle.phi = 0;
         fGunParticle.deltaPhi = 360 * degree;
      }
      fSourceType = SOURCE_TYPE_PARTICLE_GUN;
   }

   if (user_opts->Find("BEAM", beampars)) {
      double beamE0 = beampars[1] * GeV;
      double beamEpeak = beampars[2] * GeV;
      double beamEmin = ((beampars[3] > 0)? beampars[3] : 0.120) * GeV;
      double radColDist = ((beampars[4] > 0)? beampars[4] : 76.) * m;
      double colDiam = ((beampars[5] > 0)? beampars[5] : 0.0034) * m;
      double beamEmit = ((beampars[6] > 0)? beampars[6] : 2.5e-9) * m;
      double radThick = ((beampars[7] > 0)? beampars[7] : 20e-6) * m;
      double spotRMS = ((beampars[8] > 0)? beampars[8] : 5e-4) * m;
      double spotX = ((beampars[9] != 0)? beampars[9] : 0) * m;
      double spotY = ((beampars[10] != 0)? beampars[10] : 0) * m;

      if (beamE0 == 0) {
         G4cerr << "GlueXPrimaryGeneratorAction error: "
                << "BEAM card specified in control.in but required values "
                << "Emax and/or Epeak are missing, cannot continue."
                << G4endl;
         exit(-1);
      }

      fBeamEndpointEnergy = beamE0;
      fBeamPeakEnergy = beamEpeak;

      // CobremsGeneration has its own standard units that it uses:
      //  length : m
      //  angles : radians
      //  energy : GeV
      //  time   : s
      //  current: microAmps
 
      fCobremsGeneration = new CobremsGeneration(beamE0/GeV, beamEpeak/GeV);
      fCobremsGeneration->setPhotonEnergyMin(beamEmin/GeV);
      fCobremsGeneration->setCollimatorDistance(radColDist/m);
      fCobremsGeneration->setCollimatorDiameter(colDiam/m);
      fCobremsGeneration->setBeamEmittance(beamEmit/(m*radian));
      fCobremsGeneration->setTargetThickness(radThick/m);
      fCobremsGeneration->setCollimatorSpotrms(spotRMS/m);
      fPhotonBeamGenerator = new GlueXPhotonBeamGenerator(fCobremsGeneration);
      fPhotonBeamGenerator->setBeamOffset(spotX, spotY);
   }
   else {  // look up beam properties in the ccdb based on run number
      int runno = HddmOutput::getRunNo();
      extern JApplication *japp;
      if (japp == 0) {
         G4cerr << "Error in GlueXPrimaryGeneratorAction constructor - "
                << "jana global DApplication object not set, "
                << "cannot continue." << G4endl;
         exit(-1);
      }
      std::map<string, float> beam_parms;
      JCalibration *jcalib = japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
      jcalib->Get("PHOTON_BEAM/endpoint_energy", beam_parms);
      fBeamEndpointEnergy = beam_parms.at("PHOTON_BEAM_ENDPOINT_ENERGY")*GeV;
      std::map<string, float> coherent_parms;
      jcalib->Get("PHOTON_BEAM/coherent_energy", coherent_parms);
      fBeamPeakEnergy = coherent_parms.at("cohedge_energy")*GeV;

      // CobremsGeneration has its own standard units that it uses:
      //  length : m
      //  angles : radians
      //  energy : GeV
      //  time   : s
      //  current: microAmps
 
      fCobremsGeneration = new CobremsGeneration(fBeamEndpointEnergy/GeV,
                                                 fBeamPeakEnergy/GeV);
      fCobremsGeneration->setPhotonEnergyMin(fBeamEndpointEnergy/GeV * 0.25);
      fCobremsGeneration->setCollimatorDistance(76.);
      fCobremsGeneration->setCollimatorDiameter(0.005);
      fCobremsGeneration->setBeamEmittance(2e-9);
      fCobremsGeneration->setTargetThickness(50e-6);
      fCobremsGeneration->setCollimatorSpotrms(0.5e-3);
      fPhotonBeamGenerator = new GlueXPhotonBeamGenerator(fCobremsGeneration);
      std::map<string, float> beam_spot_params;
      jcalib->Get("/PHOTON_BEAM/beam_spot", beam_spot_params);
      double spotX = beam_spot_params.at("x") * cm;
      double spotY = beam_spot_params.at("y") * cm;
      fPhotonBeamGenerator->setBeamOffset(spotX, spotY);
      std::cout << "no beam card found, looking up beam properties in ccdb, "
                << "run number " << runno 
                << ": endpoint=" << fBeamEndpointEnergy/GeV << "GeV,"
                << " peak=" << fBeamPeakEnergy/GeV << "GeV" << std::endl;
   }

   std::map<int, double> bgratepars;
   std::map<int, double> bggatepars;
   if (user_opts->Find("BGRATE", bgratepars) &&
       user_opts->Find("BGGATE", bggatepars))
   {
      fBeamBackgroundRate = bgratepars[1] * 1/ns;
      fBeamBackgroundGateStart = bggatepars[1] * ns;
      fBeamBackgroundGateStop = bggatepars[2] * ns;
      if (fBeamBackgroundRate > 0 &&
          fBeamBackgroundGateStart >= fBeamBackgroundGateStop)
      {
         G4cerr << "GlueXPrimaryGeneratorAction error: "
                << "BGRATE is non-zero, but the time window specified "
                << "in BGGATE is invalid."
                << G4endl;
         exit(-1);
      }
   }

   std::map<int,double> trefsigma;
   if (user_opts->Find("trefsigma", trefsigma)) {
      fL1triggerTimeSigma = trefsigma[1] * ns;
   }
   else {
      fL1triggerTimeSigma = 10 * ns;
   }
}

GlueXPrimaryGeneratorAction::GlueXPrimaryGeneratorAction(const
                             GlueXPrimaryGeneratorAction &src)
 : GlueXPrimaryGeneratorAction()
{
   fBeamvertex = src.fBeamvertex;
   fBeamvertex_activated = src.fBeamvertex_activated;
}

GlueXPrimaryGeneratorAction &GlueXPrimaryGeneratorAction::operator=(const
                             GlueXPrimaryGeneratorAction &src)
{
   fParticleGun = new GlueXParticleGun();
   fCobremsGeneration = 0;
   fPhotonBeamGenerator = 0;
   fPrimaryGenerator = 0;

   if (fSourceType == SOURCE_TYPE_HDDM) {
      fPrimaryGenerator = new GlueXPrimaryGenerator(fHDDMistream);
   }
   else if (fSourceType == SOURCE_TYPE_PARTICLE_GUN) {
      fParticleGun->SetParticleDefinition(fGunParticle.partDef);
   }
   if (src.fPhotonBeamGenerator) {
      clone_photon_beam_generator();
   }
   fBeamvertex = src.fBeamvertex;
   fBeamvertex_activated = src.fBeamvertex_activated;
   return *this;
}

//--------------------------------------------
// ~GlueXPrimaryGeneratorAction (destructor)
//--------------------------------------------

GlueXPrimaryGeneratorAction::~GlueXPrimaryGeneratorAction()
{
   G4AutoLock barrier(&fMutex);
   fInstance.remove(this);
   if (fPrimaryGenerator)
      delete fPrimaryGenerator;
   if (fCobremsGeneration)
      delete fCobremsGeneration;
   if (fPhotonBeamGenerator)
      delete fPhotonBeamGenerator;
   delete fParticleGun;
   if (fInstance.size() == 0) {
      if (fHDDMistream)
         delete fHDDMistream;
      if (fHDDMinfile)
         delete fHDDMinfile;
   }
}

const CobremsGeneration *GlueXPrimaryGeneratorAction::GetCobremsGeneration()
{
   return fCobremsGeneration;
}

//--------------------------------------------
// GeneratePrimaries
//--------------------------------------------

void GlueXPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   switch (fSourceType) {
      case SOURCE_TYPE_HDDM:
         GeneratePrimariesHDDM(anEvent);
         break;
      case SOURCE_TYPE_COBREMS_GEN:
         GeneratePrimariesCobrems(anEvent);
         break;
      case SOURCE_TYPE_PARTICLE_GUN:
         GeneratePrimariesParticleGun(anEvent);
         break;
      default:
         G4cout << "No event source selected, cannot continue!" << G4endl;
         exit(-1);
   }
}   

//--------------------------------------------
// GeneratePrimariesParticleGun
//--------------------------------------------

void GlueXPrimaryGeneratorAction::GeneratePrimariesParticleGun(G4Event* anEvent)
{
   // GlueXUserEventInformation constructor can seed the random number
   // generator for this event, so this must happen here at the top.
   GlueXUserEventInformation *event_info = new GlueXUserEventInformation();
   anEvent->SetUserInformation(event_info);
   GlueXUserEventInformation::setWriteNoHitEvents(1);

   // Unbelievably, GEANT4's G4ParticleGun class insists on printing
   // a message whenever the momentum or energy is changed, unless
   // the other is 0. Here, we reset the particle gun energy using 
   // our own derived class. (Sheesh!!)
   fParticleGun->Reset();

   //   std::cout<<"GlueXPrimaryGeneratorAction:: GENERATE PRIMARIES PARTICLE GUN"<<std::endl;

   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   std::map<int,int> dirclutpars; 
   std::map<int,int> dircledpars; 

   // place and smear the particle gun origin
   G4ThreeVector pos(fGunParticle.pos);
   if (fGunParticle.deltaR > 0) {
      double dx, dy;
      while (true) {
         double rx = G4UniformRand() - 0.5;
         double ry = G4UniformRand() - 0.5;
         if (rx*rx + ry*ry <= 0.25) {
            dx = rx * 2 * fGunParticle.deltaR;
            dy = ry * 2 * fGunParticle.deltaR;
            break;
         }
      }
      pos += G4ThreeVector(dx, dy, 0);
   }
   if (fGunParticle.deltaZ > 0) {
      double dz = (G4UniformRand() - 0.5) * fGunParticle.deltaZ;
      pos += G4ThreeVector(0, 0, dz);
   }
   fParticleGun->SetParticlePosition(pos);

   // Assign and optionally smear the particle momentum
   double p = fGunParticle.mom;
   double thetap = fGunParticle.theta;
   double phip = fGunParticle.phi;
   if (fGunParticle.deltaMom > 0) {
      if (fGunParticle.plogOption) {
         double pmin = p - fGunParticle.deltaMom / 2;
         double pmax = p + fGunParticle.deltaMom / 2;
         pmin = (pmin > 0)? pmin : 1e-6*GeV;
         p = pmin * pow(pmax/pmin, G4UniformRand());
      }
      else {
         p += (G4UniformRand() - 0.5) * fGunParticle.deltaMom;
      }
   }
   if (fGunParticle.deltaTheta > 0) {
      if (fGunParticle.plogOption) {
         double thetamin = thetap - fGunParticle.deltaTheta / 2;
         double thetamax = thetap + fGunParticle.deltaTheta / 2;
         thetamin = (thetamin > 0)? thetamin : 1e-6*degree;
         thetap = thetamin * pow(thetamax/thetamin, G4UniformRand());
      }
      else {
         thetap += (G4UniformRand() - 0.5) * fGunParticle.deltaTheta;
      }
   }
   if (fGunParticle.deltaPhi > 0)
      phip += (G4UniformRand() - 0.5) * fGunParticle.deltaPhi;

   // Special case of Cherenkov photon gun for DIRC Look Up Tables (LUT)
   if (user_opts->Find("DIRCLUT", dirclutpars)) {

      // array of bar y-positions for LUT from JGeometry
      double y = 0.; // no shift
      double x = DIRC_LUT_X[dirclutpars[1]];
      double z = DIRC_LUT_Z;

      G4ThreeVector vec(0,0,1);
      double rand1 = G4UniformRand();
      double rand2 = G4UniformRand();
      vec.setTheta(acos(rand1));
      vec.setPhi(2*M_PI*rand2);
      vec.rotateY(M_PI/2.);
      y = DIRC_BAR_Y[dirclutpars[1]];
      if (dirclutpars[1] < 24) {
        vec.rotateY(M_PI);
      }
     
      // spread over end of bar in y and z
      y += DIRC_QZBL_DY/2.0 - DIRC_QZBL_DY*G4UniformRand();
      z += DIRC_QZBL_DZ/2.0 - DIRC_QZBL_DZ*G4UniformRand(); 

      thetap = vec.theta();
      phip = vec.phi();
      fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));
   }



   double DeltaT = 0.;
   // Special case of Cherenkov photon gun for DIRC LED generator 
   if (user_opts->Find("DIRCLED", dircledpars)){

      double x(0.),y(0.),z(0.);

      int FDTH = -1;
      vector<int> FDTHs = {}; 
      for (int par_index = 1; par_index <= 6; par_index++)
      {
	 int passed_FDTH = dircledpars[par_index];
	 if (passed_FDTH && 0 < passed_FDTH && passed_FDTH < 7)
            FDTHs.push_back(dircledpars[par_index]);
      }
      int NumFDTHs = int(FDTHs.size());
      double rand0 = G4UniformRand();
      for (int FDTH_index = 0; FDTH_index < NumFDTHs; FDTH_index++)
      {
         if (rand0 <= (FDTH_index+1)*(1./NumFDTHs))
	 {
	    FDTH = FDTHs[FDTH_index]; break;
	 }
	 else
	    continue;
      }

      switch (FDTH)
      {
         case 1:
         x = DIRC_LED_OBCN_FDTH_X;
         y = DIRC_LED_OBCN_FDTH1_Y;
         z = DIRC_LED_OBCN_FDTH_Z;
         break;

	 case 2:
         x = DIRC_LED_OBCN_FDTH_X;
         y = DIRC_LED_OBCN_FDTH2_Y;
         z = DIRC_LED_OBCN_FDTH_Z;
	 break;

	 case 3:
         x = DIRC_LED_OBCN_FDTH_X;
         y = DIRC_LED_OBCN_FDTH3_Y;
         z = DIRC_LED_OBCN_FDTH_Z;
	 break;

         case 4:
         x = DIRC_LED_OBCS_FDTH_X;
         y = DIRC_LED_OBCS_FDTH4_Y;
         z = DIRC_LED_OBCS_FDTH_Z;
         break;

	 case 5:
         x = DIRC_LED_OBCS_FDTH_X;
         y = DIRC_LED_OBCS_FDTH5_Y;
         z = DIRC_LED_OBCS_FDTH_Z;
	 break;

	 case 6:
         x = DIRC_LED_OBCS_FDTH_X;
         y = DIRC_LED_OBCS_FDTH6_Y;
         z = DIRC_LED_OBCS_FDTH_Z;
	 break;
	
	 default:
	 break;
      }

      //z -= 0.5*cm;
      double theta_range = 12.5; // in degrees
      double inclination_wrt_bars_deg = 0.;//negative->away from 3-segment mirror; positive->towards 3-segment mirror
      double angle_towards_center_deg = 0.;//bending angle for the two feedthroughs on the sides towards the center 

      if (x > 0.)
         inclination_wrt_bars_deg *= -1.;
      if (FDTH == 3 || FDTH == 4)
         angle_towards_center_deg *= -1.;
      if (FDTH == 2 || FDTH == 5)
         angle_towards_center_deg = 0.;

      G4ThreeVector vec(0,0,1);
      double rand1 = G4UniformRand();
      double rand2 = G4UniformRand();

      double costheta = -1. + rand1 * (std::cos((180.-theta_range) * M_PI / 180.) + 1.);

      double rand3 = G4UniformRand();
      double theta_to_set = acos(costheta);
      if (rand3 < 0.5)
         theta_to_set = 2*M_PI - theta_to_set;

      vec.setTheta(theta_to_set);
      vec.setPhi(2*M_PI*rand2);
      vec.rotateY(inclination_wrt_bars_deg*deg);
      vec.rotateX(angle_towards_center_deg*deg);
/*
      //For square diffuser case
      double theta_range = 25.; // in degrees
      double diffuser_X = 1.697 * cm;
      double diffuser_Y = 1.697 * cm;
      double inclination_wrt_bars_deg = -6.;

      G4ThreeVector vec(0,0,-1.);
      double rand1 = G4UniformRand();
      double rand2 = G4UniformRand();
      double a = 2.* tan(theta_range*M_PI/180.);
      vec.setX((rand1-0.5)*a);
      vec.setY((rand2-0.5)*a);
      vec.setZ(-1.);

      double diffuser_offset_x = diffuser_X/2. - diffuser_X * G4UniformRand();
      double diffuser_offset_y = diffuser_Y/2. - diffuser_Y * G4UniformRand();

      G4ThreeVector diffuser_offset_vec(diffuser_offset_x,diffuser_offset_y,0.);
      diffuser_offset_vec.rotateY(inclination_wrt_bars_deg*deg);

      x += diffuser_offset_vec.x();
      y += diffuser_offset_vec.y();
      z += diffuser_offset_vec.z();

*/

      // time smearing from LED pulse shape 
      double t_rise = 0.84;
      double t_FWHM = 1.5;

      double t_range_LED = t_rise + t_FWHM;
      double t_rand1 = G4UniformRand();
      double t_rand2 = G4UniformRand();
      double t_LEDShape = t_rand1 * t_range_LED;

      double pulse_val = -1.;
      while (t_rand2 > pulse_val)
      {
         t_rand1 = G4UniformRand();
         t_rand2 = G4UniformRand();
         t_LEDShape = t_rand1 * t_range_LED;

         if (0.<= t_LEDShape && t_LEDShape < t_rise)
            pulse_val = 1./t_rise * t_LEDShape;
         else if (t_rise <= t_LEDShape && t_LEDShape < t_FWHM)
	    pulse_val = 1.;
         else if (t_FWHM <= t_LEDShape && t_LEDShape < t_range_LED)
	    pulse_val = -1./t_rise * t_LEDShape + 1./t_rise * t_range_LED;
	 else 
	    pulse_val = 0.;

	 if (t_rand2 <= pulse_val)
	    break;
      }
      double DeltaT_LED = t_LEDShape;


      // delay from LED feedthrough cables
      double DeltaT_CableDelay = 0.;
      if (FDTH == 2 || FDTH == 5)
	 DeltaT_CableDelay = 10.;
      if (FDTH == 3 || FDTH == 6)
	 DeltaT_CableDelay = 20.;


      DeltaT = DeltaT_LED + DeltaT_CableDelay;

      G4ThreeVector pos_vec(x,y,z);
      thetap = vec.theta();
      phip = vec.phi();
      fParticleGun->SetParticlePosition(pos_vec);

   }

   G4ThreeVector mom(p * sin(thetap) * cos(phip),
                     p * sin(thetap) * sin(phip),
                     p * cos(thetap));
   fParticleGun->SetParticleMomentum(mom);
   fParticleGun->SetParticleTime(0.+DeltaT);

   // Fire the gun
   fParticleGun->GeneratePrimaryVertex(anEvent);

   // Store generated particle info so it can be written to output file
   event_info->AddPrimaryVertex(*anEvent->GetPrimaryVertex());
}

//--------------------------------------------
// GeneratePrimariesHDDM
//--------------------------------------------

void GlueXPrimaryGeneratorAction::GeneratePrimariesHDDM(G4Event* anEvent)
{
   if (fPrimaryGenerator != 0) {
      double beamDiameter = GlueXPhotonBeamGenerator::getBeamDiameter();
      double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
      double x, y, z;
      if (fBeamvertex_activated == 0)
         configure_beam_vertex();
      if (fBeamvertex_activated == 1) {
 
         // generate a random beam spot according to a gaussian model
 
         double ur = sqrt(-2 * log(G4UniformRand()));
         double phi = G4UniformRand() * 2*M_PI;
         double uz = G4UniformRand() - 0.5;
         double u[2] = {ur * cos(phi), ur * sin(phi)};
         u[0] *= fBeamvertex.sigma[0];
         u[1] *= fBeamvertex.sigma[1];
         x = u[0] * cos(fBeamvertex.alpha) + u[1] * sin(fBeamvertex.alpha);
         y =-u[0] * sin(fBeamvertex.alpha) + u[1] * cos(fBeamvertex.alpha);
         z = uz * fBeamvertex.length;
         x += fBeamvertex.x + fBeamvertex.dxdz * z;
         y += fBeamvertex.y + fBeamvertex.dydz * z;
         z += fBeamvertex.z;
      }
      else {

         // fall back on the old hard-coded cylindrical beam spot
  
         while (true) {
            x = G4UniformRand() - 0.5;
            y = G4UniformRand() - 0.5;
            if (x*x + y*y <= 0.25) {
               x *= beamDiameter;
               y *= beamDiameter;
               break;
            }
         }
         z = fTargetCenterZ + (G4UniformRand() - 0.5) * fTargetLength;
      }
      assert (fPrimaryGenerator != 0);
      assert (fPhotonBeamGenerator != 0);
      double ttag = fPhotonBeamGenerator->GenerateTriggerTime(anEvent);
      double trel = (z - fRFreferencePlaneZ) / beamVelocity;
      fPrimaryGenerator->SetParticlePosition(G4ThreeVector(x,y,z));
      fPrimaryGenerator->SetParticleTime(trel + ttag);
      fPrimaryGenerator->GeneratePrimaryVertex(anEvent);
      GlueXUserEventInformation *eventinfo;
      eventinfo = (GlueXUserEventInformation*)anEvent->GetUserInformation();
      if (eventinfo) {
         // The above-assigned vertex coordinates and time are advisory to 
         // GeneratePrimaryVertex, and may have been overridden by values
         // read from the input MC event record.
         G4PrimaryVertex *vertex = anEvent->GetPrimaryVertex();
         trel = (vertex->GetZ0() - fRFreferencePlaneZ) / beamVelocity;
         ttag = vertex->GetT0() - trel;
         double E = eventinfo->GetBeamPhotonEnergy();
         fPhotonBeamGenerator->GenerateTaggerHit(anEvent, E, ttag);
      }
      else {
         return;
      }
   }

   // Superimpose any requested background minimum-bias beam interactions

   if (fBeamBackgroundRate > 0 && fPhotonBeamGenerator != 0) {
      double t = fBeamBackgroundGateStart;
      while (true) {
         t += -log(G4UniformRand()) / fBeamBackgroundRate;
         if (t > fBeamBackgroundGateStop)
            break;
         fPhotonBeamGenerator->GenerateBeamPhoton(anEvent, t);
      }
   }

   // Add the RF sync
   fPhotonBeamGenerator->GenerateRFsync(anEvent);
}

void GlueXPrimaryGeneratorAction::GeneratePrimariesCobrems(G4Event* anEvent)
{
   if (fPhotonBeamGenerator != 0) {
      fPhotonBeamGenerator->GeneratePrimaryVertex(anEvent);
   }
}

// Convert particle types from Geant3 types to PDG scheme

int GlueXPrimaryGeneratorAction::ConvertGeant3ToPdg(int Geant3type)
{
   // This method was imported from ROOT source file TDatabasePDG.cc

   switch (Geant3type) {
      case 1   : return 22;       // photon
      case 2   : return -11;      // e+
      case 3   : return 11;       // e-
      case 4   : return 12;       // e-neutrino (NB: flavour undefined by Geant)
      case 5   : return -13;      // mu+
      case 6   : return 13;       // mu-
      case 7   : return 111;      // pi0
      case 8   : return 211;      // pi+
      case 9   : return -211;     // pi-
      case 10  : return 130;      // K long
      case 11  : return 321;      // K+
      case 12  : return -321;     // K-
      case 13  : return 2112;     // n
      case 14  : return 2212;     // p
      case 15  : return -2212;    // anti-proton
      case 16  : return 310;      // K short
      case 17  : return 221;      // eta
      case 18  : return 3122;     // Lambda
      case 19  : return 3222;     // Sigma+
      case 20  : return 3212;     // Sigma0
      case 21  : return 3112;     // Sigma-
      case 22  : return 3322;     // Xi0
      case 23  : return 3312;     // Xi-
      case 24  : return 3334;     // Omega- (PB)
      case 25  : return -2112;    // anti-neutron
      case 26  : return -3122;    // anti-Lambda
      case 27  : return -3222;    // anti-Sigma-
      case 28  : return -3212;    // anti-Sigma0
      case 29  : return -3112;    // anti-Sigma+
      case 30  : return -3322;    // anti-Xi0
      case 31  : return -3312;    // anti-Xi+
      case 32  : return -3334;    // anti-Omega+ 
      case 45  : return 1000010020; // deuteron
      case 46  : return 1000010030; // triton
      case 47  : return 1000020040; // alpha
      case 48  : return 0;        // geantino (no PDG type)
      case 49  : return 1000020030; // He3 ion
      case 50  : return -22;        // Cerenkov photon (assigned in G4.10.7 release)
      case 61  : return 1000030060;  // Li6
      case 62  : return 1000030070;  // Li7
      case 63  : return 1000040070;  // Be7
      case 64  : return 1000040090;  // Be9
      case 65  : return 1000050100;  // B10
      case 66  : return 1000050110;  // B11
      case 67  : return 1000060120;  // C12
      case 68  : return 1000070140;  // N14
      case 69  : return 1000080160;  // O16
      case 70  : return 1000090190;  // F19
      case 71  : return 1000100200;  // Ne20
      case 72  : return 1000110230;  // Na23
      case 73  : return 1000120240;  // Mg24
      case 74  : return 1000130270;  // Al27
      case 75  : return 1000140280;  // Si28
      case 76  : return 1000150310;  // P31
      case 77  : return 1000160320;  // S32
      case 78  : return 1000170350;  // Cl35
      case 79  : return 1000180360;  // Ar36
      case 80  : return 1000190390;  // K39
      case 81  : return 1000200400;  // Ca40
      case 82  : return 1000210450;  // Sc45
      case 83  : return 1000220480;  // Ti48
      case 84  : return 1000230510;  // V51
      case 85  : return 1000240520;  // Cr52
      case 86  : return 1000250550;  // Mn55
      case 87  : return 1000260560;  // Fe56
      case 88  : return 1000270590;  // Co59
      case 89  : return 1000280580;  // Ni58
      case 90  : return 1000290630;  // Cu63
      case 91  : return 1000300640;  // Zn64
      case 92  : return 1000320740;  // Ge74
      case 93  : return 1000340800;  // Se80
      case 94  : return 1000360840;  // Kr84
      case 95  : return 1000380880;  // Sr88
      case 96  : return 1000400900;  // Zr90
      case 97  : return 1000420980;  // Mo98
      case 98  : return 1000461060;  // Pd106
      case 99  : return 1000481140;  // Cd114
      case 100 : return 1000501200;  // Sn120
      case 101 : return 1000541320;  // Xe132
      case 102 : return 1000561380;  // Ba138
      case 103 : return 1000581400;  // Ce140
      case 104 : return 1000621520;  // Sm152
      case 105 : return 1000661640;  // Dy164
      case 106 : return 1000701740;  // Yb174
      case 107 : return 1000741840;  // W184
      case 108 : return 1000781940;  // Pt194
      case 109 : return 1000791970;  // Au197
      case 110 : return 1000802020;  // Hg202
      case 111 : return 1000822080;  // Pb208
      case 112 : return 1000922380;  // U238

      // These are "private" Geant3 types that were defined in hdgeant
      case 33  : return 223;      // omega(782)
      case 34  : return 333;      // phi(1020)
      case 35  : return 331;      // etaPrime(958)
      case 36  : return 0;        // unused
      case 37  : return 0;        // unused
      case 38  : return 0;        // unused
      case 39  : return 0;        // unused
      case 40  : return 0;        // unused
      case 41  : return 0;        // unused
      case 42  : return 213;      // rho(770)+
      case 43  : return -213;     // rho(770)-
      case 44  : return 113;      // rho(770)0

      case 182  : return 2224;    // Delta++
      case 183  : return 443;     // Jpsi
      case 184  : return 441;     // Eta_c
      case 185  : return 10441;   // Chi_c0
      case 186  : return 20443;   // Chi_c1
      case 187  : return 445;     // Chi_c2
      case 188  : return 100443;  // Psi2s
      case 189  : return 421;     // D0
      case 190  : return 411;     // D+
      case 191  : return 10421;   // Dstar0
      case 192  : return 10411;   // Dstar+
      case 193  : return 4022;    // Lambda_c+
      case 194  : return -421;    // anti-D0

      case 163  : return 9000111; // a0(980)
      case 164  : return 9010221; // f0(980)
      case 165  : return 313;     // K*(892)0
      case 166  : return 323;     // K*(892)+
      case 167  : return -323;    // K*(892)-
      case 168  : return -313;    // anti-K*(892)0
      case 169  : return 20323;   // K1(1400)+
      case 170  : return -20323;  // K1(1400)-
      case 171  : return 4122;    // b1(1235)+
      case 172  : return 3224;    // Sigma*(1385)+
      case 173  : return 3214;    // Sigma*(1385)0
      case 174  : return 3114;    // Sigma*(1385)-

      default  :
         G4cout << "Warning in GlueXPrimaryGeneratorAction::"
                   "ConvertGeant3ToPdg - lookup performed on unknown"
                   " Geant3 particle type " << Geant3type << ","
                   " returning 0 for the PDG particle code." << G4endl;
   }
   return 0;
}

// Convert particle types from PDG scheme to Geant3 types

int GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(int PDGtype)
{
   // Invert the table contained in ConvertGeant3ToPdg

   switch (PDGtype) {
      case          0 : return 50;    // old usage, optical photon
      case        -22 : return 50;    // optical photon
      case         22 : return 1;     // photon
      case        -11 : return 2;     // e+
      case         11 : return 3;     // e-
      case        -12 : return 4;     // anti-e-neutrino
      case         12 : return 4;     // e-neutrino
      case        -13 : return 5;     // mu+
      case         13 : return 6;     // mu-
      case        -14 : return 4;     // anti-mu-neutrino
      case         14 : return 4;     // mu-neutrino
      case        -16 : return 4;     // anti-tau-neutrino
      case         16 : return 4;     // tau-neutrino
      case        111 : return 7;     // pi0
      case        211 : return 8;     // pi+
      case       -211 : return 9;     // pi-
      case        130 : return 10;    // K long
      case        321 : return 11;    // K+
      case       -321 : return 12;    // K-
      case       2112 : return 13;    // n
      case       2212 : return 14;    // p
      case      -2212 : return 15;    // anti-proton
      case        310 : return 16;    // K short
      case        221 : return 17;    // eta
      case       3122 : return 18;    // Lambda
      case       3222 : return 19;    // Sigma+
      case       3212 : return 20;    // Sigma0
      case       3112 : return 21;    // Sigma-
      case       3322 : return 22;    // Xi0
      case       3312 : return 23;    // Xi-
      case       3334 : return 24;    // Omega-
      case      -2112 : return 25;    // anti-neutron
      case      -3122 : return 26;    // anti-Lambda
      case      -3222 : return 27;    // Sigma-
      case      -3212 : return 28;    // Sigma0
      case      -3112 : return 29;    // Sigma+
      case      -3322 : return 30;    // Xi0
      case      -3312 : return 31;    // Xi+
      case      -3334 : return 32;    // Omega+
      case 1000010020 : return 45;    // deuteron
      case 1000010030 : return 46;    // triton
      case 1000020040 : return 47;    // alpha
      case 1000020030 : return 49;    // He3
      case 1000030060 : return 61;    // Li6
      case 1000030070 : return 62;    // Li7
      case 1000040070 : return 63;    // Be7
      case 1000040090 : return 64;    // Be9
      case 1000050100 : return 65;    // B10
      case 1000050110 : return 66;    // B11
      case 1000060120 : return 67;    // C12
      case 1000070140 : return 68;    // N14
      case 1000080160 : return 69;    // O16
      case 1000090190 : return 70;    // F19
      case 1000100200 : return 71;    // Ne20
      case 1000110230 : return 72;    // Na23
      case 1000120240 : return 73;    // Mg24
      case 1000130270 : return 74;    // Al27
      case 1000140280 : return 75;    // Si28
      case 1000150310 : return 76;    // P31
      case 1000160320 : return 77;    // S32
      case 1000170350 : return 78;    // Cl35
      case 1000180360 : return 79;    // Ar36
      case 1000190390 : return 80;    // K39
      case 1000200400 : return 81;    // Ca40
      case 1000210450 : return 82;    // Sc45
      case 1000220480 : return 83;    // Ti48
      case 1000230510 : return 84;    // V51
      case 1000240520 : return 85;    // Cr52
      case 1000250550 : return 86;    // Mn55
      case 1000260560 : return 87;    // Fe56
      case 1000270590 : return 88;    // Co59
      case 1000280580 : return 89;    // Ni58
      case 1000290630 : return 90;    // Cu63
      case 1000300640 : return 91;    // Zn64
      case 1000320740 : return 92;    // Ge74
      case 1000340800 : return 93;    // Se80
      case 1000360840 : return 94;    // Kr84
      case 1000380880 : return 95;    // Sr88
      case 1000400900 : return 96;    // Zr90
      case 1000420980 : return 97;    // Mo98
      case 1000461060 : return 98;    // Pd106
      case 1000481140 : return 99;    // Cd114
      case 1000501200 : return 100;   // Sn120
      case 1000541320 : return 101;   // Xe132
      case 1000561380 : return 102;   // Ba138
      case 1000581400 : return 103;   // Ce140
      case 1000621520 : return 104;   // Sm152
      case 1000661640 : return 105;   // Dy164
      case 1000701740 : return 106;   // Yb174
      case 1000741840 : return 107;   // W184
      case 1000781940 : return 108;   // Pt194
      case 1000791970 : return 109;   // Au197
      case 1000802020 : return 110;   // Hg202
      case 1000822080 : return 111;   // Pb208
      case 1000922380 : return 112;   // U238

      // These are "private" Geant3 types that were defined in hdgeant
      case 223        : return 33;    // omega(782)
      case 333        : return 34;    // phi(1020)
      case 331        : return 35;    // etaPrime(958)
      case 213        : return 42;    // rho(770)+
      case -213       : return 43;    // rho(770)-
      case 113        : return 44;    // rho(770)0

      case 2224       : return 182;    // Delta++
      case 443        : return 183;    // Jpsi
      case 441        : return 184;    // Eta_c
      case 10441      : return 185;    // Chi_c0
      case 20443      : return 186;    // Chi_c1
      case 445        : return 187;    // Chi_c2
      case 100443     : return 188;    // Psi2s
      case 421        : return 189;    // D0
      case 411        : return 190;    // D+
      case 10421      : return 191;    // Dstar0
      case 10411      : return 192;    // Dstar+
      case 4022       : return 193;    // Lambda_c+
      case -421       : return 194;    // anti-D0

      case 9000111    : return 163;   // a0(980)
      case 9010221    : return 164;   // f0(980)
      case 313        : return 165;   // K*(892)0
      case 323        : return 166;   // K*(892)+
      case -323       : return 167;   // K*(892)-
      case -313       : return 168;   // anti-K*(892)0
      case 20323      : return 169;   // K1(1400)+
      case -20323     : return 170;   // K1(1400)-
      case 4122       : return 171;   // b1(1235)+
      case 3224       : return 172;   // Sigma*(1385)+
      case 3214       : return 173;   // Sigma*(1385)0
      case 3114       : return 174;   // Sigma*(1385)-

      default  :
         if (PDGtype < 1000000000) {
            G4cout << "Warning in GlueXPrimaryGeneratorAction::"
                      "ConvertPdgToGeant3 - lookup performed on unknown"
                      " PDG particle type " << PDGtype << ","
                      " returning 0 for the Geant3 particle code." << G4endl;
         }
   }
   return 0;
}

double GlueXPrimaryGeneratorAction::GetMassPDG(int PDGtype)
{
   return fParticleTable->FindParticle(PDGtype)->GetPDGMass();
}

double GlueXPrimaryGeneratorAction::GetMass(int Geant3Type)
{
   return GetMassPDG(ConvertGeant3ToPdg(Geant3Type));
}

G4ParticleDefinition *GlueXPrimaryGeneratorAction::GetParticle(int PDGtype)
{
   G4ParticleDefinition *p = fParticleTable->FindParticle(PDGtype);
   if (p==0) {
      if (PDGtype > 1000000000) {
         p = fParticleTable->GetIonTable()->GetIon(PDGtype);
      }
      else if (PDGtype == -22) {
         p = fParticleTable->FindParticle(22);
      }
      else {
         G4cout << "unknown particle type " << PDGtype
                << ", substituting geantino in its place!"
                << G4endl;
         p = fParticleTable->FindParticle("geantino");
      }
   }
   return p;
}

G4ParticleDefinition *GlueXPrimaryGeneratorAction::GetParticle(const G4String &name)
{
   G4ParticleDefinition *p = fParticleTable->FindParticle(name);
   if (p==0) {
      G4cout << "unknown particle type " << name
             << ", substituting geantino in its place!"
             << G4endl;
      p = fParticleTable->FindParticle("geantino");
   }
   return p;
}

void GlueXPrimaryGeneratorAction::configure_beam_vertex()
{
   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
   if (user_opts == 0) {
      G4cerr << "Error in GlueXPrimaryGeneratorAction::configure_beam_vertex"
             << " - GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }
   std::map<int, std::string> vertex_spec;
   if (user_opts->Find("VERTEX", vertex_spec)) {
      if (vertex_spec.find(1) == vertex_spec.end()) {
         G4cerr << "Error in GlueXPrimaryGeneratorAction::configure_beam_vertex"
                << " - VERTEX card found but no specification was given."
                << G4endl;
      }
   }
   else {
      fBeamvertex_activated = -1;
      return;
   }
   const char *spec = vertex_spec[1].c_str();
   if (sscanf(spec, "beam_spot(ccdb) * %lf", &fBeamvertex.length)) {
      int runno = HddmOutput::getRunNo();
      extern JApplication *japp;
      if (japp == 0) {
         G4cerr << "Error in GlueXPrimaryGeneratorAction::configure_beam_vertex"
                << " - jana global DApplication object not set, "
                << "cannot continue." << G4endl;
         exit(-1);
      }
      JCalibration *jcalib = japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
      std::map<string, float> beam_spot_params;
      jcalib->Get("/PHOTON_BEAM/beam_spot", beam_spot_params);
      fBeamvertex.x = beam_spot_params.at("x") * cm;
      fBeamvertex.y = beam_spot_params.at("y") * cm;
      fBeamvertex.z = beam_spot_params.at("z") * cm;
      fBeamvertex.var_xx = beam_spot_params.at("var_xx") * cm*cm;
      fBeamvertex.var_xy = beam_spot_params.at("var_xy") * cm*cm;
      fBeamvertex.var_yy = beam_spot_params.at("var_yy") * cm*cm;
      fBeamvertex.dxdz = beam_spot_params.at("dxdz");
      fBeamvertex.dydz = beam_spot_params.at("dydz");
      fBeamvertex.length *= cm;
   }
   else if (sscanf(spec, "beam_spot(%*[-+0-9.,]) * %lf", &fBeamvertex.length)) {
      sscanf(spec, "beam_spot(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
                   &fBeamvertex.x, &fBeamvertex.y, &fBeamvertex.z,
                   &fBeamvertex.var_xx, &fBeamvertex.var_xy,
                   &fBeamvertex.var_yy, &fBeamvertex.dxdz,
                   &fBeamvertex.dydz);
      fBeamvertex.x *= cm;
      fBeamvertex.y *= cm;
      fBeamvertex.z *= cm;
      fBeamvertex.var_xx *= cm*cm;
      fBeamvertex.var_xy *= cm*cm;
      fBeamvertex.var_yy *= cm*cm;
      fBeamvertex.length *= cm;
   }
   else {
      G4cerr << "Error in GlueXPrimaryGeneratorAction::configure_beam_vertex"
             << " - unrecognized VERTEX specification on control.in line"
             << G4endl << spec << G4endl;
      exit(-1);
   }
   double D = fBeamvertex.var_xx * fBeamvertex.var_yy -
              fBeamvertex.var_xy * fBeamvertex.var_xy;
   double A = (fBeamvertex.var_xx + fBeamvertex.var_yy)/2;
   double B = (fBeamvertex.var_xx - fBeamvertex.var_yy)/2;
   double C = fBeamvertex.var_xy;
   double evalue1 = (A + sqrt(B*B + C*C) + 1e-20) / (D + 1e-40);
   double evalue2 = (A - sqrt(B*B + C*C) + 1e-20) / (D + 1e-40);
   if (evalue1 < 0 || evalue2 < 0) {
      G4cerr << "Error in GlueXPrimaryGeneratorAction::configure_beam_vertex"
             << " - unphysical values given for the beam spot ellipse:"
             << " var_xx=" << fBeamvertex.var_xx
             << ", var_xy=" << fBeamvertex.var_xy
             << ", var_yy=" << fBeamvertex.var_yy
             << G4endl;
      exit(1);
   }
   double alpha = 0;
   if (C == 0)
      alpha = (B < 0)? 0 : M_PI/2;
   else if (B == 0)
      alpha = (C < 0)? -M_PI/4 : M_PI/4;
   else
      alpha = atan2(C,B) / 2;
   fBeamvertex.sigma[0] = sqrt(1 / evalue1);
   fBeamvertex.sigma[1] = sqrt(1 / evalue2);
   fBeamvertex.alpha = alpha;
   fBeamvertex_activated = 1;

#ifndef QUIET_CONFIGURE_BEAM_VERTEX
   G4cout << "Configured beam vertex parameters: (units of cm)" << G4endl
          << "   center x=" << fBeamvertex.x/cm
          << "   y=" << fBeamvertex.y/cm
          << "   z=" << fBeamvertex.z/cm << G4endl
          << "   sigma_x=" << sqrt(fBeamvertex.var_xx)/cm
          << "   sigma_y=" << sqrt(fBeamvertex.var_yy)/cm
          << "   sigma_xy=" << fBeamvertex.var_xy/(cm*cm) << G4endl
          << "   dxdz=" << fBeamvertex.dxdz
          << "   dydz=" << fBeamvertex.dydz
          << "   dz=" << fBeamvertex.length/cm << G4endl;
#endif
}

void GlueXPrimaryGeneratorAction::generate_beam_vertex(double v[3])
{
}

void GlueXPrimaryGeneratorAction::clone_photon_beam_generator()
{
   fCobremsGeneration = new CobremsGeneration(fBeamEndpointEnergy/GeV,
                                              fBeamPeakEnergy/GeV);
   CobremsGeneration *gen0 = (*fInstance.begin())->fCobremsGeneration;
   double beamEmin = gen0->getPhotonEnergyMin()*GeV;
   double radColDist = gen0->getCollimatorDistance()*m;
   double colDiam = gen0->getCollimatorDiameter()*m;
   double beamEmit = gen0->getBeamEmittance()*(m*radian);
   double radThick = gen0->getTargetThickness()*m;
   double spotRMS = gen0->getCollimatorSpotrms()*m;
   GlueXPhotonBeamGenerator *gen1 = (*fInstance.begin())->fPhotonBeamGenerator;
   double spotX = gen1->getBeamOffset(0);
   double spotY = gen1->getBeamOffset(1);
   fCobremsGeneration->setPhotonEnergyMin(beamEmin/GeV);
   fCobremsGeneration->setCollimatorDistance(radColDist/m);
   fCobremsGeneration->setCollimatorDiameter(colDiam/m);
   fCobremsGeneration->setBeamEmittance(beamEmit/(m*radian));
   fCobremsGeneration->setTargetThickness(radThick/m);
   fCobremsGeneration->setCollimatorSpotrms(spotRMS/m);
   fPhotonBeamGenerator = new GlueXPhotonBeamGenerator(fCobremsGeneration);
   fPhotonBeamGenerator->setBeamOffset(spotX, spotY);
}

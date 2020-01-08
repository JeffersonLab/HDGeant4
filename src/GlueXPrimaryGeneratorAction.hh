//
// GlueXPrimaryGeneratorAction class header
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread,
// but virtually all of its functions need to be serialized, so
// it maintains its own interlocks for this purpose. Resources
// are created once when the first object is instantiated, and
// destroyed once when the last object is destroyed.

#ifndef _GLUEXPRIMARYGENERATORACTION_H_
#define _GLUEXPRIMARYGENERATORACTION_H_

#include "CobremsGeneration.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "GlueXPhotonBeamGenerator.hh"
#include "GlueXPrimaryGenerator.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4AutoLock.hh"
#include "GlueXParticleGun.hh"
#include "G4Event.hh"

#include <HDDM/hddm_s.hpp>
#include <JANA/JApplication.h>

#include <fstream>

class GlueXPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
   
   enum source_type_t {
      SOURCE_TYPE_NONE,
      SOURCE_TYPE_PARTICLE_GUN,
      SOURCE_TYPE_COBREMS_GEN,
      SOURCE_TYPE_HDDM
   };

   GlueXPrimaryGeneratorAction();
   GlueXPrimaryGeneratorAction(const GlueXPrimaryGeneratorAction &src);
   GlueXPrimaryGeneratorAction &operator=(const GlueXPrimaryGeneratorAction &src);
   ~GlueXPrimaryGeneratorAction();
   
   virtual void GeneratePrimaries(G4Event* anEvent);
   void GeneratePrimariesHDDM(G4Event* anEvent);
   void GeneratePrimariesParticleGun(G4Event* anEvent);
   void GeneratePrimariesCobrems(G4Event* anEvent);

   static int ConvertGeant3ToPdg(int Geant3Type);
   static int ConvertPdgToGeant3(int PDGtype);
   static G4ParticleDefinition *GetParticle(int PDGtype);
   static G4ParticleDefinition *GetParticle(const G4String &name);
   static double GetMassPDG(int PDGtype);
   static double GetMass(int Geant3Type);
 
   static const GlueXPrimaryGeneratorAction* GetInstance();

 private:
   static int instanceCount;
   static source_type_t fSourceType;
   static std::ifstream *fHDDMinfile;
   static hddm_s::istream *fHDDMistream;
   static CobremsGeneration *fCobremsGeneration;
   static GlueXPhotonBeamGenerator *fPhotonBeamGenerator;
   static G4ParticleTable *fParticleTable;
   static GlueXParticleGun *fParticleGun;
   static GlueXPrimaryGenerator *fPrimaryGenerator;

 public:
   struct single_particle_gun_t {
      int geantType;
      int pdgType;
      G4ParticleDefinition *partDef;
      G4ThreeVector pos;
      double mom;
      double theta;
      double phi;
      double deltaR;
      double deltaZ;
      double deltaMom;
      double deltaTheta;
      double deltaPhi;
      int plogOption;
      int tlogOption;
   };

 private:
   static single_particle_gun_t fGunParticle;

   static double fBeamBackgroundRate;
   static double fBeamBackgroundGateStart;
   static double fBeamBackgroundGateStop;
   static double fL1triggerTimeSigma;
   static double fRFreferencePlaneZ;

   static int fEventCount;

   // The following parameters describe the dimensions of the target
   // that are used when generating the primary interaction vertex for
   // events from an external generator. An external event generator
   // knows nothing about the simulation geometry, so it makes sense
   // that this should be modeled in the simulation. They only apply
   // to the HDDM input source.  They are initialized to default values
   // for the GlueX liquid hydrogen target in the constructor, but
   // can be accessed/changed by the getter/setter methods below.
   static double fTargetCenterZ;
   static double fTargetLength;

 public:
   static void setTargetCenterZ(double Z) {
      fTargetCenterZ = Z;
   }
   static void setTargetLength(double L) {
      fTargetLength = L;
   }
   static double getTargetCenterZ() {
      return fTargetCenterZ;
   }
   static double getTargetLength() {
      return fTargetLength;
   }
   static void setL1triggerTimeSigma(double sigma) {
      fL1triggerTimeSigma = sigma;
   }
   static double getL1triggerTimeSigma() {
      return fL1triggerTimeSigma;
   }
   static void setRFreferencePlaneZ(double refZ) {
      fRFreferencePlaneZ = refZ;
   }
   static double getRFreferencePlaneZ() {
      return fRFreferencePlaneZ;
   }
   int getEventCount() const {
      return fEventCount;
   }

 private:
   static G4Mutex fMutex;
   static std::list<GlueXPrimaryGeneratorAction*> fInstance;

   double DIRC_LUT_X[48], DIRC_LUT_Y[48], DIRC_LUT_Z[48]; //Pre-Rotated Positions
   double BXXA_boxX[48], BXXA_boxY[48], BXXA_boxZ[48];

   G4ThreeVector DCML10_InHall, DCML11_InHall, DCML01_InHall, DCML00_InHall;
   vector<double> DCML10_rot, DCML11_rot, DCML01_rot, DCML00_rot;

   vector<double> MRAS_rot, MRAN_rot;
   vector<double> OBCS_rot, OBCN_rot;
   vector<G4ThreeVector> FDTHs_InHall;

};

#endif

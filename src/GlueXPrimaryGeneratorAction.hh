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
 
   const CobremsGeneration* GetCobremsGeneration() const;

 private:
   static std::ifstream *fHDDMinfile;
   static hddm_s::istream *fHDDMistream;
   static G4ParticleTable *fParticleTable;
   static source_type_t fSourceType;

   GlueXParticleGun *fParticleGun;
   CobremsGeneration *fCobremsGeneration;
   GlueXPhotonBeamGenerator *fPhotonBeamGenerator;
   GlueXPrimaryGenerator *fPrimaryGenerator;

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
   static double fBeamEndpointEnergy;
   static double fBeamPeakEnergy;

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

   // The following model parameters are introduced to superseded the
   // above model which assumes a uniform circular beam spot that is
   // centered on the z axis. The values for these parameters are
   // specified using the VERTEX line in control.in, either by including
   // the values directly on the VERTEX line or by a lookup reference to
   // ccdb. This more sophisticated model incorporates a Gaussian beam
   // spot with a general elliptical transverse profile, assumed to be
   // uniformly distributed along the length of the target in z. The
   // walk of the beam centroid in the transverse plane with z is also
   // allowed for by the model. Note that the length of the beam spot in
   // z is not redundant with the target thickness, as one may wish to
   // generate events with a more restricted fiducial volume than would
   // be allowed by the geometry alone. All lengths are specified in
   // standard Geant4 units.
   struct beam_spot_t {
      double x;       // centroid of beam spot in x
      double y;       // centroid of beam spot in y
      double z;       // centroid in beam spot in z
      double var_xx;  // variance of beam spot in x
      double var_xy;  // covariance of beam spot in xy
      double var_yy;  // variance of beam spot in y
      double dxdz;    // slope of beam spot center walk in xz plane
      double dydz;    // slope of beam spot center walk in yz plane
      double length;  // length of the beam spot in z
      double sigma[2];// minor,major sigmas of the gaussian ellipse
      double alpha;   // rotation angle in xy plane of the beam ellipse
   };
   struct beam_spot_t fBeamvertex;
   int fBeamvertex_activated;  // uninitialized=0, active=1, disabled=-1

   void clone_photon_beam_generator();

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
   static void setBeamBackgroundRate(double rate_GHz) {
      fBeamBackgroundRate = rate_GHz;
   }
   static double getBeamBackgroundRate() {
      return fBeamBackgroundRate;
   }
   static void setBeamBackgroundGateStart(double t) {
      fBeamBackgroundGateStart = t;
   }
   static double getBeamBackgroundGateStart() {
      return fBeamBackgroundGateStart;
   }
   static void setBeamEndpointEnergy(double E0) {
      fBeamEndpointEnergy = E0;
   }
   static double getBeamEndpointEnergy() {
      return fBeamEndpointEnergy;
   }
   static void setBeamPeakEnergy(double Epeak) {
      fBeamPeakEnergy = Epeak;
   }
   static double getBeamPeakEnergy() {
      return fBeamPeakEnergy;
   }
   void configure_beam_vertex();
   void generate_beam_vertex(double v[3]);

 private:
   static G4Mutex fMutex;
   static std::list<GlueXPrimaryGeneratorAction*> fInstance;

  //double DIRC_LUT_X[48], DIRC_LUT_Y[48], DIRC_LUT_Z[48]; //Pre-Rotated Positions
   double BXXA_boxX[48], BXXA_boxY[48], BXXA_boxZ[48];

   G4ThreeVector DCML10_InHall, DCML11_InHall, DCML01_InHall, DCML00_InHall;
   vector<double> DCML10_rot, DCML11_rot, DCML01_rot, DCML00_rot;

   vector<double> MRAS_rot, MRAN_rot;
   vector<double> OBCS_rot, OBCN_rot;
   vector<G4ThreeVector> FDTHs_InHall;
  
  static double DIRC_LUT_X[48], DIRC_LUT_Y[48], DIRC_LUT_Z[48];
  static double DIRC_QZBL_DY;
  static double DIRC_QZBL_DZ;
  static double DIRC_OWDG_DZ;
  static double DIRC_LED_OBCS_FDTH_X;
  static double DIRC_LED_OBCS_FDTH_Z;
  static double DIRC_LED_OBCN_FDTH_X;
  static double DIRC_LED_OBCN_FDTH_Z;
  static double DIRC_LED_OBCN_FDTH1_Y;
  static double DIRC_LED_OBCN_FDTH2_Y;
  static double DIRC_LED_OBCN_FDTH3_Y;
  static double DIRC_LED_OBCS_FDTH4_Y;
  static double DIRC_LED_OBCS_FDTH5_Y;
  static double DIRC_LED_OBCS_FDTH6_Y;
  
};

#endif

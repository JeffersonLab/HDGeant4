//
// GlueXUserEventInformation class header
//
// author: richard.t.jones at uconn.edu
// version: aug 1, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread.

#ifndef _GLUEXUSEREVENTINFORMATION_
#define _GLUEXUSEREVENTINFORMATION_

#include "G4VUserEventInformation.hh"
#include "G4PrimaryVertex.hh"
#include "G4ThreeVector.hh"
#include <G4TrackVector.hh>

#include <HDDM/hddm_s.hpp>

#include <map>

class BCALincidentParticle;

class GlueXUserEventInformation: public G4VUserEventInformation
{
 public:
   GlueXUserEventInformation(hddm_s::HDDM *hddmevent=NULL);
   ~GlueXUserEventInformation();

   void AddBeamParticle(int geanttype, double t0, const G4ThreeVector &pos, 
                                                  const G4ThreeVector &mom);
   void AddBeamParticle(int geanttype, double t0, const G4ThreeVector &pos, 
                                                  const G4ThreeVector &mom,
                                                  const G4ThreeVector &pol);
   void AddTargetParticle(int geanttype, double t0, const G4ThreeVector &pos, 
                                                    const G4ThreeVector &mom);
   void AddTargetParticle(int geanttype, double t0, const G4ThreeVector &pos, 
                                                    const G4ThreeVector &mom,
                                                    const G4ThreeVector &pol);
   void AddPrimaryVertex(const G4PrimaryVertex &vertex);
   void AddSecondaryVertex(const G4TrackVector &secondaries,
                                               int parentID, int mech);
   void AddMCtrajectoryPoint(const G4Step &step, int save_option);

   int GetRunNo();
   double GetBeamPhotonEnergy();
   int GetGlueXTrackID(int g4ID);
   int GetGlueXTrackID(const G4Track *track);
   void SetGlueXTrackID(int g4ID, int gluexID);
   int AssignNextGlueXTrackID(const G4Track *track = 0);
   int AssignBCALincidentID(const G4Track *track);
   const BCALincidentParticle *GetBCALincidentParticle(int incidentID);

   void SetKeepEvent(int flag) { fKeepEvent = flag; }
   int GetKeepEvent() const { return fKeepEvent; }

   static void SetStartingSeeds(const long int seeds[2]);
   void SetRandomSeeds();
   void Print() const;

   hddm_s::HDDM *getOutputRecord() {
      return fOutputRecord;
   }

   static int fWriteNoHitEvents;
   static long int *fStartingSeeds;

 protected:
   hddm_s::HDDM *fOutputRecord;
   bool fKeepEvent;
   int fNprimaries;
   int fNvertices;
   std::map<int,int> fGlueXTrackID;
   std::vector<BCALincidentParticle> fBCALincidentParticle;

 private:
   GlueXUserEventInformation(const GlueXUserEventInformation &src);
   GlueXUserEventInformation &operator=(const GlueXUserEventInformation &src);
};

class BCALincidentParticle {
 public:
   G4ThreeVector pos;
   G4ThreeVector mom;
   double E;
   int ptype;
   int trackID;
};

#endif

//
// GlueXSensitiveDetectorCDC - class header
//
// author: richard.t.jones at uconn.edu
// version: august 28, 2015
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorCDC_h
#define GlueXSensitiveDetectorCDC_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitCDCstraw.hh"
#include "GlueXHitCDCpoint.hh"

class G4Step;
class G4HCofThisEvent;

typedef std::vector<GlueXHitCDCstraw::hitinfo_t> hit_vector_t;

class GlueXSensitiveDetectorCDC : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorCDC(const G4String& name);
   GlueXSensitiveDetectorCDC(const GlueXSensitiveDetectorCDC &right);
   GlueXSensitiveDetectorCDC &operator=(const GlueXSensitiveDetectorCDC &right);
   virtual ~GlueXSensitiveDetectorCDC();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   double asic_response(double t_ns); 
   double cdc_wire_signal_mV(double t_ns, hit_vector_t &hits);
   void add_cluster(hit_vector_t &hits, GlueXHitCDCstraw::hitinfo_t &h,
                    int n_p, double t, G4ThreeVector &x);
   void polint(double *xa, double *ya, int n, double x, double *y, double *dy);

 private:
   GlueXHitsMapCDCstraw* fStrawsMap;
   GlueXHitsMapCDCpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static const double ELECTRON_CHARGE;
   static double DRIFT_SPEED;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_KEV;
   static double THRESH_MV;
   static double STRAW_RADIUS;
   static double CDC_TIME_WINDOW;
   static double GAS_GAIN;
   static double W_EFF_PER_ION;
   static int N_SECOND_PER_PRIMARY;
   static int MAX_HITS;

   static int fDrift_clusters;
   const static int CDC_DRIFT_TABLE_LEN = 78;
   static double fDrift_time[CDC_DRIFT_TABLE_LEN];
   static double fDrift_distance[CDC_DRIFT_TABLE_LEN];
   static double fBscale_par1;
   static double fBscale_par2;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif

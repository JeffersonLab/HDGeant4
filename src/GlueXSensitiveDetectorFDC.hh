//
// GlueXSensitiveDetectorFDC - class header
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorFDC_h
#define GlueXSensitiveDetectorFDC_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitFDCwire.hh"
#include "GlueXHitFDCcathode.hh"
#include "GlueXHitFDCpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorFDC : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorFDC(const G4String& name);
   GlueXSensitiveDetectorFDC(const GlueXSensitiveDetectorFDC &right);
   GlueXSensitiveDetectorFDC &operator=(const GlueXSensitiveDetectorFDC &right);
   virtual ~GlueXSensitiveDetectorFDC();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   double Ei(double x);
   double asic_response(double t_ns); 
   double fdc_wire_signal_mV(
          double t_ns, std::vector<GlueXHitFDCwire::hitinfo_t> &hits);
   double fdc_cathode_signal_mV(
          double t_ns, std::vector<GlueXHitFDCcathode::hitinfo_t> &hits);
   int add_anode_hit(std::vector<GlueXHitFDCwire::hitinfo_t> &hits, 
                     GlueXHitFDCwire::hitinfo_t &wirehit,
                     int layer, 
                     double xwire,
                     G4ThreeVector &xglobal, 
                     G4ThreeVector &xlocal,
                     double dE, double t, double &tdrift);
   void add_cathode_hit(GlueXHitFDCwire::hitinfo_t &wirehit, int packageNo,
                        double xwire, double yavalanche, double tdrift,
                        int n_p, int chamber, int module, int layer,
                        int global_wire_number);
   void polint(double *xa, double *ya, int n, double x, double *y, double *dy);
   int locate(double *xx, int n, double x);

 private:
   GlueXHitsMapFDCwire* fWiresMap;
   GlueXHitsMapFDCcathode* fCathodesMap;
   GlueXHitsMapFDCpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static const double ELECTRON_CHARGE;
   static double DRIFT_SPEED;
   static double ACTIVE_AREA_OUTER_RADIUS;
   static double ANODE_CATHODE_SPACING;
   static double TWO_HIT_TIME_RESOL;
   static int    WIRES_PER_PLANE;
   static double WIRE_SPACING;
   static double STRIP_SPACING;
   static double U_OF_WIRE_ONE;
   static int    STRIPS_PER_PLANE;
   static double CATHODE_ROT_ANGLE;
   static double U_OF_STRIP_ONE;
   static double STRIP_GAP;
   static double K2;
   static double STRIP_NODES;
   static double THRESH_KEV;
   static double THRESH_ANODE;
   static double THRESH_STRIPS;
   static double DIFFUSION_COEFF;
   static double FDC_TIME_WINDOW;
   static double GAS_GAIN;
   static double W_EFF_PER_ION;
   static double N_SECOND_PER_PRIMARY;
   static int MAX_HITS;

   static double LORENTZ_NR_PAR1;
   static double LORENTZ_NR_PAR2;
   static double LORENTZ_NZ_PAR1;
   static double LORENTZ_NZ_PAR2;
   static double DRIFT_RES_PARMS[3];
   static double DRIFT_FUNC_PARMS[6];
   static double DRIFT_BSCALE_PAR1;
   static double DRIFT_BSCALE_PAR2;

   static int fDrift_clusters;

   static double wire_dead_zone_radius[4];
   static double strip_dead_zone_radius[4];

   static int drift_table_len;
   static double *drift_table_t_ns;
   static double *drift_table_d_cm;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif

//
// GlueXSensitiveDetectorFDC - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016

#include "GlueXSensitiveDetectorFDC.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"
#include "HddmOutput.hh"

#include <CLHEP/Random/RandPoisson.h>
#include <Randomize.hh>

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <JANA/JApplication.h>
#include <JANA/Calibrations/JCalibrationManager.h>


#include <stdlib.h>
#include <math.h>

const double fC = 1e-15 * coulomb;
const double pC = 1e-12 * coulomb;
const double GlueXSensitiveDetectorFDC::ELECTRON_CHARGE = 1.6022e-4*fC;

// Drift speed 2.2cm/us is appropriate for a 90/10 Argon/Methane mixture
double GlueXSensitiveDetectorFDC::DRIFT_SPEED = 0.0055*cm/ns;

// Minimum hit time difference for two hits on the same wire
double GlueXSensitiveDetectorFDC::TWO_HIT_TIME_RESOL = 25*ns;

// Cutoff on the number of allowed hits per anode wire or cathode strip
int GlueXSensitiveDetectorFDC::MAX_HITS = 1000;

// Minimum energy deposition for a wire hit (keV, pC)
double GlueXSensitiveDetectorFDC::THRESH_KEV = 1.;
double GlueXSensitiveDetectorFDC::THRESH_ANODE = 1.*pC;
double GlueXSensitiveDetectorFDC::THRESH_STRIPS = 5.*pC;

// Straw hits are accepted from t=0 up to this maximum time
double GlueXSensitiveDetectorFDC::FDC_TIME_WINDOW = 1000*ns;

// Parameters for setting signal pulse height
double GlueXSensitiveDetectorFDC::GAS_GAIN = 8e4;

// Average number of secondary ion pairs for 40/60 Ar/CO2 mixture
double GlueXSensitiveDetectorFDC::N_SECOND_PER_PRIMARY = 1.89; 

// Average energy needed to produce an ion pair for 40/60 mixture
double GlueXSensitiveDetectorFDC::W_EFF_PER_ION = 30.2*eV;

// Geometry parameters in the FDC chambers
int GlueXSensitiveDetectorFDC::WIRES_PER_PLANE = 96;
int GlueXSensitiveDetectorFDC::STRIPS_PER_PLANE = 192;
double GlueXSensitiveDetectorFDC::ACTIVE_AREA_OUTER_RADIUS = 48.5*cm;
double GlueXSensitiveDetectorFDC::ANODE_CATHODE_SPACING = 0.5*cm;
double GlueXSensitiveDetectorFDC::WIRE_SPACING = 1.0*cm;
double GlueXSensitiveDetectorFDC::STRIP_SPACING = 0.5*cm;
double GlueXSensitiveDetectorFDC::U_OF_WIRE_ONE = -47.5*cm; //-(WIRES_PER_PLANE-1)*WIRE_SPACING/2
double GlueXSensitiveDetectorFDC::U_OF_STRIP_ONE = -47.75*cm; //-(STRIPS_PER_PLANE-1)*STRIP_SPACING/2
double GlueXSensitiveDetectorFDC::CATHODE_ROT_ANGLE = 1.309; // radians (75 degrees)
double GlueXSensitiveDetectorFDC::STRIP_GAP = 0.1*cm;
double GlueXSensitiveDetectorFDC::STRIP_NODES = 3;

// Parameters for calculating the drift time-distance relation
double GlueXSensitiveDetectorFDC::LORENTZ_NR_PAR1;
double GlueXSensitiveDetectorFDC::LORENTZ_NR_PAR2;
double GlueXSensitiveDetectorFDC::LORENTZ_NZ_PAR1;
double GlueXSensitiveDetectorFDC::LORENTZ_NZ_PAR2;
double GlueXSensitiveDetectorFDC::DRIFT_RES_PARMS[3];
double GlueXSensitiveDetectorFDC::DRIFT_FUNC_PARMS[6];
double GlueXSensitiveDetectorFDC::DRIFT_BSCALE_PAR1;
double GlueXSensitiveDetectorFDC::DRIFT_BSCALE_PAR2;

// Parameters for estimating magnetic field drift effects
double GlueXSensitiveDetectorFDC::DIFFUSION_COEFF = 1.1e-6*cm*cm/s; // cm^2/s --> 200 microns at 1 cm
double GlueXSensitiveDetectorFDC::K2 = 1.15;

// Radius of deadened region around the beam
double GlueXSensitiveDetectorFDC::wire_dead_zone_radius[4] =
                                  {3.0*cm, 3.0*cm, 3.9*cm, 3.9*cm};
double GlueXSensitiveDetectorFDC::strip_dead_zone_radius[4] =
                                  {1.3*cm, 1.3*cm, 1.3*cm, 1.3*cm};

// Drift time - distance lookup table
int GlueXSensitiveDetectorFDC::drift_table_len;
double *GlueXSensitiveDetectorFDC::drift_table_t_ns;
double *GlueXSensitiveDetectorFDC::drift_table_d_cm;

int GlueXSensitiveDetectorFDC::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorFDC::fMutex = G4MUTEX_INITIALIZER;
int GlueXSensitiveDetectorFDC::fDrift_clusters = 0;

GlueXSensitiveDetectorFDC::GlueXSensitiveDetectorFDC(const G4String& name)
 : G4VSensitiveDetector(name),
   fWiresMap(0), fCathodesMap(0), fPointsMap(0)
{
   collectionName.insert("FDCWireHitsCollection");
   collectionName.insert("FDCCathodeHitsCollection");
   collectionName.insert("FDCPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the FDC, you must delete all old
   // objects of this class and create new ones.

   G4AutoLock barrier(&fMutex);
   if (instanceCount++ == 0) {
      int runno = HddmOutput::getRunNo();
      extern JApplication *japp;
      if (japp == 0) {
         G4cerr << "Error in GlueXSensitiveDetector constructor - "
                << "jana global DApplication object not set, "
                << "cannot continue." << G4endl;
         exit(-1);
      }
      JCalibration *jcalib = japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
      std::map<string, float> fdc_parms;
      jcalib->Get("FDC/fdc_parms", fdc_parms);
      DRIFT_SPEED = fdc_parms.at("FDC_DRIFT_SPEED")*cm/ns;
      ACTIVE_AREA_OUTER_RADIUS = fdc_parms.at("FDC_ACTIVE_AREA_OUTER_RADIUS")*cm;
      ANODE_CATHODE_SPACING = fdc_parms.at("FDC_ANODE_CATHODE_SPACING")*cm;
      TWO_HIT_TIME_RESOL = fdc_parms.at("FDC_TWO_HIT_RESOL")*ns;
      WIRES_PER_PLANE = fdc_parms.at("FDC_WIRES_PER_PLANE");
      WIRE_SPACING = fdc_parms.at("FDC_WIRE_SPACING")*cm;
      STRIPS_PER_PLANE = fdc_parms.at("FDC_STRIPS_PER_PLANE");
      STRIP_SPACING = fdc_parms.at("FDC_STRIP_SPACING")*cm;
      STRIP_GAP = fdc_parms.at("FDC_STRIP_GAP")*cm;
      MAX_HITS = fdc_parms.at("FDC_MAX_HITS");
      K2 = fdc_parms.at("FDC_K2");
      STRIP_NODES = fdc_parms.at("FDC_STRIP_NODES");
      THRESH_KEV = fdc_parms.at("FDC_THRESH_KEV");
      THRESH_STRIPS = fdc_parms.at("FDC_THRESH_STRIPS");
      DIFFUSION_COEFF = fdc_parms.at("FDC_DIFFUSION_COEFF")*cm*cm/s;
      U_OF_WIRE_ONE = -(WIRES_PER_PLANE -1) * WIRE_SPACING / 2;
      U_OF_STRIP_ONE = -(STRIPS_PER_PLANE -1) * STRIP_SPACING / 2;

      // Parameters for correcting for deflection due to Lorentz force
      std::map<string, double> lorentz_parms;
      jcalib->Get("FDC/lorentz_deflection_parms", lorentz_parms);
      LORENTZ_NR_PAR1 = lorentz_parms["nr_par1"];
      LORENTZ_NR_PAR2 = lorentz_parms["nr_par2"];
      LORENTZ_NZ_PAR1 = lorentz_parms["nz_par1"];
      LORENTZ_NZ_PAR2 = lorentz_parms["nz_par2"];

      // Parameters for accounting for variation in drift distance from FDC
      std::map<string, double> drift_res_parms;
      jcalib->Get("FDC/drift_resolution_parms", drift_res_parms); 
      DRIFT_RES_PARMS[0] = drift_res_parms["p0"];   
      DRIFT_RES_PARMS[1] = drift_res_parms["p1"];
      DRIFT_RES_PARMS[2] = drift_res_parms["p2"]; 

      // Time-to-distance function parameters for FDC
      std::map<string, double> drift_func_parms;
      jcalib->Get("FDC/drift_function_parms", drift_func_parms); 
      DRIFT_FUNC_PARMS[0] = drift_func_parms["p0"];   
      DRIFT_FUNC_PARMS[1] = drift_func_parms["p1"];
      DRIFT_FUNC_PARMS[2] = drift_func_parms["p2"]; 
      DRIFT_FUNC_PARMS[3] = drift_func_parms["p3"];
      DRIFT_FUNC_PARMS[4] = 1000.;
      DRIFT_FUNC_PARMS[5] = 0.;
      std::map<string, double> drift_func_ext;
      if (jcalib->Get("FDC/drift_function_ext", drift_func_ext) == false) {
         DRIFT_FUNC_PARMS[4] = drift_func_ext["p4"]; 
         DRIFT_FUNC_PARMS[5] = drift_func_ext["p5"]; 
      }

      // Factors for taking care of B-dependence of drift time for FDC
      std::map<string, double> fdc_drift_parms;
      jcalib->Get("FDC/fdc_drift_parms", fdc_drift_parms);
      DRIFT_BSCALE_PAR1 = fdc_drift_parms["bscale_par1"];
      DRIFT_BSCALE_PAR2 = fdc_drift_parms["bscale_par2"];

      // Build a lookup table of drift time->distance for the FDC,
      // used in the code to build an efficient reverse-map function.
      drift_table_len = 1000;
      drift_table_t_ns = new double[drift_table_len];
      drift_table_d_cm = new double[drift_table_len];
      double thigh = DRIFT_FUNC_PARMS[4];
      double tstep = 0.5; //ns
	  for (int j=0; j < drift_table_len; j++) {
	     double t = j * tstep;
	     if (t < thigh) {
	        double t2 = t*t;
            drift_table_t_ns[j] = t;
	        drift_table_d_cm[j] = DRIFT_FUNC_PARMS[0] * sqrt(t) +
                                  DRIFT_FUNC_PARMS[1] * t +
	                              DRIFT_FUNC_PARMS[2] * t2 +
                                  DRIFT_FUNC_PARMS[3] * t*t2;
	     }
	     else {
	        double thigh2 = thigh * thigh;
	        drift_table_t_ns[j] = t;
            drift_table_d_cm[j] = DRIFT_FUNC_PARMS[0] * sqrt(thigh) +
	                              DRIFT_FUNC_PARMS[1] * thigh +
	                              DRIFT_FUNC_PARMS[2] * thigh2 +
	                              DRIFT_FUNC_PARMS[3] * thigh2 * thigh +
	                              DRIFT_FUNC_PARMS[5] * (t - thigh);
	     }
      }

      G4cout << "FDC: ALL parameters loaded from ccdb" << G4endl;

      // Check for "driftclusters" option in control.in

      GlueXUserOptions *opts = GlueXUserOptions::GetInstance();
      if (opts) {
         std::map<int, int> driftclusters_opts;
         if (opts->Find("driftclusters", driftclusters_opts))
            fDrift_clusters = driftclusters_opts[1];
      }
   }
}

GlueXSensitiveDetectorFDC::GlueXSensitiveDetectorFDC(
                     const GlueXSensitiveDetectorFDC &src)
 : G4VSensitiveDetector(src),
   fWiresMap(src.fWiresMap),
   fCathodesMap(src.fCathodesMap),
   fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorFDC &GlueXSensitiveDetectorFDC::operator=(const
                                         GlueXSensitiveDetectorFDC &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fWiresMap = src.fWiresMap;
   fCathodesMap = src.fCathodesMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorFDC::~GlueXSensitiveDetectorFDC() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorFDC::Initialize(G4HCofThisEvent* hce)
{
   fWiresMap = new
               GlueXHitsMapFDCwire(SensitiveDetectorName, collectionName[0]);
   fCathodesMap = new
               GlueXHitsMapFDCcathode(SensitiveDetectorName, collectionName[1]);
   fPointsMap = new
               GlueXHitsMapFDCpoint(SensitiveDetectorName, collectionName[2]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fWiresMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fCathodesMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[2]), fPointsMap);
}

G4bool GlueXSensitiveDetectorFDC::ProcessHits(G4Step* step, 
                                              G4TouchableHistory* ROhist)
{
   double dEsum = step->GetTotalEnergyDeposit();
   if (dEsum == 0)
      return false;

   G4ThreeVector pin = step->GetPreStepPoint()->GetMomentum();
   G4ThreeVector xin = step->GetPreStepPoint()->GetPosition();
   G4ThreeVector xout = step->GetPostStepPoint()->GetPosition();
   double Ein = step->GetPreStepPoint()->GetTotalEnergy();
   double tin = step->GetPreStepPoint()->GetGlobalTime();
   double tout = step->GetPostStepPoint()->GetGlobalTime();

   G4Track *track = step->GetTrack();
   int trackID = track->GetTrackID();
   const G4VTouchable* touch = step->GetPreStepPoint()->GetTouchable();
   int package = GetIdent("package", touch);
   int layer = GetIdent("layer", touch);
   if (layer == 0) {
      fprintf(stderr, "hitFDC error: FDC layer number evaluates to zero! "
              "THIS SHOULD NEVER HAPPEN! drop this particle.\n");
      return false;
   }
   else if (package == 0) {
      fprintf(stderr, "hitFDC error: FDC package number evaluates to zero! "
              "THIS SHOULD NEVER HAPPEN! drop this particle.\n");
      return false;
   }
   // Normally numeric identifiers start at 1, eg. layer, package, module
   // but if it is an index counting from zero, add the "No" suffix.
   int packNo = package - 1;
   int module = 2 * packNo + ((layer - 1) / 3) + 1;
   int chamber = (module * 10) + ((layer - 1) % 3) + 1;

#ifdef MERGE_STEPS_BEFORE_HITS_GENERATION

   // This section forces hdgeant4 to conform to the behavior of hdgeant
   // in merging together all of the steps taken by a single track as it
   // moves through a single fdc chamber (plane) before generating hits.
   // For this purpose, it borrows existing class GlueXHitFDCpoint to
   // save the information from previous tracking steps by this track
   // inside the current fdc active volume, and merging the saved info
   // into the current step when the track either stops inside the fdc
   // active volume or when it exits through an exterior boundary.
   //
   // NOTE
   // This algorithm has a blind spot for tracks that stop inside the 
   // wire layer because the track disappears without ever crossing
   // through an exterior boundary. This condition produces a warning
   // message if verboseLevel > 0, but it is an inherent inefficiency
   // in this model, and should not be considered to be a bug.

   G4int rkey = 1 << 24; // reserved key for this thread
   std::map<int,GlueXHitFDCpoint*>::iterator saved_segment =
                                    fPointsMap->GetMap()->find(rkey);
   G4String invol = step->GetPostStepPoint()->GetTouchable()
                                            ->GetVolume()->GetName();
   if (track->GetTrackStatus() == fAlive &&
       (invol.index("FDA") == 0 || invol.index("FDX") == 0))
   {
      if (saved_segment == fPointsMap->GetMap()->end()) {
         GlueXHitFDCpoint segment(chamber);
         segment.track_ = trackID;
         segment.x_cm = xin[0];     // this is borrowed storage for
         segment.y_cm = xin[1];     // local track step information,
         segment.z_cm = xin[2];     // don't worry about the units!
         segment.t_ns = tin;
         segment.px_GeV = pin[0];
         segment.py_GeV = pin[1];
         segment.pz_GeV = pin[2];
         segment.E_GeV = Ein;
         segment.dEdx_GeV_cm = dEsum;
         fPointsMap->add(rkey, segment);
         return true;
      }
      else if (saved_segment->second->track_ == trackID) {
         saved_segment->second->dEdx_GeV_cm += dEsum;
         return true;
      }
      else {
         if (verboseLevel > 0)
            fprintf(stderr, "hitFDC warning: FDC saved track segment "
                            "was lost, drop this step.\n");
         fPointsMap->GetMap()->erase(saved_segment);
         return ProcessHits(step, ROhist);
      }
   }
   else if (saved_segment != fPointsMap->GetMap()->end()) {
      if (saved_segment->second->track_ == trackID) {
         xin.setX(saved_segment->second->x_cm);
         xin.setY(saved_segment->second->y_cm);
         xin.setZ(saved_segment->second->z_cm);
         tin = saved_segment->second->t_ns;
         pin.setX(saved_segment->second->px_GeV);
         pin.setY(saved_segment->second->py_GeV);
         pin.setZ(saved_segment->second->pz_GeV);
         Ein = saved_segment->second->E_GeV;
         dEsum += saved_segment->second->dEdx_GeV_cm;
         fPointsMap->GetMap()->erase(saved_segment);
      }
      else {
         if (verboseLevel > 0)
            fprintf(stderr, "hitFDC warning: FDC saved track segment "
                            "was orphaned, drop this step.\n");
         fPointsMap->GetMap()->erase(saved_segment);
         return ProcessHits(step, ROhist);
      }
   }
   
#endif

   G4ThreeVector x = (xin + xout) / 2;
   G4ThreeVector dx = xout - xin;
   double t = (tin + tout) / 2;
   double dr = dx.mag();
   double dEdx = (dr > 1e-3*cm)? dEsum/dr : 0;

   const G4AffineTransform &local_from_global = touch->GetHistory()
                                                     ->GetTopTransform();
   G4ThreeVector xinlocal = local_from_global.TransformPoint(xin);
   G4ThreeVector xoutlocal = local_from_global.TransformPoint(xout);

   // For particles that range out inside the active volume, the
   // "out" time may sometimes be set to something enormously high.
   // This screws up the hit. Check for this case here by looking
   // at tout and making sure it is less than 1 second. If it's
   // not, then just use tin for "t".

   if (tout > 1.0*s)
      t = tin;

   double alpha = atan2(xoutlocal[0] - xinlocal[0], xoutlocal[2] - xinlocal[2]);
   double sinalpha = sin(alpha);
   double cosalpha = cos(alpha);
   G4ThreeVector xlocal = (xinlocal + xoutlocal) / 2;

   // Make a fuzzy boundary around the forward dead region 
   // by killing any track segment whose midpoint is within the boundary

   if (xlocal.perp() < wire_dead_zone_radius[packNo])
      return false;

   int wire = ceil((xlocal[0] - U_OF_WIRE_ONE) / WIRE_SPACING + 0.5);
   double xwire = U_OF_WIRE_ONE + (wire - 1) * WIRE_SPACING;
   double uwire = xinlocal[2];
   double vwire = xinlocal[0] - xwire;
   double dradius = fabs(vwire * cosalpha - uwire * sinalpha);
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);

#if VERBOSE_PRINT_POINTS
   std::cout << "fdc in u,v=" << uwire << "," << vwire
             << " out u,v=" << xoutlocal[2] << "," << xoutlocal[0] - xwire
             << " dradius=" << dradius
             << " in=" << step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()
             << " out=" << step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetName()
             << std::endl;
#endif

   // Post the hit to the points list in the
   // order of appearance in the event simulation.
 
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (trackinfo->GetGlueXHistory() == 0) {
      G4int key = (chamber << 20) + fPointsMap->entries();
      GlueXHitFDCpoint* lastPoint = (*fPointsMap)[key - 1];
      // Limit fdc truthPoints to one per chamber
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          lastPoint->chamber_ != chamber)
      {
         GlueXHitFDCpoint newPoint(chamber);
         newPoint.primary_ = (track->GetParentID() == 0);
         newPoint.track_ = trackID;
         newPoint.x_cm = xout[0]/cm;
         newPoint.y_cm = xout[1]/cm;
         newPoint.z_cm = xout[2]/cm;
         newPoint.t_ns = tout/ns;
         newPoint.px_GeV = pin[0]/GeV;
         newPoint.py_GeV = pin[1]/GeV;
         newPoint.pz_GeV = pin[2]/GeV;
         newPoint.E_GeV = Ein/GeV;
         newPoint.dradius_cm = dradius/cm;
         newPoint.dEdx_GeV_cm = dEdx/(GeV/cm);
         newPoint.ptype_G3 = g3type;
         newPoint.trackID_ = itrack;
         fPointsMap->add(key, newPoint);
      }
   }

   // Post the hit to the hits tree, ordered by wire, strip number

   if (dEsum > 0) {
      double u0 = xinlocal[0];
      double u1 = xoutlocal[0];
      int wire1 = ceil((u0 - U_OF_WIRE_ONE) / WIRE_SPACING + 0.5);
      int wire2 = ceil((u1 - U_OF_WIRE_ONE) / WIRE_SPACING + 0.5);

      // Check that wire numbers are not out of range,
      // making sure at least one wire number is valid
      if (wire1 > WIRES_PER_PLANE && wire2 > WIRES_PER_PLANE)
         return false;
      else if (wire1 < 1 && wire2 < 1)
         return false;
      wire1 = (wire1 > WIRES_PER_PLANE)? WIRES_PER_PLANE :
              (wire1 < 1)? 1 : wire1;
      wire2 = (wire2 > WIRES_PER_PLANE)? WIRES_PER_PLANE :
              (wire2 < 1)? 1 : wire2;
      int dwire = (wire1 < wire2)? 1 : -1;

      // deal with the case of tracks crossing two cells
      for (int wire = wire1; wire != wire2 + dwire; wire += dwire) {
         double xwire = U_OF_WIRE_ONE + (wire - 1) * WIRE_SPACING;
         G4ThreeVector x0;
         G4ThreeVector x1;
         double dE;
         if (wire1 == wire2) {
            dE = dEsum;
            x0 = xinlocal;
            x1 = xoutlocal;
         }
         else {
            x0[0] = xwire - 0.5 * dwire * WIRE_SPACING;
            x0[1] = xinlocal[1] + (x0[0] - xinlocal[0] + 1e-20) *
                    (xoutlocal[1] - xinlocal[1]) / 
                    (xoutlocal[0] - xinlocal[0] + 1e-20);
            x0[2] = xinlocal[2] + (x0[0] - xinlocal[0] + 1e-20) *
                    (xoutlocal[2] - xinlocal[2]) /
                    (xoutlocal[0] - xinlocal[0] + 1e-20);
            if (fabs(x0[2] - xoutlocal[2]) > fabs(xinlocal[2] - xoutlocal[2]))
               x0 = xinlocal;

            x1[0] = xwire + 0.5 * dwire * WIRE_SPACING;
            x1[1] = xinlocal[1] + (x1[0] - xinlocal[0] + 1e-20) *
                    (xoutlocal[1] - xinlocal[1]) /
                    (xoutlocal[0] - xinlocal[0] + 1e-20);
            x1[2] = xinlocal[2] + (x1[0] - xinlocal[0] + 1e-20) *
                    (xoutlocal[2] - xinlocal[2]) /
                    (xoutlocal[0] - xinlocal[0] + 1e-20);
            if (fabs(x1[2] - xinlocal[2]) > fabs(xoutlocal[2] - xinlocal[2]))
               x1 = xoutlocal;

            dE = dEsum * (x1[2] - x0[2]) / 
                         (xoutlocal[2] - xinlocal[2] + 1e-20);
         }

         int key = GlueXHitFDCwire::GetKey(chamber, wire);
         GlueXHitFDCwire *anode = (*fWiresMap)[key];
         if (anode == 0) {
            GlueXHitFDCwire newanode(chamber, wire);
            fWiresMap->add(key, newanode);
            anode = (*fWiresMap)[key];
         }

         // Add the hit to the hits vector, maintaining track time ordering,
         // re-ordering according to hit times will take place at end of event.

         int merge_hits = 0;
         std::vector<GlueXHitFDCwire::hitinfo_t>::iterator hiter;
         for (hiter = anode->hits.begin(); hiter != anode->hits.end(); ++hiter) {
            if (itrack == hiter->itrack_ && fabs(tin - hiter->t1_ns*ns) < 0.01) {
               merge_hits = 1;
               break;
            }
            else if (hiter->t0_ns*ns > tin) {
               break;
            }
         }
         if (merge_hits) {
            hiter->dE_keV += dE/keV;
            hiter->t1_ns = tout/ns;
            hiter->x1_g = xout;
            hiter->x1_l = x1;
         }
         else {
            // create new hit
            hiter = anode->hits.insert(hiter, GlueXHitFDCwire::hitinfo_t());
            hiter->dE_keV = dE/keV;
            hiter->itrack_ = itrack;
            hiter->ptype_G3 = g3type;
            hiter->t0_ns = tin/ns;
            hiter->t1_ns = tout/ns;
            hiter->x0_g = xin;
            hiter->x1_g = xout;
            hiter->x0_l = x0;
            hiter->x1_l = x1;
         }
      }
   }
   return true;
}

void GlueXSensitiveDetectorFDC::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitFDCwire*> *wires = fWiresMap->GetMap();
   std::map<int,GlueXHitFDCcathode*> *strips = fCathodesMap->GetMap();
   std::map<int,GlueXHitFDCpoint*> *points = fPointsMap->GetMap();
   if (wires->size() == 0 && strips->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitFDCwire*>::iterator witer;
   std::map<int,GlueXHitFDCcathode*>::iterator siter;
   std::map<int,GlueXHitFDCpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << wires->size() << " anode wires with hits in the FDC: "
             << G4endl;
      for (witer = wires->begin(); witer != wires->end(); ++witer)
         witer->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << strips->size() << " cathode strips with hits in the FDC: "
             << G4endl;
      for (siter = strips->begin(); siter != strips->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the FDC: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorFDC::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getForwardDCs().size() == 0)
      hitview.addForwardDCs();
   hddm_s::ForwardDC &forwardDC = hitview.getForwardDC();

   // Collect and output the wireTruthHits

   for (witer = wires->begin(); witer != wires->end(); ++witer) {

      // Merge multiple segments from a single track into one, and 
      // apply the drift time algorithm to get a single hit time for each.

      int chamber = witer->second->chamber_;
      int wire = witer->second->wire_;
      int module = chamber / 10;
      int packNo = (module - 1) / 2;
      int layer = chamber % 10;
      int glayerNo = 3*(layer - 1) + module-1;
      int global_wire_number = 96 * glayerNo + wire - 1;
      std::vector<GlueXHitFDCwire::hitinfo_t> &splits = witer->second->hits;
      std::vector<GlueXHitFDCwire::hitinfo_t> hits;
      while (splits.size() > 0) {
         for (unsigned int ih=1; ih < splits.size(); ++ih) {
            if (fabs(splits[ih].t0_ns - splits[0].t1_ns) < 0.5*ns) {
               splits[0].dE_keV += splits[ih].dE_keV;
               splits[0].t1_ns = splits[ih].t1_ns;
               splits[0].x1_g = splits[ih].x1_g;
               splits[0].x1_l = splits[ih].x1_l;
               splits.erase(splits.begin() + ih);
               --ih;
            }
         }

         // Simulate the number of primary ion pairs.
         // The total number of ion pairs depends on the energy deposition 
         // and the effective average energy to produce a pair, w_eff.
         // On average for each primary ion pair produced there are n_s_per_p 
         // secondary ion pairs produced. 
    
         double xwire = U_OF_WIRE_ONE + (wire - 1) * WIRE_SPACING;
         double dE = splits[0].dE_keV*keV;
         if (dE > THRESH_KEV*keV) {
            // Average number of primary ion pairs
            double n_p_mean = dE / W_EFF_PER_ION / (1 + N_SECOND_PER_PRIMARY);
            // number of primary ion pairs
            int n_p = CLHEP::RandPoisson::shoot(n_p_mean); 
            G4ThreeVector &x0 = splits[0].x0_l;
            G4ThreeVector &x1 = splits[0].x1_l;
       
            if (fDrift_clusters == 0) {
               G4ThreeVector dx(x1 - x0);
               double alpha = -((x0[0] - xwire) * dx[0] + x0[2] * dx[2]) /
                               (dx[0] * dx[0] + dx[2] * dx[2] + 1e-99);
               alpha = (alpha < 0)? 0 : (alpha > 1)? 1 : alpha;
               double t = (splits[0].t0_ns + splits[0].t1_ns)*ns / 2;
               G4ThreeVector x((splits[0].x0_g + splits[0].x1_g) / 2);
               G4ThreeVector xlocal(x0 + alpha * dx);
               double tdrift;
               int wire_fired = add_anode_hit(hits, splits[0], layer,
                                              xwire, x, xlocal, dE,
                                              t, tdrift);
               if (wire_fired) {
                  add_cathode_hit(splits[0], packNo, xwire, xlocal[1],
                                  tdrift, n_p, chamber, module, layer,
                                  global_wire_number);
               }
            }
            else {
               G4ThreeVector xlocal;
               G4ThreeVector x((splits[0].x0_g + splits[0].x1_g) / 2);
               double t = (splits[0].t0_ns + splits[0].t1_ns)*ns / 2;
               // Loop over the number of primary ion pairs
               for (int n=0; n < n_p; n++) {
                  // Generate a cluster at a random position
                  // along the path within the cell
                  double u = G4UniformRand();
                  xlocal = x0 + u * (x1 - x0);
                  double tdrift;
                  int wire_fired = add_anode_hit(hits, splits[0], layer,
                                                 xwire, x, xlocal, dE,
                                                 t, tdrift);
                  if (wire_fired) {
                     add_cathode_hit(splits[0], packNo, xwire, xlocal[1],
                                     tdrift, n_p, chamber, module, layer,
                                     global_wire_number);
                  }
               }
            }
         }
         splits.erase(splits.begin());
      }

      if (fDrift_clusters) {
         // store waveform data in sampled sequence with 1 ns bins
         int num_samples = (int)FDC_TIME_WINDOW;
         double *samples = new double[num_samples];
         for (int i=0; i < num_samples; i++) {
            samples[i] = fdc_wire_signal_mV(double(i), hits);
         }
 
         // take the earliest hit to identify the track parameters
         double dradius_cm = hits[0].d_cm;
         int ptype_G3 = hits[0].ptype_G3;
         int itrack = hits[0].itrack_;
         double t0_ns = hits[0].t0_ns;
         double t1_ns = hits[0].t1_ns;
         G4ThreeVector x0_g = hits[0].x0_g;
         G4ThreeVector x0_l = hits[0].x0_l;
         G4ThreeVector x1_g = hits[0].x1_g;
         G4ThreeVector x1_l = hits[0].x1_l;

         double dE_keV = 0;
         int over_threshold = 0;
         for (int i=0; i < num_samples; i++) {
            if (samples[i] > THRESH_ANODE) {
               if (!over_threshold) {
                  splits.push_back(GlueXHitFDCwire::hitinfo_t());
                  splits.back().d_cm = dradius_cm;
                  splits.back().ptype_G3 = ptype_G3;
                  splits.back().itrack_ = itrack;
                  splits.back().t_ns = double(i);
                  splits.back().t0_ns = t0_ns;
                  splits.back().t1_ns = t1_ns;
                  splits.back().x0_g = x0_g;
                  splits.back().x0_l = x0_l;
                  splits.back().x1_g = x1_g;
                  splits.back().x1_l = x1_l;
         
                  // Do an interpolation to find the time 
                  // at which the threshold was crossed.
                  double t_array[4];
                  double t_ns, t_err;
                  for (int k=0; k < 4; k++)
                     t_array[k] = i - 1 + k;
                  polint(&samples[i-1], t_array, 4,
                         THRESH_ANODE, &t_ns, &t_err);
                  splits.back().t_ns = t_ns;
                  over_threshold = 1; 
               }
               dE_keV += samples[i];
            }
            else if (over_threshold) {
               splits.back().dE_keV = dE_keV;
               over_threshold = 0;
               dE_keV = 0;
            }
         }
         delete [] samples;
      }
      else {

         // Build new reduced hit list ordered by hit time
 
         std::vector<GlueXHitFDCwire::hitinfo_t>::iterator hiter;
         for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
            int merge_hits = 0;
            std::vector<GlueXHitFDCwire::hitinfo_t>::iterator siter;
            for (siter = splits.begin(); siter != splits.end(); ++siter) {
               if (fabs(hiter->t_ns - siter->t_ns) < TWO_HIT_TIME_RESOL) {
                  merge_hits = 1;
                  break;
               }
               else if (hiter->t_ns < siter->t_ns) {
                  break;
               }
            }
            if (merge_hits) {
               if (hiter->t_ns < siter->t_ns) {
                  double dE_keV = siter->dE_keV;
                  *siter = *hiter;
                  siter->dE_keV += dE_keV;
               }
               else {
                  siter->dE_keV += hiter->dE_keV;
               }
            }
            else {
               splits.insert(siter, *hiter);
            }
         }
      }

      if (splits.size() > 0) {
         hddm_s::FdcChamberList chambers = forwardDC.getFdcChambers();
         hddm_s::FdcChamberList::iterator citer;
         for (citer = chambers.begin(); citer != chambers.end(); ++citer) {
            if (citer->getModule() == module && citer->getLayer() == layer)
               break;
         }
         if (citer == chambers.end()) {
            chambers = forwardDC.addFdcChambers(1);
            chambers(0).setModule(module);
            chambers(0).setLayer(layer);
            citer = chambers.begin();
         }
         hddm_s::FdcAnodeWireList anodes = citer->getFdcAnodeWires();
         hddm_s::FdcAnodeWireList::iterator aiter;
         for (aiter = anodes.begin(); aiter != anodes.end(); ++aiter) {
            if (aiter->getWire() == wire)
               break;
         }
         if (aiter == anodes.end()) {
            anodes = citer->addFdcAnodeWires(1);
            anodes(0).setWire(wire);
            aiter = anodes.begin();
         }
         int hitscount = splits.size();
         if (hitscount + aiter->getFdcAnodeTruthHits().size() > MAX_HITS) {
            hitscount = MAX_HITS - aiter->getFdcAnodeTruthHits().size();
            G4cerr << "GlueXSensitiveDetectorFDC::EndOfEvent warning: "
                   << "wire hit count exceeds max hit count " << MAX_HITS
                   << ", " << splits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::FdcAnodeTruthHitList thit = aiter->addFdcAnodeTruthHits(1);
            thit(0).setDE(splits[ih].dE_keV * 1e-6);
            thit(0).setT(splits[ih].t_ns);
            thit(0).setD(splits[ih].d_cm);
            thit(0).setItrack(splits[ih].itrack_);
            thit(0).setPtype(splits[ih].ptype_G3);
            thit(0).setT_unsmeared(splits[ih].t_unsmeared_ns);
         }
      }
   }

   // Collect and output the cathodeTruthHits

   for (siter = strips->begin(); siter != strips->end(); ++siter) {
      int chamber = siter->second->chamber_;
      int planeNo = siter->second->plane_;
      int stripNo = siter->second->strip_;
      int module = chamber / 10;
      int layer = chamber % 10;
      std::vector<GlueXHitFDCcathode::hitinfo_t> &hits = siter->second->hits;
      std::vector<GlueXHitFDCcathode::hitinfo_t>::iterator hiter;
      if (fDrift_clusters) {
         // store waveform data in sampled sequence with 1 ns bins
         int num_samples = (int)FDC_TIME_WINDOW;
         double *samples = new double[num_samples];
         for (int i=0; i < num_samples; i++) {
            samples[i] = fdc_cathode_signal_mV(double(i), hits);
         }

         // take the earliest hit to identify the track parameters
         int ptype_G3 = hits[0].ptype_G3;
         int itrack = hits[0].itrack_;
         double u_cm = hits[0].u_cm;
         double v_cm = hits[0].v_cm;

         hits.clear();
         int istart = 0;
         int threshold_toggle = 0;
         int FADC_BIN_SIZE = 1;
         for (int i=0; i < num_samples; i += FADC_BIN_SIZE) {
            if (samples[i] > THRESH_STRIPS) {
               if (threshold_toggle == 0) {
                  hits.push_back(GlueXHitFDCcathode::hitinfo_t());
                  hits.back().t_ns = i;
                  hits.back().ptype_G3 = ptype_G3;
                  hits.back().itrack_ = itrack;
                  hits.back().v_cm = v_cm;
                  hits.back().u_cm = u_cm;
                  istart = (i > 0)? i-1 : 0;
                  threshold_toggle=1;
               }
            }
            else if (threshold_toggle) {
               // Find the first peak
               for (int j = istart + 1; j < i - 1; j++) {
                  if (samples[j] > samples[j-1] && samples[j] > samples[j+1]) {
                     hits.back().q_fC = samples[j];
                     break;
                  }
               } 
               threshold_toggle=0; 
            }
         }
         int i = num_samples - 1;
         if (samples[i] > THRESH_STRIPS && threshold_toggle) {
            for (int j = istart + 1; j < i - 1; j++) {
               if (samples[j] > samples[j-1] && samples[j] > samples[j+1]) {
                  hits.back().q_fC = samples[j];
                  break;
               }
            }
         }
        delete [] samples;
      }
      else {
         double t_ns = -1e9;
         for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
            if (hiter == hits.begin())
               continue;
            // combine separate hits that are too close together in time
            if (fabs(hiter->t_ns - t_ns) < TWO_HIT_TIME_RESOL || 
                hiter->q_fC == 0)
            {
               // Use the time from the earlier hit but add the charge
               (hiter - 1)->q_fC += hiter->q_fC;
               hits.erase(hiter);
               hiter = hits.begin();
            }
            else {
               t_ns = hiter->t_ns;
            }
         }
      }

      if (hits.size() > 0) {
         hddm_s::FdcChamberList chambers = forwardDC.getFdcChambers();
         hddm_s::FdcChamberList::iterator citer;
         for (citer = chambers.begin(); citer != chambers.end(); ++citer) {
            if (citer->getModule() == module && citer->getLayer() == layer)
               break;
         }
         if (citer == chambers.end()) {
            chambers = forwardDC.addFdcChambers(1);
            chambers(0).setModule(module);
            chambers(0).setLayer(layer);
            citer = chambers.begin();
         }
         hddm_s::FdcCathodeStripList cathodes = citer->getFdcCathodeStrips();
         hddm_s::FdcCathodeStripList::iterator kiter;
         for (kiter = cathodes.begin(); kiter != cathodes.end(); ++kiter) {
            if (kiter->getPlane() == planeNo && kiter->getStrip() == stripNo)
               break;
         }
         if (kiter == cathodes.end()) {
            cathodes = citer->addFdcCathodeStrips(1);
            cathodes(0).setPlane(planeNo);
            cathodes(0).setStrip(stripNo);
            kiter = cathodes.begin();
         }
         int hitscount = hits.size();
         if (hitscount + kiter->getFdcCathodeTruthHits().size() > MAX_HITS) {
            hitscount = MAX_HITS - kiter->getFdcCathodeTruthHits().size();
            G4cerr << "GlueXSensitiveDetectorFDC::EndOfEvent warning: "
                   << "cathode hit count exceeds max hit count " << MAX_HITS
                   << ", " << hits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::FdcCathodeTruthHitList thit = kiter->addFdcCathodeTruthHits(1);
            thit(0).setQ(hits[ih].q_fC);
            thit(0).setT(hits[ih].t_ns);
            thit(0).setItrack(hits[ih].itrack_);
            thit(0).setPtype(hits[ih].ptype_G3);
         }
      }
   }

   // Collect and output the fdcTruthPoints

   GlueXUserEventInformation* eventinfo = (GlueXUserEventInformation*)info;
   int last_chamber = -1;
   hddm_s::FdcTruthPoint *last_point = 0;
   for (piter = points->begin(); piter != points->end(); ++piter) {
      if (piter->second->chamber_ == last_chamber &&
          piter->second->track_ == last_point->getTrack() &&
          fabs(piter->second->t_ns - last_point->getT()) < 0.1)
      {
         if (piter->second->dradius_cm < last_point->getDradius())
            last_point->setDradius(piter->second->dradius_cm);
         else
            piter->second->dradius_cm = last_point->getDradius();
         if (piter->second->t_ns > last_point->getT())
            continue;
      }
      int module = piter->second->chamber_ / 10;
      int layer = piter->second->chamber_ % 10;
      hddm_s::FdcChamberList chambers = forwardDC.getFdcChambers();
      hddm_s::FdcChamberList::iterator citer;
      for (citer = chambers.begin(); citer != chambers.end(); ++citer) {
         if (citer->getModule() == module && citer->getLayer() == layer)
            break;
      }
      if (citer == chambers.end()) {
         chambers = forwardDC.addFdcChambers(1);
         chambers(0).setModule(module);
         chambers(0).setLayer(layer);
         citer = chambers.begin();
      }
      hddm_s::FdcTruthPointList point = citer->addFdcTruthPoints(1);
      point(0).setE(piter->second->E_GeV);
      point(0).setDEdx(piter->second->dEdx_GeV_cm);
      point(0).setDradius(piter->second->dradius_cm);
      point(0).setPrimary(piter->second->primary_);
      point(0).setPtype(piter->second->ptype_G3);
      point(0).setPx(piter->second->px_GeV);
      point(0).setPy(piter->second->py_GeV);
      point(0).setPz(piter->second->pz_GeV);
      point(0).setT(piter->second->t_ns);
      point(0).setX(piter->second->x_cm);
      point(0).setY(piter->second->y_cm);
      point(0).setZ(piter->second->z_cm);
      point(0).setTrack(piter->second->track_);
      hddm_s::TrackIDList tid = point(0).addTrackIDs();
      tid(0).setItrack(piter->second->trackID_);
      int itrack = piter->second->track_;
      const GlueXUserEventInformation::parent_history_t *parent_history;
      while ((parent_history = eventinfo->GetParentHistory(itrack))) {
         hddm_s::TrackOriginList origin = point(0).addTrackOrigins(1);
	 origin(0).setItrack(parent_history->parent_id);
	 origin(0).setPtype(parent_history->g3type);
	 origin(0).setX(parent_history->x0[0]/cm);
	 origin(0).setY(parent_history->x0[1]/cm);
	 origin(0).setZ(parent_history->x0[2]/cm);
	 origin(0).setT(parent_history->t0/ns);
	 itrack = parent_history->parent_id;
      }
      last_chamber = piter->second->chamber_;
      last_point = &point(0);
   }
}

double GlueXSensitiveDetectorFDC::asic_response(double t_ns)
{
   // Simulation of the ASIC response to a pulse due to a cluster

   double par[11] = {-0.01986, 0.01802, -0.001097, 10.3, 11.72,
                     -0.03701, 35.84, 15.93, 0.006141, 80.95, 24.77};
   if (t_ns < par[3])
      return par[0] * t_ns + par[1] * t_ns*t_ns + par[2] * t_ns*t_ns*t_ns;
   else
      return (par[0] * par[3] +
              par[1] * par[3]*par[3] +
              par[2] * par[3]*par[3]*par[3]) *
             exp(-pow((t_ns - par[3])/par[4], 2)) +
             par[5] * exp(-pow((t_ns - par[6])/par[7], 2)) +
             par[8] * exp(-pow((t_ns - par[9])/par[10], 2));
}

double GlueXSensitiveDetectorFDC::fdc_wire_signal_mV(
       double t_ns, std::vector<GlueXHitFDCwire::hitinfo_t> &hits)
{
   // Simulation of signal on a wire

   double asic_gain = 0.76; // mV/fC
   double signal_mV = 0;
   std::vector<GlueXHitFDCwire::hitinfo_t>::iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      if (t_ns > hiter->t_ns) {
         double my_time_ns = t_ns - hiter->t_ns;
         signal_mV += asic_gain * hiter->dE_keV * asic_response(my_time_ns);
      }
   }
   return signal_mV;
}

double GlueXSensitiveDetectorFDC::fdc_cathode_signal_mV(
       double t_ns, std::vector<GlueXHitFDCcathode::hitinfo_t> &hits)
{
   // Simulation of signal on a cathode strip (ASIC output)

   double asic_gain = 2.3; // mv/fC
   double signal_mV = 0;
   std::vector<GlueXHitFDCcathode::hitinfo_t>::iterator hiter;
   for (hiter =hits.begin(); hiter != hits.end(); ++hiter) {
      if (t_ns > hiter->t_ns) {
         double my_time_ns = t_ns - hiter->t_ns;
         signal_mV += asic_gain * hiter->q_fC * asic_response(my_time_ns);
      }
   }
   return signal_mV;
}

void GlueXSensitiveDetectorFDC::add_cathode_hit(
                                GlueXHitFDCwire::hitinfo_t &wirehit,
                                int packageNo,
                                double xwire, 
                                double yavalanche, 
                                double tdrift,
                                int n_p,
                                int chamber, 
                                int module,
                                int layer, 
                                int global_wire_number)
{
   double q_anode;
   if (!fDrift_clusters) {
      // Total number of ion pairs.  On average for each primary ion 
      // pair produced there are n_s secondary ion pairs produced.  The
      // probability distribution is a compound poisson distribution
      // that requires generating two Poisson variables.
      double n_s_mean = n_p * N_SECOND_PER_PRIMARY;
      int n_s = CLHEP::RandPoisson::shoot(n_s_mean);
      q_anode = (n_s + n_p) * GAS_GAIN * ELECTRON_CHARGE;
   }
   else {
      // Distribute the number of secondary ionizations for this primary
      // ionization according to a Poisson distribution with mean n_s_over_p.
      // For simplicity we assume these secondary electrons and the primary
      // electron stay together as a cluster.
      int n_s = CLHEP::RandPoisson::shoot(N_SECOND_PER_PRIMARY);
      q_anode = (1 + n_s) * GAS_GAIN * ELECTRON_CHARGE;
   }

   // Mock-up of cathode strip charge distribution 
   for (int plane=1; plane < 4; plane += 2) {
      double theta = (plane == 1)?  M_PI-CATHODE_ROT_ANGLE : CATHODE_ROT_ANGLE;
      double cathode_u = -xwire * cos(theta) - yavalanche * sin(theta);
      int strip1 = ceil((cathode_u - U_OF_STRIP_ONE) / STRIP_SPACING + 0.5);
      double cathode_u1 = (strip1 - 1) * STRIP_SPACING + U_OF_STRIP_ONE;
      double delta_u = cathode_u - cathode_u1;
      for (int node = -STRIP_NODES; node <= STRIP_NODES; node++) {
         // Induce charge on the strips according to the Mathieson 
         // function tuned to results from FDC prototype
         double lambda1 = ((node - 0.5) * STRIP_SPACING +
                           STRIP_GAP / 2. - delta_u) / ANODE_CATHODE_SPACING;
         double lambda2 = ((node + 0.5) * STRIP_SPACING - 
                           STRIP_GAP / 2. - delta_u) / ANODE_CATHODE_SPACING;
         double factor = 0.25 * M_PI * K2;
         double q = 0.25 * q_anode * (tanh(factor * lambda2) - 
                                      tanh(factor * lambda1));
         int strip = strip1 + node;
         // Throw away hits on strips falling within a certain dead-zone radius
         double strip_outer_u = cathode_u1;
         strip_outer_u  += node * (STRIP_SPACING + STRIP_GAP / 2.);
         double cathode_v = -xwire * sin(theta) + yavalanche * cos(theta);
         double check_radius = sqrt(strip_outer_u * strip_outer_u +
                                    cathode_v * cathode_v);
         if (strip > 0 && strip <= STRIPS_PER_PLANE &&
             check_radius > strip_dead_zone_radius[packageNo])
         {
            int key = GlueXHitFDCcathode::GetKey(chamber, plane, strip);
            GlueXHitFDCcathode *cathode = (*fCathodesMap)[key];
            if (cathode == 0) {
               GlueXHitFDCcathode newcathode(chamber, plane, strip);
               fCathodesMap->add(key, newcathode);
               cathode = (*fCathodesMap)[key];
            }
            std::vector<GlueXHitFDCcathode::hitinfo_t>::iterator hiter;
            for (hiter = cathode->hits.begin();
                 hiter != cathode->hits.end(); ++hiter)
            {
               if (hiter->t_ns*ns > tdrift)
                  break;
            }
            // create new hit
            hiter = cathode->hits.insert(hiter, GlueXHitFDCcathode::hitinfo_t());
            hiter->t_ns = tdrift/ns;
            hiter->q_fC = q/fC;
            hiter->itrack_ = wirehit.itrack_;
            hiter->ptype_G3 = wirehit.ptype_G3;
            hiter->u_cm = xwire/cm;
            hiter->v_cm = yavalanche/cm;
         }
      } // loop over cathode strips
   } // loop over cathode planes
}

int GlueXSensitiveDetectorFDC::add_anode_hit(
                               std::vector<GlueXHitFDCwire::hitinfo_t> &hits,
                               GlueXHitFDCwire::hitinfo_t &wirehit,
                               int layer, 
                               double xwire,
                               G4ThreeVector &xglobal,
                               G4ThreeVector &xlocal,
                               double dE, 
                               double t,
                               double &tdrift)
{
   // Get the magnetic field at this cluster position        
   G4ThreeVector B = GlueXDetectorConstruction::GetInstance()
                     ->GetMagneticField(xglobal, tesla);
   double BrhoT = B.perp();
  
   // Find the angle between the wire direction and the direction of the
   // magnetic field in the x-y plane
   double wire_theta = 1.0472 * ((layer % 3) - 1);
   double wire_dir[2] = {sin(wire_theta), cos(wire_theta)};
   double phi = 0;
   if (BrhoT > 0)
      phi = acos((B[0] * wire_dir[0] + B[1] * wire_dir[1]) / BrhoT);

   // useful combinations of dx and dz
   double dx = xlocal[0] - xwire;
   double dx2 = dx * dx;
   double dx4 = dx2 * dx2;
   double dz = xlocal[2];
   double dz2 = dz * dz;
   double dz4 = dz2 * dz2;

   // Next compute the avalanche position along wire.  
   // Correct avalanche position with deflection along wire
   // due to the Lorentz force.
   double cm2 = cm * cm;
   double cm4 = cm2 * cm2;
   xlocal[1] += (LORENTZ_NR_PAR1 * B[2] * (1 + LORENTZ_NR_PAR2 * BrhoT)) * dx +
                (LORENTZ_NZ_PAR1 + LORENTZ_NZ_PAR2 * B[2]) * BrhoT * cos(phi) * xlocal[2] +
                (-0.000176 * dx * dx2 / (dz2 + 0.001*cm2));
   // Add transverse diffusion
   xlocal[1] += G4RandGauss::shoot() *
                (0.01*cm * pow((dx2 + dz2)/cm2, 0.125) + 0.0061*cm * dx2/cm2);

   // Do not use this cluster if the Lorentz force would deflect 
   // the electrons outside the active region of the detector
   if (sqrt(xlocal[1] * xlocal[1] + xwire * xwire) > ACTIVE_AREA_OUTER_RADIUS) 
      return 0;

   // Model the drift time and longitudinal diffusion as a function of 
   // position of the cluster within the cell            

#if OLD_FDC_DRIFT_TIME_MODEL
   double tdrift_unsmeared = 1086.0*ns * (1 + 0.039 * B.mag()) * dx2/cm2 +
                             1068.0*ns * dz2/cm2 + 
                             (-2.675*ns / (dz2/cm2 + 0.001) +
                               2.4e4*ns * dz2/cm2) * dx4/cm4;
#else
    double dradius = sqrt(dx2 + dz2);
    int index = locate(drift_table_d_cm, drift_table_len, dradius/cm);
    index = (index < drift_table_len - 3)? index : drift_table_len - 3;
    double *dd = &drift_table_d_cm[index];
    double tt = 0.5; //ns
    double dd10 = dd[1] - dd[0];
    double dd20 = dd[2] - dd[0];
    double dd21 = dd[2] - dd[1];
    double qa = tt*index;
    double qb = (dd20/dd10 - 2*dd10/dd20) * tt/dd21;
    double qc = (2/dd20 - 1/dd10) * tt/dd21;
    double d0 = dradius/cm - dd[0];
    double tdrift_unsmeared = qa + qb*d0 + qc*d0*d0;
#endif

   // Apply small B-field dependence on the drift time
   tdrift_unsmeared *= 1. + DRIFT_BSCALE_PAR1 + DRIFT_BSCALE_PAR2*B[2]*B[2];

   // Minimum drift time for docas near wire (very crude approximation)
   double v_max = 0.08*cm/ns; // guess for now based on Garfield, near wire 
   double tmin = dradius / v_max;

   // longitidinal diffusion, derived from Garfield calculations
   double dt = (G4RandGauss::shoot() - 0.5) *
               (39.44*ns * dx4/cm4 / (0.5 - dz2/cm2) + 
                56.00*ns * dz4/cm4 / (0.5 - dx2/cm2) +
                0.01566*ns * dx4/cm4 / (dz4/cm4 + 0.002) / (0.251 - dx2/cm2));

   double tdrift_smeared = tdrift_unsmeared + dt;
   if (tdrift_smeared < tmin) {
      tdrift_smeared = tmin;
   }

   // Avalanche time
   tdrift = t + tdrift_smeared;

   // Skip cluster if the time would go beyond readout window
   if (tdrift > FDC_TIME_WINDOW)
      return 0;

   // Record the anode hit
   std::vector<GlueXHitFDCwire::hitinfo_t>::iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      if (hiter->t_ns*ns > tdrift) {
         break;
      }
   }
   hiter = hits.insert(hiter, wirehit);
   hiter->dE_keV = dE/keV;
   hiter->t_ns = tdrift/ns;
   hiter->t_unsmeared_ns = tdrift_unsmeared/ns;
   hiter->d_cm = dradius/cm;
   return 1;
}

void GlueXSensitiveDetectorFDC::polint(double *xa, double *ya, int n,
                                       double x, double *y, double *dy)
{
   // Slightly modified versions of the "polint" routine from
   // Press, William H., Brian P. Flannery, Saul A. Teukolsky and
   // William T. Vetterling, 1986, "Numerical Recipes: The Art of
   // Scientific Computing" (Fortran), Cambrigde University Press,
   // pp. 80-82.

   double *c = NULL;
   double *d = NULL;
   double den;
   double dif;
   double dift;
   double ho;
   double hp;
   double w;

   int i;
   int m;
   int ns;

   if ((c = (double*)malloc(n*sizeof(double))) == NULL ||
       (d = (double*)malloc(n*sizeof(double))) == NULL )
   {
      fprintf(stderr, "polint error: allocating workspace\n" );
      fprintf(stderr, "polint error: setting y = 0 and dy = 1e9\n" );
      *y = 0.0;
      *dy = 1.e9;
      if (c != NULL)
         free(c);
      if (d != NULL )
         free(d);
      return;
   }

   ns = 0;
   dif = fabs(x-xa[0]);
   for (i = 0; i < n; ++i) {
       dift = fabs(x-xa[i]);
       if (dift < dif) {
           ns = i;
           dif = dift;
       }
       c[i] = ya[i];
       d[i] = ya[i];
   }
   *y = ya[ns];
   ns = ns-1;
   for (m = 0; m < n-1; ++m) {
      for (i = 0; i < n-m-1; ++i) {
         ho = xa[i]-x;
         hp = xa[i+m+1]-x;
         w = c[i+1]-d[i];
         den = ho-hp;
         if (den == 0) {
            fprintf( stderr, "polint error: den = 0\n" );
            fprintf( stderr, "polint error: setting y = 0 and dy = 1e9\n" );
            *y = 0.0;
            *dy = 1.e9;
            if (c != NULL)
               free(c);
            if (d != NULL)
               free(d);
            return;
         }
         den = w/den;
         d[i] = hp*den;
         c[i] = ho*den;
      }
      if (2*(ns+1) < n-m-1) {
         *dy = c[ns+1];
      }
      else {
         *dy = d[ns];
         ns = ns-1;
      }
      *y = (*y)+(*dy);
   }

   if (c != NULL)
      free(c);
   if (d != NULL)
      free(d);
}

int GlueXSensitiveDetectorFDC::locate(double *xx, int n, double x)
{
   // Locate a position in array xx of value x

   int j;
   int ju;
   int jm;
   int jl;
   int ascnd;

   jl = -1;
   ju = n;
   ascnd = (xx[n-1] >= xx[0]);
   while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if ((x >= xx[jm]) == ascnd)
         jl = jm;
      else
         ju = jm;
   }
   if (x == xx[0])
      j = 0;
   else if (x == xx[n-1])
      j = n-2;
   else
      j = jl; 
   return j;
}

int GlueXSensitiveDetectorFDC::GetIdent(std::string div, 
                                        const G4VTouchable *touch)
{
   const HddsG4Builder* bldr = GlueXDetectorConstruction::GetBuilder();
   std::map<std::string, std::vector<int> >::const_iterator iter;
   std::map<std::string, std::vector<int> > *identifiers;
   int max_depth = touch->GetHistoryDepth();
   for (int depth = 0; depth < max_depth; ++depth) {
      G4VPhysicalVolume *pvol = touch->GetVolume(depth);
      G4LogicalVolume *lvol = pvol->GetLogicalVolume();
      int volId = fVolumeTable[lvol];
      if (volId == 0) {
         volId = bldr->getVolumeId(lvol);
         fVolumeTable[lvol] = volId;
      }
      identifiers = &Refsys::fIdentifierTable[volId];
      if ((iter = identifiers->find(div)) != identifiers->end()) {
         int copyNum = touch->GetCopyNumber(depth);
         copyNum += (dynamic_cast<G4PVPlacement*>(pvol))? -1 : 0;
         return iter->second[copyNum];
      }
   }
   return -1;
}

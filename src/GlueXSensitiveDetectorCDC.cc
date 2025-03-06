//
// GlueXSensitiveDetectorCDC - class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 28, 2015

#include "GlueXSensitiveDetectorCDC.hh"
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

const double fC = 1e-15 * coulomb;
const double GlueXSensitiveDetectorCDC::ELECTRON_CHARGE = 1.6022e-4*fC;

// Drift speed 2.2cm/us is appropriate for a 90/10 Argon/Methane mixture
double GlueXSensitiveDetectorCDC::DRIFT_SPEED = 0.0055*cm/ns;

// Minimum hit time difference for two hits on the same wire
double GlueXSensitiveDetectorCDC::TWO_HIT_TIME_RESOL = 25*ns;

// Cutoff on the number of allowed straw hits per straw
int GlueXSensitiveDetectorCDC::MAX_HITS = 1000;

// Minimum energy deposition for a straw hit (keV and mV)
double GlueXSensitiveDetectorCDC::THRESH_KEV = 1.;
double GlueXSensitiveDetectorCDC::THRESH_MV = 1.;

// Drift distance can be somewhat larger than this in a field
double GlueXSensitiveDetectorCDC::STRAW_RADIUS = 0.776*cm;

// Straw hits are accepted from t=0 up to this maximum time
double GlueXSensitiveDetectorCDC::CDC_TIME_WINDOW = 1000*ns;

// Parameters for setting signal pulse height
double GlueXSensitiveDetectorCDC::GAS_GAIN = 1e5;

// Average number of secondary ion pairs for 50/50 Ar/CO2 mixture
double GlueXSensitiveDetectorCDC::N_SECOND_PER_PRIMARY = 1.94; 

// Average energy needed to produce an ion pair for 50/50 mixture
double GlueXSensitiveDetectorCDC::W_EFF_PER_ION = 29.5*eV;

// These shared constants/tables initialized once from ccdb at run-time.
// It is assumed that all threads share common values of all of these
// lookup tables, and as long as all threads are working on events
// belonging to the same run, this should be true.

int GlueXSensitiveDetectorCDC::instanceCount = 0;
int GlueXSensitiveDetectorCDC::fDrift_clusters = 0;
double GlueXSensitiveDetectorCDC::fDrift_time[CDC_DRIFT_TABLE_LEN];
double GlueXSensitiveDetectorCDC::fDrift_distance[CDC_DRIFT_TABLE_LEN];
double GlueXSensitiveDetectorCDC::fBscale_par1;
double GlueXSensitiveDetectorCDC::fBscale_par2;

G4Mutex GlueXSensitiveDetectorCDC::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorCDC::GlueXSensitiveDetectorCDC(const G4String& name)
 : G4VSensitiveDetector(name),
   fStrawsMap(0), fPointsMap(0)
{
   collectionName.insert("CDCStrawHitsCollection");
   collectionName.insert("CDCPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the CDC, you must delete all old
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
      std::map<string, float> cdc_parms;
      jcalib->Get("CDC/cdc_parms", cdc_parms);
      DRIFT_SPEED = cdc_parms.at("CDC_DRIFT_SPEED")*cm/ns;
      TWO_HIT_TIME_RESOL = cdc_parms.at("CDC_TWO_HIT_RESOL")*ns;
      MAX_HITS = cdc_parms.at("CDC_MAX_HITS");
      THRESH_KEV = cdc_parms.at("CDC_THRESH_KEV");

      // Decide which drift time tables to load based on
      // whether or not a magnetic field is present.
 
      G4ThreeVector B = GlueXDetectorConstruction::GetInstance()
                        ->GetMagneticField(G4ThreeVector(0, 0, 65*cm), tesla);
      if (B.mag() > 1e-3) {
         int nvalues = CDC_DRIFT_TABLE_LEN;
         std::vector< std::map<std::string, float> > values;
         jcalib->Get("CDC/cdc_drift_table", values);
         for (int k=0; k < nvalues; ++k) {
            fDrift_distance[k] = 0.01 * k; // 100 micron increments;
            fDrift_time[k] = values[k]["t"] * 1000; // from us to ns
         }

         // The /CDC/drift_parameters consists of 2 rows and 10 columns
         // The B-field scale parameters B1 and B2 are in columns 9 and 10
         // They should be the same for both rows, so we get them from the first
         std::vector< std::vector<float> > cdc_drift_parms;
         jcalib->Get("/CDC/drift_parameters", cdc_drift_parms);
         fBscale_par1 = cdc_drift_parms.at(0).at(9);
         fBscale_par2 = cdc_drift_parms.at(0).at(10);
      }
      else {
         int nvalues = CDC_DRIFT_TABLE_LEN;
         std::vector< std::map<std::string, float> > values;
         jcalib->Get("CDC/cdc_drift_table::NoBField", values);
         for (int k=0; k < nvalues; ++k) {
            fDrift_distance[k] = 0.01 * k; // 100 micron increments;
            fDrift_time[k] = values[k]["t"] * 1000; // from us to ns
         }
         fBscale_par1 = 1.;
         fBscale_par2 = 0;
      }
      G4cout << "CDC: ALL parameters loaded from ccdb" << G4endl;

      // Check for "driftclusters" option in control.in

      GlueXUserOptions *opts = GlueXUserOptions::GetInstance();
      if (opts) {
         std::map<int, int> driftclusters_opts;
         if (opts->Find("driftclusters", driftclusters_opts))
            fDrift_clusters = driftclusters_opts[1];
      }
   }
}

GlueXSensitiveDetectorCDC::GlueXSensitiveDetectorCDC(
                     const GlueXSensitiveDetectorCDC &src)
 : G4VSensitiveDetector(src),
   fStrawsMap(src.fStrawsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorCDC &GlueXSensitiveDetectorCDC::operator=(const
                                         GlueXSensitiveDetectorCDC &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fStrawsMap = src.fStrawsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorCDC::~GlueXSensitiveDetectorCDC() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorCDC::Initialize(G4HCofThisEvent* hce)
{
   fStrawsMap = new
                GlueXHitsMapCDCstraw(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
                GlueXHitsMapCDCpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fStrawsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorCDC::ProcessHits(G4Step* step, 
                                              G4TouchableHistory* ROhist)
{
   double dEsum = step->GetTotalEnergyDeposit();
   if (dEsum == 0)
      return false;

   const G4ThreeVector &pin = step->GetPreStepPoint()->GetMomentum();
   const G4ThreeVector &xin = step->GetPreStepPoint()->GetPosition();
   const G4ThreeVector &xout = step->GetPostStepPoint()->GetPosition();
   double tin = step->GetPreStepPoint()->GetGlobalTime();
   double tout = step->GetPostStepPoint()->GetGlobalTime();
   G4ThreeVector x = (xin + xout) / 2;
   G4ThreeVector dx = xout - xin;
   double t = (tin + tout) / 2;

   const G4VTouchable* touch = step->GetPreStepPoint()->GetTouchable();
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

   double drin = xinlocal.perp();
   double drout = xoutlocal.perp();

   // Find x2local = the point of closest approach to the z axis
   G4ThreeVector xlocal = (xinlocal + xoutlocal) / 2;
   G4ThreeVector dxlocal = xoutlocal - xinlocal;
   double alpha = -(xinlocal[0] * dxlocal[0] +
                    xinlocal[1] * dxlocal[1]) / (dxlocal.perp2() + 1e-30);
   G4ThreeVector xpoca = xinlocal + alpha * dxlocal; 
   double alpha01 = (alpha < 0)? 0 : (alpha > 1)? 1 : alpha;
   G4ThreeVector x2local = xinlocal + alpha01 * dxlocal; 

   // Deal with tracks exiting the ends of the straws
   if (fabs(x2local[2]) >= 75.45*cm) {
      int sign = (xoutlocal[2] > 0)? 1 : -1;
      int ring = GetIdent("ring", touch);
      if (ring <= 4 || (ring >= 13 && ring <= 16) || ring >= 25) {
         alpha = (sign * 75.45*cm - xinlocal[2]) / (dxlocal[2] + 1e-30);
         xpoca = xinlocal + alpha * dxlocal; 
         alpha01 = (alpha < 0)? 0 : (alpha > 1)? 1 : alpha;
         x2local = xinlocal + alpha01 * dxlocal;
      }
      else if (fabs(x2local[2]) >= 75.575*cm) {
         alpha = (sign * 75.575*cm - xinlocal[2]) / (dxlocal[2] + 1e-30); 
         xpoca = xinlocal + alpha * dxlocal; 
         alpha01 = (alpha < 0)? 0 : (alpha > 1)? 1 : alpha;
         x2local = xinlocal + alpha01 * dxlocal;
      }
   } 

   // Handle the case when the particle actually passes through
   // the wire volume itself. For these cases, we should set the 
   // location of the hit to be the point on the wire itself. To
   // determine if this is what is happening, we check drout to
   // see if it is very close to the wire and drin to see if it is
   // close to the tube. For the other case, when drin is close to
   // the wire, we assume it is because it is emerging from the
   // wire volume and automatically ignore those hits by returning
   // immediately.

   if (drin < 0.0050*cm) {
      return false; /* entering straw within 50 microns of wire. ignore */
   }
 
   if ((drin > (STRAW_RADIUS - 0.0200*cm) && drout < 0.0050*cm) ||
       (drin < 0.274*cm && drin > 0.234*cm && drout < 0.0050*cm))
   {
      // Either we entered within 200 microns of the straw tube and left
      // within 50 microns of the wire or we entered the stub region near the 
      // donuts at either end of the straw (the inner radius of the feedthrough 
      // region is 0.254 cm) and passed near the wire. Assume the track passed 
      // through the wire volume.

      x = xout;
      t = tout;
      xlocal = xoutlocal;
    
      // For dx, we will just assume it is twice the distance from
      // the straw to wire. Approximate the energy loss in the straw to
      // be twice the energy loss in the first half of the straw.

      dx *= 2;
      dEsum *= 2;
   }

   double dradius = xpoca.perp();
   double dr = dx.mag();
   double dEdx = (dr > 1e-3*cm)? dEsum/dr : 0;

   // Post the hit to the points list in the
   // order of appearance in the event simulation.

   int ring = GetIdent("ring", touch);
   int sector = GetIdent("sector", touch);
   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                           track->GetUserInformation();
   G4int itrack = trackinfo->GetGlueXTrackID();
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   if (trackinfo->GetGlueXHistory() == 0) {
      G4int key = fPointsMap->entries();
      GlueXHitCDCpoint* lastPoint = (*fPointsMap)[key - 1];
      // No more than one truth point per ring in the cdc
      int ring = GetIdent("ring", touch);
      if (lastPoint == 0 || lastPoint->track_ != trackID || 
          lastPoint->ring_ != ring)
      {
         GlueXHitCDCpoint newPoint;
         newPoint.ptype_G3 = g3type;
         newPoint.track_ = trackID;
         newPoint.trackID_ = itrack;
         newPoint.primary_ = (track->GetParentID() == 0);
         newPoint.t_ns = t/ns;
         newPoint.z_cm = x[2]/cm;
         newPoint.r_cm = x.perp()/cm;
         newPoint.phi_rad = x.phi();
         newPoint.dradius_cm = dradius/cm;
         newPoint.px_GeV = pin[0]/GeV;
         newPoint.py_GeV = pin[1]/GeV;
         newPoint.pz_GeV = pin[2]/GeV;
         newPoint.dEdx_GeV_cm = dEdx/(GeV/cm);
         newPoint.sector_ = sector;
         newPoint.ring_ = ring;
         fPointsMap->add(key, newPoint);
      }
   }
   
   // Post the hit to the straw hits map, ordered by straw index

   if (dEsum > 0) {
      int key = GlueXHitCDCstraw::GetKey(ring, sector);
      GlueXHitCDCstraw *straw = (*fStrawsMap)[key];
      if (straw == 0) {
         GlueXHitCDCstraw newstraw(ring, sector);
         fStrawsMap->add(key, newstraw);
         straw = (*fStrawsMap)[key];
      }

      // Add the hit to the hits vector, maintaining track time ordering,
      // re-ordering according to hit times will take place and end of event.

      int merge_hit = 0;
      hit_vector_t::iterator hiter;
      for (hiter = straw->hits.begin(); hiter != straw->hits.end(); ++hiter) {
         if (itrack == hiter->itrack_ && fabs(tin - hiter->t1_ns*ns) < 0.01) {
            merge_hit = 1;
            break;
         }
         else if (hiter->t0_ns*ns > tin) {
            break;
         }
      }
      if (merge_hit) {
         double d_cm = x2local.perp()/cm;
         if (d_cm < hiter->d_cm)
            hiter->d_cm = d_cm;
         hiter->q_fC += dEsum;
         hiter->t1_ns = tout/ns;
         hiter->x1_g = xout;
      }
      else {       // create new hit
         hiter = straw->hits.insert(hiter, GlueXHitCDCstraw::hitinfo_t());
         hiter->track_ = trackID;
         hiter->q_fC = dEsum;
         hiter->d_cm = x2local.perp()/cm;
         hiter->itrack_ = itrack;
         hiter->ptype_G3 = g3type;
         hiter->t0_ns = tin/ns;
         hiter->t1_ns = tout/ns;
         hiter->x0_g = xin;
         hiter->x1_g = xout;
      }
   }
   return true;
}

void GlueXSensitiveDetectorCDC::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int, GlueXHitCDCstraw*> *straws = fStrawsMap->GetMap();
   std::map<int, GlueXHitCDCpoint*> *points = fPointsMap->GetMap();
   if (straws->size() == 0 && points->size() == 0)
      return;
   std::map<int, GlueXHitCDCstraw*>::iterator siter;
   std::map<int, GlueXHitCDCpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << straws->size() << " straws with hits in the CDC: "
             << G4endl;
      for (siter = straws->begin(); siter != straws->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the CDC: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorCDC::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getCentralDCs().size() == 0)
      hitview.addCentralDCs();
   hddm_s::CentralDC &centralDC = hitview.getCentralDC();

   // Collect and output the cdcTruthHits

   for (siter = straws->begin(); siter != straws->end(); ++siter) {

      // Merge multiple segments from a single track into one, and 
      // apply the drift time algorithm to get a single hit time for each.

      hit_vector_t &splits = siter->second->hits;
      hit_vector_t hits;
      while (splits.size() > 0) {
         for (unsigned int ih=1; ih < splits.size(); ++ih) {
            if (fabs(splits[ih].t0_ns - splits[0].t1_ns) < 0.5*ns) {
               splits[0].t1_ns = splits[ih].t1_ns;
               splits[0].q_fC += splits[ih].q_fC;
               if (splits[ih].d_cm < splits[0].d_cm)
                  splits[0].d_cm = splits[ih].d_cm;
               splits[0].x1_g = splits[ih].x1_g;
               splits.erase(splits.begin() + ih);
               --ih;
            }
         }
         if (splits[0].q_fC/keV > THRESH_KEV) {

            // Simulate number of primary ion pairs.
            // The total number of ion pairs depends on the energy deposition
            // and the effective average energy to produce a pair.
    
            // Average number of primary ion pairs
            double t = (splits[0].t0_ns + splits[0].t1_ns)*ns / 2;
            double n_p_mean = (splits[0].q_fC / W_EFF_PER_ION) / 
                              (1 + N_SECOND_PER_PRIMARY);
            // Number of generated primary ion pairs
            int n_p = CLHEP::RandPoisson::shoot(n_p_mean);
            if (fDrift_clusters == 0) {
               add_cluster(hits, splits[0], n_p, t, splits[0].x0_g);
            }
            else {
               // Loop over the number of primary ion pairs,
               // generating a cluster at a random position
               // along the track within the straw.
               G4ThreeVector dx = splits[0].x1_g - splits[0].x0_g;
               for (int n=0; n < n_p; n++) {
                  double u = G4RandFlat::shoot();
                  G4ThreeVector x = splits[0].x0_g + u * dx;
                  add_cluster(hits, splits[0], 1, t, x);
               }
            }
         }
         splits.erase(splits.begin());
      }

      // If doing driftclusters generate a sampled waveform and 
      // analyze it to rebuild the hits list from scratch.

      if (fDrift_clusters) {
         // store waveform data in sampled sequence with 1 ns bins
         int num_samples = int(CDC_TIME_WINDOW/ns);
         double *samples = new double[num_samples];
         for (int i=0; i < num_samples; i++) {
            samples[i] = cdc_wire_signal_mV(double(i), hits);
         }

         // take the earliest hit to identify the track parameters
         double dradius_cm = hits[0].d_cm;
         int ptype_G3 = hits[0].ptype_G3;
         int track = hits[0].track_;
         int itrack = hits[0].itrack_;

         hits.clear();
         double q_mV_ns = 0.; 
         int over_threshold = 0;
         for (int i=0; i < num_samples; ++i) {
            if (samples[i] > THRESH_MV) {
               if (!over_threshold) {
                  hits.push_back(GlueXHitCDCstraw::hitinfo_t());
                  hits.back().track_ = track;
                  hits.back().t_ns = double(i);
                  hits.back().d_cm = dradius_cm;
                  hits.back().itrack_ = itrack;
                  hits.back().ptype_G3 = ptype_G3;
                  over_threshold = 1;
               }
               q_mV_ns += samples[i];
            }
            else if (over_threshold) {
               hits.back().q_fC = q_mV_ns;  // warning -- faking the units 
               over_threshold = 0;   
               q_mV_ns = 0;
            }
         }
         delete [] samples;
      }
      else {

         // New reduced hit list is ordered by hit time
 
         hit_vector_t::iterator hiter;
         for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
            int merge_hits = 0;
            hit_vector_t::iterator iter;
            for (iter = splits.begin(); iter != splits.end(); ++iter) {
               if (fabs(hiter->t_ns - iter->t_ns) < TWO_HIT_TIME_RESOL) {
                  merge_hits = 1;
                  break;
               }
               else if (hiter->t_ns < iter->t_ns) {
                  break;
               }
            }
            if (merge_hits) {
               if (hiter->t_ns < iter->t_ns) {
                  double q_fC = iter->q_fC;
                  *iter = *hiter;   
                  iter->q_fC += q_fC;
               }
               else {
                  iter->q_fC += hiter->q_fC;
               }
            }
            else {
               splits.insert(iter, *hiter);
            }
         }
      }

      if (hits.size() > 0) {
         hddm_s::CdcStrawList straw = centralDC.addCdcStraws(1);
         straw(0).setRing(siter->second->ring_);
         straw(0).setStraw(siter->second->sector_);
         int hitscount = hits.size();
         if (hitscount > MAX_HITS) {
            hitscount = MAX_HITS;
            G4cerr << "GlueXSensitiveDetectorCDC::EndOfEvent warning: "
                   << "straw hit count " << hitscount
                   << " overflows MAX_HITS=" << MAX_HITS
                   << ", " << hits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::CdcStrawTruthHitList thit = straw(0).addCdcStrawTruthHits(1);
            thit(0).setQ(hits[ih].q_fC);
            thit(0).setT(hits[ih].t_ns);
            thit(0).setD(hits[ih].d_cm);
            thit(0).setItrack(hits[ih].itrack_);
            thit(0).setPtype(hits[ih].ptype_G3);
         }
      }
   }

   // Collect and output the cdcTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::CdcTruthPointList point = centralDC.addCdcTruthPoints(1);
      point(0).setDEdx(piter->second->dEdx_GeV_cm);
      point(0).setDradius(piter->second->dradius_cm);
      point(0).setPrimary(piter->second->primary_);
      point(0).setPtype(piter->second->ptype_G3);
      point(0).setPx(piter->second->px_GeV);
      point(0).setPy(piter->second->py_GeV);
      point(0).setPz(piter->second->pz_GeV);
      point(0).setR(piter->second->r_cm);
      point(0).setPhi(piter->second->phi_rad);
      point(0).setZ(piter->second->z_cm);
      point(0).setT(piter->second->t_ns);
      point(0).setTrack(piter->second->track_);
      hddm_s::TrackIDList tid = point(0).addTrackIDs();
      tid(0).setItrack(piter->second->trackID_);
   }
}

double GlueXSensitiveDetectorCDC::asic_response(double t_ns)
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

double GlueXSensitiveDetectorCDC::cdc_wire_signal_mV(double t_ns,
                                                     hit_vector_t &hits)
{
   // Simulation of signal on a wire

   double asic_gain = 0.5; // mV/fC
   double signal_mV = 0;
   hit_vector_t::iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      if (t_ns > hiter->t_ns) {
         double my_time_ns = t_ns - hiter->t_ns;
         signal_mV += asic_gain * hiter->q_fC * asic_response(my_time_ns);
      }
  }
  return signal_mV;
}

void GlueXSensitiveDetectorCDC::add_cluster(hit_vector_t &hits,
                                            GlueXHitCDCstraw::hitinfo_t &hit,
                                            int n_p,
                                            double t, 
                                            G4ThreeVector &x)
{
   // measured charge 
   double q_fC = 0;

   // drift radius 
   double dradius_cm = hit.d_cm;

   // Find the drift time for this cluster. Drift time depends on B:
   // (dependence derived from Garfield calculations)

   G4ThreeVector B = GlueXDetectorConstruction::GetInstance()
                     ->GetMagneticField(x, tesla);
   double BmagT = B.mag();

  // Check for closeness to boundaries of the drift table

   double my_t_ns;
   double my_t_err;
   int i = (int)(dradius_cm * 100);
   if (i >= CDC_DRIFT_TABLE_LEN - 3) {
      // Do a crude linear extrapolation
      my_t_ns = fDrift_time[CDC_DRIFT_TABLE_LEN - 3] + 
                (dradius_cm - fDrift_distance[CDC_DRIFT_TABLE_LEN - 3]) *
                (fDrift_time[CDC_DRIFT_TABLE_LEN - 1] -
                 fDrift_time[CDC_DRIFT_TABLE_LEN - 3]) / 0.02;
   }
   else {
      int index = (i < 1)? 0 : i - 1;

      // Interpolate over the drift table to find 
      // an approximation for the drift time
      polint(&fDrift_distance[index],
             &fDrift_time[index], 4, dradius_cm, &my_t_ns, &my_t_err);
   }
   double tdrift_ns = my_t_ns / (fBscale_par1 + fBscale_par2 * BmagT);

   // Prevent unphysical times (drift electrons arriving 
   // at wire before particle passes the doca to the wire) 
   double v_max = 0.08; // guess for now based on Garfield, near wire 
   double tmin_ns = dradius_cm / v_max;
   if (tdrift_ns < tmin_ns) {
      tdrift_ns = tmin_ns;
   }
   double total_time = t + tdrift_ns*ns;

   // Skip cluster if the time would go beyond readout window
   if (total_time > CDC_TIME_WINDOW) {
      hit.q_fC = 0;
      hit.t_ns = 0;
      return;
   }

   if (fDrift_clusters == 0) {
      // Total number of ion pairs.  On average for each primary ion 
      // pair produced there are n_s secondary ion pairs produced.  The
      // probability distribution is a compound poisson distribution
      // that requires generating two Poisson variables.
      double n_s_mean = n_p * N_SECOND_PER_PRIMARY;
      int n_s = CLHEP::RandPoisson::shoot(n_s_mean);
      int n_t = n_s + n_p;
      q_fC = n_t * GAS_GAIN * ELECTRON_CHARGE/fC;
   }
   else {
      // Distribute the number of secondary ionizations for this primary
      // ionization according to a Poisson distribution with mean n_s_over_p.
      // For simplicity we assume these secondary electrons and the primary
      // electron stay together as a cluster.
      int n_s = CLHEP::RandPoisson::shoot(N_SECOND_PER_PRIMARY);
      // Energy deposition, equivalent to anode charge, in units of fC
      q_fC = GAS_GAIN * ELECTRON_CHARGE/fC * (1 + n_s);
   }
   hit.q_fC = q_fC;
   hit.t_ns = total_time/ns;
   hits.push_back(hit);
}

void GlueXSensitiveDetectorCDC::polint(double *xa, double *ya, int n,
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

int GlueXSensitiveDetectorCDC::GetIdent(std::string div, 
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

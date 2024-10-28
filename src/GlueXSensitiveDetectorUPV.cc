//
// GlueXSensitiveDetectorUPV - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016

#include "GlueXSensitiveDetectorUPV.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "HddmOutput.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <JANA/JApplication.h>

// Cutoff on the number of allowed hits per counter
int GlueXSensitiveDetectorUPV::MAX_HITS = 100;

// Light propagation parameters in tof bars
double GlueXSensitiveDetectorUPV::ATTENUATION_LENGTH = 150.*cm;
double GlueXSensitiveDetectorUPV::C_EFFECTIVE = 19*cm/ns;

// Minimum hit time difference for two hits on the same bar
double GlueXSensitiveDetectorUPV::TWO_HIT_TIME_RESOL = 50*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorUPV::THRESH_MEV = 5.;

int GlueXSensitiveDetectorUPV::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorUPV::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorUPV::GlueXSensitiveDetectorUPV(const G4String& name)
 : G4VSensitiveDetector(name),
   fBarHitsMap(0), fPointsMap(0)
{
   collectionName.insert("UPVBarHitsCollection");
   collectionName.insert("UPVPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the UPV, you must delete all old
   // objects of this class and create new ones.

   G4AutoLock barrier(&fMutex);
   if (instanceCount++ == 0) {
      int runno = HddmOutput::getRunNo();
      extern jana::JApplication *japp;
      if (japp == 0) {
         G4cerr << "Error in GlueXSensitiveDetector constructor - "
                << "jana global DApplication object not set, "
                << "cannot continue." << G4endl;
         exit(-1);
      }
      jana::JCalibration *jcalib = japp->GetJCalibration(runno);
      if (japp == 0) {   // dummy
         jcalib = 0;
         G4cout << "UPV: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorUPV::GlueXSensitiveDetectorUPV(
                     const GlueXSensitiveDetectorUPV &src)
 : G4VSensitiveDetector(src),
   fBarHitsMap(src.fBarHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorUPV &GlueXSensitiveDetectorUPV::operator=(const
                                         GlueXSensitiveDetectorUPV &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fBarHitsMap = src.fBarHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorUPV::~GlueXSensitiveDetectorUPV() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorUPV::Initialize(G4HCofThisEvent* hce)
{
   fBarHitsMap = new
              GlueXHitsMapUPVbar(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapUPVpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fBarHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorUPV::ProcessHits(G4Step* step, 
                                              G4TouchableHistory* ROhist)
{
   double dEsum = step->GetTotalEnergyDeposit();
   if (dEsum == 0)
      return false;

   const G4ThreeVector &pin = step->GetPreStepPoint()->GetMomentum();
   const G4ThreeVector &xin = step->GetPreStepPoint()->GetPosition();
   const G4ThreeVector &xout = step->GetPostStepPoint()->GetPosition();
   double Ein = step->GetPreStepPoint()->GetTotalEnergy();
   double tin = step->GetPreStepPoint()->GetGlobalTime();
   double tout = step->GetPostStepPoint()->GetGlobalTime();
   G4ThreeVector x = (xin + xout) / 2;
   G4ThreeVector dx = xout - xin;
   double t = (tin + tout) / 2;

   const G4VTouchable* touch = step->GetPreStepPoint()->GetTouchable();
   const G4AffineTransform &local_from_global = touch->GetHistory()
                                                     ->GetTopTransform();
   G4ThreeVector xlocal = local_from_global.TransformPoint(x);
  
   // For particles that range out inside the active volume, the
   // "out" time may sometimes be set to something enormously high.
   // This screws up the hit. Check for this case here by looking
   // at tout and making sure it is less than 1 second. If it's
   // not, then just use tin for "t".

   if (tout > 1.0*s)
      t = tin;

   // Post the hit to the points list in the
   // order of appearance in the event simulation.

   int layer = GetIdent("layer", touch);
   int row = GetIdent("row", touch);
   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   if (trackinfo->GetGlueXHistory() == 0 && Ein/MeV > THRESH_MEV) {
      G4int key = fPointsMap->entries();
      GlueXHitUPVpoint* lastPoint = (*fPointsMap)[key - 1];
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          fabs(lastPoint->t_ns - t/ns) > 0.1 ||
          fabs(lastPoint->x_cm - x[0]/cm) > 2. ||
          fabs(lastPoint->y_cm - x[1]/cm) > 2. ||
          fabs(lastPoint->z_cm - x[2]/cm) > 2.)
      {
         int pdgtype = track->GetDynamicParticle()->GetPDGcode();
         int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
         GlueXHitUPVpoint newPoint;
         newPoint.ptype_G3 = g3type;
         newPoint.track_ = trackID;
         newPoint.trackID_ = trackinfo->GetGlueXTrackID();
         newPoint.primary_ = (track->GetParentID() == 0);
         newPoint.t_ns = t/ns;
         newPoint.x_cm = x[0]/cm;
         newPoint.y_cm = x[1]/cm;
         newPoint.z_cm = x[2]/cm;
         newPoint.px_GeV = pin[0]/GeV;
         newPoint.py_GeV = pin[1]/GeV;
         newPoint.pz_GeV = pin[2]/GeV;
         newPoint.E_GeV = Ein/GeV;
         fPointsMap->add(key, newPoint);
      }
   }

   // Post the hit to the hits map, ordered by layer,row,end index

   if (dEsum > 0) {
      int key = GlueXHitUPVbar::GetKey(layer, row);
      GlueXHitUPVbar *counter = (*fBarHitsMap)[key];
      if (counter == 0) {
         GlueXHitUPVbar newcounter(layer, row);
         fBarHitsMap->add(key, newcounter);
         counter = (*fBarHitsMap)[key];
      }

      double dxleft = xlocal[0];
      double dxright = -xlocal[0];

      // Calculate time at the PMT "normalized" to the center, 
      // so a hit in the center will have time "t" at both PMTs.
      double tleft  = t + dxleft / C_EFFECTIVE;
      double tright = t + dxright / C_EFFECTIVE;

      // calculate energy seen by PM for this track step using attenuation factor
      double dEleft  = dEsum * exp(-dxleft / ATTENUATION_LENGTH);
      double dEright = dEsum * exp(-dxright / ATTENUATION_LENGTH);

      // Add the hit to the hits vector, maintaining strict time ordering

      if (dEleft > 0) {
         // add the hit on end=0 (north/left end of the bar)
         int merge_hit = 0;
         std::vector<GlueXHitUPVbar::hitinfo_t>::iterator hiter;
         for (hiter = counter->hits.begin(); hiter != counter->hits.end(); ++hiter) {
            if (hiter->end_ == 0) {
               if (fabs(hiter->t_ns*ns - tleft) < TWO_HIT_TIME_RESOL) {
                  merge_hit = 1;
                  break;
               }
               else if (hiter->t_ns*ns > tleft) {
                  break;
               }
            }
         }
         if (merge_hit) {
            // Use the time from the earlier hit but add the charge
            hiter->E_GeV += dEleft/GeV;
            if (hiter->t_ns*ns > tleft) {
               hiter->t_ns = tleft/ns;
            }
         }
         else {
            // create new hit 
            hiter = counter->hits.insert(hiter, GlueXHitUPVbar::hitinfo_t());
            hiter->end_ = 0;
            hiter->E_GeV = dEleft/GeV;
            hiter->t_ns = tleft/ns;
         }
      }

      if (dEright/MeV > 0) {
         // add the hit on end=1 (south/bottom end of the bar)
         int merge_hit = 0;
         std::vector<GlueXHitUPVbar::hitinfo_t>::iterator hiter;
         for (hiter = counter->hits.begin(); hiter != counter->hits.end(); ++hiter) {
            if (hiter->end_ == 1) {
               if (fabs(hiter->t_ns*ns - tright) < TWO_HIT_TIME_RESOL) {
                  merge_hit = 1;
                  break;
               }
               else if (hiter->t_ns*ns > tright) {
                  break;
               }
            }
         }
         if (merge_hit) {
            // Use the time from the earlier hit but add the charge
            hiter->E_GeV += dEright/GeV;
            if (hiter->t_ns*ns > tright) {
               hiter->t_ns = tright/ns;
            }
         }
         else {
            // create new hit 
            hiter = counter->hits.insert(hiter, GlueXHitUPVbar::hitinfo_t());
            hiter->end_ = 1;
            hiter->E_GeV = dEright/GeV;
            hiter->t_ns = tright/ns;
         }
      }
   }
   return true;
}

void GlueXSensitiveDetectorUPV::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitUPVbar*> *bars = fBarHitsMap->GetMap();
   std::map<int,GlueXHitUPVpoint*> *points = fPointsMap->GetMap();
   if (bars->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitUPVbar*>::iterator siter;
   std::map<int,GlueXHitUPVpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << bars->size() << " bars with hits in the UPV: "
             << G4endl;
      for (siter = bars->begin(); siter != bars->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the UPV: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorUPV::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getUpstreamEMvetos().size() == 0)
      hitview.addUpstreamEMvetos();
   hddm_s::UpstreamEMveto &upv = hitview.getUpstreamEMveto();

   // Collect and output the tofTruthHits
   for (siter = bars->begin(); siter != bars->end(); ++siter) {
      std::vector<GlueXHitUPVbar::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].E_GeV*1e3 <= THRESH_MEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::UpvPaddleList counter = upv.addUpvPaddles(1);
         counter(0).setLayer(siter->second->layer_);
         counter(0).setRow(siter->second->row_);
         // first the end=0 hits
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            if (hits[ih].end_ == 0) {
               hddm_s::UpvTruthHitList thit = counter(0).addUpvTruthHits(1);
               thit(0).setEnd(hits[ih].end_);
               thit(0).setE(hits[ih].E_GeV);
               thit(0).setT(hits[ih].t_ns);
            }
         }
         // followed by the end=1 hits
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            if (hits[ih].end_ == 1) {
               hddm_s::UpvTruthHitList thit = counter(0).addUpvTruthHits(1);
               thit(0).setEnd(hits[ih].end_);
               thit(0).setE(hits[ih].E_GeV);
               thit(0).setT(hits[ih].t_ns);
            }
         }
         int hitscount = counter(0).getUpvTruthHits().size();
         if (hitscount > 2 * MAX_HITS) {
            counter(0).deleteUpvTruthHits(-1, 2 * MAX_HITS);
            G4cerr << "GlueXSensitiveDetectorUPV::ENdOfEvent warning: "
                << "max hit count " << 2 * MAX_HITS << " exceeded, "
                << hitscount - 2 * MAX_HITS << " hits discarded."
                << G4endl;
         }
      }
   }

   // Collect and output the barTruthShowers
   GlueXUserEventInformation* eventinfo = (GlueXUserEventInformation*)info;
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::UpvTruthShowerList point = upv.addUpvTruthShowers(1);
      point(0).setPrimary(piter->second->primary_);
      point(0).setPtype(piter->second->ptype_G3);
      point(0).setPx(piter->second->px_GeV);
      point(0).setPy(piter->second->py_GeV);
      point(0).setPz(piter->second->pz_GeV);
      point(0).setE(piter->second->E_GeV);
      point(0).setX(piter->second->x_cm);
      point(0).setY(piter->second->y_cm);
      point(0).setZ(piter->second->z_cm);
      point(0).setT(piter->second->t_ns);
      point(0).setTrack(piter->second->track_);
      hddm_s::TrackIDList tid = point(0).addTrackIDs();
      tid(0).setItrack(piter->second->trackID_);
      int itrack = piter->second->track_;
      const GlueXUserEventInformation::parent_history_t *parent_history;
      while ((parent_history = eventinfo->GetParentHistory(itrack))) {
         hddm_s::TrackOriginList origin = point(0).addTrackOrigins(1);
	 origin(0).setItrack(parent_history->parent_id);
	 origin(0).setPtype(parent_history->g3type);
	 origin(0).setX(parent_history->x0[0]);
	 origin(0).setY(parent_history->x0[1]);
	 origin(0).setZ(parent_history->x0[2]);
	 origin(0).setT(parent_history->t0);
	 itrack = parent_history->parent_id;
      }
   }
}

int GlueXSensitiveDetectorUPV::GetIdent(std::string div, 
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

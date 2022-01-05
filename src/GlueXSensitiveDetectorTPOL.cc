//
// GlueXSensitiveDetectorTPOL - class implementation
//
// author: richard.t.jones at uconn.edu
// version: december 16, 2016

#include "GlueXSensitiveDetectorTPOL.hh"
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

#define sqr(x) ((x)*(x))

// Cutoff on the total number of allowed hits
int GlueXSensitiveDetectorTPOL::MAX_HITS = 100;

// Minimum hit time difference for two hits on the same wedge
double GlueXSensitiveDetectorTPOL::TWO_HIT_TIME_RESOL = 1000*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorTPOL::THRESH_MEV = 0.050;

int GlueXSensitiveDetectorTPOL::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorTPOL::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorTPOL::GlueXSensitiveDetectorTPOL(const G4String& name)
 : G4VSensitiveDetector(name),
   fHitsMap(0), fPointsMap(0)
{
   collectionName.insert("TPOLWedgeHitsCollection");
   collectionName.insert("TPOLPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the TPOL, you must delete all old
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
         G4cout << "TPOL: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorTPOL::GlueXSensitiveDetectorTPOL(
                     const GlueXSensitiveDetectorTPOL &src)
 : G4VSensitiveDetector(src),
   fHitsMap(src.fHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorTPOL &GlueXSensitiveDetectorTPOL::operator=(const
                                         GlueXSensitiveDetectorTPOL &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fHitsMap = src.fHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorTPOL::~GlueXSensitiveDetectorTPOL() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorTPOL::Initialize(G4HCofThisEvent* hce)
{
   fHitsMap = new
              GlueXHitsMapTPOLwedge(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapTPOLpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorTPOL::ProcessHits(G4Step* step, 
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

   double dr = dx.mag();
   double dEdx = (dr > 1e-3*cm)? dEsum/dr : 0;

   // Post the hit to the points list in the
   // order of appearance in the event simulation.

   int ring = 0; // GetIdent("ring", touch);
   int sector = GetIdent("sector", touch);
   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (trackinfo->GetGlueXHistory() == 0 && itrack > 0) {
      G4int key = fPointsMap->entries();
      GlueXHitTPOLpoint* lastPoint = (*fPointsMap)[key - 1];
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          fabs(lastPoint->t_ns - t/ns) > 0.1 ||
          (fabs(lastPoint->r_cm - x.perp()/cm) > 0.1 &&
           fabs(lastPoint->phi_rad - x.phi()) > 0.1) )
      {
         GlueXHitTPOLpoint newPoint;
         newPoint.ptype_G3 = g3type;
         newPoint.track_ = trackID;
         newPoint.trackID_ = itrack;
         newPoint.primary_ = (track->GetParentID() == 0);
         newPoint.phi_rad = x.phi();
         newPoint.r_cm = x.perp()/cm;
         newPoint.t_ns = t/ns;
         newPoint.px_GeV = pin[0]/GeV;
         newPoint.py_GeV = pin[1]/GeV;
         newPoint.pz_GeV = pin[2]/GeV;
         newPoint.E_GeV = Ein/GeV;
         newPoint.dEdx_GeV_cm = dEdx/(GeV/cm);
         fPointsMap->add(key, newPoint);
      }
   }

   // Post the hit to the hits map, ordered by sector,ring index

   if (dEsum > 0) {
      int key = GlueXHitTPOLwedge::GetKey(sector, ring);
      GlueXHitTPOLwedge *wedge = (*fHitsMap)[key];
      if (wedge == 0) {
         GlueXHitTPOLwedge newwedge(sector, ring);
         fHitsMap->add(key, newwedge);
         wedge = (*fHitsMap)[key];
      }

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitTPOLwedge::hitinfo_t>::iterator hiter;
      for (hiter = wedge->hits.begin(); hiter != wedge->hits.end(); ++hiter) {
         if (fabs(hiter->t_ns*ns - t) < TWO_HIT_TIME_RESOL) {
            merge_hit = 1;
            break;
         }
         else if (hiter->t_ns*ns > t) {
            break;
         }
      }
      if (merge_hit) {
         // Use the time from the earlier hit but add the charge
         hiter->dE_MeV += dEsum/MeV;
         if (hiter->t_ns*ns > t) {
            hiter->t_ns = t/ns;
            hiter->itrack_ = itrack;
            hiter->ptype_G3 = g3type;
            hiter->t0_ns = t/ns;
            hiter->r_cm = x.perp()/cm;
         }
      }
      else if ((int)wedge->hits.size() < MAX_HITS) {
         // create new hit 
         hiter = wedge->hits.insert(hiter, GlueXHitTPOLwedge::hitinfo_t());
         hiter->dE_MeV = dEsum/MeV;
         hiter->t_ns = t/ns;
         hiter->itrack_ = itrack;
         hiter->ptype_G3 = g3type;
         hiter->t0_ns = t/ns;
         hiter->r_cm = x.perp()/cm;
      }
      else {
         G4cerr << "GlueXSensitiveDetectorTPOL::ProcessHits error: "
             << "max hit count " << MAX_HITS << " exceeded, truncating!"
             << G4endl;
      }
   }
   return true;
}

void GlueXSensitiveDetectorTPOL::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitTPOLwedge*> *wedges = fHitsMap->GetMap();
   std::map<int,GlueXHitTPOLpoint*> *points = fPointsMap->GetMap();
   if (wedges->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitTPOLwedge*>::iterator siter;
   std::map<int,GlueXHitTPOLpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << wedges->size() << " wedges with hits in the TPOL: "
             << G4endl;
      for (siter = wedges->begin(); siter != wedges->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the TPOL: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorTPOL::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getTripletPolarimeters().size() == 0)
      hitview.addTripletPolarimeters();
   hddm_s::TripletPolarimeter &polarimeter = hitview.getTripletPolarimeter();

   // Collect and output the tpolTruthHits
   for (siter = wedges->begin(); siter != wedges->end(); ++siter) {
      std::vector<GlueXHitTPOLwedge::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].dE_MeV <= THRESH_MEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::TpolSectorList wedge = polarimeter.addTpolSectors(1);
         wedge(0).setSector(siter->second->sector_);
         wedge(0).setRing(siter->second->ring_);
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            hddm_s::TpolTruthHitList thit = wedge(0).addTpolTruthHits(1);
            thit(0).setDE(hits[ih].dE_MeV*MeV/GeV);
            thit(0).setT(hits[ih].t_ns);
            thit(0).setItrack(hits[ih].itrack_);
            thit(0).setPtype(hits[ih].ptype_G3);
         }
      }
   }

   // Collect and output the tpolTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::TpolTruthPointList point = polarimeter.addTpolTruthPoints(1);
      point(0).setE(piter->second->E_GeV);
      point(0).setDEdx(piter->second->dEdx_GeV_cm);
      point(0).setPrimary(piter->second->primary_);
      point(0).setPtype(piter->second->ptype_G3);
      point(0).setPx(piter->second->px_GeV);
      point(0).setPy(piter->second->py_GeV);
      point(0).setPz(piter->second->pz_GeV);
      point(0).setPhi(piter->second->phi_rad);
      point(0).setR(piter->second->r_cm);
      point(0).setT(piter->second->t_ns);
      point(0).setTrack(piter->second->track_);
      hddm_s::TrackIDList tid = point(0).addTrackIDs();
      tid(0).setItrack(piter->second->trackID_);
   }
}

int GlueXSensitiveDetectorTPOL::GetIdent(std::string div, 
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

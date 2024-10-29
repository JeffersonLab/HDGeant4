//
// GlueXSensitiveDetectorCERE - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 29, 2016

#include "GlueXSensitiveDetectorCERE.hh"
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

// Geant3-style particle index for optical photon
#define OPTICAL_PHOTON 50

// Cutoff on the number of allowed hits per section
int GlueXSensitiveDetectorCERE::MAX_HITS = 100;

// Minimum hit time difference for two hits on the same tube
double GlueXSensitiveDetectorCERE::TWO_HIT_TIME_RESOL = 50*ns;

// Minimum photoelectron count for a hit
double GlueXSensitiveDetectorCERE::THRESH_PE = 2.;

int GlueXSensitiveDetectorCERE::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorCERE::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorCERE::GlueXSensitiveDetectorCERE(const G4String& name)
 : G4VSensitiveDetector(name),
   fTubeHitsMap(0), fPointsMap(0)
{
   collectionName.insert("CERETubeHitsCollection");
   collectionName.insert("CEREPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the CERE, you must delete all old
   // objects of this class and create new ones.

   G4AutoLock tuberier(&fMutex);
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
         G4cout << "CERE: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorCERE::GlueXSensitiveDetectorCERE(
                     const GlueXSensitiveDetectorCERE &src)
 : G4VSensitiveDetector(src),
   fTubeHitsMap(src.fTubeHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock tuberier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorCERE &GlueXSensitiveDetectorCERE::operator=(const
                                         GlueXSensitiveDetectorCERE &src)
{
   G4AutoLock tuberier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fTubeHitsMap = src.fTubeHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorCERE::~GlueXSensitiveDetectorCERE() 
{
   G4AutoLock tuberier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorCERE::Initialize(G4HCofThisEvent* hce)
{
   fTubeHitsMap = new
              GlueXHitsMapCEREtube(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapCEREpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fTubeHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorCERE::ProcessHits(G4Step* step, 
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

   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (touch->GetVolume()->GetName() == "CERW") {
      if (trackinfo->GetGlueXHistory() == 0 && itrack > 0 && xin.dot(pin) > 0) {
         G4int key = fPointsMap->entries();
         GlueXHitCEREpoint* lastPoint = (*fPointsMap)[key - 1];
         if (lastPoint == 0 || lastPoint->track_ != trackID ||
             fabs(lastPoint->t_ns - t/ns) > 0.1 ||
             fabs(lastPoint->x_cm - x[0]/cm) > 2. ||
             fabs(lastPoint->y_cm - x[1]/cm) > 2. ||
             fabs(lastPoint->z_cm - x[2]/cm) > 2.)
         {
            GlueXHitCEREpoint newPoint;
            newPoint.ptype_G3 = g3type;
            newPoint.track_ = trackID;
            newPoint.trackID_ = itrack;
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
      // This sensitive detector is unique in that different volumes are used
      // to generate truth points and hits, which means that a single entry
      // to ProcessEvents can make one or the other, but not both.
      return true;
   }

   // Post the hit to the hits map, ordered by plane,tube,end index

   if (dEsum > 0) {
      int sector = GetIdent("sector", touch);
      int key = GlueXHitCEREtube::GetKey(sector);
      GlueXHitCEREtube *counter = (*fTubeHitsMap)[key];
      if (counter == 0) {
         GlueXHitCEREtube newcounter(sector);
         fTubeHitsMap->add(key, newcounter);
         counter = (*fTubeHitsMap)[key];
      }

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitCEREtube::hitinfo_t>::iterator hiter;
      for (hiter = counter->hits.begin(); hiter != counter->hits.end(); ++hiter) {
         if (fabs(hiter->t_ns*ns - t) < TWO_HIT_TIME_RESOL) {
            merge_hit = 1;
            break;
         }
         else if (hiter->t_ns*ns > t) {
            break;
         }
         if (merge_hit) {
            // Use the time from the earlier hit but add the charge
            hiter->pe_ += 1;
            if (hiter->t_ns*ns > t) {
               hiter->t_ns = t/ns;
            }
         }
         else {
            // create new hit 
            hiter = counter->hits.insert(hiter, GlueXHitCEREtube::hitinfo_t());
            hiter->pe_ = 1;
            hiter->t_ns = t/ns;
         }
      }
   }
   return true;
}

void GlueXSensitiveDetectorCERE::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitCEREtube*> *tubes = fTubeHitsMap->GetMap();
   std::map<int,GlueXHitCEREpoint*> *points = fPointsMap->GetMap();
   if (tubes->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitCEREtube*>::iterator siter;
   std::map<int,GlueXHitCEREpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << tubes->size() << " tubes with hits in the Cerenkov: "
             << G4endl;
      for (siter = tubes->begin(); siter != tubes->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the Cerenkov: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorCERE::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getCerenkovs().size() == 0)
      hitview.addCerenkovs();
   hddm_s::Cerenkov &cerenkov = hitview.getCerenkov();

   // Collect and output the cereTruthHits
   for (siter = tubes->begin(); siter != tubes->end(); ++siter) {
      std::vector<GlueXHitCEREtube::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].pe_ <= THRESH_PE) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::CereSectionList tube = cerenkov.addCereSections(1);
         tube(0).setSector(siter->second->sector_);
         int hitscount = hits.size();
         if (hitscount > MAX_HITS) {
            hitscount = MAX_HITS;
            G4cerr << "GlueXSensitiveDetectorCERE::EndOfEvent warning: "
                   << "max tube hit count " << MAX_HITS << " exceeded, "
                   << hits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::CereTruthHitList thit = tube(0).addCereTruthHits(1);
            thit(0).setPe(hits[ih].pe_);
            thit(0).setT(hits[ih].t_ns);
         }
      }
   }

   // Collect and output the cereTruthPoints
   GlueXUserEventInformation* eventinfo = (GlueXUserEventInformation*)info;
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::CereTruthPointList point = cerenkov.addCereTruthPoints(1);
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
	 origin(0).setX(parent_history->x0[0]/cm);
	 origin(0).setY(parent_history->x0[1]/cm);
	 origin(0).setZ(parent_history->x0[2]/cm);
	 origin(0).setT(parent_history->t0/ns);
	 itrack = parent_history->parent_id;
      }
   }
}

int GlueXSensitiveDetectorCERE::GetIdent(std::string div, 
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

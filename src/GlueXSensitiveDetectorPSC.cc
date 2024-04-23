//
// GlueXSensitiveDetectorPSC - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016

#include "GlueXSensitiveDetectorPSC.hh"
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
int GlueXSensitiveDetectorPSC::MAX_HITS = 100;

// Minimum hit time difference for two hits on the same paddle
double GlueXSensitiveDetectorPSC::TWO_HIT_TIME_RESOL = 25*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorPSC::THRESH_MEV = 0.010;

// Geometric parameters
int GlueXSensitiveDetectorPSC::NUM_MODULES_PER_ARM = 8;

int GlueXSensitiveDetectorPSC::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorPSC::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorPSC::GlueXSensitiveDetectorPSC(const G4String& name)
 : G4VSensitiveDetector(name),
   fCounterHitsMap(0), fPointsMap(0)
{
   collectionName.insert("PSCPaddleHitsCollection");
   collectionName.insert("PSCPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the PSC, you must delete all old
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
      if (japp == 0) {    // dummy
         jcalib = 0;
         G4cout << "PSC: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorPSC::GlueXSensitiveDetectorPSC(
                     const GlueXSensitiveDetectorPSC &src)
 : G4VSensitiveDetector(src),
   fCounterHitsMap(src.fCounterHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorPSC &GlueXSensitiveDetectorPSC::operator=(const
                                         GlueXSensitiveDetectorPSC &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fCounterHitsMap = src.fCounterHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorPSC::~GlueXSensitiveDetectorPSC() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorPSC::Initialize(G4HCofThisEvent* hce)
{
   fCounterHitsMap = new
              GlueXHitsMapPSCpaddle(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapPSCpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fCounterHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorPSC::ProcessHits(G4Step* step, 
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

   int module = GetIdent("module", touch);
   int arm = (module - 1) / NUM_MODULES_PER_ARM;
   module = (module - 1) % NUM_MODULES_PER_ARM + 1;
   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (trackinfo->GetGlueXHistory() == 0 && itrack > 0) {
      G4int key = fPointsMap->entries();
      GlueXHitPSCpoint* lastPoint = (*fPointsMap)[key - 1];
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          fabs(lastPoint->t_ns - t/ns) > 0.2 ||
          fabs(lastPoint->x_cm - x[0]/cm) > 5.0 ||
          fabs(lastPoint->y_cm - x[1]/cm) > 5.0 ||
          fabs(lastPoint->z_cm - x[2]/cm) > 5.0)
      {
         GlueXHitPSCpoint newPoint;
         newPoint.arm_ = arm;
         newPoint.module_ = module;
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
         newPoint.dEdx_GeV_cm = dEdx/(GeV/cm);
         fPointsMap->add(key, newPoint);
      }
   }

   // Post the hit to the hits map, ordered by arm,module index

   if (dEsum > 0) {
      int key = GlueXHitPSCpaddle::GetKey(arm, module);
      GlueXHitPSCpaddle *paddle = (*fCounterHitsMap)[key];
      if (paddle == 0) {
         GlueXHitPSCpaddle newpaddle(arm, module);
         fCounterHitsMap->add(key, newpaddle);
         paddle = (*fCounterHitsMap)[key];
      }

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitPSCpaddle::hitinfo_t>::iterator hiter;
      for (hiter = paddle->hits.begin(); hiter != paddle->hits.end(); ++hiter) {
         if (fabs(hiter->t_ns*ns - t) < TWO_HIT_TIME_RESOL) {
            merge_hit = 1;
            break;
         }
         else if (hiter->t_ns*ns > t) {
            break;
         }
      }
      if (merge_hit) {
         // Add the charge, do energy-weighted time averaging --disabled for bad behavior
         //hiter->t_ns = (hiter->t_ns * hiter->dE_GeV + t/ns * dEsum/GeV) /
         //              (hiter->dE_GeV + dEsum/GeV);
         //hiter->dE_GeV += dEsum/GeV;
 
         // Use the time from the earlier hit but add the charge
         hiter->dE_GeV += dEsum/GeV;
         if (hiter->t_ns*ns > t) {
             hiter->t_ns = t/ns;
         }
      }
      else {
         // create new hit 
         hiter = paddle->hits.insert(hiter, GlueXHitPSCpaddle::hitinfo_t());
         hiter->dE_GeV = dEsum/GeV;
         hiter->t_ns = t/ns;
         hiter->itrack_ = itrack;
         hiter->ptype_G3 = g3type;
      }
   }
   return true;
}

void GlueXSensitiveDetectorPSC::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitPSCpaddle*> *paddles = fCounterHitsMap->GetMap();
   std::map<int,GlueXHitPSCpoint*> *points = fPointsMap->GetMap();
   if (paddles->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitPSCpaddle*>::iterator siter;
   std::map<int,GlueXHitPSCpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << paddles->size() << " paddles with hits in the PSC: "
             << G4endl;
      for (siter = paddles->begin(); siter != paddles->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the PSC: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorPSC::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getPairSpectrometerCoarses().size() == 0)
      hitview.addPairSpectrometerCoarses();
   hddm_s::PairSpectrometerCoarse &psc = hitview.getPairSpectrometerCoarse();

   // Collect and output the PscTruthHits
   for (siter = paddles->begin(); siter != paddles->end(); ++siter) {
      std::vector<GlueXHitPSCpaddle::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].dE_GeV*1e3 <= THRESH_MEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::PscPaddleList paddle = psc.addPscPaddles(1);
         paddle(0).setArm(siter->second->arm_);
         paddle(0).setModule(siter->second->module_);
         int hitscount = hits.size();
         if (hitscount > MAX_HITS) {
            hitscount = MAX_HITS;
            G4cerr << "GlueXSensitiveDetectorPSC::EndOfEvent warning: "
                   << "max hit count " << MAX_HITS << " exceeded, "
                   << hits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::PscTruthHitList thit = paddle(0).addPscTruthHits(1);
            thit(0).setDE(hits[ih].dE_GeV);
            thit(0).setT(hits[ih].t_ns);
            thit(0).setItrack(hits[ih].itrack_);
            thit(0).setPtype(hits[ih].ptype_G3);
         }
      }
   }

   // Collect and output the pscTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::PscTruthPointList point = psc.addPscTruthPoints(1);
      point(0).setE(piter->second->E_GeV);
      point(0).setDEdx(piter->second->dEdx_GeV_cm);
      point(0).setPrimary(piter->second->primary_);
      point(0).setPtype(piter->second->ptype_G3);
      point(0).setPx(piter->second->px_GeV);
      point(0).setPy(piter->second->py_GeV);
      point(0).setPz(piter->second->pz_GeV);
      point(0).setX(piter->second->x_cm);
      point(0).setY(piter->second->y_cm);
      point(0).setZ(piter->second->z_cm);
      point(0).setT(piter->second->t_ns);
      point(0).setArm(piter->second->arm_);
      point(0).setModule(piter->second->module_);
      point(0).setTrack(piter->second->track_);
      hddm_s::TrackIDList tid = point(0).addTrackIDs();
      tid(0).setItrack(piter->second->trackID_);
   }
}

int GlueXSensitiveDetectorPSC::GetIdent(std::string div, 
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

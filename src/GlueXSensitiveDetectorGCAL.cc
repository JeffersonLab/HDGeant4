//
// GlueXSensitiveDetectorGCAL - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016

#include "GlueXSensitiveDetectorGCAL.hh"
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
#include <JANA/Calibrations/JCalibrationManager.h>


// Cutoff on the number of allowed hits per cell
int GlueXSensitiveDetectorGCAL::MAX_HITS = 100;

// Light propagation parameters in start counter
double GlueXSensitiveDetectorGCAL::ATTENUATION_LENGTH = 1e6*cm;
double GlueXSensitiveDetectorGCAL::C_EFFECTIVE = 15*cm/ns;
double GlueXSensitiveDetectorGCAL::LENGTH_OF_BLOCK = 45.*cm;

// Minimum hit time difference for two hits on the same module
double GlueXSensitiveDetectorGCAL::TWO_HIT_TIME_RESOL = 75*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorGCAL::THRESH_MEV = 30.;

int GlueXSensitiveDetectorGCAL::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorGCAL::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorGCAL::GlueXSensitiveDetectorGCAL(const G4String& name)
 : G4VSensitiveDetector(name),
   fBlockHitsMap(0), fPointsMap(0)
{
   collectionName.insert("GCALBlockHitsCollection");
   collectionName.insert("GCALPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the GCAL, you must delete all old
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
      /* JCalibration *jcalib = */ japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
      if (japp == 0) {   // dummy
         //jcalib = 0;
         G4cout << "GCAL: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorGCAL::GlueXSensitiveDetectorGCAL(
                     const GlueXSensitiveDetectorGCAL &src)
 : G4VSensitiveDetector(src),
   fBlockHitsMap(src.fBlockHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorGCAL &GlueXSensitiveDetectorGCAL::operator=(const
                                         GlueXSensitiveDetectorGCAL &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fBlockHitsMap = src.fBlockHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorGCAL::~GlueXSensitiveDetectorGCAL() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorGCAL::Initialize(G4HCofThisEvent* hce)
{
   fBlockHitsMap = new
              GlueXHitsMapGCALblock(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapGCALpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fBlockHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorGCAL::ProcessHits(G4Step* step, 
                                               G4TouchableHistory* /* ROhist */)
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

   int module = GetIdent("module", touch);
   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
      GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                             track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (track->GetCurrentStepNumber() == 1)
      trackinfo->SetGlueXHistory(3);
   if (trackinfo->GetGlueXHistory() == 0 && itrack > 0 &&
       xin.dot(pin) > 0 && Ein/MeV > THRESH_MEV)
   {
      G4int key = fPointsMap->entries();
      GlueXHitGCALpoint* lastPoint = (*fPointsMap)[key - 1];
      // Limit gcal truthPoints to one per track
      if (lastPoint == 0 || lastPoint->track_ != trackID) {
         int pdgtype = track->GetDynamicParticle()->GetPDGcode();
         int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
         GlueXHitGCALpoint newPoint;
         newPoint.ptype_G3 = g3type;
         newPoint.track_ = trackID;
         newPoint.trackID_ = itrack;
         newPoint.primary_ = (track->GetParentID() == 0);
         newPoint.t_ns = t/ns;
         newPoint.z_cm = x[2]/cm;
         newPoint.r_cm = x.perp()/cm;
         newPoint.phi_rad = x.phi();
         newPoint.px_GeV = pin[0]/GeV;
         newPoint.py_GeV = pin[1]/GeV;
         newPoint.pz_GeV = pin[2]/GeV;
         newPoint.E_GeV = Ein/GeV;
         fPointsMap->add(key, newPoint);
      }
   }

   // Post the hit to the hits map, ordered by module index

   if (dEsum > 0) {
      int key = GlueXHitGCALblock::GetKey(module);
      GlueXHitGCALblock *block = (*fBlockHitsMap)[key];
      if (block == 0) {
         GlueXHitGCALblock newblock(module);
         fBlockHitsMap->add(key, newblock);
         block = (*fBlockHitsMap)[key];
      }
      double dist = 0.5 * LENGTH_OF_BLOCK - xlocal[2];
      double dEcorr = dEsum * exp(-dist / ATTENUATION_LENGTH);
      double tcorr = t + dist / C_EFFECTIVE;

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitGCALblock::hitinfo_t>::iterator hiter;
      for (hiter = block->hits.begin(); hiter != block->hits.end(); ++hiter) {
         if (fabs(hiter->t_ns*ns - tcorr) < TWO_HIT_TIME_RESOL) {
            merge_hit = 1;
            break;
         }
         else if (hiter->t_ns*ns > tcorr) {
            break;
         }
      }
      if (merge_hit) {
         // Use the time from the earlier hit but add the charge
         hiter->E_GeV += dEcorr/GeV;
         if (hiter->t_ns*ns > tcorr) {
            hiter->t_ns = tcorr/ns;
         }
      }
      else {
         // create new hit 
         hiter = block->hits.insert(hiter, GlueXHitGCALblock::hitinfo_t());
         hiter->E_GeV = dEcorr/GeV;
         hiter->t_ns = tcorr/ns;
      }
   }
   return true;
}

void GlueXSensitiveDetectorGCAL::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitGCALblock*> *blocks = fBlockHitsMap->GetMap();
   std::map<int,GlueXHitGCALpoint*> *points = fPointsMap->GetMap();
   if (blocks->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitGCALblock*>::iterator siter;
   std::map<int,GlueXHitGCALpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << blocks->size() << " blocks with hits in the GCAL: "
             << G4endl;
      for (siter = blocks->begin(); siter != blocks->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the GCAL: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorGCAL::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getGapEMcals().size() == 0)
      hitview.addGapEMcals();
   hddm_s::GapEMcal &gcal = hitview.getGapEMcal();

   // Collect and output the gcalTruthHits
   for (siter = blocks->begin(); siter != blocks->end(); ++siter) {
      std::vector<GlueXHitGCALblock::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].E_GeV*1e3 <= THRESH_MEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::GcalCellList block = gcal.addGcalCells(1);
         block(0).setModule(siter->second->module_);
         int hitscount = hits.size();
         if (hitscount > MAX_HITS) {
            hitscount = MAX_HITS;
            G4cerr << "GlueXSensitiveDetectorGCAL::EndOfEvent warning: "
                   << "max cell hit count " << MAX_HITS << " exceeded, "
                   << hits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::GcalTruthHitList thit = block(0).addGcalTruthHits(1);
            thit(0).setE(hits[ih].E_GeV);
            thit(0).setT(hits[ih].t_ns);
         }
      }
   }

   // Collect and output the gcalTruthShowers
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::GcalTruthShowerList point = gcal.addGcalTruthShowers(1);
      point(0).setE(piter->second->E_GeV);
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

int GlueXSensitiveDetectorGCAL::GetIdent(std::string div, 
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

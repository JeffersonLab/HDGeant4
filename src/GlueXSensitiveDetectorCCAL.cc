//
// GlueXSensitiveDetectorCCAL - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016

#include "GlueXSensitiveDetectorCCAL.hh"
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

// Cutoff on the number of allowed hits per block
int GlueXSensitiveDetectorCCAL::MAX_HITS = 100;

// Geometry constants for the CCal
//int GlueXSensitiveDetectorCCAL::CENTRAL_ROW = 8;
//int GlueXSensitiveDetectorCCAL::CENTRAL_COLUMN = 8;
double GlueXSensitiveDetectorCCAL::WIDTH_OF_BLOCK = 2.0*cm;
double GlueXSensitiveDetectorCCAL::LENGTH_OF_BLOCK = 20.0*cm;

// Light propagation parameters in Compton calorimeter
double GlueXSensitiveDetectorCCAL::ATTENUATION_LENGTH = 200.*cm;
double GlueXSensitiveDetectorCCAL::C_EFFECTIVE = 13.*cm/ns;

// Minimum hit time difference for two hits on the same block
double GlueXSensitiveDetectorCCAL::TWO_HIT_TIME_RESOL = 75.*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorCCAL::THRESH_MEV = 5.;

int GlueXSensitiveDetectorCCAL::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorCCAL::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorCCAL::GlueXSensitiveDetectorCCAL(const G4String& name)
 : G4VSensitiveDetector(name),
   fBlocksMap(0), fPointsMap(0)
{
   collectionName.insert("CCALBlockHitsCollection");
   collectionName.insert("CCALPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the CCAL, you must delete all old
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
         G4cout << "CCAL: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorCCAL::GlueXSensitiveDetectorCCAL(
                     const GlueXSensitiveDetectorCCAL &src)
 : G4VSensitiveDetector(src),
   fBlocksMap(src.fBlocksMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorCCAL &GlueXSensitiveDetectorCCAL::operator=(const
                                         GlueXSensitiveDetectorCCAL &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fBlocksMap = src.fBlocksMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorCCAL::~GlueXSensitiveDetectorCCAL() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorCCAL::Initialize(G4HCofThisEvent* hce)
{
   fBlocksMap = new
              GlueXHitsMapCCALblock(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapCCALpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fBlocksMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorCCAL::ProcessHits(G4Step* step, 
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
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (track->GetCurrentStepNumber() == 1)
      trackinfo->SetGlueXHistory(4);
   if (trackinfo->GetGlueXHistory() == 0 &&
       xin.dot(pin) > 0 && Ein/MeV > THRESH_MEV)
   {
      G4int key = fPointsMap->entries();
      GlueXHitCCALpoint* lastPoint = (*fPointsMap)[key - 1];
      // Limit ccal truthPoints to one per track
      if (lastPoint == 0 || lastPoint->track_ != trackID) {
         int pdgtype = track->GetDynamicParticle()->GetPDGcode();
         int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
         GlueXHitCCALpoint newPoint;
         newPoint.ptype_G3 = g3type;
         newPoint.track_ = trackID;
         newPoint.trackID_ = itrack;
         newPoint.primary_ = (track->GetParentID() == 0);
         newPoint.t_ns = t/ns;
         newPoint.x_cm = xin[0]/cm;
         newPoint.y_cm = xin[1]/cm;
         newPoint.z_cm = xin[2]/cm;
         newPoint.px_GeV = pin[0]/GeV;
         newPoint.py_GeV = pin[1]/GeV;
         newPoint.pz_GeV = pin[2]/GeV;
         newPoint.E_GeV = Ein/GeV;
         fPointsMap->add(key, newPoint);
      }
   }

   // Post the hit to the hits map, ordered by sector index

   if (dEsum > 0) {
      int column = GetIdent("column", touch);
      int row = GetIdent("row", touch);
      int key = GlueXHitCCALblock::GetKey(column, row);
      GlueXHitCCALblock *block = (*fBlocksMap)[key];
      if (block == 0) {
         GlueXHitCCALblock newblock(column, row);
         fBlocksMap->add(key, newblock);
         block = (*fBlocksMap)[key];
      }
      double dist = 0.5 * LENGTH_OF_BLOCK - xlocal[2];
      double dEcorr = dEsum * exp(-dist / ATTENUATION_LENGTH);
      double tcorr = t + dist / C_EFFECTIVE;

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitCCALblock::hitinfo_t>::iterator hiter;
      for (hiter = block->hits.begin(); hiter != block->hits.end(); ++hiter) {
         if (fabs(hiter->t_ns*ns - tcorr) < TWO_HIT_TIME_RESOL) {
            merge_hit = 1;
            break;
         }
         else if (hiter->t_ns*ns > t) {
            break;
         }
      }
      if (merge_hit) {
         // Use the time from the earlier hit but add the energy deposition
         hiter->E_GeV += dEcorr/GeV;
         if (hiter->t_ns*ns > tcorr) {
            hiter->t_ns = tcorr/ns;
         }
      }
      else {
         // create new hit 
          hiter = block->hits.insert(hiter, GlueXHitCCALblock::hitinfo_t());
         hiter->E_GeV = dEcorr/GeV;
         hiter->t_ns = tcorr/ns;
      }
   }
   return true;
}

void GlueXSensitiveDetectorCCAL::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitCCALblock*> *blocks = fBlocksMap->GetMap();
   std::map<int,GlueXHitCCALpoint*> *points = fPointsMap->GetMap();
   if (blocks->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitCCALblock*>::iterator biter;
   std::map<int,GlueXHitCCALpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << blocks->size() << " blocks with hits in the CCAL: "
             << G4endl;
      for (biter = blocks->begin(); biter != blocks->end(); ++biter)
         biter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth showers in the CCAL: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorCCAL::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getComptonEMcals().size() == 0)
      hitview.addComptonEMcals();
   hddm_s::ComptonEMcal &comptonEMcal = hitview.getComptonEMcal();

   // Collect and output the ccalTruthHits
   for (biter = blocks->begin(); biter != blocks->end(); ++biter) {
      std::vector<GlueXHitCCALblock::hitinfo_t> &hits = biter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].E_GeV <= THRESH_MEV/1e3) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }

      if (hits.size() > 0) {
         hddm_s::CcalBlockList block = comptonEMcal.addCcalBlocks(1);
         block(0).setColumn(biter->second->column_);
         block(0).setRow(biter->second->row_);
         int hitscount = hits.size();
         if (hitscount > MAX_HITS) {
            hitscount = MAX_HITS;
            G4cerr << "GlueXSensitiveDetectorCCAL::EndOfEvent warning: "
                   << "max hit count " << MAX_HITS << " exceeded, "
                   << hits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::CcalTruthHitList thit = block(0).addCcalTruthHits(1);
            thit(0).setE(hits[ih].E_GeV);
            thit(0).setT(hits[ih].t_ns);
         }
      }
   }

   // Collect and output the fcalTruthShowers
   GlueXUserEventInformation* eventinfo = (GlueXUserEventInformation*)info;
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::CcalTruthShowerList point = comptonEMcal.addCcalTruthShowers(1);
      point(0).setE(piter->second->E_GeV);
      point(0).setPrimary(piter->second->primary_);
      point(0).setPtype(piter->second->ptype_G3);
      point(0).setPx(piter->second->px_GeV);
      point(0).setPy(piter->second->py_GeV);
      point(0).setPz(piter->second->pz_GeV);
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

int GlueXSensitiveDetectorCCAL::GetIdent(std::string div, 
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

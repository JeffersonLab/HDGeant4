//
// GlueXSensitiveDetectorPS - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016

#include "GlueXSensitiveDetectorPS.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"

#include <CLHEP/Random/RandPoisson.h>
#include <Randomize.hh>

#include "G4THitsMap.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <JANA/JApplication.h>

#include <stdio.h>
#include <malloc.h>
#include <math.h>

// Cutoff on the total number of allowed hits
int GlueXSensitiveDetectorPS::MAX_HITS = 100;

// Minimum hit time difference for two hits on the same tile
double GlueXSensitiveDetectorPS::TWO_HIT_TIME_RESOL = 25*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorPS::THRESH_MEV = 0.010;

// Geometric parameters
int GlueXSensitiveDetectorPS::NUM_COLUMNS_PER_ARM = 8;

int GlueXSensitiveDetectorPS::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorPS::fMutex = G4MUTEX_INITIALIZER;

std::map<G4LogicalVolume*, int> GlueXSensitiveDetectorPS::fVolumeTable;

GlueXSensitiveDetectorPS::GlueXSensitiveDetectorPS(const G4String& name)
 : G4VSensitiveDetector(name),
   fTileHitsMap(0), fPointsMap(0)
{
   collectionName.insert("PSTileHitsCollection");
   collectionName.insert("PSPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the PS, you must delete all old
   // objects of this class and create new ones.

   G4AutoLock barrier(&fMutex);
   if (instanceCount++ == 0) {
      extern int run_number;
      extern jana::JApplication *japp;
      if (japp == 0) {
         G4cerr << "Error in GlueXSensitiveDetector constructor - "
                << "jana global DApplication object not set, "
                << "cannot continue." << G4endl;
         exit(-1);
      }
      jana::JCalibration *jcalib = japp->GetJCalibration(run_number);
      if (japp == 0) {    // dummy
         jcalib = 0;
         G4cout << "PS: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorPS::GlueXSensitiveDetectorPS(
                     const GlueXSensitiveDetectorPS &src)
 : G4VSensitiveDetector(src),
   fTileHitsMap(src.fTileHitsMap), fPointsMap(src.fPointsMap)
{
   ++instanceCount;
}

GlueXSensitiveDetectorPS &GlueXSensitiveDetectorPS::operator=(const
                                         GlueXSensitiveDetectorPS &src)
{
   *(G4VSensitiveDetector*)this = src;
   fTileHitsMap = src.fTileHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorPS::~GlueXSensitiveDetectorPS() 
{
   --instanceCount;
}

void GlueXSensitiveDetectorPS::Initialize(G4HCofThisEvent* hce)
{
   fTileHitsMap = new 
              GlueXHitsMapPStile(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapPSpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fTileHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorPS::ProcessHits(G4Step* step, 
                                              G4TouchableHistory* unused)
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

   int column = GetIdent("column", touch);
   int arm = column / NUM_COLUMNS_PER_ARM;
   column = column % NUM_COLUMNS_PER_ARM;
   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
      GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                             track->GetUserInformation();
   if (trackinfo->GetGlueXHistory() == 0) {
      G4int key = fPointsMap->entries();
      GlueXHitPSpoint* lastPoint = (*fPointsMap)[key - 1];
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          fabs(lastPoint->t_ns - t/ns) > 0.1 ||
          fabs(lastPoint->z_cm - x[2]/cm) > 0.1)
      {
         GlueXHitPSpoint* newPoint = new GlueXHitPSpoint();
         fPointsMap->add(key, newPoint);
         int pdgtype = track->GetDynamicParticle()->GetPDGcode();
         int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
         newPoint->arm_ = arm;
         newPoint->column_ = column;
         newPoint->ptype_G3 = g3type;
         newPoint->track_ = trackID;
         newPoint->trackID_ = trackinfo->GetGlueXTrackID();
         newPoint->primary_ = (track->GetParentID() == 0);
         newPoint->t_ns = t/ns;
         newPoint->x_cm = x[0]/cm;
         newPoint->y_cm = x[1]/cm;
         newPoint->z_cm = x[2]/cm;
         newPoint->phi_rad = x.phi();
         newPoint->px_GeV = pin[0]/GeV;
         newPoint->py_GeV = pin[1]/GeV;
         newPoint->pz_GeV = pin[2]/GeV;
         newPoint->E_GeV = Ein/GeV;
         newPoint->dEdx_GeV_cm = dEdx/(GeV/cm);
      }
   }

   // Post the hit to the hits map, ordered by arm,column index

   if (dEsum > 0) {
      int key = GlueXHitPStile::GetKey(arm, column);
      GlueXHitPStile *tile = (*fTileHitsMap)[key];
      if (tile == 0) {
         tile = new GlueXHitPStile(arm, column);
         fTileHitsMap->add(key, tile);
      }

      // Add the hit to the hits vector, maintaining strict time ordering

      std::vector<GlueXHitPStile::hitinfo_t>::iterator hiter;
      for (hiter = tile->hits.begin(); hiter != tile->hits.end(); ++hiter) {
         if (fabs(hiter->t_ns*ns - t) < TWO_HIT_TIME_RESOL) {
            break;
         }
         else if (hiter->t_ns*ns > t) {
            hiter = tile->hits.insert(hiter, GlueXHitPStile::hitinfo_t());
            hiter->t_ns = 1e99;
            break;
         }
      }

      GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                             track->GetUserInformation();
      if (hiter != tile->hits.end()) {             // merge with former hit
         // Use the time from the earlier hit but add the charge
         hiter->dE_GeV += dEsum/GeV;
         if (hiter->t_ns*ns > t) {
            hiter->t_ns = t/ns;
            int pdgtype = track->GetDynamicParticle()->GetPDGcode();
            int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
            hiter->itrack_ = trackinfo->GetGlueXTrackID();
            hiter->ptype_G3 = g3type;
         }
      }
      else if ((int)tile->hits.size() < MAX_HITS)	{   // create new hit 
         GlueXHitPStile::hitinfo_t newhit;
         newhit.dE_GeV = dEsum/GeV;
         newhit.t_ns = t/ns;
         int pdgtype = track->GetDynamicParticle()->GetPDGcode();
         int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
         newhit.itrack_ = trackinfo->GetGlueXTrackID();
         newhit.ptype_G3 = g3type;
         tile->hits.push_back(newhit);
      }
      else {
         G4cerr << "GlueXSensitiveDetectorPS::ProcessHits error: "
             << "max hit count " << MAX_HITS << " exceeded, truncating!"
             << G4endl;
      }
   }
   return true;
}

void GlueXSensitiveDetectorPS::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitPStile*> *tiles = fTileHitsMap->GetMap();
   std::map<int,GlueXHitPSpoint*> *points = fPointsMap->GetMap();
   if (tiles->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitPStile*>::iterator siter;
   std::map<int,GlueXHitPSpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << tiles->size() << " tiles with hits in the PS: "
             << G4endl;
      for (siter = tiles->begin(); siter != tiles->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the PS: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorPS::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getPairSpectrometerFines().size() == 0)
      hitview.addPairSpectrometerFines();
   hddm_s::PairSpectrometerFine &psTile = hitview.getPairSpectrometerFine();

   // Collect and output the stcTruthHits
   for (siter = tiles->begin(); siter != tiles->end(); ++siter) {
      std::vector<GlueXHitPStile::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].dE_GeV*1e3 < THRESH_MEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::PsTileList tile = psTile.addPsTiles(1);
         tile(0).setArm(siter->second->arm_);
         tile(0).setColumn(siter->second->column_);
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            hddm_s::PsTruthHitList thit = tile(0).addPsTruthHits(1);
            thit(0).setDE(hits[ih].dE_GeV);
            thit(0).setT(hits[ih].t_ns);
            thit(0).setItrack(hits[ih].itrack_);
            thit(0).setPtype(hits[ih].ptype_G3);
         }
      }
   }

   // Collect and output the tileTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::PsTruthPointList point = psTile.addPsTruthPoints(1);
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
      point(0).setPhi(piter->second->phi_rad);
      point(0).setT(piter->second->t_ns);
      point(0).setArm(piter->second->arm_);
      point(0).setColumn(piter->second->column_);
      point(0).setTrack(piter->second->track_);
      hddm_s::TrackIDList tid = point(0).addTrackIDs();
      tid(0).setItrack(piter->second->trackID_);
   }
}

int GlueXSensitiveDetectorPS::GetIdent(std::string div, 
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
         if (dynamic_cast<G4PVPlacement*>(pvol))
            return iter->second[pvol->GetCopyNo() - 1];
         else
            return iter->second[pvol->GetCopyNo()];
      }
   }
   return -1;
}


//
// GlueXSensitiveDetectorFCAL - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 25, 2016

#include "GlueXSensitiveDetectorFCAL.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <JANA/JApplication.h>

// Scale factors to match reconstructed to generated
double GlueXSensitiveDetectorFCAL::SHOWER_ENERGY_SCALE_FACTOR =  1.0;
double GlueXSensitiveDetectorFCAL::MIP_ENERGY_SCALE_FACTOR = 1.35;

// Cutoff on the total number of allowed hits
int GlueXSensitiveDetectorFCAL::MAX_HITS = 100;

// Geometry constants for the FCal
double GlueXSensitiveDetectorFCAL::WIDTH_OF_BLOCK = 4.0*cm;
double GlueXSensitiveDetectorFCAL::LENGTH_OF_BLOCK = 45.0*cm;
double GlueXSensitiveDetectorFCAL::ACTIVE_RADIUS = 120.0*cm;
int GlueXSensitiveDetectorFCAL::CENTRAL_COLUMN = 29;
int GlueXSensitiveDetectorFCAL::CENTRAL_ROW = 29;

// Light propagation parameters in forward calorimeter
double GlueXSensitiveDetectorFCAL::ATTENUATION_LENGTH = 100.*cm;
double GlueXSensitiveDetectorFCAL::C_EFFECTIVE = 15.*cm/ns;

// Minimum hit time difference for two hits on the same block
double GlueXSensitiveDetectorFCAL::TWO_HIT_TIME_RESOL = 75.*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorFCAL::THRESH_MEV = 5.;

int GlueXSensitiveDetectorFCAL::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorFCAL::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorFCAL::GlueXSensitiveDetectorFCAL(const G4String& name)
 : G4VSensitiveDetector(name),
   fBlocksMap(0), fPointsMap(0)
{
   collectionName.insert("FCALBlockHitsCollection");
   collectionName.insert("FCALPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the FCAL, you must delete all old
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
      std::map<string, float> fcal_parms;
      jcalib->Get("FCAL/fcal_parms", fcal_parms);
      ATTENUATION_LENGTH = fcal_parms.at("FCAL_ATTEN_LENGTH")*cm;
      C_EFFECTIVE = fcal_parms.at("FCAL_C_EFFECTIVE")*cm/ns;
      WIDTH_OF_BLOCK = fcal_parms.at("FCAL_WIDTH_OF_BLOCK")*cm;
      LENGTH_OF_BLOCK = fcal_parms.at("FCAL_LENGTH_OF_BLOCK")*cm;
      ACTIVE_RADIUS = fcal_parms.at("FCAL_ACTIVE_RADIUS")*cm;
      CENTRAL_ROW = fcal_parms.at("FCAL_CENTRAL_ROW");
      CENTRAL_COLUMN= fcal_parms.at("FCAL_CENTRAL_COLUMN");
      TWO_HIT_TIME_RESOL = fcal_parms.at("FCAL_TWO_HIT_RESOL")*ns;
      MAX_HITS = fcal_parms.at("FCAL_MAX_HITS");
      THRESH_MEV = fcal_parms.at("FCAL_THRESH_MEV");

      G4cout << "FCAL: ALL parameters loaded from ccdb" << G4endl;
   }
}

GlueXSensitiveDetectorFCAL::GlueXSensitiveDetectorFCAL(
                     const GlueXSensitiveDetectorFCAL &src)
 : G4VSensitiveDetector(src),
   fBlocksMap(src.fBlocksMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorFCAL &GlueXSensitiveDetectorFCAL::operator=(const
                                         GlueXSensitiveDetectorFCAL &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fBlocksMap = src.fBlocksMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorFCAL::~GlueXSensitiveDetectorFCAL() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorFCAL::Initialize(G4HCofThisEvent* hce)
{
   fBlocksMap = new
              GlueXHitsMapFCALblock(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapFCALpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fBlocksMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorFCAL::ProcessHits(G4Step* step, 
                                               G4TouchableHistory* ROhist)
{
   double dEsum = step->GetTotalEnergyDeposit();
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
   if (trackinfo->GetGlueXHistory() == 0 &&
       xin.dot(pin) > 0 && Ein/MeV > THRESH_MEV)
   {
      int pdgtype = track->GetDynamicParticle()->GetPDGcode();
      int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
      GlueXHitFCALpoint newPoint;
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
      G4int key = fPointsMap->entries();
      fPointsMap->add(key, newPoint);
      trackinfo->SetGlueXHistory(2);
   }

   // Post the hit to the hits map, ordered by sector index

   if (dEsum > 0) {
      int column = GetIdent("column", touch);
      int row = GetIdent("row", touch);
      int key = GlueXHitFCALblock::GetKey(column, row);
      GlueXHitFCALblock *block = (*fBlocksMap)[key];
      if (block == 0) {
         GlueXHitFCALblock newblock(column, row);
         fBlocksMap->add(key, newblock);
         block = (*fBlocksMap)[key];
      }

      // Handle hits in the lead glass
 
      if (touch->GetVolume()->GetName() == "LGBL") {
         double dist = 0.5 * LENGTH_OF_BLOCK - xlocal[2];
         double dEcorr = dEsum * exp(-dist / ATTENUATION_LENGTH);
         double tcorr = t + dist / C_EFFECTIVE;

         // Apply effective response corrections, depending on particle type
         int pmass = track->GetDynamicParticle()->GetMass();
         if (pmass < 1 * MeV) { // must be one of e+,e-,gamma
            dEcorr *= SHOWER_ENERGY_SCALE_FACTOR;
         }
         else {
            double gamma = Ein / pmass; // nothing massless here
            dEcorr *= (gamma > 1.25)? MIP_ENERGY_SCALE_FACTOR : 0;
         }

         if (dEcorr == 0)
            return true;

         // Add the hit to the hits vector, maintaining strict time ordering

         int merge_hit = 0;
         std::vector<GlueXHitFCALblock::hitinfo_t>::iterator hiter;
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
            // Merge the time with the existing hit, add the energy deposition
            hiter->t_ns = hiter->t_ns * hiter->E_GeV + dEcorr/GeV * tcorr/ns;
            hiter->E_GeV += dEcorr/GeV;
            hiter->t_ns /= hiter->E_GeV;
         }
         else if ((int)block->hits.size() < MAX_HITS) {
            // create new hit 
            hiter = block->hits.insert(hiter, GlueXHitFCALblock::hitinfo_t());
            hiter->E_GeV = dEcorr/GeV;
            hiter->t_ns = tcorr/ns;
            hiter->dE_lightguide_GeV = 0;
            hiter->t_lightguide_ns = 0;
         }
         else {
            G4cerr << "GlueXSensitiveDetectorFCAL::ProcessHits error: "
                << "max hit count " << MAX_HITS << " exceeded, truncating!"
                << G4endl;
         }
      }

      // If not in the lead glass, it should be in the light guide

      else if (touch->GetVolume()->GetName() == "LGLG") {

         // Add the hit to the hits vector, maintaining strict time ordering

         int merge_hit = 0;
         std::vector<GlueXHitFCALblock::hitinfo_t>::iterator hiter;
         for (hiter = block->hits.begin(); hiter != block->hits.end(); ++hiter) {
            if (fabs(hiter->t_ns*ns - t) < TWO_HIT_TIME_RESOL) {
               merge_hit = 1;
               break;
            }
            else if (hiter->t_ns*ns > t) {
               break;
            }
         }
         if (merge_hit) {
            hiter->t_lightguide_ns = 
                   (hiter->t_lightguide_ns * hiter->dE_lightguide_GeV +
                    t/ns * dEsum/GeV) / (hiter->dE_lightguide_GeV + dEsum/GeV);
            hiter->dE_lightguide_GeV += dEsum/GeV;
         }
         else if ((int)block->hits.size() < MAX_HITS) {
            // create new hit 
            hiter = block->hits.insert(hiter, GlueXHitFCALblock::hitinfo_t());
            hiter->dE_lightguide_GeV = dEsum/GeV;
            hiter->t_lightguide_ns = t/ns;
            hiter->E_GeV = 0;
            hiter->t_ns = t/ns;
         }
         else {
            G4cerr << "GlueXSensitiveDetectorFCAL::ProcessHits error: "
                << "max hit count " << MAX_HITS << " exceeded, truncating!"
                << G4endl;
         }
      }
   }
   return true;
}

void GlueXSensitiveDetectorFCAL::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitFCALblock*> *blocks = fBlocksMap->GetMap();
   std::map<int,GlueXHitFCALpoint*> *points = fPointsMap->GetMap();
   if (blocks->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitFCALblock*>::iterator biter;
   std::map<int,GlueXHitFCALpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << blocks->size() << " blocks with hits in the FCAL: "
             << G4endl;
      for (biter = blocks->begin(); biter != blocks->end(); ++biter)
         biter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth showers in the FCAL: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorFCAL::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getForwardEMcals().size() == 0)
      hitview.addForwardEMcals();
   hddm_s::ForwardEMcal &forwardEMcal = hitview.getForwardEMcal();

   // Collect and output the fcalTruthHits
   for (biter = blocks->begin(); biter != blocks->end(); ++biter) {
      std::vector<GlueXHitFCALblock::hitinfo_t> &hits = biter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].E_GeV <= THRESH_MEV/1e3) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }

      int row = biter->second->row_;
      int column = biter->second->column_;
      double y0 = (row - CENTRAL_ROW) * WIDTH_OF_BLOCK;
      double x0 = (column - CENTRAL_COLUMN) * WIDTH_OF_BLOCK;
      double dist = sqrt(x0*x0 + y0*y0);

      if (dist < ACTIVE_RADIUS && hits.size() > 0) {
         hddm_s::FcalBlockList block = forwardEMcal.addFcalBlocks(1);
         block(0).setColumn(biter->second->column_);
         block(0).setRow(biter->second->row_);
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            hddm_s::FcalTruthHitList thit = block(0).addFcalTruthHits(1);
            thit(0).setE(hits[ih].E_GeV);
            thit(0).setT(hits[ih].t_ns);
            if (hits[ih].dE_lightguide_GeV > 0) {
               hddm_s::FcalTruthLightGuideList lghit = 
                       thit(0).addFcalTruthLightGuides(1);
               lghit(0).setDE(hits[ih].dE_lightguide_GeV);
               lghit(0).setT(hits[ih].t_lightguide_ns);
            }
         }
      }
   }

   // Collect and output the fcalTruthShowers
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::FcalTruthShowerList point = forwardEMcal.addFcalTruthShowers(1);
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
   }
}

int GlueXSensitiveDetectorFCAL::GetIdent(std::string div, 
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

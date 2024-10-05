//
// GlueXSensitiveDetectorSTC - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016

#include "GlueXSensitiveDetectorSTC.hh"
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
int GlueXSensitiveDetectorSTC::MAX_HITS = 100;

// Light propagation parameters in start counter
double GlueXSensitiveDetectorSTC::ATTENUATION_LENGTH = 150.*cm;
double GlueXSensitiveDetectorSTC::C_EFFECTIVE = 15*cm/ns;
double GlueXSensitiveDetectorSTC::LIGHT_GUIDE = 0.*cm;
double GlueXSensitiveDetectorSTC::ANGLE_COR = 1.054;

// Geometric parameters of STC scintillators
double GlueXSensitiveDetectorSTC::BENT_REGION = 39.465*cm;
double GlueXSensitiveDetectorSTC::STRAIGHT_LENGTH = 39.465*cm;
double GlueXSensitiveDetectorSTC::BEND_LENGTH = 3.5924*cm;
double GlueXSensitiveDetectorSTC::NOSE_LENGTH = 15.5366*cm;

// Light loss and attenuation constants
double GlueXSensitiveDetectorSTC::STRAIGHT_ATTENUATION_A[NCHANNELS];
double GlueXSensitiveDetectorSTC::STRAIGHT_ATTENUATION_B[NCHANNELS];
double GlueXSensitiveDetectorSTC::STRAIGHT_ATTENUATION_C[NCHANNELS];
double GlueXSensitiveDetectorSTC::BENDNOSE_ATTENUATION_A[NCHANNELS];
double GlueXSensitiveDetectorSTC::BENDNOSE_ATTENUATION_B[NCHANNELS];
double GlueXSensitiveDetectorSTC::BENDNOSE_ATTENUATION_C[NCHANNELS];
double GlueXSensitiveDetectorSTC::STRAIGHT_PROPAGATION_A[NCHANNELS];
double GlueXSensitiveDetectorSTC::STRAIGHT_PROPAGATION_B[NCHANNELS];
double GlueXSensitiveDetectorSTC::BEND_PROPAGATION_A[NCHANNELS];
double GlueXSensitiveDetectorSTC::BEND_PROPAGATION_B[NCHANNELS];
double GlueXSensitiveDetectorSTC::NOSE_PROPAGATION_A[NCHANNELS];
double GlueXSensitiveDetectorSTC::NOSE_PROPAGATION_B[NCHANNELS];

// Minimum hit time difference for two hits on the same paddle
double GlueXSensitiveDetectorSTC::TWO_HIT_TIME_RESOL = 25*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorSTC::THRESH_MEV = 0.150;

int GlueXSensitiveDetectorSTC::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorSTC::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorSTC::GlueXSensitiveDetectorSTC(const G4String& name)
 : G4VSensitiveDetector(name),
   fHitsMap(0), fPointsMap(0)
{
   collectionName.insert("STCPaddleHitsCollection");
   collectionName.insert("STCPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the STC, you must delete all old
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
      std::map<string, float> stc_parms;
      jcalib->Get("START_COUNTER/start_parms", stc_parms);
      ATTENUATION_LENGTH = stc_parms.at("START_ATTEN_LENGTH")*cm;
      C_EFFECTIVE = stc_parms.at("START_C_EFFECTIVE")*cm/ns;
      TWO_HIT_TIME_RESOL = stc_parms.at("START_TWO_HIT_RESOL")*ns;
      MAX_HITS = stc_parms.at("START_MAX_HITS");
      THRESH_MEV = stc_parms.at("START_THRESH_MEV");
      LIGHT_GUIDE = stc_parms.at("START_LIGHT_GUIDE")*cm;
      ANGLE_COR = stc_parms.at("START_ANGLE_COR");
      BENT_REGION = stc_parms.at("START_BENT_REGION")*cm;

      std::vector< std::map<std::string, float> > values;
      jcalib->Get("START_COUNTER/attenuation_factor", values);
      for (unsigned int k=0; k < values.size(); ++k) {
         STRAIGHT_ATTENUATION_A[k] = values[k].at("SC_STRAIGHT_ATTENUATION_A");
         STRAIGHT_ATTENUATION_B[k] = values[k].at("SC_STRAIGHT_ATTENUATION_B");
         STRAIGHT_ATTENUATION_C[k] = values[k].at("SC_STRAIGHT_ATTENUATION_C");
         BENDNOSE_ATTENUATION_A[k] = values[k].at("SC_BENDNOSE_ATTENUATION_A");
         BENDNOSE_ATTENUATION_B[k] = values[k].at("SC_BENDNOSE_ATTENUATION_B");
         BENDNOSE_ATTENUATION_C[k] = values[k].at("SC_BENDNOSE_ATTENUATION_C");
         // B factors are in units of 1/cm, A and C factors are no-care
         STRAIGHT_ATTENUATION_B[k] /= cm;
         BENDNOSE_ATTENUATION_B[k] /= cm;
      }
      jcalib->Get("START_COUNTER/propagation_time_corr", values);
      for (unsigned int k=0; k < values.size(); ++k) {
         STRAIGHT_PROPAGATION_A[k] = values[k].at("a");
         STRAIGHT_PROPAGATION_B[k] = values[k].at("b");
         BEND_PROPAGATION_A[k] = values[k].at("c");
         BEND_PROPAGATION_B[k] = values[k].at("d");
         NOSE_PROPAGATION_A[k] = values[k].at("e");
         NOSE_PROPAGATION_B[k] = values[k].at("f");
         // A factors are in units of ns, B factors are ns/cm
         STRAIGHT_PROPAGATION_A[k] *= ns;
         STRAIGHT_PROPAGATION_B[k] *= ns/cm;
         BEND_PROPAGATION_A[k] *= ns;
         BEND_PROPAGATION_B[k] *= ns/cm;
         NOSE_PROPAGATION_A[k] *= ns;
         NOSE_PROPAGATION_B[k] *= ns/cm;
      }

      G4cout << "STC: ALL parameters loaded from ccdb" << G4endl;
   }
}

GlueXSensitiveDetectorSTC::GlueXSensitiveDetectorSTC(
                     const GlueXSensitiveDetectorSTC &src)
 : G4VSensitiveDetector(src),
   fHitsMap(src.fHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorSTC &GlueXSensitiveDetectorSTC::operator=(const
                                         GlueXSensitiveDetectorSTC &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fHitsMap = src.fHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorSTC::~GlueXSensitiveDetectorSTC() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorSTC::Initialize(G4HCofThisEvent* hce)
{
   fHitsMap = new
              GlueXHitsMapSTCpaddle(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapSTCpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorSTC::ProcessHits(G4Step* step, 
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

   int sector = GetIdent("sector", touch);
   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (trackinfo->GetGlueXHistory() == 0) {
      G4int key = fPointsMap->entries();
      GlueXHitSTCpoint* lastPoint = (*fPointsMap)[key - 1];
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          fabs(lastPoint->t_ns - t/ns) > 0.1 ||
          fabs(lastPoint->z_cm - x[2]/cm) > 0.1)
      {
         GlueXHitSTCpoint newPoint;
         newPoint.ptype_G3 = g3type;
         newPoint.track_ = trackID;
         newPoint.trackID_ = itrack;
         newPoint.primary_ = (track->GetParentID() == 0);
         newPoint.sector_ = sector;
         newPoint.t_ns = t/ns;
         newPoint.z_cm = x[2]/cm;
         newPoint.r_cm = x.perp()/cm;
         newPoint.phi_rad = x.phi();
         newPoint.px_GeV = pin[0]/GeV;
         newPoint.py_GeV = pin[1]/GeV;
         newPoint.pz_GeV = pin[2]/GeV;
         newPoint.E_GeV = Ein/GeV;
         newPoint.dEdx_GeV_cm = dEdx/(GeV/cm);
         fPointsMap->add(key, newPoint);
      }
   }

   // Post the hit to the hits map, ordered by sector index

   if (dEsum > 0) {
      int key = GlueXHitSTCpaddle::GetKey(sector);
      GlueXHitSTCpaddle *paddle = (*fHitsMap)[key];
      if (paddle == 0) {
         GlueXHitSTCpaddle newpaddle(sector);
         fHitsMap->add(key, newpaddle);
         paddle = (*fHitsMap)[key];
      }

      double dbent = 0.0;
      double dpath = 0.0;
      if (xlocal[2] >= BENT_REGION){
      	 dbent = (xlocal[2] - BENT_REGION) * ANGLE_COR;
      	 dpath = BENT_REGION + dbent;
      }
      else {
      	dpath  = xlocal[2];
      }

      int sindex = sector - 1;
      double dEcorr = 9.9E+9;
      double tcorr  = 9.9E+9;
      if (xlocal[2] <= STRAIGHT_LENGTH) {
	     dEcorr = dEsum * exp(dpath * STRAIGHT_ATTENUATION_B[sindex]);
	     tcorr  = t + STRAIGHT_PROPAGATION_A[sindex] +
                      STRAIGHT_PROPAGATION_B[sindex] * dpath;

	  }
      else if (xlocal[2] <= STRAIGHT_LENGTH + BEND_LENGTH) {
         double bna_A = BENDNOSE_ATTENUATION_A[sindex];
         double bna_B = BENDNOSE_ATTENUATION_B[sindex];
         double bna_C = BENDNOSE_ATTENUATION_C[sindex];
         double sta_A = STRAIGHT_ATTENUATION_A[sindex];
	     dEcorr = dEsum * ((bna_A * exp(dpath * bna_B) + bna_C) / sta_A);
	     tcorr  = t + BEND_PROPAGATION_A[sindex] +
                      BEND_PROPAGATION_B[sindex] * dpath;
	  }
      else if (xlocal[2] <= STRAIGHT_LENGTH + BEND_LENGTH + NOSE_LENGTH) {
         double bna_A = BENDNOSE_ATTENUATION_A[sindex];
         double bna_B = BENDNOSE_ATTENUATION_B[sindex];
         double bna_C = BENDNOSE_ATTENUATION_C[sindex];
         double sta_A = STRAIGHT_ATTENUATION_A[sindex];
	     dEcorr = dEsum * ((bna_A * exp(dpath * bna_B) + bna_C) / sta_A);
	     tcorr  = t + NOSE_PROPAGATION_A[sindex] +
                      NOSE_PROPAGATION_B[sindex] * dpath;
	  }
      else {
         return false;
      }

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitSTCpaddle::hitinfo_t>::iterator hiter;
      for (hiter = paddle->hits.begin(); hiter != paddle->hits.end(); ++hiter) {
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
         hiter->dE_MeV += dEcorr/MeV;
         if (hiter->t_ns*ns > tcorr) {
            hiter->t_ns = tcorr/ns;
            hiter->itrack_ = itrack;
            hiter->ptype_G3 = g3type;
            hiter->t0_ns = t/ns;
            hiter->z_cm = x[2]/cm;
         }
      }
      else {
         // create new hit 
         hiter = paddle->hits.insert(hiter, GlueXHitSTCpaddle::hitinfo_t());
         hiter->dE_MeV = dEcorr/MeV;
         hiter->t_ns = tcorr/ns;
         hiter->itrack_ = itrack;
         hiter->ptype_G3 = g3type;
         hiter->t0_ns = t/ns;
         hiter->z_cm = x[2]/cm;
      }
   }
   return true;
}

void GlueXSensitiveDetectorSTC::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitSTCpaddle*> *paddles = fHitsMap->GetMap();
   std::map<int,GlueXHitSTCpoint*> *points = fPointsMap->GetMap();
   if (paddles->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitSTCpaddle*>::iterator siter;
   std::map<int,GlueXHitSTCpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << paddles->size() << " paddles with hits in the STC: "
             << G4endl;
      for (siter = paddles->begin(); siter != paddles->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the STC: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorSTC::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getStartCntrs().size() == 0)
      hitview.addStartCntrs();
   hddm_s::StartCntr &startCntr = hitview.getStartCntr();

   // Collect and output the stcTruthHits
   for (siter = paddles->begin(); siter != paddles->end(); ++siter) {
      std::vector<GlueXHitSTCpaddle::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].dE_MeV <= THRESH_MEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::StcPaddleList paddle = startCntr.addStcPaddles(1);
         paddle(0).setSector(siter->second->sector_);
         int hitscount = hits.size();
         if (hitscount > MAX_HITS) {
            hitscount = MAX_HITS;
            G4cerr << "GlueXSensitiveDetectorSTC::EndOfEvent warning: "
                   << "max hit count " << MAX_HITS << " exceeded, "
                   << hits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::StcTruthHitList thit = paddle(0).addStcTruthHits(1);
            thit(0).setDE(hits[ih].dE_MeV*MeV/GeV);
            thit(0).setT(hits[ih].t_ns);
            thit(0).setItrack(hits[ih].itrack_);
            thit(0).setPtype(hits[ih].ptype_G3);
         }
      }
   }

   // Collect and output the stcTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::StcTruthPointList point = startCntr.addStcTruthPoints(1);
      point(0).setE(piter->second->E_GeV);
      point(0).setDEdx(piter->second->dEdx_GeV_cm);
      point(0).setPrimary(piter->second->primary_);
      point(0).setPtype(piter->second->ptype_G3);
      point(0).setPx(piter->second->px_GeV);
      point(0).setPy(piter->second->py_GeV);
      point(0).setPz(piter->second->pz_GeV);
      point(0).setR(piter->second->r_cm);
      point(0).setPhi(piter->second->phi_rad);
      point(0).setZ(piter->second->z_cm);
      point(0).setT(piter->second->t_ns);
      point(0).setSector(piter->second->sector_);
      point(0).setTrack(piter->second->track_);
      hddm_s::TrackIDList tid = point(0).addTrackIDs();
      tid(0).setItrack(piter->second->trackID_);
   }
}

int GlueXSensitiveDetectorSTC::GetIdent(std::string div, 
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

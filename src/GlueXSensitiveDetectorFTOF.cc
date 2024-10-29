//
// GlueXSensitiveDetectorFTOF - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016

#include "GlueXSensitiveDetectorFTOF.hh"
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

// Cutoff on the number of allowed hits per bar
int GlueXSensitiveDetectorFTOF::MAX_HITS_PER_BAR = 25;

// Cutoff on the maximum time of flight
double GlueXSensitiveDetectorFTOF::MAX_TOF = 1000*ns;

// Light propagation parameters in tof bars
double GlueXSensitiveDetectorFTOF::ATTENUATION_LENGTH = 150.*cm;
double GlueXSensitiveDetectorFTOF::C_EFFECTIVE = 15*cm/ns;

// Geometric parameters of FTOF scintillators
double GlueXSensitiveDetectorFTOF::FULL_BAR_LENGTH = 252.*cm;

// Minimum hit time difference for two hits on the same bar
double GlueXSensitiveDetectorFTOF::TWO_HIT_TIME_RESOL = 25*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorFTOF::THRESH_MEV = 0.;

int GlueXSensitiveDetectorFTOF::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorFTOF::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorFTOF::GlueXSensitiveDetectorFTOF(const G4String& name)
 : G4VSensitiveDetector(name),
   fBarHitsMap(0), fPointsMap(0)
{
   collectionName.insert("FTOFBarHitsCollection");
   collectionName.insert("FTOFPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the FTOF, you must delete all old
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
      std::map<string, float> tof_parms;
      jcalib->Get("TOF/tof_parms", tof_parms);
      ATTENUATION_LENGTH = tof_parms.at("TOF_ATTEN_LENGTH")*cm;
      C_EFFECTIVE = tof_parms.at("TOF_C_EFFECTIVE")*cm/ns;
      FULL_BAR_LENGTH = tof_parms.at("TOF_PADDLE_LENGTH")*cm;
      TWO_HIT_TIME_RESOL = tof_parms.at("TOF_TWO_HIT_RESOL")*ns;
      THRESH_MEV = tof_parms.at("TOF_THRESH_MEV");
      MAX_HITS_PER_BAR = tof_parms.at("TOF_MAX_PAD_HITS");

      G4cout << "FTOF: ALL parameters loaded from ccdb" << G4endl;
   }
}

GlueXSensitiveDetectorFTOF::GlueXSensitiveDetectorFTOF(
                     const GlueXSensitiveDetectorFTOF &src)
 : G4VSensitiveDetector(src),
   fBarHitsMap(src.fBarHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorFTOF &GlueXSensitiveDetectorFTOF::operator=(const
                                         GlueXSensitiveDetectorFTOF &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fBarHitsMap = src.fBarHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorFTOF::~GlueXSensitiveDetectorFTOF() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorFTOF::Initialize(G4HCofThisEvent* hce)
{
   fBarHitsMap = new
              GlueXHitsMapFTOFbar(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapFTOFpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fBarHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorFTOF::ProcessHits(G4Step* step, 
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

   int plane = GetIdent("plane", touch);
   int column = GetIdent("column", touch);
   int barNo = GetIdent("row", touch);
   int barIndex = (column < 2)? barNo : GetIdent("paired_row", touch);
   if (barIndex < 1) {
      G4cerr << "GlueXSensitiveDetectorFTOF::ProcessHits error - "
             << "hdds geometry for FTOF is missing paired_row identifier, "
             << "cannot continue!" << G4endl;
      exit(1);
   }
   G4Track *track = step->GetTrack();
   G4int trackID = track->GetTrackID();
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (trackinfo->GetGlueXHistory() == 0 && plane == 0) {
      G4int key = fPointsMap->entries();
      GlueXHitFTOFpoint* lastPoint = (*fPointsMap)[key - 1];
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          fabs(lastPoint->t_ns - t/ns) > 0.1 ||
          fabs(lastPoint->x_cm - x[0]/cm) > 2. ||
          fabs(lastPoint->y_cm - x[1]/cm) > 2. ||
          fabs(lastPoint->z_cm - x[2]/cm) > 2.)
      {
         GlueXHitFTOFpoint newPoint;
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

   // Post the hit to the hits map, ordered by plane,bar,end index

   if (dEsum > 0) {
      int key = GlueXHitFTOFbar::GetKey(plane, barIndex);
      GlueXHitFTOFbar *counter = (*fBarHitsMap)[key];
      if (counter == 0) {
         GlueXHitFTOFbar newcounter(plane, barIndex);
         fBarHitsMap->add(key, newcounter);
         counter = (*fBarHitsMap)[key];
      }

      double dist = x[1]; // do not use local coordinate for x and y
      if (plane == 1)
         dist = x[0];
      double dxnorth = FULL_BAR_LENGTH / 2 - dist;
      double dxsouth = FULL_BAR_LENGTH / 2 + dist;

      // Calculate time at the PMT "normalized" to the center, 
      // so a hit in the center will have time "t" at both PMTs.
      double tnorth = t + dxnorth / C_EFFECTIVE;
      double tsouth = t + dxsouth / C_EFFECTIVE;
 
      // calculate energy seen by PM for this track step using attenuation factor
      double dEnorth = dEsum * exp(-dxnorth / ATTENUATION_LENGTH);
      double dEsouth = dEsum * exp(-dxsouth / ATTENUATION_LENGTH);

      // column=0 is a full bar, column=1,2 are half paddles
      if (plane == 0) {
         if (column == 1) {
            tnorth = 0;
            dEnorth = 0;
         }
         else if (column == 2) {    
            tsouth = 0;
            dEsouth =0;
         }
      }
      else {
         if (column == 2) {
            tnorth = 0;
            dEnorth = 0;
         }
         else if (column == 1) {
            tsouth = 0;
            dEsouth = 0;
         }
      }

      // Add the hit to the hits vector, maintaining strict time ordering

      if (dEnorth/MeV > 0 && tnorth < MAX_TOF) {
         // add the hit on end=0 (north/top end of the bar)
         int merge_hit = 0;
         std::vector<GlueXHitFTOFbar::hitinfo_t>::iterator hiter;
         for (hiter = counter->hits.begin(); 
              hiter != counter->hits.end(); ++hiter)
         {
            if (hiter->end_ == 0) {
               if (fabs(hiter->t_ns*ns - tnorth) < TWO_HIT_TIME_RESOL) {
                  merge_hit = 1;
                  break;
               }
               else if (hiter->t_ns*ns > tnorth) {
                  break;
               }
            }
         }

         if (merge_hit) {
            // sum the charge, do energy weighting of the time --disabled for bad behavior, rtj--
            //hiter->t_ns = (hiter->dE_GeV * hiter->t_ns +
            //               dEnorth/GeV * tnorth/ns) /
            //(hiter->dE_GeV += dEnorth/GeV);
 
            // Use the time from the earlier hit but add the charge
            hiter->dE_GeV += dEnorth/GeV;
            if (hiter->t_ns*ns > tnorth) {
               hiter->t_ns = tnorth/ns;
            }

            std::vector<GlueXHitFTOFbar::hitextra_t>::reverse_iterator xiter;
            xiter = hiter->extra.rbegin();
            if (trackID != xiter->track_ || fabs(tin/ns - xiter->t_ns) > 0.1) {
               GlueXHitFTOFbar::hitextra_t extra;
               extra.track_ = trackID;
               extra.itrack_ = itrack;
               extra.ptype_G3 = g3type;
               extra.px_GeV = pin[0]/GeV;
               extra.py_GeV = pin[1]/GeV;
               extra.pz_GeV = pin[2]/GeV;
               extra.E_GeV = Ein/GeV;
               extra.x_cm = x[0]/cm;
               extra.y_cm = x[1]/cm;
               extra.z_cm = x[2]/cm;
               extra.t_ns = tout/ns;
               extra.dist_cm = dist/cm;
               hiter->extra.push_back(extra);
            }
         }
         else {
            // create new hit 
            hiter = counter->hits.insert(hiter, GlueXHitFTOFbar::hitinfo_t());
            hiter->end_ = 0;
            hiter->t_ns = tnorth/ns;
            hiter->dE_GeV = dEnorth/GeV;
            GlueXHitFTOFbar::hitextra_t extra;
            extra.track_ = trackID;
            extra.itrack_ = itrack;
            extra.ptype_G3 = g3type;
            extra.px_GeV = pin[0]/GeV;
            extra.py_GeV = pin[1]/GeV;
            extra.pz_GeV = pin[2]/GeV;
            extra.E_GeV = Ein/GeV;
            extra.x_cm = x[0]/cm;
            extra.y_cm = x[1]/cm;
            extra.z_cm = x[2]/cm;
            extra.t_ns = tout/ns;
            extra.dist_cm = dist/cm;
            hiter->extra.push_back(extra);
         }
      }

      if (dEsouth/MeV > 0 && tsouth < MAX_TOF) {
         // add the hit on end=1 (south/bottom end of the bar)
         int merge_hit = 0;
         std::vector<GlueXHitFTOFbar::hitinfo_t>::iterator hiter;
         for (hiter = counter->hits.begin(); 
              hiter != counter->hits.end(); ++hiter)
         {
            if (hiter->end_ == 1) {
               if (fabs(hiter->t_ns*ns - tsouth) < TWO_HIT_TIME_RESOL) {
                  merge_hit = 1;
                  break;
               }
               else if (hiter->t_ns*ns > tsouth) {
                  break;
               }
            }
         }
         if (merge_hit) {
            // sum the charge, do energy weighting of the time --disabled for bad behavior, rtj--
            //hiter->t_ns = (hiter->dE_GeV * hiter->t_ns +
            //               dEsouth/GeV * tsouth/ns) /
            //(hiter->dE_GeV += dEsouth/GeV);

            // Use the time from the earlier hit but add the charge
            hiter->dE_GeV += dEsouth/GeV;
            if (hiter->t_ns*ns > tsouth) {
               hiter->t_ns = tsouth/ns;
            }

            std::vector<GlueXHitFTOFbar::hitextra_t>::reverse_iterator xiter;
            xiter = hiter->extra.rbegin();
            if (trackID != xiter->track_ || fabs(tin/ns - xiter->t_ns) > 0.1) {
               GlueXHitFTOFbar::hitextra_t extra;
               extra.track_ = trackID;
               extra.itrack_ = itrack;
               extra.ptype_G3 = g3type;
               extra.px_GeV = pin[0]/GeV;
               extra.py_GeV = pin[1]/GeV;
               extra.pz_GeV = pin[2]/GeV;
               extra.E_GeV = Ein/GeV;
               extra.x_cm = x[0]/cm;
               extra.y_cm = x[1]/cm;
               extra.z_cm = x[2]/cm;
               extra.t_ns = tout/ns;
               extra.dist_cm = dist/cm;
               hiter->extra.push_back(extra);
            }
         }
         else {
            // create new hit 
            hiter = counter->hits.insert(hiter, GlueXHitFTOFbar::hitinfo_t());
            hiter->end_ = 1;
            hiter->t_ns = tsouth/ns;
            hiter->dE_GeV = dEsouth/GeV;
            GlueXHitFTOFbar::hitextra_t extra;
            extra.track_ = trackID;
            extra.itrack_ = itrack;
            extra.ptype_G3 = g3type;
            extra.px_GeV = pin[0]/GeV;
            extra.py_GeV = pin[1]/GeV;
            extra.pz_GeV = pin[2]/GeV;
            extra.E_GeV = Ein/GeV;
            extra.x_cm = x[0]/cm;
            extra.y_cm = x[1]/cm;
            extra.z_cm = x[2]/cm;
            extra.t_ns = tout/ns;
            extra.dist_cm = dist/cm;
            hiter->extra.push_back(extra);
         }
      }
   }
   return true;
}

void GlueXSensitiveDetectorFTOF::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitFTOFbar*> *bars = fBarHitsMap->GetMap();
   std::map<int,GlueXHitFTOFpoint*> *points = fPointsMap->GetMap();
   if (bars->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitFTOFbar*>::iterator siter;
   std::map<int,GlueXHitFTOFpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << bars->size() << " bars with hits in the FTOF: "
             << G4endl;
      for (siter = bars->begin(); siter != bars->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the FTOF: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorFTOF::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getForwardTOFs().size() == 0)
      hitview.addForwardTOFs();
   hddm_s::ForwardTOF &forwardTOF = hitview.getForwardTOF();

   // Collect and output the ftofTruthHits
   for (siter = bars->begin(); siter != bars->end(); ++siter) {
      std::vector<GlueXHitFTOFbar::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].dE_GeV*1e3 <= THRESH_MEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::FtofCounterList counter = forwardTOF.addFtofCounters(1);
         counter(0).setPlane(siter->second->plane_);
         counter(0).setBar(siter->second->bar_);
         // first the end=0 hits
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            if (hits[ih].end_ == 0) {
               hddm_s::FtofTruthHitList thit = counter(0).addFtofTruthHits(1);
               thit(0).setEnd(hits[ih].end_);
               thit(0).setDE(hits[ih].dE_GeV);
               thit(0).setT(hits[ih].t_ns);
               for (int ihx=0; ihx < (int)hits[ih].extra.size(); ++ihx) {
                  hddm_s::FtofTruthExtraList xtra = thit(0).addFtofTruthExtras(1);
                  xtra(0).setItrack(hits[ih].extra[ihx].itrack_);
                  xtra(0).setPtype(hits[ih].extra[ihx].ptype_G3);
                  xtra(0).setPx(hits[ih].extra[ihx].px_GeV);
                  xtra(0).setPy(hits[ih].extra[ihx].py_GeV);
                  xtra(0).setPz(hits[ih].extra[ihx].pz_GeV);
                  xtra(0).setE(hits[ih].extra[ihx].E_GeV);
                  xtra(0).setX(hits[ih].extra[ihx].x_cm);
                  xtra(0).setY(hits[ih].extra[ihx].y_cm);
                  xtra(0).setZ(hits[ih].extra[ihx].z_cm);
                  xtra(0).setDist(hits[ih].extra[ihx].dist_cm);
               }
            }
         }
         // followed by the end=1 hits
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            if (hits[ih].end_ == 1) {
               hddm_s::FtofTruthHitList thit = counter(0).addFtofTruthHits(1);
               thit(0).setEnd(hits[ih].end_);
               thit(0).setDE(hits[ih].dE_GeV);
               thit(0).setT(hits[ih].t_ns);
               for (int ihx=0; ihx < (int)hits[ih].extra.size(); ++ihx) {
                  hddm_s::FtofTruthExtraList xtra = thit(0).addFtofTruthExtras(1);
                  xtra(0).setItrack(hits[ih].extra[ihx].itrack_);
                  xtra(0).setPtype(hits[ih].extra[ihx].ptype_G3);
                  xtra(0).setPx(hits[ih].extra[ihx].px_GeV);
                  xtra(0).setPy(hits[ih].extra[ihx].py_GeV);
                  xtra(0).setPz(hits[ih].extra[ihx].pz_GeV);
                  xtra(0).setE(hits[ih].extra[ihx].E_GeV);
                  xtra(0).setX(hits[ih].extra[ihx].x_cm);
                  xtra(0).setY(hits[ih].extra[ihx].y_cm);
                  xtra(0).setZ(hits[ih].extra[ihx].z_cm);
                  xtra(0).setDist(hits[ih].extra[ihx].dist_cm);
               }
            }
         }
         int hitscount = counter(0).getFtofTruthHits().size();
         if (hitscount > 2 * MAX_HITS_PER_BAR) {
            counter(0).deleteFtofTruthHits(-1, 2 * MAX_HITS_PER_BAR);
            G4cerr << "GlueXSensitiveDetectorFTOF::EndOfEvent warning: "
                   << "max bar hit count " << 2 * MAX_HITS_PER_BAR
                   << " exceeded, " << hitscount - 2 * MAX_HITS_PER_BAR
                   << " hits discarded." << G4endl;
         }
      }
   }

   // Collect and output the ftofTruthPoints
   GlueXUserEventInformation* eventinfo = (GlueXUserEventInformation*)info;
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::FtofTruthPointList point = forwardTOF.addFtofTruthPoints(1);
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

int GlueXSensitiveDetectorFTOF::GetIdent(std::string div, 
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

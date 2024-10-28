//
// GlueXSensitiveDetectorBCAL - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 25, 2016

#include "GlueXSensitiveDetectorBCAL.hh"
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

#define USE_ENERGY_WEIGHTED_TIMES 1

// Energy rescale factor to match reconstructed to generated
double GlueXSensitiveDetectorBCAL::SHOWER_ENERGY_SCALE_FACTOR = 1.008;

// Cutoff on the number of allowed hits per sector
int GlueXSensitiveDetectorBCAL::MAX_HITS = 100;

// Time resolution for full pulse digitization (unused at present)
double GlueXSensitiveDetectorBCAL::SIPM_TIME_BIN_WIDTH = 0.1;

// Light propagation parameters in barrel calorimeter
double GlueXSensitiveDetectorBCAL::ATTENUATION_LENGTH = 300.*cm;
double GlueXSensitiveDetectorBCAL::C_EFFECTIVE = 16.75*cm/ns;
double GlueXSensitiveDetectorBCAL::MODULE_FULL_LENGTH = 390.*cm;
double GlueXSensitiveDetectorBCAL::ATTENUATION_FULL_LENGTH = 
                     exp(-MODULE_FULL_LENGTH / ATTENUATION_LENGTH);
double GlueXSensitiveDetectorBCAL::THRESH_ATTENUATED_GEV =
                     (THRESH_MEV / 1e3) * ATTENUATION_FULL_LENGTH;

// Minimum hit time difference for two hits on the same cell
double GlueXSensitiveDetectorBCAL::TWO_HIT_TIME_RESOL = 50*ns;

// Minimum energy deposition for a hit
double GlueXSensitiveDetectorBCAL::THRESH_MEV = 1.;

int GlueXSensitiveDetectorBCAL::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorBCAL::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorBCAL::GlueXSensitiveDetectorBCAL(const G4String& name)
 : G4VSensitiveDetector(name),
   fCellsMap(0), fPointsMap(0)
{
   collectionName.insert("BCALCellHitsCollection");
   collectionName.insert("BCALPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the BCAL, you must delete all old
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
      std::map<string, float> bcal_parms;
      jcalib->Get("BCAL/bcal_parms", bcal_parms);
      THRESH_MEV = bcal_parms.at("BCAL_THRESH_MEV");
      TWO_HIT_TIME_RESOL = bcal_parms.at("BCAL_TWO_HIT_RESOL")*ns;
      MAX_HITS = bcal_parms.at("BCAL_MAX_HITS");
      std::map<string, float> mc_parms;
      jcalib->Get("BCAL/mc_parms", mc_parms);
      ATTENUATION_LENGTH = mc_parms.at("ATTEN_LENGTH")*cm;
      C_EFFECTIVE = mc_parms.at("C_EFFECTIVE")*cm/ns;
      SIPM_TIME_BIN_WIDTH = mc_parms.at("SiPM_tbin_width");

      G4cout << "BCAL: ALL parameters loaded from ccdb" << G4endl;

      ATTENUATION_FULL_LENGTH = exp(-MODULE_FULL_LENGTH /ATTENUATION_LENGTH);
      THRESH_ATTENUATED_GEV = (THRESH_MEV / 1e3) * ATTENUATION_FULL_LENGTH;
   }
}

GlueXSensitiveDetectorBCAL::GlueXSensitiveDetectorBCAL(
                     const GlueXSensitiveDetectorBCAL &src)
 : G4VSensitiveDetector(src),
   fCellsMap(src.fCellsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorBCAL &GlueXSensitiveDetectorBCAL::operator=(const
                                         GlueXSensitiveDetectorBCAL &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fCellsMap = src.fCellsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorBCAL::~GlueXSensitiveDetectorBCAL() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorBCAL::Initialize(G4HCofThisEvent* hce)
{
   fCellsMap = new
              GlueXHitsMapBCALcell(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapBCALpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fCellsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorBCAL::ProcessHits(G4Step* step, 
                                               G4TouchableHistory* RHhist)
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
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (touch->GetVolume()->GetName() == "BCL0") {
      if (track->GetCurrentStepNumber() == 1)
         trackinfo->SetGlueXHistory(1);
      if (trackinfo->GetGlueXHistory() == 0 &&
          xin[0] * pin[0] + xin[1] * pin[1] > 0 &&
          Ein > THRESH_MEV*MeV)
      {
         G4int key = fPointsMap->entries();
         GlueXHitBCALpoint* lastPoint = (*fPointsMap)[key - 1];
         // Limit bcal truthPoints to one per track
         if (lastPoint == 0 || lastPoint->track_ != trackID) {
            GlueXHitBCALpoint newPoint;
            newPoint.ptype_G3 = g3type;
            newPoint.track_ = trackID;
            newPoint.trackID_ = itrack;
            newPoint.primary_ = (track->GetParentID() == 0);
            newPoint.t_ns = t/ns;
            newPoint.z_cm = xin[2]/cm;
            newPoint.r_cm = xin.perp()/cm;
            newPoint.phi_rad = xin.phi();
            newPoint.px_GeV = pin[0]/GeV;
            newPoint.py_GeV = pin[1]/GeV;
            newPoint.pz_GeV = pin[2]/GeV;
            newPoint.E_GeV = Ein/GeV;
            fPointsMap->add(key, newPoint);
         }
 
         // The original HDGeant hits code for the BCal had a heavy-weight
         // recording system implemented to assign every hit to a particular
         // incident particle ID. This ID was unique to the BCal subsystem,
         // and not used anywhere else in the hits record. Each track that
         // entered the BCAL modules was assigned its own ID unless it was
         // within some rectangular window of an existing entry in the
         // incident_particle table, in which case it could supersede the
         // existing entry (if its energy was greater) or be ignored. When
         // BCAL hits were recorded, each was assigned the incident particle
         // id that it was closest to in the phi,z plane. There was no upper
         // limit on that distance, so the association could be considered
         // to be a weak link.
         //
         // At one time, it seems that the incident particle table was being
         // written out to the hddm record at the end of event processing
         // because there is a tag defined in the hddm record called 
         // bcalTruthIncidentParticle, but the present code in hitBCal.cc 
         // does not fill in this tag when it writes out the hits. I have
         // no idea whether this ID system is of any use or if it is an
         // artifact of some obsolete studies, but since it is being 
         // tracked and saved in the existing HDGeant bcal hits code, I
         // make some attempt to replicate that behavior here.
         // 
         // DO NOT use proximity of hits to incident particle entry points
         // in (phi,z) to associate hits to incident particles -- that is
         // a time waster and the causal linkage based on this is weak.
         // Instead, use an inheritance mechanism based on the G4 user
         // track information object. This user track info object gets
         // copied from parent track to secondaries through the entire
         // tracking workflow. The first ancestor that enters the BCal
         // gets assigned the next available BCAL incident ID, and that
         // is kept for that track and all of its descendants in the user
         // trackinfo object.
  
         trackinfo->AssignBCALincidentID(track);
      }
      return false;
   }

   // Post the hit to the hits map, ordered by sector index

   if (dEsum > 0) {
      int module = GetIdent("module", touch);
      int layer = GetIdent("layer", touch);
      int sector = GetIdent("sector", touch);
      int key = GlueXHitBCALcell::GetKey(module, layer, sector);
      GlueXHitBCALcell *cell = (*fCellsMap)[key];
      if (cell == 0) {
         GlueXHitBCALcell newcell(module, layer, sector);
         fCellsMap->add(key, newcell);
         cell = (*fCellsMap)[key];
      }

      // Add the hit to the bcal truth hits list, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitBCALcell::hitinfo_t>::iterator hiter;
      for (hiter = cell->hits.begin(); hiter != cell->hits.end(); ++hiter) {
         if (fabs(hiter->t_ns*ns - t) < TWO_HIT_TIME_RESOL) {
            merge_hit = 1;
            break;
         }
         else if (hiter->t_ns*ns > t) {
            break;
         }
      }
      if (merge_hit) {
         // Use the time from the earlier hit but add the energy deposition
#if USE_ENERGY_WEIGHTED_TIMES
         hiter->t_ns = (hiter->t_ns * hiter->E_GeV + 
                               t/ns * dEsum/GeV) /
                       (hiter->E_GeV + dEsum/GeV);
         hiter->zlocal_cm = (hiter->zlocal_cm * hiter->E_GeV +
                                    xlocal[2]/cm * dEsum/GeV) /
                            (hiter->E_GeV + dEsum/GeV);
#else
         if (hiter->t_ns*ns > t) {
            hiter->t_ns = t/ns;
            hiter->zlocal_cm = xlocal[2]/cm;
            hiter->incidentId_ = trackinfo->GetBCALincidentID();
         }
#endif
         // correction factor makes shower yields match hdgeant
         hiter->E_GeV += dEsum/GeV * SHOWER_ENERGY_SCALE_FACTOR;
      }
      else {
         // create new hit 
         hiter = cell->hits.insert(hiter, GlueXHitBCALcell::hitinfo_t());
         // correction factor makes shower yields match hdgeant
         hiter->E_GeV = dEsum/GeV * SHOWER_ENERGY_SCALE_FACTOR;
         hiter->t_ns = t/ns;
         hiter->zlocal_cm = xlocal[2]/cm;
         hiter->incidentId_ = trackinfo->GetBCALincidentID();
      }

      // Add the hit to the upstream sipm hits list, with strict time ordering

      double udist = MODULE_FULL_LENGTH/2 + xlocal[2];
      double dEup = dEsum * exp(-udist/ATTENUATION_LENGTH);
      double tup = t + udist / C_EFFECTIVE;
      merge_hit = 0;
      for (hiter = cell->hits.begin(); hiter != cell->hits.end(); ++hiter) {
         if (fabs(hiter->tup_ns*ns - tup) < TWO_HIT_TIME_RESOL) {
            merge_hit = 1;
            break;
         }
         else if (hiter->tup_ns*ns > tup) {
            break;
         }
      }
      if (merge_hit) {
         // Use the time from the earlier hit but add the energy deposition
#if USE_ENERGY_WEIGHTED_TIMES
         hiter->tup_ns = (hiter->tup_ns * hiter->Eup_GeV + 
                                 tup/ns * dEup/GeV) /
                       (hiter->Eup_GeV + dEup/GeV);
#else
         if (hiter->tup_ns*ns > tup) {
            hiter->tup_ns = tup/ns;
         }
#endif
         // correction factor makes shower yields match hdgeant
         hiter->Eup_GeV += dEup/GeV * SHOWER_ENERGY_SCALE_FACTOR;
      }
      else {
         // create new hit 
         hiter = cell->hits.insert(hiter, GlueXHitBCALcell::hitinfo_t());
         // correction factor makes shower yields match hdgeant
         hiter->Eup_GeV = dEup/GeV * SHOWER_ENERGY_SCALE_FACTOR;
         hiter->tup_ns = tup/ns;
      }

      // Add the hit to the downstream sipm hits list, with strict time ordering

      double ddist = MODULE_FULL_LENGTH/2 - xlocal[2];
      double dEdown = dEsum * exp(-ddist/ATTENUATION_LENGTH);
      double tdown = t + ddist / C_EFFECTIVE;
      merge_hit = 0;
      for (hiter = cell->hits.begin(); hiter != cell->hits.end(); ++hiter) {
         if (fabs(hiter->tdown_ns*ns - tdown) < TWO_HIT_TIME_RESOL) {
            merge_hit = 1;
            break;
         }
         else if (hiter->tdown_ns*ns > tdown) {
            break;
         }
      }
      if (merge_hit) {
         // Use the time from the earlier hit but add the energy deposition
#if USE_ENERGY_WEIGHTED_TIMES
         hiter->tdown_ns = (hiter->tdown_ns * hiter->Edown_GeV + 
                                   tdown/ns * dEdown/GeV) /
                       (hiter->Edown_GeV + dEdown/GeV);
#else
         if (hiter->tdown_ns*ns > tdown) {
            hiter->tdown_ns = tdown/ns;
         }
#endif
         // correction factor makes shower yields match hdgeant
         hiter->Edown_GeV += dEdown/GeV * SHOWER_ENERGY_SCALE_FACTOR;
      }
      else {
         // create new hit 
         hiter = cell->hits.insert(hiter, GlueXHitBCALcell::hitinfo_t());
         // correction factor makes shower yields match hdgeant
         hiter->Edown_GeV = dEdown/GeV * SHOWER_ENERGY_SCALE_FACTOR;
         hiter->tdown_ns = tdown/ns;
      }
   }
   return true;
}

void GlueXSensitiveDetectorBCAL::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitBCALcell*> *cells = fCellsMap->GetMap();
   std::map<int,GlueXHitBCALpoint*> *points = fPointsMap->GetMap();
   if (cells->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitBCALcell*>::iterator citer;
   std::map<int,GlueXHitBCALpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << cells->size() << " cells with hits in the BCAL: "
             << G4endl;
      for (citer = cells->begin(); citer != cells->end(); ++citer)
         citer->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth showers in the BCAL: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorBCAL::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getBarrelEMcals().size() == 0)
      hitview.addBarrelEMcals();
   hddm_s::BarrelEMcal &barrelEMcal = hitview.getBarrelEMcal();

   // Collect and output the bcalTruthHits
   for (citer = cells->begin(); citer != cells->end(); ++citer) {
      std::vector<GlueXHitBCALcell::hitinfo_t> &hits = citer->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].E_GeV <= THRESH_MEV/1e3) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::BcalCellList cell = barrelEMcal.addBcalCells(1);
         cell(0).setModule(citer->second->module_);
         cell(0).setLayer(citer->second->layer_);
         cell(0).setSector(citer->second->sector_);
         int hitscount = hits.size();
         if (hitscount > MAX_HITS) {
            hitscount = MAX_HITS;
            G4cerr << "GlueXSensitiveDetectorBCAL::EndOfEvent warning: "
                   << "max hit count " << MAX_HITS << " exceeded, "
                   << hits.size() - hitscount << " hits discarded."
                   << G4endl;
         }
         for (int ih=0; ih < hitscount; ++ih) {
            hddm_s::BcalTruthHitList thit = cell(0).addBcalTruthHits(1);
            thit(0).setE(hits[ih].E_GeV);
            thit(0).setT(hits[ih].t_ns);
            thit(0).setZLocal(hits[ih].zlocal_cm);
            thit(0).setIncident_id(hits[ih].incidentId_);
            if (hits[ih].Eup_GeV >= THRESH_MEV/1e3) {
               hddm_s::BcalSiPMUpHitList uhit = cell(0).addBcalSiPMUpHits(1);
               uhit(0).setE(hits[ih].Eup_GeV);
               uhit(0).setT(hits[ih].tup_ns);
            }
            if (hits[ih].Edown_GeV >= THRESH_MEV/1e3) {
               hddm_s::BcalSiPMDownHitList dhit = cell(0).addBcalSiPMDownHits(1);
               dhit(0).setE(hits[ih].Edown_GeV);
               dhit(0).setT(hits[ih].tdown_ns);
            }
         }
      }
   }

   // Collect and output the bcalTruthShowers
   GlueXUserEventInformation* eventinfo = (GlueXUserEventInformation*)info;
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::BcalTruthShowerList point = barrelEMcal.addBcalTruthShowers(1);
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

int GlueXSensitiveDetectorBCAL::GetIdent(std::string div, 
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

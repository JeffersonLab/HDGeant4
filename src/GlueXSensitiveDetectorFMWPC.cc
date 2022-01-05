//
// GlueXSensitiveDetectorFMWPC - class implementation
//
// author: richard.t.jones at uconn.edu
// version: november 29, 2016

#include "GlueXSensitiveDetectorFMWPC.hh"
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

// Cutoff on the total number of allowed hits
int GlueXSensitiveDetectorFMWPC::MAX_HITS = 100;

// Minimum hit time difference for two hits on the same wire
double GlueXSensitiveDetectorFMWPC::TWO_HIT_TIME_RESOL = 400*ns;

// Minimum photoelectron count for a hit
double GlueXSensitiveDetectorFMWPC::THRESH_KEV = 0.;

// Coordinate of wire 0, transverse to wire direction
double GlueXSensitiveDetectorFMWPC::WIRE_OFFSET = -(73.000*1.016)*cm;

// Minimum photoelectron count for a hit
double GlueXSensitiveDetectorFMWPC::WIRE_PITCH = 1.016*cm;

int GlueXSensitiveDetectorFMWPC::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorFMWPC::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorFMWPC::GlueXSensitiveDetectorFMWPC(const G4String& name)
 : G4VSensitiveDetector(name),
   fWireHitsMap(0), fPointsMap(0)
{
   collectionName.insert("FMWPCWireHitsCollection");
   collectionName.insert("FMWPCPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the FMWPC, you must delete all old
   // objects of this class and create new ones.

   G4AutoLock wirerier(&fMutex);
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
         G4cout << "FMWPC: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorFMWPC::GlueXSensitiveDetectorFMWPC(
                     const GlueXSensitiveDetectorFMWPC &src)
 : G4VSensitiveDetector(src),
   fWireHitsMap(src.fWireHitsMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock wirerier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorFMWPC &GlueXSensitiveDetectorFMWPC::operator=(const
                                         GlueXSensitiveDetectorFMWPC &src)
{
   G4AutoLock wirerier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fWireHitsMap = src.fWireHitsMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorFMWPC::~GlueXSensitiveDetectorFMWPC() 
{
   G4AutoLock wirerier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorFMWPC::Initialize(G4HCofThisEvent* hce)
{
   fWireHitsMap = new
              GlueXHitsMapFMWPCwire(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapFMWPCpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fWireHitsMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorFMWPC::ProcessHits(G4Step* step, 
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
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
                                          track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (trackinfo->GetGlueXHistory() == 0 && itrack > 0 && xin.dot(pin) > 0) {
      G4int key = fPointsMap->entries();
      GlueXHitFMWPCpoint* lastPoint = (*fPointsMap)[key - 1];
      if (lastPoint == 0 || lastPoint->track_ != trackID ||
          fabs(lastPoint->t_ns - t/ns) > 0.1 ||
          fabs(lastPoint->x_cm - x[0]/cm) > 2. ||
          fabs(lastPoint->y_cm - x[1]/cm) > 2. ||
          fabs(lastPoint->z_cm - x[2]/cm) > 2.)
      {
         int pdgtype = track->GetDynamicParticle()->GetPDGcode();
         int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
         GlueXHitFMWPCpoint newPoint;
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

   // Post the hit to the hits map, ordered by plane,wire,end index

   if (dEsum > 0) {

      // HDDS geometry does not include wires nor plane rotations. Assume 1 cm
      // wire spacing with wires at 0.5cm on either side of the beamline at
      // x,y = 0,0. Also assume odd number planes have wires in the vertical
      // direction and even numbered planes have wires in the horizontal direction.
      // Vertical wires start with wire 1 at x=-71.5 and wire 144 at x=+71.5
      // (the gas volume ends at x=+/-72.0). Horizontal wires start with wire 1
      // at y=-71.5 (i.e. closest to the ground) and wire 144 at y=+71.5
      // (i.e. closest to the sky).

      int layer = GetIdent("layer", touch);
      int wire = 0;
      if (layer % 2 != 0) {
         // Vertical wires
	wire = floor((x[0] - WIRE_OFFSET)/WIRE_PITCH);
      }
      else {
         // Horizontal wires
	wire = floor((x[1] - WIRE_OFFSET)/WIRE_PITCH);
      }
      //cout<<"MWPC: layer/wire = "<<layer<<" / "<<wire<<endl;

      if (wire < 1 || wire > 144)
         return false;
      
      int key = GlueXHitFMWPCwire::GetKey(layer, wire);
      GlueXHitFMWPCwire *counter = (*fWireHitsMap)[key];
      if (counter == 0) {
         GlueXHitFMWPCwire newwire(layer, wire);
         fWireHitsMap->add(key, newwire);
         counter = (*fWireHitsMap)[key];
      }

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitFMWPCwire::hitinfo_t>::iterator hiter;
      for (hiter = counter->hits.begin(); hiter != counter->hits.end(); ++hiter) {
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
         hiter->dE_keV += dEsum/keV;
         hiter->dx_cm += dx.mag()/cm;
         if (hiter->t_ns*ns > t) {
            hiter->t_ns = t/ns;
         }
      }
      else if ((int)counter->hits.size() < MAX_HITS) {
         // create new hit 
         hiter = counter->hits.insert(hiter, GlueXHitFMWPCwire::hitinfo_t());
         hiter->dE_keV = dEsum/keV;
         hiter->dx_cm = dx.mag()/cm;
         hiter->t_ns = t/ns;
      }
      else {
         G4cerr << "GlueXSensitiveDetectorFMWPC::ProcessHits error: "
                << "max hit count " << MAX_HITS
                << " exceeded, truncating!"
                << G4endl;
      }
   }
   return true;
}

void GlueXSensitiveDetectorFMWPC::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitFMWPCwire*> *wires = fWireHitsMap->GetMap();
   std::map<int,GlueXHitFMWPCpoint*> *points = fPointsMap->GetMap();
   if (wires->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitFMWPCwire*>::iterator siter;
   std::map<int,GlueXHitFMWPCpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << wires->size() << " wires with hits in the FMWPC: "
             << G4endl;
      for (siter = wires->begin(); siter != wires->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the FWMPC: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorFMWPC::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getForwardMWPCs().size() == 0)
      hitview.addForwardMWPCs();
   hddm_s::ForwardMWPC &fmwpc = hitview.getForwardMWPC();

   // Collect and output the fmwpcTruthHits
   for (siter = wires->begin(); siter != wires->end(); ++siter) {
      std::vector<GlueXHitFMWPCwire::hitinfo_t> &hits = siter->second->hits;
      // apply a pulse height threshold cut
      for (unsigned int ih=0; ih < hits.size(); ++ih) {
         if (hits[ih].dE_keV <= THRESH_KEV) {
            hits.erase(hits.begin() + ih);
            --ih;
         }
      }
      if (hits.size() > 0) {
         hddm_s::FmwpcChamberList wire = fmwpc.addFmwpcChambers(1);
         wire(0).setLayer(siter->second->layer_);
         wire(0).setWire(siter->second->wire_);
         for (int ih=0; ih < (int)hits.size(); ++ih) {
            hddm_s::FmwpcTruthHitList thit = wire(0).addFmwpcTruthHits(1);
            thit(0).setDE(hits[ih].dE_keV);
            thit(0).setDx(hits[ih].dx_cm);
            thit(0).setT(hits[ih].t_ns);
         }
      }
   }

   // Collect and output the fmwpcTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::FmwpcTruthPointList point = fmwpc.addFmwpcTruthPoints(1);
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
   }
}

int GlueXSensitiveDetectorFMWPC::GetIdent(std::string div, 
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

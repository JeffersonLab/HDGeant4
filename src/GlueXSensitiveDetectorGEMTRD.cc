//
// GlueXSensitiveDetectorGEMTRD - class implementation
//
// author: staylor at jlab.org
// version: february 7, 2022

#include "GlueXSensitiveDetectorGEMTRD.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "HddmOutput.hh"

#include <CLHEP/Random/RandPoisson.h>
#include <Randomize.hh>

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <JANA/JApplication.h>

// Cutoff on the total number of allowed hits and traces
int GlueXSensitiveDetectorGEMTRD::MAX_TRACES = 20;
int GlueXSensitiveDetectorGEMTRD::MAX_HITS = 100;

const double fC = 1e-15 * coulomb;
const double GlueXSensitiveDetectorGEMTRD::ELECTRON_CHARGE = 1.6022e-4*fC;

// Parameters for setting signal pulse height
double GlueXSensitiveDetectorGEMTRD::GAS_GAIN = 1e5;

// Average number of secondary ion pairs for 90/10 Xe/CO2 mixture
double GlueXSensitiveDetectorGEMTRD::N_SECOND_PER_PRIMARY = 5.6; 

// Average energy needed to produce an ion pair for 90/10 mixture
double GlueXSensitiveDetectorGEMTRD::W_EFF_PER_ION = 22.4*eV;

int GlueXSensitiveDetectorGEMTRD::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorGEMTRD::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorGEMTRD::GlueXSensitiveDetectorGEMTRD(const G4String& name)
 : G4VSensitiveDetector(name),
   fTracesMap(0), fPointsMap(0)
{
   collectionName.insert("GEMTRDTracesCollection");
   collectionName.insert("GEMTRDPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the properties of hits in the GEMTRD, you must delete all old
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
      if (japp == 0) {   // dummy
         jcalib = 0;
         G4cout << "GEMTRD: ALL parameters loaded from ccdb" << G4endl;
      }
   }
}

GlueXSensitiveDetectorGEMTRD::GlueXSensitiveDetectorGEMTRD(
                     const GlueXSensitiveDetectorGEMTRD &src)
 : G4VSensitiveDetector(src),
   fTracesMap(src.fTracesMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorGEMTRD &GlueXSensitiveDetectorGEMTRD::operator=(const
                                         GlueXSensitiveDetectorGEMTRD &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fTracesMap = src.fTracesMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorGEMTRD::~GlueXSensitiveDetectorGEMTRD() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorGEMTRD::Initialize(G4HCofThisEvent* hce)
{
   fTracesMap = new
              GlueXHitsMapGEMTRDtrace(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
              GlueXHitsMapGEMTRDpoint(SensitiveDetectorName, collectionName[1]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fTracesMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[1]), fPointsMap);
}

G4bool GlueXSensitiveDetectorGEMTRD::ProcessHits(G4Step* step, 
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
    GlueXHitGEMTRDpoint* lastPoint = (*fPointsMap)[key - 1];
    if (lastPoint == 0 || lastPoint->track_ != trackID ||
	fabs(lastPoint->t_ns - t/ns) > 0.1 ||
	fabs(lastPoint->x_cm - x[0]/cm) > 0.1 ||
	fabs(lastPoint->y_cm - x[1]/cm) > 0.1 ||
	fabs(lastPoint->z_cm - x[2]/cm) > 0.1)
      {
	int pdgtype = track->GetDynamicParticle()->GetPDGcode();
	int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
	GlueXHitGEMTRDpoint newPoint;
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
  
  // Post the trace info to the traces map, ordered by layer
  
  if (dEsum > 0) {
    int layer = GetIdent("layer", touch);
    int key = GlueXHitGEMTRDtrace::GetKey(layer);
    GlueXHitGEMTRDtrace *counter = (*fTracesMap)[key];
    if (counter == 0) {
      GlueXHitGEMTRDtrace newtrace(layer);
      fTracesMap->add(key, newtrace);
      counter = (*fTracesMap)[key];
    }
    
    // Add the trace to the traces vector
    std::vector<GlueXHitGEMTRDtrace::traceinfo_t>::iterator titer=counter->traces.end();
    if ((int)counter->traces.size() < MAX_TRACES) {
      // create new hit 
      titer = counter->traces.insert(titer, GlueXHitGEMTRDtrace::traceinfo_t());
      titer->dE_keV = dEsum/keV;
      titer->x_cm = xin[0]/cm;
      titer->y_cm = xin[1]/cm;
      titer->dxdz = dx[0]/dx[2];
      titer->dydz = dx[1]/dx[2];
      titer->t_ns = tin/ns;
    }
    else {
      G4cerr << "GlueXSensitiveDetectorGEMTRD::ProcessHits error: "
	     << "max trace count " << MAX_TRACES
	     << " exceeded, truncating!"
	     << G4endl;
    }
  }
  return true;
}

void GlueXSensitiveDetectorGEMTRD::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitGEMTRDtrace*> *traces = fTracesMap->GetMap();
   std::map<int,GlueXHitGEMTRDpoint*> *points = fPointsMap->GetMap();
   if (traces->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitGEMTRDtrace*>::iterator siter;
   std::map<int,GlueXHitGEMTRDpoint*>::iterator piter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << traces->size() << " traces in the GEMTRD: "
             << G4endl;
      for (siter = traces->begin(); siter != traces->end(); ++siter)
         siter->second->Print();

      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the GEMTRD: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorGEMTRD::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getGEMTRDs().size() == 0)
      hitview.addGEMTRDs();
   hddm_s::GEMTRD &gemtrd = hitview.getGEMTRD();

   // Collect and output the gemtrdTruthHits
   for (siter = traces->begin(); siter != traces->end(); ++siter) {
     std::vector<GlueXHitGEMTRDtrace::traceinfo_t> &traces = siter->second->traces;
     hddm_s::GemtrdChamberList chamber = gemtrd.addGemtrdChambers(1);
     chamber(0).setLayer(siter->second->layer_);
     double charge_per_ion_pair=GAS_GAIN*ELECTRON_CHARGE/fC;
     for (int it=0; it < (int)traces.size(); ++it) {
       // Position and direction at entry to gas volume
       double x0=traces[it].x_cm;
       double y0=traces[it].y_cm;
       double dxdz=traces[it].dxdz;
       double dydz=traces[it].dydz;
       // Generate the number of primary ion pairs produced in the gas    
       double dE=traces[it].dE_keV;
       double n_p_mean=dE/W_EFF_PER_ION/(1.+N_SECOND_PER_PRIMARY);
       int n_p=CLHEP::RandPoisson::shoot(n_p_mean);
       for (int ip=0;ip<n_p;ip++){
	 hddm_s::GemtrdTruthHitList thit = chamber(0).addGemtrdTruthHits(1);
	 // Generate a cluster for this primary pair
	 int n_s=CLHEP::RandPoisson::shoot(N_SECOND_PER_PRIMARY);
	 thit(0).setQ(charge_per_ion_pair*n_s);
	 // Randomly generate the drift distance for this cluster in units of 
	 // the chamber width
	 double drift_fraction=G4UniformRand();
	 // Use simple linear time-to-distance relationship for now
	 double drift_time=800*ns*drift_fraction;
	 thit(0).setT(traces[it].t_ns+drift_time);
	 // position in x and y
	 double dz=2.0*cm*(1.-drift_fraction);
	 thit(0).setX(x0+dxdz*dz);
	 thit(0).setY(y0+dydz*dz);
       }
     }
   }
   
   // Collect and output the gemtrdTruthPoints
   for (piter = points->begin(); piter != points->end(); ++piter) {
      hddm_s::GemtrdTruthPointList point = gemtrd.addGemtrdTruthPoints(1);
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

int GlueXSensitiveDetectorGEMTRD::GetIdent(std::string div, 
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

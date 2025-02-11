//
// GlueXSensitiveDetectorDIRC - class implementation
//
// author: richard.t.jones at uconn.edu
// version: november 28, 2016

#include "GlueXSensitiveDetectorDIRC.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"
#include "HddmOutput.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4TransportationManager.hh"
#include "G4ParallelWorldProcess.hh"

#include <JANA/JApplication.h>
#include <JANA/Calibrations/JCalibrationManager.h>


// Geant3-style particle index for optical photon
#define OPTICAL_PHOTON 50

// Cutoff on the total number of allowed hits
int GlueXSensitiveDetectorDIRC::MAX_HITS = 1000;
int GlueXSensitiveDetectorDIRC::MAX_PIXELS = 6912;

// Minimum hit time difference for two hits on the same tube
double GlueXSensitiveDetectorDIRC::TWO_HIT_TIME_RESOL = 50*ns;

int GlueXSensitiveDetectorDIRC::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorDIRC::fMutex = G4MUTEX_INITIALIZER;

TGraph *GlueXSensitiveDetectorDIRC::fDetEff = 0;

GlueXSensitiveDetectorDIRC::GlueXSensitiveDetectorDIRC(const G4String& name)
  : G4VSensitiveDetector(name)
{
  // The rest of this only needs to happen once, the first time an object
  // of this type is instantiated for this configuration of geometry and
  // fields. If the geometry or fields change in such a way as to modify
  // the drift-time properties of hits in the DIRC, you must delete all old
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
    JCalibration *jcalib = japp->GetService<JCalibrationManager>()->GetJCalibration(runno);
    if (japp == 0) {   // dummy
      jcalib = 0;
      G4cout << "DIRC: ALL parameters loaded from ccdb" << G4endl;
    }
  }

  GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
  std::map<int, int> dirclutpars;
  if (user_opts->Find("DIRCLUT", dirclutpars)) {
    fLutId = dirclutpars[1];
  }
  else {
    fLutId = 100;
  }
  std::map<int, int> dircledpars;
  if (user_opts->Find("DIRCLED", dircledpars)){
    fLED = true;
  }
  else {
    fLED = false;
  }

}

GlueXSensitiveDetectorDIRC::GlueXSensitiveDetectorDIRC(const GlueXSensitiveDetectorDIRC &src)
 : G4VSensitiveDetector(src)
{
  G4AutoLock barrier(&fMutex);
  ++instanceCount;
}

GlueXSensitiveDetectorDIRC &GlueXSensitiveDetectorDIRC::operator=(const
                                  GlueXSensitiveDetectorDIRC &src)
{
  G4AutoLock barrier(&fMutex);
  *(G4VSensitiveDetector*)this = src;
  return *this;
}

GlueXSensitiveDetectorDIRC::~GlueXSensitiveDetectorDIRC() 
{
  G4AutoLock barrier(&fMutex);
  --instanceCount;
}

void GlueXSensitiveDetectorDIRC::Initialize(G4HCofThisEvent* hce)
{
}

G4bool GlueXSensitiveDetectorDIRC::ProcessHits(G4Step* step, 
                                               G4TouchableHistory* ROhist)
{
  
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
  const G4TouchableHistory* touch_hist = (G4TouchableHistory*)touch;
  const G4AffineTransform &local_from_global = touch->GetHistory()->GetTopTransform();
  G4ThreeVector xlocal = local_from_global.TransformPoint(x);
  
  // For particles that range out inside the active volume, the
  // "out" time may sometimes be set to something enormously high.
  // This screws up the hit. Check for this case here by looking
  // at tout and making sure it is less than 1 second. If it's
  // not, then just use tin for "t".

  if (tout > 1.0*s) t = tin;

  // Post the hit to the points list in the
  // order of appearance in the event simulation.

  G4Track *track = step->GetTrack();
  GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
    track->GetUserInformation();
  int itrack = trackinfo->GetGlueXTrackID();
  G4String volname = touch->GetVolume()->GetName();  
  
  // radiator volume: BNNM (NN = bar number 0-47 and M is sub-bar character A-D)
  int ibar = 10*((int)volname(1,1)(0)-48)+(int)volname(2,1)(0)-48;   // this is nasty, but it works
  if (volname(0,1)(0) == 'B' && ibar >= 0 && ibar < 48) { 

    if (trackinfo->GetGlueXHistory() == 0 && itrack > 0 && xin.dot(pin) > 0) {
      int pdgtype = track->GetDynamicParticle()->GetPDGcode();
      int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
      
      GlueXHitDIRCBar barhit;
      barhit.E_GeV = Ein/GeV;
      barhit.t_ns = t/ns;
      barhit.x_cm = x[0]/cm;
      barhit.y_cm = x[1]/cm;
      barhit.z_cm = x[2]/cm;
      barhit.px_GeV = pin[0]/GeV;
      barhit.py_GeV = pin[1]/GeV;
      barhit.pz_GeV = pin[2]/GeV;
      barhit.pdg = g3type;
      barhit.bar = ibar; // from HDDS geometry
      barhit.track = itrack; // track id of the charged particle
      fHitsBar.push_back(barhit);
    }
    return true;
  }

  // wedge and mirrors volumes

  if (volname == "WM1N" || volname == "WM2N" || volname == "WM1S" || volname == "WM2S" ||
      volname == "FTMN" || volname == "FTMS" ||
      volname == "TM1N" || volname == "TM2N" || volname == "TM3N" ||
      volname == "TM1S" || volname == "TM2S" || volname == "TM3S" ||
      volname == "SM1N" || volname == "SM2N" || volname == "SM1S" || volname == "SM2S" ||
      volname == "OWDG" ||
      (volname(0,1)(0) == 'A' && volname(0,1)(1) == 'G') )
  {

    GlueXHitDIRCWob wobhit;
    wobhit.track = track->GetTrackID();

    // store normal to the closest boundary
    G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
    std::vector<G4Navigator*>::iterator iNav =
      G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();

    G4bool valid;
    G4ThreeVector localNormal = (iNav[hNavId])->GetLocalExitNormal(&valid);
    if (valid){
      int mid=-1;
      if (volname == "OWDG") {
        if (localNormal.y()<-0.999)
          mid=1;
        else if (localNormal.y()>0.999)
          mid=2;
        else if(localNormal.z()>0.999)
          mid=3;
        else if(fabs(localNormal.z()+0.86)<0.01)
          mid=4;
      }
      if (volname == "SM1N" || volname == "SM1S")
        mid = 5;
      if (volname == "SM2N" || volname == "SM2S")
        mid = 6;
      if (volname == "WM1N" || volname == "WM1S")
        mid = 7;
      if (volname == "WM2N" || volname == "WM2S")
        mid = 8;
      if (volname == "FTMN" || volname == "FTMS")
        mid = 0;
      if (volname == "TM1N" || volname == "TM1S")
        mid=91;
      if (volname == "TM2N" || volname == "TM2S")
        mid=92;
      if (volname == "TM3N" || volname == "TM3S")
        mid=93;
      if ((volname(0,1)(0) == 'A' && volname(0,1)(1) == 'G'))
        mid=100;

      if (mid!=-1) {
        G4double normalId = mid;// localNormal.x() + 10*localNormal.y() + 100*localNormal.z();
        wobhit.normalId = normalId;
        fHitsWob.push_back(wobhit);
      }
    }
    return true;
  }
  
  // PMT's pixel volume
  if (volname == "PIXV") {

      // fix propagation speed for op
      double tracklen=track->GetTrackLength()/cm;
      double en=Ein/GeV;
      double refindex= 1.43603+0.0132404*en-0.00225287*en*en+0.000500109*en*en*en;  
      double time_fixed=tracklen/(29.9792458/refindex);

#ifdef DIRC_CHECK_PROPAGATION_TIME
      double l_QZBL = 9.1 + 0.96 + 122.5; // 2 * (4*122.5) // 2 * bar length + wedge + window
      double l_EPOTEK = 0.005 + 8 * 0.005; // window+wedge glue + 6 * bar joing glue
      double l_AIR = 2 * 0.01; // air gap to mirror
      double l_H2O = tracklen - l_QZBL - l_EPOTEK - l_AIR;
      
      // hard coded propagation time for 3.5 eV OpticalPhoton
      double angle = 45/180*3.14159;
      double time_propagated = l_QZBL/(29.9792458/1.476)/cos(angle);
      time_propagated += l_EPOTEK/(29.9792458/1.616)/cos(angle);
      time_propagated += l_AIR/(29.9792458)/cos(angle);
      time_propagated += l_H2O/(29.9792458/1.343);
      G4cout<<"Propagated time = "<<time_propagated<<" and measured time = "<<t/ns<<G4endl;
#endif

      GlueXHitDIRCPmt pmthit;
      pmthit.E_GeV = Ein/GeV;
      pmthit.t_ns = t/ns;
      pmthit.t_fixed_ns = time_fixed;
      pmthit.x_cm = x[0]/cm;
      pmthit.y_cm = x[1]/cm;
      pmthit.z_cm = x[2]/cm;
      
      double box = touch_hist->GetReplicaNumber(1)-1; // [0,1]
      double pix = touch_hist->GetReplicaNumber(0)-1; // [0,6911]

      pmthit.ch = box*MAX_PIXELS + pix;

      pmthit.key_bar = -999;
      for (unsigned int i=0;i<fHitsBar.size();i++){ // get bar hit from parent track
        if (fHitsBar[i].track == track->GetParentID()) {
          pmthit.key_bar = fHitsBar[i].bar;
        }
      }    
      pmthit.track = track->GetParentID();

      int64_t pathId1 = 0;
      int64_t pathId2 = 0;
      int mid, refl=0;
      pmthit.bbrefl = false;
      for (unsigned int i=0;i<fHitsWob.size();i++) {
        if (fHitsWob[i].track == track->GetTrackID()) {
          mid =fHitsWob[i].normalId;
          if (mid>=100) {
            pmthit.bbrefl = true;
            continue;
          }
          refl++;
          if (refl <= 18) {
            if (mid<10)
              pathId1 = pathId1*10 + mid;
            else
              pathId1 = pathId1*100 + mid;
          }
          else if (refl > 18 && refl < 26) {
            if (mid<10)
              pathId2 = pathId2*10 + mid;
            else
              pathId2 = pathId2*100 + mid;
          }
        }
      }
      int64_t pathId=pathId1+pathId2;
      if (refl>19)
        pathId *=-1;
      
      pmthit.path = pathId;
      pmthit.refl = refl;
      
      if (fLutId<48) {
        //G4ThreeVector vmom = track->GetVertexMomentumDirection();
        pmthit.key_bar = fLutId;
      }
      fHitsPmt.push_back(pmthit);
  }
  return true;
}

void GlueXSensitiveDetectorDIRC::EndOfEvent(G4HCofThisEvent*)
{
  if ((fHitsBar.size() == 0 && !(fLutId<48) && !fLED) || (fHitsBar.size() == 0 && !fLED) || fHitsPmt.size() == 0 || (fHitsWob.size() == 0 && !fLED))
  {
  
    fHitsBar.clear();
    fHitsPmt.clear();
    fHitsWob.clear();
    return;
  }
  
  if (verboseLevel > 1) { 
    G4cout << G4endl
       << "--------> Hits Collection: in this event there are "
       << fHitsBar.size() << " bar hits:"
       << G4endl;
    for(unsigned int h=0; h<fHitsBar.size(); h++)
      fHitsBar[h].Print();
    
    G4cout << G4endl
       << "--------> Hits Collection: in this event there are "
       << fHitsPmt.size() << " PMT hits: "
       << G4endl;
    for(unsigned int h=0; h<fHitsPmt.size(); h++)
      fHitsPmt[h].Print();
  }

  // pack hits into ouptut hddm record
 
  G4EventManager* mgr = G4EventManager::GetEventManager();
  G4VUserEventInformation* info = mgr->GetUserInformation();
  hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
  if (record == 0) {
    G4cerr << "GlueXSensitiveDetectorDIRC::EndOfEvent error - "
           << "hits seen but no output hddm record to save them into, "
       << "cannot continue!" << G4endl;
    exit(1);
  }

  if (record->getPhysicsEvents().size() == 0) record->addPhysicsEvents();
  if (record->getHitViews().size() == 0) record->getPhysicsEvent().addHitViews();
  hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
  if (hitview.getDIRCs().size() == 0) hitview.addDIRCs();
  hddm_s::DIRC &dirc = hitview.getDIRC();

  // Collect and output the BarHits
  for(unsigned int h=0; h<fHitsBar.size(); h++){
    hddm_s::DircTruthBarHitList bhit = dirc.addDircTruthBarHits(1);
    bhit(0).setE(fHitsBar[h].E_GeV);
    bhit(0).setT(fHitsBar[h].t_ns);
    bhit(0).setX(fHitsBar[h].x_cm);
    bhit(0).setY(fHitsBar[h].y_cm);
    bhit(0).setZ(fHitsBar[h].z_cm);
    bhit(0).setPx(fHitsBar[h].px_GeV);
    bhit(0).setPy(fHitsBar[h].py_GeV);
    bhit(0).setPz(fHitsBar[h].pz_GeV);
    bhit(0).setPdg(fHitsBar[h].pdg);
    bhit(0).setBar(fHitsBar[h].bar);
    bhit(0).setTrack(fHitsBar[h].track);
  }

  // Collect and output the DircTruthPmtHit
  int hitscount = fHitsPmt.size();
  if (hitscount > MAX_HITS) {
    hitscount = MAX_HITS;
    G4cerr << "GlueXSensitiveDetectorDIRC::EndOfEvent warning: "
           << "max hit count " << MAX_HITS << " exceeded, "
           << fHitsPmt.size() - hitscount << " hits discarded."
           << G4endl;
  }
  for(int h=0; h<hitscount; h++){
    hddm_s::DircTruthPmtHitList mhit = dirc.addDircTruthPmtHits(1);
    mhit(0).setE(fHitsPmt[h].E_GeV);
    mhit(0).setT(fHitsPmt[h].t_ns);
    mhit(0).setX(fHitsPmt[h].x_cm);
    mhit(0).setY(fHitsPmt[h].y_cm);
    mhit(0).setZ(fHitsPmt[h].z_cm);
    mhit(0).setCh(fHitsPmt[h].ch);
    mhit(0).setKey_bar(fHitsPmt[h].key_bar);
#if DIRCTRUTHEXTRA
    hddm_s::DircTruthPmtHitExtraList mhitextra = mhit(0).addDircTruthPmtHitExtras(1);
    mhitextra(0).setT_fixed(fHitsPmt[h].t_fixed_ns);
    mhitextra(0).setPath(fHitsPmt[h].path);
    mhitextra(0).setRefl(fHitsPmt[h].refl);
    mhitextra(0).setBbrefl(fHitsPmt[h].bbrefl);
#endif
  }

  fHitsBar.clear();
  fHitsPmt.clear();
  fHitsWob.clear();
}

int GlueXSensitiveDetectorDIRC::GetIdent(std::string div, 
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

double GlueXSensitiveDetectorDIRC::GetDetectionEfficiency(double energy)
{
   if (fDetEff == 0)
      InitializeDetEff();
   double wavelength = CLHEP::hbarc * CLHEP::twopi / energy;
   return fDetEff->Eval(wavelength / nm);
}

void GlueXSensitiveDetectorDIRC::InitializeDetEff()
{
   G4AutoLock barrier(&fMutex);
   if (fDetEff != 0)
      return;

   // quantum efficiency for H12700
   // defined for wavelength in the range [0, 1000] nm
   double fEfficiency[1000] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.196823, 0.202016, 0.206993, 0.211762, 0.216329,
                               0.220700, 0.224884, 0.228885, 0.232711, 0.236366,
                               0.239858, 0.243192, 0.246373, 0.249408, 0.252301,
                               0.255058, 0.257683, 0.260182, 0.262560, 0.264822,
                               0.266971, 0.269013, 0.270951, 0.272791, 0.274535,
                               0.276189, 0.277756, 0.279239, 0.280643, 0.281971,
                               0.283227, 0.284413, 0.285534, 0.286592, 0.287590,
                               0.288531, 0.289419, 0.290255, 0.291043, 0.291785,
                               0.292484, 0.293142, 0.293762, 0.294345, 0.294894,
                               0.295412, 0.295899, 0.296358, 0.296791, 0.297200,
                               0.297585, 0.297950, 0.298296, 0.298624, 0.298935,
                               0.299231, 0.299514, 0.299784, 0.300043, 0.300291,
                               0.300531, 0.300763, 0.300988, 0.301206, 0.301420,
                               0.301629, 0.301835, 0.302038, 0.302238, 0.302437,
                               0.302636, 0.302833, 0.303031, 0.303230, 0.303429,
                               0.303630, 0.303832, 0.304037, 0.304244, 0.304453,
                               0.304666, 0.304881, 0.305099, 0.305321, 0.305545,
                               0.305774, 0.306005, 0.306240, 0.306478, 0.306720,
                               0.306964, 0.307212, 0.307463, 0.307717, 0.307973,
                               0.308232, 0.308493, 0.308756, 0.309022, 0.309288,
                               0.309557, 0.309826, 0.310096, 0.310367, 0.310637,
                               0.310908, 0.311178, 0.311447, 0.311715, 0.311981,
                               0.312245, 0.312507, 0.312766, 0.313022, 0.313274,
                               0.313522, 0.313766, 0.314005, 0.314238, 0.314466,
                               0.314688, 0.314903, 0.315111, 0.315311, 0.315504,
                               0.315688, 0.315863, 0.316029, 0.316185, 0.316331,
                               0.316466, 0.316590, 0.316703, 0.316804, 0.316893,
                               0.316968, 0.317031, 0.317080, 0.317114, 0.317135,
                               0.317140, 0.317130, 0.317105, 0.317063, 0.317005,
                               0.316931, 0.316839, 0.316729, 0.316602, 0.316457,
                               0.316292, 0.316109, 0.315907, 0.315686, 0.315444,
                               0.315182, 0.314900, 0.314598, 0.314274, 0.313929,
                               0.313563, 0.313175, 0.312765, 0.312333, 0.311878,
                               0.311402, 0.310902, 0.310379, 0.309834, 0.309265,
                               0.308673, 0.308057, 0.307418, 0.306755, 0.306068,
                               0.305357, 0.304622, 0.303863, 0.303080, 0.302273,
                               0.301442, 0.300586, 0.299706, 0.298802, 0.297874,
                               0.296921, 0.295944, 0.294943, 0.293918, 0.292869,
                               0.291796, 0.290699, 0.289579, 0.288434, 0.287266,
                               0.286075, 0.284860, 0.283623, 0.282362, 0.281078,
                               0.279772, 0.278443, 0.277093, 0.275720, 0.274325,
                               0.272908, 0.271470, 0.270012, 0.268532, 0.267031,
                               0.265510, 0.263970, 0.262409, 0.260829, 0.259229,
                               0.257611, 0.255974, 0.254319, 0.252646, 0.250956,
                               0.249248, 0.247524, 0.245782, 0.244025, 0.242252,
                               0.240464, 0.238661, 0.236843, 0.235011, 0.233165,
                               0.231306, 0.229433, 0.227549, 0.225652, 0.223743,
                               0.221823, 0.219892, 0.217951, 0.216000, 0.214040,
                               0.212070, 0.210092, 0.208105, 0.206111, 0.204110,
                               0.202102, 0.200087, 0.198067, 0.196041, 0.194011,
                               0.191976, 0.189937, 0.187894, 0.185849, 0.183801,
                               0.181750, 0.179699, 0.177646, 0.175592, 0.173538,
                               0.171484, 0.169431, 0.167380, 0.165330, 0.163282,
                               0.161236, 0.159194, 0.157155, 0.155120, 0.153089,
                               0.151063, 0.149042, 0.147027, 0.145018, 0.143016,
                               0.141020, 0.139032, 0.137052, 0.135080, 0.133116,
                               0.131162, 0.129217, 0.127281, 0.125356, 0.123442,
                               0.121538, 0.119645, 0.117765, 0.115896, 0.114039,
                               0.112195, 0.110364, 0.108547, 0.106743, 0.104953,
                               0.103177, 0.101415, 0.099669, 0.097937, 0.096221,
                               0.094521, 0.092836, 0.091168, 0.089515, 0.087880,
                               0.086261, 0.084659, 0.083074, 0.081506, 0.079956,
                               0.078424, 0.076910, 0.075413, 0.073935, 0.072475,
                               0.071034, 0.069611, 0.068206, 0.066821, 0.065454,
                               0.064106, 0.062776, 0.061466, 0.060175, 0.058903,
                               0.057650, 0.056416, 0.055201, 0.054005, 0.052829,
                               0.051671, 0.050532, 0.049413, 0.048312, 0.047230,
                               0.046166, 0.045122, 0.044096, 0.043088, 0.042099,
                               0.041128, 0.040175, 0.039241, 0.038324, 0.037424,
                               0.036542, 0.035678, 0.034831, 0.034000, 0.033187,
                               0.032389, 0.031609, 0.030844, 0.030096, 0.029363,
                               0.028645, 0.027943, 0.027255, 0.026583, 0.025924,
                               0.025280, 0.024650, 0.024033, 0.023430, 0.022840,
                               0.022262, 0.021697, 0.021144, 0.020603, 0.020073,
                               0.019555, 0.019048, 0.018551, 0.018065, 0.017588,
                               0.017122, 0.016665, 0.016217, 0.015777, 0.015346,
                               0.014924, 0.014509, 0.014102, 0.013702, 0.013310,
                               0.012924, 0.012544, 0.012171, 0.011803, 0.011442,
                               0.011085, 0.010734, 0.010388, 0.010047, 0.009710,
                               0.009377, 0.009048, 0.008724, 0.008403, 0.008085,
                               0.007772, 0.007461, 0.007154, 0.006850, 0.006548,
                               0.006250, 0.005955, 0.005663, 0.005374, 0.005087,
                               0.004804, 0.004523, 0.004246, 0.003972, 0.003701,
                               0.003434, 0.003171, 0.002911, 0.002656, 0.002405,
                               0.002159, 0.001918, 0.001682, 0.001453, 0.001229,
                               0.001013, 0.000803, 0.000602, 0.000409, 0.000225,
                               0.000051, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                               0.000000, 0.000000, 0.000000, 0.000000, 0.000000
                              };
   double fLambda[1000];
   for (Int_t i=0; i < 1000; i++) {
      fLambda[i] = i;
   }
   fDetEff = new TGraph(1000, fLambda, fEfficiency);
}

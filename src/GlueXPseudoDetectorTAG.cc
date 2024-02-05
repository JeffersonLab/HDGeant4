//
// class implementation for GlueXPseudoDetectorTAG
//
// author: richard.t.jones at uconn.edu
// version: october 20, 2016

#include "GlueXPseudoDetectorTAG.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXPhotonBeamGenerator.hh"

#include "G4EventManager.hh"
#include "G4SystemOfUnits.hh"

#include <JANA/JApplication.h>

// Place limits on the number of hits per counter or column per event
int GlueXPseudoDetectorTAG::HODO_MAX_HITS = 5000;
int GlueXPseudoDetectorTAG::MICRO_MAX_HITS = 5000;

double GlueXPseudoDetectorTAG::HODO_TWO_HIT_TIME_RESOL = 25.*ns;
double GlueXPseudoDetectorTAG::MICRO_TWO_HIT_TIME_RESOL = 25.*ns;
double GlueXPseudoDetectorTAG::MICRO_HIT_DE = 3.5*MeV;
double GlueXPseudoDetectorTAG::HODO_HIT_DE = 0.55*MeV;

double GlueXPseudoDetectorTAG::MICRO_LIMITS_ERANGE[2];
int GlueXPseudoDetectorTAG::MICRO_CHANNEL_NUMBER[GlueXPseudoDetectorTAG::MICRO_NCHANNELS];
double GlueXPseudoDetectorTAG::MICRO_CHANNEL_EMIN[GlueXPseudoDetectorTAG::MICRO_NCHANNELS];
double GlueXPseudoDetectorTAG::MICRO_CHANNEL_EMAX[GlueXPseudoDetectorTAG::MICRO_NCHANNELS];

double GlueXPseudoDetectorTAG::HODO_LIMITS_ERANGE[2];
int GlueXPseudoDetectorTAG::HODO_CHANNEL_NUMBER[GlueXPseudoDetectorTAG::HODO_NCHANNELS];
double GlueXPseudoDetectorTAG::HODO_CHANNEL_EMIN[GlueXPseudoDetectorTAG::HODO_NCHANNELS];
double GlueXPseudoDetectorTAG::HODO_CHANNEL_EMAX[GlueXPseudoDetectorTAG::HODO_NCHANNELS];

double GlueXPseudoDetectorTAG::TAGGER_TMIN_NS = -200.*ns;
double GlueXPseudoDetectorTAG::TAGGER_TMAX_NS = +200.*ns;

int GlueXPseudoDetectorTAG::instanceCount = 0;
G4Mutex GlueXPseudoDetectorTAG::fMutex = G4MUTEX_INITIALIZER;

GlueXPseudoDetectorTAG::GlueXPseudoDetectorTAG(int runNo)
 : fRunNo(runNo)
{
   if (runNo != 0)
      setRunNo(runNo);
   G4AutoLock barrier(&fMutex);
   instanceCount++;
}

GlueXPseudoDetectorTAG::~GlueXPseudoDetectorTAG()
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

GlueXPseudoDetectorTAG::GlueXPseudoDetectorTAG(GlueXPseudoDetectorTAG &src)
{}

GlueXPseudoDetectorTAG& GlueXPseudoDetectorTAG::operator=(GlueXPseudoDetectorTAG &src)
{
   return *this;
}

inline void GlueXPseudoDetectorTAG::setRunNo(int runno)
{
   fRunNo = runno;
   G4AutoLock barrier(&fMutex);
   extern jana::JApplication *japp;
   if (japp == 0) {
      G4cerr << "Error in GlueXPseudoDetectorTAG constructor - "
             << "jana global DApplication object not set, "
             << "cannot continue." << G4endl;
      exit(-1);
   }
   jana::JCalibration *jcalib = japp->GetJCalibration(runno);
   std::map<string, float> beam_parms;
   jcalib->Get("PHOTON_BEAM/endpoint_energy", beam_parms);
   double endpoint_energy = beam_parms.at("PHOTON_BEAM_ENDPOINT_ENERGY")*GeV;
   std::map<string, float> beam_calib;
   jcalib->Get("PHOTON_BEAM/hodoscope/endpoint_calib", beam_calib);
   double endpoint_calib = endpoint_energy;
   if (beam_calib.find("TAGGER_CALIB_ENERGY") != beam_calib.end()) {
      endpoint_calib = beam_calib.at("TAGGER_CALIB_ENERGY")*GeV;
   }
   std::vector<std::map<string, float> > micro_parms;
   jcalib->Get("PHOTON_BEAM/microscope/scaled_energy_range", micro_parms);
   double Emin = 1e9;
   double Emax = -1e9;
   for (unsigned int i=0; i < micro_parms.size(); ++i) {
      MICRO_CHANNEL_NUMBER[i] = micro_parms[i]["column"];
      MICRO_CHANNEL_EMIN[i] = micro_parms[i]["xlow"] * endpoint_calib
                              + endpoint_energy - endpoint_calib;
      MICRO_CHANNEL_EMAX[i] = micro_parms[i]["xhigh"] * endpoint_calib
                              + endpoint_energy - endpoint_calib;
      Emin = (MICRO_CHANNEL_EMIN[i] < Emin)? MICRO_CHANNEL_EMIN[i] : Emin;
      Emax = (MICRO_CHANNEL_EMAX[i] > Emax)? MICRO_CHANNEL_EMAX[i] : Emax;
   }
   MICRO_LIMITS_ERANGE[0] = Emax;
   MICRO_LIMITS_ERANGE[1] = Emin;
   std::vector<std::map<string, float> > hodo_parms;
   jcalib->Get("PHOTON_BEAM/hodoscope/scaled_energy_range", hodo_parms);
   for (unsigned int i=0; i < hodo_parms.size(); ++i) {
      HODO_CHANNEL_NUMBER[i] = hodo_parms[i]["counter"];
      HODO_CHANNEL_EMIN[i] = hodo_parms[i]["xlow"] * endpoint_calib
                              + endpoint_energy - endpoint_calib;
      HODO_CHANNEL_EMAX[i] = hodo_parms[i]["xhigh"] * endpoint_calib
                              + endpoint_energy - endpoint_calib;
      Emin = (HODO_CHANNEL_EMIN[i] < Emin)? HODO_CHANNEL_EMIN[i] : Emin;
      Emax = (HODO_CHANNEL_EMAX[i] > Emax)? HODO_CHANNEL_EMAX[i] : Emax;
   }
   HODO_LIMITS_ERANGE[0] = Emax;
   HODO_LIMITS_ERANGE[1] = Emin;

   G4cout << "TAGGER: all parameters loaded from ccdb" << G4endl;
}

int GlueXPseudoDetectorTAG::addTaggerPhoton(const G4Event *event,
                                            double energy, double time,
                                            int bg) const
{
   // look up which tagger channel is hit, if any

   int micro_channel = -1;
   int hodo_channel = -1;
   double micro_energy = energy;
   double hodo_energy = energy;

   if (energy < MICRO_LIMITS_ERANGE[0] && energy > MICRO_LIMITS_ERANGE[1]) {
      int i = MICRO_NCHANNELS * (energy - MICRO_LIMITS_ERANGE[0]) /
              (MICRO_LIMITS_ERANGE[1] - MICRO_LIMITS_ERANGE[0]);
      while (energy < MICRO_CHANNEL_EMIN[i] || MICRO_CHANNEL_EMIN[i] < 1.) 
         ++i;
      while (energy > MICRO_CHANNEL_EMAX[i] || MICRO_CHANNEL_EMAX[i] < 1.) 
         --i;
      if (energy >= MICRO_CHANNEL_EMIN[i] && energy <= MICRO_CHANNEL_EMAX[i]) {
         micro_energy = (MICRO_CHANNEL_EMIN[i] + MICRO_CHANNEL_EMAX[i]) / 2;
         micro_channel = MICRO_CHANNEL_NUMBER[i];
      }
   }
   if (energy < HODO_LIMITS_ERANGE[0] && energy > HODO_LIMITS_ERANGE[1]) {
      int i = HODO_NCHANNELS * (energy - HODO_LIMITS_ERANGE[0]) /
              (HODO_LIMITS_ERANGE[1] - HODO_LIMITS_ERANGE[0]);
      while (energy < HODO_CHANNEL_EMIN[i] || HODO_CHANNEL_EMIN[i] < 1.) 
         ++i;
      while (energy > HODO_CHANNEL_EMAX[i] || HODO_CHANNEL_EMAX[i] < 1.) 
         --i;
      if (energy >= HODO_CHANNEL_EMIN[i] && energy <= HODO_CHANNEL_EMAX[i]) {
         hodo_energy = (HODO_CHANNEL_EMIN[i] + HODO_CHANNEL_EMAX[i]) / 2;
         hodo_channel = HODO_CHANNEL_NUMBER[i];
      }
   }

   if (micro_channel < 0 && hodo_channel < 0)
      return false;

   double beam_period = GlueXPhotonBeamGenerator::getBeamBucketPeriod();
   double t = floor(time / beam_period + 0.5) * beam_period;

   // pack hit into ouptut hddm record
 
   G4VUserEventInformation* info = event->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXPseudoDetectorTAG::addTaggerPhoton error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getTaggers().size() == 0)
      hitview.addTaggers();
   hddm_s::Tagger &tagger = hitview.getTagger();

   // post microscope hit to the microscope hits tree

   if (micro_channel > -1) {
      hddm_s::MicroChannelList columns = tagger.getMicroChannels();
      hddm_s::MicroChannelList::iterator citer;
      for (citer = columns.begin(); citer != columns.end(); ++citer) {
         if (citer->getColumn() == micro_channel)
            break;
      }
      if (citer == columns.end()) {
         columns = tagger.addMicroChannels();
         citer = columns.begin();
         citer->setColumn(micro_channel);
         citer->setRow(0);
         citer->setE(micro_energy/GeV);
      }
      hddm_s::TaggerTruthHitList hits = citer->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator hiter;
      int hit = 0;
      for (hiter = hits.begin(); hiter != hits.end(); ++hiter, ++hit) {
         if (fabs(t - hiter->getT()*ns) < MICRO_TWO_HIT_TIME_RESOL)
            break;
         else if (hiter->getT()*ns > t) {
            if (citer->getTaggerTruthHits().size() < MICRO_MAX_HITS) {
               hits = citer->addTaggerTruthHits(1, hit);
               hiter = hits.begin();
               hiter->setE(energy/GeV);
               hiter->setT(1e99);
               hiter->setDE(0);
               hiter->setBg(bg);
            }
            else {
               G4cerr << "GlueXPseudoDetectorTAG::addTaggerPhoton warning: "
                      << "TAGM max hit count " << MICRO_MAX_HITS 
                      << " exceeded, discarding hit." << G4endl;
            }
            break;
         }
      }
      if (hiter != hits.end()) {              // merge with a prior hit
         if (t < hiter->getT()*ns) {
            hiter->setT(t/ns);
            hiter->setE(energy/GeV);
            hiter->setBg(bg);
         }
         hiter->setDE(hiter->getDE() + MICRO_HIT_DE/GeV);
      }
      else {                                 // make a new hit
         if (citer->getTaggerTruthHits().size() < MICRO_MAX_HITS) {
            hits = citer->addTaggerTruthHits();
            hiter = hits.begin();
            hiter->setT(t/ns);
            hiter->setE(energy/GeV);
            hiter->setDE(MICRO_HIT_DE/GeV);
            hiter->setBg(bg);
         }
         else {
            G4cerr << "GlueXPseudoDetectorTAG::addTaggerPhoton warning: "
             << "TAGM max hit count " << MICRO_MAX_HITS << " exceeded,"
             << " discarding hit." << G4endl;
         }
      }
   }
 
   // post hodoscope hit to the hodoscope hits tree

   if (hodo_channel > -1) {
      hddm_s::HodoChannelList counters = tagger.getHodoChannels();
      hddm_s::HodoChannelList::iterator citer;
      for (citer = counters.begin(); citer != counters.end(); ++citer) {
         if (citer->getCounterId() == hodo_channel)
            break;
      }
      if (citer == counters.end()) {
         counters = tagger.addHodoChannels();
         citer = counters.begin();
         citer->setCounterId(hodo_channel);
         citer->setE(hodo_energy/GeV);
      }
      hddm_s::TaggerTruthHitList hits = citer->getTaggerTruthHits();
      hddm_s::TaggerTruthHitList::iterator hiter;
      int hit = 0;
      for (hiter = hits.begin(); hiter != hits.end(); ++hiter, ++hit) {
         if (fabs(t - hiter->getT()*ns) < HODO_TWO_HIT_TIME_RESOL)
            break;
         else if (hiter->getT()*ns > t) {
            if (citer->getTaggerTruthHits().size() < HODO_MAX_HITS) {
               hits = citer->addTaggerTruthHits(1, hit);
               hiter = hits.begin();
               hiter->setE(energy/GeV);
               hiter->setT(1e99);
               hiter->setDE(0);
               hiter->setBg(bg);
            }
            else {
               G4cerr << "GlueXPseudoDetectorTAG::addTaggerPhoton warning: "
                << "TAGH max hit count " << HODO_MAX_HITS << " exceeded,"
                << " discarding hit." << G4endl;
            }
            break;
         }
      }
      if (hiter != hits.end()) {              // merge with a prior hit
         if (t < hiter->getT()*ns) {
            hiter->setT(t/ns);
            hiter->setE(energy/GeV);
            hiter->setBg(bg);
         }
         hiter->setDE(hiter->getDE() + HODO_HIT_DE/GeV);
      }
      else {                                 // make a new hit
         if (citer->getTaggerTruthHits().size() < HODO_MAX_HITS) {
            hits = citer->addTaggerTruthHits();
            hiter = hits.begin();
            hiter->setT(t/ns);
            hiter->setE(energy/GeV);
            hiter->setDE(HODO_HIT_DE/GeV);
            hiter->setBg(bg);
         }
         else {
            G4cerr << "GlueXPseudoDetectorTAG::addTaggerPhoton warning: "
             << "TAGH max hit count " << MICRO_MAX_HITS << " exceeded,"
             << " discarding hit." << G4endl;
         }
      }
   }
   return (micro_channel > -1 || hodo_channel > -1);
}

int GlueXPseudoDetectorTAG::addRFsync(const G4Event *event, double tsync) const
{
   G4VUserEventInformation* info = event->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXPseudoDetectorTAG::addRFsync error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getRFtimes().size() == 0)
      hitview.addRFtimes();
   hddm_s::RFtime &rftime = hitview.getRFtime();
   rftime.setTsync(tsync);
   return 1;
}

void GlueXPseudoDetectorTAG::Draw() const
{
   // not yet implemented
}

void GlueXPseudoDetectorTAG::Print() const
{
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   if (info != 0) {
      hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
      if (record != 0) {
         G4cout << "GlueXPseudoDetectorTAG: " << G4endl;
         hddm_s::MicroChannelList micros = record->getMicroChannels();
         hddm_s::MicroChannelList::iterator miter;
         for (miter = micros.begin(); miter != micros.end(); ++miter) {
            G4cout << "  microscope column " << miter->getColumn()
                   << ", row " << miter->getRow()
                   << ", E " << miter->getE() << " GeV:"
                   << G4endl;
            hddm_s::TaggerTruthHitList mhits = miter->getTaggerTruthHits();
            hddm_s::TaggerTruthHitList::iterator michiter;
            for (michiter = mhits.begin(); michiter != mhits.end(); ++michiter) {
               G4cout << "   dE = " << michiter->getDE() << " GeV" << G4endl
                      << "   E = " << michiter->getE() << " GeV" << G4endl
                      << "   t = " << michiter->getT() << " ns" << G4endl
                      << "   bg = " << michiter->getBg() << G4endl
                      << G4endl;
            }
         }
         hddm_s::HodoChannelList hodos = record->getHodoChannels();
         hddm_s::HodoChannelList::iterator hiter;
         for (hiter = hodos.begin(); hiter != hodos.end(); ++hiter) {
            G4cout << "  hodoscope counter " << hiter->getCounterId()
                   << ", E " << hiter->getE() << " GeV:"
                   << G4endl;
            hddm_s::TaggerTruthHitList hhits = hiter->getTaggerTruthHits();
            hddm_s::TaggerTruthHitList::iterator hodhiter;
            for (hodhiter = hhits.begin(); hodhiter != hhits.end(); ++hodhiter) {
               G4cout << "   dE = " << hodhiter->getDE() << " GeV" << G4endl
                      << "   E = " << hodhiter->getE() << " GeV" << G4endl
                      << "   t = " << hodhiter->getT() << " ns" << G4endl
                      << "   bg = " << hodhiter->getBg() << G4endl
                      << G4endl;
            }
         }
      }
   }
}

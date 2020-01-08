//
// GlueXPhysicsList class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include "GlueXPhysicsList.hh"

#include <G4ParallelWorldPhysics.hh>

#include <iomanip>   

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmExtraPhysics.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmParameters.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4OpticalProcessIndex.hh"

#if USING_DIRACXX
G4ThreadLocal GlueXBeamConversionProcess *GlueXPhysicsList::fBeamConversion = 0;
#endif

GlueXPhysicsList::GlueXPhysicsList(const GlueXDetectorConstruction *geometry,
                                   G4int verbosity)
 : G4VModularPhysicsList(), fOpticalPhysics(0)
{
   if (geometry == 0) {
      geometry = GlueXDetectorConstruction::GetInstance();
      if (geometry == 0) {
         G4cerr << "GlueXPhysicsList constructor error - "
                << "cannot construct GlueXPhysicsList until the detector "
                << "geometry has been constructed, aborting."
                << G4endl;
         exit(1);
      }
   }

   fOptions = GlueXUserOptions::GetInstance();
   if (fOptions == 0) {
      G4cerr << "Error in GlueXPhysicsList constructor - "
             << "GlueXUserOptions::GetInstance() returns null, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   G4DataQuestionaire it(photon);
   if (verbosity > 0) {
      G4cout << "<<< GlueX Physics List simulation engine: "
             << "modeled after Geant4 QGSP_FTFP_BERT 4.0"
             << G4endl;
      G4cout << G4endl;
   }

   defaultCutValue = 0.7*CLHEP::mm;  
   SetVerboseLevel(verbosity);

   // Parallel world transportation
   for (int para=1; para <= geometry->GetParallelWorldCount(); ++para) {
      G4String name = geometry->GetParallelWorldName(para);
      RegisterPhysics(new G4ParallelWorldPhysics(name, true));
   }

   // EM Physics
   RegisterPhysics( new G4EmStandardPhysics_option1(verbosity) );

   // Synchroton Radiation & GN Physics
   RegisterPhysics( new G4EmExtraPhysics(verbosity) );

   // Decays
   RegisterPhysics( new G4DecayPhysics(verbosity) );

   // Hadron Elastic scattering
   RegisterPhysics( new G4HadronElasticPhysics(verbosity) );

   // Hadron Physics
   RegisterPhysics( new G4HadronPhysicsQGSP_FTFP_BERT(verbosity));

   // Stopping Physics
   RegisterPhysics( new G4StoppingPhysics(verbosity) );

   // Ion Physics
   RegisterPhysics( new G4IonPhysics(verbosity));

   // Neutron tracking cut
   RegisterPhysics( new G4NeutronTrackingCut(verbosity));

   // Optical photons (Cerenkov only)
   int doCerenkov=0;
   std::map<int,int> ckov;
   if (fOptions->Find("CKOV", ckov)) {
      if (ckov.find(1) != ckov.end() && ckov[1] != 0)
         doCerenkov = 1;
   }
   int doAbsorption=0;
   std::map<int,int> labs;
   if (fOptions->Find("LABS", labs)) {
      if (labs.find(1) != labs.end() && labs[1] != 0)
         doAbsorption = 1;
   }
   if (doCerenkov) {
      fOpticalPhysics = new G4OpticalPhysics(verbosity);
      fOpticalPhysics->Configure(kCerenkov, 1);
      fOpticalPhysics->Configure(kScintillation, 0);
      fOpticalPhysics->Configure(kAbsorption, doAbsorption);
      fOpticalPhysics->Configure(kRayleigh, 0);
      fOpticalPhysics->Configure(kMieHG, 0);
      fOpticalPhysics->Configure(kBoundary, 1);
      fOpticalPhysics->Configure(kWLS, 0);
      RegisterPhysics(fOpticalPhysics);
   }
   else {
      fOpticalPhysics = 0;
   }
}

GlueXPhysicsList::~GlueXPhysicsList()
{}

void GlueXPhysicsList::ConstructParticle()
{
   G4VModularPhysicsList::ConstructParticle();
#if VERBOSE_PARTICLES
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4cout << "*** particle type " <<  particleName << G4endl;
   }
#endif
}

void GlueXPhysicsList::ConstructProcess()
{
   // Read special cuts from the user options file
   double tcut = 0;
   double KEcut_gamma = 0;
   double KEcut_electron = 0;
   double KEcut_neutron = 0;
   double KEcut_proton = 0;
   double KEcut_muon = 0;
   std::map<int,double> cut;
   if (fOptions->Find("CUTS", cut)) {
      if (cut.find(1) != cut.end() && cut[1] > 0)
         KEcut_gamma = cut[1]*CLHEP::GeV;
      if (cut.find(2) != cut.end() && cut[2] > 0)
         KEcut_electron = cut[2]*CLHEP::GeV;
      if (cut.find(3) != cut.end() && cut[3] > 0)
         KEcut_neutron = cut[3]*CLHEP::GeV;
      if (cut.find(4) != cut.end() && cut[4] > 0)
         KEcut_proton = cut[4]*CLHEP::GeV;
      if (cut.find(5) != cut.end() && cut[5] > 0)
         KEcut_muon = cut[5]*CLHEP::GeV;
   }
   std::map<int,double> tof;
   if (fOptions->Find("TOFMAX", tof) && tof.find(1) != tof.end()) {
      tcut = tof[1]*CLHEP::s;
   }
   else {
      tcut = 1.0e-5*CLHEP::s;
   }

   // Give cuts guidance to the electromagnetic processes
   G4EmParameters *emparams = G4EmParameters::Instance();
   if (KEcut_gamma > 0)
      emparams->SetBremsstrahlungTh(KEcut_gamma);
   if (KEcut_electron > 0)
      emparams->SetLowestElectronEnergy(KEcut_electron);
   if (KEcut_proton > 0 && KEcut_proton < KEcut_muon)
      emparams->SetLowestMuHadEnergy(KEcut_proton);
   else if (KEcut_muon > 0)
      emparams->SetLowestMuHadEnergy(KEcut_muon);

   // Construct all of the standard physics processes
   G4VModularPhysicsList::ConstructProcess();

#if USING_DIRACXX
   // Add a process for TPOL beam photon pair conversion process
   fBeamConversion = new GlueXBeamConversionProcess("TPolBeamConversion");
#endif

   // create the special cuts processes and register them
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager *mgr = particle->GetProcessManager();
      if (mgr == 0) {
         G4cerr << "Error in GlueXPhysicsList::ConstructProcess - "
                << particleName << " process manager is null, "
                << "cannot continue." << G4endl;
         exit(1);
      }
      else if (particleName == "gamma") {
#if USING_DIRACXX
         mgr->AddDiscreteProcess(fBeamConversion);
#endif
         if (KEcut_gamma > 0) {
            G4UserLimits *glimits = new G4UserLimits();
            glimits->SetUserMaxTime(tcut);
            glimits->SetUserMinEkine(KEcut_gamma);
            GlueXSpecialCuts *gcuts;
            gcuts = new GlueXSpecialCuts("GlueX tracking limits for gammas");
            gcuts->SetUserLimits(glimits);
            mgr->AddProcess(gcuts, -1, -1, 1);
         }
         else
            continue;
      }
      else if (particleName == "e-" || particleName == "e+") {
         if (KEcut_electron > 0) {
            G4UserLimits *elimits = new G4UserLimits();
            elimits->SetUserMaxTime(tcut);
            elimits->SetUserMinEkine(KEcut_electron);
            GlueXSpecialCuts *ecuts;
            ecuts = new GlueXSpecialCuts("GlueX tracking limits for electrons");
            ecuts->SetUserLimits(elimits);
            mgr->AddProcess(ecuts, -1, -1, 1);
         }
         else
            continue;
      }
      else if (particleName == "neutron") {
         if (KEcut_neutron > 0) {
            G4UserLimits *nlimits = new G4UserLimits();
            nlimits->SetUserMaxTime(tcut);
            nlimits->SetUserMinEkine(KEcut_neutron);
            GlueXSpecialCuts *ncuts;
            ncuts = new GlueXSpecialCuts("GlueX tracking limits for neutrons");
            ncuts->SetUserLimits(nlimits);
            mgr->AddProcess(ncuts, -1, -1, 1);
         }
         else
            continue;
      }
      else if (particleName == "proton" ||
               particleName == "pi+" || particleName == "pi-" ||
               particleName == "kaon+" || particleName == "kaon-")
      {
         if (KEcut_proton > 0) {
            G4UserLimits *hlimits = new G4UserLimits();
            hlimits->SetUserMaxTime(tcut);
            hlimits->SetUserMinEkine(KEcut_proton);
            GlueXSpecialCuts *hcuts;
            hcuts = new GlueXSpecialCuts("GlueX tracking limits for protons");
            hcuts->SetUserLimits(hlimits);
            mgr->AddProcess(hcuts, -1, -1, 1);
         }
         else
            continue;
      }
      else if (particleName == "mu-" || particleName == "mu+") {
         if (KEcut_muon > 0) {
            G4UserLimits *mlimits = new G4UserLimits();
            mlimits->SetUserMaxTime(tcut);
            mlimits->SetUserMinEkine(KEcut_muon);
            GlueXSpecialCuts *mcuts;
            mcuts = new GlueXSpecialCuts("GlueX tracking limits for muons");
            mcuts->SetUserLimits(mlimits);
            mgr->AddProcess(mcuts, -1, -1, 1);
         }
         else
            continue;
      }
      else if (particleName == "GenericIon") {
         // mgr->DumpInfo();
      }
      else {
         continue;
      }
      if (verboseLevel > 0) {
         G4cout << "Added user cuts for " << particleName << G4endl;
      }
   }

   // Try to limit the number of secondaries generated by physics processes
   // that will end up falling below the particle thresholds set above. In
   // no case will I ever set this to lower than 1 MeV.
   double lowE = 1.*CLHEP::keV;
   lowE = (lowE > KEcut_gamma && KEcut_gamma > 0)? KEcut_gamma : lowE;
   lowE = (lowE > KEcut_electron && KEcut_electron > 0)? KEcut_electron : lowE;
   lowE = (lowE > KEcut_neutron && KEcut_neutron > 0)? KEcut_neutron : lowE;
   lowE = (lowE > KEcut_proton && KEcut_proton > 0)? KEcut_proton : lowE;
   lowE = (lowE > KEcut_muon && KEcut_muon > 0)? KEcut_muon : lowE;
   G4ProductionCutsTable *cuts_table;
   cuts_table = G4ProductionCutsTable::GetProductionCutsTable();
   double highE = cuts_table->GetHighEdgeEnergy();
   cuts_table->SetEnergyRange(lowE, highE);
}

void GlueXPhysicsList::SetCuts()
{
   if (verboseLevel > 1) {
      G4cout << "GlueXPhysicsList::SetCuts:";
   }

   //  "G4VUserPhysicsList::SetCutsWithDefault" method sets 
   //  the default cut value for all particle types 

   SetCutsWithDefault();   

   // Update cuts based on user options
   std::map<int,double> cuts;
   if (fOptions->Find("RANGECUTS", cuts)) {
      if (cuts.find(1) != cuts.end() && cuts[1] > 0) {
         SetCutValue(cuts[1]*CLHEP::cm, "gamma");
         SetApplyCuts(1, "gamma");
      }
      if (cuts.find(2) != cuts.end() && cuts[2] > 0) {
         SetCutValue(cuts[2]*CLHEP::cm, "e-");
         SetCutValue(cuts[2]*CLHEP::cm, "e+");
         SetApplyCuts(1, "e-");
         SetApplyCuts(1, "e+");
      }
      // SetCuts for neutron not currently supported by G4
      // if (cuts.find(3) != cuts.end() && cuts[3] > 0) {
      //   SetCutValue(cuts[3]*CLHEP::cm, "neutron");
      // }
      if (cuts.find(4) != cuts.end() && cuts[4] > 0) {
         SetCutValue(cuts[4]*CLHEP::cm, "proton");
         SetApplyCuts(1, "proton");
      }
      // SetCuts for muon not currently supported by G4
      // if (cuts.find(5) != cuts.end() && cuts[5] > 0) {
      //   SetCutValue(cuts[5]*CLHEP::cm, "mu-");
      //   SetCutValue(cuts[5]*CLHEP::cm, "mu+");
      // }
   }
   std::map<int,int> flags;
   if (fOptions->Find("DRAY", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         SetCutValue(1e3*CLHEP::m, "e-");
         SetApplyCuts(1, "e-");
      }
   }

   if (verboseLevel > 0)
      G4VUserPhysicsList::DumpCutValuesTable();  
}

void GlueXPhysicsList::ListActiveProcesses()
{
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4cout << particleName << ": ApplyCuts is " 
             << ((particle->GetApplyCutsFlag())? "on" : "off")
             << G4endl;
      G4ProcessManager *pman = particle->GetProcessManager();
      G4ProcessVector *pvec = pman->GetProcessList();
      for (int nproc=0; nproc < pvec->size(); ++nproc) {
         G4cout << "     " << (*pvec)[nproc]->GetProcessName() << " : "
                << ((pman->GetProcessActivation(nproc) == 0)? "no" : "yes")
                << " ( ";
         if ((*pvec)[nproc]->isAtRestDoItIsEnabled())
            G4cout << " AtRest ";
         if ((*pvec)[nproc]->isAlongStepDoItIsEnabled())
            G4cout << " AlongStep ";
         if ((*pvec)[nproc]->isPostStepDoItIsEnabled())
            G4cout << " PostStep ";
         G4cout << ")" << G4endl;
      }
   }
}

void GlueXPhysicsList::DoMultipleScattering(G4int flag)
{
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager *pman = particle->GetProcessManager();
      G4ProcessVector *pvec = pman->GetProcessList();
      for (int nproc=0; nproc < pvec->size(); ++nproc) {
         if ((*pvec)[nproc]->GetProcessName() == "msc") {
            pman->SetProcessActivation(nproc, flag);
         }
      }
   }
}

void GlueXPhysicsList::DoBremsstrahlung(G4int flag)
{
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager *pman = particle->GetProcessManager();
      G4ProcessVector *pvec = pman->GetProcessList();
      for (int nproc=0; nproc < pvec->size(); ++nproc) {
         if ((*pvec)[nproc]->GetProcessName().contains("Brems")) {
            pman->SetProcessActivation(nproc, flag);
         }
      }
   }
}

void GlueXPhysicsList::DoComptonScattering(G4int flag)
{
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager *pman = particle->GetProcessManager();
      G4ProcessVector *pvec = pman->GetProcessList();
      for (int nproc=0; nproc < pvec->size(); ++nproc) {
         if ((*pvec)[nproc]->GetProcessName() == "compt") {
            pman->SetProcessActivation(nproc, flag);
         }
      }
   }
}

void GlueXPhysicsList::DoIonizationEnergyLoss(G4int flag)
{
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager *pman = particle->GetProcessManager();
      G4ProcessVector *pvec = pman->GetProcessList();
      for (int nproc=0; nproc < pvec->size(); ++nproc) {
         if ((*pvec)[nproc]->GetProcessName().contains("Ioni")) {
            pman->SetProcessActivation(nproc, flag);
         }
      }
   }
}

void GlueXPhysicsList::DoPairConversion(G4int flag)
{
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager *pman = particle->GetProcessManager();
      G4ProcessVector *pvec = pman->GetProcessList();
      for (int nproc=0; nproc < pvec->size(); ++nproc) {
         if ((*pvec)[nproc]->GetProcessName() == "conv") {
            pman->SetProcessActivation(nproc, flag);
         }
      }
   }
}

void GlueXPhysicsList::DoParticleDecay(G4int flag)
{
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager *pman = particle->GetProcessManager();
      G4ProcessVector *pvec = pman->GetProcessList();
      for (int nproc=0; nproc < pvec->size(); ++nproc) {
         if ((*pvec)[nproc]->GetProcessName() == "Decay") {
            pman->SetProcessActivation(nproc, flag);
         }
      }
   }
}

void GlueXPhysicsList::DoHadronicInteractions(G4int flag)
{
   GetParticleIterator()->reset();
   while ( (*GetParticleIterator())() ) {
      G4ParticleDefinition* particle = GetParticleIterator()->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager *pman = particle->GetProcessManager();
      G4ProcessVector *pvec = pman->GetProcessList();
      for (int nproc=0; nproc < pvec->size(); ++nproc) {
         if ((*pvec)[nproc]->GetProcessName().contains("Inelastic") ||
             (*pvec)[nproc]->GetProcessName().contains("hadElastic") ||
             (*pvec)[nproc]->GetProcessName().contains("CaptureAtRest") ||
             (*pvec)[nproc]->GetProcessName().contains("Nuclear"))
         {
            pman->SetProcessActivation(nproc, flag);
         }
      }
   }
}

void GlueXPhysicsList::DoDeltaRayProduction(G4int flag)
{
}

void GlueXPhysicsList::DoCerenkovRadiation(G4int flag)
{
   if (fOpticalPhysics && flag == 0) {
      G4cerr << "GlueXPhysicsList::DoCerenkovRadiation error - "
             << G4endl
             << "Cerenkov radiation was enabled at pre-init, "
             << "cannot be disabled at run time."
             << G4endl;
   }
   else if (fOpticalPhysics == 0 && flag == 1) {
      G4cerr << "GlueXPhysicsList::DoCerenkovRadiation error - "
             << G4endl
             << "Cerenkov radiation cannot be enabled at run time, "
             << G4endl
             << "please make sure this happens in the GlueXPhysics "
             << "constructor." << G4endl;
   }
}

void GlueXPhysicsList::DoOpticalAbsorption(G4int flag)
{}

void GlueXPhysicsList::SelectActiveProcesses(G4int verbosity)
{
   std::map<int,int> flags;
   if (fOptions->Find("MULS", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         DoMultipleScattering(0);
         if (verbosity > 0) {
            G4cout << "*** Multiple scattering disabled for all particles."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("BREM", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         DoBremsstrahlung(0);
         if (verbosity > 0) {
            G4cout << "*** Bremsstrahlung disabled for all particles."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("COMP", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         DoComptonScattering(0);
         if (verbosity > 0) {
            G4cout << "*** Compton scattering disabled for gammas."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("LOSS", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         DoIonizationEnergyLoss(0);
         if (verbosity > 0) {
            G4cout << "*** Ionization energy loss disabled for all particles."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("PAIR", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         DoPairConversion(0);
         if (verbosity > 0) {
            G4cout << "*** Pair conversion disabled for gammas."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("DCAY", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         DoParticleDecay(0);
         if (verbosity > 0) {
            G4cout << "*** Decays disabled for all particles."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("DRAY", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         DoDeltaRayProduction(0);
         if (verbosity > 0) {
            G4cout << "*** Delta ray production disabled for all particles."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("HADR", flags)) {
      if (flags.find(1) != flags.end() && flags[1] == 0) {
         DoHadronicInteractions(0);
         if (verbosity > 0) {
            G4cout << "*** Hadronic interactions (excludes decays) "
                   << "disabled for all particles."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("CKOV", flags)) {
      if (flags.find(1) == flags.end() || flags[1] != 0) {
         DoCerenkovRadiation(1);
         if (verbosity > 0) {
            G4cout << "*** Cerenkov radiation enabled for charged particles."
                   << G4endl;
         }
      }
   }
   if (fOptions->Find("LABS", flags)) {
      if (flags.find(1) == flags.end() || flags[1] != 0) {
         DoOpticalAbsorption(1);
         if (verbosity > 0) {
            G4cout << "*** Light absorption enabled for optical photons."
                   << G4endl;
         }
      }
   }
}

void GlueXPhysicsList::CheckProcessOrdering()
{
   std::cout << "Complete list of PostStepDoIt processes for G4OpticalPhoton,"
             << " in order of execution:" << std::endl;
   G4ProcessManager *mgr = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
   G4ProcessVector *procs = mgr->GetPostStepProcessVector(typeDoIt);
   for (int j=0; j < procs->size(); ++j) {
      std::cout << j << ": " << (*procs)[j]->GetProcessName() << std::endl;
   }
}

void GlueXPhysicsList::DoProcessReordering()
{
   // Special process ordering needed for optical photons
   G4ProcessManager *mgr = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
   G4ProcessVector *procs = mgr->GetPostStepProcessVector(typeDoIt);
   G4VProcess *paraWorld=0;
   for (int j=0; j < procs->size(); ++j) {
      std::string procname((*procs)[j]->GetProcessName());
      if (procname.substr(0, 13) == "ParallelWorld") {
         paraWorld = (*procs)[j];
      }
   }
   if (paraWorld == 0) {
      G4cerr << "ParallelWorld process not found, cannot continue!"
             << G4endl;
      exit(1);
   }
   mgr->SetProcessOrderingToSecond(paraWorld, idxPostStep);
}

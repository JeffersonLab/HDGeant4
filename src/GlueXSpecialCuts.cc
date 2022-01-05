//
// GlueXSpecialCuts class implementation
//
// author: richard.t.jones at uconn.edu
// version: january 22, 2017

#include "GlueXSpecialCuts.hh"
#include "G4TransportationProcessType.hh"

#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4LossTableManager.hh"

G4UserLimits *GlueXSpecialCuts::fUserLimits = 0;
G4Mutex GlueXSpecialCuts::fMutex = G4MUTEX_INITIALIZER;

GlueXSpecialCuts::GlueXSpecialCuts(const G4String& aName)
  : G4VProcess(aName, fUserDefined)
{
   // set Process Sub Type
   SetProcessSubType(static_cast<int>(USER_SPECIAL_CUTS));

   if (verboseLevel > 0) {
      G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

GlueXSpecialCuts::~GlueXSpecialCuts()
{}

GlueXSpecialCuts::GlueXSpecialCuts(GlueXSpecialCuts& right)
  : G4VProcess(right)
{}

G4double GlueXSpecialCuts::
PostStepGetPhysicalInteractionLength( const G4Track& aTrack,
                                            G4double previousStepSize,
                                            G4ForceCondition* condition  )
{
   // condition is set to "Not Forced"
   *condition = NotForced;

   G4double ProposedStep = DBL_MAX;
   if (fUserLimits) {

      // check max kinetic energy first
      G4double Ekine = aTrack.GetKineticEnergy();
      if (Ekine <= fUserLimits->GetUserMinEkine(aTrack))
         return 0.;

      // max track length
      ProposedStep = fUserLimits->GetUserMaxTrackLength(aTrack) -
                     aTrack.GetTrackLength();
      if (ProposedStep < 0.)
         return 0.;

      // max time limit
      G4double tlimit = fUserLimits->GetUserMaxTime(aTrack);
      if (tlimit < DBL_MAX) {
         G4double beta = aTrack.GetDynamicParticle()->GetTotalMomentum() /
	                     aTrack.GetTotalEnergy();
         G4double dTime = tlimit - aTrack.GetGlobalTime();
         G4double temp  = beta * c_light * dTime;
         if (temp < 0.) 
            return 0.;
         if (ProposedStep > temp)
            ProposedStep = temp;
      }

      // min remaining range 
      // (only for charged particle except for chargedGeantino)
      G4double Rmin = fUserLimits->GetUserMinRange(aTrack);
      if (Rmin > DBL_MIN) {
         G4ParticleDefinition* Particle = aTrack.GetDefinition();
         if (Particle->GetPDGCharge() != 0 && Particle->GetPDGMass() > 0.0) {
	        const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
	        G4double RangeNow = theLossTableManager
                                ->GetRange(Particle,Ekine,couple);
            G4double temp = RangeNow - Rmin;
            if (temp < 0.) 
               return 0.;
            if (ProposedStep > temp)
               ProposedStep = temp;
         }	 
      }
   }
   return ProposedStep;
}

G4VParticleChange*
GlueXSpecialCuts::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&  )
{

   // Stop (but do not kill) the current particle, if requested by G4UserLimits

   if (verboseLevel > 1) {
      G4cout << "    " << aTrack.GetParticleDefinition()->GetParticleName()
             << " " << aTrack.GetTrackID() << " stopped by GlueXSpecialCuts: "
             << "t = " << aTrack.GetGlobalTime() << " ns, "
             << "Ekin = " << aTrack.GetKineticEnergy() << " MeV, "
             << "total track length = " << aTrack.GetTrackLength() << " mm."
             << std::endl;
   }

   aParticleChange.Initialize(aTrack);
   aParticleChange.ProposeEnergy(0.) ;
   aParticleChange.ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy()) ;
   //aParticleChange.ProposeTrackStatus(fStopAndKill);
   return &aParticleChange;
}

const G4UserLimits *GlueXSpecialCuts::GetUserLimits() const
{
   return fUserLimits;
}

void GlueXSpecialCuts::SetUserLimits(const G4UserLimits *userlimits)
{
   G4AutoLock barrier(&fMutex);
   if (fUserLimits)
      delete fUserLimits;
   fUserLimits = new G4UserLimits(*userlimits);
}

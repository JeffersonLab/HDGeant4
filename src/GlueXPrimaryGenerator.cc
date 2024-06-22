//
// class implementation for GlueXPrimaryGenerator
//
// author: richard.t.jones at uconn.edu
// version: december 24, 2016
//
// This is the principal way that GlueX events get injected into the
// HDGeant4 simulation. An external event generator writes MC events
// into a hddm file, which is read below and translated into Geant4
// primary vertex objects to be tracked.
//

#include "GlueXPrimaryGenerator.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "HddmOutput.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"

#include <HDDM/hddm_s.hpp>

GlueXPrimaryGenerator::GlueXPrimaryGenerator(hddm_s::istream *hddm_source)
 : fHDDMistream(hddm_source)
{}

GlueXPrimaryGenerator::~GlueXPrimaryGenerator()
{}

void GlueXPrimaryGenerator::GeneratePrimaryVertex(G4Event *event)
{
   hddm_s::HDDM *hddmevent = new hddm_s::HDDM;
   while (hddmevent->getPhysicsEvents().size() == 0) {
      if (! (*fHDDMistream >> *hddmevent)) {
         event->SetEventAborted();
         G4cout << "End of file on hddm input, ending the run here." << std::endl;
         G4RunManager::GetRunManager()->AbortRun();
         return;
      }
   }

   // Override the run number on the input record
   hddm_s::PhysicsEventList pe = hddmevent->getPhysicsEvents();
   HddmOutput::setRunNo(pe(0).getRunNo());
   pe(0).setRunNo(HddmOutput::getRunNo());

   // Store generated event info so it can be written to output file
   GlueXUserEventInformation *event_info;
   event_info = new GlueXUserEventInformation(hddmevent);
   event->SetUserInformation(event_info);

   // Unpack generated event and prepare initial state for simulation
   int Nprimaries = 0;
   hddm_s::VertexList vertices = hddmevent->getVertices();
   if (vertices.size() == 0) {
      G4cout << "No vertices in input HDDM event!" << G4endl;
      event->SetEventAborted();
      return;
   }
   hddm_s::VertexList::iterator it_vertex;
   it_vertex = vertices.begin();
   event->SetEventID(it_vertex->getEventNo());
   G4ThreeVector vtx(GetParticlePosition());
   double tvtx(GetParticleTime());
   hddm_s::Origin &origin = it_vertex->getOrigin();
   double x = origin.getVx() * cm;
   double y = origin.getVy() * cm;
   double z = origin.getVz() * cm;
   double t = origin.getT() * ns;
   if (x == 0 && y == 0 && z == 0) {
      tvtx = (t == 0)? tvtx : 0;
   }
   else {
      double beamVelocity = GlueXPhotonBeamGenerator::getBeamVelocity();
      if (t == 0)
         tvtx += (z - vtx[2]) / beamVelocity;
      else
         tvtx = 0;
      vtx[0] = 0;
      vtx[1] = 0;
      vtx[2] = 0;
   }
   for (it_vertex = vertices.begin();
        it_vertex != vertices.end(); ++it_vertex)
   {
      hddm_s::Origin &origin = it_vertex->getOrigin();
      double x = origin.getVx() * cm;
      double y = origin.getVy() * cm;
      double z = origin.getVz() * cm;
      double t = origin.getT() * ns;
      x += vtx[0];
      y += vtx[1];
      z += vtx[2];
      t += tvtx;
      origin.setVx(x/cm);
      origin.setVy(y/cm);
      origin.setVz(z/cm);
      origin.setT(t/ns);
      G4ThreeVector pos(x, y, z);
      G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, t);
      hddm_s::ProductList &products = it_vertex->getProducts();
      hddm_s::ProductList::iterator it_product;
      for (it_product = products.begin();
           it_product != products.end(); ++it_product)
      {
         // ignore intermediaries in the MC record
         if (it_product->getType() <= 0)
           continue;

         int trackId = it_product->getId();
         int g3type = it_product->getType();
         int pdgtype = it_product->getPdgtype();
         G4ParticleDefinition *part;
         if (pdgtype > 0) {
            part = GlueXPrimaryGeneratorAction::GetParticle(pdgtype);
         }
         else if (g3type > 0) {
            pdgtype = GlueXPrimaryGeneratorAction::ConvertGeant3ToPdg(g3type);
            part = GlueXPrimaryGeneratorAction::GetParticle(pdgtype);
#if FORCE_PARTICLE_TYPE_CHARGED_GEANTINO
            part = GlueXPrimaryGeneratorAction::GetParticle("chargedgeantino");
#endif
         }
         else {
            G4cerr << "=== WARNING in GlueXPrimaryGenerator::"
                      "GeneratePrimaryVertex ==="
                   << G4endl
                   << "   Unknown particle found in input MC record, "
                   << "geant3 type " << g3type 
                   << ", PDG type " << pdgtype
                   << ", failing over to geantino!"
                   << G4endl;
            part = GlueXPrimaryGeneratorAction::GetParticle("geantino");
         }
         hddm_s::Momentum &momentum = it_product->getMomentum();
         double px = momentum.getPx() * GeV;
         double py = momentum.getPy() * GeV;
         double pz = momentum.getPz() * GeV;
         double Etot = momentum.getE() * GeV;
         G4PrimaryParticle *pp = new G4PrimaryParticle(part, px, py, pz, Etot);
         hddm_s::Polarization &polar = it_product->getPolarization();
         pp->SetPolarization(polar.getPx(), polar.getPy(), polar.getPz());
         vertex->SetPrimary(pp);
         event_info->SetGlueXTrackID(++Nprimaries, trackId);
         double mass = part->GetPDGMass();
         double width = part->GetPDGWidth();
         double minv = sqrt(Etot*Etot - px*px - py*py - pz*pz);
         if (fabs(minv - mass) > width * 1e3 + 3*MeV) {
            G4cerr << "=== WARNING in GlueXPrimaryGenerator::"
                      "GeneratePrimaryVertex ==="
                   << G4endl
                   << "   " << part->GetParticleName()
                   << " found in input MC record, "
                   << "geant3 type " << g3type 
                   << ", PDG type " << pdgtype
                   << " has unphysical mass: "
                   << G4endl
                   << "   expected " << G4BestUnit(mass, "Energy")
                   << ", found " << G4BestUnit(minv, "Energy")
                   << ", difference "
                   << G4BestUnit(mass - minv, "Energy")
                   << G4endl;
         }
      }
      event->AddPrimaryVertex(vertex);
   }
}

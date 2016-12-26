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
#include "G4SystemOfUnits.hh"

#include <HDDM/hddm_s.hpp>

GlueXPrimaryGenerator::GlueXPrimaryGenerator(hddm_s::istream *hddm_source)
 : fHDDMistream(hddm_source)
{}

GlueXPrimaryGenerator::~GlueXPrimaryGenerator()
{}

void GlueXPrimaryGenerator::GeneratePrimaryVertex(G4Event *event)
{
   hddm_s::HDDM *hddmevent = new hddm_s::HDDM;
   try {
      *fHDDMistream >> *hddmevent;
   }
   catch(std::exception e) {
      G4cout << e.what() << G4endl;
      event->SetEventAborted();
      return;
   }

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
   for (it_vertex = vertices.begin();
        it_vertex != vertices.end(); ++it_vertex)
   {
      event->SetEventID(it_vertex->getEventNo());
      hddm_s::Origin &origin = it_vertex->getOrigin();
      double x = origin.getVx() * cm;
      double y = origin.getVy() * cm;
      double z = origin.getVz() * cm;
      double t = origin.getT() * ns;
      if (x == 0 && y == 0 && z == 0) {
         G4ThreeVector vtx(GetParticlePosition());
         x = vtx[0];
         y = vtx[1];
         z = vtx[2];
         origin.setVx(x/cm);
         origin.setVy(y/cm);
         origin.setVz(z/cm);
      }
      if (t == 0) {
         t = GetParticleTime();
         origin.setT(t/ns);
      }
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

         int g3type = it_product->getType();
         int pdgtype = it_product->getPdgtype();
         G4ParticleDefinition *part;
         if (pdgtype > 0 && pdgtype < 999999) {
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
            G4cerr << "Unknown particle found in input MC record, "
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
         vertex->SetPrimary(new G4PrimaryParticle(part, px, py, pz, Etot));
         ++Nprimaries;
      }
      event->AddPrimaryVertex(vertex);
   }
}

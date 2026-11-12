//
// GlueXHitSTCpoint - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016

#include "GlueXHitSTCpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitSTCpoint>* GlueXHitSTCpointAllocator = 0;

GlueXHitSTCpoint::GlueXHitSTCpoint(const GlueXHitSTCpoint &src)
{
   E_GeV = src.E_GeV;
   dEdx_GeV_cm = src.dEdx_GeV_cm;
   phi_rad = src.phi_rad;
   primary_ = src.primary_;
   ptype_G3 = src.ptype_G3;
   px_GeV = src.px_GeV;
   py_GeV = src.py_GeV;
   pz_GeV = src.pz_GeV;
   r_cm = src.r_cm;
   z_cm = src.z_cm;
   t_ns = src.t_ns;
   sector_ = src.sector_;
   track_ = src.track_;
   trackID_ = src.trackID_;
}

int GlueXHitSTCpoint::operator==(const GlueXHitSTCpoint &right) const
{
   if (E_GeV       != right.E_GeV       ||
       dEdx_GeV_cm != right.dEdx_GeV_cm ||
       phi_rad     != right.phi_rad     ||
       primary_    != right.primary_    ||
       ptype_G3    != right.ptype_G3    ||
       px_GeV      != right.px_GeV      ||
       py_GeV      != right.py_GeV      ||
       pz_GeV      != right.pz_GeV      ||
       r_cm        != right.r_cm        ||
       z_cm        != right.z_cm        ||
       t_ns        != right.t_ns        ||
       sector_     != right.sector_     ||
       track_      != right.track_      ||
       trackID_    != right.trackID_    )
   {
      return 0;
   }
   return 1;
}

GlueXHitSTCpoint &GlueXHitSTCpoint::operator+=(const GlueXHitSTCpoint& /* right */)
{
   G4cerr << "Error in GlueXHitSTCpoint::operator+= - "
          << "illegal attempt to merge two TruthPoint objects in the stc!"
          << G4endl;
   return *this;
}

void GlueXHitSTCpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitSTCpoint::Print() const
{
   G4cout << "GlueXHitSTCpoint:" << G4endl
          << "   sector = " << sector_ << G4endl
          << "   track = " << track_ << G4endl
          << "   trackID = " << trackID_ << G4endl
          << "   E = " << E_GeV << " GeV" << G4endl
          << "   dEdx = " << dEdx_GeV_cm << " GeV/cm" << G4endl
          << "   phi = " << phi_rad << " rad" << G4endl
          << "   primary = " << primary_ << G4endl
          << "   ptype = " << ptype_G3 << G4endl
          << "   px = " << px_GeV << " GeV/c" << G4endl
          << "   py = " << py_GeV << " GeV/c" << G4endl
          << "   pz = " << pz_GeV << " GeV/c" << G4endl
          << "   r = " << r_cm << " cm" << G4endl
          << "   z = " << z_cm << " cm" << G4endl
          << "   t = " << t_ns << " ns" << G4endl
          << G4endl;
}

void printallhits(GlueXHitsMapSTCpoint *hitsmap)
{
   std::map<int, GlueXHitSTCpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitSTCpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

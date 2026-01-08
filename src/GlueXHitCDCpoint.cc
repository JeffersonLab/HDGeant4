//
// GlueXHitCDCpoint - class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 28, 2015

#include "GlueXHitCDCpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitCDCpoint>* GlueXHitCDCpointAllocator = 0;

GlueXHitCDCpoint::GlueXHitCDCpoint(const GlueXHitCDCpoint &src)
{
   dEdx_GeV_cm = src.dEdx_GeV_cm;
   dradius_cm = src.dradius_cm;
   phi_rad = src.phi_rad;
   primary_ = src.primary_;
   ptype_G3 = src.ptype_G3;
   px_GeV = src.px_GeV;
   py_GeV = src.py_GeV;
   pz_GeV = src.pz_GeV;
   r_cm = src.r_cm;
   z_cm = src.z_cm;
   t_ns = src.t_ns;
   track_ = src.track_;
   trackID_ = src.trackID_;
   sector_ = src.sector_;
   ring_ = src.ring_;
}

int GlueXHitCDCpoint::operator==(const GlueXHitCDCpoint &right) const
{
   if (dEdx_GeV_cm != right.dEdx_GeV_cm ||
       dradius_cm  != right.dradius_cm  ||
       phi_rad     != right.phi_rad     ||
       primary_    != right.primary_    ||
       ptype_G3    != right.ptype_G3    ||
       px_GeV      != right.px_GeV      ||
       py_GeV      != right.py_GeV      ||
       pz_GeV      != right.pz_GeV      ||
       r_cm        != right.r_cm        ||
       z_cm        != right.z_cm        ||
       t_ns        != right.t_ns        ||
       track_      != right.track_      ||
       trackID_    != right.trackID_    ||
       sector_     != right.sector_     ||
       ring_       != right.ring_       )
   {
      return 0;
   }
   return 1;
}

GlueXHitCDCpoint &GlueXHitCDCpoint::operator+=(const GlueXHitCDCpoint& /* right */)
{
   G4cerr << "Error in GlueXHitCDCpoint::operator+= - "
          << "illegal attempt to merge two TruthPoint objects in the cdc!"
          << G4endl;
   return *this;
}

void GlueXHitCDCpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitCDCpoint::Print() const
{
   G4cout << "GlueXHitCDCpoint:" << G4endl
          << "   track = " << track_ << G4endl
          << "   trackID = " << trackID_ << G4endl
          << "   dEdx = " << dEdx_GeV_cm << " GeV/cm" << G4endl
          << "   dradius = " << dradius_cm << " cm" << G4endl
          << "   phi = " << phi_rad << " rad" << G4endl
          << "   primary = " << primary_ << G4endl
          << "   ptype = " << ptype_G3 << G4endl
          << "   px = " << px_GeV << " GeV/c" << G4endl
          << "   py = " << py_GeV << " GeV/c" << G4endl
          << "   pz = " << pz_GeV << " GeV/c" << G4endl
          << "   r = " << r_cm << " cm" << G4endl
          << "   z = " << z_cm << " cm" << G4endl
          << "   t = " << t_ns << " ns" << G4endl
          << "   sector = " << sector_ << G4endl
          << "   ring = " << ring_ << G4endl
          << G4endl;
}

void printallhits(GlueXHitsMapCDCpoint *hitsmap)
{
   std::map<int, GlueXHitCDCpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCDCpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

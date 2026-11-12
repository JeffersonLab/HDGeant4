//
// GlueXHitGCALpoint - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016

#include "GlueXHitGCALpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitGCALpoint>* GlueXHitGCALpointAllocator = 0;

GlueXHitGCALpoint::GlueXHitGCALpoint(const GlueXHitGCALpoint &src)
{
   E_GeV = src.E_GeV;
   primary_ = src.primary_;
   ptype_G3 = src.ptype_G3;
   px_GeV = src.px_GeV;
   py_GeV = src.py_GeV;
   pz_GeV = src.pz_GeV;
   r_cm = src.r_cm;
   phi_rad = src.phi_rad;
   z_cm = src.z_cm;
   t_ns = src.t_ns;
   track_ = src.track_;
   trackID_ = src.trackID_;
}

int GlueXHitGCALpoint::operator==(const GlueXHitGCALpoint &right) const
{
   if (E_GeV    != right.E_GeV    ||
       primary_ != right.primary_ ||
       ptype_G3 != right.ptype_G3 ||
       px_GeV   != right.px_GeV   ||
       py_GeV   != right.py_GeV   ||
       pz_GeV   != right.pz_GeV   ||
       r_cm     != right.r_cm     ||
       phi_rad  != right.phi_rad  ||
       z_cm     != right.z_cm     ||
       t_ns     != right.t_ns     ||
       track_   != right.track_   ||
       trackID_ != right.trackID_ )
   {
      return 0;
   }
   return 1;
}

GlueXHitGCALpoint &GlueXHitGCALpoint::operator+=(const GlueXHitGCALpoint& /* right */)
{
   G4cerr << "Error in GlueXHitGCALpoint::operator+= - "
          << "illegal attempt to merge two TruthShower objects in the gcal!"
          << G4endl;
   return *this;
}

void GlueXHitGCALpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitGCALpoint::Print() const
{
   G4cout << "GlueXHitGCALpoint:" << G4endl
          << "   track = " << track_ << G4endl
          << "   trackID = " << trackID_ << G4endl
          << "   E = " << E_GeV << " GeV" << G4endl
          << "   primary = " << primary_ << G4endl
          << "   ptype = " << ptype_G3 << G4endl
          << "   px = " << px_GeV << " GeV/c" << G4endl
          << "   py = " << py_GeV << " GeV/c" << G4endl
          << "   pz = " << pz_GeV << " GeV/c" << G4endl
          << "   r = " << r_cm << " cm" << G4endl
          << "   phi = " << phi_rad << " radians" << G4endl
          << "   z = " << z_cm << " cm" << G4endl
          << "   t = " << t_ns << " ns" << G4endl
          << G4endl;
}

void printallhits(GlueXHitsMapGCALpoint *hitsmap)
{
   std::map<int, GlueXHitGCALpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitGCALpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

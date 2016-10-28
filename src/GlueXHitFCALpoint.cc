//
// GlueXHitFCALpoint - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 25, 2016

#include "GlueXHitFCALpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitFCALpoint>* GlueXHitFCALpointAllocator = 0;

int GlueXHitFCALpoint::operator==(const GlueXHitFCALpoint &right) const
{
   if (E_GeV    == right.E_GeV    &&
       primary_ == right.primary_ &&
       ptype_G3 == right.ptype_G3 &&
       px_GeV   == right.px_GeV   &&
       py_GeV   == right.py_GeV   &&
       pz_GeV   == right.pz_GeV   &&
       x_cm     == right.x_cm     &&
       y_cm     == right.y_cm     &&
       z_cm     == right.z_cm     &&
       t_ns     == right.t_ns     &&
       track_   == right.track_   &&
       trackID_ == right.trackID_ )
   {
      return 1;
   }
   return 0;
}

GlueXHitFCALpoint &GlueXHitFCALpoint::operator+=(const GlueXHitFCALpoint &right)
{
   G4cerr << "Error in GlueXHitFCALpoint::operator+= - "
          << "illegal attempt to merge two TruthShower objects in the fcal!"
          << G4endl;
   return *this;
}

void GlueXHitFCALpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitFCALpoint::Print() const
{
   G4cout << "GlueXHitFCALpoint:" << G4endl
          << "   track = " << track_ << G4endl
          << "   trackID = " << trackID_ << G4endl
          << "   E = " << E_GeV << " GeV" << G4endl
          << "   primary = " << primary_ << G4endl
          << "   ptype = " << ptype_G3 << G4endl
          << "   px = " << px_GeV << " GeV/c" << G4endl
          << "   py = " << py_GeV << " GeV/c" << G4endl
          << "   pz = " << pz_GeV << " GeV/c" << G4endl
          << "   x = " << x_cm << " cm" << G4endl
          << "   y = " << y_cm << " cm" << G4endl
          << "   z = " << z_cm << " cm" << G4endl
          << "   t = " << t_ns << " ns" << G4endl
          << G4endl;
}

void printallhits(GlueXHitsMapFCALpoint *hitsmap)
{
   std::map<int, GlueXHitFCALpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitFCALpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

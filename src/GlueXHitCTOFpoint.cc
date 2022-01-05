//
// GlueXHitCTOFpoint - class implementation
//
// author: staylor at jlab.org
// version: october 25, 2021

#include "GlueXHitCTOFpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitCTOFpoint>* GlueXHitCTOFpointAllocator = 0;

GlueXHitCTOFpoint::GlueXHitCTOFpoint(const GlueXHitCTOFpoint &src)
{
   E_GeV = src.E_GeV;
   primary_ = src.primary_;
   ptype_G3 = src.ptype_G3;
   px_GeV = src.px_GeV;
   py_GeV = src.py_GeV;
   pz_GeV = src.pz_GeV;
   x_cm = src.x_cm;
   y_cm = src.y_cm;
   z_cm = src.z_cm;
   t_ns = src.t_ns;
   track_ = src.track_;
   trackID_ = src.trackID_;
}

int GlueXHitCTOFpoint::operator==(const GlueXHitCTOFpoint &right) const
{
   if (E_GeV    != right.E_GeV    ||
       primary_ != right.primary_ ||
       ptype_G3 != right.ptype_G3 ||
       px_GeV   != right.px_GeV   ||
       py_GeV   != right.py_GeV   ||
       pz_GeV   != right.pz_GeV   ||
       x_cm     != right.x_cm     ||
       y_cm     != right.y_cm     ||
       z_cm     != right.z_cm     ||
       t_ns     != right.t_ns     ||
       track_   != right.track_   ||
       trackID_ != right.trackID_ )
   {
      return 0;
   }
   return 1;
}

GlueXHitCTOFpoint &GlueXHitCTOFpoint::operator+=(const GlueXHitCTOFpoint &right)
{
   G4cerr << "Error in GlueXHitCTOFpoint::operator+= - "
          << "illegal attempt to merge two TruthPoint objects in the fTOF!"
          << G4endl;
   return *this;
}

void GlueXHitCTOFpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitCTOFpoint::Print() const
{
   G4cout << "GlueXHitCTOFpoint:" << G4endl
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

void printallhits(GlueXHitsMapCTOFpoint *hitsmap)
{
   std::map<int, GlueXHitCTOFpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCTOFpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

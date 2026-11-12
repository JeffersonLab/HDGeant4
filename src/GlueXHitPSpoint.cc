//
// GlueXHitPSpoint - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016

#include "GlueXHitPSpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitPSpoint>* GlueXHitPSpointAllocator = 0;

GlueXHitPSpoint::GlueXHitPSpoint(const GlueXHitPSpoint &src)
{
   arm_ = src.arm_;
   column_ = src.column_;
   E_GeV = src.E_GeV;
   dEdx_GeV_cm = src.dEdx_GeV_cm;
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

int GlueXHitPSpoint::operator==(const GlueXHitPSpoint &right) const
{
   if (arm_        != right.arm_        ||
       column_     != right.column_     ||
       E_GeV       != right.E_GeV       ||
       dEdx_GeV_cm != right.dEdx_GeV_cm ||
       primary_    != right.primary_    ||
       ptype_G3    != right.ptype_G3    ||
       px_GeV      != right.px_GeV      ||
       py_GeV      != right.py_GeV      ||
       pz_GeV      != right.pz_GeV      ||
       x_cm        != right.x_cm        ||
       y_cm        != right.y_cm        ||
       z_cm        != right.z_cm        ||
       t_ns        != right.t_ns        ||
       track_      != right.track_      ||
       trackID_    != right.trackID_    )
   {
      return 0;
   }
   return 1;
}

GlueXHitPSpoint &GlueXHitPSpoint::operator+=(const GlueXHitPSpoint& /* right */)
{
   G4cerr << "Error in GlueXHitPSpoint::operator+= - "
          << "illegal attempt to merge two TruthPoint objects in the ps!"
          << G4endl;
   return *this;
}

void GlueXHitPSpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitPSpoint::Print() const
{
   G4cout << "GlueXHitPSpoint:" << G4endl
          << "   arm = " << arm_ << G4endl
          << "   column = " << column_ << G4endl
          << "   track = " << track_ << G4endl
          << "   trackID = " << trackID_ << G4endl
          << "   E = " << E_GeV << " GeV" << G4endl
          << "   dEdx = " << dEdx_GeV_cm << " GeV/cm" << G4endl
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

void printallhits(GlueXHitsMapPSpoint *hitsmap)
{
   std::map<int, GlueXHitPSpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitPSpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

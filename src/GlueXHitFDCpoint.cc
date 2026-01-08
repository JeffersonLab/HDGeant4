//
// GlueXHitFDCpoint - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016

#include "GlueXHitFDCpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitFDCpoint>* GlueXHitFDCpointAllocator = 0;

GlueXHitFDCpoint::GlueXHitFDCpoint(G4int chamber)
 : G4VHit(),
   chamber_(chamber)
{}

GlueXHitFDCpoint::GlueXHitFDCpoint(const GlueXHitFDCpoint &src)
{
   chamber_ = src.chamber_;
   E_GeV = src.E_GeV;
   dEdx_GeV_cm = src.dEdx_GeV_cm;
   dradius_cm = src.dradius_cm;
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

int GlueXHitFDCpoint::operator==(const GlueXHitFDCpoint &right) const
{
   if (E_GeV          != right.E_GeV       ||
       dEdx_GeV_cm    != right.dEdx_GeV_cm ||
       dradius_cm     != right.dradius_cm  ||
       primary_       != right.primary_    ||
       ptype_G3       != right.ptype_G3    ||
       px_GeV         != right.px_GeV      ||
       py_GeV         != right.py_GeV      ||
       pz_GeV         != right.pz_GeV      ||
       x_cm           != right.x_cm        ||
       y_cm           != right.y_cm        ||
       z_cm           != right.z_cm        ||
       t_ns           != right.t_ns        ||
       track_         != right.track_      ||
       trackID_       != right.trackID_    )
   {
      return 0;
   }
   return 1;
}

GlueXHitFDCpoint &GlueXHitFDCpoint::operator+=(const GlueXHitFDCpoint& /* right */)
{
   G4cerr << "Error in GlueXHitFDCpoint::operator+= - "
          << "illegal attempt to merge two TruthPoint objects in the fdc!"
          << G4endl;
   return *this;
}

void GlueXHitFDCpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitFDCpoint::Print() const
{
   G4cout << "GlueXHitFDCpoint:" << G4endl
          << "   track = " << track_ << G4endl
          << "   trackID = " << trackID_ << G4endl
          << "   E = " << E_GeV << " GeV" << G4endl
          << "   dEdx = " << dEdx_GeV_cm << " GeV/cm" << G4endl
          << "   dradius = " << dradius_cm << " cm" << G4endl
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

void printallhits(GlueXHitsMapFDCpoint *hitsmap)
{
   std::map<int, GlueXHitFDCpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitFDCpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

//
// GlueXHitDIRCflash - class implementation
//
// author: richard.t.jones at uconn.edu
// version: november 26, 2016
//
///////////////////////////////////////////////////////////////////
// NOTICE - This is not a real hits class for the DIRC detector.
// It was created during the DIRC conceptual design phase to measure
// the impact of the DIRC on GlueX physics channels. This needs to
// be retired as soon as the technical design is complete, and in
// its place a new class should be introduced with a name like
// GlueXHitDIRCtube that represents a single detected DIRC photon.
//
// DO NOT EDIT THIS FILE TO IMPLEMENT A REALISTIC DIRC HITS SCHEME,
// clone it instead, and give the new class a more sensible name!
// You also need to add a tag in event.xml to contain the new hit.
///////////////////////////////////////////////////////////////////

#include "GlueXHitDIRCflash.hh"

G4ThreadLocal G4Allocator<GlueXHitDIRCflash>* GlueXHitDIRCflashAllocator = 0;

GlueXHitDIRCflash::GlueXHitDIRCflash(G4int bar)
 : G4VHit(),
   bar_(bar)
{}

int GlueXHitDIRCflash::operator==(const GlueXHitDIRCflash &right) const
{
   if (bar_ !=  right.bar_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].E_GeV == right.hits[ih].E_GeV &&
          hits[ih].t_ns  == right.hits[ih].t_ns  &&
          hits[ih].x_cm  == right.hits[ih].x_cm  &&
          hits[ih].y_cm  == right.hits[ih].y_cm  &&
          hits[ih].z_cm  == right.hits[ih].z_cm)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitDIRCflash &GlueXHitDIRCflash::operator+=(const GlueXHitDIRCflash &right)
{
   if (bar_ != right.bar_) {
      G4cerr << "Error in GlueXHitDIRCflash::operator+=() - "
             << "illegal attempt to merge hits from two different flashs!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitDIRCflash::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitDIRCflash::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, GlueXHitDIRCflash::hitinfo_t());
      hiter->E_GeV = hitsrc->E_GeV;
      hiter->t_ns = hitsrc->t_ns;
      hiter->x_cm = hitsrc->x_cm;
      hiter->y_cm = hitsrc->y_cm;
      hiter->z_cm = hitsrc->z_cm;
   }
   return *this;
}

void GlueXHitDIRCflash::Draw() const
{
   // not yet implemented
}

void GlueXHitDIRCflash::Print() const
{
   G4cout << "GlueXHitDIRCflash: "
          << "   bar = " << bar_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   E = " << hiter->E_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   x = " << hiter->x_cm << " cm" << G4endl
             << "   y = " << hiter->y_cm << " cm" << G4endl
             << "   z = " << hiter->z_cm << " cm" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapDIRCflash *hitsmap)
{
   std::map<int, GlueXHitDIRCflash*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitDIRCflash*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

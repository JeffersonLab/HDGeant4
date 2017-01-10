//
// GlueXHitFDCwire - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016

#include "GlueXHitFDCwire.hh"

G4ThreadLocal G4Allocator<GlueXHitFDCwire>* GlueXHitFDCwireAllocator = 0;

GlueXHitFDCwire::GlueXHitFDCwire(G4int chamber, G4int wire)
 : G4VHit(),
   chamber_(chamber),
   wire_(wire)
{}

int GlueXHitFDCwire::operator==(const GlueXHitFDCwire &right) const
{
   if (chamber_ != right.chamber_  || wire_ !=  right.wire_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].dE_keV   == right.hits[ih].dE_keV  &&
          hits[ih].t_ns     == right.hits[ih].t_ns    &&
          hits[ih].d_cm     == right.hits[ih].d_cm    &&
          hits[ih].itrack_  == right.hits[ih].itrack_ &&
          hits[ih].t0_ns    == right.hits[ih].t0_ns   &&
          hits[ih].t1_ns    == right.hits[ih].t1_ns   &&
          hits[ih].ptype_G3 == right.hits[ih].ptype_G3)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitFDCwire &GlueXHitFDCwire::operator+=(const GlueXHitFDCwire &right)
{
   if (chamber_ != right.chamber_  || wire_ !=  right.wire_) {
      G4cerr << "Error in GlueXHitFDCwire::operator+=() - "
             << "illegal attempt to merge hits from two different wires!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitFDCwire::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitFDCwire::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, GlueXHitFDCwire::hitinfo_t());
      hiter->dE_keV = hitsrc->dE_keV;
      hiter->t_ns = hitsrc->t_ns;
      hiter->d_cm = hitsrc->d_cm;
      hiter->itrack_ = hitsrc->itrack_;
      hiter->ptype_G3 = hitsrc->ptype_G3;
      hiter->t0_ns = hitsrc->t0_ns;
      hiter->t1_ns = hitsrc->t1_ns;
   }
   return *this;
}

void GlueXHitFDCwire::Draw() const
{
   // not yet implemented
}

void GlueXHitFDCwire::Print() const
{
   G4cout << "GlueXHitFDCwire: "
          << "   chamber = " << chamber_
          << "   wire = " << wire_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   dE = " << hiter->dE_keV << " keV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   d = " << hiter->d_cm << " cm" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << "   ptype = " << hiter->ptype_G3 << G4endl
             << "   t0 = " << hiter->t0_ns << " ns" << G4endl
             << "   t1 = " << hiter->t1_ns << " ns" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapFDCwire *hitsmap)
{
   std::map<int, GlueXHitFDCwire*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitFDCwire*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

//
// GlueXHitFMWPCwire - class implementation
//
// author: richard.t.jones at uconn.edu
// version: november 28, 2016

#include "GlueXHitFMWPCwire.hh"

G4ThreadLocal G4Allocator<GlueXHitFMWPCwire>* GlueXHitFMWPCwireAllocator = 0;

GlueXHitFMWPCwire::GlueXHitFMWPCwire(G4int layer, G4int wire)
 : G4VHit(),
   layer_(layer),
   wire_(wire),
   overflow_(0)
{}

GlueXHitFMWPCwire::GlueXHitFMWPCwire(const GlueXHitFMWPCwire &src)
{
   layer_ = src.layer_;
   wire_ = src.wire_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitFMWPCwire::operator==(const GlueXHitFMWPCwire &right) const
{
   if (layer_ !=  right.layer_ || wire_ != right.wire_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].dE_keV   != right.hits[ih].dE_keV  ||
          hits[ih].t_ns     != right.hits[ih].t_ns    ||
          hits[ih].d_cm    != right.hits[ih].d_cm   ||
          hits[ih].itrack_  != right.hits[ih].itrack_ 
          )
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitFMWPCwire &GlueXHitFMWPCwire::operator+=(const GlueXHitFMWPCwire &right)
{
   if (layer_ !=  right.layer_ || wire_ != right.wire_) {
      G4cerr << "Error in GlueXHitFMWPCwire::operator+=() - "
             << "illegal attempt to merge hits from two different wires!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitFMWPCwire::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitFMWPCwire::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, *hitsrc);
   }
   overflow_ += right.overflow_;
   return *this;
}

void GlueXHitFMWPCwire::Draw() const
{
   // not yet implemented
}

void GlueXHitFMWPCwire::Print() const
{
   G4cout << "GlueXHitFMWPCwire: "
          << "   layer = " << layer_ 
          << ",  wire = " << wire_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   dE = " << hiter->dE_keV << " keV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   d = " << hiter->d_cm << " cm" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
	     << G4endl;
   }
}

void printallhits(GlueXHitsMapFMWPCwire *hitsmap)
{
   std::map<int, GlueXHitFMWPCwire*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitFMWPCwire*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

//
// GlueXHitPStile - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016

#include "GlueXHitPStile.hh"

G4ThreadLocal G4Allocator<GlueXHitPStile>* GlueXHitPStileAllocator = 0;

GlueXHitPStile::GlueXHitPStile(G4int arm, G4int column)
 : G4VHit(),
   arm_(arm),
   column_(column),
   overflow_(0)
{}

GlueXHitPStile::GlueXHitPStile(const GlueXHitPStile &src)
{
   arm_ = src.arm_;
   column_ = src.column_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitPStile::operator==(const GlueXHitPStile &right) const
{
   if (arm_ !=  right.arm_ || column_ != right.column_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].dE_GeV   != right.hits[ih].dE_GeV  ||
          hits[ih].t_ns     != right.hits[ih].t_ns    ||
          hits[ih].itrack_  != right.hits[ih].itrack_ ||
          hits[ih].ptype_G3 != right.hits[ih].ptype_G3)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitPStile &GlueXHitPStile::operator+=(const GlueXHitPStile &right)
{
   if (arm_ !=  right.arm_ || column_ != right.column_) {
      G4cerr << "Error in GlueXHitPStile::operator+=() - "
             << "illegal attempt to merge hits from two different tiles!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitPStile::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitPStile::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitPStile::Draw() const
{
   // not yet implemented
}

void GlueXHitPStile::Print() const
{
   G4cout << "GlueXHitPStile: "
          << "   arm = " << arm_
          << ",  column = " << column_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   dE = " << hiter->dE_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << "   ptype = " << hiter->ptype_G3 << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapPStile *hitsmap)
{
   std::map<int, GlueXHitPStile*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitPStile*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

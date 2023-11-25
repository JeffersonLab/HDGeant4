//
// GlueXHitGCALblock - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016

#include "GlueXHitGCALblock.hh"

G4ThreadLocal G4Allocator<GlueXHitGCALblock>* GlueXHitGCALblockAllocator = 0;

GlueXHitGCALblock::GlueXHitGCALblock(G4int module)
 : G4VHit(),
   module_(module),
   overflow_(0)
{}

GlueXHitGCALblock::GlueXHitGCALblock(const GlueXHitGCALblock &src)
{
   module_ = src.module_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitGCALblock::operator==(const GlueXHitGCALblock &right) const
{
   if (module_ !=  right.module_) {
      return 0;
   }
   else if (hits.size() != right.hits.size()) {
      return 0;
   }

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].E_GeV != right.hits[ih].E_GeV ||
          hits[ih].t_ns != right.hits[ih].t_ns)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitGCALblock &GlueXHitGCALblock::operator+=(const GlueXHitGCALblock &right)
{
   if (module_ !=  right.module_) {
      G4cerr << "Error in GlueXHitGCALblock::operator+=() - "
             << "illegal attempt to merge hits from two different blocks!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitGCALblock::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitGCALblock::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitGCALblock::Draw() const
{
   // not yet implemented
}

void GlueXHitGCALblock::Print() const
{
   G4cout << "GlueXHitGCALblock: "
          << "   module = " << module_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   E = " << hiter->E_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   zLocal = " << hiter->zlocal_cm << " cm" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapGCALblock *hitsmap)
{
   std::map<int, GlueXHitGCALblock*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitGCALblock*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

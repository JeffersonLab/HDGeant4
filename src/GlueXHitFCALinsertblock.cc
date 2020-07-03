//
// GlueXHitFCALinsertblock - class implementation
//

#include "GlueXHitFCALinsertblock.hh"

G4ThreadLocal G4Allocator<GlueXHitFCALinsertblock>* GlueXHitFCALinsertblockAllocator = 0;

GlueXHitFCALinsertblock::GlueXHitFCALinsertblock(G4int column, G4int row)
 : G4VHit(),
   column_(column),
   row_(row)
{}

GlueXHitFCALinsertblock::GlueXHitFCALinsertblock(const GlueXHitFCALinsertblock &src)
{
   column_ = src.column_;
   row_ = src.row_;
   hits = src.hits;
}

int GlueXHitFCALinsertblock::operator==(const GlueXHitFCALinsertblock &right) const
{
   if (column_ !=  right.column_ || row_ != right.row_) {
      return 0;
   }
   else if (hits.size() != right.hits.size()) {
      return 0;
   }

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].E_GeV != right.hits[ih].E_GeV ||
          hits[ih].t_ns != right.hits[ih].t_ns ||
          hits[ih].dE_lightguide_GeV != right.hits[ih].dE_lightguide_GeV ||
          hits[ih].t_lightguide_ns != right.hits[ih].t_lightguide_ns)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitFCALinsertblock &GlueXHitFCALinsertblock::operator+=(const GlueXHitFCALinsertblock &right)
{
   if (column_ !=  right.column_ || row_ != right.row_) {
      G4cerr << "Error in GlueXHitFCALinsertblock::operator+=() - "
             << "illegal attempt to merge hits from two different blocks!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitFCALinsertblock::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitFCALinsertblock::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, *hitsrc);
   }
   return *this;
}

void GlueXHitFCALinsertblock::Draw() const
{
   // not yet implemented
}

void GlueXHitFCALinsertblock::Print() const
{
   G4cout << "GlueXHitFCALinsertblock: "
          << "  column = " << column_
          << ", row = " << row_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   E = " << hiter->E_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   E(lightguide) = " 
             << hiter->dE_lightguide_GeV << " GeV" << G4endl
             << "   t(lightguide) = " 
             << hiter->t_lightguide_ns << " ns" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapFCALinsertblock *hitsmap)
{
   std::map<int, GlueXHitFCALinsertblock*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitFCALinsertblock*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

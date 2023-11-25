//
// GlueXHitCCALblock - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 28, 2016

#include "GlueXHitCCALblock.hh"

G4ThreadLocal G4Allocator<GlueXHitCCALblock>* GlueXHitCCALblockAllocator = 0;

GlueXHitCCALblock::GlueXHitCCALblock(G4int column, G4int row)
 : G4VHit(),
   column_(column),
   row_(row),
   overflow_(0)
{}

GlueXHitCCALblock::GlueXHitCCALblock(const GlueXHitCCALblock &src)
{
   column_ = src.column_;
   row_ = src.row_;
   hits = src.hits;
   overflow_ = src.overflow_;
}

int GlueXHitCCALblock::operator==(const GlueXHitCCALblock &right) const
{
   if (column_ !=  right.column_ || row_ != right.row_) {
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

GlueXHitCCALblock &GlueXHitCCALblock::operator+=(const GlueXHitCCALblock &right)
{
   if (column_ !=  right.column_ || row_ != right.row_) {
      G4cerr << "Error in GlueXHitCCALblock::operator+=() - "
             << "illegal attempt to merge hits from two different blocks!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitCCALblock::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitCCALblock::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitCCALblock::Draw() const
{
   // not yet implemented
}

void GlueXHitCCALblock::Print() const
{
   G4cout << "GlueXHitCCALblock: "
          << "  column = " << column_
          << ", row = " << row_
          << ", overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   E = " << hiter->E_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapCCALblock *hitsmap)
{
   std::map<int, GlueXHitCCALblock*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCCALblock*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

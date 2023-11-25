//
// GlueXHitUPVbar - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 29, 2016

#include "GlueXHitUPVbar.hh"

G4ThreadLocal G4Allocator<GlueXHitUPVbar>* GlueXHitUPVbarAllocator = 0;

GlueXHitUPVbar::GlueXHitUPVbar(G4int layer, G4int row)
 : G4VHit(),
   layer_(layer),
   row_(row),
   overflow_(0)
{}

GlueXHitUPVbar::GlueXHitUPVbar(const GlueXHitUPVbar &src)
{
   layer_ = src.layer_;
   row_ = src.row_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitUPVbar::operator==(const GlueXHitUPVbar &right) const
{
   if (layer_ !=  right.layer_ || row_ != right.row_ )
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].end_      != right.hits[ih].end_   ||
          hits[ih].E_GeV     != right.hits[ih].E_GeV  ||
          hits[ih].t_ns      != right.hits[ih].t_ns   ||
          hits[ih].xlocal_cm != right.hits[ih].xlocal_cm)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitUPVbar &GlueXHitUPVbar::operator+=(const GlueXHitUPVbar &right)
{
   if (layer_ !=  right.layer_ || row_ != right.row_) {
      G4cerr << "Error in GlueXHitUPVbar::operator+=() - "
             << "illegal attempt to merge hits from two different bars!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitUPVbar::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitUPVbar::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitUPVbar::Draw() const
{
   // not yet implemented
}

void GlueXHitUPVbar::Print() const
{
   G4cout << "GlueXHitUPVbar: "
          << "   layer = " << layer_ 
          << ",  row = " << row_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   end = " << hiter->end_ << G4endl
             << "   E = " << hiter->E_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   xlocal = " << hiter->xlocal_cm << " cm" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapUPVbar *hitsmap)
{
   std::map<int, GlueXHitUPVbar*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitUPVbar*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

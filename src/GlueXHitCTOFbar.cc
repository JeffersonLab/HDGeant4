//
// GlueXHitCTOFbar - class implementation
//
// author: staylor at jlab.org
// version: october 25, 2016

#include "GlueXHitCTOFbar.hh"

G4ThreadLocal G4Allocator<GlueXHitCTOFbar>* GlueXHitCTOFbarAllocator = 0;

GlueXHitCTOFbar::GlueXHitCTOFbar(G4int bar)
 : G4VHit(),
   bar_(bar)
{}

GlueXHitCTOFbar::GlueXHitCTOFbar(const GlueXHitCTOFbar &src)
{
   bar_ = src.bar_;
   hits = src.hits;
}

int GlueXHitCTOFbar::operator==(const GlueXHitCTOFbar &right) const
{
   if (bar_ != right.bar_ )
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].end_     != right.hits[ih].end_    ||
          hits[ih].dE_GeV   != right.hits[ih].dE_GeV  ||
          hits[ih].t_ns     != right.hits[ih].t_ns   
	  )
      {
         return 0;
      }

   }
   return 1;
}

GlueXHitCTOFbar &GlueXHitCTOFbar::operator+=(const GlueXHitCTOFbar &right)
{
   if (bar_ != right.bar_) {
      G4cerr << "Error in GlueXHitCTOFbar::operator+=() - "
             << "illegal attempt to merge hits from two different bars!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitCTOFbar::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitCTOFbar::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, *hitsrc);
   }
   return *this;
}

void GlueXHitCTOFbar::Draw() const
{
   // not yet implemented
}

void GlueXHitCTOFbar::Print() const
{
   G4cout << "GlueXHitCTOFbar: "
          << "   bar = " << bar_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   end = " << hiter->end_ << G4endl
             << "   dE = " << hiter->dE_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl;
      G4cout << G4endl;
   }
}

void printallhits(GlueXHitsMapCTOFbar *hitsmap)
{
   std::map<int, GlueXHitCTOFbar*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCTOFbar*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

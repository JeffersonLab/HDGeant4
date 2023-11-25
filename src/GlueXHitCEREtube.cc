//
// GlueXHitCEREtube - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 29, 2016

#include "GlueXHitCEREtube.hh"

G4ThreadLocal G4Allocator<GlueXHitCEREtube>* GlueXHitCEREtubeAllocator = 0;

GlueXHitCEREtube::GlueXHitCEREtube(G4int sector)
 : G4VHit(),
   sector_(sector),
   overflow_(0)
{}

GlueXHitCEREtube::GlueXHitCEREtube(const GlueXHitCEREtube &src)
{
   sector_ = src.sector_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitCEREtube::operator==(const GlueXHitCEREtube &right) const
{
   if (sector_ !=  right.sector_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].pe_     != right.hits[ih].pe_     ||
          hits[ih].t_ns    != right.hits[ih].t_ns    ||
          hits[ih].itrack_ != right.hits[ih].itrack_ ||
          hits[ih].z0_cm   != right.hits[ih].z0_cm)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitCEREtube &GlueXHitCEREtube::operator+=(const GlueXHitCEREtube &right)
{
   if (sector_ !=  right.sector_) {
      G4cerr << "Error in GlueXHitCEREtube::operator+=() - "
             << "illegal attempt to merge hits from two different tubes!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitCEREtube::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitCEREtube::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitCEREtube::Draw() const
{
   // not yet implemented
}

void GlueXHitCEREtube::Print() const
{
   G4cout << "GlueXHitCEREtube: "
          << "   sector = " << sector_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   pe = " << hiter->pe_ << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << "   z0 = " << hiter->z0_cm << " cm" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapCEREtube *hitsmap)
{
   std::map<int, GlueXHitCEREtube*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCEREtube*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

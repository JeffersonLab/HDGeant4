//
// GlueXHitCDCstraw - class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 28, 2015

#include "GlueXHitCDCstraw.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4Allocator<GlueXHitCDCstraw>* GlueXHitCDCstrawAllocator = 0;

GlueXHitCDCstraw::GlueXHitCDCstraw(G4int ring, G4int sector)
 : G4VHit(),
   ring_(ring),
   sector_(sector),
   overflow_(0)
{}

GlueXHitCDCstraw::GlueXHitCDCstraw(const GlueXHitCDCstraw &src)
{
   ring_ = src.ring_;
   sector_ = src.sector_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitCDCstraw::operator==(const GlueXHitCDCstraw &right) const
{
   if (ring_ != right.ring_  || sector_ !=  right.sector_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].track_   != right.hits[ih].track_  ||
          hits[ih].q_fC     != right.hits[ih].q_fC    ||
          hits[ih].t_ns     != right.hits[ih].t_ns    ||
          hits[ih].d_cm     != right.hits[ih].d_cm    ||
          hits[ih].itrack_  != right.hits[ih].itrack_ ||
          hits[ih].t0_ns    != right.hits[ih].t0_ns   ||
          hits[ih].t1_ns    != right.hits[ih].t1_ns   ||
          hits[ih].x0_g     != right.hits[ih].x0_g    ||
          hits[ih].x1_g     != right.hits[ih].x1_g    ||
          hits[ih].ptype_G3 != right.hits[ih].ptype_G3)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitCDCstraw &GlueXHitCDCstraw::operator+=(const GlueXHitCDCstraw &right)
{
   if (ring_ != right.ring_  || sector_ !=  right.sector_) {
      G4cerr << "Error in GlueXHitCDCstraw::operator+=() - "
             << "illegal attempt to merge hits from two different straws!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitCDCstraw::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitCDCstraw::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitCDCstraw::Draw() const
{
   // not yet implemented
}

void GlueXHitCDCstraw::Print() const
{
   G4cout << "GlueXHitCDCstraw: "
          << "   ring = " << ring_
          << ",  sector = " << sector_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   q = " << hiter->q_fC << " fC" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   d = " << hiter->d_cm << " cm" << G4endl
             << "   track = " << hiter->track_ << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << "   ptype = " << hiter->ptype_G3 << G4endl
             << "   t0 = " << hiter->t0_ns << " ns" << G4endl
             << "   t1 = " << hiter->t1_ns << " ns" << G4endl
             << "   x0 = (" << hiter->x0_g[0]/cm << ","
                            << hiter->x0_g[1]/cm << ","
                            << hiter->x0_g[2]/cm << ")cm"
                            << " in global coordinates" << G4endl
             << "   x1 = (" << hiter->x1_g[0]/cm << ","
                            << hiter->x1_g[1]/cm << ","
                            << hiter->x1_g[2]/cm << ")cm"
                            << " in global coordinates" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapCDCstraw *hitsmap)
{
   std::map<int, GlueXHitCDCstraw*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCDCstraw*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

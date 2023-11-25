//
// GlueXHitSTCpaddle - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016

#include "GlueXHitSTCpaddle.hh"

G4ThreadLocal G4Allocator<GlueXHitSTCpaddle>* GlueXHitSTCpaddleAllocator = 0;

GlueXHitSTCpaddle::GlueXHitSTCpaddle(G4int sector)
 : G4VHit(),
   sector_(sector),
   overflow_(0)
{}

GlueXHitSTCpaddle::GlueXHitSTCpaddle(const GlueXHitSTCpaddle &src)
{
   sector_ = src.sector_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitSTCpaddle::operator==(const GlueXHitSTCpaddle &right) const
{
   if (sector_ !=  right.sector_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].dE_MeV   != right.hits[ih].dE_MeV  ||
          hits[ih].t_ns     != right.hits[ih].t_ns    ||
          hits[ih].itrack_  != right.hits[ih].itrack_ ||
          hits[ih].t0_ns    != right.hits[ih].t0_ns   ||
          hits[ih].z_cm     != right.hits[ih].z_cm    ||
          hits[ih].ptype_G3 != right.hits[ih].ptype_G3)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitSTCpaddle &GlueXHitSTCpaddle::operator+=(const GlueXHitSTCpaddle &right)
{
   if (sector_ !=  right.sector_) {
      G4cerr << "Error in GlueXHitSTCpaddle::operator+=() - "
             << "illegal attempt to merge hits from two different paddles!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitSTCpaddle::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitSTCpaddle::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitSTCpaddle::Draw() const
{
   // not yet implemented
}

void GlueXHitSTCpaddle::Print() const
{
   G4cout << "GlueXHitSTCpaddle: "
          << "   sector = " << sector_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   dE = " << hiter->dE_MeV << " MeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << "   ptype = " << hiter->ptype_G3 << G4endl
             << "   t0 = " << hiter->t0_ns << " ns" << G4endl
             << "   z = " << hiter->z_cm << " cm" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapSTCpaddle *hitsmap)
{
   std::map<int, GlueXHitSTCpaddle*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitSTCpaddle*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

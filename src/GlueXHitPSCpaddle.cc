//
// GlueXHitPSCpaddle - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016

#include "GlueXHitPSCpaddle.hh"

G4ThreadLocal G4Allocator<GlueXHitPSCpaddle>* GlueXHitPSCpaddleAllocator = 0;

GlueXHitPSCpaddle::GlueXHitPSCpaddle(G4int arm, G4int module)
 : G4VHit(),
   arm_(arm),
   module_(module),
   overflow_(0)
{}

GlueXHitPSCpaddle::GlueXHitPSCpaddle(const GlueXHitPSCpaddle &src)
{
   arm_ = src.arm_;
   module_ = src.module_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitPSCpaddle::operator==(const GlueXHitPSCpaddle &right) const
{
   if (module_ !=  right.module_ || arm_ != right.arm_)
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

GlueXHitPSCpaddle &GlueXHitPSCpaddle::operator+=(const GlueXHitPSCpaddle &right)
{
   if (module_ !=  right.module_ || arm_ != right.arm_) {
      G4cerr << "Error in GlueXHitPSCpaddle::operator+=() - "
             << "illegal attempt to merge hits from two different paddles!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitPSCpaddle::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitPSCpaddle::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitPSCpaddle::Draw() const
{
   // not yet implemented
}

void GlueXHitPSCpaddle::Print() const
{
   G4cout << "GlueXHitPSCpaddle: "
          << "   arm = " << arm_
          << ",  module = " << module_
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

void printallhits(GlueXHitsMapPSCpaddle *hitsmap)
{
   std::map<int, GlueXHitPSCpaddle*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitPSCpaddle*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

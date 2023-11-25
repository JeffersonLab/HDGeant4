//
// GlueXHitFDCcathode - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016

#include "GlueXHitFDCcathode.hh"

G4ThreadLocal G4Allocator<GlueXHitFDCcathode>* GlueXHitFDCcathodeAllocator = 0;

GlueXHitFDCcathode::GlueXHitFDCcathode(G4int chamber, G4int plane, G4int strip)
 : G4VHit(),
   chamber_(chamber),
   plane_(plane),
   strip_(strip),
   overflow_(0)
{}

GlueXHitFDCcathode::GlueXHitFDCcathode(const GlueXHitFDCcathode &src)
{
   chamber_ = src.chamber_;
   plane_ = src.plane_;
   strip_ = src.strip_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitFDCcathode::operator==(const GlueXHitFDCcathode &right) const
{
   if (chamber_ != right.chamber_  ||
       plane_ !=  right.plane_ ||
       strip_ != right.strip_)
   {
      return 0;
   }
   else if (hits.size() != right.hits.size()) {
      return 0;
   }
   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].q_fC     != right.hits[ih].q_fC    ||
          hits[ih].t_ns     != right.hits[ih].t_ns    ||
          hits[ih].itrack_  != right.hits[ih].itrack_ ||
          hits[ih].u_cm     != right.hits[ih].u_cm    ||
          hits[ih].v_cm     != right.hits[ih].v_cm    ||
          hits[ih].ptype_G3 != right.hits[ih].ptype_G3)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitFDCcathode &GlueXHitFDCcathode::operator+=(const GlueXHitFDCcathode &right)
{
   if (chamber_ != right.chamber_  || 
       plane_ != right.plane_ ||
       strip_ != right.strip_)
   {
      G4cerr << "Error in GlueXHitFDCcathode::operator+=() - "
             << "illegal attempt to merge hits from two different strips!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitFDCcathode::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitFDCcathode::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitFDCcathode::Draw() const
{
   // not yet implemented
}

void GlueXHitFDCcathode::Print() const
{
   G4cout << "GlueXHitFDCcathode: "
          << "   chamber = " << chamber_
          << ",  plane = " << plane_
          << ",  strip = " << strip_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   q = " << hiter->q_fC << " fC" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << "   ptype = " << hiter->ptype_G3 << G4endl
             << "   u = " << hiter->u_cm << " cm" << G4endl
             << "   v = " << hiter->v_cm << " cm" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapFDCcathode *hitsmap)
{
   std::map<int, GlueXHitFDCcathode*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitFDCcathode*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

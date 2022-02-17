//
// GlueXHitGEMTRDhit - class implementation
//
// author: staylor at jlab.org
// version: february 7, 2022

#include "GlueXHitGEMTRDhit.hh"

G4ThreadLocal G4Allocator<GlueXHitGEMTRDhit>* GlueXHitGEMTRDhitAllocator = 0;

GlueXHitGEMTRDhit::GlueXHitGEMTRDhit(G4int layer)
 : G4VHit(),
   layer_(layer)
{}

GlueXHitGEMTRDhit::GlueXHitGEMTRDhit(const GlueXHitGEMTRDhit &src)
{
   layer_ = src.layer_;
   hits = src.hits;
}

int GlueXHitGEMTRDhit::operator==(const GlueXHitGEMTRDhit &right) const
{
  if (layer_ !=  right.layer_)
    return 0;
  else if (hits.size() != right.hits.size())
    return 0;
  
  for (int ih=0; ih < (int)hits.size(); ++ih) {
    if (hits[ih].q_fC     != right.hits[ih].q_fC      ||
	hits[ih].t_ns     != right.hits[ih].t_ns      ||
	hits[ih].x_cm     != right.hits[ih].x_cm      ||
	hits[ih].y_cm     != right.hits[ih].y_cm      ||
	hits[ih].d_cm     != right.hits[ih].d_cm      ||
	hits[ih].itrack_  != right.hits[ih].itrack_ 
	)
      {
	return 0;
      }
  }
  return 1;
}

GlueXHitGEMTRDhit &GlueXHitGEMTRDhit::operator+=(const GlueXHitGEMTRDhit &right)
{
   if (layer_ !=  right.layer_) {
      G4cerr << "Error in GlueXHitGEMTRDhit::operator+=() - "
             << "illegal attempt to merge hits from two different layers!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitGEMTRDhit::hitinfo_t>::iterator titer = hits.begin();
   std::vector<GlueXHitGEMTRDhit::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; titer != hits.end(); ++titer) {
         if (titer->t_ns > hitsrc->t_ns)
            break;
      }
      titer = hits.insert(titer, *hitsrc);
   }
   return *this;
}

void GlueXHitGEMTRDhit::Draw() const
{
   // not yet implemented
}

void GlueXHitGEMTRDhit::Print() const
{
   G4cout << "GlueXHitGEMTRDhit: "
          << "   layer = " << layer_ << G4endl;
   std::vector<hitinfo_t>::const_iterator titer;
   for (titer = hits.begin(); titer != hits.end(); ++titer) {
     G4cout << "   q = " << titer->q_fC << " fC" << G4endl
	    << "   t = " << titer->t_ns << " ns" << G4endl
	    << "   x = " << titer->x_cm << " cm" << G4endl
	    << "   y = " << titer->y_cm << " cm" << G4endl
	    << "   d = " << titer->d_cm << G4endl
	    << "   itrack = " << titer->itrack_ << G4endl
	    << G4endl;
   }
}

void printallhits(GlueXHitsMapGEMTRDhit *hitsmap)
{
   std::map<int, GlueXHitGEMTRDhit*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitGEMTRDhit*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

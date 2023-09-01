//
// GlueXHitITOFhit - class implementation
//
// author: staylor at jlab.org
// version: february 7, 2022

#include "GlueXHitITOFhit.hh"

G4ThreadLocal G4Allocator<GlueXHitITOFhit>* GlueXHitITOFhitAllocator = 0;

GlueXHitITOFhit::GlueXHitITOFhit(G4int id)
 : G4VHit(),
   id_(id)
{}

GlueXHitITOFhit::GlueXHitITOFhit(const GlueXHitITOFhit &src)
{
  id_ = src.id_;
  hits = src.hits;
}

int GlueXHitITOFhit::operator==(const GlueXHitITOFhit &right) const
{ 
  if (id_ != right.id_ )
      return 0;
  else if (hits.size() != right.hits.size())
    return 0;
  
  for (int ih=0; ih < (int)hits.size(); ++ih) {
    if (hits[ih].dE_GeV     != right.hits[ih].dE_GeV      ||
	hits[ih].t_ns     != right.hits[ih].t_ns      ||
	hits[ih].x_cm     != right.hits[ih].x_cm      ||
	hits[ih].y_cm     != right.hits[ih].y_cm
	)
      {
	return 0;
      }
  }
  return 1;
}

GlueXHitITOFhit &GlueXHitITOFhit::operator+=(const GlueXHitITOFhit &right)
{  
  if (id_ != right.id_) {
      G4cerr << "Error in GlueXHitCTOFbar::operator+=() - "
             << "illegal attempt to merge hits from two different ids!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitITOFhit::hitinfo_t>::iterator titer = hits.begin();
   std::vector<GlueXHitITOFhit::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; titer != hits.end(); ++titer) {
         if (titer->t_ns > hitsrc->t_ns)
            break;
      }
      titer = hits.insert(titer, *hitsrc);
   }
   return *this;
}

void GlueXHitITOFhit::Draw() const
{
   // not yet implemented
}

void GlueXHitITOFhit::Print() const
{
   G4cout << "GlueXHitITOFhit: " 
	  << "   id = " << id_  << G4endl;
   std::vector<hitinfo_t>::const_iterator titer;
   for (titer = hits.begin(); titer != hits.end(); ++titer) {
     G4cout << "   dE = " << titer->dE_GeV << " GeV" << G4endl
	    << "   t = " << titer->t_ns << " ns" << G4endl
	    << "   x = " << titer->x_cm << " cm" << G4endl
	    << "   y = " << titer->y_cm << " cm" << G4endl
	    << G4endl;
   }
}

void printallhits(GlueXHitsMapITOFhit *hitsmap)
{
   std::map<int, GlueXHitITOFhit*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitITOFhit*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

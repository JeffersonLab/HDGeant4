//
// GlueXHitBCALcell - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 24, 2016

#include "GlueXHitBCALcell.hh"

G4ThreadLocal G4Allocator<GlueXHitBCALcell>* GlueXHitBCALcellAllocator = 0;

GlueXHitBCALcell::GlueXHitBCALcell(G4int module, G4int layer, G4int sector)
 : G4VHit(),
   module_(module),
   layer_(layer),
   sector_(sector)
{}

int GlueXHitBCALcell::operator==(const GlueXHitBCALcell &right) const
{
   if (module_ !=  right.module_ || layer_ != right.layer_ || 
       sector_ !=  right.sector_)
   {
      return 0;
   }
   else if (hits.size() != right.hits.size()) {
      return 0;
   }

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].E_GeV     == right.hits[ih].E_GeV     &&
          hits[ih].zlocal_cm == right.hits[ih].zlocal_cm &&
          hits[ih].t_ns      == right.hits[ih].t_ns      &&
          hits[ih].itrack_   == right.hits[ih].itrack_ )
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitBCALcell &GlueXHitBCALcell::operator+=(const GlueXHitBCALcell &right)
{
   if (module_ !=  right.module_ || layer_ != right.layer_ || 
       sector_ !=  right.sector_)
   {
      G4cerr << "Error in GlueXHitBCALcell::operator+=() - "
             << "illegal attempt to merge hits from two different cells!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitBCALcell::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitBCALcell::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, GlueXHitBCALcell::hitinfo_t());
      hiter->E_GeV = hitsrc->E_GeV;
      hiter->zlocal_cm = hitsrc->zlocal_cm;
      hiter->t_ns = hitsrc->t_ns;
      hiter->itrack_ = hitsrc->itrack_;
   }
   return *this;
}

void GlueXHitBCALcell::Draw() const
{
   // not yet implemented
}

void GlueXHitBCALcell::Print() const
{
   G4cout << "GlueXHitBCALcell: "
          << "  module = " << module_
          << ", layer = " << layer_
          << ", sector = " << sector_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   E = " << hiter->E_GeV << " GeV" << G4endl
             << "   zlocal = " << hiter->t_ns << " cm" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapBCALcell *hitsmap)
{
   std::map<int, GlueXHitBCALcell*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitBCALcell*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

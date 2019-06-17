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

GlueXHitBCALcell::GlueXHitBCALcell(const GlueXHitBCALcell &src)
{
   module_ = src.module_;
   layer_ = src.layer_;
   sector_ = src.sector_;
   hits = src.hits;
}

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
      if (hits[ih].E_GeV       != right.hits[ih].E_GeV     ||
          hits[ih].zlocal_cm   != right.hits[ih].zlocal_cm ||
          hits[ih].t_ns        != right.hits[ih].t_ns      ||
          hits[ih].Eup_GeV     != right.hits[ih].Eup_GeV   ||
          hits[ih].Edown_GeV   != right.hits[ih].Edown_GeV ||
          hits[ih].tup_ns      != right.hits[ih].tup_ns    ||
          hits[ih].tdown_ns    != right.hits[ih].tdown_ns  ||
          hits[ih].incidentId_ != right.hits[ih].incidentId_ )
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
      hiter = hits.insert(hiter, *hitsrc);
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
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   zlocal = " << hiter->t_ns << " cm" << G4endl
             << "   incidentId = " << hiter->incidentId_ << G4endl
             << "   E(upstream end) = " << hiter->Eup_GeV << " GeV" << G4endl
             << "   E(downstream end) = " << hiter->Edown_GeV << " GeV" << G4endl
             << "   t(upstream end) = " << hiter->tup_ns << " ns" << G4endl
             << "   t(downstream end) = " << hiter->tdown_ns << " ns" << G4endl
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

/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       GlueX and KLF Collaborations                    * 
*                                                                         *                                                                                                                               
* Author: The GlueX and KLF Collaborations                                *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/

#include "GlueXHitCGEMhole.hh"

G4ThreadLocal G4Allocator<GlueXHitCGEMhole>* GlueXHitCGEMholeAllocator = 0;

GlueXHitCGEMhole::GlueXHitCGEMhole(G4int layer, G4int hole)
 : G4VHit(),
   layer_(layer),
   hole_(hole)
{}

GlueXHitCGEMhole::GlueXHitCGEMhole(const GlueXHitCGEMhole &src)
{
   layer_ = src.layer_;
   hole_ = src.hole_;
   hits = src.hits;
}

int GlueXHitCGEMhole::operator==(const GlueXHitCGEMhole &right) const
{
   if (layer_ !=  right.layer_ || hole_ != right.hole_)
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].dE_keV   != right.hits[ih].dE_keV  ||
          hits[ih].t_ns     != right.hits[ih].t_ns    ||
          hits[ih].dx_cm    != right.hits[ih].dx_cm   ||
          hits[ih].itrack_  != right.hits[ih].itrack_ ||
          hits[ih].z0_cm    != right.hits[ih].z0_cm)
      {
         return 0;
      }
   }
   return 1;
}

GlueXHitCGEMhole &GlueXHitCGEMhole::operator+=(const GlueXHitCGEMhole &right)
{
   if (layer_ !=  right.layer_ || hole_ != right.hole_) {
      G4cerr << "Error in GlueXHitCGEMhole::operator+=() - "
             << "illegal attempt to merge hits from two different holes!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitCGEMhole::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitCGEMhole::hitinfo_t>::const_iterator hitsrc;
   for (hitsrc = right.hits.begin(); hitsrc != right.hits.end(); ++hitsrc) {
      for (; hiter != hits.end(); ++hiter) {
         if (hiter->t_ns > hitsrc->t_ns)
            break;
      }
      hiter = hits.insert(hiter, *hitsrc);
   }
   return *this;
}

void GlueXHitCGEMhole::Draw() const
{
   // not yet implemented
}

void GlueXHitCGEMhole::Print() const
{
   G4cout << "GlueXHitCGEMhole: "
          << "   layer = " << layer_ << ", hole = " << hole_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   dE = " << hiter->dE_keV << " keV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   dx = " << hiter->dx_cm << " cm" << G4endl
             << "   itrack = " << hiter->itrack_ << G4endl
             << "   z0 = " << hiter->z0_cm << " cm" << G4endl
             << G4endl;
   }
}

void printallhits(GlueXHitsMapCGEMhole *hitsmap)
{
   std::map<int, GlueXHitCGEMhole*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCGEMhole*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

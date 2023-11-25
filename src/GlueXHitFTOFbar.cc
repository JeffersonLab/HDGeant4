//
// GlueXHitFTOFbar - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 26, 2016

#include "GlueXHitFTOFbar.hh"

G4ThreadLocal G4Allocator<GlueXHitFTOFbar>* GlueXHitFTOFbarAllocator = 0;

GlueXHitFTOFbar::GlueXHitFTOFbar(G4int plane, G4int bar)
 : G4VHit(),
   plane_(plane),
   bar_(bar),
   overflow_(0)
{}

GlueXHitFTOFbar::GlueXHitFTOFbar(const GlueXHitFTOFbar &src)
{
   plane_ = src.plane_;
   bar_ = src.bar_;
   overflow_ = src.overflow_;
   hits = src.hits;
}

int GlueXHitFTOFbar::operator==(const GlueXHitFTOFbar &right) const
{
   if (plane_ !=  right.plane_ || bar_ != right.bar_ )
      return 0;
   else if (hits.size() != right.hits.size())
      return 0;

   for (int ih=0; ih < (int)hits.size(); ++ih) {
      if (hits[ih].end_     != right.hits[ih].end_    ||
          hits[ih].dE_GeV   != right.hits[ih].dE_GeV  ||
          hits[ih].t_ns     != right.hits[ih].t_ns    ||
          hits[ih].extra.size() != right.hits[ih].extra.size())
      {
         return 0;
      }
      std::vector<hitextra_t>::const_iterator leftextra;
      std::vector<hitextra_t>::const_iterator rightextra;
      for (leftextra = hits[ih].extra.begin();
           leftextra != hits[ih].extra.end();
           ++leftextra, ++rightextra)
      {
         if (leftextra->track_   != rightextra->track_  ||
             leftextra->itrack_  != rightextra->itrack_ ||
             leftextra->ptype_G3 != rightextra->ptype_G3 ||
             leftextra->px_GeV   != rightextra->px_GeV  ||
             leftextra->py_GeV   != rightextra->py_GeV  ||
             leftextra->pz_GeV   != rightextra->pz_GeV  ||
             leftextra->E_GeV    != rightextra->E_GeV   ||
             leftextra->x_cm     != rightextra->x_cm    ||
             leftextra->y_cm     != rightextra->y_cm    ||
             leftextra->z_cm     != rightextra->z_cm    ||
             leftextra->t_ns     != rightextra->t_ns    ||
             leftextra->dist_cm  != rightextra->dist_cm )
         {
            return 0;
         }
      }
   }
   return 1;
}

GlueXHitFTOFbar &GlueXHitFTOFbar::operator+=(const GlueXHitFTOFbar &right)
{
   if (plane_ !=  right.plane_ || bar_ != right.bar_) {
      G4cerr << "Error in GlueXHitFTOFbar::operator+=() - "
             << "illegal attempt to merge hits from two different bars!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitFTOFbar::hitinfo_t>::iterator hiter = hits.begin();
   std::vector<GlueXHitFTOFbar::hitinfo_t>::const_iterator hitsrc;
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

void GlueXHitFTOFbar::Draw() const
{
   // not yet implemented
}

void GlueXHitFTOFbar::Print() const
{
   G4cout << "GlueXHitFTOFbar: "
          << "   plane = " << plane_ 
          << ",  bar = " << bar_
          << ",  overflow = " << overflow_ << G4endl;
   std::vector<hitinfo_t>::const_iterator hiter;
   for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
      G4cout << "   end = " << hiter->end_ << G4endl
             << "   dE = " << hiter->dE_GeV << " GeV" << G4endl
             << "   t = " << hiter->t_ns << " ns" << G4endl
             << "   extra hit truth information:" << G4endl;
      std::vector<hitextra_t>::const_iterator xiter;
      for (xiter = hiter->extra.begin(); xiter != hiter->extra.end(); ++xiter) {
         G4cout << "      track = " << xiter->track_ << G4endl
                << "      itrack = " << xiter->itrack_ << G4endl
                << "      ptype = " << xiter->ptype_G3 << G4endl
                << "      px = " << xiter->px_GeV << " GeV/c" << G4endl
                << "      py = " << xiter->py_GeV << " GeV/c" << G4endl
                << "      pz = " << xiter->pz_GeV << " GeV/c" << G4endl
                << "      E = " << xiter->E_GeV << " GeV" << G4endl
                << "      x = " << xiter->x_cm << " cm" << G4endl
                << "      y = " << xiter->y_cm << " cm" << G4endl
                << "      z = " << xiter->z_cm << " cm" << G4endl
                << "      t = " << xiter->t_ns << " cm" << G4endl
                << "      dist = " << xiter->dist_cm << " cm" << G4endl;
      }
      G4cout << G4endl;
   }
}

void printallhits(GlueXHitsMapFTOFbar *hitsmap)
{
   std::map<int, GlueXHitFTOFbar*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitFTOFbar*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

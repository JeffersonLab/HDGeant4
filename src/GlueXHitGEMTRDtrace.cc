//
// GlueXHitGEMTRDtrace - class implementation
//
// author: staylor at jlab.org
// version: february 7, 2022

#include "GlueXHitGEMTRDtrace.hh"

G4ThreadLocal G4Allocator<GlueXHitGEMTRDtrace>* GlueXHitGEMTRDtraceAllocator = 0;

GlueXHitGEMTRDtrace::GlueXHitGEMTRDtrace(G4int layer)
 : G4VHit(),
   layer_(layer)
{}

GlueXHitGEMTRDtrace::GlueXHitGEMTRDtrace(const GlueXHitGEMTRDtrace &src)
{
   layer_ = src.layer_;
   traces = src.traces;
}

int GlueXHitGEMTRDtrace::operator==(const GlueXHitGEMTRDtrace &right) const
{
   if (layer_ !=  right.layer_)
      return 0;
   else if (traces.size() != right.traces.size())
      return 0;

   for (int ih=0; ih < (int)traces.size(); ++ih) {
      if (traces[ih].dE_keV   != right.traces[ih].dE_keV  ||
          traces[ih].t_ns     != right.traces[ih].t_ns    ||
          traces[ih].x_cm     != right.traces[ih].x_cm    ||
	  traces[ih].y_cm     != right.traces[ih].y_cm    ||
	  traces[ih].dxdz     != right.traces[ih].dxdz    ||
	  traces[ih].dydz     != right.traces[ih].dydz    ||
          traces[ih].itrack_  != right.traces[ih].itrack_ 
	  )
	{
         return 0;
      }
   }
   return 1;
}

GlueXHitGEMTRDtrace &GlueXHitGEMTRDtrace::operator+=(const GlueXHitGEMTRDtrace &right)
{
   if (layer_ !=  right.layer_) {
      G4cerr << "Error in GlueXHitGEMTRDtrace::operator+=() - "
             << "illegal attempt to merge traces from two different layers!"
             << G4endl;
      return *this;
   }
   std::vector<GlueXHitGEMTRDtrace::traceinfo_t>::iterator titer = traces.begin();
   std::vector<GlueXHitGEMTRDtrace::traceinfo_t>::const_iterator tracesrc;
   for (tracesrc = right.traces.begin(); tracesrc != right.traces.end(); ++tracesrc) {
      for (; titer != traces.end(); ++titer) {
         if (titer->t_ns > tracesrc->t_ns)
            break;
      }
      titer = traces.insert(titer, *tracesrc);
   }
   return *this;
}

void GlueXHitGEMTRDtrace::Draw() const
{
   // not yet implemented
}

void GlueXHitGEMTRDtrace::Print() const
{
   G4cout << "GlueXHitGEMTRDtrace: "
          << "   layer = " << layer_ << G4endl;
   std::vector<traceinfo_t>::const_iterator titer;
   for (titer = traces.begin(); titer != traces.end(); ++titer) {
      G4cout << "   dE = " << titer->dE_keV << " keV" << G4endl
             << "   t = " << titer->t_ns << " ns" << G4endl
             << "   x = " << titer->x_cm << " cm" << G4endl
	     << "   y = " << titer->y_cm << " cm" << G4endl
	     << "   dx/dz = " << titer->dxdz << G4endl
	     << "   dy/dz = " << titer->dydz << G4endl
             << "   itrack = " << titer->itrack_ << G4endl
	     << G4endl;
   }
}

void printallhits(GlueXHitsMapGEMTRDtrace *tracesmap)
{
   std::map<int, GlueXHitGEMTRDtrace*> *map = tracesmap->GetMap();
   std::map<int, GlueXHitGEMTRDtrace*>::const_iterator iter;
   G4cout << "G4THitsMap " << tracesmap->GetName() << " with " << tracesmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
      G4cout << "  key=" << iter->first << " ";
      iter->second->Print();
   }
}

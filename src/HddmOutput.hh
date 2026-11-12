//
// HddmOutput - class header
//
// author: richard.t.jones at uconn.edu
// version: september 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state, but it
// is thread-safe in that its methods can be called concurrently
// from several different threads without conflicts.

#ifndef _HDDMOUTPUT_
#define _HDDMOUTPUT_

#include "G4AutoLock.hh"
#include "G4ios.hh"

#include <HDDM/hddm_s.hpp>
#include <fstream>

class HddmOutput
{
 public:
   HddmOutput(const std::string &filename);
   ~HddmOutput();

   static void WriteOutputHDDM(hddm_s::HDDM &record);
   static int getRunNo();
   static int getEventNo();
   static int incrementEventNo();
   static void setRunNo(int runno);
   static void setEventNo(int eventno);

 protected:
   HddmOutput(HddmOutput &src) = delete;
   HddmOutput& operator=(HddmOutput &src) = delete;

   static int fRunNo;
   static int fEventNo;
   static std::ofstream *fHDDMoutfile;
   static hddm_s::ostream *fHDDMostream;

 private:
   static int instanceCount;
   static G4Mutex fMutex;
};

inline int HddmOutput::getRunNo()
{
   return fRunNo;
}

inline int HddmOutput::getEventNo()
{
   return fEventNo;
}

inline void HddmOutput::setEventNo(int eventno)
{
   G4AutoLock barrier(&fMutex);
   fEventNo = eventno;
}

#endif

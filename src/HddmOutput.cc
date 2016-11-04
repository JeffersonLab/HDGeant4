//
// class implementation for HddmOutput
//
// author: richard.t.jones at uconn.edu
// version: september 28, 2016

#include "HddmOutput.hh"

int HddmOutput::fRunNo = 0;
int HddmOutput::fEventNo = 0;
std::ofstream *HddmOutput::fHDDMoutfile = 0;
hddm_s::ostream *HddmOutput::fHDDMostream = 0;

GlueXPseudoDetectorTAG HddmOutput::fTagger(0);

int HddmOutput::instanceCount = 0;
G4Mutex HddmOutput::fMutex = G4MUTEX_INITIALIZER;

HddmOutput::HddmOutput(const std::string &filename)
{
   G4AutoLock barrier(&fMutex);
   if (instanceCount++ > 0) {
      G4cerr << "Error in HddmOutput constructor - "
             << "attempt to declare than one HddmOutput object at a time, "
             << "cannot continue." << G4endl;
      exit(-1);
   }

   if (fHDDMostream != 0) {
      delete fHDDMostream;
      delete fHDDMoutfile;
   }
   fHDDMoutfile = new std::ofstream(filename);
   fHDDMostream = new hddm_s::ostream(*fHDDMoutfile);
}

HddmOutput::~HddmOutput()
{
   if (fHDDMostream != 0) {
      delete fHDDMostream;
      delete fHDDMoutfile;
      fHDDMostream = 0;
      fHDDMoutfile = 0;
   }
   --instanceCount;
}

HddmOutput::HddmOutput(HddmOutput &src)
{}

HddmOutput& HddmOutput::operator=(HddmOutput &src)
{
   return *this;
}

void HddmOutput::setRunNo(int runno)
{
   fRunNo = runno;
   fTagger.setRunNo(runno);
}

void HddmOutput::WriteOutputHDDM(hddm_s::HDDM &record)
{
   if (fHDDMostream != 0) {
      *fHDDMostream << record;
   }
}

int HddmOutput::incrementEventNo()
{
   G4AutoLock barrier(&fMutex);
   return ++fEventNo;
}


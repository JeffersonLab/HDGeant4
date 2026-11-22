//
// GlueXTimer class implementation
//
// author: richard.t.jones at uconn.edu
// version: january 30, 2017
 
#include "GlueXTimer.hh"
#include "G4ios.hh"

#include <sstream>

G4Mutex GlueXTimer::fMutex = G4MUTEX_INITIALIZER;
std::map<std::string, GlueXTimer*> GlueXTimer::fInstance;

GlueXTimer::GlueXTimer(std::string name)
{
   if (!sysconf(_POSIX_THREAD_CPUTIME)) {
      G4cerr << "GlueXTimer constructor error - this unix kernel "
             << "does not support POSIX_THREAD_CPUTIME, "
             << "cannot continue." << G4endl;
      exit(1);
   }

   for (int i=0; i < MAX_THREADS; ++i) {
      fIsStarted[i] = 0;
      fIsRunning[i] = 0;
      fRunningTotal[i] = 0;
   }

   G4AutoLock barrier(&fMutex);
   std::stringstream unique_name;
   unique_name << name;
   int n=0;
   while (fInstance.find(unique_name.str()) != fInstance.end()) {
      unique_name.str("");
      unique_name << name << "(" << ++n << ")";
   }
   fName = unique_name.str();
   fInstance[fName] = this;
}

GlueXTimer::~GlueXTimer()
{
   G4AutoLock barrier(&fMutex);
   fInstance.erase(fName);
}

void GlueXTimer::Reset() {
   int thread = G4Threading::G4GetThreadId() + 1;
   if (fIsRunning[thread]) {
#ifndef CLOCK_THREAD_CPUTIME_ID
      int err = -1;
#else
      int err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &fClock[thread]);
#endif
      assert(err == 0);
      (void)err;
      fTimeLastGo[thread] = fClock[thread].tv_sec + fClock[thread].tv_nsec*1e-9;
   }
   fRunningTotal[thread] = 0;
}

void GlueXTimer::Start() {
   int thread = G4Threading::G4GetThreadId() + 1;
   if (fIsStarted[thread]) {
      G4cerr << "GlueXTimer::Start warning - this timer is already started "
             << "in this thread, start request is ignored."
             << G4endl;
   }
   else {
#ifndef CLOCK_THREAD_CPUTIME_ID
      int err = -1;
#else
      int err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &fClock[thread]);
#endif
      assert(err == 0);
      (void)err;
      fTimeLastGo[thread] = fClock[thread].tv_sec + fClock[thread].tv_nsec*1e-9;
      fRunningTotal[thread] = 0;
      fIsStarted[thread] = 1;
      fIsRunning[thread] = 1;
   }
}

void GlueXTimer::Stop() {
   int thread = G4Threading::G4GetThreadId() + 1;
   if (fIsRunning[thread]) {
#ifndef CLOCK_THREAD_CPUTIME_ID
      int err = -1;
#else
      int err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &fClock[thread]);
#endif
      assert(err == 0);
      (void)err;
      double now = fClock[thread].tv_sec + fClock[thread].tv_nsec*1e-9;
      fRunningTotal[thread] += now - fTimeLastGo[thread];
      fIsRunning[thread] = 0;
      fIsStarted[thread] = 0;
   }
   else if (fIsStarted[thread]) {
      fIsStarted[thread] = 0;
   }
   else {
      G4cerr << "GlueXTimer::Stop warning - this timer not started "
             << "in this thread, stop request is ignored."
             << G4endl;
   }
}

void GlueXTimer::Suspend() {
   int thread = G4Threading::G4GetThreadId() + 1;
   if (fIsRunning[thread]) {
#ifndef CLOCK_THREAD_CPUTIME_ID
      int err = -1;
#else
      int err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &fClock[thread]);
#endif
      assert(err == 0);
      (void)err;
      double now = fClock[thread].tv_sec + fClock[thread].tv_nsec*1e-9;
      fRunningTotal[thread] += now - fTimeLastGo[thread];
      fIsRunning[thread] = 0;
   }
   else {
      G4cerr << "GlueXTimer::Suspend warning - this timer is not running "
             << "in this thread, suspend request is ignored."
             << G4endl;
   }
}

void GlueXTimer::Resume() {
   int thread = G4Threading::G4GetThreadId() + 1;
   if (fIsStarted[thread] && !fIsRunning[thread]) {
#ifndef CLOCK_THREAD_CPUTIME_ID
      int err = -1;
#else
      int err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &fClock[thread]);
#endif
      assert(err == 0);
      (void)err;
      fTimeLastGo[thread] = fClock[thread].tv_sec + fClock[thread].tv_nsec*1e-9;
      fIsRunning[thread] = 1;
   }
   else {
      G4cerr << "GlueXTimer::Resume warning - this timer not suspended "
             << "in this thread, stop request is ignored."
             << G4endl;
   }
}

G4double GlueXTimer::GetProcessTotal() const
{
   G4AutoLock barrier(&fMutex);
   double cputotal = 0;
   for (int i=0; i < MAX_THREADS; ++i) {
      cputotal += fRunningTotal[i];
   }
   return cputotal;
}

void GlueXTimer::PrintAll()
{
   G4AutoLock barrier(&fMutex);
   if (fInstance.size() > 0) {
      G4cout << "GlueXTimer report for all timers:" << G4endl;
      std::map<std::string, GlueXTimer*>::iterator iter;
      for (iter = fInstance.begin(); iter != fInstance.end(); ++iter) {
         double cputotal = 0;
         int threadtotal = 0;
         for (int i=0; i < MAX_THREADS; i++) {
            if (iter->second->fRunningTotal[i] > 0) {
               cputotal += iter->second->fRunningTotal[i];
               threadtotal += 1;
            }
         }
         G4cout << "   " << iter->first << ": "
                << cputotal << " seconds total from "
		        << threadtotal << " threads." << G4endl;
      }
   }
}

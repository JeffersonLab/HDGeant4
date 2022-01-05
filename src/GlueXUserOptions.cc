//
// GlueXUserOptions - class implementation
//
// author: richard.t.jones at uconn.edu
// version: may 28, 2015
//

#define APP_NAME "GlueXUserOptions"

#include <GlueXUserOptions.hh>
#include <G4ios.hh>
#include <unistd.h>
#include <iostream>
#include <fstream>

G4Mutex GlueXUserOptions::fMutex = G4MUTEX_INITIALIZER;
std::list<GlueXUserOptions*> GlueXUserOptions::fInstance;

GlueXUserOptions::GlueXUserOptions()
{
   // There is nothing wrong with constructing multiple instances
   // of GlueXUserOptions, but the normal usage would be as follows.
   //    1) main() creates GlueXUserOptions at startup and holds
   //       a reference to it until the application exits;
   //    2) main() loads the user options from input sources, eg.
   //       through a call to ReadControl_in();
   //    3) application library components who need to query the
   //       GlueXUserOptions object for configuration settings should
   //       get a pointer to the global object through a call to static
   //       method GetInstance() and then use that to look up values.
   // This is a thread-safe implementation of a singleton constructor.

   G4AutoLock barrier(&fMutex);
   fInstance.push_back(this);
}

GlueXUserOptions::~GlueXUserOptions()
{
   // This is a thread-safe implementation of a singleton destructor.

   G4AutoLock barrier(&fMutex);
   fInstance.remove(this);
}

GlueXUserOptions::GlueXUserOptions(const GlueXUserOptions &src)
{
   // This is a thread-safe implementation of a copy constructor

   G4AutoLock barrier(&fMutex);
   fKeyValue = src.fKeyValue;
   fInstance.push_back(this);
}

GlueXUserOptions &GlueXUserOptions::operator=(const GlueXUserOptions &src)
{
   // This is a thread-safe implementation of an assignment operator

   G4AutoLock barrier(&fMutex);
   fKeyValue = src.fKeyValue;
   return *this;
}

GlueXUserOptions *GlueXUserOptions::GetInstance()
{
   // Generally one only needs a single instance of this object
   // per process, and this static method lets any component in the
   // application obtain the primary instance, if any. If none has
   // yet been constructed, it returns zero.

   G4AutoLock barrier(&fMutex);
   if (fInstance.size() > 0)
      return *fInstance.begin();
   else
      return 0;
}

int GlueXUserOptions::ReadControl_in(const char *ctrlin)
{
   // Opens a geant3-style (FFREAD) input control file (default name
   // control.in) and stores its contents as key/value pairs.

   std::ifstream fin(ctrlin);
   if (!fin.good()) {
      G4cerr << "Error in GlueXUserOptions::ReadControl_in: "
             << "unable to open file " << ctrlin
             << G4endl;
      return 0;
   }

   G4AutoLock barrier(&fMutex);
  
   int ncards = 0;
   char card[9999];
   while (fin.getline(card, 9999).good()) {
      std::string cards(card);
      size_t strt = cards.find_first_not_of(" ");
      if (strt == cards.npos) {
         continue;
      }
      size_t stop = cards.substr(strt).find_first_of(" ");
      if (stop == cards.npos) {
         std::string key(askey(cards.substr(strt)));
         fKeyValue[key] = ".TRUE.";
         continue;
      }
      else {
         std::string key(askey(cards.substr(strt, stop)));
         if (key == "c" || key == "C") {
            continue;
         }
         size_t args = strt + stop;
         size_t argo = cards.substr(args).find_first_not_of(" ");
         if (argo == cards.npos) {
            fKeyValue[key] = ".TRUE.";
         }
         else {
            fKeyValue[key] = cards.substr(args + argo);
         }
      }
      ++ncards;
   }
   return ncards;
}

int GlueXUserOptions::Find(const char *name, 
                           std::map<int, std::string> &value) const
{
   // Look up name in the options key/value table and, if found, break
   // up the value string into tokens and return them as a map from
   // field number to string value. Return value is 1 (found) or 0
   // (not found).

   std::map<std::string, std::string>::const_iterator item =
                                       fKeyValue.find(askey(name));
   if (item == fKeyValue.end()) {
      return 0;
   }

   std::string args(item->second);
   value.clear();
   size_t p=0;
   int narg=1;
   while (p != args.npos) {
      size_t pfin;
      pfin = args.substr(p).find_first_not_of(" ");
      if (pfin == args.npos) {
         break;
      }
      p += pfin;
      pfin = args.substr(p).find_first_not_of("0123456789");
      if (pfin != args.npos && pfin > 0 && args[p + pfin] == '=') {
         std::string substr(args.substr(p, pfin));
         narg = atoi(substr.c_str());
         p += pfin + 1;
         pfin = args.substr(p).find_first_not_of("0123456789");
      }
      int nrep=1;
      if (pfin != args.npos && pfin > 0 && args[p + pfin] == '*') {
         std::string substr(args.substr(p, pfin));
         nrep = atoi(substr.c_str());
         p += pfin + 1;
      }
      if (args[p] == '\'') {
         pfin = args.substr(p + 1).find_first_of("'");
         for (int rep=0; rep < nrep; ++rep, ++narg) {
            value[narg] = args.substr(p + 1, pfin);
         }
      }
      else {
         pfin = args.substr(p).find_first_of(" ");
         for (int rep=0; rep < nrep; ++rep, ++narg) {
            value[narg] = args.substr(p, pfin);
         }
      }
      if (pfin == args.npos) {
         break;
      }
      else {
         p += pfin + 1;
      }
   }
   return 1;
}

int GlueXUserOptions::Find(const char *name, 
                           std::map<int, double> &value) const
{
   // Look up name in the options key/value table and, if found, break
   // up the value into tokens and return them as a map from field number
   // to double. If the field does not contain a textual representation
   // of double then it is not saved in the map. Return value is 1 (found)
   // or 0 (not found).

   std::map<std::string, std::string>::const_iterator item =
                                       fKeyValue.find(askey(name));
   if (item == fKeyValue.end()) {
      return 0;
   }

   std::string args(item->second);
   value.clear();
   size_t p=0;
   int narg=1;
   while (p != args.npos) {
      size_t pfin;
      pfin = args.substr(p).find_first_not_of(" ");
      if (pfin == args.npos) {
         break;
      }
      p += pfin;
      pfin = args.substr(p).find_first_not_of("0123456789");
      if (pfin != args.npos && pfin > 0 && args[p + pfin] == '=') {
         std::string substr(args.substr(p, pfin));
         narg = atoi(substr.c_str());
         p = pfin + 1;
         pfin = args.substr(p).find_first_not_of("0123456789");
      }
      int nrep=1;
      if (pfin != args.npos && pfin > 0 && args[p + pfin] == '*') {
         std::string substr(args.substr(p, pfin));
         nrep = atoi(substr.c_str());
         p += pfin + 1;
      }
      char *end;
      std::string subarg(args.substr(p));
      const char *start = subarg.c_str();
      double dval = strtod(start, &end);
      if (end > start) {
         for (int rep=0; rep < nrep; ++rep, ++narg) {
            value[narg] = dval;
         }
         p += end - start;
      }
      else {
         pfin = args.substr(p).find_first_of(" ");
         p = (pfin == args.npos)? pfin : p + pfin;
      }
   }
   return 1;
}

int GlueXUserOptions::Find(const char *name,
                           std::map<int, int> &value) const
{
   // Look up name in the options key/value table and, if found, break
   // up the value into tokens and return them as a map from field number
   // to int. If the field does not contain a textual representation
   // of int then it is not saved in the map. Return value is 1 (found)
   // or 0 (not found).

   std::map<std::string, std::string>::const_iterator item =
                                       fKeyValue.find(askey(name));
   if (item == fKeyValue.end()) {
      return 0;
   }

   std::string args(item->second);
   value.clear();
   size_t p=0;
   int narg=1;
   while (p != args.npos) {
      size_t pfin;
      pfin = args.substr(p).find_first_not_of(" ");
      if (pfin == args.npos) {
         break;
      }
      p += pfin;
      pfin = args.substr(p).find_first_not_of("0123456789");
      if (pfin != args.npos && pfin > 0 && args[p + pfin] == '=') {
         std::string substr(args.substr(p, pfin));
         narg = atoi(substr.c_str());
         p += pfin + 1;
         pfin = args.substr(p).find_first_not_of("0123456789");
      }
      int nrep=1;
      if (pfin != args.npos && pfin > 0 && args[p + pfin] == '*') {
         std::string substr(args.substr(p, pfin));
         nrep = atoi(substr.c_str());
         p += pfin + 1;
      }
      char *end;
      std::string subarg(args.substr(p));
      const char *start = subarg.c_str();
      double ival = strtol(start, &end, 0);
      if (end != start) {
         for (int rep=0; rep < nrep; ++rep, ++narg) {
            value[narg] = ival;
         }
         p += end - start;
      }
      else {
         pfin = args.substr(p).find_first_of(" ");
         p = (pfin == args.npos)? pfin : p + pfin;
      }
   }
   return 1;
}

std::string GlueXUserOptions::askey(const std::string name) const
{
   // Do the appropriate (uppercase) character mappings to make
   // a lookup name in the options file into a unique string.

   std::string key(name);
   for (std::string::iterator p = key.begin(); p != key.end(); ++p)
       *p = std::toupper(*p);
   return key;
}

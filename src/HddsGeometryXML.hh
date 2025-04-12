//
// HddsGeometry - class header
//
// author: richard.t.jones at uconn.edu
// version: november 21, 2017
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state.

#ifndef _HDDSGEOMETRYXML_
#define _HDDSGEOMETRYXML_

#include <JANA/Geometry/JGeometryXML.h>

class HddsGeometryXML : public JGeometryXML
{
 public:
   HddsGeometryXML(std::string url, int run) : JGeometryXML(url, run, "") {}
   ~HddsGeometryXML() {}

   DOMDocument *getDocument() { return doc; }

 private:
   HddsGeometryXML(const HddsGeometryXML &src);
   HddsGeometryXML &operator=(const HddsGeometryXML &src);
};

#endif

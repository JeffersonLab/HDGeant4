//
// GlueXDetectorConstruction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "shared", ie. has no thread-local state.

#ifndef GlueXDetectorConstruction_h
#define GlueXDetectorConstruction_h 1

#include <list>

#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "G4VUserDetectorConstruction.hh"
#include <G4VUserParallelWorld.hh>
#include <G4ThreeVector.hh>
#include <HddsG4Builder.hh>
#include <HddsGeometryXML.hh>
#include <GlueXMagneticField.hh>

#include "G4SystemOfUnits.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class GlueXDetectorMessenger;

class GlueXDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
     GlueXDetectorConstruction(G4String hddsFile="");
     GlueXDetectorConstruction(const GlueXDetectorConstruction &src);
     GlueXDetectorConstruction &operator=(const GlueXDetectorConstruction &src);
     virtual ~GlueXDetectorConstruction();

  public:
     virtual G4VPhysicalVolume* Construct();
     virtual void ConstructSDandField();
     virtual void CloneF();

     void SetUniformField(G4double field_T);
     void SetMaxStep (G4double step_mm);     

     G4double GetUniformField(G4double unit) const {
        return fUniformField * tesla / unit;
     }
     G4double GetMaxStep(G4double unit) const {
        return fMaxStep * mm / unit;
     }

     int GetParallelWorldCount() const;
     G4String GetParallelWorldName(int paraIndex) const;
     G4LogicalVolume* GetParallelWorldVolume(int paraIndex) const;

     G4ThreeVector GetMagneticField(G4ThreeVector pos, double unit) const;

     static const GlueXDetectorConstruction* GetInstance();
     static const HddsG4Builder* GetBuilder();

  private:
     G4double fMaxStep;		// maximum step size for tracking
     G4double fUniformField;  	// optional uniform field (for testing)
     G4MagneticField* fpMagneticField; // pointer to the field manager
     GlueXDetectorMessenger* fpDetectorMessenger;  // pointer to the Messenger
     HddsG4Builder fHddsBuilder; // hdds translator object instance

     static G4Mutex fMutex;
     static std::list<GlueXDetectorConstruction*> fInstance;

  protected:
     HddsGeometryXML *fGeometryXML;
};

class GlueXParallelWorld : public G4VUserParallelWorld
{
 public:
   GlueXParallelWorld(G4String worldName, G4LogicalVolume *topVolume)
    : G4VUserParallelWorld(worldName), fTopVolume(topVolume)
   {}
   GlueXParallelWorld(const GlueXParallelWorld &src);
   GlueXParallelWorld &operator=(const GlueXParallelWorld &src);
   ~GlueXParallelWorld();

   virtual void Construct(); 

 private:
   G4LogicalVolume *fTopVolume;
};

#endif

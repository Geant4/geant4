///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Able.hh
// Description: CCalG4Able will provide the basic functionallity of a CCal
//              "detector" to become a Geant4 Detector.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalG4Able_h
#define CCalG4Able_h 1

#include "g4std/vector"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

#include "CCalVisualisable.hh"

#include "g4std/iostream"

//Forward declartion for the CCalG4AbleTable typedef
class CCalG4Able;

//A table to hold a list of pointers to CMS Detectors
typedef  G4std::vector<CCalG4Able*> CCalG4AbleTable;

////////////////////
//At last the class
class CCalG4Able {

  friend G4std::ostream& operator<<(G4std::ostream&, const CCalG4Able&);

public:
  //////////////////////////////////////////////////////////////////////////
  //Constructor with a name
  CCalG4Able(G4String name);
  //////////////////////////////////////////////////////////////////////////
  //Destructor
  virtual ~CCalG4Able() {delete detPhysicalVolume;}

  //////////////////////////////////////////////////////////////////////////
  //Get/Set methods

  //Method to retrieve the physical volume of the detector.
  G4VPhysicalVolume* PhysicalVolume(G4VPhysicalVolume*);

  //Set the type of a logical volume to select its vis parameters.
  void setVisType(CCalVisualisable::visType, G4LogicalVolume*);


  //Sensitivity related
  void setSensitivity(G4bool sens=true) {sensitivity=sens;}
  G4bool isSensitive() const {return sensitivity;}

  //Name
  G4String G4Name() const {return g4ableName;}
  void setName(const G4String& name) {g4ableName = name;}


  //Comparison operator needed for CCalG4AbleTable. Compares phys. volumes
  G4bool operator==(const CCalG4Able& right) const;
 
protected:
  //A method that allows to add a new CCalG4Able inside this one.
  void AddCCalG4Able(CCalG4Able*);

  //This method actually constructs the volume. Pure virtual.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*) = 0;

  //Constructs the sensitive detectors and associates them to the corresponding
  //logical volumes
  virtual void constructSensitive() {}

protected:
  //////////////////////////////////////////////////////////////////////////
  // Data Members
  G4VPhysicalVolume *detPhysicalVolume; //The G4PhysVolume or 0 if not created
  CCalG4AbleTable theG4DetectorsInside; //CCalG4Able* daughters of this det.

  G4String g4ableName;

  G4bool sensitivity;                   //true if sentive, false if not.

  //Visualisation in G4 related Variables
  CCalVisualisable visProperties;       //Visualisation properties.
  G4VisAttributes* g4VisAtt[CCalVisualisable::TotalVisTypes];

};

#endif

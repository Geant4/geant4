///////////////////////////////////////////////////////////////////////////////
// File: G4Able.hh
// Description: G4Able will provide the basic functionallity of a CMS
//              "detector" to become a Geant4 Detector.
// Date: 03/98 
// Modified: 03/98 I.G.
//           11/98 I.G. => Added sensitivity features.
//           16/12/99 I.G. RWTPtrOrderedVector ->G4RWTPtrOrderedVector
//                         Needs change to STL!!!
//           27/03/00 S.B. Moved to G4Geometry and in SCRAM
//           17/11/01 P.A. G4RWTPtrOrderedVector ->STL vector
///////////////////////////////////////////////////////////////////////////////
#ifndef G4Able_h
#define G4Able_h 1

#include <vector>
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

#include "Visualisable.hh"

#include <iostream>

//Forward declartion for the G4AbleTable typedef
class G4Able;

//A table to hold a list of pointers to CMS Detectors
typedef  vector<G4Able*> G4AbleTable;

////////////////////
//At last the class
class G4Able {

  friend ostream& operator<<(ostream&, const G4Able&);

public:
  //////////////////////////////////////////////////////////////////////////
  //Constructor with a name
  G4Able(G4String name);
  //////////////////////////////////////////////////////////////////////////
  //Destructor
  virtual ~G4Able() {delete detPhysicalVolume;}

  //////////////////////////////////////////////////////////////////////////
  //Get/Set methods

  //Method to retrieve the physical volume of the detector.
  G4VPhysicalVolume* PhysicalVolume(G4VPhysicalVolume*);

  //Set the type of a logical volume to select its vis parameters.
  void setVisType(Visualisable::visType, G4LogicalVolume*);


  //Sensitivity related
  void setSensitivity(G4bool sens=true) {sensitivity=sens;}
  G4bool isSensitive() const {return sensitivity;}

  //Name
  G4String G4Name() const {return g4ableName;}
  void setName(const G4String& name) {g4ableName = name;}


  //Comparison operator needed for G4AbleTable. Compares phys. volumes
  G4bool operator==(const G4Able& right) const;
 
protected:
  //A method that allows to add a new g4able inside this one.
  void AddG4Able(G4Able*);

  //This method actually constructs the volume. Pure virtual.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*) = 0;

  //Constructs the sensitive detectors and associates them to the corresponding
  //logical volumes
  virtual void constructSensitive() {}

protected:
  //////////////////////////////////////////////////////////////////////////
  // Data Members
  G4VPhysicalVolume *detPhysicalVolume; //The G4PhysVolume or 0 if not created
  G4AbleTable theG4DetectorsInside;     //G4Able* daughters of this det.

  G4String g4ableName;

  G4bool sensitivity;                   //true if sentive, false if not.

  //Visualisation in G4 related Variables
  Visualisable visProperties;           //Visualisation properties.
  G4VisAttributes* g4VisAtt[Visualisable::TotalVisTypes];

};

#endif

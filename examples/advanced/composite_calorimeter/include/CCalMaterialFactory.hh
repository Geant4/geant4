///////////////////////////////////////////////////////////////////////////////
// File: CCalMaterialFactory.hh
// Description: CCalMaterialFactory is a singleton class to get from a file
//              the information to build CMS Materials. A G4Material is built
//              only if asked. A parallel structure with the file information
//              is keeped in a CCalMaterial pointer vector (CCalMaterialTable).
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalMaterialFactory_h
#define CCalMaterialFactory_h 1

#include <fstream>

#include "CCalMaterial.hh"
#include "CCalAMaterial.hh"

#include "G4MaterialTable.hh"
#include "G4ElementTable.hh"

typedef vector<CCalMaterial*>  CCalMaterialTable;
typedef vector<CCalAMaterial*> CCalAMaterialTable;

class CCalMaterialFactory {
public:
  enum MatDescription {byWeight, byVolume, byAtomic};

  ~CCalMaterialFactory();
  
  static CCalMaterialFactory* getInstance(const G4String&, const G4String&);
  static CCalMaterialFactory* getInstance(const G4String&);
  static CCalMaterialFactory* getInstance();

  G4Material* findMaterial(const G4String& ) const;
  G4Element*  findElement (const G4String& ) const;

  //Adds a new CCal Material to the list and returns a G4Material.
  G4Element*  addElement  (const G4String&, const G4String&, G4double,
                           G4double, G4double );
  //Adds a new CCal Material to the list and returns a G4Material.
  G4Material* addMaterial(const G4String& nam, G4double density, 
			  G4int nconst, G4String mats[], G4double prop[],
			  MatDescription md=byWeight);

  void readElements (const G4String&);
  void readMaterials(const G4String&);

protected:
  void readElements(ifstream&);
  void readMaterials(ifstream&);

private:
  //Constructor
  CCalMaterialFactory();

  G4Material*    findG4Material   (const G4String& ) const;
  CCalMaterial*  findCCalMaterial (const G4String& ) const;
  CCalAMaterial* findCCalAMaterial(const G4String& ) const;

  //Adds a CCalMaterial to the list. Used by readMaterials and addMaterial.
  CCalMaterial* addCCalMaterial(const G4String& nam, G4double density, 
				G4int nconst, G4String mats[], G4double prop[],
				MatDescription md=byWeight);
  

private:
  static CCalMaterialFactory* instance;

  static G4String elementfile; //File with the elements.
  static G4String mixturefile; //File with the materials.

  CCalMaterialTable  theCCalMaterials;  //Where the CCalMaterials are stored
  CCalAMaterialTable theCCalAMaterials; //Where the CCalAMaterials are stored

};

#endif

///////////////////////////////////////////////////////////////////////////////
// File: CMSMaterialFactory.hh
// Date: 03/98 
// Last modification: 13/03/98 I.G.
//                    20/03/98 S.B.
//                    23/03/98 I.G.
//                    16/12/99 I.G. RWTPtrOrderedVector ->G4RWTPtrOrderedVector
//                                  Needs change to STL!!!!
//                    27/03/00 S.B. Get rid of readName (in OSCAR)
//           17/11/01 P.A. G4RWTPtrOrderedVector ->STL vector
// Description: CMSMaterialFactory is a singleton class to get from a file
//              the information to build CMS Materials. A G4Material is built
//              only if asked. A parallel structure with the file information
//              is keeped in a CMSMaterial pointer vector (CMSMaterialTable).
///////////////////////////////////////////////////////////////////////////////
#ifndef CMSMaterialFactory_h
#define CMSMaterialFactory_h 1

#include <fstream>

#include "CMSMaterial.hh"
#include "CMSAMaterial.hh"

#include "G4MaterialTable.hh"
#include "G4ElementTable.hh"

typedef vector<CMSMaterial*>  CMSMaterialTable;
typedef vector<CMSAMaterial*> CMSAMaterialTable;

class CMSMaterialFactory {
public:
  enum MatDescription {byWeight, byVolume, byAtomic};

  ~CMSMaterialFactory();
  
  static CMSMaterialFactory* getInstance(const G4String&, const G4String&);
  static CMSMaterialFactory* getInstance(const G4String&);
  static CMSMaterialFactory* getInstance();

  G4Material* findMaterial(const G4String& ) const;
  G4Element*  findElement (const G4String& ) const;

  //Adds a new CMS Material to the list and returns a G4Material.
  G4Element*  addElement  (const G4String&, const G4String&, G4double,
                           G4double, G4double );
  //Adds a new CMS Material to the list and returns a G4Material.
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
  CMSMaterialFactory();

  G4Material*   findG4Material  (const G4String& ) const;
  CMSMaterial*  findCMSMaterial (const G4String& ) const;
  CMSAMaterial* findCMSAMaterial(const G4String& ) const;

  //Adds a CMSMaterial to the list. Used by readMaterials and addMaterial.
  CMSMaterial* addCMSMaterial(const G4String& nam, G4double density, 
			      G4int nconst, G4String mats[], G4double prop[],
			      MatDescription md=byWeight);
  

private:
  static CMSMaterialFactory* instance;

  static G4String elementfile; //File with the elements.
  static G4String mixturefile; //File with the materials.

  CMSMaterialTable theCMSMaterials;  //Where the CMSMaterials are stored
  CMSAMaterialTable theCMSAMaterials; //Where the CMSAMaterials are stored

};

#endif

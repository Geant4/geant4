#ifndef G4tgbMaterialMgr_h
#define G4tgbMaterialMgr_h
using namespace std;
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbMaterialMgr
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Singleton class to manage the building of transient materials, as well as the construction of the corresponding G4Material's
 */
#include "G4tgbIsotope.hh"
#include "G4tgbElement.hh"
#include "G4tgbMaterial.hh"

#include "G4tgrIsotope.hh"
#include "G4tgrElement.hh"
#include "G4tgrElement.hh"
#include "G4tgrMaterial.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"

typedef map< G4String, G4tgbIsotope* > mstgbisot;
typedef map< G4String, G4tgbElement* > mstgbelem;
typedef map< G4String, G4tgbMaterial* > mstgbmate;
typedef map< G4String, G4Isotope* > msg4isot;
typedef map< G4String, G4Element* > msg4elem;
typedef map< G4String, G4Material* > msg4mate;

class G4tgbMaterialMgr 
{
public:
  //-  enum MatDescription {byWeight, byNoAtoms, byVolume};

  ~G4tgbMaterialMgr();
  
  //! get only instance (it it does not exists, create it)
  static G4tgbMaterialMgr* GetInstance();

  //! copy the G4tgrIsotopes into G4tgbIsotopes
  void CopyIsotopes();
  //! copy the G4tgrElements into G4tgbElements
  void CopyElements();
  //! copy the G4tgrMaterials into G4tgbMaterials
  void CopyMaterials();

  //! Look for an G4Isotope that has to exists (if not found create it from the corresponding G4tgbIsotope)
  G4Isotope* FindOrBuildG4Isotope(const G4String & name);
  //! Look for an G4Isotope and if not found return 0
  G4Isotope* FindG4Isotope(const G4String & name) const;  //! Look for an G4tgbMaterial that has to exists (if not found exit)
  //! Look for an G4Isotope and if not found return 0
  G4tgbIsotope* FindG4tgbIsotope(const G4String& name, G4bool bMustExist = 0) const;

  //! Look for an G4Element that has to exists (if not found create it from the corresponding G4tgbElement)
  G4Element* FindOrBuildG4Element(const G4String & name);
  //! Look for an G4Element and if not found return 0
  G4Element* FindG4Element(const G4String & name) const;  //! Look for an G4tgbMaterial that has to exists (if not found exit)
  //! Look for an G4Element and if not found return 0
  G4tgbElement* FindG4tgbElement(const G4String& name, G4bool bMustExist = 0) const;

  //! Look for an G4Material that has to exists (if not found create it from the corresponding G4tgbMaterial)
  G4Material* FindOrBuildG4Material(const G4String& name);
  //! Look for an G4Material and if not found return 0
  G4Material* FindG4Material(const G4String& name) const;
  //! Look for an G4tgbMaterial and if not found return 0
  G4tgbMaterial* FindG4tgbMaterial(const G4String& name, G4bool bMustExist = 0) const;

 public:
  const msg4isot GetG4IsotopeList() const { 
    return theG4Isotopes;
  }
  const msg4elem GetG4ElementList() const { 
    return theG4Elements;
  }
  const msg4mate GetG4MaterialList() const {
    return theG4Materials;
  }

 private:
  //Constructor
  G4tgbMaterialMgr(){};
  
 private:
  static G4tgbMaterialMgr* theInstance;

  //! List of all tgbIsotopes created
  mstgbisot theG4tgbIsotopes;
  //! List of all tgbElements created
  mstgbelem theG4tgbElements;
  //! List of all G4tgbMaterials created
  mstgbmate theG4tgbMaterials;
  //! container of all G4Isotopes created
  msg4isot theG4Isotopes;  
  //! container of all G4Elements created
  msg4elem theG4Elements;  
  //! container of all G4Materials created
  msg4mate theG4Materials;  

};

#endif


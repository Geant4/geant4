#ifndef G4tgrMaterialFactory_h
#define G4tgrMaterialFactory_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrMaterialFactory
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Singleton class to manage the building of transient materials,
 */
#include "G4tgrIsotope.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrMaterial.hh"
#include "G4tgrMaterialSimple.hh"
#include "G4tgrMaterialMixture.hh"
#include <map>

typedef map< G4String, G4tgrIsotope* > mstgrisot;
typedef map< G4String, G4tgrElement* > mstgrelem;
typedef map< G4String, G4tgrMaterial* > mstgrmate;

class G4tgrMaterialFactory 
{
public:
  ~G4tgrMaterialFactory();
  
  //! get only instance (it it does not exists, create it)
  static G4tgrMaterialFactory* GetInstance();

  //! Build a G4tgrIsotope
  G4tgrIsotope* AddIsotope( const vector<G4String>& wl );

  //! Build a G4tgrElementSimple
  G4tgrElementSimple* AddElementSimple( const vector<G4String>& wl );
  //! Build a G4tgrElementFromIsotopes
  G4tgrElementFromIsotopes* AddElementFromIsotopes( const vector<G4String>& wl );
                                    
  //! Build a G4tgrMaterialSimple and add it to theMaterials
  G4tgrMaterialSimple* AddMaterialSimple( const vector<G4String>& wl );

  //! Build a G4tgrMaterialByWeight or G4tgrMaterialByNoAtoms or G4tgrMaterialByVolume and add it to theMaterials
  G4tgrMaterialMixture* AddMaterialMixture( const vector<G4String>& wl, const G4String& mixtType );

  //! Look for a G4tgrIsotope and if not found return 0
  G4tgrIsotope* FindIsotope(const G4String& name) const;

  //! Look for an G4tgrElement and if not found return 0
  G4tgrElement* FindElement(const G4String& name) const;

  //! Look for an G4tgrMaterial and if not found return 0
  G4tgrMaterial* FindMaterial(const G4String& name) const;

  //! dump detailed list of isotopes
  void DumpIsotopeList() const;
  //! dump detailed list of elements
  void DumpElementList() const;
  //! dump detailed list of materials
  void DumpMaterialList() const;

 public:
  const mstgrisot& GetIsotopeList() const {return theG4tgrIsotopes;}
  const mstgrelem& GetElementList() const {return theG4tgrElements;}
  const mstgrmate& GetMaterialList() const {return theG4tgrMaterials;}

 private:
  //Constructor
  G4tgrMaterialFactory(){};

  void ErrorAlreadyExists(const G4String& object,  const vector<G4String>& wl, const G4bool bNoRepeating = TRUE );

 private:
  static G4tgrMaterialFactory* theInstance;

  //! List of all G4tgrIsotopes created
  mstgrisot theG4tgrIsotopes;
  //! List of all G4tgrElements created
  mstgrelem theG4tgrElements;
  //! List of all G4tgrMaterials created
  mstgrmate theG4tgrMaterials;

};

#endif


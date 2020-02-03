//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalMaterialFactory.cc
// Description: CCalMaterialFactory is a factory class to vuild G4Material 
//              from CCalMaterial and CCalAmaterial
///////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <stdlib.h>

#include "CCalMaterialFactory.hh"
#include "CCalutils.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"

//#define ddebug
//#define debug

typedef CCalMaterial*  ptrCCalMaterial;
typedef CCalAMaterial* ptrCCalAMaterial;


CCalMaterialFactory * CCalMaterialFactory::instance = 0;
G4String CCalMaterialFactory::elementfile = "";
G4String CCalMaterialFactory::mixturefile = "";


CCalMaterialFactory* CCalMaterialFactory::getInstance(const G4String& matfile,
                                                      const G4String& mixfile){
  if ((matfile=="" || matfile==elementfile) &&
      (mixfile=="" || mixfile==mixturefile))
    return getInstance();
  else if ((matfile != "" && elementfile != "" && matfile != elementfile) ||
           (mixfile != "" && mixturefile != "" && mixfile != mixturefile)) {
    G4cerr << "ERROR: Trying to get materials from " << matfile << " and " 
         << mixfile << " while previously were retrieved from " 
         << elementfile << " and " << mixturefile << "." << G4endl;
    return 0;
  } else {
    if (elementfile == "") 
      elementfile=matfile;
    if (mixturefile == "") 
      mixturefile=mixfile;
    return getInstance();
  }
}


CCalMaterialFactory* CCalMaterialFactory::getInstance(const G4String& matfile){
  return getInstance(matfile,matfile);
}


CCalMaterialFactory* CCalMaterialFactory::getInstance(){
  if (elementfile=="" || mixturefile=="") {
    G4cerr << "ERROR: You haven't defined files to be used for materials in "
         << "CCalMaterialFactory::getInstance(const G4String&,const G4String&)"
         << G4endl;
    return 0;
  }

  if (instance==0) {
    instance = new CCalMaterialFactory;
    return instance;
  }
  else
    return instance;
}


CCalMaterialFactory::~CCalMaterialFactory(){
  CCalMaterialTable::iterator ite;
  for(ite = theCCalMaterials.begin(); ite != theCCalMaterials.end(); ite++ ){
    delete *ite;
  }
  theCCalMaterials.clear();
  CCalAMaterialTable::iterator itea;
  for(itea = theCCalAMaterials.begin(); itea != theCCalAMaterials.end(); 
      itea++ ){
    delete *itea;
  }
  theCCalAMaterials.clear();
}

  
G4Material* CCalMaterialFactory::findMaterial(const G4String & mat) const {
  G4Material* theMat=findG4Material(mat);

  if (theMat) {
#ifdef ddebug
    G4cout << "Material " << mat << " already defined. Returning previous "
         << "instance." << G4endl;
#endif
    return theMat;
  } else { 
    CCalMaterial* CCalmat=findCCalMaterial(mat);
    if (CCalmat){
      G4Material* G4Mat = new G4Material(CCalmat->Name(),
                                         CCalmat->Density()*g/cm3, 
                                         CCalmat->NElements());
      for(G4int i=0; i<CCalmat->NElements(); i++) {
        G4Element* elem = findElement(CCalmat->Element(i));
        if (!elem) {
          G4ExceptionDescription ed;
          ed << "   Could not build material " << mat << "." << G4endl;
          G4Exception("CCalMaterialFactory::findMaterial()","ccal001",
                      FatalException,ed);
        }
        G4Mat->AddElement(elem, CCalmat->Weight(i));
      }
#ifdef ddebug
    G4cout << "Material " << mat << " has been built successfully." << G4endl;
#endif
      return G4Mat;
    } else {
      G4cerr << "ERROR: Material " << mat << " not found in CCal database!!!" 
           << G4endl;
      return 0;
    }
  }
}


G4Element* CCalMaterialFactory::findElement(const G4String & mat) const {
  const G4ElementTable  theElements = *(G4Element::GetElementTable());
  for (unsigned int i=0; i<theElements.size(); i++)
    if (theElements[i]->GetName()==mat){
#ifdef ddebug
      G4cout << "Element " << mat << " found!" << G4endl;
#endif
      return theElements[i];
    }
  return 0;
}


G4Element* CCalMaterialFactory::addElement(const G4String & name,
                                           const G4String & symbol,
                                           G4double Z, G4double A,
                                           G4double density) {

  G4Element* theEl = new G4Element(name, symbol, Z, A*g/mole);
  //Make it also as a material.
  CCalAMaterial* theMat = new CCalAMaterial(name,A,density);
  theCCalAMaterials.push_back(theMat);

#ifdef ddebug
  G4cout << "Element " << name << " created!" << G4endl;
#endif
  return theEl;
}


G4Material* CCalMaterialFactory::addMaterial(const G4String& name,
                                             G4double density,
                                             G4int nconst,
                                             G4String mats[],
                                             G4double prop[],
                                             MatDescription md){
  addCCalMaterial(name, density, nconst, mats, prop, md);
  return findMaterial(name);
}


void CCalMaterialFactory::readElements(const G4String& matfile) {

  G4String path = "NULL";
  if (std::getenv("CCAL_GLOBALPATH"))
    path = std::getenv("CCAL_GLOBALPATH");

  G4cout << " ==> Opening file " << matfile << " to read elements..." << G4endl;
  std::ifstream is;
  G4bool ok = openGeomFile(is, path, matfile);
  if (!ok) {
    G4cerr << "ERROR: Could not open file " << matfile << G4endl;
    return;
  }

  // Find *DO GMAT
  findDO(is, G4String("GMAT"));
  
  readElements(is);
  
  is.close();
}


void CCalMaterialFactory::readMaterials(const G4String& matfile) {

  G4String path = "NULL";
  if (std::getenv("CCAL_GLOBALPATH"))
    path = std::getenv("CCAL_GLOBALPATH");

  G4cout << " ==> Opening file " << matfile << " to read materials..." << G4endl;
  std::ifstream is;
  bool ok = openGeomFile(is, path, matfile);
  if (!ok) {
    G4cerr << "ERROR: Could not open file " << matfile << G4endl;
    return;
  }

  // Find *DO GMIX
  findDO(is, G4String("GMIX"));

  readMaterials(is);

  is.close();
}


//===========================================================================
// Protected & private methods ==============================================


G4Material* CCalMaterialFactory::findG4Material(const G4String & mat) const {
  const G4MaterialTable theG4Materials = *(G4Material::GetMaterialTable());
  for (unsigned int i=0; i<theG4Materials.size(); i++) {
    if (theG4Materials[i]->GetName()==mat){
      return theG4Materials[i];
    }
  }
  return 0;
}


CCalMaterial* CCalMaterialFactory::findCCalMaterial(const G4String & mat) 
  const {
  for (unsigned int i=0; i<theCCalMaterials.size(); i++)
    if (theCCalMaterials[i]->Name()==mat){
#ifdef ddebug
      G4cout << "CCalMaterial " << mat << " found!" << G4endl;
#endif
      return theCCalMaterials[i];
    }
  return (CCalMaterial*) findCCalAMaterial(mat);
}


CCalAMaterial* CCalMaterialFactory::findCCalAMaterial(const G4String & mat) 
  const {
  for (unsigned int i=0; i<theCCalAMaterials.size(); i++)
    if (theCCalAMaterials[i]->Name()==mat){
#ifdef ddebug
      G4cout << "CCalMaterial " << mat << " found!" << G4endl;
#endif
      return theCCalAMaterials[i];
    }
  return 0;
}


CCalMaterial* CCalMaterialFactory::addCCalMaterial(const G4String& name, 
                                                   G4double density,
                                                   G4int nconst,
                                                   G4String mats[], 
                                                   G4double prop[],
                                                   MatDescription md){
  ptrCCalMaterial* matcol=0;
  ptrCCalAMaterial* amatcol=0;

  if (md==byAtomic)
    amatcol = new ptrCCalAMaterial[nconst];
  else
    matcol = new ptrCCalMaterial[nconst];

  for (G4int i=0; i<nconst; i++){
    if (md==byAtomic) {
      CCalAMaterial* amat = findCCalAMaterial(mats[i]);
      if (amat)
        amatcol[i]=amat;
      else {
        G4cerr << "ERROR: Trying to build" << name << " out of unknown " 
             << mats[i] << "." << G4endl 
             << "Skiping this material!" << G4endl;
        delete[] amatcol;
        return 0;
      }
    } //by Atomic fractions
    else {
      CCalMaterial* mat = findCCalMaterial(mats[i]);
      if (mat)
        matcol[i]=mat;
      else {
        G4cerr << "ERROR: Trying to build" <<name << " out of unknown " 
             << mats[i] << "." << G4endl 
             << "Skiping this material!" << G4endl;
        delete[] matcol;
        return 0;
      }
    }
  } //for

  //Let's do the CCalMaterial!
  if (md==byAtomic) {
    CCalAMaterial* amaterial = new CCalAMaterial(name, density, nconst, 
                                                 amatcol, prop);
    delete[] amatcol;
    theCCalAMaterials.push_back(amaterial);
#ifdef ddebug
    G4cout << *amaterial << G4endl;
#endif
    return amaterial;
  } else {
    CCalMaterial::FractionType ft;
    if (md == byWeight)
      ft=CCalMaterial::FTWeight;
    else
      ft=CCalMaterial::FTVolume;
    CCalMaterial* material = new CCalMaterial(name, density, nconst, 
                                            matcol, prop, ft);
    delete[] matcol;
    theCCalMaterials.push_back(material);
#ifdef ddebug
    G4cout << *material << G4endl;
#endif
    return material;
  }  
}


void CCalMaterialFactory::readElements(std::ifstream& is){
  G4String name, symbol;
  
  G4cout << "     ==> Reading elements... " << G4endl;
#ifdef debug
  G4cout << "       Element    \tsymbol\tA\tZ\tdensity\tX_0          abs_l"<< G4endl;
#endif
  //There should come the list of materials. #. Defines a comment
  //*DO defines the beguining of the Mixes block.

  readName(is,name);
  while (name != "*ENDDO") {
    //It should be an element definition
    G4double A, Z, density;
    is >> symbol >> A >> Z >> density >> jump;    
#ifdef debug
    G4cout << "       " << name << "    \t" << symbol << "\t" 
         << A << "\t" << Z << "\t" << density << G4endl;
#endif    
    addElement(name, symbol, Z, A, density);
    readName(is,name);
  };
  G4cout << "     " << G4Element::GetElementTable()->size() 
       << " elements read from file" << G4endl << G4endl;
}


void CCalMaterialFactory::readMaterials(std::ifstream& is){
  G4String name, matname;

  G4cout << "     ==> Reading materials... " << G4endl;

  //Take into account the special case of vacuum...
#ifdef debug
  G4cout <<"       \"Vacuum\"" << G4endl;
#endif
  G4double density  = universe_mean_density;   //from PhysicalConstants.h
  G4double pressure = 1.E-19*pascal;
  G4double temperature = 0.1*kelvin;
  new G4Material("Vacuum", /*Z=*/ 1., /*A=*/ 1.01*g/mole,
                 density, kStateGas, temperature, pressure);

  //There should come the list of materials. #. Defines a comment
  //*ENDDO defines the block.
  readName(is,name);
  while (name != "*ENDDO") {
    //It should be a material definition
    matname=name;
    G4int nElem;
    G4double dens;
    is >> nElem >> dens >> jump;

#ifdef debug
    G4cout <<"       " << matname
           << " made of " << nElem 
           << " elements. Density=" << dens
           << G4endl;
#endif

    G4int absnelem = std::abs(nElem);

    G4String* mats    = new G4String[absnelem];
    G4double* weights = new G4double[absnelem];
    
    G4double prop;
    for(int i=0; i<absnelem; i++) {
      readName(is, name);
      is >> prop >> jump;
      mats[i]=name;
      weights[i]=std::abs(prop);
    } //for...
    MatDescription md;
    if (nElem>0 && prop<0)
      md = byAtomic;
    else if (nElem>0)
      md = byWeight;
    else
      md = byVolume;

    addCCalMaterial(matname, dens, absnelem, mats, weights, md);
    delete[] mats;
    delete[] weights;
    
    readName(is,name);
  };  //while

  G4cout << "     " << theCCalMaterials.size() << " materials read from " 
       << mixturefile << G4endl << G4endl;
}


CCalMaterialFactory::CCalMaterialFactory() {
  readElements (elementfile);
  readMaterials(mixturefile);
}

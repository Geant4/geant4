///////////////////////////////////////////////////////////////////////////////
// File: CCalMaterialFactory.cc
// Description: CCalMaterialFactory is a factory class to vuild G4Material 
//              from CCalMaterial and CCalAmaterial
///////////////////////////////////////////////////////////////////////////////
#include "CCalMaterialFactory.hh"
#include "CCalutils.hh"
#include <fstream.h>
#include <stdlib.h>

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
    cerr << "ERROR: Trying to get materials from " << matfile << " and " 
         << mixfile << " while previously were retrieved from " 
         << elementfile << " and " << mixturefile << "." << endl;
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
    cerr << "ERROR: You haven't defined files to be used for materials in "
	 << "CCalMaterialFactory::getInstance(const G4String&,const G4String&)"
	 << endl;
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
    cout << "Material " << mat << " already defined. Returning previous "
	 << "instance." << endl;
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
	  cerr << "       Could not build material " << mat << "." << endl;
	  exit(-10);
	}
	G4Mat->AddElement(elem, CCalmat->Weight(i));
      }
#ifdef ddebug
    cout << "Material " << mat << " has been built successfully." << endl;
#endif
      return G4Mat;
    } else {
      cerr << "ERROR: Material " << mat << " not found in CCal database!!!" 
	   << endl;
      return 0;
    }
  }
}

G4Element* CCalMaterialFactory::findElement(const G4String & mat) const {
  const G4ElementTable  theElements = *(G4Element::GetElementTable());
  for (unsigned int i=0; i<theElements.size(); i++)
    if (theElements[i]->GetName()==mat){
#ifdef ddebug
      cout << "Element " << mat << " found!" << endl;
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
  cout << "Element " << name << " created!" << endl;
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

  G4String path = getenv("OSCARGLOBALPATH");
  cout << " ==> Opening file " << matfile << " to read elements..." << endl;
  ifstream is;
  bool ok = openGeomFile(is, path, matfile);
  if (!ok) {
    cerr << "ERROR: Could not open file " << matfile << endl;
    return;
  }

  // Find *DO GMAT
  findDO(is, G4String("GMAT"));
  
  readElements(is);
  
  is.close();
}

void CCalMaterialFactory::readMaterials(const G4String& matfile) {

  G4String path = getenv("OSCARGLOBALPATH");
  cout << " ==> Opening file " << matfile << " to read materials..." << endl;
  ifstream is;
  bool ok = openGeomFile(is, path, matfile);
  if (!ok) {
    cerr << "ERROR: Could not open file " << matfile << endl;
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
      cout << "CCalMaterial " << mat << " found!" << endl;
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
      cout << "CCalMaterial " << mat << " found!" << endl;
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
	cerr << "ERROR: Trying to build" << name << " out of unknown " 
	     << mats[i] << "." << endl 
	     << "Skiping this material!" << endl;
	delete[] amatcol;
	return 0;
      }
    } //by Atomic fractions
    else {
      CCalMaterial* mat = findCCalMaterial(mats[i]);
      if (mat)
	matcol[i]=mat;
      else {
	cerr << "ERROR: Trying to build" <<name << " out of unknown " 
	     << mats[i] << "." << endl 
	     << "Skiping this material!" << endl;
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
    cout << *amaterial << endl;
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
    cout << *material << endl;
#endif
    return material;
  }  
}


void CCalMaterialFactory::readElements(ifstream& is){
  G4String name, symbol;
  
  cout << "     ==> Reading elements... " << endl;
#ifdef debug
  cout << "       Element    \tsymbol\tA\tZ\tdensity\tX_0          abs_l"<< endl;
#endif
  //There should come the list of materials. #. Defines a comment
  //*DO defines the beguining of the Mixes block.

  readName(is,name);
  while (name != "*ENDDO") {
    //It should be an element definition
    G4double A, Z, density;
    is >> symbol >> A >> Z >> density >> jump;
    
#ifdef debug
    cout << "       " << name << "    \t" << symbol << "\t" 
	 << A << "\t" << Z << "\t" << density << endl;
#endif
    
    addElement(name, symbol, Z, A, density);
    
    readName(is,name);
  };
  cout << "     " << G4Element::GetElementTable()->size() 
       << " elements read from file" << endl << endl;
}


void CCalMaterialFactory::readMaterials(ifstream& is){
  G4String name, matname;

  cout << "     ==> Reading materials... " << endl;

  //Take into account the special case of vacuum...
#ifdef debug
  cout <<"       \"Vacuum\"" << endl;
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
    G4double density;
    is >> nElem >> density >> jump;

#ifdef debug
    cout <<"       " << matname
	 << " made of " << nElem 
	 << " elements. Density=" << density 
	 << endl;
#endif

    G4int absnelem = abs(nElem);

    G4String* mats    = new G4String[absnelem];
    G4double* weights = new G4double[absnelem];
    
    G4double prop;
    for(int i=0; i<absnelem; i++) {
      readName(is, name);
      
      is >> prop >> jump;

      mats[i]=name;
      weights[i]=abs(prop);
    } //for...
    MatDescription md;
    if (nElem>0 && prop<0)
      md = byAtomic;
    else if (nElem>0)
      md = byWeight;
    else
      md = byVolume;

    addCCalMaterial(matname, density, absnelem, mats, weights, md);
    delete[] mats;
    delete[] weights;
    
    readName(is,name);
  };  //while


  cout << "     " << theCCalMaterials.size() << " materials read from " 
       << mixturefile << endl << endl;
}

CCalMaterialFactory::CCalMaterialFactory() {
  readElements (elementfile);
  readMaterials(mixturefile);
}

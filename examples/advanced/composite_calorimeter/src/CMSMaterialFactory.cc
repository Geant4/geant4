///////////////////////////////////////////////////////////////////////////////
// File: CMSMaterialFactory.cc
// Date: 03/98 I. Gonzalez
// Modifications: 13/03/98 I.G.
//                27/03/00 S.B. Inside OSCAR
///////////////////////////////////////////////////////////////////////////////
#include "CMSMaterialFactory.hh"
#include "utils.hh"
#include <fstream.h>
#include <stdlib.h>

#include "G4Material.hh"

//#define ddebug
//#define debug

typedef CMSMaterial* ptrCMSMaterial;
typedef CMSAMaterial* ptrCMSAMaterial;


CMSMaterialFactory * CMSMaterialFactory::instance = 0;
G4String CMSMaterialFactory::elementfile = "";
G4String CMSMaterialFactory::mixturefile = "";

CMSMaterialFactory* CMSMaterialFactory::getInstance(const G4String & matfile,
                                                    const G4String & mixfile){
  if ((matfile=="" || matfile==elementfile) &&
      (mixfile=="" || mixfile==mixturefile))
    return getInstance();
  else if ((matfile != "" && elementfile != "" && matfile != elementfile) ||
           (mixfile != "" && mixturefile != "" && mixfile != mixturefile)) {
    cerr << "ERROR: Trying to get materials from " << matfile << " and " 
         << mixfile << " while previously were retrieved from " 
         << elementfile << " and " << mixturefile << "." << endl;
    return 0;
  }
  else {
    if (elementfile == "") 
      elementfile=matfile;
    if (mixturefile == "") 
      mixturefile=mixfile;
    return getInstance();
  }
}

CMSMaterialFactory* CMSMaterialFactory::getInstance(const G4String & matfile){
    return getInstance(matfile,matfile);
}

CMSMaterialFactory* CMSMaterialFactory::getInstance(){
  if (elementfile=="" || mixturefile=="") {
    cerr << "ERROR: You haven't defined files to be used for materials in CMSMaterialFactory::getInstance(const G4String&, const G4String&)" << endl;
    return 0;
  }

  if (instance==0) {
    instance = new CMSMaterialFactory;
    return instance;
  }
  else
    return instance;
}

CMSMaterialFactory::~CMSMaterialFactory(){
  CMSMaterialTable::iterator ite;
  for(ite = theCMSMaterials.begin(); ite != theCMSMaterials.end(); ite++ ){
    delete *ite;
  }
  theCMSMaterials.clear();
  CMSAMaterialTable::iterator itea;
  for(itea = theCMSAMaterials.begin(); itea != theCMSAMaterials.end(); itea++ ){
    delete *itea;
  }
  theCMSAMaterials.clear();
}
  
G4Material* CMSMaterialFactory::findMaterial(const G4String & mat) const {
  G4Material* theMat=findG4Material(mat);

  if (theMat) {
#ifdef ddebug
    cout << "Material " << mat << " already defined. Returning previous instance." << endl;
#endif
    return theMat;
  }
  else { 
    CMSMaterial* CMSmat=findCMSMaterial(mat);
    if (CMSmat){
      G4Material* G4Mat = new G4Material(CMSmat->Name(),
					 CMSmat->Density()*g/cm3, 
					 CMSmat->NElements());
      for(G4int i=0; i<CMSmat->NElements(); i++) {
	G4Element* elem = findElement(CMSmat->Element(i));
	if (!elem) {
	  cerr << "       Could not build material " << mat << "." << endl;
	  exit(-10);
	}
	G4Mat->AddElement(elem, CMSmat->Weight(i));
      }
#ifdef ddebug
    cout << "Material " << mat << " has been built successfully." << endl;
#endif
      return G4Mat;
    }
    else {
      cerr << "ERROR: Material " << mat << " not found in CMS database!!!" << endl;
      return 0;
    }
  }
}

G4Element* CMSMaterialFactory::findElement(const G4String & mat) const {
  const G4ElementTable  theElements = *(G4Element::GetElementTable());
  for (int i=0; i<theElements.size(); i++)
    if (theElements[i]->GetName()==mat){
#ifdef ddebug
      cout << "Element " << mat << " found!" << endl;
#endif
      return theElements[i];
    }
  return 0;
}

G4Element* CMSMaterialFactory::addElement(const G4String & name,
                                          const G4String & symbol,
                                          G4double Z, G4double A,
                                          G4double density) {

  G4Element* theEl = new G4Element(name, symbol, Z, A*g/mole);
  //Make it also as a material.
  CMSAMaterial* theMat = new CMSAMaterial(name,A,density);
  theCMSAMaterials.push_back(theMat);

#ifdef ddebug
    cout << "Element " << name << " created!" << endl;
#endif
  return theEl;
}

G4Material* CMSMaterialFactory::addMaterial(const G4String& name,
					     G4double density,
					     G4int nconst,
					     G4String mats[],
					     G4double prop[],
					     MatDescription md){
  addCMSMaterial(name, density, nconst, mats, prop, md);
  return findMaterial(name);
}

void CMSMaterialFactory::readElements(const G4String& matfile) {

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

void CMSMaterialFactory::readMaterials(const G4String& matfile) {

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

G4Material* CMSMaterialFactory::findG4Material(const G4String & mat) const {
  const G4MaterialTable theG4Materials = *(G4Material::GetMaterialTable());
  for (int i=0; i<theG4Materials.size(); i++) {
    if (theG4Materials[i]->GetName()==mat){
      return theG4Materials[i];
    }
  }
  return 0;
}

CMSMaterial* CMSMaterialFactory::findCMSMaterial(const G4String & mat) const {
  for (int i=0; i<theCMSMaterials.size(); i++)
    if (theCMSMaterials[i]->Name()==mat){
#ifdef ddebug
      cout << "CMSMaterial " << mat << " found!" << endl;
#endif
      return theCMSMaterials[i];
    }
  return (CMSMaterial*) findCMSAMaterial(mat);
}

CMSAMaterial* CMSMaterialFactory::findCMSAMaterial(const G4String & mat) const {
  for (int i=0; i<theCMSAMaterials.size(); i++)
    if (theCMSAMaterials[i]->Name()==mat){
#ifdef ddebug
      cout << "CMSMaterial " << mat << " found!" << endl;
#endif
      return theCMSAMaterials[i];
    }
  return 0;
}




CMSMaterial* CMSMaterialFactory::addCMSMaterial(const G4String& name, 
						G4double density,
						G4int nconst,
						G4String mats[], 
						G4double prop[],
						MatDescription md){
  ptrCMSMaterial* matcol=0;
  ptrCMSAMaterial* amatcol=0;

  if (md==byAtomic)
    amatcol = new ptrCMSAMaterial[nconst];
  else
    matcol = new ptrCMSMaterial[nconst];
    
  for (G4int i=0; i<nconst; i++){
    if (md==byAtomic) {
      CMSAMaterial* amat = findCMSAMaterial(mats[i]);
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
      CMSMaterial* mat = findCMSMaterial(mats[i]);
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

  //Let's do the CMSMaterial!
  if (md==byAtomic) {
    CMSAMaterial* amaterial = new CMSAMaterial(name, density, nconst, 
					       amatcol, prop);
    delete[] amatcol;
    theCMSAMaterials.push_back(amaterial);
#ifdef ddebug
    cout << *amaterial << endl;
#endif
    return amaterial;
  }  
  else {
    CMSMaterial::FractionType ft;
    if (md == byWeight)
      ft=CMSMaterial::FTWeight;
    else
      ft=CMSMaterial::FTVolume;
    CMSMaterial* material = new CMSMaterial(name, density, nconst, 
					    matcol, prop, ft);
    delete[] matcol;
    theCMSMaterials.push_back(material);
#ifdef ddebug
    cout << *material << endl;
#endif
    return material;
  }  
}



void CMSMaterialFactory::readElements(ifstream& is){
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
    
    G4Element* theEl = addElement(name, symbol, Z, A, density);
    
    readName(is,name);
  };
  cout << "     " << G4Element::GetElementTable()->size() 
       << " elements read from file" << endl << endl;
}
void CMSMaterialFactory::readMaterials(ifstream& is){
  G4String name, matname;

  cout << "     ==> Reading materials... " << endl;

  //Take into account the special case of vacuum...
#ifdef debug
  cout <<"       \"Vacuum\"" << endl;
#endif
  G4double density  = universe_mean_density;   //from PhysicalConstants.h
  G4double pressure = 1.E-19*pascal;
  G4double temperature = 0.1*kelvin;
  G4Material *vacuum = new G4Material("Vacuum",
				      /*Z=*/ 1.,
				      /*A=*/ 1.01*g/mole,
				      density,
				      kStateGas,
				      temperature,
				      pressure);

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

    addCMSMaterial(matname, density, absnelem, mats, weights, md);
    delete[] mats;
    delete[] weights;
    
    readName(is,name);
  };  //while


  cout << "     " << theCMSMaterials.size() << " materials read from " 
       << mixturefile << endl << endl;
}

CMSMaterialFactory::CMSMaterialFactory() {
  readElements (elementfile);
  readMaterials(mixturefile);
}

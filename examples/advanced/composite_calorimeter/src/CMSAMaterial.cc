///////////////////////////////////////////////////////////////////////////////
// File: CMSAMaterial.cc
// Date: 16/03/98 I. Gonzalez
// Modifications: 31/08/98 I.G. -> G4type moved to type for int, double.
///////////////////////////////////////////////////////////////////////////////
#include "CMSAMaterial.hh"

CMSAMaterial::CMSAMaterial(G4String mat, G4double dens, int nconst, 
			   CMSAMaterial** constituents, G4double* weights) {
  name=mat;
  nElem=0;
  int i=0;
  for (i=0; i<nconst; i++)
    nElem += (constituents[i]->NElements());

  theElements = new G4String[nElem];
  theWeights  = new double[nElem];

  double factor;
  int nelem=0;
  for (i=0; i<nconst; i++) {
    factor=constituents[i]->Aeff();
    for (int j=0; j<constituents[i]->NElements(); j++) {
      theElements[nelem] = constituents[i]->Element(j);
      theWeights[nelem]  = constituents[i]->Weight(j)* weights[i] * factor;
      nelem++;
    }
  }

  if (dens>0) 
    density=dens;
  else //Let's compute density
    computeDensity(nconst,(CMSMaterial**)constituents, weights, FTVolume);

  computeAeff(nconst, constituents, weights);
  closeMaterial();
}

CMSAMaterial::CMSAMaterial(G4String elemat, double Aeff, double dens) {
  name=elemat;
  density=dens;
  nElem=1;
  theElements = new G4String[nElem];
  theWeights  = new double[nElem];
  
  theElements[0] = elemat;
  theWeights[0]  = 1;

  aEff=Aeff;
}

CMSAMaterial::~CMSAMaterial() {
  //The base class destructor is called?
}

CMSAMaterial::CMSAMaterial(const CMSAMaterial& mat){
  name    = mat.name;
  density = mat.density;
  nElem   = mat.nElem;
  theElements = new G4String[nElem];
  theWeights  = new double[nElem];
  for (int i=0; i<nElem; i++){
    theElements[i]=mat.theElements[i];
    theWeights[i]=mat.theWeights[i];
  }
}

CMSAMaterial& CMSAMaterial::operator=(const CMSAMaterial& mat){
  if(theElements)
    delete[] theElements;
  if(theWeights)
    delete[] theWeights;

  name=mat.name;
  density=mat.density;
  nElem=mat.nElem;
  aEff=mat.aEff;
  
  theElements = new G4String[nElem];
  theWeights  = new double[nElem];
  for (int i=0; i<nElem; i++){
    theElements[i]=mat.theElements[i];
    theWeights[i]=mat.theWeights[i];
  }
  return *this;
}

void CMSAMaterial::computeAeff(int nconst, 
			       CMSAMaterial** constituents, 
			       double* weights){
  aEff=0;
  for (int i=0; i<nconst; i++)
    aEff += weights[i] * constituents[i]->Aeff();
}

ostream& operator<<(ostream& os, const CMSAMaterial& mat) {
  os << mat.name << endl;
  os << "Density= " << mat.density << " g/cm3. Number of Elements: "
     << mat.nElem 
     << ". Aeff= " << mat.aEff << endl;
  for (int i=0; i<mat.nElem; i++)
    os << '\t' << mat.theElements[i] << '\t' << mat.theWeights[i] << endl;
  return os;
}

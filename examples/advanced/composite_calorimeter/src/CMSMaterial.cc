///////////////////////////////////////////////////////////////////////////////
// File: CMSMaterial.cc
// Date: 12/03/98 I. Gonzalez
// Modifications: 31/08/98 I.G. -> G4type moved to type for int, double.
///////////////////////////////////////////////////////////////////////////////
#include "CMSMaterial.hh"

//#define debug

CMSMaterial::CMSMaterial(G4String mat, double dens, int nconst, 
			 CMSMaterial** constituents, double* weights,
			 FractionType ft):
  name(mat), density(dens)
{
  nElem = 0;
  
  int i=0;
  for (i=0; i<nconst; i++)
    nElem += constituents[i]->NElements();

  theElements = new G4String[nElem];
  theWeights  = new double[nElem];

  double factor;
  int nelem=0;
  for (i=0; i<nconst; i++) {
    if (ft==FTWeight)
      factor=1.0;
    else
      factor=constituents[i]->Density();
    for (int j=0; j<constituents[i]->NElements(); j++) {
      theElements[nelem] = constituents[i]->Element(j);
      theWeights[nelem]  = constituents[i]->Weight(j)* weights[i] * factor;
      nelem++;
    }
  }

  if (density<0) { //Let's compute density
    computeDensity(nconst, constituents, weights, ft);
  }
  closeMaterial();
}

CMSMaterial::CMSMaterial(const CMSMaterial& mat):
  name(mat.name), density(mat.density), nElem(mat.nElem)
{
  theElements = new G4String[nElem];
  theWeights  = new double[nElem];
  for (int i=0; i<nElem; i++){
    theElements[i]=mat.theElements[i];
    theWeights[i]=mat.theWeights[i];
  }
}

CMSMaterial::~CMSMaterial() {
  if (theElements)
    delete[] theElements;
  if (theWeights)
    delete[] theWeights;
}

void CMSMaterial::computeDensity(int nconst, 
				 CMSMaterial** constituents, double* weights, 
				 FractionType ft){
  double mass=0;
  double volume=0;
  for (int i=0; i<nconst; i++) {
    if (ft==FTWeight) {
      mass+=weights[i];
      volume+=(weights[i]/constituents[i]->Density());
    }
    else { //by volume
      mass+=(weights[i]*constituents[i]->Density());
      volume+=weights[i];
    }
  }
  density=mass/volume;
}

CMSMaterial& CMSMaterial::operator=(const CMSMaterial& mat){
  if(theElements)
    delete[] theElements;
  if(theWeights)
    delete[] theWeights;

  name=mat.name;
  density=mat.density;
  nElem=mat.nElem;
  
  theElements = new G4String[nElem];
  theWeights  = new double[nElem];
  for (int i=0; i<nElem; i++){
    theElements[i]=mat.theElements[i];
    theWeights[i]=mat.theWeights[i];
  }
  return *this;
}

G4bool CMSMaterial::operator==(const CMSMaterial& mat) const{
  return (name==mat.name);
}

G4bool CMSMaterial::operator!=(const CMSMaterial& mat) const{
  return (name!=mat.name);
}

void CMSMaterial::closeMaterial() {
  int trueConst=0;

  double norm=0;

  for (int i=0; i<nElem; i++) {
    norm+=theWeights[i];
    if (theElements[i]!="") {
      trueConst++;
      for (int j=i+1; j<nElem; j++) {
	if(theElements[i]==theElements[j]){
	  theWeights[i]+=theWeights[j];
	  theElements[j]="";
	}
      }//for j
    } //if
  }//for i

  if (trueConst != nElem) {
    G4String* newConst = new G4String[trueConst];
    double* newWeight = new double[trueConst];
    
    int newi=0;
    for(int i=0; i<nElem; i++){
      if (theElements[i]!="") {
	newConst[newi]  = theElements[i];
	newWeight[newi] = theWeights[i]/norm;
	newi++;
      }
    }

#ifdef debug    
    cout << "\tGoing from " << nElem <<" constituents to " << trueConst << endl;
#endif
    nElem=trueConst;
    
    delete[] theElements;
    delete[] theWeights;

    theElements=newConst;
    theWeights=newWeight;
  }
  else { //Let's normalize the weights
    for (int i=0; i<nElem; i++)
      theWeights[i] = theWeights[i]/norm;
  }
}

ostream& operator<<(ostream& os, const CMSMaterial& mat) {
  os << mat.name << endl;
  os << "Density= " << mat.density << " g/cm3. Number of Elements: "
     << mat.nElem << endl;
  for (int i=0; i<mat.nElem; i++)
    os << '\t' << mat.theElements[i] << '\t' << mat.theWeights[i] << endl;
  return os;
}

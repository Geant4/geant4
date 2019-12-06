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
// File: CCalMaterial.cc
// Description: CCalMaterial transient class to store information from 
//              material table (database) which is used to make a G4Material
///////////////////////////////////////////////////////////////////////////////
#include "CCalMaterial.hh"

//#define debug

CCalMaterial::CCalMaterial(G4String mat, G4double dens, G4int nconst, 
                           CCalMaterial** constituents, G4double* weights,
                           FractionType ft): name(mat), density(dens) {
  nElem = 0;
  
  G4int i=0;
  for (i=0; i<nconst; i++)
    nElem += constituents[i]->NElements();

  theElements = new G4String[nElem];
  theWeights  = new G4double[nElem];

  G4double factor;
  G4int nelem=0;
  for (i=0; i<nconst; i++) {
    if (ft==FTWeight)
      factor=1.0;
    else
      factor=constituents[i]->Density();
    for (G4int j=0; j<constituents[i]->NElements(); j++) {
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

CCalMaterial::CCalMaterial(const CCalMaterial& mat):
  name(mat.name), density(mat.density), nElem(mat.nElem) {
  theElements = new G4String[nElem];
  theWeights  = new G4double[nElem];
  for (G4int i=0; i<nElem; i++){
    theElements[i]=mat.theElements[i];
    theWeights[i]=mat.theWeights[i];
  }
}

CCalMaterial::~CCalMaterial() {
  if (theElements)
    delete[] theElements;
  if (theWeights)
    delete[] theWeights;
}

void CCalMaterial::computeDensity(G4int nconst, 
                                  CCalMaterial** constituents, 
                                  G4double* weights, FractionType ft) {
  G4double mass=0;
  G4double volume=0;
  for (G4int i=0; i<nconst; i++) {
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

CCalMaterial& CCalMaterial::operator=(const CCalMaterial& mat){
  if(theElements)
    delete[] theElements;
  if(theWeights)
    delete[] theWeights;

  name=mat.name;
  density=mat.density;
  nElem=mat.nElem;
  
  theElements = new G4String[nElem];
  theWeights  = new G4double[nElem];
  for (G4int i=0; i<nElem; i++){
    theElements[i]=mat.theElements[i];
    theWeights[i]=mat.theWeights[i];
  }
  return *this;
}

G4bool CCalMaterial::operator==(const CCalMaterial& mat) const{
  return (name==mat.name);
}

G4bool CCalMaterial::operator!=(const CCalMaterial& mat) const{
  return (name!=mat.name);
}

void CCalMaterial::closeMaterial() {
  G4int trueConst=0;

  G4double norm=0;

  for (G4int i=0; i<nElem; i++) {
    norm+=theWeights[i];
    if (theElements[i]!="") {
      trueConst++;
      for (G4int j=i+1; j<nElem; j++) {
        if(theElements[i]==theElements[j]){
          theWeights[i]+=theWeights[j];
          theElements[j]="";
        }
      }//for j
    } //if
  }//for i

  if (trueConst != nElem) {
    G4String* newConst = new G4String[trueConst];
    G4double* newWeight = new G4double[trueConst];
    
    G4int newi=0;
    for(G4int i=0; i<nElem; i++){
      if (theElements[i]!="") {
        newConst[newi]  = theElements[i];
        newWeight[newi] = theWeights[i]/norm;
        newi++;
      }
    }

#ifdef debug    
    G4cout << "\tGoing from " << nElem <<" constituents to " << trueConst << G4endl;
#endif
    nElem=trueConst;
    
    delete[] theElements;
    delete[] theWeights;

    theElements=newConst;
    theWeights=newWeight;
  }
  else { //Let's normalize the weights
    for (G4int i=0; i<nElem; i++)
      theWeights[i] = theWeights[i]/norm;
  }
}

std::ostream& operator<<(std::ostream& os, const CCalMaterial& mat) {
  os << mat.name << G4endl;
  os << "Density= " << mat.density << " g/cm3. Number of Elements: "
     << mat.nElem << G4endl;
  for (G4int i=0; i<mat.nElem; i++)
    os << '\t' << mat.theElements[i] << '\t' << mat.theWeights[i] << G4endl;
  return os;
}

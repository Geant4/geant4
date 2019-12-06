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
// File: CCalAMaterial.cc
// Description: Specialised class to store information to make G4Material 
//              from atomic proportion
///////////////////////////////////////////////////////////////////////////////
#include "CCalAMaterial.hh"

CCalAMaterial::CCalAMaterial(G4String mat, G4double dens, G4int nconst, 
                             CCalAMaterial** constituents, G4double* weights) {
  name=mat;
  nElem=0;
  G4int i=0;
  for (i=0; i<nconst; i++)
    nElem += (constituents[i]->NElements());

  theElements = new G4String[nElem];
  theWeights  = new G4double[nElem];

  G4double factor;
  G4int nelem=0;
  for (i=0; i<nconst; i++) {
    factor=constituents[i]->Aeff();
    for (G4int j=0; j<constituents[i]->NElements(); j++) {
      theElements[nelem] = constituents[i]->Element(j);
      theWeights[nelem]  = constituents[i]->Weight(j)* weights[i] * factor;
      nelem++;
    }
  }

  if (dens>0) 
    density=dens;
  else //Let's compute density
    computeDensity(nconst,(CCalMaterial**)constituents, weights, FTVolume);

  computeAeff(nconst, constituents, weights);
  closeMaterial();
}

CCalAMaterial::CCalAMaterial(G4String elemat, G4double eff, G4double dens) {
  name=elemat;
  density=dens;
  nElem=1;
  theElements = new G4String[nElem];
  theWeights  = new G4double[nElem];
  
  theElements[0] = elemat;
  theWeights[0]  = 1;

  aEff=eff;
}

CCalAMaterial::~CCalAMaterial() {
  //The base class destructor is called?
}

CCalAMaterial::CCalAMaterial(const CCalAMaterial& mat) 
  : CCalMaterial( mat ) {
  name    = mat.name;
  density = mat.density;
  nElem   = mat.nElem;
  theElements = new G4String[nElem];
  theWeights  = new G4double[nElem];
  for (G4int i=0; i<nElem; i++){
    theElements[i]=mat.theElements[i];
    theWeights[i]=mat.theWeights[i];
  }
}

CCalAMaterial& CCalAMaterial::operator=(const CCalAMaterial& mat){
  if(theElements)
    delete[] theElements;
  if(theWeights)
    delete[] theWeights;

  name=mat.name;
  density=mat.density;
  nElem=mat.nElem;
  aEff=mat.aEff;
  
  theElements = new G4String[nElem];
  theWeights  = new G4double[nElem];
  for (G4int i=0; i<nElem; i++){
    theElements[i]=mat.theElements[i];
    theWeights[i]=mat.theWeights[i];
  }
  return *this;
}

void CCalAMaterial::computeAeff(G4int nconst, 
                                CCalAMaterial** constituents, 
                                G4double* weights){
  aEff=0;
  for (G4int i=0; i<nconst; i++)
    aEff += weights[i] * constituents[i]->Aeff();
}

std::ostream& operator<<(std::ostream& os, const CCalAMaterial& mat) {
  os << mat.name << G4endl;
  os << "Density= " << mat.density << " g/cm3. Number of Elements: "
     << mat.nElem 
     << ". Aeff= " << mat.aEff << G4endl;
  for (G4int i=0; i<mat.nElem; i++)
    os << '\t' << mat.theElements[i] << '\t' << mat.theWeights[i] << G4endl;
  return os;
}

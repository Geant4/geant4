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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// 
/// \file PDBmolecule.cc
/// \brief Implementation file for PDBmolecule class

#include "PDBmolecule.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Molecule::Molecule():
fMolName(""),fMolNum(0),fMinGlobZ(0),fMaxGlobZ(0),
fMinGlobX(0),fMaxGlobX(0),fMinGlobY(0),fMaxGlobY(0),
fCenterX(0),fCenterY(0),fCenterZ(0),fDistCenterMax(0),fNbResidue(0),
fpNext(0),fpFirst(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Molecule::Molecule(const std::string& mN,int mNum)
{
  fMolName=mN;  //Molecule name
  fMolNum=mNum; //Molecule number
  fMinGlobZ=0;
  fMaxGlobZ=0;
  fMinGlobX=0;
  fMaxGlobX=0;
  fMinGlobY=0;
  fMaxGlobY=0;
  fCenterX=0;
  fCenterY=0;
  fCenterZ=0;
  fDistCenterMax=0;
  fNbResidue=0;
  fpNext=0;
  fpFirst=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Molecule *Molecule::GetNext()
{
  return fpNext;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Residue *Molecule::GetFirst()
{
  return fpFirst;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Molecule::GetID()
{
  return fMolNum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Molecule::SetNext(Molecule *moleculeNext)
{
  fpNext=moleculeNext;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Molecule::SetFirst(Residue *resFirst)
{
  fpFirst=resFirst;
}


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
// --------------------------------------------------------------
// Authors: E. Delage
// november 2013
// --------------------------------------------------------------
//
//
/// \file PDBmolecule.hh
/// \brief Definition of the PDBmolecule class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MOLECULE_H
#define MOLECULE_H

#include "PDBresidue.hh"

class Residue;

//! Molecule Class
/*!
 * This Class define Molecule model ... 
 */
class Molecule
{
public:
  //! First constructor
  Molecule();
  //! Second constructor
  Molecule(const std::string& resName,int mNum);
  //! Destructor
  ~Molecule() {};

  //! information about molecule (not implemented)
  //void PrintInfo();
  //! Get the next molecule
  Molecule *GetNext();
  //! Get the first Residue
  Residue *GetFirst();
  //! Get number Molecule
  int GetID();
  //! Set the next Molecule
  void SetNext(Molecule *);
  //! Set the first Residue
  void SetFirst(Residue *);

  std::string fMolName;   //!< Molecule name
  int fMolNum;       //!< Molecule number

  double fMinGlobZ;   //Cylinder length => min Z
  double fMaxGlobZ;
  double fMinGlobX;   //Radius => min X
  double fMaxGlobX;
  double fMinGlobY;   //=> min Y
  double fMaxGlobY;

  int fCenterX;      //!< "X center" of this Molecule (for rotation...)
  int fCenterY;      //!< "Y center" of this Molecule (for rotation...)
  int fCenterZ;//!< "Z center" of this Molecule (for rotation...)
  int fDistCenterMax;//!< dist from center to most away most of the molecule
  int fNbResidue;        //!< Number of residue inside the molecule

private:
  Molecule *fpNext;//!< Header of the next Molecule (usage before vector)
  Residue *fpFirst;//!< Header of the first Residue (usage before vector)
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

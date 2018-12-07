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
/// \file PDBresidue.hh
/// \brief Definition of the PDBresidue class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RESIDUE_H
#define RESIDUE_H

#include "PDBatom.hh"

class Atom;

//! Residue Class
/*!
 * This Class define residue model ... 
 */
class Residue
{
public:
  //! Constructor
  Residue();
  //! Second constructor
  Residue(const std::string& resName,int resSeq);
  //!Destructor
  ~Residue() {};

  //! information about molecule (not implemented)
  void Affiche();
  //! Get the next residue
  Residue *GetNext();
  //! Get the first atom
  Atom *GetFirst();
  //! Get the number of the residue
  int GetID();
  //! Set the next residue
  void SetNext(Residue *);
  //! Set the first Atom of the residue
  void SetFirst(Atom *);

  std::string fResName;       //!< Residue name
  int fResSeq;           //!< Residue sequence number
  bool fVisible;         //!< Whether Residue is visible or not
  bool fSelected;     //!< Whether Residue is selected (Highlight) or not
  int fCenterX;
  int fCenterY;
  int fCenterZ;
  int fNbAtom;       //!< Number of atoms into a residue (usage before vector)

private:
  //! Residue header for next residue  (usage before vector)
  Residue *fpNext;
  //! Atom header for first atom of residue (usage before vector)
  Atom *fpFirst;
};

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

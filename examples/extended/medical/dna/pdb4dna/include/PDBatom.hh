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
/// \file PDBatom.hh
/// \brief Definition of the Atom class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ATOM_H
#define ATOM_H

#include <string>

//! Atom Class
/*!
 * This Class define Atom model ... 
 */
class Atom
{
public:
  //! constructor with initialization
  Atom(int serial, 
       const std::string& name, 
       const std::string& resName,
       int numInRes,
       int resSeq,
       double xInit,
       double yInit,
       double zInit,
       double radius,
       double occupancy, 
       double tempFactor, 
       const std::string& element);
       
  //! Empty destructor
  ~Atom()
  {
  };

  //! Returns the next Atom
  Atom *GetNext();
  //! Return the X position for the Atom
  double GetX();
  //! Return the Y position for the Atom
  double GetY();
  //! Return the Z position for the Atom
  double GetZ();
  //! Return the Atom's ID
  int GetID();
  //! Return name of the atom
  const std::string& GetName();
  //! Return name of the element
  const std::string& GetElementName();
  //! Return name of the atom
  double GetVanDerWaalsRadius();
  //! Set the next atom
  void SetNext(Atom *);

  int fSerial;       //!< its serial number
  int fNumInRes;     //!< its number in residue sequence
  std::string fName;      //!< Atom name
  std::string fResName;   //!< Residue name
  int fResSeq;       //!< Residue sequence number
  double fX;          //!< X orthogonal coordinates in Angstroms
  double fY;          //!< Y orthogonal coordinates in Angstroms
  double fZ;          //!< Z orthogonal coordinates in Angstroms
  double fVdwRadius;  // Vand der Waals Radius in Angstrom
  double fOccupancy;  //!< Occupancy for the Atom
  std::string fElement;   //!< Element symbol extracted from 'atom name'
  double fTempFactor; //!< Temperature factor for the Atom

private:
  Atom * fpNext;       //!< Pointer to the next Atom
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

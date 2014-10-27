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
// $Id$
//
/// \file PDBbarycenter.hh
/// \brief Definition of the Barycenter class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BARY_H
#define BARY_H

//! Molecule Class
/*!
 * This Class define Molecule model ... 
 */
class Barycenter
{
public:
  //! First constructor
  Barycenter();
  //! Second constructor
  Barycenter(int bNum,float x,float y, float z,//Nucleotide bar. coordinates
      float Bx,float By, float Bz, //Base bar. coordinates
      float Sx,float Sy, float Sz, //Sugar bar. coordinates
      float Px,float Py, float Pz);//Phosphate bar. coordinates
  //! Destructor
  ~Barycenter() {};

  //! Get the next Barycenter
  Barycenter *GetNext();
  //! Get the first
  //Residue *GetFirst();
  //! Get number Barycenter
  int GetID();
  //! Set the next Barycenter
  void SetNext(Barycenter *);
  //! Set the distance between atom i and nucleotide barycenter
  void SetDistance(int i, float);
  //! Get the distance between atom i and nucleotide barycenter
  float GetDistance(int i);
  //! Set the distance between the farther atom and nucleotide barycenter
  void SetRadius(float );
  //! Get the distance between the farther atom and nucleotide barycenter
  float GetRadius();

  int fBaryNum;//!< Barycenter number
  float fDistanceTab[33];//!< distance table [0..32] (11 hydrogens!)
  float fRadius;

  float fCenterX;            //!< "X coordinate" of this nucelotide Barycenter
  float fCenterY;            //!< "Y coordinate" of this nucelotide Barycenter
  float fCenterZ;            //!< "Z coordinate" of this nucelotide Barycenter

  float fCenterBaseX;        //!< "X coordinate" of this Base Barycenter
  float fCenterBaseY;        //!< "Y coordinate" of this Base Barycenter
  float fCenterBaseZ;        //!< "Z coordinate" of this Base Barycenter

  float fCenterSugarX;       //!< "X coordinate" of this Sugar Barycenter
  float fCenterSugarY;       //!< "Y coordinate" of this Sugar Barycenter
  float fCenterSugarZ;       //!< "Z coordinate" of this Sugar Barycenter

  float fCenterPhosphateX;   //!< "X coordinate" of this Phosphate Barycenter
  float fCenterPhosphateY;   //!< "Y coordinate" of this Phosphate Barycenter
  float fCenterPhosphateZ;   //!< "Z coordinate" of this Phosphate Barycenter

private:
  Barycenter *fpNext;    //!< Header of the next Molecule (usage before vector)
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/include/G4VGlauberDataSet.hh
/// \brief Definition of the G4VGlauberDataSet class
//
// $Id: G4VGlauberDataSet.hh 81932 2014-06-06 15:39:45Z gcosmo $
//

#ifndef G4VGlauberDataSet_h
#define G4VGlauberDataSet_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4VGlauberDataSet.hh
//
// Version:             0.B
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
//
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class G4VGlauberDataSet
{
public:

  G4VGlauberDataSet ();
  virtual ~G4VGlauberDataSet ();

  G4VGlauberDataSet (const G4VGlauberDataSet &right);
  const G4VGlauberDataSet& operator=(G4VGlauberDataSet &right);
    
  inline G4int GetVerboseLevel () const;
  inline void SetVerboseLevel (const G4int i);

  G4int GetAP () const;
  G4int GetZP () const;
  G4int GetAT () const;
  G4int GetZT () const;

  G4int GetGlauberDataSetType () const;
    
  void SetArrayPointer (const G4int i);
  inline G4double *GetArrayPointerN (const G4double ppn = 0.0);
  inline G4double *GetArrayPointerM (const G4double ppn = 0.0);
    
  virtual std::ofstream & WriteDataToFile (std::ofstream &File) const;
  virtual std::ifstream & ReadDataFromFile (std::ifstream &File);

public:

  G4double      rproj;
  G4double      rtarg;
  G4double      bstep;
  G4double      bmax;
  G4int         AP;
  G4int         ZP;
  G4int         AT;
  G4int         ZT;

  G4double    * baseArrayPtrn;
  G4double    * baseArrayPtrm;
  G4double    * arrayPtrn;
  G4double    * arrayPtrm;

  G4int         glauberDataSetType;

  G4int         maxArray;
  G4int         maxig;
    
  G4int         verboseLevel; 
    
  friend std::ofstream & operator << (std::ofstream &File, 
                                      const G4VGlauberDataSet &q);
  friend std::ifstream & operator >> (std::ifstream &File, 
                                      G4VGlauberDataSet &q);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void G4VGlauberDataSet::SetVerboseLevel (const G4int verboseLevel1)
{verboseLevel = verboseLevel1;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline G4int G4VGlauberDataSet::GetVerboseLevel () const
{return verboseLevel;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

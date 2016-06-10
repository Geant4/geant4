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
/// \file hadronic/Hadr02/include/G4ParamType1GlaubAADataSet.hh
/// \brief Definition of the G4ParamType1GlaubAADataSet class
//
// $Id: G4ParamType1GlaubAADataSet.hh 81932 2014-06-06 15:39:45Z gcosmo $
//

#ifndef G4ParamType1GlaubAADataSet_h
#define G4ParamType1GlaubAADataSet_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4ParamType1GlaubAADataSet.hh
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

#include "G4FullGlaubAADataSet.hh"
#include "G4GlaubAADataSet.hh"
#include "G4Type1GlauberParameterisation.hh"

#include "G4DPMJET2_5Interface.hh"

#include "globals.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class G4ParamType1GlaubAADataSet : public G4GlaubAADataSet,
                                   public G4Type1GlauberParameterisation 
{
public:

  G4ParamType1GlaubAADataSet ();
  G4ParamType1GlaubAADataSet (const G4int AP1, const G4int AT1);
  G4ParamType1GlaubAADataSet (G4FullGlaubAADataSet *fullGlauberDataSet);
  ~G4ParamType1GlaubAADataSet ();
  G4ParamType1GlaubAADataSet (const G4ParamType1GlaubAADataSet &right);

  const 
  G4ParamType1GlaubAADataSet& operator=(G4ParamType1GlaubAADataSet &right);

  G4bool CreateGlauberData (const G4int AP1, const G4int AT1);
  G4bool CreateGlauberData (G4FullGlaubAADataSet *fullGlauberDataSet);

  inline G4double GetValueN (const G4double f, const G4double ppn = 0.0) const;
  inline G4double GetValueM (const G4double f, const G4double ppn = 0.0) const;
    
  std::ofstream & WriteDataToFile (std::ofstream &File) const;
  std::ifstream & ReadDataFromFile (std::ifstream &File);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// GetValueN
//
inline 
G4double G4ParamType1GlaubAADataSet::GetValueN (const G4double f,
                                                const G4double ppn) const
{return GetParameterisedValueN (f, ppn);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// GetValueM
//
inline 
G4double G4ParamType1GlaubAADataSet::GetValueM (const G4double f,
                                                const G4double ppn) const
{return GetParameterisedValueM (f, ppn);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

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
/// \file hadronic/Hadr02/include/G4Type1GlauberParameterisation.hh
/// \brief Definition of the G4Type1GlauberParameterisation class
//
#ifndef G4Type1GlauberParameterisation_h
#define G4Type1GlauberParameterisation_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4Type1GlauberParameterisation.hh
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
#include <fstream>
////////////////////////////////////////////////////////////////////////////////
//
class G4Type1GlauberParameterisation
{
  public:
    G4Type1GlauberParameterisation ();
   ~G4Type1GlauberParameterisation ();

    G4double GetFitParameters (const G4double *bsite, G4double *p) const;

    G4double GetParameterisedValueN (const G4double f, const G4double ppn) const;
    G4double GetParameterisedValueM (const G4double f, const G4double ppn) const;

  public:
    G4int    maxArrayp;
    G4int    maxigp;

    G4double paramn[24][10];
    G4double paramm[24][10];
    
    G4double mun1[24];
    G4double mun2[24];
    G4double cn[24];
    G4double mum1[24];
    G4double mum2[24];
    G4double cm[24];

    G4double limit1;
    G4double limit2;
    G4double limit3;
    G4double limit4;
};
////////////////////////////////////////////////////////////////////////////////
//
#endif

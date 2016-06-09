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
/// \file hadronic/Hadr02/include/G4DPMJET2_5CrossSectionParamSet.hh
/// \brief Definition of the G4DPMJET2_5CrossSectionParamSet class
//
#ifndef G4DPMJET2_5CrossSectionParamSet_h
#define G4DPMJET2_5CrossSectionParamSet_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4DPMJET2_5CrossSectionParamSet.hh
//
// Version:             0.A
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


///////////////////////////////////////////////////////////////////////////////
//
class G4DPMJET2_5CrossSectionParamSet
{
  public:
    G4DPMJET2_5CrossSectionParamSet () {c[0]=0.0; c[1]=0.0; c[2]=0.0;}
    G4DPMJET2_5CrossSectionParamSet (const G4double c0, const G4double c1,
      const G4double c2) {c[0]=c0; c[1]=c1; c[2]=c2;}
    ~G4DPMJET2_5CrossSectionParamSet () {};
    
    G4double operator[](const G4int i) const {return c[i];}
    G4DPMJET2_5CrossSectionParamSet &operator=(G4DPMJET2_5CrossSectionParamSet &right);
    G4DPMJET2_5CrossSectionParamSet operator=(G4DPMJET2_5CrossSectionParamSet right)
      {c[0] = right.c[0]; c[1] = right.c[1]; c[2] = right.c[2]; return *this;};
  
  protected:
    G4double c[3];
};
#endif

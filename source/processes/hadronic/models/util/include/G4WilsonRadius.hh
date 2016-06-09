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
// * contract 17191/03/NL/LvH (Aurora Programme).                     *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef G4WilsonRadius_h
#define G4WilsonRadius_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4WilsonRadus.hh
//
// Version:		B.1
// Date:		15/04/04
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
// Provides functions to calculate the nuclear radius as used in the
// abrasion-ablation model (see Wilson, Tripathi, Cucinotta et al, NASA
// Technical Paper 3533, October 1995), and to determine total nuclear-nuclear
// cross-sections.  The Algorithm comes from Tripathi, Cucinnota and Wilson,
// Nucl Instrum Meth Phys Res B, 152, pp425-431, 1999.
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"

////////////////////////////////////////////////////////////////////////////////
//
class G4WilsonRadius
{
  public:
    G4WilsonRadius ();
    ~G4WilsonRadius ();
  
    G4double GetWilsonRMSRadius (G4double A);
    G4double GetWilsonRadius (G4double A);

  private:
    G4double r0sq;
    G4double factor;
    G4double third;
};
////////////////////////////////////////////////////////////////////////////////
//
#endif

#ifndef G4NuclearAbrasionGeometry_h
#define G4NuclearAbrasionGeometry_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4NuclearAbrasionGeometry.hh
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
// 18 November 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation  is  the  intellectual   property  of *
// * QinetiQ.   Rights  to use  this  software  are  granted  to  the *
// * European Space Agency  under  the  provisions  of  Clause 39  of *
// * the    Standard    Conditions   for  Contract.   ESA    contract *
// * 17191/03/NL/LvH  provides  the  rights to distribute source code *
// * through  the  Geant4 Collaboration Agreement for general  public *
// * use.  Some elements of this source code may be derived from code *
// * developed   by   other   member   institutes   of   the   Geant4 *
// * Collaboration, and the provision of  their code under the Geant4 *
// * Collaboration is acknowledged.                                   *
// *                                                                  *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"

////////////////////////////////////////////////////////////////////////////////
//
class G4NuclearAbrasionGeometry
{
  public:
    G4NuclearAbrasionGeometry (G4double AP, G4double AT, G4double r);
    ~G4NuclearAbrasionGeometry ();
    void SetPeripheralThreshold (G4double);
    G4double GetPeripheralThreshold ();
    G4double F ();
    G4double P ();
    G4double GetExcitationEnergyOfProjectile ();
    G4double GetExcitationEnergyOfTarget ();

  private:
    G4double AP;
    G4double AT;
    G4double rP;
    G4double rT;
    G4double r;
    G4double n;
    G4double b;
    G4double m;

    G4double Q;
    G4double S;
    G4double T;
    G4double R;
    G4double U;
    
    G4double rth;
    G4double B;
};
#endif

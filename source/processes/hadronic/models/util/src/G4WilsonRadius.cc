// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4WilsonRadius.cc
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
//
////////////////////////////////////////////////////////////////////////////////
//
#include "G4WilsonRadius.hh"

////////////////////////////////////////////////////////////////////////////////
//
G4WilsonRadius::G4WilsonRadius ()
{
  G4double r0 = 0.84*fermi;
  r0sq        = r0 * r0;
  factor      = sqrt(5.0/3.0) * fermi;
  third       = 1.0 / 3.0;
}
////////////////////////////////////////////////////////////////////////////////
//
G4WilsonRadius::~G4WilsonRadius ()
{;}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4WilsonRadius::GetWilsonRMSRadius (G4double A)
{
  G4double radius;
  if (A > 26.0)
    radius = factor * (0.84*pow(A,third) + 0.55);
  else
  {
    G4double r[27] = {0.0, 0.85,  2.095, 1.976, 1.671, 1.986,
                           2.57,  2.41,  2.23,  2.519, 2.45,
                           2.42,  2.471, 2.440, 2.58,  2.611,
                           2.730, 2.662, 2.727, 2.9,   3.040,
                           2.867, 2.969, 2.94,  3.075, 3.11,
                           3.06};
    radius = factor * r[(G4int) (A+0.4)];
  }
  return radius;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4WilsonRadius::GetWilsonRadius (G4double A)
{
  G4double r  = GetWilsonRMSRadius(A);
  return  1.29*sqrt(r*r-r0sq);
}
////////////////////////////////////////////////////////////////////////////////
//

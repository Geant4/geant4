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
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the	      *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme). 		      *
// *                                                                  *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4NuclearAbrasionGeometry.cc
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
// 4 June 2004, J.P. Wellisch, CERN, Switzerland
// resolving technical portability issues.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4NuclearAbrasionGeometry.hh"
#include "G4WilsonRadius.hh"
////////////////////////////////////////////////////////////////////////////////
//
G4NuclearAbrasionGeometry::G4NuclearAbrasionGeometry (G4double AP1,
  G4double AT1, G4double r1)
{
//
//
// Initialise variables for interaction geometry.
//
  G4WilsonRadius aR;
  AP = AP1;
  AT = AT1;
  rP = aR.GetWilsonRadius(AP);
  rT = aR.GetWilsonRadius(AT);
  r  = r1;
  n  = rP / (rP + rT);
  b  = r / (rP + rT);
  m  = rT / rP;
  Q  = (1.0 - b)/n;
  S  = Q * Q;
  T  = S * Q;
  R  = sqrt(m*n);
  U  = 1.0/m - 2.0;
//
//
// Initialise the threshold radius-ratio at which interactions are considered
// peripheral or central.
//  
  rth = 2.0/3.0;
  B   = 10.0 * MeV;
}
////////////////////////////////////////////////////////////////////////////////
//
G4NuclearAbrasionGeometry::~G4NuclearAbrasionGeometry ()
{;}
////////////////////////////////////////////////////////////////////////////////
//
void G4NuclearAbrasionGeometry::SetPeripheralThreshold (G4double rth1)
  {if (rth1 > 0.0 && rth1 <= 1.0) rth = rth1;}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4NuclearAbrasionGeometry::GetPeripheralThreshold ()
  {return rth;}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4NuclearAbrasionGeometry::P ()
{
//
//
// Initialise the value for P, then determine the actual value depending upon
// whether the projectile is larger or smaller than the target and these radii
// in relation to the impact parameter.
//
  G4double P = 0.0;

  if (rT > rP)
  {
    if (rT-rP<=r && r<=rT+rP) P = 0.125*R*U*S - 0.125*(0.5*R*U+1.0)*T;
    else                      P = -1.0;
  }
  else
  {
    if (rP-rT<=r && r<=rP+rT) P = 0.125*R*U*S - 0.125*(0.5*sqrt(n/m)*U-
      (sqrt(1.0-m*m)/n - 1.0)*sqrt((2.0-m)/pow(m,5.0)))*T;
    else                      P = (sqrt(1.0-m*m)/n-1.0)*sqrt(1.0-b*b/n/n);
  }

  if (!(P <= 1.0 && P>= -1.0))
  {
    if (P > 1.0) P =  1.0;
    else         P = -1.0;
  }
  return P;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4NuclearAbrasionGeometry::F ()
{
//
//
// Initialise the value for F, then determine the actual value depending upon
// whether the projectile is larger or smaller than the target and these radii
// in relation to the impact parameter.
//
  G4double F = 0.0;

  if (rT > rP)
  {
    if (rT-rP<=r && r<=rT+rP) F = 0.75*R*S - 0.125*(3.0*R-1.0)*T;
    else                      F = 1.0;
  }
  else
  {
    if (rP-rT<=r && r<=rP+rT) F = 0.75*R*S - 0.125*(3.0*sqrt(n/m)-
      (1.0-pow(1.0-m*m,3.0/2.0))*sqrt(1.0-pow(1.0-m,2.0))/pow(m,3.0))*T;
    else                      F = (1.0-pow(1.0-m*m,3.0/2.0))*sqrt(1.0-b*b/n/n);
  }

  if (!(F <= 1.0 && F>= 0.0))
  {
    if (F > 1.0) F = 1.0;
    else         F = 0.0;
  }
  return F;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4NuclearAbrasionGeometry::GetExcitationEnergyOfProjectile ()
{
  G4double F1 = F();
  G4double P1 = P();
  G4double Es = 0.0;

  Es = 0.95 * MeV * 4.0 * pi * rP*rP/fermi/fermi *
       (1.0+P1-pow(1.0-F1,2.0/3.0));
//  if (rT < rP && r < rP-rT)
  if ((r-rP)/rT < rth)
  {
    G4double omega = 0.0;
    if      (AP < 12.0)  omega = 1500.0;
    else if (AP <= 16.0) omega = 1500.0 - 320.0*(AP-12.0);
    Es *= 1.0 + F1*(5.0+omega*F1*F1);
  }
  
  if (Es < 0.0) 
    Es = 0.0;
  else if (Es > B * AP)
    Es = B * AP;
  return Es;
}
//////////////////////////////////////////////////////////////////////////////////////
//
G4double G4NuclearAbrasionGeometry::GetExcitationEnergyOfTarget ()
{
//
//
// This member function declares a new G4NuclearAbrasionGeometry object but with the
// projectile and target exchanged to determine the values for F and P.  Determination
// of the excess surface area and excitation energy is as above.
//
  G4NuclearAbrasionGeometry *revAbrasionGeometry =
    new G4NuclearAbrasionGeometry(AT, AP, r);
  G4double F1 = revAbrasionGeometry->F();
  G4double P1 = revAbrasionGeometry->P();
  G4double Es = 0.0;

  Es = 0.95 * MeV * 4.0 * pi * rT*rT/fermi/fermi *
       (1.0+P1-pow(1.0-F1,2.0/3.0));
//  if (rP < rT && r < rT-rP)
  if ((r-rT)/rP < rth)
  {
    G4double omega = 0.0;
    if      (AT < 12.0)  omega = 1500.0;
    else if (AT <= 16.0) omega = 1500.0 - 320.0*(AT-12.0);
    Es *= 1.0 + F1*(5.0+omega*F1*F1);
  }
  
  if (Es < 0.0)
    Es = 0.0;
  else if (Es > B * AT)
    Es = B * AT;
  return Es;
}
////////////////////////////////////////////////////////////////////////////////
//

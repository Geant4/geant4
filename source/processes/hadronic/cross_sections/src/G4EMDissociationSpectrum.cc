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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4EMDissociationSpectrum.cc
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
////////////////////////////////////////////////////////////////////////////////
//
#include "G4EMDissociationSpectrum.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"

////////////////////////////////////////////////////////////////////////////////
//
G4EMDissociationSpectrum::G4EMDissociationSpectrum ()
{
  bessel = new G4Bessel();
}
////////////////////////////////////////////////////////////////////////////////
//
G4EMDissociationSpectrum::~G4EMDissociationSpectrum ()
{
  delete bessel;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4EMDissociationSpectrum::GetGeneralE1Spectrum
  (G4double Eg, G4double b, G4double bmin)
{
  G4double b2 = b*b;
  G4double gg  = 1.0/std::sqrt(1.0-b2);
  G4double xi = Eg * bmin / gg / b / hbarc;
  G4double K0 = bessel->K0(xi);
  G4double K1 = bessel->K1(xi);
  G4double n  = 2.0 * fine_structure_const / pi / b2 / Eg *
    (xi*K0*K1 - xi*xi*b2/2.0*(K1*K1-K0*K0));
  return n;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4EMDissociationSpectrum::GetGeneralE2Spectrum
  (G4double Eg, G4double b, G4double bmin)
{
  G4double b2 = b * b;
  G4double b4 = b2 * b2;
  G4double gg  = 1.0/std::sqrt(1.0-b2);
  G4double xi = Eg * bmin / gg / b / hbarc;
  G4double K0 = bessel->K0(xi);
  G4double K1 = bessel->K1(xi);
  G4double n  = 2.0 * fine_structure_const / pi / b4 / Eg *
    (2.0*(1.0-b2)*K1*K1 + xi*G4Pow::GetInstance()->powA((2.0-b2),2.0)*K0*K1 -
     xi*xi*b4/2.0*(K1*K1-K0*K0));
  return n;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4EMDissociationSpectrum::GetClosestApproach
  (const G4double AP, const G4double ZP, G4double AT, G4double ZT, G4double b)
{
  G4double bsq     = b * b;
  G4double gg       = 1.0/std::sqrt(1-bsq);
  G4double AProot3 = G4Pow::GetInstance()->powA(AP,1.0/3.0);
  G4double ATroot3 = G4Pow::GetInstance()->powA(AT,1.0/3.0);
  G4double bc      = 1.34 * fermi * (AProot3+ATroot3 - 0.75 *(1.0/AProot3+1.0/ATroot3));
//  G4double a0      = ZP * ZT * classic_electr_radius/bsq;
  G4double a0      = ZP * ZT * elm_coupling / (AT*AP*amu_c2/(AT+AP)) / bsq;
  G4double bmin    = 1.25 * bc + halfpi*a0/gg;
  return bmin;
}
////////////////////////////////////////////////////////////////////////////////
//

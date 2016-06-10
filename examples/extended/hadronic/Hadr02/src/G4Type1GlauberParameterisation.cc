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
/// \file hadronic/Hadr02/src/G4Type1GlauberParameterisation.cc
/// \brief Implementation of the G4Type1GlauberParameterisation class
//
// $Id: G4Type1GlauberParameterisation.cc 81932 2014-06-06 15:39:45Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4Type1GlauberParameterisation.cc
//
// Version:             0.B
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#ifdef G4_USE_DPMJET

#include "G4Type1GlauberParameterisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

//
//
// The following function performs a polynomial least squares fit to double
// precision data.  It forms part of the CERN mathlib.  This is a really
// beyond standard policy and will eventually be changed to a C++ function.
//
extern "C" {void dlsqpm_ (int*,double*,double*,int*,double*,double*,int*);}

////////////////////////////////////////////////////////////////////////////////
//
G4Type1GlauberParameterisation::G4Type1GlauberParameterisation ()
{
  limit1    = 0.1;
  limit2    = 0.60;
  limit3    = 0.95;
  limit4    = 0.9999;
  
  maxArrayp = 200;
  maxigp    = 24;
  
  for (G4int ig=0; ig<maxigp; ig++) {
    for (G4int j=0; j<10; j++) {
      paramn[ig][j] = 0.0;
      paramn[ig][j] = 0.0;
    }
    mun1[ig] = 0.0;
    mun2[ig] = 0.0;
    cn[ig]   = 0.0;
    mum1[ig] = 0.0;
    mum2[ig] = 0.0;
    cm[ig]   = 0.0;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
G4Type1GlauberParameterisation::~G4Type1GlauberParameterisation ()
{;}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4Type1GlauberParameterisation::GetFitParameters
  (const G4double *bsite, G4double *p) const
{
//
//
// Initialise parameters.
//
  G4int pt1 = -1;
  G4int pt2 = -1;
  G4int pt3 = -1;
  G4int pt4 = -1;
//
//
// Locate transitions between different parts of the curve where different
// fits are used.
//
  G4double ib[200];
  G4double lnib[200];
  G4double lnbsite[200];
  for (G4int i=0; i<maxArrayp; i++) {
    ib[i]   = (G4double) i;
    lnib[i] = std::log(ib[i]);
    if (bsite[i] < 1.0E-10) lnbsite[i] = -23.02585093;
    else                    lnbsite[i] = std::log(bsite[i]);
    if (pt1 == -1) {
      if (bsite[i] >= limit1) pt1 = i;
    }
    else if (pt2 == -1) {
      if (bsite[i] >= limit2) pt2 = i;
    }
    else if (pt3 == -1) {
      if (bsite[i] >= limit3) pt3 = i;
    }
    else if (pt4 == -1) {
      if (bsite[i] >= limit4) pt4 = i;
    }
  }
//
//
// First section determines the power-law fits for the low and intermediate b
// values;
//
  G4int deltan = pt2 - pt1;
  G4int M      = 1;
  G4double a1[2];
  G4double sd;
  G4int ifail;
  dlsqpm_ (&deltan,&lnbsite[pt1],&lnib[pt1],&M,&a1[0],&sd,&ifail);
  
  G4double c1 = a1[0];
  G4double M1 = a1[1];

//  G4double M2 = (lnib[3] - lnib[2]) / (lnbsite[3] - lnbsite[2]);
//  G4double c2 = lnib[2] - M2*lnbsite[2];
  deltan = 3;
  dlsqpm_ (&deltan,&lnbsite[1],&lnib[1],&M,&a1[0],&sd,&ifail);
  
  G4double c2 = a1[0];
  G4double M2 = a1[1];
  
  p[0] = std::exp(c2);
  p[1] = M2;
  p[2] = std::exp(c1);
  p[3] = M1;
  if (std::abs(M2-M1) > 1.0E-10) {
    p[4] = exp(-(c2-c1)/(M2-M1));
  }
  else {
    p[4] = limit2 / 2.0;
  }
//
//
// This next bit solves for gamma to determine the inflection at high-b values.
// The algorthM used is EXTREEEEEMELY crude .... but practical and robust.
// It's a linear search.
//
  G4double delta = 1.0E+99;
  G4double gam   = 0.0;
  for (G4int ig = 120; ig < 1000; ig++) {
    G4double DELTA = 0.0;
    G4double GAMMA = (G4double) ig / 100.0;
    G4double EXPON = p[3] / GAMMA;
    for (G4int i = pt2; i<pt3; i++) {
      G4double f  = bsite[i];
      G4double B  = p[2] * std::pow(f,p[3]) / 
        std::pow((1.0-std::pow(f,GAMMA)),EXPON);
      G4double epsilon = std::abs((ib[i]-B)/ib[i]);
      if (epsilon > DELTA) DELTA = epsilon;
    }
    if (DELTA < delta) {
      gam   = GAMMA;
      delta = DELTA;
    }
  }
  p[5] = gam;
//
//
// For the final part of the curve, we use a cubic polynomial 
// fit to -ln(1-bsite)
// versus ib.  This does seem to work quite well.
//
  G4double phi[200];
  for (G4int i = pt3; i<pt4; i++) {
    phi[i] = -std::log(1.0 - bsite[i]);
  }
  deltan = pt4-pt3;
  M      = 3;
  G4double a2[4];
  
  dlsqpm_ (&deltan,&phi[pt3],&ib[pt3],&M,&a2[0],&sd,&ifail);
  
  p[6] = a2[0];
  p[7] = a2[1];
  p[8] = a2[2];
  p[9] = a2[3];
  
  return 0.0;
}
////////////////////////////////////////////////////////////////////////////////
//
// GetValueN
//
G4double 
G4Type1GlauberParameterisation::GetParameterisedValueN(const G4double f,
                                                       const G4double ppn) const
{
  G4int ig = 0;
  if (ppn < 1.0E-10) {
    return 0;
  }
  else {
    ig = G4int(2.0*std::log10(ppn)) - 2;
  }
  if      (ig < 0)  ig = 0;
  else if (ig > 23) ig = 23;
    
  G4double v = 0.0;
  if (f <= paramn[ig][4]) {
    v = paramn[ig][0] * std::pow(f,paramn[ig][1]);
  }
  else if (f <= limit3) {
    v = paramn[ig][2] * std::pow(f,paramn[ig][3]) /
      std::pow((1.0 - std::pow(f,paramn[ig][5])),paramn[ig][3]/paramn[ig][5]);
  }
  else {
    G4double l = -std::log(1.0-f);
    v = paramn[ig][6] +
        paramn[ig][7]*l +
        paramn[ig][8]*l*l +
        paramn[ig][9]*std::pow(l,3.0);
  }

  return v;
}
////////////////////////////////////////////////////////////////////////////////
//
// GetValueM
//
G4double 
G4Type1GlauberParameterisation::GetParameterisedValueM(const G4double f,
                                                       const G4double ppn) const
{
  G4int ig = 0;
  if (ppn < 1.0E-10) {
    return 0;
  }
  else {
    ig = G4int(2.0*std::log10(ppn)) - 2;
  }
  if      (ig < 0)  ig = 0;
  else if (ig > 23) ig = 23;
  
  G4double v = 0.0;
  if (f <= paramm[ig][4]) {
    v = paramm[ig][0] * std::pow(f,paramm[ig][1]);
  }
  else if (f <= limit3) {
    v = paramm[ig][2] * std::pow(f,paramm[ig][3]) /
      std::pow((1.0 - std::pow(f,paramm[ig][5])),paramm[ig][3]/paramm[ig][5]);
  }
  else {
    G4double l = -std::log(1.0-f);
    v = paramn[ig][6] +
        paramn[ig][7]*l +
        paramn[ig][8]*l*l +
        paramn[ig][9]*std::pow(l,3.0);
  }
  
  return v;
}
////////////////////////////////////////////////////////////////////////////////
//
// GetInverseValueN
//
/*G4double G4Type1GlauberParameterisation::GetInverseValueN (const G4int b,
  const G4double ppn) const
{
  G4int ig = 0;
  if (ppn < 1.0E-10) {
    return 0;
  }
  else {
    ig = G4int(2.0*std::log10(ppn)) - 2;
  }
  if (ig > 23) ig = 23; 
  
  return GetInverseValueN(b,ig);
}*/
////////////////////////////////////////////////////////////////////////////////
//
// GetInverseValueM
//
/*G4double G4Type1GlauberParameterisation::GetInverseValueM (const G4int b,
  const G4double ppn) const
{
  G4int ig = 0;
  if (ppn < 1.0E-10) {
    return 0;
  }
  else {
    ig = G4int(2.0*std::log10(ppn)) - 2;
  }
  if (ig > 23) ig = 23; 
  
  return GetInverseValueN(b,ig);
}*/
////////////////////////////////////////////////////////////////////////////////
//
// GetInverseValueN
//
/*G4double G4Type1GlauberParameterisation::GetInverseValueN (const G4int b,
  const G4int ig) const
{
  G4double v = 0.0;
  G4double x = (G4double) b;
  if (b <= 1) {
    v = 0.0;
  }
  else {
    G4double f1 = 1.0 + std::pow(paramn[ig][0]/x,mun1[ig]);
    G4double f2 = 1.0 + std::pow(cn[ig]/x,mun2[ig]);
    v = std::pow(f1*f2,1.0/paramn[ig][5]);
  }
  
  return v;
}*/
////////////////////////////////////////////////////////////////////////////////
//
// GetInverseValueM
//
/*G4double G4Type1GlauberParameterisation::GetInverseValueM (const G4int b,
  const G4int ig) const
{
  G4double v = 0.0;
  G4double x = (G4double) b;
  if (b <= 1) {
    v = 0.0;
  }
  else {
    G4double f1 = 1.0 + std::pow(paramm[ig][0]/x,mum1[ig]);
    G4double f2 = 1.0 + std::pow(cm[ig]/x,mum2[ig]);
    v = std::pow(f1*f2,1.0/paramm[ig][5]);
  }
  
  return v;
}*/
////////////////////////////////////////////////////////////////////////////////
//
#endif

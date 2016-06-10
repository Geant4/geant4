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
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4Pow.cc 69381 2013-05-02 09:58:14Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4Pow
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.05.2009
//
// Modifications:
// 08.01.2011 V.Ivanchenko extended maxZ from 256 to 512
//
// -------------------------------------------------------------------

#include "G4Pow.hh"

G4Pow* G4Pow::fpInstance = 0;

// -------------------------------------------------------------------

G4Pow* G4Pow::GetInstance()
{
  if (fpInstance == 0)
  {
    fpInstance = new G4Pow;
  }
  return fpInstance;
}

// -------------------------------------------------------------------

G4Pow::G4Pow()
  : onethird(1.0/3.0), minA(0.5000001), maxA(255.5)
{
  const G4int maxZ = 512; 
  const G4int maxZfact = 170; 

  pz13.resize(maxZ,0.0);
  lz.resize(maxZ,0.0);
  fact.resize(maxZfact,0.0);
  logfact.resize(maxZ,0.0);

  G4double f = 1.0;
  G4double logf = 0.0;
  fact[0] = 1.0;

  for(G4int i=1; i<maxZ; ++i)
  {
    G4double x  = G4double(i);
    pz13[i] = std::pow(x,onethird);
    lz[i]   = std::log(x);
    if(i < maxZfact)
    { 
      f *= x; 
      fact[i] = f;
    }
    logf += lz[i];
    logfact[i] = logf;
  }
}

// -------------------------------------------------------------------

G4Pow::~G4Pow()
{
  delete fpInstance; fpInstance = 0;
}

// -------------------------------------------------------------------

G4double G4Pow::powN(G4double x, G4int n)
{
  if(std::abs(n) > 8) { return std::pow(x, G4double(n)); }
  G4double res = 1.0;
  if(n >= 0) { for(G4int i=0; i<n; ++i) { res *= x; } }
  else if((n < 0) && (x != 0.0))
  {
    G4double y = 1.0/x;
    G4int nn = -n;
    for(G4int i=0; i<nn; ++i) { res *= y; }
  }
  return res;
}

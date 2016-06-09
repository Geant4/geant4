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
// $Id: G4Pow.cc,v 1.2 2009/07/03 11:35:06 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
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
// -------------------------------------------------------------------

#include "G4Pow.hh"

G4Pow* G4Pow::fInstance = 0;

// -------------------------------------------------------------------

G4Pow* G4Pow::GetInstance()
{
  if (fInstance == 0)
  {
    static G4Pow manager;
    fInstance = &manager;
  }
  return fInstance;
}

// -------------------------------------------------------------------

G4Pow::G4Pow()
{
  const G4int maxZ = 256; 

  pz13.resize(maxZ,0.0);
  lz.resize(maxZ,0.0);
  fact.resize(maxZ,0.0);

  onethird = 1.0/3.0;
  G4double f = 1.0;

  for(G4int i=1; i<maxZ; i++)
  {
    G4double x  = G4double(i);
    pz13[i] = std::pow(x,onethird);
    lz[i]   = std::log(x);
    f      *= x;
    fact[i] = f;
  }
  fact[0] = 1.0;
}

// -------------------------------------------------------------------

G4Pow::~G4Pow()
{}

// -------------------------------------------------------------------

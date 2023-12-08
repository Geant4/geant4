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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4BetaSpectrumSampler.cc                                          //
//  Author: D.H. Wright                                                       //
//  Date:   21 November 2022                                                  //
//  Description: samples a spectrum which is a piece-wise linear function of  //
//               energy. The CDF is calculated by trapezoidal integration     //
//               and within bins the line y = mx + b is sampled.              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4BetaSpectrumSampler.hh"
#include "Randomize.hh"

G4double G4BetaSpectrumSampler::shoot(const G4int npoints, const G4double* aCDF,
                                      const G4double estep)
{
  G4double prob = aCDF[npoints - 1]*G4UniformRand();
  G4int i = 0;
  for (; i<npoints; ++i) { if (prob <= aCDF[i]) { break; } }
  const G4double p1 = (i > 0) ? aCDF[i - 1] : aCDF[0];
  const G4double p2 = aCDF[i];
  const G4double delta = p2 - p1;
  const G4double x = (delta > 0.0) ? estep*i - estep*(p2 - prob)/delta : estep*i;
  return x;
}


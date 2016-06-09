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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// J. M. Quesada (12 October 2009) new implementation of Gamma function in configuration weight 
// J. M. Quesada (09 March 2010) Kappa is set to 6.
// 01.04.2011 General cleanup by V.Ivanchenko: integer Z and A, constructor
// 23.04.2011 V.Ivanchenko: make this class to be a simple container and no physics

#include "G4FermiConfiguration.hh"

G4FermiConfiguration::G4FermiConfiguration(const std::vector<const G4VFermiFragment*>& v)
{
  Configuration = v;
  totalA = totalZ = 0;
  totalMass = 0.0;
  size_t nn = v.size();
  for(size_t i=0; i<nn; ++i) {
    totalA += v[i]->GetA();
    totalZ += v[i]->GetZ();
    totalMass += v[i]->GetTotalEnergy();
  }
}

G4FermiConfiguration::~G4FermiConfiguration()
{}

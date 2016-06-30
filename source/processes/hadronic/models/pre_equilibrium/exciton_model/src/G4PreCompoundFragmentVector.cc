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
// $Id: G4PreCompoundFragmentVector.cc 96603 2016-04-25 13:29:51Z gcosmo $
//
// Hadronic Process: Nuclear Preequilibrium
// by V. Lara 
//
// Modified:
// 27.08.2010 V.Ivanchenko moved constructor and destructor to source, 
//            simplify run time computations making inlined 
//

#include "G4PreCompoundFragmentVector.hh"

G4PreCompoundFragmentVector::G4PreCompoundFragmentVector(pcfvector * avector) 
  : theChannels(nullptr), nChannels(0) 
{
  SetVector(avector);
}

G4PreCompoundFragmentVector::~G4PreCompoundFragmentVector()
{}

void G4PreCompoundFragmentVector::SetVector(pcfvector * avector)
{
  theChannels = avector;
  if(theChannels) {
    nChannels = theChannels->size();
    probabilities.resize(nChannels);
  }
}

//for inverse cross section choice
void G4PreCompoundFragmentVector::SetOPTxs(G4int opt)
{    
  for (G4int i=0; i< nChannels; ++i) { 
    (*theChannels)[i]->SetOPTxs(opt);
  }
}

//for superimposed Coulomb Barrier for inverse  cross sections 
void G4PreCompoundFragmentVector::UseSICB(G4bool use)
{    
  for (G4int i=0; i< nChannels; ++i) { 
    (*theChannels)[i]->UseSICB(use);
  }
}


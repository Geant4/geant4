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
// $Id: G4VPreCompoundTransitions.cc 68028 2013-03-13 13:48:15Z gcosmo $
//
// V.Ivanchenko 20.08.2010
//
 
#include "G4VPreCompoundTransitions.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4VPreCompoundTransitions::G4VPreCompoundTransitions()
  :useNGB(false),useCEMtr(false),  
   TransitionProb1(0.0), TransitionProb2(0.0), TransitionProb3(0.0)
{}

G4VPreCompoundTransitions::~G4VPreCompoundTransitions() 
{}

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
//
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4VCellScorer.cc
//
// ----------------------------------------------------------------------

#include "G4VCellScorer.hh"
#include "G4ios.hh"

G4VCellScorer::G4VCellScorer()
{
  static G4bool warn=true;
  if (warn)
  {
    G4cout
         << "--------------------------------------------------------" << G4endl
         << "WARNING: Class  <G4VCellScorer>  is  now obsolete.  It |" << G4endl
         << "         will be removed starting from the next Geant4 |" << G4endl
         << "         major release.  Please, consider switching to |" << G4endl
         << "         general purpose scoring functionality.        |" << G4endl
         << "--------------------------------------------------------"
         << G4endl;
    warn = false;
  }
}

G4VCellScorer::~G4VCellScorer()
{}

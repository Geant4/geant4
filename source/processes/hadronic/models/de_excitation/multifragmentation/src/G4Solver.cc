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
// $Id: G4Solver.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#include "G4Solver.hh"

template <class Function> 
G4Solver<Function>::G4Solver(const G4Solver & right)
{
    MaxIter = right.MaxIter;
    tolerance = right.tolerance;
    a = right.a;
    b = right.b;
    root = right.root;
}

// operators
template <class Function> 
G4Solver<Function> & G4Solver<Function>::operator=(const G4Solver & right)
{
    MaxIter = right.MaxIter;
    tolerance = right.tolerance;
    a = right.a;
    b = right.b;
    root = right.root;	
    return *this;
}

template <class Function> 
G4bool G4Solver<Function>::operator==(const G4Solver & right) const
{
    if (this == &right) return true;
    else return false;
}

template <class Function> 
G4bool G4Solver<Function>::operator!=(const G4Solver & right) const
{
    return !operator==(right);
}
	




//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Solver.cc,v 1.2 2005/06/04 13:27:48 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
	




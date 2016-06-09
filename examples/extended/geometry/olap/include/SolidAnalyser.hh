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
/// \file geometry/olap/include/SolidAnalyser.hh
/// \brief Definition of the SolidAnalyser class
//
//
// $Id$
//
// 
// --------------------------------------------------------------
// SolidAnalyser
//
// Sigleton providing Geant4 solids specifications in a generic way.
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef SolidAnalyser_h
#define SolidAnalyser_h

#include <vector>
#include <algorithm>

#include "globals.hh"

// supported Geant4 solids:
class G4VSolid;
class G4Box;
class G4Cons;
class G4Polycone;
class G4Polyhedra;
class G4Trap;
class G4Trd;
class G4Tubs;

class SolidAnalyser
{

public:

   static SolidAnalyser * GetSolidAnalyser();
   
   //void Reset();
   
   // user method to retrieve information
   G4int GetParam(const G4VSolid *,
                  std::vector<std::pair<G4String,G4double> > & ) const;
   
   ~SolidAnalyser();

protected:   
   G4int  GetParam(const G4Box *,
                   std::vector<std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Cons *,
                   std::vector<std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Polycone  *,
                   std::vector<std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Polyhedra *,
                   std::vector<std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Trap *,
                   std::vector<std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Trd *,
                   std::vector<std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Tubs *,
                   std::vector<std::pair<G4String,G4double> > & ) const;
   
private:

   SolidAnalyser();
   
   G4int NotImplemented(const G4VSolid *,
                       std::vector<std::pair<G4String,G4double> > & ) const;

   static SolidAnalyser * theInstance;
};

std::ostream & operator<<(std::ostream& flux,
                            std::vector<std::pair<G4String,G4double> >& v);

#endif

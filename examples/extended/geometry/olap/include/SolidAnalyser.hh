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
// $Id: SolidAnalyser.hh,v 1.1 2002-06-04 07:40:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "g4std/vector"
#include "g4std/algorithm"

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
                  G4std::vector<G4std::pair<G4String,G4double> > & ) const;
   
   ~SolidAnalyser();

protected:   
   G4int  GetParam(const G4Box *,
                   G4std::vector<G4std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Cons *,
                   G4std::vector<G4std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Polycone  *,
                   G4std::vector<G4std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Polyhedra *,
                   G4std::vector<G4std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Trap *,
                   G4std::vector<G4std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Trd *,
                   G4std::vector<G4std::pair<G4String,G4double> > & ) const;
   G4int  GetParam(const G4Tubs *,
                   G4std::vector<G4std::pair<G4String,G4double> > & ) const;
   
private:

   SolidAnalyser();
   
   G4int NotImplemented(const G4VSolid *,
                       G4std::vector<G4std::pair<G4String,G4double> > & ) const;

   static SolidAnalyser * theInstance;
};

G4std::ostream & operator<<(G4std::ostream& flux,
                            G4std::vector<G4std::pair<G4String,G4double> >& v);

#endif

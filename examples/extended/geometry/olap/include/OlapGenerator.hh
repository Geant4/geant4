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
// $Id: OlapGenerator.hh,v 1.1 2002-06-04 07:40:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapGenerator
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OlapGenerator_h
#define OlapGenerator_h

#include "g4std/vector"
#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"

class G4VisExtent;

class OlapGrid 
{
public:
   OlapGrid();
   void Reset() { for (G4int i=0;i<3;i++) count[i]=0; axis=0; }
   void Next();
   G4std::vector<G4int> count; // 
   G4std::vector<G4int> gsize; // gridsize
   G4int axis, eventsPerRun;
};


class OlapGenerator : public  G4VUserPrimaryGeneratorAction
{

public:
   OlapGenerator();
   ~OlapGenerator();

   void GeneratePrimaries(G4Event*);
   
   void SetExtent(const G4VisExtent&);
   
   void SetExtent(G4double);
   
   void SetAutoIncrement(G4bool b) { autoinc = b; }
   
   void SetGrid(G4int,G4int,G4int);
   
   void Reset();
   
   G4int GetAxis() { return grid.axis; };
   
   OlapGrid grid; // basically a vector of <int>, size==3
   G4std::vector<G4double> ext;
   G4ThreeVector posAB, posBA;
   G4bool autoinc;
};
#endif

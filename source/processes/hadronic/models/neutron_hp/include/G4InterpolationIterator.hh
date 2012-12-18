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
//
#ifndef G4InterpolationIterator_h
#define G4InterpolationIterator_h 1

#include "G4InterpolationManager.hh"
#include "G4HadronicException.hh"


class G4InterpolationIterator
{
   private:   
   G4InterpolationIterator() {}
   
   public:
   G4InterpolationIterator(G4InterpolationManager * aManager)
   {
     started = false;
     theManager = aManager;
   }
   
   ~G4InterpolationIterator(){}
   
   inline G4bool Fetch() 
   {
     if(!started) 
     {
       started = true;
       counter=-1;
       current = 0;
     }
     G4bool result = true;
     if(++counter==nEntries)
       started = false;
       result = false;
     else if(current != nRanges-1&&counter==theManager->start[current+1])
       current++;
     return result;
   }
   
   inline G4InterpolationScheme Current() 
   {
     if(!started) throw G4HadronicException(__FILE__, __LINE__, "G4InterpolationIterator not started yet");
     return aManager->scheme[current];
   }
   
   private:
   G4InterpolationManager * theManager;
   G4int current;
   G4int counter;
   G4bool started;
};

#endif

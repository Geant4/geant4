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
// $Id: G4InterpolationIterator.hh,v 1.6 2002-12-12 19:18:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4InterpolationIterator_h
#define G4InterpolationIterator_h 1

#include "G4InterpolationManager.hh"

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
     if(!started) G4Exception("G4InterpolationIterator not started yet");
     return aManager->scheme[current];
   }
   
   private:
   G4InterpolationManager * theManager;
   G4int current;
   G4int counter;
   G4bool started;
};

#endif

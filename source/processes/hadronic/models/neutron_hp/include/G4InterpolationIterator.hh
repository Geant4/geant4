// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4InterpolationIterator.hh,v 1.2 1999-06-29 18:43:42 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Here only for historical reasons.....HPW

#ifndef G4InterpolationIterator_h
#define G4InterpolationIterator_h 1

#include "G4InterpolationManager.hh"

class G4InterpolationIterator
{
   private:   
   G4InterpolationIterator();
   
   public:
   G4InterpolationIterator(G4InterpolationManager * aManager);
   
   ~G4InterpolationIterator();
   
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
     {
       started = false;
       result = false;
     }
     else if(current != nRanges-1&&counter==theManager->start[current+1])
       current++;
     return result;
   }
   
   inline G4InterpolationScheme Current() 
   {
     if(!started) G4Exception("G4InterpolationIterator not started yet");
     return theManager->scheme[current];
   }
   
   private:
   G4InterpolationManager * theManager;
   G4int current;
   G4int nEntries;
   G4int nRanges;
   G4int counter;
   G4bool started;
};

#endif

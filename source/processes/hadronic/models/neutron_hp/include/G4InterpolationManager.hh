// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4InterpolationManager.hh,v 1.5 1999-12-15 14:53:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4InterpolationManager_h
#define G4InterpolationManager_h 1

#include "globals.hh"
#include "G4InterpolationScheme.hh"
#include "G4ios.hh"
#include "g4std/fstream"

class G4InterpolationManager
{
   public:
   
   friend class G4InterpolationIterator;
   
   G4InterpolationManager()
   {
     nRanges = 1;
     start = new G4int[1];
     start[0] = 0;
     range = new G4int[1];
     range [0] = 100000;
     scheme = new G4InterpolationScheme[1];
     scheme[0] = LINLIN;
     nEntries = 0;
   }
   
   ~G4InterpolationManager()
   {
     if(start!=NULL) delete [] start;
     if(range!=NULL) delete [] range;
     if(scheme!=NULL) delete [] scheme;
   }
   
   G4InterpolationManager & operator= (const G4InterpolationManager & aManager)
   {
     if(this != &aManager)
     {
       nRanges = aManager.nRanges;
       nEntries = aManager.nEntries;
       if(scheme!=NULL) delete [] scheme;
       if(start!=NULL) delete [] start;
       if(range!=NULL) delete [] range;
       scheme = new G4InterpolationScheme[nRanges];
       start = new G4int[nRanges];
       range = new G4int[nRanges];
       for(G4int i=0; i<nRanges; i++)
       {
	 scheme[i]=aManager.scheme[i];
	 start[i]=aManager.start[i];
	 range[i]=aManager.range[i];
       }
     }
     return *this;
   }
   
   inline void Init(G4int aScheme, G4int aRange)
   {
     nRanges = 1;
     start[0] = 0;
     range [0] = aRange;
     scheme[0] = MakeScheme(aScheme);
     nEntries = aRange;
   }
   inline void Init(G4InterpolationScheme aScheme, G4int aRange)
   {
     nRanges = 1;
     start[0] = 0;
     range [0] = aRange;
     scheme[0] = aScheme;
     nEntries = aRange;
   }
   
   inline void Init(G4std::ifstream & aDataFile)
   {
     delete [] start;
     delete [] range;
     delete [] scheme;
     aDataFile >> nRanges;
     start = new G4int[nRanges];
     range = new G4int[nRanges];
     scheme = new G4InterpolationScheme[nRanges];
     start[0] = 0;
     G4int it;
     for(G4int i=0; i<nRanges; i++)
     {
       aDataFile>>range[i];
       if(i!=0) start[i] = start[i-1]+range[i-1];
       aDataFile>>it;
       scheme[i] = MakeScheme(it);
     }
     nEntries = start[nRanges-1]+range[nRanges-1];
   }
   
   G4InterpolationScheme MakeScheme(G4int it);
   
   inline G4InterpolationScheme GetScheme(G4int index) const
   {
     G4int it = 0;
     for(G4int i=1; i<nRanges; i++)
     {
       if(index<start[i]) break;
       it = i;
     }
     return scheme[it];
   }
   
   inline void CleanUp()
   {
     nRanges = 0;
     nEntries = 0;
   }
   
   inline G4InterpolationScheme GetInverseScheme(G4int index)
   {
     G4InterpolationScheme result = GetScheme(index);
     if(result == HISTO) 
     {
       result = RANDOM;
     }
     else if(result == LINLOG)  
     {
       result = LOGLIN;
     }
     else if(result == LOGLIN)  
     {
       result = LINLOG;
     }
     else if(result == CHISTO) 
     {
       result = CRANDOM;
     }
     else if(result == CLINLOG) 
     {
       result = CLOGLIN;
     }
     else if(result == CLOGLIN) 
     {
       result = CLINLOG;
     }
     else if(result == UHISTO) 
     {
       result = URANDOM;
     }
     else if(result == ULINLOG) 
     {
       result = ULOGLIN;
     }
     else if(result == ULOGLIN) 
     {
       result = ULINLOG;
     }
     return result;
   }
   
   void AppendScheme(G4int aPoint, const G4InterpolationScheme & aScheme);
   
   private:
   
   G4int nRanges;
   G4InterpolationScheme * scheme;
   G4int * start;
   G4int * range;  
   G4int nEntries;
    
};
#endif

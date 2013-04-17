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
#ifndef G4InterpolationManager_h
#define G4InterpolationManager_h 1

#include "globals.hh"
#include "G4InterpolationScheme.hh"
#include "G4ios.hh"
#include <fstream>
#include "G4HadronicException.hh"

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
     if(start!=0) delete [] start;
     if(range!=0) delete [] range;
     if(scheme!=0) delete [] scheme;
   }
   
   G4InterpolationManager & operator= (const G4InterpolationManager & aManager)
   {
     if(this != &aManager)
     {
       nRanges = aManager.nRanges;
       nEntries = aManager.nEntries;
       if(scheme!=0) delete [] scheme;
       if(start!=0) delete [] start;
       if(range!=0) delete [] range;
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
   
   inline void Init(std::istream & aDataFile)
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
       //***************************************************************
       //EMendoza -> there is a bug here.
       /*
       if(i!=0) start[i] = start[i-1]+range[i-1];
       */
       //***************************************************************
       if(i!=0) start[i] = range[i-1];
       //***************************************************************
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

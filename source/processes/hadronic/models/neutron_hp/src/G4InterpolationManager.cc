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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4InterpolationManager.hh"

   G4InterpolationScheme G4InterpolationManager::MakeScheme(G4int it)
   {
     G4InterpolationScheme result(LINLIN);
     switch(it)
     {
      case 1:
        result = HISTO;
        break;
      case 2:
        result = LINLIN;
        break;
      case 3:
        result = LINLOG;
        break;
      case 4:
        result = LOGLIN;
        break;
      case 5:
        result = LOGLOG;
        break;
      case 11:
        result = CHISTO;
        break;
      case 12:
        result = CLINLIN;
        break;
      case 13:
        result = CLINLOG;
        break;
      case 14:
        result = CLOGLIN;
        break;
      case 15:
        result = CLOGLOG;
        break;
      case 21:
        result = UHISTO;
        break;
      case 22:
        result = ULINLIN;
        break;
      case 23:
        result = ULINLOG;
        break;
      case 24:
        result = ULOGLIN;
        break;
      case 25:
        result = ULOGLOG;
        break;
      default:
        G4Exception("G4InterpolationManager: unknown interpolation scheme");
        break;        
     }
     return result;
   }

   void G4InterpolationManager::AppendScheme(G4int aPoint, const G4InterpolationScheme & aScheme)
   {
     if(aPoint!=nEntries) 
     {
       G4cout <<"G4InterpolationManager::AppendScheme - "<<aPoint<<" "<<nEntries<<G4endl;
       G4Exception("Wrong usage of G4InterpolationManager::AppendScheme");
     }
     if(nEntries==0)
     {
       nEntries = 1;
       nRanges = 1;
       start[0] = 0;
       range [0] = 1;
       scheme[0] = aScheme;
     }
     else if(aScheme==scheme[nRanges-1])
     {
       ++range[nRanges-1];
       nEntries++;
     }
     else
     {
       nEntries++;
       nRanges++;
       G4int i;
       G4int * buffer = new G4int[nRanges];
       G4int * buffer1 = new G4int[nRanges];
       G4InterpolationScheme* buff2 = new G4InterpolationScheme[nRanges];
       for(i=0; i<nRanges-1; i++) 
       {
         buffer[i]  = start[i];
         buffer1[i] = range[i];
         buff2[i]   = scheme[i];
       }
       delete [] start;
       delete [] range;
       delete [] scheme;
       start  = buffer;
       range  = buffer1;
       scheme = buff2;
       start[nRanges-1]  = start[nRanges-2]+range[nRanges-2];
       range[nRanges-1]  = 1;
       scheme[nRanges-1] = aScheme;
     }
   }

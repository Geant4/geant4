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
// $Id: G4HEKaonZeroShortInelastic.hh,v 1.9 2002-12-12 19:17:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Gheisha High Energy model class -- header file
// H. Fesefeldt, RWTH Aachen 23-October-1996
// Last modified: 10-December-1996

// A prototype of the Gheisha High Energy collision model.

#ifndef G4HEKaonZeroShortInelastic_h
#define G4HEKaonZeroShortInelastic_h 1

#include "G4HEKaonZeroInelastic.hh"
#include "G4HEAntiKaonZeroInelastic.hh"

class G4HEKaonZeroShortInelastic : public G4HEInelastic  
{
 public: 
        G4HEKaonZeroShortInelastic() 
           {
              theMinEnergy =  20*GeV;
              theMaxEnergy = 10*TeV;
              MAXPART      = 2048;
              verboseLevel = 0; 
           }

        ~G4HEKaonZeroShortInelastic(){ };
         
        G4int vecLength;

        void SetMaxNumberOfSecondaries(G4int maxnumber)
             { MAXPART = maxnumber;};
        void SetVerboseLevel(G4int verbose)
             { verboseLevel = verbose;};
        G4int GetNumberOfSecondaries()
             { return vecLength;};           

        G4VParticleChange * ApplyYourself( const G4Track &aTrack, G4Nucleus &targetNucleus );

};
#endif                     
                                         


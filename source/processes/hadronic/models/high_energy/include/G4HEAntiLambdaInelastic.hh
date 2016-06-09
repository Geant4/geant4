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
// $Id: G4HEAntiLambdaInelastic.hh,v 1.10 2003/07/01 15:42:21 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
//
// G4 Gheisha High Energy model class -- header file
// H. Fesefeldt, RWTH Aachen 23-October-1996
// Last modified: 10-December-1996

// A prototype of the Gheisha High Energy collision model.

#ifndef G4HEAntiLambdaInelastic_h
#define G4HEAntiLambdaInelastic_h 1

#include "G4HEInelastic.hh"

class G4HEAntiLambdaInelastic : public G4HEInelastic  
{
 public: 
        G4HEAntiLambdaInelastic() : G4HEInelastic()
           {
              theMinEnergy =  20*GeV;
              theMaxEnergy = 10*TeV;
              MAXPART      = 2048;
              verboseLevel = 0; 
           }

        ~G4HEAntiLambdaInelastic(){ };
         
        G4int vecLength;
        
        G4HadFinalState * ApplyYourself( const G4HadProjectile &aTrack, G4Nucleus &targetNucleus );

        G4int  GetNumberOfSecondaries()
               { return vecLength; }         

        void   FirstIntInCasAntiLambda(G4bool &inElastic, const G4double availableEnergy,
                                       G4HEVector pv[],
                                       G4int &vecLen, 
                                       G4HEVector incidentParticle,
                                       G4HEVector targetParticle,
                                       const G4double atomicWeight);
};
#endif                     
                                         


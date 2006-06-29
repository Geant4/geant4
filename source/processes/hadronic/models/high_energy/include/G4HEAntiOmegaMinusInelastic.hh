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
// $Id: G4HEAntiOmegaMinusInelastic.hh,v 1.13 2006-06-29 20:28:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Gheisha High Energy model class -- header file
// H. Fesefeldt, RWTH Aachen 23-October-1996
// Last modified: 10-December-1996

// A prototype of the Gheisha High Energy collision model.

#ifndef G4HEAntiOmegaMinusInelastic_h
#define G4HEAntiOmegaMinusInelastic_h 1

// Class description:
// High energy parameterized model for anti-Omega- inelastic scattering.  This
// class is responsible for producing the final state of the interaction and
// is typically valid for incident anti-Omega- energies above 20 GeV.  This
// physics may be invoked by registering an instance of the class with
// G4AntiOmegaMinusInelasticProcess in the user's physics list.
//
// This class is derived from G4HEInelastic which in turn is derived from
// G4HadronicInteraction.

// Class Description - End

#include "G4HEInelastic.hh"

class G4HEAntiOmegaMinusInelastic : public G4HEInelastic  
{
 public:  // with description
        G4HEAntiOmegaMinusInelastic() : G4HEInelastic()
           {
              theMinEnergy =  20*GeV;
              theMaxEnergy = 10*TeV;
              MAXPART      = 2048;
              verboseLevel = 0; 
           }

        ~G4HEAntiOmegaMinusInelastic(){ };
         
        G4int vecLength;
        
        G4HadFinalState * ApplyYourself( const G4HadProjectile &aTrack, G4Nucleus &targetNucleus );

        G4int  GetNumberOfSecondaries()
               { return vecLength; }         

        void   FirstIntInCasAntiOmegaMinus(G4bool &inElastic, const G4double availableEnergy,
                                           G4HEVector pv[],
                                           G4int &vecLen, 
                                           G4HEVector incidentParticle,
                                           G4HEVector targetParticle,
                                           const G4double atomicWeight);
};
#endif                     
                                         


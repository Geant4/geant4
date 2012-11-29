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
// $Id$
//
//
// G4 Gheisha High Energy model class -- header file
// H. Fesefeldt, RWTH Aachen 23-October-1996
// A prototype of the Gheisha High Energy collision model.

#ifndef G4HESigmaZeroInelastic_h
#define G4HESigmaZeroInelastic_h 1

// Class description:
// High energy parameterized model for Sigma0 inelastic scattering.  This
// class is responsible for producing the final state of the interaction and
// is typically valid for incident Sigma0 energies above 20 GeV.  There is
// currently no corresponding inelastic process to which this model can 
// be assigned.
//
// This class is derived from G4HEInelastic which in turn is derived from
// G4HadronicInteraction.

// Class Description - End

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4HELambdaInelastic.hh"

class G4HESigmaZeroInelastic : public G4HEInelastic  
{
  public:  // with description 
    G4HESigmaZeroInelastic() : G4HEInelastic("G4HESigmaZeroInelastic")
    {
      theMinEnergy = 20*CLHEP::GeV;
      theMaxEnergy = 10*CLHEP::TeV;
      MAXPART      = 2048;
      verboseLevel = 0; 
      G4cout << "WARNING: model G4HESigmaZeroInelastic is being deprecated and will\n"
             << "disappear in Geant4 version 10.0"  << G4endl; 
    }

    ~G4HESigmaZeroInelastic() {};

    virtual void ModelDescription(std::ostream&) const;
         
    G4int verboseLevel;
    G4int MAXPART;
    G4int vecLength;

    void SetMaxNumberOfSecondaries(G4int maxnumber)
        { MAXPART = maxnumber; };
    void SetVerboseLevel(G4int verbose)
        { verboseLevel = verbose;};
 
    G4HadFinalState* ApplyYourself(const G4HadProjectile &aTrack, G4Nucleus &targetNucleus);

    G4int GetNumberOfSecondaries()
         { return vecLength; }         
};
#endif                     
                                         


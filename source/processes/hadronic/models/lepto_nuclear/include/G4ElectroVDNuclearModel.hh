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
// Author:  D.H. Wright
// Date:    1 May 2012
//
// Description: model for electron and positron interaction with nuclei
//              using the equivalent photon spectrum.  A real gamma is 
//              produced from the virtual photon spectrum and is then 
//              interacted hadronically by the Bertini cascade at low
//              energies.  At high energies the gamma is converted to a
//              pi0 and interacted with the FTFP model.

#ifndef G4ElectroVDNuclearModel_h
#define G4ElectroVDNuclearModel_h 1

#include "G4HadronicInteraction.hh"

class G4CascadeInterface;
class G4TheoFSGenerator;
class G4LundStringFragmentation;
class G4ExcitedStringDecay;
class G4ElectroNuclearCrossSection;
class G4VCrossSectionDataSet;

class G4ElectroVDNuclearModel : public G4HadronicInteraction
{
  public: 

    G4ElectroVDNuclearModel();
    ~G4ElectroVDNuclearModel();
    
    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                   G4Nucleus& aTargetNucleus);

    virtual void ModelDescription(std::ostream& outFile) const;

  private:

    G4DynamicParticle* CalculateEMVertex(const G4HadProjectile& aTrack,
                                         G4Nucleus& aTargetNucleus);

    void CalculateHadronicVertex(G4DynamicParticle* incident,
                                 G4Nucleus& target);

    G4double leptonKE;
    G4double photonEnergy;
    G4double photonQ2;

    G4ElectroNuclearCrossSection* electroXS;
    G4VCrossSectionDataSet* gammaXS;

    G4TheoFSGenerator* ftfp;
    G4LundStringFragmentation* theFragmentation;
    G4ExcitedStringDecay* theStringDecay;
    G4CascadeInterface* bert;

    G4int secID;  // Creator model ID for the secondaries created by this model
};

#endif

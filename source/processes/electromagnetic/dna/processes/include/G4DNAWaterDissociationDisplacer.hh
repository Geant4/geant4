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
// $Id: G4DNAWaterDissociationDisplacer.hh 94010 2015-11-05 10:08:33Z gcosmo $
//
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 


#ifndef G4DNAWaterDissociationDisplacer_h
#define G4DNAWaterDissociationDisplacer_h 1

#include "G4VMolecularDissociationDisplacer.hh"
#include "G4DNARevertProbability.hh"

class G4DNAWaterDissociationDisplacer : public G4VMolecularDecayDisplacer
{
public :
    G4DNAWaterDissociationDisplacer();
    virtual ~G4DNAWaterDissociationDisplacer() ;

    virtual std::vector<G4ThreeVector> GetProductsDisplacement(const G4MolecularDissociationChannel*) const;
    virtual G4ThreeVector GetMotherMoleculeDisplacement(const G4MolecularDissociationChannel*) const;
    G4ThreeVector radialDistributionOfElectron() const;
    G4ThreeVector radialDistributionOfProducts(G4double r_rms) const;
    static G4double ElectronProbaDistribution(G4double r);

#if defined G4EM_ALLOC_EXPORT
    G4DLLEXPORT static const DisplacementType Ionisation_DissociationDecay;
    G4DLLEXPORT static const DisplacementType A1B1_DissociationDecay;
    G4DLLEXPORT static const DisplacementType B1A1_DissociationDecay;
    G4DLLEXPORT static const DisplacementType AutoIonisation;
    G4DLLEXPORT static const DisplacementType DissociativeAttachment;
#else
    G4DLLIMPORT static const DisplacementType Ionisation_DissociationDecay;
    G4DLLIMPORT static const DisplacementType A1B1_DissociationDecay;
    G4DLLIMPORT static const DisplacementType B1A1_DissociationDecay;
    G4DLLIMPORT static const DisplacementType AutoIonisation;
    G4DLLIMPORT static const DisplacementType DissociativeAttachment;
#endif

private :
    std::function<G4double(G4double)> fProba1DFunction;
    std::vector<G4double> fElectronThermalization;
    G4DNARevertProbability fFastElectronDistrib;
};
#endif


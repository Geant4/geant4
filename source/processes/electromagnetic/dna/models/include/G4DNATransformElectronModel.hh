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
// $Id: G4DNATransformElectronModel.hh 85244 2014-10-27 08:24:13Z gcosmo $
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


#ifndef G4DNATransformElectronModel_h
#define G4DNATransformElectronModel_h 1

#include "G4VEmModel.hh"


/**
  * When an electron reaches the highest energy domain of G4DNATransformElectronModel,
  * it is then automatically converted into a solvated electron without thermalization
  * displacement (assumed to be already thermalized).
  */

class G4DNATransformElectronModel: public G4VEmModel
{
public :
    G4DNATransformElectronModel(const G4ParticleDefinition* p = 0,
                                const G4String& nam = "DNATransformElectronModel");
    virtual ~G4DNATransformElectronModel();

    virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

    virtual G4double CrossSectionPerVolume(  const G4Material* material,
            const G4ParticleDefinition* p,
            G4double ekin,
            G4double emin,
            G4double emax);

    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                   const G4MaterialCutsCouple*,
                                   const G4DynamicParticle*,
                                   G4double tmin,
                                   G4double maxEnergy);

    inline void SetVerbose(int);

    inline void SetEpsilonEnergy(G4double);

protected:

    G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

    // Water density table
    const std::vector<G4double>* fpWaterDensity;

    G4bool fIsInitialised;
    G4int fVerboseLevel;
    G4double fEpsilon;

    G4DNATransformElectronModel & operator=(const  G4DNATransformElectronModel &right);
    G4DNATransformElectronModel(const  G4DNATransformElectronModel&);

};

inline void G4DNATransformElectronModel::SetVerbose(int flag)
{
    fVerboseLevel = flag ;
}

inline void G4DNATransformElectronModel::SetEpsilonEnergy(G4double eps)
{
    fEpsilon = eps ;
}

#endif

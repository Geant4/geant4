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
// $Id: $
//
// Author:  D.H. Wright
// Date:    1 May 2012
//
// Description: model for electron and positron interaction with nuclei
//              using the equivalent photon spectrum.  A real gamma is 
//              produced from the virtual photon spectrum and is then 
//              interacted hadronically by the Bertini cascade at low
//              energies.  At high energies the gamma is treated as a 
//              pi0 and interacted with the nucleus using the FTFP model.
//              The electro- and photo-nuclear cross sections of
//              M. Kossov are used to generate the virtual photon
//              spectrum.
//

#include "G4ElectroVDNuclearModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ElectroNuclearCrossSection.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CascadeInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"

#include "G4HadFinalState.hh"
#include "G4HadronicInteractionRegistry.hh"

G4ElectroVDNuclearModel::G4ElectroVDNuclearModel()
 : G4HadronicInteraction("G4ElectroVDNuclearModel"),
   leptonKE(0.0), photonEnergy(0.0), photonQ2(0.0)
{
  SetMinEnergy(0.0);
  SetMaxEnergy(1*PeV);
  electroXS = 
    (G4ElectroNuclearCrossSection*)G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4ElectroNuclearCrossSection::Default_Name());
  gammaXS = 
    (G4PhotoNuclearCrossSection*)G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4PhotoNuclearCrossSection::Default_Name());

  // reuse existing pre-compound model
  G4GeneratorPrecompoundInterface* precoInterface 
    = new G4GeneratorPrecompoundInterface();
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4VPreCompoundModel* pre = static_cast<G4VPreCompoundModel*>(p);
  if(!pre) { pre = new G4PreCompoundModel(); }
  precoInterface->SetDeExcitation(pre);

  // string model
  ftfp = new G4TheoFSGenerator();
  ftfp->SetTransport(precoInterface);
  theFragmentation = new G4LundStringFragmentation();
  theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  G4FTFModel* theStringModel = new G4FTFModel();
  theStringModel->SetFragmentationModel(theStringDecay);
  ftfp->SetHighEnergyGenerator(theStringModel);
    
  // Build Bertini model
  bert = new G4CascadeInterface();
}

G4ElectroVDNuclearModel::~G4ElectroVDNuclearModel()
{
  delete theFragmentation;
  delete theStringDecay;
}
    
void G4ElectroVDNuclearModel::ModelDescription(std::ostream& outFile) const 
{
  outFile << "G4ElectroVDNuclearModel handles the inelastic scattering\n"
          << "of e- and e+ from nuclei using the equivalent photon\n"
          << "approximation in which the incoming lepton generates a\n"
          << "virtual photon at the electromagnetic vertex, and the\n"
          << "virtual photon is converted to a real photon.  At low\n"
          << "energies, the photon interacts directly with the nucleus\n"
          << "using the Bertini cascade.  At high energies the photon\n"
          << "is converted to a pi0 which interacts using the FTFP\n"
          << "model.  The electro- and gamma-nuclear cross sections of\n"
          << "M. Kossov are used to generate the virtual photon spectrum\n";
}


G4HadFinalState*
G4ElectroVDNuclearModel::ApplyYourself(const G4HadProjectile& aTrack,
                                       G4Nucleus& targetNucleus)
{
    // Set up default particle change (just returns initial state)
    theParticleChange.Clear();
    theParticleChange.SetStatusChange(isAlive);
    leptonKE = aTrack.GetKineticEnergy();
    theParticleChange.SetEnergyChange(leptonKE);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit() );
    
    // Set up sanity checks for real photon production
    G4DynamicParticle lepton(aTrack.GetDefinition(), aTrack.Get4Momentum() );
    
    // Need to call GetElementCrossSection before calling GetEquivalentPhotonEnergy.
    G4Material* mat = 0;
    G4int targZ = targetNucleus.GetZ_asInt();
    electroXS->GetElementCrossSection(&lepton, targZ, mat);
    
    photonEnergy = electroXS->GetEquivalentPhotonEnergy();
    // Photon energy cannot exceed lepton energy
    if (photonEnergy < leptonKE) {
        photonQ2 = electroXS->GetEquivalentPhotonQ2(photonEnergy);
        G4double dM = G4Proton::Proton()->GetPDGMass() + G4Neutron::Neutron()->GetPDGMass();
        // Photon
        if (photonEnergy > photonQ2/dM) {
            // Produce recoil lepton and transferred photon
            G4DynamicParticle* transferredPhoton = CalculateEMVertex(aTrack, targetNucleus);
            // Interact gamma with nucleus
            if (transferredPhoton) CalculateHadronicVertex(transferredPhoton, targetNucleus);
        }
    }
    return &theParticleChange;
}


G4DynamicParticle*
G4ElectroVDNuclearModel::CalculateEMVertex(const G4HadProjectile& aTrack,
                                           G4Nucleus& targetNucleus)
{
  G4DynamicParticle photon(G4Gamma::Gamma(), photonEnergy,
                           G4ThreeVector(0.,0.,1.) );

  // Get gamma cross section at Q**2 = 0 (real gamma)
  G4int targZ = targetNucleus.GetZ_asInt();
  G4Material* mat = 0;
  G4double sigNu =
    gammaXS->GetElementCrossSection(&photon, targZ, mat);

  // Change real gamma energy to equivalent energy and get cross section at that energy 
  G4double dM = G4Proton::Proton()->GetPDGMass() + G4Neutron::Neutron()->GetPDGMass();
  photon.SetKineticEnergy(photonEnergy - photonQ2/dM);      
  G4double sigK =
    gammaXS->GetElementCrossSection(&photon, targZ, mat);
  G4double rndFraction = electroXS->GetVirtualFactor(photonEnergy, photonQ2);

  // No gamma produced, return null ptr
  if (sigNu*G4UniformRand() > sigK*rndFraction) return 0;

  // Scatter the lepton
  G4double mProj = aTrack.GetDefinition()->GetPDGMass();
  G4double mProj2 = mProj*mProj;
  G4double iniE = leptonKE + mProj;               // Total energy of incident lepton
  G4double finE = iniE - photonEnergy;            // Total energy of scattered lepton
  theParticleChange.SetEnergyChange(finE-mProj);
  G4double iniP = std::sqrt(iniE*iniE-mProj2);    // Incident lepton momentum
  G4double finP = std::sqrt(finE*finE-mProj2);    // Scattered lepton momentum
  G4double cost = (iniE*finE - mProj2 - photonQ2/2.)/iniP/finP;  // cos(theta) from Q**2
  if (cost > 1.) cost= 1.;
  if (cost < -1.) cost=-1.;
  G4double sint = std::sqrt(1.-cost*cost);

  G4ThreeVector dir = aTrack.Get4Momentum().vect().unit();
  G4ThreeVector ortx = dir.orthogonal().unit();   // Ortho-normal to scattering plane
  G4ThreeVector orty = dir.cross(ortx);           // Third unit vector
  G4double phi = twopi*G4UniformRand();
  G4double sinx = sint*std::sin(phi);
  G4double siny = sint*std::cos(phi);
  G4ThreeVector findir = cost*dir+sinx*ortx+siny*orty;
  theParticleChange.SetMomentumChange(findir);    // change lepton direction

  // Create a gamma with momentum equal to momentum transfer
  G4ThreeVector photonMomentum = iniP*dir - finP*findir;
  G4DynamicParticle* gamma = new G4DynamicParticle(G4Gamma::Gamma(),
                                                   photonEnergy, photonMomentum);
  return gamma;
}


void
G4ElectroVDNuclearModel::CalculateHadronicVertex(G4DynamicParticle* incident,
                                                 G4Nucleus& target)
{
  G4HadFinalState* hfs = 0;
  G4double gammaE = incident->GetTotalEnergy();

  if (gammaE < 10*GeV) {
    G4HadProjectile projectile(*incident);
    hfs = bert->ApplyYourself(projectile, target);
  } else {
    // At high energies convert incident gamma to a pion
    G4double piMass = G4PionZero::PionZero()->GetPDGMass();
    G4double piMom = std::sqrt(gammaE*gammaE - piMass*piMass);
    G4ThreeVector piMomentum(incident->GetMomentumDirection() );
    piMomentum *= piMom;
    G4DynamicParticle theHadron(G4PionZero::PionZero(), piMomentum);
    G4HadProjectile projectile(theHadron);
    hfs = ftfp->ApplyYourself(projectile, target);
  }

  delete incident;

  // Copy secondaries from sub-model to model
  theParticleChange.AddSecondaries(hfs);
}


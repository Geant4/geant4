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
//

#include "Tst25DetectorConstruction.hh"
#include "Tst25RunAction.hh"
#include "Tst25PrimaryGeneratorAction.hh"
#include "Tst25PhysicsList.hh"
#include "Tst25SteppingAction.hh"
#include "Tst25StackingAction.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

// #include "../src/G4HadronKineticModel.cc"
// #include "../src/G4Absorber.cc"
// #include "../src/G4AntiProtonField.cc"
// #include "../src/G4ConstDistribution.cc"
// #include "../src/G4ExpDistribution.cc"
// #include "../src/G4FieldPropagation.cc"
// #include "../src/G4GeneratorPrecompoundInterface.cc"
// #include "../src/G4IsotropicDistribution.cc"
// #include "../src/G4KM_NuclearFermiDensity.cc"
// #include "../src/G4KM_NuclearShellModelDensity.cc"
// #include "../src/G4KM_NucleonEqRhs.cc"
// #include "../src/G4KM_OpticalEqRhs.cc"
// #include "../src/G4KaonMinusField.cc"
// #include "../src/G4KaonPlusField.cc"
// #include "../src/G4KaonZeroField.cc"
// #include "../src/G4NeutronField.cc"
// #include "../src/G4PionMinusField.cc"
// #include "../src/G4PionPlusField.cc"
// #include "../src/G4PionZeroField.cc"
// #include "../src/G4ProtonField.cc"
// #include "../src/G4RKFieldIntegrator.cc"
// #include "../src/G4RKPropagation.cc"
// #include "../src/G4ScattererStub.cc"
// #include "../src/G4SigmaMinusField.cc"
// #include "../src/G4SigmaPlusField.cc"
// #include "../src/G4SigmaZeroField.cc"
// #include "../src/G4VFieldPropagation.cc"
// #include "../src/G4VNuclearField.cc"
// #include "../src/G4XNull.cc"
// 
// 
// // scattering
// #include "../src/G4AngularDistribution.cc"
// #include "../src/G4AngularDistributionNP.cc"
// #include "../src/G4AngularDistributionPP.cc"
// #include "../src/G4BaryonPartialWidth.cc"
// #include "../src/G4BaryonWidth.cc"
// #include "../src/G4Clebsch.cc"
// #include "../src/G4CollisionComposite.cc"
// #include "../src/G4CollisionDeltaN.cc"
// #include "../src/G4CollisionDeltaNToNN.cc"
// #include "../src/G4CollisionDeltaStarN.cc"
// #include "../src/G4CollisionDeltaStarNToNN.cc"
// #include "../src/G4CollisionInitialState.cc"
// #include "../src/G4CollisionManager.cc"
// #include "../src/G4CollisionMesonBaryon.cc"
// #include "../src/G4CollisionMesonBaryonElastic.cc"
// #include "../src/G4CollisionMesonBaryonToResonance.cc"
// #include "../src/G4CollisionNN.cc"
// #include "../src/G4CollisionNNElastic.cc"
// #include "../src/G4CollisionNNToDeltaDelta.cc"
// #include "../src/G4CollisionNNToDeltaDelta1600.cc"
// #include "../src/G4CollisionNNToDeltaDelta1620.cc"
// #include "../src/G4CollisionNNToDeltaDelta1700.cc"
// #include "../src/G4CollisionNNToDeltaDelta1900.cc"
// #include "../src/G4CollisionNNToDeltaDelta1905.cc"
// #include "../src/G4CollisionNNToDeltaDelta1910.cc"
// #include "../src/G4CollisionNNToDeltaDelta1920.cc"
// #include "../src/G4CollisionNNToDeltaDelta1930.cc"
// #include "../src/G4CollisionNNToDeltaDelta1950.cc"
// #include "../src/G4CollisionNNToDeltaDeltastar.cc"
// #include "../src/G4CollisionNNToDeltaNstar.cc"
// #include "../src/G4CollisionNNToNDelta.cc"
// #include "../src/G4CollisionNNToNDelta1600.cc"
// #include "../src/G4CollisionNNToNDelta1620.cc"
// #include "../src/G4CollisionNNToNDelta1700.cc"
// #include "../src/G4CollisionNNToNDelta1900.cc"
// #include "../src/G4CollisionNNToNDelta1905.cc"
// #include "../src/G4CollisionNNToNDelta1910.cc"
// #include "../src/G4CollisionNNToNDelta1920.cc"
// #include "../src/G4CollisionNNToNDelta1930.cc"
// #include "../src/G4CollisionNNToNDelta1950.cc"
// #include "../src/G4CollisionNNToNDeltastar.cc"
// #include "../src/G4CollisionNNToNNstar.cc"
// #include "../src/G4CollisionNNToStarStar.cc"
// #include "../src/G4CollisionNStarN.cc"
// #include "../src/G4CollisionNStarNToNN.cc"
// #include "../src/G4CollisionPN.cc"
// #include "../src/G4CollisionPtr.cc"
// #include "../src/G4CollisionnpElastic.cc"
// #include "../src/G4ConcreteDeltaStarNToNN.cc"
// #include "../src/G4ConcreteMesonBaryonToResonance.cc"
// #include "../src/G4ConcreteNNToDeltaDelta.cc"
// #include "../src/G4ConcreteNNToDeltaDeltastar.cc"
// #include "../src/G4ConcreteNNToDeltaNstar.cc"
// #include "../src/G4ConcreteNNToNDelta.cc"
// #include "../src/G4ConcreteNNToNDeltaStar.cc"
// #include "../src/G4ConcreteNNToNNStar.cc"
// #include "../src/G4ConcreteNNTwoBodyResonance.cc"
// #include "../src/G4ConcreteNStarNToNN.cc"
// #include "../src/G4CrossSectionComposite.cc"
// #include "../src/G4CrossSectionPatch.cc"
// #include "../src/G4CrossSectionSourcePtr.cc"
// #include "../src/G4DetailedBalancePhaseSpaceIntegral.cc"
// #include "../src/G4GeneralElasticCollision.cc"
// #include "../src/G4LowEXsection.cc"
// #include "../src/G4MesonPartialWidth.cc"
// #include "../src/G4MesonWidth.cc"
// #include "../src/G4PartialWidthTable.cc"
// #include "../src/G4ParticleTypeConverter.cc"
// #include "../src/G4ResonanceDecay.cc"
// #include "../src/G4ResonanceNames.cc"
// #include "../src/G4Scatterer.cc"
// #include "../src/G4VAngularDistribution.cc"
// #include "../src/G4VAnnihilationCollision.cc"
// #include "../src/G4VCollision.cc"
// #include "../src/G4VCrossSectionSource.cc"
// #include "../src/G4VElasticCollision.cc"
// #include "../src/G4VScatterer.cc"
// #include "../src/G4VScatteringCollision.cc"
// #include "../src/G4VXResonance.cc"
// #include "../src/G4XAnnihilationChannel.cc"
// #include "../src/G4XAqmElastic.cc"
// #include "../src/G4XAqmInelastic.cc"
// #include "../src/G4XAqmInelasticRes.cc"
// #include "../src/G4XAqmInelasticString.cc"
// #include "../src/G4XAqmTotal.cc"
// #include "../src/G4XDeltaDeltaTable.cc"
// #include "../src/G4XDeltaDeltaToNDelta.cc"
// #include "../src/G4XDeltaDeltaToNN.cc"
// #include "../src/G4XDeltaDeltastarTable.cc"
// #include "../src/G4XDeltaNstarTable.cc"
// #include "../src/G4XMesonBaryonElastic.cc"
// #include "../src/G4XMesonBaryonOneS.cc"
// #include "../src/G4XMesonBaryonTwoS.cc"
// #include "../src/G4XMesonMesonElastic.cc"
// #include "../src/G4XNDeltaTable.cc"
// #include "../src/G4XNDeltaTo.cc"
// #include "../src/G4XNDeltaToDeltaDelta.cc"
// #include "../src/G4XNDeltaToNN.cc"
// #include "../src/G4XNDeltastarTable.cc"
// #include "../src/G4XNNElastic.cc"
// #include "../src/G4XNNElasticLowE.cc"
// #include "../src/G4XNNTotal.cc"
// #include "../src/G4XNNTotalLowE.cc"
// #include "../src/G4XNNstarTable.cc"
// #include "../src/G4XPDGElastic.cc"
// #include "../src/G4XPDGTotal.cc"
// #include "../src/G4XResonance.cc"
// #include "../src/G4XnpElastic.cc"
// #include "../src/G4XnpElasticLowE.cc"
// #include "../src/G4XnpTotal.cc"
// #include "../src/G4XnpTotalLowE.cc"
// #include "../src/G4XpimNTotal.cc"
// #include "../src/G4XpipNTotal.cc"
// //util
// #include "../src/G4ExcitedString.cc"
// #include "../src/G4Fancy3DNucleus.cc"
// #include "../src/G4FermiMomentum.cc"
// #include "../src/G4Fragment.cc"
// #include "../src/G4GeneralPhaseSpaceDecay.cc"
// #include "../src/G4KineticTrack.cc"
// #include "../src/G4KineticTrackVector.cc"
// #include "../src/G4NuclearFermiDensity.cc"
// #include "../src/G4NuclearShellModelDensity.cc"
// #include "../src/G4Nucleon.cc"
// #include "../src/G4Parton.cc"
// #include "../src/G4SampleResonance.cc"
// //management
// #include "../src/G4V3DNucleus.cc"
// #include "../src/G4VAnnihilator.cc"
// #include "../src/G4VCrossSectionBase.cc"
// #include "../src/G4VElasticScatterer.cc"
// #include "../src/G4VHighEnergyGenerator.cc"
// #include "../src/G4VInElasticScatterer.cc"
// #include "../src/G4VIntraNuclearTransportModel.cc"
// #include "../src/G4VKineticNucleon.cc"
// #include "../src/G4VNuclearDensity.cc"
// #include "../src/G4VPreCompoundModel.cc"
// //hadutil
// #include "../src/G4ReactionProduct.cc"
// // magnetic field
// #include "../src/G4CashKarpRKF45.cc"
// #include "../src/G4ChordFinder.cc"
// #include "../src/G4ClassicalRK4.cc"
// #include "../src/G4DELPHIMagField.cc"
// #include "../src/G4EqMagElectricField.cc"
// #include "../src/G4EquationOfMotion.cc"
// #include "../src/G4ExplicitEuler.cc"
// #include "../src/G4FieldManager.cc"
// #include "../src/G4FieldTrack.cc"
// #include "../src/G4HarmonicPolMagField.cc"
// #include "../src/G4HelixExplicitEuler.cc"
// #include "../src/G4HelixHeum.cc"
// #include "../src/G4HelixImplicitEuler.cc"
// #include "../src/G4HelixSimpleRunge.cc"
// #include "../src/G4ImplicitEuler.cc"
// #include "../src/G4LineCurrentMagField.cc"
// #include "../src/G4LineSection.cc"
// #include "../src/G4MagErrorStepper.cc"
// #include "../src/G4MagHelicalStepper.cc"
// #include "../src/G4MagIntegratorDriver.cc"
// #include "../src/G4MagIntegratorStepper.cc"
// #include "../src/G4Mag_EqRhs.cc"
// #include "../src/G4Mag_SpinEqRhs.cc"
// #include "../src/G4Mag_UsualEqRhs.cc"
// #include "../src/G4QuadrupoleMagField.cc"
// #include "../src/G4RKG3_Stepper.cc"
// #include "../src/G4SimpleHeum.cc"
// #include "../src/G4SimpleRunge.cc"
// #include "../src/G4UniformElectricField.cc"
// #include "../src/G4UniformMagField.cc"
// //pre-compound
// #include "../src/G4PreCompoundEmission.cc"
// #include "../src/G4PreCompoundFragmentVector.cc"
// #include "../src/G4PreCompoundModel.cc"
// #include "../src/G4PreCompoundParameters.cc"
// #include "../src/G4PreCompoundTransitions.cc"
// #include "../src/G4VPreCompoundFragment.cc"
// #include "../src/G4VPreCompoundIon.cc"
// #include "../src/G4VPreCompoundNucleon.cc"

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new Tst25DetectorConstruction);
  runManager->SetUserInitialization(new Tst25PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst25RunAction);
  runManager->SetUserAction(new Tst25PrimaryGeneratorAction);
  runManager->SetUserAction(new Tst25StackingAction);

  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;
  return 0;
}


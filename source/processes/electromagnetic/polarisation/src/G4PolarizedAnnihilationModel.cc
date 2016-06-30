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
// $Id: G4PolarizedAnnihilationModel.cc 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedAnnihilationModel
//
// Author:        Andreas Schaelicke
//
// Creation date: 01.05.2005
//
// Modifications:
// 18-07-06 use newly calculated cross sections (P. Starovoitov)
// 21-08-06 update interface (A. Schaelicke)
// 17-11-06 add protection agaist e+ zero energy PostStep (V.Ivanchenko)
// 10-07-07 copied Initialise() method from G4eeToTwoGammaModel to provide a  
//          local ParticleChangeForGamma object and reduce overhead 
//          in SampleSecondaries()  (A. Schaelicke)
//
//
// Class Description:
//
// Implementation of polarized gamma Annihilation scattering on free electron
// 

// -------------------------------------------------------------------
#include "G4PolarizedAnnihilationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4StokesVector.hh"
#include "G4PolarizedAnnihilationCrossSection.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4TrackStatus.hh"
#include "G4Gamma.hh"

G4PolarizedAnnihilationModel::G4PolarizedAnnihilationModel(const G4ParticleDefinition* p, 
							   const G4String& nam)
  : G4eeToTwoGammaModel(p,nam),
    crossSectionCalculator(nullptr),
    verboseLevel(0),
    gParticleChange(nullptr),
    gIsInitialised(false)
{
  crossSectionCalculator=new G4PolarizedAnnihilationCrossSection();
}

G4PolarizedAnnihilationModel::~G4PolarizedAnnihilationModel()
{
  if (crossSectionCalculator) delete crossSectionCalculator;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedAnnihilationModel::Initialise(const G4ParticleDefinition*,
                                     const G4DataVector&)
{
  //  G4eeToTwoGammaModel::Initialise(part,dv);
  if(gIsInitialised) return;
  gParticleChange = GetParticleChangeForGamma();
  gIsInitialised = true;
}

G4double G4PolarizedAnnihilationModel::ComputeCrossSectionPerElectron(
                                const G4ParticleDefinition* pd,
                                      G4double kinEnergy, 
                                      G4double cut,
                                      G4double emax)
{
  G4double xs = G4eeToTwoGammaModel::ComputeCrossSectionPerElectron(pd,kinEnergy,
								cut,emax);

  G4double polzz = theBeamPolarization.z()*theTargetPolarization.z();
  G4double poltt = theBeamPolarization.x()*theTargetPolarization.x() 
                 + theBeamPolarization.y()*theTargetPolarization.y();
  if (polzz!=0 || poltt!=0) {
    G4double xval,lasym,tasym;
    ComputeAsymmetriesPerElectron(kinEnergy,xval,lasym,tasym);
    xs*=(1.+polzz*lasym+poltt*tasym);
  }

  return xs;
}

void G4PolarizedAnnihilationModel::ComputeAsymmetriesPerElectron(G4double ene,
					       G4double & valueX,
					       G4double & valueA,
					       G4double & valueT)
{
  // *** calculate asymmetries
  G4double gam = 1. + ene/electron_mass_c2;
  G4double xs0=crossSectionCalculator->TotalXSection(0.,1.,gam,
					       G4StokesVector::ZERO,
					       G4StokesVector::ZERO);
  G4double xsA=crossSectionCalculator->TotalXSection(0.,1.,gam,
					       G4StokesVector::P3,
					       G4StokesVector::P3);
  G4double xsT1=crossSectionCalculator->TotalXSection(0.,1.,gam,
					       G4StokesVector::P1,
					       G4StokesVector::P1);
  G4double xsT2=crossSectionCalculator->TotalXSection(0.,1.,gam,
					       G4StokesVector::P2,
					       G4StokesVector::P2);
  G4double xsT=0.5*(xsT1+xsT2);
  
  valueX=xs0;
  valueA=xsA/xs0-1.;
  valueT=xsT/xs0-1.;
  //  G4cout<<valueX<<"\t"<<valueA<<"\t"<<valueT<<"   energy = "<<gam<<G4endl;
  if ( (valueA < -1) || (1 < valueA)) {
    G4cout<< " ERROR PolarizedAnnihilationPS::ComputeAsymmetries \n";
    G4cout<< " something wrong in total cross section calculation (valueA)\n";
    G4cout<<"*********** LONG "<<valueX<<"\t"<<valueA<<"\t"<<valueT<<"   energy = "<<gam<<G4endl;
  }
  if ( (valueT < -1) || (1 < valueT)) {
    G4cout<< " ERROR PolarizedAnnihilationPS::ComputeAsymmetries \n";
    G4cout<< " something wrong in total cross section calculation (valueT)\n";
    G4cout<<"****** TRAN "<<valueX<<"\t"<<valueA<<"\t"<<valueT<<"   energy = "<<gam<<G4endl;
  }
}


void G4PolarizedAnnihilationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						     const G4MaterialCutsCouple* /*couple*/,
						     const G4DynamicParticle* dp,
						     G4double /*tmin*/,
						     G4double /*maxEnergy*/) 
{
//   G4ParticleChangeForGamma*  gParticleChange 
//     = dynamic_cast<G4ParticleChangeForGamma*>(pParticleChange);
  const G4Track * aTrack = gParticleChange->GetCurrentTrack();

  // kill primary 
  gParticleChange->SetProposedKineticEnergy(0.);
  gParticleChange->ProposeTrackStatus(fStopAndKill);

  // V.Ivanchenko add protection against zero kin energy
  G4double PositKinEnergy = dp->GetKineticEnergy();

  if(PositKinEnergy < DBL_MIN) {

    G4double cosTeta = 2.*G4UniformRand()-1.;
    G4double sinTeta = std::sqrt((1.0 - cosTeta)*(1.0 + cosTeta));
    G4double phi     = twopi * G4UniformRand();
    G4ThreeVector dir(sinTeta*std::cos(phi), sinTeta*std::sin(phi), cosTeta);
    fvect->push_back( new G4DynamicParticle(G4Gamma::Gamma(), dir, electron_mass_c2));
    fvect->push_back( new G4DynamicParticle(G4Gamma::Gamma(),-dir, electron_mass_c2));
    return;
  }

  // *** obtain and save target and beam polarization ***
  G4PolarizationManager * polarizationManager = G4PolarizationManager::GetInstance();

  // obtain polarization of the beam
  theBeamPolarization = aTrack->GetPolarization();

  // obtain polarization of the media
  G4VPhysicalVolume*  aPVolume  = aTrack->GetVolume();
  G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();
  const G4bool targetIsPolarized = polarizationManager->IsPolarized(aLVolume);
  theTargetPolarization = polarizationManager->GetVolumePolarization(aLVolume);

  if (verboseLevel >= 1) {
    G4cout << "G4PolarizedComptonModel::SampleSecondaries in "
           <<  aLVolume->GetName() << G4endl;
  }

  // transfer target electron polarization in frame of positron
  if (targetIsPolarized)
      theTargetPolarization.rotateUz(dp->GetMomentumDirection());
  
  G4ParticleMomentum PositDirection = dp->GetMomentumDirection();

  // polar asymmetry:
  G4double polarization = theBeamPolarization.p3()*theTargetPolarization.p3();

  G4double gamam1 = PositKinEnergy/electron_mass_c2;
  G4double gama   = gamam1+1. , gamap1 = gamam1+2.;
  G4double sqgrate = std::sqrt(gamam1/gamap1)/2. , sqg2m1 = std::sqrt(gamam1*gamap1);

  // limits of the energy sampling
  G4double epsilmin = 0.5 - sqgrate , epsilmax = 0.5 + sqgrate;
  G4double epsilqot = epsilmax/epsilmin;
  
  //
  // sample the energy rate of the created gammas 
  // note: for polarized partices, the actual dicing strategy 
  //       will depend on the energy, and the degree of polarization !!
  //
  G4double epsil;
  G4double gmax=1. + std::fabs(polarization); // crude estimate

  //G4bool check_range=true;

  crossSectionCalculator->Initialize(epsilmin, gama, 0.,  theBeamPolarization, theTargetPolarization);
  if (crossSectionCalculator->DiceEpsilon()<0) {
    G4cout<<"ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
	  <<"epsilmin DiceRoutine not appropriate ! "<<crossSectionCalculator->DiceEpsilon()<<G4endl;
    //check_range=false;
  }

  crossSectionCalculator->Initialize(epsilmax, gama, 0.,  theBeamPolarization, theTargetPolarization);
  if (crossSectionCalculator->DiceEpsilon()<0) {
    G4cout<<"ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
	  <<"epsilmax DiceRoutine not appropriate ! "<<crossSectionCalculator->DiceEpsilon()<<G4endl;
    //check_range=false;
  }

  G4int ncount=0;
  G4double trejectmax=0.;
  G4double treject;


  do {
    // 
    epsil = epsilmin*std::pow(epsilqot,G4UniformRand());

    crossSectionCalculator->Initialize(epsil, gama, 0., theBeamPolarization, theTargetPolarization,1);

    treject = crossSectionCalculator->DiceEpsilon(); 
    treject*=epsil;

    if (treject>gmax  || treject<0.) 
      G4cout<<"ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
	    <<" eps ("<<epsil<<") rejection does not work properly: "<<treject<<G4endl;
    ++ncount;
    if (treject>trejectmax) trejectmax=treject;
    if (ncount>1000) {
      G4cout<<"WARNING  in PolarizedAnnihilationPS::PostStepDoIt\n"
	    <<"eps dicing very inefficient ="<<trejectmax/gmax
	    <<", "<<treject/gmax<<".  For secondary energy = "<<epsil<<"   "<<ncount<<G4endl;
      break;
    }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while( treject < gmax*G4UniformRand() );

  //
  // scattered Gamma angles. ( Z - axis along the parent positron)
  //
   
  G4double cost = (epsil*gamap1-1.)/(epsil*sqg2m1);
  G4double sint = std::sqrt((1.+cost)*(1.-cost));
  G4double phi  = 0.;
  G4double   beamTrans = std::sqrt(sqr(theBeamPolarization.p1()) + sqr(theBeamPolarization.p2()));
  G4double targetTrans = std::sqrt(sqr(theTargetPolarization.p1()) + sqr(theTargetPolarization.p2()));

  //  G4cout<<"phi dicing START"<<G4endl;
  do{
    phi  = twopi * G4UniformRand();
    crossSectionCalculator->Initialize(epsil, gama, 0., theBeamPolarization, theTargetPolarization,2);

    G4double gdiced =crossSectionCalculator->getVar(0);
    gdiced += crossSectionCalculator->getVar(3)*theBeamPolarization.p3()*theTargetPolarization.p3();
    gdiced += 1.*(std::fabs(crossSectionCalculator->getVar(1)) 
		  + std::fabs(crossSectionCalculator->getVar(2)))*beamTrans*targetTrans;
    gdiced += 1.*std::fabs(crossSectionCalculator->getVar(4))
      *(std::fabs(theBeamPolarization.p3())*targetTrans + std::fabs(theTargetPolarization.p3())*beamTrans);

    G4double gdist = crossSectionCalculator->getVar(0);
    gdist += crossSectionCalculator->getVar(3)*theBeamPolarization.p3()*theTargetPolarization.p3();
    gdist += crossSectionCalculator->getVar(1)*(std::cos(phi)*theBeamPolarization.p1() 
						+ std::sin(phi)*theBeamPolarization.p2())
                                              *(std::cos(phi)*theTargetPolarization.p1() 
						+ std::sin(phi)*theTargetPolarization.p2());
    gdist += crossSectionCalculator->getVar(2)*(std::cos(phi)*theBeamPolarization.p2() 
						- std::sin(phi)*theBeamPolarization.p1())
                                              *(std::cos(phi)*theTargetPolarization.p2() 
						- std::sin(phi)*theTargetPolarization.p1());
    gdist += crossSectionCalculator->getVar(4)
      *(std::cos(phi)*theBeamPolarization.p3()*theTargetPolarization.p1()
	+ std::cos(phi)*theBeamPolarization.p1()*theTargetPolarization.p3() 
	+ std::sin(phi)*theBeamPolarization.p3()*theTargetPolarization.p2() 
	+ std::sin(phi)*theBeamPolarization.p2()*theTargetPolarization.p3());

    treject = gdist/gdiced;
    //G4cout<<" treject = "<<treject<<" at phi = "<<phi<<G4endl;
     if (treject>1.+1.e-10 || treject<0){
       G4cout<<"!!!ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
	     <<" phi rejection does not work properly: "<<treject<<G4endl;
       G4cout<<" gdiced = "<<gdiced<<G4endl;
       G4cout<<" gdist = "<<gdist<<G4endl;
       G4cout<<" epsil = "<<epsil<<G4endl;
     }
     
     if (treject<1.e-3) {
       G4cout<<"!!!ERROR in PolarizedAnnihilationPS::PostStepDoIt\n"
	    <<" phi rejection does not work properly: "<<treject<<"\n";
       G4cout<<" gdiced="<<gdiced<<"   gdist="<<gdist<<"\n";
       G4cout<<" epsil = "<<epsil<<G4endl;
     }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while( treject < G4UniformRand() );
  //  G4cout<<"phi dicing END"<<G4endl;

  G4double dirx = sint*std::cos(phi) , diry = sint*std::sin(phi) , dirz = cost;

  //
  // kinematic of the created pair
  //
  G4double TotalAvailableEnergy = PositKinEnergy + 2*electron_mass_c2;
  G4double Phot1Energy = epsil*TotalAvailableEnergy;
  G4double Phot2Energy =(1.-epsil)*TotalAvailableEnergy;

  // *** prepare calculation of polarization transfer ***
  G4ThreeVector Phot1Direction (dirx, diry, dirz);

  // get interaction frame
  G4ThreeVector  nInteractionFrame = 
    G4PolarizationHelper::GetFrame(PositDirection,Phot1Direction);
     
  // define proper in-plane and out-of-plane component of initial spins
  theBeamPolarization.InvRotateAz(nInteractionFrame,PositDirection);
  theTargetPolarization.InvRotateAz(nInteractionFrame,PositDirection);

  // calculate spin transfere matrix

  crossSectionCalculator->Initialize(epsil,gama,phi,theBeamPolarization,theTargetPolarization,2);

  // **********************************************************************

  Phot1Direction.rotateUz(PositDirection);   
  // create G4DynamicParticle object for the particle1  
  G4DynamicParticle* aParticle1= new G4DynamicParticle (G4Gamma::Gamma(),
							Phot1Direction, Phot1Energy);
  finalGamma1Polarization=crossSectionCalculator->GetPol2();
  G4double n1=finalGamma1Polarization.mag2();
  if (n1>1) {
    G4cout<<"ERROR: PolarizedAnnihilation Polarization Vector at epsil = "
	  <<epsil<<" is too large!!! \n"
	  <<"annihi pol1= "<<finalGamma1Polarization<<", ("<<n1<<")\n";
    finalGamma1Polarization*=1./std::sqrt(n1);
  }

  // define polarization of first final state photon
  finalGamma1Polarization.SetPhoton();
  finalGamma1Polarization.RotateAz(nInteractionFrame,Phot1Direction);
  aParticle1->SetPolarization(finalGamma1Polarization.p1(),
			      finalGamma1Polarization.p2(),
			      finalGamma1Polarization.p3());

  fvect->push_back(aParticle1);


  // **********************************************************************

  G4double Eratio= Phot1Energy/Phot2Energy;
  G4double PositP= std::sqrt(PositKinEnergy*(PositKinEnergy+2.*electron_mass_c2));
  G4ThreeVector Phot2Direction (-dirx*Eratio, -diry*Eratio,
				(PositP-dirz*Phot1Energy)/Phot2Energy); 
  Phot2Direction.rotateUz(PositDirection); 
  // create G4DynamicParticle object for the particle2 
  G4DynamicParticle* aParticle2= new G4DynamicParticle (G4Gamma::Gamma(),
							Phot2Direction, Phot2Energy);

  // define polarization of second final state photon
  finalGamma2Polarization=crossSectionCalculator->GetPol3();
  G4double n2=finalGamma2Polarization.mag2();
  if (n2>1) {
    G4cout<<"ERROR: PolarizedAnnihilation Polarization Vector at epsil = "<<epsil<<" is too large!!! \n";
    G4cout<<"annihi pol2= "<<finalGamma2Polarization<<", ("<<n2<<")\n";
    
    finalGamma2Polarization*=1./std::sqrt(n2);
  }
  finalGamma2Polarization.SetPhoton();
  finalGamma2Polarization.RotateAz(nInteractionFrame,Phot2Direction);
  aParticle2->SetPolarization(finalGamma2Polarization.p1(),
			      finalGamma2Polarization.p2(),
			      finalGamma2Polarization.p3());

  fvect->push_back(aParticle2);
}

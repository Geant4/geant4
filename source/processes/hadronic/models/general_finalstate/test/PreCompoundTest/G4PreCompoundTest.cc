#include "G4ios.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>

#include "globals.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"


#include "PCTProjectile.hh"
#include "PCTTarget.hh"
#include "PCTCompositeNucleus.hh"
#include "PCTCrossSection.hh"
#include "PCTOptions.hh"
#include "PCTWriter.hh"
#include "PCTBinaryCascade.hh"

struct DeleteFragment
{
    template<typename T>
    void operator()(const T* ptr) const
	{
	    delete ptr;
	}
};


int main()
{

//     G4ParticleDefinition *theGamma = G4Gamma::GammaDefinition();
//     G4ParticleDefinition *theElectron = G4Electron::ElectronDefinition();
//     G4ParticleDefinition *theNeutron = G4Neutron::NeutronDefinition();
//     G4ParticleDefinition *theProton = G4Proton::ProtonDefinition();   
//     G4ParticleDefinition *theDeuteron = G4Deuteron::DeuteronDefinition();
//     G4ParticleDefinition *theTriton = G4Triton::TritonDefinition();
//     G4ParticleDefinition *theHelium3 = G4He3::He3Definition();
//     G4ParticleDefinition *theAlpha = G4Alpha::AlphaDefinition();
//     G4ParticleDefinition *theIon = G4GenericIon::GenericIonDefinition();
//     theProton->SetCuts(1.0);
//     theGamma->SetCuts(1.0);
//     theElectron->SetCuts(1.0);
//     theNeutron->SetCuts(1.0);
//     theDeuteron->SetCuts(1.0);
//     theTriton->SetCuts(1.0);
//     theHelium3->SetCuts(1.0);
//     theAlpha->SetCuts(1.0);
//     theIon->SetCuts(1.0);
    
    // --------
    // Greeting
    // --------
    G4cout << "\nPreCompound Model Test program\n";
    
    // -------------------------
    // Ask information for setup
    // -------------------------
    PCTOptions Opt;
    Opt.Initialize();
    
    // --------------
    // The Projectile
    // --------------
    PCTProjectile theProjectile("proton");
    theProjectile.SetDirection(Opt.GetProjectileDirection());
    theProjectile.SetKineticEnergy(Opt.GetProjectileKineticEnergy());

    // ------------------
    // The Target Nucleus
    // ------------------
    PCTTarget theTarget(Opt);

    // -----------------------------
    // Prepare the Composite Nucleus
    // -----------------------------
    PCTCompositeNucleus CompositeNucleus(&theProjectile,&theTarget,
					 Opt.GetNumberOfHoles(),
					 Opt.GetNumberOfParticles(),
					 Opt.GetNumberOfCharged());

    // ---------------------------
    // Calculate the Cross Section
    // ---------------------------
    PCTCrossSection * theCrossSection = new PCTCrossSection(&theProjectile,&theTarget);
    G4double theXS = theCrossSection->GetCrossSection();


    // -----------------
    // Prepare the Model
    // -----------------
    G4ExcitationHandler * theExcitationHandlerPtr = new G4ExcitationHandler();
    G4Evaporation * theEvaporationModelPtr = new G4Evaporation();
    theExcitationHandlerPtr->SetEvaporation(theEvaporationModelPtr);
  
    G4PreCompoundModel * thePreCompoundModelPtr =
      new G4PreCompoundModel(theExcitationHandlerPtr);
    
    if (Opt.GetPreeqEmissionMode() == HETC_emission_mode) thePreCompoundModelPtr->UseHETCEmission();
    if (Opt.GetPreeqTransitionMode() == GNASH_transition_mode) thePreCompoundModelPtr->UseGNASHTransition();

    if (Opt.GetEvapMode() == GEM_mode)  theEvaporationModelPtr->SetGEMChannel();
    if (Opt.UsingFermi()) theExcitationHandlerPtr->SetMaxAandZForFermiBreakUp(13,7);
    else theExcitationHandlerPtr->SetMaxAandZForFermiBreakUp(1,1);


    PCTBinaryCascade * theINCptr(0);
    if (Opt.UsingINC()) 
      {
	theINCptr = new PCTBinaryCascade();
	theINCptr->SetDeExcitation(thePreCompoundModelPtr);
      }

    // -----------------
    // Prepare the file
    // -----------------
    PCTWriter theFile;
    theFile.OpenFile(Opt.GetFilename());
    theFile.WriteHeader(CompositeNucleus,theXS);
    
    std::cout << "Iteration: " << std::setw(10) << 0;
    std::cout.flush();  
    for (G4int i = 0; i < Opt.GetNumberOfIterations(); i++) 
      {
	std::cout << "\b\b\b\b\b\b\b\b\b\b" << std::setw(10) << i+1;	
	std::cout.flush();
// 	std::cout << "########################\n";
// 	std::cout << "# Event n.: "<< std::setw(10) << i+1 << " #\n";
// 	std::cout << "########################\n";
	G4ReactionProductVector * theReactionProductVectorPtr(0);
	// DeExcite the nucleus 
	if (Opt.UsingINC())
	  {
	    const G4Fragment * target = CompositeNucleus.GetNewCNucleus();
	    G4int tA = G4int(target->GetA()) - CompositeNucleus.GetProjectile()->GetA();
	    G4int tZ = G4int(target->GetZ()) - CompositeNucleus.GetProjectile()->GetZ();
	    theReactionProductVectorPtr = theINCptr->DeExcite(CompositeNucleus.GetProjectile(),tA,tZ);
	    if (!theReactionProductVectorPtr)
	      {
		theReactionProductVectorPtr = thePreCompoundModelPtr->DeExcite(*target);
	      }
	  }
	else 
	  {
	    theReactionProductVectorPtr = thePreCompoundModelPtr->DeExcite(*(CompositeNucleus.GetNewCNucleus()));
	  }
	theFile.WriteReaction(i+1,theReactionProductVectorPtr,CompositeNucleus.GetLastCNucleus());
	
	std::for_each(theReactionProductVectorPtr->begin(),theReactionProductVectorPtr->end(),
		      DeleteFragment());
	delete theReactionProductVectorPtr;
	
      }	      
    std::cout << std::endl;
    theFile.CloseFile();      

    //    delete thePreCompoundModelPtr;
    delete theExcitationHandlerPtr;
    return 0;
}

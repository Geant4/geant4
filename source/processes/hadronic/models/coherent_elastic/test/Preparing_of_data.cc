//  The main program for "G4ElasticHadrNucleusHE" class

#include <vector>
#include <utility>

#include "G4DynamicParticle.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Nucleus.hh"
#include "G4IonConstructor.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "fstream"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4Lambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiLambda.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4AntiOmegaMinus.hh"

 int main()

 {
	G4IonConstructor Ion;
	Ion.ConstructParticle();

        // G4double CrossSection, ElasticCrossSec, Q2, dE, Tkin,
        //          InelasticCrossSec, Energy, Momentum;
        G4double Tkin, Momentum;

        G4double        px = 0;
        G4double        py = 0;
        G4double        pz = 1000;
        G4ThreeVector   inVector(px, py, pz);
        G4ThreeVector   outVector, aPosition(0., 0., 0.);
	std::vector<G4ParticleDefinition *> setOfParticles;
	setOfParticles.push_back(G4Proton::ProtonDefinition()); // nucleon
	setOfParticles.push_back(G4PionPlus::PionPlusDefinition());
	setOfParticles.push_back(G4PionMinus::PionMinusDefinition());
	
	setOfParticles.push_back(G4KaonPlus::KaonPlusDefinition());
	setOfParticles.push_back(G4KaonMinus::KaonMinusDefinition());
	// setOfParticles.push_back(G4KaonZero::KaonZeroDefinition());
	// setOfParticles.push_back(G4AntiKaonZero::AntiKaonZeroDefinition());
	
	setOfParticles.push_back(G4SigmaPlus::SigmaPlusDefinition());
	setOfParticles.push_back(G4SigmaMinus::SigmaMinusDefinition());
	setOfParticles.push_back(G4AntiSigmaPlus::AntiSigmaPlusDefinition());
	setOfParticles.push_back(G4AntiSigmaMinus::AntiSigmaMinusDefinition());

	setOfParticles.push_back(G4Lambda::LambdaDefinition());
	setOfParticles.push_back(G4AntiLambda::AntiLambdaDefinition());

	setOfParticles.push_back(G4XiMinus::XiMinusDefinition());
	setOfParticles.push_back(G4XiZero::XiZeroDefinition());
	setOfParticles.push_back(G4AntiXiMinus::AntiXiMinusDefinition());
	setOfParticles.push_back(G4AntiXiZero::AntiXiZeroDefinition());
	
	setOfParticles.push_back(G4OmegaMinus::OmegaMinusDefinition());
	setOfParticles.push_back(G4AntiOmegaMinus::AntiOmegaMinusDefinition());
	
	std::vector<double> setOfNuclei;
	for(int i=4; i<208; i++) setOfNuclei.push_back(i);
	setOfNuclei.push_back(209);
	setOfNuclei.push_back(210);
	setOfNuclei.push_back(211);
	setOfNuclei.push_back(212);
	setOfNuclei.push_back(233);
	setOfNuclei.push_back(234);
	setOfNuclei.push_back(235);
	setOfNuclei.push_back(236);
	setOfNuclei.push_back(237);
	setOfNuclei.push_back(238);


        G4DynamicParticle *  aParticle = new G4DynamicParticle;
        G4VParticleChange *  aChange   = new G4VParticleChange;
        G4Nucleus            aNucleus;

// const    G4double N = 90, Z = 45;

        // G4double aTime = 0.1;
        for(size_t ip=0; ip<setOfParticles.size(); ip++)
	{
	  aParticle->SetDefinition(setOfParticles[ip]);
  //  ----------  The preparing files for n nuclei -----------
          for(size_t ik=0; ik<setOfNuclei.size(); ik++)
          {
            const G4double   N = setOfNuclei[ik];
            const G4double   Z = N/2.;
            aNucleus.SetParameters( N, Z);
            //   The preparation of file for total Energy from 1.5 GeV to 1000 GeV
            G4ElasticHadrNucleusHE    aElasticRandom(setOfParticles[ip], &aNucleus);
	    G4cout<<G4endl<< " The array is created !!! "<<G4endl;
	  }   //  Nucleus
        }

//  -----------------  The end of preparing -----------------

/*

//  =========================================================

//  -----------  The sampling of 500000 events  -------------  
//       -----------  for last nucleus  ----------------
//       -----------     at 100 GeV     ----------------
       G4Track *    secTrack;

       for(G4int i1=1; i1<=2; i1+=5)
       {
         G4ElasticHadrNucleusHE  aElasticRandom(aParticle, pNucl,
                                  2000., 10000., 5);
                 Momentum = i1*100000;
                 Tkin     = sqrt(Momentum*Momentum+938.27*938.27)-938.27;

                 inVector.setZ(Momentum);
                 aParticle->SetMomentum(inVector);
        G4Track  aTrack   = G4Track( aParticle, aTime, aPosition);

                 inVector =  aParticle->GetMomentum();

            std::ofstream TestFile("q2from4.dat", std::ios::out);
            TestFile.precision(9);
            TestFile.setf(std::ios::scientific);

            for(G4int i2=1; i2<500001; i2++)
            {
              Q2 = aElasticRandom.RandomElastic1(aParticle,  pNucl);
              TestFile<<"  "<<Q2<<G4endl;
             }  //  i2

      G4cout<<" end cicle of Q2 "<<G4endl;
           }      //  i1
*/

      return 0;

  }           


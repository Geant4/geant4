//  The main program for "G4ElasticHadrNucleusHE" class
//  Using data preparing in the directory 'Elastic'

#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4SigmaPlus.hh"
#include "G4Lambda.hh"
#include "G4SigmaMinus.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Nucleus.hh"
#include "G4IonConstructor.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "fstream"

 int main()

 {
	G4IonConstructor Ion;
	Ion.ConstructParticle();

        G4double CrossSection, ElasticCrossSec, Q2, dE, Tkin,
                 InelasticCrossSec, Energy, Momentum;

        G4double        px = 0;
        G4double        py = 0;
        G4double        pz = 1000;
        G4ThreeVector   inVector(px, py, pz);
        G4ThreeVector   outVector, aPosition(0., 0., 0.);
        G4double        setOfNuclei[8]=
               {12.0, 16.0, 28.0, 40.0, 58.0, 64.0, 90.0, 208.0}; 

        G4PionPlus        *  aPionP    = G4PionPlus::
                                            PionPlusDefinition();
        G4PionMinus       *  aPionM    = G4PionMinus::
                                            PionMinusDefinition();
        G4KaonPlus        *  aKaonP    = G4KaonPlus::
                                            KaonPlusDefinition();
        G4KaonMinus       *  aKaonM    = G4KaonMinus::
                                            KaonMinusDefinition();
        G4Proton          *  aProton   = G4Proton::Proton();
        G4AntiProton      *  aProtonA  = G4AntiProton::AntiProton();
        G4Lambda          *  aLambda   = G4Lambda::Lambda();
//        G4SigmaPlus       *  aSigmaP   = G4SigmaPlus::
//                                            SigmaPLusDefinition();
        G4SigmaMinus      *  aSigmaM   = G4SigmaMinus::
                                            SigmaMinusDefinition();

        G4DynamicParticle *  aParticle = new G4DynamicParticle;
        G4VParticleChange *  aChange   = new G4VParticleChange;
        G4Nucleus            aNucleus, * pNucl;

// const    G4double N = 90, Z = 45;

        G4double aTime = 0.1;
        aParticle->SetDefinition(aSigmaM);

//  ------  The using of data files for some  nucleus ------
//  -----------  The sampling of 50000 events  -------------  
//       -----------  for last nucleus  ----------------
//       -----------     at 100 GeV     ----------------

        for(G4int ik=5; ik<6; ik++)
        {
const G4double   N = setOfNuclei[ik];
const G4double   Z = N/2.0;

        aNucleus.SetParameters( N, Z);
        pNucl   = &aNucleus;


        G4double aTime = 0.1;

        Momentum = 100000;
        Tkin     = sqrt(Momentum*Momentum+938.27*938.27)-938.27;

        inVector.setZ(Momentum);
        aParticle->SetMomentum(inVector);

        G4ElasticHadrNucleusHE  aElasticRandom(aParticle, pNucl, 5);

    G4Track  aTrack   = G4Track( aParticle, aTime, aPosition);

            inVector =  aParticle->GetMomentum();

            std::ofstream TestFile("q2for64.dat", std::ios::out);
            TestFile.precision(9);
            TestFile.setf(std::ios::scientific);

            for(G4int i2=1; i2<50001; i2++)
            {
              Q2 = aElasticRandom.RandomElastic1(aParticle,  pNucl);
              TestFile<<"  "<<Q2<<G4endl;
            }  //  i2

        }   //  Nucleus

      return 0;
  }           


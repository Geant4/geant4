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
// $Id: G4QFragmentation.cc,v 1.16 2009-07-24 16:37:03 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//                 History: 
//     Created by Mikhail Kossov, October 2006
//     CHIPS QGS fragmentation class 
//     For comparison similar member functions can be found in the G4 classes:
//     G4QGSParticipants
//     G4QGSModels
//     G4ExcitedStringDecay
// -----------------------------------------------------------------------------
// Short description: CHIPS QG string fragmentation class
// Rhe key member function is Scatter, making the interaction (see G4QCollision)
// -----------------------------------------------------------------------------

//#define debug
//#define edebug
//#define pdebug
//#define sdebug
//#define ppdebug

#include "G4QFragmentation.hh"

// Promoting model parameters from local variables class properties @@(? M.K.)

// Definition of static parameters
G4int    G4QFragmentation::nCutMax=7; 
G4double G4QFragmentation::theNucleonRadius=0.;           // M.K.? Was 1.5*fermi
// Parameters of diffractional fragmentation
G4double G4QFragmentation::widthOfPtSquare=-0.72*GeV*GeV; // pt -width2 forStringExcitation

G4QFragmentation::G4QFragmentation(){}

G4QFragmentation::~G4QFragmentation() {}

G4QHadronVector* G4QFragmentation::Breeder(const G4QNucleus &aNucleus,
                                           const G4QHadron  &aPrimary)
{ // This is the main member function, which returns the resulting vector of hadrons
  static const G4double  eps = 0.001;                              // Tolerance in MeV
  static const G4double  mProt = G4Proton::Proton()->GetPDGMass(); // Mass of proton
  static const G4double  mPi0  = G4PionZero::PionZero()->GetPDGMass(); // Mass of Pi0
  G4QStringVector* strings=0;
  G4QHadron aProjectile = aPrimary;       // As a primary is const
  G4LorentzRotation toZ;                  // Lorentz Transformation to the projectileSystem
  G4LorentzVector proj4M=aProjectile.Get4Momentum(); // Projectile 4-momentum in LS
#ifdef debug
  G4double tgMass=aNucleus.GetGSMass();   // Ground state mass of the nucleus
  G4cout<<"G4QFragmentation::Breeder:Called,p4M="<<proj4M<<",A="<<aNucleus<<tgMass<<G4endl;
#endif
  G4int tZ=aNucleus.GetZ();
  G4int tN=aNucleus.GetN();
  G4int tPDG=90000000+tZ*1000+tN;
  toZ.rotateZ(-proj4M.phi());
  toZ.rotateY(-proj4M.theta());
  G4LorentzVector zProj4M=toZ*proj4M;     // Proj 4-momentum in LS rotated to Z axis
  aProjectile.Set4Momentum(zProj4M);      // Now the projectile moves along Z axis
#ifdef edebug
  G4int totChg=aProjectile.GetCharge()+tZ;// Charge of the projectile+target for the CHECK
  G4int totBaN=aProjectile.GetBaryonNumber()+tZ+tN;// Baryon Number of Proj+Targ for CHECK
  G4LorentzVector tgLS4M(0.,0.,0.,tgMass);// Target 4-momentum in LS
  G4LorentzVector totLS4M=proj4M+tgLS4M;  // Total 4-momentum in LS
  G4LorentzVector totZLS4M=zProj4M+tgLS4M;// Total 4-momentum in LS with momentum along Z
  G4cout<<"-EMC-G4QFragmentation::Breeder: tLS4M="<<totLS4M<<", tZLS4M="<<totZLS4M<<G4endl;
  // === From nere all consideration is made in the rotated LS frame (proj is along Z) ===
#endif
  G4LorentzRotation toLab(toZ.inverse()); // Lorentz Transfornation "ZLS"->LS (at the end)
  G4int pPDG=aProjectile.GetPDGCode();    // The PDG code of the projectile
  G4double projM=aProjectile.GetMass();   // Mass of the projectile
  G4QInteractionVector theInteractions;   // A vector of interactions (taken from the Body)
  G4QPartonPairVector  thePartonPairs;    // The parton pairs (taken from the Body)
  G4int                ModelMode=SOFT;    // The current model type (taken from the Body)
  G4QNucleus* theNucleus = new G4QNucleus(tZ,tN); // Create a copyInTheHip toMoveFromLStoCM
  G4int attempts = 0;
  G4int maxAttempts=20;                   // @@ A parameter !!
  while(!strings)
  {
    if (attempts++ > maxAttempts ) 
    {
      G4cerr<<"***G4QFragmentation::Breeder: "<<attempts<<" to create a string ( > max="
            <<maxAttempts<<") --> try to increase maxAttempts"<<G4endl;
      G4Exception("G4QFragmentation::Breeder:","72",FatalException,"StringCreation");
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Breeder: The attempt # "<<attempts<<G4endl;
#endif
    theNucleus->InitByPDG(tPDG);          // Reinit the Nucleus for the new Attempt
#ifdef pdebug
    G4cout<<"G4QFragmentation::Breeder: Nucl4Mom="<<theNucleus->Get4Momentum()<<G4endl;
#endif
    theNucleus->Init3D();                 // 3D-initialisation(nucleons) of theNucleusClone
#ifdef edebug
    G4LorentzVector sum1=theNucleus->GetNucleons4Momentum();// Sum ofNucleons4M inRotatedLS
    G4cout<<"-EMC-G4QFragmentation::Breeder: attempt#"<<attempts<<", Nuc4M="<<sum1<<G4endl;
#endif
    G4ThreeVector theCurrentVelocity(0.,0.,0.);         // "CM" velosity (temporary)
    // @@ "target Nucleon" == "Proton at rest" case (M.K. ?)
    G4double nCons = 1;                                 // 1 or baryonNum of the Projectile
    G4int projAbsB=std::abs(aProjectile.GetBaryonNumber());// Fragment/Baryon (Meson/AntiB)
    if(projAbsB>1) nCons = projAbsB;                    // @@ what if this is a meson ?
    G4LorentzVector proj4M = aProjectile.Get4Momentum();// 4-mom of the projectile hadron
    G4double pz_per_projectile = proj4M.pz()/nCons;     // Momentum per nucleon (hadron?)
    // @@ use M_A/A instead of mProt ------------ M.K.
    G4double e_per_projectile = proj4M.e()/nCons+mProt; // @@ Add MassOfTargetProtonAtRest
    G4double vz = pz_per_projectile/e_per_projectile;   // Speed of the "oneNuclenCM"
#ifdef pdebug
    G4cout<<"G4QFragmentation::Breeder: Projectile4M="<<proj4M<<", Vz="<<vz<<", nC="
          <<nCons<<", pE="<<e_per_projectile<<G4endl;
#endif
    theCurrentVelocity.setZ(vz);                        // CM (w/r to one nucleon) velosity
    if(theNucleus) theNucleus->DoLorentzBoost(-theCurrentVelocity);// BoostTgNucleus to"CM"
#ifdef edebug
    G4LorentzVector sum2=theNucleus->GetNucleons4Momentum();// Sum ofNucleons4M inRotatedCM
    G4cout<<"-EMC-G4QFragmentation::Breeder: AfterBoost, v="<<vz<<", Nuc4M="<<sum2<<G4endl;
#endif
    G4LorentzVector cmProjMom = proj4M;                 // Copy the original proj4M in LS
    cmProjMom.boost(-theCurrentVelocity);               // Bring the LS proj4Mom to "CM"
    G4double eKin = cmProjMom.e()-projM;                // Primary kinetic energy (MeV!)
    // @@ The maxCuts can improve the performance at low energies
    G4int maxCuts = std::min( 7, std::max( 1, static_cast<int>(std::log(eKin/GeV)) ) );
#ifdef pdebug
    G4cout<<"G4QFragmentation::Breeder: Proj4MInCM="<<cmProjMom<<", pPDG="<<pPDG<<G4endl;
#endif
    //
    // >>>>>>>>>> Find collisions meeting collision conditions
    //
    G4QHadron* cmProjectile = new G4QHadron(pPDG,cmProjMom); // HipCopy of the CMProjectile
    G4QPomeron theProbability(pPDG);                    // the PDG must be a data member
    G4double outerRadius = theNucleus->GetOuterRadius();// Get the nucleus frontiers
#ifdef pdebug
    G4cout<<"G4QFragmentation::Breeder: OuterRad="<<outerRadius<<",mCut="<<maxCuts<<G4endl;
#endif
    // Check the reaction threshold 
    theNucleus->StartLoop();                            // Prepare Loop ovdefr nucleons
    G4QHadron* pNucleon = theNucleus->GetNextNucleon(); // Get the next nucleon to try
    G4double s = (cmProjMom + pNucleon->Get4Momentum()).mag2(); // Squared CM Energy
    G4double ThresholdMass = projM + pNucleon->GetMass(); // @@ Nucleon can be virtual...?
    ModelMode = SOFT;                                   // NOT-Diffractive hadronization
    if (s < sqr(ThresholdMass))                         // At ThP=0 is impossible(virtNucl)
    {
      G4cerr<<"***G4QFragmentation::Breeder: ThrM="<<ThresholdMass<<" > sqrt(s)="
            <<std::sqrt(s)<<G4endl;
      G4Exception("G4QFragmentation::Breeder:","72",FatalException,"LowEnergy");
    }
    if (s < sqr(ThresholdMass))                        // --> Only diffractive interaction
    {
#ifdef debug
      G4cout<<"G4QFragmentation::Breeder:*OnlyDiffraction*ThM="<<ThresholdMass<<">sqrt(s)="
            <<std::sqrt(s)<<" -> only Diffraction is possible"<<G4endl; // @@ Dif toQuasmon
#endif
      ModelMode = DIFFRACTIVE;
    }
    // Clean up all previous interactions and reset the counters
    std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
    theInteractions.clear();
    G4int totalCuts = 0;
    G4double impactUsed = 0.;

    while(!theInteractions.size()) // Till interactions are created inside the while LOOP
    {
#ifdef debug
      G4cout<<"G4QFragmentation::Breeder: *** Enter the interaction LOOP ***"<<G4endl;
#endif
      // choose random impact parameter
      std::pair<G4double, G4double> theImpactParameter;
      theImpactParameter = theNucleus->ChooseImpactXandY(outerRadius+theNucleonRadius);
      G4double impactX = theImpactParameter.first; 
      G4double impactY = theImpactParameter.second;
    
      // LOOP over nuclei of the target nucleus to select collisions
      theNucleus->StartLoop();                          // To get the same nucleon
      G4int nucleonCount = 0;                           // debug
      //G4QFragmentation_NPart = 0;                       // ? M.K.
      // @@ Get Next nucleon should update energies of the nucleons to conserv Energy !!
      while( (pNucleon = theNucleus->GetNextNucleon()) && totalCuts < maxCuts)
      {
        nucleonCount++; // for debug only (?)
#ifdef sdebug
        G4cout<<"G4QFragmentation::Breeder:LOOP overNucleons, totCuts="<<totalCuts<<G4endl;
#endif
        // Needs to be moved to Probability class @@@
        G4double s = (cmProjMom + pNucleon->Get4Momentum()).mag2();
        G4double Distance2 = sqr(impactX - pNucleon->GetPosition().x()) +
                             sqr(impactY - pNucleon->GetPosition().y());
        G4double Probability = theProbability.GetInelasticProbability(s, Distance2);// INEL
        // test for inelastic collision
        G4double rndNumber = G4UniformRand();           // For the printing purpose
        // ModelMode = DIFFRACTIVE;
#ifdef sdebug
        G4cout<<"G4QFragmentation::Breeder: NLOOP prob="<<Probability<<", rndm="<<rndNumber
              <<", d="<<std::sqrt(Distance2)<<G4endl;
#endif
        if (Probability > rndNumber) // Inelastic (diffractive or soft) interaction (JOB)
        {
          G4QHadron* aTarget = new G4QHadron(*pNucleon);// Copy selected nucleon for String
#ifdef edebug
          G4cout<<">EMC>G4QFragmentation::Breeder: Target Nucleon is filled, 4M/PDG="
                <<aTarget->Get4Momentum()<<aTarget->GetPDGCode()<<G4endl;
#endif
          // Now the energy of the nucleons must be updated in CMS
          theNucleus->DoLorentzBoost(theCurrentVelocity);// Boost theResNucleus toRotatedLS
          theNucleus->SubtractNucleon(pNucleon);         // Pointer to the used nucleon
          theNucleus->DoLorentzBoost(-theCurrentVelocity);// Boost theResNucleus back to CM
          if((theProbability.GetDiffractiveProbability(s,Distance2)/Probability >
              G4UniformRand() && ModelMode==SOFT ) || ModelMode==DIFFRACTIVE)
          { 
            // ------------->>>> diffractive interaction @@ IsSingleDiffractive called once
            if(IsSingleDiffractive()) ExciteSingDiffParticipants(cmProjectile, aTarget);
            else                          ExciteDiffParticipants(cmProjectile, aTarget);
            G4QInteraction* anInteraction = new G4QInteraction(cmProjectile);
            anInteraction->SetTarget(aTarget); 
            theInteractions.push_back(anInteraction);   //--> now theInteractions not empty
            anInteraction->SetNumberOfDiffractiveCollisions(1); // Why not increment? M.K.?
            // @@ Why not breake the NLOOP, if only one diffractive can happend?
            totalCuts++;                               // UpdateOfNucleons in't necessary
#ifdef debug
            G4cout<<"G4QFragmentation::Breeder:NLOOP DiffInteract, tC="<<totalCuts<<G4endl;
#endif
          }
          else
          {
            // -------------->>>>> nondiffractive = soft interaction
            // sample nCut+1 (cut Pomerons) pairs of strings can be produced
            G4int nCut;                                // Result in a chosen number of cuts
            G4double* running = new G4double[nCutMax]; // @@ This limits the max cuts
            for(nCut = 0; nCut < nCutMax; nCut++)      // Calculates multiCut probabilities
            {
              running[nCut]= theProbability.GetCutPomeronProbability(s, Distance2, nCut+1);
              if(nCut) running[nCut] += running[nCut-1];// Sum up with the previous one
            }
            G4double random = running[nCutMax-1]*G4UniformRand();
            for(nCut = 0; nCut < nCutMax; nCut++) if(running[nCut] > random) break;
            delete [] running;
#ifdef debug
            G4cout<<"G4QFragmentation::Breeder:NLOOP-Soft Chosen nCut="<<nCut<<G4endl;
#endif
            // @@ If nCut>0 interaction with a few nucleons is possible
            // @@ nCut is found with big efforts and now nCut=0 ?? M.K. ?? !!
            //nCut = 0; // @@ in original code ?? @@
            aTarget->IncrementCollisionCount(nCut+1); // @@ What about multyNucleon target?
            cmProjectile->IncrementCollisionCount(nCut+1);
            G4QInteraction* anInteraction = new G4QInteraction(cmProjectile);
            anInteraction->SetTarget(aTarget);
            anInteraction->SetNumberOfSoftCollisions(nCut+1);
            theInteractions.push_back(anInteraction);
            totalCuts += nCut+1;
#ifdef debug
            G4cout<<"G4QFragmentation::Breeder:NLOOP SoftInteract, tC="<<totalCuts<<G4endl;
#endif
            impactUsed=Distance2;
          }
        }
      }
      // When nucleon count is incremented, the LOOP stops, so nucleonCount==1 always!
#ifdef debug
      G4cout<<"G4QFragmentation::Breeder: NUCLEONCOUNT="<<nucleonCount<<G4endl;
#endif
    }
    //#ifdef debug
    G4cout<<"G4QFragmentation::Breeder: CUTDEBUG="<<totalCuts<<", ImpactParam="
          <<impactUsed<<", #ofInter="<<theInteractions.size()<<G4endl;
    //#endif
    //
    // ------------------ now build the parton pairs for the strings ------------------
    //
    for(unsigned i=0; i<theInteractions.size(); i++)
    {
      theInteractions[i]->SplitHadrons();
#ifdef edebug
      G4QHadron* projH=theInteractions[i]->GetProjectile(); // Projectile of theInteraction
      G4QHadron* targH=theInteractions[i]->GetTarget();     // Target of the Interaction
      G4LorentzVector pSP(0.,0.,0.,0.);                // Sum of parton's 4mom's for proj
      G4LorentzVector tSP(0.,0.,0.,0.);                // Sum of parton's 4mom's for proj
      std::list<G4QParton*> projCP=projH->GetColor();  // Pointers to proj Color-partons
      std::list<G4QParton*> projAC=projH->GetAntiColor();// PointersTo projAntiColorPartons
      std::list<G4QParton*> targCP=targH->GetColor();  // Pointers to targ Color-partons
      std::list<G4QParton*> targAC=targH->GetAntiColor();// PointersTo targAntiColorPartons
      std::list<G4QParton*>::iterator picp = projCP.begin();
      std::list<G4QParton*>::iterator pecp = projCP.end();
      std::list<G4QParton*>::iterator piac = projAC.begin();
      std::list<G4QParton*>::iterator peac = projAC.end();
      std::list<G4QParton*>::iterator ticp = targCP.begin();
      std::list<G4QParton*>::iterator tecp = targCP.end();
      std::list<G4QParton*>::iterator tiac = targAC.begin();
      std::list<G4QParton*>::iterator teac = targAC.end();
      for(; picp!=pecp&& piac!=peac&& ticp!=tecp&& tiac!=teac; ++picp,++piac,++ticp,++tiac)
      {
        pSP+=(*picp)->Get4Momentum();
        pSP+=(*piac)->Get4Momentum();
        tSP+=(*ticp)->Get4Momentum();
        tSP+=(*tiac)->Get4Momentum();
      }
      G4cout<<"-EMC-G4QFragmentation::Breeder: Interaction#"<<i<<",dP4M="
            <<projH->Get4Momentum()-pSP<<",dT4M="<<targH->Get4Momentum()-tSP<<G4endl;
#endif
    }  
    // 
    // >>>>>>>> make soft collisions (ordering is vital)
    //
    G4QInteractionVector::iterator i;
#ifdef debug
    G4cout<<"G4QFragmentation::Breeder: Creation ofSoftCollisionPartonPair STARTS"<<G4endl;
#endif
    for(i = theInteractions.begin(); i != theInteractions.end(); i++)   
    {
      G4QInteraction* anIniteraction = *i;
      G4QPartonPair*  aPair=0;
      G4int nSoftCollisions = anIniteraction->GetNumberOfSoftCollisions();
#ifdef debug
      G4cout<<"G4QFragmentation::Breeder: #0f SOFT collisions ="<<nSoftCollisions<<G4endl;
#endif
      if (nSoftCollisions)
      { 
        G4QHadron* pProjectile = anIniteraction->GetProjectile();
        G4QHadron* pTarget     = anIniteraction->GetTarget();
        for (G4int j = 0; j < nSoftCollisions; j++)
        {
          aPair = new G4QPartonPair(pTarget->GetNextParton(),
                                    pProjectile->GetNextAntiParton(),
                                    G4QPartonPair::SOFT, G4QPartonPair::TARGET);
          thePartonPairs.push_back(aPair); // A target pair (Why TAGRET?)
          aPair = new G4QPartonPair(pProjectile->GetNextParton(),
                                    pTarget->GetNextAntiParton(),
                                    G4QPartonPair::SOFT, G4QPartonPair::PROJECTILE);
          thePartonPairs.push_back(aPair); // A projectile pair (Why Projectile?)
#ifdef debug
          G4cout<<"G4QFragmentation::Breeder: SOFT, 2 parton pairs're filled"<<G4endl;
#endif
        }  
        delete *i;
        i=theInteractions.erase(i);       // Soft interactions are converted & erased
        i--;
      }
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Breeder: -> Parton pairs for SOFT strings are made"<<G4endl;
#endif  
		  //
    // >>>>>>>>>>>>>>> make the rest as the diffractive interactions
    //
    for(unsigned i = 0; i < theInteractions.size(); i++) // Interactions are reduced bySoft
    {
      // The double or single diffraction is defined by the presonce of proj/targ partons
      G4QInteraction* anIniteraction = theInteractions[i];
      G4QPartonPair* aPartonPair;
#ifdef debug
      G4cout<<"G4QFragmentation::Breeder: CreationOfDiffractivePartonPairs, i="<<i<<G4endl;
#endif
      // the projectile diffraction parton pair is created first
      G4QHadron* aProjectile = anIniteraction->GetProjectile();
      G4QParton* aParton = aProjectile->GetNextParton();
      if (aParton)
      {
        aPartonPair=new G4QPartonPair(aParton, aProjectile->GetNextAntiParton(), 
                                      G4QPartonPair::DIFFRACTIVE,
                                      G4QPartonPair::PROJECTILE);
        thePartonPairs.push_back(aPartonPair);
#ifdef debug
        G4cout<<"G4QFragmentation::Breeder: proj Diffractive PartonPair is filled"<<G4endl;
#endif
      }
      // then the target diffraction parton pair is created
      G4QHadron* aTarget = anIniteraction->GetTarget();
      aParton = aTarget->GetNextParton();
      if (aParton)
      {
        aPartonPair = new G4QPartonPair(aParton, aTarget->GetNextAntiParton(), 
                                        G4QPartonPair::DIFFRACTIVE, G4QPartonPair::TARGET);
        thePartonPairs.push_back(aPartonPair);
#ifdef debug
        G4cout<<"G4QFragmentation::Breeder: targ Diffractive PartonPair is filled"<<G4endl;
#endif
      }
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Breeder:DiffractivePartonPairs are created"<<G4endl;
#endif  
    //
    // >>>>>>>>>>>>>> clean-up  Interactions and cmProjectile, if necessary
    //
    std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
    theInteractions.clear();
    delete cmProjectile;
#ifdef debug
    G4cout<<"G4QFragmentation::Breeder: Temporary objects are cleaned up"<<G4endl;
#endif  
    // This function prepares theBoost for transformation of secondaries to LS (-ProjRot!)
    theNucleus->DoLorentzBoost(theCurrentVelocity);// Boost theResidualNucleus to RotatedLS
    // @@ Nucleus isn't completely in LS, it needs the toZ (-ProjRot) rotation to consE/M
#ifdef pdebug
    G4cout<<"G4QFragmentation::Breeder: >>>>>> Strings are created "<<G4endl;
#endif
    G4QPartonPair* aPair;
    strings = new G4QStringVector;
    G4QString* aString=0;
    while(thePartonPairs.size()) // @@ At present noDifference in stringBuild (? M.K.)
    {
      aPair = thePartonPairs.back();           // Get the parton pair
      thePartonPairs.pop_back();               // Clean up the parton pair in the vector
#ifdef debug
      G4cout<<"G4QFragmentation::Breeder: StringType="<<aPair->GetCollisionType()<<G4endl;
#endif
      if (aPair->GetCollisionType() == G4QPartonPair::DIFFRACTIVE)
      {
#ifdef debug
        G4cout<<"G4QFragmentation::Breeder:Prepared for creationOfDiffExcitation"<<G4endl;
#endif
        aString = BuildString(aPair);           // Just to shorten the new
#ifdef debug
        G4cout<<"G4QFragmentation::Breeder:DifString4M="<<aString->Get4Momentum()<<G4endl;
#endif
      }
      else
      {
#ifdef debug
        G4cout<<"G4QFragmentation::Breeder: Prepared for creationOfTheSoftString"<<G4endl;
#endif
        aString = BuildString(aPair);           // Just to shorten the new
#ifdef debug
        G4cout<<"G4QFragmentation::Breeder:SoftString4M="<<aString->Get4Momentum()<<G4endl;
#endif
      }
      aString->Boost(theCurrentVelocity);       // ! Strings are moved to ZLS when pushed !
      strings->push_back(aString);
      delete aPair;
    } // End of the String Creation LOOP
#ifdef edebug
    G4LorentzVector sum=theNucleus->Get4Momentum();// Nucleus 4Mom in rotatedLS
    G4int rChg=totChg-theNucleus->GetZ();
    G4int rBaN=totBaN-theNucleus->GetA();
    G4int nStrings=strings->size();
    G4cout<<"-EMC-G4QFragmentation::Breeder:#ofString="<<nStrings<<",tNuc4M="<<sum<<G4endl;
    for(G4int i=0; i<nStrings; i++)
    {
      G4LorentzVector strI4M=(*strings)[i]->Get4Momentum();
      sum+=strI4M;
      G4int sChg=(*strings)[i]->GetCharge();
      G4int sBaN=(*strings)[i]->GetBaryonNumber();
      G4int LPDG=(*strings)[i]->GetLeftParton()->GetPDGCode();
      G4int RPDG=(*strings)[i]->GetRightParton()->GetPDGCode();
      G4QContent LQC=(*strings)[i]->GetLeftParton()->GetQC();
      G4QContent RQC=(*strings)[i]->GetRightParton()->GetQC();
      rChg-=sChg;
      rBaN-=sBaN;
      G4cout<<"-EMC-G4QFragmentation::Breeder: String#"<<i<<", 4M="<<strI4M<<",LPDG="<<LPDG
            <<LQC<<",RPDG="<<RPDG<<RQC<<", Ch="<<sChg<<", BN="<<sBaN<<G4endl;
    }
    G4cout<<"-EMC-G4QFragm::Breed: r4M="<<sum-totZLS4M<<",rC="<<rChg<<",rB="<<rBaN<<G4endl;
#endif
  } // End of the LOOP of the strings creation
#ifdef debug
  G4cout<<"G4QFragmentation::Breeder: BeforeRotation #0fStrings="<<strings->size()<<G4endl;
#endif
  //
  // ------------------ At this point the strings must be created -----------------
  //
  for(unsigned astring=0; astring < strings->size(); astring++)
          (*strings)[astring]->LorentzRotate(toLab); // Recove Z-direction in LS ("LS"->LS)
  theNucleus->DoLorentzRotation(toLab); // Recove Z-direction in LS ("LS"->LS) for rNucleus
  // Now everything is in LS system
#ifdef edebug
  // @@ Should be tested for the primary projectile not along Z !
  G4LorentzVector sum=theNucleus->Get4Momentum();    // Nucleus 4Mom in LS
  G4int rChg=totChg-theNucleus->GetZ();
  G4int rBaN=totBaN-theNucleus->GetA();
  G4int nStrings=strings->size();
  G4cout<<"-EMCLS-G4QFragmentation::Breeder:#ofS="<<nStrings<<",tNuc4M(E=M)="<<sum<<G4endl;
  for(G4int i=0; i<nStrings; i++)
  {
    G4LorentzVector strI4M=(*strings)[i]->Get4Momentum();
    sum+=strI4M;
    G4int sChg=(*strings)[i]->GetCharge();
    rChg-=sChg;
    G4int sBaN=(*strings)[i]->GetBaryonNumber();
    rBaN-=sBaN;
    G4cout<<"-EMCLS-G4QFragm::Breeder:String#"<<i<<",4M="<<strI4M<<strI4M.m()<<",Charge="
          <<sChg<<",BaryN="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFragm::Breeder:r4M="<<sum-totLS4M<<",rC="<<rChg<<",rB="<<rBaN<<G4endl;
#endif
  //
  // --- Strings are created, but we should get rid of too light strings (Mmin+MPi0) -----
  //
  G4QStringVector::iterator ist;
  G4int problem=0;                                   // 0="no problem", incremented by ASIS
  for(ist = strings->begin(); ist < strings->end(); ist++)
  {
    G4bool bad=true;
    G4LorentzVector cS4M=(*ist)->Get4Momentum();
    G4double cSM2=cS4M.m2();                         // Squared mass of the String
    G4QParton* cLeft=(*ist)->GetLeftParton();
    G4QParton* cRight=(*ist)->GetRightParton();
    G4int cLT=cLeft->GetType();
    G4int cRT=cRight->GetType();
    G4int cLPDG=cLeft->GetPDGCode();
    G4int cRPDG=cRight->GetPDGCode();
    G4int aLPDG=0;
    G4int aRPDG=0;
    if     (cLPDG > 7) aLPDG=  cLPDG/100;
    else if(cLPDG <-7) aLPDG=(-cLPDG)/100;
    if     (cRPDG > 7) aRPDG=  cRPDG/100;
    else if(cRPDG <-7) aRPDG=(-cRPDG)/100;
    G4int L1=0;
    G4int L2=0;
    if(aLPDG)
    {
      L1=aLPDG/10;
      L2=aLPDG%10;
    }
    G4int R1=0;
    G4int R2=0;
    if(aRPDG)
    {
      R1=aRPDG/10;
      R2=aRPDG%10;
    }
#ifdef debug
    G4double cSM=0.;
    if(cSM2>0.) cSM=std::sqrt(cSM2);
    G4cout<<"G4QFragmentation::Breed:NeedsFusion? cLPDG="<<cLPDG<<",cRPDG="<<cRPDG<<",cM2="
          <<cSM2<<",cM="<<cSM<<G4endl;
#endif
    if(cSM2>0.)                                      // Mass can be calculated
    {
      G4bool single=true;
      G4double miM=0.;                               // Proto of the Min HadronString Mass
      if(cLT==2 && cRT==2)
      {
        if(L1!=R1 && L1!=R2 && L2!=R1 && L2!=R2)     // Unreducable DiQ-aDiQ
        {
          single=false;
          G4QPDGCode tmp;
          std::pair<G4int,G4int> paB=tmp.MakeTwoBaryons(L1, L2, R1, R2);
          miM=G4QPDGCode(paB.first).GetMass()+G4QPDGCode(paB.second).GetMass();
        }
      }
	if(single) miM=G4QPDGCode((*ist)->GetQC().GetSPDGCode()).GetMass() + mPi0;//MinHSMass
#ifdef debug
      G4cout<<"G4QFrag::Breed:*IsItGood? realM="<<std::sqrt(cSM2)<<" > GSM="<<miM<<G4endl;
#endif
      if(std::sqrt(cSM2) > miM) bad=false;           // String is OK
    }
    if(bad)                                          // String should be merged with others
    {
#ifdef debug
      G4cout<<"G4QFrag::Breed:TryFuse,L1="<<L1<<",L2="<<L2<<",R1="<<R1<<",R2="<<R2<<G4endl;
#endif
      G4int cST=cLT+cRT;
      G4double excess=-DBL_MAX;                      // The value to be maximized excess M
      G4double maxiM2=-DBL_MAX;                      // The value to be maximized M2
      G4QStringVector::iterator sst;                 // Selected partner string
      G4QStringVector::iterator pst;
      G4int sLPDG=0;                                 // selectedLeft (like inStringPartner)
      G4int sRPDG=0;                                 // selectedRight(like inStringPartner)
      G4int sOrd=0;                                  // selected Order of PartonFusion
      G4bool minC=true;                              // for the case when M2<0
      if(cSM2>0.) minC=false;                        // If M2>0 already don'tSearchFor M2>0
      for(pst = strings->begin(); pst < strings->end(); pst++) if(pst != ist)
      {
        G4LorentzVector pS4M=(*pst)->Get4Momentum()+cS4M; // Summed 4-momentum
        G4int nLPDG=0;                               // new Left (like in theStringPartner)
        G4int nRPDG=0;                               // new Right(like in theStringPartner)
        G4double pExcess=-DBL_MAX;                   // Prototype of the excess
        G4double pSM2=pS4M.m2();                     // Squared mass of the Fused Strings
#ifdef debug
        G4cout<<"->G4QFragm::Breeder: sum4M="<<pS4M<<", M2="<<pSM2<<G4endl;
#endif
        //if(pSM2>0.)                                  // The partner can be a candidate
        //{
        G4QParton* pLeft=(*pst)->GetLeftParton();
        G4QParton* pRight=(*pst)->GetRightParton();
        G4int pLT=pLeft->GetType();
        G4int pRT=pRight->GetType();
        G4int pLPDG=pLeft->GetPDGCode();
        G4int pRPDG=pRight->GetPDGCode();
        G4int LPDG=0;
        G4int RPDG=0;
        if     (pLPDG > 7) LPDG=  pLPDG/100;
        else if(pLPDG <-7) LPDG=(-pLPDG)/100;
        if     (pRPDG > 7) RPDG=  pRPDG/100;
        else if(pRPDG <-7) RPDG=(-pRPDG)/100;
        G4int pL1=0;
        G4int pL2=0;
        if(LPDG)
        {
          pL1=LPDG/10;
          pL2=LPDG%10;
        }
        G4int pR1=0;
        G4int pR2=0;
        if(RPDG)
        {
          pR1=RPDG/10;
          pR2=RPDG%10;
        }
        G4int pST=pLT+pRT;
#ifdef debug
        G4cout<<"G4QFragm::Breeder: Partner/w pLPDG="<<pLPDG<<", pRPDG="<<pRPDG<<", pM2="
              <<pSM2<<G4endl;
#endif
        // Different fromCompactAlrorithm ofStringFusionAfterDecay (no DiQaDiQ reduction)
        G4int tf=0;                                // Type combination flag
        G4int af=0;                                // Annihilatio combination flag
        if     (cST==2 && pST==2)                  // QaQ + QaQ = DiQaDiQ (always)
        {
          tf=1;
          af=1;
        }
        else if(cST==2 && pST==3)                  // QaQ + QDiQ/aQaDiQ = QDiQ/aQaDiQ (s)
        {
          tf=2;
          if     (pLPDG > 7 &&
                  ( (cLPDG<0 && (-cLPDG==pL1 || -cLPDG==pL2) ) ||
                    (cRPDG<0 && (-cRPDG==pL1 || -cRPDG==pL2) )
                  )
                 ) af=1;
          else if(pRPDG > 7 &&
                  ( (cLPDG<0 && (-cLPDG==pR1 || -cLPDG==pR2) ) ||
                    (cRPDG<0 && (-cRPDG==pR1 || -cRPDG==pR2) )
                  )
                 ) af=2;
          else if(pLPDG <-7 &&
                  ( (cLPDG>0 && ( cLPDG==pL1 || cLPDG==pL2) ) ||
                    (cRPDG>0 && ( cRPDG==pL1 || cRPDG==pL2) )
                  )
                 ) af=3;
          else if(pRPDG <-7 &&
                  ( (cLPDG>0 && ( cLPDG==pR1 || cLPDG==pR2) ) ||
                    (cRPDG>0 && ( cRPDG==pR1 || cRPDG==pR2) )
                  )
                 ) af=4;
#ifdef debug
          else G4cout<<"G4QFragmentation::Breeder:2(QaQ+QDiQ/aQaDiQ) Can't fuse"<<G4endl;
#endif
        }
        else if(cST==3 && pST==2)                  // QDiQ/aQaDiQ + QaQ = QDiQ/aQaDiQ (s)
        {
          tf=3;
          if     (cLPDG > 7 &&
                  ( (pLPDG<0 && (-pLPDG==L1 || -pLPDG==L2) ) ||
                    (pRPDG<0 && (-pRPDG==L1 || -pRPDG==L2) )
                  )
                 ) af=1;
          else if(cRPDG > 7 &&
                  ( (pLPDG<0 && (-pLPDG==R1 || -pLPDG==R2) ) ||
                    (pRPDG<0 && (-pRPDG==R1 || -pRPDG==R2) )
                  )
                 ) af=2;
          else if(cLPDG <-7 &&
                  ( (pLPDG>0 && ( pLPDG==L1 || pLPDG==L2) ) ||
                    (pRPDG>0 && ( pRPDG==L1 || pRPDG==L2) )
                  )
                 ) af=3;
          else if(cRPDG <-7 &&
                  ( (pLPDG>0 && ( pLPDG==R1 || pLPDG==R2) ) ||
                    (pRPDG>0 && ( pRPDG==R1 || pRPDG==R2) )
                  )
                 ) af=4;
#ifdef debug
          else G4cout<<"G4QFragmentation::Breeder:3(QDiQ/aQaDiQ+QaQ) Can't fuse"<<G4endl;
#endif
        }
        else if(cST==2 && pST==4)                  // QaQ + aDiQDiQ = QaQ (double)
        {
          tf=4;
          if     (pLPDG > 7) // pRPDG <-7
          {
            if     ( (-cLPDG==pL1 || -cLPDG==pL2) && (cRPDG==pR1 || cRPDG==pR2) ) af=1;
            else if( (-cRPDG==pL1 || -cRPDG==pL2) && (cLPDG==pR1 || cLPDG==pR2) ) af=2;
          }
          else if(pRPDG > 7) // pLPDG <-7
          {
            if     ( (-cRPDG==pR1 || -cRPDG==pR2) && (cLPDG==pL1 || cLPDG==pL2) ) af=3;
            else if( (-cLPDG==pR1 || -cLPDG==pR2) && (cRPDG==pL1 || cRPDG==pL2) ) af=4;
          }
#ifdef debug
          else G4cout<<"-G4QFragmentation::Breeder: 4 (QaQ+aQDiQDiQ) Can't fuse"<<G4endl;
#endif
        }
        else if(cST==4 && pST==2)                  // aDiQDiQ + QaQ = QaQ (double)
        {
          tf=5;
          if     (cLPDG > 7) // cRPDG<-7
          {
            if     ( (-pLPDG==L1 || -pLPDG==L2) && (pRPDG==R1 || pRPDG==R2) ) af=1;
            else if( (-pRPDG==L1 || -pRPDG==L2) && (pLPDG==R1 || pLPDG==R2) ) af=2;
          }
          else if(cRPDG > 7) // cLPDG<-7
          {
            if     ( (-pRPDG==R1 || -pRPDG==R2) && (pLPDG==L1 || pLPDG==L2) ) af=3;
            else if( (-pLPDG==R1 || -pLPDG==R2) && (pRPDG==L1 || pRPDG==L2) ) af=4;
          }
#ifdef debug
          else G4cout<<"-G4QFragmentation::Breeder: 5 (aQDiQDiQ+QaQ) Can't fuse"<<G4endl;
#endif
        }
        else if(cST==3 && pST==3)                  // QDiQ + aQaDiQ = QaQ (double)
        {
          tf=6;
          if(pLPDG > 7)
          {
            if     (cLPDG<-7 && (-cRPDG==pL1 || -cRPDG==pL2) && (pRPDG==L1 || pRPDG==L2))
              af=1;
            else if(cRPDG<-7 && (-cLPDG==pL1 || -cLPDG==pL2) && (pRPDG==R1 || pRPDG==R2))
              af=2;
          }
          else if(pRPDG > 7)
          {
            if     (cLPDG<-7 && (-cRPDG==pR1 || -cRPDG==pR2) && (pLPDG==L1 || pLPDG==L2))
              af=3;
            else if(cRPDG<-7 && (-cLPDG==pR1 || -cLPDG==pR2) && (pLPDG==R1 || pLPDG==R2))
              af=4;
          }
          else if(cLPDG > 7)
          {
            if     (pLPDG<-7 && (-pRPDG==L1 || -pRPDG==L2) && (cRPDG==pL1 || cRPDG==pL2))
              af=5;
            else if(pRPDG<-7 && (-pLPDG==L1 || -pLPDG==L2) && (cRPDG==pR1 || cRPDG==pR2))
              af=6;
          }
          else if(cRPDG > 7)
          {
            if     (pLPDG<-7 && (-pRPDG==R1 || -pRPDG==R2) && (cLPDG==pL1 || cLPDG==pL2))
              af=7;
            else if(pRPDG<-7 && (-pLPDG==R1 || -pLPDG==R2) && (cLPDG==pR1 || cLPDG==pR2))
              af=8;
          }
#ifdef debug
          else G4cout<<"-G4QFragmentation::Breeder: 6 (QDiQ+aQaDiQ) Can't fuse"<<G4endl;
#endif
        }
#ifdef debug
        G4cout<<"G4QFrag::Breed: ***Possibility***, tf="<<tf<<", af="<<af<<G4endl;
#endif
        if(tf && af)
        {
          // Strings can be fused, update the max excess and remember usefull switches
          G4int order=0;                           // LL/RR (1) or LR/RL (2) PartonFusion
          switch (tf)
          {
            case 1: // ------------------------------------> QaQ + QaQ = DiQaDiQ (always)
              if     (cLPDG > 0 && pLPDG > 0)
              {
                order= 1;
                if     (cLPDG > pLPDG) nLPDG=cLPDG*1000+pLPDG*100+1;
                else if(cLPDG < pLPDG) nLPDG=pLPDG*1000+cLPDG*100+1;
                else                   nLPDG=pLPDG*1000+cLPDG*100+3;
                if     (cRPDG < pRPDG) nRPDG=cRPDG*1000+pRPDG*100-1;
                else if(cRPDG > pRPDG) nRPDG=pRPDG*1000+cRPDG*100-1;
                else                   nRPDG=pRPDG*1000+cRPDG*100-3;
              }
              else if(cLPDG < 0 && pLPDG < 0)
              {
                order= 1;
                if     (cRPDG > pRPDG) nRPDG=cRPDG*1000+pRPDG*100+1;
                else if(cRPDG < pRPDG) nRPDG=pRPDG*1000+cRPDG*100+1;
                else                   nRPDG=pRPDG*1000+cRPDG*100+3;
                if     (cLPDG < pLPDG) nLPDG=cLPDG*1000+pLPDG*100-1;
                else if(cLPDG > pLPDG) nLPDG=pLPDG*1000+cLPDG*100-1;
                else                   nLPDG=pLPDG*1000+cLPDG*100-3;
              }
              else if(cRPDG > 0 && pLPDG > 0)
              {
                order=-1;
                if     (cRPDG > pLPDG) nLPDG=cRPDG*1000+pLPDG*100+1;
                else if(cRPDG < pLPDG) nLPDG=pLPDG*1000+cRPDG*100+1;
                else                   nLPDG=pLPDG*1000+cRPDG*100+3;
                if     (cLPDG < pRPDG) nRPDG=cLPDG*1000+pRPDG*100-1;
                else if(cLPDG > pRPDG) nRPDG=pRPDG*1000+cLPDG*100-1;
                else                   nRPDG=pRPDG*1000+cLPDG*100-3;
              }
              else if(cRPDG < 0 && pLPDG < 0)
              {
                order=-1;
                if     (cLPDG > pRPDG) nRPDG=cLPDG*1000+pRPDG*100+1;
                else if(cLPDG < pRPDG) nRPDG=pRPDG*1000+cLPDG*100+1;
                else                   nRPDG=pRPDG*1000+cLPDG*100+3;
                if     (cRPDG < pLPDG) nLPDG=cRPDG*1000+pLPDG*100-1;
                else if(cRPDG > pLPDG) nLPDG=pLPDG*1000+cRPDG*100-1;
                else                   nLPDG=pLPDG*1000+cRPDG*100-3;
              }
              break;
            case 2: // ------------------------> QaQ + QDiQ/aQaDiQ = QDiQ/aQaDiQ (single)
             switch (af)
             {
               case 1: // ....................... pLPDG > 7
                 if(cLPDG < 0)
                 {
                   order= 1;
                   if     (cRPDG > pRPDG) nRPDG=cRPDG*1000+pRPDG*100+1;
                   else if(cRPDG < pRPDG) nRPDG=pRPDG*1000+cRPDG*100+1;
                   else                   nRPDG=pRPDG*1000+cRPDG*100+3;
                   if  (-cLPDG == pL1)    nLPDG=pL2;
                   else                   nLPDG=pL1; // -cLPDG == pL2
                 }
                 else // cRPDG < 0
                 {
                   order=-1;
                   if     (cLPDG > pRPDG) nRPDG=cLPDG*1000+pRPDG*100+1;
                   else if(cLPDG < pRPDG) nRPDG=pRPDG*1000+cLPDG*100+1;
                   else                   nRPDG=pRPDG*1000+cLPDG*100+3;
                   if  (-cRPDG == pL1)    nLPDG=pL2;
                   else                   nLPDG=pL1; // -cRPDG == pL2
                 }
                 break;
               case 2: // ....................... pRPDG > 7
                 if(cLPDG < 0)
                 {
                   order=-1;
                   if     (cRPDG > pLPDG) nLPDG=cRPDG*1000+pLPDG*100+1;
                   else if(cRPDG < pLPDG) nLPDG=pLPDG*1000+cRPDG*100+1;
                   else                   nLPDG=pLPDG*1000+cRPDG*100+3;
                   if  (-cLPDG == pR1)    nRPDG=pR2;
                   else                   nRPDG=pR1; // -cLPDG == pR2
                 }
                 else // cRPDG < 0
                 {
                   order= 1;
                   if     (cLPDG > pLPDG) nLPDG=cLPDG*1000+pLPDG*100+1;
                   else if(cLPDG < pLPDG) nLPDG=pLPDG*1000+cLPDG*100+1;
                   else                   nLPDG=pLPDG*1000+cLPDG*100+3;
                   if  (-cRPDG == pR1)    nRPDG=pR2;
                   else                   nRPDG=pR1; // -cRPDG == pR2
                 }
                 break;
               case 3: // ....................... pLPDG <-7
                 if(cLPDG > 0)
                 {
                   order= 1;
                   if     (cRPDG < pRPDG) nRPDG=cRPDG*1000+pRPDG*100-1;
                   else if(cRPDG > pRPDG) nRPDG=pRPDG*1000+cRPDG*100-1;
                   else                   nRPDG=pRPDG*1000+cRPDG*100-3;
                   if  ( cLPDG == pL1)    nLPDG=-pL2;
                   else                   nLPDG=-pL1; // cLPDG == pL2
                 }
                 else // cRPDG > 0
                 {
                   order=-1;
                   if     (cLPDG < pRPDG) nRPDG=cLPDG*1000+pRPDG*100-1;
                   else if(cLPDG > pRPDG) nRPDG=pRPDG*1000+cLPDG*100-1;
                   else                   nRPDG=pRPDG*1000+cLPDG*100-3;
                   if  ( cRPDG == pL1)    nLPDG=-pL2;
                   else                   nLPDG=-pL1; // cRPDG == pL2
                 }
                 break;
               case 4: // ....................... pLRDG <-7
                 if(cLPDG > 0)
                 {
                   order=-1;
                   if     (cRPDG < pLPDG) nLPDG=cRPDG*1000+pLPDG*100-1;
                   else if(cRPDG > pLPDG) nLPDG=pRPDG*1000+cLPDG*100-1;
                   else                   nLPDG=pLPDG*1000+cRPDG*100-3;
                   if  ( cLPDG == pR1)    nRPDG=-pR2;
                   else                   nRPDG=-pR1; // cLPDG == pR2
                 }
                 else // cRPDG > 0
                 {
                   order= 1;
                   if     (cLPDG < pLPDG) nLPDG=cLPDG*1000+pLPDG*100-1;
                   else if(cLPDG > pLPDG) nLPDG=pLPDG*1000+cLPDG*100-1;
                   else                   nLPDG=pLPDG*1000+cLPDG*100-3;
                   if  ( cRPDG == pR1)    nRPDG=-pR2;
                   else                   nRPDG=-pR1; // cRPDG == pR2
                 }
                 break;
             }
             break;
            case 3: // ------------------------> QDiQ/aQaDiQ + QaQ = QDiQ/aQaDiQ (single)
             switch (af)
             {
               case 1: // ....................... cLPDG > 7
                 if(pLPDG < 0)
                 {
                   order= 1;
                   if     (pRPDG > cRPDG) nRPDG=pRPDG*1000+cRPDG*100+1;
                   else if(pRPDG < cRPDG) nRPDG=cRPDG*1000+pRPDG*100+1;
                   else                   nRPDG=cRPDG*1000+pRPDG*100+3;
                   if  (-pLPDG == L1)     nLPDG=L2;
                   else                   nLPDG=L1; // -pLPDG == L2
                 }
                 else // pRPDG < 0
                 {
                   order=-1;
                   if     (pLPDG > cRPDG) nLPDG=pLPDG*1000+cRPDG*100+1;
                   else if(pLPDG < cRPDG) nLPDG=cRPDG*1000+pLPDG*100+1;
                   else                   nLPDG=cRPDG*1000+pLPDG*100+3;
                   if  (-pRPDG == L1)     nRPDG=L2;
                   else                   nRPDG=L1; // -pRPDG == L2
                 }
                 break;
               case 2: // ....................... cRPDG > 7
                 if(pLPDG < 0)
                 {
                   order=-1;
                   if     (pRPDG > cLPDG) nRPDG=pRPDG*1000+cLPDG*100+1;
                   else if(pRPDG < cLPDG) nRPDG=cLPDG*1000+pRPDG*100+1;
                   else                   nRPDG=cLPDG*1000+pRPDG*100+3;
                   if  (-pLPDG == R1)     nLPDG=R2;
                   else                   nLPDG=R1; // -pLPDG == R2
                 }
                 else // pRPDG < 0
                 {
                   order= 1;
                   if     (pLPDG > cLPDG) nLPDG=pLPDG*1000+cLPDG*100+1;
                   else if(pLPDG < cLPDG) nLPDG=cLPDG*1000+pLPDG*100+1;
                   else                   nLPDG=cLPDG*1000+pLPDG*100+3;
                   if  (-pRPDG == R1)     nRPDG=R2;
                   else                   nRPDG=R1; // -pRPDG == R2
                 }
                 break;
               case 3: // ....................... cLPDG <-7
                 if(pLPDG > 0)
                 {
                   order= 1;
                   if     (pRPDG < cRPDG) nRPDG=pRPDG*1000+cRPDG*100-1;
                   else if(pRPDG > cRPDG) nRPDG=cRPDG*1000+pRPDG*100-1;
                   else                   nRPDG=cRPDG*1000+pRPDG*100-3;
                   if  ( pLPDG == L1)     nLPDG=-L2;
                   else                   nLPDG=-L1; // pLPDG == L2
                 }
                 else // pRPDG > 0
                 {
                   order=-1;
                   if     (pLPDG < cRPDG) nLPDG=pLPDG*1000+cRPDG*100-1;
                   else if(pLPDG > cRPDG) nLPDG=cRPDG*1000+pLPDG*100-1;
                   else                   nLPDG=cRPDG*1000+pLPDG*100-3;
                   if  ( pRPDG == L1)     nRPDG=-L2;
                   else                   nRPDG=-L1; // pRPDG == L2
                 }
                 break;
               case 4: // ....................... cLRDG <-7
                 if(pLPDG > 0)
                 {
                   order=-1;
                   if     (pRPDG < cLPDG) nRPDG=pRPDG*1000+cLPDG*100-1;
                   else if(pRPDG > cLPDG) nRPDG=cRPDG*1000+pLPDG*100-1;
                   else                   nRPDG=cLPDG*1000+pRPDG*100-3;
                   if  ( pLPDG == R1)     nLPDG=-R2;
                   else                   nLPDG=-R1; // pLPDG == R2
                 }
                 else // pRPDG > 0
                 {
                   order= 1;
                   if     (pLPDG < cLPDG) nLPDG=pLPDG*1000+cLPDG*100-1;
                   else if(pLPDG > cLPDG) nLPDG=cLPDG*1000+pLPDG*100-1;
                   else                   nLPDG=cLPDG*1000+pLPDG*100-3;
                   if  ( pRPDG == R1)     nRPDG=-R2;
                   else                   nRPDG=-R1; // pRPDG == R2
                 }
                 break;
             }
             break;
            case 4: // ------------------------------------> QaQ + aDiQDiQ = QaQ (double)
             switch (af)
             {
               case 1: // ....................... pLPDG > 7 && pRPDG <-7 === LL/RR
                 order= 1;
                 if(-cLPDG == pL1) nLPDG= pL2;
                 else              nLPDG= pL1;
                 if( cRPDG == pR1) nRPDG=-pR2;
                 else              nRPDG=-pR1;
                 break;
               case 2: // ...................... pLPDG > 7 && pRPDG <-7 === LR/RL
                 order=-1;
                 if(-cRPDG == pL1) nLPDG= pL2;
                 else              nLPDG= pL1;
                 if( cLPDG == pR1) nRPDG=-pR2;
                 else              nRPDG=-pR1;
                 break;
               case 3: // ...................... pRPDG > 7 && pLPDG <-7 === LL/RR
                 order= 1;
                 if( cLPDG == pL1) nLPDG=-pL2;
                 else              nLPDG=-pL1;
                 if(-cRPDG == pR1) nRPDG= pR2;
                 else              nRPDG= pR1;
                 break;
               case 4: // ...................... pRPDG > 7 && pLPDG <-7 === LR/RL
                 order=-1;
                 if( cRPDG == pL1) nLPDG=-pL2;
                 else              nLPDG=-pL1;
                 if(-cLPDG == pR1) nRPDG= pR2;
                 else              nRPDG= pR1;
                 break;
             }
             break;
            case 5: // ------------------------------------> aDiQDiQ + QaQ = QaQ (double)
             switch (af)
             {
               case 1: // ...................... cLPDG > 7 && cRPDG <-7 === LL/RR
                 order= 1;
                 if(-pLPDG == L1) nLPDG= L2;
                 else             nLPDG= L1;
                 if( pRPDG == R1) nRPDG=-R2;
                 else             nRPDG=-R1;
                 break;
               case 2: // ...................... cLPDG > 7 && cRPDG <-7 === LR/RL
                 order=-1;
                 if(-pRPDG == L1) nRPDG= L2;
                 else             nRPDG= L1;
                 if( pLPDG == R1) nLPDG=-R2;
                 else             nLPDG=-R1;
                 break;
               case 3: // ...................... cRPDG > 7 && cLPDG <-7 === LL/RR
                 order= 1;
                 if( pLPDG == L1) nLPDG=-L2;
                 else             nLPDG=-L1;
                 if(-pRPDG == R1) nRPDG= R2;
                 else             nRPDG= R1;
                 break;
               case 4: // ...................... cRPDG > 7 && cLPDG <-7 === LR/RL
                 order=-1;
                 if( pRPDG == L1) nRPDG=-L2;
                 else             nRPDG=-L1;
                 if(-pLPDG == R1) nLPDG= R2;
                 else             nLPDG= R1;
                 break;
             }
             break;
            case 6: // ------------------------------------> QDiQ + aQaDiQ = QaQ (double)
             switch (af)
             {
               case 1:
                 order=-1;
                 if(-cRPDG == pL1) nLPDG= pL2;
                 else              nLPDG= pL1;
                 if( pRPDG ==  L1) nRPDG= -L2;
                 else              nRPDG= -L1;
                 break;
               case 2:
                 order= 1;
                 if(-cLPDG == pL1) nLPDG= pL2;
                 else              nLPDG= pL1;
                 if( pRPDG ==  R1) nRPDG= -R2;
                 else              nRPDG= -R1;
                 break;
               case 3:
                 order= 1;
                 if(-cRPDG == pR1) nRPDG= pR2;
                 else              nRPDG= pR1;
                 if( pLPDG ==  L1) nLPDG= -L2;
                 else              nLPDG= -L1;
                 break;
               case 4:
                 order=-1;
                 if(-cLPDG == pR1) nRPDG= pR2;
                 else              nRPDG= pR1;
                 if( pLPDG ==  R1) nLPDG= -R2;
                 else              nLPDG= -R1;
                 break;
               case 5:
                 order=-1;
                 if(-pRPDG ==  L1) nRPDG=  L2;
                 else              nRPDG=  L1;
                 if( cRPDG == pL1) nLPDG=-pL2;
                 else              nLPDG=-pL1;
                 break;
               case 6:
                 order= 1;
                 if(-pLPDG ==  L1) nLPDG=  L2;
                 else              nLPDG=  L1;
                 if( cRPDG == pR1) nRPDG=-pR2;
                 else              nRPDG=-pR1;
                 break;
               case 7:
                 order= 1;
                 if(-pRPDG ==  R1) nRPDG=  R2;
                 else              nRPDG=  R1;
                 if( cLPDG == pL1) nLPDG=-pL2;
                 else              nLPDG=-pL1;
                 break;
               case 8:
                 order=-1;
                 if(-pLPDG ==  R1) nLPDG=  R2;
                 else              nLPDG=  R1;
                 if( cLPDG == pR1) nRPDG=-pR2;
                 else              nRPDG=-pR1;
                 break;
             }
             break;
          }
          if(!order) G4cerr<<"-Warning-G4QFrag::Breed: t="<<tf<<", a="<<af<<", cL="<<cLPDG
                           <<", cR="<<cRPDG<<", pL="<<pLPDG<<", pR="<<pRPDG<<G4endl;
          else
          {
            // With theNewHypotheticalPartons the min mass must be calculated & compared
            G4int LT=1;
            if(std::abs(nLPDG) > 7) ++LT;
            G4int RT=1;
            if(std::abs(nRPDG) > 7) ++RT;
            G4double minM=0.;
            G4bool sing=true;
            if(cLT==2 && cRT==2)
            {
              G4int aLPDG=0;
              G4int aRPDG=0;
              if(cLPDG>0)
              {
                aLPDG=nLPDG/100;
                aRPDG=(-nRPDG)/100;
              }
              else //cRPDG>0
              {
                aRPDG=nRPDG/100;
                aLPDG=(-nLPDG)/100;
              }
              G4int nL1=aLPDG/10;
              G4int nL2=aLPDG%10;
              G4int nR1=aRPDG/10;
              G4int nR2=aRPDG%10;
              if(nL1!=nR1 && nL1!=nR2 && nL2!=nR1 && nL2!=nR2) // Unreducable DiQ-aDiQ
              {
                sing=false;
                G4QPDGCode tmp;
                std::pair<G4int,G4int> pB=tmp.MakeTwoBaryons(nL1, nL2, nR1, nR2);
                minM=G4QPDGCode(pB.first).GetMass()+G4QPDGCode(pB.second).GetMass();
              }
            }
            if(sing)
            {
              G4QContent newStQC(std::make_pair(nLPDG,nRPDG)); // NewString QuarkContent
              G4int minPDG=newStQC.GetSPDGCode(); // PDG of the Lightest Hadron=String
              minM=G4QPDGCode(minPDG).GetMass() + mPi0; // Min SingleHadron=String Mass
            }
            // Compare this mass
            G4bool win=false;
            G4double    pSM=0.;
            if(pSM2>0.) pSM=std::sqrt(pSM2);
            if(minC && pSM2 > maxiM2)             // Up to now any positive mass is good
            {
              maxiM2=pSM2;
              win=true;
            }
            else if(!minC || pSM > minM)
            {
              pExcess=pSM-minM;
              if(minC || pExcess > excess)
              {
                minC=false;
                excess=pExcess;
                win=true;
              }
            }
            if(win)
            {
              sst=pst;
              sLPDG=nLPDG;
              sRPDG=nRPDG;
              sOrd=order;
            }
          } // End of IF(new partons are created)
        } // End of IF(compatible partons)
        //} // End of positive squared mass of the fused string
      } // End of the LOOP over the possible partners (with the exclusive if for itself)
      if(sOrd)                                       // The best pStringCandidate was found
      {
        G4LorentzVector cL4M=cLeft->Get4Momentum();
        G4LorentzVector cR4M=cRight->Get4Momentum();
        G4QParton* pLeft=(*sst)->GetLeftParton();
        G4QParton* pRight=(*sst)->GetRightParton();
        G4LorentzVector pL4M=pLeft->Get4Momentum();
        G4LorentzVector pR4M=pRight->Get4Momentum();
#ifdef debug
        G4cout<<"G4QFragmentation::Breed:cS4M="<<cS4M<<" fused/w pS4M="<<pL4M+pR4M<<G4endl;
#endif
        if(sOrd>0)
        {
          pL4M+=cL4M;
          pR4M+=cR4M;
        }
        else
        {
          pL4M+=cR4M;
          pR4M+=cL4M;
        }
        pLeft->SetPDGCode(sLPDG);
        pLeft->Set4Momentum(pL4M);
        pRight->SetPDGCode(sRPDG);
        pRight->Set4Momentum(pR4M);
        delete (*ist);
        strings->erase(ist);
        ist--;
#ifdef debug
        G4LorentzVector ss4M=pL4M+pR4M;
        G4cout<<"G4QFragmentation::Breeder: Created,4M="<<ss4M<<", m2="<<ss4M.m2()<<G4endl;
#endif
      } // End of the IF(the best partnerString candidate was found)
      else
      {
#ifdef debug
        G4cout<<"-Warning-G4QFrag::Breed:S4M="<<cS4M<<", M2="<<cSM2<<" Leave ASIS"<<G4endl;
#endif
        ++problem;
      }
    }
  }
#ifdef edebug
  // This print has meaning only if something appear between it and the StringFragmLOOP
  G4LorentzVector t4M=theNucleus->Get4Momentum();    // Nucleus 4Mom in LS
  G4int rC=totChg-theNucleus->GetZ();
  G4int rB=totBaN-theNucleus->GetA();
  G4int nStr=strings->size();
  G4cout<<"-EMCLS-G4QFr::Breed: AfterSUPPRESION #ofS="<<nStr<<",tNuc4M(E=M)="<<sum<<G4endl;
  for(G4int i=0; i<nStr; i++)
  {
    G4LorentzVector strI4M=(*strings)[i]->Get4Momentum();
    t4M+=strI4M;
    G4int sChg=(*strings)[i]->GetCharge();
    rC-=sChg;
    G4int sBaN=(*strings)[i]->GetBaryonNumber();
    rB-=sBaN;
    G4cout<<"-EMCLS-G4QFragm::Breeder: St#"<<i<<", 4M="<<strI4M<<", M="<<strI4M.m()<<", C="
          <<sChg<<", B="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFragm::Breeder:r4M="<<t4M-totLS4M<<",rC="<<rC<<",rB="<<rB<<G4endl;
#endif
  //
  // --- If a problem is foreseen then the DiQaDiQ strings should be reduced if possible --
  //
  if(problem)
  {
    G4int nOfStr=strings->size();
#ifdef debug
    G4cout<<"G4QFragmentation::Breeder: SecurityDiQaDiQReduction, #OfStr="<<nOfStr<<G4endl;
#endif
    for (G4int astring=0; astring < nOfStr; astring++)
    {
      G4QString* curString=(*strings)[astring];
      G4QParton* cLeft=curString->GetLeftParton();
      G4QParton* cRight=curString->GetRightParton();
      G4int LT=cLeft->GetType();
      G4int RT=cRight->GetType();
      if(LT==2 && RT==2)
      {
#ifdef debug
        G4int sPDG=cLeft->GetPDGCode();
        G4int nPDG=cRight->GetPDGCode();
        G4cout<<"G4QFragmentation::Breeder:TrySelfRedString,L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
        if( cLeft->ReduceDiQADiQ(cLeft, cRight) ) // DiQ-aDiQ pair was successfully reduced
        {
#ifdef debug
          sPDG=cLeft->GetPDGCode();
          nPDG=cRight->GetPDGCode();
          G4cout<<"+G4QFragm::Breed:#"<<astring<<" Reduced, L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
        }
#ifdef debug
        else G4cout<<"--*--G4QFragm::Breed:#"<<astring<<" DQ-aDQ reduction Failed"<<G4endl;
#endif
      } // End of the found DiQ/aDiQ pair
    }
  }
  //
  // ------------ At this point the strings are fragmenting to hadrons in LS -------------
  //
  G4QHadronVector* theResult = new G4QHadronVector;
  G4int nOfStr=strings->size();
#ifdef debug
  G4cout<<"G4QFragmentation::Breeder: BeforeFragmentation, #OfStr="<<nOfStr<<G4endl;
#endif
  for (G4int astring=0; astring < nOfStr; astring++)
  {
#ifdef edebug
    G4LorentzVector sum=theNucleus->Get4Momentum();  // Nucleus 4Mom in LS
    G4int rChg=totChg-theNucleus->GetZ();
    G4int rBaN=totBaN-theNucleus->GetA();
    G4int nOfHadr=theResult->size();
    G4cout<<"-EMCLS-G4QFragmentation::Breeder:#ofSt="<<nOfStr<<",#ofHad="<<nOfHadr<<G4endl;
    for(G4int i=astring; i<nOfStr; i++)
    {
      G4LorentzVector strI4M=(*strings)[i]->Get4Momentum();
      sum+=strI4M;
      G4int sChg=(*strings)[i]->GetCharge();
      rChg-=sChg;
      G4int sBaN=(*strings)[i]->GetBaryonNumber();
      rBaN-=sBaN;
      G4cout<<"-EMCLS-G4QF::Breed:S#"<<i<<",4M="<<strI4M<<",C="<<sChg<<",B="<<sBaN<<G4endl;
    }
    for(G4int i=0; i<nOfHadr; i++)
    {
      G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
      sum+=hI4M;
      G4int hChg=(*theResult)[i]->GetCharge();
      rChg-=hChg;
      G4int hBaN=(*theResult)[i]->GetBaryonNumber();
      rBaN-=hBaN;
      G4cout<<"-EMCLS-G4QFr::Breed: H#"<<i<<",4M="<<hI4M<<",C="<<hChg<<",B="<<hBaN<<G4endl;
    }
    G4cout<<"....-EMCLS-G4QFrag::Br:r4M="<<sum-totLS4M<<",rC="<<rChg<<",rB="<<rBaN<<G4endl;
#endif
    G4QString* curString=(*strings)[astring];
    if(!curString->GetDirection()) continue;  // Historic for the dead strings: DoesNotWork
#ifdef edebug
    G4int curStrChg = curString->GetCharge();
    G4int curStrBaN = curString->GetBaryonNumber();
#endif
    G4LorentzVector curString4M = curString->Get4Momentum();
#ifdef debug
    G4cout<<"====>G4QFragmentation::Breeder: String#"<<astring<<",s4M/m="<<curString4M
          <<curString4M.m()<<", LPart="<<curString->GetLeftParton()->GetPDGCode()
          <<", RPart="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
#endif
    G4QHadronVector* theHadrons = 0;           // Prototype of theStringFragmentationOUTPUT
    theHadrons=curString->FragmentString(true);// !! Fragmenting the String !!
    if (!theHadrons)                           // The string can not be fragmented
    {
      // First try to correct the diQ-antiDiQ strings, converting them to Q-antiQ
      G4QParton* cLeft=curString->GetLeftParton();
      G4QParton* cRight=curString->GetRightParton();
      G4int sPDG=cLeft->GetPDGCode();
      G4int nPDG=cRight->GetPDGCode();
      G4int LT=cLeft->GetType();
      G4int RT=cRight->GetType();
      G4int LS=LT+RT;
      if(LT==2 && RT==2)
      {
#ifdef debug
        G4cout<<"G4QFragmentation::Breeder:TryReduceString, L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
        if( cLeft->ReduceDiQADiQ(cLeft, cRight) ) // DiQ-aDiQ pair was successfully reduced
        {
          LT=1;
          RT=1;
          LS=2;
          sPDG=cLeft->GetPDGCode();
          nPDG=cRight->GetPDGCode();
#ifdef debug
          G4cout<<"G4QFragmentation::Breeder:AfterReduction,L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
          theHadrons=curString->FragmentString(true);//!! Try to fragment the new String !!
          cLeft=curString->GetLeftParton();
          cRight=curString->GetRightParton();
#ifdef debug
          G4cout<<"G4QFrag::Breed:L="<<cLeft->Get4Momentum()<<",R="<<cRight->Get4Momentum()
                <<G4endl;
#endif
        }
#ifdef debug
        else G4cout<<"^G4QFragmentation::Breeder: DQ-aDQ reduction to Q-aQ Failed"<<G4endl;
#endif
      } // End of the SelfReduction
#ifdef debug
      G4cout<<"G4QFrag::Breed:AfterRedAttempt, theH="<<theHadrons<<", L4M="
            <<cLeft->Get4Momentum()<<", R4M="<<cRight->Get4Momentum()<<G4endl;
#endif
      unsigned next=astring+1;                 // The next string position
      if (!theHadrons)                         // The string can not be fragmented
      {
        G4int fusionDONE=0; // StringFusion didn't happen (1=Fuse L+L/R+R, -1=Fuse L+R/R+L)
        if(next < strings->size())             // TheString isn't theLastString can fuse
        {
          G4int fustr=0;                       // The found partner index (never can be 0)
          G4int swap=0;                        // working interger for swapping parton PDG
          G4double Vmin=DBL_MAX;               // Prototype of the found Velocity Distance 
          G4int dPDG=nPDG;
          G4int qPDG=sPDG;
          if(dPDG<-99 || (dPDG>0&&dPDG<7) || qPDG>99 || (qPDG<0 && qPDG>-7))
          {
            swap=qPDG;
            qPDG=dPDG;
            dPDG=swap;
          }
          if(dPDG>99) dPDG/=100;
          if(qPDG<-99) qPDG=-(-qPDG)/100;
#ifdef debug
          G4cout<<"G4QFrag::Breed:TryFuseStringS, q="<<qPDG<<", a="<<dPDG<<", n="<<next
                <<G4endl;
#endif
          G4ThreeVector curV=curString4M.vect()/curString4M.e();
          G4int reduce=0;                      // a#of reduced Q-aQ pairs
          G4int restr=0;                       // To use beyon the LOOP for printing
          G4int MPS=0;                         // PLS for the selected string
          for (restr=next; restr < nOfStr; restr++)
          {
            G4QString* reString=(*strings)[restr];
            G4QParton* Left=reString->GetLeftParton();
            G4QParton* Right=reString->GetRightParton();
            G4int uPDG=Left->GetPDGCode();
            G4int mPDG=Right->GetPDGCode();
            G4int PLT =Left->GetType();
            G4int PRT =Right->GetType();
            G4int aPDG=mPDG;
            G4int rPDG=uPDG;
            if(aPDG<-99 || (aPDG>0 && aPDG<7) || rPDG>99 || (rPDG<0 && rPDG>-7))
            {
              swap=rPDG;
              rPDG=aPDG;
              aPDG=swap;
            }
            if(aPDG > 99) aPDG/=100;
            if(rPDG <-99) rPDG=-(-rPDG)/100;
            // Try to reduce two DQ-aDQ strings
#ifdef debug
            G4cout<<"G4QFragm::Breed: TryReduce #"<<restr<<",q="<<rPDG<<",a="<<aPDG<<G4endl;
#endif
            if(LT==2 && RT==2 && PLT==2 && PRT==2)    // Have a chance for the reduction
            {
              G4int cQ1=(-qPDG)/10;
              G4int cQ2=(-qPDG)%10;
              G4int cA1=dPDG/10;
              G4int cA2=dPDG%10;
              G4int pQ1=(-rPDG)/10;
              G4int pQ2=(-rPDG)%10;
              G4int pA1=aPDG/10;
              G4int pA2=aPDG%10;
#ifdef debug
		  G4cout<<"G4QFragment::Breeder: cQ="<<cQ1<<","<<cQ2<<", cA="<<cA1<<","<<cA2
                    <<", pQ="<<pQ1<<","<<pQ2<<", pA="<<pA1<<","<<pA2<<G4endl;
#endif
              G4bool iQA = (cA1==pQ1 || cA1==pQ2 || cA2==pQ1 || cA2==pQ2);
              G4bool iAQ = (cQ1==pA1 || cQ1==pA2 || cQ2==pA1 || cQ2==pA2);
              if(iQA) reduce++;
              if(iAQ) reduce++;
              if  (reduce==2)                  // Two quark pairs can be reduced
              {
                if(sPDG>0 && uPDG<0)           // LL/RR Reduction
                {
                  std::pair<G4int,G4int> resLL=ReducePair(sPDG/100, (-uPDG)/100);
                  G4int newCL=resLL.first;
                  G4int newPL=resLL.second;
                  if(!newCL || !newPL)
                  {
                    G4cerr<<"*G4QFragmentation::Breeder:CL="<<newCL<<",PL="<<newPL<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-LL");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair((-nPDG)/100, mPDG/100);
                  G4int newCR=resRR.first;
                  G4int newPR=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QFragmentation::Breeder:CR="<<newCR<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-RR");
                  }
                  cLeft->SetPDGCode(newCL);    // Reset the left quark of curString
                  cRight->SetPDGCode(-newCR);  // Reset the right quark of curString
                  Left->SetPDGCode(-newPL);    // Reset the left quark of reString
                  Right->SetPDGCode(newPR);    // Reset the right quark of reString
                  break;                       // Break out of the reString internal LOOP
                }
                else if(sPDG>0)                // LR Reduction
                {
                  std::pair<G4int,G4int> resLL=ReducePair(sPDG/100, (-mPDG)/100);
                  G4int newCL=resLL.first;
                  G4int newPR=resLL.second;
                  if(!newCL || !newPR)
                  {
                    G4cerr<<"*G4QFragmentation::Breeder:CL="<<newCL<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-LR");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair((-nPDG)/100, uPDG/100);
                  G4int newCR=resRR.first;
                  G4int newPL=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QFragmentation::Breeder:CR="<<newCR<<",PL="<<newPL<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-LR");
                  }
                  cLeft->SetPDGCode(newCL);    // Reset the left quark of curString
                  cRight->SetPDGCode(-newCR);  // Reset the right quark of curString
                  Left->SetPDGCode(newPL);     // Reset the left quark of reString
                  Right->SetPDGCode(-newPR);   // Reset the right quark of reString
                  break;                       // Break out of the reString internal LOOP
                }
                else                           // RL Reduction
                {
                  std::pair<G4int,G4int> resLL=ReducePair((-sPDG)/100, mPDG/100);
                  G4int newCL=resLL.first;
                  G4int newPR=resLL.second;
                  if(!newCL || !newPR)
                  {
                    G4cerr<<"*G4QFragmentation::Breeder:CL="<<newCL<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-RL");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair(nPDG/100, (-uPDG)/100);
                  G4int newCR=resRR.first;
                  G4int newPL=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QFragmentation::Breeder:CR="<<newCR<<",PL="<<newPL<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-RL");
                  }
                  cLeft->SetPDGCode(-newCL);   // Reset the left quark of curString
                  cRight->SetPDGCode(newCR);   // Reset the right quark of curString
                  Left->SetPDGCode(-newPL);    // Reset the left quark of reString
                  Right->SetPDGCode(newPR);    // Reset the right quark of reString
                  break;                       // Break out of the reString internal LOOP
                }
              } // End of IF(possible reduction)
            }
            // Check StringsCanBeFused: all QaQ+QaQ(22), single QaQ+QDiQ/aQaDtQ(23/32),
            //                          double QaQ+DiQaDiQ(24/42), QDiQ+aDiQaQ(34/43)
#ifdef debug
            G4cout<<"G4QFragm::Breed:TryFuse/w #"<<restr<<",q="<<rPDG<<",a="<<aPDG<<G4endl;
#endif
            G4int PLS=PLT+PRT;
            if( (LS==2 && PLS==2) ||                           // QaQ+QaQ always to DiQaDiQ
                ( ( (LS==2 && PLS==3) || (LS==3 && PLS==2) ) &&// QaQ w QDiQ/aQaDiQ(single)
                  ( (aPDG> 7 && (-dPDG==aPDG/10   || -dPDG==aPDG%10) )   || // cAQ w DQ
                    (dPDG> 7 && (-aPDG==dPDG/10   || -aPDG==dPDG%10) )   || // AQ w cDQ
                    (rPDG<-7 && (qPDG==(-rPDG)/10 || qPDG==(-rPDG)%10) ) || // cQ w ADQ
                    (qPDG<-7 && (rPDG==(-qPDG)/10 || rPDG==(-qPDG)%10) ) || // Q w cADQ
                    (aPDG< 0 && -aPDG==qPDG) || (dPDG< 0 && -dPDG==rPDG)    // aQ w Q 
                  )
                )                 ||                           // -----------------------
                ( ( LS==2 && PLS==4                          &&// QaQ w DiQaDiQ (double)
                    (aPDG> 7 && (-dPDG == aPDG/10 || -dPDG == aPDG%10) ) &&
                    (rPDG<-7 && (qPDG==(-rPDG)/10 || qPDG==(-rPDG)%10) )
                  )       ||
                  ( LS==4 && PLS==2                          &&// DiQaDiQ w QaQ (double)
                    (dPDG> 7 && (-aPDG == dPDG/10 || -aPDG == dPDG%10) ) &&
                    (qPDG<-7 && (rPDG==(-qPDG)/10 || rPDG==(-qPDG)%10) )
                  )
                )                 ||                           // -----------------------
                ( LS==3 && PLS==3                            &&// QDiQ w aDiQaQ (double)
                  ( (aPDG> 7 && (-dPDG == aPDG/10 || -dPDG == aPDG%10) &&
                     qPDG<-7 && (rPDG==(-qPDG)/10 || rPDG==(-qPDG)%10)
                    )       ||
                    (dPDG> 7 && (-aPDG == dPDG/10 || -aPDG == dPDG%10) &&
                     rPDG<-7 && (dPDG==(-rPDG)/10 || dPDG==(-rPDG)%10) 
                    )
                  )
                )
              )
            {
              G4LorentzVector reString4M = reString->Get4Momentum();
              G4ThreeVector reV = reString4M.vect()/reString4M.e();
              G4double dV=(curV-reV).mag2();   // Squared difference between relVelocities
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder: StringCand#"<<restr<<", q="<<rPDG<<", a="
                    <<aPDG<<", L="<<uPDG<<", R="<<mPDG<<",dV="<<dV<<G4endl;
#endif
              if(dV < Vmin)
              {
                Vmin=dV;
                fustr=restr;
                MPS=PLS;
              }
            }
          }
          if(reduce==2)                        // String mutual reduction happened
          {
#ifdef debug
            G4cout<<"-G4QFragmentation::Breeder:Reduced #"<<astring<<" & #"<<restr<<G4endl;
#endif
            astring--;                         // String was QCreduced using another String
            continue;                          // Repeat fragmentation of the same string
          }
          if(fustr)                            // The partner was found -> fuse strings
          {
#ifdef debug
            G4cout<<"G4QFragm::Breeder: StPartner#"<<fustr<<",LT="<<LT<<",RT="<<RT<<G4endl;
#endif
            G4QString* fuString=(*strings)[fustr];
            G4QParton* Left=fuString->GetLeftParton();
            G4QParton* Right=fuString->GetRightParton();
            G4int uPDG=Left->GetPDGCode();
            G4int mPDG=Right->GetPDGCode();
            G4int rPDG=uPDG;
            G4int aPDG=mPDG;
            if(aPDG<-99 || (aPDG>0 && aPDG<7) || rPDG>99 || (rPDG<0 && rPDG>-7))
            {
              swap=rPDG;
              rPDG=aPDG;
              aPDG=swap;
            }
            if(aPDG > 99) aPDG/=100;
            if(rPDG <-99) rPDG=-(-rPDG)/100;
            // Check that the strings can fuse (no anti-diquarks assumed)
            G4LorentzVector resL4M(0.,0.,0.,0.);
            G4LorentzVector resR4M(0.,0.,0.,0.);
            G4LorentzVector L4M=cLeft->Get4Momentum();
            G4LorentzVector R4M=cRight->Get4Momentum();
#ifdef debug
            G4cout<<"G4QFragmentation::Breeder:BeforeFuDir,sL="<<sPDG<<",nR="<<nPDG<<",uL="
                  <<uPDG<<",mR="<<mPDG<<",L4M="<<L4M<<",R4M="<<R4M<<G4endl;
#endif
            G4LorentzVector PL4M=Left->Get4Momentum();
            G4LorentzVector PR4M=Right->Get4Momentum();
            fusionDONE=AnnihilationOrder(LS, MPS, uPDG, mPDG, sPDG, nPDG);
            if     (fusionDONE>0)
            {
              if(fusionDONE>1)                             // Anihilation algorithm
              {
                if     ( (uPDG<0 || nPDG<0) && -uPDG==nPDG ) Left->SetPDGCode(sPDG);
                else if( (mPDG<0 || sPDG<0) && -mPDG==sPDG ) Right->SetPDGCode(nPDG);
                else if( (uPDG<0 || sPDG<0) && -uPDG==sPDG ) Left->SetPDGCode(nPDG);
                else if( (mPDG<0 || nPDG<0) && -mPDG==nPDG ) Right->SetPDGCode(sPDG);
              }
              {
                Left->SetPDGCode(SumPartonPDG(uPDG, sPDG));
                Right->SetPDGCode(SumPartonPDG(mPDG, nPDG));
              }
              Left->Set4Momentum(L4M+PL4M);
              Right->Set4Momentum(R4M+PR4M);
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder:LL/RR s4M="<<fuString->Get4Momentum()
                    <<",S="<<L4M+PL4M+R4M+PR4M<<", L="<<Left->Get4Momentum()<<", R="
                    <<Right->Get4Momentum()<<G4endl;
#endif
            }
            else if(fusionDONE<0)
            {
              Left->SetPDGCode(SumPartonPDG(uPDG, nPDG));
              Left->Set4Momentum(L4M+PR4M);
              Right->SetPDGCode(SumPartonPDG(mPDG, sPDG));
              Right->Set4Momentum(R4M+PL4M);
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder:LR/RL s4M="<<fuString->Get4Momentum()
                    <<",S="<<L4M+PL4M+R4M+PR4M<<", L="<<Left->Get4Momentum()<<", R="
                    <<Right->Get4Momentum()<<G4endl;
#endif
            }
#ifdef debug
            else G4cout<<"-Warning-G4QFragmentation::Breeder: WrongStringFusion"<<G4endl;
#endif
#ifdef edebug
            G4cout<<"#EMC#G4QFragmentation::Breeder:StringFused,F="<<fusionDONE<<",L="<<L4M
                  <<",R="<<R4M<<",pL="<<PL4M<<",pR="<<PR4M<<",nL="<<Left->Get4Momentum()
                  <<",nR="<<Right->Get4Momentum()<<",S="<<fuString->Get4Momentum()<<G4endl;
#endif
            if(fusionDONE)
            {
#ifdef debug
              G4cout<<"###G4QFragmentation::Breeder: Str#"<<astring<<" fused/w Str#"<<fustr
                    <<"->"<<fuString->Get4Momentum()<<fuString->Get4Momentum().m()<<G4endl;
#endif
              continue;                          // @@ killing of the curString is excluded
            }
          }
          else
          {
            G4cerr<<"**G4QFragmentation::Breeder:*NoPart*M="<<curString->Get4Momentum().m()
                  <<", F="<<fusionDONE<<", LPDG="<<curString->GetLeftParton()->GetPDGCode()
                  <<", RPDG="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
            // @@ Temporary exception for the future solution
            //G4Exception("G4QFragmentation::Breeder:","72",FatalException,"NoStrgsFused");
          }
        }
        if(!fusionDONE)                          // Fusion of theString failed, try hadrons
        {
          G4int nHadr=theResult->size();         // #of collected resulting hadrons upToNow
#ifdef debug
          G4cout<<"+++4QFragmentation::Breeder:*TryHadr* M="<<curString->Get4Momentum().m()
                <<", nH="<<nHadr<<", LPDG="<<curString->GetLeftParton()->GetPDGCode()
                <<", RPDG="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
#endif
          // The only solution now is to try fusion with the secondary hadron
          while( (nHadr=theResult->size()) && !theHadrons)
          {
#ifdef edebug
            for(G4int i=0; i<nHadr; i++)
            {
              G4LorentzVector h4M=(*theResult)[i]->Get4Momentum();
              G4int hPDG=(*theResult)[i]->GetPDGCode();
              G4QContent hQC=(*theResult)[i]->GetQC();
              G4cout<<"-EMC-G4QFrag::Breed:H#"<<i<<",4M="<<h4M<<hQC<<",PDG="<<hPDG<<G4endl;
            }
#endif
            G4int fusDONE=0;                      // String+Hadron fusion didn't happen
            G4int fuhad=-1;                       // The found hadron index
            G4int newPDG=0;                       // PDG ofTheParton afterMerging with Hadr
            G4int secPDG=0;                       // Second PDG oParton afterMerging w/Hadr
            G4double maM2=-DBL_MAX;               // Prototype of the max ResultingStringM2
            G4LorentzVector selH4M(0.,0.,0.,0.);  // 4-mom of the selected hadron
            G4QHadron* selHP=0;                   // Pointer to the used hadron for erasing
            G4QString* cString=(*strings)[astring];// Must be the last string by definition
            G4LorentzVector cString4M = cString->Get4Momentum();
            G4QParton* cLeft=cString->GetLeftParton();
            G4QParton* cRight=cString->GetRightParton();
            G4int sumT=cLeft->GetType()+cRight->GetType();
            G4int sPDG=cLeft->GetPDGCode();
            G4int nPDG=cRight->GetPDGCode();
            G4int cMax=0;                         // Both or only one cuark can merge
            for (G4int reh=0; reh < nHadr; reh++)
            {
              G4QHadron* curHadr=(*theResult)[reh];
              G4QContent curQC=curHadr->GetQC();  // Quark content of the hadron
              G4int partPDG1=curQC.AddParton(sPDG);
              G4int partPDG2=curQC.AddParton(nPDG);
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder: Hadron # "<<reh<<", QC="<<curQC
                      <<", P1="<<partPDG1<<", P2="<<partPDG2<<G4endl;
#endif
              if(partPDG1 || partPDG2)            // Hadron can merge at least w/one parton
              {
                G4int cCur=1;
                if(sumT>3 && partPDG1 && partPDG2) cCur=2;
                G4LorentzVector curHadr4M = curHadr->Get4Momentum();
                G4double M2=(cString4M+curHadr4M).m2();// SquaredMass of theResultingString
#ifdef debug
                G4cout<<"G4QFragmentation::Breeder:*IN*Hadron#"<<reh<<",M2="<<M2<<G4endl;
#endif
                if( (sumT<4 || cCur>=cMax) && M2 > maM2) // DeepAnnihilation for DiQ-aDiQ 
                {
                  maM2=M2;
                  fuhad=reh;
                  if(partPDG1)
                  {
                    fusDONE= 1;
                    newPDG=partPDG1;
                    if(partPDG2)
                    {
                      secPDG=partPDG2;
                      cMax=cCur;
                    }
                  }
                  else
                  {
                    fusDONE=-1;
                    newPDG=partPDG2;
                  }
#ifdef debug
                  G4cout<<"G4QFrag::Br:*Selected*,P1="<<partPDG1<<",P2="<<partPDG2<<G4endl;
#endif
                  selH4M=curHadr4M;
                  selHP=curHadr;
                } // End of IF(update selection)
              } // End of IF(HadronCanMergeWithTheString)
            } // End of the LOOP over Hadrons
#ifdef debug
            G4cout<<"G4QFragmentation::Breeder: fuh="<<fuhad<<",fus="<<fusDONE<<G4endl;
#endif
            if(fuhad>-1)                          // The partner was found -> fuse H&S
            {
              if     (fusDONE>0)                  // Fuse hadron with the left parton
              {
                cLeft->SetPDGCode(newPDG);
                G4LorentzVector newLeft=cLeft->Get4Momentum()+selH4M;
                cLeft->Set4Momentum(newLeft);
                if(secPDG && cMax>1)
                {
#ifdef debug
                  G4cout<<"G4QFragm::Br:TryReduce,nPDG="<<newPDG<<",sPDG="<<secPDG<<G4endl;
#endif
                  cLeft->ReduceDiQADiQ(cLeft, cRight);
                }
#ifdef debug
                G4cout<<"G4QFragmentation::Breeder: Left, s4M="<<curString->Get4Momentum()
                      <<", L4M="<<newLeft<<", R4M="<<cRight->Get4Momentum()<<", h4M="
                      <<selH4M<<", newL="<<newPDG<<", oldR="<<cRight->GetPDGCode()<<G4endl;
#endif
              }
              else if(fusDONE<0)                  // Fuse hadron with the right parton
              {
                cRight->SetPDGCode(newPDG);
                G4LorentzVector newRight=cRight->Get4Momentum()+selH4M;
                cRight->Set4Momentum(newRight);
#ifdef debug
                G4cout<<"G4QFragmentation::Breeder: Right, s4M="<<curString->Get4Momentum()
                      <<", L4M="<<cLeft->Get4Momentum()<<", R4M="<<newRight<<", h4M="
                      <<selH4M<<", newR="<<newPDG<<", oldL="<<cLeft->GetPDGCode()<<G4endl;
#endif
              }
#ifdef debug
              else G4cout<<"-G4QFragmentation::Breeder: Wrong String+HadronFusion"<<G4endl;
#endif
#ifdef debug
              if(fusDONE) G4cout<<"####G4QFragmentation::Breeder: String #"<<astring
                                <<" is fused with Hadron #"<<fuhad
                                <<", new4Mom="<<curString->Get4Momentum()
                                <<", M2="<<curString->Get4Momentum().m2()
                                <<", QC="<<curString->GetQC()<<G4endl;
#endif
            }
            else
            {
              G4cerr<<"**G4QFragmentation::Breeder:*NoH*,M="<<curString->Get4Momentum().m()
                    <<", LPDG="<<curString->GetLeftParton()->GetPDGCode()
                    <<", RPDG="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
              // @@ Temporary exception for the future solution
              //G4Exception("G4QFragmentation::Breeder:","72",FatalException,"SHNotFused");
              break;                           // Breake the While LOOP
            } // End of the namespace where both Fusion and reduction have failed
            // The fused hadron must be excluded from theResult
#ifdef debug
            G4cout<<"G4QFragmentation::Breeder: before HR, nH="<<theResult->size()<<G4endl;
            G4int icon=0;                              // Loop counter
#endif
            G4QHadronVector::iterator ih;
            G4bool found=false;
            for(ih = theResult->begin(); ih != theResult->end(); ih++)
            {
#ifdef debug
              G4cout<<"G4QFrag::Breeder:#"<<icon<<", i="<<(*ih)<<", sH="<<selHP<<G4endl;
              icon++;
#endif
              if((*ih)==selHP)
              {
#ifdef debug
                G4cout<<"G4QFragm::Breed: *HadronFound*,PDG="<<selHP->GetPDGCode()<<G4endl;
#endif
#ifdef edebug
                G4LorentzVector p4M=selHP->Get4Momentum();
                G4int           Chg=selHP->GetCharge();
                G4int           BaN=selHP->GetBaryonNumber();
                curString4M+=p4M;
                curStrChg+=Chg;
                curStrBaN+=BaN;
                G4cout<<"-EMC->>>>G4QFragmentation::Breeder: S+=H, 4M="<<curString4M<<",M="
                      <<curString4M.m()<<", Charge="<<curStrChg<<", B="<<curStrBaN<<G4endl;
#endif
                delete selHP;                          // delete the Hadron
                theResult->erase(ih);                  // erase the Hadron from theResult
                found=true;
                break;                                 // beak the LOOP over hadrons
              }
            } // End of the LOOP over hadrons
#ifdef debug
            if(!found) G4cout<<"*G4QFragmentation::Breeder:nH="<<theResult->size()<<G4endl;
#endif
            // New attempt of the string decay
            theHadrons=curString->FragmentString(true);//! Try to fragment the new String !
#ifdef debug
            G4cout<<"G4QFrag::Breeder: tH="<<theHadrons<<",nH="<<theResult->size()<<G4endl;
#endif
          } // End of the while LOOP over the fusion with hadrons
#ifdef debug
          G4cout<<"*G4QFragmentation::Breeder: *CanTryToDecay?* nH="<<theHadrons<<", next="
                <<next<<" =? nS="<<strings->size()<<", nR="<<theResult->size()<<G4endl;
#endif
          if(!theHadrons && next == strings->size() && !(theResult->size()))// TryToDecay
          {
            G4QContent miQC=curString->GetQC();        // QContent of the Lightest Hadron
            G4int miPDG=miQC.GetSPDGCode();            // PDG of the Lightest Hadron
            if(miPDG == 10)                        // ==> Decay Hadron-Chipolino
            {
              G4QChipolino QCh(miQC);              // define theTotalResidual as aChipolino
              G4QPDGCode   h1QPDG=QCh.GetQPDG1();  // QPDG of the first hadron
              G4double     h1M   =h1QPDG.GetMass();// Mass of the first hadron
              G4QPDGCode   h2QPDG=QCh.GetQPDG2();  // QPDG of the second hadron 
              G4double     h2M   =h2QPDG.GetMass();// Mass of the second hadron
              G4double     ttM   =curString4M.m(); // Real Mass of the Chipolino
              if(h1M+h2M<ttM+eps)                  // Two particles decay of Chipolino
              {
                G4LorentzVector h14M(0.,0.,0.,h1M);
                G4LorentzVector h24M(0.,0.,0.,h2M);
                if(std::fabs(ttM-h1M-h2M)<=eps)
                {
                  G4double part1=h1M/(h1M+h2M);
                  h14M=part1*curString4M;
                  h24M=curString4M-h14M;
                }
                else
                {
                  if(!G4QHadron(curString4M).DecayIn2(h14M,h24M))
                  {
                    G4cerr<<"***G4QFragmentation::Breeder: tM="<<ttM<<"->h1="<<h1QPDG<<"("
                          <<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")="<<h1M+h2M<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ChiDec");
                  }
                }
                G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
                theResult->push_back(h1H);         // (delete equivalent)  
#ifdef debug
                G4LorentzVector f4M=h1H->Get4Momentum();
                G4int           fPD=h1H->GetPDGCode();
                G4int           fCg=h1H->GetCharge();
                G4int           fBN=h1H->GetBaryonNumber();
                G4cout<<"-EMC->>>G4QFragmentation::Breeder:St=Hadr ChiPro1 is filled, f4M="
                      <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
                G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
                theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
                G4LorentzVector s4M=h2H->Get4Momentum();
                G4int           sPD=h2H->GetPDGCode();
                G4int           sCg=h2H->GetCharge();
                G4int           sBN=h2H->GetBaryonNumber();
                G4cout<<"-EMC->>>G4QFragmentation::Breeder:St=Hadr ChiPro2 is filled, s4M="
                      <<s4M<<", sPDG="<<sPD<<", sCg="<<sCg<<", sBN="<<sBN<<G4endl;
#endif
#ifdef edebug
                G4cout<<"-EMC-..Chi..G4QFragmentation::Breeder:DecCHECK,Ch4M="<<curString4M
                      <<", ChQC="<<miQC<<", d4M="<<curString4M-h14M-h24M<<G4endl;
#endif
                break;                               // Go out of the main StringDecay LOOP
              }
            }
            else if(miPDG)                                 // Decay Hadron as a Resonans
            {
              if     (miPDG>0 &&   miPDG%10 < 3) miPDG+=2; // Convert to Resonans
              else if(miPDG<0 && (-miPDG)%10< 3) miPDG-=2; // Convert to antiResonans
              G4Quasmon Quasm;
              G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
              G4QHadronVector* tmpQHadVec=Quasm.DecayQHadron(sHad); // It deletes sHad
#ifdef pdebug
              G4cout<<"G4QFragmentation::Breeder: DecLast,nH="<<tmpQHadVec->size()<<G4endl;
#endif
              for(unsigned aH=0; aH < tmpQHadVec->size(); aH++)
              {
                theResult->push_back((*tmpQHadVec)[aH]);// TheDecProductOfHadronDecIsFilled
#ifdef debug
                G4QHadron*   prodH =(*tmpQHadVec)[aH];
                G4LorentzVector p4M=prodH->Get4Momentum();
                G4int           PDG=prodH->GetPDGCode();
                G4int           Chg=prodH->GetCharge();
                G4int           BaN=prodH->GetBaryonNumber();
                G4cout<<"-EMC->>>G4QFragmentation::Breeder:Str=Had,pH#"<<aH<<" filled, 4M="
                    <<p4M<<", PDG="<<PDG<<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
              }
              tmpQHadVec->clear();
              delete tmpQHadVec;  // Who calls DecayQHadron is responsible for clear&delete
              break;                               // Go out of the main String Decay LOOP
            }
          } // End of the DecayOfTheLast
        } // End of IF(String-Hadron fusion)
      } // End of IF(NO_Hadrons) for String-String and String-Hadron fusion
      // The last hope is to CORREC the string, using other strings (ForwardInLOOP)
#ifdef debug
      G4cout<<"G4QFragmentation::Breeder: theH="<<theHadrons<<"?=0, next="<<next<<G4endl;
#endif
      if(!theHadrons && next < strings->size())       // ForwardInLOOP strings exist
      {
        // @@ string can be not convertable to one hadron (2,0.0,0,2,0) - To be improved
        G4QContent miQC=curString->GetQC(); // QContent of the Lightest Hadron
        G4int miPDG=miQC.GetSPDGCode();// PDG of the Lightest Hadron
        G4double miM=0.;               // Prototype of the Mass of the Cur LightestHadron
        if(miPDG!=10) miM=G4QPDGCode(miPDG).GetMass(); // Mass of theCurLightestOneHadron
        else
        {
          G4QChipolino QCh(miQC);      // define the TotalString as a Chipolino
          miM=QCh.GetQPDG1().GetMass()+QCh.GetQPDG2().GetMass();//MinMass of theStringChipo
        }
        G4double cM=curString4M.m();   // Actual mass of the Cur String
        if(std::fabs(cM-miM) < eps)    // Convert to hadron(2 hadrons) without calculations
        {
          if(miPDG!=10)
          {
            G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
            theResult->push_back(sHad);// Fill the curString as a hadron
#ifdef debug
            G4cout<<">>>G4QFragmentation::Breeder:S->H="<<miPDG<<curString4M<<G4endl;
#endif
          }
          else
          {
            G4QChipolino QCh(miQC);               // define the TotalResidual as aChipolino
            G4QPDGCode   h1QPDG=QCh.GetQPDG1();   // QPDG of the first hadron
            G4double     h1M   =h1QPDG.GetMass(); // Mass of the first hadron
            G4QPDGCode   h2QPDG=QCh.GetQPDG2();   // QPDG of the second hadron 
            G4double     h2M   =h2QPDG.GetMass(); // Mass of the second hadron
            G4double     pt1   =h1M/(h1M+h2M);    // 4-mom part of the first hadron
            G4LorentzVector h14M=pt1*curString4M; // 4-mom of the first hadron
            G4LorentzVector h24M=curString4M-h14M;// 4-mom of the second hadron
            G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
            theResult->push_back(h1H);            // (delete equivalent)  
#ifdef debug
            G4LorentzVector f4M=h1H->Get4Momentum();
            G4int           fPD=h1H->GetPDGCode();
            G4int           fCg=h1H->GetCharge();
            G4int           fBN=h1H->GetBaryonNumber();
            G4cout<<"-EMC->>>G4QFragmentation::Breeder: Str=2HadrAR Prod-F is filled, f4M="
                  <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
            G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
            theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
            G4LorentzVector s4M=h2H->Get4Momentum();
            G4int           sPD=h2H->GetPDGCode();
            G4int           sCg=h2H->GetCharge();
            G4int           sBN=h2H->GetBaryonNumber();
            G4cout<<"-EMC->>>G4QFragmentation::Breeder: Str=2HadrAR Prod-S is filled, s4M="
                  <<s4M<<", sPDG="<<sPD<<", sCg="<<sCg<<", sBN="<<sBN<<G4endl;
#endif
          }
          continue;                    // Continue the LOOP over the curStrings
        }
        else                           // Try to recover (+/-) to min Mass
        {
          G4ThreeVector cP=curString4M.vect(); // Momentum of the curString
          G4double      cE=curString4M.e();    // Energy of the curString
          G4ThreeVector curV=cP/cE;    // curRelativeVelocity
          G4double miM2=miM*miM;
          G4double cM2=cM*cM;
          G4int restr=0;               // To use beyon the LOOP for printing
          G4int fustr=0;               // Selected String-partner (0 = NotFound)
          G4double selX=0.;            // Selected value of x
          G4double maD=-DBL_MAX;       // Maximum Free Mass
          G4double Vmin=DBL_MAX;       // min Velocity Distance
          G4LorentzVector s4M(0.,0.,0.,0.); // Selected 4-mom of the hadron
          for (restr=next; restr < nOfStr; restr++)
          {
            G4QString* pString=(*strings)[restr];
            G4LorentzVector p4M=pString->Get4Momentum();
            G4ThreeVector pP=p4M.vect(); // Momentum of the partnerString
            G4double      pE=p4M.e();    // Energy of the partnerString
            G4ThreeVector pV=pP/pE;    // curRelativeVelocity
            G4double D2=cE*pE-cP.dot(pP); 
            G4double pM2=p4M.m2();
            G4double dM4=pM2*(miM2-cM2);
            G4double x=(std::sqrt(D2*D2+dM4)-D2)/pM2; // what we should get from p
            if(x > 0)                  // We are getting x part of p4M
            {
              G4int pPDG=pString->GetQC().GetSPDGCode(); // PDG ofTheLightestHadron
              G4double pM=G4QPDGCode(pPDG).GetMass(); // Mass of the LightestHadron
              G4double rM=std::sqrt(pM2); // Real mass of the string-partner
              G4double delta=(1.-x)*rM-pM;// @@ Minimum CM disterbance measure
              if(delta > maD)
              {
                maD=delta;
                fustr=restr;
                selX=x;
                s4M=p4M;
              }
            }
            else                       // We are adding to p4M, so use RelVelocity
            {
              G4double dV=(curV-pV).mag2();// SquaredDifferenceBetweenRelVelocities
              if(dV < Vmin)
              {
                Vmin=dV;
                fustr=restr;
                selX=x;
                s4M=p4M;
              }
            }
          } // End of the LOOP over string-partners for Correction
#ifdef edebug
          G4LorentzVector sum4M=s4M+curString4M;
#endif
          G4QString* pString=(*strings)[fustr];
          curString4M+=selX*s4M;
          if(std::abs(miPDG)%10 > 2)                  // Decay Hadron-Resonans
          {
            G4Quasmon Quasm;
            G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
            G4QHadronVector* tmpQHadVec=Quasm.DecayQHadron(sHad); // It deletes sHad
#ifdef pdebug
            G4cout<<"G4QFragmentation::Breeder: DecStrHad,nH="<<tmpQHadVec->size()<<G4endl;
#endif
            for(unsigned aH=0; aH < tmpQHadVec->size(); aH++)
            {
              theResult->push_back((*tmpQHadVec)[aH]);// TheDecayProductofTheHadronIsFilled
#ifdef debug
              G4QHadron*   prodH =(*tmpQHadVec)[aH];
              G4LorentzVector p4M=prodH->Get4Momentum();
              G4int           PDG=prodH->GetPDGCode();
              G4int           Chg=prodH->GetCharge();
              G4int           BaN=prodH->GetBaryonNumber();
              G4cout<<"-EMC->>>G4QFragmentation::Breeder:Str=Hadr PrH#"<<aH<<" filled, 4M="
                    <<p4M<<", PDG="<<PDG<<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
            }
            tmpQHadVec->clear();
            delete tmpQHadVec;  // Who calls DecayQHadron is responsible for clear & delete
          }
          else if(miPDG == 10)                   // ==> Decay Hadron-Chipolino
          {
            G4QChipolino QCh(miQC);              // define the TotalResidual as a Chipolino
            G4QPDGCode   h1QPDG=QCh.GetQPDG1();  // QPDG of the first hadron
            G4double     h1M   =h1QPDG.GetMass();// Mass of the first hadron
            G4QPDGCode   h2QPDG=QCh.GetQPDG2();  // QPDG of the second hadron 
            G4double     h2M   =h2QPDG.GetMass();// Mass of the second hadron
            G4double     ttM   =curString4M.m(); // Real Mass of the Chipolino
            if(h1M+h2M<ttM+eps)                  // Two particles decay of Chipolino
            {
              G4LorentzVector h14M(0.,0.,0.,h1M);
              G4LorentzVector h24M(0.,0.,0.,h2M);
              if(std::fabs(ttM-h1M-h2M)<=eps)
              {
                G4double part1=h1M/(h1M+h2M);
                h14M=part1*curString4M;
                h24M=curString4M-h14M;
              }
              else
              {
                if(!G4QHadron(curString4M).DecayIn2(h14M,h24M))
                {
                  G4cerr<<"***G4QFragmentation::Breeder: tM="<<ttM<<"->h1="<<h1QPDG<<"("
                        <<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")="<<h1M+h2M<<G4endl;
                  G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ChipoDec");
                }
              }
              G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
              theResult->push_back(h1H);         // (delete equivalent)  
#ifdef debug
              G4LorentzVector f4M=h1H->Get4Momentum();
              G4int           fPD=h1H->GetPDGCode();
              G4int           fCg=h1H->GetCharge();
              G4int           fBN=h1H->GetBaryonNumber();
              G4cout<<"-EMC->>>G4QFragmentation::Breeder:Str=Hadr Prod-F is filled, f4M="
                    <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
              G4LorentzVector s4M=h2H->Get4Momentum();
              G4int           sPD=h2H->GetPDGCode();
              G4int           sCg=h2H->GetCharge();
              G4int           sBN=h2H->GetBaryonNumber();
              G4cout<<"-EMC->>>G4QFragmentation::Breeder:Str=Hadr Prod-S is filled, s4M="
                    <<s4M<<", sPDG="<<sPD<<", sCg="<<sCg<<", sBN="<<sBN<<G4endl;
#endif
#ifdef edebug
              G4cout<<"-EMC-..Chipo..G4QFragmentation::Breeder:DecCHECK,Ch4M="<<curString4M
                    <<", ChQC="<<miQC<<", d4M="<<curString4M-h14M-h24M<<G4endl;
#endif
            }
            else
            {
              G4cerr<<"***G4QFragm::Breeder: tM="<<ttM<<miQC<<"->h1="<<h1QPDG<<"(" <<h1M
                    <<")+h2="<<h1QPDG<<"("<<h2M<<") = "<<h1M+h2M<<G4endl;
              G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ChipoDecMass");
            }
          }
          else
          {
            G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
            theResult->push_back(sHad);             // The original string-hadron is filled
#ifdef debug
            G4cout<<"-EMC->>>G4QFragmentation::Breeder: Str=Hadr Filled, 4M="<<curString4M
                  <<", PDG="<<miPDG<<G4endl;
#endif
          }
          G4double corF=1-selX;
          G4QParton* Left=pString->GetLeftParton();
          G4QParton* Right=pString->GetRightParton();
          Left->Set4Momentum(corF*Left->Get4Momentum());
          Right->Set4Momentum(corF*Right->Get4Momentum());
#ifdef edebug
          G4cout<<"-EMC-...Cor...G4QFragmentation::Breeder:CorCHECK Sum="<<sum4M
                <<" =? "<<curString4M+pString->Get4Momentum()<<", M="<<miM<<" ?= "
                <<curString4M.m()<<G4endl;
#endif
#ifdef debug
          G4cout<<">>>G4QFragmentation::Breeder:*Corrected* String->Hadr="<<miPDG
                <<curString4M<<" by String #"<<restr<<G4endl;
#endif
          continue;                            // Continue the LOOP over the curStrings
        } // Try to recover String+String correction-algorithm End
      } // End of IF(Can try to correct by String-String)
#ifdef debug
      else G4cerr<<"***G4QFragmentation::Breeder: **No SSCorrection**,next="<<next<<G4endl;
#endif
      G4int nofRH=theResult->size();           // #of resulting Hadrons
      if(!theHadrons && nofRH)                 // Hadrons are existing for Correction
      {
#ifdef debug
        G4cerr<<"!G4QFragmentation::Breeder:Can try SHCor, nH="<<theResult->size()<<G4endl;
#endif
        // @@ string can be not convertable to one hadron (2,0.0,0,2,0) - To be improved
        G4QContent miQC=curString->GetQC();    // QContent of the Lightest Hadron
        G4int miPDG=miQC.GetSPDGCode();        // PDG of the Lightest Hadron
        G4double miM=0.;                       // Prototype of Mass of theCurLightestHadron
        if(miPDG!=10) miM=G4QPDGCode(miPDG).GetMass(); // Mass of theCurLightestOneHadron
        else
        {
          G4QChipolino QCh(miQC);              // define the TotalString as a Chipolino
          miM=QCh.GetQPDG1().GetMass()+QCh.GetQPDG2().GetMass();//MinMass of theStringChipo
        }
        G4double cM=curString4M.m();           // Actual mass of the Cur String
        G4ThreeVector cP=curString4M.vect();   // Momentum of the curString
        G4double      cE=curString4M.e();      // Energy of the curString
        G4ThreeVector curV=cP/cE;              // curRelativeVelocity
        G4int reha=0;                          // Hadron # to use beyon the LOOP
        G4int fuha=0;                          // Selected Hadron-partner (0 = NotFound)
        G4double Vmin=DBL_MAX;                 // min Velocity Distance
        G4LorentzVector s4M(0.,0.,0.,0.);      // Selected 4-mom of the hadron+string
        for (reha=next; reha < nofRH; reha++)
        {
          G4QHadron* pHadron=(*theResult)[reha];
          G4LorentzVector p4M=pHadron->Get4Momentum();
          G4double         pM=p4M.m();         // Mass of the Hadron
          G4LorentzVector t4M=p4M+curString4M; // Total momentum
          if(t4M.m2() > sqr(pM+cM))
          {
            G4ThreeVector pP=p4M.vect();       // Momentum of the partnerString
            G4double      pE=p4M.e();          // Energy of the partnerString
            G4ThreeVector pV=pP/pE;            // curRelativeVelocity
            G4double dV=(curV-pV).mag2();      // SquaredDifferenceBetweenRelVelocities
            if(dV < Vmin)
            {
              Vmin=dV;
              fuha=reha;
              s4M=t4M;
            }
          }
        } // End of the LOOP over string-partners for Correction
        // @@@@@@@@@@@@@@@ G4QHadron* pHadron=(*theResult)[fuha]; // Necessary for update
        // @@@@@@@@@@@@@@@@@@@@ Decay the compound in minString & Hadron
        if(std::abs(miPDG)%10 > 2)                  // Decay Hadron-Resonans
        {
          G4Quasmon Quasm;
          G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
          G4QHadronVector* tmpQHadVec=Quasm.DecayQHadron(sHad); // It deletes sHad
#ifdef pdebug
          G4cout<<"G4QFragmentation::Breeder: DecStrHad,nH="<<tmpQHadVec->size()<<G4endl;
#endif
          for(unsigned aH=0; aH < tmpQHadVec->size(); aH++)
          {
            theResult->push_back((*tmpQHadVec)[aH]);// TheDecayProductofTheHadronIsFilled
#ifdef debug
            G4QHadron*   prodH =(*tmpQHadVec)[aH];
            G4LorentzVector p4M=prodH->Get4Momentum();
            G4int           PDG=prodH->GetPDGCode();
            G4int           Chg=prodH->GetCharge();
            G4int           BaN=prodH->GetBaryonNumber();
            G4cout<<"-EMC->>>G4QFragmentation::Breeder:Str=Hadr PrH#"<<aH<<" filled, 4M="
                  <<p4M<<", PDG="<<PDG<<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
          }
          tmpQHadVec->clear();
          delete tmpQHadVec;  // Who calls DecayQHadron is responsible for clear & delete
        }
        else if(miPDG == 10)                   // ==> Decay Hadron-Chipolino
        {
          G4QChipolino QCh(miQC);              // define the TotalResidual as a Chipolino
          G4QPDGCode   h1QPDG=QCh.GetQPDG1();  // QPDG of the first hadron
          G4double     h1M   =h1QPDG.GetMass();// Mass of the first hadron
          G4QPDGCode   h2QPDG=QCh.GetQPDG2();  // QPDG of the second hadron 
          G4double     h2M   =h2QPDG.GetMass();// Mass of the second hadron
          G4double     ttM   =curString4M.m(); // Real Mass of the Chipolino
          if(h1M+h2M<ttM+eps)                  // Two particles decay of Chipolino
          {
            G4LorentzVector h14M(0.,0.,0.,h1M);
            G4LorentzVector h24M(0.,0.,0.,h2M);
            if(std::fabs(ttM-h1M-h2M)<=eps)
            {
              G4double part1=h1M/(h1M+h2M);
              h14M=part1*curString4M;
              h24M=curString4M-h14M;
            }
            else
            {
              if(!G4QHadron(curString4M).DecayIn2(h14M,h24M))
              {
                G4cerr<<"***G4QFragmentation::Breeder: tM="<<ttM<<"->h1="<<h1QPDG<<"("
                      <<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")="<<h1M+h2M<<G4endl;
                G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ChipoDec");
              }
            }
            G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
            theResult->push_back(h1H);         // (delete equivalent)  
#ifdef debug
            G4LorentzVector f4M=h1H->Get4Momentum();
            G4int           fPD=h1H->GetPDGCode();
            G4int           fCg=h1H->GetCharge();
            G4int           fBN=h1H->GetBaryonNumber();
            G4cout<<"-EMC->>>G4QFragmentation::Breeder:Str=Hadr Prod-F is filled, f4M="
                  <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
            G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
            theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
            G4LorentzVector s4M=h2H->Get4Momentum();
            G4int           sPD=h2H->GetPDGCode();
            G4int           sCg=h2H->GetCharge();
            G4int           sBN=h2H->GetBaryonNumber();
            G4cout<<"-EMC->>>G4QFragmentation::Breeder:Str=Hadr Prod-S is filled, s4M="
                  <<s4M<<", sPDG="<<sPD<<", sCg="<<sCg<<", sBN="<<sBN<<G4endl;
#endif
#ifdef edebug
            G4cout<<"-EMC-..Chipo..G4QFragmentation::Breeder:DecCHECK,Ch4M="<<curString4M
                  <<", ChQC="<<miQC<<", d4M="<<curString4M-h14M-h24M<<G4endl;
#endif
          }
          else
          {
            G4cerr<<"***G4QFragm::Breeder: tM="<<ttM<<miQC<<"->h1="<<h1QPDG<<"(" <<h1M
                  <<")+h2="<<h1QPDG<<"("<<h2M<<") = "<<h1M+h2M<<G4endl;
            G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ChipoDecMass");
          }
        }
        else
        {
          G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
          theResult->push_back(sHad);             // The original string-hadron is filled
#ifdef debug
          G4cout<<"-EMC->>>G4QFragmentation::Breeder: Str=Hadr Filled, 4M="<<curString4M
                <<", PDG="<<miPDG<<G4endl;
#endif
        }
#ifdef edebug
        G4cout<<"-EMC-...Cor...G4QFragmentation::Breeder:CorCHECK Sum="<<s4M
              <<" =? "<<curString4M+pHadron->Get4Momentum()<<", M="<<miM<<" ?= "
              <<curString4M.m()<<G4endl;
#endif
#ifdef debug
        G4cout<<">>>G4QFragmentation::Breeder:*Corrected* String->Hadr="<<miPDG
              <<curString4M<<" by String #"<<reha<<G4endl;
#endif
        continue;                    // Continue the LOOP over the curStrings
      }
    } // End of IF(NO_Hadrons) = Problem solution namespace
    G4Quasmon tmpQ;                                 // @@ an issue of Q to decay resonances
    for(unsigned aTrack=0; aTrack<theHadrons->size(); aTrack++)
    {
      G4QHadron* curHadron=(*theHadrons)[aTrack];
      G4int hPDG=curHadron->GetPDGCode();
#ifdef edebug
      G4LorentzVector curH4M=curHadron->Get4Momentum();
      G4int           curHCh=curHadron->GetCharge();
      G4int           curHBN=curHadron->GetBaryonNumber();
#endif
#ifdef debug
      G4cout<<">>>>>>>>G4QFragmentation::Breeder:S#"<<astring<<",H#"<<aTrack<<",PDG="<<hPDG
            <<",4M="<<curHadron->Get4Momentum()<<G4endl;
#endif
      if(std::abs(hPDG)%10 > 2)
      {
        G4QHadronVector* tmpQHadVec=tmpQ.DecayQHadron(curHadron); // It deletes curHadron
#ifdef pdebug
        G4cout<<"G4QFragmentation::Breeder:-DECAY'S DONE-,nH="<<tmpQHadVec->size()<<G4endl;
#endif
        //G4int tmpS=tmpQHadVec->size(); // "The elegant method" (tested) is commented
        //theResult->resize(tmpS+theResult->size()); // Resize theQHadrons length
        //copy(tmpQHadVec->begin(), tmpQHadVec->end(), theResult->end()-tmpS);
        for(unsigned aH=0; aH < tmpQHadVec->size(); aH++)
        {
          theResult->push_back((*tmpQHadVec)[aH]);// TheDecayProduct of TheHadron is filled
#ifdef edebug
          G4QHadron*   prodH =(*tmpQHadVec)[aH];
          G4LorentzVector p4M=prodH->Get4Momentum();
          G4int           PDG=prodH->GetPDGCode();
          G4int           Chg=prodH->GetCharge();
          G4int           BaN=prodH->GetBaryonNumber();
          curH4M-=p4M;
          curString4M-=p4M;
          curStrChg-=Chg;
          curStrBaN-=BaN;
          curHCh-=Chg;
          curHBN-=BaN;
          G4cout<<"-EMC->>>>G4QFragmentation::Breeder:Str*Filled, 4M="<<p4M<<", PDG="<<PDG
                <<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
        }
#ifdef edebug
        G4cout<<"-EMC-.G4QFr::Br:Dec,r4M="<<curH4M<<",rC="<<curHCh<<",rB="<<curHBN<<G4endl;
#endif
        tmpQHadVec->clear();
        delete tmpQHadVec;  // Who calls DecayQHadron is responsible for clear & delete
      }
      else
      {
        theResult->push_back(curHadron);        // The original hadron is filled
#ifdef edebug
        curString4M-=curH4M;
        G4int curCh=curHadron->GetCharge();
        G4int curBN=curHadron->GetBaryonNumber();
        curStrChg-=curCh;
        curStrBaN-=curBN;
        G4cout<<"-EMC->>>>>>G4QFragmentation::Breeder: curH filled 4M="<<curH4M<<",PDG="
              <<curHadron->GetPDGCode()<<", Chg="<<curCh<<", BaN="<<curBN<<G4endl;
#endif
      }
    }
    G4LorentzVector r4M=theNucleus->GetNucleons4Momentum();// Sum of N4M's in RotatedLS=LS
    // clean up (the issues are filled to theResult)
    delete theHadrons;
#ifdef edebug
    G4cout<<"-EMC-.........G4QFragmentation::Breeder: StringDecay CHECK, r4M="<<curString4M
          <<", rChg="<<curStrChg<<", rBaN="<<curStrBaN<<G4endl;
#endif
  } // End of the main LOOP over decaying strings
  G4LorentzVector r4M=theNucleus->Get4Momentum(); // Nucleus 4-momentum in LS
  G4int rPDG=theNucleus->GetPDG();
  G4QHadron* resNuc = new G4QHadron(rPDG,r4M);
  theResult->push_back(resNuc);                          // Fill the residual nucleus
#ifdef edebug
  // @@ Should be tested for the primary projectile not along Z !
  G4LorentzVector s4M(0.,0.,0.,0.); // Sum of the Result in LS
  G4int rCh=totChg;
  G4int rBN=totBaN;
  G4int nHadr=theResult->size();
  G4cout<<"-EMCLS-G4QFragmentation::Breeder:#ofHadr="<<nHadr<<",rN="<<r4M.m()<<"="
        <<G4QNucleus(rPDG).GetGSMass()<<G4endl;
  for(G4int i=0; i<nHadr; i++)
  {
    G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
    s4M+=hI4M;
    G4int hChg=(*theResult)[i]->GetCharge();
    rCh-=hChg;
    G4int hBaN=(*theResult)[i]->GetBaryonNumber();
    rBN-=hBaN;
    G4cout<<"-EMCLS-G4QFragmentation::Breeder: H#"<<i<<", 4M="<<hI4M<<", PDG="
          <<(*theResult)[i]->GetPDGCode()<<", C="<<hChg<<", B="<<hBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFragm::Breed: LS r4M="<<s4M-totLS4M<<",rC="<<rCh<<",rB="<<rBN<<G4endl;
#endif
  //
  // >>>>>>>>>>>>>> clean-up used Nucleons, Strings. and theNucleus
  //
  delete theNucleus;
  std::for_each(strings->begin(), strings->end(), DeleteQString() );
  delete strings;

  return theResult; // This is the resulting hadron vector with the resNucleus (the answer)
}

// Excite double diffractive string
G4bool G4QFragmentation::ExciteDiffParticipants(G4QHadron* projectile,
                                                G4QHadron* target) const
{
  G4LorentzVector Pprojectile=projectile->Get4Momentum();
  G4double Mprojectile=projectile->GetMass();
  G4double Mprojectile2=Mprojectile*Mprojectile;
  G4LorentzVector Ptarget=target->Get4Momentum();
  G4double Mtarget=target->GetMass();
  G4double Mtarget2=Mtarget*Mtarget;
#ifdef debug
  G4cout<<"G4QFragm::ExciteDiffPartici:Ep="<<Pprojectile.e()<<",Et="<<Ptarget.e()<<G4endl;
#endif
  // Transform momenta to cms and then rotate parallel to z axis;
  G4LorentzVector Psum=Pprojectile+Ptarget;
  G4LorentzRotation toCms(-Psum.boostVector()); // Boost Rotation to CMS
  G4LorentzVector Ptmp=toCms*Pprojectile;
  if(Ptmp.pz()<=0.) // "String" moving backwards in CMS, abort collision !! ? M.K.
  {
#ifdef debug
    G4cout<<"G4QFragmentation::ExciteDiffParticipants: *1* abort Collision!! *1*"<<G4endl;
#endif
    return false; 
  }         
  toCms.rotateZ(-Ptmp.phi());
  toCms.rotateY(-Ptmp.theta());
#ifdef debug
  G4cout<<"G4QFragment::ExciteDiffParticipantts: Be4Boost Pproj="<<Pprojectile<<", Ptarg="
        <<Ptarget<<G4endl;
#endif
  G4LorentzRotation toLab(toCms.inverse()); // Boost Rotation to LabSys (LS)
  Pprojectile.transform(toCms);
  Ptarget.transform(toCms);
#ifdef debug
  G4cout<< "G4QFragment::ExciteDiffParticipantts: AfterBoost Pproj="<<Pprojectile<<"Ptarg="
        <<Ptarget<<", cms4M="<<Pprojectile+Ptarget<<G4endl;
  G4cout<<"G4QFragment::ExciteDiffParticipantts: ProjX+="<<Pprojectile.plus()<<", ProjX-="
        <<Pprojectile.minus()<<", TargX+="<< Ptarget.plus()<<", TargX-="<<Ptarget.minus()
        <<G4endl;
#endif
  G4LorentzVector Qmomentum(0.,0.,0.,0.);
  G4int whilecount=0;
  do
  {
    //  Generate pt  
    G4double maxPtSquare=sqr(Ptarget.pz());
    if(whilecount++>=500 && whilecount%100==0) // @@ M.K. Hardwired limits 
#ifdef debug
    G4cout<<"G4QFragmentation::ExciteDiffParticipantts: can loop, loopCount="<<whilecount
          <<", maxPtSquare="<<maxPtSquare<<G4endl;
#endif
    if(whilecount>1000)                        // @@ M.K. Hardwired limits 
    {
#ifdef debug
      G4cout<<"G4QFragmentation::ExciteDiffParticipants: *2* abort Loop!! *2*"<<G4endl;
#endif
      return false;    //  Ignore this interaction 
    }
    Qmomentum=G4LorentzVector(GaussianPt(widthOfPtSquare,maxPtSquare),0);
#ifdef debug
    G4cout<<"G4QFragment::ExciteDiffParticipants: generated Pt="<<Qmomentum<<", ProjPt="
          <<Pprojectile+Qmomentum<<", TargPt="<<Ptarget-Qmomentum<<G4endl;
#endif
    //  Momentum transfer
    G4double Xmin=0.;
    G4double Xmax=1.;
    G4double Xplus =ChooseX(Xmin,Xmax);
    G4double Xminus=ChooseX(Xmin,Xmax);
#ifdef debug
    G4cout<<"G4QFragment::ExciteDiffParticip: X-plus="<<Xplus<<",X-minus="<<Xminus<<G4endl;
#endif
    G4double pt2=Qmomentum.vect().mag2();
    G4double Qplus =-pt2/Xminus/Ptarget.minus();
    G4double Qminus= pt2/Xplus /Pprojectile.plus();
    Qmomentum.setPz((Qplus-Qminus)/2);
    Qmomentum.setE( (Qplus+Qminus)/2);
#ifdef debug
    G4cout<<"G4QFragment::ExciteDiffParticip: Qplus="<<Qplus<<", Qminus="<<Qminus<<", pt2="
          <<pt2<<", Qmomentum="<<Qmomentum<<", ProjM="<<(Pprojectile+Qmomentum).mag()
          <<", TargM="<<(Ptarget-Qmomentum).mag()<<G4endl;
#endif
  } while((Pprojectile+Qmomentum).mag2()<=Mprojectile2 ||
          (Ptarget-Qmomentum).mag2()<=Mtarget2);
  Pprojectile += Qmomentum;
  Ptarget     -= Qmomentum;
#ifdef debug
  G4cout<<"G4QFragment::ExciteDiffParticipan: Proj(Q)="<<Pprojectile<<", Targ(Q)="<<Ptarget
        <<", Proj(back)="<<toLab*Pprojectile<<", Targ(bac)="<< toLab*Ptarget << G4endl;
#endif
  // Transform back and update SplitableHadron Participant.
  Pprojectile.transform(toLab);
  Ptarget.transform(toLab);
#ifdef debug
  G4cout<< "G4QFragmentation::ExciteDiffParticipants: TargetMass="<<Ptarget.mag()<<G4endl;
#endif
  target->Set4Momentum(Ptarget);  
#ifdef debug
  G4cout<<"G4QFragment::ExciteDiffParticipants:ProjectileMass="<<Pprojectile.mag()<<G4endl;
#endif
  projectile->Set4Momentum(Pprojectile);
  return true;
} // End of ExciteDiffParticipants


// Excite single diffractive string
G4bool G4QFragmentation::ExciteSingDiffParticipants(G4QHadron* projectile,
                                                    G4QHadron* target) const
{
  G4LorentzVector Pprojectile=projectile->Get4Momentum();
  G4double Mprojectile=projectile->GetMass();
  G4double Mprojectile2=Mprojectile*Mprojectile;
  G4LorentzVector Ptarget=target->Get4Momentum();
  G4double Mtarget=target->GetMass();
  G4double Mtarget2=Mtarget*Mtarget;
#ifdef debug
  G4cout<<"G4QFragm::ExSingDiffPartici:Ep="<<Pprojectile.e()<<",Et="<<Ptarget.e()<<G4endl;
#endif
  G4bool KeepProjectile= G4UniformRand() > 0.5;
  // Reset minMass of the non diffractive particle to its value, (minus for rounding...)
  if(KeepProjectile ) 
  {
#ifdef debug
    G4cout<<"--1/2--G4QFragmentation::ExSingDiffParticipants: Projectile is fixed"<<G4endl;
#endif
    Mprojectile2 = projectile->GetMass2()*(1.-perCent); // Isn't it too big reduction? M.K.
  }
  else
  {
#ifdef debug
    G4cout<<"---1/2---G4QFragmentation::ExSingDiffParticipants: Target is fixed"<<G4endl;
#endif
    Mtarget2 = target->GetMass2()*(1.-perCent); // @@ Isn't it too big reduction? M.K.
  }
  // @@ From this point it repeats the Diffractional excitation (? Use flag ?)
  // Transform momenta to cms and then rotate parallel to z axis;
  G4LorentzVector Psum=Pprojectile+Ptarget;
  G4LorentzRotation toCms(-Psum.boostVector()); // Boost Rotation to CMS
  G4LorentzVector Ptmp=toCms*Pprojectile;
  if(Ptmp.pz()<=0.) // "String" moving backwards in CMS, abort collision !! ? M.K.
  {
#ifdef debug
    G4cout<<"G4QFragment::ExciteSingDiffParticipants: *1* abort Collision!! *1*"<<G4endl;
#endif
    return false; 
  }         
  toCms.rotateZ(-Ptmp.phi());
  toCms.rotateY(-Ptmp.theta());
#ifdef debug
  G4cout<<"G4QFragm::ExciteSingDiffParticipantts: Be4Boost Pproj="<<Pprojectile<<",Ptarg="
        <<Ptarget<<G4endl;
#endif
  G4LorentzRotation toLab(toCms.inverse()); // Boost Rotation to LabSys (LS)
  Pprojectile.transform(toCms);
  Ptarget.transform(toCms);
#ifdef debug
  G4cout<< "G4QFragment::ExciteDiffParticipantts: AfterBoost Pproj="<<Pprojectile<<"Ptarg="
        <<Ptarget<<", cms4M="<<Pprojectile+Ptarget<<G4endl;

  G4cout<<"G4QFragment::ExciteDiffParticipantts: ProjX+="<<Pprojectile.plus()<<", ProjX-="
        <<Pprojectile.minus()<<", TargX+="<< Ptarget.plus()<<", TargX-="<<Ptarget.minus()
        <<G4endl;
#endif
  G4LorentzVector Qmomentum(0.,0.,0.,0.);
  G4int whilecount=0;
  do
  {
    //  Generate pt  
    G4double maxPtSquare=sqr(Ptarget.pz());
    if(whilecount++>=500 && whilecount%100==0) // @@ M.K. Hardwired limits 
#ifdef debug
    G4cout<<"G4QFragment::ExciteSingDiffParticipantts: can loop, loopCount="<<whilecount
          <<", maxPtSquare="<<maxPtSquare<<G4endl;
#endif
    if(whilecount>1000)                        // @@ M.K. Hardwired limits 
    {
#ifdef debug
      G4cout<<"G4QFragmentation::ExciteSingDiffParticipants: *2* abort Loop!! *2*"<<G4endl;
#endif
      return false;    //  Ignore this interaction 
    }
    Qmomentum=G4LorentzVector(GaussianPt(widthOfPtSquare,maxPtSquare),0);
#ifdef debug
    G4cout<<"G4QFragm::ExciteSingDiffParticipants: generated Pt="<<Qmomentum<<", ProjPt="
          <<Pprojectile+Qmomentum<<", TargPt="<<Ptarget-Qmomentum<<G4endl;
#endif
    //  Momentum transfer
    G4double Xmin=0.;
    G4double Xmax=1.;
    G4double Xplus =ChooseX(Xmin,Xmax);
    G4double Xminus=ChooseX(Xmin,Xmax);
#ifdef debug
    G4cout<<"G4QFragm::ExciteSingDiffPartici:X-plus="<<Xplus<<",X-minus="<<Xminus<<G4endl;
#endif
    G4double pt2=G4ThreeVector(Qmomentum.vect()).mag2();
    G4double Qplus =-pt2/Xminus/Ptarget.minus();
    G4double Qminus= pt2/Xplus /Pprojectile.plus();
    if (KeepProjectile)
      Qminus=(projectile->GetMass2()+pt2)/(Pprojectile.plus()+Qplus) - Pprojectile.minus();
    else Qplus=Ptarget.plus() - (target->GetMass2()+pt2)/(Ptarget.minus()-Qminus);  
    Qmomentum.setPz((Qplus-Qminus)/2);
    Qmomentum.setE( (Qplus+Qminus)/2);
#ifdef debug
    G4cout<<"G4QFragm::ExciteDiffParticip: Qplus="<<Qplus<<", Qminus="<<Qminus<<", pt2="
          <<pt2<<", Qmomentum="<<Qmomentum<<", ProjM="<<(Pprojectile+Qmomentum).mag()
          <<", TargM="<<(Ptarget-Qmomentum).mag()<<G4endl;
#endif
    // while is different from the Double Diffractive Excitation (@@ !)
    //} while((Pprojectile+Qmomentum).mag2()<= Mprojectile2 ||
    //        (Ptarget-Qmomentum).mag2()<=Mtarget2);
  } while((Ptarget-Qmomentum).mag2()<=Mtarget2 ||
          (Pprojectile+Qmomentum).mag2()<=Mprojectile2 ||
          (Ptarget-Qmomentum).e() < 0. || (Pprojectile+Qmomentum).e() < 0.);
  Pprojectile += Qmomentum;
  Ptarget     -= Qmomentum;
#ifdef debug
  G4cout<<"G4QFragmentation::ExciteSingDiffParticipan: Proj(Q)="<<Pprojectile<<"(E="
        <<Pprojectile.e()<<"), Targ(Q)="<<Ptarget<<"(E="<<Ptarget.e()
        <<"), Proj(back)="<<toLab*Pprojectile<<", Targ(bac)="<< toLab*Ptarget << G4endl;
#endif
  // Transform back and update SplitableHadron Participant.
  Pprojectile.transform(toLab);
  Ptarget.transform(toLab);
#ifdef debug
  G4cout<< "G4QFragm::ExciteSingleDiffParticipants: TargetMass="<<Ptarget.mag()<<G4endl;
#endif
  target->Set4Momentum(Ptarget);  
#ifdef debug
  G4cout<<"G4QFragm::ExciteSingleParticipants:ProjectileMass="<<Pprojectile.mag()<<G4endl;
#endif
  projectile->Set4Momentum(Pprojectile);
  return true;
} // End of ExciteSingleDiffParticipants

void G4QFragmentation::SetParameters(G4int nCM, G4double radNuc, G4double SigPt)
{//  =============================================================================
  nCutMax            = nCM;            // max number of pomeron cuts
  theNucleonRadius   = radNuc;         // effective radius of the nucleon inside Nucleus
  widthOfPtSquare    = -2*SigPt*SigPt; // width^2 of pt for string excitation
}

G4double G4QFragmentation::ChooseX(G4double Xmin, G4double Xmax) const
{
  // choose an x between Xmin and Xmax with P(x) ~ 1/x @@ M.K. -> 1/sqrt(x)
  //G4double range=Xmax-Xmin;
  if(Xmax == Xmin) return Xmin;
  if( Xmin < 0. || Xmax < Xmin) 
  {
    G4cerr<<"***G4QFragmentation::ChooseX: Xmin="<<Xmin<<", Xmax="<<Xmax<< G4endl;
    G4Exception("G4QFragmentation::ChooseX:","72",FatalException,"Bad X or X-Range");
  }
  //G4double x;
  //do {x=Xmin+G4UniformRand()*range;} while ( Xmin/x < G4UniformRand() );
  G4double sxi=std::sqrt(Xmin);
  G4double x=sqr(sxi+G4UniformRand()*(std::sqrt(Xmax)-sxi));
#ifdef debug
  G4cout<<"G4QFragmentation::ChooseX: DiffractiveX="<<x<<G4endl;
#endif
  return x;
} // End of ChooseX

// Pt distribution @@ one can use 1/(1+A*Pt^2)^B
G4ThreeVector G4QFragmentation::GaussianPt(G4double widthSq, G4double maxPtSquare) const
{
  G4double pt2=0.;
  do{pt2=widthSq*std::log(G4UniformRand());} while (pt2>maxPtSquare);
  pt2=std::sqrt(pt2);
  G4double phi=G4UniformRand()*twopi;
  return G4ThreeVector(pt2*std::cos(phi),pt2*std::sin(phi),0.);    
} // End of GaussianPt

G4int G4QFragmentation::SumPartonPDG(G4int PDG1, G4int PDG2) const
{
  if      (PDG1 < 7 && PDG1 > 0 && PDG2 < 7 && PDG2 > 0) // Sum up two Q in DiQ (S=0)
  {
    if(PDG1 > PDG2) return PDG1*1000+PDG2*100+1;
    else            return PDG2*1000+PDG1*100+1;
  }
  else if (PDG1 >-7 && PDG1 < 0 && PDG2 >-7 && PDG2 < 0) // Sum up two AQ in ADiQ (S=0)
  {
    if(-PDG1 > -PDG2) return PDG1*1000+PDG2*100-1;
    else              return PDG2*1000+PDG1*100-1;
  }
  else if (PDG1 <-99 && PDG2 < 7 && PDG2 > 0) // Sum up Q and ADiQ in AQ
  {
    G4int PDG=-PDG1/100;
    if(PDG2==PDG/10) return -PDG%10;
    if(PDG2==PDG%10) return -PDG/10;
    else
    {
      G4cerr<<"***4QFragmentation::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
      G4Exception("G4QFragmentation::SumPartonPDG:","72",FatalException,"Q&ADiQ notMatch");
    }
  }
  else if (PDG2 <-99 && PDG1 < 7 && PDG1 > 0) // Sum up ADiQ and Q in AQ
  {
    G4int PDG=-PDG2/100;
    if(PDG1==PDG/10) return -PDG%10;
    if(PDG1==PDG%10) return -PDG/10;
    else
    {
      G4cerr<<"***4QFragmentation::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
      G4Exception("G4QFragmentation::SumPartonPDG:","72",FatalException,"ADiQ&Q notMatch");
    }
  }
  else if (PDG1 > 99 && PDG2 >-7 && PDG2 < 0) // Sum up DiQ and AQ in Q
  {
    G4int PDG=PDG1/100;
    if(PDG2==-PDG/10) return PDG%10;
    if(PDG2==-PDG%10) return PDG/10;
    else
    {
      G4cerr<<"***4QFragmentation::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
      G4Exception("G4QFragmentation::SumPartonPDG:","72",FatalException,"DiQ&AQ notMatch");
    }
  }
  else if (PDG2 > 99 && PDG1 >-7 && PDG1 < 0) // Sum up AQ and DiQ in Q
  {
    G4int PDG=PDG2/100;
    if(PDG1==-PDG/10) return PDG%10;
    if(PDG1==-PDG%10) return PDG/10;
    else
    {
      G4cerr<<"***G4QFragmentation::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
      G4Exception("G4QFragmentation::SumPartonPDG:","72",FatalException,"AQ&DiQ notMatch");
    }
  }
  else
  {
    G4cerr<<"***G4QFragmentation::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
    G4Exception("G4QFragmentation::SumPartonPDG:","72",FatalException,"Can'tSumUpPartons");
  }
  return 0;
} // End of SumPartonPDG

// Reduces quark pairs (unsigned 2 digits) to quark singles (unsigned)
std::pair<G4int,G4int> G4QFragmentation::ReducePair(G4int P1, G4int P2) const
{
#ifdef debug
  G4cout<<"G4QFragmentation::ReducePair: **Called** P1="<<P1<<", P2="<<P2<<G4endl;
#endif
  G4int P11=P1/10;
  G4int P12=P1%10;
  G4int P21=P2/10;
  G4int P22=P2%10;
  if     (P11==P21) return std::make_pair(P12,P22);
  else if(P11==P22) return std::make_pair(P12,P21);
  else if(P12==P21) return std::make_pair(P11,P22);
  else if(P12==P22) return std::make_pair(P11,P21);
#ifdef debug
  G4cout<<"-Warning-G4QFragmentation::ReducePair: Reduction Failed"<<G4endl;
#endif
  return std::make_pair(0,0);                         // Reduction failed
}

// Select LL/RR (1) or LR/RL (-1) annihilation order (0, if the annihilation is impossible)
G4int G4QFragmentation::AnnihilationOrder(G4int LS, G4int MPS, G4int uPDG, G4int mPDG,
                                          G4int sPDG, G4int nPDG) // ^L^^ Partner^^R^
{//                                             ^L^^ Curent ^^R^
  G4int Ord=0;
  //       Curent   Partner
  if      (LS==2 && MPS==2 )                 // Fuse 2 Q-aQ strings to DiQ-aDiQ
  {
#ifdef debug
    G4cout<<"G4QFragmentation::AnnihilationOrder: QaQ/QaQ->DiQ-aDiQ, uPDG="<<uPDG<<G4endl;
#endif
    if     ( (uPDG>0 && sPDG>0 && mPDG<0 && nPDG<0) || 
             (uPDG<0 && sPDG<0 && mPDG>0 && nPDG>0) ) Ord= 1; // LL/RR
    else if( (uPDG>0 && nPDG>0 && mPDG<0 && sPDG<0) || 
             (uPDG<0 && nPDG<0 && mPDG>0 && sPDG>0) ) Ord=-1; // LR/RL
    else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 22 fusion, L="<<uPDG
               <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
  }
  else if ( LS==2 && MPS==3 )
  {
    if     (uPDG > 7)                                // Fuse pLDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihOrder: pLDiQ, sPDG="<<sPDG<<", nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<0 && (-sPDG==uPDG/1000 || -sPDG==(uPDG/100)%10) ) Ord= 1; // LL/RR
      else if( nPDG<0 && (-nPDG==uPDG/1000 || -nPDG==(uPDG/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong pLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (mPDG > 7)                               // Fuse pRDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihOrder: pRDiQ, sPDG="<<sPDG<<", nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<0 && (-sPDG==mPDG/1000 || -sPDG==(mPDG/100)%10) ) Ord=-1; // LR/RL
      else if( nPDG<0 && (-nPDG==mPDG/1000 || -nPDG==(mPDG/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong pRDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (uPDG <-7)                               // Fuse pLaDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihOrder: pLaDiQ, sPDG="<<sPDG<<", nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG>0 && (sPDG==(-uPDG)/1000 || sPDG==((-uPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( nPDG>0 && (nPDG==(-uPDG)/1000 || nPDG==((-uPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong pLaDiQ, L="<<uPDG
                 <<", R="<<mPDG<<", cL="<<sPDG<<", cR="<<nPDG<<G4endl;
    }
    else if (mPDG <-7)                              // Fuse pRaDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihOrder: pRaDiQ, sPDG="<<sPDG<<", nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG>0 && (sPDG==(-mPDG)/1000 || sPDG==((-mPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( nPDG>0 && (nPDG==(-mPDG)/1000 || nPDG==((-mPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong pRaDiQ, L="<<uPDG
                 <<", R="<<mPDG<<", cL="<<sPDG<<", cR="<<nPDG<<G4endl;
    }
    else if( (sPDG<0 && (-sPDG==mPDG || -sPDG==uPDG) ) ||
             (nPDG<0 && (-nPDG==mPDG || -nPDG==uPDG) ) ) Ord= 2; // @@ Annihilation fusion
#ifdef debug
    else G4cout<<"-Warning-G4QFragmentation::AnnihilatOrder: Wrong 23StringFusion"<<G4endl;
    G4cout<<"G4QFragmentation::AnnihilationOrder: Ord="<<Ord<<",sPDG="<<sPDG<<",nPDG="
          <<nPDG<<", uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
  }
  else if ( LS==3 && MPS==2 )
  {
    if     (sPDG > 7)                                // Fuse cLDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihOrder: cLDiQ, uPDG="<<uPDG<<", mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<0 && (-uPDG==sPDG/1000 || -uPDG==(sPDG/100)%10) ) Ord= 1; // LL/RR
      else if( mPDG<0 && (-mPDG==sPDG/1000 || -mPDG==(sPDG/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong cLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (nPDG > 7)                               // Fuse cRDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihOrder: cRDiQ, uPDG="<<uPDG<<", mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<0 && (-uPDG==nPDG/1000 || -uPDG==(nPDG/100)%10) ) Ord=-1; // LR/RL
      else if( mPDG<0 && (-mPDG==nPDG/1000 || -mPDG==(nPDG/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong cRDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (sPDG <-7)                               // Fuse cLaDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihOrder: cLaDiQ, uPDG="<<uPDG<<", mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG>0 && (uPDG==(-sPDG)/1000 || uPDG==((-sPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( mPDG>0 && (mPDG==(-sPDG)/1000 || mPDG==((-sPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong cLaDiQ, L="<<uPDG
                 <<", R="<<mPDG<<", cL="<<sPDG<<", cR="<<nPDG<<G4endl;
    }
    else if (nPDG <-7)                              // Fuse cRaDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihOrder: cRaDiQ, uPDG="<<uPDG<<", mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG>0 && (uPDG==(-nPDG)/1000 || uPDG==((-nPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( mPDG>0 && (mPDG==(-nPDG)/1000 || mPDG==((-nPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong cRaDiQ, L="<<uPDG
                 <<", R="<<mPDG<<", cL="<<sPDG<<", cR="<<nPDG<<G4endl;
    }
    else if( (uPDG<0 && (-uPDG==sPDG || -uPDG==nPDG) ) ||
             (mPDG<0 && (-mPDG==sPDG || -mPDG==nPDG) ) ) Ord=2; // @@ Annihilation fusion
#ifdef debug
    else G4cout<<"-Warning-G4QFragmentation::AnnihilatOrder: Wrong 32StringFusion"<<G4endl;
    G4cout<<"G4QFragmentation::AnnihilationOrder: Ord="<<Ord<<",sPDG="<<sPDG<<",nPDG="
          <<nPDG<<", uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
  }
  else if ( (LS==2 && MPS==4) || (LS==4 && MPS==2) )
  {
    if     (uPDG > 7)  // Annihilate pLDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihilOrder:pLDiQ, sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<0 && (-sPDG==uPDG/1000 || -sPDG==(uPDG/100)%10) &&
               (nPDG==(-mPDG)/1000 || nPDG==((-mPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( nPDG<0 && (-nPDG==uPDG/1000 || -nPDG==(uPDG/100)%10) &&
               (sPDG==(-mPDG)/1000 || sPDG==((-mPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 24 pLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (mPDG >7) // Annihilate pRDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihilOrder:PRDiQ, sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<0 && (-sPDG==mPDG/1000 || -sPDG==(mPDG/100)%10) &&
               (nPDG==(-uPDG)/1000 || nPDG==((-uPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( nPDG<0 && (-nPDG==mPDG/1000 || -nPDG==(mPDG/100)%10) &&
               (sPDG==(-uPDG)/1000 || sPDG==((-uPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 24 pLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (sPDG > 7)   // Annihilate cLDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihilOrder:cLDiQ, uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<0 && (-uPDG==sPDG/1000 || -uPDG==(sPDG/100)%10) &&
               (mPDG==(-nPDG)/1000 || mPDG==((-nPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( mPDG<0 && (-mPDG==sPDG/1000 || -mPDG==(sPDG/100)%10) &&
               (uPDG==(-nPDG)/1000 || uPDG==((-nPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 24 cLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (nPDG > 7)   // Annihilate cRDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihilOrder:cRDiQ, uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<0 && (-uPDG==nPDG/1000 || -uPDG==(nPDG/100)%10) &&
               (mPDG==(-sPDG)/1000 || mPDG==((-sPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( mPDG<0 && (-mPDG==nPDG/1000 || -mPDG==(nPDG/100)%10) &&
               (uPDG==(-sPDG)/1000 || uPDG==((-sPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 24 cRDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
#ifdef debug
    else G4cout<<"-Warning-G4QFragmentation::AnnihilOrder: Wrong 24 StringFusion"<<G4endl;
    G4cout<<"G4QFragmentation::AnnihilationOrder: Ord="<<Ord<<",sPDG="<<sPDG<<",nPDG="
          <<nPDG<<", uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
  }
  else if ( LS==3 && MPS==3 )
  {
    if     (uPDG > 7)  // Annihilate pLDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihilOrder: pLDiQ, sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<-7 && (-nPDG==uPDG/1000 || -nPDG==(uPDG/100)%10) &&
               (mPDG==(-sPDG)/1000 || mPDG==((-sPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( nPDG<-7 && (-sPDG==uPDG/1000 || -sPDG==(uPDG/100)%10) &&
               (mPDG==(-nPDG)/1000 || mPDG==((-nPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 33 pLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if(mPDG > 7)  // Annihilate pRDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihilOrder: pRDiQ, sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<-7 && (-nPDG==mPDG/1000 || -nPDG==(mPDG/100)%10) &&
               (uPDG==(-sPDG)/1000 || uPDG==((-sPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( nPDG<-7 && (-sPDG==mPDG/1000 || -sPDG==(mPDG/100)%10) &&
               (uPDG==(-nPDG)/1000 || uPDG==((-nPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 33 pRDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if(sPDG > 7)  // Annihilate cLDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihilOrder: cLDiQ, uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<-7 && (-mPDG==sPDG/1000 || -mPDG==(sPDG/100)%10) &&
               (nPDG==(-uPDG)/1000 || nPDG==((-uPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( mPDG<-7 && (-uPDG==sPDG/1000 || -uPDG==(sPDG/100)%10) &&
               (nPDG==(-mPDG)/1000 || nPDG==((-mPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 33 cLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if(nPDG > 7)  // Annihilate cRDiQ
    {
#ifdef debug
      G4cout<<"G4QFragmentation::AnnihilOrder: cRDiQ, uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<-7 && (-mPDG==nPDG/1000 || -mPDG==(nPDG/100)%10) &&
               (nPDG==(-uPDG)/1000 || nPDG==((-uPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( mPDG<-7 && (-uPDG==nPDG/1000 || -sPDG==(nPDG/100)%10) &&
               (sPDG==(-mPDG)/1000 || sPDG==((-mPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QFragmentation::AnnihilationOrder: Wrong 33 cRDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
#ifdef debug
    else G4cout<<"-Warning-G4QFragmentation::AnnihilOrder: Wrong 33 StringFusion"<<G4endl;
    G4cout<<"G4QFragmentation::AnnihilationOrder: Ord="<<Ord<<",sPDG="<<sPDG<<",nPDG="
          <<nPDG<<", uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
  }
  return Ord;
}

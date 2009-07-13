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
// $Id: G4QFragmentation.cc,v 1.13 2009-07-13 08:59:22 mkossov Exp $
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
//#define ppdebug

#include "G4QFragmentation.hh"

// Promoting model parameters from local variables class properties @@(? M.K.)

// Definition of static parameters
G4int    G4QFragmentation::nCutMax=7; 
G4double G4QFragmentation::ThresholdParameter = 0.*MeV;// Make it appliedFrom E=0.
G4double G4QFragmentation::QGSMThershold=0.*GeV;            // for E=0, (was 3.*GeV) (?)
G4double G4QFragmentation::theNucleonRadius=1.5*fermi;      // M.K.?
// Parameters of diffractional fragmentation
G4double G4QFragmentation::widthOfPtSquare=-0.72*GeV*GeV; // pt -width2 forStringExcitation
G4double G4QFragmentation::minExtraMass=250.*MeV;// minimum excitation mass 
G4double G4QFragmentation::minmass=250.*MeV;     // mean pion transverse mass for Xmin 

G4QFragmentation::G4QFragmentation(){}

G4QFragmentation::~G4QFragmentation() {}

G4QHadronVector* G4QFragmentation::Scatter(const G4QNucleus &aNucleus,
                                           const G4QHadron  &aPrimary)
{ // This is the main member function, which returns the resulting vector of hadrons
  static const G4double  mProt = G4Proton::Proton()->GetPDGMass(); // Mass of proton
  G4QStringVector* strings=0;
  G4QHadron aProjectile = aPrimary;       // As a primary is const
  G4LorentzRotation toZ;                  // Lorentz Transformation to the projectileSystem
  G4LorentzVector proj4M=aProjectile.Get4Momentum(); // Projectile 4-momentum in LS
#ifdef debug
  G4double tgMass=aNucleus.GetGSMass();   // Ground state mass of the nucleus
  G4cout<<"G4QFragmentation::Scatter:Called,p4M="<<proj4M<<",A="<<aNucleus<<tgMass<<G4endl;
#endif
  G4int tZ=aNucleus.GetZ();
  G4int tN=aNucleus.GetN();
  G4int tPDG=90000000+tZ*1000+tN;
  toZ.rotateZ(-1*proj4M.phi());
  toZ.rotateY(-1*proj4M.theta());
  G4LorentzVector zProj4M=toZ*proj4M;     // Proj 4-momentum in LS rotated to Z axis
  aProjectile.Set4Momentum(zProj4M);      // Now the projectile moves along Z exis
#ifdef edebug
  G4int totChg=aProjectile.GetCharge()+tZ;// Charge of the progectile+target for the CHECK
  G4LorentzVector tgLS4M(0.,0.,0.,tgMass);// Target 4-momentum in LS
  G4LorentzVector totLS4M=proj4M+tgLS4M;  // Total 4-momentum in LS
  G4LorentzVector totZLS4M=zProj4M+tgLS4M;// Total 4-momentum in LS with momentum along Z
  G4cout<<"-EMC-G4QFragmentation::Scatter: tLS4M="<<totLS4M<<", tZLS4M="<<totZLS4M<<G4endl;
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
      G4cerr<<"***G4QFragmentation::Scatter: "<<attempts<<" to create a string ( > max="
            <<maxAttempts<<") --> try to increase maxAttempts"<<G4endl;
      G4Exception("G4QFragmentation::Scatter:","72",FatalException,"StringCreation");
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: The attempt # "<<attempts<<G4endl;
#endif
    theNucleus->InitByPDG(tPDG);          // Reinit the Nucleus for the new Attempt
#ifdef pdebug
    G4cout<<"G4QFragmentation::Scatter: Nucl4Mom="<<theNucleus->Get4Momentum()<<G4endl;
#endif
    theNucleus->Init3D();                 // 3D-initialisation(nucleons) of theNucleusClone
#ifdef edebug
    G4LorentzVector sum1=theNucleus->GetNucleons4Momentum();// Sum ofNucleons4M inRotatedLS
    G4cout<<"-EMC-G4QFragmentation::Scatter: attempt#"<<attempts<<", Nuc4M="<<sum1<<G4endl;
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
    G4cout<<"G4QFragmentation::Scatter: Projectile4M="<<proj4M<<", Vz="<<vz<<", nC="
          <<nCons<<", pE="<<e_per_projectile<<G4endl;
#endif
    theCurrentVelocity.setZ(vz);                        // CM (w/r to one nucleon) velosity
    if(theNucleus) theNucleus->DoLorentzBoost(-theCurrentVelocity);// BoostTgNucleus to"CM"
#ifdef edebug
    G4LorentzVector sum2=theNucleus->GetNucleons4Momentum();// Sum ofNucleons4M inRotatedCM
    G4cout<<"-EMC-G4QFragmentation::Scatter: AfterBoost, v="<<vz<<", Nuc4M="<<sum2<<G4endl;
#endif
    G4LorentzVector cmProjMom = proj4M;                 // Copy the original proj4M in LS
    cmProjMom.boost(-theCurrentVelocity);               // Bring the LS proj4Mom to "CM"
    G4double eKin = cmProjMom.e()-projM;                // Primary kinetic energy (MeV!)
    G4double maxCuts=1.5*eKin;                          // Important only for LowEnergies
#ifdef pdebug
    G4cout<<"G4QFragmentation::Scatter: Proj4MInCM="<<cmProjMom<<", pPDG="<<pPDG<<G4endl;
#endif
    //
    // >>>>>>>>>> Find collisions meeting collision conditions
    //
    G4QHadron* cmProjectile = new G4QHadron(pPDG,cmProjMom); // HipCopy of the CMProjectile
    G4QPomeron theProbability(pPDG);                    // the PDG must be a data member
    G4double outerRadius = theNucleus->GetOuterRadius();// Get the nucleus frontiers
#ifdef pdebug
    G4cout<<"G4QFragmentation::Scatter: OuterRad="<<outerRadius<<",mCut="<<maxCuts<<G4endl;
#endif
    // Check the reaction threshold 
    theNucleus->StartLoop();                            // Prepare Loop ovdefr nucleons
    G4QHadron* pNucleon = theNucleus->GetNextNucleon(); // Get the next nucleon to try
    G4double s = (cmProjMom + pNucleon->Get4Momentum()).mag2(); // Squared CM Energy
    G4double ThresholdMass = projM + pNucleon->GetMass(); // @@ Nucleon can be virtual...?
    ModelMode = SOFT;                                   // NOT-Diffractive hadronization
    if (s < sqr(ThresholdMass + ThresholdParameter))    // At ThP=0 is impossible(virtNucl)
    {
      G4cerr<<"***G4QFragmentation::Scatter: ThrM="<<ThresholdMass<<" + ThrPa="
            <<ThresholdParameter<<" = "<<ThresholdMass+ThresholdParameter<<" > sqrt(s)="
            <<std::sqrt(s)<<G4endl;
      G4Exception("G4QFragmentation::Scatter:","72",FatalException,"LowEnergy");
    }
    if (s < sqr(ThresholdMass + QGSMThershold))         // --> Only diffractive interaction
    {
#ifdef debug
      G4cout<<"G4QFragmentation::Scatter: *OnlyDiffraction* ThrM="<<ThresholdMass
            <<" + ThrQGS="<<QGSMThershold<<" = "<<ThresholdMass+QGSMThershold<<" >sqrt(s)="
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
      G4cout<<"G4QFragmentation::Scatter: *** Enter the interaction LOOP ***"<<G4endl;
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
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter:LOOP overNucleons, totCuts="<<totalCuts<<G4endl;
#endif
        // Needs to be moved to Probability class @@@
        G4double s = (cmProjMom + pNucleon->Get4Momentum()).mag2();
        G4double Distance2 = sqr(impactX - pNucleon->GetPosition().x()) +
                             sqr(impactY - pNucleon->GetPosition().y());
        G4double Probability = theProbability.GetInelasticProbability(s, Distance2);// INEL
        // test for inelastic collision
        G4double rndNumber = G4UniformRand();           // For the printing purpose
        // ModelMode = DIFFRACTIVE;
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter: NLOOP prob="<<Probability<<", rndm="<<rndNumber
              <<", d="<<std::sqrt(Distance2)<<G4endl;
#endif
        if (Probability > rndNumber) // Inelastic (diffractive or soft) interaction (JOB)
        {
          G4QHadron* aTarget = new G4QHadron(*pNucleon);// Copy selected nucleon for String
#ifdef edebug
          G4cout<<">EMC>G4QFragmentation::Scatter: Target Nucleon is filled, 4M/PDG="
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
            totalCuts += 1;                             // UpdateOfNucleons in't necessary
#ifdef debug
            G4cout<<"G4QFragmentation::Scatter:NLOOP DiffInteract, tC="<<totalCuts<<G4endl;
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
            G4cout<<"G4QFragmentation::Scatter:NLOOP-Soft Chosen nCut="<<nCut<<G4endl;
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
            G4cout<<"G4QFragmentation::Scatter:NLOOP SoftInteract, tC="<<totalCuts<<G4endl;
#endif
            impactUsed=Distance2;
          }
        }
      }
      // When nucleon count is incremented, the LOOP stops, so nucleonCount==1 always!
#ifdef debug
      G4cout<<"G4QFragmentation::Scatter: NUCLEONCOUNT="<<nucleonCount<<G4endl;
#endif
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: CUTDEBUG="<<totalCuts<<", ImpactParam="
          <<impactUsed<<", #ofInter="<<theInteractions.size()<<G4endl;
#endif
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
      G4cout<<"-EMC-G4QFragmentation::Scatter: Interaction#"<<i<<",dP4M="
            <<projH->Get4Momentum()-pSP<<",dT4M="<<targH->Get4Momentum()-tSP<<G4endl;
#endif
    }  
    // 
    // >>>>>>>> make soft collisions (ordering is vital)
    //
    G4QInteractionVector::iterator i;
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: Creation ofSoftCollisionPartonPair STARTS"<<G4endl;
#endif
    for(i = theInteractions.begin(); i != theInteractions.end(); i++)   
    {
      G4QInteraction* anIniteraction = *i;
      G4QPartonPair*  aPair=0;
      G4int nSoftCollisions = anIniteraction->GetNumberOfSoftCollisions();
#ifdef debug
      G4cout<<"G4QFragmentation::Scatter: #0f SOFT collisions ="<<nSoftCollisions<<G4endl;
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
          G4cout<<"G4QFragmentation::Scatter: SOFT, 2 parton pairs're filled"<<G4endl;
#endif
        }  
        delete *i;
        i=theInteractions.erase(i);       // Soft interactions are converted & erased
        i--;
      }
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: -> Parton pairs for SOFT strings are made"<<G4endl;
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
      G4cout<<"G4QFragmentation::Scatter: CreationOfDiffractivePartonPairs, i="<<i<<G4endl;
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
        G4cout<<"G4QFragmentation::Scatter: proj Diffractive PartonPair is filled"<<G4endl;
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
        G4cout<<"G4QFragmentation::Scatter: targ Diffractive PartonPair is filled"<<G4endl;
#endif
      }
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter:DiffractivePartonPairs are created"<<G4endl;
#endif  
    //
    // >>>>>>>>>>>>>> clean-up  Interactions and cmProjectile, if necessary
    //
    std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
    theInteractions.clear();
    delete cmProjectile;
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: Temporary objects are cleaned up"<<G4endl;
#endif  
    // This function prepares theBoost for transformation of secondaries to LS (-ProjRot!)
    theNucleus->DoLorentzBoost(theCurrentVelocity);// Boost theResidualNucleus to RotatedLS
    // @@ Nucleus isn't completely in LS, it needs the toZ (-ProjRot) rotation to consE/M
#ifdef pdebug
    G4cout<<"G4QFragmentation::Scatter: >>>>>> Strings are created "<<G4endl;
#endif
    G4QPartonPair* aPair;
    strings = new G4QStringVector;
    G4QString* aString=0;
    while(thePartonPairs.size()) // @@ At present noDifference in stringBuild (? M.K.)
    {
      aPair = thePartonPairs.back();           // Get the parton pair
      thePartonPairs.pop_back();               // Clean up the parton pair in the vector
#ifdef debug
      G4cout<<"G4QFragmentation::Scatter: StringType="<<aPair->GetCollisionType()<<G4endl;
#endif
      if (aPair->GetCollisionType() == G4QPartonPair::DIFFRACTIVE)
      {
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter:Prepared for creationOfDiffExcitation"<<G4endl;
#endif
        aString = BuildString(aPair);           // Just to shorten the new
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter:DifString4M="<<aString->Get4Momentum()<<G4endl;
#endif
      }
      else
      {
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter: Prepared for creationOfTheSoftString"<<G4endl;
#endif
        aString = BuildString(aPair);           // Just to shorten the new
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter:SoftString4M="<<aString->Get4Momentum()<<G4endl;
#endif
      }
      aString->Boost(theCurrentVelocity);       // ! Strings are moved to ZLS when pushed !
      strings->push_back(aString);
      delete aPair;
    } // End of the String Creation LOOP
#ifdef edebug
    G4LorentzVector sum=theNucleus->Get4Momentum();// Nucleus 4Mom in rotatedLS
    G4int rChg=totChg-theNucleus->GetZ();
    G4int nStrings=strings->size();
    G4cout<<"-EMC-G4QFragmentation::Scatter:#ofString="<<nStrings<<",tNuc4M="<<sum<<G4endl;
    for(G4int i=0; i<nStrings; i++)
    {
      G4LorentzVector strI4M=(*strings)[i]->Get4Momentum();
      sum+=strI4M;
      G4int sChg=(*strings)[i]->GetCharge();
      G4int LPDG=(*strings)[i]->GetLeftParton()->GetPDGCode();
      G4int RPDG=(*strings)[i]->GetRightParton()->GetPDGCode();
      G4QContent LQC=(*strings)[i]->GetLeftParton()->GetQC();
      G4QContent RQC=(*strings)[i]->GetRightParton()->GetQC();
      rChg-=sChg;
      G4cout<<"-EMC-G4QFragmentation::Scatter: String#"<<i<<", 4M="<<strI4M<<",LPDG="<<LPDG
            <<LQC<<",RPDG="<<RPDG<<RQC<<",Ch="<<sChg<<G4endl;
    }
    G4cout<<"-EMC-G4QFragmentation::Scatter: res4M="<<sum-totZLS4M<<", rCh="<<rChg<<G4endl;
#endif
  } // End of the LOOP of the strings creation
#ifdef debug
  G4cout<<"G4QFragmentation::Scatter: BeforeRotation #0fStrings="<<strings->size()<<G4endl;
#endif
  //
  // ------------------ At this point the strings must be created -----------------
  //
  for(unsigned astring=0; astring < strings->size(); astring++)
          (*strings)[astring]->LorentzRotate(toLab); // Recove Z-direction in LS ("LS"->LS)
  // Now everything is in LS system
#ifdef edebug
  // @@ Should be tested for the primary projectile not along Z !
  G4LorentzVector sum=theNucleus->Get4Momentum(); // Nucleus 4Mom in LS
  G4int rChg=totChg-theNucleus->GetZ();
  G4int nStrings=strings->size();
  G4cout<<"-EMCLS-G4QFragmentation::Scatter:#ofS="<<nStrings<<",tNuc4M(E=M)="<<sum<<G4endl;
  for(G4int i=0; i<nStrings; i++)
  {
    G4LorentzVector strI4M=(*strings)[i]->Get4Momentum();
    sum+=strI4M;
    G4int sChg=(*strings)[i]->GetCharge();
    rChg-=sChg;
    G4cout<<"-EMCLS-G4QFragmentation::Scatter:S#"<<i<<",4M="<<strI4M<<",Ch="<<sChg<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFragmentation::Scatter: res4M="<<sum-totLS4M<<",resCh="<<rChg<<G4endl;
#endif
  //
  // ------------ At this point the strings are fragmenting to hadrons in LS -------------
  //
  G4QHadronVector* theResult = new G4QHadronVector;
  G4LorentzVector KTsum(0.,0.,0.,0.);
  G4int nOfStr=strings->size();
#ifdef debug
  G4cout<<"G4QFragmentation::Scatter: BeforeFragmentation, #OfStr="<<nOfStr<<G4endl;
#endif
  for (G4int astring=0; astring < nOfStr; astring++)
  {
    G4QString* curString=(*strings)[astring];
    if(!curString->GetDirection()) continue;
#ifdef edebug
    G4int curStrChg = curString->GetCharge();
#endif
    G4LorentzVector curString4M = curString->Get4Momentum();
    KTsum+= curString4M;                       // @@ ?
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: String#"<<astring<<",s4M/m="<<curString4M
          <<curString4M.m()<<", LPart="<<curString->GetLeftParton()->GetPDGCode()
          <<", RPart="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
#endif
    if(!(KTsum.e()<1.) && !(KTsum.e()>-1.))    // NAN check
    {
      G4cerr<<"***G4QFragmentation::Scatter: KTsum="<<KTsum<<G4endl;
      G4Exception("G4QFragmentation::Scatter:","72",FatalException,"NANin3Vector");
    }
    G4QHadronVector* theHadrons = 0;           // Prototype of theStringFragmentationOUTPUT
    theHadrons=curString->FragmentString(true);// !! Fragmentin the String !!
    if (!theHadrons)                           // The string can not be fragmented
    {
      // First try to correct the diQ-antiDiQ strings, converting them to Q-antiQ
      G4QParton* cLeft=curString->GetLeftParton();
      G4QParton* cRight=curString->GetRightParton();
      G4int sPDG=cLeft->GetPDGCode();
      G4int nPDG=cRight->GetPDGCode();
      G4int LT=cLeft->GetType();
      G4int RT=cRight->GetType();
      G4int fusionDONE=0; // String fusion didn't happend (1=Fuse L+L/R+R, -1=Fuse L+R/R+L)
      if(LT==2 && RT==2)
      {
        G4bool LisADQ=false;                   // Left parton is not an antiDiQ (default)
        G4int qPDG=sPDG;
        if(qPDG<-99)
        {
          LisADQ=true;
          qPDG=(-qPDG)/100;
        }
        else qPDG/=100;
        G4bool RisADQ=false;                   // Left parton is not an antiDiQ (default)
        G4int dPDG=nPDG;
        if(dPDG<-99)
        {
          RisADQ=true;
          dPDG=(-dPDG)/100;
        }
        else dPDG/=100;
        G4int L1=qPDG/10;
        G4int L2=qPDG%10;
        G4int R1=dPDG/10;
        G4int R2=dPDG%10;
        if(L1==R1 || L1==R2 || L2==R1 || L2==R2) // Annihilation condition
        {
          sPDG=L1;                             // Default for L2==R2
          nPDG=R1;                             // Default for L2==R2
          if     (L1==R1)
          {
            sPDG=L2;
            nPDG=R2;
          }
          else if(L1==R2)
          {
            sPDG=L2;
            nPDG=R1;
          }
          else if(L2==R1)
          {
            sPDG=L1;
            nPDG=R2;
          }
          if(LisADQ < 0) sPDG=-sPDG;
          if(RisADQ < 0) nPDG=-nPDG;
          LT=1;
          RT=1;
          cLeft->SetPDGCode(sPDG);             // Reset the left quark
          cRight->SetPDGCode(nPDG);            // Reset the right quark
          theHadrons=curString->FragmentString(true);//!! Try to fragment the new String !!
        }
      }
      if (!theHadrons)                         // The string can not be fragmented
      {
        G4int next=astring+1;                  // The next string position
        if(next < nOfStr)                      // The string is not the last string
        {
          G4int fustr=0;                       // The found partner index (never can be 0)
          G4int swap=0;                        // working interger for swapping parton PDG
          G4double Vmin=DBL_MAX;               // Prototype of the found Velocity Distance 
          if(nPDG<-99 || (nPDG>0&&nPDG<7) || sPDG>99 || (sPDG<0 && sPDG>-7))
          {
            swap=sPDG;
            sPDG=nPDG;
            nPDG=swap;
          }
          G4int dPDG=nPDG;
          if(dPDG>99) dPDG/=100;
          G4int qPDG=sPDG;
          if(qPDG<-99) qPDG=(-qPDG)/100;
#ifdef debug
          G4cout<<"^^G4QFragmentation::Scatter:No Hadrons, q="<<qPDG<<", a="<<dPDG<<G4endl;
#endif
          G4ThreeVector curV=curString4M.vect()/curString4M.e();
          for (G4int restr=next; restr < nOfStr; restr++)
          {
            G4QString* reString=(*strings)[restr];
            G4QParton* Left=reString->GetLeftParton();
            G4QParton* Right=reString->GetRightParton();
            G4int uPDG=Left->GetPDGCode();
            G4int mPDG=Right->GetPDGCode();
            G4int PLT =Left->GetType();
            G4int PRT =Right->GetType();
            if(mPDG<-99 || (mPDG>0&&mPDG<7) || uPDG>99 || (uPDG<0 && uPDG>-7))
            {
              swap=uPDG;
              uPDG=mPDG;
              mPDG=swap;
            }
            G4int aPDG=mPDG;
            if(aPDG>99) aPDG/=100;
            G4int rPDG=uPDG;
            if(rPDG<-99) rPDG=(-rPDG)/100;
            // Check that the strings can fuse (no anti-diquarks assumed)
            if( (LT==1 && RT==1 && PLT==1 && PRT==1)           ||
                (aPDG>0 && (dPDG==-aPDG/10 || dPDG==-aPDG%10)) ||
                (dPDG>0 && (aPDG==-dPDG/10 || aPDG==-dPDG%10)) ||
                (rPDG>7 && (qPDG== rPDG/10 || qPDG== rPDG%10)) ||
                (qPDG>7 && (rPDG== qPDG/10 || rPDG== qPDG%10)) )
            {
              G4LorentzVector reString4M = reString->Get4Momentum();
              G4ThreeVector reV = reString4M.vect()/reString4M.e();
              G4double dV=(curV-reV).mag2();   // Squared difference between relVelocities
#ifdef debug
              G4cout<<"G4QFragmentation::Scatter:Cand#"<<restr<<", q="<<rPDG<<", a="<<aPDG
                    <<",dV="<<dV<<G4endl;
#endif
              if(dV < Vmin)
              {
                Vmin=dV;
                fustr=restr;
              }
            }
          }
          if(fustr)                            // The partner was found -> fuse strings
          {
            G4QString* fuString=(*strings)[fustr];
            G4QParton* Left=fuString->GetLeftParton();
            G4QParton* Right=fuString->GetRightParton();
            G4int uPDG=Left->GetPDGCode();
            G4int mPDG=Right->GetPDGCode();
            G4int PLT =Left->GetType();
            G4int PRT =Right->GetType();
            if(mPDG<-99 || (mPDG>0&&mPDG<7) || uPDG>99 || (uPDG<0 && uPDG>-7))
            {
              swap=uPDG;
              uPDG=mPDG;
              mPDG=swap;
            }
            G4int rPDG=uPDG;
            if(rPDG<-99) rPDG=(-rPDG)/100;
            G4int aPDG=mPDG;
            if(aPDG>99) aPDG/=100;
            // Check that the strings can fuse (no anti-diquarks assumed)
            G4LorentzVector resL4M(0.,0.,0.,0.);
            G4LorentzVector resR4M(0.,0.,0.,0.);
            G4LorentzVector L4M=cLeft->Get4Momentum();
            G4LorentzVector R4M=cRight->Get4Momentum();
            G4LorentzVector PL4M=Left->Get4Momentum();
            G4LorentzVector PR4M=Right->Get4Momentum();
            if      (LT==1 && RT==1 && PLT==1 && PRT==1) // Fuse 2 Q-aQ strings to DiQ-aDiQ
            {
              if     ( (uPDG>0 && sPDG>0 && mPDG<0 && nPDG<0) || 
                       (uPDG<0 && sPDG<0 && mPDG>0 && nPDG>0) ) fusionDONE=1; // LL/RR
              else if( (uPDG>0 && nPDG>0 && mPDG<0 && sPDG<0) || 
                       (uPDG<0 && nPDG<0 && mPDG>0 && sPDG>0) ) fusionDONE=-1; // LR/RL
              else G4cerr<<"-Warning-G4QFragmentation::Scatter: Wrong QQ-fusion, L="<<uPDG
                         <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
            }
            else if (aPDG > 7 && (dPDG == -aPDG/10 || dPDG == -aPDG%10))  // Fuse cAQ+pDiQ
            {
              if     ( (uPDG>7 && sPDG>-7 && sPDG<0) ||
                       (mPDG>7 && nPDG>-7 && nPDG<0) ) fusionDONE=1; // LL/RR
              else if( (uPDG>7 && nPDG>-7 && nPDG<0) ||
                       (mPDG>7 && sPDG>-7 && sPDG<0) ) fusionDONE=-1; // LR/RL
              else G4cerr<<"-Warning-G4QFragmentation::Scatter: Wrong cAQ+pDiQ, L="<<uPDG
                         <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
            }
            else if (dPDG > 7 && (aPDG == -dPDG/10 || aPDG == -dPDG%10))  // Fuse cDiQ+pAQ
            {
              if     ( (sPDG>7 && uPDG>-7 && uPDG<0) ||
                       (nPDG>7 && mPDG>-7 && mPDG<0) ) fusionDONE=1; // LL/RR
              else if( (sPDG>7 && mPDG>-7 && mPDG<0) ||
                       (nPDG>7 && uPDG>-7 && uPDG<0) ) fusionDONE=-1; // LR/RL
              else G4cerr<<"-Warning-G4QFragmentation::Scatter: Wrong cDiQ+pAQ, L="<<uPDG
                         <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
            }
            else if (rPDG > 7 && (qPDG == rPDG/10 || qPDG == rPDG%10))    // Fuse cQ+pADiQ
            {
              if     ( (uPDG<-7 && sPDG<7 && sPDG>0) ||
                       (mPDG>-7 && nPDG<7 && nPDG>0) ) fusionDONE=1; // LL/RR
              else if( (uPDG<-7 && nPDG<7 && nPDG>0) ||
                       (mPDG<-7 && sPDG<7 && sPDG>0) ) fusionDONE=-1; // LR/RL
              else G4cerr<<"-Warning-G4QFragmentation::Scatter: Wrong cQ+pADiQ, L="<<uPDG
                         <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
            }
            else if (qPDG > 7 && (rPDG == qPDG/10 || rPDG == qPDG%10))    // Fuse cADiQ+pQ
            {
              if     ( (sPDG<-7 && uPDG<7 && uPDG>0) ||
                       (nPDG>-7 && mPDG<7 && mPDG>0) ) fusionDONE=1; // LL/RR
              else if( (sPDG<-7 && mPDG<7 && mPDG>0) ||
                       (nPDG<-7 && uPDG<7 && uPDG>0) ) fusionDONE=-1; // LR/RL
              else G4cerr<<"-Warning-G4QFragmentation::Scatter: Wrong cADiQ+pQ, L="<<uPDG
                         <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
            }
#ifdef debug
            else G4cout<<"-Warning-G4QFragmentation::Scatter: WrongStringFusion"<<G4endl;
#endif
            if     (fusionDONE>0)
            {
              Left->SetPDGCode(SumPartonPDG(uPDG, sPDG));
              Left->Set4Momentum(L4M+PL4M);
              Right->SetPDGCode(SumPartonPDG(mPDG, nPDG));
              Right->Set4Momentum(R4M+PR4M);
#ifdef debug
              G4cout<<"G4QFragmentation::Scatter:LL/RR s4M="<<fuString->Get4Momentum()
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
              G4cout<<"G4QFragmentation::Scatter:LR/RL s4M="<<fuString->Get4Momentum()
                    <<",S="<<L4M+PL4M+R4M+PR4M<<", L="<<Left->Get4Momentum()<<", R="
                    <<Right->Get4Momentum()<<G4endl;
#endif
            }
#ifdef debug
            else G4cout<<"G4QFragmentation::Scatter: WrongStringFusion"<<G4endl;
#endif
#ifdef debug
            if(fusionDONE) G4cout<<"####G4QFragmentation::Scatter: String #"<<astring
                                 <<" is fused with #"<<fustr<<", new4Mom="
                                 <<fuString->Get4Momentum()<<G4endl;
#endif
#ifdef edebug
            G4cout<<"#EMC#G4QFragmentation::Scatter:StringFused,F="<<fusionDONE<<",L="<<L4M
                  <<",R="<<R4M<<",pL="<<PL4M<<",pR="<<PR4M<<",nL="<<Left->Get4Momentum()
                  <<",nR="<<Right->Get4Momentum()<<",S="<<fuString->Get4Momentum()<<G4endl;
#endif
          }
          else
          {
            G4cerr<<"**G4QFragmentation::Scatter:*NoPart*M="<<curString->Get4Momentum().m()
                  <<", LPDG="<<curString->GetLeftParton()->GetPDGCode()
                  <<", RPDG="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
            // @@ Temporary exception for the future solution
            //G4Exception("G4QFragmentation::Scatter:","72",FatalException,"NoStrgsFused");
          }
        }
        else                                     // The string is the last string
        {
          G4cerr<<"***4QFragmentation::Scatter:*TheLast* M="<<curString->Get4Momentum().m()
                <<", LPDG="<<curString->GetLeftParton()->GetPDGCode()
                <<", RPDG="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
          // @@ Temporary exception for the future solution
          //G4Exception("G4QFragmentation::Scatter:","72",FatalException,"LastStrgFusion");
        }
        if(fusionDONE) continue;
      } // End of IF(NO_Hadrons)
      if(fusionDONE)
      {
        curString->KillString();
#ifdef debug
        G4cout<<"####G4QFragmentation::Scatter: String #"<<astring<<" is killed !"<<G4endl;
#endif
      }
    }
    G4Quasmon tmpQ;                                 // @@ an issue of Q to decay resonances
    for(unsigned aTrack=0; aTrack<theHadrons->size(); aTrack++)
    {
      G4QHadron* curHadron=(*theHadrons)[aTrack];
      G4int hPDG=curHadron->GetPDGCode();
#ifdef edebug
      G4LorentzVector curH4M=curHadron->Get4Momentum();
      G4int           curHCh=curHadron->GetCharge();
#endif
#ifdef debug
						G4cout<<">>>>>>>>G4QFragmentation::Scatter:S#"<<astring<<",H#"<<aTrack<<",PDG="<<hPDG
            <<",4M="<<curHadron->Get4Momentum()<<G4endl;
#endif
      if(std::abs(hPDG)%10 > 2)
      {
        G4QHadronVector* tmpQHadVec=tmpQ.DecayQHadron(curHadron);
#ifdef pdebug
        G4cout<<"G4QFragmentation::Scatter:-DECAY'S DONE-,nH="<<tmpQHadVec->size()<<G4endl;
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
          curH4M-=p4M;
          curString4M-=p4M;
          curStrChg-=Chg;
          curHCh-=Chg;
          G4cout<<"-EMC->>>>G4QFragmentation::Scatter:Str*Filled, 4M="<<p4M<<", PDG="<<PDG
                <<", Chg="<<Chg<<G4endl;
#endif
        }
#ifdef edebug
        G4cout<<"-EMC-.G4QFragmentation::Scatter:Dec,r4M="<<curH4M<<",rC="<<curHCh<<G4endl;
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
        curStrChg-=curCh;
        G4cout<<"-EMC->>>>>>G4QFragmentation::Scatter: curH filled 4M="<<curH4M<<",PDG="
              <<curHadron->GetPDGCode()<<", Chg="<<curCh<<G4endl;
#endif
      }
    }
    G4LorentzVector r4M=theNucleus->GetNucleons4Momentum();// Sum of N4M's in RotatedLS=LS
    // clean up (the issues are filled to theResult)
    delete theHadrons;
#ifdef edebug
    G4cout<<"-EMC-.........G4QFragmentation::Scatter: StringDecay CHECK, r4M="<<curString4M
          <<", rChg="<<curStrChg<<G4endl;
#endif
  }
  G4LorentzVector r4M=theNucleus->Get4Momentum(); // Nucleus 4-momentum in LS
  G4int rPDG=theNucleus->GetPDG();
  G4QHadron* resNuc = new G4QHadron(rPDG,r4M);
  theResult->push_back(resNuc);                          // Fill the residual nucleus
#ifdef edebug
  // @@ Should be tested for the primary projectile not along Z !
  G4LorentzVector s4M(0.,0.,0.,0.); // Sum of the Result in LS
  G4int rCh=totChg;
  G4int nHadr=theResult->size();
  G4cout<<"-EMCLS-G4QFragmentation::Scatter:#ofHadr="<<nHadr<<",rN="<<r4M.m()<<"="
        <<G4QNucleus(rPDG).GetGSMass()<<G4endl;
  for(G4int i=0; i<nHadr; i++)
  {
    G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
    s4M+=hI4M;
    G4int hChg=(*theResult)[i]->GetCharge();
    rCh-=hChg;
    G4cout<<"-EMCLS-G4QFragmentation::Scatter: H#"<<i<<",4M="<<hI4M<<",Chg="<<hChg<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFragmentation::Scatter: LS r4M="<<s4M-totLS4M<<", rCh="<<rCh<<G4endl;
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
  G4double Mprojectile=projectile->GetMass() + minExtraMass;
  G4double Mprojectile2=Mprojectile*Mprojectile;
  G4LorentzVector Ptarget=target->Get4Momentum();
  G4double Mtarget=target->GetMass() + minExtraMass;
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
    G4double Xmin = minmass/(Pprojectile.e() + Ptarget.e());
    G4double Xmax=1.;
    G4double Xplus =ChooseX(Xmin,Xmax);
    G4double Xminus=ChooseX(Xmin,Xmax);
#ifdef debug
    G4cout<<"G4QFragment::ExciteDiffParticip: X-plus="<<Xplus<<",X-minus="<<Xminus<<G4endl;
#endif
    G4double pt2=G4ThreeVector(Qmomentum.vect()).mag2();
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
  G4double Mprojectile=projectile->GetMass() + minExtraMass;
  G4double Mprojectile2=Mprojectile*Mprojectile;
  G4LorentzVector Ptarget=target->Get4Momentum();
  G4double Mtarget=target->GetMass() + minExtraMass;
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
    G4double Xmin = minmass/(Pprojectile.e() + Ptarget.e());
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

void G4QFragmentation::SetParameters(G4int nCM, G4double thresh, G4double QGSMth,
                           G4double radNuc, G4double SigPt, G4double extraM, G4double minM)
{//  =============================================================================
  nCutMax            = nCM;            // max number of pomeron cuts
  ThresholdParameter = thresh;         // internal threshold
  QGSMThershold      = QGSMth;         // QGSM threshold
  theNucleonRadius   = radNuc;         // effective radius of the nucleon inside Nucleus
  widthOfPtSquare    = -2*SigPt*SigPt; // width^2 of pt for string excitation
  minExtraMass       = extraM;         // minimum excitation mass 
  minmass            = minM;           // mean pion transverse mass; used for Xmin 
}

G4double G4QFragmentation::ChooseX(G4double Xmin, G4double Xmax) const
{
// choose an x between Xmin and Xmax with P(x) ~ 1/x
//  to be improved...
  G4double range=Xmax-Xmin;
  if( Xmin<= 0. || range <=0.) 
  {
    G4cerr<<"***G4QFragmentation::ChooseX: Xmin="<<Xmin<<", Xmax="<<Xmax<< G4endl;
    G4Exception("G4QFragmentation::ChooseX:","72",FatalException,"BadXRange");
  }
  G4double x;
  do {x=Xmin+G4UniformRand()*range;} while ( Xmin/x < G4UniformRand() );
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
      G4cerr<<"***4QFragmentation::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
      G4Exception("G4QFragmentation::SumPartonPDG:","72",FatalException,"AQ&DiQ notMatch");
    }
  }
  else
  {
    G4cerr<<"***4QFragmentation::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
    G4Exception("G4QFragmentation::SumPartonPDG:","72",FatalException,"Can'tSumUpPartons");
  }
  return 0;
} // End of SumPartonPDG

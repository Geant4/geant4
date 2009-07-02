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
// $Id: G4QFragmentation.cc,v 1.9 2009-07-02 07:17:09 mkossov Exp $
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

#include "globals.hh"
#include "G4QFragmentation.hh"
#include "G4LorentzVector.hh"
#include <utility>

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

G4QFragmentation::G4QFragmentation()
{
  // Construct Shortlived particles (needed after the 2006 Particles revolution)
  G4ShortLivedConstructor ShortLived; // @@ CHIPS decays can be used 
  ShortLived.ConstructParticle();     // @@ Not necessary in CHIPS
}

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
  G4cout<<"G4QFragmentation::Scatter: Called, proj4M="<<proj4M<<", Nuc="<<aNucleus<<G4endl;
#endif
  G4int tZ=aNucleus.GetZ();
  G4int tN=aNucleus.GetN();
  G4int tPDG=90000000+tZ*1000+tN;
  toZ.rotateZ(-1*proj4M.phi());
  toZ.rotateY(-1*proj4M.theta());
  G4LorentzVector zProj4M=toZ*proj4M;     // Proj 4-momentum in LS rotated to Z axis
  aProjectile.Set4Momentum(zProj4M);      // Now the projectile moves along Z exis
#ifdef edebug
  G4double tgMass=aNucleus.GetGSMass();
  G4LorentzVector tgLS4M(0.,0.,0.,tgMass);// Target 4-momentum in LS
  G4LorentzVector totLS4M=proj4M+tgLS4M;  // Total 4-momentum in LS
  G4LorentzVector totZLS4M=zProj4M+tgLS4M;// Total 4-momentum in LS with momentum along Z
  G4cout<<"-EMC-G4QFragmentation::Scatter: tLS4M="<<totLS4M<<", tZLS4M="<<totZLS4M<<G4endl;
  // === From nere all consideration is made in the rotated LS frame (proj is along Z) ===
#endif
  G4LorentzRotation toLab(toZ.inverse()); // Lorentz Transfornation back to LS (at the end)
  G4int pPDG=aProjectile.GetPDGCode();    // The PDG code of the projectile
  G4double projM=aProjectile.GetMass();   // Mass of the projectile
  G4QInteractionVector theInteractions;   // A vector of interactions (taken from the Body)
  G4QHadronVector      theTargets;        // Selevcted target-nucleons (taken from theBody)
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
    G4cout<<"G4QFragmentation::Scatter: Proj4MInCM="<<cmProjMom<<", Cut="<<maxCuts<<G4endl;
#endif
    //
    // >>>>>>>>>> Find collisions meeting collision conditions
    //
    G4QHadron* cmProjectile = new G4QHadron(pPDG,cmProjMom); // HipCopy of the CMProjectile
    G4QPomeron theProbability(pPDG);                    // the PDG must be a data member
    G4double outerRadius = theNucleus->GetOuterRadius();// Get the nucleus frontiers
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: OuterRadius = "<<outerRadius<<G4endl;
#endif
    // Check the reaction threshold 
    theNucleus->StartLoop();                            // Prepare Loop ovdefr nucleons
    G4QHadron* pNucleon = theNucleus->GetNextNucleon(); // Get the next nucleon to try
    G4double s = (cmProjMom + pNucleon->Get4Momentum()).mag2(); // Squaed CM Energy
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
      theNucleus->StartLoop();
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
        if (Probability > rndNumber) // Inelastic (diffractive or soft) interaction
        {
          G4QHadron* aTarget = new G4QHadron(*pNucleon);// Copy selected nucleon for String
          // @@ Instead of theTargets update theNucleus (M.K.)
          theTargets.push_back(aTarget);
#ifdef edebug
          G4cout<<">EMC>G4QFragmentation::Scatter: Target Nucleon is filled"<<G4endl;
#endif
          if((theProbability.GetDiffractiveProbability(s,Distance2)/Probability >
              G4UniformRand() && ModelMode==SOFT ) || ModelMode==DIFFRACTIVE)
          { 
            // ---------------->>>> diffractive interaction
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
    for(unsigned i=0; i<theInteractions.size(); i++) theInteractions[i]->SplitHadrons();
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: Parton pairs are builded"<<G4endl;
#endif  
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
          thePartonPairs.push_back(aPair); // A target pair (@@Chrck 4M)
          aPair = new G4QPartonPair(pProjectile->GetNextParton(),
                                    pTarget->GetNextAntiParton(),
                                    G4QPartonPair::SOFT, G4QPartonPair::PROJECTILE);
          thePartonPairs.push_back(aPair); // A projectile pair (@@Chrck 4M)
#ifdef debug
          G4cout<<"G4QFragmentation::Scatter: SOFT, 2 parton pairs're filled"<<G4endl;
#endif
        }  
        delete *i;
        i=theInteractions.erase(i);
        i--;
      }
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: Strings are hadronized"<<G4endl;
#endif  
		  //
    // >>>>>>>>>>>>>>> make the rest as the diffractive fragmentation
    //
    for(unsigned i = 0; i < theInteractions.size(); i++) // Interactions are reduced bySoft
    {
      G4QInteraction* anIniteraction = theInteractions[i];
      G4QHadron* aProjectile = anIniteraction->GetProjectile();
      G4QParton* aParton = aProjectile->GetNextParton();
      G4QPartonPair* aPartonPair;
#ifdef debug
      G4cout<<"G4QFragmentation::Scatter: Creation ofDiffractivePartonPair STARTS"<<G4endl;
#endif
      // the projectile diffraction first
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
      // then the target diffraction
      G4QHadron* aTarget = anIniteraction->GetTarget();
      aParton = aTarget->GetNextParton();
      if (aParton)
      {
        aPartonPair = new G4QPartonPair(aParton, aTarget->GetNextAntiParton(), 
                                        G4QPartonPair::DIFFRACTIVE, G4QPartonPair::TARGET);
        thePartonPairs.push_back(aPartonPair);
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter: proj Diffractive PartonPair is filled"<<G4endl;
#endif
      }
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter:DiffractiveExcitations are made"<<G4endl;
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
    theNucleus->DoLorentzBoost(theCurrentVelocity);// Boost theResidualNucleusToLS
    // @@ Nucleus isn't completely in LS, it needs the toZ (-ProjRot) rotation to consE/M
#ifdef pdebug
    G4cout<<"G4QFragmentation::Scatter: >>> >>> Strings are created "<<G4endl;
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
        aString = BuildString(aPair);           // @@ ?
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter:DifString4M="<<aString->Get4Momentum()<<G4endl;
#endif
      }
      else
      {
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter: Prepared for creationOfTheSoftString"<<G4endl;
#endif
        aString = BuildString(aPair);           // @@ ?
#ifdef debug
        G4cout<<"G4QFragmentation::Scatter:SoftString4M="<<aString->Get4Momentum()<<G4endl;
#endif
      }
      aString->Boost(theCurrentVelocity);
      strings->push_back(aString);
      delete aPair;
    } // End of the String Creation LOOP
#ifdef edebug
    G4LorentzVector sum=theNucleus->GetNucleons4Momentum(); // Sum of N4Mom's in rotated LS
    G4int nStrings=strings->size();
    G4cout<<"-EMC-G4QFragmentation::Scatter:#ofS="<<nStrings<<", tgNuc4M="<<sum<<G4endl;
    for(G4int i=0; i<nStrings; i++)
    {
      G4LorentzVector strI4M=(*strings)[i]->Get4Momentum();
      sum+=strI4M;
      G4cout<<"-EMC-G4QFragmentation::Scatter:Str#"<<i<<",4M="<<strI4M<<",S="<<sum<<G4endl;
    }
    G4int nTgNuc = theTargets.size();
    for(G4int i=0; i<nTgNuc; i++)
    {
      G4LorentzVector tgI4M=theTargets[i]->Get4Momentum();
      sum-=tgI4M;
      G4cout<<"-EMC-G4QFragmentation::Scatter:Tg#"<<i<<", 4M="<<tgI4M<<", S="<<sum<<G4endl;
    }
    G4cout<<"-EMC-G4QFragmentation::Scatter:#ofN="<<nTgNuc<<", R4M="<<sum-totZLS4M<<G4endl;
#endif
  } // End of the LOOP of the strings creation
  //
  // ------------------ At this point the strings must be created -----------------
  //
  G4double stringEnergy=0.;                    // To sum up the total string energy
  for( unsigned astring=0; astring < strings->size(); astring++)
  {
    // Make energy sum and bring strings to LS, models generate them along z
    G4QString* curString=(*strings)[astring];
    stringEnergy += curString->GetLeftParton()->Get4Momentum().t();
    stringEnergy += curString->GetRightParton()->Get4Momentum().t();
    curString->LorentzRotate(toLab); // @@ Remnber about the projectile Rotation
  }
#ifdef debug
  G4cout<<"G4QFragmentation::Scatter: Total energy of strings = "<<stringEnergy<<G4endl;
#endif
  //
  // ------------ At this point the strings are fragmenting to hadrons -------------
  //
  G4QHadronVector* theResult = new G4QHadronVector;
  G4LorentzVector KTsum(0.,0.,0.);
  G4bool NeedEnergyCorrector=false;            // @@ M.K.?
#ifdef debug
  G4cout<<"G4QFragmentation::Scatter: BefFragmentation, nString="<<strings->size()<<G4endl;
#endif
  for( unsigned astring=0; astring < strings->size(); astring++)
  {
    G4QString* curString=(*strings)[astring];
    G4LorentzVector curString4M = curString->Get4Momentum();
    KTsum+= curString4M;
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: String#"<<astring<<",sum4M="<<KTsum<<G4endl;
#endif
    if(!(KTsum.e()<1.) && !(KTsum.e()>-1.))    // NAN check
    {
      G4cerr<<"***G4QFragmentation::Scatter: KTsum="<<KTsum<<G4endl;
      G4Exception("G4QFragmentation::Scatter:","72",FatalException,"NANin3Vector");
    }
    G4QHadronVector* theHadrons = 0;           // A prototype of the string frag. output
    if( curString->IsExcited() ) theHadrons=curString->FragmentString(true); // FragmString
    else
    {
      //theHadrons = new G4QHadronVector;
      //theHadrons->push_back(curString->GetAsQHadron()); //@@ NotImplem
    }    
    if (!theHadrons)
    {
      G4cerr<<"-Warning-G4QFragmentation::Scatter: No Hadrons produced"<<G4endl;
      continue;
    }
    G4LorentzVector KTsum1(0.,0.,0.,0.);
    for(unsigned aTrack=0; aTrack<theHadrons->size(); aTrack++)
    {
      G4QHadron* curHadron=(*theHadrons)[aTrack];
      theResult->push_back(curHadron);
      KTsum1+= curHadron->Get4Momentum();
    }
 
    if(std::abs(1.-curString4M.e()/KTsum1.e()) > perMillion) NeedEnergyCorrector=true;

    //      clean up (the issues are filled to theResult)
#ifdef debug
    G4cout<<"G4QFragmentation::Scatter: Done, EnCor="<<NeedEnergyCorrector<<G4endl;
#endif
    delete theHadrons;
  }
#ifdef debug
  G4cout<<"G4QFragmentation::Scatter: String4mom="<<KTsum<<G4endl; 
#endif
  if(NeedEnergyCorrector) EnergyAndMomentumCorrector(theResult, KTsum);
#ifdef debug
  G4cout<<"G4QFragmentation::Scatter:***Done***, after Fragm NS="<<strings->size()<<G4endl;
#endif
  //
  // >>>>>>>>>>>>>> clean-up used Nucleons, Strings. and theNucleus
  //
  std::for_each(theTargets.begin(), theTargets.end(), DeleteQHadron());
  theTargets.clear();
  delete theNucleus;
  std::for_each(strings->begin(), strings->end(), DeleteQString() );
  delete strings;

  return theResult; // This is the resulting hadron vector (the answer)
}

// @@ Shoul be updated or deleted, if the energy-momentum nonconsertvation does not exist
G4bool G4QFragmentation::EnergyAndMomentumCorrector(G4QHadronVector* Output,
                                                    G4LorentzVector& TotalCollisionMom)   
{
  const int    nAttemptScale = 500;
  const double ErrLimit = 1.E-5;
  if (Output->empty()) return TRUE;
  G4LorentzVector SumMom;
  G4double        SumMass = 0;     
  G4double        TotalCollisionMass = TotalCollisionMom.m();
  if( !(TotalCollisionMass<1) && !(TotalCollisionMass>-1) )
  {
    G4cerr<<"***G4QFragmentation::EnergyAndMomentumCorrect:M="<<TotalCollisionMass<<G4endl;
    G4Exception("G4QFragmentation::EnergyAndMomentumCorr:","72",FatalException,"NAN_totM");
  }
  // Calculate sum hadron 4-momenta and summing hadron mass
  unsigned int cHadron;
  for(cHadron = 0; cHadron < Output->size(); cHadron++)
  {
    SumMom  += Output->operator[](cHadron)->Get4Momentum();
    if( !(SumMom<1) && !(SumMom>-1) )
    {
      G4cerr<<"***G4QFragmentation::EnergyAndMomentumCorrector: SumMom="<<SumMom<<G4endl;
      G4Exception("G4QFragmentation::EnergyAndMomentumCorr:","72",FatalException,"NANMom");
    }
    SumMass += (*Output)[cHadron]->GetMass();
    if(!(SumMass<1) && !(SumMass>-1))
    {
      G4cerr<<"***G4QFragmentation::EnergyAndMomentumCorrector: SumMass="<<SumMass<<G4endl;
      G4Exception("G4QFragmentation::EnergyAndMomentumCor:","72",FatalException,"NANMass");
    }
  }
  // Cannot correct a single particle
  if(Output->size() < 2) return FALSE;
  if (SumMass > TotalCollisionMass) return FALSE;
  SumMass = SumMom.m2();
  if (SumMass < 0) return FALSE;
  SumMass = std::sqrt(SumMass);
  // Compute c.m.s. hadron velocity and boost KTV to hadron c.m.s.
  G4ThreeVector Beta = -SumMom.boostVector();
  G4int nOut=Output->size();
  if(nOut) for(G4int o=0; o<nOut; o++) (*Output)[o]->Boost(Beta);
  // Scale total c.m.s. hadron energy (hadron system mass).
  // It should be equal interaction mass
  G4double Scale = 1;
  G4int cAttempt = 0;
  G4double Sum = 0;
  G4bool success = false;
  for(cAttempt = 0; cAttempt < nAttemptScale; cAttempt++)
  {
    Sum = 0;
    for(cHadron = 0; cHadron < Output->size(); cHadron++)
    {
      G4LorentzVector HadronMom = Output->operator[](cHadron)->Get4Momentum();
      HadronMom.setVect(Scale*HadronMom.vect());
      G4double E = std::sqrt(HadronMom.vect().mag2() + sqr((*Output)[cHadron]->GetMass()));
      HadronMom.setE(E);
      Output->operator[](cHadron)->Set4Momentum(HadronMom);
      Sum += E;
    }   
    Scale = TotalCollisionMass/Sum;    
    if (Scale - 1 <= ErrLimit) 
    {
      success = true;
      break;
    }
#ifdef debug
    G4cout<<"G4QFragmentation::EnergyAndMomentumCorrector: Scale-1="<<Scale-1<<", TotM="
          <<TotalCollisionMass<<", Sum="<<Sum<<G4endl;
#endif
  }
#ifdef debug
  if(!success)
  {
    G4cout<<"***G4QFragmentation::EnergyAndMomentumCorrector: Scale #1 at end of loop M="
          <<TotalCollisionMass<<", S"<<Sum<<", Sc="<<Scale
          <<" Increase number of attempts or increase the ErrLimit parameter"<<G4endl;
    //G4Exception("G4QFragmentation::EnergyAndMomCor:","72",FatalException,"NoCorrection");
  }
#endif
  // Compute c.m.s. interaction velocity and KTV back boost
  Beta = TotalCollisionMom.boostVector(); 
  nOut=Output->size();
  if(nOut) for(G4int o=0; o<nOut; o++) (*Output)[o]->Boost(Beta);
  return TRUE;
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
  G4cout<< "G4QFragm::ExciteSingDiffParticipants: TargetMass="<<Ptarget.mag()<<G4endl;
#endif
  target->Set4Momentum(Ptarget);  
#ifdef debug
  G4cout<<"G4QFragm::ExciteDiffParticipants:ProjectileMass="<<Pprojectile.mag()<<G4endl;
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
  G4double pt2; do{pt2=widthSq*std::log(G4UniformRand());} while (pt2>maxPtSquare);
  pt2=std::sqrt(pt2);
  G4double phi=G4UniformRand()*twopi;
  return G4ThreeVector(pt2*std::cos(phi),pt2*std::sin(phi),0.);    
} // End of GaussianPt

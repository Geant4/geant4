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
// $Id: G4QIonIonCollision.cc,v 1.9 2010-06-19 07:46:44 mkossov Exp $
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

#include "G4QIonIonCollision.hh"

// Promoting model parameters from local variables class properties @@(? M.K.)

// Definition of static parameters
G4int    G4QIonIonCollision::nCutMax=7; 
G4double G4QIonIonCollision::stringTension=1.*GeV/fermi;
G4double G4QIonIonCollision::tubeDensity  =1./fermi;
// Parameters of diffractional fragmentation
G4double G4QIonIonCollision::widthOfPtSquare=-0.75*GeV*GeV; // ptWidth2 forStringExcitation

G4QIonIonCollision::G4QIonIonCollision(G4QNucleus &pNucleus, const G4QNucleus &tNucleus)
{
  static const G4double mProt= G4QPDGCode(2212).GetMass(); // Mass of proton
  static const G4double mPi0= G4QPDGCode(111).GetMass();   // Mass of Pi0
  theWorld= G4QCHIPSWorld::Get();         // Get a pointer to the CHIPS World
  theResult = new G4QHadronVector;        // Define theResultingHadronVector
  G4bool stringsInitted=false;            // Strings are initiated
  G4LorentzRotation toZ;                  // Lorentz Transformation to the projectileSystem
  G4LorentzVector proj4M=pNucleus.Get4Momentum(); // Projectile nucleus 4-momentum in LS
  //G4double projM=pNucleus.GetGSMass();    // Ground state mass of the projectile nucleus
  G4double targM=tNucleus.GetGSMass();    // Ground state mass of the target nucleus
#ifdef edebug
  G4cout<<"G4QIn::Constr:*Called*,pA="<<pNucleus<<proj4M<<",tA="<<tNucleus<<targM<<G4endl;
#endif
  G4int pZ=pNucleus.GetZ();
  G4int pN=pNucleus.GetN();
  G4int pA=pZ+pN;
  G4int pPDG=90000000+pZ*1000+pN;
  G4int tZ=tNucleus.GetZ();
  G4int tN=tNucleus.GetN();
  G4int tA=tZ+tN;
  G4int tPDG=90000000+tZ*1000+tN;
  toZ.rotateZ(-proj4M.phi());
  toZ.rotateY(-proj4M.theta());
  G4LorentzVector zProj4M=toZ*proj4M;     // Proj 4-momentum in LS rotated to Z axis
  pNucleus.Set4Momentum(zProj4M);         // Now the projectile nucleus moves along Z axis
#ifdef edebug
  G4int totChg=pZ+tZ;                     // Charge of the projectile+target for the CHECK
  G4int totBaN=pA+tA;                     // Baryon Number of Proj+Targ for CHECK
  G4LorentzVector tgLS4M(0.,0.,0.,targM); // Target 4-momentum in LS
  G4LorentzVector totLS4M=proj4M+tgLS4M;  // Total 4-momentum in LS
  G4LorentzVector totZLS4M=zProj4M+tgLS4M;// Total 4-momentum in LS with momentum along Z
  G4cout<<"-EMC-G4QIonIonCollision::Constr: tLS4M="<<totLS4M<<",tZLS4M="<<totZLS4M<<G4endl;
  // === From nere all consideration is made in the rotated LS frame (proj is along Z) ===
#endif
  G4LorentzRotation toLab(toZ.inverse()); // Lorentz Transfornation "ZLS"->LS (at the end)
  G4QInteractionVector theInteractions;   // A vector of interactions (taken from the Body)
  G4QPartonPairVector  thePartonPairs;    // The parton pairs (taken from the Body)
  G4int                ModelMode=SOFT;    // The current model type (taken from the Body)
  theTargNucleus=G4QNucleus(tZ,tN);       // Create theTargNucleus to Move From LS to CM
  theTargNucleus.InitByPDG(tPDG);         // Reinit the Nucleus for the new Attempt?
  theTargNucleus.Init3D();                // 3D-initialisation(nucleons)ofTheTGNucleusClone
#ifdef debug
  G4cout<<"G4QIonIonCollision::Constr:TargNuclMom="<<theTargNucleus.Get4Momentum()<<G4endl;
#endif
  theProjNucleus=G4QNucleus(pZ,pN);       // Create theProjNucleus to Move From LS to CM
  theProjNucleus.InitByPDG(pPDG);         // Reinit the Nucleus for the new Attempt?
  theProjNucleus.Init3D();                // 3D-initialisation(nucleons)ofThePrNucleusClone
#ifdef debug
  G4cout<<"G4QIonIonCollision::Constr:ProjNuclMom="<<theProjNucleus.Get4Momentum()<<G4endl;
#endif
#ifdef edebug
  G4LorentzVector sumP1=theProjNucleus.GetNucleons4Momentum();// Sum ofPrNucleons4M inRotLS
  G4LorentzVector sumT1=theTargNucleus.GetNucleons4Momentum();// Sum ofTgNucleons4M inRotLS
  G4cout<<"-EMC-G4QIonIonCollision::Construct: pNuc4M="<<sumP1<<", tNuc4M="<<sumT1<<G4endl;
#endif
  G4ThreeVector theCMVelocity(0.,0.,0.);             // Prototype of the "CM" velocity
  G4ThreeVector theALVelocity(0.,0.,0.);             // Prototype of "CMAntiLab" velocity
  // Very important! Here the projectile 4M is recalculated in the ZLS (previously in LS)
  proj4M = pNucleus.Get4Momentum();                  // 4-mom of theProjectileHadron inZLS
  G4double pz_per_projectile = proj4M.pz()/pA;       // Momentum per nucleon
  G4double e_per_projectile = proj4M.e()/pA+targM/tA;// Add MassOfTargetNucleon
  G4double vz = pz_per_projectile/e_per_projectile;  // Speed of the "oneNuclenCM"
#ifdef debug
  G4cout<<"G4QIonIonCollision::Construct: (КщеЯ)Projectile4M="<<proj4M<<", Vz="<<vz
        <<", pA="<<pA<<", pE="<<e_per_projectile<<G4endl;
#endif
  theCMVelocity.setZ(vz);                            // CM (w/r to one nucleon) velosity
  theProjNucleus.DoLorentzBoost(-theCMVelocity);     // BoostProjNucleus to "CM"
  theTargNucleus.DoLorentzBoost(-theCMVelocity);     // BoostProjNucleus to "CM"
  G4LorentzVector prN4M=theProjNucleus.Get4Momentum();
  G4double avz=prN4M.pz()/prN4M.e();                 // Speed of AntiLabSys in CMS
  theALVelocity.setZ(avz);                           // CM (w/r to one nucleon) velosity
#ifdef edebug
  G4LorentzVector sumP2=theProjNucleus.GetNucleons4Momentum();// Sum of Nucleons4M in RotCM
  G4LorentzVector sumT2=theTargNucleus.GetNucleons4Momentum();// Sum of Nucleons4M in RotCM
  G4cout<<"-EMC-G4QIonIonCollision::Construct: Boosted, vCM="<<vz<<", vAL="<<avz<<", tN4M="
        <<sumT2<<", pN4M="<<sumP2<<G4endl;
#endif
  G4int maxCuts = 7;                                 // @@ Can be reduced for low energies
  //
  // >>>>>>>>>> Find collisions meeting collision conditions and the NN interaction XS
  //
  G4double outerRadius = theProjNucleus.GetOuterRadius()+theTargNucleus.GetOuterRadius();
  G4QProbability theProbability(2212);               // *** proj/targ nucleons ***
  // Clean up all previous interactions and reset the counters
#ifdef debug
  G4cout<<"G4QIonIonCollision::Construct: theIntSize="<<theInteractions.size()<<G4endl;
#endif
  G4int attCnt=-1;
  G4int maxAtt=27;
  // From here a loop over nucleons of the projectile Nucleus (@@ Z-sorting isn't implem!!)
  //
  G4QHadron* pNucleon; // Proto of the Projectile Nucleon pointer
  G4QHadron* tNucleon; // Proto of the Target Nucleon pointer
  G4QNucleus curProjNucleus;
  G4QNucleus curTargNucleus;
  while(!theInteractions.size() && ++attCnt < maxAtt)// TillAtLeastOneInteractionIsCreated
  {
    // std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
    // Here we need to clean up the ProjNucleon and the TargNucleon in the Interactions !
    G4int nint=theInteractions.size();               // For the 1st attempt should be zero
    for(G4int ni=0; ni < nint; ++ni)
    {
      G4QInteraction* curInt = theInteractions[ni];
      delete curInt->GetProjectile();
      delete curInt->GetTarget();
      delete curInt;
    }
    theInteractions.clear();
    // Now we need to make a copy of the projectile and the target nuclei with 3D info
    if(attCnt)                                       // This is not theFirstAttempt->Clean
    {
      curProjNucleus.DeleteNucleons();
      curTargNucleus.DeleteNucleons();
    }
    curProjNucleus = theProjNucleus;
    curTargNucleus = theTargNucleus;
    // choose random impact parameter
    std::pair<G4double, G4double> theImpactParameter;
    theImpactParameter = curProjNucleus.ChooseImpactXandY(outerRadius);
    G4double impactX = theImpactParameter.first;
    G4double impactY = theImpactParameter.second;
#ifdef debug
    G4cout<<"G4QIonIonCollision::Con:At#="<<attCnt<<",X="<<impactX<<",Y="<<impactY<<G4endl;
#endif
    curProjNucleus.StartLoop();                      // Prepare Loop ovder nucleons
    G4int pnCnt=0;
    G4int pnA=curProjNucleus.GetA();
#ifdef debu
    G4cout<<"G4QIonIonCollision::Constr: Before the WHILE over pNucleons,pA="<<pnA<<G4endl;
#endif
    while((pNucleon=curProjNucleus.GetNextNucleon()) && pnCnt < pnA) // @@ can be for LOOP?
    {
      ++pnCnt;
      G4QInteractionVector curInteractions;          // A temporary vector of interactions
      G4LorentzVector pNuc4M=pNucleon->Get4Momentum();
      G4ThreeVector pNucR=pNucleon->GetPosition();   // Position of the pNucleon WRTo pNucl
      G4double pImpX = impactX+pNucR.x();
      G4double pImpY = impactY+pNucR.y();
      G4double prEn=proj4M.e();                      // For mesons
      G4int proB=pNucleus.GetBaryonNumber();
      if     (proB>0) prEn-=pNucleus.GetMass();      // For baryons
      else if(proB<0) prEn+=mProt;                   // For anti-baryons
#ifdef debug
      G4double impactUsed = 0.;
      G4cout<<"G4QIonIonCollision::Construct: estimated energy, prEn="<<prEn<<G4endl;
#endif
      G4int totalCuts = 0;
      // @@ This is a fake (random) s calculation @@ can be done inside the TARG-while
      G4int tnA=curTargNucleus.GetA();
      G4double pImp=std::sqrt(pImpX*pImpX+pImpY*pImpY);   
      G4double eflen=curTargNucleus.GetThickness(pImp); // EffectiveLength
      maxEn=eflen*stringTension;                     // maxAbsorbedEnergy in IndUnits=MeV
      maxNuc=static_cast<int>(eflen*tubeDensity+0.5);// max #0f involved nuclear nucleons
#ifdef debug
      G4cout<<"G4QIonIonCollision::Con:pE="<<prEn<<" < mE="<<maxEn<<",mN="<<maxNuc<<G4endl;
#endif
      if(prEn < maxEn)                               // Create DIN interaction and go out
      {
        G4QHadron* aProjectile = new G4QHadron(*pNucleon);// Copy selected PrNuc for String
        curTargNucleus.StartLoop();                  // Initialize newSelection forNucleons
        tNucleon=curTargNucleus.GetNextNucleon();    // Select a nucleon
        G4QHadron* aTarget = new G4QHadron(*tNucleon);// Copy selected nucleon for String
        G4QInteraction* anInteraction = new G4QInteraction(aProjectile);
        anInteraction->SetTarget(aTarget); 
        anInteraction->SetNumberOfDINRCollisions(1); // Consider this interaction as DINR
        curInteractions.push_back(anInteraction);    //--> now the Interaction is not empty
        curProjNucleus.DoLorentzBoost(-theALVelocity);// Boost theResPrNucleus toRotAntiLab
        curProjNucleus.SubtractNucleon(pNucleon);    // Pointer to the used ProjNucleon
        curProjNucleus.DoLorentzBoost(theALVelocity);// Boost theResProjNucleus back to CM
        curTargNucleus.DoLorentzBoost(theCMVelocity);// Boost theResTgNucleus to RotatedLS
        curTargNucleus.SubtractNucleon(tNucleon);    // Pointer to the used TargNucleon
        curTargNucleus.DoLorentzBoost(-theCMVelocity);// Boost theResTargNucleus back to CM
#ifdef debug
        G4cout<<"G4QIonIonCollision::Construct: DINR interaction is created"<<G4endl;
#endif
        break;                                       // Break the WHILE of the pNucleons
      }
      // LOOP over nuclei of the target nucleus to select collisions
      curTargNucleus.StartLoop();                    // To get the same nucleon
      G4int tnCnt = 0;
#ifdef debu
      G4cout<<"G4QIonIonCollision::Constr: Before WHILE over tNucleons, tA="<<tnA<<G4endl;
#endif
      while((tNucleon=curTargNucleus.GetNextNucleon()) && tnCnt<tnA && totalCuts<maxCuts)
      {
        ++tnCnt;
        G4LorentzVector tNuc4M=tNucleon->Get4Momentum();
#ifdef debug
        G4cout<<"G4QIonIonCollision::Constr: OuterR="<<outerRadius<<", mC="<<maxCuts
              <<", pA="<<curProjNucleus<<", tA="<<curTargNucleus<<G4endl;
#endif
        // Check the reaction threshold 
        G4double s = (tNuc4M + pNuc4M).mag2();         // Squared CM Energy of compound
        G4double ThresholdMass = pNucleon->GetMass() + tNucleon->GetMass();
#ifdef debug
        G4cout<<"G4QInel::Constr: s="<<s<<", ThreshM="<<sqr(ThresholdMass)<<G4endl;
#endif
        ModelMode = SOFT;                              // NOT-Diffractive hadronization
        if (s < 0.)                                    // At ThP=0 is impossible(virtNucl)
        {
          G4cerr<<"*G4QInelast::Constr:s="<<s<<",pN4M="<<pNuc4M<<",tN4M="<<tNuc4M<<G4endl;
          G4Exception("G4QIonIonCollision::Construct:","72",FatalException,"LowEn(NegS)");
        }
        if (s < sqr(ThresholdMass))                    // --> Only diffractive interaction
        {
#ifdef debug
          G4cout<<"G4QIonIonCollision::Constr: ***OnlyDiffraction***, ThM="<<ThresholdMass
                <<">sqrt(s)="<<std::sqrt(s)<<" -> only Diffraction is possible"<<G4endl;
#endif
          ModelMode = DIFFRACTIVE;
        }
        G4ThreeVector tNucR=tNucleon->GetPosition(); // Position of the tNucleon WRTo tNucl
        G4double dImpX = pImpX-tNucR.x();
        G4double dImpY = pImpY-tNucR.y();
        G4double Distance2=dImpX*dImpX+dImpY*dImpY;
#ifdef sdebug
        G4cout<<"G4QIonIonCollision::Construct: s="<<s<<", D2="<<Distance2<<G4endl;
#endif
        // Needs to be moved to Probability class @@
        if(s<=10000.)
        {
          G4cout<<"-Warning-G4QIonIonCollision::Construct: s < .01 GeV^2, p4M="
                <<pNucleon->Get4Momentum()<<",t4M="<<tNucleon->Get4Momentum()<<G4endl;
          continue;                                  // skip the rest of the targetNucleons
        }
        G4double Probability = theProbability.GetPomInelProbability(s, Distance2);// P_INEL
        // test for inelastic collision
#ifdef sdebug
        G4cout<<"G4QIonIonCollision::Construct: Probubility="<<Probability<<G4endl;
#endif
        G4double rndNumber = G4UniformRand();        // For the printing purpose
        // ModelMode = DIFFRACTIVE;
#ifdef sdebug
        G4cout<<"G4QIonIonCollision::Construct: NLOOP prob="<<Probability<<", rndm="
              <<rndNumber<<", d="<<std::sqrt(Distance2)<<G4endl;
#endif
        if (Probability > rndNumber) // Inelastic (diffractive or soft) interaction (JOB)
        {
          G4QHadron* aProjectile = new G4QHadron(*pNucleon);// Copy selected pNuc to String
          G4QHadron* aTarget = new G4QHadron(*tNucleon);// Copy selected tNucleon to String
#ifdef edebug
          G4cout<<"--->EMC-->G4QIonIonCollision::Construct: TargNucleon is filled, 4M/PDG="
                <<aTarget->Get4Momentum()<<aTarget->GetPDGCode()<<G4endl;
#endif
          // Now the energy of the nucleons must be updated in CMS
          curProjNucleus.DoLorentzBoost(-theALVelocity);// Boost theResNucleus toRotatedLS
          curProjNucleus.SubtractNucleon(pNucleon);     // Pointer to the used nucleon
          curProjNucleus.DoLorentzBoost(theALVelocity); // Boost theResNucleus back to CM
          curTargNucleus.DoLorentzBoost(theCMVelocity); // Boost theResNucleus toRotatedLS
          curTargNucleus.SubtractNucleon(tNucleon);     // Pointer to the used nucleon
          curTargNucleus.DoLorentzBoost(-theCMVelocity);// Boost theResNucleus back to CM
          if((theProbability.GetPomDiffProbability(s,Distance2)/Probability >
              G4UniformRand() && ModelMode==SOFT ) || ModelMode==DIFFRACTIVE)
          { 
            // ------------->>>> diffractive interaction @@ IsSingleDiffractive called once
            if(IsSingleDiffractive()) ExciteSingDiffParticipants(aProjectile, aTarget);
            else                          ExciteDiffParticipants(aProjectile, aTarget);
            G4QInteraction* anInteraction = new G4QInteraction(aProjectile);
            anInteraction->SetTarget(aTarget); 
            anInteraction->SetNumberOfDiffractiveCollisions(1); // Why not increment? M.K.?
            curInteractions.push_back(anInteraction);//--> now theInteractions not empty
            // @@ Why not breake the NLOOP, if only one diffractive can happend?
            totalCuts++;                             // UpdateOfNucleons in't necessary
#ifdef debug
            G4cout<<"G4QIonIonCollision::Constr:NLOOP DiffInteract,tC="<<totalCuts<<G4endl;
#endif
          }
          else
          {
            // -------------->>>>> nondiffractive = soft interaction
            // sample nCut+1 (cut Pomerons) pairs of strings can be produced
            G4int nCut;                              // Result in a chosen number of cuts
            G4double* running = new G4double[nCutMax];// @@ This limits the max cuts
            for(nCut = 0; nCut < nCutMax; nCut++)    // Calculates multiCut probabilities
            {
              running[nCut]= theProbability.GetCutPomProbability(s, Distance2, nCut+1);
              if(nCut) running[nCut] += running[nCut-1];// Sum up with the previous one
            }
            G4double random = running[nCutMax-1]*G4UniformRand();
            for(nCut = 0; nCut < nCutMax; nCut++) if(running[nCut] > random) break; // tNuc
            delete [] running;
#ifdef debug
            G4cout<<"G4QIonIonCollision::Construct: NLOOP-Soft Chosen nCut="<<nCut<<G4endl;
#endif
            // @@ If nCut>0 interaction with a few nucleons is possible
            // @@ nCut is found with big efforts and now nCut=0 ?? M.K. ?? !!
            //nCut = 0; // @@ in original code ?? @@
            aTarget->IncrementCollisionCount(nCut+1); // @@ What about multyNucleon target?
            aProjectile->IncrementCollisionCount(nCut+1);
            G4QInteraction* anInteraction = new G4QInteraction(aProjectile);
            anInteraction->SetTarget(aTarget);
            anInteraction->SetNumberOfSoftCollisions(nCut+1);
            curInteractions.push_back(anInteraction); // Now curInteractions are not empty
            totalCuts += nCut+1;
#ifdef debug
            G4cout<<"G4QIonIonCollision::Constr:NLOOP SoftInteract,tC="<<totalCuts<<G4endl;
            impactUsed=Distance2;
#endif
          }
        }// Probability selection
      } // End of While over target nucleons
      G4int nCurInt=curInteractions.size();
      for(G4int ni=0; ni < nCurInt; ++ni) theInteractions.push_back(curInteractions[ni]);
      curInteractions.clear();                      // Probably, not necessary...
    } // End of WHILE over projectile nucleons
  } // End of WHILE over attempts to find at least one nucleus-nucleus interaction
  theProjNucleus.DeleteNucleons();
  theTargNucleus.DeleteNucleons();
  theProjNucleus = curProjNucleus;
  theTargNucleus = curTargNucleus;
  curProjNucleus.DeleteNucleons();
  curTargNucleus.DeleteNucleons();
  G4int nInt=theInteractions.size();
#ifdef debug
  G4cout<<"G4QIonIonCollision::Constr: #ofInteractions = "<<nInt<<", #OfDINR = "
        <<theInteractions[0]->GetNumberOfDINRCollisions()<<G4endl;
#endif
  if(!nInt || (nInt==1 && theInteractions[0]->GetNumberOfDINRCollisions()==1)) // QFreeInel
  {
    G4QHadron* aTarget;
    G4QHadron* aProjectile;
    if(nInt)                                         // Take Targ/Proj from the Interaction
    {
     	aTarget=theInteractions[0]->GetTarget();
	     aProjectile=theInteractions[0]->GetProjectile();
      delete theInteractions[0];
      theInteractions.clear();
    }
    else                                             // Create a new target nucleon (?)
    {
      theProjNucleus.StartLoop();                    // To get the same nucleon from projec
      pNucleon=theProjNucleus.GetNextNucleon();      // Get theNucleon to create projectNuc
      aProjectile = new G4QHadron(*pNucleon);        // Copy selected pNucleon for String
      theProjNucleus.DoLorentzBoost(-theALVelocity); // Boost theResProjNucleus toRotatedLS
      theProjNucleus.SubtractNucleon(pNucleon);      // Pointer to SelProjNucleon to delete
      theProjNucleus.DoLorentzBoost(theALVelocity);  // Boost theResProNucleus back to CMS
      theTargNucleus.StartLoop();                    // To get the same nucleon from target
      tNucleon=theTargNucleus.GetNextNucleon();      // Get theNucleon to create targetNucl
      aTarget = new G4QHadron(*tNucleon);            // Copy selected tNucleon for String
      theTargNucleus.DoLorentzBoost(theCMVelocity);  // Boost theResTargNucleus toRotatedLS
      theTargNucleus.SubtractNucleon(tNucleon);      // Pointer to SelTargNucleon to delete
      theTargNucleus.DoLorentzBoost(-theCMVelocity); // Boost theResTargNucleus back to CMS
    }
    G4QContent QQC=aTarget->GetQC()+aProjectile->GetQC(); // QContent of the compound
    G4LorentzVector Q4M=aTarget->Get4Momentum()+aProjectile->Get4Momentum(); // 4-mom of Q
    delete aTarget;
    delete aProjectile;
    // @@ 4-Mom should be converted to LS for the TargQuasmon and to AL for the ProjQuasmon
    Q4M.boost(theCMVelocity);
    Q4M=toLab*Q4M;
    G4Quasmon* stringQuasmon = new G4Quasmon(QQC, Q4M);
    theQuasmons.push_back(stringQuasmon);
    theTargNucleus.DoLorentzBoost(theCMVelocity);   // BoostTheResidualNucleus toRotatedLS
    theTargNucleus.DoLorentzRotation(toLab);        // Recove Z-direction in LS ("LS"->LS)
    return;
  }
  //
  // ------------------ now build the parton pairs for the strings ------------------
  //
  for(G4int i=0; i<nInt; i++)
  {
    theInteractions[i]->SplitHadrons();
#ifdef edebug
    G4QHadron* projH=theInteractions[i]->GetProjectile(); // Projectile of theInteraction
    G4QHadron* targH=theInteractions[i]->GetTarget();     // Target of the Interaction
    G4LorentzVector pSP(0.,0.,0.,0.);               // Sum of parton's 4mom's for proj
    G4LorentzVector tSP(0.,0.,0.,0.);               // Sum of parton's 4mom's for proj
    std::list<G4QParton*> projCP=projH->GetColor(); // Pointers to proj Color-partons
    std::list<G4QParton*> projAC=projH->GetAntiColor();// PointersTo projAntiColorPartons
    std::list<G4QParton*> targCP=targH->GetColor(); // Pointers to targ Color-partons
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
    G4cout<<"-EMC-G4QIonIonCollision::Construct: Interaction#"<<i<<",dP4M="
          <<projH->Get4Momentum()-pSP<<",dT4M="<<targH->Get4Momentum()-tSP<<G4endl;
#endif
  }  
  // 
  // >>>>>>>> make soft collisions (ordering is vital)
  //
  G4QInteractionVector::iterator it;
#ifdef debug
  G4cout<<"G4QIonIonCollision::Constr: Creation ofSoftCollisionPartonPair STARTS"<<G4endl;
#endif
  G4bool rep=true;
  while(rep && theInteractions.size())
  {
   for(it = theInteractions.begin(); it != theInteractions.end(); ++it)   
   {
    G4QInteraction* anIniteraction = *it;
    G4QPartonPair*  aPair=0;
    G4int nSoftCollisions = anIniteraction->GetNumberOfSoftCollisions();
#ifdef debug
    G4cout<<"G4QIonIonCollision::Construct: #0f SoftCollisions ="<<nSoftCollisions<<G4endl;
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
        G4cout<<"--->G4QIonIonCollision::Constr: SOFT, 2 parton pairs are filled"<<G4endl;
#endif
      }  
      delete *it;
      it=theInteractions.erase(it);      // SoftInteractions're converted&erased
      if( it != theInteractions.begin() )// To avoid going below begin() (for Windows)
      {
        it--;
        rep=false;
      }
      else
      {
        rep=true;
        break;
      }
    }
    else rep=false;
   }
  }
#ifdef debug
  G4cout<<"G4QIonIonCollision::Constr: -> Parton pairs for SOFT strings are made"<<G4endl;
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
    G4cout<<"G4QIonIonCollision::Constr: CreationOfDiffractivePartonPairs, i="<<i<<G4endl;
#endif
    // the projectile diffraction parton pair is created first
    G4QHadron* aProjectile = anIniteraction->GetProjectile();
    G4QParton* aParton = aProjectile->GetNextParton();
    if (aParton)
    {
      aPartonPair = new G4QPartonPair(aParton, aProjectile->GetNextAntiParton(), 
                                      G4QPartonPair::DIFFRACTIVE,
                                      G4QPartonPair::PROJECTILE);
      thePartonPairs.push_back(aPartonPair);
#ifdef debug
      G4cout<<"G4QIonIonCollision::Constr: proj Diffractive PartonPair is filled"<<G4endl;
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
      G4cout<<"G4QIonIonCollision::Constr: target Diffractive PartonPair's filled"<<G4endl;
#endif
    }
  }
#ifdef debug
  G4cout<<"G4QIonIonCollision::Construct: DiffractivePartonPairs are created"<<G4endl;
#endif  
  //
  // >>>>>>>>>>>>>> clean-up  Interactions, if necessary
  //
  std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
  theInteractions.clear();
#ifdef debug
  G4cout<<"G4QIonIonCollision::Construct: Temporary objects are cleaned up"<<G4endl;
#endif  
  // This function prepares theBoost for transformation of secondaries to LS (-ProjRot!)
  theProjNucleus.DoLorentzBoost(theCMVelocity);// Boost theResidualProjNucleus to RotatedLS
  theTargNucleus.DoLorentzBoost(theCMVelocity);// Boost theResidualTargNucleus to RotatedLS
  // @@ Nucleus isn't completely in LS, it needs the toZ (-ProjRot) rotation to consE/M
#ifdef debug
  G4cout<<">>>>>>>>>>>>G4QIonIonCollision::Construct: >>>>>> Strings are created "<<G4endl;
#endif
  G4QPartonPair* aPair;
  G4QString* aString=0;
  while(thePartonPairs.size()) // @@ At present noDifference in stringBuild (? M.K.)
  {
    aPair = thePartonPairs.back();           // Get the parton pair
    thePartonPairs.pop_back();               // Clean up thePartonPairPointer in the Vector
#ifdef debug
    G4cout<<"G4QIonIonCollision::Construct:StringType="<<aPair->GetCollisionType()<<G4endl;
#endif
    aString= new G4QString(aPair);
#ifdef debug
    G4cout<<"G4QIonIonCollision::Construct:NewString4M="<<aString->Get4Momentum()<<G4endl;
#endif
    aString->Boost(theCMVelocity);       // ! Strings are moved to ZLS when pushed !
    strings.push_back(aString);
    stringsInitted=true;
    delete aPair;
  } // End of the String Creation LOOP
#ifdef edebug
  G4LorentzVector sum=theProjNucleus.Get4Momentum()+theTargNucleus.Get4Momentum(); // in LS
  G4int rChg=totChg-theProjNucleus.GetZ()-theTargNucleus.GetZ();
  G4int rBaN=totBaN-theProjNucleus.GetA()-theTargNucleus.GetA();
  G4int nStrings=strings.size();
  G4cout<<"-EMC-G4QIonIonCollision::Constr:#ofString="<<nStrings<<",tNuc4M="<<sum<<G4endl;
  for(G4int i=0; i<nStrings; i++)
  {
    G4QString* prString=strings[i];
    G4LorentzVector strI4M=prString->Get4Momentum();
    sum+=strI4M;
    G4int      sChg=prString->GetCharge();
    G4int      sBaN=prString->GetBaryonNumber();
    G4int      LPDG=prString->GetLeftParton()->GetPDGCode();
    G4int      RPDG=prString->GetRightParton()->GetPDGCode();
    G4QContent LQC =prString->GetLeftParton()->GetQC();
    G4QContent RQC =prString->GetRightParton()->GetQC();
    rChg-=sChg;
    rBaN-=sBaN;
    G4cout<<"-EMC-G4QIonIonCollision::Construct: String#"<<i<<", 4M="<<strI4M<<", LPDG="
          <<LPDG<<LQC<<",RPDG="<<RPDG<<RQC<<", Ch="<<sChg<<", BN="<<sBaN<<G4endl;
  }
  G4cout<<"-EMC-G4QInel::Constr: r4M="<<sum-totZLS4M<<", rC="<<rChg<<", rB="<<rBaN<<G4endl;
#endif
  if(!stringsInitted)
  {
    G4cerr<<"*****G4QIonIonCollision::Construct:**** No strings are created *****"<<G4endl;
    G4Exception("G4QIonIonCollision::Construct:","72",FatalException,"No Strings created");
  }
#ifdef debug
  G4cout<<"G4QIonIonCollision::Constr:BeforeRotation, #0fStrings="<<strings.size()<<G4endl;
#endif
  //
  // ---------------- At this point the strings must be already created in "LS" -----------
  //
  for(unsigned astring=0; astring < strings.size(); astring++)
            strings[astring]->LorentzRotate(toLab); // Recove Z-direction in LS ("LS"->LS)
  theProjNucleus.DoLorentzRotation(toLab); // Recove Z-dir in LS ("LS"->LS) for ProjNucleus
  theTargNucleus.DoLorentzRotation(toLab); // Recove Z-dir in LS ("LS"->LS) for TargNucleus
  // Now everything is in LS system
#ifdef edebug
  G4LorentzVector sm=theProjNucleus.Get4Momentum()+theTargNucleus.Get4Momentum();// NucInLS
  G4int rCg=totChg-theProjNucleus.GetZ()-theTargNucleus.GetZ();
  G4int rBC=totBaN-theProjNucleus.GetA()-theTargNucleus.GetA();
  G4int nStrs=strings.size();
  G4cout<<"-EMCLS-G4QIonIonCollision::Constr:#ofS="<<nStrings<<",Nuc4M(E=M)="<<sum<<G4endl;
  for(G4int i=0; i<nStrs; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    sm+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    rCg-=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    rBC-=sBaN;
    G4cout<<"-EMCLS-G4QInel::Construct:String#"<<i<<",4M="<<strI4M<<strI4M.m()<<",Charge="
          <<sChg<<",BaryN="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-...G4QInel::Constr: r4M="<<sm-totLS4M<<", rC="<<rCg<<",rB="<<rBC<<G4endl;
#endif
  //
  // --- Strings are created, but we should try to get rid of negative mass strings -----
  //
  SwapPartons();
  //
  // --- Strings are created, but we should get rid of too light strings (Mmin+MPi0) -----
  //
  G4int problem=0;                                   // 0="no problem", incremented by ASIS
  G4QStringVector::iterator ist;
  G4bool con=true;
  while(con && strings.size())
  {
   for(ist = strings.begin(); ist < strings.end(); ++ist)
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
    G4double cSM=cSM2;
    if(cSM2>0.) cSM=std::sqrt(cSM2);
#ifdef debug
    G4cout<<"G4QIonIonCollision::Construct: Needs Fusion? cLPDG="<<cLPDG<<",cRPDG="<<cRPDG
          <<",cM(cM2 if neg)="<<cSM<<G4endl;
#endif
    if(cSM>0.)                                       // Mass can be calculated
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
	//if(single) miM=G4QPDGCode((*ist)->GetQC().GetSPDGCode()).GetMass();//MinHSMass
#ifdef debug
      G4cout<<"G4QInel::Const:*IsItGood? realM="<<std::sqrt(cSM2)<<" > GSM="<<miM<<G4endl;
#endif
      if(std::sqrt(cSM2) > miM) bad=false;           // String is OK
    }
    if(bad)                                          // String should be merged with others
    {
#ifdef debug
      G4cout<<"G4QInel::Const:TryFuse,L1="<<L1<<",L2="<<L2<<",R1="<<R1<<",R2="<<R2<<G4endl;
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
      for(pst = strings.begin(); pst < strings.end(); pst++) if(pst != ist)
      {
        G4LorentzVector pS4M=(*pst)->Get4Momentum()+cS4M; // Summed 4-momentum
        G4int nLPDG=0;                               // new Left (like in theStringPartner)
        G4int nRPDG=0;                               // new Right(like in theStringPartner)
        G4double pExcess=-DBL_MAX;                   // Prototype of the excess
        G4double pSM2=pS4M.m2();                     // Squared mass of the Fused Strings
#ifdef debug
        G4cout<<"->G4QIonIonCollision::Construct: sum4M="<<pS4M<<", M2="<<pSM2<<G4endl;
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
        G4cout<<"G4QInel::Construct: Partner/w pLPDG="<<pLPDG<<", pRPDG="<<pRPDG<<", pM2="
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
                  ( (cLPDG<0 && (-cLPDG==pL1 || -cLPDG==pL2 || -cLPDG==pRPDG) ) ||
                    (cRPDG<0 && (-cRPDG==pL1 || -cRPDG==pL2 || -cRPDG==pRPDG) )
                  )
                 ) af=1;
          else if(pRPDG > 7 &&
                  ( (cLPDG<0 && (-cLPDG==pR1 || -cLPDG==pR2 || -cLPDG==pLPDG) ) ||
                    (cRPDG<0 && (-cRPDG==pR1 || -cRPDG==pR2 || -cRPDG==pLPDG) )
                  )
                 ) af=2;
          else if(pLPDG <-7 &&
                  ( (cLPDG>0 && ( cLPDG==pL1 || cLPDG==pL2 || cLPDG==-pRPDG) ) ||
                    (cRPDG>0 && ( cRPDG==pL1 || cRPDG==pL2 || cRPDG==-pRPDG) )
                  )
                 ) af=3;
          else if(pRPDG <-7 &&
                  ( (cLPDG>0 && ( cLPDG==pR1 || cLPDG==pR2 || cLPDG==-pLPDG) ) ||
                    (cRPDG>0 && ( cRPDG==pR1 || cRPDG==pR2 || cRPDG==-pLPDG) )
                  )
                 ) af=4;
#ifdef debug
          else G4cout<<"G4QIonIonCollision::Constr:2(QaQ+QDiQ/aQaDiQ) Can't fuse"<<G4endl;
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
          else G4cout<<"G4QIonIonCollision::Constr:3(QDiQ/aQaDiQ+QaQ) Can't fuse"<<G4endl;
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
          else G4cout<<"-G4QIonIonCollision::Constr: 4 (QaQ+aQDiQDiQ) Can't fuse"<<G4endl;
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
          else G4cout<<"-G4QIonIonCollision::Constr: 5 (aQDiQDiQ+QaQ) Can't fuse"<<G4endl;
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
          else G4cout<<"-G4QIonIonCollision::Constr: 6 (QDiQ+aQaDiQ) Can't fuse"<<G4endl;
#endif
        }
#ifdef debug
        G4cout<<"G4QIonIonCollision::Constr: **Possibility**, tf="<<tf<<",af="<<af<<G4endl;
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
                   if(-cLPDG==pRPDG)
                   {
                     nLPDG=pLPDG;
                     nRPDG=cRPDG;
                   }
                   else
                   {
                     if     (cRPDG > pRPDG) nRPDG=cRPDG*1000+pRPDG*100+1;
                     else if(cRPDG < pRPDG) nRPDG=pRPDG*1000+cRPDG*100+1;
                     else                   nRPDG=pRPDG*1000+cRPDG*100+3;
                     if  (-cLPDG == pL1)    nLPDG=pL2;
                     else                   nLPDG=pL1; // -cLPDG == pL2
                   }
                 }
                 else // cRPDG < 0
                 {
                   order=-1;
                   if(-cRPDG==pRPDG)
                   {
                     nLPDG=pLPDG;
                     nRPDG=cLPDG;
                   }
                   else
                   {
                     if     (cLPDG > pRPDG) nRPDG=cLPDG*1000+pRPDG*100+1;
                     else if(cLPDG < pRPDG) nRPDG=pRPDG*1000+cLPDG*100+1;
                     else                   nRPDG=pRPDG*1000+cLPDG*100+3;
                     if  (-cRPDG == pL1)    nLPDG=pL2;
                     else                   nLPDG=pL1; // -cRPDG == pL2
                   }
                 }
                 break;
               case 2: // ....................... pRPDG > 7
                 if(cLPDG < 0)
                 {
                   order=-1;
                   if(-cLPDG==pLPDG)
                   {
                     nLPDG=cRPDG;
                     nRPDG=pRPDG;
                   }
                   else
                   {
                     if     (cRPDG > pLPDG) nLPDG=cRPDG*1000+pLPDG*100+1;
                     else if(cRPDG < pLPDG) nLPDG=pLPDG*1000+cRPDG*100+1;
                     else                   nLPDG=pLPDG*1000+cRPDG*100+3;
                     if  (-cLPDG == pR1)    nRPDG=pR2;
                     else                   nRPDG=pR1; // -cLPDG == pR2
                   }
                 }
                 else // cRPDG < 0
                 {
                   order= 1;
                   if(-cRPDG==pLPDG)
                   {
                     nLPDG=cLPDG;
                     nRPDG=pRPDG;
                   }
                   else
                   {
                     if     (cLPDG > pLPDG) nLPDG=cLPDG*1000+pLPDG*100+1;
                     else if(cLPDG < pLPDG) nLPDG=pLPDG*1000+cLPDG*100+1;
                     else                   nLPDG=pLPDG*1000+cLPDG*100+3;
                     if  (-cRPDG == pR1)    nRPDG=pR2;
                     else                   nRPDG=pR1; // -cRPDG == pR2
                   }
                 }
                 break;
               case 3: // ....................... pLPDG <-7
                 if(cLPDG > 0)
                 {
                   order= 1;
                   if(cLPDG==-pRPDG)
                   {
                     nLPDG=pLPDG;
                     nRPDG=cRPDG;
                   }
                   else
                   {
                     if     (cRPDG < pRPDG) nRPDG=cRPDG*1000+pRPDG*100-1;
                     else if(cRPDG > pRPDG) nRPDG=pRPDG*1000+cRPDG*100-1;
                     else                   nRPDG=pRPDG*1000+cRPDG*100-3;
                     if  ( cLPDG == pL1)    nLPDG=-pL2;
                     else                   nLPDG=-pL1; // cLPDG == pL2
                   }
                 }
                 else // cRPDG > 0
                 {
                   order=-1;
                   if(cRPDG==-pRPDG)
                   {
                     nLPDG=pLPDG;
                     nRPDG=cLPDG;
                   }
                   else
                   {
                     if     (cLPDG < pRPDG) nRPDG=cLPDG*1000+pRPDG*100-1;
                     else if(cLPDG > pRPDG) nRPDG=pRPDG*1000+cLPDG*100-1;
                     else                   nRPDG=pRPDG*1000+cLPDG*100-3;
                     if  ( cRPDG == pL1)    nLPDG=-pL2;
                     else                   nLPDG=-pL1; // cRPDG == pL2
                   }
                 }
                 break;
               case 4: // ....................... pRPDG <-7
                 if(cLPDG > 0)
                 {
                   order=-1;
                   if(cLPDG==-pLPDG)
                   {
                     nLPDG=cRPDG;
                     nRPDG=pRPDG;
                   }
                   else
                   {
                     if     (cRPDG < pLPDG) nLPDG=cRPDG*1000+pLPDG*100-1;
                     else if(cRPDG > pLPDG) nLPDG=pLPDG*1000+cRPDG*100-1;
                     else                   nLPDG=pLPDG*1000+cRPDG*100-3;
                     if  ( cLPDG == pR1)    nRPDG=-pR2;
                     else                   nRPDG=-pR1; // cLPDG == pR2
                   }
                 }
                 else // cRPDG > 0
                 {
                   order= 1;
                   if(cRPDG==-pLPDG)
                   {
                     nLPDG=cLPDG;
                     nRPDG=pRPDG;
                   }
                   else
                   {
                     if     (cLPDG < pLPDG) nLPDG=cLPDG*1000+pLPDG*100-1;
                     else if(cLPDG > pLPDG) nLPDG=pLPDG*1000+cLPDG*100-1;
                     else                   nLPDG=pLPDG*1000+cLPDG*100-3;
                     if  ( cRPDG == pR1)    nRPDG=-pR2;
                     else                   nRPDG=-pR1; // cRPDG == pR2
                   }
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
               case 3: // ....................... cLPDG <-7 (cRPDG <0)
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
               case 4: // ....................... cRPDG <-7 (cLPDG <0)
                 if(pLPDG > 0)                       // pRPDG & cLPDG are anti-quarks
                 {
                   order=-1;
                   if     (pRPDG < cLPDG) nRPDG=pRPDG*1000+cLPDG*100-1;
                   else if(pRPDG > cLPDG) nRPDG=cLPDG*1000+pRPDG*100-1;
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
          if(!order) G4cerr<<"-Warning-G4QInel::Constr: t="<<tf<<", a="<<af<<", cL="<<cLPDG
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
#ifdef debug
                G4cout<<"G4QIonIonCollis::Constr:aLPDG="<<aLPDG<<", aRPDG="<<aRPDG<<G4endl;
#endif
                sing=false;
                G4QPDGCode tmp;
                std::pair<G4int,G4int> pB=tmp.MakeTwoBaryons(nL1, nL2, nR1, nR2);
                minM=G4QPDGCode(pB.first).GetMass()+G4QPDGCode(pB.second).GetMass();
              }
            }
            if(sing)
            {
              std::pair<G4int,G4int> newPair = std::make_pair(nLPDG,nRPDG);
              G4QContent newStQC(newPair);        // NewString QuarkContent
#ifdef debug
              G4cout<<"G4QIn::Con: LPDG="<<nLPDG<<",RPDG="<<nRPDG<<",QC="<<newStQC<<G4endl;
#endif
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
        G4cout<<"G4QIonIonCollis::Constr:cS4M="<<cS4M<<" fused/w pS4M="<<pL4M+pR4M<<G4endl;
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
        strings.erase(ist);
#ifdef debug
        G4LorentzVector ss4M=pL4M+pR4M;
        G4cout<<"G4QIonIonCollision::Constr:Created,4M="<<ss4M<<",m2="<<ss4M.m2()<<G4endl;
#endif
        if( ist != strings.begin() ) // To avoid going below begin() (for Windows)
        {
          ist--;
          con=false;
#ifdef debug
          G4cout<<"G4QIonIonCollision::Construct: *** IST Decremented ***"<<G4endl;
#endif
        }
        else
        {
          con=true;
#ifdef debug
          G4cout<<"G4QIonIonCollision::Construct: *** IST Begin ***"<<G4endl;
#endif
          break;
        }
      } // End of the IF(the best partnerString candidate was found)
      else
      {
#ifdef debug
        G4cout<<"-Warning-G4QInel::Const: S4M="<<cS4M<<",M2="<<cSM2<<" Leave ASIS"<<G4endl;
#endif
        ++problem;
        con=false;
      }
    }
    else con=false;
   } // End of loop over ist iterator
#ifdef debug
   G4cout<<"G4QIonIonCollision::Construct: *** IST While *** , con="<<con<<G4endl;
#endif
  } // End of "con" while 
#ifdef edebug
  // This print has meaning only if something appear between it and the StringFragmLOOP
  G4LorentzVector t4M=theProjNucleus.Get4Momentum()+theTargNucleus.Get4Momentum();//NucInLS
  G4int rC=totChg-theProjNucleus.GetZ()-theTargNucleus.GetZ();
  G4int rB=totBaN-theProjNucleus.GetA()-theTargNucleus.GetA();
  G4int nStr=strings.size();
  G4cout<<"-EMCLS-G4QIn::Const: AfterSUPPRESION #ofS="<<nStr<<",tNuc4M(E=M)="<<sum<<G4endl;
  for(G4int i=0; i<nStr; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    t4M+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    rC-=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    rB-=sBaN;
    G4cout<<"-EMCLS-G4QIonIonCollision::Construct: St#"<<i<<", 4M="<<strI4M<<", M="
          <<strI4M.m()<<", C="<<sChg<<", B="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QInel::Construct: r4M="<<t4M-totLS4M<<", rC="<<rC<<", rB="<<rB<<G4endl;
#endif
  //
  // --- If a problem is foreseen then the DiQaDiQ strings should be reduced if possible --
  //
#ifdef debug
    G4cout<<"G4QIonIonCollision::Construct: problem="<<problem<<G4endl;
#endif
  if(problem)
  {
    G4int nOfStr=strings.size();
#ifdef debug
    G4cout<<"G4QIonIonCollision::Constr:SecurityDiQaDiQReduction, #OfStr="<<nOfStr<<G4endl;
#endif
    for (G4int astring=0; astring < nOfStr; astring++)
    {
      G4QString* curString=strings[astring];
      G4QParton* cLeft=curString->GetLeftParton();
      G4QParton* cRight=curString->GetRightParton();
      G4int LT=cLeft->GetType();
      G4int RT=cRight->GetType();
      G4int sPDG=cLeft->GetPDGCode();
      G4int nPDG=cRight->GetPDGCode();
      if(LT==2 && RT==2)
      {
#ifdef debug
        G4cout<<"G4QIonIonCollis::Constr:TrySelfReduString,L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
        if( cLeft->ReduceDiQADiQ(cLeft, cRight) ) // DiQ-aDiQ pair was successfully reduced
        {
          sPDG=cLeft->GetPDGCode();
          nPDG=cRight->GetPDGCode();
#ifdef debug
          G4cout<<"+G4QInel::Const:#"<<astring<<" Reduced, L="<<sPDG<<", R="<<nPDG<<G4endl;
#endif
        }
#ifdef debug
        else G4cout<<"--*--G4QInel::Const: #"<<astring<<" DQ-aDQ reduction Failed"<<G4endl;
#endif
      } // End of the found DiQ/aDiQ pair
      else if(sPDG==3 && nPDG==-3)
      {
        sPDG= 1;
        nPDG=-1;
        cLeft->SetPDGCode(sPDG);
        cRight->SetPDGCode(nPDG);
      }
      else if(sPDG==-3 && nPDG==3)
      {
        sPDG=-1;
        nPDG= 1;
        cLeft->SetPDGCode(sPDG);
        cRight->SetPDGCode(nPDG);
      }
    }
    SwapPartons();
  } // End of IF(problem)
#ifdef edebug
  G4LorentzVector u4M=theProjNucleus.Get4Momentum()+theTargNucleus.Get4Momentum();//NucInLS
  G4int rCh=totChg-theProjNucleus.GetZ()-theTargNucleus.GetZ();
  G4int rBa=totBaN-theProjNucleus.GetA()-theTargNucleus.GetA();
  G4int nStri=strings.size();
  G4cout<<"-EMCLS-G4QIn::Const: FinalConstruct, #ofSt="<<nStri<<",tN4M(E=M)="<<t4M<<G4endl;
  for(G4int i=0; i<nStri; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    u4M+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    rCh-=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    rBa-=sBaN;
    G4cout<<"-EMCLS-G4QIonIonCollision::Construct: St#"<<i<<", 4M="<<strI4M<<", M="
          <<strI4M.m()<<", C="<<sChg<<", B="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QInel::Construct: r4M="<<u4M-totLS4M<<",rC="<<rCh<<",rB="<<rBa<<G4endl;
#endif
}

G4QIonIonCollision::~G4QIonIonCollision()
{
  std::for_each(strings.begin(), strings.end(), DeleteQString() );
}

G4QHadronVector* G4QIonIonCollision::Fragment()
{ // This is the member function fragmenting Strings & Quasmons (in nuclear matter)
#ifdef debug
  G4cout<<"*******>G4QIonIonCollision::Fragment: ***Called***, Res="<<theResult<<G4endl;
#endif
  G4int striNum=strings.size();                             // Find out if there're strings
  G4int hadrNum=theResult->size();                          // Find out if there're hadrons
#ifdef edebug
  G4int nQm=theQuasmons.size();
  G4LorentzVector totLS4M=theProjNucleus.Get4Momentum()+theTargNucleus.Get4Momentum(); //LS
  G4int totChg=theProjNucleus.GetZ()+theTargNucleus.GetZ();
  G4int totBaN=theTargNucleus.GetA()+theTargNucleus.GetA();
  G4cout<<"-EMCLS-G4QIn::Fragment: CHECKRecovery, #ofS="<<striNum<<", Nuc4M(E=M)="<<totLS4M
        <<",#Q="<<nQm<<",#H="<<hadrNum<<G4endl;
  for(G4int i=0; i < striNum; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    totLS4M+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    totChg+=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    totBaN+=sBaN;
    G4cout<<"-EMCLS-G4QIonIonCollision::Fragment:String#"<<i<<", 4M="<<strI4M<<", M="
          <<strI4M.m()<<", C="<<sChg<<", B="<<sBaN<<G4endl;
  }
  for(G4int i=0; i < nQm; i++)
  {
    G4LorentzVector hI4M=theQuasmons[i]->Get4Momentum();
    totLS4M+=hI4M;
    G4int hChg=theQuasmons[i]->GetCharge();
    totChg+=hChg;
    G4int hBaN=theQuasmons[i]->GetBaryonNumber();
    totBaN+=hBaN;
    G4cout<<"-EMCLS-G4QIonIonCollision::Fragment: Quasmon#"<<i<<", 4M="<<hI4M<<", C="<<hChg
          <<", B="<<hBaN<<G4endl;
  }
  for(G4int i=0; i < hadrNum; i++)
  {
    G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
    totLS4M+=hI4M;
    G4int hChg=(*theResult)[i]->GetCharge();
    totChg+=hChg;
    G4int hBaN=(*theResult)[i]->GetBaryonNumber();
    totBaN+=hBaN;
    G4cout<<"-EMCLS-G4QIn::Fragment:H#"<<i<<",4M="<<hI4M<<",C="<<hChg<<",B="<<hBaN<<G4endl;
  }
#endif
#ifdef debug
  G4cout<<"***>G4QIonIonCollision::Fragm: #OfStr="<<striNum<<", #OfRes="<<hadrNum<<G4endl;
#endif
  if(!striNum && hadrNum)                                   // Quasi-elastic or decoupled p
  {
#ifdef debug
    G4cout<<"***>G4QIonIonCollision::Fragm:**Quasi-Elastic**, #OfResult="<<hadrNum<<G4endl;
#endif
    return theResult;
  }
  else if(striNum) Breeder();                               // Strings fragmentation
  else                                                      // No strings, make HadrNucleus
  {
    if(hadrNum)
    {
      for(G4int ih=0; ih<hadrNum; ih++) delete (*theResult)[ih];
      theResult->clear();
    }
    G4LorentzVector rp4M=theProjNucleus.Get4Momentum();     // Nucleus 4-momentum in LS
    G4int rpPDG=theProjNucleus.GetPDG();                    // Nuclear PDG
    G4QHadron* resPNuc = new G4QHadron(rpPDG,rp4M);         // Nucleus -> Hadron
    theResult->push_back(resPNuc);                          // Fill the residual nucleus
    G4LorentzVector rt4M=theTargNucleus.Get4Momentum();     // Nucleus 4-momentum in LS
    G4int rtPDG=theTargNucleus.GetPDG();                    // Nuclear PDG
    G4QHadron* resTNuc = new G4QHadron(rtPDG,rt4M);         // Nucleus -> Hadron
    theResult->push_back(resTNuc);                          // Fill the residual nucleus
  }
  G4int nQuas=theQuasmons.size();                           // Size of the Quasmon OUTPUT
  G4int theRS=theResult->size();                            // Size of Hadron Output by now
#ifdef debug
  G4cout<<"***>G4QIonIonCollision::Fragm:beforeEnv, #OfQ="<<nQuas<<",#OfR="<<theRS<<G4endl;
#endif
  if(nQuas && theRS)
  {

    G4QHadron* resNuc = (*theResult)[theRS-1];              // Pointer to Residual Nucleus
    G4LorentzVector resNuc4M = resNuc->Get4Momentum();      // 4-Momentum of the Nucleuz
    G4int           resNucPDG= resNuc->GetPDGCode();        // PDG Code of the Nucleus
    G4QNucleus      theEnv(resNucPDG);                      // NucleusHadron->NucleusAtRest
    delete resNuc;                                          // Delete resNucleus as aHadron
    theResult->pop_back();                                  // Exclude the nucleus from HV
    --theRS;                                                // Reduce the OUTPUT by theNucl
#ifdef pdebug
    G4cout<<"G4QIonIonCollision::Fragm:#OfRemainingHadron="<<theRS<<", A="<<theEnv<<G4endl;
#endif
    // Now we need to be sure that the compound nucleus is heavier than the Ground State
    for(G4int j=theRS-1; j>-2; --j)                         // try to reach M_compound>M_GS
    {
      G4LorentzVector qsum4M=resNuc4M;                      // Proto compound 4-momentum
      G4QContent qsumQC=theEnv.GetQCZNS();                  // Proto compound Quark Content
#ifdef pdebug
      G4cout<<"G4QIonIonCollis::Fragm: rN4M"<<qsum4M<<qsum4M.m()<<",rNQC="<<qsumQC<<G4endl;
#endif
      G4Quasmon* firstQ=0;                                  // Prototype of theFirstQuasmon
      G4LorentzVector first4M;                              // Proto of the FirstQuasmon 4M
      G4QContent firstQC;                                   // Proto of the FirstQuasmon QC
      for(G4int i=0; i<nQuas; ++i)                          // LOOP over Quasmons
      {
        G4Quasmon* curQuasm=theQuasmons[i];                 // current Quasmon
        G4LorentzVector cur4M=curQuasm->Get4Momentum();     // 4-Mom of the Quasmon
        G4QContent curQC=curQuasm->GetQC();                 // Quark Content of the Quasmon
        qsum4M+=cur4M;                                      // Add quasmon's 4-momentum
        qsumQC+=curQC;                                      // Add quasmon's Quark Content
#ifdef pdebug
        G4cout<<"G4QIonIonCollis::Fragm: Q#"<<i<<", QQC="<<curQC<<", sQC="<<qsumQC<<G4endl;
#endif
        if(!i)                                              // Remember 1-st for correction
        {
          firstQ =curQuasm;
          first4M=cur4M;
          firstQC=curQC;
        }
      }
      G4int miPDG=qsumQC.GetSPDGCode();                     // PDG of minM of hadron/fragm.
      G4double gsM=0.;                                      // Proto minM of had/frag forQC
      if(miPDG == 10)
      {
        G4QChipolino QCh(qsumQC);                           // define TotNuc as a Chipolino
        gsM=QCh.GetQPDG1().GetMass()+QCh.GetQPDG2().GetMass(); // Sum of Hadron Masses
        //gsM=theWorld->GetQParticle(QCh.GetQPDG1())->MinMassOfFragm() +
        //    theWorld->GetQParticle(QCh.GetQPDG2())->MinMassOfFragm();
      }
      // @@ it is not clear, why it does not work ?!
      //else if(miPDG>80000000)                             // Compound Nucleus
      //{
      //  G4QNucleus rtN(qsumQC);                           // CreatePseudoNucl for totComp
      //  gsM=rtN.GetGSMass();                              // MinMass of residQ+(Env-ParC)
      //}
      else if(miPDG < 80000000 && std::abs(miPDG)%10 > 2)
                           gsM=theWorld->GetQParticle(G4QPDGCode(miPDG))->MinMassOfFragm();
      else gsM=G4QPDGCode(miPDG).GetMass();                 // minM of hadron/fragm. for QC
      G4double reM=qsum4M.m();                              // real mass of the compound
#ifdef pdebug
      G4cout<<"G4QIonIonCollision::Fragm: PDG="<<miPDG<<", rM="<<reM<<",GSM="<<gsM<<G4endl;
#endif
      if(reM > gsM) break;                                  // CHIPS can be called
      if(j > -1)                                            // Can try to add hadrons to Q0
      {
        G4QHadron* cH = (*theResult)[j];                    // Pointer to the last Hadron
        G4LorentzVector h4M = cH->Get4Momentum();           // 4-Momentum of the Hadron
        G4QContent      hQC = cH->GetQC();                  // QC of the Hadron
        firstQ->Set4Momentum(first4M+h4M);                  // Update the Quasmon's 4-Mom
        firstQ->SetQC(firstQC+hQC);                         // Update the Quasmon's QCont
        delete cH;                                          // Delete the Hadron
        theResult->pop_back();                              // Exclude the hadron from HV
#ifdef pdebug
        G4cout<<"G4QInel::Fragm: H#"<<j<<", hQC="<<hQC<<",hPDG="<<cH->GetPDGCode()<<G4endl;
#endif
      }
      else
      {
        G4cerr<<"***G4QIonIonCollision::Fra:PDG="<<miPDG<<",M="<<reM<<",GSM="<<gsM<<G4endl;
        G4Exception("G4QIonIonCollision::Fragment:","27",FatalException,"Can'tRecoverGSM");
      }
    }
    G4double nucE=resNuc4M.e();                             // Total energy of the nuclEnv
    if(nucE<1.E-12) nucE=0.;                                // Computer accuracy safety
    G4ThreeVector   nucVel(0.,0.,0.);                       // Proto of the NucleusVelocity
    G4QHadronVector* output=0;                              // NucleusFragmentation Hadrons
    G4QEnvironment* pan= new G4QEnvironment(theEnv);        // ---> DELETED --->----------+
#ifdef pdebug
    G4cout<<"G4QIonIonCollision::Fragm: nucE="<<nucE<<", nQ="<<nQuas<<G4endl; //          |
#endif
    if(nucE) nucVel=resNuc4M.vect()/nucE;                   // The NucleusVelocity        |
    for(G4int i=0; i<nQuas; ++i)                            // LOOP over Quasmons         |
    {                                                       //                            |
      G4Quasmon* curQuasm=theQuasmons[i];                   // current Quasmon            |
      if(nucE) curQuasm->Boost(-nucVel);                    // Boost it to CMS of Nucleus |
      pan->AddQuasmon(curQuasm);                            // Fill the predefined Quasmon|
#ifdef pdebug
      G4LorentzVector cQ4M=curQuasm->Get4Momentum();        // Just for printing          |
      G4cout<<"G4QIonIonCollision::Fragment: Quasmon# "<<i<<" added, 4M="<<cQ4M<<G4endl;//|
#endif
    }                                                       //                            |
    try                                                     //                            |
    {                                                       //                            |
      delete output;                                        //                            |
      output = pan->Fragment();// DESTROYED after theHadrons are transferred to theResult |
    }                                                       //                          | |
    catch (G4QException& error)                             //                          | |
    {                                                       //                          | |
      G4cerr<<"***G4QIonIonCollision::Fragment: G4QE Exception is catched"<<G4endl; //  | |
      G4Exception("G4QIonIonCollision::Fragment:","27",FatalException,"CHIPSCrash");//  | |
    }                                                       //                          | |
    delete pan;                              // Delete the Nuclear Environment <-----<--+-+
    if(output)                               // Output exists                           |
    {                                        //                                         |
      G4int nOut=output->size();             // #ofHadrons in the Nuclear Fragmentation |
      for(G4int j=0; j<nOut; j++)            // LOOP over Hadrons transferring to LS    |
      {                                      //                                         |
        G4QHadron* curHadron=(*output)[j];   // Hadron from the nucleus fragmentation   |
        if(nucE) curHadron->Boost(nucVel);   // Boost it back to Laboratory System      |
        theResult->push_back(curHadron);     // Transfer it to the result               |
      }                                      //                                         |
      delete output;                         // Delete the OUTPUT <-----<-----<-----<---+
    }
  }
  else if(!striNum) G4cout<<"-Warning-G4QIonIonCollision::Fragment:NothingWasDone"<<G4endl;
#ifdef debug
  G4cout<<"====>G4QIonIonCollision::Fragment: Final #OfResult="<<theResult->size()<<G4endl;
#endif
    G4int nQ =theQuasmons.size();
    if(nQ) theQuasmons.clear();                              // @@ Not necesary ?
#ifdef edebug
    G4LorentzVector f4M(0.,0.,0.,0.);                        // Sum of the Result in LS
    G4int fCh=totChg;
    G4int fBN=totBaN;
    G4int nHd=theResult->size();
    G4cout<<"-EMCLS-G4QIonIonCollision::Fragment: #ofHadr="<<nHd<<",#OfQuasm="<<nQ<<G4endl;
    for(G4int i=0; i<nHd; i++)
    {
      G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
      f4M+=hI4M;
      G4int hChg=(*theResult)[i]->GetCharge();
      fCh-=hChg;
      G4int hBaN=(*theResult)[i]->GetBaryonNumber();
      fBN-=hBaN;
      G4cout<<"-EMCLS-G4QIonIonCollision::Fragment: Hadron#"<<i<<", 4M="<<hI4M<<", PDG="
            <<(*theResult)[i]->GetPDGCode()<<", C="<<hChg<<", B="<<hBaN<<G4endl;
    }
    G4cout<<"-EMCLS-G4QInel::Fragm:, r4M="<<f4M-totLS4M<<", rC="<<fCh<<",rB="<<fBN<<G4endl;
#endif
  return theResult;
} // End of fragmentation

void G4QIonIonCollision::Breeder()
{ // This is the member function, which returns the resulting vector of Hadrons & Quasmons
  static const G4double  eps = 0.001;                              // Tolerance in MeV
  //
  // ------------ At this point the strings are fragmenting to hadrons in LS -------------
  //
#ifdef edebug
  G4LorentzVector totLS4M=theProjNucleus.Get4Momentum()+theTargNucleus.Get4Momentum(); //LS
  G4int totChg=theProjNucleus.GetZ()+theTargNucleus.GetZ();
  G4int totBaN=theProjNucleus.GetA()+theTargNucleus.GetA();
  G4int nStri=strings.size();
  G4cout<<"-EMCLS-G4QIn::Breed: CHECKRecovery #ofS="<<nStri<<",N4M(E=M)="<<totLS4M<<G4endl;
  for(G4int i=0; i<nStri; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    totLS4M+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    totChg+=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    totBaN+=sBaN;
    G4cout<<"-EMCLS-G4QIonIonCollision::Breeder: St#"<<i<<", 4M="<<strI4M<<", M="
          <<strI4M.m()<<", C="<<sChg<<", B="<<sBaN<<G4endl;
  }
#endif
  G4int nOfStr=strings.size();
#ifdef debug
  G4cout<<"G4QIonIonCollision::Breeder: BeforeFragmentation, #OfStr="<<nOfStr<<G4endl;
#endif
  G4LorentzVector ft4M(0.,0.,0.,0.);
  G4QContent      ftQC(0,0,0,0,0,0);
  G4bool          ftBad=false;
  for(G4int i=0; i < nOfStr; ++i)
  {
    G4LorentzVector pS4M=strings[i]->Get4Momentum(); // String 4-momentum
    ft4M+=pS4M;
    G4QContent pSQC=strings[i]->GetQC();             // String Quark Content
    ftQC+=pSQC;
    if(pS4M.m2() < 0.) ftBad=true;
#ifdef debug
    G4cout<<">G4QIonIonCollision::Breed:1stTest,S#"<<i<<",4M="<<pS4M<<",QC="<<pSQC<<G4endl;
#endif
  }
  if(ftBad)
  {
    G4Quasmon* stringQuasmon = new G4Quasmon(ftQC, ft4M);
#ifdef debug
    G4cout<<"->G4QIonIonCollision::Breed:*TotQ*,QC="<<ftQC<<",4M="<<ft4M<<ft4M.m()<<G4endl;
#endif
    theQuasmons.push_back(stringQuasmon);
    G4LorentzVector rp4M=theProjNucleus.Get4Momentum(); // Nucleus 4-momentum in LS
    G4int rpPDG=theProjNucleus.GetPDG();
    G4QHadron* resPNuc = new G4QHadron(rpPDG,rp4M);
    theResult->push_back(resPNuc);                  // Fill the residual projectile nucleus
    G4LorentzVector rt4M=theTargNucleus.Get4Momentum(); // Nucleus 4-momentum in LS
    G4int rtPDG=theTargNucleus.GetPDG();
    G4QHadron* resTNuc = new G4QHadron(rtPDG,rt4M);
    theResult->push_back(resTNuc);                  // Fill the residual target nucleus
    return;
  }
  for (G4int astring=0; astring < nOfStr; astring++)
  {
#ifdef edebug
    G4LorentzVector sum=theProjNucleus.Get4Momentum()+theTargNucleus.Get4Momentum(); //inLS
    G4int rChg=totChg-theProjNucleus.GetZ()-theTargNucleus.GetZ();
    G4int rBaN=totBaN-theProjNucleus.GetA()-theTargNucleus.GetA();
    G4int nOfHadr=theResult->size();
    G4cout<<"-EMCLS-G4QIonIonCollision::Breed:#ofSt="<<nOfStr<<",#ofHad="<<nOfHadr<<G4endl;
    for(G4int i=astring; i<nOfStr; i++)
    {
      G4LorentzVector strI4M=strings[i]->Get4Momentum();
      sum+=strI4M;
      G4int sChg=strings[i]->GetCharge();
      rChg-=sChg;
      G4int sBaN=strings[i]->GetBaryonNumber();
      rBaN-=sBaN;
      G4cout<<"-EMCLS-G4QI::Breed:S#"<<i<<",4M="<<strI4M<<",C="<<sChg<<",B="<<sBaN<<G4endl;
    }
    for(G4int i=0; i<nOfHadr; i++)
    {
      G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
      sum+=hI4M;
      G4int hChg=(*theResult)[i]->GetCharge();
      rChg-=hChg;
      G4int hBaN=(*theResult)[i]->GetBaryonNumber();
      rBaN-=hBaN;
      G4cout<<"-EMCLS-G4QIn::Breed: H#"<<i<<",4M="<<hI4M<<",C="<<hChg<<",B="<<hBaN<<G4endl;
    }
    G4cout<<"....-EMCLS-G4QInel::Br:r4M="<<sum-totLS4M<<",rC="<<rChg<<",rB="<<rBaN<<G4endl;
#endif
    G4QString* curString=strings[astring];
    if(!curString->GetDirection()) continue;  // Historic for the dead strings: DoesNotWork
#ifdef edebug
    G4int curStrChg = curString->GetCharge();
    G4int curStrBaN = curString->GetBaryonNumber();
#endif
    G4LorentzVector curString4M = curString->Get4Momentum();
#ifdef debug
    G4cout<<"====>G4QIonIonCollision::Breeder: String#"<<astring<<",s4M/m="<<curString4M
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
        G4cout<<"G4QIonIonCollision::Breed:TryReduceString, L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
        if( cLeft->ReduceDiQADiQ(cLeft, cRight) ) // DiQ-aDiQ pair was successfully reduced
        {
          LT=1;
          RT=1;
          LS=2;
          sPDG=cLeft->GetPDGCode();
          nPDG=cRight->GetPDGCode();
#ifdef debug
          G4cout<<"G4QIonIonCollision::Breed:AfterReduction,L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
          theHadrons=curString->FragmentString(true);//!! Try to fragment the new String !!
          cLeft=curString->GetLeftParton();
          cRight=curString->GetRightParton();
#ifdef debug
          G4cout<<"G4QInel::Breed:L="<<cLeft->Get4Momentum()<<",R="<<cRight->Get4Momentum()
                <<G4endl;
#endif
        }
#ifdef debug
        else G4cout<<"^G4QIonIonCollision::Breed: DQ-aDQ reduction to Q-aQ Failed"<<G4endl;
#endif
      } // End of the SelfReduction
#ifdef debug
      G4cout<<"G4QIonIonCollision::Breed:AfterRedAttempt, theH="<<theHadrons<<", L4M="
            <<cLeft->Get4Momentum()<<", R4M="<<cRight->Get4Momentum()<<G4endl;
#endif
      unsigned next=astring+1;                 // The next string position
      if (!theHadrons)                         // The string can not be fragmented
      {
        G4int fusionDONE=0; // StringFusion didn't happen (1=Fuse L+L/R+R, -1=Fuse L+R/R+L)
        if(next < strings.size())             // TheString isn't theLastString can fuse
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
          G4cout<<"G4QIonIonCollision::Breed:TryFuseStringS, q="<<qPDG<<", a="<<dPDG
                <<", n="<<next<<G4endl;
#endif
          G4ThreeVector curV=curString4M.vect()/curString4M.e();
          G4int reduce=0;                      // a#of reduced Q-aQ pairs
          G4int restr=0;                       // To use beyon the LOOP for printing
          G4int MPS=0;                         // PLS for the selected string
          for (restr=next; restr < nOfStr; restr++)
          {
            reduce=0;
            G4QString* reString=strings[restr];
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
            G4cout<<"G4Qnel::Breed: TryReduce #"<<restr<<", q="<<rPDG<<",a="<<aPDG<<G4endl;
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
		  G4cout<<"G4QIonIonCollision::Breeder: cQ="<<cQ1<<","<<cQ2<<", cA="<<cA1<<","
                    <<cA2<<", pQ="<<pQ1<<","<<pQ2<<", pA="<<pA1<<","<<pA2<<G4endl;
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
                    G4cerr<<"*G4QIonIonCollision::Breed:CL="<<newCL<<",PL="<<newPL<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"2LL-");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair((-nPDG)/100, mPDG/100);
                  G4int newCR=resRR.first;
                  G4int newPR=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QIonIonCollision::Breed:CR="<<newCR<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"2RR-");
                  }
                  cLeft->SetPDGCode(newCL);    // Reset the left quark of curString
                  cRight->SetPDGCode(-newCR);  // Reset the right quark of curString
                  Left->SetPDGCode(-newPL);    // Reset the left quark of reString
                  Right->SetPDGCode(newPR);    // Reset the right quark of reString
                  break;                       // Break out of the reString internal LOOP
                }
                else if(sPDG<0 && uPDG>0)           // LL/RR Reduction
                {
                  std::pair<G4int,G4int> resLL=ReducePair((-sPDG)/100, uPDG/100);
                  G4int newCL=resLL.first;
                  G4int newPL=resLL.second;
                  if(!newCL || !newPL)
                  {
                    G4cerr<<"*G4QIonIonCollision::Breed:CL="<<newCL<<",PL="<<newPL<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"2LL+");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair(nPDG/100, (-mPDG)/100);
                  G4int newCR=resRR.first;
                  G4int newPR=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QIonIonCollision::Breed:CR="<<newCR<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"2RR+");
                  }
                  cLeft->SetPDGCode(-newCL);   // Reset the left quark of curString
                  cRight->SetPDGCode(newCR);   // Reset the right quark of curString
                  Left->SetPDGCode(newPL);     // Reset the left quark of reString
                  Right->SetPDGCode(-newPR);   // Reset the right quark of reString
                  break;                       // Break out of the reString internal LOOP
                }
                else if(sPDG>0 && mPDG<0)                // LR Reduction
                {
                  std::pair<G4int,G4int> resLL=ReducePair(sPDG/100, (-mPDG)/100);
                  G4int newCL=resLL.first;
                  G4int newPR=resLL.second;
                  if(!newCL || !newPR)
                  {
                    G4cerr<<"*G4QIonIonCollision::Breed:CL="<<newCL<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"2-LR");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair((-nPDG)/100, uPDG/100);
                  G4int newCR=resRR.first;
                  G4int newPL=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QIonIonCollision::Breed:CR="<<newCR<<",PL="<<newPL<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"2-LR");
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
                    G4cerr<<"*G4QIonIonCollision::Breed:CL="<<newCL<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"2-RL");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair(nPDG/100, (-uPDG)/100);
                  G4int newCR=resRR.first;
                  G4int newPL=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QIonIonCollision::Breed:CR="<<newCR<<",PL="<<newPL<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"2-RL");
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
            G4cout<<"G4QInel::Breed: TryFuse/w #"<<restr<<",q="<<rPDG<<",a="<<aPDG<<G4endl;
#endif
            G4int PLS=PLT+PRT;
            if( (LS==2 && PLS==2) ||                           // QaQ+QaQ always to DiQaDiQ
                ( ( (LS==2 && PLS==3) || (LS==3 && PLS==2) ) &&// QaQ w QDiQ/aQaDiQ(single)
                  ( (aPDG> 7 && (-dPDG==aPDG/10   || -dPDG==aPDG%10) )   || // cAQ w DQ
                    (dPDG> 7 && (-aPDG==dPDG/10   || -aPDG==dPDG%10) )   || // AQ w cDQ
                    (rPDG<-7 && (qPDG==(-rPDG)/10 || qPDG==(-rPDG)%10) ) || // cQ w ADQ
                    (qPDG<-7 && (rPDG==(-qPDG)/10 || rPDG==(-qPDG)%10) )    // Q w cADQ
                    //|| (aPDG< 0 && -aPDG==qPDG) || (dPDG< 0 && -dPDG==rPDG) // aQ w Q 
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
              G4cout<<"G4QIonIonCollision::Breeder: StringCand#"<<restr<<", q="<<rPDG
                    <<", a="<<aPDG<<", L="<<uPDG<<", R="<<mPDG<<",dV="<<dV<<G4endl;
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
            G4cout<<"-G4QIonIonCollision::Breed:Reduced #"<<astring<<" & #"<<restr<<G4endl;
#endif
            astring--;                         // String was QCreduced using another String
            continue;                          // Repeat fragmentation of the same string
          }
          if(fustr)                            // The partner was found -> fuse strings
          {
#ifdef debug
            G4cout<<"G4QInel::Breeder: StPartner#"<<fustr<<", LT="<<LT<<",RT="<<RT<<G4endl;
#endif
            G4QString* fuString=strings[fustr];
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
            G4cout<<"G4QIonIonCollision::Breeder:BeforeFuDir,sL="<<sPDG<<",nR="<<nPDG
                  <<",uL="<<uPDG<<",mR="<<mPDG<<",L4M="<<L4M<<",R4M="<<R4M<<G4endl;
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
              G4cout<<"G4QIonIonCollision::Breeder:LL/RR s4M="<<fuString->Get4Momentum()
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
              G4cout<<"G4QIonIonCollision::Breeder:LR/RL s4M="<<fuString->Get4Momentum()
                    <<",S="<<L4M+PL4M+R4M+PR4M<<", L="<<Left->Get4Momentum()<<", R="
                    <<Right->Get4Momentum()<<G4endl;
#endif
            }
#ifdef debug
            else G4cout<<"-Warning-G4QIonIonCollision::Breeder: WrongStringFusion"<<G4endl;
#endif
#ifdef edebug
            G4cout<<"#EMC#G4QIonIonCollision::Breed:StringFused,F="<<fusionDONE<<",L="<<L4M
                  <<",R="<<R4M<<",pL="<<PL4M<<",pR="<<PR4M<<",nL="<<Left->Get4Momentum()
                  <<",nR="<<Right->Get4Momentum()<<",S="<<fuString->Get4Momentum()<<G4endl;
#endif
            if(fusionDONE)
            {
#ifdef debug
              G4cout<<"###G4QIonIonCollision::Breed: Str#"<<astring<<" fused/w Str#"<<fustr
                    <<"->"<<fuString->Get4Momentum()<<fuString->Get4Momentum().m()<<G4endl;
#endif
              continue;                          // @@ killing of the curString is excluded
            }
          }
#ifdef debug
          else
          {

            G4cerr<<"**G4QIonIonCollision::Breed:*NoPart*M="<<curString->Get4Momentum().m()
                  <<", F="<<fusionDONE<<", LPDG="<<curString->GetLeftParton()->GetPDGCode()
                  <<", RPDG="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
          }
#endif
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
              G4cout<<"-EMC-G4QInel::Breed:H#"<<i<<",4M="<<h4M<<hQC<<",PDG="<<hPDG<<G4endl;
            }
#endif
            G4int fusDONE=0;                      // String+Hadron fusion didn't happen
            G4int fuhad=-1;                       // The found hadron index
            G4int newPDG=0;                       // PDG ofTheParton afterMerging with Hadr
            G4int secPDG=0;                       // Second PDG oParton afterMerging w/Hadr
            G4double maM2=-DBL_MAX;               // Prototype of the max ResultingStringM2
            G4LorentzVector selH4M(0.,0.,0.,0.);  // 4-mom of the selected hadron
            G4QHadron* selHP=0;                   // Pointer to the used hadron for erasing
            G4QString* cString=strings[astring];  // Must be the last string by definition
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
              G4cout<<"G4QIonIonCollision::Breeder: Hadron # "<<reh<<", QC="<<curQC
                      <<", P1="<<partPDG1<<", P2="<<partPDG2<<G4endl;
#endif
              if(partPDG1 || partPDG2)            // Hadron can merge at least w/one parton
              {
                G4int cCur=1;
                if(sumT>3 && partPDG1 && partPDG2) cCur=2;
                G4LorentzVector curHadr4M = curHadr->Get4Momentum();
                G4double M2=(cString4M+curHadr4M).m2();// SquaredMass of theResultingString
#ifdef debug
                G4cout<<"G4QIonIonCollision::Breeder:*IN*Hadron#"<<reh<<",M2="<<M2<<G4endl;
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
                  G4cout<<"G4QInel::Br:*Selected*,P1="<<partPDG1<<",P2="<<partPDG2<<G4endl;
#endif
                  selH4M=curHadr4M;
                  selHP=curHadr;
                } // End of IF(update selection)
              } // End of IF(HadronCanMergeWithTheString)
            } // End of the LOOP over Hadrons
#ifdef debug
            G4cout<<"G4QIonIonCollision::Breeder: fuh="<<fuhad<<",fus="<<fusDONE<<G4endl;
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
                  G4cout<<"G4QInel::Br:TryReduce, nPDG="<<newPDG<<",sPDG="<<secPDG<<G4endl;
#endif
                  cLeft->ReduceDiQADiQ(cLeft, cRight);
                }
#ifdef debug
                G4cout<<"G4QIonIonCollision::Breed: Left, s4M="<<curString->Get4Momentum()
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
                G4cout<<"G4QIonIonCollision::Breed: Right, s4M="<<curString->Get4Momentum()
                      <<", L4M="<<cLeft->Get4Momentum()<<", R4M="<<newRight<<", h4M="
                      <<selH4M<<", newR="<<newPDG<<", oldL="<<cLeft->GetPDGCode()<<G4endl;
#endif
              }
#ifdef debug
              else G4cout<<"-G4QIonIonCollision::Breed: Wrong String+HadronFusion"<<G4endl;
#endif
#ifdef debug
              if(fusDONE) G4cout<<"####G4QIonIonCollision::Breeder: String #"<<astring
                                <<" is fused with Hadron #"<<fuhad
                                <<", new4Mom="<<curString->Get4Momentum()
                                <<", M2="<<curString->Get4Momentum().m2()
                                <<", QC="<<curString->GetQC()<<G4endl;
#endif
            }
            else
            {
#ifdef debug
              G4cerr<<"**G4QIonIonCollision::Breed:*NoH*,M="<<curString->Get4Momentum().m()
                    <<", LPDG="<<curString->GetLeftParton()->GetPDGCode()
                    <<", RPDG="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
              // @@ Temporary exception for the future solution
              //G4Exception("G4QIonIonCollision::Breed:","72",FatalException,"SHNotFused");
#endif
              break;                           // Breake the While LOOP
            } // End of the namespace where both Fusion and reduction have failed
            // The fused hadron must be excluded from theResult
#ifdef debug
            G4cout<<"G4QIonIonCollision::Breed: before HR, nH="<<theResult->size()<<G4endl;
            G4int icon=0;                              // Loop counter
#endif
            G4QHadronVector::iterator ih;
            G4bool found=false;
            for(ih = theResult->begin(); ih != theResult->end(); ih++)
            {
#ifdef debug
              G4cout<<"G4QInelast::Breeder:#"<<icon<<", i="<<(*ih)<<", sH="<<selHP<<G4endl;
              icon++;
#endif
              if((*ih)==selHP)
              {
#ifdef debug
                G4cout<<"G4QInel::Breed: *HadronFound*, PDG="<<selHP->GetPDGCode()<<G4endl;
#endif
                G4LorentzVector p4M=selHP->Get4Momentum();
                curString4M+=p4M;
#ifdef edebug
                G4int Chg=selHP->GetCharge();
                G4int BaN=selHP->GetBaryonNumber();
                curStrChg+=Chg;
                curStrBaN+=BaN;
                G4cout<<"-EMC->>>>G4QIonIonCollision::Breed: S+=H, 4M="<<curString4M<<",M="
                      <<curString4M.m()<<", Charge="<<curStrChg<<", B="<<curStrBaN<<G4endl;
#endif
                delete selHP;                          // delete the Hadron
                theResult->erase(ih);                  // erase the Hadron from theResult
                found=true;
                break;                                 // beak the LOOP over hadrons
              }
            } // End of the LOOP over hadrons
#ifdef debug
            if(!found) G4cout<<"*G4QIonIonCollision::Breed:nH="<<theResult->size()<<G4endl;
#endif
            // New attempt of the string decay
            theHadrons=curString->FragmentString(true);//! Try to fragment the new String !
#ifdef debug
            G4cout<<"G4QInel::Breeder: tH="<<theHadrons<<",nH="<<theResult->size()<<G4endl;
#endif
          } // End of the while LOOP over the fusion with hadrons
#ifdef debug
          G4cout<<"*G4QIonIonCollision::Breed: *CanTryToDecay?* nH="<<theHadrons<<", next="
                <<next<<" =? nS="<<strings.size()<<", nR="<<theResult->size()<<G4endl;
#endif
          if(!theHadrons && next == strings.size() && !(theResult->size()))// TryToDecay
          {
            G4QContent miQC=curString->GetQC();    // QContent of the Lightest Hadron
            G4int miPDG=miQC.GetSPDGCode();         // PDG of the Lightest Hadron
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
                    G4cerr<<"***G4QIonIonCollision::Breeder: tM="<<ttM<<"->h1="<<h1QPDG
                          <<"("<<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")="<<h1M+h2M<<G4endl;
                    G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"QDec");
                  }
                }
                G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
                theResult->push_back(h1H);         // (delete equivalent)  
#ifdef debug
                G4LorentzVector f4M=h1H->Get4Momentum();
                G4int           fPD=h1H->GetPDGCode();
                G4int           fCg=h1H->GetCharge();
                G4int           fBN=h1H->GetBaryonNumber();
                G4cout<<"-EMC->>G4QIonIonCollision::Breed:String=HadrChiPro1's filled,f4M="
                      <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
                G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
                theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
                G4LorentzVector s4M=h2H->Get4Momentum();
                G4int           sPD=h2H->GetPDGCode();
                G4int           sCg=h2H->GetCharge();
                G4int           sBN=h2H->GetBaryonNumber();
                G4cout<<"-EMC->>G4QIonIonCollision::Breed:String=HadrChiPro2's filled,s4M="
                      <<s4M<<", sPDG="<<sPD<<", sCg="<<sCg<<", sBN="<<sBN<<G4endl;
#endif
#ifdef edebug
                G4cout<<"-EMC-..Chi..G4QIonIonCollision::Breeder: DecayCHECK, Ch4M="
                      <<curString4M<<", d4M="<<curString4M-h14M-h24M<<G4endl;
#endif
                break;                               // Go out of the main StringDecay LOOP
              }
              else
              {
                G4Quasmon* stringQuasmon = new G4Quasmon(miQC, curString4M);// String->Quas
                theQuasmons.push_back(stringQuasmon);
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
              G4int tmpN=tmpQHadVec->size();
#ifdef debug
              G4cout<<"G4QIonIonCollision::Breeder: Decay the Last, Res#H="<<tmpN<<G4endl;
#endif
              if(tmpN>1)
              {
                for(G4int aH=0; aH < tmpN; aH++)
                {
                  theResult->push_back((*tmpQHadVec)[aH]);//TheDecayProdOfHadronDecIsFilled
#ifdef debug
                  G4QHadron*   prodH =(*tmpQHadVec)[aH];
                  G4LorentzVector p4M=prodH->Get4Momentum();
                  G4int           PDG=prodH->GetPDGCode();
                  G4int           Chg=prodH->GetCharge();
                  G4int           BaN=prodH->GetBaryonNumber();
                  G4cout<<"-EMC->>G4QIonIonCollis::Breed:String=Hadr,H#"<<aH<<" filled,4M="
                        <<p4M<<", PDG="<<PDG<<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
                }
              }
              else
              {
                G4Quasmon* stringQuasmon = new G4Quasmon(miQC, curString4M);// String->Quas
#ifdef debug
                G4cout<<"G4QIonIonCollision::Breeder:==> to Quasm="<<miQC<<curString4M
                      <<", pNuc="<<theProjNucleus<<theProjNucleus.Get4Momentum()<<", tNuc="
                      <<theTargNucleus<<theTargNucleus.Get4Momentum()<<", NString="
                      <<strings.size()<<", nR="<<theResult->size()<<", nQ="
                      <<theQuasmons.size()<<G4endl;
#endif
                theQuasmons.push_back(stringQuasmon);
                delete sHad;
                tmpQHadVec->clear();
                delete tmpQHadVec;  // WhoCallsDecayQHadron is responsible for clear&delete
                break;                               // Go out of the main StringDecay LOOP
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
      G4cout<<"G4QIonIonCollision::Breeder: theH="<<theHadrons<<"?=0, next="<<next<<G4endl;
#endif
      if(!theHadrons && next < strings.size())       // ForwardInLOOP strings exist
      {
        // @@ string can be not convertable to one hadron (2,0.0,0,2,0) - To be improved
        G4QContent miQC=curString->GetQC(); // QContent of the Lightest Hadron
        G4int miPDG=miQC.GetSPDGCode();// PDG of the Lightest Hadron
#ifdef debug
        G4cout<<">>>G4QIonIonCollision::Breeder: SQC="<<miQC<<", miSPDG="<<miPDG<<G4endl;
#endif
        G4double miM=0.;               // Prototype of the Mass of the Cur LightestHadron
        if(miPDG!=10) miM=G4QPDGCode(miPDG).GetMass(); // Mass of theCurLightestOneHadron
        else
        {
          G4QChipolino QCh(miQC);      // define the TotalString as a Chipolino
          miM=QCh.GetQPDG1().GetMass()+QCh.GetQPDG2().GetMass();//MinMass of theStringChipo
        }
        G4double cM2=curString4M.m2(); // Actual squared mass of the Cur String
#ifdef debug
        G4cout<<">>>G4QIonIonCollision::Breeder: minMass="<<miM<<", realM2="<<cM2<<G4endl;
#endif
        G4double   cM=0.;
        if(cM2>0.)
        {
          cM=std::sqrt(cM2);
          if(std::fabs(cM-miM) < eps)    // Convert to hadron(2 hadrons) w/o calculations
          {
            if(miPDG!=10)
            {
              G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
              theResult->push_back(sHad);// Fill the curString as a hadron
#ifdef debug
              G4cout<<">>>G4QIonIonCollision::Breeder:S->H="<<miPDG<<curString4M<<G4endl;
#endif
            }
            else
            {
              G4QChipolino QCh(miQC);               // define TotalResidual as a Chipolino
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
              G4cout<<"-EMC->>G4QIonIonCollision::Breed:Str=2HadrAR Prod-F is filled, f4M="
                    <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
              G4LorentzVector s4M=h2H->Get4Momentum();
              G4int           sPD=h2H->GetPDGCode();
              G4int           sCg=h2H->GetCharge();
              G4int           sBN=h2H->GetBaryonNumber();
              G4cout<<"-EMC->>G4QIonIonCollision::Breed:Str=2HadrAR Prod-S is filled, s4M="
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
            G4int restr=0;            // To use beyon the LOOP for printing
            G4int fustr=0;               // Selected String-partner (0 = NotFound)
            G4double selX=0.;            // Selected value of x
            G4double maD=-DBL_MAX;       // Maximum Free Mass
            G4double Vmin=DBL_MAX;       // min Velocity Distance
            G4LorentzVector s4M(0.,0.,0.,0.); // Selected 4-mom of the hadron
#ifdef debug
            G4cout<<"G4QIonIonCollision::Breeder: TryRecover, cV="<<curV<<G4endl;
#endif
            nOfStr=strings.size();
            for(restr=next; restr < nOfStr; ++restr) if(restr != astring)
            {
              G4QString* pString=strings[restr];
              G4LorentzVector p4M=pString->Get4Momentum();
              G4ThreeVector pP=p4M.vect();  // Momentum of the partnerString
              G4double      pE=p4M.e();     // Energy of the partnerString
              G4double D2=cE*pE-cP.dot(pP); 
              G4double pM2=p4M.m2();
              G4double dM4=pM2*(miM2-cM2);
              G4double D=D2*D2+dM4;
              G4double x=-1.;                // Bad preexpectation 
              if(D >= 0. && pM2>.01) x=(std::sqrt(D)-D2)/pM2; // what we should get from p
#ifdef debug
		  else G4cout<<"G4QIonIonCollis::Breed:D="<<D<<",D2="<<D2<<",dM="<<dM4<<G4endl;
              G4cout<<"G4QIonIonCollision::Breed: pM2="<<pM2<<",D2="<<D2<<",x="<<x<<G4endl;
#endif
              if(x > 0. && x < 1.)          // We are getting x part of p4M
              {
                G4QContent pQC=pString->GetQC(); // Quark Content of The Partner
                G4int pPDG=pQC.GetSPDGCode();// PDG of The Lightest Hadron for the Partner
                G4double pM=0.;             // Mass of the LightestHadron
                if(pPDG==10)
                {
                  G4QChipolino QCh(pQC);    // define the TotalString as a Chipolino
                  pM=QCh.GetQPDG1().GetMass()+QCh.GetQPDG2().GetMass();// Mass of Chipolino
                }
                else pM=G4QPDGCode(pPDG).GetMass();// Mass of theLightestHadron for Partner
                G4double rM=std::sqrt(pM2); // Real mass of the string-partner
                G4double delta=(1.-x)*rM-pM;// @@ Minimum CM disterbance measure
                if(delta > 0. && delta > maD)
                {
                  maD=delta;
#ifdef debug
                  G4cout<<"G4QIonIonCollision::Breed: Subtr,S#"<<restr<<",d="<<maD<<G4endl;
#endif
                  fustr=restr;
                  selX=x;
                  s4M=p4M;
                }
              }
              else if(x <= 0.)               // We are adding to p4M, so use RelVelocity
              {
                G4ThreeVector pV=pP/pE;      // curRelativeVelocity
                G4double dV=(curV-pV).mag2();// SquaredDifferenceBetweenRelVelocities
                if(dV < Vmin)
                {
#ifdef debug
                  G4cout<<"G4QIonIonCollision::Breed: FreeAdd,S#"<<restr<<",x="<<x<<G4endl;
#endif
                  Vmin=dV;
                  fustr=restr;
                  selX=x;
                  s4M=p4M;
                }
              }
#ifdef debug
              G4cout<<"G4QIonIonCollision::Breed:EndOfLOOP r="<<restr<<"<"<<nOfStr<<G4endl;
#endif
            } // End of the LOOP over string-partners for Correction
#ifdef debug
              G4cout<<"G4QIonIonCollision::Breeder: AfterLOOP fustr="<<fustr<<G4endl;
#endif
            if(fustr)
            {
#ifdef edebug
              G4LorentzVector sum4M=s4M+curString4M;
              G4cout<<"G4QIonIonCollision::Breeder: Found Sum4M="<<sum4M<<G4endl;
#endif
              G4QString* pString=strings[fustr];
              curString4M+=selX*s4M;
              if(std::abs(miPDG)%10 > 2)                  // Decay String-Hadron-Resonance
              {
                G4Quasmon Quasm;
                G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
                G4QHadronVector* tmpQHadVec=Quasm.DecayQHadron(sHad); // It deletes sHad
#ifdef debug
                G4cout<<"G4QIonIonCollision::Breed:DecStH,nH="<<tmpQHadVec->size()<<G4endl;
#endif
                for(unsigned aH=0; aH < tmpQHadVec->size(); aH++)
                {
                  theResult->push_back((*tmpQHadVec)[aH]);//TheDecayProdOfHadron is filled
#ifdef debug
                  G4QHadron*   prodH =(*tmpQHadVec)[aH];
                  G4LorentzVector p4M=prodH->Get4Momentum();
                  G4int           PDG=prodH->GetPDGCode();
                  G4int           Chg=prodH->GetCharge();
                  G4int           BaN=prodH->GetBaryonNumber();
                  G4cout<<"-EMC->>G4QIonIonCollision::Breed:St=Had,pH#"<<aH<<" filled, 4M="
                        <<p4M<<", PDG="<<PDG<<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
                }
                tmpQHadVec->clear();
                delete tmpQHadVec;  // Who calls DecayQHadron is responsibleRorClear&Delete
              }
              else if(miPDG == 10)                   // ==> Decay Hadron-Chipolino
              {
                G4QChipolino QCh(miQC);              // define theTotalResid as aChipolino
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
                      G4cerr<<"***G4QIonIonCollision::Breeder: tM="<<ttM<<"->h1="<<h1QPDG
                            <<"("<<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")="<<h1M+h2M<<G4endl;
                      G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"CD");
                    }
                  }
                  G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
                  theResult->push_back(h1H);        // (delete equivalent)  
#ifdef debug
                  G4LorentzVector f4M=h1H->Get4Momentum();
                  G4int           fPD=h1H->GetPDGCode();
                  G4int           fCg=h1H->GetCharge();
                  G4int           fBN=h1H->GetBaryonNumber();
                  G4cout<<"-EMC->>>G4QIonIonCollision::Breed:Str=Hadr Prod-F's filled,f4M="
                        <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
                  G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
                  theResult->push_back(h2H);        // (delete equivalent)  
#ifdef debug
                  G4LorentzVector s4M=h2H->Get4Momentum();
                  G4int           sPD=h2H->GetPDGCode();
                  G4int           sCg=h2H->GetCharge();
                  G4int           sBN=h2H->GetBaryonNumber();
                  G4cout<<"-EMC->>>G4QIonIonCollision::Breed:Str=Hadr Prod-S's filled,s4M="
                        <<s4M<<", sPDG="<<sPD<<", sCg="<<sCg<<", sBN="<<sBN<<G4endl;
#endif
#ifdef edebug
                  G4cout<<"-EMC-Chipo.G4QIonIonCollision::Breed:DecCHECK,c4M="<<curString4M
                        <<", ChQC="<<miQC<<", d4M="<<curString4M-h14M-h24M<<G4endl;
#endif
                }
                else
                {
                  G4cerr<<"***G4QInel::Breeder: tM="<<ttM<<miQC<<"->h1="<<h1QPDG<<"(" <<h1M
                        <<")+h2="<<h1QPDG<<"("<<h2M<<") = "<<h1M+h2M<<G4endl;
                  G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"ChiDec");
                }
              }
              else
              {
                G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
                theResult->push_back(sHad);         // The original string-hadron is filled
#ifdef debug
                G4cout<<"-EMC->>>G4QIonIonCollision::Breeder:Str=Hadr Filled, 4M="
                      <<curString4M<<", PDG="<<miPDG<<G4endl;
#endif
              }
              G4double corF=1-selX;
              G4QParton* Left=pString->GetLeftParton();
              G4QParton* Right=pString->GetRightParton();
              Left->Set4Momentum(corF*Left->Get4Momentum());
              Right->Set4Momentum(corF*Right->Get4Momentum());
#ifdef edebug
              G4cout<<"-EMC-...Cor...G4QIonIonCollision::Breeder:CorCHECK Sum="<<sum4M
                    <<" =? "<<curString4M+pString->Get4Momentum()<<", M="<<miM<<" ?= "
                    <<curString4M.m()<<G4endl;
#endif
#ifdef debug
              G4cout<<">>>G4QIonIonCollision::Breeder:*Corrected* String->Hadr="<<miPDG
                    <<curString4M<<" by String #"<<fustr<<G4endl;
#endif
              continue;                            // Continue the LOOP over the curStrings
            } // End of Found combination for the String on string correction
          } // End of the Try-to-recover String+String correction algorithm
        } // End of IF(CM2>0.)
      } // End of IF(Can try to correct by String-String)
#ifdef debug
      else G4cerr<<"***G4QIonIonCollision::Breed:**No SSCorrection**, next="<<next<<G4endl;
#endif
      // ------------ At this point we can reduce the 3/-3 meson to 1/-1 meson ------------
      G4QParton* lpcS=curString->GetLeftParton();
      G4QParton* rpcS=curString->GetRightParton();
      G4int lPDGcS=lpcS->GetPDGCode();
      G4int rPDGcS=rpcS->GetPDGCode();
      if     (lPDGcS==3 && rPDGcS==-3)
      {
        lpcS->SetPDGCode( 1);
        rpcS->SetPDGCode(-1);
      }
      else if(lPDGcS==-3 && rPDGcS==3)
      {
        lpcS->SetPDGCode(-1);
        rpcS->SetPDGCode( 1);
      }
      // -------- Now the only hope is Correction, using the already prodused Hadrons -----
      G4int nofRH=theResult->size();            // #of resulting Hadrons
#ifdef debug
      G4cout<<"G4QIonIonCollision::Breeder: theH="<<theHadrons<<", #OfH="<<nofRH<<G4endl;
#endif
      if(!theHadrons && nofRH)                  // Hadrons are existing for SH Correction
      {
#ifdef debug
        G4cerr<<"!G4QIonIonCollision::Breeder:CanTrySHCor, nH="<<theResult->size()<<G4endl;
#endif
        // @@ string can be not convertable to one hadron (2,0.0,0,2,0) - To be improved
        G4QContent miQC=curString->GetQC();     // QContent of the Lightest Hadron
        G4int miPDG=miQC.GetSPDGCode();         // PDG of the Lightest Hadron
        G4double miM=0.;                        // Prototype ofMass of theCurLightestHadron
        if(miPDG==10)                           // Mass of the Cur Lightest ChipolinoHadron
        {
          G4QChipolino QCh(miQC);               // define the TotalString as a Chipolino
          miM=QCh.GetQPDG1().GetMass()+QCh.GetQPDG2().GetMass();//MinMass of theStringChipo
        }
        else miM=G4QPDGCode(miPDG).GetMass();   // Mass of the Cur Lightest Hadron
        G4double spM=0.;                        // Mass of the selected Hadron-Partner
        G4ThreeVector cP=curString4M.vect();    // Momentum of the curString
        G4double      cE=curString4M.e();       // Energy of the curString
        G4ThreeVector curV=cP/cE;               // curRelativeVelocity
        G4int reha=0;                           // Hadron # to use beyon the LOOP
        G4int fuha=0;                           // Selected Hadron-partner (0 = NotFound)
        G4double dMmin=DBL_MAX;                 // min Excess of the mass
        G4LorentzVector s4M(0.,0.,0.,0.);       // Selected 4-mom of the Hadron+String
        G4double sM=0.;                         // Selected Mass of the Hadron+String
        for (reha=next; reha < nofRH; reha++)   // LOOP over already collected hadrons
        {
          G4QHadron* pHadron=(*theResult)[reha];// Pointer to the current Partner-Hadron
          G4LorentzVector p4M=pHadron->Get4Momentum();
          G4double         pM=p4M.m();          // Mass of the Partner-Hadron
          G4LorentzVector t4M=p4M+curString4M;  // Total momentum of the compound
          G4double        tM2=t4M.m2();         // Squared total mass of the compound
          if(tM2 >= sqr(pM+miM+eps))            // Condition of possible correction
          {
            G4double tM=std::sqrt(tM2);         // Mass of the Hadron+String compound
            G4double dM=tM-pM-miM;              // Excess of the compound mass
            if(dM < dMmin)
            {
              dMmin=dM;
              fuha=reha;
              spM=pM;
              s4M=t4M;
              sM=tM;
            }
          }
        } // End of the LOOP over string-partners for Correction
        if(fuha)                                // The hadron-partner was found
        { 
          G4QHadron* pHadron=(*theResult)[fuha];// Necessary for update
          G4LorentzVector mi4M(0.,0.,0.,miM);   // Prototype of the new String=Hadron
          if(miM+spM<sM+eps)                    // Decay into CorrectedString+theSameHadron
          {
            G4LorentzVector sp4M(0.,0.,0.,spM);
            if(std::fabs(sM-miM-spM)<=eps)
            {
              G4double part1=miM/(miM+spM);
              mi4M=part1*s4M;
              sp4M=s4M-mi4M;
            }
            else
            {
              if(!G4QHadron(s4M).DecayIn2(mi4M,sp4M))
              {
                G4cerr<<"***G4QIonIonCollision::Breeder: *SH*, tM="<<sM<<"->h1=("<<miPDG
                      <<")"<<miM<<" + h2="<<spM<<" = "<<miM+spM<<G4endl;
                G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"SHChiDec");
              }
            }
            pHadron->Set4Momentum(sp4M);
#ifdef debug
            G4cout<<"-EMC->...>G4QIonIonCollision::Breed: H# "<<fuha<<" is updated, new4M="
                  <<sp4M<<G4endl;
#endif
          }
          else
          {
            G4cerr<<"***G4QInel::Breeder: HS Failed, tM="<<sM<<"->h1M=("<<miPDG<<")"<<miM
                  <<"+h2M="<<spM<<" = "<<miM+spM<<G4endl;
            G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"HSChipoliDec");
          }
          if(std::abs(miPDG)%10 > 2)                  // Decay Hadron-Resonans
          {
            G4Quasmon Quasm;
            G4QHadron* sHad = new G4QHadron(miPDG,mi4M);
            G4QHadronVector* tmpQHadVec=Quasm.DecayQHadron(sHad); // It deletes sHad
#ifdef debug
            G4cout<<"G4QInelast::Breeder: *HS* DecStrHad, nH="<<tmpQHadVec->size()<<G4endl;
#endif
            for(unsigned aH=0; aH < tmpQHadVec->size(); aH++)
            {
              theResult->push_back((*tmpQHadVec)[aH]);// TheDecayProductOfTheHadronIsFilled
#ifdef debug
              G4QHadron*   prodH =(*tmpQHadVec)[aH];
              G4LorentzVector p4M=prodH->Get4Momentum();
              G4int           PDG=prodH->GetPDGCode();
              G4int           Chg=prodH->GetCharge();
              G4int           BaN=prodH->GetBaryonNumber();
              G4cout<<"-EMC->>>G4QIonIonCollision::Breed:Str+Hadr PrH#"<<aH<<" filled, 4M="
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
            if(h1M+h2M<miM+eps)                  // Two particles decay of Chipolino
            {
              G4LorentzVector h14M(0.,0.,0.,h1M);
              G4LorentzVector h24M(0.,0.,0.,h2M);
              if(std::fabs(ttM-h1M-h2M)<=eps)
              {
                G4double part1=h1M/(h1M+h2M);
                h14M=part1*mi4M;
                h24M=mi4M-h14M;
              }
              else
              {
                if(!G4QHadron(mi4M).DecayIn2(h14M,h24M))
                {
                  G4cerr<<"***G4QIonIonCollision::Breeder: HS tM="<<ttM<<"->h1="<<h1QPDG
                        <<"("<<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")="<<h1M+h2M<<G4endl;
                  G4Exception("G4QIonIonCollision::Breeder:","72",FatalException,"ChiDec");
                }
              }
              G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
              theResult->push_back(h1H);         // (delete equivalent)  
#ifdef debug
              G4LorentzVector f4M=h1H->Get4Momentum();
              G4int           fPD=h1H->GetPDGCode();
              G4int           fCg=h1H->GetCharge();
              G4int           fBN=h1H->GetBaryonNumber();
              G4cout<<"-EMC->>G4QIonIonCollision::Breed: CorStrHadr Prod-1 is filled, f4M="
                    <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
              G4LorentzVector n4M=h2H->Get4Momentum();
              G4int           nPD=h2H->GetPDGCode();
              G4int           nCg=h2H->GetCharge();
              G4int           nBN=h2H->GetBaryonNumber();
              G4cout<<"-EMC->>>G4QIonIonCollision::Breed:CorStrHadr Prod-2 is filled, n4M="
                    <<n4M<<", nPDG="<<nPD<<", nCg="<<nCg<<", nBN="<<nBN<<G4endl;
#endif
#ifdef edebug
              G4cout<<"-EMC-...HS-Chipo...G4QIonIonCollision::Breeder:DecCHECK, Ch4M="
                    <<mi4M<<", ChQC="<<miQC<<", d4M="<<mi4M-h14M-h24M<<G4endl;
#endif
            }
          }
          else
          {
            G4QHadron* sHad = new G4QHadron(miPDG, mi4M);
            theResult->push_back(sHad);          // The original string=hadron is filled
#ifdef debug
            G4cout<<">>>>>>G4QIonIonCollision::Breeder: CorStr=Hadr is Filled, 4M="
                  <<curString4M<<", StPDG="<<miPDG<<G4endl;
#endif
          }
#ifdef edebug
          G4cout<<"-EMC-...Cor...G4QIonIonCollision::Breeder:StHadCor CHECK Sum="<<s4M
                <<" =? "<<mi4M+pHadron->Get4Momentum()<<G4endl;
#endif
#ifdef debug
          G4cout<<">>>G4QIonIonCollision::Breeder:*Corrected* String+Hadr="<<miPDG
                <<mi4M<<" by Hadron #"<<reha<<G4endl;
#endif
          continue;                    // Continue the LOOP over the curStrings
        }
        else
        {
#ifdef debug
          G4cout<<"-EMC->>>G4QIonIonCollision::Breeder: Str+Hadr Failed, 4M="<<curString4M
                <<", PDG="<<miPDG<<G4endl;
#endif
        }
        // @@@ convert string to Quasmon with curString4M
        G4QContent curStringQC=curString->GetQC();
        G4Quasmon* stringQuasmon = new G4Quasmon(curStringQC, curString4M);
        theQuasmons.push_back(stringQuasmon);
        continue;                      // Continue the LOOP over the curStrings
      } // End of IF(Can try the String-Hadron correction
    } // End of IF(NO_Hadrons) = Problem solution namespace
    G4Quasmon tmpQ;                                 // @@ an issue of Q to decay resonances
    G4int nHfin=0;
    if(theHadrons) nHfin=theHadrons->size();
    else // !! Sum Up all strings and convert them in a Quasmon (Exception for development)
    {
      G4LorentzVector ss4M(0.,0.,0.,0.);
      G4QContent      ssQC(0,0,0,0,0,0);
      G4int tnSt=strings.size();
      for(G4int i=astring; i < tnSt; ++i)
      {
        G4LorentzVector pS4M=strings[i]->Get4Momentum(); // String 4-momentum
        ss4M+=pS4M;
        G4QContent pSQC=strings[i]->GetQC();             // String Quark Content
        ssQC+=pSQC;
#ifdef debug
        G4cout<<"====>G4QIonIonCollision::Breed:S#"<<i<<",4M="<<pS4M<<",QC="<<pSQC<<G4endl;
#endif
      }
#ifdef debug
      G4cout<<"==>G4QIonIonCollision::Breed:AllStrings are summed up in a Quasmon"<<G4endl;
#endif
      G4Quasmon* stringQuasmon = new G4Quasmon(ssQC, ss4M);
      theQuasmons.push_back(stringQuasmon);
      break;                                   // break the LOOP over Strings
    }
#ifdef debug
    G4cout<<"G4QIonIonCollision::Breeder: Trying to decay hadrons #ofHRes="<<nHfin<<G4endl;
#endif
    for(G4int aTrack=0; aTrack<nHfin; aTrack++)
    {
      G4QHadron* curHadron=(*theHadrons)[aTrack];
      G4int hPDG=curHadron->GetPDGCode();
#ifdef edebug
      G4LorentzVector curH4M=curHadron->Get4Momentum();
      G4int           curHCh=curHadron->GetCharge();
      G4int           curHBN=curHadron->GetBaryonNumber();
#endif
#ifdef debug
      G4cout<<">>>>>>>>G4QIonIonCollision::Breeder:S#"<<astring<<",H#"<<aTrack<<",PDG="
            <<hPDG<<",4M="<<curHadron->Get4Momentum()<<G4endl;
#endif
      if(std::abs(hPDG)%10 > 2)
      {
        G4QHadronVector* tmpQHadVec=tmpQ.DecayQHadron(curHadron); // It deletes curHadron
#ifdef debug
        G4cout<<"G4QIonIonCollision::Breed:-DECAY'S DONE-,nH="<<tmpQHadVec->size()<<G4endl;
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
          G4cout<<"-EMC->>>>G4QIonIonCollision::Breed:Str*Filled, 4M="<<p4M<<", PDG="<<PDG
                <<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
        }
#ifdef edebug
        G4cout<<"-EMC-.G4QIn::Br:Dec,r4M="<<curH4M<<",rC="<<curHCh<<",rB="<<curHBN<<G4endl;
#endif
        tmpQHadVec->clear();
        delete tmpQHadVec;  // Who calls DecayQHadron is responsible for clear & delete
      }
      else                                      // Chipolino is not checked here
      {
        theResult->push_back(curHadron);        // The original hadron is filled
#ifdef edebug
        curString4M-=curH4M;
        G4int curCh=curHadron->GetCharge();
        G4int curBN=curHadron->GetBaryonNumber();
        curStrChg-=curCh;
        curStrBaN-=curBN;
        G4cout<<"-EMC->>>>>>G4QIonIonCollision::Breeder: curH filled 4M="<<curH4M<<",PDG="
              <<curHadron->GetPDGCode()<<", Chg="<<curCh<<", BaN="<<curBN<<G4endl;
#endif
      }
    }
    // clean up (the issues are filled to theResult)
    if(theHadrons) delete theHadrons;
#ifdef edebug
    G4cout<<"-EMC-.......G4QIonIonCollision::Breeder: StringDecay CHECK, r4M="<<curString4M
          <<", rChg="<<curStrChg<<", rBaN="<<curStrBaN<<G4endl;
#endif
  } // End of the main LOOP over decaying strings
  G4LorentzVector rp4M=theProjNucleus.Get4Momentum(); // projNucleus 4-momentum in LS
  G4int rpPDG=theProjNucleus.GetPDG();
  G4QHadron* resPNuc = new G4QHadron(rpPDG,rp4M);
  theResult->push_back(resPNuc);                      // Fill the residual nucleus
  G4LorentzVector rt4M=theTargNucleus.Get4Momentum(); // Nucleus 4-momentum in LS
  G4int rtPDG=theTargNucleus.GetPDG();
  G4QHadron* resTNuc = new G4QHadron(rtPDG,rt4M);
  theResult->push_back(resTNuc);                      // Fill the residual nucleus
#ifdef edebug
  G4LorentzVector s4M(0.,0.,0.,0.);                   // Sum of the Result in LS
  G4int rCh=totChg;
  G4int rBN=totBaN;
  G4int nHadr=theResult->size();
  G4int nQuasm=theQuasmons.size();
  G4cout<<"-EMCLS-G4QInel::Breeder: #ofHadr="<<nHadr<<", #OfQ="<<nQuasm<<", rPA="<<rp4M.m()
        <<"="<<G4QNucleus(rpPDG).GetGSMass()<<", rTA="<<rt4M.m()<<"="
        <<G4QNucleus(rtPDG).GetGSMass()<<G4endl;
  for(G4int i=0; i<nHadr; i++)
  {
    G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
    s4M+=hI4M;
    G4int hChg=(*theResult)[i]->GetCharge();
    rCh-=hChg;
    G4int hBaN=(*theResult)[i]->GetBaryonNumber();
    rBN-=hBaN;
    G4cout<<"-EMCLS-G4QIonIonCollision::Breeder: Hadron#"<<i<<", 4M="<<hI4M<<", PDG="
          <<(*theResult)[i]->GetPDGCode()<<", C="<<hChg<<", B="<<hBaN<<G4endl;
  }
  for(G4int i=0; i<nQuasm; i++)
  {
    G4LorentzVector hI4M=theQuasmons[i]->Get4Momentum();
    s4M+=hI4M;
    G4int hChg=theQuasmons[i]->GetCharge();
    rCh-=hChg;
    G4int hBaN=theQuasmons[i]->GetBaryonNumber();
    rBN-=hBaN;
    G4cout<<"-EMCLS-G4QIonIonCollision::Breeder: Quasmon#"<<i<<", 4M="<<hI4M<<", C="<<hChg
          <<", B="<<hBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QInel::Breed: LS r4M="<<s4M-totLS4M<<",rC="<<rCh<<",rB="<<rBN<<G4endl;
#endif
  // Now we need to coolect particles for creation of a Quasmon @@ improve !!
  G4int nRes=theResult->size();
  if(nRes>2)
  {
    G4QHadronVector::iterator ih;
    G4QHadronVector::iterator nih;
    G4QHadronVector::iterator mih;
    G4QHadronVector::iterator lst=theResult->end();
    lst--;
    G4double minMesEn=DBL_MAX;
    G4double minBarEn=DBL_MAX;
    G4bool nfound=false;
    G4bool mfound=false;
    for(ih = theResult->begin(); ih < theResult->end(); ++ih) if(ih != lst)
    {
      G4LorentzVector h4M=(*ih)->Get4Momentum();
      G4int          hPDG=(*ih)->GetPDGCode();
#ifdef debug
      G4cout<<"%->G4QIonIonCollision::Breeder: TRY hPDG="<<hPDG<<", h4M="<<h4M<<G4endl;
#endif
      if(hPDG>1111 && hPDG<3333)
      {
        G4double bE=h4M.e()-(*ih)->GetMass();
        if(bE < minBarEn)
        {
          minBarEn=bE;
          nih=ih;
          nfound=true;
        }
      }
      else if(hPDG>-1111)
      {
        G4double mE=h4M.e();
        if(mE < minMesEn)
        {
          minMesEn=mE;
          mih=ih;
          mfound=true;
        }
      }
    }
    if(nfound && mfound && minMesEn+minBarEn < maxEn)
    {
      G4QHadron* resNuc = (*theResult)[nRes-1];              // ResidualNucleus is theLastH
      theResult->pop_back();                                 // Fill the QHad-nucleus later
      G4LorentzVector q4M=(*nih)->Get4Momentum()+(*mih)->Get4Momentum();
      G4QContent qQC=(*nih)->GetQC()+(*mih)->GetQC();
      maxEn -= minBarEn+minMesEn;                            // Reduce the energy limit
#ifdef debug
      G4cout<<"%->G4QIonIonCollision::Breeder:Exclude,4M="<<(*nih)->Get4Momentum()
            <<", & 4m="<<(*mih)->Get4Momentum()<<G4endl;
#endif
      delete (*nih);
      delete (*mih);
      if(nih > mih)                                          // Exclude nucleon first
      {
        theResult->erase(nih);                               // erase Baryon from theResult
        theResult->erase(mih);                               // erase Meson from theResult
      }
      else                                                   // Exclude meson first
      {
        theResult->erase(mih);                               // erase Baryon from theResult
        theResult->erase(nih);                               // erase Meson from theResult
      }
#ifdef debug
      G4cout<<"%->G4QI::Breed: BeforeLOOP, dE="<<maxEn<<", nR="<<theResult->size()<<G4endl;
#endif
      if(maxEn > 0.)                                         // More hadrons to absorbe
      {
        for(ih = theResult->begin(); ih < theResult->end(); ih++)
        {
          G4LorentzVector h4M=(*ih)->Get4Momentum();
          G4int          hPDG=(*ih)->GetPDGCode();
          G4double dE=0.;
          if     (hPDG> 1111 && hPDG<3333) dE=h4M.e()-(*ih)->GetMass(); // Baryons
          else if(hPDG>-1111) dE=h4M.e();                    // Mesons
          else continue;                                     // Don't consider anti-baryons
          if(dE < maxEn)
          {
            maxEn-=dE;
            q4M+=h4M;
            qQC+=(*ih)->GetQC();
#ifdef debug
            G4cout<<"%->G4QIonIonCollision::Breed:Exclude,4M="<<h4M<<",dE="<<maxEn<<G4endl;
#endif
            delete (*ih);
            theResult->erase(ih);                            // erase Hadron from theResult
            --ih;
          }
        }
      }
      G4Quasmon* softQuasmon = new G4Quasmon(qQC, q4M);      // SoftPart -> Quasmon
#ifdef debug
      G4cout<<"%->G4QIonIonCollision::Breed:QuasmonIsFilled,4M="<<q4M<<",QC="<<qQC<<G4endl;
#endif
      theQuasmons.push_back(softQuasmon);
      theResult->push_back(resNuc);
    }
#ifdef edebug
    G4LorentzVector f4M(0.,0.,0.,0.);                        // Sum of the Result in LS
    G4int fCh=totChg;
    G4int fBN=totBaN;
    G4int nHd=theResult->size();
    G4int nQm=theQuasmons.size();
    G4cout<<"-EMCLS-G4QIonIonCollision::Breeder:#ofHadr="<<nHd<<", #OfQ="<<nQm
          <<",rPA="<<rt4M.m()<<"="<<G4QNucleus(rtPDG).GetGSMass()
          <<",rTA="<<rt4M.m()<<"="<<G4QNucleus(rtPDG).GetGSMass()<<G4endl;
    for(G4int i=0; i<nHd; i++)
    {
      G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
      f4M+=hI4M;
      G4int hChg=(*theResult)[i]->GetCharge();
      fCh-=hChg;
      G4int hBaN=(*theResult)[i]->GetBaryonNumber();
      fBN-=hBaN;
      G4cout<<"-EMCLS-G4QIonIonCollision::Breeder: Hadron#"<<i<<", 4M="<<hI4M<<", PDG="
            <<(*theResult)[i]->GetPDGCode()<<", C="<<hChg<<", B="<<hBaN<<G4endl;
    }
    for(G4int i=0; i<nQm; i++)
    {
      G4LorentzVector hI4M=theQuasmons[i]->Get4Momentum();
      f4M+=hI4M;
      G4int hChg=theQuasmons[i]->GetCharge();
      fCh-=hChg;
      G4int hBaN=theQuasmons[i]->GetBaryonNumber();
      fBN-=hBaN;
      G4cout<<"-EMCLS-G4QIonIonCollision::Breeder: Quasmon#"<<i<<", 4M="<<hI4M<<", C="
            <<hChg<<", B="<<hBaN<<G4endl;
    }
    G4cout<<"-EMCLS-G4QInel::Breed:, r4M="<<f4M-totLS4M<<", rC="<<fCh<<",rB="<<fBN<<G4endl;
#endif
  } // End of the soft Quasmon Creation
  return;
} // End of Breeder

// Excite double diffractive string
G4bool G4QIonIonCollision::ExciteDiffParticipants(G4QHadron* projectile,
                                                G4QHadron* target) const
{
  G4LorentzVector Pprojectile=projectile->Get4Momentum();
  G4double Mprojectile=projectile->GetMass();
  G4double Mprojectile2=Mprojectile*Mprojectile;
  G4LorentzVector Ptarget=target->Get4Momentum();
  G4double Mtarget=target->GetMass();
  G4double Mtarget2=Mtarget*Mtarget;
#ifdef debug
  G4cout<<"G4QInel::ExciteDiffPartici: Ep="<<Pprojectile.e()<<", Et="<<Ptarget.e()<<G4endl;
#endif
  // Transform momenta to cms and then rotate parallel to z axis;
  G4LorentzVector Psum=Pprojectile+Ptarget;
  G4LorentzRotation toCms(-Psum.boostVector()); // Boost Rotation to CMS
  G4LorentzVector Ptmp=toCms*Pprojectile;
  if(Ptmp.pz()<=0.) // "String" moving backwards in CMS, abort collision !! ? M.K.
  {
#ifdef debug
    G4cout<<"G4QIonIonCollision::ExciteDiffParticipants:*1* abort Collision!! *1*"<<G4endl;
#endif
    return false; 
  }         
  toCms.rotateZ(-Ptmp.phi());
  toCms.rotateY(-Ptmp.theta());
#ifdef debug
  G4cout<<"G4QIonIonCollision::ExciteDiffParticipantts:BeforBoost Pproj="<<Pprojectile
        <<",Ptarg="<<Ptarget<<G4endl;
#endif
  G4LorentzRotation toLab(toCms.inverse()); // Boost Rotation to LabSys (LS)
  Pprojectile.transform(toCms);
  Ptarget.transform(toCms);
#ifdef debug
  G4cout<< "G4QInelast::ExciteDiffParticipantts: AfterBoost Pproj="<<Pprojectile<<",Ptarg="
        <<Ptarget<<", cms4M="<<Pprojectile+Ptarget<<G4endl;
  G4cout<<"G4QIonIonCollis::ExciteDiffParticipants:ProjX+="<<Pprojectile.plus()<<",ProjX-="
        <<Pprojectile.minus()<<", tX+="<< Ptarget.plus()<<",tX-="<<Ptarget.minus()<<G4endl;
#endif
  G4LorentzVector Qmomentum(0.,0.,0.,0.);
  G4int whilecount=0;
#ifdef debug
  G4cout<<"G4QIonIonCollision::ExciteDiffParticipants: Before DO"<<G4endl;
#endif
  do
  {
    //  Generate pt  
    G4double maxPtSquare=sqr(Ptarget.pz());
#ifdef debug
    G4cout<<"G4QIonIonCollision::ExciteDiffParticipants: maxPtSq="<<maxPtSquare<<G4endl;
    if(whilecount++>=500 && whilecount%100==0) // @@ M.K. Hardwired limits 
    G4cout<<"G4QIonIonCollision::ExciteDiffParticipantts: can loop, loopCount="<<whilecount
          <<", maxPtSquare="<<maxPtSquare<<G4endl;
#endif
    if(whilecount>1000)                        // @@ M.K. Hardwired limits 
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::ExciteDiffParticipants: *2* abort Loop!! *2*"<<G4endl;
#endif
      return false;    //  Ignore this interaction 
    }
    Qmomentum=G4LorentzVector(GaussianPt(widthOfPtSquare,maxPtSquare),0);
#ifdef debug
    G4cout<<"G4QIonIonCollision::ExciteDiffParticipants: generated Pt="<<Qmomentum
          <<", ProjPt="<<Pprojectile+Qmomentum<<", TargPt="<<Ptarget-Qmomentum<<G4endl;
#endif
    //  Momentum transfer
    G4double Xmin=0.;
    G4double Xmax=1.;
    G4double Xplus =ChooseX(Xmin,Xmax);
    G4double Xminus=ChooseX(Xmin,Xmax);
#ifdef debug
    G4cout<<"G4QIonIonCollision::ExciteDiffParticipant:X+="<<Xplus<<",X-="<<Xminus<<G4endl;
#endif
    G4double pt2=Qmomentum.vect().mag2();
    G4double Qplus =-pt2/Xminus/Ptarget.minus();
    G4double Qminus= pt2/Xplus /Pprojectile.plus();
    Qmomentum.setPz((Qplus-Qminus)/2);
    Qmomentum.setE( (Qplus+Qminus)/2);
#ifdef debug
    G4cout<<"G4QInelast::ExciteDiffParticip: Qplus="<<Qplus<<", Qminus="<<Qminus<<", pt2="
          <<pt2<<", Qmomentum="<<Qmomentum<<", ProjM="<<(Pprojectile+Qmomentum).mag()
          <<", TargM="<<(Ptarget-Qmomentum).mag()<<G4endl;
#endif
  } while((Pprojectile+Qmomentum).mag2()<=Mprojectile2 ||
          (Ptarget-Qmomentum).mag2()<=Mtarget2);
  Pprojectile += Qmomentum;
  Ptarget     -= Qmomentum;
#ifdef debug
  G4cout<<"G4QInelast::ExciteDiffParticipan: Proj(Q)="<<Pprojectile<<", Targ(Q)="<<Ptarget
        <<", Proj(back)="<<toLab*Pprojectile<<", Targ(bac)="<< toLab*Ptarget << G4endl;
#endif
  // Transform back and update SplitableHadron Participant.
  Pprojectile.transform(toLab);
  Ptarget.transform(toLab);
#ifdef debug
  G4cout<< "G4QIonIonCollision::ExciteDiffParticipants:TargetMass="<<Ptarget.mag()<<G4endl;
#endif
  target->Set4Momentum(Ptarget);  
#ifdef debug
  G4cout<<"G4QIonIonCollision::ExciteDiffParticipant:ProjMass="<<Pprojectile.mag()<<G4endl;
#endif
  projectile->Set4Momentum(Pprojectile);
  return true;
} // End of ExciteDiffParticipants


// Excite single diffractive string
G4bool G4QIonIonCollision::ExciteSingDiffParticipants(G4QHadron* projectile,
                                                    G4QHadron* target) const
{
  G4LorentzVector Pprojectile=projectile->Get4Momentum();
  G4double Mprojectile=projectile->GetMass();
  G4double Mprojectile2=Mprojectile*Mprojectile;
  G4LorentzVector Ptarget=target->Get4Momentum();
  G4double Mtarget=target->GetMass();
  G4double Mtarget2=Mtarget*Mtarget;
#ifdef debug
  G4cout<<"G4QInel::ExSingDiffPartici: Ep="<<Pprojectile.e()<<", Et="<<Ptarget.e()<<G4endl;
#endif
  G4bool KeepProjectile= G4UniformRand() > 0.5;
  // Reset minMass of the non diffractive particle to its value, (minus for rounding...)
  if(KeepProjectile ) 
  {
#ifdef debug
    G4cout<<"-1/2-G4QIonIonCollision::ExSingDiffParticipants: Projectile is fixed"<<G4endl;
#endif
    Mprojectile2 = projectile->GetMass2()*(1.-perCent); // Isn't it too big reduction? M.K.
  }
  else
  {
#ifdef debug
    G4cout<<"---1/2---G4QIonIonCollision::ExSingDiffParticipants: Target is fixed"<<G4endl;
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
    G4cout<<"G4QIonIonCollision::ExciteSingDiffParticipants: *1*abortCollision*1*"<<G4endl;
#endif
    return false; 
  }         
  toCms.rotateZ(-Ptmp.phi());
  toCms.rotateY(-Ptmp.theta());
#ifdef debug
  G4cout<<"G4QInel::ExciteSingDiffParticipantts: Be4Boost Pproj="<<Pprojectile<<", Ptarg="
        <<Ptarget<<G4endl;
#endif
  G4LorentzRotation toLab(toCms.inverse()); // Boost Rotation to LabSys (LS)
  Pprojectile.transform(toCms);
  Ptarget.transform(toCms);
#ifdef debug
  G4cout<< "G4QInelast::ExciteDiffParticipantts: AfterBoost Pproj="<<Pprojectile<<"Ptarg="
        <<Ptarget<<", cms4M="<<Pprojectile+Ptarget<<G4endl;

  G4cout<<"G4QInelast::ExciteDiffParticipantts: ProjX+="<<Pprojectile.plus()<<", ProjX-="
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
    G4cout<<"G4QIonIonCollision::ExciteSingDiffParticipantts: can loop, loopCount="
          <<whilecount<<", maxPtSquare="<<maxPtSquare<<G4endl;
#endif
    if(whilecount>1000)                        // @@ M.K. Hardwired limits 
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::ExciteSingDiffParticipants:*2* abortLoop!! *2*"<<G4endl;
#endif
      return false;    //  Ignore this interaction 
    }
    Qmomentum=G4LorentzVector(GaussianPt(widthOfPtSquare,maxPtSquare),0);
#ifdef debug
    G4cout<<"G4QInelast::ExciteSingDiffParticipants: generated Pt="<<Qmomentum
          <<", ProjPt="<<Pprojectile+Qmomentum<<", TargPt="<<Ptarget-Qmomentum<<G4endl;
#endif
    //  Momentum transfer
    G4double Xmin=0.;
    G4double Xmax=1.;
    G4double Xplus =ChooseX(Xmin,Xmax);
    G4double Xminus=ChooseX(Xmin,Xmax);
#ifdef debug
    G4cout<<"G4QInel::ExciteSingDiffPartici: X-plus="<<Xplus<<", X-minus="<<Xminus<<G4endl;
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
    G4cout<<"G4QInel::ExciteDiffParticip: Qplus="<<Qplus<<", Qminus="<<Qminus<<", pt2="
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
  G4cout<<"G4QIonIonCollision::ExciteSingDiffParticipan: Proj(Q)="<<Pprojectile<<"(E="
        <<Pprojectile.e()<<"), Targ(Q)="<<Ptarget<<"(E="<<Ptarget.e()
        <<"), Proj(back)="<<toLab*Pprojectile<<", Targ(bac)="<< toLab*Ptarget << G4endl;
#endif
  // Transform back and update SplitableHadron Participant.
  Pprojectile.transform(toLab);
  Ptarget.transform(toLab);
#ifdef debug
  G4cout<< "G4QIonIonCollision::ExciteSingleDiffParticipants: TgM="<<Ptarget.mag()<<G4endl;
#endif
  target->Set4Momentum(Ptarget);  
#ifdef debug
  G4cout<<"G4QIonIonCollision::ExciteSingleParticipants:ProjM="<<Pprojectile.mag()<<G4endl;
#endif
  projectile->Set4Momentum(Pprojectile);
  return true;
} // End of ExciteSingleDiffParticipants

void G4QIonIonCollision::SetParameters(G4int nC,G4double stT, G4double tbD, G4double SigPt)
{//  ======================================================================================
  nCutMax            = nC;             // max number of pomeron cuts
  stringTension      = stT;            // string tension for absorbed energy
  tubeDensity        = tbD;            // Flux Tube Density of nuclear nucleons
  widthOfPtSquare    = -2*SigPt*SigPt; // width^2 of pt for string excitation
}

G4double G4QIonIonCollision::ChooseX(G4double Xmin, G4double Xmax) const
{
  // choose an x between Xmin and Xmax with P(x) ~ 1/x @@ M.K. -> 1/sqrt(x)
  //G4double range=Xmax-Xmin;
  if(Xmax == Xmin) return Xmin;
  if( Xmin < 0. || Xmax < Xmin) 
  {
    G4cerr<<"***G4QIonIonCollision::ChooseX: Xmin="<<Xmin<<", Xmax="<<Xmax<< G4endl;
    G4Exception("G4QIonIonCollision::ChooseX:","72",FatalException,"Bad X or X-Range");
  }
  //G4double x;
  //do {x=Xmin+G4UniformRand()*range;} while ( Xmin/x < G4UniformRand() );
  G4double sxi=std::sqrt(Xmin);
  G4double x=sqr(sxi+G4UniformRand()*(std::sqrt(Xmax)-sxi));
#ifdef debug
  G4cout<<"G4QIonIonCollision::ChooseX: DiffractiveX="<<x<<G4endl;
#endif
  return x;
} // End of ChooseX

// Add CHIPS exponential Pt distribution (see Fragmentation)
G4ThreeVector G4QIonIonCollision::GaussianPt(G4double widthSq, G4double maxPtSquare) const
{
#ifdef debug
  G4cout<<"G4QIonIonCollision::GaussianPt: w^2="<<widthSq<<",maxPt2="<<maxPtSquare<<G4endl;
#endif
  G4double pt2=0.;
  G4double rm=maxPtSquare/widthSq;                      // Negative
  if(rm>-.01) pt2=widthSq*(std::sqrt(1.-G4UniformRand()*(1.-sqr(1.+rm)))-1.);
  else        pt2=widthSq*std::log(1.-G4UniformRand()*(1.-std::exp(rm)));
  pt2=std::sqrt(pt2);                                   // It is not pt2, but pt now
  G4double phi=G4UniformRand()*twopi;
  return G4ThreeVector(pt2*std::cos(phi),pt2*std::sin(phi),0.);    
} // End of GaussianPt

G4int G4QIonIonCollision::SumPartonPDG(G4int PDG1, G4int PDG2) const
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
      G4Exception("G4QIonIonCollision::SumPartonPDG:","72",FatalException,"Q&ADiQnoMatch");
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
      G4Exception("G4QIonIonCollision::SumPartonPDG:","72",FatalException,"ADiQ&QnoMatch");
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
      G4Exception("G4QIonIonCollision::SumPartonPDG:","72",FatalException,"DiQ&AQnoMatch");
    }
  }
  else if (PDG2 > 99 && PDG1 >-7 && PDG1 < 0) // Sum up AQ and DiQ in Q
  {
    G4int PDG=PDG2/100;
    if(PDG1==-PDG/10) return PDG%10;
    if(PDG1==-PDG%10) return PDG/10;
    else
    {
      G4cerr<<"***G4QIonIonCollision::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
      G4Exception("G4QIonIonCollision::SumPartonPDG:","72",FatalException,"AQ&DiQnoMatch");
    }
  }
  else
  {
    G4cerr<<"***G4QIonIonCollision::SumPartonPDG: PDG1="<<PDG1<<", PDG2="<<PDG2<<G4endl;
    G4Exception("G4QIonIonCollision::SumPartonPDG:","72",FatalException,"Can'tSumUpParts");
  }
  return 0;
} // End of SumPartonPDG

// Reduces quark pairs (unsigned 2 digits) to quark singles (unsigned)
std::pair<G4int,G4int> G4QIonIonCollision::ReducePair(G4int P1, G4int P2) const
{
#ifdef debug
  G4cout<<"G4QIonIonCollision::ReducePair: **Called** P1="<<P1<<", P2="<<P2<<G4endl;
#endif
  G4int P11=P1/10;
  G4int P12=P1%10;
  G4int P21=P2/10;
  G4int P22=P2%10;
  if     (P11==P21) return std::make_pair(P12,P22);
  else if(P11==P22) return std::make_pair(P12,P21);
  else if(P12==P21) return std::make_pair(P11,P22);
  else if(P12==P22) return std::make_pair(P11,P21);
  //#ifdef debug
  G4cout<<"-Warning-G4QIonIonCollision::ReducePair:**Failed**,P1="<<P1<<",P2="<<P2<<G4endl;
  //#endif
  return std::make_pair(0,0);                         // Reduction failed
}

// Select LL/RR (1) or LR/RL (-1) annihilation order (0, if the annihilation is impossible)
G4int G4QIonIonCollision::AnnihilationOrder(G4int LS, G4int MPS, G4int uPDG, G4int mPDG,
                                          G4int sPDG, G4int nPDG) // ^L^^ Partner^^R^
{//                                             ^L^^ Curent ^^R^
  G4int Ord=0;
  //       Curent   Partner
  if      (LS==2 && MPS==2 )                 // Fuse 2 Q-aQ strings to DiQ-aDiQ
  {
#ifdef debug
    G4cout<<"G4QIonIonCollision::AnnihilationOrder:QaQ/QaQ->DiQ-aDiQ, uPDG="<<uPDG<<G4endl;
#endif
    if     ( (uPDG>0 && sPDG>0 && mPDG<0 && nPDG<0) || 
             (uPDG<0 && sPDG<0 && mPDG>0 && nPDG>0) ) Ord= 1; // LL/RR
    else if( (uPDG>0 && nPDG>0 && mPDG<0 && sPDG<0) || 
             (uPDG<0 && nPDG<0 && mPDG>0 && sPDG>0) ) Ord=-1; // LR/RL
    else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong 22 fusion, L="
               <<uPDG<<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
  }
  else if ( LS==2 && MPS==3 )
  {
    if     (uPDG > 7)                                // Fuse pLDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihOrder:pLDiQ, sPDG="<<sPDG<<", nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<0 && (-sPDG==uPDG/1000 || -sPDG==(uPDG/100)%10) ) Ord= 1; // LL/RR
      else if( nPDG<0 && (-nPDG==uPDG/1000 || -nPDG==(uPDG/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong pLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (mPDG > 7)                               // Fuse pRDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihOrder:pRDiQ, sPDG="<<sPDG<<", nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<0 && (-sPDG==mPDG/1000 || -sPDG==(mPDG/100)%10) ) Ord=-1; // LR/RL
      else if( nPDG<0 && (-nPDG==mPDG/1000 || -nPDG==(mPDG/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong pRDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (uPDG <-7)                               // Fuse pLaDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihOrder:pLaDiQ, sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG>0 && (sPDG==(-uPDG)/1000 || sPDG==((-uPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( nPDG>0 && (nPDG==(-uPDG)/1000 || nPDG==((-uPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong pLaDiQ, L="<<uPDG
                 <<", R="<<mPDG<<", cL="<<sPDG<<", cR="<<nPDG<<G4endl;
    }
    else if (mPDG <-7)                              // Fuse pRaDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihOrder:pRaDiQ, sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG>0 && (sPDG==(-mPDG)/1000 || sPDG==((-mPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( nPDG>0 && (nPDG==(-mPDG)/1000 || nPDG==((-mPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong pRaDiQ, L="<<uPDG
                 <<", R="<<mPDG<<", cL="<<sPDG<<", cR="<<nPDG<<G4endl;
    }
    else if( (sPDG<0 && (-sPDG==mPDG || -sPDG==uPDG) ) ||
             (nPDG<0 && (-nPDG==mPDG || -nPDG==uPDG) ) ) Ord= 2; // @@ Annihilation fusion
#ifdef debug
    else G4cout<<"-Warning-G4QIonIonCollision::AnnihilatOrder:Wrong23StringFusion"<<G4endl;
    G4cout<<"G4QIonIonCollision::AnnihilationOrder: Ord="<<Ord<<",sPDG="<<sPDG<<",nPDG="
          <<nPDG<<", uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
  }
  else if ( LS==3 && MPS==2 )
  {
    if     (sPDG > 7)                                // Fuse cLDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihOrder:cLDiQ, uPDG="<<uPDG<<", mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<0 && (-uPDG==sPDG/1000 || -uPDG==(sPDG/100)%10) ) Ord= 1; // LL/RR
      else if( mPDG<0 && (-mPDG==sPDG/1000 || -mPDG==(sPDG/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong cLDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (nPDG > 7)                               // Fuse cRDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihOrder:cRDiQ, uPDG="<<uPDG<<", mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<0 && (-uPDG==nPDG/1000 || -uPDG==(nPDG/100)%10) ) Ord=-1; // LR/RL
      else if( mPDG<0 && (-mPDG==nPDG/1000 || -mPDG==(nPDG/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong cRDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (sPDG <-7)                               // Fuse cLaDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihOrder:cLaDiQ, uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG>0 && (uPDG==(-sPDG)/1000 || uPDG==((-sPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( mPDG>0 && (mPDG==(-sPDG)/1000 || mPDG==((-sPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong cLaDiQ, L="<<uPDG
                 <<", R="<<mPDG<<", cL="<<sPDG<<", cR="<<nPDG<<G4endl;
    }
    else if (nPDG <-7)                              // Fuse cRaDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihOrder:cRaDiQ, uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG>0 && (uPDG==(-nPDG)/1000 || uPDG==((-nPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( mPDG>0 && (mPDG==(-nPDG)/1000 || mPDG==((-nPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong cRaDiQ, L="<<uPDG
                 <<", R="<<mPDG<<", cL="<<sPDG<<", cR="<<nPDG<<G4endl;
    }
    else if( (uPDG<0 && (-uPDG==sPDG || -uPDG==nPDG) ) ||
             (mPDG<0 && (-mPDG==sPDG || -mPDG==nPDG) ) ) Ord=2; // @@ Annihilation fusion
#ifdef debug
    else G4cout<<"-Warning-G4QIonIonCollision::AnnihilatOrder:Wrong32StringFusion"<<G4endl;
    G4cout<<"G4QIonIonCollision::AnnihilationOrder: Ord="<<Ord<<",sPDG="<<sPDG<<",nPDG="
          <<nPDG<<", uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
  }
  else if ( (LS==2 && MPS==4) || (LS==4 && MPS==2) )
  {
    if     (uPDG > 7)  // Annihilate pLDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihilOrder:pLDiQ,sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<0 && (-sPDG==uPDG/1000 || -sPDG==(uPDG/100)%10) &&
               (nPDG==(-mPDG)/1000 || nPDG==((-mPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( nPDG<0 && (-nPDG==uPDG/1000 || -nPDG==(uPDG/100)%10) &&
               (sPDG==(-mPDG)/1000 || sPDG==((-mPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong 24 pLDiQ, L="
                 <<uPDG<<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (mPDG >7) // Annihilate pRDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihilOrder:PRDiQ,sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<0 && (-sPDG==mPDG/1000 || -sPDG==(mPDG/100)%10) &&
               (nPDG==(-uPDG)/1000 || nPDG==((-uPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( nPDG<0 && (-nPDG==mPDG/1000 || -nPDG==(mPDG/100)%10) &&
               (sPDG==(-uPDG)/1000 || sPDG==((-uPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong 24 pLDiQ, L="
                 <<uPDG<<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (sPDG > 7)   // Annihilate cLDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihilOrder:cLDiQ,uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<0 && (-uPDG==sPDG/1000 || -uPDG==(sPDG/100)%10) &&
               (mPDG==(-nPDG)/1000 || mPDG==((-nPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( mPDG<0 && (-mPDG==sPDG/1000 || -mPDG==(sPDG/100)%10) &&
               (uPDG==(-nPDG)/1000 || uPDG==((-nPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong 24 cLDiQ, L="
                 <<uPDG<<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if (nPDG > 7)   // Annihilate cRDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihilOrder:cRDiQ,uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<0 && (-uPDG==nPDG/1000 || -uPDG==(nPDG/100)%10) &&
               (mPDG==(-sPDG)/1000 || mPDG==((-sPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( mPDG<0 && (-mPDG==nPDG/1000 || -mPDG==(nPDG/100)%10) &&
               (uPDG==(-sPDG)/1000 || uPDG==((-sPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong 24 cRDiQ, L="
                 <<uPDG<<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
#ifdef debug
    else G4cout<<"-Warning-G4QIonIonCollision::AnnihilOrder: Wrong 24StringFusion"<<G4endl;
    G4cout<<"G4QIonIonCollision::AnnihilationOrder: Ord="<<Ord<<",sPDG="<<sPDG<<",nPDG="
          <<nPDG<<", uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
  }
  else if ( LS==3 && MPS==3 )
  {
    if     (uPDG > 7)  // Annihilate pLDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihilOrder:pLDiQ,sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<-7 && (-nPDG==uPDG/1000 || -nPDG==(uPDG/100)%10) &&
               (mPDG==(-sPDG)/1000 || mPDG==((-sPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( nPDG<-7 && (-sPDG==uPDG/1000 || -sPDG==(uPDG/100)%10) &&
               (mPDG==(-nPDG)/1000 || mPDG==((-nPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong 33 pLDiQ, L="
                 <<uPDG<<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if(mPDG > 7)  // Annihilate pRDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihilOrder:pRDiQ,sPDG="<<sPDG<<",nPDG="<<nPDG<<G4endl;
#endif
      if     ( sPDG<-7 && (-nPDG==mPDG/1000 || -nPDG==(mPDG/100)%10) &&
               (uPDG==(-sPDG)/1000 || uPDG==((-sPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( nPDG<-7 && (-sPDG==mPDG/1000 || -sPDG==(mPDG/100)%10) &&
               (uPDG==(-nPDG)/1000 || uPDG==((-nPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong33pRDiQ, L="<<uPDG
                 <<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if(sPDG > 7)  // Annihilate cLDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihilOrder:cLDiQ,uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<-7 && (-mPDG==sPDG/1000 || -mPDG==(sPDG/100)%10) &&
               (nPDG==(-uPDG)/1000 || nPDG==((-uPDG)/100)%10) ) Ord=-1; // LR/RL
      else if( mPDG<-7 && (-uPDG==sPDG/1000 || -uPDG==(sPDG/100)%10) &&
               (nPDG==(-mPDG)/1000 || nPDG==((-mPDG)/100)%10) ) Ord= 1; // LL/RR
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong 33 cLDiQ, L="
                 <<uPDG<<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
    else if(nPDG > 7)  // Annihilate cRDiQ
    {
#ifdef debug
      G4cout<<"G4QIonIonCollision::AnnihilOrder:cRDiQ,uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
      if     ( uPDG<-7 && (-mPDG==nPDG/1000 || -mPDG==(nPDG/100)%10) &&
               (nPDG==(-uPDG)/1000 || nPDG==((-uPDG)/100)%10) ) Ord= 1; // LL/RR
      else if( mPDG<-7 && (-uPDG==nPDG/1000 || -sPDG==(nPDG/100)%10) &&
               (sPDG==(-mPDG)/1000 || sPDG==((-mPDG)/100)%10) ) Ord=-1; // LR/RL
      else G4cerr<<"-Warning-G4QIonIonCollision::AnnihilationOrder: Wrong 33 cRDiQ, L="
                 <<uPDG<<",R="<<mPDG<<",cL="<<sPDG<<",cR="<<nPDG<<G4endl;
    }
#ifdef debug
    else G4cout<<"-Warning-G4QIonIonCollision::AnnihilOrder: Wrong 33StringFusion"<<G4endl;
    G4cout<<"G4QIonIonCollision::AnnihilationOrder: Ord="<<Ord<<",sPDG="<<sPDG<<",nPDG="
          <<nPDG<<", uPDG="<<uPDG<<",mPDG="<<mPDG<<G4endl;
#endif
  }
  return Ord;
}

void G4QIonIonCollision::SwapPartons() // Swap string partons, if a string has negative M2
{
  static const G4double baryM=800.;                  // Mass excess for baryons
  G4QStringVector::iterator ist;
  for(ist = strings.begin(); ist < strings.end(); ist++)
  {
    G4LorentzVector cS4M=(*ist)->Get4Momentum();
    G4double cSM2=cS4M.m2();                         // Squared mass of the String
    if(cSM2<0.1)                                     // Parton Swapping is needed
    {
      G4QParton* cLeft=(*ist)->GetLeftParton();      // Current String Left Parton 
      G4QParton* cRight=(*ist)->GetRightParton();    // Current String Right Parton 
      G4int cLPDG=cLeft->GetPDGCode();
      G4int cRPDG=cRight->GetPDGCode();
      G4LorentzVector cL4M=cLeft->Get4Momentum();
      G4LorentzVector cR4M=cRight->Get4Momentum();
      G4int cLT=cLeft->GetType();
      G4int cRT=cRight->GetType();
      G4QStringVector::iterator sst;                 // Selected partner string
      G4QStringVector::iterator pst;                 // LOOP iterator
      G4double maxM=-DBL_MAX;                        // Swapping providing the max mass
      G4int    sDir=0;                               // Selected direction of swapping
      for(pst = strings.begin(); pst < strings.end(); pst++) if(pst != ist)
      {
        G4QParton* pLeft=(*pst)->GetLeftParton();    // Partner String Left Parton 
        G4QParton* pRight=(*pst)->GetRightParton();  // Partner String Right Parton 
        G4int pLPDG=pLeft->GetPDGCode();
        G4int pRPDG=pRight->GetPDGCode();
        G4LorentzVector pL4M=pLeft->Get4Momentum();
        G4LorentzVector pR4M=pRight->Get4Momentum();
        G4int pLT=pLeft->GetType();
        G4int pRT=pRight->GetType();
        G4double LM=0.;
        G4double RM=0.;
        if(((cLPDG<-7 || (cLPDG>0 && cLPDG< 7) ) && (pLPDG<-7 || (pLPDG>0 && pLPDG< 7) ))||
           ((cLPDG> 7 || (cLPDG<0 && cLPDG>-7) ) && (pLPDG> 7 || (pLPDG<0 && cLPDG>-7) )))
        {
          G4double pLM2=(cL4M+pR4M).m2();                      // new partner M2
          G4double cLM2=(cR4M+pL4M).m2();                      // new partner M2
          if(pLM2>0. && cLM2>0.)
          {
            G4double pLM=std::sqrt(pLM2);
            if(cLT+pRT==3) pLM-=baryM;
            G4double cLM=std::sqrt(cLM2);
            if(cRT+pLT==3) cLM-=baryM;
            LM=std::min(pLM2,cLM2);
          }
        }
        if(((cRPDG<-7 || (cRPDG>0 && cRPDG< 7) ) && (pRPDG<-7 || (pRPDG>0 && pRPDG< 7) ))||
           ((cRPDG> 7 || (cRPDG<0 && cRPDG>-7) ) && (pRPDG> 7 || (pRPDG<0 && cRPDG>-7) )) )
        {
          G4double pRM2=(cR4M+pL4M).m2();                      // new partner M2
          G4double cRM2=(cL4M+pR4M).m2();                      // new partner M2
          if(pRM2>0. && cRM2>0.)
          {
            G4double pRM=std::sqrt(pRM2);
            if(cRT+pLT==3) pRM-=baryM;
            G4double cRM=std::sqrt(cRM2);
            if(cLT+pRT==3) cRM-=baryM;
            RM=std::min(pRM,cRM);
          }
        }
        G4int dir=0;
        G4double sM=0.;
        if( LM && LM > RM )
        {
          dir= 1;
          sM=LM;
        }
        else if(RM)
        {
          dir=-1;
          sM=RM;
        }
        if(sM > maxM)
        {
          sst=pst;
          maxM=sM;
          sDir=dir;
        }
      }
      if(sDir)
      {
        G4QParton* pLeft=(*sst)->GetLeftParton();    // Partner String Left Parton 
        G4QParton* pRight=(*sst)->GetRightParton();  // Partner String Right Parton 
        G4QParton* swap=pLeft;                       // Prototype remember the partner Left
        if(sDir>0)                                   // swap left partons
        {
          (*sst)->SetLeftParton(cLeft);
          (*ist)->SetLeftParton(swap);
        }
        else
        {
          swap=pRight;
          (*sst)->SetRightParton(cRight);
          (*ist)->SetRightParton(swap);
        }
      }
#ifdef debug
      else G4cout<<"***G4QIonIonCollision::SwapPartons:**Failed**,cLPDG="<<cLPDG<<",cRPDG="
                 <<cRPDG<<",-->cM2="<<cSM2<<G4endl;
#endif
      
    }
  }
}

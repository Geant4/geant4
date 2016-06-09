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
// $Id$
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// Promoting model parameters from local variables class properties @@(? M.K.)

// Definition of static parameters
G4int    G4QFragmentation::nCutMax=27;
//
//G4double G4QFragmentation::stringTension=GeV/fermi; // Minimum value
G4double G4QFragmentation::stringTension=2*GeV/fermi; // ST of intranuclear string
//G4double G4QFragmentation::stringTension=2.7*GeV/fermi; // Maximum value
//
//G4double G4QFragmentation::tubeDensity  =0./fermi;  // Nucleons per fermi (Min Value)
G4double G4QFragmentation::tubeDensity  =.5/fermi;  // Nucleons per fermi
//G4double G4QFragmentation::tubeDensity  =1./fermi;  // Nucleons per fermi (Max Value)
// Parameters of diffractional fragmentation (was .72*)
G4double G4QFragmentation::widthOfPtSquare=-GeV*GeV;// pt -width2 forStringExcitation

G4QFragmentation::G4QFragmentation(const G4QNucleus &aNucleus, const G4QHadron &aPrimary)
{
  static const G4double mProt= G4QPDGCode(2212).GetMass(); // Mass of proton
  static const G4double mProt2= mProt*mProt;               // SquaredMass of proton
  static const G4double mPi0= G4QPDGCode(111).GetMass();   // Mass of Pi0
  static const G4double thresh= 1000;                      // Diffraction threshold (MeV)
  theWorld= G4QCHIPSWorld::Get();                          // Pointer to CHIPS World
  theQuasiElastic=G4QuasiFreeRatios::GetPointer();         // Pointer to CHIPS quasiElastic
  theDiffraction=G4QDiffractionRatio::GetPointer();        // Pointer to CHIPS Diffraction
  theResult = new G4QHadronVector;        // Define theResultingHadronVector
  G4bool stringsInitted=false;            // Strings are initiated
  G4QHadron aProjectile = aPrimary;       // As a primary is const
  G4LorentzRotation toZ;                  // Lorentz Transformation to the projectileSystem
  G4LorentzVector proj4M=aProjectile.Get4Momentum(); // Projectile 4-momentum in LS
#ifdef edebug
  G4LorentzVector targ4M=aProjectile.Get4Momentum(); // Target 4-momentum in LS
  G4double tgMass=aNucleus.GetGSMass();   // Ground state mass of the nucleus
  G4double cM=0.;
  G4double cs=(proj4M+targ4M).mag2();   // s of the compound
  if(cs > 0.) cM=std::sqrt(cs);
  G4cout<<"G4QFragmentation::Construct: *Called*, p4M="<<proj4M<<", A="<<aNucleus<<tgMass
        <<",M="<<cM<<",s="<<cs<<",t4M="<<targ4M<<G4endl;
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
  G4cout<<"-EMC-G4QFragmentation::Construct:tLS4M="<<totLS4M<<",tZLS4M="<<totZLS4M<<G4endl;
  // === From nere all consideration is made in the rotated LS frame (proj is along Z) ===
#endif
  G4LorentzRotation toLab(toZ.inverse()); // Lorentz Transfornation "ZLS"->LS (at the end)
  G4int pPDG=aProjectile.GetPDGCode();    // The PDG code of the projectile
  G4double projM=aProjectile.GetMass();   // Mass of the projectile
  G4QInteractionVector theInteractions;   // A vector of interactions (taken from the Body)
  G4QPartonPairVector  thePartonPairs;    // The parton pairs (taken from the Body)
  G4int                ModelMode=SOFT;    // The current model type (taken from the Body)
  theNucleus=G4QNucleus(tZ,tN);           // Create theNucleus (Body) to Move From LS to CM
  theNucleus.InitByPDG(tPDG);             // Reinit the Nucleus for the new Attempt
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: Nucleus4Mom="<<theNucleus.Get4Momentum()<<G4endl;
#endif
  theNucleus.Init3D();                    // 3D-initialisation(nucleons) of theNucleusClone
  // Now we can make the Quasi-Elastic (@@ Better to select a nucleon from the perifery)
  std::pair<G4double,G4double> ratios=std::make_pair(0.,0.);
  G4int apPDG=std::abs(pPDG);
  if(apPDG>99)                            // Diffraction or Quasi-elastic 
  {
   G4double pMom=proj4M.vect().mag();                  // proj.Momentum in MeV (indepUnits)
   ratios = theQuasiElastic->GetRatios(pMom, pPDG, tZ, tN);
   G4double qeRat = ratios.first*ratios.second;        // quasi-elastic part [qe/in]
   G4double difRat= theDiffraction->GetRatio(pMom, pPDG, tZ, tN); // diffrPart [d/(in-qe)]
   //if(qeRat < 1.) difRat /= (1.-qeRat)*.5;             // Double for Top/Bottom @@ ?
   if(qeRat<1. && proj4M.e()>thresh) difRat /= 1.-qeRat; // Both Top & Bottom
   else difRat=0.;                                     // Close diffraction for low Energy
   difRat += qeRat;                                    // for the diffraction selection
   G4double rnd=G4UniformRand();
   if(rnd < qeRat)                                     // --> Make Quasi-Elastic reaction
   {
    theNucleus.StartLoop();                            // Prepare Loop ovder nucleons
    G4QHadron* pNucleon = theNucleus.GetNextNucleon(); // Get the next nucleon to try
    G4LorentzVector pN4M=pNucleon->Get4Momentum();     // Selected Nucleon 4-momentum
    G4int pNPDG=pNucleon->GetPDGCode();                // Selected Nucleon PDG code
#ifdef debug
    G4cout<<">QE>G4QFragmentation::Construct:TryQE,R="<<ratios.first*ratios.second<<",N4M="
          <<pN4M<<",NPDG="<<pNPDG<<G4endl;
#endif
    std::pair<G4LorentzVector,G4LorentzVector> QEout=theQuasiElastic->Scatter(pNPDG,pN4M,
                                                                              pPDG,proj4M);
    G4bool CoulB = true;                               // No Under Coulomb Barrier flag
    G4double CB=theNucleus.CoulombBarrier(1, 1);       // Coulomb barrier for protons
    if( (pNPDG==2212 && QEout.first.e()-mProt < CB) ||
        (pPDG==2212 && QEout.second.e()-mProt < CB) ) CoulB = false; // at least one UCB
    if(QEout.first.e() > 0 && CoulB)                   // ==> Successful scattering
    {
      G4QHadron* qfN = new G4QHadron(pNPDG,QEout.first);
      theResult->push_back(qfN);                       // Scattered Quasi-free nucleon  
      G4QHadron* qfP = new G4QHadron(pPDG,QEout.second);
      theResult->push_back(qfP);                       // Scattered Projectile  
      theNucleus.SubtractNucleon(pNucleon);            // Exclude the used nucleon from Nuc
      G4LorentzVector r4M=theNucleus.Get4Momentum();   // Nucleus 4-momentum in LS
      G4int rPDG=theNucleus.GetPDG();                  // Nuclear PDG
      G4QHadron* resNuc = new G4QHadron(rPDG,r4M);     // The residual Nucleus
      theResult->push_back(resNuc);                    // Fill the residual nucleus
#ifdef debug
      G4cout<<"-->QE-->G4QFragmentation::Construct:QEDone, N4M="<<QEout.first<<", p4M="
            <<QEout.second<<G4endl;
#endif
      return;                                          // The Quasielastic is the only act
    }
   } // End of quasi-elastic reaction
   else if(rnd < difRat)                               // --> Make diffractive reaction
   {
#ifdef debug
     G4cout<<"-->Dif-->G4QFragmentation::Construct: qe="<<qeRat<<", dif="<<difRat-qeRat
           <<",P="<<proj4M.vect().mag()<<", tZ="<<tZ<<", tN="<<tN<<G4endl;
#endif
     G4QHadronVector* out=0;
     //if(G4UniformRand()>0.5) out=theDiffraction->TargFragment(pPDG, proj4M, tZ, tN);
     //else                    out=theDiffraction->ProjFragment(pPDG, proj4M, tZ, tN);
     //theDiffraction->TargFragment(pPDG, proj4M, tZ, tN); // !! is not debugged !!
     out=theDiffraction->ProjFragment(pPDG, proj4M, tZ, tN);
     G4int nSec=out->size();                           // #of secondaries in diffReaction
     if(nSec>1) for(G4int i=0; i<nSec; i++) theResult->push_back((*out)[i]);
     else if(nSec>0)
     {
#ifdef debug
       G4cout<<"-Warning-G4QFragmentation::Construct: 1 secondary in Diffractionp 4M="
             <<proj4M<<", s4M="<<(*out)[0]->Get4Momentum()<<G4endl;
#endif
       delete (*out)[0];
     }
     out->clear();                                     // Clean up the output QHadronVector
     delete out;                                       // Delete the output QHadronVector
     if(nSec>1) return;
   } // End of diffraction
  }
#ifdef edebug
  G4LorentzVector sum1=theNucleus.GetNucleons4Momentum(); // Sum ofNucleons4M inRotatedLS
  G4cout<<"-EMC-G4QFragmentation::Construct: Nuc4M="<<sum1<<G4endl;
#endif
  G4ThreeVector theCurrentVelocity(0.,0.,0.);        // "CM" velosity (temporary)
  // @@ "target Nucleon" == "Proton at rest" case (M.K. ?)
  G4double nCons = 1;                                // 1 or baryonNum of the Projectile
  G4int projAbsB=std::abs(aProjectile.GetBaryonNumber());// Fragment/Baryon (Meson/AntiB)
  if(projAbsB>1) nCons = projAbsB;                   // @@ what if this is a meson ?
  // Very important! Here the projectile 4M is recalculated in the ZLS (previously in LS)
  proj4M = aProjectile.Get4Momentum();               // 4-mom of theProjectileHadron inZLS
  G4double pz_per_projectile = proj4M.pz()/nCons;    // Momentum per nucleon (hadron?)
  // @@ use M_A/A instead of mProt ------------ M.K.
  G4double e_per_projectile = proj4M.e()/nCons+mProt;// @@ Add MassOfTargetProtonAtRest
  G4double vz = pz_per_projectile/e_per_projectile;  // Speed of the "oneNuclenCM"
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: Projectile4M="<<proj4M<<", Vz="<<vz<<", nC="
        <<nCons<<", pE="<<e_per_projectile<<G4endl;
#endif
  theCurrentVelocity.setZ(vz);                       // CM (w/r to one nucleon) velosity
  theNucleus.DoLorentzBoost(-theCurrentVelocity);    // BoostTgNucleus to"CM"
#ifdef edebug
  G4LorentzVector sum2=theNucleus.GetNucleons4Momentum();// Sum of Nucleons 4M in RotatedCM
  G4cout<<"-EMC-G4QFragmentation::Construct: AfterBoost, v="<<vz<<", Nuc4M="<<sum2<<G4endl;
#endif
  G4LorentzVector cmProjMom = proj4M;                // Copy the original proj4M in LS
  cmProjMom.boost(-theCurrentVelocity);              // Bring the LS proj4Mom to "CM"
  G4double kE=cmProjMom.e()-projM;
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: kE="<<kE<<G4endl;
#endif
  G4int maxCt=1;
  if(kE > 720.) maxCt=static_cast<int>(std::log(kE/270.)); // 270 MeV !
  // @@ The maxCuts can improve the performance at low energies
  //G4int maxCuts = 7;
  G4int maxCuts=std::min( 7 , std::max(1, maxCt) );
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: Proj4MInCM="<<cmProjMom<<", pPDG="<<pPDG<<G4endl;
#endif
  //
  // ---------->> Find collisions meeting collision conditions
  //
  G4QHadron* cmProjectile = new G4QHadron(pPDG,cmProjMom); // HipCopy of the CMProjectile
  // @@ Do not forget to delete the probability! 
  G4QProbability theProbability(pPDG);               // thePDG must be a data member
  G4double outerRadius = theNucleus.GetOuterRadius();// Get the nucleus frontiers
#ifdef debug
  G4cout<<"G4QFrag::Constr:OutR="<<outerRadius<<",mC="<<maxCuts<<",A="<<theNucleus<<G4endl;
#endif
  G4QHadron* pNucleon=0;
  // Check the reaction threshold 
  G4int theNA=theNucleus.GetA();
  G4LorentzVector pNuc4M=theNucleus.Get4Momentum()/theNA;
  G4double s_value = (cmProjMom + pNuc4M).mag2();          // Squared CM Energy of compound
  G4double ThresholdMass = projM + theNucleus.GetGSMass()/theNA;
#ifdef debug
  G4cout<<"G4QFrag::Construc: p4M="<<cmProjMom<<", tgN4M="<<pNuc4M<<", s="<<s_value<<", ThreshM="
        <<ThresholdMass<<G4endl;
#endif
  ModelMode = SOFT;                                  // NOT-Diffractive hadronization
  if (s_value < 0.)                                  // At ThP=0 is impossible(virtNucl)
  {
    G4cerr<<"***G4QFragmentation::Construct: s="<<s_value<<", pN4M="<<pNuc4M<<G4endl;
    G4Exception("G4QFragmentation::Construct:","72",FatalException,"LowEnergy(NegativeS)");
  }
  else if(s_value < mProt2)
  {
    theNucleus.StartLoop();                          // To get the same nucleon
    G4QHadron* aTarget=0;                            // Prototype of the target
    G4QHadron* pProjectile=0;                        // Prototype of the projectile
    G4QHadron* bNucleon=0;                           // Prototype of the best nucleon
    G4double   maxS=0.;                              // Maximum s found
    while((bNucleon=theNucleus.GetNextNucleon()))    // Loop over all nuclei to get theBest
    {
      G4LorentzVector cp4M=bNucleon->Get4Momentum(); // 4-mom of the current nucleon
      G4double cs=(cmProjMom + cp4M).mag2();         // Squared CM Energy of compound
      if(cs > maxS)                                  // Remember nucleon with the biggest s
      {
        maxS=cs;
        pNucleon=bNucleon;
      }
    }
    aTarget = new G4QHadron(*pNucleon);              // Copy selected nucleon for String
    pProjectile =cmProjectile;
    theNucleus.DoLorentzBoost(theCurrentVelocity);   // Boost theResNucleus toRotatedLS
    theNucleus.SubtractNucleon(pNucleon);            // Pointer to theUsedNucleon to delete
    theNucleus.DoLorentzBoost(-theCurrentVelocity);  // Boost theResNucleus back to CM
    G4QContent QQC=aTarget->GetQC()+pProjectile->GetQC(); // QContent of the compound
    G4LorentzVector Q4M=aTarget->Get4Momentum()+pProjectile->Get4Momentum(); // 4-mom of Q
    delete aTarget;
    delete pProjectile;
    //if(maxNuc>1)                                     // Absorb moreNucleons to theQuasmon
    //{
    //  for(G4int i=1; i<maxNuc; ++i)
    //  {
    //    pNucleon=theNucleus.GetNextNucleon();        // Get the next nucleon
    //    QQC+=pNucleon->GetQC();                      // Add it to the Quasmon
    //    Q4M+=pNucleon->Get4Momentum();
    //    theNucleus.DoLorentzBoost(theCurrentVelocity); // Boost theResNucleus toRotatedLS
    //    theNucleus.SubtractNucleon(pNucleon);        // Exclude the used nucleon from Nuc
    //    theNucleus.DoLorentzBoost(-theCurrentVelocity);// Boost theResNucleus back to CM
    //  }
    //}
    // 4-Mom should be converted to LS
    Q4M.boost(theCurrentVelocity);
    Q4M=toLab*Q4M;
    G4Quasmon* stringQuasmon = new G4Quasmon(QQC, Q4M);
    theQuasmons.push_back(stringQuasmon);
    theNucleus.DoLorentzBoost(theCurrentVelocity);   // BoostTheResidualNucleus toRotatedLS
    theNucleus.DoLorentzRotation(toLab);// Recove Z-direction in LS ("LS"->LS) for rNucleus
    return;
  }
  if (s_value < sqr(ThresholdMass))                    // --> Only diffractive interaction
  {
#ifdef debug
    G4cout<<"G4QFragmentation::Construct:*OnlyDiffraction*ThM="<<ThresholdMass<<">sqrt(s)="
          <<std::sqrt(s_value)<<" -> only Diffraction is possible"<<G4endl; // @@ Dif toQuasmon
#endif
    ModelMode = DIFFRACTIVE;
  }
  // Clean up all previous interactions and reset the counters
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: theIntSize="<<theInteractions.size()<<G4endl;
#endif
  std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
  theInteractions.clear();
  G4int totalCuts = 0;
  G4int attCnt=0;
  //G4int maxAtt=227;
  G4int maxAtt=27;
  G4double prEn=proj4M.e();                           // For mesons
  G4int proB=aProjectile.GetBaryonNumber();
  if     (proB>0) prEn-=aProjectile.GetMass();        // For baryons
  else if(proB<0) prEn+=mProt;                        // For anti-baryons
#ifdef debug
  G4double impactUsed = 0.;
  G4cout<<"G4QFragmentation::Construct: estimated energy, prEn="<<prEn<<G4endl;
#endif
  while(!theInteractions.size() && ++attCnt < maxAtt) // Till Interaction is created
  {
#ifdef debug
    G4cout<<"G4QFragmentation::Construct: *EnterTheInteractionLOOP*, att#"<<attCnt<<G4endl;
#endif
    // choose random impact parameter
    std::pair<G4double, G4double> theImpactParameter;
    theImpactParameter = theNucleus.ChooseImpactXandY(outerRadius);
    G4double impactX = theImpactParameter.first; 
    G4double impactY = theImpactParameter.second;
#ifdef debug
    G4cout<<"G4QFragmentation::Construct: Impact Par X="<<impactX<<", Y="<<impactY<<G4endl;
#endif
    G4double impar=std::sqrt(impactX*impactX+impactY*impactY);   
    G4int nA=theNucleus.GetA();
    G4double eflen=theNucleus.GetThickness(impar);   // EffectiveLength
    maxEn=eflen*stringTension;                       // max absorbed energy in IndUnits=MeV
    maxNuc=static_cast<int>(eflen*tubeDensity+0.5);  // max #0f involved nuclear nucleons
#ifdef debug
    G4cout<<"G4QFragment::Construct: pE="<<prEn<<" <? mE="<<maxEn<<", mN="<<maxNuc<<G4endl;
#endif
    if(prEn < maxEn)                                 // Create DIN interaction and go out
    {
      theNucleus.StartLoop();                        // Initialize newSelection forNucleons
      pNucleon=theNucleus.GetNextNucleon();          // Select a nucleon
      G4QHadron* aTarget = new G4QHadron(*pNucleon); // Copy selected nucleon for String
      G4QInteraction* anInteraction = new G4QInteraction(cmProjectile);
      anInteraction->SetTarget(aTarget); 
      anInteraction->SetNumberOfDINRCollisions(1);   // Consider this interaction as DINR
      theInteractions.push_back(anInteraction);      //--> now theInteractions not empty
      theNucleus.DoLorentzBoost(theCurrentVelocity); // Boost theResNucleus toRotatedLS
      theNucleus.SubtractNucleon(pNucleon);          // Pointer to the used nucleon
      theNucleus.DoLorentzBoost(-theCurrentVelocity);// Boost theResNucleus back to CM
#ifdef debug
      G4cout<<"G4QFragmentation::Construct: DINR interaction is created"<<G4endl;
#endif
      break;                                         // Break the WHILE of interactions
    }
    // LOOP over nuclei of the target nucleus to select collisions
    theNucleus.StartLoop();                          // To get the same nucleon
    G4int nucleonCount = 0;
#ifdef debug
    G4cout<<"G4QFragment::Construct:BeforeWhileOveNuc, A="<<nA<<",p4M="<<cmProjMom<<G4endl;
#endif
    while( (pNucleon=theNucleus.GetNextNucleon()) && nucleonCount<nA && totalCuts<maxCuts)
    {
      ++nucleonCount;
      // Needs to be moved to Probability class @@@
      s_value = (cmProjMom + pNucleon->Get4Momentum()).mag2();
#ifdef debug
      G4cout<<"G4QFrag::Constr:N# "<<nucleonCount<<", s="<<s_value<<", tgN4M="
            <<pNucleon->Get4Momentum()<<G4endl;
#endif
      if(s_value<=10000.)
      {
#ifdef debug
        G4cout<<"G4QFragmentation::Construct: SKIP, s<.01 GeV^2, p4M="<<cmProjMom
              <<",t4M="<<pNucleon->Get4Momentum()<<G4endl;
#endif
        continue;
      }
#ifdef sdebug
      G4cout<<"G4QFragmentation::Construct:LOOPovNuc,nC="<<nucleonCount<<", s="<<s_value<<G4endl;
      G4cout<<"G4QFragmentation::Construct:LOOPovNuc, R="<<pNucleon->GetPosition()<<G4endl;
#endif
      G4double Distance2 = sqr(impactX - pNucleon->GetPosition().x()) +
                           sqr(impactY - pNucleon->GetPosition().y());
#ifdef sdebug
      G4cout<<"G4QFragmentation::Construct: s="<<s_value<<", D2="<<Distance2<<G4endl;
#endif
      G4double Probability = theProbability.GetPomInelProbability(s_value, Distance2); // PomINEL
      // test for inelastic collision
#ifdef sdebug
      G4cout<<"G4QFragmentation::Construct: Probubility="<<Probability<<G4endl;
#endif
      G4double rndNumber = G4UniformRand();           // For the printing purpose
      // ModelMode = DIFFRACTIVE;
#ifdef sdebug
      G4cout<<"G4QFragmentation::Construct: NLOOP prob="<<Probability<<", rndm="<<rndNumber
            <<", d="<<std::sqrt(Distance2)<<G4endl;
#endif
      if (Probability > rndNumber) // Inelastic (diffractive or soft) interaction (JOB)
      {
        G4QHadron* aTarget = new G4QHadron(*pNucleon);// Copy for String (ValgrindComplain)
#ifdef edebug
        G4cout<<"--->EMC-->G4QFragmentation::Construct: Target Nucleon is filled, 4M/PDG="
              <<aTarget->Get4Momentum()<<aTarget->GetPDGCode()<<G4endl;
#endif
        // Now the energy of the nucleons must be updated in CMS
        theNucleus.DoLorentzBoost(theCurrentVelocity);// Boost theResNucleus toRotatedLS
        theNucleus.SubtractNucleon(pNucleon);         // Pointer to the used nucleon
        theNucleus.DoLorentzBoost(-theCurrentVelocity);// Boost theResNucleus back to CM
        if((theProbability.GetPomDiffProbability(s_value,Distance2)/Probability >
            G4UniformRand() && ModelMode==SOFT ) || ModelMode==DIFFRACTIVE)
        { 
          // --------------->> diffractive interaction @@ IsSingleDiffractive called once
          if(IsSingleDiffractive()) ExciteSingDiffParticipants(cmProjectile, aTarget);
          else                          ExciteDiffParticipants(cmProjectile, aTarget);
          G4QInteraction* anInteraction = new G4QInteraction(cmProjectile);
          anInteraction->SetTarget(aTarget); 
          anInteraction->SetNumberOfDiffractiveCollisions(1); // Why not increment? M.K.?
          theInteractions.push_back(anInteraction);   //--> now theInteractions not empty
          // @@ Why not breake the NLOOP, if only one diffractive can happend?
          totalCuts++;                               // UpdateOfNucleons in't necessary
#ifdef debug
          G4cout<<"G4QFragmentation::Construct:NLOOP DiffInteract, tC="<<totalCuts<<G4endl;
#endif
        }
        else
        {
          // ---------------->> nondiffractive = soft interaction
          // sample nCut+1 (cut Pomerons) pairs of strings can be produced
          G4int nCut;                                // Result in a chosen number of cuts
          G4double* running = new G4double[nCutMax]; // @@ This limits the max cuts
          for(nCut = 0; nCut < nCutMax; nCut++)      // Calculates multiCut probabilities
          {
            running[nCut]= theProbability.GetCutPomProbability(s_value, Distance2, nCut+1);
            if(nCut) running[nCut] += running[nCut-1];// Sum up with the previous one
          }
          G4double random = running[nCutMax-1]*G4UniformRand();
          for(nCut = 0; nCut < nCutMax; nCut++) if(running[nCut] > random) break;
          delete [] running;
#ifdef debug
          G4cout<<"G4QFragmentation::Construct: NLOOP-Soft Chosen nCut="<<nCut<<G4endl;
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
          G4cout<<"G4QFragmentation::Construct:NLOOP SoftInteract, tC="<<totalCuts<<G4endl;
          impactUsed=Distance2;
#endif
        }
      }
    } // End of While over nucleons
    // When nucleon count is incremented, the LOOP stops, so nucleonCount==1 always!
#ifdef debug
    G4cout<<"G4QFragmentation::Construct: NUCLEONCOUNT="<<nucleonCount<<G4endl;
#endif
  }
  G4int nInt=theInteractions.size();
#ifdef debug
  G4cout<<"G4QFrag::Con:CUT="<<totalCuts<<",ImpPar="<<impactUsed<<",#ofInt="<<nInt<<G4endl;
#endif
  // --- Use this line ---
  if(!nInt || (nInt==1 && theInteractions[0]->GetNumberOfDINRCollisions()==1)) // @@ ?? @@
  {
    G4QHadron* aTarget=0;
    G4QHadron* pProjectile=0;
    if(nInt)                                         // Take Targ/Proj from the Interaction
    {
	aTarget=theInteractions[0]->GetTarget();
	pProjectile=theInteractions[0]->GetProjectile();
      delete theInteractions[0];
      theInteractions.clear();
    }
    else                                             // Create a new target nucleon
    {
      theNucleus.StartLoop();                        // To get the same nucleon
      pNucleon=theNucleus.GetNextNucleon();          // Get the nucleon to create
      aTarget = new G4QHadron(*pNucleon);            // Copy selected nucleon for String
      pProjectile =cmProjectile;
      theNucleus.DoLorentzBoost(theCurrentVelocity); // Boost theResNucleus toRotatedLS
      theNucleus.SubtractNucleon(pNucleon);          // Pointer to theUsedNucleon to delete
      theNucleus.DoLorentzBoost(-theCurrentVelocity);// Boost theResNucleus back to CM
    }
    G4QContent QQC=aTarget->GetQC()+pProjectile->GetQC(); // QContent of the compound
    G4LorentzVector Q4M=aTarget->Get4Momentum()+pProjectile->Get4Momentum(); // 4-mom of Q
    delete aTarget;
    delete pProjectile;
    //if(maxNuc>1)                                     // Absorb moreNucleons to theQuasmon
    //{
    //  for(G4int i=1; i<maxNuc; ++i)
    //  {
    //    pNucleon=theNucleus.GetNextNucleon();        // Get the next nucleon
    //    QQC+=pNucleon->GetQC();                      // Add it to the Quasmon
    //    Q4M+=pNucleon->Get4Momentum();
    //    theNucleus.DoLorentzBoost(theCurrentVelocity); // Boost theResNucleus toRotatedLS
    //    theNucleus.SubtractNucleon(pNucleon);        // Exclude the used nucleon from Nuc
    //    theNucleus.DoLorentzBoost(-theCurrentVelocity);// Boost theResNucleus back to CM
    //  }
    //}
    // 4-Mom should be converted to LS
    Q4M.boost(theCurrentVelocity);
    Q4M=toLab*Q4M;
    G4Quasmon* stringQuasmon = new G4Quasmon(QQC, Q4M);
    theQuasmons.push_back(stringQuasmon);
    theNucleus.DoLorentzBoost(theCurrentVelocity);   // BoostTheResidualNucleus toRotatedLS
    theNucleus.DoLorentzRotation(toLab);// Recove Z-direction in LS ("LS"->LS) for rNucleus
    return;
  }
  //
  // ------------------ now build the parton pairs for the strings ------------------
  //
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: Before PartPairCreation nInt="<<nInt<<G4endl;
#endif
  for(G4int i=0; i<nInt; i++)
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
    G4cout<<"-EMC-G4QFragmentation::Construct: Interaction#"<<i<<",dP4M="
          <<projH->Get4Momentum()-pSP<<",dT4M="<<targH->Get4Momentum()-tSP<<G4endl;
#endif
  }  
  // 
  // ------->> make soft collisions (ordering is vital)
  //
  G4QInteractionVector::iterator it;
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: Creation ofSoftCollisionPartonPair STARTS"<<G4endl;
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
    G4cout<<"G4QFragmentation::Construct: #0f SOFT collisions ="<<nSoftCollisions<<G4endl;
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
        G4cout<<"--->G4QFragmentation::Construct: SOFT, 2 parton pairs are filled"<<G4endl;
#endif
      }  
      delete *it;
      it=theInteractions.erase(it);      // Soft interactions are converted & erased
      if( it != theInteractions.begin() )// To avoid going below begin() (for Windows)
      {
        it--;
        rep=false;
#ifdef debug
        G4cout<<"G4QFragmentation::Construct: *** Decremented ***"<<G4endl;
#endif
      }
      else
      {
        rep=true;
#ifdef debug
        G4cout<<"G4QFragmentation::Construct: *** IT Begin ***"<<G4endl;
#endif
        break;
      }
    }
    else rep=false;
#ifdef debug
    G4cout<<"G4QFragmentation::Construct: #0fSC="<<nSoftCollisions<<", r="<<rep<<G4endl;
#endif
   }
#ifdef debug
   G4cout<<"G4QFragmentation::Construct: *** IT While *** , r="<<rep<<G4endl;
#endif
  }
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: -> Parton pairs for SOFT strings are made"<<G4endl;
#endif  
  //
  // ---------->> make the rest as the diffractive interactions
  //
  for(unsigned i = 0; i < theInteractions.size(); i++) // Interactions are reduced bySoft
  {
    // The double or single diffraction is defined by the presonce of proj/targ partons
    G4QInteraction* anIniteraction = theInteractions[i];
    G4QPartonPair* aPartonPair;
#ifdef debug
    G4cout<<"G4QFragmentation::Construct: CreationOfDiffractivePartonPairs, i="<<i<<G4endl;
#endif
    // the projectile diffraction parton pair is created first
    G4QHadron* pProjectile = anIniteraction->GetProjectile();
    G4QParton* aParton = pProjectile->GetNextParton();
    if (aParton)
    {
      aPartonPair = new G4QPartonPair(aParton, pProjectile->GetNextAntiParton(), 
                                      G4QPartonPair::DIFFRACTIVE,
                                      G4QPartonPair::PROJECTILE);
      thePartonPairs.push_back(aPartonPair);
#ifdef debug
      G4cout<<"G4QFragmentation::Construct: proj Diffractive PartonPair is filled"<<G4endl;
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
      G4cout<<"G4QFragmentation::Constr: target Diffractive PartonPair is filled"<<G4endl;
#endif
    }
  }
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: DiffractivePartonPairs are created"<<G4endl;
#endif  
  //
  // ---------->> clean-up  Interactions and cmProjectile, if necessary
  //
  std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
  theInteractions.clear();
  delete cmProjectile;
#ifdef debug
  G4cout<<"G4QFragmentation::Construct: Temporary objects are cleaned up"<<G4endl;
#endif  
  // This function prepares theBoost for transformation of secondaries to LS (-ProjRot!)
  theNucleus.DoLorentzBoost(theCurrentVelocity);// Boost theResidualNucleus to RotatedLS
  // @@ Nucleus isn't completely in LS, it needs the toZ (-ProjRot) rotation to consE/M
#ifdef debug
  G4cout<<"--------->>G4QFragmentation::Construct: ------->> Strings are created "<<G4endl;
#endif
  G4QPartonPair* aPair;
  G4QString* aString=0;
  while(thePartonPairs.size()) // @@ At present noDifference in stringBuild (? M.K.)
  {
    aPair = thePartonPairs.back();           // Get the parton pair
    thePartonPairs.pop_back();               // Clean up thePartonPairPointer in the Vector
#ifdef debug
    G4cout<<"G4QFragmentation::Construct: StringType="<<aPair->GetCollisionType()<<G4endl;
#endif
    aString= new G4QString(aPair);
#ifdef debug
    G4cout<<"G4QFragmentation::Construct:NewString4M="<<aString->Get4Momentum()<<G4endl;
#endif
    aString->Boost(theCurrentVelocity);       // ! Strings are moved to ZLS when pushed !
    strings.push_back(aString);
    stringsInitted=true;
    delete aPair;
  } // End of the String Creation LOOP
#ifdef edebug
  G4LorentzVector sum=theNucleus.Get4Momentum();// Nucleus 4Mom in rotatedLS
  G4int rChg=totChg-theNucleus.GetZ();
  G4int rBaN=totBaN-theNucleus.GetA();
  G4int nStrings=strings.size();
  G4cout<<"-EMC-G4QFragmentation::Construct:#ofString="<<nStrings<<",tNuc4M="<<sum<<G4endl;
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
    G4cout<<"-EMC-G4QFragmentation::Construct: String#"<<i<<", 4M="<<strI4M<<",LPDG="<<LPDG
          <<LQC<<",RPDG="<<RPDG<<RQC<<", Ch="<<sChg<<", BN="<<sBaN<<G4endl;
  }
  G4cout<<"-EMC-G4QFragm::Constr: r4M="<<sum-totZLS4M<<",rC="<<rChg<<",rB="<<rBaN<<G4endl;
#endif
  if(!stringsInitted)
  {
    G4cerr<<"******G4QFragmentation::Construct:***** No strings are created *****"<<G4endl;
    G4Exception("G4QFragmentation::Construct:","72",FatalException,"NoStrings're created");
  }
#ifdef debug
  G4cout<<"G4QFragmentation::Constr: BeforeRotation, #0fStrings="<<strings.size()<<G4endl;
#endif
  //
  // ---------------- At this point the strings must be already created in "LS" -----------
  //
  for(unsigned astring=0; astring < strings.size(); astring++)
            strings[astring]->LorentzRotate(toLab); // Recove Z-direction in LS ("LS"->LS)
  theNucleus.DoLorentzRotation(toLab); // Recove Z-direction in LS ("LS"->LS) for rNucleus
  // Now everything is in LS system
#ifdef edebug
  G4LorentzVector sm=theNucleus.Get4Momentum();    // Nucleus 4Mom in LS
  G4int rCg=totChg-theNucleus.GetZ();
  G4int rBC=totBaN-theNucleus.GetA();
  G4int nStrs=strings.size();
  G4cout<<"-EMCLS-G4QFragmentation::Constr: #ofS="<<nStrings<<",tNuc4M(E=M)="<<sum<<G4endl;
  for(G4int i=0; i<nStrs; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    sm+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    rCg-=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    rBC-=sBaN;
    G4cout<<"-EMCLS-G4QFragm::Construct:String#"<<i<<",4M="<<strI4M<<strI4M.m()<<",Charge="
          <<sChg<<",BaryN="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-...G4QFragm::Constr:r4M="<<sm-totLS4M<<",rC="<<rCg<<",rB="<<rBC<<G4endl;
#endif
  //
  // --- Strings are created, but we should try to get rid of negative mass strings -----
  //
  SwapPartons();
#ifdef edebug
  sm=theNucleus.Get4Momentum();    // Nucleus 4Mom in LS
  rCg=totChg-theNucleus.GetZ();
  rBC=totBaN-theNucleus.GetA();
  nStrs=strings.size();
  G4cout<<"-EMCLS-G4QFrag::Constr:AfterSwap #ofS="<<nStrings<<",tNuc4M(E=M)="<<sum<<G4endl;
  for(G4int i=0; i<nStrs; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    sm+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    rCg-=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    rBC-=sBaN;
    G4cout<<"-EMCLS-G4QFragm::Construct:String#"<<i<<",4M="<<strI4M<<strI4M.m()<<",Charge="
          <<sChg<<",BaryN="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-...G4QFragm::Constr:r4M="<<sm-totLS4M<<",rC="<<rCg<<",rB="<<rBC<<G4endl;
#endif
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
    G4cout<<"G4QFrag::Constr:NeedsFusion? cLPDG="<<cLPDG<<",cRPDG="<<cRPDG<<",cM(cM2If<0)="
          <<cSM<<",c4M"<<cS4M<<G4endl;
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
      G4cout<<"G4QFrag::Const:*IsItGood? realM="<<std::sqrt(cSM2)<<" > GSM="<<miM<<G4endl;
#endif
      if(std::sqrt(cSM2) > miM) bad=false;           // String is OK
    }
    if(bad)                                          // String should be merged with others
    {
#ifdef debug
      G4cout<<"G4QFrag::Const:TryFuse,L1="<<L1<<",L2="<<L2<<",R1="<<R1<<",R2="<<R2<<G4endl;
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
        G4LorentzVector sS4M=(*pst)->Get4Momentum(); // Partner's 4-momentum
        G4LorentzVector pS4M=sS4M+cS4M;              // Summed 4-momentum
        G4int nLPDG=0;                               // new Left (like in theStringPartner)
        G4int nRPDG=0;                               // new Right(like in theStringPartner)
        G4double pExcess=-DBL_MAX;                   // Prototype of the excess
        G4double pSM2=pS4M.m2();                     // Squared mass of the Fused Strings
#ifdef debug
        G4cout<<"->G4QFragm::Construct: sum4M="<<pS4M<<",M2="<<pSM2<<",p4M="<<sS4M<<G4endl;
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
        G4cout<<"G4QFragm::Construct: Partner/w pLPDG="<<pLPDG<<", pRPDG="<<pRPDG<<", pM2="
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
          else G4cout<<"G4QFragmentation::Construct:2(QaQ+QDiQ/aQaDiQ) Can't fuse"<<G4endl;
#endif
        }
        else if(cST==3 && pST==2)                  // QDiQ/aQaDiQ + QaQ = QDiQ/aQaDiQ (s)
        {
          tf=3;
          if     (cLPDG > 7 &&
                  ( (pLPDG<0 && (-pLPDG==L1 || -pLPDG==L2 || -pLPDG==cRPDG) ) ||
                    (pRPDG<0 && (-pRPDG==L1 || -pRPDG==L2 || -pRPDG==cRPDG) )
                  )
                 ) af=1;
          else if(cRPDG > 7 &&
                  ( (pLPDG<0 && (-pLPDG==R1 || -pLPDG==R2 || -pLPDG==cLPDG) ) ||
                    (pRPDG<0 && (-pRPDG==R1 || -pRPDG==R2 || -pRPDG==cLPDG) )
                  )
                 ) af=2;
          else if(cLPDG <-7 &&
                  ( (pLPDG>0 && ( pLPDG==L1 || pLPDG==L2 || pLPDG==-cRPDG) ) ||
                    (pRPDG>0 && ( pRPDG==L1 || pRPDG==L2 || pRPDG==-cRPDG) )
                  )
                 ) af=3;
          else if(cRPDG <-7 &&
                  ( (pLPDG>0 && ( pLPDG==R1 || pLPDG==R2 || pLPDG==-cLPDG) ) ||
                    (pRPDG>0 && ( pRPDG==R1 || pRPDG==R2 || pRPDG==-cLPDG) )
                  )
                 ) af=4;
#ifdef debug
          else G4cout<<"G4QFragmentation::Construct:3(QDiQ/aQaDiQ+QaQ) Can't fuse"<<G4endl;
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
          else G4cout<<"-G4QFragmentation::Construct: 4 (QaQ+aQDiQDiQ) Can't fuse"<<G4endl;
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
          else G4cout<<"-G4QFragmentation::Construct: 5 (aQDiQDiQ+QaQ) Can't fuse"<<G4endl;
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
          else G4cout<<"-G4QFragmentation::Construct: 6 (QDiQ+aQaDiQ) Can't fuse"<<G4endl;
#endif
        }
#ifdef debug
        G4cout<<"G4QFragmentation::Const: ***Possibility***, tf="<<tf<<", af="<<af<<G4endl;
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
                   if(-pLPDG==cRPDG)
                   {
                     nLPDG=cLPDG;
                     nRPDG=pRPDG;
                   }
                   else
                   {
                     if     (pRPDG > cRPDG) nRPDG=pRPDG*1000+cRPDG*100+1;
                     else if(pRPDG < cRPDG) nRPDG=cRPDG*1000+pRPDG*100+1;
                     else                   nRPDG=cRPDG*1000+pRPDG*100+3;
                     if  (-pLPDG == L1)     nLPDG=L2;
                     else                   nLPDG=L1; // -pLPDG == L2
                   }
                 }
                 else // pRPDG < 0
                 {
                   order=-1;
                   if(-pRPDG==cRPDG)
                   {
                     nLPDG=cLPDG;
                     nRPDG=pLPDG;
                   }
                   else
                   {
                     if     (pLPDG > cRPDG) nLPDG=pLPDG*1000+cRPDG*100+1;
                     else if(pLPDG < cRPDG) nLPDG=cRPDG*1000+pLPDG*100+1;
                     else                   nLPDG=cRPDG*1000+pLPDG*100+3;
                     if  (-pRPDG == L1)     nRPDG=L2;
                     else                   nRPDG=L1; // -pRPDG == L2
                   }
                 }
                 break;
               case 2: // ....................... cRPDG > 7
                 if(pLPDG < 0)
                 {
                   order=-1;
                   if(-pLPDG==cLPDG)
                   {
                     nLPDG=pRPDG;
                     nRPDG=cRPDG;
                   }
                   else
                   {
                     if     (pRPDG > cLPDG) nRPDG=pRPDG*1000+cLPDG*100+1;
                     else if(pRPDG < cLPDG) nRPDG=cLPDG*1000+pRPDG*100+1;
                     else                   nRPDG=cLPDG*1000+pRPDG*100+3;
                     if  (-pLPDG == R1)     nLPDG=R2;
                     else                   nLPDG=R1; // -pLPDG == R2
                   }
                 }
                 else // pRPDG < 0
                 {
                   order= 1;
                   if(-pRPDG==cLPDG)
                   {
                     nLPDG=pLPDG;
                     nRPDG=cRPDG;
                   }
                   else
                   {
                     if     (pLPDG > cLPDG) nLPDG=pLPDG*1000+cLPDG*100+1;
                     else if(pLPDG < cLPDG) nLPDG=cLPDG*1000+pLPDG*100+1;
                     else                   nLPDG=cLPDG*1000+pLPDG*100+3;
                     if  (-pRPDG == R1)     nRPDG=R2;
                     else                   nRPDG=R1; // -pRPDG == R2
                   }
                 }
                 break;
               case 3: // ....................... cLPDG <-7 (cRPDG <0)
                 if(pLPDG > 0)
                 {
                   order= 1;
                   if(pLPDG==-cRPDG)
                   {
                     nLPDG=cLPDG;
                     nRPDG=pRPDG;
                   }
                   else
                   {
                     if     (pRPDG < cRPDG) nRPDG=pRPDG*1000+cRPDG*100-1;
                     else if(pRPDG > cRPDG) nRPDG=cRPDG*1000+pRPDG*100-1;
                     else                   nRPDG=cRPDG*1000+pRPDG*100-3;
                     if  ( pLPDG == L1)     nLPDG=-L2;
                     else                   nLPDG=-L1; // pLPDG == L2
                   }
                 }
                 else // pRPDG > 0
                 {
                   order=-1;
                   if(pRPDG==-cRPDG)
                   {
                     nLPDG=cLPDG;
                     nRPDG=pLPDG;
                   }
                   else
                   {
                     if     (pLPDG < cRPDG) nLPDG=pLPDG*1000+cRPDG*100-1;
                     else if(pLPDG > cRPDG) nLPDG=cRPDG*1000+pLPDG*100-1;
                     else                   nLPDG=cRPDG*1000+pLPDG*100-3;
                     if  ( pRPDG == L1)     nRPDG=-L2;
                     else                   nRPDG=-L1; // pRPDG == L2
                   }
                 }
                 break;
               case 4: // ....................... cRPDG <-7 (cLPDG <0)
                 if(pLPDG > 0)                       // pRPDG & cLPDG are anti-quarks
                 {
                   order=-1;
                   if(pLPDG==-cLPDG)
                   {
                     nLPDG=pRPDG;
                     nRPDG=cRPDG;
                   }
                   else
                   {
                     if     (pRPDG < cLPDG) nRPDG=pRPDG*1000+cLPDG*100-1;
                     else if(pRPDG > cLPDG) nRPDG=cLPDG*1000+pRPDG*100-1;
                     else                   nRPDG=cLPDG*1000+pRPDG*100-3;
                     if  ( pLPDG == R1)     nLPDG=-R2;
                     else                   nLPDG=-R1; // pLPDG == R2
                   }
                 }
                 else // pRPDG > 0
                 {
                   order= 1;
                   if(pRPDG==-cLPDG)
                   {
                     nLPDG=pLPDG;
                     nRPDG=cRPDG;
                   }
                   else
                   {
                     if     (pLPDG < cLPDG) nLPDG=pLPDG*1000+cLPDG*100-1;
                     else if(pLPDG > cLPDG) nLPDG=cLPDG*1000+pLPDG*100-1;
                     else                   nLPDG=cLPDG*1000+pLPDG*100-3;
                     if  ( pRPDG == R1)     nRPDG=-R2;
                     else                   nRPDG=-R1; // pRPDG == R2
                   }
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
          if(!order) G4cerr<<"-Warning-G4QFrag::Constr: t="<<tf<<", a="<<af<<", cL="<<cLPDG
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
              aLPDG=0;
              aRPDG=0;
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
                G4cout<<"G4QFragmentation::Const:aLPDG="<<aLPDG<<", aRPDG="<<aRPDG<<G4endl;
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
              G4cout<<"G4QFr::Con: LPDG="<<nLPDG<<",RPDG="<<nRPDG<<",QC="<<newStQC<<G4endl;
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
        G4cout<<"G4QFragmentation::Const:cS4M="<<cS4M<<" fused/w pS4M="<<pL4M+pR4M<<G4endl;
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
        G4cout<<"G4QFragmentation::Construct:Created,4M="<<ss4M<<",m2="<<ss4M.m2()<<G4endl;
#endif
        if( ist != strings.begin() ) // To avoid going below begin() (for Windows)
        {
          ist--;
          con=false;
#ifdef debug
          G4cout<<"G4QFragmentation::Construct: *** IST Decremented ***"<<G4endl;
#endif
        }
        else
        {
          con=true;
#ifdef debug
          G4cout<<"G4QFragmentation::Construct: *** IST Begin ***"<<G4endl;
#endif
          break;
        }
      } // End of the IF(the best partnerString candidate was found)
      else
      {
#ifdef debug
        G4cout<<"-Warning-G4QFragm::Const:S4M="<<cS4M<<",M2="<<cSM2<<" Leave ASIS"<<G4endl;
#endif
        ++problem;
        con=false;
      }
    } // End of IF
    else con=false;
   } // End of loop over ist iterator
#ifdef debug
   G4cout<<"G4QFragmentation::Construct: *** IST While *** , con="<<con<<G4endl;
#endif
  } // End of "con" while 
#ifdef edebug
  // This print has meaning only if something appear between it and the StringFragmLOOP
  G4LorentzVector t4M=theNucleus.Get4Momentum();    // Nucleus 4Mom in LS
  G4int rC=totChg-theNucleus.GetZ();
  G4int rB=totBaN-theNucleus.GetA();
  G4int nStr=strings.size();
  G4cout<<"-EMCLS-G4QFr::Const: AfterSUPPRESION #ofS="<<nStr<<",tNuc4M(E=M)="<<sum<<G4endl;
  for(G4int i=0; i<nStr; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    t4M+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    rC-=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    rB-=sBaN;
    G4cout<<"-EMCLS-G4QFragm::Construct: St#"<<i<<", 4M="<<strI4M<<", M="<<strI4M.m()
          <<", C="<<sChg<<", B="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFragm::Construct:r4M="<<t4M-totLS4M<<",rC="<<rC<<",rB="<<rB<<G4endl;
#endif
  //
  // --- If a problem is foreseen then the DiQaDiQ strings should be reduced if possible --
  //
#ifdef debug
    G4cout<<"G4QFragmentation::Construct: problem="<<problem<<G4endl;
#endif
  if(problem)
  {
    G4int nOfStr=strings.size();
#ifdef debug
    G4cout<<"G4QFragmentation::Construct:SecurityDiQaDiQReduction,#OfStr="<<nOfStr<<G4endl;
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
        G4cout<<"G4QFragmentation::Constr:TrySelfReduString,L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
        if( cLeft->ReduceDiQADiQ(cLeft, cRight) ) // DiQ-aDiQ pair was successfully reduced
        {
          sPDG=cLeft->GetPDGCode();
          nPDG=cRight->GetPDGCode();
#ifdef debug
          G4cout<<"+G4QFragm::Const:#"<<astring<<" Reduced, L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
        }
#ifdef debug
        else G4cout<<"--*--G4QFragm::Const:#"<<astring<<" DQ-aDQ reduction Failed"<<G4endl;
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
  G4LorentzVector u4M=theNucleus.Get4Momentum();    // Nucleus 4Mom in LS
  G4int rCh=totChg-theNucleus.GetZ();
  G4int rBa=totBaN-theNucleus.GetA();
  G4int nStri=strings.size();
  G4cout<<"-EMCLS-G4QFr::Const: FinalConstruct, #ofSt="<<nStri<<",tN4M(E=M)="<<t4M<<G4endl;
  for(G4int i=0; i<nStri; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    u4M+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    rCh-=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    rBa-=sBaN;
    G4cout<<"-EMCLS-G4QFragm::Construct: St#"<<i<<", 4M="<<strI4M<<", M="<<strI4M.m()
          <<", C="<<sChg<<", B="<<sBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFragm::Construct:r4M="<<u4M-totLS4M<<",rC="<<rCh<<",rB="<<rBa<<G4endl;
#endif
} // End of the Constructer

G4QFragmentation::~G4QFragmentation()
{
  std::for_each(strings.begin(), strings.end(), DeleteQString() );
}

G4QHadronVector* G4QFragmentation::Fragment()
{ // This is the member function fragmenting Strings & Quasmons (in nuclear matter)
  static const G4QPDGCode nQPDG(2112);
  static const G4double   mProt = G4QPDGCode(2212).GetMass(); // Mass of proton
  static const G4double   mNeut = G4QPDGCode(2112).GetMass(); // Mass of neutron
  static const G4double   mPiCh = G4QPDGCode(211).GetMass();  // Mass of chgdPion
  static const G4double   mPiZr = G4QPDGCode(111).GetMass();  // Mass of neutrPion
  //static const G4double   mHe3 = G4QPDGCode(2112).GetNuclMass(2,1,0);
  static const G4LorentzVector  nul4M(0.,0.,0.,0.);          // Zero (vacuum) 4M
  //static const G4double eps=0.003;
#ifdef debug
  G4cout<<"*******>G4QFragmentation::Fragment: ***Called***, Res="<<theResult<<G4endl;
#endif
  G4int striNum=strings.size();                             // Find out if there're strings
  G4int hadrNum=theResult->size();                          // Find out if there're hadrons
#ifdef edebug
  G4int nQm=theQuasmons.size();
  G4LorentzVector totLS4M=theNucleus.Get4Momentum();        // Nucleus 4Mom in LS
  G4int totChg=theNucleus.GetZ();
  G4int totBaN=theNucleus.GetA();
  G4cout<<"-EMCLS-G4QF::Fragment: CHECKRecovery, #ofS="<<striNum<<", #Nuc4M(E=M)="<<totLS4M
        <<",#Q="<<nQm<<",#H="<<hadrNum<<G4endl;
  for(G4int i=0; i < striNum; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    totLS4M+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    totChg+=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    totBaN+=sBaN;
    G4cout<<"-EMCLS-G4QFragm::Fragment: String#"<<i<<", 4M="<<strI4M<<", M="<<strI4M.m()
          <<", C="<<sChg<<", B="<<sBaN<<G4endl;
  }
  for(G4int i=0; i < nQm; i++)
  {
    G4LorentzVector hI4M=theQuasmons[i]->Get4Momentum();
    totLS4M+=hI4M;
    G4int hChg=theQuasmons[i]->GetCharge();
    totChg+=hChg;
    G4int hBaN=theQuasmons[i]->GetBaryonNumber();
    totBaN+=hBaN;
    G4cout<<"-EMCLS-G4QFragmentation::Fragment: Quasmon#"<<i<<", 4M="<<hI4M<<", C="<<hChg
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
    G4cout<<"-EMCLS-G4QFr::Fragment:H#"<<i<<",4M="<<hI4M<<",C="<<hChg<<",B="<<hBaN<<G4endl;
  }
#endif
#ifdef debug
  G4cout<<"***>G4QFragmentation::Fragment: #OfStr="<<striNum<<", #OfRes="<<hadrNum<<G4endl;
#endif
  if(!striNum && hadrNum)                                   // Quasi-elastic or decoupled p
  {
#ifdef debug
    G4cout<<"***>G4QFragmentation::Fragment:**Quasi-Elastic**,#OfResult="<<hadrNum<<G4endl;
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
    G4LorentzVector r4M=theNucleus.Get4Momentum();          // Nucleus 4-momentum in LS
    G4int rPDG=theNucleus.GetPDG();                         // Nuclear PDG
    G4QHadron* resNuc = new G4QHadron(rPDG,r4M);            // Nucleus -> Hadron
    theResult->push_back(resNuc);                           // Fill the residual nucleus
  }
  G4int nQuas=theQuasmons.size();                           // Size of the Quasmon OUTPUT
  G4int theRS=theResult->size();                            // Size of Hadron Output by now
#ifdef debug
  G4cout<<"***>G4QFragmentation::Fragment:beforeEnv,#OfQ="<<nQuas<<",#OfR="<<theRS<<G4endl;
#endif
  if(nQuas && theRS)
  {
    G4QHadron* resNuc = (*theResult)[theRS-1];              // Pointer to Residual Nucleus
    G4LorentzVector resNuc4M = resNuc->Get4Momentum();      // 4-Momentum of the Nucleuz
    G4int           resNucPDG= resNuc->GetPDGCode();        // PDG Code of the Nucleus
    if(resNucPDG==90000000 || resNuc4M.m2()<800000.)        // m_N^2 = 880000 MeV^2
    {
      resNuc4M=G4LorentzVector(0.,0.,0.,0.);
      if(resNucPDG == 90000000) resNuc->Set4Momentum(resNuc4M);
    }
#ifdef edebug
    G4int rnChg=resNuc->GetCharge();
    G4int rnBaN=resNuc->GetBaryonNumber();
#endif
    G4QNucleus      theEnv(resNucPDG);                      // NucleusHadron->NucleusAtRest
    delete resNuc;                                          // Delete resNucleus as aHadron
    theResult->pop_back();                                  // Exclude the nucleus from HV
    --theRS;                                                // Reduce the OUTPUT by theNucl
#ifdef debug
    G4cout<<"G4QFragmentation::Fragment:#OfRemainingHadron="<<theRS<<",A="<<theEnv<<G4endl;
#endif
    // Now we need to be sure that the compound nucleus is heavier than the Ground State
    for(G4int j=theRS-1; j>-2; --j)                         // try to reach M_compound>M_GS
    {
      G4LorentzVector qsum4M=resNuc4M;                      // Proto compound 4-momentum
      G4QContent qsumQC=theEnv.GetQCZNS();                  // Proto compound Quark Content
#ifdef debug
      G4cout<<"G4QFragmentation::Fragm:rN4M"<<qsum4M<<qsum4M.m()<<",rNQC="<<qsumQC<<G4endl;
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
#ifdef debug
        G4cout<<"G4QFr::Fr:Q#"<<i<<",Q4M="<<cur4M<<",QQC="<<curQC<<",sQC="<<qsumQC<<G4endl;
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
#ifdef debug
      G4cout<<"G4QFragmentation::Fragment: QC="<<qsumQC<<",PDG="<<miPDG<<G4endl;
#endif
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
#ifdef debug
      G4cout<<"G4QFragmentation::Fragment: PDG="<<miPDG<<",rM="<<reM<<",GSM="<<gsM<<G4endl;
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
#ifdef debug
        G4cout<<"G4QFragm::Fragm: H#"<<j<<",hQC="<<hQC<<",hPDG="<<cH->GetPDGCode()<<G4endl;
#endif
      }
      else
      {
        G4cerr<<"*G4QFr::Fr:PDG="<<miPDG<<",M="<<reM<<",GSM="<<gsM<<",QC="<<qsumQC<<G4endl;
        G4Exception("G4QFragmentation::Fragment:","27",FatalException,"Can't recover GSM");
      }
    }
    G4double nucE=resNuc4M.e();                             // Total energy of the nuclEnv
    if(nucE < 1.E-12) nucE=0.;                              // Computer accuracy safety
    else if(resNucPDG==22 && nQuas==1)                      // NuclEnv for nQ=1 is a photon
    {
      G4Quasmon* aQuasm=theQuasmons[0];                     // the only Quasmon
      aQuasm->Set4Momentum(aQuasm->Get4Momentum()+resNuc4M);// add the gammaEnv to the Q
      nucE=0.;
    }
    G4ThreeVector   nucVel(0.,0.,0.);                       // Proto of the NucleusVelocity
    G4QHadronVector* output=0;                              // NucleusFragmentation Hadrons
    G4QEnvironment* pan= new G4QEnvironment(theEnv);        // ---> DELETED --->----------+
#ifdef debug
    G4cout<<"G4QFragm::Fragm: nucE="<<nucE<<",nQ="<<nQuas<<G4endl; //                     |
#endif
    if(nucE) nucVel=resNuc4M.vect()/nucE;                   // The NucleusVelocity        |
#ifdef edebug
    G4LorentzVector sq4M=resNuc4M-totLS4M;                  // 4-mom deficit              |
    G4int           sqCg=rnChg-totChg;                      // Charge deficit             |
    G4int           sqBN=rnBaN-totBaN;                      // Baryon number deficit      |
#endif
    for(G4int i=0; i<nQuas; ++i)                            // LOOP over Quasmons         |
    {                                                       //                            |
      G4Quasmon* curQuasm=theQuasmons[i];                   // current Quasmon            |
#ifdef debug
      if(nucE) G4cout<<"G4QFr::Fr:V="<<nucVel<<",Q="<<curQuasm->Get4Momentum()<<",R=" //  |
                     <<resNuc4M<<resNucPDG<<G4endl;         //                            |
#endif
      if(nucE) curQuasm->Boost(-nucVel);                    // Boost it to CMS of Nucleus |
      pan->AddQuasmon(curQuasm);                            // Fill the predefined Quasmon|
#ifdef edebug
      G4LorentzVector cQ4M=curQuasm->Get4Momentum();        // Just for printing          |
      G4cout<<"G4QFragmentation::Fragment: Quasmon# "<<i<<" added, 4M="<<cQ4M<<G4endl; // |
      sq4M+=cQ4M;                                           // Sum up total 4Mom          |
      sqCg+=curQuasm->GetCharge();                          // Sum up the Charge          |
      sqBN+=curQuasm->GetBaryonNumber();                    // Sum up the baryon number   |
#endif
    }                                                       //                            |
#ifdef edebug
    G4cout<<"-EMCLS-G4QFrag::Fragm: r4M="<<sq4M<<", rC="<<sqCg<<", rB="<<sqBN<<G4endl; // |
#endif
    try                                                     //                            |
    {                                                       //                            |
#ifdef debug
    G4cout<<"G4QFrag::Fragm: *** Before Del Output ***"<<G4endl; //                       |
#endif
      delete output;                                        //                            |
#ifdef debug
    G4cout<<"G4QFrag::Fragm: *** After Del Output ***"<<G4endl; //                        |
#endif
      output = pan->Fragment();// DESTROYED after theHadrons are transferred to theResult |
    }                                                       //                          | |
    catch (G4QException& error)                             //                          | |
    {                                                       //                          | |
      G4cerr<<"***G4QFragmentation::Fragment: G4QE Exception is catched"<<G4endl; //    | |
      // G4Exception("G4QFragmentation::Fragment:","27",FatalException,"CHIPSCrash");//    | |
      G4Exception("G4QFragmentation::Fragment()", "HAD_CHPS_0027",
                  FatalException, "CHIPSCrash");
    }                                                       //                          | |
#ifdef debug
    G4cout<<"G4QFrag::Fragm: *** Before Del Pan ***"<<G4endl; //                        | |
#endif
    delete pan;                              // Delete the Nuclear Environment <-----<--+-*
#ifdef debug
    G4cout<<"G4QFrag::Fragm: *** After Del Pan ***"<<G4endl; //                         |
#endif
    if(output)                               // Output exists                           |
    {                                        //                                         |
      G4int nOut=output->size();             // #ofHadrons in the Nuclear Fragmentation |
      for(G4int j=0; j<nOut; j++)            // LOOP over Hadrons transferring to LS    |
      {                                      //                                         |
        G4QHadron* curHadr=(*output)[j];     // Hadron from the nucleus fragmentation   |
        if(nucE) curHadr->Boost(nucVel);     // Boost it back to Laboratory System      |
        G4int hPDG=curHadr->GetPDGCode();    // PDGC of the hadron                      |
        G4LorentzVector h4M=curHadr->Get4Momentum(); // 4-mom of the hadron             |
        if((hPDG!=90000000 && hPDG!=22) || h4M!=nul4M) theResult->push_back(curHadr); //|
        else delete curHadr;                 // delete zero-photons                     |
      }                                      //                                         |
      delete output;                         // Delete the OUTPUT <-----<-----<-----<---+
    }
  }
  else if(!striNum) G4cout<<"-Warning-G4QFragmentation::Fragment:Nothing was done"<<G4endl;
#ifdef debug
  G4cout<<"=--=>G4QFragmentation::Fragment: Final #OfResult="<<theResult->size()<<G4endl;
#endif
  G4int nQ =theQuasmons.size();
  if(nQ) theQuasmons.clear();                              // @@ Not necesary ?
  G4int nHd=theResult->size();
#ifdef edebug
  G4LorentzVector f4M(0.,0.,0.,0.);                        // Sum of the Result in LS
  G4int fCh=totChg;
  G4int fBN=totBaN;
  G4cout<<"-EMCLS-G4QFragmentation::Fragment: #ofHadr="<<nHd<<", #OfQuasm="<<nQ<<G4endl;
  for(G4int i=0; i<nHd; i++)
  {
    G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
    f4M+=hI4M;
    G4int hChg=(*theResult)[i]->GetCharge();
    fCh-=hChg;
    G4int hBaN=(*theResult)[i]->GetBaryonNumber();
    fBN-=hBaN;
    G4cout<<"-EMCLS-G4QFragmentation::Fragment: Hadron#"<<i<<", 4M="<<hI4M<<", PDG="
          <<(*theResult)[i]->GetPDGCode()<<", C="<<hChg<<", B="<<hBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFrag::Fragm: r4M="<<f4M-totLS4M<<", rC="<<fCh<<", rB="<<fBN<<G4endl;
#endif
  //G4QHadron* resNuc = theResult->back();              // Pointer to the Residual Nucleus
  G4QHadron* resNuc = (*theResult)[nHd-1];              // Pointer to the Residual Nucleus
  G4int rnBn = resNuc->GetBaryonNumber();
  G4int rnCg = resNuc->GetCharge();
  if(rnBn==1 && (rnCg==-2 || rnCg==3 || rnCg==-1 || rnCg==2)) // E/Delta decay
  {
    G4LorentzVector tot4M=resNuc->Get4Momentum();       // 4-mom to be split
    G4int nPDG=2212;                                    // Proton as a default
    G4int mPDG=211;                                     // PiPlus as a default
    G4double nM=mProt;                                  // Proton mass as a default
    if(rnCg<0)
    {
      nPDG=2112;
      mPDG=-211;
      nM=mNeut;
    }
    G4LorentzVector m14M(0.,0.,0.,mPiCh);
    G4LorentzVector n4M(0.,0.,0.,nM);
    if(rnCg==-2 || rnCg==3)                             // Decay In 3
    {
      G4LorentzVector m24M(0.,0.,0.,mPiCh);
      if(!G4QHadron(tot4M).DecayIn3(m14M,m24M,n4M))
      {
        G4cerr<<"***G4QFrag::Frag: tM="<<tot4M.m()<<" -> m1="<<mPiCh<<" + m2="<<mPiCh
              <<" + nM="<<nM<<" = "<<2*mPiCh+nM<<G4endl;
        G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ImpossibleDecayIn3");
      }
      theResult->pop_back();
      delete resNuc;
      G4QHadron* m1H = new G4QHadron(mPDG,m14M);
      theResult->push_back(m1H);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn3, M1="<<mPDG<<m14M<<G4endl;
#endif
      G4QHadron* m2H = new G4QHadron(mPDG,m24M);
      theResult->push_back(m2H);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn3, M2="<<mPDG<<m24M<<G4endl;
#endif
      G4QHadron* nH = new G4QHadron(nPDG,n4M);
      theResult->push_back(nH);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn3, Nucleon="<<nPDG<<n4M<<G4endl;
#endif
    }
    else                                            // Decay in 2
    {
      if(!G4QHadron(tot4M).DecayIn2(m14M,n4M))
      {
        G4cerr<<"***G4QFrag::Frag: tM="<<tot4M.m()<<" -> m1="<<mPiCh
              <<" + nM="<<nM<<" = "<<mPiCh+nM<<G4endl;
        G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ImpossibleDecayIn2");
      }
      theResult->pop_back();
      delete resNuc;
      G4QHadron* m1H = new G4QHadron(mPDG,m14M);
      theResult->push_back(m1H);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn2, M1="<<mPDG<<m14M<<G4endl;
#endif
      G4QHadron* nH = new G4QHadron(nPDG,n4M);
      theResult->push_back(nH);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn2, Nucleon="<<nPDG<<n4M<<G4endl;
#endif
    }
  }
  if(rnBn==2)                                       // Di-baryon
  {
    if(!rnCg)                                       // Di-neutron pair
    {
      G4LorentzVector tot4M=resNuc->Get4Momentum(); // 4-mom to be split
      G4LorentzVector n14M(0.,0.,0.,mNeut);
      G4LorentzVector n24M(0.,0.,0.,mNeut);
      if(!G4QHadron(tot4M).DecayIn2(n14M,n24M))
      {
        G4cerr<<"***G4QFrag::Frag: tM="<<tot4M.m()<<" -> n*2="<<2*mNeut<<G4endl;
        G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ImpossibleDecay-2n");
      }
      theResult->pop_back();
      delete resNuc;
      G4QHadron* n1H = new G4QHadron(2112,n14M);
      theResult->push_back(n1H);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn2, Neutron1="<<n14M<<G4endl;
#endif
      G4QHadron* n2H = new G4QHadron(2112,n24M);
      theResult->push_back(n2H);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn2, Neutron2="<<n24M<<G4endl;
#endif
    }
    else if(rnCg==2)                                // Di-proton pair
    {
      G4LorentzVector tot4M=resNuc->Get4Momentum(); // 4-mom to be split
      G4LorentzVector n14M(0.,0.,0.,mProt);
      G4LorentzVector n24M(0.,0.,0.,mProt);
      if(!G4QHadron(tot4M).DecayIn2(n14M,n24M))
      {
        G4cerr<<"***G4QFrag::Frag: tM="<<tot4M.m()<<" -> n*2="<<2*mProt<<G4endl;
        G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ImpossibleDecay-2p");
      }
      theResult->pop_back();
      delete resNuc;
      G4QHadron* n1H = new G4QHadron(2212,n14M);
      theResult->push_back(n1H);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn2, Proton1="<<n14M<<G4endl;
#endif
      G4QHadron* n2H = new G4QHadron(2212,n24M);
      theResult->push_back(n2H);
#ifdef debug
      G4cout<<"G4QFragment::Fragment:DecayIn2, Proton2="<<n24M<<G4endl;
#endif
    }
  } // End of the residual dibaryon decay
  // Now we should check and correct the final state dibaryons (NN-pairs) and 90000000->22
  nHd=theResult->size();
  G4int maxChg=0;                              // max Charge of the hadron found
#ifdef debug
  G4int maxBN=0;                               // max Baryon Number of the hadron found
#endif
  G4QContent maxQC(0,0,0,0,0,0);               // QC for maxChgH for particle UndCoulBar QC
  for(G4int i=0; i<nHd; ++i)
  {
    G4int found=0;
    G4QHadron* cHadr = (*theResult)[i];
    G4int hPDG= cHadr->GetPDGCode();
    if(hPDG==90000000 || hPDG==22)
    {
      G4QHadron* curHadr=(*theResult)[i];
      G4LorentzVector curh4M=curHadr->Get4Momentum();
      if     ( curh4M.e() > 0.) curHadr->SetPDGCode(22);
      else if( curh4M == nul4M) // Kill such a creature
      {
        G4QHadron* theLast = (*theResult)[nHd-1];
        if(theLast != curHadr) // Copy the Last to the current hadron
        {
          curHadr->Set4Momentum(theLast->Get4Momentum()); //4-Mom of CurHadr
          G4QPDGCode lQP=theLast->GetQPDG();
          if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP);
          else curHadr->SetQC(theLast->GetQC());
        }
        theResult->pop_back(); // theLastQHadron is excluded from OUTPUT
        --nHd;
        delete theLast;        //*!!When kill, delete theLastQHadr as an Instance!*
        if(i == nHd-1) break;  // @@ Why it was anyway break ??
      }
    }
    //else if(hPDG==2212 || hPDG==2112) // @@ Why this isotopic correction is necessary ?? 
    else if(2 > 3) // @@ The isotopic exchange (correction) is closed for acceleration
    {
      for(G4int j=i+1; j<nHd; ++j)
      {
        G4int pPDG=(*theResult)[j]->GetPDGCode();
        if(hPDG==pPDG)                       // The pp or nn pair is found
        {
          G4LorentzVector h4M=(*theResult)[i]->Get4Momentum();
          G4LorentzVector p4M=(*theResult)[j]->Get4Momentum();
          G4LorentzVector d4M=h4M+p4M;
          G4double E=d4M.m();                // Proto of tot CM energy (@@ was .mag() ??)
          if(hPDG==2212) E -= mProt+mProt;   // Reduction to tot kin energy in CM
          else           E -= mNeut+mNeut;
          if(E < 140. && G4UniformRand() < .6)// A close pair was found @@ Par 140. @@
          {
            G4int          piPDG= 211;       // Pi+ default for nn pairs
            if(hPDG==2212) piPDG=-211;       // Pi- for pp pairs
            for(G4int k=0; k<nHd; ++k)
            {
              G4int mPDG=(*theResult)[k]->GetPDGCode();
              // @@ if the isotopic exchange is made to increase Pi0, then only piPDG
              // @@ if the isotopic exchange is made to reduce Pi0, then only pi0
              if(mPDG==111 || mPDG==piPDG)   // Appropriate for correction pion is found
              {
                G4LorentzVector m4M=(*theResult)[k]->Get4Momentum();// Must meson be close?
                G4double mN=mProt;           // Final nucleon after charge exchange (nn)
                G4int  nPDG=2212;            // Default for (nn)
                G4int  tPDG=-211;            // Proto Pion after charge exchange from Pi0
                if(hPDG==2212)               // (pp)
                {
                  mN=mNeut;
                  nPDG=2112;
                  tPDG= 211;
                }
                G4double mPi=mPiZr;          // Pion after the charge exchange from Pi+/-
                G4int   sPDG=111;
                if(mPDG==111)
                {
                  mPi=mPiCh;                 // Pion after the charge exchange from Pi0
                  sPDG=tPDG;
                }
                //G4cout<<"G4QFrag::Frag: H="<<hPDG<<", P="<<pPDG<<", M="<<mPDG<<", N="
                //      <<nPDG<<", S="<<sPDG<<G4endl;
                G4double D=mPi+mN;           // The same for both identical nucleons
                G4LorentzVector t4M=m4M+h4M; // meson+ 1st nicleon
                G4LorentzVector n4M=h4M;
                G4double D2=D*D;
                G4double S=t4M.m2();         //  (@@ was .mag2() ??)
                if(S > D2)         found= 1; // 1st nucleon correction can be done
                else                         // Try the 2nd nucleon
                {
                  t4M=m4M+p4M;               // meson+ 1st nicleon
                  n4M=p4M;
                  S=t4M.m2();                //  (@@ was .mag2() ??)
                  if(S > D2)      found=-1;  // 2nd nucleon correction can be done
                }
                if(found)                    // Isotopic Correction
                {
                  G4ThreeVector tV=t4M.vect()/t4M.e();
                  //G4cout<<"G4QFragment::Fragment: Before 4M/M2="<<m4M<<m4M.m2()<<G4endl;
                  m4M.boost(-tV);            // boost the meson to piN CM
                  //G4cout<<"G4QFragment::Fragment: After 4M/M2="<<m4M<<m4M.m2()<<G4endl;
                  n4M.boost(-tV);            // boost the nucleon to piN CM
                  G4double mPi2=mPi*mPi;
                  G4double mN2=mN*mN;
                  G4double C=S-mPi2-mN2;
                  G4double p2=(C*C/4.-mPi2*mN2)/S;
                  if(p2 < 0.) G4cout<<"-Warning-G4QFragment::Fragment: P2="<<p2<<G4endl;
                  G4double pc2=m4M.vect().mag2();
                  //G4double nc2=n4M.vect().mag2();
                  G4double r=1.;
                  if(pc2 < .00000000000001)
                          G4cout<<"-Warning-G4QFragment::Fragment: PC2="<<pc2<<m4M<<G4endl;
                  else r=std::sqrt(p2/pc2);
                  m4M.setV(r*m4M.vect());     // conservs the pion direction (not close!)
                  m4M.setE(std::sqrt(mPi2+p2));
                  //G4cout<<"G4QFragment::Fragment: Changed 4M/M2="<<m4M<<m4M.m2()<<", pc2="
                  //      <<pc2<<", nc2="<<nc2<<G4endl;
                  n4M.setV(r*n4M.vect());
                  n4M.setE(std::sqrt(mN2+p2));
                  m4M.boost(tV);               // Boost the meson back to the Lab system
                  n4M.boost(tV);               // Boost the nucleon back to the Lab system
                  (*theResult)[k]->SetPDGCode(sPDG);
                  (*theResult)[k]->Set4Momentum(m4M);
                  if(found > 0)                // Nucleon correction
                  {
                    (*theResult)[i]->SetPDGCode(nPDG);
                    (*theResult)[i]->Set4Momentum(n4M);
                  }
                  else
                  {
                    (*theResult)[j]->SetPDGCode(nPDG);
                    (*theResult)[j]->Set4Momentum(n4M);
                  }
                  break;                       // Break the pion LOOP
                }
              }
            } // End of the pion LOOP
            if(found) break;                   // Break the nucleon partner LOOP
          } // End of Par 140. IF
        } // End of the identical nucleon IF
      } // End of the nucleon partner LOOP
    } // End of cur=nucleon IF (now closed)
    // Here we can find a hadron with the maximum charge = the residual nuclear environment
    G4int hChg=cHadr->GetCharge();
    if(hChg > maxChg)
    {
      maxChg = hChg;
      maxQC = cHadr->GetQC();
#ifdef debug
      maxBN = cHadr->GetBaryonNumber();
#endif
    }
  } // End of the primary hadron LOOP
  G4QNucleus ResNucEnv(maxQC); // vacuum if not found (check maxChg & maxBN when used)
#ifdef debug
  G4cout<<"G4QFragmentation::Fra: ResNucEnv with A="<<maxBN<<", Z="<<maxChg<<G4endl;
#endif
  // --- The photon && UCB suppressor ---
  G4LorentzVector sum(0.,0.,0.,0.);            // total 4-mom of the found gammas
  G4QContent sumQC(0,0,0,0,0,0);               // aquired positive particle UndCuBar QC
  G4int sumCount=0;                            // Counter of the found gammas
  G4int nHadr=theResult->size();               // #of hadrons in the output so far
  G4bool frag=false;                           // presence of fragments (A>1)
  if(nHadr>2) for(unsigned f=0; f<theResult->size(); f++) //Check that there's a fragment
  {
    G4int fBN=(*theResult)[f]->GetBaryonNumber(); // Baryon number of the fragment
#ifdef debug
    G4int fPDG=(*theResult)[f]->GetPDGCode();  // PDG code of the possible fragment
    G4LorentzVector fLV=(*theResult)[f]->Get4Momentum(); // 4Mom of the possible fragment
    G4cout<<"G4QFragmentation::Fra:"<<f<<",PDG="<<fPDG<<",fBN="<<fBN<<",f4M="<<fLV<<G4endl;
#endif
    if(fBN>1)                                  // At least one fragment (A>1) is found
    {
      frag=true;
      break;
    }
  }
#ifdef debug
  G4cout<<"G4QFrag::Frag:=>Before Gamma Suppression<=, nH="<<nHadr<<",frag="<<frag<<G4endl;
#endif
  if(nHadr>2 && frag) for(G4int h=nHadr-1; h>=0; h--)//Collect gammas & kill DecayedHadrons
  {
    G4QHadron* curHadr = (*theResult)[h];      // Get a pointer to the current Hadron
    G4int   hF = curHadr->GetNFragments();     // This is historic ... (was decayed flag)
    G4int hPDG = curHadr->GetPDGCode();        // PDG of the hadron
    G4LorentzVector h4M=curHadr->Get4Momentum();// 4Mom of the hadron
    if(hPDG==89999003||hPDG==90002999)G4cout<<"-Warning-G4QFr::Fr:nD-/pD++="<<hPDG<<G4endl;
#ifdef debug
    G4cout<<"G4QFragmentation::Fragm: h#"<<h<<", hPDG="<<hPDG<<", hNFrag="<<hF<<G4endl;
#endif
    G4int hChg = curHadr->GetCharge();         // Charge of the hadron
    G4bool UCB = false;                        // Not UCB yet
    if(hChg > 0 && hPDG!=321)                  // All positive but not K+
    {
      G4int hBN = curHadr->GetBaryonNumber();  // Baryon Number of the hadron
      G4double hCB=ResNucEnv.CoulombBarrier(hChg,hBN); // Coulomb barrier
      if(h4M.e()-h4M.m() < hCB) UCB = true;    // The hadron should be absorbed
    }
    if(hF || hPDG==22 || UCB)                  // Must be absorbed (decayed, photon, UCB)
    {
      G4int last=theResult->size()-1;
      G4QHadron* theLast = (*theResult)[last]; //Get Ptr to the Last Hadron
      if(hPDG==22 || UCB)                      // Absorb if this is gamma
      {
        sum+=h4M;                              // Add 4Mom of hadron to the "sum"
        sumCount++;
        if(UCB) sumQC+=curHadr->GetQC();       // Collect the absorbed QC
#ifdef debug
        G4cout<<"G4QFragmentation::Frag: gam4M="<<h4M<<" is added to s4M="<<sum<<G4endl;
#endif
      }
      nHadr = static_cast<G4int>(theResult->size())-1;
      if(h < last)                             // Need swap theCurHadron with the Last
      {
        curHadr->SetNFragments(0);
        curHadr->Set4Momentum(theLast->Get4Momentum());
        G4QPDGCode lQP=theLast->GetQPDG();     // The QPDG of the last
        if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of LastHadr
        else curHadr->SetQC(theLast->GetQC()); // CurHadrPDG instead of LastHadrPDG
#ifdef debug
        G4cout<<"G4QFragmentation::Fragment: Exchange with the last is done"<<G4endl;
#endif
      }
      theResult->pop_back();                   // theLastQHadron is excluded from theResult
      delete theLast;//!!When kill,DON'T forget to delete theLastQHadron as an instance!!
#ifdef debug
      G4cout<<"G4QFragmentation::Fragment: The last is compessed"<<G4endl;
#endif
    }
  }
#ifdef debug
  G4cout<<"G4QFragment::Frag: nH="<<nHadr<<"="<<theResult->size()<<", sum="<<sum<<G4endl;
#endif
  if(nHadr > 1) for(unsigned hdr=0; hdr<theResult->size()-1; hdr++)//Ord:theBigestIsTheLast
  {
    G4QHadron* curHadr = (*theResult)[hdr];    // Get a pointer to the current Hadron
#ifdef debug
    G4cout<<"G4QFrag::Frag: h#"<<hdr<<"<"<<nHadr<<", hPDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
    G4QHadron* theLast = (*theResult)[theResult->size()-1]; //Get Ptr to the Last Hadron
    G4int hB           = curHadr->GetBaryonNumber();
    G4int lB           = theLast->GetBaryonNumber();
#ifdef debug
    G4cout<<"G4QFra::Fra:hBN="<<hB<<"<lBN="<<lB<<",lstPDG="<<theLast->GetPDGCode()<<G4endl;
#endif
    if(lB < hB)                                // Must be swapped
    {
      G4QPDGCode   hQPDG = curHadr->GetQPDG();
      G4LorentzVector h4m= curHadr->Get4Momentum();
      curHadr->Set4Momentum(theLast->Get4Momentum());
      G4QPDGCode lQP=theLast->GetQPDG();       // The QPDG of the last
      if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of LastHadr
      else curHadr->SetQC(theLast->GetQC());   // CurHadrPDG instead of LastHadrPDG
      theLast->Set4Momentum(h4m);
      theLast->SetQPDG(hQPDG);
    }
  }
  nHadr=theResult->size(); // --> At this point the last hadron is the biggest nucleus
  if(sumCount)
  {
    G4QHadron* theLast = (*theResult)[nHadr-1];// Get a pointer to the Last Hadron
    G4int nucEnvBN=theLast->GetBaryonNumber(); // Initial A of residualNuclearEnvironment
    if ( nucEnvBN > 0 )                       // "Absorb phot/UCB & evaporate/decay" case
    {
      G4QHadron* theNew  = new G4QHadron(theLast); // Make New Hadron of the Last Hadron
#ifdef debug
      G4cout<<"G4QFra::Fra:BeforeLastSub,n="<<nHadr<<",PDG="<<theNew->GetPDGCode()<<G4endl;
#endif
      theResult->pop_back();                   // the last QHadron is excluded from OUTPUT
      delete theLast;//*!When kill,DON'T forget to delete theLastQHadron as an instance!*
      nHadr--;                                 // TheLastHadron only virtually exists now
      G4QContent newQC=theNew->GetQC();        // QContent of the fragment=ResNuclEnv
      G4LorentzVector new4M=theNew->Get4Momentum(); // 4-mom of the fragment=ResNuclEnv
#ifdef debug
      G4cout<<"G4QFra::Fra:gSum4M="<<sum<<" is added to "<<new4M<<", QC="<<newQC<<G4endl;
#endif
      G4LorentzVector exRes4M = new4M + sum;   //Icrease 4Mom of theLast by sum 4Mon
      G4QContent exResQC = newQC + sumQC;      //Icrease QCont of theLast by sumQC
      theNew->Set4Momentum(exRes4M);
      theNew->SetQC(exResQC);
#ifdef debug
      G4cout<<"G4QFra::Fra:BeforeEvap, n="<<nHadr<<", nPDG="<<theNew->GetPDGCode()<<G4endl;
#endif
      EvaporateResidual(theNew); // Try to evaporate theNucl. (del. equiv.)
      nHadr=theResult->size();
    } // End of "the last is the nucleus" case
    else G4cout<<"-Warning-G4QFragmentation::Fragment:RA="<<nucEnvBN<<",E/M cons?"<<G4endl;
  } // End of "There are gammas to suppress"
  // End of the Gamma Suppression
  return theResult;
} // End of fragmentation

void G4QFragmentation::Breeder()
{ // This is the member function, which returns the resulting vector of Hadrons & Quasmons
  static const G4double  eps = 0.001;                              // Tolerance in MeV
  //static const G4QContent vacQC(0,0,0,0,0,0);
  static const G4LorentzVector vac4M(0.,0.,0.,0.);
  //
  // ------------ At this point the strings are fragmenting to hadrons in LS -------------
  //
#ifdef edebug
  G4LorentzVector totLS4M=theNucleus.Get4Momentum();    // Nucleus 4Mom in LS
  G4int totChg=theNucleus.GetZ();
  G4int totBaN=theNucleus.GetA();
  G4int nStri=strings.size();
  G4cout<<"-EMCLS-G4QFr::Breed: CHECKRecovery #ofS="<<nStri<<",N4M(E=M)="<<totLS4M<<G4endl;
  for(G4int i=0; i<nStri; i++)
  {
    G4LorentzVector strI4M=strings[i]->Get4Momentum();
    totLS4M+=strI4M;
    G4int sChg=strings[i]->GetCharge();
    totChg+=sChg;
    G4int sBaN=strings[i]->GetBaryonNumber();
    totBaN+=sBaN;
    G4cout<<"-EMCLS-G4QFragm::Breeder: St#"<<i<<", 4M="<<strI4M<<", M="<<strI4M.m()
          <<", C="<<sChg<<", B="<<sBaN<<G4endl;
  }
#endif
  G4int nOfStr=strings.size();
#ifdef debug
  G4cout<<"G4QFragmentation::Breeder: BeforeFragmentation, #OfStr="<<nOfStr<<G4endl;
#endif
  G4LorentzVector ft4M(0.,0.,0.,0.);
  G4QContent      ftQC(0,0,0,0,0,0);
  G4bool          ftBad=false;
  for(G4int i=0; i < nOfStr; ++i)
  {
    G4QString* crStr=strings[i];
    G4LorentzVector pS4M=crStr->Get4Momentum();     // String 4-momentum
    ft4M+=pS4M;
    G4QContent pSQC=crStr->GetQC();                 // String Quark Content
    ftQC+=pSQC;
    if(pS4M.m2() < 0.) ftBad=true;
#ifdef debug
    G4cout<<">G4QFrag::Br:1stTest,S#"<<i<<",P="<<crStr<<",4M="<<pS4M<<",QC="<<pSQC<<G4endl;
#endif
  }
  if(ftBad)
  {
    G4Quasmon* stringQuasmon = new G4Quasmon(ftQC, ft4M);
#ifdef debug
    G4cout<<"->G4QFragmentation::Breeder:*TotQ*,QC="<<ftQC<<",4M="<<ft4M<<ft4M.m()<<G4endl;
#endif
    theQuasmons.push_back(stringQuasmon);
    G4LorentzVector r4M=theNucleus.Get4Momentum();  // Nucleus 4-momentum in LS
    G4int rPDG=theNucleus.GetPDG();
    G4QHadron* resNuc = new G4QHadron(rPDG,r4M);
    theResult->push_back(resNuc);                   // Fill the residual nucleus
    return;
  }
  for (G4int astring=0; astring < nOfStr; astring++)
  {
#ifdef edebug
    G4LorentzVector sum=theNucleus.Get4Momentum();  // Nucleus 4Mom in LS
    G4int rChg=totChg-theNucleus.GetZ();
    G4int rBaN=totBaN-theNucleus.GetA();
    G4int nOfHadr=theResult->size();
    G4cout<<"-EMCLS-G4QFragmentation::Breeder:#ofSt="<<nOfStr<<",#ofHad="<<nOfHadr<<G4endl;
    for(G4int i=astring; i<nOfStr; i++)
    {
      G4LorentzVector strI4M=strings[i]->Get4Momentum();
      sum+=strI4M;
      G4int sChg=strings[i]->GetCharge();
      rChg-=sChg;
      G4int sBaN=strings[i]->GetBaryonNumber();
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
    G4QString* curString=strings[astring];
    if(!curString->GetDirection()) continue;  // Historic for the dead strings: DoesNotWork
    G4int curStrChg = curString->GetCharge();
    G4int curStrBaN = curString->GetBaryonNumber();
    G4LorentzVector curString4M = curString->Get4Momentum();
#ifdef debug
    G4cout<<"=--=>G4QFragmentation::Breeder: String#"<<astring<<",s4M/m="<<curString4M
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
        if(next < strings.size())              // TheString isn't theLastString can fuse
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
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-LL-");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair((-nPDG)/100, mPDG/100);
                  G4int newCR=resRR.first;
                  G4int newPR=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QFragmentation::Breeder:CR="<<newCR<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-RR-");
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
                    G4cerr<<"*G4QFragmentation::Breeder:CL="<<newCL<<",PL="<<newPL<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-LL+");
                  }
                  std::pair<G4int,G4int> resRR=ReducePair(nPDG/100, (-mPDG)/100);
                  G4int newCR=resRR.first;
                  G4int newPR=resRR.second;
                  if(!newCR || !newPR)
                  {
                    G4cerr<<"*G4QFragmentation::Breeder:CR="<<newCR<<",PR="<<newPR<<G4endl;
                    G4Exception("G4QFragmentation::Breeder:","72",FatalException,"2-RR+");
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
#ifdef debug
          else
          {

            G4cout<<"**G4QFragmentation::Breeder:*NoPart*M="<<curString->Get4Momentum().m()
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
            G4QString* cString=strings[astring];  // Must be the last string by definition
            G4LorentzVector cString4M = cString->Get4Momentum();
            cLeft=cString->GetLeftParton();
            cRight=cString->GetRightParton();
            G4int sumT=cLeft->GetType()+cRight->GetType();
            sPDG=cLeft->GetPDGCode();
            nPDG=cRight->GetPDGCode();
            G4int cMax=0;                         // Both or only one cuark can merge
            for (G4int reh=0; reh < nHadr; reh++)
            {
              G4QHadron* curHadr=(*theResult)[reh];
              G4int curPDG=curHadr->GetPDGCode(); // PDGCode of the hadron
              G4QContent curQC=curHadr->GetQC();  // Quark content of the hadron
              if(curPDG==331 && sPDG!=3 && nPDG!=3 && sPDG!=-3 && nPDG!=-3) // eta' red
              {
                if(sPDG==2 || sPDG==-2 || nPDG==2 || nPDG==-2)
                                                             curQC=G4QContent(0,1,0,0,1,0);
                else                                         curQC=G4QContent(1,0,0,1,0,0);
              }
              else if(curPDG==221 && sPDG!=2 && nPDG!=2 && sPDG!=-2 && nPDG!=-2) // eta
                                                             curQC=G4QContent(1,0,0,1,0,0);
              else if(curPDG==111 && sPDG!=1 && nPDG!=1 && sPDG!=-1 && nPDG!=-1) // eta
                                                             curQC=G4QContent(0,1,0,0,1,0);
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
#ifdef debug
              G4cout<<"**G4QFragmentation::Breeder:*NoH*,M="<<curString->Get4Momentum().m()
                    <<", LPDG="<<curString->GetLeftParton()->GetPDGCode()
                    <<", RPDG="<<curString->GetRightParton()->GetPDGCode()<<G4endl;
              // @@ Temporary exception for the future solution
              //G4Exception("G4QFragmentation::Breeder:","72",FatalException,"SHNotFused");
#endif
              break;                           // Breake the While LOOP
            } // End of the namespace where both Fusion and reduction have failed
            // The fused hadron must be excluded from theResult
#ifdef debug
            G4cout<<"G4QFragmentation::Breeder: before HR, nH="<<theResult->size()<<G4endl;
            G4int icon=0;                              // Loop counter
#endif
            G4QHadronVector::iterator ih;
#ifdef debug
            G4bool found=false;
#endif
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
                G4LorentzVector p4M=selHP->Get4Momentum();
                //
                curString4M+=p4M;
                G4int Chg=selHP->GetCharge();
                G4int BaN=selHP->GetBaryonNumber();
                curStrChg+=Chg;
                curStrBaN+=BaN;
#ifdef edebug
                G4cout<<"-EMC->>G4QFragmentation::Breeder: S+=H, 4M="<<curString4M<<", M="
                      <<curString4M.m()<<", Charge="<<curStrChg<<", B="<<curStrBaN<<G4endl;
#endif
                delete selHP;                          // delete the Hadron
                theResult->erase(ih);                  // erase the Hadron from theResult
#ifdef debug
                found=true;
#endif
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
                G4cout<<"-EMC->>G4QFragment::Breeder: String=Hadr ChiPro1 is filled, f4M="
                      <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
                G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
                theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
                G4LorentzVector s4M=h2H->Get4Momentum();
                G4int           sPD=h2H->GetPDGCode();
                G4int           sCg=h2H->GetCharge();
                G4int           sBN=h2H->GetBaryonNumber();
                G4cout<<"-EMC->>G4QFragmentation::Breeder: String=Hadr ChiPro2 is filled, s4M="
                      <<s4M<<", sPDG="<<sPD<<", sCg="<<sCg<<", sBN="<<sBN<<G4endl;
#endif
#ifdef edebug
                G4cout<<"-EMC-..Chi..G4QFragmentation::Breeder: DecayCHECK, Ch4M="
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
              G4cout<<"G4QFragmentation::Breeder: Decay the Last, Res#H="<<tmpN<<G4endl;
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
                  G4cout<<"-EMC->>G4QFragment::Breeder:String=Hadr,H#"<<aH<<" filled, 4M="
                        <<p4M<<", PDG="<<PDG<<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
                }
              }
              else
              {
                G4Quasmon* stringQuasmon = new G4Quasmon(miQC, curString4M);// String->Quas
#ifdef debug
                G4cout<<"G4QFragmentat::Breeder:==> to Quasm="<<miQC<<curString4M<<", Nuc="
                      <<theNucleus<<theNucleus.Get4Momentum()<<", NString="<<strings.size()
                      <<", nR="<<theResult->size()<<", nQ="<<theQuasmons.size()<<G4endl;
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
      G4cout<<"G4QFragmentation::Breeder: theH="<<theHadrons<<"?=0, next="<<next<<G4endl;
#endif
      if(!theHadrons && next < strings.size())       // ForwardInLOOP strings exist
      {
        // @@ string can be not convertable to one hadron (2,0.0,0,2,0) - To be improved
        G4QContent miQC=curString->GetQC(); // QContent of the Lightest Hadron
        G4int miPDG=miQC.GetSPDGCode();// PDG of the Lightest Hadron
#ifdef debug
        G4cout<<"---->>G4QFragmentation::Breeder: SQC="<<miQC<<", miSPDG="<<miPDG<<G4endl;
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
        G4cout<<"---->>G4QFragmentation::Breeder: minMass="<<miM<<", realM2="<<cM2<<G4endl;
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
              G4cout<<"----->>G4QFragmentation::Breeder:S->H="<<miPDG<<curString4M<<G4endl;
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
              G4cout<<"-EMC->>G4QFragmentation::Breeder:Str=2HadrAR Prod-F is filled, f4M="
                    <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
              G4LorentzVector s4M=h2H->Get4Momentum();
              G4int           sPD=h2H->GetPDGCode();
              G4int           sCg=h2H->GetCharge();
              G4int           sBN=h2H->GetBaryonNumber();
              G4cout<<"-EMC->>G4QFragmentation::Breeder:Str=2HadrAR Prod-S is filled, s4M="
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
            G4int restr=0;               // To use beyon the LOOP for printing
            G4int fustr=0;               // Selected String-partner (0 = NotFound)
            G4double selX=0.;            // Selected value of x
            G4double maD=-DBL_MAX;       // Maximum Free Mass
            G4double Vmin=DBL_MAX;       // min Velocity Distance
            G4LorentzVector s4M(0.,0.,0.,0.); // Selected 4-mom of the hadron
#ifdef debug
            G4cout<<"G4QFr::Breed:TryRecover,V="<<curV<<",cM2="<<cM2<<",miM="<<miM<<G4endl;
#endif
            nOfStr=strings.size();
            for(restr=next; restr < nOfStr; ++restr) if(restr != astring)
            {
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder: rS="<<restr<<", nS="<<nOfStr<<G4endl;
#endif
              G4QString* pString=strings[restr];
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder: pString="<<pString<<G4endl;
#endif
              G4LorentzVector p4M=pString->Get4Momentum();
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder: p4M="<<p4M<<G4endl;
#endif
              G4ThreeVector pP=p4M.vect();  // Momentum of the partnerString
              G4double      pE=p4M.e();     // Energy of the partnerString
              G4double D2=cE*pE-cP.dot(pP); 
              G4double pM2=p4M.m2();
#ifdef debug
              G4cout<<"G4QFrag::Breeder: pM2="<<pM2<<",miM2="<<miM2<<",cM2="<<cM2<<G4endl;
#endif
              G4double dM4=pM2*(miM2-cM2);
              G4double D=D2*D2+dM4;
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder: D="<<D<<",dM4="<<dM4<<",D2="<<D2<<G4endl;
#endif
              G4double x=-1.;               // Bad preexpectation 
              if(D > 0. && pM2>.01) x=(std::sqrt(D)-D2)/pM2; // what we should get from p
#ifdef debug
              else G4cout<<"G4QFragment::Breeder: D="<<D<<",D2="<<D2<<",pM4="<<dM4<<G4endl;
              G4cout<<"G4QFragmentation::Breeder: pM2="<<pM2<<",D2="<<D2<<",x="<<x<<G4endl;
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
                  G4cout<<"G4QFragmentation::Breeder: Subtr,S#"<<restr<<",d="<<maD<<G4endl;
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
                  G4cout<<"G4QFragmentation::Breeder: FreeAdd,S#"<<restr<<",x="<<x<<G4endl;
#endif
                  Vmin=dV;
                  fustr=restr;
                  selX=x;
                  s4M=p4M;
                }
              }
#ifdef debug
              G4cout<<"G4QFragmentation::Breeder:EndOfLOOP r="<<restr<<"<"<<nOfStr<<G4endl;
#endif
            } // End of the LOOP over string-partners for Correction
#ifdef debug
            G4cout<<"G4QFragmentation::Breeder: AfterLOOP fustr="<<fustr<<G4endl;
#endif
            if(fustr)
            {
#ifdef edebug
              G4LorentzVector sum4M=s4M+curString4M;
              G4cout<<"G4QFragmentation::Breeder: Found Sum4M="<<sum4M<<G4endl;
#endif
              G4QString* pString=strings[fustr];
              curString4M+=selX*s4M;
              if(std::abs(miPDG)%10 > 2)                  // Decay String-Hadron-Resonance
              {
                G4Quasmon Quasm;
                G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
                G4QHadronVector* tmpQHadVec=Quasm.DecayQHadron(sHad); // It deletes sHad
#ifdef debug
                G4cout<<"G4QFragmentation::Breeder:DecStH,nH="<<tmpQHadVec->size()<<G4endl;
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
                  G4cout<<"-EMC->>G4QFragmentation::Breeder:St=Had,pH#"<<aH<<" filled, 4M="
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
                      G4cerr<<"***G4QFragmentation::Breeder: tM="<<ttM<<"->h1="<<h1QPDG
                            <<"("<<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")="<<h1M+h2M<<G4endl;
                      G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ChDe");
                    }
                  }
                  G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
                  theResult->push_back(h1H);        // (delete equivalent)  
#ifdef debug
                  G4LorentzVector f4M=h1H->Get4Momentum();
                  G4int           fPD=h1H->GetPDGCode();
                  G4int           fCg=h1H->GetCharge();
                  G4int           fBN=h1H->GetBaryonNumber();
                  G4cout<<"-EMC->>G4QFragmentation::Breeder:Str=Hadr Prod-F's filled, f4M="
                        <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
                  G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
                  theResult->push_back(h2H);        // (delete equivalent)  
#ifdef debug
                  G4LorentzVector s4M=h2H->Get4Momentum();
                  G4int           sPD=h2H->GetPDGCode();
                  G4int           sCg=h2H->GetCharge();
                  G4int           sBN=h2H->GetBaryonNumber();
                  G4cout<<"-EMC->>G4QFragmentation::Breeder:Str=Hadr Prod-S's filled, s4M="
                        <<s4M<<", sPDG="<<sPD<<", sCg="<<sCg<<", sBN="<<sBN<<G4endl;
#endif
#ifdef edebug
                  G4cout<<"-EMC-Chipo.G4QFragmentation::Breeder:DecCHECK,c4M="<<curString4M
                        <<", ChQC="<<miQC<<", d4M="<<curString4M-h14M-h24M<<G4endl;
#endif
                }
                else
                {
                  G4cerr<<"***G4QFragm::Breeder:tM="<<ttM<<miQC<<"->h1="<<h1QPDG<<"(" <<h1M
                        <<")+h2="<<h1QPDG<<"("<<h2M<<") = "<<h1M+h2M<<G4endl;
                  G4Exception("G4QFragmentation::Breeder:","72",FatalException,"ChiDecMa");
                }
              }
              else
              {
                G4QHadron* sHad = new G4QHadron(miPDG,curString4M);
                theResult->push_back(sHad);         // The original string-hadron is filled
#ifdef debug
                G4cout<<"-EMC->>G4QFragmentation::Breeder:Str=Hadr Filled, 4M="
                      <<curString4M<<", PDG="<<miPDG<<G4endl;
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
              G4cout<<"---->>G4QFragmentation::Breeder:*Corrected* String->Hadr="<<miPDG
                    <<curString4M<<" by String #"<<fustr<<G4endl;
#endif
              continue;                            // Continue the LOOP over the curStrings
            } // End of Found combination for the String on string correction
          } // End of the Try-to-recover String+String correction algorithm
        } // End of IF(CM2>0.)
      } // End of IF(Can try to correct by String-String)
#ifdef debug
      else G4cout<<"***G4QFragmentation::Breeder: **No SSCorrection**,next="<<next<<G4endl;
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
      G4cout<<"G4QFragmentation::Breeder: theH="<<theHadrons<<", #OfH="<<nofRH<<G4endl;
#endif
      if(!theHadrons && nofRH)                  // Hadrons are existing for SH Correction
      {
#ifdef debug
        G4cout<<"!G4QFragmentation::Breeder:Can try SHCor, nH="<<theResult->size()<<G4endl;
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
        G4int fuha=-1;                          // Selected Hadron-partner (0 = NotFound)
        G4double dMmin=DBL_MAX;                 // min Excess of the mass
        G4LorentzVector s4M(0.,0.,0.,0.);       // Selected 4-mom of the Hadron+String
        G4double sM=0.;                         // Selected Mass of the Hadron+String
        for (reha=0; reha < nofRH; reha++)      // LOOP over already collected hadrons
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
#ifdef debug
          else G4cout<<"G4QFragmentation::Breeder:H# "<<reha<<",tM="<<std::sqrt(tM2)<<" < "
                     <<" mS="<<miM<<" + mH="<<pM<<" = "<<pM+miM<<G4endl;
#endif
        } // End of the LOOP over string-partners for Correction
#ifdef debug
        G4cout<<"G4QFragment::Breeder: fuha="<<fuha<<", dMmin="<<dMmin<<G4endl;
#endif
        if(fuha>-1)                             // The hadron-partner was found
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
                G4cerr<<"***G4QFragmentation::Breeder: *SH*, tM="<<sM<<"->h1=("<<miPDG<<")"
                      <<miM<<" + h2="<<spM<<" = "<<miM+spM<<G4endl;
                G4Exception("G4QFragmentation::Breeder:","72",FatalException,"SHChipoDec");
              }
            }
            pHadron->Set4Momentum(sp4M);
#ifdef debug
            G4cout<<"-EMC->...>G4QFragmentation::Breeder: H# "<<fuha<<" is updated, new4M="
                  <<sp4M<<G4endl;
#endif
          }
          else
          {
            G4cerr<<"***G4QFragm::Breeder: HS Failed, tM="<<sM<<"->h1M=("<<miPDG<<")"<<miM
                  <<"+h2M="<<spM<<" = "<<miM+spM<<G4endl;
            G4Exception("G4QFragmentation::Breeder:","72",FatalException,"HSChipoDecMass");
          }
          if(std::abs(miPDG)%10 > 2)                  // Decay Hadron-Resonans
          {
            G4Quasmon Quasm;
            G4QHadron* sHad = new G4QHadron(miPDG,mi4M);
            G4QHadronVector* tmpQHadVec=Quasm.DecayQHadron(sHad); // It deletes sHad
#ifdef debug
            G4cout<<"G4QFragment::Breeder:*HS* DecStrHad, nH="<<tmpQHadVec->size()<<G4endl;
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
              G4cout<<"-EMC->>G4QFragmentation::Breeder: Str+Hadr PrH#"<<aH<<" filled, 4M="
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
                  G4cerr<<"***G4QFragmentation::Breeder: HS tM="<<ttM<<"->h1="<<h1QPDG<<"("
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
              G4cout<<"-EMC->>G4QFragmentation::Breeder: CorStrHadr Prod-1 is filled, f4M="
                    <<f4M<<", fPDG="<<fPD<<", fCg="<<fCg<<", fBN="<<fBN<<G4endl;
#endif
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theResult->push_back(h2H);         // (delete equivalent)  
#ifdef debug
              G4LorentzVector n4M=h2H->Get4Momentum();
              G4int           nPD=h2H->GetPDGCode();
              G4int           nCg=h2H->GetCharge();
              G4int           nBN=h2H->GetBaryonNumber();
              G4cout<<"-EMC->>G4QFragmentation::Breeder: CorStrHadr Prod-2 is filled, n4M="
                    <<n4M<<", nPDG="<<nPD<<", nCg="<<nCg<<", nBN="<<nBN<<G4endl;
#endif
#ifdef edebug
              G4cout<<"-EMC-...HS-Chipo...G4QFragmentation::Breeder:DecCHECK, Ch4M="<<mi4M
                    <<", ChQC="<<miQC<<", d4M="<<mi4M-h14M-h24M<<G4endl;
#endif
            }
          }
          else
          {
            G4QHadron* sHad = new G4QHadron(miPDG, mi4M);
            theResult->push_back(sHad);          // The original string=hadron is filled
#ifdef debug
            G4cout<<"----->>G4QFragmentation::Breeder: CorStr=Hadr is Filled, 4M="
                  <<curString4M<<", StPDG="<<miPDG<<G4endl;
#endif
          }
#ifdef edebug
          G4cout<<"-EMC-...Cor...G4QFragmentation::Breeder:StHadCor CHECK Sum="<<s4M
                <<" =? "<<mi4M+pHadron->Get4Momentum()<<G4endl;
#endif
#ifdef debug
          G4cout<<"------->>G4QFragmentation::Breeder: *Corrected* String+Hadr="<<miPDG
                <<mi4M<<" by Hadron #"<<reha<<G4endl;
#endif
          continue;                    // Continue the LOOP over the curStrings
        }
        else
        {
#ifdef debug
          G4cout<<"G4QFragmentation::Breeder: Str+Hadr Failed, 4M="<<curString4M
                <<", PDG="<<miPDG<<" -> Now try to recover the string as a hadron"<<G4endl;
#endif
          //for (reha=0; reha < nofRH; reha++)      // LOOP over already collected hadrons
          //{
          //  G4QHadron* pHadron=(*theResult)[reha];// Pointer to the CurrentPartnerHadron
          //  G4LorentzVector p4M=pHadron->Get4Momentum();
          //  G4double         pM=p4M.m();          // Mass of the Partner-Hadron
          //  G4LorentzVector t4M=p4M+curString4M;  // Total momentum of the compound
          //  G4double        tM2=t4M.m2();         // Squared total mass of the compound
          //  if(tM2 >= sqr(pM+miM+eps))            // Condition of possible correction
          //  {
          //    G4double tM=std::sqrt(tM2);         // Mass of the Hadron+String compound
          //    G4double dM=tM-pM-miM;              // Excess of the compound mass
          //    if(dM < dMmin)
          //    {
          //      dMmin=dM;
          //      fuha=reha;
          //      spM=pM;
          //      s4M=t4M;
          //      sM=tM;
          //    }
          //  }
          //} // End of the LOOP over string-partners for Correction
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
        G4cout<<"=--=>G4QFragmentation::Breeder:S#"<<i<<",4M="<<pS4M<<",QC="<<pSQC<<G4endl;
#endif
      }
#ifdef debug
      G4cout<<"==>G4QFragmentation::Breeder:AllStrings are summed up in a Quasmon"<<G4endl;
#endif
      G4Quasmon* stringQuasmon = new G4Quasmon(ssQC, ss4M);
      theQuasmons.push_back(stringQuasmon);
      break;                                   // break the LOOP over Strings
    }
#ifdef debug
    G4cout<<"G4QFragmentation::Breeder: Trying to decay hadrons #ofHRes="<<nHfin<<G4endl;
#endif
    for(G4int aTrack=0; aTrack<nHfin; aTrack++)
    {
      G4QHadron* curHadron=(*theHadrons)[aTrack];
      G4int hPDG=curHadron->GetPDGCode();
      G4LorentzVector curH4M=curHadron->Get4Momentum();
      G4int           curHCh=curHadron->GetCharge();
      G4int           curHBN=curHadron->GetBaryonNumber();
#ifdef debug
      G4cout<<"----->>G4QFragmentation::Breeder:S#"<<astring<<",H#"<<aTrack<<",PDG="<<hPDG
            <<",4M="<<curHadron->Get4Momentum()<<G4endl;
#endif
      if(std::abs(hPDG)%10 > 2)
      {
        G4QHadronVector* tmpQHadVec=tmpQ.DecayQHadron(curHadron); // It deletes curHadron
#ifdef debug
        G4cout<<"G4QFragmentation::Breeder:-DECAY'S DONE-,nH="<<tmpQHadVec->size()<<G4endl;
#endif
        for(unsigned aH=0; aH < tmpQHadVec->size(); aH++)
        {
          theResult->push_back((*tmpQHadVec)[aH]);// TheDecayProduct of TheHadron is filled
          //
          G4QHadron*   prodH =(*tmpQHadVec)[aH];
          G4LorentzVector p4M=prodH->Get4Momentum();
          G4int           Chg=prodH->GetCharge();
          G4int           BaN=prodH->GetBaryonNumber();
          curString4M-=p4M;
          curStrChg-=Chg;
          curStrBaN-=BaN;
          curH4M-=p4M;
          curHCh-=Chg;
          curHBN-=BaN;
#ifdef edebug
          G4int           PDG=prodH->GetPDGCode();
          G4cout<<"-EMC->>G4QFragmentation::Breeder:String*Filled, 4M="<<p4M<<", PDG="<<PDG
                <<", Chg="<<Chg<<", BaN="<<BaN<<G4endl;
#endif
        }
#ifdef edebug
        G4cout<<"-EMC-.G4QFr::Br:Dec,r4M="<<curH4M<<",rC="<<curHCh<<",rB="<<curHBN<<G4endl;
#endif
        tmpQHadVec->clear();
        delete tmpQHadVec;  // Who calls DecayQHadron is responsible for clear & delete
      }
      else                                      // Chipolino is not checked here
      {
        theResult->push_back(curHadron);        // The original hadron is filled
        //
        curString4M-=curH4M;
        G4int curCh=curHadron->GetCharge();
        G4int curBN=curHadron->GetBaryonNumber();
        curStrChg-=curCh;
        curStrBaN-=curBN;
#ifdef edebug
        G4cout<<"-EMC->>-->>G4QFragmentation::Breeder: curH filled 4M="<<curH4M<<",PDG="
              <<curHadron->GetPDGCode()<<", Chg="<<curCh<<", BaN="<<curBN<<G4endl;
#endif
      }
    }
    // clean up (the issues are filled to theResult)
    if(theHadrons) delete theHadrons;
#ifdef edebug
    G4cout<<"-EMC-.........G4QFragmentation::Breeder: StringDecay CHECK, r4M="<<curString4M
          <<", rChg="<<curStrChg<<", rBaN="<<curStrBaN<<G4endl;
#endif
    // Trap with the debugging warning --- Starts ---
    if(curStrChg || curStrBaN || curString4M.t() > eps || std::fabs(curString4M.x()) > eps
       || std::fabs(curString4M.y()) > eps || std::fabs(curString4M.z()) > eps )
    {
      G4double dEn=curString4M.t();
      G4double dPx=curString4M.x();
      G4double dPy=curString4M.y();
      G4double dPz=curString4M.z();
      G4int nHadr=theResult->size();
      G4double hEn=0.;
      G4double hPx=0.;
      G4double hPy=0.;
      G4double hPz=0.;
      G4int    hCh=0;
      G4int    hBN=0;
      G4double mEn=0.;
      G4double mPx=0.;
      G4double mPy=0.;
      G4double mPz=0.;
      G4int    mCh=0;
      G4int    mBN=0;
      for(G4int i=0; i<nHadr; i++)
      {
        mEn=hEn; // Previous hadron
        mPx=hPx;
        mPy=hPy;
        mPz=hPz;
        mCh=hCh;
        mBN=hBN;
        G4QHadron* curHadr = (*theResult)[i];
        G4LorentzVector hI4M = curHadr->Get4Momentum();
        hEn=hI4M.t();
        hPx=hI4M.x();
        hPy=hI4M.y();
        hPz=hI4M.z();
        hCh=curHadr->GetCharge();
        hBN=curHadr->GetBaryonNumber();
        G4cout<<"G4QFragmentation::Breeder: H#"<<i<<", d4M="<<curString4M+hI4M
              <<",dCh="<<hCh+curStrChg<<",dBN="<<hBN+curStrBaN<<G4endl;
        if( !(hCh+curStrChg) && !(hBN+curStrBaN) && std::fabs(dEn+hEn)<eps &&
            std::fabs(dPx+hPx)<eps && std::fabs(dPy+hPy)<eps && std::fabs(dPz+hPz)<eps )
        {
          G4cout<<"G4QFragmentation::Breeder: ***Cured*** Redundent Hadron # "<<i<<G4endl;
          G4QHadron* theLast = (*theResult)[nHadr-1];
          curHadr->Set4Momentum(theLast->Get4Momentum()); //4-Mom of CurHadr
          G4QPDGCode lQP=theLast->GetQPDG();
          if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP);
          else curHadr->SetQC(theLast->GetQC());
          theResult->pop_back(); // theLastQHadron is excluded from OUTPUT
          delete theLast;        //*!!When kill, delete theLastQHadr as an Instance!*
          break;
        }
        if( !(hCh+mCh+curStrChg) && !(hBN+mBN+curStrBaN) && std::fabs(dEn+hEn+mEn)<eps &&
            std::fabs(dPx+hPx+mPx)<eps && std::fabs(dPy+hPy+mPy)<eps &&
            std::fabs(dPz+hPz+mPz)<eps && i>0)
        { 
          G4cout<<"G4QFragmentation::Breeder:***Cured*** Redundent 2Hadrons i="<<i<<G4endl;
          G4QHadron* preHadr = (*theResult)[i-1];
          G4QHadron* theLast = (*theResult)[nHadr-1];
          if(i < nHadr-1)        // Only cur can overlap with the two last hadrons
          {                      // Put the last to the previous
            preHadr->Set4Momentum(theLast->Get4Momentum()); // must be 4-Mom of preHadr
            G4QPDGCode lQP=theLast->GetQPDG();
            if(lQP.GetPDGCode()!=10) preHadr->SetQPDG(lQP);
            else preHadr->SetQC(theLast->GetQC());
          }
          theResult->pop_back(); // theLastQHadron's excluded from OUTPUT(even if Cur=Last)
          delete theLast;        //*!!When kill, delete theLastQHadr as an Instance!*
          theLast = (*theResult)[nHadr-2]; // nHadr is not changed -> so it's LastButOne
          if(i < nHadr-2)        // The two current and the two Last are not overlaped
          {                      // Put the last but one to the current
            curHadr->Set4Momentum(theLast->Get4Momentum()); // must be 4-Mom of curHadr
            G4QPDGCode lQP=theLast->GetQPDG();
            if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP);
            else curHadr->SetQC(theLast->GetQC());
          }
          theResult->pop_back(); // theLastQHadron's excluded from OUTPUT(even for overlap)
          delete theLast;        //*!!When kill, delete theLastQHadr as an Instance!*
          nHadr=theResult->size(); // Just a precaution... should be nHadr-2
          break;
        }
        // If the redundent particle decay in 3 FS hadrons -> the corresponding Improvement
        G4cout<<"*Warning*G4QFragmentation::Breeder: Nonconservation isn't cured!"<<G4endl;
      }
    }
    // Trap with the debugging warning ^^^  Ends  ^^^
  } // End of the main LOOP over decaying strings
  G4LorentzVector r4M=theNucleus.Get4Momentum(); // Nucleus 4-momentum in LS
  G4int rPDG=theNucleus.GetPDG();
  G4QHadron* resNuc = new G4QHadron(rPDG,r4M);
  theResult->push_back(resNuc);                          // Fill the residual nucleus
#ifdef edebug
  G4LorentzVector s4M(0.,0.,0.,0.); // Sum of the Result in LS
  G4int rCh=totChg;
  G4int rBN=totBaN;
  G4int nHadr=theResult->size();
  G4int nQuasm=theQuasmons.size();
  G4cout<<"-EMCLS-G4QFragmentation::Breeder:#ofHadr="<<nHadr<<", #OfQuasm="<<nQuasm<<",rN="
        <<r4M.m()<<"="<<G4QNucleus(rPDG).GetGSMass()<<G4endl;
  for(G4int i=0; i<nHadr; i++)
  {
    G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
    s4M+=hI4M;
    G4int hChg=(*theResult)[i]->GetCharge();
    rCh-=hChg;
    G4int hBaN=(*theResult)[i]->GetBaryonNumber();
    rBN-=hBaN;
    G4cout<<"-EMCLS-G4QFragmentation::Breeder:(1) Hadron#"<<i<<", 4M="<<hI4M<<", PDG="
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
    G4cout<<"-EMCLS-G4QFragmentation::Breeder:(1) Quasmon#"<<i<<", 4M="<<hI4M<<", C="<<hChg
          <<", B="<<hBaN<<G4endl;
  }
  G4cout<<"-EMCLS-G4QFragm::Breed: LS r4M="<<s4M-totLS4M<<",rC="<<rCh<<",rB="<<rBN<<G4endl;
#endif
  // Now we need to coolect particles for creation of a Quasmon @@ improve !!
  G4int nRes=theResult->size();
#ifdef ppdebug
  G4cout<<"G4QFragmentation::Breeder: Strings4M="<<ft4M<<", nRes="<<nRes<<G4endl;
#endif
  G4ThreeVector LS3Mom=ft4M.v();
  G4ThreeVector LSDir=LS3Mom.unit();
  if(nRes > 2 && maxEn > 0.)
  {
    std::list<std::pair<G4double, G4QHadron*> > theSorted;         // Output
    std::list<std::pair<G4double, G4QHadron*> >::iterator current; // Input
    for(G4int secondary = 0; secondary<nRes-1; ++secondary)
    {
      G4QHadron*      ih    =theResult->operator[](secondary);
      G4LorentzVector h4M   =ih->Get4Momentum();
      G4double        hM2   =ih->GetMass2();
      G4ThreeVector   h3M   =h4M.v();
      G4double        toSort=DBL_MAX;
      if(hM2>0.00001) toSort=h4M.e()+h3M.dot(LSDir)/std::sqrt(hM2);// monotonic as rapidity
#ifdef ppdebug
      G4cout<<"G4QFragmentation::Breeder:#"<<secondary<<",M2="<<hM2<<",s="<<toSort<<G4endl;
#endif
      std::pair<G4double, G4QHadron*> it;
      it.first      = toSort;
      it.second     = ih;
      G4bool inserted = false;
      for(current = theSorted.begin(); current!=theSorted.end(); ++current)
      {
        if((*current).first > toSort)        // The current is smaller then existing
        {
          theSorted.insert(current, it);     // It shifts the others up
          inserted = true;
          break;
        }
      }
      if(!inserted) theSorted.push_back(it); // It is bigger than any previous
    }
    theResult->clear();                      // Clear and refill theResult by StriHardPart
    G4LorentzVector q4M(0.,0.,0.,0.);
    G4QContent qQC(0,0,0,0,0,0);
    for(current = theSorted.begin(); current!=theSorted.end(); ++current)
    {
      G4QHadron*       ih= (*current).second;
      G4LorentzVector h4M= ih->Get4Momentum();
      G4int          hPDG= ih->GetPDGCode();
      G4double         dE= 0.;
      G4bool       tested=true; 
      if     (hPDG> 1111 && hPDG< 3335) dE=h4M.e()-ih->GetMass(); // Baryons
      else if(hPDG>-1111 && hPDG<1111 && hPDG!=22) dE=h4M.e();    // Mesons (Photons Hard)
      //else if(hPDG<-1111 && hPDG>-3335) dE=h4M.e()+ih->GetMass(); // Antiaryons Don'tUse
      else tested=false;                                          // Skip other
#ifdef ppdebug
      G4cout<<"G4QFragmentation::Breeder:dE="<<dE<<",mE="<<maxEn<<",t="<<tested<<G4endl;
#endif
      if(tested && dE < maxEn)
      {
        maxEn-=dE;
        q4M+=h4M;
        qQC+=ih->GetQC();
#ifdef debug
        G4cout<<"%->G4QFragmentation::Breeder:Exclude,4M="<<h4M<<",dE="<<maxEn<<G4endl;
#endif
        delete ih;
      }
      else  theResult->push_back(ih);                      // Delete equivalent
    } // End of Loop over sorted pairs
    G4Quasmon* softQuasmon = new G4Quasmon(qQC, q4M);      // SoftPart -> Quasmon
#ifdef debug
    G4cout<<"%->G4QFragmentation::Breeder:QuasmonIsFilled,4M="<<q4M<<",QC="<<qQC<<G4endl;
#endif
    if(q4M != vac4M) theQuasmons.push_back(softQuasmon);
    else delete softQuasmon;
    theResult->push_back(resNuc);
#ifdef edebug
    G4LorentzVector f4M(0.,0.,0.,0.);                      // Sum of the Result in LS
    G4int fCh=totChg;
    G4int fBN=totBaN;
    G4int nHd=theResult->size();
    G4int nQm=theQuasmons.size();
    G4cout<<"-EMCLS-G4QFragmentation::Breeder:#ofHadr="<<nHd<<", #OfQuasm="<<nQm<<",rN="
          <<r4M.m()<<"="<<G4QNucleus(rPDG).GetGSMass()<<G4endl;
    for(G4int i=0; i<nHd; i++)
    {
      G4LorentzVector hI4M=(*theResult)[i]->Get4Momentum();
      f4M+=hI4M;
      G4int hChg=(*theResult)[i]->GetCharge();
      fCh-=hChg;
      G4int hBaN=(*theResult)[i]->GetBaryonNumber();
      fBN-=hBaN;
      G4cout<<"-EMCLS-G4QFragmentation::Breeder:(2) Hadron#"<<i<<", 4M="<<hI4M<<", PDG="
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
      G4cout<<"-EMCLS-G4QFragmentation::Breeder:(2) Quasmon#"<<i<<", 4M="<<hI4M<<", C="
            <<hChg<<", B="<<hBaN<<G4endl;
    }
    G4cout<<"-EMCLS-G4QFragm::Breed:, r4M="<<f4M-totLS4M<<",rC="<<fCh<<",rB="<<fBN<<G4endl;
#endif
  } // End of the soft Quasmon Creation
  return;
} // End of Breeder

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
  G4cout<<"G4QFragment::ExciteDiffParticipantts:BeforBoost Pproj="<<Pprojectile<<", Ptarg="
        <<Ptarget<<G4endl;
#endif
  G4LorentzRotation toLab(toCms.inverse()); // Boost Rotation to LabSys (LS)
  Pprojectile.transform(toCms);
  Ptarget.transform(toCms);
#ifdef debug
  G4cout<< "G4QFragment::ExciteDiffParticipantts: AfterBoost Pproj="<<Pprojectile<<"Ptarg="
        <<Ptarget<<", cms4M="<<Pprojectile+Ptarget<<G4endl;
  G4cout<<"G4QFragment::ExciteDiffParticipants: ProjX+="<<Pprojectile.plus()<<", ProjX-="
        <<Pprojectile.minus()<<", TargX+="<< Ptarget.plus()<<", TargX-="<<Ptarget.minus()
        <<G4endl;
#endif
  G4LorentzVector Qmomentum(0.,0.,0.,0.);
  G4int whilecount=0;
#ifdef debug
  G4cout<<"G4QFragmentation::ExciteDiffParticipants: Before DO"<<G4endl;
#endif
  do
  {
    //  Generate pt  
    G4double maxPtSquare=sqr(Ptarget.pz());
#ifdef debug
    G4cout<<"G4QFragmentation::ExciteDiffParticipants: maxPtSq="<<maxPtSquare<<G4endl;
    if(whilecount++>=500 && whilecount%100==0) // @@ M.K. Hardwired limits 
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

void G4QFragmentation::SetParameters(G4int nC, G4double stT,  G4double tbD, G4double SigPt)
{
  nCutMax            = nC;             // max number of pomeron cuts
  stringTension      = stT;            // string tension for absorbed energy
  tubeDensity        = tbD;            // Flux Tube Density of nuclear nucleons
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
#ifdef debug
  G4cout<<"G4QFragmentation::GaussianPt:width2="<<widthSq<<",maxPt2="<<maxPtSquare<<G4endl;
#endif
  G4double pt2=0.;
  G4double rm=maxPtSquare/widthSq;                      // Negative
  if(rm>-.01) pt2=widthSq*(std::sqrt(1.-G4UniformRand()*(1.-sqr(1.+rm)))-1.);
  else        pt2=widthSq*std::log(1.-G4UniformRand()*(1.-std::exp(rm)));
  pt2=std::sqrt(pt2);                                   // It is not pt2, but pt now
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
  //#ifdef debug
  G4cout<<"-Warning-G4QFragmentation::ReducePair:**Failed**, P1="<<P1<<", P2="<<P2<<G4endl;
  //#endif
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

void G4QFragmentation::SwapPartons() // Swap string partons, if a string has negative M2
{
  static const G4double baryM=800.;                  // Mass excess for baryons
  G4QStringVector::iterator ist;
  for(ist = strings.begin(); ist < strings.end(); ist++)
  {
    G4QParton* cLeft=(*ist)->GetLeftParton();        // Current String Left Parton 
    G4QParton* cRight=(*ist)->GetRightParton();      // Current String Right Parton 
    G4LorentzVector cL4M=cLeft->Get4Momentum();
    G4LorentzVector cR4M=cRight->Get4Momentum();
    G4LorentzVector cS4M=cL4M+cR4M;
    G4double cSM2=cS4M.m2();                         // Squared mass of the String
    if(std::fabs(cSM2)<.01)                          // Small correction
    {
      G4double dM2=.001-cSM2;
      G4double E=cS4M.e();
      G4double dE=std::sqrt(E*E+dM2)-E;
      G4double LE=cL4M.e();
      G4double RE=cR4M.e();
      if(LE<RE) cLeft->Set4Momentum( G4LorentzVector(LE+dE) );
      else      cRight->Set4Momentum( G4LorentzVector(RE+dE) );
      cSM2=.001;                                     // Correction
    }
    if(cSM2<0.011)                                   // Parton Swapping is needed
    {
      G4int cLPDG=cLeft->GetPDGCode();
      G4int cRPDG=cRight->GetPDGCode();
      G4int cLT=cLeft->GetType();
      G4int cRT=cRight->GetType();
      G4QStringVector::iterator sst;                 // Selected partner string
      G4QStringVector::iterator pst;                 // LOOP iterator
      G4double maxM=-DBL_MAX;                        // Swapping providing the max mass
      G4int    sDir=0;                               // Selected direction of swapping
#ifdef debug
      G4cout<<"G4QFragmentation::SwapPartons: M2=="<<cSM2<<", 4M="<<cS4M<<",LPDG="<<cLPDG
            <<",RPDG="<<cRPDG<<G4endl;
#endif
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
#ifdef debug
        G4cout<<"G4QFragmentation::SwapPartons: p4M="<<cS4M<<",pM2="<<cS4M.m2()<<",LPDG="
              <<pLPDG<<",RPDG="<<pRPDG<<G4endl;
#endif
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
      else G4cout<<"***G4QFragmentation::SwapPartons:**Failed**,cLPDG="<<cLPDG<<",cRPDG="
                 <<cRPDG<<",-->cM2="<<cSM2<<G4endl;
#endif
    }
  }
}

//Evaporate Residual Nucleus/Fragment (@@Modified function from G4QEnvironment can be used)
void G4QFragmentation::EvaporateResidual(G4QHadron* qH)
{
  static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mNeut = G4QPDGCode(2112).GetMass();
  static const G4double mProt = G4QPDGCode(2212).GetMass();
  static const G4double mAlPr = mAlph+mProt;
  static const G4double mAlNt = mAlph+mNeut;
  static const G4double dProt = mProt+mProt;
  static const G4double dNeut = mNeut+mNeut;
  static const G4double dAlph = mAlph+mAlph;
  static const G4double eps=.003;
  G4QEnvironment envir(theNucleus);
  G4int thePDG = qH->GetPDGCode();           // Get PDG code of the Residual Nucleus
  G4int theBN  = qH->GetBaryonNumber();      // A (Baryon number of the nucleus)
  G4QContent  theQC  = qH->GetQC();          // Quark Content of the hadron
  G4int theS=theQC.GetStrangeness();         // S (Strangeness of the nucleus)
#ifdef debug
  G4cout<<"G4QFragment::EvaRes:-Called- PDG="<<thePDG<<",4M="<<qH->Get4Momentum()
        <<",QC="<<theQC<<", BN="<<theBN<<G4endl;
#endif
  if(thePDG==10)
  {
#ifdef debug
    G4cout<<"G4QFragment::EvaRes: Cgipolino QC="<<theQC<<qH->Get4Momentum()<<G4endl;
#endif
    G4QContent   chQC=qH->GetQC();           // Quark content of the Hadron-Chipolino
    G4QChipolino QCh(chQC);                  // Define a Chipolino instance for the Hadron
    G4LorentzVector ch4M=qH->Get4Momentum(); // 4Mom of the Hadron-Chipolino
    G4QPDGCode h1QPDG=QCh.GetQPDG1();        // QPDG of the first hadron
    G4double   h1M   =h1QPDG.GetMass();      // Mass of the first hadron
    G4QPDGCode h2QPDG=QCh.GetQPDG2();        // QPDG of the second hadron
    G4double   h2M   =h2QPDG.GetMass();      // Mass of the second hadron
    G4double   chM2  =ch4M.m2();             // Squared Mass of the Chipolino
    if( sqr(h1M+h2M) < chM2 )                // Decay is possible
    {
      G4LorentzVector h14M(0.,0.,0.,h1M);
      G4LorentzVector h24M(0.,0.,0.,h2M);
      if(!G4QHadron(ch4M).DecayIn2(h14M,h24M))
      {
        G4cerr<<"***G4QFrag::EvaporateResid: CM="<<std::sqrt(chM2)<<" -> h1="<<h1QPDG<<"("
              <<h1M<<") + h2="<<h1QPDG<<"("<<h2M<<") = "<<h1M+h2M<<" **Failed**"<<G4endl;
        // throw G4QException("*G4QFragmentation::EvaporateResidual:QChipolino DecIn2 error");
        G4Exception("G4QFragmentation::EvaporateResidual()", "HAD_CHPS_0000",
                    FatalException, "QChipolino DecIn2 error");
      }
      delete qH;                             // Kill the primary Chipolino
      G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
      theResult->push_back(h1H);             // (delete equivalent)
#ifdef debug
      G4cout<<"G4QFragm::EvaporateResidual: Chipolino -> H1="<<h1QPDG<<h14M<<G4endl;
#endif
      qH = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
      theResult->push_back(qH);              // (delete equivalent)
#ifdef debug
      G4cout<<"G4QE::EvaporateResidual: Chipolino -> H2="<<h2QPDG<<h24M<<G4endl;
#endif
    }
    else
    {
      G4cerr<<"***G4QFragment::EvaporateResid: Chipolino="<<qH->GetQC()<<qH->Get4Momentum()
            <<", chipoM="<<std::sqrt(chM2)<<" < m1="<<h1M<<"("<<h1QPDG<<") + m2="<<h2M
            <<"("<<h2QPDG<<") = "<<h1M+h2M<<G4endl;
      // throw G4QException("G4QFragmentation::EvaporateResidual: LowMassChipolino in Input");
      G4Exception("G4QFragmentation::EvaporateResidual()", "HAD_CHPS_0001",
                  FatalException, "LowMassChipolino in Input");
    }
    return;
  }
  else if(theS<0)                            // Antistrange nucleus
  {
#ifdef debug
    G4cout<<"G4QFragment::EvaRes: AntistrangeNucleus="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
    envir.DecayAntistrange(qH, theResult);   // (delete equivalent)
    return;
  }
  else if(theBN==1)
  {
#ifdef debug
    G4cout<<"G4QFragmentation::EvaporateResid:Baryon="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
    envir.DecayBaryon(qH, theResult);        // (delete equivalent)
    return;
  }
  else if(!theBN) // @@ In future it is usefull to add the MesonExcitationDecay (?!)
  {
#ifdef debug
    G4LorentzVector mesLV=qH->Get4Momentum();
    G4cout<<"G4QFragmentation::EvaRes:(!)Meson(!) PDG="<<thePDG<<",4M="<<mesLV<<mesLV.m()
          <<",QC="<<qH->GetQC()<<",MPDG="<<G4QPDGCode(thePDG).GetMass()<<G4endl;
#endif
    envir.DecayMeson(qH, theResult);         // @@ To be written
    return;
  }
  G4int theC=theQC.GetCharge();              // P
#ifdef debug
  G4cout<<"G4QFragment::EvaRes: qH.Charge = "<<theC<<G4endl;
#endif
  if(!thePDG) thePDG = theQC.GetSPDGCode();  // If there is no PDG code, get it from QC
  if( thePDG == 10 && theBN > 0 ) thePDG=theQC.GetZNSPDGCode();
  if(theS>0) thePDG-=theS*999999;            // @@ May hide hypernuclear problems (G4) ! @@
#ifdef debug
  G4cout<<"G4QFragment::EvaRes: S="<<theS<<", newPDG="<<thePDG<<G4endl;
#endif
  G4double totGSM = G4QNucleus(thePDG).GetGSMass();// TheGroundStMass of theTotalResNucleus
#ifdef debug
  G4cout<<"G4QFragment::EvaRes: totGSM="<<totGSM<<G4endl;
#endif
  if(theBN==2)
  {
    if(!theC)        totGSM=dNeut;           // nn, nL, LL
    else if(theC==2) totGSM=dProt;           // pp
    else             totGSM=mDeut;           // np, Lp
  }
  else if(theBN==5)
  {
    if     (theC==3) totGSM=mAlPr;           // effective "Alph+p"
    else if(theC==2) totGSM=mAlNt;           // effective "Alph+n"
  }
  else if(theBN==8)   totGSM=dAlph;          // effective "Be8"
  // @@ Should be more (else if) for bigger A=theBN
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of theTotalResidNucleus
  G4double    totMass = q4M.m();             // Get theRealMass of theTotalResidNucleus
#ifdef debug
    G4cout<<"G4QFragment::EvaRes: Excitation = "<<totMass-totGSM<<G4endl;
#endif
  if(std::fabs(totMass-totGSM) < eps)
  {
    theResult->push_back(qH);               // fill As It Is
  }
  else if(totMass > totGSM)
  {
#ifdef debug
    G4cout<<"G4QFragment::EvaRes: try Evaporate Nucleus PDG="<<thePDG<<G4endl;
#endif
    theNucleus.EvaporateNucleus(qH, theResult);
#ifdef debug
    G4cout<<"G4QFragment::EvaRes: ** Evaporation is done **"<<G4endl;
#endif
    //delete qH;
    qH=0;
  }
  else                                       // Correction must be done
  {
#ifdef debug
    G4cout<<"-War-G4QFr::EvaRes:*Must correct* "<<theQC<<q4M<<totMass<<"<"<<totGSM<<G4endl;
#endif
  }
#ifdef qdebug
  if (qH)
  {
    G4cout<<"-W-G4QFragmentation::EvaporateResid:*Deleted*,PDG="<<qH->GetPDGCode()<<G4endl;
    delete qH;
  }
#endif
  return;
} // End of EvaporateResidual

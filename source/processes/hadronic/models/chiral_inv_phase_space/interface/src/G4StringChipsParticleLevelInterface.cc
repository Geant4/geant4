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
// Short description: Interface of QGSC to CHIPS (SoftHadrons) 
// ----------------------------------------------------------------
//

//#define debug
//#define trapdebug
//#define pdebug
//#define ppdebug

#include <utility>
#include <list>
#include <vector>

#include "G4StringChipsParticleLevelInterface.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4KineticTrackVector.hh"
#include "G4Nucleon.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4LorentzRotation.hh"
#include "G4HadronicException.hh"
// #define CHIPSdebug
// #define CHIPSdebug_1

#ifdef hdebug_SCPLI
const G4int G4StringChipsParticleLevelInterface::nbh=200;
G4double    G4StringChipsParticleLevelInterface::bhmax=20.;
G4double    G4StringChipsParticleLevelInterface::ehmax=20.;
G4double    G4StringChipsParticleLevelInterface::bhdb=0.;
G4double    G4StringChipsParticleLevelInterface::ehde=0.;
G4double    G4StringChipsParticleLevelInterface::toth=0.;
G4int       G4StringChipsParticleLevelInterface::bover=0;
G4int       G4StringChipsParticleLevelInterface::eover=0;
G4int*      G4StringChipsParticleLevelInterface::bhis =
                                       new G4int[G4StringChipsParticleLevelInterface::nbh];
G4int*      G4StringChipsParticleLevelInterface::ehis =
                                       new G4int[G4StringChipsParticleLevelInterface::nbh];
#endif

G4StringChipsParticleLevelInterface::G4StringChipsParticleLevelInterface()
{
#ifdef debug
  G4cout<<"G4StringChipsParticleLevelInterface::Constructor is called"<<G4endl;
#endif
  //theEnergyLossPerFermi = 1.*GeV;
  theEnergyLossPerFermi = 1.5*GeV;
  nop = 152;                               // clusters (A<6)
  fractionOfSingleQuasiFreeNucleons = 0.5; // It is A-dependent (C=.85, U=.40) - M.K.
  fractionOfPairedQuasiFreeNucleons = 0.05;
  clusteringCoefficient = 5.;
  temperature = 180.;
  halfTheStrangenessOfSee = 0.3; // = s/d = s/u
  etaToEtaPrime = 0.3;
  fusionToExchange = 1.;
  //theInnerCoreDensityCut = 50.;
  theInnerCoreDensityCut = 70.;
  
  if(getenv("ChipsParameterTuning"))
  {
    G4cout << "Please enter the energy loss per fermi in GeV"<<G4endl;
    G4cin >> theEnergyLossPerFermi;
    theEnergyLossPerFermi *= GeV;
    G4cout << "Please enter nop"<<G4endl;
    G4cin >> nop;
    G4cout << "Please enter the fractionOfSingleQuasiFreeNucleons"<<G4endl;
    G4cin >> fractionOfSingleQuasiFreeNucleons;
    G4cout << "Please enter the fractionOfPairedQuasiFreeNucleons"<<G4endl;
    G4cin >> fractionOfPairedQuasiFreeNucleons;
    G4cout << "Please enter the clusteringCoefficient"<<G4endl;
    G4cin >> clusteringCoefficient;
    G4cout << "Please enter the temperature"<<G4endl;
    G4cin >> temperature;
    G4cout << "Please enter halfTheStrangenessOfSee"<<G4endl;
    G4cin >> halfTheStrangenessOfSee;
    G4cout << "Please enter the etaToEtaPrime"<<G4endl;
    G4cin >> etaToEtaPrime;
    G4cout << "Please enter the fusionToExchange"<<G4endl;
    G4cin >> fusionToExchange;
    G4cout << "Please enter the cut-off for calculating the nuclear radius in percent"<<G4endl;
    G4cin >> theInnerCoreDensityCut;
  }
}

G4HadFinalState* G4StringChipsParticleLevelInterface::
ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
#ifdef debug
  G4cout<<"G4StringChipsParticleLevelInterface::ApplyYourself is called"<<G4endl;
#endif
  return theModel.ApplyYourself(aTrack, theNucleus);
}

G4ReactionProductVector* G4StringChipsParticleLevelInterface::
Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
  static const G4double mProt=G4Proton::Proton()->GetPDGMass();
  static const G4double mNeut=G4Neutron::Neutron()->GetPDGMass();
  static const G4double mLamb=G4Lambda::Lambda()->GetPDGMass();
  static const G4double mKChg=G4KaonPlus::KaonPlus()->GetPDGMass();
  static const G4double mKZer=G4KaonZero::KaonZero()->GetPDGMass();
  static const G4double mPiCh=G4PionMinus::PionMinus()->GetPDGMass();
  static const G4int pcl=4; // clusterization parameter for Energy Flow
  static const G4QContent ProtQC(1,2,0,0,0,0);
  static const G4QContent NeutQC(2,1,0,0,0,0);
  static const G4QContent LambQC(1,1,1,0,0,0);
  static const G4QContent KPlsQC(0,1,0,0,0,1);
  static const G4QContent KMinQC(0,0,1,0,1,0);
  static const G4QContent AKZrQC(1,0,0,0,0,1);
  static const G4QContent KZerQC(1,0,0,0,0,1);
  static const G4QContent PiMiQC(1,0,0,0,1,0);
  static const G4QContent PiPlQC(0,1,0,1,0,0);
#ifdef debug
  G4cout<<"G4StringChipsParticleLevelInterface::Propagate is called"<<G4endl;
#endif
  // Protection for non physical conditions
  
  if(theSecondaries->size() == 1) 
  {
    G4ReactionProductVector* theFastResult = new G4ReactionProductVector;
    G4ReactionProduct* theFastSec;
    theFastSec = new G4ReactionProduct((*theSecondaries)[0]->GetDefinition());
    G4LorentzVector current4Mom = (*theSecondaries)[0]->Get4Momentum();
    theFastSec->SetTotalEnergy(current4Mom.t());
    theFastSec->SetMomentum(current4Mom.vect());
    theFastResult->push_back(theFastSec);
    return theFastResult;
    //throw G4HadronicException(__FILE__,__LINE__,
    //       "G4StringChipsParticleLevelInterface: Only one particle from String models!");
  }
  
  // target properties needed in constructor of quasmon, and for boosting to
  // target rest frame
  // remove all nucleons already involved in STRING interaction, to make the ResidualTarget
  theNucleus->StartLoop();
  G4Nucleon * aNucleon;
  G4int resA = 0;
  G4int resZ = 0;
  G4ThreeVector hitMomentum(0,0,0);
  G4double hitMass = 0;
  unsigned int hitCount = 0;
  while((aNucleon = theNucleus->GetNextNucleon()))
  {
    if(!aNucleon->AreYouHit())
    {
      resA++;                                                  // Collect A of the ResidNuc
      resZ+=G4int (aNucleon->GetDefinition()->GetPDGCharge()); // Collect Z of the ResidNuc
    }
    else
    {
      hitMomentum += aNucleon->GetMomentum().vect();           // Sum 3-mom of StringHadr's
      hitMass += aNucleon->GetMomentum().m();                  // Sum masses of StringHadrs
      hitCount ++;                                             // Calculate STRING hadrons
    }
  }
  G4int targetPDGCode = 90000000 + 1000*resZ + (resA-resZ);    // PDG of theResidualNucleus
  G4double targetMass=mNeut;
  if (!resZ)                                                   // Nucleus of only neutrons
  {
    if (resA>1) targetMass*=resA;
  }
  else targetMass=G4ParticleTable::GetParticleTable()->FindIon(resZ,resA,0,resZ)
                                                                            ->GetPDGMass();
  G4double targetEnergy = std::sqrt(hitMomentum.mag2()+targetMass*targetMass);
  // !! @@ Target should be at rest: hitMomentum=(0,0,0) @@ !! M.K. (go to this system)
  G4LorentzVector targ4Mom(-1.*hitMomentum, targetEnergy);
  
  // Calculate the mean energy lost
  std::pair<G4double, G4double> theImpact = theNucleus->RefetchImpactXandY();
  G4double impactX = theImpact.first;
  G4double impactY = theImpact.second;
  G4double impactPar2 = impactX*impactX + impactY*impactY;
  G4double radius2 = theNucleus->GetNuclearRadius(theInnerCoreDensityCut*perCent);
  //G4double radius2 = theNucleus->GetNuclearRadius(theInnerCoreDensityCut*perCent);
  radius2 *= radius2;
  G4double pathlength = 0.;
#ifdef ppdebug
  G4cout<<"G4StringChipsParticleLevelInterface::Propagate: r="<<std::sqrt(radius2)/fermi
        <<", b="<<std::sqrt(impactPar2)/fermi<<", R="<<theNucleus->GetOuterRadius()/fermi
        <<", b/r="<<std::sqrt(impactPar2/radius2)<<G4endl; 
#endif
  if(radius2 - impactPar2>0) pathlength = 2.*std::sqrt(radius2 - impactPar2);
  //pathlength = 27.;                         // *** Only temporary *** Always CHIPS
#ifdef hdebug_SCPLI
  toth+=1.;                                 // increment total number of measurements
  G4double bfm=std::sqrt(impactPar2)/fermi; // impact parameter
  G4double efm=pathlength/fermi;            // energy absorption length
  G4int    nbi=static_cast<G4int>(bfm/bhdb);
  G4int    nei=static_cast<G4int>(efm/ehde);
  if(nbi<nbh) bhis[nbi]++;
  else        bover++;
  if(nei<nbh) ehis[nei]++;
  else        eover++;
#endif
  G4double theEnergyLostInFragmentation = theEnergyLossPerFermi*pathlength/fermi;
  
  // now select all particles in range
  std::list<std::pair<G4double, G4KineticTrack *> > theSorted;         // Output
  std::list<std::pair<G4double, G4KineticTrack *> >::iterator current; // Input
  for(unsigned int secondary = 0; secondary<theSecondaries->size(); secondary++)
  {
    G4LorentzVector a4Mom = theSecondaries->operator[](secondary)->Get4Momentum();
#ifdef CHIPSdebug
    G4cout<<"G4StringChipsParticleLevelInterface::Propagate: ALL STRING particles "
          << theSecondaries->operator[](secondary)->GetDefinition()->GetPDGCharge()<<" "
          << theSecondaries->operator[](secondary)->GetDefinition()->GetPDGEncoding()<<" "
          << a4Mom <<G4endl; 
#endif
#ifdef pdebug
    G4cout<<"G4StringChipsParticleLevelInterface::Propagate: in C="
          <<theSecondaries->operator[](secondary)->GetDefinition()->GetPDGCharge()<<",PDG="
          <<theSecondaries->operator[](secondary)->GetDefinition()->GetPDGEncoding()
          <<",4M="<<a4Mom<<", current nS="<<theSorted.size()<<G4endl; 
#endif
    G4double toSort = a4Mom.rapidity();          // Rapidity is used for the ordering (?!)
    std::pair<G4double, G4KineticTrack *> it;
    it.first = toSort;
    it.second = theSecondaries->operator[](secondary);
    G4bool inserted = false;
    for(current = theSorted.begin(); current!=theSorted.end(); current++)
    {
      if((*current).first > toSort)        // The current is smaller then existing
      {
        theSorted.insert(current, it);     // It shifts the others up
        inserted = true;
        break;
      }
    }
    if(!inserted) theSorted.push_back(it); // It is bigger than any previous (@@Check this)
  }
  
  G4LorentzVector proj4Mom(0.,0.,0.,0.);
  G4int nD  = 0;
  G4int nU  = 0;
  G4int nS  = 0;
  G4int nAD = 0;
  G4int nAU = 0;
  G4int nAS = 0;
  std::list<std::pair<G4double,G4KineticTrack*> >::iterator firstEscape=theSorted.begin();
  G4double runningEnergy = 0;
  G4int particleCount = 0;
  G4LorentzVector theLow = (*(theSorted.begin())).second->Get4Momentum();
  G4LorentzVector theHigh;

#ifdef CHIPSdebug
  G4cout<<"G4StringChipsParticleLevelInterface::Propagate: CHIPS ENERGY LOST "
        <<theEnergyLostInFragmentation<<". Sorted rapidities event start"<<G4endl;
#endif
#ifdef pdebug
  G4cout<<"G4StringChipsParticleLevelInterface::Propagate: total CHIPS energy = "
        <<theEnergyLostInFragmentation<<". Start rapidity sorting nS="<<theSorted.size()
        <<G4endl;
#endif

  G4QHadronVector projHV;
  std::vector<G4QContent> theContents;
  std::vector<G4LorentzVector*> theMomenta;
  G4ReactionProductVector* theResult = new G4ReactionProductVector;
  G4ReactionProduct* theSec;
  G4KineticTrackVector* secondaries;
  G4KineticTrackVector* secsec;
#ifdef pdebug
  G4cout<<"G4StringChipsParticleLevelInterface::Propagate: Absorption nS="
        <<theSorted.size()<<G4endl; 
#endif
  G4bool EscapeExists = false;
  for(current = theSorted.begin(); current!=theSorted.end(); current++)
  {
#ifdef pdebug
    G4cout<<"G4StringChipsParticleLevelInterface::Propagate: nq="
          <<(*current).second->GetDefinition()->GetQuarkContent(3)<<", naq="
          <<(*current).second->GetDefinition()->GetAntiQuarkContent(3)<<", PDG="
          <<(*current).second->GetDefinition()->GetPDGEncoding()<<",4M="
          <<(*current).second->Get4Momentum()<<G4endl; 
#endif
    firstEscape = current;              // Remember to make decays for the rest
    G4KineticTrack* aResult = (*current).second;
    // This is an old (H.P.) solution, which makes an error in En/Mom conservation
    //
    // @@ Now it does not include strange particle for the absorption in nuclei (?!) M.K.
    //if((*current).second->GetDefinition()->GetQuarkContent(3)!=0 || 
    //   (*current).second->GetDefinition()->GetAntiQuarkContent(3) !=0) // Strange quarks
    //{
    //  G4ParticleDefinition* pdef = aResult->GetDefinition();
    //  secondaries = NULL;
    //  if ( pdef->GetPDGWidth() > 0 && pdef->GetPDGLifeTime() < 5E-17*s )
    //     secondaries = aResult->Decay(); // @@ Decay of only strange resonances (?!) M.K.
    //  if ( secondaries == NULL )           // No decay
    //  {
    //    theSec = new G4ReactionProduct(aResult->GetDefinition());
    //    G4LorentzVector current4Mom = aResult->Get4Momentum();
    //    current4Mom.boost(targ4Mom.boostVector()); // boost from the targetAtRes system
    //    theSec->SetTotalEnergy(current4Mom.t());
    //    theSec->SetMomentum(current4Mom.vect());
    //    theResult->push_back(theSec);
    //  } 
    //  else                                 // The decay happened
    //  {
    //    for (unsigned int aSecondary=0; aSecondary<secondaries->size(); aSecondary++)
    //    {
    //      theSec = 
    //           new G4ReactionProduct(secondaries->operator[](aSecondary)->GetDefinition());
    //      G4LorentzVector current4Mom =secondaries->operator[](aSecondary)->Get4Momentum();
    //      current4Mom.boost(targ4Mom.boostVector());
    //      theSec->SetTotalEnergy(current4Mom.t());
    //      theSec->SetMomentum(current4Mom.vect());
    //      theResult->push_back(theSec);
    //    }
    //    std::for_each(secondaries->begin(), secondaries->end(), DeleteKineticTrack());
    //    delete secondaries;
    //  }
    //}
    //
    //runningEnergy += (*current).second->Get4Momentum().t();
    //if((*current).second->GetDefinition() == G4Proton::Proton())
    //                                   runningEnergy-=G4Proton::Proton()->GetPDGMass();
    //if((*current).second->GetDefinition() == G4Neutron::Neutron())
    //                                   runningEnergy-=G4Neutron::Neutron()->GetPDGMass();
    //if((*current).second->GetDefinition() == G4Lambda::Lambda())
    //                                   runningEnergy-=G4Lambda::Lambda()->GetPDGMass();
    //
    // New solution starts from here (M.Kossov March 2006) [Strange particles included]
    runningEnergy += aResult->Get4Momentum().t();
    G4double charge=aResult->GetDefinition()->GetPDGCharge(); // Charge of the particle
    G4int strang=aResult->GetDefinition()->GetQuarkContent(3);// Its strangeness
    G4int baryn=aResult->GetDefinition()->GetBaryonNumber();  // Its baryon number
    if     (baryn>0 && charge>0 && strang<1) runningEnergy-=mProt;  // For positive baryons
    else if(baryn>0 && strang<1) runningEnergy-=mNeut;              // For neut/neg baryons
    else if(baryn>0) runningEnergy-=mLamb;                          // For strange baryons
    else if(baryn<0) runningEnergy+=mProt;                          // For anti-particles
    // ------------ End of the new solution
#ifdef CHIPSdebug
    G4cout<<"G4StringChipsParticleLevelInterface::Propagate: sorted rapidities "
                                    <<(*current).second->Get4Momentum().rapidity()<<G4endl;
#endif

#ifdef pdebug
    G4cout<<"G4StringChipsParticleLevelInterface::Propagate: E="<<runningEnergy<<", EL="
                                                    <<theEnergyLostInFragmentation<<G4endl;
#endif
    if(runningEnergy > theEnergyLostInFragmentation)
    {
      EscapeExists = true;
      break;
    }
#ifdef CHIPSdebug
    G4cout <<"G4StringChipsParticleLevelInterface::Propagate: ABSORBED STRING particles "
           <<(*current).second->GetDefinition()->GetPDGCharge()<<" "
           << (*current).second->GetDefinition()->GetPDGEncoding()<<" "
           << (*current).second->Get4Momentum() <<G4endl; 
#endif
#ifdef pdebug
    G4cout<<"G4StringChipsParticleLevelInterface::Propagate:C="
          <<current->second->GetDefinition()->GetPDGCharge()<<", PDG="
          <<current->second->GetDefinition()->GetPDGEncoding()<<", 4M="
          <<current->second->Get4Momentum()<<G4endl; 
#endif

    // projectile 4-momentum in target rest frame needed in constructor of QHadron
    particleCount++;
    theHigh = (*current).second->Get4Momentum(); 
    proj4Mom = (*current).second->Get4Momentum(); 
    proj4Mom.boost(-1.*targ4Mom.boostVector());   // Back to the system of nucleusAtRest
    nD = (*current).second->GetDefinition()->GetQuarkContent(1);
    nU = (*current).second->GetDefinition()->GetQuarkContent(2);
    nS = (*current).second->GetDefinition()->GetQuarkContent(3);
    nAD = (*current).second->GetDefinition()->GetAntiQuarkContent(1);
    nAU = (*current).second->GetDefinition()->GetAntiQuarkContent(2);
    nAS = (*current).second->GetDefinition()->GetAntiQuarkContent(3);
    G4QContent aProjectile(nD, nU, nS, nAD, nAU, nAS);

#ifdef CHIPSdebug_1
    G4cout <<G4endl;
    G4cout <<"G4StringChipsParticleLevelInterface::Propagate: Quark content: d="<<nD
           <<", u="<<nU<<", s="<<nS<< "Anti-quark content: anit-d="<<nAD<<", anti-u="<<nAU
           <<", anti-s="<<nAS<<". G4QContent is constructed"<<endl;
#endif

    theContents.push_back(aProjectile);
    G4LorentzVector* aVec = new G4LorentzVector((1./MeV)*proj4Mom); // @@ MeV is basic

#ifdef CHIPSdebug_1
    G4cout<<"G4StringChipsParticleLevelInterface::Propagate: projectile momentum = "
          <<*aVec<<G4endl;
    G4cout << G4endl;
#endif
   
    theMomenta.push_back(aVec);
  }
  std::vector<G4QContent> theFinalContents;
  std::vector<G4LorentzVector*> theFinalMomenta;

  // Multiquasmon case: each particle creates a quasmon
  //for(unsigned int hp = 0; hp<theContents.size(); hp++)
  //{
  //  G4QHadron* aHadron = new G4QHadron(theContents[hp], *(theMomenta[hp]) );
  //  projHV.push_back(aHadron);
  //}
  // Energy flow: one Quasmon for each B>0 collection ----------
  G4QContent EnFlowQC(0,0,0,0,0,0);
  G4LorentzVector EnFlow4M(0.,0.,0.,0.);
  //G4bool empty=true;
  G4int  barys=0;
  G4int  stras=0;
  G4int  chars=0;
  for(G4int hp = theContents.size()-1; hp>=0; hp--)
  {
     G4QContent curQC=theContents[hp];
     G4int baryn = curQC.GetBaryonNumber();
     G4int stran = curQC.GetStrangeness();
     G4int charg = curQC.GetCharge();
     EnFlowQC += curQC;                        // Keep collecting energy flow
     EnFlow4M += *(theMomenta[hp]);
     barys += baryn;                           // Collected baryon number
     stras += stran;                           // Collected strangeness
     chars += charg;                           // Collected charge
     //empty = false;
  }
  if(barys>pcl)                               // Split in two or more parts (to survive!)
  {
     G4int nprt=(barys-1)/pcl+1;               // Number of parts (pcl=4: 2:5-8,3:9-12...)
     G4int curb=barys;
     while (nprt>0)
     {
        nprt--;                                 // One part is going to be created
        G4int brnm=pcl;                         // Baryon number of splitting part
        curb-=brnm;                             // The residual baryon number
        G4double prtM=0.;                       // The resulting GS mass of the part
        G4double resM=0.;                       // The resulting GS mass of the residual
        G4QContent prtQC(0,0,0,0,0,0);          // The resulting Quark Content of the part
        G4int strm=0;                           // Max strangeness per part (stras=0)
        if(stras>0) strm=(stras-1)/nprt+1;      // Max strangeness per part (stras>0)
        else if(stras<0) strm=(stras+1)/nprt-1; // Max strangeness per part (stras<0)
        G4int chgm=0;                           // Max charge per part (chars=0)
        if(stras>0) chgm=(chars-1)/nprt+1;      // Max strangeness per part (chars>0)
        else if(stras<0) chgm=(chars+1)/nprt-1; // Max strangeness per part (chars<0)
        // ---> calculate proposed separated part
        //@@ Convert it to a CHIPS function (Which class? G4QH::Conctruct?)
        if(!strm)                               // --> The total strangness = 0 (n/p/pi-)
        {
          if(chgm<0)                            // (n/pi-)
          {
            prtM=(-chgm)*mPiCh+brnm*mNeut;
            prtQC=(-chgm)*PiMiQC+brnm*NeutQC;
          }
          else                                  // (n/p)
          {
            prtM=chgm*mProt+(brnm-chgm)*mNeut;
            prtQC=chgm*ProtQC+(brnm-chgm)*NeutQC;
          }
        }
        else if(strm>=brnm)                     // ---> BigPositiveStrangeness(L/Pi+/K0/K-)
        {
          G4int stmb=strm-brnm;
          if(chgm<0)                            // (L/K-/K0)
          {
            prtM=(-chgm)*mKChg+brnm*mLamb+std::abs(stmb+chgm)*mKZer;
            prtQC=(-chgm)*KMinQC+brnm*LambQC;
            if(stmb>-chgm) prtQC+=(stmb+chgm)*KZerQC;
            else if(stmb<-chgm) prtQC+=(-stmb-chgm)*AKZrQC;
          }
          else                                  // (L/K0/pi+)
          {
            prtM=chgm*mPiCh+(strm-brnm)*mKZer+brnm*mLamb;
            prtQC=chgm*PiPlQC+(strm-brnm)*KZerQC+brnm*LambQC;
          }
        }
        else if(strm>0)                         // ---> PositiveStrangeness<B (L/n/p/Pi+-)
        {
          G4int bmst=brnm-strm;
          if(chgm<0)                            // (L/n/Pi-)
          {
            prtM=(-chgm)*mPiCh+strm*mLamb+bmst*mNeut;
            prtQC=(-chgm)*PiMiQC+strm*LambQC+bmst*NeutQC;
          }
          else if(chgm>=bmst)                   // (L/p/Pi+)
          {
            prtM=(chgm-bmst)*mPiCh+strm*mLamb+bmst*mProt;
            prtQC=(chgm-bmst)*PiPlQC+strm*LambQC+bmst*ProtQC;
          }
          else                                  // ch<bmst (L/p/n)
          {
            prtM=chgm*mProt+strm*mLamb+(bmst-chgm)*mNeut;
            prtQC=chgm*ProtQC+strm*LambQC+(bmst-chgm)*NeutQC;
          }
        }
        else                                    // ---> NegativeStrangeness (N/K+/aK0/Pi-)
        {
          G4int bmst=brnm-strm;
          if(chgm>=bmst)                        // (K+/p/Pi+)
          {
            prtM=(-strm)*mKChg+brnm*mProt+(chgm-bmst)*mPiCh;
            prtQC=(-strm)*KPlsQC+brnm*ProtQC+(chgm-bmst)*PiPlQC;
          }
          else if(chgm>=-strm)                  // (K+/p/n)
          {
            prtM=(-strm)*mKChg+chgm*mProt+(brnm-chgm)*mNeut;
            prtQC=(-strm)*KPlsQC+chgm*ProtQC+(brnm-chgm)*NeutQC;
          }
          else if(chgm>=0)                      // (K+/aK0/n)
          {
            prtM=chgm*mKChg+(-chgm-strm)*mKZer+brnm*mNeut;
            prtQC=chgm*KPlsQC+(-chgm-strm)*AKZrQC+brnm*NeutQC;
          }
          else                                  // ch<0 (aK0/n/Pi-)
          {
            prtM=(-strm)*mKChg+(-chgm)*mPiCh+brnm*mNeut;
            prtQC=(-strm)*KPlsQC+(-chgm)*PiMiQC+brnm*NeutQC;
          }
        }
        EnFlowQC-=prtQC;
        chgm=chars-chgm;                        // Just to keep the same notation
        strm=stras-strm;
        brnm=curb;
        if(!strm)                               // --> The total strangness = 0 (n/p/pi-)
        {
          if(chgm<0) resM=(-chgm)*mPiCh+brnm*mNeut;
          else       resM=chgm*mProt+(brnm-chgm)*mNeut;
        }
        else if(strm>=brnm)                     // ---> BigPositiveStrangeness(L/Pi+/K0/K-)
        {
          G4int stmb=strm-brnm;
          if(chgm<0) resM=(-chgm)*mKChg+brnm*mLamb+std::abs(stmb+chgm)*mKZer;
          else       resM=chgm*mPiCh+(strm-brnm)*mKZer+brnm*mLamb;
        }
        else if(strm>0)                         // ---> PositiveStrangeness<B (L/n/p/Pi+-)
        {
          G4int bmst=brnm-strm;
          if     (chgm<0)     resM=(-chgm)*mPiCh+strm*mLamb+bmst*mNeut;
          else if(chgm>=bmst) resM=(chgm-bmst)*mPiCh+strm*mLamb+bmst*mProt;
          else                resM=chgm*mProt+strm*mLamb+(bmst-chgm)*mNeut;
        }
        else                                    // ---> NegativeStrangeness (N/K+/aK0/Pi-)
        {
          G4int bmst=brnm-strm;
          if     (chgm>=bmst)  resM=(-strm)*mKChg+brnm*mProt+(chgm-bmst)*mPiCh;
          else if(chgm>=-strm) resM=(-strm)*mKChg+chgm*mProt+(brnm-chgm)*mNeut;
          else if(chgm>=0)     resM=chgm*mKChg+(-chgm-strm)*mKZer+brnm*mNeut;
          else                 resM=(-strm)*mKChg+(-chgm)*mPiCh+brnm*mNeut;
        }
        G4LorentzVector prt4M=(prtM/(prtM+resM))*EnFlow4M;
        EnFlow4M-=prt4M;
        EnFlowQC-=prtQC;
        G4QHadron* aHadron = new G4QHadron(prtQC, prt4M);
        projHV.push_back(aHadron);
        if(nprt==1)
        {
          G4QHadron* fHadron = new G4QHadron(EnFlowQC, EnFlow4M);
          projHV.push_back(fHadron);
          nprt=0;
        }
#ifdef debug
        G4cout<<"G4StringChipsParticleLevelInterface::Propagate: nprt="<<nprt<<G4endl;
#endif
     } // End of WHILE
  }
  else
  {
     G4QHadron* aHadron = new G4QHadron(EnFlowQC, EnFlow4M);
     projHV.push_back(aHadron);
  }

  // construct the quasmon
  size_t i;
  for (i=0; i<theFinalMomenta.size(); i++) delete theFinalMomenta[i];
  for (i=0; i<theMomenta.size();      i++) delete theMomenta[i];
  theFinalMomenta.clear();
  theMomenta.clear();

  G4QNucleus::SetParameters(fractionOfSingleQuasiFreeNucleons,
                            fractionOfPairedQuasiFreeNucleons,
                            clusteringCoefficient, 
                            fusionToExchange);
  G4Quasmon::SetParameters(temperature, halfTheStrangenessOfSee, etaToEtaPrime);

#ifdef CHIPSdebug
  G4cout<<"G4StringChipsParticleLevelInterface::Propagate: G4QNucleus parameters "
        <<fractionOfSingleQuasiFreeNucleons<<" "<<fractionOfPairedQuasiFreeNucleons
        <<" "<<clusteringCoefficient<<G4endl;
  G4cout<<"G4Quasmon parameters "<<temperature<<" "<<halfTheStrangenessOfSee<<" "
        <<etaToEtaPrime << G4endl;
  G4cout<<"The Target PDG code = "<<targetPDGCode<<G4endl;
  G4cout<<"The projectile momentum = "<<1./MeV*proj4Mom<<G4endl;
  G4cout<<"The target momentum = "<<1./MeV*targ4Mom<<G4endl;
#endif

  // now call chips with this info in place
  G4QHadronVector* output = 0;
  if (particleCount!=0 && resA!=0)
  {
    //  G4QCHIPSWorld aWorld(nop);              // Create CHIPS World of nop particles
    G4QCHIPSWorld::Get()->GetParticles(nop);
    G4QEnvironment* pan= new G4QEnvironment(projHV, targetPDGCode);
#ifdef pdebug
      G4cout<<"G4StringChipsParticleLevelInterface::Propagate: CHIPS fragmentation, rA="
            <<resA<<", #AbsPt="<<particleCount<<G4endl; 
#endif
    try
    {
      output = pan->Fragment();                 // The main fragmentation member function
    }
    catch(G4HadronicException& aR)
    {
      G4cerr << "Exception thrown of G4StringChipsParticleLevelInterface "<<G4endl;
      G4cerr << " targetPDGCode = "<< targetPDGCode <<G4endl;
      G4cerr << " The projectile momentum = "<<1./MeV*proj4Mom<<G4endl;
      G4cerr << " The target momentum = "<<1./MeV*targ4Mom<<G4endl<<G4endl;
      G4cerr << " Dumping the information in the pojectile list"<<G4endl;
      for(i=0; i< projHV.size(); i++)
      {
        G4cerr <<"  Incoming 4-momentum and PDG code of "<<i<<"'th hadron: "
        <<" "<< projHV[i]->Get4Momentum()<<" "<<projHV[i]->GetPDGCode()<<G4endl;
      }
      throw;
    }
    delete pan;
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4StringChipsParticleLevelInterface::Propagate: NO CHIPS fragmentation, rA="
          <<resA<<", #AbsPt="<<particleCount<<G4endl; 
#endif
    output = new G4QHadronVector;
  }   

  // clean up impinging particles
  std::for_each(projHV.begin(), projHV.end(), DeleteQHadron());
  projHV.clear();
   
  // Fill the result.
#ifdef CHIPSdebug
  G4cout << "NEXT EVENT, EscapeExists="<<EscapeExists<<G4endl;
#endif

  // first decay and add all escaping particles.
  if (EscapeExists) for (current = firstEscape; current!=theSorted.end(); current++)
  {
    G4KineticTrack* aResult = (*current).second;
    G4ParticleDefinition* pdef=aResult->GetDefinition();
    secondaries = NULL;
    //if(pdef->GetPDGWidth() > 0 && pdef->GetPDGLifeTime() < 5E-17*s ) // HPW version
    if ( pdef->IsShortLived() )
    {
      secondaries = aResult->Decay();
      for (unsigned int aSecondary=0; aSecondary<secondaries->size(); aSecondary++)
      {
        G4KineticTrack* bResult=secondaries->operator[](aSecondary);
        G4ParticleDefinition* sdef=bResult->GetDefinition();
        if ( sdef->IsShortLived() )
        {
          secsec = bResult->Decay();
          for (unsigned int bSecondary=0; bSecondary<secsec->size(); bSecondary++)
          {
             G4KineticTrack* cResult=secsec->operator[](bSecondary);
             G4ParticleDefinition* cdef=cResult->GetDefinition();
             theSec = new G4ReactionProduct(cdef);
             G4LorentzVector cur4Mom = cResult->Get4Momentum();
             cur4Mom.boost(targ4Mom.boostVector());
             theSec->SetTotalEnergy(cur4Mom.t());
             theSec->SetMomentum(cur4Mom.vect());
#ifdef trapdebug
             if(cdef->GetPDGEncoding()==113) G4cout
                <<"G4StringChipsParticleLevelInterface::Propagate: *Rho0* QGS dec2 PDG="
                <<cdef->GetPDGEncoding()<<",4M="<<cur4Mom<<", grandparPDG= "
                <<pdef->GetPDGEncoding()<<", parPDG= "<<sdef->GetPDGEncoding()<<G4endl; 
#endif
#ifdef pdebug
             G4cout<<"G4StringChipsParticleLevelInterface::Propagate: *OUT* QGS dec2 PDG="
              <<sdef->GetPDGEncoding()<<",4M="<<cur4Mom<<G4endl; 
#endif
             theResult->push_back(theSec);
          }
          std::for_each(secsec->begin(), secsec->end(), DeleteKineticTrack());
          delete secsec;
        }
        else
        {
          theSec = new G4ReactionProduct(sdef);
          G4LorentzVector current4Mom = bResult->Get4Momentum();
          current4Mom.boost(targ4Mom.boostVector());
          theSec->SetTotalEnergy(current4Mom.t());
          theSec->SetMomentum(current4Mom.vect());
#ifdef trapdebug
          if(sdef->GetPDGEncoding()==113)
            G4cout<<"G4StringChipsParticleLevelInterface::Propagate:*Rho0* QGS decay PDG="
                  <<sdef->GetPDGEncoding()<<",4M="<<current4Mom<<", parentPDG= "
                  <<pdef->GetPDGEncoding()<<G4endl; 
        //throw G4HadronicException(__FILE__,__LINE__,
        //                          "G4StringChipsParticleLevelInterface: Rho0 is found!");
#endif
#ifdef pdebug
          G4cout<<"G4StringChipsParticleLevelInterface::Propagate: *OUT* QGS decay PDG="
                <<sdef->GetPDGEncoding()<<",4M="<<current4Mom<<G4endl; 
#endif
          theResult->push_back(theSec);
        }
      }
      std::for_each(secondaries->begin(), secondaries->end(), DeleteKineticTrack());
      delete secondaries;
    }
    else
    {
      theSec = new G4ReactionProduct(aResult->GetDefinition());
      G4LorentzVector current4Mom = aResult->Get4Momentum();
      current4Mom.boost(targ4Mom.boostVector());
      theSec->SetTotalEnergy(current4Mom.t());
      theSec->SetMomentum(current4Mom.vect());
#ifdef trapdebug
      if(aResult->GetDefinition()->GetPDGEncoding()==113)
        G4cout<<"G4StringChipsParticleLevelInterface::Propagate: *OUT* QGS stable PDG="
              <<aResult->GetDefinition()->GetPDGEncoding()<<",4M="<<current4Mom<<G4endl; 
#endif
#ifdef pdebug
      G4cout<<"G4StringChipsParticleLevelInterface::Propagate: *OUT* QGS stable PDG="
            <<aResult->GetDefinition()->GetPDGEncoding()<<",4M="<<current4Mom<<G4endl; 
#endif
      theResult->push_back(theSec);
    } 
  }
  std::for_each(theSecondaries->begin(), theSecondaries->end(), DeleteKineticTrack());
  delete theSecondaries;
    
  // now add the quasmon output
  G4int maxParticle=output->size();
#ifdef CHIPSdebug
  G4cout << "Number of particles from string"<<theResult->size()<<G4endl;
  G4cout << "Number of particles from chips"<<maxParticle<<G4endl;
#endif
#ifdef pdebug
  G4cout << "Number of particles from QGS="<<theResult->size()<<G4endl;
  G4cout << "Number of particles from CHIPS="<<maxParticle<<G4endl;
#endif
  if(maxParticle) for(G4int particle = 0; particle < maxParticle; particle++)
  {
    if(output->operator[](particle)->GetNFragments() != 0) 
    {
      delete output->operator[](particle);
      continue;
    }
    G4int pdgCode = output->operator[](particle)->GetPDGCode();


#ifdef CHIPSdebug
    G4cerr << "PDG code of chips particle = "<<pdgCode<<G4endl;
#endif

    G4ParticleDefinition * theDefinition;
    // Note that I still have to take care of strange nuclei
    // For this I need the mass calculation, and a changed interface
    // for ion-table ==> work for Hisaya @@@@@@@
    // Then I can sort out the pdgCode. I also need a decau process 
    // for strange nuclei; may be another chips interface
    if(pdgCode>90000000) 
    {
      G4int aZ = (pdgCode-90000000)/1000;
      if (aZ>1000) aZ=aZ%1000;  // patch for strange nuclei, to be repaired @@@@
      G4int anN = pdgCode-90000000-1000*aZ;
      if(anN>1000) anN=anN%1000; // patch for strange nuclei, to be repaired @@@@

      if(pdgCode==90000999) theDefinition = G4PionPlus::PionPlusDefinition();
      else if(pdgCode==89999001) theDefinition = G4PionMinus::PionMinusDefinition();
      else if(pdgCode==90999999) theDefinition = G4KaonZero::KaonZeroDefinition();
      else if(pdgCode==90999000) theDefinition = G4KaonMinus::KaonMinusDefinition();
      else if(pdgCode==89001000) theDefinition = G4KaonPlus::KaonPlusDefinition();
      else if(pdgCode==89000001) theDefinition = G4AntiKaonZero::AntiKaonZeroDefinition();
      else if(pdgCode==91000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==92000000) theDefinition = G4Lambda::LambdaDefinition(); //NLambd?
      else if(pdgCode==93000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==94000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==95000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==96000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==97000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==98000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(aZ == 0 && anN == 1) theDefinition = G4Neutron::Neutron();
      else theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,anN+aZ,0,aZ);
    }    
    else theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(pdgCode);

#ifdef CHIPSdebug
    G4cout << "Particle code produced = "<< pdgCode <<G4endl;
#endif

    if(theDefinition)
    {
      theSec = new G4ReactionProduct(theDefinition);
      G4LorentzVector current4Mom = output->operator[](particle)->Get4Momentum();
      current4Mom.boost(targ4Mom.boostVector());
      theSec->SetTotalEnergy(current4Mom.t());
      theSec->SetMomentum(current4Mom.vect());
#ifdef pdebug
      G4cout<<"G4StringChipsParticleLevelInterface::Propagate: *OUT* CHIPS PDG="
              <<theDefinition->GetPDGEncoding()<<",4M="<<current4Mom<<G4endl; 
#endif
      theResult->push_back(theSec);
    }
#ifdef pdebug
    else
    {
      G4cerr <<"G4StringChipsParticleLevelInterface::Propagate: WARNING"<<G4endl;
      G4cerr << "Getting unknown pdgCode from chips in ParticleLevelInterface"<<G4endl;
      G4cerr << "skipping particle with pdgCode = "<<pdgCode<<G4endl<<G4endl;
    }
#endif
    
#ifdef CHIPSdebug
    G4cout <<"CHIPS particles "<<theDefinition->GetPDGCharge()<<" "
           << theDefinition->GetPDGEncoding()<<" "
    << output->operator[](particle)->Get4Momentum() <<G4endl; 
#endif

    delete output->operator[](particle);
  }
  else
  {
    if(resA>0)
    {
      G4ParticleDefinition* theDefinition = G4Neutron::Neutron();
      if(resA==1) // The residual nucleus at rest must be added to conserve BaryN & Charge
      {
        if(resZ == 1) theDefinition = G4Proton::Proton();
      }
      else theDefinition = G4ParticleTable::GetParticleTable()->FindIon(resZ,resA,0,resZ);
      theSec = new G4ReactionProduct(theDefinition);
      theSec->SetTotalEnergy(theDefinition->GetPDGMass());
      theSec->SetMomentum(G4ThreeVector(0.,0.,0.));
      theResult->push_back(theSec);
      if(!resZ && resA>0) for(G4int ni=1; ni<resA; ni++) 
      {
        theSec = new G4ReactionProduct(theDefinition);
        theSec->SetTotalEnergy(theDefinition->GetPDGMass());
        theSec->SetMomentum(G4ThreeVector(0.,0.,0.));
        theResult->push_back(theSec);
      }
    }
  }
  delete output;

#ifdef CHIPSdebug
  G4cout << "Number of particles"<<theResult->size()<<G4endl;
  G4cout << G4endl;
  G4cout << "QUASMON preparation info "
         << 1./MeV*proj4Mom<<" "
  << 1./MeV*targ4Mom<<" "
  << nD<<" "<<nU<<" "<<nS<<" "<<nAD<<" "<<nAU<<" "<<nAS<<" "
  << hitCount<<" "
  << particleCount<<" "
  << theLow<<" "
  << theHigh<<" "
  << G4endl;
#endif

  return theResult;
} 

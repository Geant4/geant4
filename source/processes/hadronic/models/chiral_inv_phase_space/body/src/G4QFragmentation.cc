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
// $Id: G4QFragmentation.cc,v 1.3 2007/05/02 14:59:55 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//                 History: 
//     Created by Mikhail Kossov, October 2006
//     CHIPS QGS fragmentation class 
//     For comparison mirror member functions are taken from G4 classes:
//     G4QGSParticipants
//     G4QGSModels
//     G4ExcitedStringDecay
// -----------------------------------------------------------------------------
//

//#define debug
//#define pdebug
//#define ppdebug

#include "globals.hh"
#include "G4QFragmentation.hh"
#include "G4LorentzVector.hh"
#include <utility>

// Promoting model parameters from local variables class properties @@(? M.K.)

// Definition of static parameters
G4int    G4QFragmentation::nCutMax=7; 
G4double G4QFragmentation::ThersholdParameter=450.*MeV; 
G4double G4QFragmentation::QGSMThershold=3.*GeV; 
G4double G4QFragmentation::theNucleonRadius=1.5*fermi;
 // Parameters of diffractional fragmentation
G4double G4QFragmentation::widthOfPtSquare=-0.72*GeV*GeV; // pt -width2 forStringExcitation
G4double G4QFragmentation::minExtraMass=250.*MeV;// minimum excitation mass 
G4double G4QFragmentation::minmass=250.*MeV;	    // mean pion transverse mass for Xmin 

G4QFragmentation::G4QFragmentation() : theNucleus(0)
{
  // Construct Shortlived particles (needed after the 2006 Particles revolution)
	 G4ShortLivedConstructor ShortLived;
	 ShortLived.ConstructParticle();
}

G4QFragmentation::~G4QFragmentation() {if(theNucleus) delete theNucleus;}

G4QHadronVector* G4QFragmentation::Scatter(const G4QNucleus &aNucleus,
                                           const G4QHadron &aPrimary)
{  
  G4QStringVector* strings=0;

  G4QHadron thePrimary = aPrimary;
  
  G4LorentzRotation toZ;
  G4LorentzVector Ptmp=thePrimary.Get4Momentum();
  toZ.rotateZ(-1*Ptmp.phi());
  toZ.rotateY(-1*Ptmp.theta());
  thePrimary.Set4Momentum(toZ*Ptmp);
  G4LorentzRotation toLab(toZ.inverse());

  G4int attempts = 0, maxAttempts=20;
  while(!strings)
  {
  	 if (attempts++ > maxAttempts ) 
  	 {
				 	G4cerr<<"***G4QFragmentation::Scatter: "<<attempts<<" to create a string ( > max="
            <<maxAttempts<<") --> try to increase maxAttempts"<<G4endl;
      G4Exception("G4QFragmentation::Scatter:","72",FatalException,"StringCreation");
  	 }
	   InitModel(aNucleus,thePrimary);
  	 strings = GetStrings();
  }
  
  G4QHadronVector* theResult = 0;
  G4double stringEnergy(0);
  for( unsigned astring=0; astring < strings->size(); astring++)
  {
    //    rotate string to lab frame, models have it aligned to z
    stringEnergy += (*strings)[astring]->GetLeftParton()->Get4Momentum().t();
    stringEnergy += (*strings)[astring]->GetRightParton()->Get4Momentum().t();
    (*strings)[astring]->LorentzRotate(toLab);
  }
#ifdef debug
  G4cout<<"G4QFragmentation::Scatter: Total string energy = "<<stringEnergy<<G4endl;
#endif
  theResult = FragmentStrings(strings);
  std::for_each(strings->begin(), strings->end(), DeleteQString() );
  delete strings;

  return theResult;
}

void G4QFragmentation::InitModel(const G4QNucleus& aNucleus,const G4QHadron& aProjectile)
{
  static const G4double 	mProt = G4Proton::Proton()->GetPDGMass(); // Mass of proton
  Init(aNucleus.GetN(),aNucleus.GetZ());
  theCurrentVelocity.setX(0);    
  theCurrentVelocity.setY(0); 
  // @@ "target Nucleon" == "Proton at rest" (M.K. ?)
  G4double nCons = 1;                         // 1 or baryon number of the projectile
  G4int projAbsB=std::abs(aProjectile.GetBaryonNumber()); // Baryon/anti-baryon
  if(projAbsB>1) nCons = projAbsB;
  G4LorentzVector proj4M = aProjectile.Get4Momentum();
  G4double pz_per_projectile = proj4M.pz()/nCons;
  G4double e_per_projectile = proj4M.e()/nCons+mProt; // @@ Add mass of TargetProtonAtRest
  //G4double e_per_projectile = aProjectile.Get4Momentum()*aProjectile.Get4Momentum();
  //e_per_projectile += proj4M.vect()*proj4M.vect();
  //e_per_projectile /=nCons*nCons;
	 //e_per_projectile = std::sqrt(e_per_projectile);
	 //e_per_projectile += G4Proton::Proton()->GetPDGMass();//@@Add mass of TargetProtonAtRest
  G4double vz = pz_per_projectile/e_per_projectile;
#ifdef debug
		G4cout<<"G4QFragmentation::Init: Projectile4M="<<proj4M<<", vz="<<vz<<G4endl;
#endif
  theCurrentVelocity.setZ(vz);
  DoLorentzBoost(-theCurrentVelocity); 
  G4LorentzVector Mom = proj4M;
  Mom.boost(-theCurrentVelocity);
  G4QHadron theProjectile(aProjectile.GetQPDG(),Mom);
#ifdef debug
  G4cout<<"G4QFragmentation::Init: PreInteractionMomentum"<<Mom<<G4endl;
#endif
  BuildInteractions(theProjectile);
  GetWoundedNucleus()->DoLorentzBoost(theCurrentVelocity);
}

G4QStringVector* G4QFragmentation::GetStrings()
{
  G4QPartonPair* aPair;
  G4QStringVector* theStrings = new G4QStringVector;
  G4QString* aString=0;
  while((aPair = GetNextPartonPair())) // @@ At present no difference in stringBuild ? M.K.
  {
    if (aPair->GetCollisionType() == G4QPartonPair::DIFFRACTIVE)
    {
      aString = BuildString(aPair);           // @@ ?
#ifdef debug
      G4cout<<"G4QFragmentation::GetStrings:DifString4M="<<aString->Get4Momentum()<<G4endl;
#endif
    }
    else
    {
      aString = BuildString(aPair);           // @@ ?
#ifdef debug
      G4cout<<"G4QFragmentation::GetStrings:SftString4M="<<aString->Get4Momentum()<<G4endl;
#endif
    }
    aString->Boost(theCurrentVelocity);  
    theStrings->push_back(aString);
    delete aPair;
  }
#ifdef debug
  for(G4int i=0; i<theStrings->size(); i++) G4cout<<"G4QFragmentation::GetStrings:String #"
                                    <<i<<", 4M="<<(*theStrings)[i]->Get4Momentum()<<G4endl;
#endif
  return theStrings;
}

void G4QFragmentation::BuildInteractions(const G4QHadron &thePrimary) 
{
  
  // Find the collisions and collition conditions
  G4QHadron* aProjectile = SelectInteractions(thePrimary);

  // now build the parton pairs
  SplitHadrons();
  
  // soft collisions, ordering is vital
  PerformSoftCollisions();
   
  // the rest is diffractive
  PerformDiffractiveCollisions();
  
  // clean-up, if necessary
  std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
  theInteractions.clear();
  std::for_each(theTargets.begin(), theTargets.end(), DeleteQHadron());
  theTargets.clear();
  delete aProjectile;
}

//
G4QHadron* G4QFragmentation::SelectInteractions(const G4QHadron &thePrimary)
{
  G4QHadron* aProjectile = new G4QHadron(thePrimary);
  G4QPomeron theProbability(thePrimary.GetPDGCode()); // must be data member
  G4double outerRadius = theNucleus->GetOuterRadius();

  // Check reaction threshold 
  theNucleus->StartLoop();
  G4QHadron* pNucleon = theNucleus->GetNextNucleon();
  G4LorentzVector aPrimaryMomentum=thePrimary.Get4Momentum();
  G4double s = (aPrimaryMomentum + pNucleon->Get4Momentum()).mag2();
  G4double primaryMass=thePrimary.GetMass();
  G4double ThresholdMass = primaryMass + pNucleon->GetMass(); 
  ModelMode = SOFT;
  if (sqr(ThresholdMass + ThersholdParameter) > s)
  {
    G4cerr<<"***G4QFragmentation::SelectInteractions: ThrM="<<ThresholdMass<<" + ThrPa="
          <<ThersholdParameter<<" = "<<ThresholdMass+ThersholdParameter<<" > std::sqrt(s)="
          <<std::sqrt(s)<<G4endl;
    G4Exception("G4QFragmentation::SelectInteractions:","72",FatalException,"LowEnergy");
  }
  if (sqr(ThresholdMass + QGSMThershold) > s) // thus only diffractive in cascade!
  {
#ifdef debug
    G4cout<<"G4QFragmentation::SelectInteractions: ThrM="<<ThresholdMass<<" + ThrQGS="
          <<QGSMThershold<<" = "<<ThresholdMass+QGSMThershold<<" > std::sqrt(s)="
          <<std::sqrt(s)<<" -> only Diffraction is possible"<<G4endl; // @@ to Quasmon
#endif
    ModelMode = DIFFRACTIVE;
  }
 
  // first find the collisions HPW
  std::for_each(theInteractions.begin(),theInteractions.end(), DeleteQInteraction());
  theInteractions.clear();
  G4int totalCuts = 0;
  G4double impactUsed = 0;

  G4double eKin = aPrimaryMomentum.e()-primaryMass; // Primary kinetic energy ? GeV ? M.K.

  while(theInteractions.size() == 0)
  {
    // choose random impact parameter HPW
    std::pair<G4double, G4double> theImpactParameter;
    theImpactParameter = theNucleus->ChooseImpactXandY(outerRadius+theNucleonRadius);
    G4double impactX = theImpactParameter.first; 
    G4double impactY = theImpactParameter.second;
    
    // loop over nuclei to find collissions HPW
    theNucleus->StartLoop();
    G4int nucleonCount = 0; // debug
    //G4QFragmentation_NPart = 0;                        // ? M.K.
    while( (pNucleon = theNucleus->GetNextNucleon()) )
    {
      if(totalCuts>1.5*eKin) break;
      nucleonCount++; // debug
      // Needs to be moved to Probability class @@@
      G4double s = (aPrimaryMomentum + pNucleon->Get4Momentum()).mag2();
      G4double Distance2 = sqr(impactX - pNucleon->GetPosition().x()) +
                           sqr(impactY - pNucleon->GetPosition().y());
      G4double Probability = theProbability.GetInelasticProbability(s, Distance2);  
      // test for inelastic collision
      G4double rndNumber = G4UniformRand();
      // ModelMode = DIFFRACTIVE;
      if (Probability > rndNumber)
      {
#ifdef debug
        G4cout<<"G4QFragmentation::SelectInteractions: p="<<Probability<<", r="<<rndNumber
              <<", d="<<std::sqrt(Distance2)<<G4endl;
#endif
        G4QHadron* aTarget = new G4QHadron(*pNucleon);
        //G4QFragmentation_NPart ++;                       // ? M.K.
	       theTargets.push_back(aTarget);
 	      pNucleon=aTarget;
        if((theProbability.GetDiffractiveProbability(s,Distance2)/Probability >
            G4UniformRand() && (ModelMode==SOFT)) || ModelMode==DIFFRACTIVE)
	       { 
	         // diffractive interaction occurs
	         if(IsSingleDiffractive()) ExciteSingDiffParticipants(aProjectile, aTarget);
	         else                          ExciteDiffParticipants(aProjectile, aTarget);
          G4QInteraction* aInteraction = new G4QInteraction(aProjectile);
          aInteraction->SetTarget(aTarget); 
          theInteractions.push_back(aInteraction);
	         aInteraction->SetNumberOfDiffractiveCollisions(1);
          totalCuts += 1;
	       }
	       else
	       {
	         // nondiffractive soft interaction occurs
	         // sample nCut+1 (cut Pomerons) pairs of strings can be produced
          G4int nCut;
          G4double* running = new G4double[nCutMax];
          running[0] = 0;
          for(nCut = 0; nCut < nCutMax; nCut++)
          {
	           running[nCut] = theProbability.GetCutPomeronProbability(s, Distance2, nCut+1);
            if(nCut!=0) running[nCut] += running[nCut-1];
          }
	         G4double random = running[nCutMax-1]*G4UniformRand();
          for(nCut = 0; nCut < nCutMax; nCut++) {if(running[nCut] > random) break;}
          delete [] running;
          nCut = 0;
          aTarget->IncrementCollisionCount(nCut+1);
          aProjectile->IncrementCollisionCount(nCut+1);
          G4QInteraction* aInteraction = new G4QInteraction(aProjectile);
          aInteraction->SetTarget(aTarget);
          aInteraction->SetNumberOfSoftCollisions(nCut+1);
          theInteractions.push_back(aInteraction);
          totalCuts += nCut+1;
          impactUsed=Distance2;
        }
      }
    }
#ifdef debug
    G4cout<<"G4QFragmentation::SelectInteractions: NUCLEONCOUNT="<<nucleonCount<<G4endl;
#endif
  }
#ifdef debug
  G4cout<<"G4QFragmentation::SelectInteractions: CUTDEBUG="<<totalCuts
        <<", ImpactParameter="<<impactUsed<<G4endl;
#endif
  return aProjectile;
}

void G4QFragmentation::PerformDiffractiveCollisions()
{
  for(unsigned i = 0; i < theInteractions.size(); i++) 
  {
    G4QInteraction* anIniteraction = theInteractions[i];
    G4QHadron* aProjectile = anIniteraction->GetProjectile();
    G4QParton* aParton = aProjectile->GetNextParton();
    G4QPartonPair* aPartonPair;
    // projectile first
    if (aParton)
    {
      aPartonPair = new G4QPartonPair(aParton, aProjectile->GetNextAntiParton(), 
                                      G4QPartonPair::DIFFRACTIVE, 
                                      G4QPartonPair::PROJECTILE);
      thePartonPairs.push_back(aPartonPair);
    }
    // then target
    G4QHadron* aTarget = anIniteraction->GetTarget();
    aParton = aTarget->GetNextParton();
    if (aParton)
    {
      aPartonPair = new G4QPartonPair(aParton, aTarget->GetNextAntiParton(), 
                                      G4QPartonPair::DIFFRACTIVE, 
                                      G4QPartonPair::TARGET);
      thePartonPairs.push_back(aPartonPair);
    }
  }
}

void G4QFragmentation::PerformSoftCollisions()
{
  G4QInteractionVector::iterator i;
  for(i = theInteractions.begin(); i != theInteractions.end(); i++)   
  {
    G4QInteraction* anIniteraction = *i;
    G4QPartonPair* aPair=0;
    if (anIniteraction->GetNumberOfSoftCollisions())
    { 
      G4QHadron* pProjectile = anIniteraction->GetProjectile();
      G4QHadron* pTarget     = anIniteraction->GetTarget();
      for (G4int j = 0; j < anIniteraction->GetNumberOfSoftCollisions(); j++)
      {
        aPair= new G4QPartonPair(pTarget->GetNextParton(),pProjectile->GetNextAntiParton(),
                                 G4QPartonPair::SOFT, G4QPartonPair::TARGET);
        thePartonPairs.push_back(aPair);
        aPair= new G4QPartonPair(pProjectile->GetNextParton(),pTarget->GetNextAntiParton(),
                                G4QPartonPair::SOFT, G4QPartonPair::PROJECTILE);
        thePartonPairs.push_back(aPair);
      }  
      delete *i;
      i=theInteractions.erase(i);
      i--;
    }
  }
}

G4QPartonPair* G4QFragmentation::GetNextPartonPair()
{
  if(thePartonPairs.empty()) return 0;
  G4QPartonPair* result = thePartonPairs.back();
  thePartonPairs.pop_back();
  return result;
}

G4QHadronVector* G4QFragmentation::FragmentStrings(const G4QStringVector* theStrings)
{
  G4QHadronVector* theResult = new G4QHadronVector;

  G4LorentzVector KTsum(0.,0.,0.);
  G4bool NeedEnergyCorrector=false;
  
  for( unsigned astring=0; astring < theStrings->size(); astring++)
  {
    KTsum+= (*theStrings)[astring]->Get4Momentum();
    if(!(KTsum.e()<1.) && !(KTsum.e()>-1.))
    {
      G4cerr<<"***G4QFragmentation::FragmentStrings: KTsum="<<KTsum<<G4endl;
      G4Exception("G4QFragmentation::FragmentStrings:","72",FatalException,"NANin3Vector");
    }
    G4QHadronVector* generatedHadrons = 0;     // A prototype of the string output
    if( (*theStrings)[astring]->IsExcited() )
    {
      generatedHadrons=(*theStrings)[astring]->FragmentString(true);// Fragment QGSM String
    }
    else
    {
	     //generatedHadrons = new G4QHadronVector;
	     //generatedHadrons->push_back((*theStrings)[astring]->GetAsQHadron()); //@@ NotImplem
    }    
	   if (!generatedHadrons) 
	   {
		    G4cerr<<"G4QFragmentation::FragmentStrings: No Hadrons produced" << G4endl;
		    continue;
	   }
    G4LorentzVector KTsum1(0.,0.,0.,0.);
    for(unsigned aTrack=0; aTrack<generatedHadrons->size(); aTrack++)
	   {
		    theResult->push_back((*generatedHadrons)[aTrack]);
		    KTsum1+= (*generatedHadrons)[aTrack]->Get4Momentum();
	   }
	
	   if(std::abs((KTsum1.e()-(*theStrings)[astring]->Get4Momentum().e())/KTsum1.e())
       > perMillion) NeedEnergyCorrector=true;
    //      clean up
	   delete generatedHadrons;
  }
#ifdef debug
  G4cout<<"G4QFragmentation::FragmentStrings: String4mom="<<KTsum<<endl; 
#endif
  if(NeedEnergyCorrector) EnergyAndMomentumCorrector(theResult, KTsum);
  return theResult;
}

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
  if(!success)
  {
    G4cout<<"***G4QFragmentation::EnergyAndMomentumCorrector: Scale #1 at end of loop M="
          <<TotalCollisionMass<<", S"<<Sum<<", Sc="<<Scale
          <<" Increase number of attempts or increase ERRLIMIT"<<G4endl;
    G4Exception("G4QFragmentation::SelectInteractions:","72",FatalException,"NotCorrect");
  }
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
	 	   return false; 	  //  Ignore this interaction 
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
	 	   return false; 	  //  Ignore this interaction 
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
  ThersholdParameter = thresh;         // internal threshold
  QGSMThershold      = QGSMth;         // QGSM threshold
  theNucleonRadius   = radNuc;         // effective radius of the nucleon inside Nucleus
	 widthOfPtSquare    = -2*SigPt*SigPt;	// width^2 of pt for string excitation
	 minExtraMass       = extraM;	        // minimum excitation mass 
	 minmass            = minM;	          // mean pion transverse mass; used for Xmin 
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

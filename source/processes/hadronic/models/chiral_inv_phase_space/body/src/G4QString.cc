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
// $Id: G4QString.cc,v 1.3 2006/12/12 11:02:22 mkossov Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QString ----------------
//      by Mikhail Kossov Oct, 2006
//      class for excited string used by Parton String Models
//   For comparison mirror member functions are taken from G4 classes:
//   G4FragmentingString
//   G4ExcitedStringDecay
// ------------------------------------------------------------

//#define debug

#include "G4QString.hh"
#include <algorithm>
// Static parameters definition
G4double G4QString::MassCut=350.*MeV;     // minimum mass cut for the string
G4double G4QString::ClusterMass=150.*MeV; // minimum cluster mass for fragmentation
G4double G4QString::SigmaQT=0.5*GeV;      // quarkTransverseMomentum distribution parameter
G4double G4QString::DiquarkSuppress=0.1;  // is Diquark suppression parameter  
G4double G4QString::DiquarkBreakProb=0.1; // is Diquark breaking probability 
G4double G4QString::SmoothParam=0.9;      // QGS model parameter
G4double G4QString::StrangeSuppress=0.435;// Strangeness suppression (u:d:s=1:1:0.3 ?M.K.)
G4double G4QString::widthOfPtSquare=-0.72*GeV*GeV; // pt -width2 forStringExcitation
G4int G4QString::StringLoopInterrupt=999; // String fragmentation LOOP limit 
G4int G4QString::ClusterLoopInterrupt=500;// Cluster fragmentation LOOP limit 

G4QString::G4QString() : theDirection(0),thePosition(G4ThreeVector(0.,0.,0.)),
 hadronizer(new G4QHadronBuilder){}

G4QString::G4QString(G4QParton* Color, G4QParton* AntiColor, G4int Direction) :
  hadronizer(new G4QHadronBuilder)
{
  ExciteString(Color, AntiColor, Direction);
}

G4QString::G4QString(G4QParton* QCol,G4QParton* Gluon,G4QParton* QAntiCol,G4int Direction):
	 theDirection(Direction),thePosition(QCol->GetPosition()),hadronizer(new G4QHadronBuilder)
{
  thePartons.push_back(QCol);
  thePartons.push_back(Gluon);
  thePartons.push_back(QAntiCol);
}

G4QString::G4QString(const G4QString &right) : theDirection(right.GetDirection()),
thePosition(right.GetPosition()), hadronizer(new G4QHadronBuilder){}

G4QString::~G4QString()
 {if(thePartons.size())std::for_each(thePartons.begin(),thePartons.end(),DeleteQParton());}

G4LorentzVector G4QString::Get4Momentum() const
{
	 G4LorentzVector momentum(0.,0.,0.,0.);
	 for(unsigned i=0; i<thePartons.size(); i++) momentum += thePartons[i]->Get4Momentum();
	 return momentum;
}

void G4QString::LorentzRotate(const G4LorentzRotation & rotation)
{
	 for(unsigned i=0; i<thePartons.size(); i++)
	    thePartons[i]->Set4Momentum(rotation*thePartons[i]->Get4Momentum());
}

void G4QString::InsertParton(G4QParton* aParton, const G4QParton* addafter)
{
  G4QPartonVector::iterator insert_index;   // Begin by default (? M.K.)
	
	 if(addafter != NULL) 
	 {
	   insert_index=std::find(thePartons.begin(), thePartons.end(), addafter);
	   if (insert_index == thePartons.end())		// No object addafter in thePartons
	   {
				 	G4cerr<<"***G4QString::InsertParton: Address Parton is not found"<<G4endl;
      G4Exception("G4QString::InsertParton:","72",FatalException,"NoAddressParton");
	   }
	 }
	 thePartons.insert(insert_index+1, aParton);
} 

G4LorentzRotation G4QString::TransformToCenterOfMass()
{
	 G4LorentzRotation toCms(-1*Get4Momentum().boostVector());
	 for(unsigned i=0; i<thePartons.size(); i++)
	                     thePartons[i]->Set4Momentum(toCms*thePartons[i]->Get4Momentum());
	 return toCms;
}

G4LorentzRotation G4QString::TransformToAlignedCms()
{
	 G4LorentzVector momentum=Get4Momentum();
	 G4LorentzRotation toAlignedCms(-1*momentum.boostVector());
	 momentum= toAlignedCms* thePartons[0]->Get4Momentum();
	 toAlignedCms.rotateZ(-1*momentum.phi());
	 toAlignedCms.rotateY(-1*momentum.theta());
	 for(unsigned index=0; index<thePartons.size(); index++)
	 {
	   momentum=toAlignedCms * thePartons[index]->Get4Momentum();
	   thePartons[index]->Set4Momentum(momentum);
	 }
	 return toAlignedCms;
}

void G4QString::Boost(G4ThreeVector& Velocity)
{
  for(unsigned cParton=0; cParton<thePartons.size(); cParton++)
  {
    G4LorentzVector Mom = thePartons[cParton]->Get4Momentum();
    Mom.boost(Velocity);
    thePartons[cParton]->Set4Momentum(Mom);
  }
}

G4QParton* G4QString::GetColorParton(void) const
{
  G4QParton* start = *(thePartons.begin());
  G4QParton* end = *(thePartons.end()-1);
  G4int Encoding = start->GetPDGCode();
  if (Encoding<-1000 || (Encoding  < 1000 && Encoding > 0)) return start;
  return end; 
}

G4QParton* G4QString::GetAntiColorParton(void) const
{
  G4QParton* start = *(thePartons.begin());
  G4QParton* end = *(thePartons.end()-1);
  G4int Encoding = start->GetPDGCode();
  if (Encoding < -1000 || (Encoding  < 1000 && Encoding > 0)) return end; 
  return start; 
}

// Fill parameters
void G4QString::SetParameters(G4double mCut,G4double clustM, G4double sigQT,G4double DQSup,
    G4double DQBU, G4double smPar, G4double SSup, G4double SigPt, G4int SLmax, G4int CLmax)
{//  =============================================================================
  MassCut         = mCut;           // minimum mass cut for the string
  ClusterMass     = clustM;         // minimum cluster mass for the fragmentation
  SigmaQT         = sigQT;          // quark transverse momentum distribution parameter 
  DiquarkSuppress = DQSup;          // is Diquark suppression parameter  
  DiquarkBreakProb= DQBU;           // is Diquark breaking probability 
  SmoothParam     = smPar;          // QGS model parameter
  StrangeSuppress = SSup;           // Strangeness suppression parameter
	 widthOfPtSquare = -2*SigPt*SigPt;	// width^2 of pt for string excitation
  StringLoopInterrupt = SLmax;      // String fragmentation LOOP limit 
  ClusterLoopInterrupt= CLmax;      // Cluster fragmentation LOOP limit 
}

// Pt distribution @@ one can use 1/(1+A*Pt^2)^B
G4ThreeVector G4QString::GaussianPt(G4double widthSquare, G4double maxPtSquare) const
{
	 G4double pt2; do{pt2=widthSquare*std::log(G4UniformRand());} while (pt2>maxPtSquare);
		pt2=std::sqrt(pt2);
	 G4double phi=G4UniformRand()*twopi;
	 return G4ThreeVector(pt2*std::cos(phi),pt2*std::sin(phi),0.);    
} // End of GaussianPt

// Diffractively excite the string
void G4QString::DiffString(G4QHadron* hadron, G4bool isProjectile)
{
	 hadron->SplitUp();
	 G4QParton* start = hadron->GetNextParton();
	 if( start==NULL) 
	 {
    G4cerr<<"***G4QString::DiffString: No start parton found, nothing is done"<<G4endl;
	   return;
 	}
	 G4QParton* end = hadron->GetNextParton();
	 if( end==NULL) 
	 {
    G4cerr<<"***G4QString::DiffString: No end parton found, nothing is done"<<G4endl;
	   return;
	 }
	 if(isProjectile)	ExciteString(end, start, 1); //  1 = Projectile
	 else            	ExciteString(start, end,-1); // -1 = Target
	 SetPosition(hadron->GetPosition());
  // momenta of string ends
	 G4double ptSquared= hadron->Get4Momentum().perp2();
  G4double hmins=hadron->Get4Momentum().minus();
  G4double hplus=hadron->Get4Momentum().plus();
	 G4double transMassSquared=hplus*hmins;
 	G4double maxMomentum = std::sqrt(transMassSquared) - std::sqrt(ptSquared);
	 G4double maxAvailMomentumSquared = maxMomentum * maxMomentum;
 	G4ThreeVector pt=GaussianPt(widthOfPtSquare,maxAvailMomentumSquared);
 	G4LorentzVector Pstart(G4LorentzVector(pt,0.));
	 G4LorentzVector Pend(hadron->Get4Momentum().px(), hadron->Get4Momentum().py(), 0.);
	 Pend-=Pstart;
	 G4double tm1=hmins+(Pend.perp2()-Pstart.perp2())/hplus;
	 G4double tm2=std::sqrt( std::max(0., tm1*tm1-4*Pend.perp2()*hmins/hplus ) );
 	G4int Sign= isProjectile ? TARGET : PROJECTILE;
 	G4double endMinus  = 0.5 * (tm1 + Sign*tm2);
	 G4double startMinus= hmins - endMinus;
 	G4double startPlus = Pstart.perp2() / startMinus;
	 G4double endPlus   = hplus - startPlus;
	 Pstart.setPz(0.5*(startPlus - startMinus));
	 Pstart.setE (0.5*(startPlus + startMinus));
	 Pend.setPz  (0.5*(endPlus - endMinus));
	 Pend.setE   (0.5*(endPlus + endMinus));
	 start->Set4Momentum(Pstart);
	 end->Set4Momentum(Pend);
#ifdef debug
		G4cout<<"G4QString::DiffString: StartQ="<<start->GetPDGCode()<<start->Get4Momentum()<<"("
        <<start->Get4Momentum().mag()<<"), EndQ="<<end->GetPDGCode()<<end ->Get4Momentum()
        <<"("<<end->Get4Momentum().mag()<<"), sumOfEnds="<<Pstart+Pend<<", H4M(original)="
        <<hadron->Get4Momentum()<<G4endl;
#endif
} // End of DiffString (The string is excited diffractively)

// Excite the string by two partons
void G4QString::ExciteString(G4QParton* Color, G4QParton* AntiColor, G4int Direction)
{
		theDirection = Direction;
  thePosition  = Color->GetPosition();
  thePartons.push_back(Color);
  thePartons.push_back(AntiColor);
} // End of ExciteString

// LUND Longitudinal fragmentation
G4double G4QString::GetLundLightConeZ(G4double zmin, G4double zmax, G4int , // @@ ? M.K.
                                      G4QHadron* pHadron, G4double Px, G4double Py)
{
    static const G4double  alund = 0.7/GeV/GeV; 
    // If blund get restored, you MUST adapt the calculation of zOfMaxyf.
    //static const G4double  blund = 1;
    G4double z, yf;
    G4double Mt2 = Px*Px + Py*Py + pHadron->GetMass2();
    G4double zOfMaxyf=alund*Mt2/(alund*Mt2+1.);
    G4double maxYf=(1.-zOfMaxyf)/zOfMaxyf * std::exp(-alund*Mt2/zOfMaxyf);
    do
    {
       z = zmin + G4UniformRand()*(zmax-zmin);
       // yf = std::pow(1. - z, blund)/z*std::exp(-alund*Mt2/z);
	      yf = (1-z)/z * std::exp(-alund*Mt2/z);
    } while (G4UniformRand()*maxYf>yf); 
    return z;
} // End of GetLundLightConeZ


// QGSM Longitudinal fragmentation
G4double G4QString::GetQGSMLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,
																																						G4QHadron* , G4double, G4double) // @@ M.K.
{
  static const G4double arho=0.5;
  static const G4double aphi=0.;
  static const G4double an=-0.5;
  static const G4double ala=-0.75;
  static const G4double aksi=-1.;
  static const G4double alft=0.5;
  G4double z;    
  G4double theA(0), d2, yf;
  G4int absCode = std::abs(PartonEncoding);
  // @@ Crazy algorithm is used for simple power low...
  if (absCode < 10)                            // Quarks (@@ 9 can be a gluon)
  { 
    if(absCode == 1 || absCode == 2) theA = arho;
    else if(absCode == 3)            theA = aphi;
    else
    {
     G4cerr<<"***G4QString::GetQGSMLightConeZ: CHIPS is SU(3), quakCode="<<absCode<<G4endl;
     G4Exception("G4QString::GetQGSMLightConeZ:","72",FatalException,"WrongQuark");
    }
    do 	
    {
      z  = zmin + G4UniformRand()*(zmax - zmin);
      d2 = alft - theA;
      yf = std::pow( 1.-z, d2);
    } 
    while (G4UniformRand()>yf);
  }
  else                                     // Di-quarks (@@ Crazy codes are not checked)
  {       
    if (absCode==3101||absCode==3103||absCode==3201||absCode==3203)   d2=alft-ala-ala+arho;
    else if(absCode==1103||absCode==2101||absCode==2203||absCode==2103) d2=alft-an-an+arho;
    else                                                            d2=alft-aksi-aksi+arho;
    do  
    {
      z = zmin + G4UniformRand()*(zmax - zmin);
      yf = std::pow(1.-z, d2);
    } 
    while (G4UniformRand()>yf); 
  }
  return z;
} // End of GetQGSMLightConeZ

// Diffractively excite the string
G4QHadronVector* G4QString::FragmentString(G4bool QL)
{
  // Can no longer modify Parameters for Fragmentation.
	 // PastInitPhase=true;                                     // Now static
	
  // 	check if string has enough mass to fragment...
	 G4QHadronVector* LeftVector=LightFragmentationTest();
	 if(LeftVector) return LeftVector;
	
	 LeftVector = new G4QHadronVector;
	 G4QHadronVector* RightVector = new G4QHadronVector;

  // this should work but it's only a semi deep copy.
  // %GF	G4ExcitedString theStringInCMS(theString);
  G4QString* theStringInCMS=CPExcited();                      // must be deleted
	 G4LorentzRotation toCms=theStringInCMS->TransformToAlignedCms();
	 G4bool success=false;
	 G4bool inner_sucess=true;
	 G4int attempt=0;
	 while (!success && attempt++<StringLoopInterrupt)              //@@ It's Loop with break
	 {
		  G4QString* currentString = new G4QString(*theStringInCMS);
		  if(LeftVector->size())
    {
      std::for_each(LeftVector->begin(), LeftVector->end(), DeleteQHadron());
		    LeftVector->clear();
    }
		  if(RightVector->size())
		  {
      std::for_each(RightVector->begin(), RightVector->end(), DeleteQHadron());
		    RightVector->clear();
		  }
		  inner_sucess=true;  // set false on failure..
		  while (!StopFragmentation())
		  {  // Split current string into hadron + new string
		   	G4QString* newString=0;  // used as output from SplitUp...
			   G4QHadron* Hadron=Splitup(QL);
			   if(Hadron && IsFragmentable())
			   {
			     if(currentString->GetDecayDirection()>0) LeftVector->push_back(Hadron);
       	else RightVector->push_back(Hadron);
			     delete currentString;
			     currentString=newString;
			   }
      else
      {
			     // abandon ... start from the beginning
			     if (newString) delete newString;
			     if (Hadron)    delete Hadron;
			     inner_sucess=false;
			     break;
			   }
		  } 
		// Split current string into 2 final Hadrons
		  if( inner_sucess && SplitLast(currentString, LeftVector, RightVector) ) 	success=true;
		  delete currentString;
	 }
	 delete theStringInCMS;
	 if (!success)
	 {
		  if(RightVector->size())
		  {
      std::for_each(RightVector->begin(), RightVector->end(), DeleteQHadron());
		    RightVector->clear();
		  }
		  delete RightVector;
		  if(LeftVector->size())
    {
      std::for_each(LeftVector->begin(), LeftVector->end(), DeleteQHadron());
		    LeftVector->clear();
    }
		  return LeftVector;
	 }		
	 // Join Left and Right Vectors into LeftVector in correct order.
	 while(!RightVector->empty())
	 {
	    LeftVector->push_back(RightVector->back());
	    RightVector->erase(RightVector->end()-1);
	 }
	 delete RightVector;
	 CalculateHadronTimePosition(Get4Momentum().mag(), LeftVector);
	 G4LorentzRotation toObserverFrame(toCms.inverse());
	 for(unsigned C1 = 0; C1 < LeftVector->size(); C1++)
	 {
	   G4QHadron* Hadron = (*LeftVector)[C1];
	   G4LorentzVector Momentum = Hadron->Get4Momentum();
	   Momentum = toObserverFrame*Momentum;
	   Hadron->Set4Momentum(Momentum);
	   G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
	   Momentum = toObserverFrame*Coordinate;
	   Hadron->SetFormationTime(Momentum.e());
	   Hadron->SetPosition(GetPosition()+Momentum.vect());
	 }
	 return LeftVector;
} // End of FragmentLundString

// Creates a string, using only the end-partons of the string
G4QString* G4QString::CPExcited()
{
	G4QParton* Left = new G4QParton(GetLeftParton());
	G4QParton* Right= new G4QParton(GetRightParton());
	return new G4QString(Left,Right,GetDirection());
} // End of CPExcited

// Simple decay of the string
G4QHadronVector* G4QString::LightFragmentationTest()
{
  static const G4double MassCut = 0.35*GeV; 
  // Check string decay threshold
		
	 G4QHadronVector* result=0;  // return 0 when string exceeds the mass cut
	
	 G4QHadronPair hadrons((G4QHadron*)0, (G4QHadron*)0); // pair of hadrons
	 if(sqr(FragmentationMass(0,&hadrons)+MassCut)<Mass2()) return 0; //Par MassCut
	
	 result = new G4QHadronVector;
        
	 if(hadrons.second==0)                   // Second hadron exists
	 {
	   // Substitute string by a light hadron, Note that Energy is not conserved here! @@
	   G4ThreeVector Mom3 = Get4Momentum().vect();
	   G4LorentzVector Mom(Mom3, std::sqrt(Mom3.mag2() + hadrons.first->GetMass2()) );
    result->push_back(new G4QHadron(hadrons.first, 0, GetPosition(), Mom));
	 }
  else 
	 {
	   //... string was qq--qqbar type: Build two stable hadrons,
	   G4LorentzVector  Mom1, Mom2;
	   Sample4Momentum(&Mom1, hadrons.first->GetMass(), &Mom2, hadrons.second->GetMass(),
			                                                          Get4Momentum().mag());
	   result->push_back(new G4QHadron(hadrons.first, 0, GetPosition(), Mom1));
	   result->push_back(new G4QHadron(hadrons.second,0, GetPosition(), Mom2));
    G4ThreeVector Velocity = Get4Momentum().boostVector();
    G4int L=result->size(); if(L) for(G4int i=0; i<L; i++) (*result)[i]->Boost(Velocity);
	 }
	 return result;
} // End of LightFragmentationTest

// Calculate Fragmentation Mass (if pdefs # 0, returns two hadrons)
G4double G4QString::FragmentationMass(G4QHcreate build,	G4QHadronPair* pdefs)
{
  G4double mass;
  // Example how to use an interface to different member functions
	 if(build==0) build=&G4QHadronBuilder::BuildLowSpin; // @@ Build S Hadrons?
  G4QHadron* Hadron1 = 0;                            // @@ Not initialized
  G4QHadron* Hadron2 = 0;
  if(!FourQuarkString() )
  {
    // spin 0 meson or spin 1/2 barion will be built
    Hadron1 = (hadronizer->*build)(GetLeftParton(), GetRightParton());
    mass    = Hadron1->GetMass();
  }
  else
  {
    // string is qq--qqbar: Build two stable hadrons, with extra uubar or ddbar quark pair
	   G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
	   if (GetLeftParton()->GetPDGCode() < 0) iflc = -iflc;
	   //... theSpin = 4; spin 3/2 baryons will be built
	   Hadron1 = (hadronizer->*build)(GetLeftParton(),CreateParton(iflc));
	   Hadron2 = (hadronizer->*build)(GetRightParton(),CreateParton(-iflc));
    mass    = Hadron1->GetMass() + Hadron2->GetMass();
  }
	 if(pdefs) // need to return hadrons as well.... 
	 {
	   pdefs->first  = Hadron1;
	   pdefs->second = Hadron2;
	 } 
  return mass;
} // End of FragmentationMass

 // Checks that the string is qq-(qq)bar string
G4bool G4QString::FourQuarkString() const
{
  return    GetLeftParton()->GetParticleSubType() == "di_quark"
         && GetRightParton()->GetParticleSubType()== "di_quark";
} // End of FourQuarkString

void G4QString::CalculateHadronTimePosition(G4double StringMass, G4QHadronVector* Hadrons)
{
  // `yo-yo` formation time
  static const G4double kappa = 1.0 * GeV/fermi;
  static const G4double dkappa = kappa+kappa;
  for(unsigned c1 = 0; c1 < Hadrons->size(); c1++)
  {
    G4double SumPz = 0; 
    G4double SumE  = 0;
    for(unsigned c2 = 0; c2 < c1; c2++)
    {
      G4LorentzVector hc2M=(*Hadrons)[c2]->Get4Momentum();
      SumPz += hc2M.pz();
      SumE  += hc2M.e();   
    }
    G4QHadron* hc1=(*Hadrons)[c1];
    G4LorentzVector hc1M=hc1->Get4Momentum();
    G4double HadronE = hc1M.e();
    G4double HadronPz= hc1M.pz();
    hc1->SetFormationTime((StringMass-SumPz-SumPz+HadronE-HadronPz)/dkappa);
    hc1->SetPosition(G4ThreeVector(0,0,(StringMass-SumE-SumE-HadronE+HadronPz)/dkappa));
  } 
} // End of CalculateHadronTimePosition

void G4QString::SetLeftPartonStable()
{
     theStableParton=GetLeftParton();
     theDecayParton=GetRightParton();
     decaying=Right;
}

void G4QString::SetRightPartonStable()
{
     theStableParton=GetRightParton();
     theDecayParton=GetLeftParton();
     decaying=Left;
}

G4int G4QString::GetDecayDirection() const
{
	 if      (decaying == Left ) return +1;
	 else if (decaying == Right) return -1;
	 else
	 {
				G4cerr<<"***G4QString::GetDecayDirection: wrong DecayDirection="<<decaying<<G4endl;
    G4Exception("G4QString::GetDecayDirection:","72",FatalException,"WrongDecayDirection");
	 }
	 return 0;
}

G4ThreeVector G4QString::StablePt()
{
	 if (decaying == Left ) return Ptright;
	 else if (decaying == Right ) return Ptleft;
	 else
	 {
				G4cerr<<"***G4QString::StablePt: wrong DecayDirection="<<decaying<<G4endl;
    G4Exception("G4QString::StablePt:","72",FatalException,"WrongDecayDirection");
	 }
	 return G4ThreeVector();
}

G4ThreeVector G4QString::DecayPt()
{
	 if      (decaying == Left  ) return Ptleft;
	 else if (decaying == Right ) return Ptright;
	 else
	 {
				G4cerr<<"***G4QString::DecayPt: wrong DecayDirection="<<decaying<<G4endl;
    G4Exception("G4QString::DecayPt:","72",FatalException,"WrongDecayDirection");
	 }
	 return G4ThreeVector();
}

G4double G4QString::LightConeDecay()
{
	 if      (decaying == Left  ) return Pplus;
	 else if (decaying == Right ) return Pminus;
	 else
	 {
				G4cerr<<"***G4QString::LightConeDecay: wrong DecayDirection="<<decaying<<G4endl;
    G4Exception("G4QString::LightConeDecay:","72",FatalException,"WrongDecayDirection");
	 }
	 return 0;
}

G4LorentzVector G4QString::GetFragmentation4Mom() const
{
  G4LorentzVector momentum(Ptleft+Ptright,0.5*(Pplus+Pminus));
	 momentum.setPz(0.5*(Pplus-Pminus));
	 return momentum;
}

// Random choice of string end to use it for creating the hadron (decay)   
G4QHadron* G4QString::Splitup(G4bool QL)
{
  SideOfDecay = (G4UniformRand() < 0.5)? 1: -1;
  if(SideOfDecay<0) SetLeftPartonStable();  // Decay Right parton
  else              SetRightPartonStable(); // Decay Left parton
  G4QParton* newStringEnd;
  G4QHadron* Hadron;
  if(DecayIsQuark()) Hadron=QuarkSplitup(GetDecayParton(), newStringEnd);    // MF1
  else               Hadron= DiQuarkSplitup(GetDecayParton(), newStringEnd); // MF2
  // create new String from old, ie. keep Left and Right order, but replace decay
  G4LorentzVector* HadronMomentum=SplitEandP(Hadron, QL);                    // MF3
  if(HadronMomentum)
  {    
	   G4ThreeVector   Pos(0.,0.,0.);
	   Hadron->Set4Momentum(*HadronMomentum);
	   UpdateString(newStringEnd, HadronMomentum);
	   delete HadronMomentum;
  }      
  return Hadron;
} // End of Splitup

void G4QString::UpdateString(G4QParton* decay, const G4LorentzVector* mom)
{
	 decaying=None;
	 if(decaying == Left)
	 {
    G4QParton* theFirst = thePartons[0];
				delete theFirst;
		  theFirst = decay;
		  Ptleft  -= mom->vect();
		  Ptleft.setZ(0.);
	 }
  else if (decaying == Right)
	 {
    G4QParton* theLast = thePartons[thePartons.size()-1];
    delete theLast;
		  theLast  = decay;
		  Ptright -= mom->vect();
		  Ptright.setZ(0.);
	 }
  else
	 {
				G4cerr<<"***G4QString::UpdateString: wrong oldDecay="<<decaying<<G4endl;
    G4Exception("G4QString::UpdateString","72",FatalException,"WrongDecayDirection");
	 }
	 Pplus  -= mom->e() + mom->pz();
	 Pminus -= mom->e() - mom->pz();
} // End of UpdateString

// QL=true for QGSM and QL=false for Lund fragmentation
G4LorentzVector* G4QString::SplitEandP(G4QHadron* pHadron, G4bool QL)
{
  G4double HadronMass = pHadron->GetMass();
  // calculate and assign hadron transverse momentum component HadronPx andHadronPy
  G4ThreeVector HadronPt = SampleQuarkPt() + DecayPt();
  HadronPt.setZ(0.);
  //...  sample z to define hadron longitudinal momentum and energy
  //... but first check the available phase space
  G4double DecayQuarkMass2  = sqr(GetDecayParton()->GetMass()); // Mass of quark? M.K.
  G4double HadronMass2T = HadronMass*HadronMass + HadronPt.mag2();
  if (DecayQuarkMass2 + HadronMass2T >= SmoothParam*Mass2() )  return 0;		// restart!
  //... then compute allowed z region  z_min <= z <= z_max 
  G4double zMin = HadronMass2T/Mass2();
  G4double zMax = 1. - DecayQuarkMass2/Mass2();
  if (zMin >= zMax) return 0;		// have to start all over!  
  G4double z=0;
  if(QL) z = GetQGSMLightConeZ(zMin, zMax, GetDecayParton()->GetPDGCode(), pHadron,
                               HadronPt.x(), HadronPt.y());      
  else   z = GetLundLightConeZ(zMin, zMax, GetDecayParton()->GetPDGCode(), pHadron,
                               HadronPt.x(), HadronPt.y());      
  //... now compute hadron longitudinal momentum and energy
  // longitudinal hadron momentum component HadronPz
  G4double zl= z*LightConeDecay();
  G4double HadronE = (zl + HadronMass2T/zl)/2;
  HadronPt.setZ( GetDecayDirection() * (zl - HadronE) );
  G4LorentzVector* a4Momentum= new G4LorentzVector(HadronPt,HadronE);

  return a4Momentum;
}

G4ThreeVector G4QString::SampleQuarkPt()
{
   G4double Pt = SigmaQT * std::sqrt( -std::log(G4UniformRand()));
   G4double phi = twopi*G4UniformRand();
   return G4ThreeVector(Pt * std::cos(phi),Pt * std::sin(phi),0);
}

void G4QString::Sample4Momentum(G4LorentzVector* Mom, G4double Mass,
                         G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass)
{
  G4double m2 = Mass*Mass;
  G4double am2= AntiMass*AntiMass;
  G4double dub=InitialMass*InitialMass - m2 - am2;
  G4double r_val = dub - 4*m2*am2;
  G4double Pabs = (r_val > 0.) ? std::sqrt(r_val)/(InitialMass*InitialMass) : 0;
  //... sample unit vector       
  G4double r  = G4UniformRand();                    // @@ G4RandomDirection()
  G4double pz = 1. - r - r;                         // cos(theta)
  G4double st = std::sqrt(1. - pz * pz) * Pabs;
  G4double phi= twopi*G4UniformRand();
  G4double px = st*std::cos(phi);
  G4double py = st*std::sin(phi);
  pz *= Pabs;
  G4double p2=Pabs*Pabs;
  Mom->setPx(px); Mom->setPy(py); Mom->setPz(pz);
  Mom->setE(std::sqrt(p2 + Mass*Mass));
  AntiMom->setPx(-px); AntiMom->setPy(-py); AntiMom->setPz(-pz);
  AntiMom->setE (std::sqrt(Pabs*Pabs + AntiMass*AntiMass));
} // End of Sample4Momentum

G4bool G4QString::SplitLast(G4QString* string, G4QHadronVector* LeftVector,
                                               G4QHadronVector* RightVector)
{
  //... perform last cluster decay
  G4ThreeVector ClusterVel =string->Get4Momentum().boostVector();
  G4double ResidualMass    =string->Mass(); 
  G4double ClusterMassCut = ClusterMass;
  G4int cClusterInterrupt = 0;
  G4QHadron* LeftHadron;
  G4QHadron* RightHadron;
  do
  {
    if(cClusterInterrupt++ >= ClusterLoopInterrupt) return false;
	   G4QParton* quark = 0;
	   string->SetLeftPartonStable(); // to query quark contents..
	   if(string->DecayIsQuark() && string->StableIsQuark()) // There're quarks on clusterEnds
		    LeftHadron= QuarkSplitup(string->GetLeftParton(), quark);
	   else
    {
	     //... there is a Diquark on cluster ends
		    G4int IsParticle;
		    if(string->StableIsQuark())IsParticle=(string->GetLeftParton()->GetPDGCode()>0)?-1:1;
		    else                       IsParticle=(string->GetLeftParton()->GetPDGCode()>0)?1:-1;
      G4QPartonPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      quark = QuarkPair.GetParton2();
      LeftHadron=hadronizer->Build(QuarkPair.GetParton1(), string->GetLeftParton());
	   }
    RightHadron = hadronizer->Build(string->GetRightParton(), quark);
    //... repeat procedure, if mass of cluster is too low to produce hadrons
    //... ClusterMassCut = 0.15*GeV model parameter
	   if ( quark->GetParticleSubType()== "quark" ) ClusterMassCut = 0.;
	   else                                         ClusterMassCut = ClusterMass;
  } while(ResidualMass <= LeftHadron->GetMass() + RightHadron->GetMass() + ClusterMassCut);
  //... compute hadron momenta and energies   
  G4LorentzVector  LeftMom, RightMom;
  G4ThreeVector    Pos;
  Sample4Momentum(&LeftMom,LeftHadron->GetMass(),&RightMom,RightHadron->GetMass(),
                                                                             ResidualMass);
  LeftMom.boost(ClusterVel);
  RightMom.boost(ClusterVel);
  LeftVector->push_back(new G4QHadron(LeftHadron, 0, Pos, LeftMom));
  RightVector->push_back(new G4QHadron(RightHadron, 0, Pos, RightMom));

  return true;
} // End of SplitLast

G4QHadron*	G4QString::QuarkSplitup(G4QParton*	decay, G4QParton *&created)
{
  G4int IsParticle=(decay->GetPDGCode()>0) ? -1 : +1; // if we have a quark, we need antiquark (or diquark)
  G4QPartonPair QuarkPair = CreatePartonPair(IsParticle);
  created = QuarkPair.GetParton2();
  return hadronizer->Build(QuarkPair.GetParton1(), decay);
} // End of QuarkSplitup

//
G4QHadron* G4QString::DiQuarkSplitup(G4QParton* decay, G4QParton* &created)
{
  //... can Diquark break or not? 
  if (G4UniformRand() < DiquarkBreakProb )
  {
    //... Diquark break
    G4int stableQuarkEncoding = decay->GetPDGCode()/1000;
    G4int decayQuarkEncoding = (decay->GetPDGCode()/100)%10;
    if (G4UniformRand() < 0.5)
    {
      G4int Swap = stableQuarkEncoding;
      stableQuarkEncoding = decayQuarkEncoding;
      decayQuarkEncoding = Swap;
    }
    G4int IsParticle=(decayQuarkEncoding>0) ? -1 : +1; 
	   // if we have a quark, we need antiquark)
    G4QPartonPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
    //... Build new Diquark
    G4int QuarkEncoding=QuarkPair.GetParton2()->GetPDGCode();
    G4int i10  = std::max(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
    G4int i20  = std::min(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
    G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
    G4int NewDecayEncoding = -1*IsParticle*(i10 * 1000 + i20 * 100 + spin);
    created = CreateParton(NewDecayEncoding);
    G4QParton* decayQuark=CreateParton(decayQuarkEncoding);
    return hadronizer->Build(QuarkPair.GetParton1(), decayQuark);
  }
  else
  {
    //... Diquark does not break
    G4int IsParticle=(decay->GetPDGCode()>0) ? +1 : -1; 
    	// if we have a diquark, we need quark)
    G4QPartonPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
    created = QuarkPair.GetParton2();
    return hadronizer->Build(QuarkPair.GetParton1(), decay);
 }
} // End of DiQuarkSplitup

G4QPartonPair G4QString::CreatePartonPair(G4int NeedParticle, G4bool AllowDiquarks)
{
  //  NeedParticle = {+1 for Particle, -1 for AntiParticle}
  if(AllowDiquarks && G4UniformRand() < DiquarkSuppress)
  {
    // Create a Diquark - AntiDiquark pair , first in pair is anti to IsParticle
    G4int q1  = SampleQuarkFlavor();
    G4int q2  = SampleQuarkFlavor();
    G4int spin = (q1 != q2 && G4UniformRand() <= 0.5)? 1 : 3;
    // Convention: quark with higher PDG number is first
    G4int PDGcode = (std::max(q1,q2) * 1000 + std::min(q1,q2) * 100 + spin) * NeedParticle;
    return G4QPartonPair(CreateParton(-PDGcode),CreateParton(PDGcode));
  }
  else
  {
    // Create a Quark - AntiQuark pair, first in pair  IsParticle
    G4int PDGcode=SampleQuarkFlavor()*NeedParticle;
    return G4QPartonPair(CreateParton(PDGcode),CreateParton(-PDGcode));
  }
} // End of CreatePartonPair

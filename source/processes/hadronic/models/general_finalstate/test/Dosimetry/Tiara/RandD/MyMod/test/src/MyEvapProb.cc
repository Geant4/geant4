//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: MyEvapProb.cc,v 1.1 2003-10-08 12:32:19 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "MyEvapProb.hh"
#include "G4PairingCorrection.hh"

MyEvapProb::MyEvapProb(const G4int anA, const G4int aZ, const G4double aGamma) : 
  theA(anA),
  theZ(aZ),
  Gamma(aGamma) 
{
  theEvapLDPptr = new MyLDP;
}


MyEvapProb::MyEvapProb(const MyEvapProb &right)
{
    G4Exception("G4EvaporationProbability::copy_constructor meant to not be accessable");
}




const MyEvapProb & MyEvapProb::
operator=(const MyEvapProb &right)
{
    G4Exception("G4EvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool MyEvapProb::operator==(const MyEvapProb &right) const
{
    return false;
}

G4bool MyEvapProb::operator!=(const MyEvapProb &right) const
{
    return true;
}

G4double MyEvapProb::EmissionProbability(const G4Fragment & fragment, const G4double anEnergy)
{
    G4double EmissionProbability = 0.0;
    G4double MaximalKineticEnergy = anEnergy;
    if (MaximalKineticEnergy > 0.0 && fragment.GetExcitationEnergy() > 0.0) {
	EmissionProbability = CalcProbability(fragment,MaximalKineticEnergy);
	// 		// Next there is a loop over excited states for this channel summing probabilities
	// 		G4double SavedGamma = Gamma;
	// 		for (G4int i = 0; i < ExcitationEnergies->length(); i++) {
	// 			if (ExcitationSpins->operator()(i) < 0.1) continue;
	// 			Gamma = ExcitationSpins->operator()(i)*A; // A is the channel mass number
	// 			// substract excitation energies
	// 			MaximalKineticEnergy -= ExcitationEnergies->operator()(i);
	// 			// update probability
	// 			EmissionProbability += CalcProbability(fragment,MaximalKineticEnergy);
	// 				EmissionProbability += tmp;
	// 			}
	// 		// restore Gamma and MaximalKineticEnergy
	// 		MaximalKineticEnergy = SavedMaximalKineticEnergy;
	// 		Gamma = SavedGamma;
	// 		}
    }
    return ((EmissionProbability<0) ? 0 : EmissionProbability);
}


double MyEvapProb::ro(double U,int Prime,int mat)
{
  if(U < 1e-10) return 0;
  double sq;
    if(mat==0) sq = theEvapLDPptr->LevelDensityParameter(m_NuclA,m_NuclB,U);
  else sq = theEvapLDPptr->LevelDensityParameter(fragm.GetA(),fragm.GetZ(),U);
  double a = sq;
  double Ux;
  if(mat==0) Ux = 150/G4NucleiProperties::GetNuclearMass(m_NuclA,m_NuclB)+2.4;
  else Ux = 150/G4NucleiProperties::GetNuclearMass(fragm.GetA(),fragm.GetZ())+2.4;
  double T = 1/(sqrt(a/U)/* - 1.5/Ux*/),res=0,E0=0;
  if(U >= Ux+m_delta0){
    U = U/*-m_delta0*/;
    sq = sqrt(sq*(U));
    switch(Prime){
    case 0:
      res= exp(2*sq)/pow(a,1./4.)/pow(U,5./4.); break;
    case 1:
      res=ro(U)*(1/sqrt(U)-5./(4.*(U))); break;
    case 2:
      res= ro(U)*(1/sqrt(U)-5./(4.*(U)) + (5./(4*(U))-1/(2*sqrt(U))))/(U);
    default : break;
    }
  }
  else{
    E0 = Ux+m_delta0 - T*(log(T)-0.25*log(a)-1.25*log(Ux)+2*sqrt(a*Ux));
    if(Prime==0) res= exp((U-E0+m_delta0)/T)/T;
    else if(Prime==1) res= exp((U-E0+m_delta0)/T)/T/T;
    else if(Prime==2) res= exp((U-E0+m_delta0)/T)/T/T/T;
  }
  //  G4cout<<Ux<<" "<<Ux+m_delta0<<" "<<U<<" "<<T<<" "<<E0<<" "<<res<<G4endl;
  return res;
}

double MyEvapProb::func(double U,int Prime)
{
  G4double res=0;
  switch(Prime){
  case 0: res= ro(Estar-U)/ro(Estar,0,1)*(U+Beta);
      if(res<0){
	G4cout<<"Negative function in integral "<<Estar-U<<" "<<ro(Estar-U)<<" "<<Emax<<" "<<ro(Estar,0,1)<<" "<<(U+Beta)<<G4endl;
	return res/(U+Beta);
      }
      break;
  case 1:
    if(U+Beta>0)
      res = Beta*ro(Estar-U)/ro(Estar,0,1) - ro(Estar-U,1)/ro(Estar,0,1)*(U+Beta);
    else res = ro(Estar-U,1)/ro(Estar,0,1);
    break;
  case 2:
    if(U+Beta > 0)
      res = -2*ro(Estar-U,1)/ro(Estar,0,1)*Beta + ro(Estar-U,2)/ro(Estar,0,1)*(U+Beta);
    else res = ro(Estar-U,2)/ro(Estar,0,1);
  }
  return res;
}

double MyEvapProb::Integrator(double Min,double Max,double Left,double Right)
{
  if(Max-Min < 1e-5) return (Max-Min)*(Left+Right)*0.5;
  G4double middle1 = (Max+Min)*0.5;
  G4double middle = func(middle1);
  G4double S = (middle1-Min)*(middle+Left) + (Max-middle1)*(Right+middle);
  if(fabs((S-(Max-Min)*(Left+Right))/S) <0.005*(Max-Min)){
    return S;
  }
  return Integrator(Min,middle1,Left,middle) + Integrator(middle1,Max,middle,Right);
}

G4double MyEvapProb::BindingEnergy()
{
  // G4NucleTable* pTbl = G4ParticleTable::GetParticleTable();
  /*  double Mass = G4NucleiProperties::GetNuclearMass(fragm.GetA(),fragm.GetZ());
      Mass -= G4NucleiProperties::GetNuclearMass(fragm.GetA()-theA,fragm.GetZ()-theZ);*/
  return G4NucleiProperties::GetBindingEnergy(fragm.GetA()-theA,fragm.GetZ()-theZ)*theA;
}

G4double MyEvapProb::FindMaximum(G4double left,G4double right,G4double hint)
{
  G4double MaxValue=0;
  G4double f;
  G4double T,S=0,S1,i = (right-left);
  T = theEvapLDPptr->LevelDensityParameter(fragm.GetA()-theA,fragm.GetZ()-theZ,(left+right)*0.5);
  if(T>5){
    f = /*(sqrt(T+right+1./4.)-1./2.)/T*/(S=(T-sqrt(T*(T-5))/2))*S/T+m_delta0;
    if(f < (S=(T+sqrt(T*(T-5))/2))*S/T+m_delta0)
      f = S*S/T+m_delta0;
    T = f;
    if(T > right) T = right;
    if(T<left) T = left;
  }
  else T = right;
  if(hint!=0) T=hint;
  MaxValue = func(T);
  do{
    S = T;
    S1 = T + func(T,1)/fabs(func(T,2));
    if(S1<left)
      S1 = (right-left)*G4UniformRand()*i+left;
    else if(S1>right)
      S1 = right - (right-left)*G4UniformRand()*i;
    f = func(S1);
    if(f > MaxValue){
      T = S1;
      MaxValue = f;
    }
    else{
      if(exp((MaxValue-f)/i)*i < G4UniformRand()){
	T = S1;
	MaxValue = f;
      }
    }
    i *= 0.7;
  }while(fabs((T-S)/S) > 1e-10);
  return MaxValue;
}

G4double MyEvapProb::SampleEnergy(const G4Fragment& fragment,const G4double MaximalKineticEnergy)
{
  G4double Max,left,right,random1,f;
  fragm = fragment;
  //  fragm.SetA(m_NuclA);
  //  fragm.SetZ(m_NuclB);
  left = CalcBetaParam(fragm);
  if(left<0) left = 0;
  right = MaximalKineticEnergy;
  /*  right = fragment.GetExcitationEnergy();
  right = right - BindingEnergy()-
#ifdef __USE_BIASING__
G4PairingCorrection::GetInstance()->GetPairingCorrection(fragment.GetA(),fragment.GetZ());
#else
G4PairingCorrection::GetPairingCorrection(fragment.GetA(),fragment.GetZ());
#endif
  */
  Max = FindMaximum(left,right)*10;
  do{
    random1 = G4UniformRand();
    f = func((right-left)*random1+left);
    if(f>Max){
      G4cout<<"Error in finding maximum: on "<<random1<<" it is "<<(f-Max)<<"/"<<Max<<" percent higher"<<G4endl;
      Max = FindMaximum(left,right,random1);
      continue;
    }
  }while(f < G4UniformRand()*Max);
  return random1;
}


G4double MyEvapProb::CalcProbability(const G4Fragment & fragment, 
						   const G4double MaximalKineticEnergy)
    // Calculate integrated probability (width) for rvaporation channel
{
  fragm = fragment;
  G4double ResidualA = G4double(fragment.GetA() - theA);
  G4double ResidualZ = G4double(fragment.GetZ() - theZ);
  m_NuclA = ResidualA;
  m_NuclB = ResidualZ;
  //  fragm.SetA(ResidualA);
  // fragm.SetZ(ResidualZ);
  G4double U = fragment.GetExcitationEnergy();
  G4double A = fragment.GetA();
  G4double Z = fragment.GetZ();
	
  G4double NuclearMass = G4NucleiProperties::GetNuclearMass(theA,theZ);
  Estar = U -BindingEnergy();
  if(Estar < 0) return 0;

  G4double Alpha = CalcAlphaParam(fragment);
  G4double Beta = CalcBetaParam(fragment);

  G4double delta0 =
#ifdef __USE_BIASING__
G4PairingCorrection::GetInstance()->GetPairingCorrection(fragment.GetA(),fragment.GetZ());
#else
G4PairingCorrection::GetPairingCorrection(fragment.GetA(),fragment.GetZ());
#endif
  m_delta0 = delta0;
  G4double RN = 1.5*fermi;
  /*
  G4double Ux = 150/ResidualA + 2.5;
  G4double Ex = Ux + delta0;
  G4double LDP = theEvapLDPptr->LevelDensityParameter(ResidualA,ResidualZ,U);
  G4double T = 1/(sqrt(LDP/Ux)- 1.5/Ux);
  G4double E0 = Ex -T*(log(T)-0.25*log(LDP) -1.25*log(Ux) + 2*sqrt(LDP*Ux));
  G4double t = (Estar-CoulombBarrier(fragment))/T;
  G4double tx = Ex/T;
  G4double s = 2*sqrt(LDP*(Estar-CoulombBarrier(fragment)-delta0));
  G4double sx = 2*sqrt(LDP*(Ex-delta0));
  G4double I0 = exp(-E0/T)*(exp(t)-1);
  G4double I1 = exp(-E0/T)*T*((t-tx+1)*exp(tx)-t-1);
  G4double res;
  G4cout<<T<<" "<<E0<<" "<<t<<" "<<tx<<" "<<s<<" "<<sx<<" "<<I0<<" "<<I1<<G4endl;
  if(Estar-CoulombBarrier(fragment) < Ex){
    if(CoulombBarrier(fragment)!=0)
      res = I1+(Beta/CoulombBarrier(fragment)+CoulombBarrier(fragment))*I0;
    else
      res = I1+Beta*I0;
  }
  else{
    G4double I2 = 2*1.4142136*(1/(s*sqrt(s))+1.5/(s*s*sqrt(s))+3.75/(s*s*s*sqrt(s)) - (1/(s*sqrt(s)) + 1.5/(s*s*sqrt(s)) + 3.75/(s*s*s*sqrt(s)))*exp(sx-s));
    G4double I3 = 2./sqrt(s) + 4./(s*sqrt(s)) + 13.5/(s*s*sqrt(s)) + 60/(s*s*s*sqrt(s)) + 325.125/(s*s*s*s*sqrt(s));
    I3 -= ((s*s-sx*sx)/(sx*sqrt(sx))+(1.5*s*s + 0.5*sx*sx)/(sx*sx*sqrt(s))+ (3.75*s*s+0.25*sx*sx)/(sx*sx*sx*sqrt(sx)) +(12.875*s*s + 0.625*sx*sx)/(sx*sx*sx*sx*sqrt(sx)) + (59.0625*s*s + 0.9375*sx*sx)/(sx*sx*sx*sx*sx*sqrt(sx)) + (324.8*s*s+3.28*sx*sx)/(sx*




sx*sx*sx*sx*sx*sqrt(sx)))*exp(sx-s);
    I3 /= (1.4142136*LDP);
    G4cout<<I2<<" "<<I3<<G4endl;
    if(CoulombBarrier(fragment)!=0)
      res = I1 + I3*exp(s) + (Beta/CoulombBarrier(fragment)+CoulombBarrier(fragment))*(I0 + I2*exp(s));
    else
      res = I1 + I3*exp(s) + Beta*(I0 + I2*exp(s));
  }
  if(U<Ex){
    res /= 1/T*exp((U-E0)/T);
  }
  else{
    res /= pi/12*exp(2*sqrt(U-delta0))/(sqrt(sqrt(LDP*(U-delta0)))*(U-delta0));
  }
  G4cout<<res<<G4endl;
  G4cout<<Alpha/12*(RN*RN*pi*pow(theA,2./3.))*Gamma*NuclearMass*res/(pi*pi*hbar_Planck*hbar_Planck)<<G4endl;
  return Alpha/12*(RN*RN*pi*pow(theA,2./3.))*Gamma*NuclearMass*res/(pi*pi*hbar_Planck*hbar_Planck);
  */
  Emax = MaximalKineticEnergy;
  G4double Coef = Alpha*RN*RN*pow(ResidualA,2./3.)*(Gamma)*NuclearMass/(8*pi*pi*pi/*hbar_Planck*hbar_Planck*/);
  G4double Width;
  if(GetZ()==0)
    Width = Coef*Integrator(0,Estar-delta0,func(0),func(Estar-delta0));
  else{
    Width = Coef*Integrator(fabs(Beta),Estar-delta0,func(fabs(Beta)),func(Estar-delta0));
  }
  return ((Width<0.) ? 0. : Width);
}



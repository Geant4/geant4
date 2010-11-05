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
// $Id: G4GEMProbability.cc,v 1.15 2010-11-05 14:43:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------
//
// Geant4 class G4GEMProbability
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept 2001) 
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept 2001)
//
// J. M. Quesada : several fixes in total GEM width 
// J. M. Quesada 14/07/2009 bug fixed in total emission width (hbarc)
// J. M. Quesada (September 2009) several fixes:
//       -level density parameter of residual calculated at its right excitation energy.
//       -InitialLeveldensity calculated according to the right conditions of the 
//        initial excited nucleus.
// J. M. Quesada 19/04/2010 fix in  emission probability calculation.
// V.Ivanchenko  20/04/2010 added usage of G4Pow and use more safe computation
// V.Ivanchenko 18/05/2010 trying to speedup the most slow method
//            by usage of G4Pow, integer Z and A; moved constructor, 
//            destructor and virtual functions to source

#include "G4GEMProbability.hh"
#include "G4PairingCorrection.hh"
#include "G4Pow.hh"
#include "G4IonTable.hh"

G4GEMProbability:: G4GEMProbability() : 
  theA(0), theZ(0), Spin(0.0), theCoulombBarrierPtr(0),
  ExcitationEnergies(0), ExcitationSpins(0), ExcitationLifetimes(0), Normalization(1.0)
{
  theEvapLDPptr = new G4EvaporationLevelDensityParameter;
  fG4pow = G4Pow::GetInstance(); 
  fPairCorr = G4PairingCorrection::GetInstance();
}

G4GEMProbability:: G4GEMProbability(G4int anA, G4int aZ, G4double aSpin) : 
  theA(anA), theZ(aZ), Spin(aSpin), theCoulombBarrierPtr(0),
  ExcitationEnergies(0), ExcitationSpins(0), ExcitationLifetimes(0), Normalization(1.0)
{
  theEvapLDPptr = new G4EvaporationLevelDensityParameter;
  fG4pow = G4Pow::GetInstance(); 
  fPairCorr = G4PairingCorrection::GetInstance();
}
    
G4GEMProbability::~G4GEMProbability()
{
  delete theEvapLDPptr;
}

G4double G4GEMProbability::CalcAlphaParam(const G4Fragment & ) const 
{
  return 1.0;
}

G4double G4GEMProbability::CalcBetaParam(const G4Fragment & ) const 
{
  return 1.0;
}

G4double G4GEMProbability::CCoeficient(const G4double ) const 
{
  return 0.0;
}

G4double G4GEMProbability::EmissionProbability(const G4Fragment & fragment,
                                               G4double MaximalKineticEnergy)
{
  G4double EmissionProbability = 0.0;
    
  if (MaximalKineticEnergy > 0.0 && fragment.GetExcitationEnergy() > 0.0) {
    G4double CoulombBarrier = GetCoulombBarrier(fragment);
      
    EmissionProbability = CalcProbability(fragment,MaximalKineticEnergy,CoulombBarrier);
    Normalization = EmissionProbability;
    // Next there is a loop over excited states for this channel summing probabilities
    if (ExcitationEnergies  &&  ExcitationSpins && ExcitationLifetimes) {
      G4double SavedSpin = Spin;
      for (size_t i = 0; i < ExcitationEnergies->size(); ++i) {
	Spin = ExcitationSpins->operator[](i);
	// substract excitation energies
	G4double Tmax = MaximalKineticEnergy - ExcitationEnergies->operator[](i);
	if (Tmax > 0.0) {
	  G4double width = CalcProbability(fragment,Tmax,CoulombBarrier);
	  //JMQ April 2010 added condition to prevent reported crash
	  // update probability
	  if (width > 0. && 
	      hbar_Planck*fG4pow->logZ(2) < width*ExcitationLifetimes->operator[](i)) {
	    EmissionProbability += width;
	  }
	}
      }
      // Restore Spin
      Spin = SavedSpin;
    }
  }
  return EmissionProbability;
}


G4double G4GEMProbability::CalcProbability(const G4Fragment & fragment, 
                                           G4double MaximalKineticEnergy,
                                           G4double V)

// Calculate integrated probability (width) for evaporation channel
{
  G4int A = fragment.GetA_asInt();
  G4int Z = fragment.GetZ_asInt();

  G4int ResidualA = A - theA;
  G4int ResidualZ = Z - theZ;
  G4double U = fragment.GetExcitationEnergy();
  
  G4double NuclearMass = fragment.ComputeGroundStateMass(theZ, theA);
  
  G4double Alpha = CalcAlphaParam(fragment);
  G4double Beta = CalcBetaParam(fragment);
    
  //                       ***RESIDUAL***
  //JMQ (September 2009) the following quantities refer to the RESIDUAL:
  
  G4double delta0 = fPairCorr->GetPairingCorrection(ResidualA, ResidualZ);  
  
  G4double a = theEvapLDPptr->LevelDensityParameter(ResidualA,
  						    ResidualZ,MaximalKineticEnergy+V-delta0);
  G4double Ux = (2.5 + 150.0/G4double(ResidualA))*MeV;
  G4double Ex = Ux + delta0;
  G4double T  = 1.0/(std::sqrt(a/Ux) - 1.5/Ux);
  //JMQ fixed bug in units
  G4double E0 = Ex - T*(std::log(T/MeV) - std::log(a*MeV)/4.0 
	- 1.25*std::log(Ux/MeV) + 2.0*std::sqrt(a*Ux));
  //                      ***end RESIDUAL ***
  
  //                       ***PARENT***
  //JMQ (September 2009) the following quantities refer to the PARENT:
     
  G4double deltaCN = fPairCorr->GetPairingCorrection(A, Z);				       
  G4double aCN     = theEvapLDPptr->LevelDensityParameter(A, Z, U-deltaCN);
  G4double UxCN    = (2.5 + 150.0/G4double(A))*MeV;
  G4double ExCN    = UxCN + deltaCN;
  G4double TCN     = 1.0/(std::sqrt(aCN/UxCN) - 1.5/UxCN);
  //	                 ***end PARENT***

  G4double Width;
  G4double InitialLevelDensity;
  G4double t = MaximalKineticEnergy/T;
  if ( MaximalKineticEnergy < Ex ) {
    //JMQ 190709 bug in I1 fixed (T was  missing)
    Width = (I1(t,t)*T + (Beta+V)*I0(t))/std::exp(E0/T);
    //JMQ 160909 fix:  InitialLevelDensity has been taken away 
    //(different conditions for initial CN..) 
  } else {

    //VI minor speedup
    G4double expE0T = std::exp(E0/T);
    const G4double sqrt2 = std::sqrt(2.0);

    G4double tx = Ex/T;
    G4double s = 2.0*std::sqrt(a*(MaximalKineticEnergy-delta0));
    G4double sx = 2.0*std::sqrt(a*(Ex-delta0));
    Width = I1(t,tx)*T/expE0T + I3(s,sx)*std::exp(s)/(sqrt2*a);
    // For charged particles (Beta+V) = 0 beacuse Beta = -V
    if (theZ == 0) {
      Width += (Beta+V)*(I0(tx)/expE0T + 2.0*sqrt2*I2(s,sx)*std::exp(s));
    }
  }
  
  //JMQ 14/07/2009 BIG BUG : NuclearMass is in MeV => hbarc instead of hbar_planck must be used
  //    G4double g = (2.0*Spin+1.0)*NuclearMass/(pi2* hbar_Planck*hbar_Planck);
  G4double g = (2.0*Spin+1.0)*NuclearMass/(pi2* hbarc*hbarc);
  
  //JMQ 190709 fix on Rb and  geometrical cross sections according to Furihata's paper 
  //                      (JAERI-Data/Code 2001-105, p6)
  //    G4double RN = 0.0;
  G4double Rb = 0.0;
  if (theA > 4) 
    {
      //	G4double R1 = std::pow(ResidualA,1.0/3.0);
      //	G4double R2 = std::pow(G4double(theA),1.0/3.0);
      G4double Ad = fG4pow->Z13(ResidualA);
      G4double Aj = fG4pow->Z13(theA);
      //	RN = 1.12*(R1 + R2) - 0.86*((R1+R2)/(R1*R2));
      Rb = 1.12*(Aj + Ad) - 0.86*((Aj+Ad)/(Aj*Ad))+2.85;
      Rb *= fermi;
    }
  else if (theA>1)
    {
      G4double Ad = fG4pow->Z13(ResidualA);
      G4double Aj = fG4pow->Z13(theA);
      Rb=1.5*(Aj+Ad)*fermi;
    }
  else
    {
      G4double Ad = fG4pow->Z13(ResidualA);
      Rb = 1.5*Ad*fermi;
    }
  //   G4double GeometricalXS = pi*RN*RN*std::pow(ResidualA,2./3.); 
  G4double GeometricalXS = pi*Rb*Rb; 
  //end of JMQ fix on Rb by 190709
  
  //JMQ 160909 fix: initial level density must be calculated according to the 
  // conditions at the initial compound nucleus 
  // (it has been removed from previous "if" for the residual) 

  if ( U < ExCN ) 
    {
      //JMQ fixed bug in units
      //VI moved the computation here
      G4double E0CN = ExCN - TCN*(std::log(TCN/MeV) - std::log(aCN*MeV)/4.0 
				  - 1.25*std::log(UxCN/MeV) + 2.0*std::sqrt(aCN*UxCN));
      InitialLevelDensity = (pi/12.0)*std::exp((U-E0CN)/TCN)/TCN;
    } 
  else 
    {
      //InitialLevelDensity = (pi/12.0)*std::exp(2*std::sqrt(aCN*(U-deltaCN)))/std::pow(aCN*std::pow(U-deltaCN,5.0),1.0/4.0);
      //VI speedup
      G4double x  = U-deltaCN;
      G4double x1 = std::sqrt(aCN*x);
      InitialLevelDensity = (pi/12.0)*std::exp(2*x1)/(x*std::sqrt(x1));
    }

  //JMQ 190709 BUG : pi instead of sqrt(pi) must be here according to Furihata's report:
  //    Width *= std::sqrt(pi)*g*GeometricalXS*Alpha/(12.0*InitialLevelDensity); 
  Width *= pi*g*GeometricalXS*Alpha/(12.0*InitialLevelDensity); 
   
  return Width;
}

G4double G4GEMProbability::I3(G4double s, G4double sx)
{
  G4double s2 = s*s;
  G4double sx2 = sx*sx;
  G4double S = 1.0/std::sqrt(s);
  G4double S2 = S*S;
  G4double Sx = 1.0/std::sqrt(sx);
  G4double Sx2 = Sx*Sx;
  
  G4double p1 = S *(2.0 + S2 *( 4.0 + S2 *( 13.5 + S2 *( 60.0 + S2 * 325.125 ))));
  G4double p2 = Sx*Sx2 *(
			 (s2-sx2) + Sx2 *(
					  (1.5*s2+0.5*sx2) + Sx2 *(
								   (3.75*s2+0.25*sx2) + Sx2 *(
											      (12.875*s2+0.625*sx2) + Sx2 *(
															    (59.0625*s2+0.9375*sx2) + Sx2 *(324.8*s2+3.28*sx2))))));
  
  p2 *= std::exp(sx-s);
  
  return p1-p2; 
}

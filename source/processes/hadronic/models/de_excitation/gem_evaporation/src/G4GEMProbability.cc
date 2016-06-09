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

#include "G4GEMProbability.hh"
#include "G4PairingCorrection.hh"



G4GEMProbability::G4GEMProbability(const G4GEMProbability &) : G4VEmissionProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4GEMProbability::copy_constructor meant to not be accessable");
}




const G4GEMProbability & G4GEMProbability::
operator=(const G4GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4GEMProbability::operator==(const G4GEMProbability &) const
{
  return false;
}

G4bool G4GEMProbability::operator!=(const G4GEMProbability &) const
{
  return true;
}

G4double G4GEMProbability::EmissionProbability(const G4Fragment & fragment,
                                               const G4double MaximalKineticEnergy)
{
  G4double EmissionProbability = 0.0;
    
  if (MaximalKineticEnergy > 0.0 && fragment.GetExcitationEnergy() > 0.0) {
    G4double CoulombBarrier = GetCoulombBarrier(fragment);
      
    EmissionProbability = CalcProbability(fragment,MaximalKineticEnergy,CoulombBarrier);
    Normalization = EmissionProbability;
    // Next there is a loop over excited states for this channel summing probabilities
    if (ExcitationEnergies  &&  ExcitationSpins && ExcitationLifetimes) {
      G4double SavedSpin = Spin;
      for (unsigned int i = 0; i < ExcitationEnergies->size(); i++) {
	Spin = ExcitationSpins->operator[](i);
	// substract excitation energies
	G4double Tmax = MaximalKineticEnergy - ExcitationEnergies->operator[](i);
	if (Tmax > 0.0) {
	  G4double width = CalcProbability(fragment,Tmax,CoulombBarrier);
	  // update probability
	  if (hbar_Planck*std::log(2.0)/width < ExcitationLifetimes->operator[](i)) {
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
                                           const G4double MaximalKineticEnergy,
                                           const G4double V)
// Calculate integrated probability (width) for evaporation channel
{
  G4double ResidualA = static_cast<G4double>(fragment.GetA() - theA);
  G4double ResidualZ = static_cast<G4double>(fragment.GetZ() - theZ);
  G4double U = fragment.GetExcitationEnergy();
  
  G4double NuclearMass = G4ParticleTable::GetParticleTable()->
    GetIonTable()->GetNucleusMass(theZ,theA);
  
  
  G4double Alpha = CalcAlphaParam(fragment);
  G4double Beta = CalcBetaParam(fragment);
  
  
  //                       ***RESIDUAL***
  //JMQ (September 2009) the following quantities refer to the RESIDUAL:
  
  G4double delta0 = G4PairingCorrection::GetInstance()->GetPairingCorrection( static_cast<G4int>( ResidualA ) , static_cast<G4int>( ResidualZ ) );  
  
  G4double a = theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(ResidualA),
  						    static_cast<G4int>(ResidualZ),MaximalKineticEnergy+V-delta0);
  G4double Ux = (2.5 + 150.0/ResidualA)*MeV;
  G4double Ex = Ux + delta0;
  G4double T = 1.0/(std::sqrt(a/Ux) - 1.5/Ux);
  //JMQ fixed bug in units
  G4double E0 = Ex - T*(std::log(T/MeV) - std::log(a*MeV)/4.0 - 1.25*std::log(Ux/MeV) + 2.0*std::sqrt(a*Ux));
  //                      ***end RESIDUAL ***
  
  //                       ***PARENT***
  //JMQ (September 2009) the following quantities refer to the PARENT:
     
  G4double deltaCN=G4PairingCorrection::GetInstance()->
    GetPairingCorrection(static_cast<G4int>(fragment.GetA()),static_cast<G4int>(fragment.GetZ()));				       
  G4double aCN = theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(fragment.GetA()),
						      static_cast<G4int>(fragment.GetZ()),U-deltaCN);
  G4double UxCN = (2.5 + 150.0/fragment.GetA())*MeV;
  G4double ExCN = UxCN + deltaCN;
  G4double TCN = 1.0/(std::sqrt(aCN/UxCN) - 1.5/UxCN);
  //JMQ fixed bug in units
  G4double E0CN = ExCN - TCN*(std::log(TCN/MeV) - std::log(aCN*MeV)/4.0 - 1.25*std::log(UxCN/MeV) + 2.0*std::sqrt(aCN*UxCN));
//                       ***end PARENT***

  G4double Width;
  G4double InitialLevelDensity;
  G4double t = MaximalKineticEnergy/T;
  if ( MaximalKineticEnergy < Ex ) {
    //JMQ 190709 bug in I1 fixed (T was  missing)
    Width = (I1(t,t)*T + (Beta+V)*I0(t))/std::exp(E0/T);
//JMQ 160909 fix:  InitialLevelDensity has been taken away (different conditions for initial CN..) 
  } else {
    G4double tx = Ex/T;
    G4double s = 2.0*std::sqrt(a*(MaximalKineticEnergy-delta0));
    G4double sx = 2.0*std::sqrt(a*(Ex-delta0));
    Width = I1(t,tx)*T/std::exp(E0/T) + I3(s,sx)*std::exp(s)/(std::sqrt(2.0)*a);
    // For charged particles (Beta+V) = 0 beacuse Beta = -V
    if (theZ == 0) {
      Width += (Beta+V)*(I0(tx)/std::exp(E0/T) + 2.0*std::sqrt(2.0)*I2(s,sx)*std::exp(s));
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
      G4double Ad = std::pow(ResidualA,1.0/3.0);
      G4double Aj = std::pow(G4double(theA),1.0/3.0);
      //	RN = 1.12*(R1 + R2) - 0.86*((R1+R2)/(R1*R2));
      Rb = 1.12*(Aj + Ad) - 0.86*((Aj+Ad)/(Aj*Ad))+2.85;
      Rb *= fermi;
    }
  else if (theA>1)
    {
      G4double Ad = std::pow(ResidualA,1.0/3.0);
      G4double Aj = std::pow(G4double(theA),1.0/3.0);
      Rb=1.5*(Aj+Ad)*fermi;
    }
  else
    {
      G4double Ad = std::pow(ResidualA,1.0/3.0);
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
        InitialLevelDensity = (pi/12.0)*std::exp((U-E0CN)/TCN)/TCN;
      } 
    else 
      {
        InitialLevelDensity = (pi/12.0)*std::exp(2*std::sqrt(aCN*(U-deltaCN)))/std::pow(aCN*std::pow(U-deltaCN,5.0),1.0/4.0);
      }
 //


  //JMQ 190709 BUG : pi instead of sqrt(pi) must be here according to Furihata's report:
  //    Width *= std::sqrt(pi)*g*GeometricalXS*Alpha/(12.0*InitialLevelDensity); 
  Width *= pi*g*GeometricalXS*Alpha/(12.0*InitialLevelDensity); 
   

  return Width;
}




G4double G4GEMProbability::I0(const G4double t)
{
  G4double result = (std::exp(t) - 1.0);
  return result;
}

G4double G4GEMProbability::I1(const G4double t, const G4double tx)
{
  G4double result = t - tx + 1.0;
  result *= std::exp(tx);
  result -= (t + 1.0);
  return result;
}


G4double G4GEMProbability::I2(const G4double s, const G4double sx)
{
  G4double S = 1.0/std::sqrt(s);
  G4double Sx = 1.0/std::sqrt(sx);
  
  G4double p1 = S*S*S*( 1.0 + S*S*( 1.5 + 3.75*S*S) );
  G4double p2 = Sx*Sx*Sx*( 1.0 + Sx*Sx*( 1.5 + 3.75*Sx*Sx) )*std::exp(sx-s);
  
  return p1-p2;
}

G4double G4GEMProbability::I3(const G4double s, const G4double sx)
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

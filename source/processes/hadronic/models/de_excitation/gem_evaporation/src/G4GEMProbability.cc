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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept 2001)
//


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
    //    G4double ResidualZ = static_cast<G4double>(fragment.GetZ() - theZ);
    G4double U = fragment.GetExcitationEnergy();
	
    G4double NuclearMass = G4ParticleTable::GetParticleTable()->
        GetIonTable()->GetNucleusMass(theZ,theA);

    G4double a = theEvapLDPptr->LevelDensityParameter(static_cast<G4int>(fragment.GetA()),
						      static_cast<G4int>(fragment.GetZ()),U);
    G4double delta0 = G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(fragment.GetA()),
									       static_cast<G4int>(fragment.GetZ()));
    
    G4double Alpha = CalcAlphaParam(fragment);
    G4double Beta = CalcBetaParam(fragment);
    
    
    G4double Ux = (2.5 + 150.0/ResidualA)*MeV;
    G4double Ex = Ux + delta0;
    G4double T = 1.0/(std::sqrt(a/Ux) - 1.5/Ux);
    G4double E0 = Ex - T*(std::log(T*MeV) - std::log(a/MeV)/4.0 - 1.25*std::log(Ux*MeV) + 2.0*std::sqrt(a*Ux));
    
    G4double Width;
    G4double InitialLevelDensity;
    if ( MaximalKineticEnergy < Ex ) {
        G4double t = MaximalKineticEnergy/T;
        Width = (I1(t,t) + (Beta+V)*I0(t))/std::exp(E0/T);
        InitialLevelDensity = (pi/12.0)*std::exp((U-E0)/T)/T;
    } else {
        G4double t = MaximalKineticEnergy/T;
        G4double tx = Ex/T;
        G4double s = 2.0*std::sqrt(a*(MaximalKineticEnergy-delta0));
        G4double sx = 2.0*std::sqrt(a*(Ex-delta0));
        Width = I1(t,tx)/std::exp(E0/T) + I3(s,sx)*std::exp(s)/(std::sqrt(2.0)*a);
        // For charged particles (Beta+V) = 0 beacuse Beta = -V
        if (theZ == 0) {
            Width += (Beta+V)*(I0(tx)/std::exp(E0/T) + 2.0*std::sqrt(2.0)*I2(s,sx)*std::exp(s));
        }
        InitialLevelDensity = (pi/12.0)*std::exp(2*std::sqrt(a*(U-delta0)))/std::pow(a*std::pow(U-delta0,5.0),1.0/4.0);
    }
    

    G4double g = (2.0*Spin+1.0)*NuclearMass/(pi2* hbar_Planck*hbar_Planck);

    G4double RN = 0.0;
    if (theA > 4) 
    {
	G4double R1 = std::pow(ResidualA,1.0/3.0);
	G4double R2 = std::pow(G4double(theA),1.0/3.0);
	RN = 1.12*(R1 + R2) - 0.86*((R1+R2)/(R1*R2));
	RN *= fermi;
    }
    else 
    {
	RN = 1.5*fermi;
    }
    G4double GeometricalXS = pi*RN*RN*std::pow(ResidualA,2./3.); 


    Width *= std::sqrt(pi)*g*GeometricalXS*Alpha/(12.0*InitialLevelDensity); 
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

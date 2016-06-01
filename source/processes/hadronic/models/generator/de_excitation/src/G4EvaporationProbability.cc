// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4EvaporationProbability.hh"




G4EvaporationProbability::G4EvaporationProbability(const G4EvaporationProbability &right)
{
 G4Exception("G4EvaporationProbability::copy_constructor meant to not be accessable");
}




const G4EvaporationProbability & G4EvaporationProbability::
operator=(const G4EvaporationProbability &right)
{
  G4Exception("G4EvaporationProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4EvaporationProbability::operator==(const G4EvaporationProbability &right) const
{
  return false;
}

G4bool G4EvaporationProbability::operator!=(const G4EvaporationProbability &right) const
{
  return true;
}


G4double G4EvaporationProbability::EmissionProbability(const G4Fragment & fragment, const G4double photonExcitation)
  // Calculate integrated probability (width) for rvaporation channel:
  // If fragment has A_f <= 4 it will be used Dostrovsky's
  // approximation for the inverse reaction cross section.  If
  // fragment has A_f > 4 it will be used Botvina's approximation for
  // the inverse reaction cross section.
{
  // first af all a test
  if (theChannel->GetMaximalKineticEnergy() <= 0.0 || fragment.GetExcitationEnergy() <= 0.0) return 0.0;
  
  // We take decision on which approximation we'll use.
  if (theChannel->GetA() <= 4) return DostrovskyApproximation(fragment.GetA(),fragment.GetExcitationEnergy());
  else return BotvinaApproximation(fragment.GetA(),fragment.GetExcitationEnergy());
}



G4double G4EvaporationProbability::DostrovskyApproximation(const G4int A, const G4double U)
  // Width for evaporation channel with Dostrovsky's approximation for
  // inverse cross section.  
{
  G4double SystemEntropy = 2.0*sqrt((theChannel->GetLevelDensityParameter()/(1./MeV)) *
				    A * U/MeV);
 // r0 -> Absorption Radius R=r0 A^(1/3)
 G4double r0 = 2.173*(1.0+0.006103*theChannel->GetZ()*theChannel->GetResidualZ())/
   (1.0+0.009443*theChannel->GetZ()*theChannel->GetResidualZ());   
 // compute the integrated probability of evaporation channel

 G4double RN = 1.5;
 G4double CC;
 if (theChannel->GetA() == 1) CC = 0.2; // neutron, proton
 // deuterium, triton, alpha, 
 else if ((theChannel->GetZ() == 1 && (theChannel->GetA() == 2 || theChannel->GetA() == 3)) ||
	  (theChannel->GetZ() == 2 && (theChannel->GetA() == 3 || theChannel->GetA() == 4))) 
   CC = 0.1;
 // He5, He6, Li5, Li6, ...., O17, O18
 else CC = pow(G4double(theChannel->GetA())/G4double(theChannel->GetResidualA()),2.0/3.0); 


 G4double ALFA;
 G4double BETA;
 if (theChannel->GetZ() == 0) {  // neutron 
   ALFA = 0.76+2.2/pow(theChannel->GetResidualA(),1.0/3.0);
   BETA = (2.12/pow(theChannel->GetResidualA(),2.0/3.0) - 0.05)/ALFA;
 } else {
   ALFA = 1.0 + CC;
   BETA = 0.0;
 }


 G4double Q1 = (theChannel->GetLevelDensityParameter()/(1./MeV)) * theChannel->GetResidualA();
 G4double Q2 = Q1*theChannel->GetMaximalKineticEnergy()/MeV;
 G4double Q3 = (theChannel->GetGamma()*pow(theChannel->GetResidualA(),2.0/3.0))*(ALFA/(Q1*Q1))*
   (G4double(theChannel->GetResidualA())/G4double(theChannel->GetResidualA()+theChannel->GetA()))*
   (pi*RN*RN)/(2.0*41.5*pi2);
 G4double Q4 = (2.0*BETA*Q1-3.0)/2.0 + Q2;
 G4double Q5 = (2.0*BETA*Q1-3.0)*(sqrt(Q2)-0.5)+2.0*Q2;

 G4double PEX1;
 if (SystemEntropy > 160.0) PEX1 = 0.0;
 else PEX1 = Q4*exp(-SystemEntropy);

 G4double PP2 = SystemEntropy - 2.0*sqrt(Q2);

 G4double PEX2;
 if (PP2 > 160.0) PEX2 = 0.0;
 else PEX2 = Q5*exp(-PP2);

 return Q3*(PEX1+PEX2);
}



G4double G4EvaporationProbability::BotvinaApproximation(const G4int A, const G4double U)
  // Width for evaporation channel with Botvina's approximation for
  // inverse cross section.
{

  G4double SystemEntropy = 2.0*sqrt((theChannel->GetLevelDensityParameter()/(1./MeV))*
				    A * U/MeV);
  // r0 -> Absorption Radius R=r0 A^(1/3)
  G4double r0 = 2.173*(1.0+0.006103*theChannel->GetZ()*theChannel->GetResidualZ())/
    (1.0+0.009443*theChannel->GetZ()*theChannel->GetResidualZ());   

  // compute the integrated probability of evaporation channel
  G4double DALF = 0.869+9.91/theChannel->GetResidualZ();

  G4double KinPlusCoul = theChannel->GetMaximalKineticEnergy()/MeV + theChannel->GetCoulombBarrier()/MeV;
  if (KinPlusCoul <= theChannel->GetCoulombBarrier()/(5.0*MeV) || theChannel->GetZ() == 0) return 0.0;
    
  G4double CC = pow(G4double(theChannel->GetA())/G4double(theChannel->GetResidualA()),2.0/3.0);

  G4double ALFA = 1.0+CC;

  G4double Q1 = (theChannel->GetLevelDensityParameter()/(1./MeV)) * theChannel->GetResidualA();

  G4double Q3 = theChannel->GetGamma() * pow(theChannel->GetResidualA(),2.0/3.0) * (ALFA/(Q1*Q1))
    *(G4double(theChannel->GetResidualA())/G4double(theChannel->GetResidualA()+theChannel->GetA()))*
    ((pi*r0*r0)/(2.0*41.5*pi2));

  G4double TempMKE = theChannel->GetMaximalKineticEnergy()/MeV - 1.0;
  G4double prob3 = 0.0;
  if (TempMKE > 0.0) {
    G4double Q2 = Q1*TempMKE;
    G4double Q4 = (2.0*Q1-3.0)/2.0 + Q2;
    G4double Q5 = (2.0*Q1-3.0)*(sqrt(Q2)-0.5)+2.0*Q2;
    prob3 = Q3*(Q4*exp(-SystemEntropy)+Q5*exp(2.0*sqrt(Q2)-SystemEntropy));
  }
  G4double EX = theChannel->GetCoulombBarrier()/MeV + 1.0;
  G4double EM = KinPlusCoul-Q1/(DALF*DALF);
  G4double prob4 = 0.0;
  G4double SQ = 0.0;
  G4double CSI = 0.0;
  G4double F1CSI = 0.0;
  if (EM >= EX) { 
    SQ = sqrt(Q1/(KinPlusCoul-EX));
    G4double F1X = DALF-SQ;
    G4double F2X = SQ/(2.0*(KinPlusCoul-EX));
    if (F1X >= (0.5*F2X)) {
      CSI = 0.693/F1X;
      G4double SQCSI = sqrt(Q1/(KinPlusCoul-EX+CSI));
      G4double F2CSI = SQCSI/(2.0*(KinPlusCoul-EX+CSI));
      prob4 = Q3*2.0*Q1*Q1*(1.0/F1X)*
	exp(2.0*sqrt(Q1*(KinPlusCoul-EX))-SystemEntropy-F2CSI*CSI*CSI/2.0);
    } else {
      CSI = 0.48/sqrt(0.5*F2X);
      F1CSI = DALF-sqrt(Q1/(KinPlusCoul-EX+CSI));
      prob4 = Q3*Q1*Q1*sqrt(6.2832/F2X)*
	exp(2.0*sqrt(Q1*(KinPlusCoul-EX))-SystemEntropy-F1CSI*CSI);
    }
  } else if (EM < EX && EM > theChannel->GetCoulombBarrier()/(5.0*MeV)) { 
    SQ = sqrt(Q1/(KinPlusCoul-EM));
    G4double F2M = SQ/(2.0*(KinPlusCoul-EM));
    CSI = 0.48/sqrt(0.5*F2M);
    F1CSI = DALF - sqrt(Q1/(KinPlusCoul-EM+CSI));
    prob4 = Q3*2.0*Q1*Q1*sqrt(6.2832/F2M)*
      exp(DALF*(EM-EX)+2.0*sqrt(Q1*(KinPlusCoul-EM))-SystemEntropy-F1CSI*CSI);
  } else return prob3;
  return prob3+prob4;
}


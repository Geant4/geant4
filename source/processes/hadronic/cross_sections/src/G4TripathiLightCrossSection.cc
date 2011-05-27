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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme).                     *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TripathiLightCrossSection.cc
//
// Version:		B.1
// Date:		15/04/04
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// J. M. Quesada 24 November 2010 bug fixed in X_m 
// (according to eq. 14 in
//     R.K. Tripathi et al. Nucl. Instr. and Meth. in Phys. Res. B 155 (1999) 349-356)
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4TripathiLightCrossSection.hh"
#include "G4WilsonRadius.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"


G4TripathiLightCrossSection::G4TripathiLightCrossSection ()
 : G4VCrossSectionDataSet("TripathiLightIons")
{
  // Constructor only needs to instantiate the object which provides functions
  // to calculate the nuclear radius, and some other constants used to
  // calculate cross-sections.

  theWilsonRadius = new G4WilsonRadius();
  r_0             = 1.1  * fermi;
  third           = 1.0/3.0;

  // The following variable is set to true if
  // G4TripathiLightCrossSection::GetCrossSection is going to be called from
  // within G4TripathiLightCrossSection::GetCrossSection to check whether the
  // cross-section is behaviing anomalously in the low-energy region.

  lowEnergyCheck = false;
}
///////////////////////////////////////////////////////////////////////////////
//
G4TripathiLightCrossSection::~G4TripathiLightCrossSection ()
{
//
//
// Destructor just needs to delete the pointer to the G4WilsonRadius object.
//
  delete theWilsonRadius;
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool G4TripathiLightCrossSection::IsApplicable
  (const G4DynamicParticle* theProjectile, const G4Element* theTarget)
{
  G4int Z = G4lrint(theTarget->GetZ());
  G4int A = G4lrint(theTarget->GetN());
  return IsIsoApplicable(theProjectile, Z, A);
}


G4bool
G4TripathiLightCrossSection::IsIsoApplicable(const G4DynamicParticle* theProjectile,
                                             G4int ZZ, G4int AA)
{
  G4bool result = false;
  const G4double AT = AA;
  const G4double ZT = ZZ;
  const G4double ZP = theProjectile->GetDefinition()->GetPDGCharge();
  const G4double AP = theProjectile->GetDefinition()->GetBaryonNumber();
  if (theProjectile->GetKineticEnergy()/
      theProjectile->GetDefinition()->GetBaryonNumber()<10.0*GeV &&
     ((AT==1 && ZT==1) || (AP==1 && ZP==1) ||
      (AT==1 && ZT==0) || (AP==1 && ZP==0) ||
      (AT==2 && ZT==1) || (AP==2 && ZP==1) ||
      (AT==3 && ZT==2) || (AP==3 && ZP==2) ||
      (AT==4 && ZT==2) || (AP==4 && ZP==2))) result = true;
  return result;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double
G4TripathiLightCrossSection::GetZandACrossSection(const G4DynamicParticle* theProjectile,
                                                  G4int ZZ, G4int AA, G4double /*theTemperature*/)
{
  // Initialise the result.
  G4double result = 0.0;

  // Get details of the projectile and target (nucleon number, atomic number,
  // kinetic enery and energy/nucleon.

  const G4double AT = AA;
  const G4double ZT = ZZ;
  const G4double EA = theProjectile->GetKineticEnergy()/MeV;
  const G4double AP = theProjectile->GetDefinition()->GetBaryonNumber();
  const G4double ZP = theProjectile->GetDefinition()->GetPDGCharge();
  G4double E  = EA / AP;

  // Determine target mass and energy within the centre-of-mass frame.

  G4double mT = G4ParticleTable::GetParticleTable()
                ->GetIonTable()
                ->GetIonMass(static_cast<G4int>(ZT), static_cast<G4int>(AT));
  G4LorentzVector pT(0.0, 0.0, 0.0, mT);
  G4LorentzVector pP(theProjectile->Get4Momentum());
  pT = pT + pP;
  G4double E_cm = (pT.mag()-mT-pP.m())/MeV;
  if (E_cm < DBL_MIN) return 0.;

  // Determine nuclear radii.  Note that the r_p and r_T are defined differently
  // from Wilson et al.

  G4WilsonRadius theWilsonNuclearRadius;
  G4double r_rms_p = theWilsonRadius->GetWilsonRMSRadius(AP);
  G4double r_rms_t = theWilsonRadius->GetWilsonRMSRadius(AT);

  G4double r_p = 1.29*r_rms_p;
  G4double r_t = 1.29*r_rms_t;

  G4double Radius = (r_p + r_t)/fermi + 1.2*(std::pow(AT, third) + std::pow(AP, third))/
    std::pow(E_cm, third);

  G4double B = 1.44 * ZP * ZT / Radius;

  // Now determine other parameters associated with the parametric
  // formula, depending upon the projectile and target.

  G4double T1 = 0.0;
  G4double D  = 0.0;
  G4double G  = 0.0;

  if ((AT==1 && ZT==1) || (AP==1 && ZP==1)) {
    T1 = 23.0;
    D  = 1.85 + 0.16/(1+std::exp((500.0-E)/200.0));

  } else if ((AT==1 && ZT==0) || (AP==1 && ZP==0)) {
    T1 = 18.0;
    D  = 1.85 + 0.16/(1+std::exp((500.0-E)/200.0));

  } else if ((AT==2 && ZT==1) || (AP==2 && ZP==1)) {
    T1 = 23.0;
    D  = 1.65 + 0.1/(1+std::exp((500.0-E)/200.0));

  } else if ((AT==3 && ZT==2) || (AP==3 && ZP==2)) {
    T1 = 40.0;
    D  = 1.55;

  } else if (AP==4 && ZP==2) {
    if      (AT==4 && ZT==2) {T1 = 40.0; G = 300.0;}
    else if (ZT==4)          {T1 = 25.0; G = 300.0;}
    else if (ZT==7)          {T1 = 40.0; G = 500.0;}
    else if (ZT==13)         {T1 = 25.0; G = 300.0;}
    else if (ZT==26)         {T1 = 40.0; G = 300.0;}
    else                     {T1 = 40.0; G = 75.0;}
    D = 2.77 - 8.0E-3*AT + 1.8E-5*AT*AT-0.8/(1.0+std::exp((250.0-E)/G));
  }
  else if (AT==4 && ZT==2)
  {
    if      (AP==4 && ZP==2) {T1 = 40.0; G = 300.0;}
    else if (ZP==4)          {T1 = 25.0; G = 300.0;}
    else if (ZP==7)          {T1 = 40.0; G = 500.0;}
    else if (ZP==13)         {T1 = 25.0; G = 300.0;}
    else if (ZP==26)         {T1 = 40.0; G = 300.0;}
    else                     {T1 = 40.0; G = 75.0;}
    D = 2.77 - 8.0E-3*AP + 1.8E-5*AP*AP-0.8/(1.0+std::exp((250.0-E)/G));
  }

  // C_E, S, deltaE, X1, S_L and X_m correspond directly with the original
  // formulae of Tripathi et al in his report.

  G4double C_E = D*(1.0-std::exp(-E/T1)) -
                 0.292*std::exp(-E/792.0)*std::cos(0.229*std::pow(E,0.453));

  G4double S = std::pow(AP,third)*std::pow(AT,third)/(std::pow(AP,third) + std::pow(AT,third));

  G4double deltaE = 0.0;
  G4double X1     = 0.0;
  if (AT >= AP)
  {
    deltaE = 1.85*S + 0.16*S/std::pow(E_cm,third) - C_E + 0.91*(AT-2.0*ZT)*ZP/AT/AP;
    X1     = 2.83 - 3.1E-2*AT + 1.7E-4*AT*AT;
  }
  else
  {
    deltaE = 1.85*S + 0.16*S/std::pow(E_cm,third) - C_E + 0.91*(AP-2.0*ZP)*ZT/AT/AP;
    X1     = 2.83 - 3.1E-2*AP + 1.7E-4*AP*AP;
  }
  G4double S_L = 1.2 + 1.6*(1.0-std::exp(-E/15.0));
  //JMQ 241110 bug fixed 
  //  G4double X_m = 1.0 - X1*std::exp(-E/X1*S_L);
  G4double X_m = 1.0 - X1*std::exp(-E/(X1*S_L));

  // R_c is also highly dependent upon the A and Z of the projectile and
  // target.

  G4double R_c = 1.0;
  if (AP==1 && ZP==1)
  {
    if      (AT==2 && ZT==1) R_c = 13.5;
    else if (AT==3 && ZT==2) R_c = 21.0;
    else if (AT==4 && ZT==2) R_c = 27.0;
    else if (ZT==3)          R_c = 2.2;
  }
  else if (AT==1 && ZT==1)
  {
    if      (AP==2 && ZP==1) R_c = 13.5;
    else if (AP==3 && ZP==2) R_c = 21.0;
    else if (AP==4 && ZP==2) R_c = 27.0;
    else if (ZP==3)          R_c = 2.2;
  }
  else if (AP==2 && ZP==1)
  {
    if       (AT==2 && ZT==1) R_c = 13.5;
    else if (AT==4 && ZT==2)  R_c = 13.5;
    else if (AT==12 && ZT==6) R_c = 6.0;
  }
  else if (AT==2 && ZT==1)
  {
    if       (AP==2 && ZP==1) R_c = 13.5;
    else if (AP==4 && ZP==2)  R_c = 13.5;
    else if (AP==12 && ZP==6) R_c = 6.0;
  }
  else if ((AP==4 && ZP==2 && (ZT==73 || ZT==79)) ||
           (AT==4 && ZT==2 && (ZP==73 || ZP==79))) R_c = 0.6;

  // Find the total cross-section.  Check that it's value is positive, and if
  // the energy is less that 10 MeV/nuc, find out if the cross-section is
  // increasing with decreasing energy.  If so this is a sign that the function
  // is behaving badly at low energies, and the cross-section should be
  // set to zero.

  result = pi * r_0*r_0 *
           std::pow((std::pow(AT,third) + std::pow(AP,third) + deltaE),2.0) *
           (1.0 - R_c*B/E_cm) * X_m;
  if (!lowEnergyCheck)
  {
    if (result < 0.0)
      result = 0.0;
    else if (E < 6.0*MeV)
    {
      G4double f  = 0.95;
      G4DynamicParticle slowerProjectile = *theProjectile;
      slowerProjectile.SetKineticEnergy(f * EA * MeV);
      // G4TripathiLightCrossSection theTripathiLightCrossSection; // MHM 20090824 Not needed
	  // theTripathiLightCrossSection.SetLowEnergyCheck(true);
	  G4bool savelowenergy=lowEnergyCheck;
	  SetLowEnergyCheck(true);
      G4double resultp =
        GetZandACrossSection(&slowerProjectile, ZZ, AA, 0.0);
	  SetLowEnergyCheck(savelowenergy);
      if (resultp >result) result = 0.0;
    }
  }

  return result;
}


G4double G4TripathiLightCrossSection::GetCrossSection
  (const G4DynamicParticle* theProjectile, const G4Element* theTarget,
  G4double theTemperature)
{
  G4int nIso = theTarget->GetNumberOfIsotopes();
  G4double xsection = 0;
     
  if (nIso) {
    G4double sig;
    G4IsotopeVector* isoVector = theTarget->GetIsotopeVector();
    G4double* abundVector = theTarget->GetRelativeAbundanceVector();
    G4int ZZ;
    G4int AA;
     
    for (G4int i = 0; i < nIso; i++) {
      ZZ = (*isoVector)[i]->GetZ();
      AA = (*isoVector)[i]->GetN();
      sig = GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
      xsection += sig*abundVector[i];
    }
   
  } else {
    G4int ZZ = G4lrint(theTarget->GetZ());
    G4int AA = G4lrint(theTarget->GetN());
    xsection = GetZandACrossSection(theProjectile, ZZ, AA, theTemperature);
  }
    
  return xsection;
}


///////////////////////////////////////////////////////////////////////////////
//
void G4TripathiLightCrossSection::SetLowEnergyCheck (G4bool aLowEnergyCheck)
{
  lowEnergyCheck = aLowEnergyCheck;
}


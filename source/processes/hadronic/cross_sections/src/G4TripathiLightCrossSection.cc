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
// 24 November 2010 J. M. Quesada bug fixed in X_m 
// (according to eq. 14 in
//     R.K. Tripathi et al. Nucl. Instr. and Meth. in Phys. Res. B 155 (1999) 349-356)
//
// 19 Aug 2011 V.Ivanchenko move to new design and make x-section per element
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4TripathiLightCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4WilsonRadius.hh"
#include "G4NucleiProperties.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"
#include "G4Pow.hh"

G4TripathiLightCrossSection::G4TripathiLightCrossSection ()
 : G4VCrossSectionDataSet("TripathiLightIons")
{
  // Constructor only needs to instantiate the object which provides functions
  // to calculate the nuclear radius, and some other constants used to
  // calculate cross-sections.

  theWilsonRadius = new G4WilsonRadius();
  r_0             = 1.1  * fermi;

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
  // Destructor just needs to delete the pointer to the G4WilsonRadius object.
  //
  delete theWilsonRadius;
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool
G4TripathiLightCrossSection::IsElementApplicable(const G4DynamicParticle* theProjectile,
						 G4int ZT, const G4Material*)
{
  G4bool result = false;
  G4int AT = G4lrint(G4NistManager::Instance()->GetAtomicMassAmu(ZT));
  G4int ZP = G4lrint(theProjectile->GetDefinition()->GetPDGCharge()/eplus);
  G4int AP = theProjectile->GetDefinition()->GetBaryonNumber();
  if (theProjectile->GetKineticEnergy()/AP < 10.0*GeV &&
      ((AT==1 && ZT==1) || (AP==1 && ZP==1) ||
       (AT==1 && ZT==0) || (AP==1 && ZP==0) ||
       (AT==2 && ZT==1) || (AP==2 && ZP==1) ||
       (AT==3 && ZT==2) || (AP==3 && ZP==2) ||
       (AT==4 && ZT==2) || (AP==4 && ZP==2))) { result = true; }
  return result;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double
G4TripathiLightCrossSection::GetElementCrossSection(const G4DynamicParticle* theProjectile,
						    G4int ZT, const G4Material*)
{
  // Initialise the result.
  G4double result = 0.0;

  // Get details of the projectile and target (nucleon number, atomic number,
  // kinetic enery and energy/nucleon.

  G4double xAT= G4NistManager::Instance()->GetAtomicMassAmu(ZT);
  G4int    AT = G4lrint(xAT);
  G4double EA = theProjectile->GetKineticEnergy()/MeV;
  G4int    AP = theProjectile->GetDefinition()->GetBaryonNumber();
  G4double xAP= G4double(AP);
  G4double ZP = G4lrint(theProjectile->GetDefinition()->GetPDGCharge()/eplus);
  G4double E  = EA / xAP;

  G4Pow* g4pow = G4Pow::GetInstance();
  
  G4double AT13 = g4pow->Z13(AT);
  G4double AP13 = g4pow->Z13(AP);

  // Determine target mass and energy within the centre-of-mass frame.

  G4double mT = G4NucleiProperties::GetNuclearMass(AT, ZT);
  G4LorentzVector pT(0.0, 0.0, 0.0, mT);
  G4LorentzVector pP(theProjectile->Get4Momentum());
  pT += pP;
  G4double E_cm = (pT.mag()-mT-pP.m())/MeV;

  //G4cout << G4endl;
  //G4cout << "### EA= " << EA << " ZT= " << ZT << " AT= " << AT 
  //	 << "  ZP= " << ZP << " AP= " << AP << " E_cm= " << E_cm 
  //	 << " Elim= " << (0.8 + 0.04*ZT)*xAP << G4endl;

  if (E_cm <= 0.0) { return 0.; }
  if (E_cm <= (0.8 + 0.04*ZT)*xAP && !lowEnergyCheck) { return 0.; }
  
  G4double E_cm13 = g4pow->A13(E_cm);

  // Determine nuclear radii.  Note that the r_p and r_T are defined differently
  // from Wilson et al.

  G4double r_rms_p = theWilsonRadius->GetWilsonRMSRadius(xAP);
  G4double r_rms_t = theWilsonRadius->GetWilsonRMSRadius(xAT);

  G4double r_p = 1.29*r_rms_p;
  G4double r_t = 1.29*r_rms_t;

  G4double Radius = (r_p + r_t)/fermi + 1.2*(AT13 + AP13)/E_cm13;

  G4double B = 1.44 * ZP * ZT / Radius;

  // Now determine other parameters associated with the parametric
  // formula, depending upon the projectile and target.

  G4double T1 = 0.0;
  G4double D  = 0.0;
  G4double G  = 0.0;

  if ((AT==1 && ZT==1) || (AP==1 && ZP==1)) {
    T1 = 23.0;
    D  = 1.85 + 0.16/(1+G4Exp((500.0-E)/200.0));

  } else if ((AT==1 && ZT==0) || (AP==1 && ZP==0)) {
    T1 = 18.0;
    D  = 1.85 + 0.16/(1+G4Exp((500.0-E)/200.0));

  } else if ((AT==2 && ZT==1) || (AP==2 && ZP==1)) {
    T1 = 23.0;
    D  = 1.65 + 0.1/(1+G4Exp((500.0-E)/200.0));

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
    D = 2.77 - 8.0E-3*AT + 1.8E-5*AT*AT-0.8/(1.0+G4Exp((250.0-E)/G));
  }
  else if (AT==4 && ZT==2) {
    if      (AP==4 && ZP==2) {T1 = 40.0; G = 300.0;}
    else if (ZP==4)          {T1 = 25.0; G = 300.0;}
    else if (ZP==7)          {T1 = 40.0; G = 500.0;}
    else if (ZP==13)         {T1 = 25.0; G = 300.0;}
    else if (ZP==26)         {T1 = 40.0; G = 300.0;}
    else                     {T1 = 40.0; G = 75.0;}
    D = 2.77 - 8.0E-3*AP + 1.8E-5*AP*AP-0.8/(1.0+G4Exp((250.0-E)/G));
  }

  // C_E, S, deltaE, X1, S_L and X_m correspond directly with the original
  // formulae of Tripathi et al in his report.
  //G4cout << "E= " << E << " T1= " << T1 << "  AP= " << AP << " ZP= " << ZP 
  //	 << " AT= " << AT << " ZT= " << ZT << G4endl;
  G4double C_E = D*(1.0-G4Exp(-E/T1)) -
                 0.292*G4Exp(-E/792.0)*std::cos(0.229*G4Pow::GetInstance()->powA(E,0.453));

  G4double S = AP13*AT13/(AP13 + AT13);

  G4double deltaE = 0.0;
  G4double X1     = 0.0;
  if (AT >= AP)
  {
    deltaE = 1.85*S + 0.16*S/E_cm13 - C_E + 0.91*(AT-2*ZT)*ZP/(xAT*xAP);
    X1     = 2.83 - 3.1E-2*AT + 1.7E-4*AT*AT;
  }
  else
  {
    deltaE = 1.85*S + 0.16*S/E_cm13 - C_E + 0.91*(AP-2*ZP)*ZT/(xAT*xAP);
    X1     = 2.83 - 3.1E-2*AP + 1.7E-4*AP*AP;
  }
  G4double S_L = 1.2 + 1.6*(1.0-G4Exp(-E/15.0));
  //JMQ 241110 bug fixed 
  G4double X_m = 1.0 - X1*G4Exp(-E/(X1*S_L));

  //G4cout << "deltaE= " << deltaE << "  X1= " << X1 << " S_L= " << S_L << " X_m= " << X_m << G4endl;

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

  G4double xr = r_0*(AT13 + AP13 + deltaE);
  result = pi * xr * xr * (1.0 - R_c*B/E_cm) * X_m;
  //G4cout << "       result= " << result << " E= " << E << "  check= "<< lowEnergyCheck << G4endl;
  if (result < 0.0) {
    result = 0.0;

  } else if (!lowEnergyCheck && E < 6.0) {
    G4double f  = 0.95;
    G4DynamicParticle slowerProjectile = *theProjectile;
    slowerProjectile.SetKineticEnergy(f * EA * MeV);

    G4bool savelowenergy = lowEnergyCheck;
    SetLowEnergyCheck(true);
    G4double resultp = GetElementCrossSection(&slowerProjectile, ZT);
    SetLowEnergyCheck(savelowenergy);
    //G4cout << "           newres= " << resultp << " f*EA= " << f*EA << G4endl;  
    if (resultp > result) { result = 0.0; }
  }

  return result;
}


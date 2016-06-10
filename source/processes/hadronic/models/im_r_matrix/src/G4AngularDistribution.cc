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
// hpw: done, but low quality at present.

#include "globals.hh"
#include "G4Log.hh"
#include "G4SystemOfUnits.hh"
#include "G4AngularDistribution.hh"
#include "Randomize.hh"

G4AngularDistribution::G4AngularDistribution(G4bool symmetrize)
  : sym(symmetrize)
{
  // The following are parameters of the model - not to be confused with the PDG values!

      mSigma = 0.55; 
      cmSigma = 1.20; 
      gSigma  = 9.4;

      mOmega = 0.783; 
      cmOmega = 0.808; 
      gOmega = 10.95;

      mPion = 0.138; 
      cmPion = 0.51; 
      gPion  = 7.27;

      mNucleon = 0.938;   

      // Definition of constants for pion-Term (no s-dependence)

      m42   = 4. * mNucleon * mNucleon;
      mPion2  = mPion * mPion;           
      cmPion2 = cmPion * cmPion;
      dPion1 = cmPion2-mPion2;
      dPion2 = dPion1 * dPion1;
      cm6gp = 1.5 * (cmPion2*cmPion2*cmPion2) * (gPion*gPion*gPion*gPion) * m42 * m42 / dPion2;

      cPion_3 = -(cm6gp/3.);
      cPion_2 = -(cm6gp * mPion2/dPion1);
      cPion_1 = -(cm6gp * mPion2 * (2. * cmPion2 + mPion2) / dPion2);
      cPion_m = -(cm6gp * cmPion2 * mPion2 / dPion2);
      cPion_L = -(cm6gp * 2. * cmPion2 * mPion2 * (cmPion2 + mPion2) / dPion2 / dPion1);
      cPion_0 = -(cPion_3 + cPion_2 + cPion_1 + cPion_m);
    
     // Definition of constants for sigma-Term (no s-dependence)

      G4double gSigmaSq = gSigma * gSigma; 

      mSigma2  = mSigma * mSigma;
      cmSigma2 = cmSigma * cmSigma;
      cmSigma4 = cmSigma2 * cmSigma2;
      cmSigma6 = cmSigma2 * cmSigma4;
      dSigma1 = m42 - cmSigma2;
      dSigma2 = m42 - mSigma2;
      dSigma3 = cmSigma2 - mSigma2;

      G4double dSigma1Sq = dSigma1 * dSigma1;
      G4double dSigma2Sq = dSigma2 * dSigma2;
      G4double dSigma3Sq = dSigma3 * dSigma3;

      cm2gs = 0.5 * cmSigma2 * gSigmaSq*gSigmaSq / dSigma3Sq;     


      cSigma_3 = -(cm2gs * dSigma1Sq / 3.);
      cSigma_2 = -(cm2gs * cmSigma2 * dSigma1 * dSigma2 / dSigma3);
      cSigma_1 = -(cm2gs * cmSigma4 * (2. * dSigma1 + dSigma2) * dSigma2 / dSigma3Sq);
      cSigma_m = -(cm2gs * cmSigma6 * dSigma2Sq / mSigma2 / dSigma3Sq);
      cSigma_L = -(cm2gs * cmSigma6 * dSigma2 * (dSigma1 + dSigma2) * 2. / (dSigma3 * dSigma3Sq));
      cSigma_0 = -(cSigma_3 + cSigma_2 + cSigma_1 + cSigma_m);

      // Definition of constants for omega-Term

      G4double gOmegaSq = gOmega * gOmega;

      mOmega2  = mOmega * mOmega;
      cmOmega2 = cmOmega * cmOmega;
      cmOmega4 = cmOmega2 * cmOmega2;
      cmOmega6 = cmOmega2 * cmOmega4;
      dOmega1 = m42 - cmOmega2;
      dOmega2 = m42 - mOmega2;
      dOmega3 = cmOmega2 - mOmega2;
      sOmega1 = cmOmega2 + mOmega2;

      G4double dOmega3Sq = dOmega3 * dOmega3;

      cm2go = 0.5 * cmOmega2 * gOmegaSq * gOmegaSq / dOmega3Sq;

      cOmega_3 =  cm2go / 3.;
      cOmega_2 = -(cm2go * cmOmega2 / dOmega3);
      cOmega_1 =  cm2go * cmOmega4 / dOmega3Sq;
      cOmega_m =  cm2go * cmOmega6 / (dOmega3Sq * mOmega2);
      cOmega_L = -(cm2go * cmOmega6 * 4. / (dOmega3 * dOmega3Sq));

      // Definition of constants for mix-Term

      G4double fac1Tmp = (gSigma * gOmega * cmSigma2 * cmOmega2);
      fac1 = -(fac1Tmp *  fac1Tmp * m42);  
      dMix1 = cmOmega2 - cmSigma2;
      dMix2 = cmOmega2 - mSigma2;
      dMix3 = cmSigma2 - mOmega2;

      G4double dMix1Sq = dMix1 * dMix1;
      G4double dMix2Sq = dMix2 * dMix2;
      G4double dMix3Sq = dMix3 * dMix3;

      cMix_o1 =    fac1 / (cmOmega2  * dMix1Sq * dMix2 * dOmega3);
      cMix_s1 =    fac1 / (cmSigma2  * dMix1Sq * dMix3 * dSigma3);
      cMix_Omega = fac1 / (dOmega3Sq * dMix3Sq * (mOmega2 - mSigma2));
      cMix_sm =    fac1 / (dSigma3Sq * dMix2Sq * (mSigma2 - mOmega2)); 
      fac2 = (-fac1) / (dMix1*dMix1Sq * dOmega3Sq * dMix2Sq);
      fac3 = (-fac1) / (dMix1*dMix1Sq * dSigma3Sq * dMix3Sq); 
      
      cMix_oLc = fac2 * (3. * cmOmega2*cmOmega4            - cmOmega4 * cmSigma2 
		       - 2. * cmOmega4 * mOmega2           - 2. * cmOmega4 * mSigma2 
		       + cmOmega2 * mOmega2 * mSigma2      + cmSigma2 * mOmega2 * mSigma2 
		       - 4. * cmOmega4 * m42               + 2. * cmOmega2 * cmSigma2 * m42 
		       + 3. * cmOmega2 * mOmega2 * m42     - cmSigma2 * mOmega2 * m42 
		       + 3. * cmOmega2 * mSigma2 * m42     - cmSigma2 * mSigma2 * m42 
		       - 2. * mOmega2 * mSigma2 * m42);
      
      cMix_oLs = fac2 * (8. * cmOmega4                     - 4. * cmOmega2 * cmSigma2 
			 - 6. * cmOmega2 * mOmega2         + 2. * cmSigma2 * mOmega2 
			 - 6. * cmOmega2 * mSigma2         + 2. * cmSigma2 * mSigma2 
			 + 4. * mOmega2 * mSigma2);
      
      cMix_sLc = fac3 * (cmOmega2 * cmSigma4               - 3. * cmSigma6    
			 + 2. * cmSigma4 * mOmega2         + 2. * cmSigma4 * mSigma2 
			 - cmOmega2 * mOmega2 * mSigma2    - cmSigma2 * mOmega2 * mSigma2 
			 - 2. * cmOmega2 * cmSigma2 * m42  + 4. * cmSigma4 * m42 
			 + cmOmega2 * mOmega2 * m42        - 3. * cmSigma2 * mOmega2 * m42 
			 + cmOmega2 * mSigma2 * m42        - 3. * cmSigma2 * mSigma2 * m42 
			 + 2. * mOmega2 * mSigma2 * m42);

      cMix_sLs = fac3 * (4. * cmOmega2 * cmSigma2          - 8. * cmSigma4
		       - 2. * cmOmega2 * mOmega2           + 6. * cmSigma2 * mOmega2 
		       - 2. * cmOmega2 * mSigma2           + 6. * cmSigma2 * mSigma2 
		       - 4. * mOmega2 * mSigma2);
}


G4AngularDistribution::~
G4AngularDistribution()
{ }


G4double G4AngularDistribution::CosTheta(G4double S, G4double m_1, G4double m_2) const
{
   G4double random = G4UniformRand();
   G4double dCosTheta = 2.;
   G4double cosTheta = -1.;

   // For jmax=12 the accuracy is better than 0.1 degree 
   G4int jMax = 12;

   for (G4int j = 1; j <= jMax; ++j)
   {
      // Accuracy is 2^-jmax 
      dCosTheta *= 0.5;
      G4double cosTh = cosTheta + dCosTheta;
      if(DifferentialCrossSection(S, m_1, m_2, cosTh) <= random) cosTheta = cosTh;
    }

   // Randomize in final interval in order to avoid discrete angles 
   cosTheta += G4UniformRand() * dCosTheta;


   if (cosTheta > 1. || cosTheta < -1.)
     throw G4HadronicException(__FILE__, __LINE__, "G4AngularDistribution::CosTheta - std::cos(theta) outside allowed range");

   return cosTheta;
}


G4double G4AngularDistribution::DifferentialCrossSection(G4double sIn, G4double m_1, G4double m_2, 
							 G4double cosTheta) const
{
// local calculus is in GeV, ie. normalize input
  sIn = sIn/sqr(GeV)+m42/2.;
  m_1  = m_1/GeV;
  m_2  = m_2/GeV;
//  G4cout << "Here we go"<<sIn << " "<<m1 << " " << m2 <<" " m42<< G4endl;
// scaling from masses other than p,p.
  G4double S = sIn - (m_1+m_2) * (m_1+m_2) + m42;
  G4double tMax = S - m42;
  G4double tp = 0.5 * (cosTheta + 1.) * tMax; 
  G4double twoS = 2. * S;
  
  // Define s-dependent stuff for omega-Term
  G4double brak1 = (twoS-m42) * (twoS-m42);
  G4double bOmega_3 = cOmega_3 * (-2. * cmOmega4 - 2. * cmOmega2 * twoS - brak1);
  G4double bOmega_2 = cOmega_2 * ( 2. * cmOmega2 * mOmega2 + sOmega1 * twoS + brak1);
  G4double bOmega_1 = cOmega_1 * (-4. * cmOmega2 * mOmega2 
			 - 2. * mOmega2*mOmega2 
			 - 2. * (cmOmega2 + 2 * mOmega2) * twoS 
			 - 3. * brak1);
  G4double bOmega_m = cOmega_m * (-2. * mOmega2*mOmega2 - 2. * mOmega2 * twoS - brak1);
  G4double bOmega_L = cOmega_L * (sOmega1 * mOmega2 + (cmOmega2 + 3. * mOmega2) * S + brak1);
  G4double bOmega_0 = -(bOmega_3 + bOmega_2 + bOmega_1 + bOmega_m);
  
  // Define s-dependent stuff for mix-Term            
  G4double bMix_o1 = cMix_o1 * (dOmega1 - twoS);
  G4double bMix_s1 = cMix_s1 * (dSigma1 - twoS);
  G4double bMix_Omega = cMix_Omega * (dOmega2 - twoS);
  G4double bMix_sm = cMix_sm * (dSigma2 - twoS);
  G4double bMix_oL = cMix_oLc + cMix_oLs * S;
  G4double bMix_sL = cMix_sLc + cMix_sLs * S;
  
  G4double t1_Pion = 1. / (1. + tMax / cmPion2);
  G4double t2_Pion = 1. + tMax / mPion2;
  G4double t1_Sigma = 1. / (1. + tMax / cmSigma2);
  G4double t2_Sigma = 1. + tMax / mSigma2;
  G4double t1_Omega = 1. / (1. + tMax / cmOmega2);
  G4double t2_Omega = 1. + tMax / mOmega2;
  
  G4double norm = Cross(t1_Pion, t1_Sigma, t1_Omega,
			t2_Pion, t2_Sigma, t2_Omega,
			bMix_o1, bMix_s1, bMix_Omega,
			bMix_sm, bMix_oL, bMix_sL,
			bOmega_0, bOmega_1, bOmega_2,
			bOmega_3, bOmega_m, bOmega_L);
  
  t1_Pion = 1. / (1. + tp / cmPion2);
  t2_Pion = 1. + tp / mPion2;
  t1_Sigma = 1. / (1. + tp / cmSigma2);
  t2_Sigma = 1. + tp / mSigma2;
  t1_Omega = 1. / (1. + tp / cmOmega2);
  t2_Omega = 1. + tp / mOmega2;
  
  G4double dSigma;
  if (sym) 
    { 
      G4double to;
      norm = 2. * norm;
      to = tMax - tp;
      G4double t3_Pion = 1. / (1. + to / cmPion2);
      G4double t4_Pion = 1. + to / mPion2;
      G4double t3_Sigma = 1. / (1. + to / cmSigma2);
      G4double t4_Sigma = 1. + to / mSigma2;
      G4double t3_Omega = 1. / (1. + to / cmOmega2);
      G4double t4_Omega = 1. + to / mOmega2;
      
      dSigma = ( Cross(t1_Pion, t1_Sigma, t1_Omega,
		       t2_Pion,t2_Sigma, t2_Omega,
		       bMix_o1, bMix_s1, bMix_Omega,
		       bMix_sm, bMix_oL, bMix_sL,
		       bOmega_0, bOmega_1, bOmega_2,
		       bOmega_3, bOmega_m, bOmega_L) - 
		 Cross(t3_Pion,t3_Sigma, t3_Omega,
		       t4_Pion, t4_Sigma, t4_Omega,
		       bMix_o1, bMix_s1, bMix_Omega,
		       bMix_sm, bMix_oL, bMix_sL,
		       bOmega_0, bOmega_1, bOmega_2,
		       bOmega_3, bOmega_m, bOmega_L) ) 
	/ norm + 0.5;
    }
  else
    {
      dSigma = Cross(t1_Pion, t1_Sigma, t1_Omega, 
		     t2_Pion, t2_Sigma, t2_Omega,
		     bMix_o1, bMix_s1, bMix_Omega,
		     bMix_sm, bMix_oL, bMix_sL,
		     bOmega_0, bOmega_1, bOmega_2,
		     bOmega_3, bOmega_m, bOmega_L) 
	/ norm;
    }
  
  return dSigma;
}


G4double G4AngularDistribution::Cross(G4double tpPion,
				      G4double tpSigma,
				      G4double tpOmega,
				      G4double tmPion,
				      G4double tmSigma,
				      G4double tmOmega,
				      G4double bMix_o1,
				      G4double bMix_s1,
				      G4double bMix_Omega,
				      G4double bMix_sm,
				      G4double bMix_oL,
				      G4double bMix_sL,
				      G4double bOmega_0,
				      G4double bOmega_1,
				      G4double bOmega_2,
				      G4double bOmega_3,
				      G4double bOmega_m,
				      G4double bOmega_L) const
{
  G4double cross = 0;
     //  Pion
    cross += ((cPion_3 * tpPion  + cPion_2)  * tpPion  + cPion_1)  * tpPion  + cPion_m/tmPion   + cPion_0  +  cPion_L * G4Log(tpPion*tmPion);
//    G4cout << "cross1 "<< cross<<G4endl;
    //  Sigma
    cross += ((cSigma_3 * tpSigma + cSigma_2) * tpSigma + cSigma_1) * tpSigma + cSigma_m/tmSigma + cSigma_0 + cSigma_L * G4Log(tpSigma*tmSigma);
//    G4cout << "cross2 "<< cross<<G4endl;
    // Omega
    cross += ((bOmega_3 * tpOmega + bOmega_2) * tpOmega + bOmega_1) * tpOmega + bOmega_m/tmOmega + bOmega_0 + bOmega_L * G4Log(tpOmega*tmOmega)
    // Mix
    +  bMix_o1 * (tpOmega - 1.)
    +  bMix_s1 * (tpSigma - 1.)
    +  bMix_Omega * G4Log(tmOmega)
    +  bMix_sm * G4Log(tmSigma)
    +  bMix_oL * G4Log(tpOmega)
    +  bMix_sL * G4Log(tpSigma);
/*      G4cout << "cross3 "<< cross<<" "
             <<bMix_o1<<" "
             <<bMix_s1<<" "
             <<bMix_Omega<<" "
             <<bMix_sm<<" "
             <<bMix_oL<<" "
             <<bMix_sL<<" "
             <<tpOmega<<" "
             <<tpSigma<<" "
             <<tmOmega<<" "
             <<tmSigma<<" "
             <<tpOmega<<" "
             <<tpSigma
             <<G4endl;
*/
  return cross;

}

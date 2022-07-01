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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// This file contains the source code of various functions all related to the
// calculation of microroughness.
//
// see A. Steyerl, Z. Physik 254 (1972) 169.
//
// A description of the functions can be found in the corresponding header file
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UCNMicroRoughnessHelper.hh"

#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4UCNMicroRoughnessHelper* G4UCNMicroRoughnessHelper::fpInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor

G4UCNMicroRoughnessHelper::G4UCNMicroRoughnessHelper() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UCNMicroRoughnessHelper::~G4UCNMicroRoughnessHelper()
{
  fpInstance = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UCNMicroRoughnessHelper* G4UCNMicroRoughnessHelper::GetInstance()
{
  if(fpInstance == nullptr)
  {
    fpInstance = new G4UCNMicroRoughnessHelper;
  }
  return fpInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4UCNMicroRoughnessHelper::S2(G4double costheta2, G4double klk2) const
{
  // case 1: radicand is positive,
  // case 2: radicand is negative, cf. p. 174 of the Steyerl paper

  if(costheta2 >= klk2)
  {
    return 4 * costheta2 /
           (2 * costheta2 - klk2 +
            2 * std::sqrt(costheta2 * (costheta2 - klk2)));
  }
  else
  {
    return std::norm(2 * std::sqrt(costheta2) /
                     (std::sqrt(costheta2) +
                      std::sqrt(std::complex<G4double>(costheta2 - klk2))));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4UCNMicroRoughnessHelper::SS2(G4double costhetas2, G4double klks2) const
{
  // cf. p. 175 of the Steyerl paper

  return 4*costhetas2/
                   (2*costhetas2+klks2+2*std::sqrt(costhetas2*(costhetas2+klks2)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UCNMicroRoughnessHelper::Fmu(G4double k2, G4double thetai,
                                        G4double thetao, G4double phio,
                                        G4double b2, G4double w2,
                                        G4double AngCut) const
{
  G4double mu_squared;

  // Checks if the distribution is peaked around the specular direction

  if((std::fabs(thetai - thetao) < AngCut) && (std::fabs(phio) < AngCut))
  {
    mu_squared=0.;
  }
  else
  {
    // cf. p. 177 of the Steyerl paper

    G4double sinthetai = std::sin(thetai);
    G4double sinthetao = std::sin(thetao);
    mu_squared         = k2 * (sinthetai * sinthetai + sinthetao * sinthetao -
                       2. * sinthetai * sinthetao * std::cos(phio));
    }

  // cf. p. 177 of the Steyerl paper

  return b2*w2/twopi*std::exp(-mu_squared*w2/2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4UCNMicroRoughnessHelper::FmuS(G4double k, G4double kS,
                                G4double thetai, G4double thetaSo,
                                G4double phiSo,
                                G4double b2, G4double w2,
                                G4double AngCut, G4double thetarefract) const
{
  G4double mu_squared;

  // Checks if the distribution is peaked around the direction of
  // unperturbed refraction

  if((std::fabs(thetarefract - thetaSo) < AngCut) &&
     (std::fabs(phiSo) < AngCut))
  {
    mu_squared=0.;
  }
  else
  {
    G4double sinthetai  = std::sin(thetai);
    G4double sinthetaSo = std::sin(thetaSo);

    // cf. p. 177 of the Steyerl paper
    mu_squared = k * k * sinthetai * sinthetai +
                 kS * kS * sinthetaSo * sinthetaSo -
                 2. * k * kS * sinthetai * sinthetaSo * std::cos(phiSo);
    }

  // cf. p. 177 of the Steyerl paper

  return b2*w2/twopi*std::exp(-mu_squared*w2/2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4UCNMicroRoughnessHelper::IntIplus(G4double E, G4double fermipot,
                                    G4double theta_i, G4int AngNoTheta,
                                    G4int AngNoPhi, G4double b2,
                                    G4double w2, G4double* max,
                                    G4double AngCut) const
{
  *max=0.;

  // max_theta_o saves the theta-position of the max probability,
  // the previous value is saved in a_max_theta_o

  G4double a_max_theta_o, max_theta_o=theta_i, a_max_phi_o, max_phi_o=0.;

  // max_phi_o saves the phi-position of the max probability,
  // the previous value is saved in a_max_phi_o

  // Definition of the stepsizes in theta_o and phi_o

  G4double theta_o;
  G4double phi_o;
  G4double Intens;
  G4double ang_steptheta=90.*degree/(AngNoTheta-1);
  G4double ang_stepphi=360.*degree/(AngNoPhi-1);
  G4double costheta_i=std::cos(theta_i);
  G4double costheta_i_squared=costheta_i*costheta_i;

  // (k_l/k)^2
  G4double kl4d4=neutron_mass_c2/hbarc_squared*neutron_mass_c2/
                                 hbarc_squared*fermipot*fermipot;

  // (k_l/k)^2
  G4double klk2=fermipot/E;

  G4double costheta_o_squared;

  // k^2
  G4double k2=2*neutron_mass_c2*E/hbarc_squared;

  G4double wkeit=0.;

  // Loop through theta_o
 
  for (theta_o=0.*degree; theta_o<=90.*degree+1e-6; theta_o+=ang_steptheta)
    {
      costheta_o_squared=std::cos(theta_o)*std::cos(theta_o);

      // Loop through phi_o

      for (phi_o=-180.*degree; phi_o<=180.*degree+1e-6; phi_o+=ang_stepphi)
	{
          //calculates probability for a certain theta_o,phi_o pair

	  Intens=kl4d4/costheta_i*S2(costheta_i_squared,klk2)*
                 S2(costheta_o_squared,klk2)*
                 Fmu(k2,theta_i,theta_o,phi_o,b2,w2,AngCut)*std::sin(theta_o);

          //G4cout << "S2:  " << S2(costheta_i_squared,klk2) << G4endl;
          //G4cout << "S2:  " << S2(costheta_o_squared,klk2) << G4endl;
          //G4cout << "Fmu: " << Fmu(k2,theta_i,theta_o,phi_o,b2,w2,AngCut)*sin(theta_o) << G4endl;
          // checks if the new probability is larger than the
          // maximum found so far

          if (Intens>*max)
	    {
	      *max=Intens;
	      max_theta_o=theta_o;
	      max_phi_o=phi_o;
	    }

          // Adds the newly found probability to the integral probability

	  wkeit+=Intens*ang_steptheta*ang_stepphi;
	}
    }

  // Fine-Iteration to find maximum of distribution
  // only if the energy is not zero

  if (E>1e-10*eV)
    {

  // Break-condition for refining

  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while ((ang_stepphi>=AngCut*AngCut) || (ang_steptheta>=AngCut*AngCut))
    {
      a_max_theta_o=max_theta_o;
      a_max_phi_o=max_phi_o;
      ang_stepphi /= 2.;
      ang_steptheta /= 2.;

      //G4cout << ang_stepphi/degree << ", "
      //       << ang_steptheta/degree << ","
      //       << AngCut/degree << G4endl;

      for (theta_o=a_max_theta_o-ang_steptheta;
           theta_o<=a_max_theta_o-ang_steptheta+1e-6;
           theta_o+=ang_steptheta)
	{
	  //G4cout << "theta_o: " << theta_o/degree << G4endl;
	  costheta_o_squared=std::cos(theta_o)*std::cos(theta_o);
	  for (phi_o=a_max_phi_o-ang_stepphi;
               phi_o<=a_max_phi_o+ang_stepphi+1e-6;
               phi_o+=ang_stepphi)
	    {
	      //G4cout << "phi_o: " << phi_o/degree << G4endl;
	      Intens=kl4d4/costheta_i*S2(costheta_i_squared, klk2)*
                     S2(costheta_o_squared,klk2)*
                     Fmu(k2,theta_i,theta_o,phi_o,b2,w2,AngCut)*std::sin(theta_o);
	      if (Intens>*max)
		{
		  *max=Intens;
		  max_theta_o=theta_o;
		  max_phi_o=phi_o;
		}
	    }
	}
    }
    }
  return wkeit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UCNMicroRoughnessHelper::IntIminus(G4double E, G4double fermipot,
                                              G4double theta_i,
                                              G4int AngNoTheta,
                                              G4int AngNoPhi, G4double b2,
                                              G4double w2, G4double* max,
                                              G4double AngCut) const
{
  G4double a_max_thetas_o, max_thetas_o = theta_i;
  G4double a_max_phis_o, max_phis_o = 0.;
  G4double thetas_o;
  G4double phis_o;
  G4double IntensS;
  G4double ang_steptheta=180.*degree/(AngNoTheta-1);
  G4double ang_stepphi=180.*degree/(AngNoPhi-1);
  G4double costheta_i=std::cos(theta_i);
  G4double costheta_i_squared=costheta_i*costheta_i;

  *max = 0.;
  G4double wkeit=0.;

  if(E < fermipot)
  {
    return wkeit;
  }

  //k_l^4/4
  G4double kl4d4=neutron_mass_c2/hbarc_squared*neutron_mass_c2/
                 hbarc_squared*fermipot*fermipot;
  // (k_l/k)^2
  G4double klk2=fermipot/E;

  // (k_l/k')^2
  G4double klks2=fermipot/(E-fermipot);

  // k'/k
  G4double ksdk=std::sqrt((E-fermipot)/E);

  G4double costhetas_o_squared;

  // k
  G4double k=std::sqrt(2*neutron_mass_c2*E/hbarc_squared);

  // k'
  G4double kS=ksdk*k;

  for (thetas_o=0.*degree; thetas_o<=90.*degree+1e-6; thetas_o+=ang_steptheta)
    {
      costhetas_o_squared=std::cos(thetas_o)*std::cos(thetas_o);

      for (phis_o=-180.*degree; phis_o<=180.*degree+1e-6; phis_o+=ang_stepphi)
	{
          //cf. Steyerl-paper p. 176, eq. 21
	  if (costhetas_o_squared>=-klks2) {

            //calculates probability for a certain theta'_o, phi'_o pair

            G4double thetarefract = thetas_o;
            if(std::fabs(std::sin(theta_i) / ksdk) <= 1.)
            {
              thetarefract = std::asin(std::sin(theta_i) / ksdk);
            }

      IntensS = kl4d4/costheta_i*ksdk*S2(costheta_i_squared, klk2)*
                      SS2(costhetas_o_squared,klks2)*
                      FmuS(k,kS,theta_i,thetas_o,phis_o,b2,w2,AngCut,thetarefract)*
                      std::sin(thetas_o);
	  } else {
	    IntensS=0.;
          }
            // checks if the new probability is larger than
            // the maximum found so far
	  if (IntensS>*max)
	    {
	      *max=IntensS;
	    }
	  wkeit+=IntensS*ang_steptheta*ang_stepphi;
	}
    }

  // Fine-Iteration to find maximum of distribution

  if (E>1e-10*eV)
    {

  // Break-condition for refining

  while (ang_stepphi>=AngCut*AngCut || ang_steptheta>=AngCut*AngCut)
    {
      a_max_thetas_o=max_thetas_o;
      a_max_phis_o=max_phis_o;
      ang_stepphi /= 2.;
      ang_steptheta /= 2.;
      //G4cout << ang_stepphi/degree << ", " << ang_steptheta/degree 
      //       << ", " << AngCut/degree << G4endl;
      for (thetas_o=a_max_thetas_o-ang_steptheta;
           thetas_o<=a_max_thetas_o-ang_steptheta+1e-6;
           thetas_o+=ang_steptheta)
	{
	  costhetas_o_squared=std::cos(thetas_o)*std::cos(thetas_o);
	  for (phis_o=a_max_phis_o-ang_stepphi;
               phis_o<=a_max_phis_o+ang_stepphi+1e-6;
               phis_o+=ang_stepphi)
	    {
              G4double thetarefract = thetas_o;
              if(std::fabs(std::sin(theta_i) / ksdk) <= 1.)
              {
                thetarefract = std::asin(std::sin(theta_i) / ksdk);
              }

        IntensS=kl4d4/costheta_i*ksdk*S2(costheta_i_squared, klk2)*
                      SS2(costhetas_o_squared,klks2)*
                      FmuS(k,kS,theta_i,thetas_o,phis_o,b2,w2,AngCut,thetarefract)*
                      std::sin(thetas_o);
	      if (IntensS>*max)
		{
		  *max=IntensS;
		  max_thetas_o=thetas_o;
		  max_phis_o=phis_o;
		}
	    }
	}
    }
    }
  return wkeit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UCNMicroRoughnessHelper::ProbIplus(G4double E, G4double fermipot,
                                              G4double theta_i,
                                              G4double theta_o,
                                              G4double phi_o,
                                              G4double b, G4double w,
                                              G4double AngCut) const
{
  //k_l^4/4
  G4double kl4d4=neutron_mass_c2/hbarc_squared*neutron_mass_c2/
                 hbarc_squared*fermipot*fermipot;

  // (k_l/k)^2
  G4double klk2=fermipot/E;

  G4double costheta_i=std::cos(theta_i); 
  G4double costheta_o=std::cos(theta_o);

  // eq. 20 on page 176 in the steyerl paper

  return kl4d4/costheta_i*S2(costheta_i*costheta_i, klk2)*
         S2(costheta_o*costheta_o,klk2)*
   Fmu(2*neutron_mass_c2*E/hbarc_squared,theta_i,theta_o,phi_o,b*b,w*w,AngCut)*
   std::sin(theta_o);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UCNMicroRoughnessHelper::ProbIminus(G4double E, G4double fermipot,
                                               G4double theta_i,
                                               G4double thetas_o,
                                               G4double phis_o, G4double b,
                                               G4double w, G4double AngCut) const
{
  //k_l^4/4
  G4double kl4d4=neutron_mass_c2/hbarc_squared*neutron_mass_c2/
                 hbarc_squared*fermipot*fermipot;
  // (k_l/k)^2
  G4double klk2=fermipot/E;

  // (k_l/k')^2
  G4double klks2=fermipot/(E-fermipot);

  if (E < fermipot) {
     G4cout << " ProbIminus E < fermipot " << G4endl;
     return 0.;
  }

  // k'/k
  G4double ksdk=std::sqrt((E-fermipot)/E);

  // k
  G4double k=std::sqrt(2*neutron_mass_c2*E/hbarc_squared);

  // k'/k
  G4double kS=ksdk*k;

  G4double costheta_i=std::cos(theta_i);
  G4double costhetas_o=std::cos(thetas_o);

  // eq. 20 on page 176 in the steyerl paper

  G4double thetarefract = thetas_o;
  if(std::fabs(std::sin(theta_i) / ksdk) <= 1.)
  {
    thetarefract = std::asin(std::sin(theta_i) / ksdk);
  }

  return kl4d4/costheta_i*ksdk*S2(costheta_i*costheta_i, klk2)*
         SS2(costhetas_o*costhetas_o,klks2)*
         FmuS(k,kS,theta_i,thetas_o,phis_o,b*b,w*w,AngCut,thetarefract)*
         std::sin(thetas_o);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

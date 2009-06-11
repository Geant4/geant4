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
// $Id: G4FinalStateElasticChampion.cc,v 1.10 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------


#include "G4FinalStateElasticChampion.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4FinalStateElasticChampion::G4FinalStateElasticChampion()
{
  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  
  lowEnergyLimit = 8.23* eV; // SI : i/o of 7. * eV;
  highEnergyLimit = 10. * keV;
  
  G4double scaleFactor = 1e-16*cm*cm;

  char *path = getenv("G4LEDATA");
 
  if (!path)
    G4Exception("G4FinalStateElasticChampion: constructor: G4LEDATA environment variable not set");

  if (electronDef != 0)
  {
    std::ostringstream eFullFileName;
    eFullFileName << path << "/dna/sigmadiff_elastic_e_champion.dat";
    std::ifstream eDiffCrossSection(eFullFileName.str().c_str());
     
    if (!eDiffCrossSection) G4Exception("G4FinalStateElasticChampion: constructor: error opening electron DATA FILE");
      
    eTdummyVec.push_back(0.);

    while(!eDiffCrossSection.eof())
    {
	  double tDummy;
	  double eDummy;
	  eDiffCrossSection>>tDummy>>eDummy;

	  // SI : mandatory eVecm initialization
          if (tDummy != eTdummyVec.back()) 
          { 
            eTdummyVec.push_back(tDummy); 
            eVecm[tDummy].push_back(0.);
          }
	  
          eDiffCrossSection>>eDiffCrossSectionData[0][tDummy][eDummy];

	  // SI : only if not end of file reached !
          if (!eDiffCrossSection.eof()) eDiffCrossSectionData[0][tDummy][eDummy]*=scaleFactor;
	  
          if (eDummy != eVecm[tDummy].back()) eVecm[tDummy].push_back(eDummy);
          
    }

  }
  else
  {
      G4Exception("G4FinalStateElastiChampion : constructor: electron is not defined");
  }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4FinalStateElastiChampion is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4FinalStateElasticChampion::~G4FinalStateElasticChampion()
{ 
  eVecm.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4FinalStateProduct& G4FinalStateElasticChampion::GenerateFinalState(const G4Track& track, const G4Step& )
{
  product.Clear();

  G4double k = track.GetDynamicParticle()->GetKineticEnergy();

  if (k>= lowEnergyLimit && k < highEnergyLimit)
  {  
    G4double cosTheta = RandomizeCosTheta(k);
  
    G4double phi = 2. * pi * G4UniformRand();

    G4ThreeVector zVers = track.GetDynamicParticle()->GetMomentumDirection();
    G4ThreeVector xVers = zVers.orthogonal();
    G4ThreeVector yVers = zVers.cross(xVers);

    G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
    G4double yDir = xDir;
    xDir *= std::cos(phi);
    yDir *= std::sin(phi);

    G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));

    product.ModifyPrimaryParticle(zPrimeVers,k);
  }

  if (k<lowEnergyLimit)
  {
    product.KillPrimaryParticle();
  }
  
  return product;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4FinalStateElasticChampion::DifferentialCrossSection
  (G4ParticleDefinition * particleDefinition, G4double k, G4double theta) 							  
{
  G4double sigma = 0.;
  G4double valueT1 = 0;
  G4double valueT2 = 0;
  G4double valueE21 = 0;
  G4double valueE22 = 0;
  G4double valueE12 = 0;
  G4double valueE11 = 0;
  G4double xs11 = 0;   
  G4double xs12 = 0; 
  G4double xs21 = 0; 
  G4double xs22 = 0; 

  //SI : ensure the correct computation of cross section at the 180*deg limit
  if (theta==180.) theta=theta-1e-9;

  if (particleDefinition == G4Electron::ElectronDefinition()) 
  {
    std::vector<double>::iterator t2 = std::upper_bound(eTdummyVec.begin(),eTdummyVec.end(), k);
    std::vector<double>::iterator t1 = t2-1;
 
    std::vector<double>::iterator e12 = std::upper_bound(eVecm[(*t1)].begin(),eVecm[(*t1)].end(), theta);
    std::vector<double>::iterator e11 = e12-1;
	  
    std::vector<double>::iterator e22 = std::upper_bound(eVecm[(*t2)].begin(),eVecm[(*t2)].end(), theta);
    std::vector<double>::iterator e21 = e22-1;
	  	
    valueT1  =*t1;
    valueT2  =*t2;
    valueE21 =*e21;
    valueE22 =*e22;
    valueE12 =*e12;
    valueE11 =*e11;

    xs11 = eDiffCrossSectionData[0][valueT1][valueE11];
    xs12 = eDiffCrossSectionData[0][valueT1][valueE12];
    xs21 = eDiffCrossSectionData[0][valueT2][valueE21];
    xs22 = eDiffCrossSectionData[0][valueT2][valueE22];
  }
     
  G4double xsProduct = xs11 * xs12 * xs21 * xs22;
  
  if (xs11==0 || xs12==0 ||xs21==0 ||xs22==0) return (0.);
     
  if (xsProduct != 0.)
  {
    sigma = QuadInterpolator(  valueE11, valueE12, 
    			       valueE21, valueE22, 
			       xs11, xs12, 
			       xs21, xs22, 
			       valueT1, valueT2, 
			       k, theta );
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4FinalStateElasticChampion::LinLogInterpolate(G4double e1, 
						        G4double e2, 
						        G4double e, 
						        G4double xs1, 
						        G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = std::exp(d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4FinalStateElasticChampion::LogLogInterpolate(G4double e1, 
						        G4double e2, 
						        G4double e, 
						        G4double xs1, 
						        G4double xs2)
{
  G4double a = (std::log10(xs2)-std::log10(xs1)) / (std::log10(e2)-std::log10(e1));
  G4double b = std::log10(xs2) - a*std::log10(e2);
  G4double sigma = a*std::log10(e) + b;
  G4double value = (std::pow(10.,sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4FinalStateElasticChampion::QuadInterpolator(G4double e11, G4double e12, 
						       G4double e21, G4double e22, 
						       G4double xs11, G4double xs12, 
						       G4double xs21, G4double xs22, 
						       G4double t1, G4double t2, 
						       G4double t, G4double e)
{
// Log-Log
/*
  G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
*/

// Lin-Log
  G4double interpolatedvalue1 = LinLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4FinalStateElasticChampion::RandomizeCosTheta(G4double k) 
{

 G4double cosTheta = -1;
 G4double cumul = 0;
 G4double value = 0;
 
 // Number of integration steps in the 0-180 deg range
 G4int iMax=180;
 
 G4double random = G4UniformRand();
 
 // Cumulate differential cross section
 for (G4int i=0; i<iMax; i++) 
 {
   cumul = cumul + DifferentialCrossSection(G4Electron::ElectronDefinition(),k/eV,G4double(i)*180./(iMax-1));
 }
 
 // Select theta angle
 for (G4int i=0; i<iMax; i++) 
 {
   value = value + DifferentialCrossSection(G4Electron::ElectronDefinition(),k/eV,G4double(i)*180./(iMax-1)) / cumul;
   if (random < value) 
   {
     cosTheta=std::cos( deg*G4double(i)*180./(iMax-1) ); 
     break;
   } 
 } 

 return cosTheta;
}

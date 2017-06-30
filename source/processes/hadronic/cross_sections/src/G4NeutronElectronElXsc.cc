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
// 16.05.17   V. Grichine first implementation


#include "G4NeutronElectronElXsc.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"
// #include "G4Integrator.hh"

#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTable.hh"

#include "G4Neutron.hh"
#include "G4Electron.hh"


using namespace std;
using namespace CLHEP;

G4NeutronElectronElXsc::G4NeutronElectronElXsc()
 : G4VCrossSectionDataSet("NuElectronCcXsc")
{
  // neutron magneton squared

  fM   = neutron_mass_c2; // neutron mass
  fM2  = fM*fM;
  fme  = electron_mass_c2;
  fme2 = fme*fme;
  fMv2 = 0.7056*GeV*GeV;
  fee = fme;
  fee2 = fee*fee;
  fAm = 0.001;

  fCofXsc  = pi*fine_structure_const*fine_structure_const*hbarc*hbarc;
  fCofXsc  *= 3.6481;  // neutron Fm^2(0)
  fCofXsc /= fM*fM;

  // G4cout<<"fCofXsc = "<<fCofXsc*GeV/cm2<<" cm2/GeV"<<G4endl;

  // G4cout<<"hbarc = "<<hbarc/MeV/fermi<<" MeV*fermi"<<G4endl;

  fCutEnergy = 0.; // default value

  fEnergyBin = 200;
  fMinEnergy = 1.*MeV;
  fMaxEnergy = 10000.*GeV;

  fEnergyXscVector = new G4PhysicsLogVector(fMinEnergy, fMaxEnergy, fEnergyBin);

  for( G4int iTkin = 0; iTkin < fEnergyBin; iTkin++)  fEnergyXscVector->PutValue(iTkin, fXscArray[iTkin]*microbarn); 

  fBiasingFactor = 1.;

  // Initialise();
}

G4NeutronElectronElXsc::~G4NeutronElectronElXsc() 
{
  if( fEnergyXscVector ) 
  {
    delete fEnergyXscVector;
    fEnergyXscVector = 0;
  }
}

//////////////////////////////////////////////////////
//
// For neutrons in the precalculated energy interval

G4bool 
G4NeutronElectronElXsc::IsElementApplicable( const G4DynamicParticle* aPart, G4int, const G4Material*)
{
  G4bool result  = false;
  G4String pName = aPart->GetDefinition()->GetParticleName();
  G4double Tkin = aPart->GetKineticEnergy();  

  if( pName == "neutron" && 
      Tkin >= fMinEnergy && 
      Tkin <= fMaxEnergy    ) result = true;
  
  return result;
}

//////////////////////////////////////////////////

void  G4NeutronElectronElXsc::Initialise()
{
  G4int iTkin;
  G4double Tkin, rosxsc, xsc, delta, err=1.e-5;
  const G4ThreeVector mDir = G4ThreeVector(0.,0.,1.);
  const G4ParticleDefinition* pD = G4Neutron::Neutron();
  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");

  G4DynamicParticle dP;

  for( iTkin = 0; iTkin < fEnergyBin; iTkin++)
  {
    Tkin  = fEnergyXscVector->GetLowEdgeEnergy(iTkin);
    dP = G4DynamicParticle(pD, mDir, Tkin);   
    rosxsc = GetRosenbluthXsc(&dP, 1, mat);
    fEnergyXscVector->PutValue(iTkin, rosxsc);  // xscV.PutValue(evt, rosxsc); //
    xsc= fEnergyXscVector->Value(Tkin); // xsc= xscV.GetLowEdgeEnergy(evt); //
    delta = 0.5*std::abs( (rosxsc-xsc) )/(rosxsc+xsc);
    if(delta > err) G4cout<<Tkin/GeV<<" GeV, rosxsc = "<<rosxsc/microbarn<<"umb, v-xsc = "<<xsc<<" umb"<<G4endl;
  }
  return;
}

////////////////////////////////////////////////////

G4double G4NeutronElectronElXsc::
GetElementCrossSection(const G4DynamicParticle* aPart, G4int ZZ,  
		       const G4Material*) 
{
  G4double result = 0., Tkin;

  Tkin = aPart->GetKineticEnergy();

  result = fEnergyXscVector->Value(Tkin);

  result *= ZZ;  // incoherent sum over  all element electrons

  result *= fBiasingFactor;

  return result;
}

////////////////////////////////////////////////////
//
// Integration of the Rosenbluth differential xsc

G4double G4NeutronElectronElXsc::
GetRosenbluthXsc(const G4DynamicParticle* aPart, G4int ZZ,  
		       const G4Material*) 
{
  G4double result = 0., momentum;

  fee      = aPart->GetTotalEnergy()*fme/fM;
  fee2     = fee*fee;
  momentum = sqrt( fee2 - fme2 );
  fAm      = CalculateAm(momentum);

  // G4Integrator<G4NeutronElectronElXsc, G4double(G4NeutronElectronElXsc::*)(G4double)> integral;

  // result = integral.Legendre96( this, &G4NeutronElectronElXsc::XscIntegrand, 0., 1. );

  result *= fCofXsc;
 
  result *= ZZ;  // incoherent sum over  all element electrons

  return result;
}

/////////////////////////////////////////
//
// Rosenbluth relation in the neutron rest frame, 
// x = sin^2(theta/2), theta is the electron scattering angle
// Magnetic form factor in the dipole approximation.

G4double G4NeutronElectronElXsc::XscIntegrand(G4double x)
{
  G4double result = 1., q2, znq2, znf, znf2, znf4;

  znq2 = 1. + 2.*fee*x/fM;

  q2 = 4.*fee2*x/znq2;

  znf  = 1 + q2/fMv2;
  znf2 = znf*znf;
  znf4 = znf2*znf2;

  result /= ( x + fAm )*znq2*znq2*znf4; 

  result *= ( 1 - x )/( 1 + q2/4./fM2 ) + 2.*x;

  return result;
}

//////////////////////////////////////////////////////////
//
// Rosenbluth xsc in microbarn from 1*MeV to 10*Tev, 200 points

const G4double G4NeutronElectronElXsc::fXscArray[200] = {                            
1.52681, 1.54903, 1.57123, 1.59341, 1.61556, 1.63769, 1.6598, 1.68189, 1.70396, 
1.72601, 1.74805, 1.77007, 1.79208, 1.81407, 1.83605, 1.85801, 1.87997, 1.90192, 
1.92385, 1.94578, 1.96771, 1.98962, 2.01154, 2.03345, 2.05535, 2.07725, 2.09915, 
2.12105, 2.14295, 2.16485, 2.18675, 2.20865, 2.23055, 2.25244, 2.27433, 2.29621, 
2.31807, 2.33992, 2.36173, 2.38351, 2.40524, 2.42691, 2.4485, 2.47, 2.49138, 
2.51262, 2.53369, 2.55457, 2.57524, 2.59565, 2.61577, 2.63559, 2.65505, 2.67414, 
2.69281, 2.71104, 2.72881, 2.74607, 2.76282, 2.77903, 2.79467, 2.80974, 2.82422, 
2.83811, 2.85139, 2.86408, 2.87616, 2.88764, 2.89854, 2.90885, 2.91859, 2.92777, 
2.93641, 2.94453, 2.95213, 2.95924, 2.96588, 2.97207, 2.97782, 2.98316, 2.98811, 
2.99268, 2.9969, 3.00078, 3.00435, 3.00761, 3.01059, 3.01331, 3.01578, 3.01801, 
3.02003, 3.02185, 3.02347, 3.02491, 3.02619, 3.02732, 3.0283, 3.02915, 3.02988, 
3.03049, 3.03099, 3.03139, 3.03169, 3.03191, 3.03203, 3.03208, 3.03205, 3.03195, 
3.03177, 3.03152, 3.0312, 3.03081, 3.03034, 3.0298, 3.02919, 3.02849, 3.02771, 
3.02684, 3.02588, 3.02482, 3.02365, 3.02237, 3.02097, 3.01943, 3.01775, 3.0159, 
3.01389, 3.01169, 3.00929, 3.00666, 3.00379, 3.00065, 2.99722, 2.99347, 2.98936, 
2.98487, 2.97996, 2.97459, 2.9687, 2.96226, 2.9552, 2.94748, 2.93903, 2.92977, 
2.91965, 2.90858, 2.89649, 2.88329, 2.86889, 2.85321, 2.83615, 2.81764, 2.7976, 
2.77594, 2.7526, 2.72754, 2.70071, 2.67209, 2.64171, 2.60957, 2.57575, 2.54031, 
2.50336, 2.46504, 2.42548, 2.38484, 2.34328, 2.30099, 2.2581, 2.21478, 2.17115, 
2.12735, 2.08345, 2.03954, 1.99569, 1.95191, 1.90825, 1.86471, 1.82129, 1.77799, 
1.7348, 1.69171, 1.64869, 1.60575, 1.56286, 1.52, 1.47718, 1.43437, 1.39157, 
1.34877, 1.30596, 1.26314, 1.22031, 1.17746, 1.13459, 1.0917, 1.04879, 1.00585, 
0.962892, 0.919908 };

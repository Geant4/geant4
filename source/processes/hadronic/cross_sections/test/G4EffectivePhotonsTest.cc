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
#include "G4ios.hh"
#include "G4Electron.hh"
#include "G4Element.hh"
#include "G4ElectroNuclearCrossSection.hh"

int main()
{
  G4ElectroNuclearCrossSection eACrossSection;
  const G4int nN=6;                          // A#of plots (nuclei in the test)
  ofstream file0("h1_epl.out", ios::out);
  file0.setf( ios::scientific, ios::floatfield );
  ofstream file1("h1_eph.out", ios::out);
  file1.setf( ios::scientific, ios::floatfield );
  ofstream file2("cu_epl.out", ios::out);
  file2.setf( ios::scientific, ios::floatfield );
  ofstream file3("cu_eph.out", ios::out);
  file3.setf( ios::scientific, ios::floatfield );
  G4Element* theElement[nN]={new G4Element("Hydrogen", "H", 1, 1.*g/mole),
							 new G4Element("Deuterium", "D", 1, 2.*g/mole),
							 new G4Element("Carbon", "C", 6, 12.*g/mole),
							 new G4Element("Aluminum", "Al", 13, 27.*g/mole),
							 new G4Element("Copper", "Cu", 29, 64.*g/mole),
							 new G4Element("Lead", "Pb", 82, 208.*g/mole)};
  G4ParticleDefinition* theParticleDefinition = G4Electron::ElectronDefinition();
  G4DynamicParticle* theDynamicParticle;
  for(G4int n=0; n<4; n++)
  {
    G4double ekin=1000.;
    if(n==3) ekin=10000.;
    else if(n==1) ekin=10000000.;
    G4int   nel=0;
    if(n>1) nel=3;
    theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
											   G4ParticleMomentum(1.,0.,0.), ekin*MeV);
	G4double sig = eACrossSection.GetCrossSection(theDynamicParticle,theElement[nel])/millibarn;
	G4cout<<"n="<<n<<", Cross Section="<<sig<<G4endl;
	delete theDynamicParticle;
	G4int prob[100];
	for(G4int k=0; k<100; k++) prob[k]=0;
	for(G4int i=0; i<10000000; i++)
	{
	  G4double phe=eACrossSection.GetEffectivePhotonEnergy();
      //if(phe<100.)G4Exception("Too low photon energy");
	  G4int ind=static_cast<int>(phe/10.);
	  if(n==3) ind/=10;
	  else if(n==1) ind/=10000;
	  prob[ind]++;
	  if(!(i%1000000)) G4cout<<i<<" equivalent photons are generated"<<G4endl;
	}
	//G4cout<<"Now printing"<<G4endl;
	for(G4int j=0; j<100; j++)
	{
	  //G4cout<<"j="<<j<<G4endl;
	  if     (!n)   file0<<j*10    <<" "<<prob[j]<<G4endl;
	  else if(n==1) file1<<j*100000<<" "<<prob[j]<<G4endl;
	  else if(n==2) file2<<j*10    <<" "<<prob[j]<<G4endl;
	  else if(n==3) file3<<j*100   <<" "<<prob[j]<<G4endl;
	}
  }
  return EXIT_SUCCESS;
}







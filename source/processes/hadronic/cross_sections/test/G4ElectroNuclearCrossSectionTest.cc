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
#define debug
#define qedebug

//#include "G4UIterminal.hh"
#include "G4ios.hh"
//#include "G4ProcessManager.hh"
//#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Element.hh"
#include "G4ElectroNuclearCrossSection.hh"

int main()
{
  G4ElectroNuclearCrossSection eACrossSection;
  const G4int nN=6;                          // A#of plots (nuclei in the test)
  const G4int nC=227;                        // A#of points in the plots (energy values)
#ifdef qedebug
  const G4int nP=10;                         // A#of points in nu/Q2 for each eE (energy values)
#endif
  std::ofstream fileH1("h1_e.out", std::ios::out);
  fileH1.setf( std::ios::scientific, std::ios::floatfield );
  std::ofstream fileH2("h2_e.out", std::ios::out);
  fileH2.setf( std::ios::scientific, std::ios::floatfield );
  std::ofstream fileC12("c12_e.out", std::ios::out);
  fileC12.setf( std::ios::scientific, std::ios::floatfield );
  std::ofstream fileAl27("al27_e.out", std::ios::out);
  fileAl27.setf( std::ios::scientific, std::ios::floatfield );
  std::ofstream fileCu("cu_e.out", std::ios::out);
  fileCu.setf( std::ios::scientific, std::ios::floatfield );
  std::ofstream filePb("pb_e.out", std::ios::out);
  filePb.setf( std::ios::scientific, std::ios::floatfield );
  G4double low[nN]={137. , 2.22, 12.5 , 8.3 , 7.  , 7.35};
  G4double high[nN]={2.e7, 2.e7, 2.e7, 2.e7, 2.e7, 2.e7};
  G4int zA[nN]  ={1, 1,  6, 13, 29,  82};
  G4int aA[nN]  ={1, 2, 12, 27, 64, 208};
  G4Element* theElement[nN]={new G4Element("Hydrogen",  "H" , zA[0], aA[0]*g/mole),
							 new G4Element("Deuterium", "D" , zA[1], aA[1]*g/mole),
							 new G4Element("Carbon",    "C" , zA[2], aA[2]*g/mole),
							 new G4Element("Aluminum",  "Al", zA[3], aA[3]*g/mole),
							 new G4Element("Copper",    "Cu", zA[4], aA[4]*g/mole),
							 new G4Element("Lead",      "Pb", zA[5], aA[5]*g/mole)};
  G4ParticleDefinition* theParticleDefinition = G4Electron::ElectronDefinition();
  G4DynamicParticle* theDynamicParticle;
  G4double sig=0.;
  for(G4int n=0; n<nN; n++)                 // >>>> LOOP over Elements
  {
	G4double lekin = std::log(low[n]);
	G4double dlekin= std::exp((std::log(std::log(high[n]))-std::log(lekin))/(nC-1));
	lekin /= dlekin;
#ifdef debug
	G4cout<<"G4ElNucCSTest:>>>>n="<<n<<",l="<<low[n]<<",h="<<high[n]<<",n="<<nC<<",d="<<dlekin<<G4endl;
#endif
	for(G4int ll=0; ll<nC; ll++)            // >>>> LOOP over electron energies with the log-step
	{
	  lekin*=dlekin;
      G4double ekin=std::exp(lekin);
	  theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
												 G4ParticleMomentum(1.,0.,0.), ekin*MeV);
	  sig = eACrossSection.GetCrossSection(theDynamicParticle,theElement[n])/millibarn;
	  delete theDynamicParticle;
#ifdef debug
      G4cout<<"G4ElNucCSTest:****>>> Z="<<zA[n]<<", A="<<aA[n]<<", eE="<<ekin<<", sigma="<<sig<<G4endl;
#endif
#ifdef qedebug
	  if(sig>0.)for(G4int qe=0; qe<nP; qe++) // >>>> LOOP over equiv. photon energies for particular eE
	  {
	    G4double phE = eACrossSection.GetEquivalentPhotonEnergy();
	    G4double phQ2= eACrossSection.GetEquivalentPhotonQ2(phE);
	    G4double x = phQ2/1878./phE;
        G4cout<<"G4ElNucCSTest:.......eE="<<ekin<<", phE="<<phE<<", phQ2="<<phQ2<<", x="<<x<<G4endl;
	  } // End of the equivalent photon LOOP
#endif
      // Fill the result for the electro-nuclear cross section plot
      if(!n) fileH1<<ekin<<" "<<sig<<G4endl;
      else if(n==1) fileH2<<ekin<<" "<<sig<<G4endl;
      else if(n==2) fileC12<<ekin<<" "<<sig<<G4endl;
      else if(n==3) fileAl27<<ekin<<" "<<sig<<G4endl;
      else if(n==4) fileCu<<ekin<<" "<<sig<<G4endl;
      else if(n==5) filePb<<ekin<<" "<<sig<<G4endl;
    } // End of the LOOP over electron energies
  } // End of the LOOP over elements
  return EXIT_SUCCESS;
}

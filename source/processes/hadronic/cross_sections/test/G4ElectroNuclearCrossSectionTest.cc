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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#define debug
#define qedebug

#include "G4ios.hh"
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
  std::ofstream fileH1("h1_e.out", ios::out);
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
  G4double high[nN]={4.e7, 4.e7, 4.e7, 4.e7, 4.e7, 4.e7};
  G4Element* theElement[nN]={new G4Element("Hydrogen", "H", 1, 1.*g/mole),
							 new G4Element("Deuterium", "D", 1, 2.*g/mole),
							 new G4Element("Carbon", "C", 6, 12.*g/mole),
							 new G4Element("Aluminum", "Al", 13, 27.*g/mole),
							 new G4Element("Copper", "Cu", 29, 64.*g/mole),
							 new G4Element("Lead", "Pb", 82, 208.*g/mole)};
  G4ParticleDefinition* theParticleDefinition = G4Electron::ElectronDefinition();
  G4DynamicParticle* theDynamicParticle;
  G4double sig=0.;
  for(G4int n=0; n<6; n++)
  {
	G4double lekin = log(low[n]);
	G4double dlekin= exp((log(log(high[n]))-log(lekin))/(nC-1));
	lekin /= dlekin;
#ifdef debug
	G4cout<<"G4ElNucCSTest:>>>>>>n="<<n<<",l="<<low[n]<<",h="<<high[n]<<",n="<<nC<<",d="<<dlekin<<G4endl;
#endif
	for(G4int ll=0; ll<nC; ll++)
	{
	  lekin*=dlekin;
      G4double ekin=exp(lekin);
	  theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
												 G4ParticleMomentum(1.,0.,0.), ekin*MeV);
	  sig = eACrossSection.GetCrossSection(theDynamicParticle,theElement[n])/millibarn;
	  delete theDynamicParticle;
#ifdef debug
      G4cout<<"G4ElNucCSTest: n="<<n<<", e="<<ekin<<", sigma="<<sig<<G4endl;
#endif
#ifdef qedebug
	  for(G4int qe=0; qe<nP; qe++)
	  {
	    G4double phE = eACrossSection.GetEquivalentPhotonEnergy();
	    G4double phQ2= eACrossSection.GetEquivalentPhotonQ2(phE);
	    G4double x = phQ2/1878./phE;
        G4cout<<"G4ElNucCSTest:............eE="<<ekin<<", phE="<<phE<<", phQ2="<<phQ2<<", x="<<x<<G4endl;
	  }
#endif
      if(!n) fileH1<<ekin<<" "<<sig<<G4endl;
      else if(n==1) fileH2<<ekin<<" "<<sig<<G4endl;
      else if(n==2) fileC12<<ekin<<" "<<sig<<G4endl;
      else if(n==3) fileAl27<<ekin<<" "<<sig<<G4endl;
      else if(n==4) fileCu<<ekin<<" "<<sig<<G4endl;
      else if(n==5) filePb<<ekin<<" "<<sig<<G4endl;
    } // End of the point LOOP
  } // End of the plot LOOP
  return EXIT_SUCCESS;
}

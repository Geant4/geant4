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
#include "G4ios.hh"
#include "G4Electron.hh"
#include "G4Element.hh"
#include "G4ElectroNuclearCrossSection.hh"

int main()
{
  G4ElectroNuclearCrossSection eACrossSection;
  const G4int nN=6;                          // A#of plots (nuclei in the test)
  const G4int nC=227;                        // A#of points in the plots (energy values)
  ofstream fileH1("h1_e.out", ios::out);
  fileH1.setf( ios::scientific, ios::floatfield );
  ofstream fileH2("h2_e.out", ios::out);
  fileH2.setf( ios::scientific, ios::floatfield );
  ofstream fileC12("c12_e.out", ios::out);
  fileC12.setf( ios::scientific, ios::floatfield );
  ofstream fileAl27("al27_e.out", ios::out);
  fileAl27.setf( ios::scientific, ios::floatfield );
  ofstream fileCu("cu_e.out", ios::out);
  fileCu.setf( ios::scientific, ios::floatfield );
  ofstream filePb("pb_e.out", ios::out);
  filePb.setf( ios::scientific, ios::floatfield );
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
	//G4cout<<"n="<<n<<", low="<<low[n]<<", high="<<high[n]<<", np="<<nC<<", d="<<dlekin<<G4endl;
	for(G4int ll=0; ll<nC; ll++)
	{
	  lekin*=dlekin;
      G4double ekin=exp(lekin);
	  theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
												 G4ParticleMomentum(1.,0.,0.), ekin*MeV);
	  sig = eACrossSection.GetCrossSection(theDynamicParticle,theElement[n])/millibarn;
	  delete theDynamicParticle;
      //G4cout<<"n="<<n<<", e="<<ekin<<", sigma="<<sig<<G4endl;
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

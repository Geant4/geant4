//Ho preso la struttura dal test AtomicDexecitation.cc

#include "globals.hh"
#include "G4ios.hh"
#include <vector>
#include <iostream>
#include <fstream>
#include "G4VhShellCrossSection.hh"
#include "G4hShellCrossSection.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4Proton.hh"

int main()
{ 
  G4int Z; 
  G4double incidentEnergy;
  G4double mass;
  G4double deltaEnergy;
  size_t shellNumber;
  
  G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  G4VhShellCrossSection* shell = new G4hShellCrossSection;

  //massa del protone in Kg
  //mass = 1.67262158e-27 * kg;

  G4Proton* aProtone = G4Proton::Proton();
  mass = aProtone->GetPDGMass();

  std::vector<G4double> energies; 

  energies.push_back(0.005); // 5 KeV
  energies.push_back(0.010);
  energies.push_back(0.050);
  energies.push_back(0.100);
  energies.push_back(0.500);
  energies.push_back(1.000);
  energies.push_back(5.000);
  energies.push_back(10.00);
  energies.push_back(15.00);
  energies.push_back(50.00); // 50 MeV


  G4cout << "Enter shell number: " << G4endl;
  G4cin >> shellNumber;


  //Z is the atomic number
  for (Z = 1; Z<=92; Z++)
    { 
      G4cout << "Z = " << Z << G4endl;
      
      //Cross section for each incident energy
      for (G4int k=0; k<10;k++)
	{
	  incidentEnergy = energies[k]*MeV;
	 
	  // *************************************************************** //
	  // From Grizinsky, Phys.Rev. 138,2A, A305 & A322, deltaEnergy      //
	  // is the maximum energy of the delta ray emitted. So, deltaEnergy //
	  // will be proton energy - biding energy of the selected shell.    //
	  // *************************************************************** //

	  G4double bindingEnergy = transitionManager->Shell(Z,shellNumber)->BindingEnergy();

	  deltaEnergy = incidentEnergy - bindingEnergy;

	  std::vector<G4double> CS = shell->Probabilities(Z,incidentEnergy,mass,deltaEnergy);
	  
	  //Volendo scrivere sullo schermo	  
	  G4cout << " Incident Energy : " << incidentEnergy/MeV << " MeV -- ";
	  G4cout << "Cross Section = " << CS[shellNumber]/barn <<" barn " << G4endl; 


	  //Volendo scrivere i risultati in un file.

	  //         ofstream outfile( "mytest.out");
	  //         outfile << CS[1]/barn << " Barn" <<  G4endl;
	  //         outfile.close();

	  // 	   ofstream outfile( "mytest.out", std::ios::out );
	  //         fOut.write((G4double) (CS[1]), sizeof (G4int));
	  //            //Per verificare che il file si apra.
	  //            if( ! outfile) {
	  // 	     G4cerr << "Sorry! We anable to open the file.out!" << G4endl;
	  // 	     break;
	  // 	     // return -1;
	  // 	     }
	  // 	     mytest.out  << "Cross Section = " << CS[1]/barn << "Barn" <<  G4endl;	     
	}
    } 
  delete shell;

  G4cout<<"END OF THE MAIN PROGRAM"<<G4endl;
}

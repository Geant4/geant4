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

#include "G4FermiBreakUp.hh"

#include "G4NucleiProperties.hh"




int main()
{
  G4int A, Z;
  
  std::cout << "Please enter Z : ";
  std::cin >> Z;
  std::cout << "Please enter A : ";
  std::cin >> A;
  
  G4double energy;

  std::cout << "Please enter excitation energy (MeV) : ";
  std::cin >> energy;
  

  energy *= MeV;
  
  G4int iterations;
  std::cout << "Please enter number of interations : ";
  std::cin >> iterations;
  
  G4double NuclearMass = G4NucleiProperties::GetNuclearMass(A,Z) + energy;

  
  G4LorentzVector p4(0.0,0.0,0.0,NuclearMass);
  G4Fragment aFragment(A,Z,p4);
  
  G4FermiBreakUp model;

  std::cout << "Excited Nucleus\n";
  std::cout << aFragment << "\n\n";
  
  G4LorentzVector p4null(0.0,0.0,0.0,0.0);

  //  std::ios::fmtflags old_floatfield = std::cout.flags();
  std::cout.setf(std::ios::scientific,std::ios::floatfield);
  for (G4int i = 0; i < iterations; ++i)
    {
      G4LorentzVector p4cons;
      p4cons = p4null;
      std::cout << "############ event " << i+1 << "###########################\n";
      G4FragmentVector * theResult = model.BreakItUp(aFragment);
      std::cout << "Number of fragments: " << theResult->size() << '\n';
      std::cout << "-------------------------------------------\n";
      for (G4FragmentVector::iterator i = theResult->begin(); i != theResult->end(); ++i)
        {
          std::cout << (*i);
          std::cout << "\n-------------------------------------------\n";
          p4cons += (*i)->GetMomentum();
          delete (*i);
        }
      delete theResult;
      std::cout << "        Initial momentum : " << p4/MeV << '\n';
      std::cout << "Fragments total momentum : " << p4cons/MeV << '\n';
      std::cout << "            Conservation : " << (p4-p4cons)/MeV << '\n';
    }
  return 0;
}

      

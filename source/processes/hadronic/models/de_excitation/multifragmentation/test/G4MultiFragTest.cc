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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#include "G4StatMF.hh"

#include "G4NucleiProperties.hh"




int main()
{
  G4int A, Z;
  
  std::cout << "Please enter Z : ";
  std::cin >> Z;
  std::cout << "Please enter A : ";
  std::cin >> A;
  
  G4double energy;

  std::cout << "Please enter excitation energy per nucleon (MeV) : ";
  std::cin >> energy;
  

  energy *= static_cast<G4double>(A)*MeV;
  
  G4int iterations;
  std::cout << "Please enter number of interations : ";
  std::cin >> iterations;
  
  G4double NuclearMass = G4NucleiProperties::GetNuclearMass(A,Z) + energy;

  
  G4LorentzVector p4(0.0,0.0,0.0,NuclearMass);
  G4Fragment aFragment(A,Z,p4);
  
  G4StatMF model;

  std::cout << "Excited Nucleus\n";
  std::cout << aFragment << "\n\n";
  
  G4LorentzVector p4null(0.0,0.0,0.0,0.0);
  
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
      std::cout << "        Initial momentum : " << p4 << '\n';
      std::cout << "Fragments total momentum : " << p4cons << '\n';
      std::cout << "            Conservation : " << p4-p4cons << '\n';
    }
  return 0;
}

      

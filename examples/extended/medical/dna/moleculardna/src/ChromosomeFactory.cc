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
#include "ChromosomeFactory.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "VirtualChromosome.hh"
#include "CylindricalChromosome.hh"
#include "RodChromosome.hh"
#include "SphericalChromosome.hh"
#include "EllipticalChromosome.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VirtualChromosome* ChromosomeFactory::MakeChromosome(
  const G4String& name, const std::vector<G4String>& commands)
{
  VirtualChromosome* chromosome = nullptr;
  G4String chromosome_type      = commands[0];
  if(chromosome_type == CylindricalChromosome::fShape)
  {
    // interpret the command for cylinder
    // expect cyl rad height x y z unit rx ry rz
    // rotations are in degrees
    if(commands.size() == 7)
    {
      G4double unit = G4UnitDefinition::GetValueOf(commands[6]);
      G4double radius  = std::stod(commands[1]) * unit;
      G4double hgt  = std::stod(commands[2]) * unit;
      G4ThreeVector center(std::stod(commands[3]) * unit,
                           std::stod(commands[4]) * unit,
                           std::stod(commands[5]) * unit);
      chromosome = new CylindricalChromosome(name, center, radius, hgt);
    }
    else if(commands.size() == 10)  // euler angles given
    {
      G4double unit = G4UnitDefinition::GetValueOf(commands[6]);
      G4double radius  = std::stod(commands[1]) * unit;
      G4double hgt  = std::stod(commands[2]) * unit;
      G4ThreeVector center(std::stod(commands[3]) * unit,
                           std::stod(commands[4]) * unit,
                           std::stod(commands[5]) * unit);

      // Rotations are first X, then Y, Then Z
      G4RotationMatrix rot;
      rot.rotateX(std::stod(commands[7]) * pi / 180);
      rot.rotateY(std::stod(commands[8]) * pi / 180);
      rot.rotateZ(std::stod(commands[9]) * pi / 180);

      chromosome = new CylindricalChromosome(name, center, radius, hgt, rot);
    }
    else
    {
      G4cout << "The arguments for a cylinder are:" << G4endl;
      G4cout << "1)    name cyl rad height x y z unit" << G4endl;
      G4cout << "2)    name cyl rad height x y z unit rx ry rz" << G4endl;
      G4cout << "Note that rotations are in degrees" << G4endl;
      InvalidReading(chromosome_type);
    }
  }
  else if(chromosome_type == RodChromosome::fShape)
  {
    // interpret the command for rod
    // expect cyl rad height x y z unit rx ry rz
    // rotations are in degrees
    if(commands.size() == 7)
    {
      G4double unit = G4UnitDefinition::GetValueOf(commands[6]);
      G4double radius  = std::stod(commands[1]) * unit;
      G4double hgt  = std::stod(commands[2]) * unit;
      G4ThreeVector center(std::stod(commands[3]) * unit,
                           std::stod(commands[4]) * unit,
                           std::stod(commands[5]) * unit);
      chromosome = new RodChromosome(name, center, radius, hgt);
    }
    else if(commands.size() == 10)  // euler angles given
    {
      G4double unit = G4UnitDefinition::GetValueOf(commands[6]);
      G4double radius  = std::stod(commands[1]) * unit;
      G4double hgt  = std::stod(commands[2]) * unit;
      G4ThreeVector center(std::stod(commands[3]) * unit,
                           std::stod(commands[4]) * unit,
                           std::stod(commands[5]) * unit);

      // Rotations are first X, then Y, Then Z
      G4RotationMatrix rot;
      rot.rotateX(std::stod(commands[7]) * pi / 180);
      rot.rotateY(std::stod(commands[8]) * pi / 180);
      rot.rotateZ(std::stod(commands[9]) * pi / 180);

      chromosome = new RodChromosome(name, center, radius, hgt, rot);
    }
    else
    {
      G4cout << "The arguments for a cylinder are:" << G4endl;
      G4cout << "1)    name rod rad height x y z unit" << G4endl;
      G4cout << "2)    name rod rad height x y z unit rx ry rz" << G4endl;
      G4cout << "Note that rotations are in degrees" << G4endl;
      InvalidReading(chromosome_type);
    }
  }
  else if(chromosome_type == SphericalChromosome::fShape)
  {
    // interpret the command for sphere
    // expect sphere rad x y z unit rx ry rz
    // rotations are in degrees
    if(commands.size() == 6)
    {
      G4double unit = G4UnitDefinition::GetValueOf(commands[5]);
      G4double radius  = std::stod(commands[1]) * unit;
      G4ThreeVector center(std::stod(commands[2]) * unit,
                           std::stod(commands[3]) * unit,
                           std::stod(commands[4]) * unit);
      chromosome = new SphericalChromosome(name, center, radius);
    }
    else if(commands.size() == 9)  // euler angles given
    {
      G4double unit = G4UnitDefinition::GetValueOf(commands[5]);
      G4double radius  = std::stod(commands[1]) * unit;
      G4ThreeVector center(std::stod(commands[2]) * unit,
                           std::stod(commands[3]) * unit,
                           std::stod(commands[4]) * unit);

      // Rotations are first X, then Y, Then Z
      G4RotationMatrix rot;
      rot.rotateX(std::stod(commands[6]) * pi / 180);
      rot.rotateY(std::stod(commands[7]) * pi / 180);
      rot.rotateZ(std::stod(commands[8]) * pi / 180);

      chromosome = new SphericalChromosome(name, center, radius, rot);
    }
    else
    {
      G4cout << "The arguments for a sphere are:" << G4endl;
      G4cout << "1)    name sphere rad x y z unit" << G4endl;
      G4cout << "2)    name sphere rad x y z unit rx ry rz" << G4endl;
      G4cout << "Note that rotations are in degrees" << G4endl;
      InvalidReading(chromosome_type);
    }
  }
  else if(chromosome_type == EllipticalChromosome::fShape)
  {
    // interpret the command for Ellipse
    // expect ellipse sx sy sz x y z unit rx ry rz
    // rotations are in degrees
    // sx sy sz are semi-major axes
    if(commands.size() == 8)
    {
      G4double unit = G4UnitDefinition::GetValueOf(commands[7]);
      G4double sx   = std::stod(commands[1]) * unit;
      G4double sy   = std::stod(commands[2]) * unit;
      G4double sz   = std::stod(commands[3]) * unit;
      G4ThreeVector center(std::stod(commands[4]) * unit,
                           std::stod(commands[5]) * unit,
                           std::stod(commands[6]) * unit);
      chromosome = new EllipticalChromosome(name, center, sx, sy, sz);
    }
    else if(commands.size() == 11)  // euler angles given
    {
      G4double unit = G4UnitDefinition::GetValueOf(commands[7]);
      G4double sx   = std::stod(commands[1]) * unit;
      G4double sy   = std::stod(commands[2]) * unit;
      G4double sz   = std::stod(commands[3]) * unit;
      G4ThreeVector center(std::stod(commands[4]) * unit,
                           std::stod(commands[5]) * unit,
                           std::stod(commands[6]) * unit);

      // Rotations are first X, then Y, Then Z
      G4RotationMatrix rot;
      rot.rotateX(std::stod(commands[8]) * pi / 180);
      rot.rotateY(std::stod(commands[9]) * pi / 180);
      rot.rotateZ(std::stod(commands[10]) * pi / 180);

      chromosome = new EllipticalChromosome(name, center, sx, sy, sz, rot);
    }
    else
    {
      G4cout << "The arguments for a ellipse are:" << G4endl;
      G4cout << "1)    name ellipse sx sy sz x y z unit" << G4endl;
      G4cout << "2)    name ellipse sx sy sz x y z unit rx ry rz" << G4endl;
      G4cout << "Note that rotations are in degrees" << G4endl;
      G4cout << "Note that dimensions (sx, sy, sz) are semi-major axes"
             << G4endl;
      InvalidReading(chromosome_type);
    }
  }
  else
  {
    chromosome = nullptr;
    InvalidReading(chromosome_type);
  }
  return chromosome;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChromosomeFactory::InvalidReading(const G4String& chromosome_type)
{
  G4ExceptionDescription errmsg;
  errmsg << "Chromosome type: " << chromosome_type << " is not valid" << G4endl;
  G4Exception("ChromosomeFactory::MakeChromosome", "", FatalException, errmsg);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeFactory::Test()
{
  G4cout << "------------------------------------------------------------"
         << "--------------------" << G4endl;
  G4cout << "Chromosome Test" << G4endl;
  // populate vector for tests
  std::vector<std::vector<G4String>*> tests;
  G4double r;
  for(G4int ii = 0; ii != 50; ii++)
  {
    // No rotation constructor test
    auto vec = new std::vector<G4String>();
    r        = G4UniformRand();
    if(r < 0.25)
    {
      vec->push_back("rod");
      vec->push_back(std::to_string(10 * G4UniformRand() + 2));
      vec->push_back(std::to_string(5 * G4UniformRand() + 10));
    }
    else if(r < 0.5)
    {
      vec->push_back("sphere");
      vec->push_back(std::to_string(10 * G4UniformRand() + 2));
    }
    else if(r < 0.75)
    {
      vec->push_back("cyl");
      vec->push_back(std::to_string(10 * G4UniformRand() + 2));
      vec->push_back(std::to_string(5 * G4UniformRand() + 5));
    }
    else
    {
      vec->push_back("ellipse");
      vec->push_back(std::to_string(10 * G4UniformRand() + 2));
      vec->push_back(std::to_string(5 * G4UniformRand() + 5));
      vec->push_back(std::to_string(20 * G4UniformRand() - 6));
    }

    vec->push_back(std::to_string(10 * (G4UniformRand() - 0.5)));
    vec->push_back(std::to_string(10 * (G4UniformRand() - 0.5)));
    vec->push_back(std::to_string(10 * (G4UniformRand() - 0.5)));
    vec->push_back("um");
    tests.push_back(vec);
  }
  for(G4int ii = 0; ii != 50; ++ii)
  {
    // Test with rotation
    auto vec = new std::vector<G4String>();
    r        = G4UniformRand();
    if(r < 0.25)
    {
      vec->push_back("rod");
      vec->push_back(std::to_string(10 * G4UniformRand() + 2));
      vec->push_back(std::to_string(5 * G4UniformRand() + 10));
    }
    else if(r < 0.5)
    {
      vec->push_back("sphere");
      vec->push_back(std::to_string(10 * G4UniformRand() + 2));
    }
    else if(r < 0.75)
    {
      vec->push_back("cyl");
      vec->push_back(std::to_string(10 * G4UniformRand() + 2));
      vec->push_back(std::to_string(5 * G4UniformRand() + 5));
    }
    else
    {
      vec->push_back("ellipse");
      vec->push_back(std::to_string(10 * G4UniformRand() + 2));
      vec->push_back(std::to_string(5 * G4UniformRand() + 5));
      vec->push_back(std::to_string(20 * G4UniformRand() - 6));
    }

    vec->push_back(std::to_string(10 * (G4UniformRand() - 0.5)));
    vec->push_back(std::to_string(10 * (G4UniformRand() - 0.5)));
    vec->push_back(std::to_string(10 * (G4UniformRand() - 0.5)));

    vec->push_back("um");

    vec->push_back(std::to_string(720 * (G4UniformRand() - 0.5)));
    vec->push_back(std::to_string(720 * (G4UniformRand() - 0.5)));
    vec->push_back(std::to_string(720 * (G4UniformRand() - 0.5)));

    tests.push_back(vec);
  }

  VirtualChromosome* chromo;
  for(auto& test : tests)
  {
    chromo       = this->MakeChromosome("test", (*test));
    G4int passes = 0;
    G4int n      = 1000;
    for(G4int jj = 0; jj != n; jj++)
    {
      G4ThreeVector point = chromo->RandomPointInChromosome();
      if(chromo->PointInChromosome(point))
      {
        passes++;
      }
      else
      {
        G4cout << point << G4endl;
      }
    }
    if(passes == n)
    {
      G4cout << "Chromosome Test Passed for " << chromo->GetShape()
             << " based on " << n << " test points" << G4endl;
    }
    else
    {
      G4cout << "Chromosome Test Failed for " << chromo->GetShape()
             << " based on " << n << " test points (only " << passes << " pass)"
             << G4endl;
      chromo->Print();
    }
    delete chromo;
  }

  for(auto& test : tests)
  {
    delete test;
  }
  G4cout << "------------------------------------------------------------"
         << "--------------------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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
// $Id: G4DNAWaterDissociationDisplacer.cc 94289 2015-11-11 08:33:40Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"
#include "G4H2O2.hh"
#include "Randomize.hh"
#include "G4Molecule.hh"
#include "G4MolecularConfiguration.hh"
#include "G4RandomDirection.hh"

using namespace std;

const DisplacementType G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay =
    G4VMolecularDecayDisplacer::AddDisplacement();
const DisplacementType G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay =
    G4VMolecularDecayDisplacer::AddDisplacement();
const DisplacementType G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay =
    G4VMolecularDecayDisplacer::AddDisplacement();
const DisplacementType G4DNAWaterDissociationDisplacer::AutoIonisation =
    G4VMolecularDecayDisplacer::AddDisplacement();
const DisplacementType G4DNAWaterDissociationDisplacer::DissociativeAttachment =
    G4VMolecularDecayDisplacer::AddDisplacement();

//------------------------------------------------------------------------------

G4double G4DNAWaterDissociationDisplacer::ElectronProbaDistribution(G4double r)
{
  return (2*r+1)*exp(-2*r);
}

//------------------------------------------------------------------------------

G4DNAWaterDissociationDisplacer::G4DNAWaterDissociationDisplacer() :
    G4VMolecularDecayDisplacer(),
    fFastElectronDistrib(0, 10, 0.001)
{
  fProba1DFunction =
      std::bind((G4double(*)(G4double))
                &G4DNAWaterDissociationDisplacer::ElectronProbaDistribution,
                std::placeholders::_1);

  size_t nBins = 100;

  fElectronThermalization.reserve(nBins);

  double eps = 1./((int)nBins);

  // G4cout << "eps = " << eps << G4endl;

  double proba = eps;

  fElectronThermalization.push_back(0.);

  for(size_t i = 1 ; i < nBins ; ++i)
  {
    double r = fFastElectronDistrib.Revert(proba, fProba1DFunction);
    fElectronThermalization.push_back(r);
    proba+=eps;
  }
}

//------------------------------------------------------------------------------

G4DNAWaterDissociationDisplacer::~G4DNAWaterDissociationDisplacer()
{
  ;
}

//------------------------------------------------------------------------------

G4ThreeVector
G4DNAWaterDissociationDisplacer::
GetMotherMoleculeDisplacement(const G4MolecularDissociationChannel* theDecayChannel) const
{
  G4int decayType = theDecayChannel->GetDisplacementType();

  G4double RMSMotherMoleculeDisplacement = 0;

  if (decayType == Ionisation_DissociationDecay)
  {
    RMSMotherMoleculeDisplacement = 2.0 * nanometer;
  }
  else if (decayType == A1B1_DissociationDecay)
  {
    RMSMotherMoleculeDisplacement = 0. * nanometer;
  }
  else if (decayType == B1A1_DissociationDecay)
  {
    RMSMotherMoleculeDisplacement = 0. * nanometer;
  }
  else if (decayType == AutoIonisation)
  {
    RMSMotherMoleculeDisplacement = 2.0 * nanometer;
  }
  else if (decayType == DissociativeAttachment)
  {
    RMSMotherMoleculeDisplacement = 0. * nanometer;
  }

  if (RMSMotherMoleculeDisplacement == 0)
  {
    return G4ThreeVector(0, 0, 0);
  }
  auto RandDirection = radialDistributionOfProducts(
      RMSMotherMoleculeDisplacement);

  return RandDirection;
}

//------------------------------------------------------------------------------

vector<G4ThreeVector>
G4DNAWaterDissociationDisplacer::
GetProductsDisplacement(const G4MolecularDissociationChannel* theDecayChannel) const
{
  G4int nbProducts = theDecayChannel->GetNbProducts();
  vector<G4ThreeVector> theProductDisplacementVector(nbProducts);

  typedef map<const G4MoleculeDefinition*, G4double> RMSmap;
  RMSmap theRMSmap;

  G4int decayType = theDecayChannel->GetDisplacementType();

  if (decayType == Ionisation_DissociationDecay)
  {
    if (fVerbose) G4cout << "Ionisation_DissociationDecay" << G4endl;
    G4double RdmValue = G4UniformRand();

    if(RdmValue< 0.5)
    {
      // H3O
      theRMSmap[G4H3O::Definition()] = 0.* nanometer;
      // OH
      theRMSmap[G4OH::Definition()] = 0.8* nanometer;
    }
    else
    {
      // H3O
      theRMSmap[G4H3O::Definition()] = 0.8* nanometer;
      // OH
      theRMSmap[G4OH::Definition()] = 0.* nanometer;
    }

    for(int i = 0; i < nbProducts; i++)
    {
      G4double theRMSDisplacement;
      G4MolecularConfiguration* product = theDecayChannel->GetProduct(i);
      theRMSDisplacement = theRMSmap[product->GetDefinition()];

      if(theRMSDisplacement==0)
      {
        theProductDisplacementVector[i] = G4ThreeVector();
      }
      else
      {
        auto RandDirection = radialDistributionOfProducts(theRMSDisplacement);
        theProductDisplacementVector[i] = RandDirection;
      }
    }
  }
  else if(decayType == A1B1_DissociationDecay)
  {
    if(fVerbose)
      G4cout<<"A1B1_DissociationDecay"<<G4endl;
    G4double theRMSDisplacement = 2.4 * nanometer;
    auto RandDirection =
        radialDistributionOfProducts(theRMSDisplacement);

    for(G4int i =0; i < nbProducts; i++)
    {
      G4MolecularConfiguration* product = theDecayChannel->GetProduct(i);

      if(product->GetDefinition()== G4OH::Definition())
      {
        theProductDisplacementVector[i] = -1./18.*RandDirection;
      }
      else if(product->GetDefinition() == G4Hydrogen::Definition())
      {
        theProductDisplacementVector[i] = +17./18.*RandDirection;
      }
    }
  }
  else if(decayType == B1A1_DissociationDecay)
  {
    if(fVerbose)
      G4cout<<"B1A1_DissociationDecay"<<G4endl;
    G4double theRMSDisplacement = 0.8 * nanometer;
    auto RandDirection =
        radialDistributionOfProducts(theRMSDisplacement);

    G4int NbOfOH = 0;
    for(G4int i =0; i < nbProducts; i++)
    {
      G4MolecularConfiguration* product = theDecayChannel->GetProduct(i);
      if(product->GetDefinition() == G4H2::Definition())
      {
        theProductDisplacementVector[i] = -2./18.*RandDirection;
      }
      else if(product->GetDefinition() == G4OH::Definition())
      {
        G4ThreeVector OxygenDisplacement = +16./18.*RandDirection;
        G4double OHRMSDisplacement = 1.1 * nanometer;

        auto OHDisplacement =
            radialDistributionOfProducts(OHRMSDisplacement);

        if(NbOfOH==0)
        {
          OHDisplacement = 1./2.*OHDisplacement;
        }
        else
        {
          OHDisplacement = -1./2.*OHDisplacement;
        }

        theProductDisplacementVector[i] =
            OHDisplacement + OxygenDisplacement;

        NbOfOH ++;
      }
    }
  }
  else if(decayType == AutoIonisation)
  {
    if(fVerbose)
    G4cout<<"AutoIonisation"<<G4endl;
    G4double RdmValue = G4UniformRand();

    if(RdmValue< 0.5)
    {
      // H3O
      theRMSmap[G4H3O::Definition()] = 0.* nanometer;
      // OH
      theRMSmap[G4OH::Definition()] = 0.8* nanometer;
    }
    else
    {
      // H3O
      theRMSmap[G4H3O::Definition()] = 0.8* nanometer;
      // OH
      theRMSmap[G4OH::Definition()] = 0.* nanometer;
    }

    for(G4int i =0; i < nbProducts; i++)
    {
      G4MolecularConfiguration* product = theDecayChannel->GetProduct(i);
      auto theRMSDisplacement = theRMSmap[product->GetDefinition()];

      if(theRMSDisplacement==0)
      {
        theProductDisplacementVector[i] = G4ThreeVector();
      }
      else
      {
        auto RandDirection =
            radialDistributionOfProducts(theRMSDisplacement);
        theProductDisplacementVector[i] = RandDirection;
      }
      if(product->GetDefinition() == G4Electron_aq::Definition())
      {
        theProductDisplacementVector[i]=radialDistributionOfElectron();
      }
    }
  }
  else if(decayType == DissociativeAttachment)
  {
    if(fVerbose)
    G4cout<<"DissociativeAttachment"<<G4endl;
    G4double theRMSDisplacement = 0.8 * nanometer;
    auto RandDirection =
        radialDistributionOfProducts(theRMSDisplacement);

    G4int NbOfOH = 0;
    for(G4int i =0; i < nbProducts; i++)
    {
      G4MolecularConfiguration* product = theDecayChannel->GetProduct(i);
      if(product->GetDefinition() == G4H2::Definition())
      {
        theProductDisplacementVector[i] = -2./18.*RandDirection;
      }
      else if(product->GetDefinition() == G4OH::Definition())
      {
        G4ThreeVector OxygenDisplacement = +16./18.*RandDirection;
        G4double OHRMSDisplacement = 1.1 * nanometer;

        auto OHDisplacement =
            radialDistributionOfProducts(OHRMSDisplacement);

        if(NbOfOH==0)
        {
          OHDisplacement = 1./2.*OHDisplacement;
        }
        else
        {
          OHDisplacement = -1./2.*OHDisplacement;
        }

        theProductDisplacementVector[i] = OHDisplacement + OxygenDisplacement;

        NbOfOH ++;
      }
    }
  }

  return theProductDisplacementVector;
}

//------------------------------------------------------------------------------

G4ThreeVector
G4DNAWaterDissociationDisplacer::
radialDistributionOfProducts(G4double Rrms) const
{
  static const double inverse_sqrt_3 = 1./sqrt(3.);
  double sigma = Rrms*inverse_sqrt_3;
  double x = G4RandGauss::shoot(0.,sigma);
  double y = G4RandGauss::shoot(0.,sigma);
  double z = G4RandGauss::shoot(0.,sigma);
  return G4ThreeVector(x,y,z);
}

//------------------------------------------------------------------------------

G4ThreeVector
G4DNAWaterDissociationDisplacer::radialDistributionOfElectron() const
{
  G4double rand_value = G4UniformRand();
  size_t nBins = fElectronThermalization.size();
  size_t bin = size_t(rand_value*nBins);
  return fElectronThermalization[bin]*G4RandomDirection();
}

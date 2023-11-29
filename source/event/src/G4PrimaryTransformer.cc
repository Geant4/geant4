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
// G4PrimaryTransformer class implementation
//
// Author: Makoto Asai, 1999
// --------------------------------------------------------------------

#include "G4PrimaryTransformer.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4DecayProducts.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "Randomize.hh"

G4PrimaryTransformer::G4PrimaryTransformer()
{
  particleTable = G4ParticleTable::GetParticleTable();
  CheckUnknown();
}

void G4PrimaryTransformer::CheckUnknown()
{
  unknown = particleTable->FindParticle("unknown");
  unknownParticleDefined = unknown != nullptr;
  opticalphoton = particleTable->FindParticle("opticalphoton");
  opticalphotonDefined = opticalphoton != nullptr;
}
    
G4TrackVector*
G4PrimaryTransformer::GimmePrimaries(G4Event* anEvent, G4int trackIDCounter)
{
  trackID = trackIDCounter;

  for(auto tr : TV) delete tr;
  TV.clear();

  // Loop over vertices
  //
  G4PrimaryVertex* nextVertex = anEvent->GetPrimaryVertex();
  while(nextVertex != nullptr) // Loop checking 12.28.2015 M.Asai
  { 
    GenerateTracks(nextVertex);
    nextVertex = nextVertex->GetNext();
  }
  return &TV;
}

void G4PrimaryTransformer::GenerateTracks(G4PrimaryVertex* primaryVertex)
{
  G4double X0 = primaryVertex->GetX0();
  G4double Y0 = primaryVertex->GetY0();
  G4double Z0 = primaryVertex->GetZ0();
  G4double T0 = primaryVertex->GetT0();
  G4double WV = primaryVertex->GetWeight();

#ifdef G4VERBOSE
  if(verboseLevel>2)
  { 
    primaryVertex->Print();
  }
  else if (verboseLevel==1)
  {
    G4cout << "G4PrimaryTransformer::PrimaryVertex ("
           << X0 / mm << "(mm),"
           << Y0 / mm << "(mm),"
           << Z0 / mm << "(mm),"
           << T0 / nanosecond << "(nsec))" << G4endl;
  }
#endif

  G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
  while( primaryParticle != nullptr ) // Loop checking 12.28.2015 M.Asai
  {
    GenerateSingleTrack( primaryParticle, X0, Y0, Z0, T0, WV );
    primaryParticle = primaryParticle->GetNext();
  }
}

void G4PrimaryTransformer::
GenerateSingleTrack( G4PrimaryParticle* primaryParticle,
                     G4double x0, G4double y0, G4double z0,
                     G4double t0, G4double wv)
{
  G4ParticleDefinition* partDef = GetDefinition(primaryParticle);
  if(!IsGoodForTrack(partDef))
  {  // The particle cannot be converted to G4Track, check daughters
#ifdef G4VERBOSE
    if(verboseLevel>2)
    {
      G4cout << "Primary particle (PDGcode " << primaryParticle->GetPDGcode()
             << ") --- Ignored" << G4endl;
    }
#endif 
    G4PrimaryParticle* daughter = primaryParticle->GetDaughter();
    while(daughter != nullptr) // Loop checking 12.28.2015 M.Asai
    {
      GenerateSingleTrack(daughter,x0,y0,z0,t0,wv);
      daughter = daughter->GetNext();
    }
  }
  else  // The particle is defined in GEANT4
  {
    // Create G4DynamicParticle object
#ifdef G4VERBOSE
    if(verboseLevel>1)
    {
      G4cout << "Primary particle (" << partDef->GetParticleName()
             << ") --- Transferred with momentum "
             << primaryParticle->GetMomentum()
             << G4endl;
    }
#endif
    auto* DP = 
      new G4DynamicParticle(partDef,
                            primaryParticle->GetMomentumDirection(),
                            primaryParticle->GetKineticEnergy());
    if(opticalphotonDefined && partDef==opticalphoton
       && primaryParticle->GetPolarization().mag2()==0.)
    {
      if(nWarn<10)
      {
        G4Exception("G4PrimaryTransformer::GenerateSingleTrack",
                    "ZeroPolarization", JustWarning,
                    "Polarization of the optical photon is null.\
                     Random polarization is assumed.");
        G4cerr << "This warning message is issued up to 10 times." << G4endl;
        ++nWarn;
      }

      G4double angle = G4UniformRand() * 360.0*deg;
      G4ThreeVector normal (1., 0., 0.);
      G4ThreeVector kphoton = DP->GetMomentumDirection();
      G4ThreeVector product = normal.cross(kphoton);
      G4double modul2       = product*product;

      G4ThreeVector e_perpend (0., 0., 1.);
      if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
      G4ThreeVector e_paralle    = e_perpend.cross(kphoton);

      G4ThreeVector polar = std::cos(angle)*e_paralle
                          + std::sin(angle)*e_perpend;
      DP->SetPolarization(polar.x(),polar.y(),polar.z());
    }
    else
    {
      DP->SetPolarization(primaryParticle->GetPolX(),
                          primaryParticle->GetPolY(),
                          primaryParticle->GetPolZ());
    }
    if(primaryParticle->GetProperTime()>=0.0)
    {
      DP->SetPreAssignedDecayProperTime(primaryParticle->GetProperTime());
    }

    // Set Mass if it is specified
    //
    G4double pmas = primaryParticle->GetMass();
    if(pmas>=0.) { DP->SetMass(pmas); }

    // Set Charge if it is specified
    //
    if (primaryParticle->GetCharge()<DBL_MAX)
    {
      if (partDef->GetAtomicNumber() <0)
      {
        DP->SetCharge(primaryParticle->GetCharge());
      }
      else  // ions
      {
        G4int iz = partDef->GetAtomicNumber();
        G4int iq = static_cast<G4int>(primaryParticle->GetCharge()/eplus);
        G4int n_e = iz - iq;
        if (n_e>0) DP->AddElectron(0,n_e);  
      }
    } 
    // Set decay products to the DynamicParticle
    //
    SetDecayProducts( primaryParticle, DP );

    // Set primary particle
    //
    DP->SetPrimaryParticle(primaryParticle);

    // Set PDG code if it is different from G4ParticleDefinition
    //
    if(partDef->GetPDGEncoding()==0 && primaryParticle->GetPDGcode()!=0)
    {
      DP->SetPDGcode(primaryParticle->GetPDGcode());
    }

    // Check the particle is properly constructed
    //
    if(!CheckDynamicParticle(DP))
    {
      delete DP;
      return;
    }

    // Create G4Track object
    //
    G4Track* track = new G4Track(DP,t0,G4ThreeVector(x0,y0,z0));

    // Set trackID and let primary particle know it
    //
    ++trackID;
    track->SetTrackID(trackID);
    primaryParticle->SetTrackID(trackID);

    // Set parentID to 0 as a primary particle
    //
    track->SetParentID(0);

    // Set weight ( vertex weight * particle weight )
    //
    track->SetWeight(wv*(primaryParticle->GetWeight()));

    // Store it to G4TrackVector
    //
    TV.push_back( track );
  }
}

void G4PrimaryTransformer::
SetDecayProducts(G4PrimaryParticle* mother, G4DynamicParticle* motherDP)
{
  G4PrimaryParticle* daughter = mother->GetDaughter();
  if(daughter == nullptr) return;
  auto* decayProducts
    = (G4DecayProducts*)(motherDP->GetPreAssignedDecayProducts() );
  if(decayProducts == nullptr)
  {
    decayProducts = new G4DecayProducts(*motherDP);
    motherDP->SetPreAssignedDecayProducts(decayProducts);
  }
  while(daughter != nullptr)
  {
    G4ParticleDefinition* partDef = GetDefinition(daughter);
    if(!IsGoodForTrack(partDef))
    { 
#ifdef G4VERBOSE
      if(verboseLevel>2)
      {
        G4cout << " >> Decay product (PDGcode " << daughter->GetPDGcode()
               << ") --- Ignored" << G4endl;
      }
#endif 
      SetDecayProducts(daughter,motherDP);
    }
    else
    {
#ifdef G4VERBOSE
      if(verboseLevel>1)
      {
        G4cout << " >> Decay product (" << partDef->GetParticleName()
               << ") --- Attached with momentum " << daughter->GetMomentum()
               << G4endl;
      }
#endif
      auto* DP 
        = new G4DynamicParticle(partDef,daughter->GetMomentum());
      DP->SetPrimaryParticle(daughter);

      // Decay proper time for daughter
      //
      if(daughter->GetProperTime()>=0.0)
      {
        DP->SetPreAssignedDecayProperTime(daughter->GetProperTime());
      }

      // Set Charge and Mass is specified
      //
      if (daughter->GetCharge()<DBL_MAX)
      {
        DP->SetCharge(daughter->GetCharge());
      } 
      G4double pmas = daughter->GetMass();
      if(pmas>=0.)
      {
        DP->SetMass(pmas);
      }

      // Set Polarization
      //
      DP->SetPolarization(daughter->GetPolX(),
                          daughter->GetPolY(),
                          daughter->GetPolZ());
      decayProducts->PushProducts(DP);
      SetDecayProducts(daughter,DP);

      // Check the particle is properly constructed
      //
      if(!CheckDynamicParticle(DP))
      {
        delete DP;
        return;
      }
    }
    daughter = daughter->GetNext();
  }
}

void G4PrimaryTransformer::SetUnknnownParticleDefined(G4bool vl)
{
  unknownParticleDefined = vl;
  if(unknownParticleDefined && (unknown == nullptr))
  {
    G4cerr << "unknownParticleDefined cannot be set true because" << G4endl
           << "G4UnknownParticle is not defined in the physics list." << G4endl
           << "Command ignored." << G4endl;
    unknownParticleDefined = false;
  }
}

G4bool G4PrimaryTransformer::CheckDynamicParticle(G4DynamicParticle* DP)
{
  if(IsGoodForTrack(DP->GetDefinition())) return true;
  auto* decayProducts
    = (G4DecayProducts*)(DP->GetPreAssignedDecayProducts());
  if(decayProducts != nullptr && decayProducts->entries()>0) return true;
  G4cerr << G4endl
         << "G4PrimaryTransformer: a shortlived primary particle is found"
         << G4endl
         << " without any valid decay table nor pre-assigned decay mode."
         << G4endl;
  G4Exception("G4PrimaryTransformer", "InvalidPrimary", JustWarning,
              "This primary particle will be ignored.");
  return false;
}

G4ParticleDefinition*
G4PrimaryTransformer::GetDefinition(G4PrimaryParticle* pp)
{
  G4ParticleDefinition* partDef = pp->GetG4code();
  if(partDef == nullptr)
  {
    partDef = particleTable->FindParticle(pp->GetPDGcode());
  }
  if(unknownParticleDefined && ((partDef == nullptr)||partDef->IsShortLived()))
  {
    partDef = unknown;
  }
  return partDef;
}

G4bool G4PrimaryTransformer::IsGoodForTrack(G4ParticleDefinition* pd)
{
  if(pd == nullptr)
  { return false; }
  if(!(pd->IsShortLived()))
  { return true; }
  //
  // Following two lines should be removed if the user does not want to make
  // shortlived primary particle with proper decay table to be converted into
  // a track.
  //
  if(pd->GetDecayTable() != nullptr)
  { return true; }

  return false;
}

//******************************************************************************
// PrimaryGeneratorAction.cc
//
// 1.00 JMV, LLNL, Jan-2007:  First version.
//******************************************************************************
//
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "MultipleSource.hh"

//----------------------------------------------------------------------------//
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  // Specify center and radius of HEU ball
  G4ThreeVector* center = new G4ThreeVector(0.*m, 0.*m, 0.*m);
  G4double radius = 3.97*cm;

  // Specify isotopic composition and fission rates in fissions/sec
  G4DataVector* isotopeList = new G4DataVector(2);
  isotopeList->operator[](0) = 92238;
  isotopeList->operator[](1) = 92235;

  G4DataVector* intensityList = new G4DataVector(2);
  intensityList->operator[](0) = 2.368;
  intensityList->operator[](1) = .7475;
  

  time = 0; // set to 0 initially

  SponFissIsotope* currentIsotope;
  G4SPSPosDistribution* posDist;

  totalIntensity = 0.;
  nisotopes = isotopeList->size();
  for (G4int i=0; i < nisotopes; i++) {
    currentIsotope = new SponFissIsotope(static_cast<G4int> (isotopeList->operator[](i)));
    posDist = currentIsotope->GetPosDist();
    posDist->SetPosDisType("Volume");
    posDist->SetPosDisShape("Sphere");
    posDist->SetCentreCoords(*center);
    posDist->SetRadius(radius);
    totalIntensity += intensityList->operator[](i);
    if (i == 0) fissionSource = new MultipleSource(currentIsotope, intensityList->operator[](i));
    else fissionSource->AddaSource(currentIsotope, intensityList->operator[](i));
//    fissionSource->SetVerbosity(2);
  }
}

//----------------------------------------------------------------------------//
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  fissionSource->ClearAll();
  delete fissionSource;
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  // Compute at what time is the next fission
  G4double ran = G4UniformRand();
  time -= std::log(ran/totalIntensity)*second;
//  G4cout << "next spontaneous fission will be at " << time/second << G4endl;
  for (G4int i=0; i < nisotopes; i++) {
    fissionSource->SetCurrentSourceto(i);
    fissionSource->SetParticleTime(time);
  }
  
  // Pick an isotope among the ones listed above
  // and fission it.
  fissionSource->GeneratePrimaryVertex(anEvent);
}

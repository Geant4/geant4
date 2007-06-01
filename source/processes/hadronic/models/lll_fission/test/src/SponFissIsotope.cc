//******************************************************************************
// SponFissIsotope.cc
//
// 1.00 JMV, LLNL, Jan-2007:  First version.
//******************************************************************************
//
#include "SponFissIsotope.hh"

//----------------------------------------------------------------------------//
SponFissIsotope::SponFissIsotope()
{
}

//----------------------------------------------------------------------------//
SponFissIsotope::SponFissIsotope(G4int iso)
{
  neutron_definition = G4Neutron::Neutron();
  photon_definition = G4Gamma::Gamma();

  // verbosity
  verbosityLevel = 0;

  isotope = iso;
}

//----------------------------------------------------------------------------//
SponFissIsotope::~SponFissIsotope()
{
}

//----------------------------------------------------------------------------//
void SponFissIsotope::GeneratePrimaryVertex(G4Event* anEvent)
{ 
  // Generate a spontaneous fission using the fission library and emit
  // the neutrons and gamma-rays

  G4double time = GetParticleTime()/second;
  genspfissevt_(&isotope, &time);
  
  G4int nPrompt, gPrompt;
  nPrompt = getnnu_();
  gPrompt = getpnu_();

  if(verbosityLevel > 1) {
    G4cout << " nPrompt: " << nPrompt << G4endl
           << " gPrompt: " << gPrompt << G4endl;
  }

  // Position
  posDist = GetPosDist();
  G4ThreeVector particle_position = posDist->GenerateOne();

  // create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position, 0.);

  G4double mom, momx, momy, momz, eng;

  if(verbosityLevel >= 2)
    G4cout << "Creating primaries and assigning to vertex" << G4endl;

  G4DynamicParticle* it;
  // Build neutrons
  for(G4int i=0; i<nPrompt; i++)
  {
    it = new G4DynamicParticle();
    it->SetDefinition(neutron_definition);
    eng = getneng_(&i)*MeV;
    it->SetKineticEnergy(eng);
    mom = it->GetTotalMomentum();

    momx = mom*getndircosu_(&i);
    momy = mom*getndircosv_(&i);
    momz = mom*getndircosw_(&i);

    G4PrimaryParticle* particle = new G4PrimaryParticle(
                             neutron_definition,
                             momx, momy, momz,
                             eng);
    particle->SetMass(neutron_definition->GetPDGMass());
    particle->SetCharge(neutron_definition->GetPDGCharge());
    particle->SetPolarization(particle_polarization.x(),
                              particle_polarization.y(),
                              particle_polarization.z());

    if(verbosityLevel > 1){
      G4cout << "Particle name: "<<particle->GetG4code()->GetParticleName() << G4endl;
      G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
      G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
    }
    vertex->SetPrimary(particle);
  }
  // Build gammas
  for(G4int i=0; i<gPrompt; i++)
  {
    it = new G4DynamicParticle();
    it->SetDefinition(photon_definition);
    eng = getpeng_(&i)*MeV;
    it->SetKineticEnergy(eng);
    mom = it->GetTotalMomentum();

    momx = mom*getpdircosu_(&i);
    momy = mom*getpdircosv_(&i);
    momz = mom*getpdircosw_(&i);

    G4PrimaryParticle* particle = new G4PrimaryParticle(
                             photon_definition,
                             momx, momy, momz,
                             eng);
    particle->SetMass(photon_definition->GetPDGMass());
    particle->SetCharge(photon_definition->GetPDGCharge());
    particle->SetPolarization(particle_polarization.x(),
                              particle_polarization.y(),
                              particle_polarization.z());

    if(verbosityLevel > 1){
      G4cout << "Particle name: "<<particle->GetG4code()->GetParticleName() << G4endl;
      G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
      G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
    }

    vertex->SetPrimary(particle);
  }
  vertex->SetT0(time*second);
//  G4cout << "         Time: "<<vertex->GetT0()/second << G4endl;

  anEvent->AddPrimaryVertex( vertex );
  if(verbosityLevel > 1)
    G4cout << " Primary Vetex generated !"<< G4endl;
}

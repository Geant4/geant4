// $Id: MedicalBeam.cc,v 1.1 2006-02-27 09:52:54 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   MedicalBeam.cc
//
//                                         2005 Q
// ====================================================================
#include "MedicalBeam.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryVertex.hh"

using namespace CLHEP;

#include <cmath>

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////
MedicalBeam::MedicalBeam()
  : particle(0), 
    kineticE(1.*MeV),
    sourcePosition(G4ThreeVector()),
    SSD(1.*m), 
    fieldShape(MedicalBeam::SQUARE),
    fieldR(10.*cm)
//////////////////////////
{
  fieldXY[0]= fieldXY[1]= 10.*cm;
}


///////////////////////////
MedicalBeam::~MedicalBeam()
///////////////////////////
{
}


////////////////////////////////////////////////////////
G4ThreeVector MedicalBeam::GenerateBeamDirection() const
////////////////////////////////////////////////////////
{
  // uniform distribution in a limitted solid angle
  G4double dr;
  if(fieldShape == MedicalBeam::SQUARE) {
    dr= std::sqrt(sqr(fieldXY[0]/2.)+sqr(fieldXY[1]/2.));
  } else {
    dr= fieldR;
  }

  G4double sin0= dr/SSD;
  G4double cos0= std::sqrt(1.-sqrt(sin0));

  G4double dcos, dsin, dphi, z;

  G4double x= DBL_MAX;
  G4double y= DBL_MAX;

  G4double xmax, ymax;
  if(fieldShape == MedicalBeam::SQUARE) {
    xmax= fieldXY[0]/2./SSD;
    ymax= fieldXY[1]/2./SSD;
  } else {
    xmax= ymax= DBL_MAX-1.;
  }

  while(! (std::abs(x)< xmax && std::abs(y)< ymax) ) {
    dcos= RandFlat::shoot(cos0, 1.);
    dsin= std::sqrt(1.-sqr(dcos));
    dphi= RandFlat::shoot(0., twopi);

    x= std::cos(dphi)*dsin*dcos;
    y= std::sin(dphi)*dsin*dcos;
  }
  z= dcos;

  return G4ThreeVector(x,y,z);
}


/////////////////////////////////////////////////////
void MedicalBeam::GeneratePrimaries(G4Event* anEvent)
/////////////////////////////////////////////////////
{
  if(particle==0) return;

  // create a new vertex
  G4PrimaryVertex* vertex= new G4PrimaryVertex(sourcePosition, 0.*ns);

  // momentum
  G4double mass= particle-> GetPDGMass();
  G4double p= std::sqrt(sqr(mass+kineticE)-sqr(mass));
  G4ThreeVector pmon= p*GenerateBeamDirection();
  G4PrimaryParticle* primary= new G4PrimaryParticle(particle, 
						    pmon.x(), 
						    pmon.y(), 
						    pmon.z());
  // set primary to vertex
  vertex-> SetPrimary(primary);

  // set vertex to event
  anEvent-> AddPrimaryVertex(vertex);
}


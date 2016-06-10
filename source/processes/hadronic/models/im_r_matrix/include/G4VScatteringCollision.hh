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
// $Id: G4VScatteringCollision.hh,v 1.3 2006-06-29 20:36:14 gunter Exp $ //
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4VElasticCollision
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4VSCATTERINGCOLLISION_HH
#define G4VSCATTERINGCOLLISION_HH

#include "globals.hh"
#include "G4Log.hh"
#include "G4VCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"

class G4KineticTrack;


class G4VScatteringCollision : public G4VCollision
{

public:

  G4VScatteringCollision();
  void establish_G4MT_TLS_G4VScatteringCollision();
  virtual ~G4VScatteringCollision();

  G4bool operator==(const G4VScatteringCollision &right) const;
  G4bool operator!=(const G4VScatteringCollision &right) const;

  virtual G4KineticTrackVector* FinalState(const G4KineticTrack& trk1, 
					      const G4KineticTrack& trk2) const;
  virtual const G4VAngularDistribution* GetAngularDistribution() const
  {
    return theAngularDistribution;
  }

private:
  G4VScatteringCollision(const G4VScatteringCollision &);
  G4VScatteringCollision & operator= (const G4VScatteringCollision &);

protected:

  virtual const std::vector<const G4ParticleDefinition*> & GetOutgoingParticles() const = 0;

private:  

  double BrWigInt0(const double x, const double gamma, const double m0) const
    { return 2.0*gamma*std::atan( 2.0 * (x-m0)/ gamma  ); }

  G4double BrWigInt1(const G4double x, const G4double gamma, const G4double m0) const
    { return 0.5*gamma*gamma*G4Log( (x-m0)*(x-m0)+gamma*gamma/4.0 ) + m0*BrWigInt0(x,gamma,m0); }

  double BrWigInv(const double x, const double gamma, const double m0) const
    { return 0.5*gamma*std::tan( 0.5*x/gamma )+m0; }
  
  double SampleResonanceMass(const double poleMass, 
			     const double width,
			     const double minMass,
			     const double maxMass) const;
private:

 G4VAngularDistribution * theAngularDistribution;

};    




#endif

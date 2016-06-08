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
//
// $Id: G4VAnnihilationCollision.hh,v 1.2 2002/12/12 19:17:42 gunter Exp $ //
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

#ifndef G4VANNIHILATIONCOLLISION_HH
#define G4VANNIHILATIONCOLLISION_HH

#include "globals.hh"
#include "G4VCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"

class G4KineticTrack;


class G4VAnnihilationCollision : public G4VCollision
{

public:

  G4VAnnihilationCollision();

  virtual ~G4VAnnihilationCollision();

  G4bool operator==(const G4VAnnihilationCollision &right) const;
  G4bool operator!=(const G4VAnnihilationCollision &right) const;

  virtual G4KineticTrackVector* FinalState(const G4KineticTrack& trk1, 
					      const G4KineticTrack& trk2) const;
  virtual const G4VAngularDistribution* GetAngularDistribution() const
  {
    return theAngularDistribution;
  }

protected:

  virtual const G4ParticleDefinition* GetOutgoingParticle(const G4KineticTrack& trk1, 
							  const G4KineticTrack& trk2) const = 0;

private:  

  double BrWigInt0(const double x, const double gamma, const double m0) const
    { return 2.0*gamma*atan( 2.0 * (x-m0)/ gamma  ); }

  G4double BrWigInt1(const G4double x, const G4double gamma, const G4double m0) const
    { return 0.5*gamma*gamma*log( (x-m0)*(x-m0)+gamma*gamma/4.0 ) + m0*BrWigInt0(x,gamma,m0); }

  double BrWigInv(const double x, const double gamma, const double m0) const
    { return 0.5*gamma*tan( 0.5*x/gamma )+m0; }
  
  double SampleResonanceMass(const double poleMass, 
			     const double width,
			     const double minMass,
			     const double maxMass) const;
private:

 G4VAngularDistribution * theAngularDistribution;

};    




#endif

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
// $Id: G4QuasiElasticChannel.hh,v 1.2 2008-04-08 15:41:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Author : Gunter Folger March 2007
// Class Description
// Final state production model for theoretical models of hadron inelastic
// quasi elastic scattering in geant4;
// Class Description - End
//
// Modified:

#ifndef G4QuasiElasticChannel_h
#define G4QuasiElasticChannel_h

#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4KineticTrackVector.hh"

#include "G4QuasiFreeRatios.hh"

#include "G4Fancy3DNucleus.hh"

class G4QuasiElasticChannel
{
  public:
	G4QuasiElasticChannel();
	~G4QuasiElasticChannel();
	
	G4double GetFraction(G4Nucleus &theNucleus,
				const G4DynamicParticle & thePrimary);

	G4KineticTrackVector * Scatter(G4Nucleus &theNucleus,
					const G4DynamicParticle & thePrimary);
					
  private:
        G4QuasiElasticChannel(const G4QuasiElasticChannel &);
	const G4QuasiElasticChannel & operator=(const G4QuasiElasticChannel &);
	int operator==(const G4QuasiElasticChannel &) const;
	int operator!=(const G4QuasiElasticChannel &) const;

   private:
   	G4QuasiFreeRatios * theQuasiElastic;
	G4Fancy3DNucleus  the3DNucleus;
};

#endif

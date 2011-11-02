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
/*
 * =============================================================================
 *
 *       Filename:  CexmcFakeCrossSectionData.hh
 *
 *    Description:  fake cross section data (to be never really applied)
 *
 *        Version:  1.0
 *        Created:  18.09.2011 21:47:10
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_FAKE_CROSS_SECTION_DATA_HH
#define CEXMC_FAKE_CROSS_SECTION_DATA_HH

#include <G4VCrossSectionDataSet.hh>

class  G4DynamicParticle;
class  G4Element;
class  G4ParticleDefinition;


class  CexmcFakeCrossSectionData : public G4VCrossSectionDataSet
{
    public:
        G4bool    IsApplicable( const G4DynamicParticle *, const G4Element * )
        {
            return false;
        }

        G4double  GetCrossSection( const G4DynamicParticle *,
                                   const G4Element *, G4double )
        {
            return 0;
        }

        void      BuildPhysicsTable( const G4ParticleDefinition & )
        {}

        void      DumpPhysicsTable( const G4ParticleDefinition & )
        {}
};


#endif


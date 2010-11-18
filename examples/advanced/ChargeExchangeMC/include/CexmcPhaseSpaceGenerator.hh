/*
 * =============================================================================
 *
 *       Filename:  CexmcPhaseSpaceGenerator.hh
 *
 *    Description:  phase space generator interface
 *
 *        Version:  1.0
 *        Created:  08.09.2010 13:42:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PHASE_SPACE_GENERATOR_HH
#define CEXMC_PHASE_SPACE_GENERATOR_HH

#include <vector>
#include <G4Types.hh>
#include <G4LorentzVector.hh>


struct  CexmcPhaseSpaceOutVectorElement
{
    CexmcPhaseSpaceOutVectorElement( G4LorentzVector *  lVec, G4double  mass ) :
        lVec( lVec ), mass( mass )
    {}

    G4LorentzVector *  lVec;

    G4double           mass;
};


typedef std::vector< const G4LorentzVector * >  CexmcPhaseSpaceInVector;

typedef std::vector< CexmcPhaseSpaceOutVectorElement >
                                                CexmcPhaseSpaceOutVector;


class  CexmcPhaseSpaceGenerator
{
    public:
        CexmcPhaseSpaceGenerator();

        virtual ~CexmcPhaseSpaceGenerator();

    public:
        virtual G4bool    CheckKinematics( void );

        virtual G4double  Generate( void ) = 0;

    public:
        void  SetParticles( const CexmcPhaseSpaceInVector &  inVec_,
                            const CexmcPhaseSpaceOutVector &  outVec_ );

        void  SetFermiEnergyDependence( G4bool  on = true );

    protected:
        virtual void    ParticleChangeHook( void );

        virtual void    FermiEnergyDepStatusChangeHook( void );

    protected:
        CexmcPhaseSpaceInVector   inVec;

        CexmcPhaseSpaceOutVector  outVec;

        G4bool                    fermiEnergyDepIsOn;

        G4double                  totalEnergy;

        G4double                  totalMass;
};


#endif


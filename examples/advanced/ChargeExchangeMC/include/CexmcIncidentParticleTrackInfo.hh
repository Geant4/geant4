/*
 * =============================================================================
 *
 *       Filename:  CexmcIncidentParticleTrackInfo.hh
 *
 *    Description:  incident particle track info
 *
 *        Version:  1.0
 *        Created:  18.05.2010 13:04:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_INCIDENT_PARTICLE_TRACK_INFO_HH
#define CEXMC_INCIDENT_PARTICLE_TRACK_INFO_HH

#include "CexmcTrackInfo.hh"


class  CexmcIncidentParticleTrackInfo : public CexmcTrackInfo
{
    public:
        explicit CexmcIncidentParticleTrackInfo( CexmcTrackType  trackType =
                                                            CexmcInsipidTrack );

    public:
        G4int     GetTypeInfo( void ) const;

    public:
        G4double  GetCurrentTrackLengthInTarget( void ) const;

        void      AddTrackLengthInTarget( G4double  value );

        void      SetNeedsTrackLengthResampling( G4bool  on = true );

        G4double  GetFinalTrackLengthInTarget( void ) const;

        void      SetFinalTrackLengthInTarget( G4double  value );

        void      ResetCurrentTrackLengthInTarget( void );

        G4bool    NeedsTrackLengthResampling( void ) const;

        G4bool    IsStudiedProcessActivated( void ) const;

        void      ActivateStudiedProcess( G4bool  on = true );

    private:
        G4double  currentTrackLengthInTarget;

        G4double  finalTrackLengthInTarget;

        G4bool    isStudiedProcessActivated;

        G4bool    needsTrackLengthResampling;
};


inline G4double  CexmcIncidentParticleTrackInfo::GetCurrentTrackLengthInTarget(
                                                                    void ) const
{
    return currentTrackLengthInTarget;
}


inline void  CexmcIncidentParticleTrackInfo::AddTrackLengthInTarget(
                                                            G4double  value )
{
    currentTrackLengthInTarget += value;
}


inline void  CexmcIncidentParticleTrackInfo::SetNeedsTrackLengthResampling(
                                                                    G4bool  on )
{
    needsTrackLengthResampling = on;
}


inline G4double  CexmcIncidentParticleTrackInfo::GetFinalTrackLengthInTarget(
                                                                    void ) const
{
    return finalTrackLengthInTarget;
}


inline void  CexmcIncidentParticleTrackInfo::SetFinalTrackLengthInTarget(
                                                            G4double  value )
{
    finalTrackLengthInTarget = value;
}


inline void  CexmcIncidentParticleTrackInfo::ResetCurrentTrackLengthInTarget(
                                                                        void )
{
    currentTrackLengthInTarget = 0.;
}


inline G4bool  CexmcIncidentParticleTrackInfo::NeedsTrackLengthResampling(
                                                                    void ) const
{
    return needsTrackLengthResampling;
}


inline G4bool  CexmcIncidentParticleTrackInfo::IsStudiedProcessActivated(
                                                                    void ) const
{
    return isStudiedProcessActivated;
}


inline void  CexmcIncidentParticleTrackInfo::ActivateStudiedProcess(
                                                                G4bool  on )
{
    isStudiedProcessActivated = on;
}


#endif


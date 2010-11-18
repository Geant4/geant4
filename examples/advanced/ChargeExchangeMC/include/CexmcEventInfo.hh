/*
 * =============================================================================
 *
 *       Filename:  CexmcEventInfo.hh
 *
 *    Description:  event information passed to run manager
 *
 *        Version:  1.0
 *        Created:  04.12.2009 14:47:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_EVENT_INFO_HH
#define CEXMC_EVENT_INFO_HH

#include <G4Types.hh>
#include "G4VUserEventInformation.hh"


class  CexmcEventInfo : public G4VUserEventInformation
{
    public:
        CexmcEventInfo( G4bool  edTriggerIsOk, G4bool  tpTriggerIsOk,
                        G4bool  reconstructionIsOk );

    public:
        void  Print( void ) const;

    public:
        G4bool  EdTriggerIsOk( void ) const;

        G4bool  TpTriggerIsOk( void ) const;

        G4bool  ReconstructionIsOk( void ) const;

    private:
        G4bool  edTriggerIsOk;

        G4bool  tpTriggerIsOk;

        G4bool  reconstructionIsOk;
};


inline G4bool  CexmcEventInfo::EdTriggerIsOk( void ) const
{
    return edTriggerIsOk;
}


inline G4bool  CexmcEventInfo::TpTriggerIsOk( void ) const
{
    return tpTriggerIsOk;
}


inline G4bool  CexmcEventInfo::ReconstructionIsOk( void ) const
{
    return reconstructionIsOk;
}


#endif


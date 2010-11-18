/*
 * =============================================================================
 *
 *       Filename:  CexmcRunManager.hh
 *
 *    Description:  run manager
 *
 *        Version:  1.0
 *        Created:  03.11.2009 20:17:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_RUN_MANAGER_HH
#define CEXMC_RUN_MANAGER_HH

#include <set>
#include <limits>
#include <boost/archive/binary_oarchive.hpp>
#include <G4RunManager.hh>
#include "CexmcRunSObject.hh"
#include "CexmcException.hh"
#include "CexmcCommon.hh"

class  CexmcRunManagerMessenger;
class  CexmcPhysicsManager;
class  CexmcEventFastSObject;
#ifdef CEXMC_USE_CUSTOM_FILTER
class  CexmcCustomFilterEval;
#endif


typedef std::set< CexmcOutputDataType >  CexmcOutputDataTypeSet;


class  CexmcRunManager : public G4RunManager
{
    public:
        explicit CexmcRunManager( const G4String &  projectId = "",
                                  const G4String &  rProject = "",
                                  G4bool  overrideExistingProject = false );

        virtual ~CexmcRunManager();

    public:
        void  SetPhysicsManager( CexmcPhysicsManager *  physicsManager_ );

        void  SetProductionModelType(
                            CexmcProductionModelType  productionModelType_ );

        void  SetGdmlFileName( const G4String &  gdmlFileName_ );

        void  SetGdmlFileValidation( G4bool  on = true );

        void  SetGuiMacroName( const G4String &  guiMacroName_ );

#ifdef CEXMC_USE_CUSTOM_FILTER
        void  SetCustomFilter( const G4String &  cfFileName_ );
#endif

        void  SetEventCountPolicy( CexmcEventCountPolicy  value );

        void  SetEventDataVerboseLevel( CexmcEventDataVerboseLevel  value );

        void  ReadProject( void );

        void  SaveProject( void );

        void  PrintReadRunData( void ) const;

        void  ReadAndPrintEventsData( void ) const;

        void  PrintReadData( const CexmcOutputDataTypeSet &  outputData ) const;

        void  ReplayEvents( G4int  nEvents = 0 );

        void  SeekTo( G4int  eventNmb = 1 );

        void  EnableLiveHistograms( G4bool  on = true );

        void  SkipInteractionsWithoutEDTonWrite( G4bool  on = true );

    public:
        CexmcPhysicsManager *     GetPhysicsManager( void );

        CexmcProductionModelType  GetProductionModelType( void ) const;

        G4String                  GetGdmlFileName( void ) const;

        G4bool                    ShouldGdmlFileBeValidated( void ) const;

        G4String                  GetGuiMacroName( void ) const;

        G4bool                    ProjectIsSaved( void ) const;

        G4bool                    ProjectIsRead( void ) const;

        G4String                  GetProjectsDir( void ) const;

        G4String                  GetProjectId( void ) const;

        boost::archive::binary_oarchive *  GetEventsArchive( void ) const;

        boost::archive::binary_oarchive *  GetFastEventsArchive( void ) const;

        G4bool                    AreLiveHistogramsEnabled( void ) const;

        CexmcEventDataVerboseLevel  GetEventDataVerboseLevel( void ) const;

        void                      BeamParticleChangeHook( void );

    protected:
        void  DoEventLoop( G4int  nEvent, const char *  macroFile,
                           G4int  nSelect );

    private:
        void  DoCommonEventLoop( G4int  nEvent, const G4String &  cmd,
                                 G4int  nSelect );

        void  DoReadEventLoop( G4int  nEvent );

        void  SaveCurrentTPTEvent( const CexmcEventFastSObject &  evFastSObject,
                                   const CexmcAngularRangeList &  angularRanges,
                                   G4bool  writeToDatabase );

    private:
        void  ReadPreinitProjectData( void );

    private:
        CexmcBasePhysicsUsed        basePhysicsUsed;

        CexmcProductionModelType    productionModelType;

        G4String                    gdmlFileName;

        G4bool                      shouldGdmlFileBeValidated;

        G4bool                      zipGdmlFile;

        G4String                    projectsDir;

        G4String                    projectId;

        G4String                    rProject;

        G4String                    guiMacroName;

        G4String                    cfFileName;

#ifdef CEXMC_USE_CUSTOM_FILTER
        CexmcCustomFilterEval *     customFilter;
#endif

        CexmcEventCountPolicy       eventCountPolicy;

        G4bool                      areLiveHistogramsEnabled;

        G4bool                      skipInteractionsWithoutEDTonWrite;

        CexmcEventDataVerboseLevel  evDataVerboseLevel;

        CexmcEventDataVerboseLevel  rEvDataVerboseLevel;

        CexmcRunSObject             sObject;

    private:
        G4int                       numberOfEventsProcessed;

        G4int                       numberOfEventsProcessedEffective;

        G4int                       curEventRead;

    private:
        boost::archive::binary_oarchive *  eventsArchive;

        boost::archive::binary_oarchive *  fastEventsArchive;

    private:
        CexmcPhysicsManager *       physicsManager;

    private:
        CexmcRunManagerMessenger *  messenger;
};


inline void  CexmcRunManager::SetPhysicsManager(
                                        CexmcPhysicsManager *  physicsManager_ )
{
    physicsManager = physicsManager_;
}


inline void  CexmcRunManager::SetProductionModelType(
                                CexmcProductionModelType  productionModelType_ )
{
    if ( ProjectIsRead() )
        throw CexmcException( CexmcCmdIsNotAllowed );

    productionModelType = productionModelType_;
}


inline void  CexmcRunManager::SetGdmlFileName( const G4String &  gdmlFileName_ )
{
    if ( ProjectIsRead() )
        throw CexmcException( CexmcCmdIsNotAllowed );

    gdmlFileName = gdmlFileName_;
}


inline void  CexmcRunManager::SetGdmlFileValidation( G4bool  on )
{
    shouldGdmlFileBeValidated = on;
}


inline void  CexmcRunManager::SetGuiMacroName( const G4String &  guiMacroName_ )
{
    guiMacroName = guiMacroName_;
}


inline void  CexmcRunManager::SetEventCountPolicy(
                                            CexmcEventCountPolicy  value )
{
    if ( ProjectIsRead() )
        throw CexmcException( CexmcCmdIsNotAllowed );

    eventCountPolicy = value;
}


inline void  CexmcRunManager::SetEventDataVerboseLevel(
                                            CexmcEventDataVerboseLevel  value )
{
    if ( ProjectIsRead() && value > rEvDataVerboseLevel )
        throw CexmcException( CexmcPoorEventData );

    evDataVerboseLevel = value;
}


inline CexmcPhysicsManager *  CexmcRunManager::GetPhysicsManager( void )
{
    return physicsManager;
}


inline CexmcProductionModelType
                        CexmcRunManager::GetProductionModelType( void ) const
{
    return productionModelType;
}


inline G4String  CexmcRunManager::GetGdmlFileName( void ) const
{
    return gdmlFileName;
}


inline G4bool  CexmcRunManager::ShouldGdmlFileBeValidated( void ) const
{
    return shouldGdmlFileBeValidated;
}


inline G4String  CexmcRunManager::GetGuiMacroName( void ) const
{
    return guiMacroName;
}


inline G4bool  CexmcRunManager::ProjectIsSaved( void ) const
{
    return projectId != "";
}


inline G4bool  CexmcRunManager::ProjectIsRead( void ) const
{
    return rProject != "";
}


inline G4String  CexmcRunManager::GetProjectsDir( void ) const
{
    return projectsDir;
}


inline G4String  CexmcRunManager::GetProjectId( void ) const
{
    return projectId;
}


inline boost::archive::binary_oarchive *  CexmcRunManager::GetEventsArchive(
                                                                    void ) const
{
    return eventsArchive;
}


inline boost::archive::binary_oarchive *  CexmcRunManager::GetFastEventsArchive(
                                                                    void ) const
{
    return fastEventsArchive;
}


inline void  CexmcRunManager::ReplayEvents( G4int  nEvents )
{
    if ( ! ProjectIsRead() )
        return;

    if ( nEvents == 0 )
        nEvents = std::numeric_limits< G4int >::max();

    BeamOn( nEvents );
}


inline void  CexmcRunManager::SeekTo( G4int  eventNmb )
{
    if ( ! ProjectIsRead() )
        return;

    curEventRead = eventNmb;
}


inline void  CexmcRunManager::EnableLiveHistograms( G4bool  on )
{
    areLiveHistogramsEnabled = on;
}


inline G4bool  CexmcRunManager::AreLiveHistogramsEnabled( void ) const
{
    return areLiveHistogramsEnabled;
}


inline CexmcEventDataVerboseLevel  CexmcRunManager::GetEventDataVerboseLevel(
                                                                    void ) const
{
    return evDataVerboseLevel;
}


inline void  CexmcRunManager::SkipInteractionsWithoutEDTonWrite( G4bool  on )
{
    skipInteractionsWithoutEDTonWrite = on;
}


#endif


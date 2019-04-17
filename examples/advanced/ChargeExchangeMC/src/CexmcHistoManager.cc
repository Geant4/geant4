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
 * ============================================================================
 *
 *       Filename:  CexmcHistoManager.cc
 *
 *    Description:  histograming manager (singleton)
 *
 *        Version:  1.0
 *        Created:  26.11.2009 21:00:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifdef CEXMC_USE_ROOT

#include <iostream>
#include <iomanip>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TObject.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TString.h>
#ifdef CEXMC_USE_ROOTQT
#include <TCanvas.h>
#include <TList.h>
#include <QApplication>
#include <QFont>
#include <QMenu>
#include <G4UIQt.hh>
#include "CexmcMessenger.hh"
#endif
#include <G4LogicalVolume.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4SystemOfUnits.hh>
#include "CexmcHistoManager.hh"
#include "CexmcHistoManagerMessenger.hh"
#include "CexmcProductionModel.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcRunManager.hh"
#include "CexmcSetup.hh"
#include "CexmcException.hh"
#include "CexmcHistoWidget.hh"


namespace
{
    const G4double  CexmcHistoBeamMomentumMin( 0.0 * GeV );
    const G4double  CexmcHistoBeamMomentumMax( 1.0 * GeV );
    const G4double  CexmcHistoBeamMomentumResolution( 0.5 * MeV );
    const G4double  CexmcHistoTPResolution( 0.1 * cm );
    const G4double  CexmcHistoTPSafetyArea( 1.0 * cm );
    const G4double  CexmcHistoMassResolution( 1.0 * MeV );
    const G4double  CexmcHistoEnergyMax( 1.0 * GeV );
    const G4double  CexmcHistoEnergyResolution( 1.0 * MeV );
    const G4double  CexmcHistoMissEnergyMin( -0.1 * GeV );
    const G4double  CexmcHistoMissEnergyMax( 0.2 * GeV );
    const G4double  CexmcHistoMissEnergyResolution( 0.2 * MeV );
    const G4double  CexmcHistoAngularResolution( 0.5 );
    const G4double  CexmcHistoAngularCResolution( 0.001 );
#ifdef CEXMC_USE_ROOTQT
    const G4int     CexmcHistoCanvasWidth( 800 );
    const G4int     CexmcHistoCanvasHeight( 600 );
#endif
    const G4String  CexmcHistoDirectoryHandle( "histograms" );
    const G4String  CexmcHistoDirectoryTitle( "Histograms" );
}


CexmcHistoManager *  CexmcHistoManager::instance( NULL );


CexmcHistoManager *  CexmcHistoManager::Instance( void )
{
    if ( instance == NULL )
        instance = new CexmcHistoManager;

    return instance;
}


void  CexmcHistoManager::Destroy( void )
{
    delete instance;
    instance = NULL;
}


CexmcHistoManager::CexmcHistoManager() : outFile( NULL ),
    isInitialized( false ), opName( "" ), nopName( "" ), opMass( 0. ),
    nopMass( 0. ), verboseLevel( 0 ),
#ifdef CEXMC_USE_ROOTQT
    rootCanvas( NULL ), areLiveHistogramsEnabled( false ),
    isHistoMenuInitialized( false ), drawOptions1D( "" ), drawOptions2D( "" ),
    drawOptions3D( "" ), histoMenuHandle( "" ), histoMenuLabel( "" ),
#endif
    messenger( NULL )
{
    for ( int  i( 0 ); i < CexmcHistoType_SIZE; ++i )
    {
        histos.insert( CexmcHistoPair( CexmcHistoType( i ),
                                       CexmcHistoVector() ) );
    }

    messenger = new CexmcHistoManagerMessenger( this );
}


CexmcHistoManager::~CexmcHistoManager()
{
    if ( outFile )
    {
        outFile->Write();
        outFile->Close();
    }

    /* all histograms will be deleted by outFile destructor! */
    delete outFile;
#ifdef CEXMC_USE_ROOTQT
    delete rootCanvas;
#endif
    delete messenger;
}


void  CexmcHistoManager::AddHisto( const CexmcHistoData &  data,
                                   const CexmcAngularRange &  aRange )
{
    G4String  fullName( data.name );
    G4String  fullTitle( data.title );
    G4String  rangeTypeLabel;
    G4String  triggerTypeLabel;
    G4String  decorTriggerTypeLabel;

    if ( data.isARHisto )
    {
        if ( data.isARRec )
        {
            fullName += "_arrec";
            rangeTypeLabel = "rec";
        }
        else
        {
            fullName += "_arreal";
            rangeTypeLabel = "real";
        }
    }

    switch ( data.triggerType )
    {
    case CexmcTPT :
        fullName += "_tpt";
        decorTriggerTypeLabel = " --tpt--";
        fullTitle += decorTriggerTypeLabel;
        triggerTypeLabel = "tpt";
        break;
    case CexmcEDT :
        fullName += "_edt";
        decorTriggerTypeLabel = " --edt--";
        fullTitle += decorTriggerTypeLabel;
        triggerTypeLabel = "edt";
        break;
    case CexmcRT :
        fullName += "_rt";
        decorTriggerTypeLabel = " --rt--";
        fullTitle += decorTriggerTypeLabel;
        triggerTypeLabel = "rt";
        break;
    default :
        break;
    }

    CexmcHistosMap::iterator  found( histos.find( data.type ) );

    if ( found == histos.end() )
        throw CexmcException( CexmcWeirdException );

    CexmcHistoVector &  histoVector( found->second );

    if ( data.isARHisto )
    {
        G4bool  dirOk( false );

        dirOk = gDirectory->Get( fullName ) != NULL;

        if ( ! dirOk )
            dirOk = ( gDirectory->mkdir( fullName, fullTitle ) != NULL );

        if ( dirOk )
            gDirectory->cd( fullName );

        std::ostringstream  histoName;
        std::ostringstream  histoTitle;
        histoName << data.name << "_r" << aRange.index + 1 << rangeTypeLabel <<
                "_" << triggerTypeLabel;
        histoTitle << data.title << " {range " << aRange.index + 1 <<
                rangeTypeLabel << " [" << std::fixed <<
                std::setprecision( 4 ) << aRange.top << ", " <<
                aRange.bottom << ")}" << decorTriggerTypeLabel;
        CreateHisto( histoVector, data.impl, histoName.str(), histoTitle.str(),
                     data.axes );

        if ( dirOk )
            gDirectory->cd( ".." );
    }
    else
    {
        CreateHisto( histoVector, data.impl, fullName, fullTitle, data.axes );
    }
}


void  CexmcHistoManager::CreateHisto( CexmcHistoVector &  histoVector,
                        CexmcHistoImpl  histoImpl, const G4String &  name,
                        const G4String &  title, const CexmcHistoAxes &  axes )
{
    TH1 *  histo( NULL );

    switch ( histoImpl )
    {
    case Cexmc_TH1F :
        histo = new TH1F( name, title, axes.at( 0 ).nBins,
                          axes.at( 0 ).nBinsMin, axes.at( 0 ).nBinsMax );
        break;
    case Cexmc_TH2F :
        histo = new TH2F( name, title, axes.at( 0 ).nBins,
                          axes.at( 0 ).nBinsMin, axes.at( 0 ).nBinsMax,
                          axes.at( 1 ).nBins, axes.at( 1 ).nBinsMin,
                          axes.at( 1 ).nBinsMax );
        break;
    case Cexmc_TH3F :
        histo = new TH3F( name, title, axes.at( 0 ).nBins,
                          axes.at( 0 ).nBinsMin, axes.at( 0 ).nBinsMax,
                          axes.at( 1 ).nBins, axes.at( 1 ).nBinsMin,
                          axes.at( 1 ).nBinsMax, axes.at( 2 ).nBins,
                          axes.at( 2 ).nBinsMin, axes.at( 2 ).nBinsMax );
        break;
    default :
        break;
    }

    if ( histo )
        histoVector.push_back( histo );
}


void  CexmcHistoManager::Initialize( void )
{
    if ( isInitialized )
        return;

    CexmcRunManager *       runManager( static_cast< CexmcRunManager * >(
                                            G4RunManager::GetRunManager() ) );
    CexmcPhysicsManager *   physicsManager( runManager->GetPhysicsManager() );

    if ( ! physicsManager )
        throw CexmcException( CexmcWeirdException );

    CexmcProductionModel *  productionModel( physicsManager->
                                                        GetProductionModel() );

    if ( ! productionModel )
        throw CexmcException( CexmcWeirdException );

    G4ParticleDefinition *  outputParticle(
                                productionModel->GetOutputParticle() );
    G4ParticleDefinition *  nucleusOutputParticle(
                                productionModel->GetNucleusOutputParticle() );

    if ( ! outputParticle || ! nucleusOutputParticle )
        throw CexmcException( CexmcIncompleteProductionModel );

    opName = outputParticle->GetParticleName();
    nopName = nucleusOutputParticle->GetParticleName();
    opMass = outputParticle->GetPDGMass();
    nopMass = nucleusOutputParticle->GetPDGMass();

    G4String                  title;
    Int_t                     nBinsX;
    Int_t                     nBinsY;
    Double_t                  nBinsMinX;
    Double_t                  nBinsMaxX;
    Double_t                  nBinsMinY;
    Double_t                  nBinsMaxY;
    CexmcHistoAxes            axes;

    if ( runManager->ProjectIsSaved() )
    {
        G4String  projectsDir( runManager->GetProjectsDir() );
        G4String  resultsFile( projectsDir + "/" + runManager->GetProjectId() +
                               ".root" );
        outFile = new TFile( resultsFile, "recreate" );
    }
    else
    {
        outFile = new TDirectoryFile( CexmcHistoDirectoryHandle,
                                      CexmcHistoDirectoryTitle );
        gDirectory->cd( CexmcHistoDirectoryHandle );
    }

    if ( ! outFile )
        throw CexmcException( CexmcWeirdException );

    const CexmcSetup *  setup( static_cast< const CexmcSetup * >(
                                runManager->GetUserDetectorConstruction() ) );
    if ( ! setup )
        throw CexmcException( CexmcWeirdException );

    const G4LogicalVolume *  lVolume( setup->GetVolume( CexmcSetup::Monitor ) );

    if ( ! lVolume )
        throw CexmcException( CexmcIncompatibleGeometry );

    nBinsMinX = CexmcHistoBeamMomentumMin;
    nBinsMaxX = CexmcHistoBeamMomentumMax;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) /
                    CexmcHistoBeamMomentumResolution );
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcMomentumBP_TPT_Histo, Cexmc_TH1F, false,
          false, CexmcTPT, "mombp", "Beam momentum at the monitor", axes ) );
    AddHisto( CexmcHistoData( CexmcMomentumBP_RT_Histo, Cexmc_TH1F, false,
          false, CexmcRT, "mombp", "Beam momentum at the monitor", axes ) );
    if ( verboseLevel > 0 )
    {
        AddHisto( CexmcHistoData( CexmcMomentumIP_TPT_Histo, Cexmc_TH1F, false,
            false, CexmcTPT, "momip", "Momentum of the incident particle",
            axes ) );
    }

    G4Box *   box( dynamic_cast< G4Box * >( lVolume->GetSolid() ) );

    if ( ! box )
        throw CexmcException( CexmcIncompatibleGeometry );

    G4double  width( box->GetXHalfLength() * 2 );
    G4double  height( box->GetYHalfLength() * 2 );
    G4double  halfWidth( width / 2 + CexmcHistoTPSafetyArea );
    G4double  halfHeight( height / 2 + CexmcHistoTPSafetyArea );

    nBinsX = Int_t( halfWidth * 2 / CexmcHistoTPResolution );
    nBinsY = Int_t( halfHeight * 2 / CexmcHistoTPResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, -halfWidth, halfWidth ) );
    axes.push_back( CexmcHistoAxisData( nBinsY, -halfHeight, halfHeight ) );
    AddHisto( CexmcHistoData( CexmcTPInMonitor_TPT_Histo, Cexmc_TH2F, false,
        false, CexmcTPT, "tpmon", "Track points (mon)", axes ) );

    lVolume = setup->GetVolume( CexmcSetup::Target );
    G4Tubs *  tube( dynamic_cast< G4Tubs * >( lVolume->GetSolid() ) );

    if ( ! tube )
        throw CexmcException( CexmcIncompatibleGeometry );

    G4double  radius( tube->GetOuterRadius() );
    height = tube->GetZHalfLength() * 2;
    halfWidth = radius + CexmcHistoTPSafetyArea;
    halfHeight = height / 2 + CexmcHistoTPSafetyArea;

    nBinsX = Int_t( halfWidth * 2 / CexmcHistoTPResolution );
    nBinsY = Int_t( halfWidth * 2 / CexmcHistoTPResolution );
    Int_t  nBinsZ( Int_t( halfHeight * 2 / CexmcHistoTPResolution ) );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, -halfWidth, halfWidth ) );
    axes.push_back( CexmcHistoAxisData( nBinsY, -halfWidth, halfWidth ) );
    axes.push_back( CexmcHistoAxisData( nBinsZ, -halfHeight, halfHeight ) );
    AddHisto( CexmcHistoData( CexmcTPInTarget_TPT_Histo, Cexmc_TH3F, false,
        false, CexmcTPT, "tptar", "Track points (tar)", axes ) );
    AddHisto( CexmcHistoData( CexmcTPInTarget_RT_Histo, Cexmc_TH3F, false,
        false, CexmcRT, "tptar", "Track points (tar)", axes ) );

    title = "Reconstructed masses (" + nopName + " vs. " + opName + ")";
    nBinsMinX = opMass / 2;
    nBinsMaxX = opMass + opMass / 2;
    nBinsMinY = nopMass / 2;
    nBinsMaxY = nopMass + nopMass / 2;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoMassResolution );
    nBinsY = Int_t( ( nBinsMaxY - nBinsMinY ) / CexmcHistoMassResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    axes.push_back( CexmcHistoAxisData( nBinsY, nBinsMinY, nBinsMaxY ) );
    AddHisto( CexmcHistoData( CexmcRecMasses_EDT_Histo, Cexmc_TH2F, false,
        false, CexmcEDT, "recmasses", title, axes ) );
    AddHisto( CexmcHistoData( CexmcRecMasses_RT_Histo, Cexmc_TH2F, false,
        false, CexmcRT, "recmasses", title, axes ) );

    nBinsMinX = 0.;
    nBinsMaxX = CexmcHistoEnergyMax;
    nBinsMinY = 0.;
    nBinsMaxY = CexmcHistoEnergyMax;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoEnergyResolution );
    nBinsY = Int_t( ( nBinsMaxY - nBinsMinY ) / CexmcHistoEnergyResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    axes.push_back( CexmcHistoAxisData( nBinsY, nBinsMinY, nBinsMaxY ) );
    AddHisto( CexmcHistoData( CexmcAbsorbedEnergy_EDT_Histo, Cexmc_TH2F, false,
        false, CexmcEDT, "ae", "Absorbed energy (rc vs. lc)", axes ) );
    AddHisto( CexmcHistoData( CexmcAbsorbedEnergy_RT_Histo, Cexmc_TH2F, false,
        false, CexmcRT, "ae", "Absorbed energy (rc vs. lc)", axes ) );

    SetupARHistos( runManager->GetPhysicsManager()->GetProductionModel()->
                   GetAngularRanges() );

    isInitialized = true;
}


void  CexmcHistoManager::SetupARHistos( const CexmcAngularRangeList &  aRanges )
{
    TIter      objs( gDirectory->GetList() );
    TObject *  obj( NULL );

    while ( ( obj = ( TObject * )objs() ) )
    {
        TString   name( obj->GetName() );

        if ( obj->IsFolder() && ( name.Contains( TString( "_arreal_" ) ) ||
                                  name.Contains( TString( "_arrec_" ) ) ) )
        {
            gDirectory->cd( name );
            gDirectory->DeleteAll();
            gDirectory->cd( ".." );
        }
    }

    for ( CexmcHistosMap::iterator  k( histos.begin() ); k != histos.end();
                                                                        ++k )
    {
        if ( k->second.empty() )
            continue;

        if ( k->first >= CexmcHistoType_ARReal_START &&
             k->first <= CexmcHistoType_ARReal_END )
        {
            k->second.clear();
        }
    }

    for ( CexmcAngularRangeList::const_iterator  k( aRanges.begin() );
                                                    k != aRanges.end(); ++k )
    {
        AddARHistos( *k );
    }
}


void  CexmcHistoManager::AddARHistos( const CexmcAngularRange &  aRange )
{
    G4String                  title;
    Int_t                     nBinsX;
    Double_t                  nBinsMinX;
    Double_t                  nBinsMaxX;
    CexmcHistoAxes            axes;

    title = "Reconstructed mass of " + opName;
    nBinsMinX = opMass / 2;
    nBinsMaxX = opMass + opMass / 2;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoMassResolution );
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcRecMassOP_ARReal_RT_Histo, Cexmc_TH1F, true,
        false, CexmcRT, "recmassop", title, axes ), aRange );

    title = "Reconstructed mass of " + nopName;
    nBinsMinX = nopMass / 2;
    nBinsMaxX = nopMass + nopMass / 2;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoMassResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcRecMassNOP_ARReal_RT_Histo, Cexmc_TH1F, true,
        false, CexmcRT, "recmassnop", title, axes ), aRange );

    G4RunManager *      runManager( G4RunManager::GetRunManager() );
    const CexmcSetup *  setup( static_cast< const CexmcSetup * >(
                                runManager->GetUserDetectorConstruction() ) );
    if ( ! setup )
        throw CexmcException( CexmcWeirdException );

    const G4LogicalVolume *  lVolume( setup->GetVolume(
                                                    CexmcSetup::Calorimeter ) );

    G4Box *   box( dynamic_cast< G4Box * >( lVolume->GetSolid() ) );

    if ( ! box )
        throw CexmcException( CexmcIncompatibleGeometry );

    G4double  width( box->GetXHalfLength() * 2 );
    G4double  height( box->GetYHalfLength() * 2 );
    G4double  halfWidth( width / 2 + CexmcHistoTPSafetyArea );
    G4double  halfHeight( height / 2 + CexmcHistoTPSafetyArea );

    nBinsX = Int_t( halfWidth * 2 / CexmcHistoTPResolution );
    Int_t  nBinsY( Int_t( halfHeight * 2 / CexmcHistoTPResolution ) );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, -halfWidth, halfWidth ) );
    axes.push_back( CexmcHistoAxisData( nBinsY, -halfHeight, halfHeight ) );

    /* looks like there is no possibility to draw descending xaxis in ROOT,
     * so imagine that you look at calorimeters from behind, i.e. your face to
     * the beam */
    AddHisto( CexmcHistoData( CexmcOPDPAtLeftCalorimeter_ARReal_EDT_Histo,
        Cexmc_TH2F, true, false, CexmcEDT, "opdpcl",
        "Gamma position on the surface (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcOPDPAtRightCalorimeter_ARReal_EDT_Histo,
        Cexmc_TH2F, true, false, CexmcEDT, "opdpcr",
        "Gamma position on the surface (rc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcOPDPAtLeftCalorimeter_ARReal_RT_Histo,
        Cexmc_TH2F, true, false, CexmcRT, "opdpcl",
        "Gamma position on the surface (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcOPDPAtRightCalorimeter_ARReal_RT_Histo,
        Cexmc_TH2F, true, false, CexmcRT, "opdpcr",
        "Gamma position on the surface (rc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcRecOPDPAtLeftCalorimeter_ARReal_EDT_Histo,
        Cexmc_TH2F, true, false, CexmcEDT, "recopdpcl",
        "Reconstructed gamma position on the surface (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcRecOPDPAtRightCalorimeter_ARReal_EDT_Histo,
        Cexmc_TH2F, true, false, CexmcEDT, "recopdpcr",
        "Reconstructed gamma position on the surface (rc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcRecOPDPAtLeftCalorimeter_ARReal_RT_Histo,
        Cexmc_TH2F, true, false, CexmcRT, "recopdpcl",
        "Reconstructed gamma position on the surface (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcRecOPDPAtRightCalorimeter_ARReal_RT_Histo,
        Cexmc_TH2F, true, false, CexmcRT, "recopdpcr",
        "Reconstructed gamma position on the surface (rc)", axes ), aRange );

    nBinsMinX = 0.;
    nBinsMaxX = CexmcHistoEnergyMax;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoEnergyResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcKinEnAtLeftCalorimeter_ARReal_TPT_Histo,
        Cexmc_TH1F, true, false, CexmcTPT, "kecl",
        "Kinetic energy of gamma (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcKinEnAtRightCalorimeter_ARReal_TPT_Histo,
        Cexmc_TH1F, true, false, CexmcTPT, "kecr",
        "Kinetic energy of gamma (rc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcKinEnAtLeftCalorimeter_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "kecl",
        "Kinetic energy of gamma (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcKinEnAtRightCalorimeter_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "kecr",
        "Kinetic energy of gamma (rc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcAbsEnInLeftCalorimeter_ARReal_EDT_Histo,
        Cexmc_TH1F, true, false, CexmcEDT, "aecl",
        "Absorbed energy (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcAbsEnInRightCalorimeter_ARReal_EDT_Histo,
        Cexmc_TH1F, true, false, CexmcEDT, "aecr",
        "Absorbed energy (rc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcAbsEnInLeftCalorimeter_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "aecl",
        "Absorbed energy (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcAbsEnInRightCalorimeter_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "aecr",
        "Absorbed energy (rc)", axes ), aRange );

    nBinsMinX = CexmcHistoMissEnergyMin;
    nBinsMaxX = CexmcHistoMissEnergyMax;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) /
                    CexmcHistoMissEnergyResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcMissEnFromLeftCalorimeter_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "mecl",
        "Missing energy (lc)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcMissEnFromRightCalorimeter_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "mecr",
        "Missing energy (rc)", axes ), aRange );

    title = "Kinetic energy of newborn " + opName + " (lab)";
    nBinsMinX = 0.;
    nBinsMaxX = CexmcHistoEnergyMax;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoEnergyResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcKinEnOP_LAB_ARReal_TPT_Histo,
        Cexmc_TH1F, true, false, CexmcTPT, "keop_lab", title, axes ), aRange );
    AddHisto( CexmcHistoData( CexmcKinEnOP_LAB_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "keop_lab", title, axes ), aRange );

    title = "Angle of newborn " + opName + " (scm)";
    nBinsMinX = -1.0;
    nBinsMaxX = 1.0;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoAngularCResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcAngleOP_SCM_ARReal_TPT_Histo,
        Cexmc_TH1F, true, false, CexmcTPT, "aop_scm", title, axes ), aRange );
    AddHisto( CexmcHistoData( CexmcAngleOP_SCM_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "aop_scm", title, axes ), aRange );

    title = "Reconstruced angle of newborn " + opName + " (scm)";
    AddHisto( CexmcHistoData( CexmcRecAngleOP_SCM_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "recaop_scm", title, axes ), aRange );

    title = "Real - reconstruced angle of newborn " + opName + " (scm)";
    AddHisto( CexmcHistoData( CexmcDiffAngleOP_SCM_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "diffaop_scm", title, axes ),
        aRange );

    nBinsMinX = 0.;
    nBinsMaxX = 360.;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoAngularResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcOpenAngle_ARReal_TPT_Histo,
        Cexmc_TH1F, true, false, CexmcTPT, "oa",
        "Open angle between the gammas", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcOpenAngle_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "oa",
        "Open angle between the gammas", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcRecOpenAngle_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "recoa",
        "Reconstructed open angle between the gammas", axes ), aRange );

    nBinsMinX = -180.;
    nBinsMaxX = 180.;
    nBinsX = Int_t( ( nBinsMaxX - nBinsMinX ) / CexmcHistoAngularResolution );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, nBinsMinX, nBinsMaxX ) );
    AddHisto( CexmcHistoData( CexmcDiffOpenAngle_ARReal_RT_Histo,
        Cexmc_TH1F, true, false, CexmcRT, "diffoa",
        "Real - reconstructed open angle between the gammas", axes ), aRange );

    lVolume = setup->GetVolume( CexmcSetup::Target );
    G4Tubs *  tube( dynamic_cast< G4Tubs * >( lVolume->GetSolid() ) );

    if ( ! tube )
        throw CexmcException( CexmcIncompatibleGeometry );

    G4double  radius( tube->GetOuterRadius() );
    height = tube->GetZHalfLength() * 2;
    halfWidth = radius + CexmcHistoTPSafetyArea;
    halfHeight = height / 2 + CexmcHistoTPSafetyArea;

    nBinsX = Int_t( halfWidth * 2 / CexmcHistoTPResolution );
    nBinsY = Int_t( halfWidth * 2 / CexmcHistoTPResolution );
    Int_t  nBinsZ( Int_t( halfHeight * 2 / CexmcHistoTPResolution ) );
    axes.clear();
    axes.push_back( CexmcHistoAxisData( nBinsX, -halfWidth, halfWidth ) );
    axes.push_back( CexmcHistoAxisData( nBinsY, -halfWidth, halfWidth ) );
    axes.push_back( CexmcHistoAxisData( nBinsZ, -halfHeight, halfHeight ) );
    AddHisto( CexmcHistoData( CexmcTPInTarget_ARReal_TPT_Histo, Cexmc_TH3F,
        true, false, CexmcTPT, "tptar", "Track points (tar)", axes ), aRange );
    AddHisto( CexmcHistoData( CexmcTPInTarget_ARReal_RT_Histo, Cexmc_TH3F,
        true, false, CexmcRT, "tptar", "Track points (tar)", axes ), aRange );
}


void  CexmcHistoManager::Add( CexmcHistoType  histoType, unsigned int  index,
                              G4double  x )
{
    CexmcHistosMap::iterator  found( histos.find( histoType ) );
    if ( found == histos.end() || histos[ histoType ].size() <= index )
        throw CexmcException( CexmcWeirdException );

    histos[ histoType ][ index ]->Fill( x );
}


void  CexmcHistoManager::Add( CexmcHistoType  histoType, unsigned int  index,
                              G4double  x, G4double  y )
{
    CexmcHistosMap::iterator  found( histos.find( histoType ) );
    if ( found == histos.end() || histos[ histoType ].size() <= index )
        throw CexmcException( CexmcWeirdException );

    /* no cast needed because TH1 has virtual method
     * Fill( Double_t, Double_t ) */
    histos[ histoType ][ index ]->Fill( x, y );
}


void  CexmcHistoManager::Add( CexmcHistoType  histoType, unsigned int  index,
                              G4double  x, G4double  y, G4double  z )
{
    CexmcHistosMap::iterator  found( histos.find( histoType ) );
    if ( found == histos.end() || histos[ histoType ].size() <= index )
        throw CexmcException( CexmcWeirdException );

    /* cast needed because TH1 does not have virtual method
     * Fill( Double_t, Double_t, Double_t ) */
    TH3 *  histo( static_cast< TH3 * >( histos[ histoType ][ index ] ) );

    histo->Fill( x, y, z );
}


void  CexmcHistoManager::Add( CexmcHistoType  histoType, unsigned int  index,
                              G4int  binX, G4int  binY, G4double  value )
{
    CexmcHistosMap::iterator  found( histos.find( histoType ) );
    if ( found == histos.end() || histos[ histoType ].size() <= index )
        throw CexmcException( CexmcWeirdException );

    ++binX;
    ++binY;
    Double_t  curValue( histos[ histoType ][ index ]->GetBinContent(
                                                                binX, binY ) );
    histos[ histoType ][ index ]->SetBinContent( binX, binY,
                                                     curValue + value / GeV );
}


void  CexmcHistoManager::List( void ) const
{
    /* BEWARE: list will be printed on stdout */
    gDirectory->ls();
}


void  CexmcHistoManager::Print( const G4String &  value )
{
    TObject *  histo( gDirectory->FindObjectAny( value.c_str() ) );

    if ( ! histo )
    {
        G4cout << "Histogram '" << value << "' was not found" << G4endl;
        return;
    }

    /* BEWARE: histo will be printed on stdout */
    histo->Print( "range" );
}


#ifdef CEXMC_USE_ROOTQT

void  CexmcHistoManager::Draw( const G4String &  histoName,
                               const G4String &  histoDrawOptions )
{
    if ( ! areLiveHistogramsEnabled )
    {
        G4cout << "Live histograms option is disabled" << G4endl;
        return;
    }

    TObject *  histo( gDirectory->FindObjectAny( histoName ) );

    if ( ! histo )
    {
        G4cout << "Histogram '" << histoName << "' was not found" << G4endl;
        return;
    }

    if ( ! rootCanvas )
    {
        /* save default application font because rootCanvas will break it */
        QFont  defaultAppFont( QApplication::font() );
        rootCanvas = new CexmcHistoWidget;
        QApplication::setFont( defaultAppFont );
        rootCanvas->resize( CexmcHistoCanvasWidth, CexmcHistoCanvasHeight );
        rootCanvas->GetCanvas()->cd();
    }

    histo->Draw( histoDrawOptions );
    rootCanvas->show();
    rootCanvas->GetCanvas()->Update();
}


void  CexmcHistoManager::EnableLiveHistograms( G4UIsession *  session,
                                               G4bool  on )
{
    areLiveHistogramsEnabled = on;

    if ( ! on || ! isInitialized )
        return;

    G4UIQt *  qtSession( dynamic_cast< G4UIQt * >( session ) );

    if ( ! qtSession )
        return;

    if ( ! histoMenuHandle.empty() && ! isHistoMenuInitialized )
    {
        qtSession->AddMenu( histoMenuHandle, histoMenuLabel );
        BuildMenuTree( qtSession, histoMenuHandle, gDirectory->GetList() );
        isHistoMenuInitialized = true;
    }
}


void  CexmcHistoManager::BuildMenuTree( G4UIQt *  session,
                                        const G4String &  menu, TList *  ls )
{
    TIter      objs( ls );
    TObject *  obj( NULL );

    while ( ( obj = ( TObject * )objs() ) )
    {
        G4String  name( obj->GetName() );
        G4String  title( obj->GetTitle() );

        if ( obj->IsFolder() )
        {
            AddSubmenu( session, menu, name, title );
            BuildMenuTree( session, name, ( ( TDirectory * )obj )->GetList() );
        }
        else
        {
            G4String  options( name );

            do
            {
                if ( obj->InheritsFrom( TH3::Class() ) &&
                     ! drawOptions3D.empty() )
                {
                    title = G4String( "3: " ) + title;
                    options += G4String( " " ) + drawOptions3D;
                    break;
                }
                if ( obj->InheritsFrom( TH2::Class() ) &&
                     ! drawOptions2D.empty() )
                {
                    title = G4String( "2: " ) + title;
                    options += G4String( " " ) + drawOptions2D;
                    break;
                }
                if ( obj->InheritsFrom( TH1::Class() ) &&
                     ! drawOptions1D.empty() )
                {
                    options += G4String( " " ) + drawOptions1D;
                    break;
                }
            } while ( false );

            G4String  cmd( CexmcMessenger::histoDirName + "draw " + options );
            session->AddButton( menu, title.c_str(), cmd );
        }
    }
}


void  CexmcHistoManager::AddSubmenu( G4UIQt *  session,
                                     const G4String &  parent,
                                     const G4String &  name,
                                     const G4String &  label )
{
  QMenu *  menu( new QMenu( label.c_str() ) );
  QMenu *  parentMenu( ( QMenu * )session->GetInteractor( parent ) );

  parentMenu->addMenu( menu );
  session->AddInteractor( name, ( G4Interactor )menu );
}

#endif

#endif


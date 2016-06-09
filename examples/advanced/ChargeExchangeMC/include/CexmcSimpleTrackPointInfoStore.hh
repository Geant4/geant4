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
 *       Filename:  CexmcSimpleTrackPointInfoStore.hh
 *
 *    Description:  serialization helper for track point info objects
 *
 *        Version:  1.0
 *        Created:  31.12.2009 13:55:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_SIMPLE_TRACK_POINT_INFO_STORE_HH
#define CEXMC_SIMPLE_TRACK_POINT_INFO_STORE_HH

#ifdef CEXMC_USE_PERSISTENCY

#include <boost/serialization/access.hpp>
#include "CexmcSimpleThreeVectorStore.hh"
#include "CexmcCommon.hh"

class  CexmcTrackPointInfo;


class  CexmcSimpleTrackPointInfoStore
{
    friend class  boost::serialization::access;
#ifdef CEXMC_USE_CUSTOM_FILTER
    friend class  CexmcASTEval;
#endif

    public:
        CexmcSimpleTrackPointInfoStore();

        CexmcSimpleTrackPointInfoStore( const CexmcTrackPointInfo &  tpInfo );

    public:
        operator CexmcTrackPointInfo() const;

    private:
        template  < typename  Archive >
        void  serialize( Archive &  archive, const unsigned int  version );

    private:
        CexmcSimpleThreeVectorStore  positionLocal;

        CexmcSimpleThreeVectorStore  positionWorld;

        CexmcSimpleThreeVectorStore  directionLocal;

        CexmcSimpleThreeVectorStore  directionWorld;

        G4double                     momentumAmp;

        G4int                        particlePDGEncoding;

        G4int                        trackId;

        CexmcTrackType               trackType;
};


template  < typename  Archive >
void  CexmcSimpleTrackPointInfoStore::serialize( Archive &  archive,
                                                 const unsigned int )
{
    archive & positionLocal;
    archive & positionWorld;
    archive & directionLocal;
    archive & directionWorld;
    archive & momentumAmp;
    archive & particlePDGEncoding;
    archive & trackId;
    archive & trackType;
}

#endif

#endif


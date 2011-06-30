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
# <<BEGIN-copyright>>
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# <<END-copyright>>
*/
#ifndef G4GIDI_target_h_included
#define G4GIDI_target_h_included 1

#include <vector>
#include <string>

#include <tpia_target.h>

typedef struct crossSectionData_s crossSectionData;
typedef struct channelID_s channelID;
typedef struct G4GIDI_Product_s G4GIDI_Product;

struct crossSectionData_s {
    int start, end;
    std::vector<double> crossSection;
};

struct channelID_s {
    std::string ID;
};

//typedef struct G4GIDI_Product_s {
struct G4GIDI_Product_s {
    int A, Z, m;
    double kineticEnergy, px, py, pz;
};

class G4GIDI_target {

    public:
        void init( const char *fileName );
        std::string equalProbableBinSampleMethod;
        int nElasticIndices, nCaptureIndices, nFissionIndices, nOthersIndices;
        int *elasticIndices, *captureIndices, *fissionIndices, *othersIndices;

    public:
        GIDI::statusMessageReporting smr;
        std::string name;
        std::string sourceFilename;
        double mass;
        GIDI::tpia_target *target;

        G4GIDI_target( const char *fileName );       
        G4GIDI_target( std::string &fileName );       
        ~G4GIDI_target( );

        std::string *getName( void );
        std::string *getFilename( void );
        int getZ( void );
        int getA( void );
        int getM( void );
        double getMass( void );
        int getTemperatures( double *temperatures );
        int readTemperature( int index );
        std::string getEqualProbableBinSampleMethod( void );
        int setEqualProbableBinSampleMethod( std::string method );

        int getNumberOfChannels( void );
        int getNumberOfProductionChannels( void );
        std::vector<channelID> *getChannelIDs( void );
        std::vector<channelID> *getProductionChannelIDs( void );
        std::vector<channelID> *getChannelIDs2( GIDI::tpia_channel **channels, int n );

        std::vector<double> *getEnergyGridAtTIndex( int index );

        double getTotalCrossSectionAtE( double e_in, double temperature );
        double getElasticCrossSectionAtE( double e_in, double temperature );
        double getCaptureCrossSectionAtE( double e_in, double temperature );
        double getFissionCrossSectionAtE( double e_in, double temperature );
        double getOthersCrossSectionAtE( double e_in, double temperature );
        double sumChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature );
        int sampleChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature, double (*rng)( void * ), void *rngState );

        double getElasticFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        std::vector<G4GIDI_Product> *getCaptureFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        std::vector<G4GIDI_Product> *getFissionFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        std::vector<G4GIDI_Product> *getOthersFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        std::vector<G4GIDI_Product> *getFinalState( int nIndices, int *indices, double e_in, double temperature, double (*rng)( void * ), void *rngState );
};

#endif      // End of G4GIDI_target_h_included

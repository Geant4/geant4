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
#ifndef GIDI4GEANT_target_h_included
#define GIDI4GEANT_target_h_included 1

#include <vector>
#include <string>
//using namespace std;

#include <tpia_target.h>

typedef struct crossSectionData_s crossSectionData;
typedef struct channelID_s channelID;
typedef struct GIDI4GEANT_Product_s GIDI4GEANT_Product;

struct crossSectionData_s {
    int start, end;
    vector<double> crossSection;
};

struct channelID_s {
    string ID;
};

//typedef struct GIDI4GEANT_Product_s {
struct GIDI4GEANT_Product_s {
    int A, Z, m;
    double kineticEnergy, px, py, pz;
};

class GIDI4GEANT_target {

    public:
        void init( const char *fileName );
        string equalProbableBinSampleMethod;
        int nElasticIndices, nCaptureIndices, nFissionIndices, nOthersIndices;
        int *elasticIndices, *captureIndices, *fissionIndices, *othersIndices;

    public:
        statusMessageReporting smr;
        string name;
        string sourceFilename;
        double mass;
        tpia_target *target;

        GIDI4GEANT_target( const char *fileName );       
        GIDI4GEANT_target( string &fileName );       
        ~GIDI4GEANT_target( );

        string *getName( void );
        string *getFilename( void );
        int getZ( void );
        int getA( void );
        int getM( void );
        double getMass( void );
        int getTemperatures( double *temperatures );
        int readTemperature( int index );
        string getEqualProbableBinSampleMethod( void );
        int setEqualProbableBinSampleMethod( string method );

        int getNumberOfChannels( void );
        int getNumberOfProductionChannels( void );
        vector<channelID> *getChannelIDs( void );
        vector<channelID> *getProductionChannelIDs( void );
        vector<channelID> *getChannelIDs2( tpia_channel **channels, int n );

        vector<double> *getEnergyGridAtTIndex( int index );

        double getTotalCrossSectionAtE( double e_in, double temperature );
        double getElasticCrossSectionAtE( double e_in, double temperature );
        double getCaptureCrossSectionAtE( double e_in, double temperature );
        double getFissionCrossSectionAtE( double e_in, double temperature );
        double getOthersCrossSectionAtE( double e_in, double temperature );
        double sumChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature );
        int sampleChannelCrossSectionAtE( int nIndices, int *indices, double e_in, double temperature, double (*rng)( void * ), void *rngState );

        double getElasticFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        vector<GIDI4GEANT_Product> *getCaptureFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        vector<GIDI4GEANT_Product> *getFissionFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        vector<GIDI4GEANT_Product> *getOthersFinalState( double e_in, double temperature, double (*rng)( void * ), void *rngState );
        vector<GIDI4GEANT_Product> *getFinalState( int nIndices, int *indices, double e_in, double temperature, double (*rng)( void * ), void *rngState );
};

#endif      // End of GIDI4GEANT_target_h_included

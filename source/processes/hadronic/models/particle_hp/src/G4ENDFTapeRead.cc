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
 * File:   G4ENDFTapeRead.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on September 6, 2011, 10:01 AM
 */

#include <fstream>
#include <map>
#include <vector>

#include "globals.hh"
#include "G4ParticleHPManager.hh"

#include "G4ENDFTapeRead.hh"
#include "G4ENDFYieldDataContainer.hh"
#include "G4FFGDebuggingMacros.hh"
#include "G4FFGDefaultValues.hh"
#include "G4FFGEnumerations.hh"
#include "G4TableTemplate.hh"

G4ENDFTapeRead::
G4ENDFTapeRead( G4String FileLocation,
                G4String FileName,
                G4FFGEnumerations::YieldType WhichYield,
                G4FFGEnumerations::FissionCause /*WhichCause*/ )
:   /* Cause_(WhichCause),*/
    Verbosity_(G4FFGDefaultValues::Verbosity),
    YieldType_(WhichYield)
{
    // Initialize the class
    Initialize(FileLocation + FileName);
}

G4ENDFTapeRead::
G4ENDFTapeRead( G4String FileLocation,
                G4String FileName,
                G4FFGEnumerations::YieldType WhichYield,
                G4FFGEnumerations::FissionCause /*WhichCause*/,
                G4int Verbosity )
:   /*Cause_(WhichCause),*/
    Verbosity_(Verbosity),
    YieldType_(WhichYield)
{
    // Initialize the class
    Initialize(FileLocation + FileName);
}

G4ENDFTapeRead::
G4ENDFTapeRead( std::istringstream& dataStream,
                G4FFGEnumerations::YieldType WhichYield,
                G4FFGEnumerations::FissionCause /*WhichCause*/,
                G4int Verbosity )
:   /*Cause_(WhichCause),*/
    Verbosity_(Verbosity),
    YieldType_(WhichYield)
{
    // Initialize the class
    Initialize(dataStream);
}

void G4ENDFTapeRead::
Initialize( G4String dataFile )
{
    std::istringstream dataStream(std::ios::in);
    G4ParticleHPManager::GetInstance()->GetDataStream(dataFile, dataStream);

    Initialize(dataStream);
}

void G4ENDFTapeRead::
Initialize( std::istringstream& dataStream )
{    
G4FFG_FUNCTIONENTER__

    EnergyGroups_ = 0;
    EnergyGroupValues_ = NULL;

    YieldContainerTable_ = new G4TableTemplate< G4ENDFYieldDataContainer >;

    try
    {
        ReadInData(dataStream);
    } catch (std::exception & e)
    {
        delete YieldContainerTable_;

        G4FFG_FUNCTIONLEAVE__
        throw e;
    }

G4FFG_FUNCTIONLEAVE__
}

G4double* G4ENDFTapeRead::
G4GetEnergyGroupValues( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return EnergyGroupValues_;
}

G4int G4ENDFTapeRead::
G4GetNumberOfEnergyGroups( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return EnergyGroups_;
}

G4int G4ENDFTapeRead::
G4GetNumberOfFissionProducts( void )
{
G4FFG_FUNCTIONENTER__

    G4int NumberOfElements = (G4int)YieldContainerTable_->G4GetNumberOfElements();

G4FFG_FUNCTIONLEAVE__
    return NumberOfElements;
}

G4ENDFYieldDataContainer* G4ENDFTapeRead::
G4GetYield( G4int WhichYield )
{
G4FFG_DATA_FUNCTIONENTER__

    G4ENDFYieldDataContainer* Container = NULL;
    if(WhichYield >= 0 && WhichYield < YieldContainerTable_->G4GetNumberOfElements())
    {
        Container = YieldContainerTable_->G4GetContainer(WhichYield);
    }

G4FFG_DATA_FUNCTIONLEAVE__
    return Container;
}

void G4ENDFTapeRead::
G4SetVerbosity(G4int WhatVerbosity)
{
G4FFG_FUNCTIONENTER__

    this->Verbosity_ = WhatVerbosity;

G4FFG_FUNCTIONLEAVE__
}

void G4ENDFTapeRead::
ReadInData( std::istringstream& dataStream )
{
G4FFG_FUNCTIONENTER__

    // Check if the file exists
    if(!dataStream.good())
    {
        G4Exception("G4ENDFTapeRead::ReadInData()",
                    "Illegal file name",
                    JustWarning,
                    "Fission product data not available");

        // TODO create/use a specialized exception
        G4FFG_FUNCTIONLEAVE__
        throw std::exception();
    }

    // Code to read in from a pure ENDF data tape.
    // Commented out while pre-formatted Geant4 ENDF data is being used
//    G4int CurrentEnergyGroup = -1;
//    std::vector< G4double > NewDoubleVector;
//    std::vector< G4double > EnergyPoints;
//    std::vector< G4int > Product;
//    std::vector< G4FFGEnumerations::MetaState > MetaState;
//    std::vector< std::vector< G4double > > Yield;
//    std::vector< std::vector< G4double > > Error;
//    G4String DataBlock;
//    std::size_t InsertExponent;
//    G4double Parts[6];
//    G4double dummy;
//    G4int MAT;
//    G4int MF;
//    G4int MT;
//    G4int LN;
//    G4int Block;
//    G4int EmptyProduct;
//    G4int Location;
//    G4int ItemCounter = 0;
//    G4int FirstLineInEnergyGroup = 0;
//    G4int LastLineInEnergyGroup = 0;
//    G4bool FoundEnergyGroup = false;
//    G4bool FoundPID = false;
//
//    while(getline(DataFile, Temp))
//    {
//        // Format the string so that it can be interpreted correctly
//        DataBlock = Temp.substr(0, 66);
//        Temp = Temp.substr(66);
//        InsertExponent = 0;
//        while((InsertExponent = DataBlock.find_first_of("-+", InsertExponent)) != G4String::npos)
//        {
//            DataBlock.insert(InsertExponent, 1, 'e');
//            InsertExponent += 2;
//        }
//        sscanf(DataBlock.c_str(), "%11le %11le %11le %11le %11le %11le",
//            &Parts[0], &Parts[1], &Parts[2], &Parts[3], &Parts[4], &Parts[5]);
//        sscanf(Temp.substr(0, 4).c_str(), "%i", &MAT);
//        sscanf(Temp.substr(4, 2).c_str(), "%i", &MF);
//        sscanf(Temp.substr(6, 3).c_str(), "%i", &MT);
//        sscanf(Temp.substr(9).c_str(), "%i", &LN);
//
//        if(MT == YieldType_)
//        {
//            if(LN == 1)
//            {
//                if(FoundPID != true)
//                {
//                    // The first line of an ENDF section for MT = 454 or 459
//                    // always contains the parent PID
//                    // This section can potentially be expanded to check and
//                    // verify that it is the correct nucleus
//                    FoundPID = true;
//
//                    continue;
//                }
//            } else if(FoundPID == true && FoundEnergyGroup == false)
//            {
//                // Skip this line if it is not the energy definition line
//                if(Parts[1] != 0 || Parts[3] != 0)
//                {
//                    continue;
//                }
//
//                // The first block is the incident neutron energy
//                // information.
//                // Check to make sure that it is spontaneous or neutron
//                // induced.
//                if(Cause_ == G4FFGEnumerations::NEUTRON_INDUCED)
//                {
//                    if(Parts[0] != 0)
//                    {
//                        FoundEnergyGroup = true;
//                    }
//                } else if(Cause_ == G4FFGEnumerations::SPONTANEOUS)
//                {
//                    if(Parts[0] == 0)
//                    {
//                        FoundEnergyGroup = true;
//                    }
//                } else
//                { // Maybe more fission causes here if added later
//                    FoundEnergyGroup = false;
//                }
//
//                if(FoundEnergyGroup == true)
//                {
//                    // Convert to eV
//                    Parts[0] *= eV;
//
//                    // Calculate the parameters
//                    FirstLineInEnergyGroup = LN;
//                    LastLineInEnergyGroup = FirstLineInEnergyGroup +
//                        ceil(Parts[4] / 6.0);
//                    ItemCounter = 0;
//                    EmptyProduct = 0;
//
//                    // Initialize the data storage
//                    CurrentEnergyGroup++;
//                    EnergyPoints.push_back(Parts[0]);
//                    Yield.push_back(NewDoubleVector);
//                    Yield.back().resize(Product.size(), 0);
//                    Error.push_back(NewDoubleVector);
//                    Error.back().resize(Product.size(), 0);
//
//                    continue;
//                }
//            }
//
//            if(LN > FirstLineInEnergyGroup && LN <= LastLineInEnergyGroup)
//            {
//                for(Block = 0; Block < 6; Block++)
//                {
//                    if(EmptyProduct > 0)
//                    {
//                        EmptyProduct--;
//
//                        continue;
//                    }
//                    switch(ItemCounter % 4)
//                    {
//                        case 0:
//                            // Determine if the block is empty
//                            if(Parts[Block] == 0)
//                            {
//                                EmptyProduct = 3;
//
//                                continue;
//                            }
//
//                            // Determine if this product is already loaded
//                            for(Location = 0; Location < (signed)Product.size(); Location++)
//                            {
//                                if(Parts[Block] == Product.at(Location) &&
//                                   Parts[Block + 1] == MetaState.at(Location))
//                                {
//                                    break;
//                                }
//                            }
//
//                            // The product hasn't been created yet
//                            // Add it and initialize the other vectors
//                            if(Location == (signed)Product.size())
//                            {
//                                Product.push_back(Parts[Block]);
//                                MetaState.push_back((G4FFGEnumerations::MetaState)Parts[Block + 1]);
//                                Yield.at(CurrentEnergyGroup).push_back(0.0);
//                                Error.at(CurrentEnergyGroup).push_back(0.0);
//                            }
//                            break;
//
//                        case 2:
//                            Yield.at(CurrentEnergyGroup).at(Location) = Parts[Block];
//                            break;
//
//                        case 3:
//                            Error.at(CurrentEnergyGroup).at(Location) = Parts[Block];
//                            break;
//                    }
//
//                    ItemCounter++;
//                }
//            }
//
//            if (LN == LastLineInEnergyGroup)
//            {
//                FoundEnergyGroup = false;
//            }
//        }
//    }
//
//    G4ENDFYieldDataContainer* NewDataContainer;
//    EnergyGroups_ = EnergyPoints.size();
//    EnergyGroupValues_ = new G4double[EnergyGroups_];
//    G4int NewProduct;
//    G4FFGEnumerations::MetaState NewMetaState;
//    G4double* NewYield = new G4double[EnergyGroups_];
//    G4double* NewError = new G4double[EnergyGroups_];
//
//    for(G4int i = 0; i < EnergyGroups_; i++)
//    {
//        // Load the energy values
//        EnergyGroupValues_[i] = EnergyPoints.at(i);
//
//        // Make all the vectors the same size
//        Yield[i].resize(maxSize, 0.0);
//        Error[i].resize(maxSize, 0.0);
//    }
//
//    // Load the data into the yield table
//    for(ItemCounter = 0; ItemCounter < (signed)Product.size(); ItemCounter++)
//    {
//        NewProduct = Product.at(ItemCounter);
//        NewMetaState = MetaState.at(ItemCounter);
//
//        for(CurrentEnergyGroup = 0; CurrentEnergyGroup < EnergyGroups_; CurrentEnergyGroup++)
//        {
//            NewYield[CurrentEnergyGroup] = Yield.at(CurrentEnergyGroup).at(ItemCounter);
//            NewYield[CurrentEnergyGroup] = Error.at(CurrentEnergyGroup).at(ItemCounter);
//        }
//
//        NewDataContainer = YieldContainerTable_->G4GetNewContainer(EnergyGroups_ + 1);
//        NewDataContainer->SetProduct(NewProduct);
//        NewDataContainer->SetMetaState(NewMetaState);
//        NewDataContainer->SetYieldProbability(NewYield);
//        NewDataContainer->SetYieldError(NewError);
//    }
//
//    delete[] NewYield;
//    delete[] NewError;

    G4int MT;
    G4bool correctMT;
    G4int MF;
    G4double dummy;
    G4int blockCount;
    G4int currentEnergy = 0;
    G4double incidentEnergy;
    G4int itemCount;
    // TODO correctly implement the interpolation in the fission product yield
    G4int interpolation;
    G4int isotope;
    G4int metastate;
    G4int identifier;
    G4double yield;
    // "error" is included here in the event that errors are included in the future
    G4double error = 0.0;
    G4int maxSize = 0;

    std::vector< G4double > projectileEnergies;
    std::map< const G4int, std::pair< std::vector< G4double >, std::vector< G4double > > > intermediateData;
    std::map< const G4int, std::pair< std::vector< G4double >, std::vector< G4double > > >::iterator dataIterator;

    while(dataStream.good()) // Loop checking, 11.05.2015, T. Koi
    {
        dataStream >> MT >> MF >> dummy >> blockCount;

        correctMT = MT == YieldType_;

        for(G4int b = 0; b < blockCount; ++b)
        {
            dataStream >> incidentEnergy >> itemCount >> interpolation;
            maxSize = maxSize >= itemCount ? maxSize : itemCount;

            if(correctMT)
            {
                // Load in the energy of the projectile
                projectileEnergies.push_back(incidentEnergy);
                currentEnergy = G4int(projectileEnergies.size() - 1);
            } else
            {
                // !!! Do not break since we need to parse through the !!!
                // !!! entire data file for all possible energies      !!!
            }

            for(G4int i = 0; i < itemCount; ++i)
            {
                dataStream >> isotope >> metastate >> yield;

                if(correctMT)
                {
                    identifier = isotope * 10 + metastate;

                    dataIterator = intermediateData.insert(std::make_pair(
                            identifier,
                            std::make_pair(
                                    std::vector< G4double >(projectileEnergies.size(), 0.0),
                                    std::vector< G4double >(projectileEnergies.size(), 0.0)))).first;

                    if(dataIterator->second.first.size() < projectileEnergies.size())
                    {
                        dataIterator->second.first.resize(projectileEnergies.size());
                        dataIterator->second.second.resize(projectileEnergies.size());
                    }

                    dataIterator->second.first[currentEnergy] = yield;
                    dataIterator->second.second[currentEnergy] = error;
                } else
                {
                    // !!! Do not break since we need to parse through the !!!
                    // !!! entire data file for all possible energies      !!!
                }
            }
        }
    }

    G4ENDFYieldDataContainer* NewDataContainer;
    EnergyGroups_ = (G4int)projectileEnergies.size();
    EnergyGroupValues_ = new G4double[EnergyGroups_];
    G4int NewProduct;
    G4FFGEnumerations::MetaState NewMetaState;
    G4double* NewYield = new G4double[EnergyGroups_];
    G4double* NewError = new G4double[EnergyGroups_];

    for(G4int energyGroup = 0; energyGroup < EnergyGroups_; energyGroup++)
    {
        // Load the energy values
        EnergyGroupValues_[energyGroup] = projectileEnergies[energyGroup];
    }

    // Load the data into the yield table
    for(dataIterator = intermediateData.begin(); dataIterator != intermediateData.end(); ++dataIterator)
    {
        identifier = dataIterator->first;
        metastate = identifier % 10;
        switch(metastate)
        {
        case 1:
            NewMetaState = G4FFGEnumerations::META_1;
            break;

        case 2:
            NewMetaState = G4FFGEnumerations::META_2;
            break;

        default:
            G4Exception("G4ENDFTapeRead::ReadInData()",
                        "Unsupported state",
                        JustWarning,
                        "Unsupported metastable state supplied in fission yield data. Defaulting to the ground state");
            // Fall through
        case 0:
            NewMetaState = G4FFGEnumerations::GROUND_STATE;
            break;
        }
        NewProduct = (identifier - metastate) / 10;

        for(G4int energyGroup = 0; energyGroup < EnergyGroups_; energyGroup++)
        {
            if(energyGroup < (signed)dataIterator->second.first.size())
            {
                yield = dataIterator->second.first[energyGroup];
                error = dataIterator->second.second[energyGroup];
            } else
            {
                yield = 0.0;
                error = 0.0;
            }

            NewYield[energyGroup] = yield;
            NewError[energyGroup] = error;
        }

        NewDataContainer = YieldContainerTable_->G4GetNewContainer(EnergyGroups_);
        NewDataContainer->SetProduct(NewProduct);
        NewDataContainer->SetMetaState(NewMetaState);
        NewDataContainer->SetYieldProbability(NewYield);
        NewDataContainer->SetYieldError(NewError);
    }

    delete[] NewYield;
    delete[] NewError;

G4FFG_FUNCTIONLEAVE__
}

G4ENDFTapeRead::
~G4ENDFTapeRead( void )
{
G4FFG_FUNCTIONENTER__

    delete[] EnergyGroupValues_;
    delete YieldContainerTable_;

G4FFG_FUNCTIONLEAVE__
}


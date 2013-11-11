/*
 * File:   G4ShiftedGaussian.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on July 20, 2011, 11:55 AM
 */
 
#include <utility>
 
#include "globals.hh"

#include "G4FFGDebuggingMacros.hh"
#include "G4FFGDefaultValues.hh"
#include "G4ShiftedGaussian.hh"

G4ShiftedGaussian::
G4ShiftedGaussian( void )
{
    // Set the default verbosity
    Verbosity_ = G4FFGDefaultValues::Verbosity;
    
    // Initialize the class
    Initialize();
}

G4ShiftedGaussian::
G4ShiftedGaussian( G4int Verbosity )
{
    // Set the default verbosity
    Verbosity_ = Verbosity;
    
    // Initialize the class
    Initialize();
}

void G4ShiftedGaussian::
Initialize( void )
{
G4FFG_FUNCTIONENTER__

    // Nothing here

G4FFG_FUNCTIONLEAVE__
}

G4double G4ShiftedGaussian::
G4FindShiftedMean( G4double RequestedMean,
                   G4double RequestedStdDev)
{
G4FFG_SAMPLING_FUNCTIONENTER__

    G4int VectorSize = ShiftedMean_.size();

    for(G4int i = 0; i < VectorSize; i++)
    {
        if(ShiftedMean_[i].first.first == RequestedMean)
        {
            if(ShiftedMean_[i].first.second == RequestedStdDev)
            {
G4FFG_SAMPLING_FUNCTIONLEAVE__
                return ShiftedMean_[i].second;
            }
        }
    }

G4FFG_SAMPLING_FUNCTIONLEAVE__
    return 0;
}

void G4ShiftedGaussian::
G4InsertShiftedMean( G4double ShiftedMean,
                     G4double RequestedMean,
                     G4double RequestedStdDev)
{
G4FFG_SAMPLING_FUNCTIONENTER__

    ShiftedMean_.push_back(
        std::make_pair(
            std::make_pair(
                RequestedMean,
                RequestedStdDev),
            ShiftedMean
        )
    );

G4FFG_SAMPLING_FUNCTIONLEAVE__
    return;
}

void G4ShiftedGaussian::
G4SetVerbosity(G4int Verbosity)
{
G4FFG_SAMPLING_FUNCTIONENTER__

    Verbosity_ = Verbosity;

G4FFG_SAMPLING_FUNCTIONLEAVE__
}

G4ShiftedGaussian::
~G4ShiftedGaussian()
{
G4FFG_FUNCTIONENTER__

    // Nothing here!
G4FFG_FUNCTIONLEAVE__
}


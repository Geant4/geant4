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
 * File:   G4FPYSamplingOps.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on June 30, 2011, 11:10 AM
 */

#include <iostream>

#include <CLHEP/Random/Stat.h>
#include "Randomize.hh"
#include "globals.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

#include "G4FFGDebuggingMacros.hh"
#include "G4FFGDefaultValues.hh"
#include "G4FFGEnumerations.hh"
#include "G4FPYSamplingOps.hh"
#include "G4ShiftedGaussian.hh"
#include "G4WattFissionSpectrumValues.hh"

G4FPYSamplingOps::
G4FPYSamplingOps( void )
{
    // Set the default verbosity
    Verbosity_ = G4FFGDefaultValues::Verbosity;
    
    // Initialize the class
    Initialize();
}

G4FPYSamplingOps::
G4FPYSamplingOps( G4int Verbosity )
{
    // Set the default verbosity
    Verbosity_ = Verbosity;
    
    // Initialize the class
    Initialize();
}

void G4FPYSamplingOps::
Initialize( void )
{
G4FFG_FUNCTIONENTER__

    // Get the pointer to the random number generator
    //RandomEngine_ = CLHEP::HepRandom::getTheEngine();
    RandomEngine_ = G4Random::getTheEngine();

    // Initialize the data members
    ShiftedGaussianValues_ = new G4ShiftedGaussian;
    Mean_ = 0;
    StdDev_ = 0;
    NextGaussianIsStoredInMemory_ = FALSE;
    GaussianOne_ = 0;
    GaussianTwo_ = 0;
    Tolerance_ = 0.000001;
    WattConstants_ = new WattSpectrumConstants;
    WattConstants_->Product = 0;

G4FFG_FUNCTIONLEAVE__
}

G4double G4FPYSamplingOps::
G4SampleGaussian( G4double Mean,
                  G4double StdDev )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    // Determine if the parameters have changed
    G4bool ParametersChanged = (Mean_ != Mean || StdDev_ != StdDev);
    if(ParametersChanged == TRUE)
    {
        // Set the new values if the parameters have changed
        NextGaussianIsStoredInMemory_ = FALSE;

        Mean_ = Mean;
        StdDev_ = StdDev;
    }
    
    G4double Sample = SampleGaussian();
    
G4FFG_SAMPLING_FUNCTIONLEAVE__
    return Sample;
}

G4double G4FPYSamplingOps::
G4SampleGaussian( G4double Mean,
                  G4double StdDev,
                  G4FFGEnumerations::GaussianRange Range )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    if(Range == G4FFGEnumerations::ALL)
    {
        // Call the overloaded function
        G4double Sample = G4SampleGaussian(Mean, StdDev);
        
G4FFG_SAMPLING_FUNCTIONLEAVE__
        return Sample;
    }

    // Determine if the parameters have changed
    G4bool ParametersChanged = (Mean_ != Mean ||
                                StdDev_ != StdDev);
    if(ParametersChanged == TRUE)
    {
        if(Mean <= 0)
        {
            std::ostringstream Temp;
            Temp << "Mean value of " << Mean << " out of range";
            G4Exception("G4FPYGaussianOps::G4SampleIntegerGaussian()",
            Temp.str().c_str(),
            JustWarning,
            "A value of '0' will be used instead.");

G4FFG_SAMPLING_FUNCTIONLEAVE__
            return 0;
        }

        // Set the new values if the parameters have changed and then perform
        // the shift
        Mean_ = Mean;
        StdDev_ = StdDev;

        ShiftParameters(G4FFGEnumerations::DOUBLE);
    }

    // Sample the Gaussian distribution
    G4double Rand;
    do
    {
        Rand = SampleGaussian();
    } while(Rand < 0); // Loop checking, 11.05.2015, T. Koi

G4FFG_SAMPLING_FUNCTIONLEAVE__
    return Rand;
}

G4int G4FPYSamplingOps::
G4SampleIntegerGaussian( G4double Mean,
                         G4double StdDev )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    // Determine if the parameters have changed
    G4bool ParametersChanged = (Mean_ != Mean || StdDev_ != StdDev);
    if(ParametersChanged == TRUE)
    {
        // Set the new values if the parameters have changed
        NextGaussianIsStoredInMemory_ = FALSE;

        Mean_ = Mean;
        StdDev_ = StdDev;
    }

    // Return the sample integer value
    G4int Sample = (G4int)(std::floor(SampleGaussian()));

G4FFG_SAMPLING_FUNCTIONLEAVE__
    return Sample;
}

G4int G4FPYSamplingOps::
G4SampleIntegerGaussian( G4double Mean,
                         G4double StdDev,
                         G4FFGEnumerations::GaussianRange Range )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    if(Range == G4FFGEnumerations::ALL)
    {
        // Call the overloaded function
        G4int Sample = G4SampleIntegerGaussian(Mean, StdDev);

G4FFG_SAMPLING_FUNCTIONLEAVE__
        return Sample;
    } else
    {
        // ParameterShift() locks up if the mean is less than 1.
        std::ostringstream Temp;
        if(Mean < 1)
        {
        //    Temp << "Mean value of " << Mean << " out of range";
        //    G4Exception("G4FPYGaussianOps::G4SampleIntegerGaussian()",
        //                Temp.str().c_str(),
        //                JustWarning,
        //                "A value of '0' will be used instead.");

        //    return 0;
        }
        
        if(Mean / StdDev < 2)
        {
            //Temp << "Non-ideal conditions:\tMean:" << Mean << "\tStdDev: "
            //        << StdDev;
            //G4Exception("G4FPYGaussianOps::G4SampleIntegerGaussian()",
            //            Temp.str().c_str(),
            //            JustWarning,
            //            "Results not guaranteed: Consider tightening the standard deviation");
        }

        // Determine if the parameters have changed
        G4bool ParametersChanged = (Mean_ != Mean ||
                                    StdDev_ != StdDev);
        if(ParametersChanged == TRUE)
        {
            // Set the new values if the parameters have changed and then perform
            // the shift
            Mean_ = Mean;
            StdDev_ = StdDev;

            ShiftParameters(G4FFGEnumerations::INT);
        }

        G4int RandInt;
        // Sample the Gaussian distribution - only non-negative values are
        // accepted
        do
        {
            RandInt = (G4int)floor(SampleGaussian());
        } while (RandInt < 0); // Loop checking, 11.05.2015, T. Koi

G4FFG_SAMPLING_FUNCTIONLEAVE__
        return RandInt;
    }
}

G4double G4FPYSamplingOps::
G4SampleUniform( void )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    G4double Sample = RandomEngine_->flat();
    
G4FFG_SAMPLING_FUNCTIONLEAVE__
    return Sample;
}

G4double G4FPYSamplingOps::
G4SampleUniform( G4double Lower,
                 G4double Upper )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    // Calculate the difference
    G4double Difference = Upper - Lower;

    // Scale appropriately and return the value
    G4double Sample = G4SampleUniform() * Difference + Lower;

G4FFG_SAMPLING_FUNCTIONLEAVE__
    return Sample;
}

G4double G4FPYSamplingOps::
G4SampleWatt( G4int WhatIsotope,
              G4FFGEnumerations::FissionCause WhatCause,
              G4double WhatEnergy )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    // Determine if the parameters are different
    //TK modified 131108 
    //if(WattConstants_->Product != WhatIsotope
    if(WattConstants_->Product != WhatIsotope/10
            || WattConstants_->Cause != WhatCause
            || WattConstants_->Energy!= WhatEnergy )
    {
        // If the parameters are different or have not yet been defined then we
        // need to re-evaluate the constants
        //TK modified 131108 
        //WattConstants_->Product = WhatIsotope;
        WattConstants_->Product = WhatIsotope/10;
        WattConstants_->Cause = WhatCause;
        WattConstants_->Energy = WhatEnergy;

        EvaluateWattConstants();
    }

    G4double X = -G4Log(G4SampleUniform());
    G4double Y = -G4Log(G4SampleUniform());
    G4int icounter=0;
    G4int icounter_max=1024;
    while(G4Pow::GetInstance()->powN(Y - WattConstants_->M*(X + 1), 2)
            > WattConstants_->B * WattConstants_->L * X) // Loop checking, 11.05.2015, T. Koi
    {
        icounter++;
        if ( icounter > icounter_max ) {
           G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
           break;
        }
        X = -G4Log(G4SampleUniform());
        Y = -G4Log(G4SampleUniform());
    }

G4FFG_SAMPLING_FUNCTIONLEAVE__
    return WattConstants_->L * X;
}

void G4FPYSamplingOps::
G4SetVerbosity(G4int Verbosity)
{
G4FFG_SAMPLING_FUNCTIONENTER__

    Verbosity_ = Verbosity;
    
    ShiftedGaussianValues_->G4SetVerbosity(Verbosity_);

G4FFG_SAMPLING_FUNCTIONLEAVE__
}

G4bool G4FPYSamplingOps::
CheckAndSetParameters( void )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    G4double ShiftedMean = ShiftedGaussianValues_->G4FindShiftedMean(Mean_, StdDev_);
    if(ShiftedMean == 0)
    {
G4FFG_SAMPLING_FUNCTIONLEAVE__
        return FALSE;
    }

    Mean_ = ShiftedMean;

G4FFG_SAMPLING_FUNCTIONLEAVE__
    return TRUE;
}

void G4FPYSamplingOps::
EvaluateWattConstants( void )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    G4double A, K;
    A = K = 0;
    // Use the default values if IsotopeIndex is not reset
    G4int IsotopeIndex = 0;

    if(WattConstants_->Cause == G4FFGEnumerations::SPONTANEOUS)
    {
        // Determine if the isotope requested exists in the lookup tables
        for(G4int i = 0; SpontaneousWattIsotopesIndex[i] != -1; i++)
        {
            if(SpontaneousWattIsotopesIndex[i] ==
                    WattConstants_->Product)
            {
                IsotopeIndex = i;

                break;
            }
        }

        // Get A and B
        A = SpontaneousWattConstants[IsotopeIndex][0];
        WattConstants_->B = SpontaneousWattConstants[IsotopeIndex][1];
    } else if (WattConstants_->Cause == G4FFGEnumerations::NEUTRON_INDUCED)
    {
        // Determine if the isotope requested exists in the lookup tables
        for(G4int i = 0; NeutronInducedWattIsotopesIndex[i] != -1; i++)
        {
            if(NeutronInducedWattIsotopesIndex[i] == WattConstants_->Product)
            {
                IsotopeIndex = i;
                break;
            }
        }

        // Determine the Watt fission spectrum constants based on the energy of
        // the incident neutron
        if(WattConstants_->Energy == G4FFGDefaultValues::ThermalNeutronEnergy)
        {
            A = NeutronInducedWattConstants[IsotopeIndex][0][0];
            WattConstants_->B = NeutronInducedWattConstants[IsotopeIndex][0][1];
        } else if (WattConstants_->Energy > 14.0 * CLHEP::MeV)
        {
            G4Exception("G4FPYSamplingOps::G4SampleWatt()",
                        "Incident neutron energy above 14 MeV requested.",
                        JustWarning,
                        "Using Watt fission constants for 14 Mev.");

            A = NeutronInducedWattConstants[IsotopeIndex][2][0];
            WattConstants_->B = NeutronInducedWattConstants[IsotopeIndex][2][1];
        } else
        {
            G4int EnergyIndex = 0;
            G4double EnergyDifference = 0;
            G4double RangeDifference, ConstantDifference;

            for(G4int i = 1; IncidentEnergyBins[i] != -1; i++)
            {
                if(WattConstants_->Energy <= IncidentEnergyBins[i])
                {
                    EnergyIndex = i;
                    EnergyDifference = IncidentEnergyBins[EnergyIndex] - WattConstants_->Energy;
                    if(EnergyDifference != 0)
                    {
                        std::ostringstream Temp;
                        Temp << "Incident neutron energy of ";
                        Temp << WattConstants_->Energy << " MeV is not ";
                        Temp << "explicitly listed in the data tables";
//                        G4Exception("G4FPYSamplingOps::G4SampleWatt()",
//                                    Temp.str().c_str(),
//                                    JustWarning,
//                                    "Using linear interpolation.");
                    }
                    break;
                }
            }

            RangeDifference = IncidentEnergyBins[EnergyIndex] - IncidentEnergyBins[EnergyIndex - 1];

            // Interpolate the value for 'a'
            ConstantDifference =
                NeutronInducedWattConstants[IsotopeIndex][EnergyIndex][0] -
                NeutronInducedWattConstants[IsotopeIndex]
                                           [EnergyIndex - 1][0];
            A = (EnergyDifference / RangeDifference) * ConstantDifference +
                NeutronInducedWattConstants[IsotopeIndex]
                                           [EnergyIndex - 1][0];

            // Interpolate the value for 'b'
            ConstantDifference =
                NeutronInducedWattConstants[IsotopeIndex][EnergyIndex][1] -
                NeutronInducedWattConstants[IsotopeIndex]
                                           [EnergyIndex - 1][1];
            WattConstants_->B =
                (EnergyDifference / RangeDifference) * ConstantDifference +
                NeutronInducedWattConstants[IsotopeIndex]
                                           [EnergyIndex - 1][1];
        }
    } else
    {
        // Produce an error since an unsupported fission type was requested and
        // abort the current run
        G4String Temp = "Watt fission spectra data not available for ";
        if(WattConstants_->Cause == G4FFGEnumerations::PROTON_INDUCED)
        {
            Temp += "proton induced fission.";
        } else if(WattConstants_->Cause == G4FFGEnumerations::GAMMA_INDUCED)
        {
            Temp += "gamma induced fission.";
        } else
        {
            Temp += "!Warning! unknown cause.";
        }
        G4Exception("G4FPYSamplingOps::G4SampleWatt()",
                    Temp,
                    RunMustBeAborted,
                    "Fission events will not be sampled in this run.");
    }

    // Calculate the constants
    K = 1 + (WattConstants_->B / (8.0 * A));
    WattConstants_->L = (K + G4Pow::GetInstance()->powA((K * K - 1), 0.5)) / A;
    WattConstants_->M = A * WattConstants_->L - 1;

G4FFG_SAMPLING_FUNCTIONLEAVE__
}

G4double G4FPYSamplingOps::
SampleGaussian( void )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    if(NextGaussianIsStoredInMemory_ == TRUE)
    {
        NextGaussianIsStoredInMemory_ = FALSE;

G4FFG_SAMPLING_FUNCTIONLEAVE__
        return GaussianTwo_;
    }

    // Define the variables to be used
    G4double Radius;
    G4double MappingFactor;

    // Sample from the unit circle (21.4% rejection probability)
    do
    {
        GaussianOne_ = 2.0 * G4SampleUniform() - 1.0;
        GaussianTwo_ = 2.0 * G4SampleUniform() - 1.0;
        Radius = GaussianOne_*GaussianOne_ + GaussianTwo_*GaussianTwo_;
    } while (Radius > 1.0); // Loop checking, 11.05.2015, T. Koi

    // Translate the values to Gaussian space
    MappingFactor = std::sqrt(-2.0*G4Log(Radius)/Radius) * StdDev_;
    GaussianOne_ = Mean_ + GaussianOne_*MappingFactor;
    GaussianTwo_ = Mean_ + GaussianTwo_*MappingFactor;

    // Set the flag that a value is now stored in memory
    NextGaussianIsStoredInMemory_ = TRUE;

G4FFG_SAMPLING_FUNCTIONLEAVE__
    return GaussianOne_;
}

void G4FPYSamplingOps::
ShiftParameters( G4FFGEnumerations::GaussianReturnType Type )
{
G4FFG_SAMPLING_FUNCTIONENTER__

    // Set the flag that any second value stored is no longer valid
    NextGaussianIsStoredInMemory_ = FALSE;

    // Check if the requested parameters have already been calculated
    if(CheckAndSetParameters() == TRUE)
    {
G4FFG_SAMPLING_FUNCTIONLEAVE__
        return;
    }

    // If the requested type is INT, then perform an iterative refinement of the
    // mean that is going to be sampled
    if(Type == G4FFGEnumerations::INT)
    {
        // Return if the mean is greater than 7 standard deviations away from 0
        // since there is essentially 0 probability that a sampled number will
        // be less than 0
        if(Mean_ > 7 * StdDev_)
        {
G4FFG_SAMPLING_FUNCTIONLEAVE__
            return;
        }
        // Variables that contain the area and weighted area information for
        // calculating the statistical mean of the Gaussian distribution when
        // converted to a step function
        G4double ErfContainer, AdjustedErfContainer, Container;

        // Variables that store lower and upper bounds information
        G4double LowErf, HighErf;

        // Values for determining the convergence of the solution
        G4double AdjMean = Mean_;
        G4double Delta = 1.0;
        G4bool HalfDelta = false;
        G4bool ToleranceCheck = false;


        // Normalization constant for use with erf()
        const G4double Normalization = StdDev_ * std::sqrt(2.0);

        // Determine the upper limit of the estimates
        const G4int UpperLimit = (G4int) std::ceil(Mean_ + 9 * StdDev_);

        // Calculate the integral of the Gaussian distribution

        G4int icounter=0;
        G4int icounter_max=1024;
        do
        {
            icounter++;
            if ( icounter > icounter_max ) {
	       G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
               break;
            }
            // Set the variables
            ErfContainer = 0;
            AdjustedErfContainer = 0;

            // Calculate the area and weighted area
            for(G4int i = 0; i <= UpperLimit; i++)
            {
                // Determine the lower and upper bounds
                LowErf = ((AdjMean - i) / Normalization);
                HighErf = ((AdjMean - (i + 1.0)) / Normalization);

                // Correct the bounds for how they lie on the x-axis with
                // respect to the mean
                if(LowErf <= 0)
                {
                    LowErf *= -1;
                    HighErf *= -1;
                    //Container = (erf(HighErf) - erf(LowErf))/2.0;
                    #if defined WIN32-VC
			Container = (CLHEP::HepStat::erf(HighErf) - CLHEP::HepStat::erf(LowErf))/2.0;
                    #else
                        Container = (erf(HighErf) - erf(LowErf))/2.0;
                    #endif              
                } else if (HighErf < 0)
                {
                    HighErf *= -1;

                    //Container = (erf(HighErf) + erf(LowErf))/2.0;
                    #if defined WIN32-VC
		        Container = (CLHEP::HepStat::erf(HighErf) + CLHEP::HepStat::erf(LowErf))/2.0;
                    #else
                        Container = (erf(HighErf) + erf(LowErf))/2.0;
                    #endif
                } else
                {
                    //Container = (erf(LowErf) - erf(HighErf))/2.0;
                    #if defined WIN32-VC
			Container = (CLHEP::HepStat::erf(LowErf) - CLHEP::HepStat::erf(HighErf))/2.0;
                    #else
                        Container = (erf(LowErf) - erf(HighErf))/2.0;
                    #endif
                }

                #if defined WIN32-VC
                //TK modified to avoid problem caused by QNAN
		    if ( Container != Container) Container = 0;
                #endif

                // Add up the weighted area
                ErfContainer += Container;
                AdjustedErfContainer += Container * i;
            }

            // Calculate the statistical mean
            Container = AdjustedErfContainer / ErfContainer;

            // Is it close enough to what we want?
            ToleranceCheck = (std::fabs(Mean_ - Container) < Tolerance_);
            if(ToleranceCheck == TRUE)
            {
                break;
            }

            // Determine the step size
            if(HalfDelta == TRUE)
            {
                Delta /= 2;
            }

            // Step in the appropriate direction
            if(Container > Mean_)
            {
                AdjMean -= Delta;
            } else
            {
                HalfDelta = TRUE;
                AdjMean += Delta;
            }
        } while(TRUE); // Loop checking, 11.05.2015, T. Koi

        ShiftedGaussianValues_->G4InsertShiftedMean(AdjMean, Mean_, StdDev_);
        Mean_ = AdjMean;
    } else if(Mean_ / 7 < StdDev_)
    {
        // If the requested type is double, then just re-define the standard
        // deviation appropriately - chances are approximately 2.56E-12 that
        // the value will be negative using this sampling scheme
        StdDev_ = Mean_ / 7;
    }
    
G4FFG_SAMPLING_FUNCTIONLEAVE__
}

G4FPYSamplingOps::~G4FPYSamplingOps( void )
{
G4FFG_FUNCTIONENTER__

    delete ShiftedGaussianValues_;
    delete WattConstants_;
    
G4FFG_FUNCTIONLEAVE__
}

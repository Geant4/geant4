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
 * File:   G4FissionProductYieldDist.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on May 11, 2011, 12:04 PM
 */

/* * * * * * * * * * * * * * * *   References   * * * * * * * * * * * * * * * *
 *                                                                            *
 *  1.  "Systematics of fission fragment total kinetic energy release",       *
 *      V. E. Viola, K. Kwiatkowski, and M. Walker, Physical Review C, 31.4,  *
 *      April 1985                                                            *
 *  2.  "Reactor Handbook", United States Atomic Energy Commission,           *
 *      III.A:Physics, 1962                                                   *
 *  3.  "Properties of the Alpha Particles Emitted in the Spontaneous Fission *
 *      of Cf252", Z. Fraenkel and S. G. Thompson, Physical Review Letters,   *
 *      13.14, October 1964                                                   *
 *                                                                            *
 * * * * * * * * * * * * * * * *   References   * * * * * * * * * * * * * * * */


//#include <ios>
//#include <iostream>

#include "G4Alpha.hh"
#include "G4Gamma.hh"
#include "G4Ions.hh"
#include "G4Neutron.hh"
//#include "G4NeutronHPNames.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "G4Exp.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadFinalState.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ReactionProduct.hh"

#include "G4ArrayOps.hh"
#include "G4ENDFTapeRead.hh"
#include "G4ENDFYieldDataContainer.hh"
#include "G4FFGDebuggingMacros.hh"
#include "G4FFGDefaultValues.hh"
#include "G4FFGEnumerations.hh"
#include "G4FPYNubarValues.hh"
#include "G4FPYSamplingOps.hh"
#include "G4FPYTreeStructures.hh"
#include "G4FissionProductYieldDist.hh"
#include "G4Pow.hh"

using CLHEP::pi;

#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
G4Mutex G4FissionProductYieldDist::fissprodMutex = G4MUTEX_INITIALIZER;
#endif

G4FissionProductYieldDist::
G4FissionProductYieldDist( G4int WhichIsotope,
                           G4FFGEnumerations::MetaState WhichMetaState,
                           G4FFGEnumerations::FissionCause WhichCause,
                           G4FFGEnumerations::YieldType WhichYieldType,
                           std::istringstream& dataStream )
:   Isotope_(WhichIsotope),
    MetaState_(WhichMetaState),
    Cause_(WhichCause),
    YieldType_(WhichYieldType),
    Verbosity_(G4FFGDefaultValues::Verbosity)
{
    G4FFG_FUNCTIONENTER__

    try
    {
        // Initialize the class
        Initialize(dataStream);
    } catch (std::exception& e)
    {
        G4FFG_FUNCTIONLEAVE__
        throw e;
    }

    G4FFG_FUNCTIONLEAVE__
}

G4FissionProductYieldDist::
G4FissionProductYieldDist(G4int WhichIsotope,
                          G4FFGEnumerations::MetaState WhichMetaState,
                          G4FFGEnumerations::FissionCause WhichCause,
                          G4FFGEnumerations::YieldType WhichYieldType,
                          G4int Verbosity,
                          std::istringstream& dataStream)
 : Isotope_(WhichIsotope), MetaState_(WhichMetaState), Cause_(WhichCause),
   YieldType_(WhichYieldType), Verbosity_(Verbosity)
{
  G4FFG_FUNCTIONENTER__

  try
  {
    // Initialize the class
    Initialize(dataStream);
  } catch (std::exception& e)
  {
    G4FFG_FUNCTIONLEAVE__
    throw e;
  }

  G4FFG_FUNCTIONLEAVE__
}


void G4FissionProductYieldDist::Initialize(std::istringstream& dataStream)
{
G4FFG_FUNCTIONENTER__

  IncidentEnergy_ = 0.0;
  TernaryProbability_ = 0;
  AlphaProduction_ = 0;
  SetNubar();

    // Set miscellaneous variables
    AlphaDefinition_ = static_cast<G4Ions*>(G4Alpha::Definition());
    NeutronDefinition_ = static_cast<G4Ions*>(G4Neutron::Definition());
    GammaDefinition_ = G4Gamma::Definition();
    SmallestZ_ = SmallestA_ = LargestZ_ = LargestA_ = NULL;

    // Construct G4NeutronHPNames: provides access to the element names
    ElementNames_ = new G4ParticleHPNames;
    // Get the pointer to G4ParticleTable: stores all G4Ions
    IonTable_ = G4IonTable::GetIonTable();
    // Construct the pointer to the random engine
    // TODO Make G4FPSamplingOps a singleton so that only one instance is used across all classes
    RandomEngine_ = new G4FPYSamplingOps;

    try
    {
        // Read in and sort the probability data
        ENDFData_ = new G4ENDFTapeRead(dataStream,
                                       YieldType_,
                                       Cause_,
                                       Verbosity_);
//        ENDFData_ = new G4ENDFTapeRead(MakeDirectoryName(),
//                                       MakeFileName(Isotope_, MetaState_),
//                                       YieldType_,
//                                       Cause_);
        YieldEnergyGroups_ = ENDFData_->G4GetNumberOfEnergyGroups();
        DataTotal_ = new G4double[YieldEnergyGroups_];
        MaintainNormalizedData_ = new G4double[YieldEnergyGroups_];
        YieldEnergies_ = new G4double[YieldEnergyGroups_];
        G4ArrayOps::Copy(YieldEnergyGroups_, YieldEnergies_, ENDFData_->G4GetEnergyGroupValues());
        MakeTrees();
        ReadProbabilities();
    } catch (std::exception& e)
    {
        delete ElementNames_;
        delete RandomEngine_;

        G4FFG_FUNCTIONLEAVE__
        throw e;
    }

G4FFG_FUNCTIONLEAVE__
}

G4DynamicParticleVector* G4FissionProductYieldDist::G4GetFission(void)
{
G4FFG_FUNCTIONENTER__

#ifdef G4MULTITHREADED
  G4AutoLock lk(&G4FissionProductYieldDist::fissprodMutex);
#endif

  // Check to see if the user has set the alpha production to a somewhat
  // reasonable level
  CheckAlphaSanity();

  // Generate the new G4DynamicParticle pointers to identify key locations in
  // the G4DynamicParticle chain that will be passed to the G4FissionEvent
  G4ReactionProduct* FirstDaughter = NULL;
  G4ReactionProduct* SecondDaughter = NULL;
  std::vector< G4ReactionProduct* >* Alphas = new std::vector< G4ReactionProduct* >;
  std::vector< G4ReactionProduct* >* Neutrons = new std::vector< G4ReactionProduct* >;
  std::vector< G4ReactionProduct* >* Gammas = new std::vector< G4ReactionProduct* >;

  // Generate all the nucleonic fission products
  // How many nucleons do we have to work with?
  //TK modified 131108 
  const G4int ParentA = (Isotope_/10) % 1000;
  const G4int ParentZ = ((Isotope_/10) - ParentA) / 1000;
  RemainingA_ = ParentA;
  RemainingZ_ = ParentZ;

  // Don't forget the extra nucleons depending on the fission cause
  switch(Cause_)
  {
    case G4FFGEnumerations::NEUTRON_INDUCED:
       ++RemainingA_;
       break;

    case G4FFGEnumerations::PROTON_INDUCED:
       ++RemainingZ_;
       break;

    case G4FFGEnumerations::GAMMA_INDUCED:
    case G4FFGEnumerations::SPONTANEOUS:
    default:
    // Nothing to do here
    break;
  }

  // Ternary fission can be set by the user. Thus, it is necessary to
  // sample the alpha particle first and the first daughter product
  // second. See the discussion in
  // G4FissionProductYieldDist::G4GetFissionProduct() for more information
  // as to why the fission events are sampled this way.
  GenerateAlphas(Alphas);

  // Generate the first daughter product
  FirstDaughter = new G4ReactionProduct(GetFissionProduct());
  RemainingA_ -= FirstDaughter->GetDefinition()->GetAtomicMass();
  RemainingZ_ -= FirstDaughter->GetDefinition()->GetAtomicNumber();
  if (Verbosity_ & G4FFGEnumerations::DAUGHTER_INFO) {
    G4FFG_SPACING__
    G4FFG_LOCATION__
            
    G4cout << " -- First daughter product sampled" << G4endl;
    G4FFG_SPACING__
    G4cout << "  Name:       " << FirstDaughter->GetDefinition()->GetParticleName() << G4endl;
    G4FFG_SPACING__
    G4cout << "  Z:          " << FirstDaughter->GetDefinition()->GetAtomicNumber() << G4endl;
    G4FFG_SPACING__
    G4cout << "  A:          " << FirstDaughter->GetDefinition()->GetAtomicMass() << G4endl;
    G4FFG_SPACING__
    G4cout << "  Meta State: " << (FirstDaughter->GetDefinition()->GetPDGEncoding() % 10) << G4endl;
  }

  GenerateNeutrons(Neutrons);

  // Now that all the nucleonic particles have been generated, we can
  // calculate the composition of the second daughter product.
  G4int NewIsotope = RemainingZ_ * 1000 + RemainingA_;
  SecondDaughter = new G4ReactionProduct(GetParticleDefinition(NewIsotope, G4FFGEnumerations::GROUND_STATE));
  if (Verbosity_ & G4FFGEnumerations::DAUGHTER_INFO) {
    G4FFG_SPACING__
    G4FFG_LOCATION__
            
    G4cout << " -- Second daughter product sampled" << G4endl;
    G4FFG_SPACING__
    G4cout << "  Name:       " << SecondDaughter->GetDefinition()->GetParticleName() << G4endl;
    G4FFG_SPACING__
    G4cout << "  Z:          " << SecondDaughter->GetDefinition()->GetAtomicNumber() << G4endl;
    G4FFG_SPACING__
    G4cout << "  A:          " << SecondDaughter->GetDefinition()->GetAtomicMass() << G4endl;
    G4FFG_SPACING__
    G4cout << "  Meta State: " << (SecondDaughter->GetDefinition()->GetPDGEncoding() % 10) << G4endl;
  }

    // Calculate how much kinetic energy will be available
        // 195 to 205 MeV are available in a fission reaction, but about 20 MeV
        // are from delayed sources. We are concerned only with prompt sources,
        // so sample a Gaussian distribution about 20 MeV and subtract the
        // result from the total available energy. Also, the energy of fission
        // neutrinos is neglected. Fission neutrinos would add ~11 MeV
        // additional energy to the fission. (Ref 2)
        // Finally, add in the kinetic energy contribution of the fission
        // inducing particle, if any.
        const G4double TotalKE = RandomEngine_->G4SampleUniform(195.0, 205.0) * MeV
                                 - RandomEngine_->G4SampleGaussian(20.0, 3.0) * MeV
                                 + IncidentEnergy_;
        RemainingEnergy_ = TotalKE;

    // Calculate the energies of the alpha particles and neutrons
        // Once again, since the alpha particles are user defined, we must
        // sample their kinetic energy first. SampleAlphaEnergies() returns the
        // amount of energy consumed by the alpha particles, so remove the total
        // energy alloted to the alpha particles from the available energy
        SampleAlphaEnergies(Alphas);

        // Second, the neutrons are sampled from the Watt fission spectrum.
        SampleNeutronEnergies(Neutrons);

    // Calculate the total energy available to the daughter products
        // A Gaussian distribution about the average calculated energy with
        // a standard deviation of 1.5 MeV (Ref. 2) is used. Since the energy
        // distribution is dependant on the alpha particle generation and the
        // Watt fission sampling for neutrons, we only have the left-over energy
        // to work with for the fission daughter products.
        G4double FragmentsKE=0.;
        G4int icounter=0;
        G4int icounter_max=1024;
        do
        {
           icounter++;
           if ( icounter > icounter_max ) {
	      G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
              break;
           }
            FragmentsKE = RandomEngine_->G4SampleGaussian(RemainingEnergy_, 1.5 *MeV);
        } while(FragmentsKE > RemainingEnergy_); // Loop checking, 11.05.2015, T. Koi

        // Make sure that we don't produce any sub-gamma photons
        if((RemainingEnergy_ - FragmentsKE) / (100 * keV) < 1.0)
        {
            FragmentsKE = RemainingEnergy_;
        }

        // This energy has now been allotted to the fission fragments.
        // Subtract FragmentsKE from the total available fission energy.
        RemainingEnergy_ -= FragmentsKE;

    // Sample the energies of the gamma rays
        // Assign the remainder, if any, of the energy to the gamma rays
        SampleGammaEnergies(Gammas);

    // Prepare to balance the momenta of the system
        // We will need these for sampling the angles and balancing the momenta
        // of all the fission products
        G4double NumeratorSqrt, NumeratorOther, Denominator;
        G4double Theta, Phi, Magnitude;
        G4ThreeVector Direction;
        G4ParticleMomentum ResultantVector(0, 0, 0);

        if(Alphas->size() > 0)
        {
        // Sample the angles of the alpha particles and neutrons, then calculate
        // the total moment contribution to the system
            // The average angle of the alpha particles with respect to the
            // light fragment is dependent on the ratio of the kinetic energies.
            // This equation was determined by performing a fit on data from
            // Ref. 3 using the website:
            // http://soft.arquimedex.com/linear_regression.php
            //
            // RatioOfKE    Angle (rad)     Angle (degrees)
            // 1.05         1.257           72
            // 1.155        1.361           78
            // 1.28         1.414           81
            // 1.5          1.518           87
            // 1.75         1.606           92
            // 1.9          1.623           93
            // 2.2          1.728           99
            // This equation generates the angle in radians. If the RatioOfKE is
            // greater than 2.25 the angle defaults to 1.3963 rad (100 degrees)
            G4double MassRatio = FirstDaughter->GetDefinition()->GetPDGMass() / SecondDaughter->GetDefinition()->GetPDGMass();
            
            // Invert the mass ratio if the first daughter product is the lighter fragment
            if(MassRatio < 1)
            {
                MassRatio = 1 / MassRatio;
            }
            
            // The empirical equation is valid for mass ratios up to 2.75
            if(MassRatio > 2.75)
            {
                MassRatio = 2.75;
            }
            const G4double MeanAlphaAngle = 0.3644 * MassRatio * MassRatio * MassRatio
                             	 	 	 	- 1.9766 * MassRatio * MassRatio
                             	 	 	 	+ 3.8207 * MassRatio
                             	 	 	 	- 0.9917;

            // Sample the directions of the alpha particles with respect to the
            // light fragment. For the moment we will assume that the light
            // fragment is traveling along the z-axis in the positive direction.
            const G4double MeanAlphaAngleStdDev = 0.0523598776;
            G4double PlusMinus;

            for(unsigned int i = 0; i < Alphas->size(); ++i)
            {
                PlusMinus = std::acos(RandomEngine_->G4SampleGaussian(0, MeanAlphaAngleStdDev)) - (pi / 2);
                Theta = MeanAlphaAngle + PlusMinus;
                if(Theta < 0)
                {
                    Theta = 0.0 - Theta;
                } else if(Theta > pi)
                {
                    Theta = (2 * pi - Theta);
                }
                Phi = RandomEngine_->G4SampleUniform(-pi, pi);

                Direction.setRThetaPhi(1.0, Theta, Phi);
                Magnitude = std::sqrt(2 * Alphas->at(i)->GetKineticEnergy() * Alphas->at(i)->GetDefinition()->GetPDGMass());
                Alphas->at(i)->SetMomentum(Direction * Magnitude);
                ResultantVector += Alphas->at(i)->GetMomentum();
            }
        }

        // Sample the directions of the neutrons.
        if(Neutrons->size() != 0)
        {
        	for(unsigned int i = 0; i < Neutrons->size(); ++i)
        	{
        	    Theta = std::acos(RandomEngine_->G4SampleUniform(-1, 1));
                Phi = RandomEngine_->G4SampleUniform(-pi, pi);

                Direction.setRThetaPhi(1.0, Theta, Phi);
                Magnitude = std::sqrt(2 * Neutrons->at(i)->GetKineticEnergy() * Neutrons->at(i)->GetDefinition()->GetPDGMass());
                Neutrons->at(i)->SetMomentum(Direction * Magnitude);
                ResultantVector += Neutrons->at(i)->GetMomentum();
            }
        }

        // Sample the directions of the gamma rays
        if(Gammas->size() != 0)
		{
        	for(unsigned int i = 0; i < Gammas->size(); ++i)
            {
                Theta = std::acos(RandomEngine_->G4SampleUniform(-1, 1));
                Phi = RandomEngine_->G4SampleUniform(-pi, pi);

                Direction.setRThetaPhi(1.0, Theta, Phi);
                Magnitude = Gammas->at(i)->GetKineticEnergy() / CLHEP::c_light;
                Gammas->at(i)->SetMomentum(Direction * Magnitude);
                ResultantVector += Gammas->at(i)->GetMomentum();
            }
        }

    // Calculate the momenta of the two daughter products
        G4ReactionProduct* LightFragment;
        G4ReactionProduct* HeavyFragment;
        G4ThreeVector LightFragmentDirection;
        G4ThreeVector HeavyFragmentDirection;
        G4double ResultantX, ResultantY, ResultantZ;
        ResultantX = ResultantVector.getX();
        ResultantY = ResultantVector.getY();
        ResultantZ = ResultantVector.getZ();

        if(FirstDaughter->GetDefinition()->GetPDGMass() < SecondDaughter->GetDefinition()->GetPDGMass())
        {
            LightFragment = FirstDaughter;
            HeavyFragment = SecondDaughter;
        } else
        {
            LightFragment = SecondDaughter;
            HeavyFragment = FirstDaughter;
        }
        const G4double LightFragmentMass = LightFragment->GetDefinition()->GetPDGMass();
        const G4double HeavyFragmentMass = HeavyFragment->GetDefinition()->GetPDGMass();

        LightFragmentDirection.setRThetaPhi(1.0, 0, 0);

        // Fit the momenta of the daughter products to the resultant vector of
        // the remaining fission products. This will be done in the Cartesian
        // coordinate system, not spherical. This is done using the following
        // table listing the system momenta and the corresponding equations:
        //              X               Y               Z
        //
        //      A       0               0               P
        //
        //      B       -R_x            -R_y            -P - R_z
        //
        //      R       R_x             R_y             R_z
        //
        // v = sqrt(2*m*k)  ->  k = v^2/(2*m)
        // tk = k_A + k_B
        // k_L = P^2/(2*m_L)
        // k_H = ((-R_x)^2 + (-R_y)^2 + (-P - R_z)^2)/(2*m_H)
        // where:
        // P: momentum of the light daughter product
        // R: the remaining fission products' resultant vector
        // v: momentum
        // m: mass
        // k: kinetic energy
        // tk: total kinetic energy available to the daughter products
        //
        // Below is the solved form for P, with the solution generated using
        // the WolframAlpha website:
        // http://www.wolframalpha.com/input/?i=
        // solve+((-x)^2+%2B+(-y)^2+%2B+(-P-z)^2)%2F(2*B)+%2B+L^2%2F(2*A)+%3D+k+
        // for+P
        //
        //
        // nsqrt = sqrt(m_L*(m_L*(2*m_H*tk - R_x^2 - R_y^2) +
        //                   m_H*(2*m_H*tk - R_x^2 - R_y^2 - R_z^2))
        NumeratorSqrt = std::sqrt(LightFragmentMass *
                            (LightFragmentMass * (2 * HeavyFragmentMass * FragmentsKE
                                                  - ResultantX * ResultantX
                                                  - ResultantY * ResultantY)
                             + HeavyFragmentMass * (2 * HeavyFragmentMass * FragmentsKE
                                                    - ResultantX * ResultantX
                                                    - ResultantY * ResultantY
                                                    - ResultantZ - ResultantZ)));

        // nother = m_L*R_z
        NumeratorOther = LightFragmentMass * ResultantZ;

        // denom = m_L + m_H
        Denominator = LightFragmentMass + HeavyFragmentMass;

        // P = (nsqrt + nother) / denom
        const G4double LightFragmentMomentum = (NumeratorSqrt + NumeratorOther) / Denominator;
        const G4double HeavyFragmentMomentum = std::sqrt(ResultantX * ResultantX
                                                    + ResultantY * ResultantY
                                                    + G4Pow::GetInstance()->powN(LightFragmentMomentum + ResultantZ, 2));

        // Finally! We now have everything we need for the daughter products
        LightFragment->SetMomentum(LightFragmentDirection * LightFragmentMomentum);
        HeavyFragmentDirection.setX(-ResultantX);
        HeavyFragmentDirection.setY(-ResultantY);
        HeavyFragmentDirection.setZ((-LightFragmentMomentum - ResultantZ));
        // Don't forget to normalize the vector to the unit sphere
        HeavyFragmentDirection.setR(1.0);
        HeavyFragment->SetMomentum(HeavyFragmentDirection * HeavyFragmentMomentum);

        if(Verbosity_ & (G4FFGEnumerations::DAUGHTER_INFO | G4FFGEnumerations::MOMENTUM_INFO))
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            G4cout << " -- Daugher product momenta finalized" << G4endl;
            G4FFG_SPACING__
        }

    // Load all the particles into a contiguous set
        //TK modifed 131108
        //G4DynamicParticleVector* FissionProducts = new G4DynamicParticleVector(2 + Alphas->size() + Neutrons->size() + Gammas->size());
        G4DynamicParticleVector* FissionProducts = new G4DynamicParticleVector();
        // Load the fission fragments
        FissionProducts->push_back(MakeG4DynamicParticle(LightFragment));
        FissionProducts->push_back(MakeG4DynamicParticle(HeavyFragment));
        // Load the neutrons
        for(unsigned int i = 0; i < Neutrons->size(); i++)
        {
            FissionProducts->push_back(MakeG4DynamicParticle(Neutrons->at(i)));
        }
        // Load the gammas
        for(unsigned int i = 0; i < Gammas->size(); i++)
        {
            FissionProducts->push_back(MakeG4DynamicParticle(Gammas->at(i)));
        }
        // Load the alphas
        for(unsigned int i = 0; i < Alphas->size(); i++)
        {
            FissionProducts->push_back(MakeG4DynamicParticle(Alphas->at(i)));
        }

    // Rotate the system to a random location so that the light fission fragment
    // is not always traveling along the positive z-axis
        // Sample Theta and Phi.
        G4ThreeVector RotationAxis;

        Theta = std::acos(RandomEngine_->G4SampleUniform(-1, 1));
        Phi = RandomEngine_->G4SampleUniform(-pi, pi);
        RotationAxis.setRThetaPhi(1.0, Theta, Phi);

        // We will also check the net momenta
        ResultantVector.set(0.0, 0.0, 0.0);
        for(unsigned int i = 0; i < FissionProducts->size(); i++)
        {
        	Direction = FissionProducts->at(i)->GetMomentumDirection();
        	Direction.rotateUz(RotationAxis);
        	FissionProducts->at(i)->SetMomentumDirection(Direction);
            ResultantVector += FissionProducts->at(i)->GetMomentum();
        }

        // Warn if the sum momenta of the system is not within a reasonable
        // tolerance
        G4double PossibleImbalance = ResultantVector.mag();
        if(PossibleImbalance > 0.01)
        {
            std::ostringstream Temp;
            Temp << "Momenta imbalance of ";
            Temp << PossibleImbalance / (MeV / CLHEP::c_light);
            Temp << " MeV/c in the system";
            G4Exception("G4FissionProductYieldDist::G4GetFission()",
                        Temp.str().c_str(),
                        JustWarning,
                        "Results may not be valid");
        }

    // Clean up
    delete FirstDaughter;
    delete SecondDaughter;
    G4ArrayOps::DeleteVectorOfPointers(*Neutrons);
    G4ArrayOps::DeleteVectorOfPointers(*Gammas);
    G4ArrayOps::DeleteVectorOfPointers(*Alphas);

G4FFG_FUNCTIONLEAVE__
    return FissionProducts;
}

G4Ions* G4FissionProductYieldDist::
G4GetFissionProduct( void )
{
G4FFG_FUNCTIONENTER__

    G4Ions* Product = FindParticle(RandomEngine_->G4SampleUniform());

G4FFG_FUNCTIONLEAVE__
    return Product;
}

void G4FissionProductYieldDist::
G4SetAlphaProduction(G4double WhatAlphaProduction)
{
G4FFG_FUNCTIONENTER__

    AlphaProduction_ = WhatAlphaProduction;

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
G4SetEnergy( G4double WhatIncidentEnergy )
{
G4FFG_FUNCTIONENTER__

    if(Cause_ != G4FFGEnumerations::SPONTANEOUS)
    {
        IncidentEnergy_ = WhatIncidentEnergy;
    } else
    {
        IncidentEnergy_ = 0 * GeV;
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
G4SetTernaryProbability( G4double WhatTernaryProbability )
{
G4FFG_FUNCTIONENTER__

    TernaryProbability_ = WhatTernaryProbability;

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
G4SetVerbosity(G4int Verbosity)
{
G4FFG_FUNCTIONENTER__

    Verbosity_ = Verbosity;
    
    ENDFData_->G4SetVerbosity(Verbosity_);
    RandomEngine_->G4SetVerbosity(Verbosity_);

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
CheckAlphaSanity( void )
{
G4FFG_FUNCTIONENTER__

    // This provides comfortable breathing room at 16 MeV per alpha
    if(AlphaProduction_ > 10)
    {
        AlphaProduction_ = 10;
    } else if(AlphaProduction_ < -7)
    {
        AlphaProduction_ = -7;
    }

G4FFG_FUNCTIONLEAVE__
}

G4Ions* G4FissionProductYieldDist::
FindParticle( G4double RandomParticle )
{
G4FFG_FUNCTIONENTER__

    // Determine which energy group is currently in use
    G4bool isExact = false;
    G4bool lowerExists = false;
    G4bool higherExists = false;
    G4int energyGroup;
    for(energyGroup = 0; energyGroup < YieldEnergyGroups_; energyGroup++)
    {
        if(IncidentEnergy_ == YieldEnergies_[energyGroup])
        {
            isExact = true;
            break;
        }

        if(energyGroup == 0 && IncidentEnergy_ < YieldEnergies_[energyGroup])
        {
            // Break if the energy is less than the lowest energy
            higherExists = true;
            break;
        } else if(energyGroup == YieldEnergyGroups_ - 1)
        {
            // The energy is greater than any values in the yield data.
            lowerExists = true;
            break;
        } else
        {
            // Break if the energy is less than the lowest energy
            if(IncidentEnergy_ > YieldEnergies_[energyGroup])
            {
                energyGroup--;
                lowerExists = true;
                higherExists = true;
                break;
            }
        }
    }

    // Determine which particle it is
    G4Ions* FoundParticle = NULL;
    if(isExact || YieldEnergyGroups_ == 1)
    {
        // Determine which tree contains the random value
        G4int tree;
        for(tree = 0; tree < TreeCount_; tree++)
        {
            // Break if a tree is identified as containing the random particle
            if(RandomParticle <= Trees_[tree].ProbabilityRangeEnd[energyGroup])
            {
                break;
            }
        }
        ProbabilityBranch* Branch = Trees_[tree].Trunk;

        // Iteratively traverse the tree until the particle addressed by the random
        // variable is found
        G4bool RangeIsSmaller;
        G4bool RangeIsGreater;
        while((RangeIsSmaller = (RandomParticle < Branch->ProbabilityRangeBottom[energyGroup]))
               || (RangeIsGreater = (RandomParticle > Branch->ProbabilityRangeTop[energyGroup])))
        // Loop checking, 11.05.2015, T. Koi
        {
            if(RangeIsSmaller)
            {
                Branch = Branch->Left;
            } else {
                Branch = Branch->Right;
            }
        }

        FoundParticle = Branch->Particle;
    } else if(lowerExists && higherExists)
    {
        // We need to do some interpolation
        FoundParticle = FindParticleInterpolation(RandomParticle, energyGroup);
    } else
    {
        // We need to do some extrapolation
        FoundParticle = FindParticleExtrapolation(RandomParticle, lowerExists);
    }

    // Return the particle
G4FFG_FUNCTIONLEAVE__
    return FoundParticle;
}

G4Ions* G4FissionProductYieldDist::
FindParticleExtrapolation( G4double RandomParticle,
                           G4bool LowerEnergyGroupExists )
{
G4FFG_FUNCTIONENTER__

    G4Ions* FoundParticle = NULL;
    G4int NearestEnergy;
    G4int NextNearestEnergy;

    // Check to see if we are extrapolating above or below the data set    
    if(LowerEnergyGroupExists == true)
    {
        NearestEnergy = YieldEnergyGroups_ - 1;
        NextNearestEnergy = NearestEnergy - 1;
    } else
    {
        NearestEnergy = 0;
        NextNearestEnergy = 1;
    }

    for(G4int Tree = 0; Tree < TreeCount_ && FoundParticle == NULL; Tree++)
    {
        FoundParticle = FindParticleBranchSearch(Trees_[Tree].Trunk,
                                                 RandomParticle,
                                                 NearestEnergy,
                                                 NextNearestEnergy);
    }

G4FFG_FUNCTIONLEAVE__
    return FoundParticle;
}

G4Ions* G4FissionProductYieldDist::
FindParticleInterpolation( G4double RandomParticle,
                           G4int LowerEnergyGroup )
{
G4FFG_FUNCTIONENTER__

    G4Ions* FoundParticle = NULL;
    G4int HigherEnergyGroup = LowerEnergyGroup + 1;

    for(G4int Tree = 0; Tree < TreeCount_ && FoundParticle == NULL; Tree++)
    {
        FoundParticle = FindParticleBranchSearch(Trees_[Tree].Trunk,
                                                 RandomParticle,
                                                 LowerEnergyGroup,
                                                 HigherEnergyGroup);
    }

G4FFG_FUNCTIONLEAVE__
    return FoundParticle;
}

G4Ions* G4FissionProductYieldDist::
FindParticleBranchSearch( ProbabilityBranch* Branch,
                          G4double RandomParticle,
                          G4int EnergyGroup1,
                          G4int EnergyGroup2 )
{
G4FFG_RECURSIVE_FUNCTIONENTER__

    G4Ions* Particle;

    // Verify that the branch exists
    if(Branch == NULL)
    {
        Particle = NULL;
    } else if(EnergyGroup1 >= Branch->IncidentEnergiesCount
              || EnergyGroup2 >= Branch->IncidentEnergiesCount
              || EnergyGroup1 == EnergyGroup2
              || Branch->IncidentEnergies[EnergyGroup1] == Branch->IncidentEnergies[EnergyGroup2])
    {
        // Set NULL if any invalid conditions exist
        Particle = NULL;
    } else
    {
        // Everything check out - proceed
        G4Ions* FoundParticle = NULL;
        G4double Intercept;
        G4double Slope;
        G4double RangeAtIncidentEnergy;
        G4double Denominator = Branch->IncidentEnergies[EnergyGroup1] - Branch->IncidentEnergies[EnergyGroup2];
        
        // Calculate the lower probability bounds
        Slope = (Branch->ProbabilityRangeBottom[EnergyGroup1] - Branch->ProbabilityRangeBottom[EnergyGroup2]) / Denominator;
        Intercept = Branch->ProbabilityRangeBottom[EnergyGroup1] - Slope * Branch->IncidentEnergies[EnergyGroup1];
        RangeAtIncidentEnergy = Slope * IncidentEnergy_ + Intercept;

        // Go right if the particle is below the probability bounds
        if(RandomParticle < RangeAtIncidentEnergy)
        {
            FoundParticle = FindParticleBranchSearch(Branch->Left,
                                                     RandomParticle,
                                                     EnergyGroup1,
                                                     EnergyGroup2);
        } else
        {
            // Calculate the upper probability bounds
            Slope = (Branch->ProbabilityRangeTop[EnergyGroup1] - Branch->ProbabilityRangeTop[EnergyGroup2]) / Denominator;
            Intercept = Branch->ProbabilityRangeTop[EnergyGroup1] - Slope * Branch->IncidentEnergies[EnergyGroup1];
            RangeAtIncidentEnergy = Slope * IncidentEnergy_ + Intercept;

            // Go left if the particle is above the probability bounds
            if(RandomParticle > RangeAtIncidentEnergy)
            {
                FoundParticle = FindParticleBranchSearch(Branch->Right,
                                                         RandomParticle,
                                                         EnergyGroup1,
                                                         EnergyGroup2);
            } else
            {
                // If the particle is bounded then we found it!
                FoundParticle = Branch->Particle;
            }
        }
        
        Particle = FoundParticle;
    }

G4FFG_RECURSIVE_FUNCTIONLEAVE__
    return Particle;
}

void G4FissionProductYieldDist::
GenerateAlphas( std::vector< G4ReactionProduct* >* Alphas )
{
G4FFG_FUNCTIONENTER__

    // Throw the dice to determine if ternary fission occurs
    G4bool MakeAlphas = RandomEngine_->G4SampleUniform() <= TernaryProbability_;
    if(MakeAlphas)
    {
        G4int NumberOfAlphasToProduce;

        // Determine how many alpha particles to produce for the ternary fission
        if(AlphaProduction_ < 0)
        {
            NumberOfAlphasToProduce = RandomEngine_->G4SampleIntegerGaussian(AlphaProduction_ * -1,
                                                                             1,
                                                                             G4FFGEnumerations::POSITIVE);
        } else
        {
             NumberOfAlphasToProduce = (G4int)AlphaProduction_;
        }

        //TK modifed 131108
        //Alphas->resize(NumberOfAlphasToProduce);
        for(int i = 0; i < NumberOfAlphasToProduce; i++)
        {
            // Set the G4Ions as an alpha particle
            Alphas->push_back(new G4ReactionProduct(AlphaDefinition_));

            // Remove 4 nucleons (2 protons and 2 neutrons) for each alpha added
            RemainingZ_ -= 2;
            RemainingA_ -= 4;
        }
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
GenerateNeutrons( std::vector< G4ReactionProduct* >* Neutrons )
{
G4FFG_FUNCTIONENTER__

    G4int NeutronProduction;
    NeutronProduction = RandomEngine_->G4SampleIntegerGaussian(Nubar_, NubarWidth_, G4FFGEnumerations::POSITIVE);

    //TK modifed 131108
    //Neutrons->resize(NeutronProduction);
    for(int i = 0; i < NeutronProduction; i++)
    {
        // Define the fragment as a neutron
        Neutrons->push_back(new G4ReactionProduct(NeutronDefinition_));

        // Remove 1 nucleon for each neutron added
        RemainingA_--;
    }

G4FFG_FUNCTIONLEAVE__
}

G4Ions* G4FissionProductYieldDist::
GetParticleDefinition( G4int Product,
                       //TK modified 131108
                       //G4FFGEnumerations::MetaState MetaState )
                       G4FFGEnumerations::MetaState /*MetaState*/ )
{
G4FFG_DATA_FUNCTIONENTER__

    G4Ions* Temp;

    // Break Product down into its A and Z components
    G4int A = Product % 1000; // Extract A
    G4int Z = (Product - A) / 1000; // Extract Z

    // Check to see if the particle is registered using the PDG code
    // TODO Add metastable state when supported by G4IonTable::GetIon()
    Temp = static_cast<G4Ions*>(IonTable_->GetIon(Z, A));
    
    // Removed in favor of the G4IonTable::GetIon() method
//    // Register the particle if it does not exist
//    if(Temp == NULL)
//    {
//        // Define the particle properties
//        G4String Name = MakeIsotopeName(Product, MetaState);
//        // Calculate the rest mass using a function already in Geant4
//        G4double Mass = G4NucleiProperties::
//                        GetNuclearMass((double)A, (double)Z );
//        G4double Charge = Z*eplus;
//        G4int BaryonNum = A;
//        G4bool Stable = TRUE;
//
//        // I am unsure about the following properties:
//        //     2*Spin, Parity, C-conjugation, 2*Isospin, 2*Isospin3, G-parity.
//        // Perhaps is would be a good idea to have a physicist familiar with
//        // Geant4 nomenclature to review and correct these properties.
//        Temp = new G4Ions (
//            // Name           Mass           Width          Charge
//               Name,          Mass,          0.0,           Charge,
//
//            // 2*Spin         Parity         C-conjugation  2*Isospin
//               0,             1,             0,             0,
//
//            // 2*Isospin3     G-parity       Type           Lepton number
//               0,             0,             "nucleus",     0,
//
//            // Baryon number  PDG encoding   Stable         Lifetime
//               BaryonNum,     PDGCode,       Stable,        -1,
//
//            // Decay table    Shortlived     SubType        Anti_encoding
//               NULL,          FALSE,        "generic",     0,
//
//            // Excitation
//               0.0);
//        Temp->SetPDGMagneticMoment(0.0);
//
//        // Declare that there is no anti-particle
//        Temp->SetAntiPDGEncoding(0);
//
//        // Define the processes to use in transporting the particles
//        std::ostringstream osAdd;
//        osAdd << "/run/particle/addProcManager " << Name;
//        G4String cmdAdd = osAdd.str();
//
//        // set /control/verbose 0
//        G4int tempVerboseLevel = G4UImanager::GetUIpointer()->GetVerboseLevel();
//        G4UImanager::GetUIpointer()->SetVerboseLevel(0);
//
//        // issue /run/particle/addProcManage
//        G4UImanager::GetUIpointer()->ApplyCommand(cmdAdd);
//
//        // retrieve  /control/verbose
//        G4UImanager::GetUIpointer()->SetVerboseLevel(tempVerboseLevel);
//    }

G4FFG_DATA_FUNCTIONLEAVE__
    return Temp;
}

G4String G4FissionProductYieldDist::
MakeDirectoryName( void )
{
G4FFG_FUNCTIONENTER__

    // Generate the file location starting in the Geant4 data directory
    std::ostringstream DirectoryName;
    DirectoryName << G4FindDataDir("G4NEUTRONHPDATA") << G4FFGDefaultValues::ENDFFissionDataLocation;

    // Return the directory structure
G4FFG_FUNCTIONLEAVE__
    return DirectoryName.str();
}

G4String G4FissionProductYieldDist::
MakeFileName( G4int Isotope,
              G4FFGEnumerations::MetaState MetaState )
{
G4FFG_FUNCTIONENTER__


    // Create the unique identifying name for the particle
    std::ostringstream FileName;

    // Determine if a leading 0 is needed (ZZZAAA or 0ZZAAA)
    if(Isotope < 100000)
    {
        FileName << "0";
    }

    // Add the name of the element and the extension
    FileName << MakeIsotopeName(Isotope, MetaState) << ".fpy";

G4FFG_FUNCTIONLEAVE__
    return FileName.str();
}

G4DynamicParticle* G4FissionProductYieldDist::
MakeG4DynamicParticle( G4ReactionProduct* ReactionProduct )
{
G4FFG_DATA_FUNCTIONENTER__

	G4DynamicParticle* DynamicParticle = new G4DynamicParticle(ReactionProduct->GetDefinition(), ReactionProduct->GetMomentum());

G4FFG_DATA_FUNCTIONLEAVE__
	return DynamicParticle;
}

G4String G4FissionProductYieldDist::
MakeIsotopeName( G4int Isotope,
                 G4FFGEnumerations::MetaState MetaState )
{
G4FFG_DATA_FUNCTIONENTER__

    // Break Product down into its A and Z components
    G4int A = Isotope % 1000;
    G4int Z = (Isotope - A) / 1000;

    // Create the unique identifying name for the particle
    std::ostringstream IsotopeName;

    IsotopeName << Z << "_" << A;
    
    // If it is metastable then append "m" to the name
    if(MetaState != G4FFGEnumerations::GROUND_STATE)
    {
        IsotopeName << "m";
        
        // If it is a second isomeric state then append "2" to the name
        if(MetaState == G4FFGEnumerations::META_2)
        {
            IsotopeName << "2";
        }
    }
    // Add the name of the element and the extension
    IsotopeName << "_" << ElementNames_->theString[Z - 1];

G4FFG_DATA_FUNCTIONLEAVE__
    return IsotopeName.str();
}

void G4FissionProductYieldDist::
MakeTrees( void )
{
G4FFG_FUNCTIONENTER__

    // Allocate the space
    // We will make each tree a binary search
    // The maximum number of iterations required to find a single fission product
    // based on it's probability is defined by the following:
    //      x = number of fission products
    //      Trees       = T(x)  = ceil( ln(x) )
    //      Rows/Tree   = R(x)  = ceil(( sqrt( (8 * x / T(x)) + 1) - 1) / 2)
    //      Maximum     = M(x)  = T(x) + R(x)
    //      Results: x    =>    M(x)
    //               10         5
    //               100        10
    //               1000       25
    //               10000      54
    //               100000     140
    TreeCount_ = (G4int)ceil((G4double)log((G4double)ENDFData_->G4GetNumberOfFissionProducts()));
    Trees_ = new ProbabilityTree[TreeCount_];

    // Initialize the range of each node
    for(G4int i = 0; i < TreeCount_; i++)
    {
        Trees_[i].ProbabilityRangeEnd = new G4double[YieldEnergyGroups_];
        Trees_[i].Trunk = NULL;
        Trees_[i].BranchCount = 0;
        Trees_[i].IsEnd = FALSE;
    }
    // Mark the last tree as the ending tree
    Trees_[TreeCount_ - 1].IsEnd = TRUE;

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
ReadProbabilities( void )
{
G4FFG_DATA_FUNCTIONENTER__

    G4int ProductCount = ENDFData_->G4GetNumberOfFissionProducts();
    BranchCount_ = 0;
    G4ArrayOps::Set(YieldEnergyGroups_, DataTotal_, 0.0);

    // Loop through all the products
    for(G4int i = 0; i < ProductCount; i++)
    {
        // Acquire the data and sort it
        SortProbability(ENDFData_->G4GetYield(i));
    }

    // Generate the true normalization factor, since round-off errors may result
    // in non-singular normalization of the data files. Also, reset DataTotal_
    // since it is used by Renormalize() to set the probability segments.
    G4ArrayOps::Divide(YieldEnergyGroups_, MaintainNormalizedData_, 1.0, DataTotal_);
    G4ArrayOps::Set(YieldEnergyGroups_, DataTotal_, 0.0);

    // Go through all the trees one at a time
    for(G4int i = 0; i < TreeCount_; i++)
    {
        Renormalize(Trees_[i].Trunk);
        // Set the max range of the tree to DataTotal
        G4ArrayOps::Copy(YieldEnergyGroups_, Trees_[i].ProbabilityRangeEnd, DataTotal_);
    }

G4FFG_DATA_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
Renormalize( ProbabilityBranch* Branch)
{
G4FFG_RECURSIVE_FUNCTIONENTER__

    // Check to see if Branch exists. Branch will be a null pointer if it
    // doesn't exist
    if(Branch != NULL)
    {
        // Call the lower branch to set the probability segment first, since it
        // supposed to have a lower probability segment that this node
        Renormalize(Branch->Left);

        // Set this node as the next sequential probability segment
        G4ArrayOps::Copy(YieldEnergyGroups_, Branch->ProbabilityRangeBottom, DataTotal_);
        G4ArrayOps::Multiply(YieldEnergyGroups_, Branch->ProbabilityRangeTop, MaintainNormalizedData_);
        G4ArrayOps::Add(YieldEnergyGroups_, Branch->ProbabilityRangeTop, DataTotal_);
        G4ArrayOps::Copy(YieldEnergyGroups_, DataTotal_, Branch->ProbabilityRangeTop);

        // Now call the upper branch to set those probability segments
        Renormalize(Branch->Right);
    }

G4FFG_RECURSIVE_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
SampleAlphaEnergies( std::vector< G4ReactionProduct* >* Alphas )
{
G4FFG_FUNCTIONENTER__

    // The condition of sampling more energy from the fission products than is
    // alloted is statistically unfavorable, but it could still happen. The
    // do-while loop prevents such an occurrence from happening
    G4double MeanAlphaEnergy = 16.0;
    G4double TotalAlphaEnergy;

    do
    {
        G4double AlphaEnergy;
        TotalAlphaEnergy = 0;

        // Walk through the alpha particles one at a time and sample each's
        // energy
        for(unsigned int i = 0; i < Alphas->size(); i++)
        {
            AlphaEnergy = RandomEngine_->G4SampleGaussian(MeanAlphaEnergy,
                                                          2.35,
                                                          G4FFGEnumerations::POSITIVE) * MeV;
            // Assign the energy to the alpha particle
            Alphas->at(i)->SetKineticEnergy(AlphaEnergy);

            // Add up the total amount of kinetic energy consumed.
            TotalAlphaEnergy += AlphaEnergy;
        }

        // If true, decrement the mean alpha energy by 0.1 and try again.
        MeanAlphaEnergy -= 0.1;
    } while(TotalAlphaEnergy >= RemainingEnergy_); // Loop checking, 11.05.2015, T. Koi

    // Subtract the total amount of energy that was assigned.
    RemainingEnergy_ -= TotalAlphaEnergy;

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
SampleGammaEnergies( std::vector< G4ReactionProduct* >* Gammas )
{
G4FFG_FUNCTIONENTER__

	// Make sure that there is energy to assign to the gamma rays
    if(RemainingEnergy_ != 0)
    {
        G4double SampleEnergy;

        // Sample from RemainingEnergy until it is all gone. Also,
        // RemainingEnergy should not be smaller than
        // G4FFGDefaultValues::MeanGammaEnergy. This will prevent the
        // sampling of a fractional portion of the Gaussian distribution
        // in an attempt to find a new gamma ray energy.
        G4int icounter=0;
        G4int icounter_max=1024;
        while(RemainingEnergy_ >= G4FFGDefaultValues::MeanGammaEnergy ) // Loop checking, 11.05.2015, T. Koi
        {
           icounter++;
           if ( icounter > icounter_max ) {
	      G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
              break;
           }
            SampleEnergy = RandomEngine_->
                    G4SampleGaussian(G4FFGDefaultValues::MeanGammaEnergy, 1.0 * MeV, G4FFGEnumerations::POSITIVE);
            // Make sure that we didn't sample more energy than was available
            if(SampleEnergy <= RemainingEnergy_)
            {
                // If this energy assignment would leave less energy than the
                // 'intrinsic' minimal energy of a gamma ray then just assign
                // all of the remaining energy
                if(RemainingEnergy_ - SampleEnergy < 100 * keV)
                {
                    SampleEnergy = RemainingEnergy_;
                }

                // Create the new particle
                Gammas->push_back(new G4ReactionProduct());

                // Set the properties
                Gammas->back()->SetDefinition(GammaDefinition_);
                Gammas->back()->SetTotalEnergy(SampleEnergy);

                // Calculate how much is left
                RemainingEnergy_ -= SampleEnergy;
            }
        }

        // If there is anything left over, the energy must be above 100 keV but
        // less than G4FFGDefaultValues::MeanGammaEnergy. Arbitrarily assign
        // RemainingEnergy to a new particle
        if(RemainingEnergy_ > 0)
        {
        	SampleEnergy = RemainingEnergy_;
            Gammas->push_back(new G4ReactionProduct());

            // Set the properties
            Gammas->back()->SetDefinition(GammaDefinition_);
            Gammas->back()->SetTotalEnergy(SampleEnergy);

            // Calculate how much is left
            RemainingEnergy_ -= SampleEnergy;
        }
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
SampleNeutronEnergies( std::vector< G4ReactionProduct* >* Neutrons )
{
G4FFG_FUNCTIONENTER__

    // The condition of sampling more energy from the fission products than is
    // alloted is statistically unfavorable, but it could still happen. The
    // do-while loop prevents such an occurrence from happening
    G4double TotalNeutronEnergy=0.;
    G4double NeutronEnergy=0.;

    // Make sure that we don't sample more energy than is available
    G4int icounter=0;
    G4int icounter_max=1024;
    do
    {
       icounter++;
       if ( icounter > icounter_max ) {
          G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
           break;
       }
        TotalNeutronEnergy = 0;

        // Walk through the neutrons one at a time and sample the energies.
        // The gamma rays have not yet been sampled, so the last neutron will
        // have a NULL value for NextFragment
        for(unsigned int i = 0; i < Neutrons->size(); i++)
        {
            // Assign the energy to the neutron
        	NeutronEnergy = RandomEngine_->G4SampleWatt(Isotope_, Cause_, IncidentEnergy_);
        	Neutrons->at(i)->SetKineticEnergy(NeutronEnergy);

            // Add up the total amount of kinetic energy consumed.
            TotalNeutronEnergy +=NeutronEnergy;
        }
    } while (TotalNeutronEnergy > RemainingEnergy_); // Loop checking, 11.05.2015, T. Koi

    // Subtract the total amount of energy that was assigned.
    RemainingEnergy_ -= TotalNeutronEnergy;

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
SetNubar( void )
{
G4FFG_FUNCTIONENTER__

    G4int* WhichNubar;
    G4int* NubarWidth;
    G4double XFactor, BFactor;

    if(Cause_ == G4FFGEnumerations::SPONTANEOUS)
    {
        WhichNubar = const_cast<G4int*>(&SpontaneousNubar_[0][0]);
        NubarWidth = const_cast<G4int*>(&SpontaneousNubarWidth_[0][0]);
    } else
    {
        WhichNubar = const_cast<G4int*>(&NeutronInducedNubar_[0][0]);
        NubarWidth = const_cast<G4int*>(&NeutronInducedNubarWidth_[0][0]);
    }

    XFactor = G4Pow::GetInstance()->powA(10.0, -13.0);
    BFactor = G4Pow::GetInstance()->powA(10.0, -4.0);
    Nubar_ = *(WhichNubar + 1) * IncidentEnergy_ * XFactor
             + *(WhichNubar + 2) * BFactor;
    while(*WhichNubar != -1) // Loop checking, 11.05.2015, T. Koi
    {
        if(*WhichNubar == Isotope_)
        {
            Nubar_ = *(WhichNubar + 1) * IncidentEnergy_ * XFactor
                     + *(WhichNubar + 2) * BFactor;

            break;
        }
        WhichNubar += 3;
    }

    XFactor = G4Pow::GetInstance()->powN((G4double)10, -6);
    NubarWidth_ = *(NubarWidth + 1) * XFactor;
    while(*WhichNubar != -1) // Loop checking, 11.05.2015, T. Koi
    {
        if(*WhichNubar == Isotope_)
        {
            NubarWidth_ = *(NubarWidth + 1) * XFactor;

            break;
        }
        WhichNubar += 2;
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
SortProbability( G4ENDFYieldDataContainer* YieldData )
{
G4FFG_DATA_FUNCTIONENTER__

    // Initialize the new branch
    ProbabilityBranch* NewBranch = new ProbabilityBranch;
    NewBranch->IncidentEnergiesCount = YieldEnergyGroups_;
    NewBranch->Left = NULL;
    NewBranch->Right = NULL;
    NewBranch->Particle = GetParticleDefinition(YieldData->GetProduct(), YieldData->GetMetaState());
    NewBranch->IncidentEnergies = new G4double[YieldEnergyGroups_];
    NewBranch->ProbabilityRangeTop = new G4double[YieldEnergyGroups_];
    NewBranch->ProbabilityRangeBottom = new G4double[YieldEnergyGroups_];
    G4ArrayOps::Copy(YieldEnergyGroups_, NewBranch->ProbabilityRangeTop, YieldData->GetYieldProbability());
    G4ArrayOps::Copy(YieldEnergyGroups_, NewBranch->IncidentEnergies, YieldEnergies_);
    G4ArrayOps::Add(YieldEnergyGroups_, DataTotal_, YieldData->GetYieldProbability());

    // Check to see if the this is the smallest/largest particle. First, check
    // to see if this is the first particle in the system
    if(SmallestZ_ == NULL)
    {
        SmallestZ_ = SmallestA_ = LargestZ_ = LargestA_ = NewBranch->Particle;
    } else
    {
        G4bool IsSmallerZ = NewBranch->Particle->GetAtomicNumber() < SmallestZ_->GetAtomicNumber();
        G4bool IsSmallerA = NewBranch->Particle->GetAtomicMass() < SmallestA_->GetAtomicMass();
        G4bool IsLargerZ = NewBranch->Particle->GetAtomicNumber() > LargestZ_->GetAtomicNumber();
        G4bool IsLargerA = NewBranch->Particle->GetAtomicMass() > LargestA_->GetAtomicMass();

        if(IsSmallerZ)
        {
            SmallestZ_ = NewBranch->Particle;
        }
        
        if(IsLargerZ)
        {
            LargestA_ = NewBranch->Particle;
        }
        
        if(IsSmallerA)
        {
            SmallestA_ = NewBranch->Particle;
        }
        
        if(IsLargerA)
        {
            LargestA_ = NewBranch->Particle;
        }
    }

    // Place the new branch
    // Determine which tree the new branch goes into
    G4int WhichTree = (G4int)floor((G4double)(BranchCount_ % TreeCount_));
    ProbabilityBranch** WhichBranch = &(Trees_[WhichTree].Trunk);
    Trees_[WhichTree].BranchCount++;

    // Search for the position
    // Determine where the branch goes
    G4int BranchPosition = (G4int)floor((G4double)(BranchCount_ / TreeCount_)) + 1;

    // Run through the tree until the end branch is reached
    while(BranchPosition > 1) // Loop checking, 11.05.2015, T. Koi
    {
        if(BranchPosition & 1)
        {
            // If the 1's bit is on then move to the next 'right' branch
            WhichBranch = &((*WhichBranch)->Right);
        } else
        {
            // If the 1's bit is off then move to the next 'down' branch
            WhichBranch = &((*WhichBranch)->Left);
        }
        
        BranchPosition >>= 1;
    }

    *WhichBranch = NewBranch;
    BranchCount_++;

G4FFG_DATA_FUNCTIONLEAVE__
}

G4FissionProductYieldDist::
~G4FissionProductYieldDist( void )
{
G4FFG_FUNCTIONENTER__

    // Burn each tree, one by one
    G4int WhichTree = 0;
    while(Trees_[WhichTree].IsEnd != TRUE) // Loop checking, 11.05.2015, T. Koi
    {
        BurnTree(Trees_[WhichTree].Trunk);
        delete Trees_[WhichTree].Trunk;
        delete[] Trees_[WhichTree].ProbabilityRangeEnd;
        WhichTree++;
    }

    // Delete each dynamically allocated variable
    delete ENDFData_;
    delete[] Trees_;
    delete[] DataTotal_;
    delete[] MaintainNormalizedData_;
    delete ElementNames_;
    delete RandomEngine_;

G4FFG_FUNCTIONLEAVE__
}

void G4FissionProductYieldDist::
BurnTree( ProbabilityBranch* Branch )
{
G4FFG_RECURSIVE_FUNCTIONENTER__

    // Check to see it Branch exists. Branch will be a null pointer if it
    // doesn't exist
    if(Branch)
    {
        // Burn down before you burn up
        BurnTree(Branch->Left);
        delete Branch->Left;
        BurnTree(Branch->Right);
        delete Branch->Right;

        delete[] Branch->IncidentEnergies;
        delete[] Branch->ProbabilityRangeTop;
        delete[] Branch->ProbabilityRangeBottom;
    }

G4FFG_RECURSIVE_FUNCTIONLEAVE__
}


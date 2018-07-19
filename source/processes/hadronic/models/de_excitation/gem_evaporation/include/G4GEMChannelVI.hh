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
// $Id: G4GEMChannelVI.hh 98577 2016-07-25 13:05:12Z vnivanch $
//
// GEM de-excitation model
// by V. Ivanchenko (July 2016)
//

#ifndef G4GEMChannelVI_h
#define G4GEMChannelVI_h 1

#include "G4VEvaporationChannel.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4NucleiProperties.hh"
#include "G4VCoulombBarrier.hh"
#include "G4Exp.hh"

class G4Pow;
class G4PairingCorrection;
class G4VCoulombBarrier;
class G4LevelManager;
class G4NuclearLevelData;

const G4int NPOINTSGEM = 10;

class G4GEMChannelVI : public G4VEvaporationChannel
{
public:

  explicit G4GEMChannelVI(G4int theA, G4int theZ);

  virtual ~G4GEMChannelVI();
    
  virtual void Initialise() final;

  virtual G4double GetEmissionProbability(G4Fragment* theNucleus) final;

  virtual G4Fragment* EmittedFragment(G4Fragment* theNucleus) final;

  virtual void Dump() const;

private: 

  G4double IntegratedProbability(G4double exc);

  G4double ProbabilityDistributionFunction(G4double exc, G4double resExc);

  G4double FindLevel(const G4LevelManager*, 
		     G4double exc, G4double exclim);

  inline G4double I0(G4double t);
  inline G4double I1(G4double t, G4double tx);
  inline G4double I2(G4double s0, G4double sx);
  G4double I3(G4double s0, G4double sx);

  G4GEMChannelVI(const G4GEMChannelVI & right) = delete;  
  const G4GEMChannelVI & operator=(const G4GEMChannelVI & right) = delete;
  G4bool operator==(const G4GEMChannelVI & right) const = delete;
  G4bool operator!=(const G4GEMChannelVI & right) const = delete;

  G4Pow* fG4pow;
    
  const G4VCoulombBarrier* cBarrier;

  const G4PairingCorrection* pairingCorrection;

  const G4LevelManager* levelManager;

  G4NuclearLevelData* nData;
  
  G4int A;
  G4int Z;
  G4int resA;
  G4int resZ;
  G4int fragA;
  G4int fragZ;
  G4int nWarn;

  G4double massGround;
  G4double maxLevelE;
  G4double Z13;
  G4double A13;

  G4double massFrag;
  G4double eCBarrier;
  G4double resMassGround;
  G4double maxKinEnergy;
  G4double resZ13;
  G4double resA13;
  G4double delta0;
  G4double delta1;

  G4double alphaP;
  G4double betaP;

  G4double maxExc;
  G4double maxProb;
  G4double coeff;
  G4double levelDensity;

  static const G4double ws[NPOINTSGEM];
  static const G4double xs[NPOINTSGEM];

};

inline G4double G4GEMChannelVI::I0(G4double t)
{
  return G4Exp(t) - 1.0;
}

inline G4double G4GEMChannelVI::I1(G4double t, G4double tx)
{
  return (t - tx + 1.0)*G4Exp(tx) - t - 1.0;
}


inline G4double G4GEMChannelVI::I2(G4double s0, G4double sx)
{
  G4double S = 1.0/std::sqrt(s0);
  G4double Sx = 1.0/std::sqrt(sx);
  
  G4double p1 = S*S*S*( 1.0 + S*S*( 1.5 + 3.75*S*S) );
  G4double p2 = Sx*Sx*Sx*( 1.0 + Sx*Sx*( 1.5 + 3.75*Sx*Sx) )*G4Exp(sx-s0);
  
  return p1-p2;
}

#endif

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
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation,
//      21-5-98 V.Grichine
//      28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation
//      04.03.05, V.Grichine: get local field interface
//      18-05-06 H. Burkhardt: Energy spectrum from function rather than table
//
///////////////////////////////////////////////////////////////////////////

#include "G4SynchrotronRadiation.hh"

#include "G4DipBustGenerator.hh"
#include "G4Electron.hh"
#include "G4EmProcessSubType.hh"
#include "G4Log.hh"
#include "G4LossTableManager.hh"
#include "G4Gamma.hh"
#include "G4PhysicalConstants.hh"
#include "G4PropagatorInField.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicsModelCatalog.hh"

///////////////////////////////////////////////////////////////////////
//  Constructor
G4SynchrotronRadiation::G4SynchrotronRadiation(const G4String& processName,
                                               G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
  , theGamma(G4Gamma::Gamma())
{
  G4TransportationManager* transportMgr =
    G4TransportationManager::GetTransportationManager();

  fFieldPropagator = transportMgr->GetPropagatorInField();

  secID = G4PhysicsModelCatalog::GetModelID("model_SynRad");
  SetProcessSubType(fSynchrotronRadiation);
  verboseLevel = 1;
  FirstTime    = true;
  FirstTime1   = true;
  genAngle     = nullptr;
  SetAngularGenerator(new G4DipBustGenerator());
  theManager = G4LossTableManager::Instance();
  theManager->Register(this);
}

/////////////////////////////////////////////////////////////////////////
// Destructor
G4SynchrotronRadiation::~G4SynchrotronRadiation()
{
  delete genAngle;
  theManager->DeRegister(this);
}

/////////////////////////////// METHODS /////////////////////////////////

void G4SynchrotronRadiation::SetAngularGenerator(G4VEmAngularDistribution* p)
{
  if(p != genAngle)
  {
    delete genAngle;
    genAngle = p;
  }
}

G4bool G4SynchrotronRadiation::IsApplicable(
  const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0.0 && !particle.IsShortLived());
}

/////////////////////////////////////////////////////////////////////////
// Production of synchrotron X-ray photon
// Geant4 internal units.
G4double G4SynchrotronRadiation::GetMeanFreePath(const G4Track& trackData,
                                                 G4double,
                                                 G4ForceCondition* condition)
{
  // gives the MeanFreePath in Geant4 internal units
  G4double MeanFreePath = DBL_MAX;

  const G4DynamicParticle* aDynamicParticle = trackData.GetDynamicParticle();

  *condition = NotForced;

  G4double gamma =
    aDynamicParticle->GetTotalEnergy() / aDynamicParticle->GetMass();

  G4double particleCharge = aDynamicParticle->GetDefinition()->GetPDGCharge();

  if(gamma < 1.0e3 || 0.0 == particleCharge)
  {
    MeanFreePath = DBL_MAX;
  }
  else
  {
    G4ThreeVector FieldValue;
    const G4Field* pField   = nullptr;
    G4bool fieldExertsForce = false;

    G4FieldManager* fieldMgr =
      fFieldPropagator->FindAndSetFieldManager(trackData.GetVolume());

    if(fieldMgr != nullptr)
    {
      // If the field manager has no field, there is no field !
      fieldExertsForce = (fieldMgr->GetDetectorField() != nullptr);
    }

    if(fieldExertsForce)
    {
      pField                     = fieldMgr->GetDetectorField();
      G4ThreeVector globPosition = trackData.GetPosition();

      G4double globPosVec[4], FieldValueVec[6];

      globPosVec[0] = globPosition.x();
      globPosVec[1] = globPosition.y();
      globPosVec[2] = globPosition.z();
      globPosVec[3] = trackData.GetGlobalTime();

      pField->GetFieldValue(globPosVec, FieldValueVec);

      FieldValue =
        G4ThreeVector(FieldValueVec[0], FieldValueVec[1], FieldValueVec[2]);

      G4ThreeVector unitMomentum = aDynamicParticle->GetMomentumDirection();
      G4ThreeVector unitMcrossB  = FieldValue.cross(unitMomentum);
      G4double perpB             = unitMcrossB.mag();

      static const G4double fLambdaConst =
        std::sqrt(3.0) * eplus / (2.5 * fine_structure_const * c_light);

      if(perpB > 0.0)
      {
        MeanFreePath = fLambdaConst *
                       aDynamicParticle->GetDefinition()->GetPDGMass() /
                       (perpB * particleCharge * particleCharge);
      }
      if(verboseLevel > 0 && FirstTime)
      {
        G4cout << "G4SynchrotronRadiation::GetMeanFreePath "
               << " for particle "
               << aDynamicParticle->GetDefinition()->GetParticleName() << ":"
               << '\n'
               << "  MeanFreePath = " << G4BestUnit(MeanFreePath, "Length")
               << G4endl;
        if(verboseLevel > 1)
        {
          G4ThreeVector pvec = aDynamicParticle->GetMomentum();
          G4double Btot      = FieldValue.getR();
          G4double ptot      = pvec.getR();
          G4double rho       = ptot / (MeV * c_light * Btot);
          // full bending radius
          G4double Theta = unitMomentum.theta(FieldValue);
          // angle between particle and field
          G4cout << "  B = " << Btot / tesla << " Tesla"
                 << "  perpB = " << perpB / tesla << " Tesla"
                 << "  Theta = " << Theta
                 << " std::sin(Theta)=" << std::sin(Theta) << '\n'
                 << "  ptot  = " << G4BestUnit(ptot, "Energy")
                 << "  rho   = " << G4BestUnit(rho, "Length") << G4endl;
        }
        FirstTime = false;
      }
    }
  }
  return MeanFreePath;
}

///////////////////////////////////////////////////////////////////////////////
G4VParticleChange* G4SynchrotronRadiation::PostStepDoIt(
  const G4Track& trackData, const G4Step& stepData)

{
  aParticleChange.Initialize(trackData);

  const G4DynamicParticle* aDynamicParticle = trackData.GetDynamicParticle();

  G4double gamma = aDynamicParticle->GetTotalEnergy() /
                   (aDynamicParticle->GetDefinition()->GetPDGMass());

  G4double particleCharge = aDynamicParticle->GetDefinition()->GetPDGCharge();
  if(gamma <= 1.0e3 || 0.0 == particleCharge)
  {
    return G4VDiscreteProcess::PostStepDoIt(trackData, stepData);
  }

  G4ThreeVector FieldValue;
  const G4Field* pField = nullptr;

  G4bool fieldExertsForce = false;
  G4FieldManager* fieldMgr =
    fFieldPropagator->FindAndSetFieldManager(trackData.GetVolume());

  if(fieldMgr != nullptr)
  {
    // If the field manager has no field, there is no field !
    fieldExertsForce = (fieldMgr->GetDetectorField() != nullptr);
  }

  if(fieldExertsForce)
  {
    pField                     = fieldMgr->GetDetectorField();
    G4ThreeVector globPosition = trackData.GetPosition();
    G4double globPosVec[4], FieldValueVec[6];
    globPosVec[0] = globPosition.x();
    globPosVec[1] = globPosition.y();
    globPosVec[2] = globPosition.z();
    globPosVec[3] = trackData.GetGlobalTime();

    pField->GetFieldValue(globPosVec, FieldValueVec);
    FieldValue =
      G4ThreeVector(FieldValueVec[0], FieldValueVec[1], FieldValueVec[2]);

    G4ThreeVector unitMomentum = aDynamicParticle->GetMomentumDirection();
    G4ThreeVector unitMcrossB  = FieldValue.cross(unitMomentum);
    G4double perpB             = unitMcrossB.mag();
    if(perpB > 0.0)
    {
      // M-C of synchrotron photon energy
      G4double energyOfSR = GetRandomEnergySR(
        gamma, perpB, aDynamicParticle->GetDefinition()->GetPDGMass());

      // check against insufficient energy
      if(energyOfSR <= 0.0)
      {
        return G4VDiscreteProcess::PostStepDoIt(trackData, stepData);
      }
      G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();
      G4ThreeVector gammaDirection =
        genAngle->SampleDirection(aDynamicParticle, energyOfSR, 1, nullptr);

      G4ThreeVector gammaPolarization = FieldValue.cross(gammaDirection);
      gammaPolarization               = gammaPolarization.unit();

      // create G4DynamicParticle object for the SR photon
      auto aGamma =
        new G4DynamicParticle(theGamma, gammaDirection, energyOfSR);
      aGamma->SetPolarization(gammaPolarization.x(), gammaPolarization.y(),
                              gammaPolarization.z());

      aParticleChange.SetNumberOfSecondaries(1);

      // Update the incident particle
      G4double newKinEnergy = kineticEnergy - energyOfSR;

      if(newKinEnergy > 0.)
      {
        aParticleChange.ProposeEnergy(newKinEnergy);
      }
      else
      {
        aParticleChange.ProposeEnergy(0.);
      }

      // Create the G4Track
      G4Track* aSecondaryTrack = new G4Track(aGamma, trackData.GetGlobalTime(), trackData.GetPosition());
      aSecondaryTrack->SetTouchableHandle(stepData.GetPostStepPoint()->GetTouchableHandle());
      aSecondaryTrack->SetParentID(trackData.GetTrackID());
      aSecondaryTrack->SetCreatorModelID(secID);
      aParticleChange.AddSecondary(aSecondaryTrack);

    }
  }
  return G4VDiscreteProcess::PostStepDoIt(trackData, stepData);
}

///////////////////////////////////////////////////////////////////////////////
G4double G4SynchrotronRadiation::InvSynFracInt(G4double x)
// direct generation
{
  // from 0 to 0.7
  static constexpr G4double aa1           = 0;
  static constexpr G4double aa2           = 0.7;
  static constexpr G4int ncheb1           = 27;
  static constexpr G4double cheb1[ncheb1] = {
    1.22371665676046468821,     0.108956475422163837267,
    0.0383328524358594396134,   0.00759138369340257753721,
    0.00205712048644963340914,  0.000497810783280019308661,
    0.000130743691810302187818, 0.0000338168760220395409734,
    8.97049680900520817728e-6,  2.38685472794452241466e-6,
    6.41923109149104165049e-7,  1.73549898982749277843e-7,
    4.72145949240790029153e-8,  1.29039866111999149636e-8,
    3.5422080787089834182e-9,   9.7594757336403784905e-10,
    2.6979510184976065731e-10,  7.480422622550977077e-11,
    2.079598176402699913e-11,   5.79533622220841193e-12,
    1.61856011449276096e-12,    4.529450993473807e-13,
    1.2698603951096606e-13,     3.566117394511206e-14,
    1.00301587494091e-14,       2.82515346447219e-15,
    7.9680747949792e-16
  };
  //   from 0.7 to 0.9132260271183847
  static constexpr G4double aa3           = 0.9132260271183847;
  static constexpr G4int ncheb2           = 27;
  static constexpr G4double cheb2[ncheb2] = {
    1.1139496701107756,     0.3523967429328067,     0.0713849171926623,
    0.01475818043595387,    0.003381255637322462,   0.0008228057599452224,
    0.00020785506681254216, 0.00005390169253706556, 0.000014250571923902464,
    3.823880733161044e-6,   1.0381966089136036e-6,  2.8457557457837253e-7,
    7.86223332179956e-8,    2.1866609342508474e-8,  6.116186259857143e-9,
    1.7191233618437565e-9,  4.852755117740807e-10,  1.3749966961763457e-10,
    3.908961987062447e-11,  1.1146253766895824e-11, 3.1868887323415814e-12,
    9.134319791300977e-13,  2.6211077371181566e-13, 7.588643377757906e-14,
    2.1528376972619e-14,    6.030906040404772e-15,  1.9549163926819867e-15
  };
  // Chebyshev with exp/log  scale
  // a = -Log[1 - SynFracInt[1]]; b = -Log[1 - SynFracInt[7]];
  static constexpr G4double aa4           = 2.4444485538746025480;
  static constexpr G4double aa5           = 9.3830728608909477079;
  static constexpr G4int ncheb3           = 28;
  static constexpr G4double cheb3[ncheb3] = {
    1.2292683840435586977,        0.160353449247864455879,
    -0.0353559911947559448721,    0.00776901561223573936985,
    -0.00165886451971685133259,   0.000335719118906954279467,
    -0.0000617184951079161143187, 9.23534039743246708256e-6,
    -6.06747198795168022842e-7,   -3.07934045961999778094e-7,
    1.98818772614682367781e-7,    -8.13909971567720135413e-8,
    2.84298174969641838618e-8,    -9.12829766621316063548e-9,
    2.77713868004820551077e-9,    -8.13032767247834023165e-10,
    2.31128525568385247392e-10,   -6.41796873254200220876e-11,
    1.74815310473323361543e-11,   -4.68653536933392363045e-12,
    1.24016595805520752748e-12,   -3.24839432979935522159e-13,
    8.44601465226513952994e-14,   -2.18647276044246803998e-14,
    5.65407548745690689978e-15,   -1.46553625917463067508e-15,
    3.82059606377570462276e-16,   -1.00457896653436912508e-16
  };
  static constexpr G4double aa6           = 33.122936966163038145;
  static constexpr G4int ncheb4           = 27;
  static constexpr G4double cheb4[ncheb4] = {
    1.69342658227676741765,      0.0742766400841232319225,
    -0.019337880608635717358,    0.00516065527473364110491,
    -0.00139342012990307729473,  0.000378549864052022522193,
    -0.000103167085583785340215, 0.0000281543441271412178337,
    -7.68409742018258198651e-6,  2.09543221890204537392e-6,
    -5.70493140367526282946e-7,  1.54961164548564906446e-7,
    -4.19665599629607704794e-8,  1.13239680054166507038e-8,
    -3.04223563379021441863e-9,  8.13073745977562957997e-10,
    -2.15969415476814981374e-10, 5.69472105972525594811e-11,
    -1.48844799572430829499e-11, 3.84901514438304484973e-12,
    -9.82222575944247161834e-13, 2.46468329208292208183e-13,
    -6.04953826265982691612e-14, 1.44055805710671611984e-14,
    -3.28200813577388740722e-15, 6.96566359173765367675e-16,
    -1.294122794852896275e-16
  };

  if(x < aa2)
    return x * x * x * Chebyshev(aa1, aa2, cheb1, ncheb1, x);
  else if(x < aa3)
    return Chebyshev(aa2, aa3, cheb2, ncheb2, x);
  else if(x < 1 - 0.0000841363)
  {
    G4double y = -G4Log(1 - x);
    return y * Chebyshev(aa4, aa5, cheb3, ncheb3, y);
  }
  else
  {
    G4double y = -G4Log(1 - x);
    return y * Chebyshev(aa5, aa6, cheb4, ncheb4, y);
  }
}

G4double G4SynchrotronRadiation::GetRandomEnergySR(G4double gamma,
                                                   G4double perpB,
                                                   G4double mass_c2)
{
  static const G4double fEnergyConst =
    1.5 * c_light * c_light * eplus * hbar_Planck;
  G4double Ecr = fEnergyConst * gamma * gamma * perpB / mass_c2;

  if(verboseLevel > 0 && FirstTime1)
  {
    // mean and rms of photon energy
    G4double Emean = 8. / (15. * std::sqrt(3.)) * Ecr;
    G4double E_rms = std::sqrt(211. / 675.) * Ecr;
    G4long prec     = G4cout.precision();
    G4cout << "G4SynchrotronRadiation::GetRandomEnergySR :" << '\n'
           << std::setprecision(4) << "  Ecr   = " << G4BestUnit(Ecr, "Energy")
           << '\n'
           << "  Emean = " << G4BestUnit(Emean, "Energy") << '\n'
           << "  E_rms = " << G4BestUnit(E_rms, "Energy") << G4endl;
    FirstTime1 = false;
    G4cout.precision(prec);
  }

  G4double energySR = Ecr * InvSynFracInt(G4UniformRand());
  return energySR;
}

///////////////////////////////////////////////////////////////////////////////
void G4SynchrotronRadiation::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if(0 < verboseLevel && &part == G4Electron::Electron())
    ProcessDescription(G4cout);
  // same for all particles, print only for one (electron)
}

///////////////////////////////////////////////////////////////////////////////
void G4SynchrotronRadiation::ProcessDescription(std::ostream& out) const
{
  out << GetProcessName()
      << ":  Incoherent Synchrotron Radiation\n"
         "Good description for long magnets at all energies.\n";
}

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
//------------------ G4XrayReflection physics process -----------------------
//
// History:
// 14-09-23 H. Burkhardt, initial implementation
//
// --------------------------------------------------------------------------

#include "G4XrayReflection.hh"

#include "G4EmProcessSubType.hh"
#include "G4Exp.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4VDiscreteProcess.hh"

G4double G4XrayReflection::fSurfaceRoughness = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4XrayReflection::G4XrayReflection(const G4String& processName,
                                   G4ProcessType type)  //  Constructor
  : G4VDiscreteProcess(processName, type)
{
  SetProcessSubType(fGammaReflection);
  SaveHenkeDataAsMaterialProperty();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4XrayReflection::IsApplicable(const G4ParticleDefinition& particle)
{
  return (&particle == G4Gamma::Gamma());  // apply only to gamma
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4XrayReflection::Reflectivity(const G4double GamEner, const G4double SinIncidentAngle,
                                        const G4Material* theMat) const
{
  G4double theReflectivity = 0;
  const G4MaterialPropertiesTable* theMatProp = theMat->GetMaterialPropertiesTable();
  if (SinIncidentAngle < 0.9 && theMatProp != nullptr)
  {  // avoid perpendicular refl. at straight entry  and require
    // data available
    G4MaterialPropertyVector* RealIndex = theMatProp->GetProperty(kREALRINDEX);
    G4MaterialPropertyVector* ImagIndex = theMatProp->GetProperty(kIMAGINARYRINDEX);
    if (nullptr == RealIndex || nullptr == ImagIndex) { return theReflectivity; }
    const G4double delta = RealIndex->Value(GamEner);
    const G4double beta = ImagIndex->Value(GamEner);
    const G4double sin2 = std::pow(SinIncidentAngle, 2);
    const G4double rho2 =
      0.5 * (sin2 - 2 * delta + std::sqrt(std::pow(sin2 - 2 * delta, 2) + 4 * beta * beta));
    const G4double rho = std::sqrt(rho2);
    const G4double Refl_sigma = (rho2 * std::pow(SinIncidentAngle - rho, 2) + std::pow(beta, 2))
                                / (rho2 * std::pow(SinIncidentAngle + rho, 2) + std::pow(beta, 2));
    const G4double coscot = std::sqrt(1 - sin2) / SinIncidentAngle;
    const G4double pi_over_sigma = (rho2 * std::pow(rho - coscot, 2) + std::pow(beta, 2))
                                   / (rho2 * std::pow(rho + coscot, 2) + std::pow(beta, 2));
    const G4double Refl_pi = Refl_sigma * pi_over_sigma;
    theReflectivity = 0.5 * (Refl_sigma + Refl_pi);  // unpolarized
    G4double RoughAtten = 1;
    if (fSurfaceRoughness > 0) {
      G4double kiz = SinIncidentAngle * GamEner / CLHEP::hbarc;
      G4double kjz = SinIncidentAngle * (1 - delta) * GamEner / CLHEP::hbarc;
      RoughAtten = G4Exp(-2 * kiz * kjz * fSurfaceRoughness * fSurfaceRoughness);  // Nevot–Croce
      theReflectivity *= RoughAtten;
    }
    if (GetVerboseLevel() > 1)
      G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
             << std::right << std::setw(4) << __LINE__ << " GamEner=" << GamEner
             << " fSurfaceRoughness=" << G4BestUnit(fSurfaceRoughness, "Length")
             << " RoughAtten=" << RoughAtten << " SinIncidentAngle=" << SinIncidentAngle
             << " delta=" << delta << " beta=" << beta << " Refl_sigma=" << Refl_sigma
             << " Refl_pi=" << Refl_pi << " theReflectivity=" << theReflectivity << G4endl;
  }
  return theReflectivity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4XrayReflection::GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                                           G4ForceCondition* condition)
{
  *condition = NotForced;
  G4double GamEner = aTrack.GetDynamicParticle()->GetTotalEnergy();
  if (GamEner < 30. * eV || GamEner > 30. * keV)
    return DBL_MAX;  // do nothing below and above the limits

  if (GetVerboseLevel() > 2)
    G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
           << std::right << std::setw(4) << __LINE__ << " GamEner=" << GamEner / keV
           << " keV previousStepSize=" << previousStepSize
           << " TrackLength=" << aTrack.GetTrackLength() << " StepLength=" << aTrack.GetStepLength()
           << G4endl;

  G4double MeanFreePath = DBL_MAX;  // by default no reflection
  G4VPhysicalVolume* Volume = aTrack.GetVolume();
  if (fLastVolume && Volume != fLastVolume && aTrack.GetTrackLength() > 0) {  // at a boundary
    const G4Material* theLastMat = fLastVolume->GetLogicalVolume()->GetMaterial();
    const G4Material* theMat = Volume->GetLogicalVolume()->GetMaterial();

    G4double last_density = theLastMat->GetDensity();
    G4double density = theMat->GetDensity();
    if (density > last_density) {  // density has increased
      G4Navigator* theNavigator =
        G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
      G4bool valid = false;
      G4ThreeVector theSurfaceNormal =
        theNavigator->GetGlobalExitNormal(aTrack.GetPosition(), &valid);
      if (valid) fSurfaceNormal = theSurfaceNormal;
      G4double SinIncidentAngle =
        aTrack.GetDynamicParticle()->GetMomentumDirection() * fSurfaceNormal;
      if (G4UniformRand() < Reflectivity(GamEner, SinIncidentAngle, theMat)) {
        MeanFreePath = 0;
      }
      G4ThreeVector Position = aTrack.GetPosition();  // only for info
      const G4VSolid* LastSolid_Volume = fLastVolume->GetLogicalVolume()->GetSolid();  // for info
      if (GetVerboseLevel() > 1 && MeanFreePath == 0)
        G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
               << std::right << std::setw(4) << __LINE__
               << " trigger reflection SinIncidentAngle=" << SinIncidentAngle
               << " at z=" << Position.getZ() / meter << " m" << G4endl;
      else if (GetVerboseLevel() > 2)
        G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
               << std::right << std::setw(4) << __LINE__ << " volume has changed "
               << " last logical volume name =" << fLastVolume->GetLogicalVolume()->GetName()
               << " last logical volume material name =" << theLastMat->GetName()
               << " last density=" << last_density << " part/cm3 ? "
               << " logical volume name =" << Volume->GetLogicalVolume()->GetName()
               << " logical volume material name =" << theMat->GetName() << " density=" << density
               << " part/cm3 ? "
               << " LastSolid_Volume->Inside(Position)=" << LastSolid_Volume->Inside(Position)
               << " sin(IncidentAngle)=" << SinIncidentAngle << " MeanFreePath=" << MeanFreePath
               << G4endl;
    }
  }
  fLastVolume = Volume;
  return MeanFreePath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4XrayReflection::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);  // copy the current position to the changed particle
  G4ThreeVector PhotDir = aTrack.GetDynamicParticle()->GetMomentumDirection();
  G4ThreeVector para_part = (PhotDir * fSurfaceNormal) * fSurfaceNormal;
  G4ThreeVector photon_reflected = PhotDir - 2 * para_part;  // invert the parallel component
  if (GetVerboseLevel() > 1)
    G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
           << std::right << std::setw(4) << __LINE__ << " fSurfaceNormal=" << fSurfaceNormal
           << " StepLength=" << aStep.GetStepLength() << " PhotDir=" << PhotDir
           << " photon_reflected=" << photon_reflected << " para_part=" << para_part
           << " aParticleChange.GetTrackStatus()=" << aParticleChange.GetTrackStatus()
           << G4endl;

  aParticleChange.ProposeTrackStatus(
    fStopAndKill);  // needed when working with primary gamma to get rid of
  // primary
  auto ReflectedPhoton = new G4DynamicParticle(G4Gamma::Gamma(), photon_reflected,
                                               aTrack.GetDynamicParticle()->GetTotalEnergy());
  aParticleChange.AddSecondary(ReflectedPhoton);
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4XrayReflection::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  ProcessDescription(G4cout);

  if (GetVerboseLevel() > 2)
    G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
           << std::right << std::setw(4) << __LINE__
           << " is gamma=" << (&part == G4Gamma::Definition())
           << " fSurfaceRoughness=" << G4BestUnit(fSurfaceRoughness, "Length") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4XrayReflection::ProcessDescription(std::ostream& out) const
{
  if (G4Threading::IsMasterThread())
    out << '\n' << GetProcessName() << ": Gamma specular reflection for energies > 30 eV.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4XrayReflection::ReadHenkeXrayData(std::string ElName, std::vector<G4double>& Ephot,
                                          std::vector<G4double>& f1, std::vector<G4double>& f2)
{
  std::transform(ElName.begin(), ElName.end(), ElName.begin(),
                 ::tolower);  // henke_physical_reference uses lower case filanames
  const G4String DataDir = G4EmParameters::Instance()->GetDirLEDATA() + "/XRayReflection_data/";
  const G4String InpFname = DataDir + ElName + ".nff";
  std::ifstream infile(InpFname);
  if (!infile.is_open()) {
    G4cout << "ReadHenkeXrayReflData " << InpFname << " not found" << G4endl;
    return 1;  // failure
  }
  std::vector<std::string> VarName(3);
  infile >> VarName[0] >> VarName[1] >> VarName[2];
  if (GetVerboseLevel())
    G4cout << "ReadHenkeXrayData variable names " << VarName[0] << " " << VarName[1] << " "
           << VarName[2] << G4endl;
  G4double E_eV_i, f1_i, f2_i;
  Ephot.resize(0);
  f1.resize(0);
  f2.resize(0);
  for (;;) {
    infile >> E_eV_i >> f1_i >> f2_i;
    if (infile.eof()) break;
    Ephot.push_back(E_eV_i * eV);
    f1.push_back(f1_i);
    f2.push_back(f2_i);
  }
  return 0;  // success
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4XrayReflection::SaveHenkeDataAsMaterialProperty()
{
  // loop through the material table and load set up MaterialPropertiesTable
  // with Henke data used to calculate the reflection
  auto materialTable = G4Material::GetMaterialTable();
  for (auto a_material : *materialTable) {
    auto N = a_material->GetTotNbOfAtomsPerVolume();
    if (GetVerboseLevel() > 2)
      if (GetVerboseLevel() > 2)
        G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
               << std::right << std::setw(4) << __LINE__ << " " << a_material->GetName()
               << " NbOfAtomsPerVolume()=" << N
               << " NumberOfElements()=" << a_material->GetNumberOfElements() << G4endl;
    // calculate the reflectivity from input data. Implemented for dense
    // materials of a single element
    if (a_material->GetNumberOfElements() == 1 && a_material->GetDensity() > 1) {
      G4double factor = N * CLHEP::classic_electr_radius / CLHEP::twopi;
      std::vector<G4double> Ephot, f1, f2;
      const G4Element* theElement = a_material->GetElement(0);
      G4int iret = ReadHenkeXrayData(theElement->GetName(), Ephot, f1, f2);
      if (iret) {
        if (GetVerboseLevel() > 2)
          G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
                 << std::right << std::setw(4) << __LINE__ << " no Henke data found for "
                 << a_material->GetName() << " " << theElement->GetName() << G4endl;
      }
      else {
        std::vector<G4double> RealIndex(Ephot.size()), ImagIndex(Ephot.size());
        for (std::size_t i = 0; i < Ephot.size(); ++i) {
          G4double lambda = CLHEP::twopi * CLHEP::hbarc / Ephot[i];
          G4double lambda_sqr = lambda * lambda;
          RealIndex[i] = fmax(0, factor * lambda_sqr * f1[i]);  // delta or 1-RealIndex
          ImagIndex[i] = factor * lambda_sqr * f2[i];  // beta  or  -ImagIndex
          if (GetVerboseLevel() > 2)
            G4cout << "Ephot=" << std::setw(10) << Ephot[i] / eV << " eV delta=" << std::setw(10)
                   << RealIndex[i] << " beta=" << std::setw(10) << ImagIndex[i] << G4endl;
        }  // photon energy
        G4MaterialPropertiesTable* proptab = a_material->GetMaterialPropertiesTable();
        if(proptab == nullptr) {
          proptab = new G4MaterialPropertiesTable();
          a_material->SetMaterialPropertiesTable(proptab);
        }
        proptab->AddProperty("REALRINDEX", Ephot, RealIndex);  // 1-RealIndex
        proptab->AddProperty("IMAGINARYRINDEX", Ephot, ImagIndex);
        if (GetVerboseLevel() > 2)
          G4cout << std::left << std::setw(12) << __FILE__ << " " << __FUNCTION__ << " line "
                 << std::right << std::setw(4) << __LINE__ << " " << a_material->GetName()
                 << " " << theElement->GetName()
                 << " reflection data saved in PropertiesTable" << G4endl;
      }  // data found
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4XrayReflection::SetSurfaceRoughness(const G4double value)
{
  fSurfaceRoughness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

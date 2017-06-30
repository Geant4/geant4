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
// $Id: G4WentzelVIModel.cc 104802 2017-06-19 07:11:40Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4WentzelVIModel
//
// Author:      V.Ivanchenko 
//
// Creation date: 09.04.2008 from G4MuMscModel
//
// Modifications:
// 27-05-2010 V.Ivanchenko added G4WentzelOKandVIxSection class to
//              compute cross sections and sample scattering angle
//
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// G.Wentzel, Z. Phys. 40 (1927) 590.
// H.W.Lewis, Phys Rev 78 (1950) 526.
// J.M. Fernandez-Varea et al., NIM B73 (1993) 447.
// L.Urban, CERN-OPEN-2006-077.

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WentzelVIModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmParameters.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4WentzelVIModel::G4WentzelVIModel(G4bool comb, const G4String& nam) 
  : G4VMscModel(nam),
    ssFactor(1.05),
    invssFactor(1.0),
    currentCouple(nullptr),
    cosThetaMin(1.0),
    cosThetaMax(-1.0),
    fSecondMoments(nullptr),
    idx2(0),
    numlimit(0.1),
    singleScatteringMode(false),
    isCombined(comb),
    useSecondMoment(false)
{
  SetSingleScatteringFactor(1.25);
  invsqrt12 = 1./sqrt(12.);
  tlimitminfix = 1.e-6*mm;
  lowEnergyLimit = 1.0*eV;
  particle = 0;
  nelments = 5;
  xsecn.resize(nelments);
  prob.resize(nelments);
  wokvi = nullptr;
  fixedCut = -1.0;

  minNCollisions = 10;

  preKinEnergy = effKinEnergy = tPathLength = zPathLength = lambdaeff 
    = currentRange = xtsec = cosTetMaxNuc = 0.0;
  currentMaterialIndex = 0;

  fParticleChange = nullptr;
  currentCuts = nullptr;
  currentMaterial = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WentzelVIModel::~G4WentzelVIModel()
{
  delete wokvi;
  if(fSecondMoments && IsMaster()) {
    delete fSecondMoments;
    fSecondMoments = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::Initialise(const G4ParticleDefinition* p,
                                  const G4DataVector& cuts)
{
  if(!wokvi) { wokvi = new G4WentzelOKandVIxSection(isCombined); }

  // reset parameters
  SetupParticle(p);
  currentRange = 0.0;

  if(isCombined) {
    G4double tet = PolarAngleLimit();
    if(tet <= 0.0)           { cosThetaMax = 1.0; }
    else if(tet < CLHEP::pi) { cosThetaMax = cos(tet); }
  }
  //G4cout << "G4WentzelVIModel::Initialise " << p->GetParticleName() 
  //	 << " " << this << " " << wokvi << G4endl;

  wokvi->Initialise(p, cosThetaMax);
  /*  
  G4cout << "G4WentzelVIModel: " << particle->GetParticleName()
         << "  1-cos(ThetaLimit)= " << 1 - cosThetaMax 
         << " SingScatFactor= " << ssFactor
         << G4endl;
  */
  currentCuts = &cuts;

  // set values of some data members
  fParticleChange = GetParticleChangeForMSC(p);

  // build second moment table only if transport table is build
  G4PhysicsTable* table = GetCrossSectionTable();
  if(useSecondMoment && IsMaster() && table) {

    //G4cout << "### G4WentzelVIModel::Initialise: build 2nd moment table "
    //           << table << G4endl;
    fSecondMoments =  
      G4PhysicsTableHelper::PreparePhysicsTable(fSecondMoments);
    // Access to materials
    const G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();

    G4bool splineFlag = true;
    G4PhysicsVector* aVector = nullptr;
    G4PhysicsVector* bVector = nullptr;
    G4double emin = std::max(LowEnergyLimit(), LowEnergyActivationLimit());
    G4double emax = std::min(HighEnergyLimit(), HighEnergyActivationLimit());
    if(emin < emax) {
      size_t n = G4EmParameters::Instance()->NumberOfBinsPerDecade()
        *G4lrint(std::log10(emax/emin));
      if(n < 3) { n = 3; }

      for(size_t i=0; i<numOfCouples; ++i) {

        //G4cout<< "i= " << i << " Flag=  " << fSecondMoments->GetFlag(i) 
        //      << G4endl;
        if(fSecondMoments->GetFlag(i)) {
          DefineMaterial(theCoupleTable->GetMaterialCutsCouple(i));
       
          delete (*fSecondMoments)[i];
          if(!aVector) { 
            aVector = new G4PhysicsLogVector(emin, emax, n);
            bVector = aVector;
          } else {
            bVector = new G4PhysicsVector(*aVector);
          }
          for(size_t j=0; j<n; ++j) {
            G4double e = bVector->Energy(j); 
            bVector->PutValue(j, ComputeSecondMoment(p, e)*e*e);
          }
          if(splineFlag) { bVector->FillSecondDerivatives(); }
          (*fSecondMoments)[i] = bVector;  
        }
      }
    } 
    //G4cout << *fSecondMoments << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::InitialiseLocal(const G4ParticleDefinition*, 
                                       G4VEmModel* masterModel)
{
  fSecondMoments = static_cast<G4WentzelVIModel*>(masterModel)
    ->GetSecondMomentTable(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4WentzelVIModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
    SetCurrentCouple(cup); 
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex(); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* p,
                             G4double kinEnergy,
                             G4double Z, G4double,
                             G4double cutEnergy, G4double)
{
  G4double cross = 0.0;
  if(p != particle) { SetupParticle(p); }
  if(kinEnergy < lowEnergyLimit) { return cross; }
  if(!CurrentCouple()) {
    G4Exception("G4WentzelVIModel::ComputeCrossSectionPerAtom", "em0011",
                FatalException, " G4MaterialCutsCouple is not defined");
    return 0.0;
  }
  DefineMaterial(CurrentCouple());
  cosTetMaxNuc = wokvi->SetupKinematic(kinEnergy, currentMaterial);
  if(cosTetMaxNuc < 1.0) {
    G4double cut  = (0.0 < fixedCut) ? fixedCut : cutEnergy;
    G4double cost = wokvi->SetupTarget(G4lrint(Z), cut);
    cross = wokvi->ComputeTransportCrossSectionPerAtom(cost);
    /*
    if(p->GetParticleName() == "e-")      
    G4cout << "G4WentzelVIModel::CS: Z= " << G4int(Z) << " e(MeV)= "<<kinEnergy 
           << " 1-cosN= " << 1 - cosTetMaxNuc << " cross(bn)= " << cross/barn
           << " " << particle->GetParticleName() << G4endl;
    */
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::StartTracking(G4Track* track)
{
  SetupParticle(track->GetDynamicParticle()->GetDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeTruePathLengthLimit(
                             const G4Track& track,
                             G4double& currentMinimalStep)
{
  G4double tlimit = currentMinimalStep;
  const G4DynamicParticle* dp = track.GetDynamicParticle();
  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus = sp->GetStepStatus();
  singleScatteringMode = false;

  //G4cout << "G4WentzelVIModel::ComputeTruePathLengthLimit stepStatus= " 
  //         << stepStatus << "  " << track.GetDefinition()->GetParticleName() 
  //         << G4endl;

  // initialisation for each step, lambda may be computed from scratch
  preKinEnergy = dp->GetKineticEnergy();
  effKinEnergy = preKinEnergy;
  DefineMaterial(track.GetMaterialCutsCouple());
  lambdaeff = GetTransportMeanFreePath(particle,preKinEnergy);
  currentRange = GetRange(particle,preKinEnergy,currentCouple);
  cosTetMaxNuc = wokvi->SetupKinematic(preKinEnergy, currentMaterial);
  
  //G4cout << "lambdaeff= " << lambdaeff << " Range= " << currentRange
  // << " tlimit= " << tlimit << " 1-cost= " << 1 - cosTetMaxNuc << G4endl;
  
  // extra check for abnormal situation
  // this check needed to run MSC with eIoni and eBrem inactivated
  if(tlimit > currentRange) { tlimit = currentRange; }

  // stop here if small range particle
  if(tlimit < tlimitminfix) { 
    return ConvertTrueToGeom(tlimit, currentMinimalStep); 
  }

  // pre step
  G4double presafety = sp->GetSafety();
  // far from geometry boundary
  if(currentRange < presafety) {
    return ConvertTrueToGeom(tlimit, currentMinimalStep);
  }

  // compute presafety again if presafety <= 0 and no boundary
  // i.e. when it is needed for optimization purposes
  if(stepStatus != fGeomBoundary && presafety < tlimitminfix) {
    presafety = ComputeSafety(sp->GetPosition(), tlimit); 
    if(currentRange < presafety) {
      return ConvertTrueToGeom(tlimit, currentMinimalStep);
    }
  }
  /*       
  G4cout << "e(MeV)= " << preKinEnergy/MeV
         << "  " << particle->GetParticleName() 
         << " CurLimit(mm)= " << tlimit/mm <<" safety(mm)= " << presafety/mm
         << " R(mm)= " <<currentRange/mm
         << " L0(mm^-1)= " << lambdaeff*mm 
         << G4endl;
  */
  // natural limit for high energy
  G4double rlimit = std::max(facrange*currentRange, 
                             (1.0 - cosTetMaxNuc)*lambdaeff*invssFactor);

  // low-energy e-
  if(cosThetaMax > cosTetMaxNuc) {
    rlimit = std::min(rlimit, facsafety*presafety);
  }
   
  // cut correction
  G4double rcut = currentCouple->GetProductionCuts()->GetProductionCut(1);
  //G4cout << "rcut= " << rcut << " rlimit= " << rlimit << " presafety= " 
  // << presafety << " 1-cosThetaMax= " <<1-cosThetaMax 
  //<< " 1-cosTetMaxNuc= " << 1-cosTetMaxNuc << G4endl;
  if(rcut > rlimit) { rlimit = std::min(rlimit, rcut*sqrt(rlimit/rcut)); }
 
  tlimit = std::min(tlimit, rlimit);
  tlimit = std::max(tlimit, tlimitminfix);

  // step limit in infinite media
  tlimit = std::min(tlimit, 50*currentMaterial->GetRadlen()/facgeom);

  //compute geomlimit and force few steps within a volume
  if (steppingAlgorithm == fUseDistanceToBoundary 
      && stepStatus == fGeomBoundary) {

    G4double geomlimit = ComputeGeomLimit(track, presafety, currentRange);
    tlimit = std::min(tlimit, geomlimit/facgeom);
  } 
  /*         
  G4cout << particle->GetParticleName() << " e= " << preKinEnergy
         << " L0= " << lambdaeff << " R= " << currentRange
         << " tlimit= " << tlimit  
           << " currentMinimalStep= " << currentMinimalStep << G4endl;
  */
  return ConvertTrueToGeom(tlimit, currentMinimalStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeGeomPathLength(G4double truelength)
{
  zPathLength = tPathLength = truelength;

  // small step use only single scattering
  cosThetaMin = 1.0;
  ComputeTransportXSectionPerVolume(cosThetaMin);
  //G4cout << "xtsec= " << xtsec << "  Nav= " 
  //         << zPathLength*xtsec << G4endl;
  if(0.0 >= lambdaeff || G4int(zPathLength*xtsec) < minNCollisions) {
    singleScatteringMode = true;
    lambdaeff = DBL_MAX;

  } else {
    //G4cout << "ComputeGeomPathLength: tLength= " << tPathLength
    //           << " Leff= " << lambdaeff << G4endl; 
    // small step
    if(tPathLength < numlimit*lambdaeff) {
      G4double tau = tPathLength/lambdaeff;
      zPathLength *= (1.0 - 0.5*tau + tau*tau/6.0);

      // medium step
    } else {
      G4double e1 = 0.0;
      if(currentRange > tPathLength) {
        e1 = GetEnergy(particle,currentRange-tPathLength,currentCouple);
      }
      effKinEnergy = 0.5*(e1 + preKinEnergy);
      cosTetMaxNuc = wokvi->SetupKinematic(effKinEnergy, currentMaterial);
      lambdaeff = GetTransportMeanFreePath(particle, effKinEnergy);
      //G4cout << " tLength= "<< tPathLength<< " Leff= " << lambdaeff << G4endl;
      zPathLength = lambdaeff;
      if(tPathLength*numlimit < lambdaeff) {
        zPathLength *= (1.0 - G4Exp(-tPathLength/lambdaeff));
      }
    }
  }
  //G4cout << "Comp.geom: zLength= "<<zPathLength<<" tLength= "
  //         << tPathLength<< " Leff= " << lambdaeff << G4endl;
  return zPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeTrueStepLength(G4double geomStepLength)
{
  // initialisation of single scattering x-section
  /*
  G4cout << "ComputeTrueStepLength: Step= " << geomStepLength 
         << "  geomL= " << zPathLength
         << "  Lambda= " <<  lambdaeff 
           << " 1-cosThetaMaxNuc= " << 1 - cosTetMaxNuc << G4endl;
  */
  if(singleScatteringMode) {
    zPathLength = tPathLength = geomStepLength;

  } else {

    // step defined by transportation
    // change both geom and true step lengths 
    if(geomStepLength < zPathLength) { 

      // single scattering
      if(G4int(geomStepLength*xtsec) < minNCollisions) {
        zPathLength = tPathLength = geomStepLength;
        lambdaeff = DBL_MAX;
        singleScatteringMode = true;

        // multiple scattering
      } else {
        // small step
        if(geomStepLength < numlimit*lambdaeff) {
          G4double tau = geomStepLength/lambdaeff;
          tPathLength = geomStepLength*(1.0 + 0.5*tau + tau*tau/3.0); 

          // energy correction for a big step
        } else {
          tPathLength *= geomStepLength/zPathLength;
          G4double e1 = 0.0;
          if(currentRange > tPathLength) {
            e1 = GetEnergy(particle,currentRange-tPathLength,currentCouple);
          }
          effKinEnergy = 0.5*(e1 + preKinEnergy);
          cosTetMaxNuc = wokvi->SetupKinematic(effKinEnergy, currentMaterial);
          lambdaeff = GetTransportMeanFreePath(particle, effKinEnergy);
          G4double tau = geomStepLength/lambdaeff;

          if(tau < 0.999999) { tPathLength = -lambdaeff*G4Log(1.0 - tau); } 
          else               { tPathLength = currentRange; }
        }
        zPathLength = geomStepLength;
      }
    }
  }
  // check of step length
  // define threshold angle between single and multiple scattering 
  if(!singleScatteringMode) {
    cosThetaMin -= ssFactor*tPathLength/lambdaeff; 
    xtsec = 0.0;

    // recompute transport cross section - do not change energy
    // anymore - cannot be applied for big steps
    if(cosThetaMin > cosTetMaxNuc) {
      // new computation
      G4double cross = ComputeTransportXSectionPerVolume(cosThetaMin);
      //G4cout << "%%%% cross= " << cross << "  xtsec= " << xtsec 
      //           << " 1-cosTMin= " << 1.0 - cosThetaMin << G4endl;
      if(cross <= 0.0) {
        singleScatteringMode = true;
        tPathLength = zPathLength; 
        lambdaeff = DBL_MAX;
        cosThetaMin = 1.0;
      } else if(xtsec > 0.0) {
        
        lambdaeff = 1./cross; 
        G4double tau = zPathLength*cross;
        if(tau < numlimit) { 
          tPathLength = zPathLength*(1.0 + 0.5*tau + tau*tau/3.0); 
        } else if(tau < 0.999999) { 
          tPathLength = -lambdaeff*G4Log(1.0 - tau); 
        } else { 
          tPathLength = currentRange;
        }
      }
    } 
  }
  tPathLength = std::min(tPathLength, currentRange);
  /*      
  G4cout <<"Comp.true: zLength= "<<zPathLength<<" tLength= "<<tPathLength
         <<" Leff(mm)= "<<lambdaeff/mm<<" sig0(1/mm)= " << xtsec <<G4endl;
  G4cout << particle->GetParticleName() << " 1-cosThetaMin= " << 1-cosThetaMin
         << " 1-cosTetMaxNuc= " << 1-cosTetMaxNuc 
         << " e(MeV)= " << preKinEnergy/MeV << "  "  
         << " SSmode= " << singleScatteringMode << G4endl;
  */
  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector& 
G4WentzelVIModel::SampleScattering(const G4ThreeVector& oldDirection,
                                   G4double /*safety*/)
{
  fDisplacement.set(0.0,0.0,0.0);
  //G4cout << "!##! G4WentzelVIModel::SampleScattering for " 
  //         << particle->GetParticleName() << G4endl;

  // ignore scattering for zero step length and energy below the limit
  if(preKinEnergy < lowEnergyLimit || tPathLength <= 0.0) 
    { return fDisplacement; }
  
  G4double invlambda = 0.0;
  if(lambdaeff < DBL_MAX) { invlambda = 0.5/lambdaeff; }

  // use average kinetic energy over the step
  G4double cut = (*currentCuts)[currentMaterialIndex];
  if(fixedCut > 0.0) { cut = fixedCut; }
  /*  
  G4cout <<"SampleScat: E0(MeV)= "<< preKinEnergy/MeV
           << " Leff= " << lambdaeff <<" sig0(1/mm)= " << xtsec 
          << " xmsc= " <<  tPathLength*invlambda 
         << " safety= " << safety << G4endl;
  */
  // step limit due msc
  G4int nMscSteps = 1;
  G4double x0 = tPathLength;
  G4double z0 = x0*invlambda;
  //G4double zzz = 0.0;
  G4double prob2 = 0.0;

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  // large scattering angle case - two step approach
  if(!singleScatteringMode) {
    static const G4double zzmin = 0.05;
    if(useSecondMoment) { 
      G4double z1 = invlambda*invlambda;
      G4double z2 = SecondMoment(particle, currentCouple, effKinEnergy);
      prob2 = (z2 - z1)/(1.5*z1 - z2);
    }
    //    if(z0 > zzmin && safety > tlimitminfix) { 
    if(z0 > zzmin) { 
      x0 *= 0.5; 
      z0 *= 0.5;
      nMscSteps = 2;
    } 
    //if(z0 > zzmin) { zzz = G4Exp(-1.0/z0); }
    G4double zzz = 0.0;
    if(z0 > zzmin) { 
      zzz = G4Exp(-1.0/z0); 
      z0 += zzz; 
      prob2 *= (1 + zzz);
    }
    prob2 /= (1 + prob2);
  } 

  // step limit due to single scattering
  G4double x1 = 2*tPathLength;
  if(0.0 < xtsec) { x1 = -G4Log(rndmEngine->flat())/xtsec; }

  // no scattering case
  if(singleScatteringMode && x1 > tPathLength)  
    { return fDisplacement; }

  const G4ElementVector* theElementVector = 
    currentMaterial->GetElementVector();
  G4int nelm = currentMaterial->GetNumberOfElements();

  // geometry
  G4double sint, cost, phi;
  G4ThreeVector temp(0.0,0.0,1.0);

  // current position and direction relative to the end point
  // because of magnetic field geometry is computed relatively to the 
  // end point of the step 
  G4ThreeVector dir(0.0,0.0,1.0);
  fDisplacement.set(0.0,0.0,-zPathLength);

  G4double mscfac = zPathLength/tPathLength;

  // start a loop 
  G4double x2 = x0;
  G4double step, z;
  G4bool singleScat;
  /*   
    G4cout << "Start of the loop x1(mm)= " << x1 << "  x2(mm)= " << x2 
    << " 1-cost1= " << 1 - cosThetaMin << " SSmode= " << singleScatteringMode 
           << " xtsec= " << xtsec << " Nst= "  << nMscSteps << G4endl;
  */
  do {

    //G4cout << "# x1(mm)= "<< x1<< " x2(mm)= "<< x2 << G4endl;
    // single scattering case
    if(singleScatteringMode && x1 > x2) { 
      fDisplacement += x2*mscfac*dir;
      break; 
    }

    // what is next single of multiple?
    if(x1 <= x2) { 
      step = x1;
      singleScat = true;
    } else {
      step = x2;
      singleScat = false;
    }

    //G4cout << "# step(mm)= "<< step<< "  singlScat= "<< singleScat << G4endl;

    // new position
    fDisplacement += step*mscfac*dir;

    if(singleScat) {

      // select element
      G4int i = 0;
      if(nelm > 1) {
        G4double qsec = rndmEngine->flat()*xtsec;
        for (; i<nelm; ++i) { if(xsecn[i] >= qsec) { break; } }
      }
      G4double cosTetM = 
        wokvi->SetupTarget((*theElementVector)[i]->GetZasInt(), cut);
      //G4cout << "!!! " << cosThetaMin << "  " << cosTetM << "  " 
      //     << prob[i] << G4endl;
      temp = wokvi->SampleSingleScattering(cosThetaMin, cosTetM, prob[i]);

      // direction is changed
      temp.rotateUz(dir);
      dir = temp;
      //G4cout << dir << G4endl;

      // new proposed step length
      x2 -= step; 
      x1  = -G4Log(rndmEngine->flat())/xtsec; 

    // multiple scattering
    } else { 
      --nMscSteps;
      x1 -= step;
      x2  = x0;

      // sample z in interval 0 - 1
      G4bool isFirst = true;
      if(prob2 > 0.0 && rndmEngine->flat() < prob2) { isFirst = false; } 
      do {
        //z = -z0*G4Log(1.0 - (1.0 - zzz)*rndmEngine->flat());
        if(isFirst) { z = -G4Log(rndmEngine->flat()); }
        else        { z = G4RandGamma::shoot(rndmEngine, 2.0, 2.0); }
        z *= z0;
        // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
      } while(z > 1.0);

      cost = 1.0 - 2.0*z/*factCM*/;
      if(cost > 1.0)       { cost = 1.0; }
      else if(cost < -1.0) { cost =-1.0; }
      sint = sqrt((1.0 - cost)*(1.0 + cost));
      phi  = twopi*rndmEngine->flat();
      G4double vx1 = sint*cos(phi);
      G4double vy1 = sint*sin(phi);

      // lateral displacement  
      if (latDisplasment) {
        G4double rms = invsqrt12*sqrt(2*z0);
        G4double r  = x0*mscfac;
        G4double dx = r*(0.5*vx1 + rms*G4RandGauss::shoot(rndmEngine,0.0,1.0));
        G4double dy = r*(0.5*vy1 + rms*G4RandGauss::shoot(rndmEngine,0.0,1.0));
        G4double d  = r*r - dx*dx - dy*dy;

        // change position
        if(d >= 0.0)  { 
          temp.set(dx,dy,sqrt(d) - r);
          temp.rotateUz(dir); 
          fDisplacement += temp;
        }
      }
      // change direction
      temp.set(vx1,vy1,cost);
      temp.rotateUz(dir);
      dir = temp;
    }
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (0 < nMscSteps);
    
  dir.rotateUz(oldDirection);

  //G4cout<<"G4WentzelVIModel sampling is done 1-cost= "<< 1.-dir.z()<<G4endl;
  // end of sampling -------------------------------

  fParticleChange->ProposeMomentumDirection(dir);

  // lateral displacement  
  fDisplacement.rotateUz(oldDirection);

  /*
         G4cout << " r(mm)= " << fDisplacement.mag() 
                << " safety= " << safety
                << " trueStep(mm)= " << tPathLength
                << " geomStep(mm)= " << zPathLength
                << " x= " << fDisplacement.x() 
                << " y= " << fDisplacement.y() 
                << " z= " << fDisplacement.z()
                << G4endl;
  */

  //G4cout<< "G4WentzelVIModel::SampleScattering end NewDir= " << dir<< G4endl;
  return fDisplacement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeTransportXSectionPerVolume(G4double cosTheta)
{
  // prepare recomputation of x-sections
  const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    currentMaterial->GetVecNbOfAtomsPerVolume();
  G4int nelm = currentMaterial->GetNumberOfElements();
  if(nelm > nelments) {
    nelments = nelm;
    xsecn.resize(nelm);
    prob.resize(nelm);
  }

  // check consistency
  xtsec = 0.0;
  if(cosTetMaxNuc >= cosTheta) { return 0.0; }

  G4double cut = (*currentCuts)[currentMaterialIndex];
  if(fixedCut > 0.0) { cut = fixedCut; }

  // loop over elements
  G4double xs = 0.0;
  for (G4int i=0; i<nelm; ++i) {
    G4double costm = 
      wokvi->SetupTarget((*theElementVector)[i]->GetZasInt(), cut);
    G4double density = theAtomNumDensityVector[i];

    G4double esec = 0.0;
    if(costm < cosTheta) {  

      // recompute the transport x-section
      if(1.0 > cosTheta) {
        xs += density*wokvi->ComputeTransportCrossSectionPerAtom(cosTheta);
      }
      // recompute the total x-section
      G4double nucsec = wokvi->ComputeNuclearCrossSection(cosTheta, costm);
      esec = wokvi->ComputeElectronCrossSection(cosTheta, costm);
      nucsec += esec;
      if(nucsec > 0.0) { esec /= nucsec; }
      xtsec += nucsec*density;
    }
    xsecn[i] = xtsec;
    prob[i]  = esec;
    //G4cout << i << "  xs= " << xs << " xtsec= " << xtsec 
    //       << " 1-cosTheta= " << 1-cosTheta 
    //           << " 1-cosTetMaxNuc2= " <<1-cosTetMaxNuc2<< G4endl;
  }
  
  //G4cout << "ComputeXS result:  xsec(1/mm)= " << xs 
  //         << " txsec(1/mm)= " << xtsec <<G4endl; 
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel:: ComputeSecondMoment(const G4ParticleDefinition* p,
                                                G4double kinEnergy)
{
  G4double xs = 0.0;

  SetupParticle(p);
  cosTetMaxNuc = wokvi->SetupKinematic(kinEnergy, currentMaterial);

  if(cosTetMaxNuc >= 1.0) { return xs; }

  const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    currentMaterial->GetVecNbOfAtomsPerVolume();
  G4int nelm = currentMaterial->GetNumberOfElements();

  G4double cut = (*currentCuts)[currentMaterialIndex];
  if(fixedCut > 0.0) { cut = fixedCut; }

  // loop over elements
  for (G4int i=0; i<nelm; ++i) {
    G4double costm = 
      wokvi->SetupTarget((*theElementVector)[i]->GetZasInt(), cut);
    xs += theAtomNumDensityVector[i]
      *wokvi->ComputeSecondTransportMoment(costm);
  }
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4WentzelVIModel::SetSingleScatteringFactor(G4double val)
{
  if(val > 0.05) {
    ssFactor = val;
    invssFactor = 1.0/(val - 0.05);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

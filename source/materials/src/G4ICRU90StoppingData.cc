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
//
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:    G4ICRU90StoppingData
//
// Description:  Data on electroninc stopping power from ICRU 90
//
// Author:       Lucas Norberto Burigo
//
// Creation date: 03.09.2018
//
// Modifications: 25.09.2018 V.Ivanchenko adopted for material sub-library
// 
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ICRU90StoppingData.hh" 

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ICRU90StoppingData::G4ICRU90StoppingData() : isInitialized(false)
{
  // 1st initialisation 
  for(size_t i=0; i<nvectors; ++i) { 
    materials[i]    = nullptr; 
    sdata_proton[i] = nullptr; 
    sdata_alpha[i]  = nullptr;
  }
  FillData();

  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ICRU90StoppingData::~G4ICRU90StoppingData()
{
  for(size_t i=0; i<nvectors; ++i) { 
    delete sdata_proton[i]; 
    delete sdata_alpha[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4ICRU90StoppingData::Initialise()
{
  if(isInitialized) { return; }
  // this method may be called several times during initialisation
  G4int nmat = (G4int)G4Material::GetNumberOfMaterials();
  if(nmat == (G4int)nvectors) { return; }

  static const G4String nameNIST_ICRU90[3] = 
    {"G4_AIR","G4_WATER","G4_GRAPHITE"};

  // loop via material list to add extra data
  for(G4int i=0; i<nmat; ++i) {
    const G4Material* mat = (*(G4Material::GetMaterialTable()))[i];

    G4bool isThere = false;
    for(auto& material : materials)
    {
      if(mat == material)
      {
        isThere = true;
        break;
      }
    }
    if(!isThere) {
      // check list of NIST materials
      G4String mname = mat->GetName();
      for(G4int j=0; j<nvectors; ++j) {
        if(mname == nameNIST_ICRU90[j]) {
          materials[j] = mat;
          break;
	}
      }
    }
    isInitialized = ((materials[0] != nullptr) && (materials[1] != nullptr) &&
                     (materials[2] != nullptr));
    if(isInitialized) { return; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU90StoppingData::GetElectronicDEDXforProton(
         const G4Material* mat, G4double kinEnergy) const
{
  G4int idx = GetIndex(mat);
  return (idx < 0) ? 0.0 : GetDEDX(sdata_proton[idx], kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU90StoppingData::GetElectronicDEDXforAlpha(
         const G4Material* mat, G4double scaledKinEnergy) const
{
  G4int idx = GetIndex(mat);
  return (idx < 0) ? 0.0 : GetDEDX(sdata_alpha[idx], scaledKinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ICRU90StoppingData::FillData() 
{
  G4double T0_proton[57] = {  0.0010, 0.001500, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0080, 0.010, 0.0150, 0.020, 0.030, 0.040, 0.050, 0.060, 0.080, 0.10, 0.150, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.50, 2.00, 3.00, 4.00, 5.00, 6.00, 8.00, 10.00, 15.00, 20.00, 30.00, 40.00, 50.00, 60.00, 80.00, 100.00, 150.00, 200.00, 300.00, 400.00, 500.00, 600.00, 800.00, 1000.00, 1500.00, 2000.00, 3000.00, 4000.00, 5000.00, 6000.00, 8000.00, 10000.00  };

  G4double T0_alpha[49] = {  0.0010, 0.001500, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0080, 0.010, 0.0150, 0.020, 0.030, 0.040, 0.050, 0.060, 0.080, 0.10, 0.150, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.50, 2.00, 3.00, 4.00, 5.00, 6.00, 8.00, 10.00, 15.00, 20.00, 30.00, 40.00, 50.00, 60.00, 80.00, 100.00, 150.00, 200.00, 300.00, 400.00, 500.00, 600.00, 800.00, 1000.00};

  static const G4float e0_proton[57] = {  119.70f, 146.70f, 169.30f, 207.40f, 239.50f, 267.80f, 293.30f, 338.70f, 378.70f, 450.40f, 506.70f, 590.50f, 648.30f, 687.70f, 713.20f, 734.10f, 729.00f, 667.20f, 592.20f, 476.50f, 401.20f, 349.80f, 312.10f, 258.70f, 222.70f, 168.20f, 137.00f, 101.70f, 81.920f, 69.050f, 59.940f, 47.810f, 40.040f, 28.920f, 22.930f, 16.520f, 13.120f, 10.980f, 9.5140f, 7.6180f, 6.4410f, 4.8150f, 3.9750f, 3.1170f, 2.6860f, 2.4310f, 2.2650f, 2.0690f, 1.9620f, 1.850f, 1.820f, 1.8280f, 1.8610f, 1.8980f, 1.9340f, 1.9980f, 2.0520f  };

  static const G4float e0_alpha[49] = {  87.50f, 108.60f, 126.70f, 157.30f, 183.50f, 206.70f, 227.90f, 265.90f, 299.60f, 372.30f, 434.30f, 539.50f, 629.00f, 708.10f, 779.80f, 906.80f, 1018.00f, 1247.00f, 1429.00f, 1693.00f, 1861.00f, 1961.00f, 2008.00f, 2002.00f, 1922.00f, 1626.00f, 1381.00f, 1071.00f, 885.80f, 760.60f, 669.60f, 545.30f, 463.40f, 342.30f, 274.70f, 200.10f, 159.30f, 133.20f, 115.00f, 91.190f, 76.140f, 54.930f, 43.690f, 31.840f, 25.630f, 21.790f, 19.170f, 15.830f, 13.790f  };

  static const G4float e1_proton[57] = {  133.70f, 163.80f, 189.10f, 231.60f, 267.50f, 299.00f, 327.60f, 378.20f, 422.90f, 503.60f, 567.30f, 662.80f, 729.00f, 774.00f, 802.60f, 824.10f, 814.50f, 736.00f, 658.50f, 543.50f, 464.30f, 406.50f, 362.40f, 299.70f, 257.40f, 193.40f, 156.90f, 116.00f, 93.190f, 78.420f, 68.010f, 54.170f, 45.320f, 32.690f, 25.890f, 18.640f, 14.790f, 12.380f, 10.720f, 8.5780f, 7.250f, 5.4170f, 4.470f, 3.5040f, 3.0180f, 2.7310f, 2.5440f, 2.3230f, 2.2030f, 2.0650f, 2.0170f, 1.9980f, 2.010f, 2.0290f, 2.050f, 2.090f, 2.1240f  };

  static const G4float e1_alpha[49] = {  98.910f, 122.80f, 143.10f, 177.60f, 206.90f, 233.00f, 256.80f, 299.30f, 337.00f, 418.10f, 487.10f, 603.90f, 703.00f, 790.50f, 869.60f, 1009.00f, 1131.00f, 1383.00f, 1582.00f, 1873.00f, 2062.00f, 2178.00f, 2240.00f, 2256.00f, 2190.00f, 1877.00f, 1599.00f, 1239.00f, 1021.00f, 874.70f, 768.70f, 623.90f, 529.00f, 389.40f, 311.80f, 226.70f, 180.20f, 150.60f, 130.00f, 103.00f, 85.930f, 61.940f, 49.230f, 35.860f, 28.850f, 24.520f, 21.560f, 17.80f, 15.50f  };

  static const G4float e2_proton[57] = {  118.50f, 145.10f, 167.60f, 205.30f, 237.00f, 265.00f, 290.30f, 335.20f, 374.80f, 435.30f, 481.10f, 550.30f, 611.00f, 663.00f, 697.60f, 726.70f, 726.60f, 671.20f, 598.60f, 483.30f, 407.10f, 354.60f, 315.90f, 262.40f, 226.10f, 171.10f, 139.40f, 103.50f, 83.260f, 70.10f, 60.810f, 48.450f, 40.550f, 29.250f, 23.170f, 16.680f, 13.230f, 11.070f, 9.5910f, 7.6750f, 6.4860f, 4.8430f, 3.9940f, 3.1250f, 2.6870f, 2.4260f, 2.2550f, 2.050f, 1.9360f, 1.8060f, 1.760f, 1.7440f, 1.7560f, 1.7750f, 1.7960f, 1.8340f, 1.8680f  };

  static const G4float e2_alpha[49] = {  192.30f, 228.90f, 259.00f, 308.30f, 348.90f, 384.00f, 415.30f, 469.90f, 517.00f, 615.10f, 695.50f, 826.20f, 932.70f, 1024.00f, 1104.00f, 1240.00f, 1354.00f, 1574.00f, 1731.00f, 1929.00f, 2027.00f, 2063.00f, 2060.00f, 1993.00f, 1891.00f, 1615.00f, 1387.00f, 1085.00f, 898.20f, 772.00f, 680.40f, 554.40f, 471.20f, 347.80f, 278.80f, 202.80f, 161.20f, 134.80f, 116.30f, 92.140f, 76.890f, 55.430f, 44.050f, 32.090f, 25.810f, 21.930f, 19.280f, 15.910f, 13.840f  };

  sdata_proton[0] = AddData(57, T0_proton, e0_proton);
  sdata_proton[1] = AddData(57, T0_proton, e1_proton);
  sdata_proton[2] = AddData(57, T0_proton, e2_proton);

  sdata_alpha[0] = AddData(49, T0_alpha, e0_alpha);
  sdata_alpha[1] = AddData(49, T0_alpha, e1_alpha);
  sdata_alpha[2] = AddData(49, T0_alpha, e2_alpha);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsFreeVector* G4ICRU90StoppingData::AddData(G4int n, const G4double* e, 
						   const G4float* dedx)
{
  static const G4double fac = CLHEP::MeV*CLHEP::cm2/CLHEP::g;

  auto* data = new G4PhysicsFreeVector(n, e[0], e[n - 1], true);
  for(G4int i=0; i<n; ++i) { 
    data->PutValues(i, e[i]*CLHEP::MeV, ((G4double)dedx[i])*fac); 
  }
  data->FillSecondDerivatives();
  return data;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

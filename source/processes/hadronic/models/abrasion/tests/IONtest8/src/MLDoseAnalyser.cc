////////////////////////////////////////////////////////////////////////////////
//
#include "MLDoseAnalyser.hh"
#include <iomanip>

#include "G4RunManager.hh"
#ifndef USEHBOOK
#include "MLVersion.hh"
#endif

#include "G4UnitsTable.hh"
#include <iomanip>
////////////////////////////////////////////////////////////////////////////////
//
MLDoseAnalyser::MLDoseAnalyser (MLAnalysisManager* anMan,
  MLGeometryConstruction * det)
{
//
//
// Set default conditions.
//
  doseUnit  = "rad";
  Nb        = -9999;
  factor.clear();
  factorSqr.clear();
  ADose.clear();
  ADoseSqr.clear();
  dose.clear();
  doseErr.clear();
  doseLayers.clear();
//
//
// Set Analysis Manager.
//
  analysisManager = anMan;
//
//
// Get the geometry.
//
  geometry         = det;
  materialsManager = geometry->GetMLMaterial();
}
////////////////////////////////////////////////////////////////////////////////
//
MLDoseAnalyser::~MLDoseAnalyser ()
{}
////////////////////////////////////////////////////////////////////////////////
//
void MLDoseAnalyser::AddDoseLayer (G4int i)
{
  G4int NbOfLayers = geometry->GetNbOfLayers();
  if (i <= 0 || i > NbOfLayers) {
    G4cerr <<G4endl;
    G4cerr <<"AddDoseLayer : Layer number " <<i
           <<" must be > 0 and <= " <<NbOfLayers <<"." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  }else if (std::find(doseLayers.begin(),doseLayers.end(),i-1)
    != doseLayers.end()) {
    G4cerr <<G4endl;
    G4cerr <<"AddDOseLayer : Layer number " <<i << " already selected."
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
    
  } else {
    doseLayers.push_back(i-1);
    
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLDoseAnalyser::DeleteDoseLayer (G4int i)
{
  // search the material by its name
  //
  G4int NbOfLayers = geometry->GetNbOfLayers();
  if (i < 0 || i > NbOfLayers) {
    G4cerr <<G4endl;
    G4cerr <<"DeleteDoseLayer: layer number " <<i
           <<" out of range." <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else  if ( i == 0) {
    doseLayers.clear();

  } else if (std::find(doseLayers.begin(),doseLayers.end(),i-1)
    == doseLayers.end()) {
    G4cerr <<G4endl;
    G4cerr <<"DeleteDoseLayer : " <<i <<" layer has not been selected!"
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    doseLayers.erase(std::find(doseLayers.begin(),doseLayers.end(),i-1));

  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLDoseAnalyser::ListDoseLayer ()
{
  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
  G4cout <<"Dose analysis will be carried out for the following layers:"
         <<G4endl;

  if (doseLayers.size() == 0) {
    G4cout <<G4endl;
    G4cout <<"No layer has been selected for ionising dose analysis."
           <<G4endl;
  } else {
    for (size_t k=0; k<doseLayers.size(); k++) {
      G4int i = doseLayers[k];
      G4cout <<"  Layer No. " <<std::setw(3) <<i+1 <<" "
             <<std::setw(20) <<materialsManager->GetMaterial
               (geometry->GetLayerMaterialIndex(i))->GetName() <<": "
             <<std::setw(6)
             <<G4BestUnit(geometry->GetLayerThickness(i),"Length")
             <<G4endl;
    }
  }
  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLDoseAnalyser::DeleteLayer (G4int i)
{
  size_t j = 0;
  while (j < doseLayers.size()) {
    if (doseLayers[j] == i) {
      doseLayers.erase(doseLayers.begin()+j);
    } else if (doseLayers[j] > i) {
      doseLayers[j]--;
      j++;
    } else {
      j++;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLDoseAnalyser::BeginOfRunAction (G4double nFactor)
{
  G4int iLayer = -1;
  factor.clear();
  factorSqr.clear();
  ADose.clear();
  ADoseSqr.clear();
  dose.clear();
  doseErr.clear();

  for (size_t i=0; i<doseLayers.size(); i++) {
    iLayer = doseLayers[i];
    factor.push_back(nFactor * cm2);
    if (doseUnit ==  "eV") {
      factor[i] /= eV;
    } else if (doseUnit ==  "keV") {
      factor[i] /= keV;
    } else if (doseUnit ==  "MeV") {
      factor[i] /= MeV;
    } else if (doseUnit ==  "GeV") {
      factor[i] /= GeV;
    } else if (doseUnit ==  "TeV") {
      factor[i] /= TeV;
    } else if (doseUnit ==  "PeV") {
      factor[i] /= PeV;
    } else if (doseUnit ==  "rad") {
      factor[i] *= 100.* 1.602e-19 / eV / (geometry->GetLogicalLayer(iLayer)
        ->GetMaterial()->GetDensity()/(kg/cm3));
      if (geometry->GetShape() == SLAB) {
        factor[i] /= (geometry->GetLayerThickness(iLayer)/cm);
      } else {
        G4double rtot = geometry->GetLayerRadius(0)/cm; // the outer layer radius
        G4double rout = geometry->GetLayerRadius(iLayer)/cm;
        G4double rin  = geometry->GetLayerRadius(iLayer+1)/cm;
        factor[i]    *= 4.*rtot*rtot/(4./3.*(pow(rout,3.)-pow(rin,3.)));
      }
    } else if (doseUnit ==  "Gy") {
      factor[i] *= 1.602e-19 / eV / (geometry->GetLogicalLayer(iLayer)
        ->GetMaterial()->GetDensity()/(kg/cm3));
      if (geometry->GetShape() == SLAB) {
        factor[i] /= (geometry->GetLayerThickness(iLayer)/cm);
      } else {
        G4double rtot = geometry->GetLayerRadius(0)/cm; // the outer layer radius
        G4double rout = geometry->GetLayerRadius(iLayer)/cm;
        G4double rin  = geometry->GetLayerRadius(iLayer+1)/cm;
        factor[i]    *= 4.*rtot*rtot/(4./3.*(pow(rout,3.)-pow(rin,3.)));
      }
    }

    factorSqr.push_back(factor[i] * factor[i]);

    ADose.push_back(0.);
    ADoseSqr.push_back(0.);
    dose.push_back(-1.);
    doseErr.push_back(-1.);
  }
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLDoseAnalyser::GetNBlocks ()
{
  return Nb;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLDoseAnalyser::GetAddedNBlocks (G4int i)
{
  if (doseLayers.size() == 0) {
    Nb = i;
  } else {
    Nb = i + 1;
  }
  return Nb;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLDoseAnalyser::AddToDoseDetector (G4int id, G4double ener)
{
  ADose[id]    += ener;
  ADoseSqr[id] += ener*ener;
}

////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
void  MLDoseAnalyser::NormaliseOutput (G4double nevent)
{
  for (size_t i=0; i < doseLayers.size(); i++) {
    dose[i]    = factor[i] * ADose[i] / nevent;
    if (nevent > 1) {
      doseErr[i] = sqrt((factorSqr[i]*ADoseSqr[i]/nevent - dose[i]*dose[i])
        / (nevent-1));
    } else {
      doseErr[i] = -1.0;
    }
  }
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
CSVofstream & operator << (CSVofstream &CSVFile, MLDoseAnalyser &q)
{
//
//
// Initialise some variables.
//
  G4RunManager *runManager = G4RunManager::GetRunManager();
  MLGeometryConstruction *geometry =
    (MLGeometryConstruction*)(runManager->GetUserDetectorConstruction());
  G4int j = 0;
//
//
// Check whether any layers were selected for PHS analysis - if not don't
// output anything.
//
  if (q.GetNbOfDLayers() == 0) return CSVFile;
#ifdef SPENVIS
  G4int Nc =  4;          // Number of comment lines.
#else
  G4int Nc =  2;          // Number of comment lines.
#endif
  G4int Nm =  0;          // Number of meta-variable lines.
  G4int Na =  0;          // Number of annotation lines.
  G4int Nv =  5;          // Number of variable lines.
  G4int Nd =  5;          // Number of data columns.
  G4int Nl = q.GetNbOfDLayers(); // Number of data lines.
  G4int Nh = Nc + Nm + Nv + 1;
                         // Total number of header records(incl. current one.
//
//
// Output header record for CSV file.
//
  CSVFile <<"'*'"
          <<std::setw(6) <<Nh <<","
          <<std::setw(6) <<Nc <<","
          <<std::setw(6) <<Nm <<","
          <<std::setw(6) <<Na <<","
          <<std::setw(6) <<Nv <<","
          <<std::setw(6) <<Nd <<","
          <<std::setw(6) <<Nl <<","
          <<std::setw(6) <<q.GetNBlocks()-1 <<G4endl;
//
//
// Output comment lines - basically the version number of MULASSIS and of
// SPENVIS if applicable.
//
  MLVersion codeVersion = MLVersion();
  CSVFile <<"'DOSE ANALYSIS'" <<G4endl;
  CSVFile <<"'MULASSIS " <<codeVersion.version <<"'" <<G4endl;
#ifdef SPENVIS
  CSVFile <<"'SPENVIS " <<codeVersion.SPENVISver <<"'" <<G4endl;
  CSVFile <<"'" <<codeVersion.implementation <<"'" <<G4endl;
#endif
//
//
// Output meta-variable lines.  (CURRENTLY NONE)
//

//
//
// Output annotation lines.  (CURRENTLY NONE)
//

//
//
// Output variable lines.
//
  CSVFile <<"'Layer',   ,  1,'Layer number'"
          <<G4endl;
  CSVFile <<"'Thickness',  'cm',  1,'Thickness of layer'"
          <<G4endl;
  CSVFile <<"'Density',  'g/cm3',  1,'Density of layer'"
          <<G4endl;
  CSVFile <<"'Dose','" <<q.GetDoseUnit()
          <<"',  1,'Dose/energy deposition'" <<G4endl;
  CSVFile <<"'Error','" <<q.GetDoseUnit()
          <<"',  1,'Error dose/energy deposition'" <<G4endl;
//
//
// Loop through layers in geometry for which dose information is to be output.
//
  for (G4int i = 0; i < q.GetNbOfDLayers(); i++) {
    j = q.GetDLayerIdx(i);
    CSVFile.setf (std::ios::scientific, std::ios::floatfield );
    CSVFile <<std::setw(6) <<j+1 <<","
            <<std::setw(12) <<(geometry->GetLayerThickness(j))/cm <<","
            <<std::setw(12) <<(geometry->GetLogicalLayer(j)->GetMaterial()
              ->GetDensity())/(g/cm3) <<","
            <<std::setw(12) <<q.GetDose(i) <<","
            <<std::setw(12) <<q.GetDoseErr(i)
            <<G4endl;
  }

  if (q.GetNBlocks() > 1) {
    CSVFile <<"'End of Block'" <<G4endl;
  } else {
    CSVFile <<"'End of File'" <<G4endl;
  }

  return CSVFile;
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
RPTofstream & operator << (RPTofstream &RPTFile, MLDoseAnalyser &q)
{
//
//
// Initialise some variables.
//
  G4RunManager *runManager = G4RunManager::GetRunManager();
  MLGeometryConstruction *geometry =
    (MLGeometryConstruction*)(runManager->GetUserDetectorConstruction());
  G4int j = 0;
//
//
// Output header.
//
  RPTFile <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<"Dose Analysis:" <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<G4endl;
  if (q.GetNbOfDLayers() == 0) {
    RPTFile <<"No dose analysis output was requested" <<G4endl;
    RPTFile <<G4endl;
    return RPTFile;
  }
  RPTFile <<"Dose units = " <<q.GetDoseUnit() <<G4endl;
  RPTFile <<G4endl;
  RPTFile <<"  Layer  Thickness     Density      Dose       Error" <<G4endl;
  for (G4int i = 0; i < q.GetNbOfDLayers(); i++) {
    j = q.GetDLayerIdx(i);
    RPTFile.setf(std::ios::fixed);
    RPTFile <<std::setw(7) <<j+1 <<"  ";
    RPTFile.outG4BestUnit(G4BestUnit(geometry
      ->GetLayerThickness(j),"Length"),12);
    RPTFile <<"  ";
    RPTFile.outG4BestUnit(G4BestUnit(geometry->GetLogicalLayer(j)
      ->GetMaterial()->GetDensity(),"Volumic Mass"),12);
    RPTFile <<" ";

    RPTFile.setf(std::ios::scientific, std::ios::floatfield);
    RPTFile.precision(4);
    RPTFile <<std::setw(10) <<q.GetDose(i) <<" "
            <<std::setw(10) <<q.GetDoseErr(i)
            <<G4endl;
  }

  return RPTFile;
}
#endif
////////////////////////////////////////////////////////////////////////////////

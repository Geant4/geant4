////////////////////////////////////////////////////////////////////////////////
//
#include "MLNIELAnalyser.hh"

#include "globals.hh"
#include <iomanip>
#include <vector>

#include "MLGeometryConstruction.hh"
#include "MLFluenceAnalyser.hh"
#include "MLAnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#ifndef USEHBOOK
#include "MLVersion.hh"
#endif
////////////////////////////////////////////////////////////////////////////////
//
MLNIELAnalyser::MLNIELAnalyser (MLAnalysisManager* anMan,
  MLGeometryConstruction * det)
{
//
//
// Set default conditions.
//
  Nb       = -9999;
  doseUnit = "rad";
  nielType = CERN;
  nielFunc = new MLNIELFunction();
  factor.clear();
  factorSqr.clear();
  totalniel.clear();
  totalErr.clear();
  ATotal.clear();
  ATotalSqr.clear();
  protonniel.clear();
  protonErr.clear();
  AProton.clear();
  AProtonSqr.clear();
  EvtProtonniel.clear();
  neutronniel.clear();
  neutronErr.clear();
  ANeutron.clear();
  ANeutronSqr.clear();
  EvtNeutronniel.clear();
  electronniel.clear();
  electronErr.clear();
  AElectron.clear();
  AElectronSqr.clear();
  EvtElectronniel.clear();
  pionniel.clear();
  pionErr.clear();
  APion.clear();
  APionSqr.clear();
  EvtPionniel.clear();
  nielLayers.clear();
//
//
// Set Analysis Manager and Fluence Analyser.
//
  analysisManager = anMan;
  fluenceAnalyser = analysisManager->GetFluenceAnalyser();
//
//
// Get the geometry.
//
  geometry = det;
}
////////////////////////////////////////////////////////////////////////////////
//
MLNIELAnalyser::~MLNIELAnalyser ()
{
  delete nielFunc;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::AddNielLayer (G4int ival)
{
  G4int j = ival-1;
  G4int k = 0;
  if (std::find(nielLayers.begin(),nielLayers.end(),j) != nielLayers.end()) {
    G4cerr <<"AddNielLayer: the boundary specified is already selected."
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    while (k < geometry->GetNbOfFLayers()) {
      if (geometry->GetFLayerIdx(k) == j) break;
      k++;
    }

    if (k < geometry->GetNbOfFLayers()) {
      nielLayers.push_back(j);
    } else {
      G4cerr <<"AddNielLayer: the boundary specified is not a fluence detector."
             <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::DeleteNielLayer (G4int ival)
{
  G4int j = ival-1;
  if (j < 0) {
    nielLayers.clear();
  } else {

    size_t k = 0;
    while (k < nielLayers.size()) {
      if (nielLayers[k] == j) break;
      k++;
    }

    if (k < nielLayers.size()) {
      nielLayers.erase(nielLayers.begin()+k);
    } else {
      G4cerr <<"DeleteNielLayer: the boundary specified is not in the list."
             <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::ListNielLayer ()
{
  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
  G4cout << "Niel analysis will be carried out for the following layers:"
         <<G4endl;

  G4int i = 0;
  if (nielLayers.size() == 0) {
    G4cout <<G4endl;
    G4cout <<"    no layer has been selected !" <<G4endl;
  } else {
    for (size_t k=0; k<nielLayers.size(); k++) {
      i = nielLayers[k];
      if (i == geometry->GetNbOfLayers()) {
        G4cout <<"  Boundary before layer No. " <<std::setw(3) <<i+1
               <<G4endl;
      } else {
        G4cout <<"  Boundary before layer No. " <<std::setw(3) <<i+1 <<"  "
               <<std::setw(20) <<geometry->GetMLMaterial()->GetMaterial
                  (geometry->GetLayerMaterialIndex(i))->GetName() <<": "
               <<std::setw(6)
               <<G4BestUnit(geometry->GetLayerThickness(i),"Length")
               <<G4endl;
      }
    }
  }

  G4cout <<G4endl;
  G4cout <<"-------------------------------------------------------------"
         <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::DeleteLayer (G4int i)
{
  size_t j = 0;
  while (j < nielLayers.size()) {
    if (nielLayers[j] == i) {
      nielLayers.erase(nielLayers.begin()+j);
    } else if (nielLayers[j] > i) {
      nielLayers[j]--;
      j++;
    } else {
      j++;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::SelectNielFunction (G4String sval)
{
  if (sval == "cern" || sval == "CERN") {
    nielType = CERN;
  } else if (sval == "jpl" || sval == "JPL") {
    nielType = JPL;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::BeginOfRunAction (G4double nFactor)
{
  G4int iLayer = -1;

  factor.clear();
  factorSqr.clear();

  totalniel.clear();
  totalErr.clear();
  ATotal.clear();
  ATotalSqr.clear();
  protonniel.clear();
  protonErr.clear();
  AProton.clear();
  AProtonSqr.clear();
  EvtProtonniel.clear();
  neutronniel.clear();
  neutronErr.clear();
  ANeutron.clear();
  ANeutronSqr.clear();
  EvtNeutronniel.clear();
  electronniel.clear();
  electronErr.clear();
  AElectron.clear();
  AElectronSqr.clear();
  EvtElectronniel.clear();
  pionniel.clear();
  pionErr.clear();
  APion.clear();
  APionSqr.clear();
  EvtPionniel.clear();
//
//
// Note that MLNIELFunction generates the NIEL in MeV/(mg/cm2) therefore need
// change to kilograms.  Also normalise according to the incident fluence.
//
  for (size_t i=0; i<nielLayers.size(); i++) {
    iLayer = nielLayers[i];
    factor.push_back(nFactor * cm2);
    if (doseUnit == "MeV/g") {
      factor[i] *= gram / milligram;
    } else if (doseUnit == "rad") {
      factor[i] *= 100. * 1.602e-19 / eV * kilogram / milligram;
    } else if (doseUnit == "Gy") {
      factor[i] *= 1.602e-19 / eV * kilogram / milligram;
    }
    if (geometry->GetShape() == SPHERE) {
      G4double rtot = geometry->GetLayerRadius(0); // the outer layer radius
      G4double rout = geometry->GetLayerRadius(iLayer);
      factor[i]    *= rtot*rtot/rout/rout;
    }

    factorSqr.push_back(factor[i] * factor[i]);

    EvtProtonniel.push_back(0.);
    EvtNeutronniel.push_back(0.);
    EvtElectronniel.push_back(0.);
    EvtPionniel.push_back(0.);

    ATotal.push_back(0.);
    ATotalSqr.push_back(0.);
    AProton.push_back(0.);
    AProtonSqr.push_back(0.);
    ANeutron.push_back(0.);
    ANeutronSqr.push_back(0.);
    AElectron.push_back(0.);
    AElectronSqr.push_back(0.);
    APion.push_back(0.);
    APionSqr.push_back(0.);

    totalniel.push_back(-1.);
    totalErr.push_back(-1.);
    protonniel.push_back(-1.);
    protonErr.push_back(-1.);
    neutronniel.push_back(-1.);
    neutronErr.push_back(-1.);
    electronniel.push_back(-1.);
    electronErr.push_back(-1.);
    pionniel.push_back(-1.);
    pionErr.push_back(-1.);

  }
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLNIELAnalyser::GetNBlocks ()
{
  return Nb;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLNIELAnalyser::GetAddedNBlocks (G4int i)
{
  if (nielLayers.size() == 0) {
    Nb = i;
  } else {
    Nb = i + 1;
  }
  return Nb;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::BeginOfEventAction()
{
  for (size_t i=0; i<nielLayers.size(); i++) {
    EvtProtonniel[i]   = 0.;
    if (nielType == CERN) {
      EvtNeutronniel[i]  = 0.;
      EvtPionniel[i]     = 0.;
      EvtElectronniel[i] = 0.;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::AddToNIELDetector (G4int id, G4int part, G4double theta,
  G4double energy, G4double weight)
{
  G4double cosT = fabs(cos(theta*rad));
  if (cosT < 1.0e-4) cosT = 1.0e-4;
  G4double w = weight/cosT;
//
//
// Note that the units of NIEL generated by GetNielFunction are MeV/(mg/cm2).
// Renormalisation is taken account of using the normalisation factors
// "factor[i]" established in the BeginOfRunAction.
//
  if (part == 2212) {// proton
    EvtProtonniel[id] += nielFunc->GetNielFunction(0,nielType,energy/MeV) * w;
  } else if (nielType == CERN) {
    if (part == 2112) { // neutron
      EvtNeutronniel[id] +=
        nielFunc->GetNielFunction(1,nielType,energy/MeV) * w;
    } else if (part == 211 || part == -211) { // pion
      EvtPionniel[id] += nielFunc->GetNielFunction(2,nielType,energy/MeV) * w;
    } else if (part == 11) { // electron
      EvtElectronniel[id] +=
        nielFunc->GetNielFunction(3,nielType,energy/MeV) * w;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLNIELAnalyser::EndOfEventAction ()
{
  G4double total  = 0.;

  for (size_t i=0; i<nielLayers.size(); i++) {
    AProton[i]    += EvtProtonniel[i];
    AProtonSqr[i] += EvtProtonniel[i] * EvtProtonniel[i];
    total          = EvtProtonniel[i];

    if (nielType == CERN) {
      ANeutron[i]     += EvtNeutronniel[i];
      ANeutronSqr[i]  += EvtNeutronniel[i] * EvtNeutronniel[i];
      APion[i]        += EvtPionniel[i];
      APionSqr[i]     += EvtPionniel[i] * EvtPionniel[i];
      AElectron[i]    += EvtElectronniel[i];
      AElectronSqr[i] += EvtElectronniel[i] * EvtElectronniel[i];
      total           += EvtNeutronniel[i] + EvtPionniel[i] +
                         EvtElectronniel[i];
    }
    ATotal[i]    += total;
    ATotalSqr[i] += total * total;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
void MLNIELAnalyser::NormaliseOutput (G4double nevent)
{
  for (size_t i=0; i < nielLayers.size(); i++) {
    protonniel[i] = factor[i] * AProton[i] / nevent;
    if (nevent > 1) {
      protonErr[i] = sqrt((factorSqr[i]*AProtonSqr[i]/nevent -
                     protonniel[i]*protonniel[i]) / (nevent-1));
    } else {
      protonErr[i] = -1.0;
    }

    if (nielType == CERN) {
      neutronniel[i] = factor[i] * ANeutron[i] / nevent;
      if (nevent > 1) {
        neutronErr[i] = sqrt((factorSqr[i]*ANeutronSqr[i]/nevent -
                        neutronniel[i]*neutronniel[i]) / (nevent-1));
      } else {
        neutronErr[i] = -1.0;
      }
      pionniel[i] = factor[i] * APion[i] / nevent;
      if (nevent > 1) {
        pionErr[i] = sqrt((factorSqr[i]*APionSqr[i]/nevent -
                     pionniel[i]*pionniel[i]) / (nevent-1));
      } else {
        pionErr[i] = -1.0;
      }
      electronniel[i] = factor[i] * AElectron[i] / nevent;
      if (nevent > 1) {
        electronErr[i] = sqrt((factorSqr[i]*AElectronSqr[i]/nevent -
                     electronniel[i]*electronniel[i]) / (nevent-1));
      } else {
        electronErr[i] = -1.0;
      }
    }
    totalniel[i] = factor[i] * ATotal[i] / nevent;
    if (nevent > 1) {
      totalErr[i] = sqrt((factorSqr[i]*ATotalSqr[i]/nevent -
                     totalniel[i]*totalniel[i]) / (nevent-1));
    } else {
      totalErr[i] = -1.0;
    }
  }
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
CSVofstream & operator << (CSVofstream &CSVFile, MLNIELAnalyser &q)
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
// Check whether any layers were selected for NIEL analysis - if not don't
// output anything.
//
  if (q.GetNbOfNielLayers() == 0) return CSVFile;
#ifdef SPENVIS
  G4int Nc =  4; // Number of comment lines.
#else
  G4int Nc =  2; // Number of comment lines.
#endif
  G4int Nm =  0; // Number of meta-variable lines.
  G4int Na =  0; // Number of annotation lines.
  G4int Nv =  0; // Number of variable lines.
  G4int Nd =  0; // Number of data columns.
  if (q.GetNielType() == JPL) {
    Nv =  5;
    Nd =  5;
  } else if (q.GetNielType() == CERN) {
    Nv =  13;
    Nd =  13;
  }
  G4int Nl = q.GetNbOfNielLayers(); // Number of data lines.
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
  CSVFile <<"'NIEL ANALYSIS'" <<G4endl;
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
  if (q.GetNielType() == JPL) {
//
//
// JPL coefficients only permit calculation of proton NIEL.
//
    CSVFile <<"'ProNIEL', '" <<q.GetDoseUnit() <<"',  1,'Proton NIEL'"
            <<G4endl;
    CSVFile <<"'ErrorProNIEL', '" <<q.GetDoseUnit()
            <<"',  1,'Error in Proton NIEL'" <<G4endl;
    for (size_t i = 0; i < (size_t) q.GetNbOfNielLayers(); i++) {
      j = q.GetNLayerIdx(i);
      CSVFile.setf (std::ios::fixed);
      CSVFile <<std::setw(6) <<j+1 <<",";
      CSVFile.setf (std::ios::scientific, std::ios::floatfield );
      if (j == geometry->GetNbOfLayers()) {
        CSVFile <<"         0.0,"
                <<"         0.0,"
                <<std::setw(12) <<q.GetProton(i)      <<","
                <<std::setw(12) <<q.GetProtonErr(i)
                <<G4endl;
      } else {
        CSVFile <<std::setw(12) <<(geometry->GetLayerThickness(j))/cm <<","
                <<std::setw(12) <<(geometry->GetLogicalLayer(j)->GetMaterial()
                                 ->GetDensity()/(g/cm3)) <<","
                <<std::setw(12) <<q.GetProton(i)      <<","
                <<std::setw(12) <<q.GetProtonErr(i)
                <<G4endl;
      }
    }
  } else if (q.GetNielType() == CERN) {
//
//
// CERN ROSE NIEL coefficients are more detailed.  Output total, proton,
// neutron, charged-pion and electron NIEL.
//
    CSVFile <<"'TotalNIEL', '" <<q.GetDoseUnit() <<"',  1, 'Total NIEL'"
            <<G4endl;
    CSVFile <<"'ErrorTotalNIEL', '" <<q.GetDoseUnit()
            <<"',  1,'Error in Total NIEL'" <<G4endl;
    CSVFile <<"'ProNIEL', '" <<q.GetDoseUnit() <<"',  1, 'Proton NIEL'"
            <<G4endl;
    CSVFile <<"'ErrorProNIEL', '" <<q.GetDoseUnit()
            <<"',  1,'Error in Proton NIEL'" <<G4endl;
    CSVFile <<"'NeutNIEL', '" <<q.GetDoseUnit() <<"',  1, 'Neutron NIEL'"
            <<G4endl;
    CSVFile <<"'ErrorNeutNIEL', '" <<q.GetDoseUnit()
            <<"',  1,'Error in Neutron NIEL'" <<G4endl;
    CSVFile <<"'PionNIEL', '" <<q.GetDoseUnit() <<"',  1, 'Charged-pion NIEL'"
            <<G4endl;
    CSVFile <<"'ErrorPionNIEL', '" <<q.GetDoseUnit()
            <<"',  1,'Error in charged-pion NIEL'" <<G4endl;
    CSVFile <<"'EleNIEL', '" <<q.GetDoseUnit() <<"',  1, 'Electron NIEL'"
            <<G4endl;
    CSVFile <<"'ErrorEleNIEL', '" <<q.GetDoseUnit()
            <<"',  1,'Error in Electron NIEL'" <<G4endl;
    for (size_t i = 0; i < (size_t) q.GetNbOfNielLayers(); i++) {
      j = q.GetNLayerIdx(i);
      CSVFile.setf (std::ios::fixed);
      CSVFile <<std::setw(6) <<j+1 <<",";
      CSVFile.setf (std::ios::scientific, std::ios::floatfield);
      if (j == geometry->GetNbOfLayers()) {
        CSVFile <<"         0.0,"
                <<"         0.0,"
                <<std::setw(12) <<q.GetTotal(i)       <<","
                <<std::setw(12) <<q.GetTotalErr(i)    <<","
                <<std::setw(12) <<q.GetProton(i)      <<","
                <<std::setw(12) <<q.GetProtonErr(i)   <<","
                <<std::setw(12) <<q.GetNeutron(i)     <<","
                <<std::setw(12) <<q.GetNeutronErr(i)  <<","
                <<std::setw(12) <<q.GetPion(i)        <<","
                <<std::setw(12) <<q.GetPionErr(i)     <<","
                <<std::setw(12) <<q.GetElectron(i)    <<","
                <<std::setw(12) <<q.GetElectronErr(i)
                <<G4endl;
      } else {
        CSVFile <<std::setw(12) <<(geometry->GetLayerThickness(j))/cm <<","
                <<std::setw(12) <<(geometry->GetLogicalLayer(j)->GetMaterial()
                               ->GetDensity())/(g/cm3) <<","
                <<std::setw(12) <<q.GetTotal(i)       <<","
                <<std::setw(12) <<q.GetTotalErr(i)    <<","
                <<std::setw(12) <<q.GetProton(i)      <<","
                <<std::setw(12) <<q.GetProtonErr(i)   <<","
                <<std::setw(12) <<q.GetNeutron(i)     <<","
                <<std::setw(12) <<q.GetNeutronErr(i)  <<","
                <<std::setw(12) <<q.GetPion(i)        <<","
                <<std::setw(12) <<q.GetPionErr(i)     <<","
                <<std::setw(12) <<q.GetElectron(i)    <<","
                <<std::setw(12) <<q.GetElectronErr(i)
                <<G4endl;
        }
    }
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
RPTofstream & operator << (RPTofstream &RPTFile, MLNIELAnalyser &q)
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
  RPTFile <<"NIEL Analysis:" <<G4endl;
  RPTFile <<"-------------------------------------------------------------"
          <<G4endl;
  RPTFile <<G4endl;
  if (q.GetNbOfNielLayers() == 0) {
    RPTFile <<"No NIEL output was requested" <<G4endl;
    RPTFile <<G4endl;
    return RPTFile;
  }
  if (q.GetNielType() == JPL) {
//
//
// Output only the proton NIEL derived from JPL data.
//
    RPTFile <<"(Proton NIEL based on JPL NIEL coefficients)" <<G4endl;
    RPTFile <<"NIEL units = " <<q.GetDoseUnit() <<G4endl;
    RPTFile <<G4endl;
    RPTFile <<"Boundary  Thickness   Density    NIEL       Error" <<G4endl;
    for (G4int i = 0; i < q.GetNbOfNielLayers(); i++) {
      j = q.GetNLayerIdx(i);
      RPTFile.setf(std::ios::fixed);
      if (j == geometry->GetNbOfLayers()) {
	//        RPTFile <<"       1  0.0           0.0          ";
	RPTFile <<std::setw(8) << j+1 << "  0.0           0.0    "; 
      } else {
        RPTFile <<std::setw(8) <<j+1 <<" ";
        RPTFile.outG4BestUnit(G4BestUnit(geometry
          ->GetLayerThickness(j),"Length"),10);
        RPTFile <<" ";
        RPTFile.outG4BestUnit(G4BestUnit(geometry->GetLogicalLayer(j)
          ->GetMaterial()->GetDensity(),"Volumic Mass"),10);
        RPTFile <<" ";
      }
      RPTFile.setf(std::ios::scientific, std::ios::floatfield);
      RPTFile.precision(4);
      RPTFile <<std::setw(10) <<q.GetProton(i) <<" "
              <<std::setw(10) <<q.GetProtonErr(i)
              <<G4endl;
    }
    RPTFile <<G4endl;
  } else {
//
//
// Calculate all the NIEL values based on CERN ROSE data.
//
    RPTFile <<"(Based on CERN NIEL coefficients)" <<G4endl;
    RPTFile <<"NIEL units = " <<q.GetDoseUnit() <<G4endl;
    RPTFile <<G4endl;
    RPTFile <<"Boundary  Thickness     Density      Total      Total      Proton     Proton     Neutron    Neutron    Pion       Pion       Electron   Electron" <<G4endl;
    RPTFile <<"                                     NIEL       Error      NIEL       Error      NIEL       Error      NIEL       Error      NIEL       Error" <<G4endl;
    for (G4int i = 0; i < q.GetNbOfNielLayers(); i++) {
      j = q.GetNLayerIdx(i);
      RPTFile.setf(std::ios::fixed);
      if (j == geometry->GetNbOfLayers()) {
	//        RPTFile <<"       1  0.0           0.0          ";
        RPTFile <<std::setw(8) << j+1 << "  0.0           0.0          "; 
      } else {
        RPTFile <<std::setw(8) <<j+1 <<"  ";
        RPTFile.outG4BestUnit(G4BestUnit(geometry
          ->GetLayerThickness(j),"Length"),12);
        RPTFile <<"  ";
        RPTFile.outG4BestUnit(G4BestUnit(geometry->GetLogicalLayer(j)
          ->GetMaterial()->GetDensity(),"Volumic Mass"),12);
        RPTFile <<" ";
      }
      RPTFile.setf(std::ios::scientific, std::ios::floatfield);
      RPTFile.precision(4);
      RPTFile <<std::setw(10) <<q.GetTotal(i) <<" "
              <<std::setw(10) <<q.GetTotalErr(i) <<" "
              <<std::setw(10) <<q.GetProton(i) <<" "
              <<std::setw(10) <<q.GetProtonErr(i) <<" "
              <<std::setw(10) <<q.GetNeutron(i) <<" "
              <<std::setw(10) <<q.GetNeutronErr(i) <<" "
              <<std::setw(10) <<q.GetPion(i) <<" "
              <<std::setw(10) <<q.GetPionErr(i) <<" "
              <<std::setw(10) <<q.GetElectron(i) <<" "
              <<std::setw(10) <<q.GetElectronErr(i)
              <<G4endl;
    }
  }

  return RPTFile;
}
#endif
////////////////////////////////////////////////////////////////////////////////

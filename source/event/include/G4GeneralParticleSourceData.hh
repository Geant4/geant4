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
// G4GeneralParticleSourceData
//
// Class Description:
//
// This class uses the singleton pattern to create a single copy of the data
// needed for the G4GPS class. As yet, only the largest parts have been split
// off.
//
// Thread Safety considerations:
// Singleton creation (G4GeneralParticleSourceData::Instance()) is thread safe.
// Getters are thread-safe if no other thread is adding/deleting a source or
// normalizing. However, please note that the setters are usually accessed via
// messenger. This should be instantiated in master thread.
//
// For convenience Lock and Unlock methods are provided that can be used to
// serialize calls.
// For example:
//    gpsdata=G4GeneralParticleSourceData::Instance();
//    gpsdata->Lock();
//    gpsdata->AddASource(1.0);
//    gpsdata->Unlock();

// Author: Andrew Green, 20.03.2014
// --------------------------------------------------------------------
#ifndef G4GPS_DATA_HH
#define G4GPS_DATA_HH 1

#include "G4SingleParticleSource.hh"
#include "G4Threading.hh"

class G4GeneralParticleSourceData
{
  public:

    static G4GeneralParticleSourceData* Instance();
    
    void AddASource(G4double intensity);
    void DeleteASource(G4int idx);
    void ClearSources();
    
    void IntensityNormalise();
    
    inline G4bool Normalised() const
      { return normalised; }
    
    G4SingleParticleSource* GetCurrentSource(G4int idx);
    inline G4SingleParticleSource* GetCurrentSource() const
      { return currentSource; }

    inline G4int GetSourceVectorSize() const
      { return G4int(sourceVector.size()); }
    inline G4int GetIntensityVectorSize() const
      { return G4int(sourceIntensity.size()); }
    inline G4double GetIntensity(G4int idx) const
      { return sourceIntensity.at(idx); }
    inline G4double GetSourceProbability(G4int idx) const
      { return sourceProbability.at(idx); }
    
    void SetCurrentSourceIntensity(G4double);

    inline void SetFlatSampling(G4bool fSamp)
      { flat_sampling = fSamp; }
    inline G4bool GetFlatSampling() const
      { return flat_sampling; }

    inline void SetMultipleVertex(G4bool flag)
      { multiple_vertex = flag; }
    inline G4bool GetMultipleVertex() const
      { return multiple_vertex; }

    inline G4int GetCurrentSourceIdx() const
      { return currentSourceIdx; }

    void SetVerbosityAllSources(G4int vl);

    void Lock();
    void Unlock();
      //Lock/Unlock shared mutex
    
  private:

    G4GeneralParticleSourceData();
   ~G4GeneralParticleSourceData(); 
    
  private:
    
    std::vector<G4SingleParticleSource*> sourceVector;
    std::vector <G4double> sourceIntensity;
    std::vector <G4double> sourceProbability;
    
    G4bool multiple_vertex = false;
    G4bool flat_sampling = false;
    G4bool normalised = false;
    
    G4int currentSourceIdx = 0;
    G4SingleParticleSource* currentSource = nullptr;
    G4Mutex mutex;
};

#endif

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
// G4FermiBreakUp alternative de-excitation model
// by A. Novikov (January 2025)
//

#ifndef G4FERMISINGLETON_HH
#define G4FERMISINGLETON_HH

#include <memory>

#include "G4FermiLogger.hh"

namespace fbu
{
template<typename T>
class G4FermiSingleton
{
  public:
    G4FermiSingleton()
    {
      if (FERMI_UNLIKELY(instance_ == nullptr)) {
        instance_ = std::make_unique<T>();
      }
    }

    template<typename... Args>
    G4FermiSingleton(Args&&... args)
    {
      Reset(std::forward<Args>(args)...);
    }

    G4FermiSingleton(T* ptr) { Reset(ptr); }

    template<typename... Args>
    static void Reset(Args&&... args)
    {
      instance_ = std::make_unique<T>(std::forward<Args>(args)...);
    }

    static void Reset(T* ptr) { instance_.reset(ptr); }

    static T& Instance()
    {
      return *G4FermiSingleton();
    }

    T& operator*() { return *instance_; }

    const T& operator*() const { return *instance_; }

    T* operator->() { return instance_.get(); }

    const T& operator->() const { return *instance_; }

  private:
    static inline std::unique_ptr<T> instance_;
};

}  // namespace fbu

#endif  // G4FERMISINGLETON_HH

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
//  ---------------------------------------------------------------
//  GEANT4 template class header
//  Class Description:
//      This class delegates which singletons a task references
//  ---------------------------------------------------------------
//  Author: Jonathan Madsen
//  ---------------------------------------------------------------

#include "G4AutoLock.hh"
#include "G4Threading.hh"

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <thread>
#include <tuple>
#include <type_traits>
#include <vector>

namespace G4Traits
{
  template <typename T>
  struct TaskSingletonKey
  {
    using type         = G4int;
    using compare_type = std::less<type>;
  };
}  // namespace G4Traits

//----------------------------------------------------------------------------//

/// \class G4TaskSingletonEvaluator
/// \brief This structure must be specialized and use overloads to the
/// constructor
///
template <typename T>
struct G4TaskSingletonEvaluator;

//----------------------------------------------------------------------------//

template <typename T>
class G4TaskSingletonDelegator;

//----------------------------------------------------------------------------//

template <typename T>
class G4TaskSingletonData
{
  using key_type     = typename G4Traits::TaskSingletonKey<T>::type;
  using compare_type = typename G4Traits::TaskSingletonKey<T>::compare_type;

  friend struct G4TaskSingletonEvaluator<T>;
  friend class G4TaskSingletonDelegator<T>;
  using this_type = G4TaskSingletonData<T>;

 private:
  static std::unique_ptr<this_type>& GetInstance()
  {
    static auto _instance = std::unique_ptr<this_type>(new this_type);
    return _instance;
  }

 private:
  std::map<key_type, T*, compare_type> m_data;
};

//----------------------------------------------------------------------------//

template <typename T>
struct G4TaskSingletonEvaluator
{
  using key_type  = typename G4Traits::TaskSingletonKey<T>::type;
  using data_type = G4TaskSingletonData<T>;

  template <typename... Args>
  G4TaskSingletonEvaluator(key_type&, Args&&...)
  {
    throw std::runtime_error("Not specialized!");
  }
};

//----------------------------------------------------------------------------//

template <typename T>
class G4TaskSingletonDelegator
{
 public:
  using pointer        = T*;
  using evaluator_type = G4TaskSingletonEvaluator<T>;
  using data_type      = G4TaskSingletonData<T>;
  using key_type       = typename G4Traits::TaskSingletonKey<T>;

  template <typename... Args>
  static void Configure(Args&&... args)
  {
    auto& _data = data_type::GetInstance();
    evaluator_type _eval(std::forward<Args>(args)...);
    auto _ptr = _data->m_data.at(_eval.index);
    _eval.modifier(_ptr);
    T::SetInstance(_data->m_data.at(_key));
  }
};

//----------------------------------------------------------------------------//

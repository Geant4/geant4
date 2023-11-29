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
// G4ConvergenceTester class implementation
//
// Author: Koi, Tatsumi (SLAC/SCCS)
// --------------------------------------------------------------------

#include "G4ConvergenceTester.hh"
#include "G4AutoLock.hh"
#include <iomanip>

namespace
{
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
}

G4ConvergenceTester::G4ConvergenceTester(const G4String& theName)
  : name(theName)
{
  nonzero_histories.clear();
  largest_scores.clear();
  largest_scores.push_back(0.0);

  history_grid.resize(noBinOfHistory, 0);
  mean_history.resize(noBinOfHistory, 0.0);
  var_history.resize(noBinOfHistory, 0.0);
  sd_history.resize(noBinOfHistory, 0.0);
  r_history.resize(noBinOfHistory, 0.0);
  vov_history.resize(noBinOfHistory, 0.0);
  fom_history.resize(noBinOfHistory, 0.0);
  shift_history.resize(noBinOfHistory, 0.0);
  e_history.resize(noBinOfHistory, 0.0);
  r2eff_history.resize(noBinOfHistory, 0.0);
  r2int_history.resize(noBinOfHistory, 0.0);

  timer = new G4Timer();
  timer->Start();
  cpu_time.clear();
  cpu_time.push_back(0.0);
}

G4ConvergenceTester::~G4ConvergenceTester()
{
  delete timer;
}

void G4ConvergenceTester::AddScore(G4double x)
{
  G4AutoLock l(&aMutex);

  timer->Stop();
  cpu_time.push_back(timer->GetSystemElapsed() + timer->GetUserElapsed());

  if(x < 0.0)
  {
    std::ostringstream message;
    message << "Expecting zero or positive number as inputs,\n"
            << "but received a negative number.";
    G4Exception("G4ConvergenceTester::AddScore()", "Warning",
                JustWarning, message);
  }

  if(x == 0.0)
  {
  }
  else
  {
    nonzero_histories.insert(std::pair<G4int, G4double>(n, x));
    if(x > largest_scores.back())
    {
      // Following search should become faster if begin from bottom.
      for(auto it = largest_scores.begin(); it != largest_scores.end(); ++it)
      {
        if(x > *it)
        {
          largest_scores.insert(it, x);
          break;
        }
      }

      if(largest_scores.size() > 201)
      {
        largest_scores.pop_back();
      }
    }
    sum += x;
  }

  // Data has been added so statistics have now been updated to new values
  statsAreUpdated = false;
  ++n;
  l.unlock();
  return;
}

void G4ConvergenceTester::calStat()
{
  efficiency = G4double(nonzero_histories.size()) / n;

  mean = sum / n;

  G4double sum_x2 = 0.0;
  var             = 0.0;
  shift           = 0.0;
  vov             = 0.0;

  G4double xi;
  for(const auto& nonzero_historie : nonzero_histories)
  {
    xi = nonzero_historie.second;
    sum_x2 += xi * xi;
    var += (xi - mean) * (xi - mean);
    shift += (xi - mean) * (xi - mean) * (xi - mean);
    vov += (xi - mean) * (xi - mean) * (xi - mean) * (xi - mean);
  }

  var += (n - nonzero_histories.size()) * mean * mean;
  shift += (n - nonzero_histories.size()) * mean * mean * mean * (-1);
  vov += (n - nonzero_histories.size()) * mean * mean * mean * mean;

  if(var != 0.0)
  {
    vov = vov / (var * var) - 1.0 / n;

    var = var / (n - 1);

    sd = std::sqrt(var);

    r = sd / mean / std::sqrt(G4double(n));

    r2eff = (1 - efficiency) / (efficiency * n);
    r2int = sum_x2 / (sum * sum) - 1 / (efficiency * n);

    shift = shift / (2 * var * n);

    fom = 1 / (r * r) / cpu_time.back();
  }

  // Find Largest History
  // G4double largest = 0.0;
  largest                        = 0.0;
  largest_score_happened         = 0;
  G4double spend_time_of_largest = 0.0;
  for(const auto& nonzero_historie : nonzero_histories)
  {
    if(std::abs(nonzero_historie.second) > largest)
    {
      largest                = nonzero_historie.second;
      largest_score_happened = nonzero_historie.first;
      spend_time_of_largest =
        cpu_time[nonzero_historie.first + 1] - cpu_time[nonzero_historie.first];
    }
  }

  mean_1  = 0.0;
  var_1   = 0.0;
  shift_1 = 0.0;
  vov_1   = 0.0;
  sd_1    = 0.0;
  r_1     = 0.0;
  vov_1   = 0.0;

  mean_1 = (sum + largest) / (n + 1);

  for(const auto& nonzero_historie : nonzero_histories)
  {
    xi = nonzero_historie.second;
    var_1 += (xi - mean_1) * (xi - mean_1);
    shift_1 += (xi - mean_1) * (xi - mean_1) * (xi - mean_1);
    vov_1 += (xi - mean_1) * (xi - mean_1) * (xi - mean_1) * (xi - mean_1);
  }
  xi = largest;
  var_1 += (xi - mean_1) * (xi - mean_1);
  shift_1 += (xi - mean_1) * (xi - mean_1) * (xi - mean_1);
  vov_1 += (xi - mean_1) * (xi - mean_1) * (xi - mean_1) * (xi - mean_1);

  var_1 += (n - nonzero_histories.size()) * mean_1 * mean_1;

  if(var_1 != 0.0)
  {
    shift_1 += (n - nonzero_histories.size()) * mean_1 * mean_1 * mean_1 * (-1);
    vov_1 += (n - nonzero_histories.size()) * mean_1 * mean_1 * mean_1 * mean_1;

    vov_1 = vov_1 / (var_1 * var_1) - 1.0 / (n + 1);

    var_1 = var_1 / n;

    sd_1 = std::sqrt(var_1);

    r_1 = sd_1 / mean_1 / std::sqrt(G4double(n + 1));

    shift_1 = shift_1 / (2 * var_1 * (n + 1));

    fom_1 = 1 / (r * r) / (cpu_time.back() + spend_time_of_largest);
  }

  if(nonzero_histories.size() < 500)
  {
    calcSLOPE = false;
  }
  else
  {
    G4int i = G4int(nonzero_histories.size());

    // 5% criterion
    G4int j = G4int(i * 0.05);
    while(G4int(largest_scores.size()) > j)
    {
      largest_scores.pop_back();
    }
    calc_slope_fit(largest_scores);
  }

  calc_grid_point_of_history();
  calc_stat_history();

  // statistics have been calculated and this function does not need
  // to be called again until data has been added
  statsAreUpdated = true;
}

void G4ConvergenceTester::calc_grid_point_of_history()
{
  // histroy_grid [ 0,,,15 ]
  // history_grid [0] 1/16 ,,, history_grid [15] 16/16
  // if number of event is x then history_grid [15] become x-1.
  // 16 -> noBinOfHisotry

  for(G4int i = 1; i <= noBinOfHistory; ++i)
  {
    history_grid[i - 1] = G4int(n / (G4double(noBinOfHistory)) * i - 0.1);
  }
}

void G4ConvergenceTester::calc_stat_history()
{
  if(history_grid[0] == 0)
  {
    showHistory = false;
    return;
  }

  for(G4int i = 0; i < noBinOfHistory; ++i)
  {
    G4int ith = history_grid[i];

    G4int nonzero_till_ith = 0;
    G4double xi;
    G4double mean_till_ith = 0.0;

    for(const auto& itr : nonzero_histories)
    {
      if(itr.first <= ith)
      {
        xi = itr.second;
        mean_till_ith += xi;
        ++nonzero_till_ith;
      }
    }

    if(nonzero_till_ith == 0)
    {
      continue;
    }

    mean_till_ith   = mean_till_ith / (ith + 1);
    mean_history[i] = mean_till_ith;

    G4double sum_x2_till_ith = 0.0;
    G4double var_till_ith    = 0.0;
    G4double vov_till_ith    = 0.0;
    G4double shift_till_ith  = 0.0;

    for(const auto& itr : nonzero_histories)
    {
      if(itr.first <= ith)
      {
        xi = itr.second;
        sum_x2_till_ith += std::pow(xi, 2.0);
        var_till_ith += std::pow(xi - mean_till_ith, 2.0);
        shift_till_ith += std::pow(xi - mean_till_ith, 3.0);
        vov_till_ith += std::pow(xi - mean_till_ith, 4.0);
      }
    }

    var_till_ith +=
      ((ith + 1) - nonzero_till_ith) * std::pow(mean_till_ith, 2.0);
    vov_till_ith +=
      ((ith + 1) - nonzero_till_ith) * std::pow(mean_till_ith, 4.0);

    G4double sum_till_ith = mean_till_ith * (ith + 1);

    if(!(std::fabs(var_till_ith) > 0.0))
    {
      continue;
    }
    if(!(std::fabs(mean_till_ith) > 0.0))
    {
      continue;
    }
    if(!(std::fabs(sum_till_ith) > 0.0))
    {
      continue;
    }

    vov_till_ith = vov_till_ith / std::pow(var_till_ith, 2.0) - 1.0 / (ith + 1);
    vov_history[i] = vov_till_ith;

    var_till_ith   = var_till_ith / (ith + 1 - 1);
    var_history[i] = var_till_ith;
    sd_history[i]  = std::sqrt(var_till_ith);
    r_history[i] =
      std::sqrt(var_till_ith) / mean_till_ith / std::sqrt(1.0 * (ith + 1));

    if(std::fabs(cpu_time[ith]) > 0.0 && std::fabs(r_history[i]) > 0.0)
    {
      fom_history[i] = 1.0 / std::pow(r_history[i], 2.0) / cpu_time[ith];
    }
    else
    {
      fom_history[i] = 0.0;
    }

    shift_till_ith +=
      ((ith + 1) - nonzero_till_ith) * std::pow(mean_till_ith, 3.0) * (-1.0);
    shift_till_ith   = shift_till_ith / (2 * var_till_ith * (ith + 1));
    shift_history[i] = shift_till_ith;

    e_history[i] = 1.0 * nonzero_till_ith / (ith + 1);
    if(std::fabs(e_history[i]) > 0.0)
    {
      r2eff_history[i] = (1 - e_history[i]) / (e_history[i] * (ith + 1));

      r2int_history[i] = (sum_x2_till_ith) / std::pow(sum_till_ith, 2.0) -
                         1 / (e_history[i] * (ith + 1));
    }
  }
}

void G4ConvergenceTester::ShowResult(std::ostream& out)
{
  // if data has been added since the last computation of the statistical values
  // (not statsAreUpdated) call calStat to recompute the statistical values
  if(!statsAreUpdated)
  {
    calStat();
  }

  out << std::setprecision(6);

  out << G4endl;
  out << "G4ConvergenceTester Output Result of " << name << G4endl;
  out << std::setw(20) << "EFFICIENCY = " << std::setw(13) << efficiency
      << G4endl;
  out << std::setw(20) << "MEAN = " << std::setw(13) << mean << G4endl;
  out << std::setw(20) << "VAR = " << std::setw(13) << var << G4endl;
  out << std::setw(20) << "SD = " << std::setw(13) << sd << G4endl;
  out << std::setw(20) << "R = " << std::setw(13) << r << G4endl;
  out << std::setw(20) << "SHIFT = " << std::setw(13) << shift << G4endl;
  out << std::setw(20) << "VOV = " << std::setw(13) << vov << G4endl;
  out << std::setw(20) << "FOM = " << std::setw(13) << fom << G4endl;

  out << std::setw(20) << "THE LARGEST SCORE = " << std::setw(13) << largest
      << " and it happened at " << largest_score_happened << "th event"
      << G4endl;
  if(mean != 0)
  {
    out << std::setw(20) << "Affected Mean = " << std::setw(13) << mean_1
        << " and its ratio to original is " << mean_1 / mean << G4endl;
  }
  else
  {
    out << std::setw(20) << "Affected Mean = " << std::setw(13) << mean_1
        << G4endl;
  }
  if(var != 0)
  {
    out << std::setw(20) << "Affected VAR = " << std::setw(13) << var_1
        << " and its ratio to original is " << var_1 / var << G4endl;
  }
  else
  {
    out << std::setw(20) << "Affected VAR = " << std::setw(13) << var_1
        << G4endl;
  }
  if(r != 0)
  {
    out << std::setw(20) << "Affected R = " << std::setw(13) << r_1
        << " and its ratio to original is " << r_1 / r << G4endl;
  }
  else
  {
    out << std::setw(20) << "Affected R = " << std::setw(13) << r_1 << G4endl;
  }
  if(shift != 0)
  {
    out << std::setw(20) << "Affected SHIFT = " << std::setw(13) << shift_1
        << " and its ratio to original is " << shift_1 / shift << G4endl;
  }
  else
  {
    out << std::setw(20) << "Affected SHIFT = " << std::setw(13) << shift_1
        << G4endl;
  }
  if(fom != 0)
  {
    out << std::setw(20) << "Affected FOM = " << std::setw(13) << fom_1
        << " and its ratio to original is " << fom_1 / fom << G4endl;
  }
  else
  {
    out << std::setw(20) << "Affected FOM = " << std::setw(13) << fom_1
        << G4endl;
  }

  if(!showHistory)
  {
    out << "Number of events of this run is too small to do convergence tests."
        << G4endl;
    return;
  }

  check_stat_history(out);

  // check SLOPE and output result
  if(calcSLOPE)
  {
    if(slope >= 3)
    {
      noPass++;
      out << "SLOPE is large enough" << G4endl;
    }
    else
    {
      out << "SLOPE is not large enough" << G4endl;
    }
  }
  else
  {
    out << "Number of non zero history too small to calculate SLOPE" << G4endl;
  }

  out << "This result passes " << noPass << " / " << noTotal
      << " Convergence Test." << G4endl;
  out << G4endl;
}

void G4ConvergenceTester::ShowHistory(std::ostream& out)
{
  if(!showHistory)
  {
    out << "Number of events of this run is too small to show history."
        << G4endl;
    return;
  }

  out << std::setprecision(6);

  out << G4endl;
  out << "G4ConvergenceTester Output History of " << name << G4endl;
  out << "i/" << noBinOfHistory << " till_ith      mean" << std::setw(13)
      << "var" << std::setw(13) << "sd" << std::setw(13) << "r" << std::setw(13)
      << "vov" << std::setw(13) << "fom" << std::setw(13) << "shift"
      << std::setw(13) << "e" << std::setw(13) << "r2eff" << std::setw(13)
      << "r2int" << G4endl;
  for(G4int i = 1; i <= noBinOfHistory; i++)
  {
    out << std::setw(4) << i << " " << std::setw(5) << history_grid[i - 1]
        << std::setw(13) << mean_history[i - 1] << std::setw(13)
        << var_history[i - 1] << std::setw(13) << sd_history[i - 1]
        << std::setw(13) << r_history[i - 1] << std::setw(13)
        << vov_history[i - 1] << std::setw(13) << fom_history[i - 1]
        << std::setw(13) << shift_history[i - 1] << std::setw(13)
        << e_history[i - 1] << std::setw(13) << r2eff_history[i - 1]
        << std::setw(13) << r2int_history[i - 1] << G4endl;
  }
}

void G4ConvergenceTester::check_stat_history(std::ostream& out)
{
  // 1 sigma rejection for null hypothesis

  std::vector<G4double> first_ally;
  std::vector<G4double> second_ally;

  // use 2nd half of hisories
  std::size_t N = mean_history.size() / 2;
  std::size_t i;

  G4double pearson_r;
  G4double t;

  first_ally.resize(N);
  second_ally.resize(N);

  G4double sum_of_var =
    std::accumulate(var_history.begin(), var_history.end(), 0.0);
  if(sum_of_var == 0.0)
  {
    out << "Variances in all historical grids are zero." << G4endl;
    out << "Terminating checking behavior of statistics numbers." << G4endl;
    return;
  }

  // Mean

  for(i = 0; i < N; ++i)
  {
    first_ally[i]  = history_grid[N + i];
    second_ally[i] = mean_history[N + i];
  }

  pearson_r = calc_Pearson_r((G4int)N, first_ally, second_ally);
  t         = pearson_r * std::sqrt((N - 2) / (1 - pearson_r * pearson_r));

  if(t < 0.429318)  // Student t of (Degree of freedom = N-2 )
  {
    out << "MEAN distribution is  RANDOM" << G4endl;
    noPass++;
  }
  else
  {
    out << "MEAN distribution is not RANDOM" << G4endl;
  }

  // R

  for(i = 0; i < N; ++i)
  {
    first_ally[i]  = 1.0 / std::sqrt(G4double(history_grid[N + i]));
    second_ally[i] = r_history[N + i];
  }

  pearson_r = calc_Pearson_r(G4int(N), first_ally, second_ally);
  t         = pearson_r * std::sqrt((N - 2) / (1 - pearson_r * pearson_r));

  if(t > 1.090546)
  {
    out << "r follows 1/std::sqrt(N)" << G4endl;
    noPass++;
  }
  else
  {
    out << "r does not follow 1/std::sqrt(N)" << G4endl;
  }

  if(is_monotonically_decrease(second_ally))
  {
    out << "r is monotonically decrease " << G4endl;
  }
  else
  {
    out << "r is NOT monotonically decrease " << G4endl;
  }

  if(r_history.back() < 0.1)
  {
    out << "r is less than 0.1. r = " << r_history.back() << G4endl;
    noPass++;
  }
  else
  {
    out << "r is NOT less than 0.1. r = " << r_history.back() << G4endl;
  }

  // VOV
  for(i = 0; i < N; ++i)
  {
    first_ally[i]  = 1.0 / history_grid[N + i];
    second_ally[i] = vov_history[N + i];
  }

  pearson_r = calc_Pearson_r(G4int(N), first_ally, second_ally);
  t         = pearson_r * std::sqrt((N - 2) / (1 - pearson_r * pearson_r));

  if(t > 1.090546)
  {
    out << "VOV follows 1/std::sqrt(N)" << G4endl;
    noPass++;
  }
  else
  {
    out << "VOV does not follow 1/std::sqrt(N)" << G4endl;
  }

  if(is_monotonically_decrease(second_ally))
  {
    out << "VOV is monotonically decrease " << G4endl;
  }
  else
  {
    out << "VOV is NOT monotonically decrease " << G4endl;
  }

  // FOM

  for(i = 0; i < N; ++i)
  {
    first_ally[i]  = history_grid[N + i];
    second_ally[i] = fom_history[N + i];
  }

  pearson_r = calc_Pearson_r(G4int(N), first_ally, second_ally);
  t         = pearson_r * std::sqrt((N - 2) / (1 - pearson_r * pearson_r));

  if(t < 0.429318)
  {
    out << "FOM distribution is RANDOM" << G4endl;
    noPass++;
  }
  else
  {
    out << "FOM distribution is not RANDOM" << G4endl;
  }
}

G4double G4ConvergenceTester::calc_Pearson_r(G4int N,
                                             std::vector<G4double> first_ally,
                                             std::vector<G4double> second_ally)
{
  G4double first_mean  = 0.0;
  G4double second_mean = 0.0;

  G4int i;
  for(i = 0; i < N; i++)
  {
    first_mean += first_ally[i];
    second_mean += second_ally[i];
  }
  first_mean  = first_mean / N;
  second_mean = second_mean / N;

  G4double a = 0.0;
  for(i = 0; i < N; ++i)
  {
    a += (first_ally[i] - first_mean) * (second_ally[i] - second_mean);
  }

  G4double b1 = 0.0;
  G4double b2 = 0.0;
  for(i = 0; i < N; ++i)
  {
    b1 += (first_ally[i] - first_mean) * (first_ally[i] - first_mean);
    b2 += (second_ally[i] - second_mean) * (second_ally[i] - second_mean);
  }

  G4double rds = a / std::sqrt(b1 * b2);

  return rds;
}

G4bool G4ConvergenceTester::is_monotonically_decrease(
  const std::vector<G4double>& ally)
{
  for(auto it = ally.cbegin(); it != ally.cend() - 1; ++it)
  {
    if(*it < *(it + 1))
    {
      return FALSE;
    }
  }

  ++noPass;
  return TRUE;
}

void G4ConvergenceTester::calc_slope_fit(const std::vector<G4double>&)
{
  // create PDF bins
  G4double max = largest_scores.front();
  G4int last   = G4int(largest_scores.size());
  G4double min = 0.0;
  if(largest_scores.back() != 0)
  {
    min = largest_scores.back();
  }
  else
  {
    min  = largest_scores[last - 1];
    last = last - 1;
  }

  if(max * 0.99 < min)
  {
    // upper limit is assumed to have been reached
    slope = 10.0;
    return;
  }

  std::vector<G4double> pdf_grid;

  pdf_grid.resize(noBinOfPDF + 1);  // no grid  = no bins + 1
  pdf_grid[0]          = max;
  pdf_grid[noBinOfPDF] = min;
  G4double log10_max   = std::log10(max);
  G4double log10_min   = std::log10(min);
  G4double log10_delta = log10_max - log10_min;
  for(G4int i = 1; i < noBinOfPDF; ++i)
  {
    pdf_grid[i] = std::pow(10.0, log10_max - log10_delta / 10.0 * (i));
  }

  std::vector<G4double> pdf;
  pdf.resize(noBinOfPDF);

  for(G4int j = 0; j < last; ++j)
  {
    for(G4int i = 0; i < 11; ++i)
    {
      if(largest_scores[j] >= pdf_grid[i + 1])
      {
        pdf[i] += 1.0 / (pdf_grid[i] - pdf_grid[i + 1]) / n;
        break;
      }
    }
  }

  f_xi.resize(noBinOfPDF);
  f_yi.resize(noBinOfPDF);
  for(G4int i = 0; i < noBinOfPDF; ++i)
  {
    f_xi[i] = (pdf_grid[i] + pdf_grid[i + 1]) / 2;
    f_yi[i] = pdf[i];
  }

  // number of variables ( a and k )
  minimizer = new G4SimplexDownhill<G4ConvergenceTester>(this, 2);
  // G4double minimum =  minimizer->GetMinimum();
  std::vector<G4double> mp = minimizer->GetMinimumPoint();
  G4double k               = mp[1];

  // G4cout << "SLOPE " << 1/mp[1]+1 << G4endl;
  // G4cout << "SLOPE  a " << mp[0] << G4endl;
  // G4cout << "SLOPE  k " << mp[1] << G4endl;
  // G4cout << "SLOPE  minimum " << minimizer->GetMinimum() << G4endl;

  slope = 1 / mp[1] + 1;
  if(k < 1.0 / 9)  // Please look Pareto distribution with "sigma=a" and "k"
  {
    slope = 10;
  }
  if(slope > 10)
  {
    slope = 10;
  }
}

G4double G4ConvergenceTester::slope_fitting_function(std::vector<G4double> x)
{
  G4double a = x[0];
  G4double k = x[1];

  if(a <= 0)
  {
    return 3.402823466e+38;  // FLOAT_MAX
  }
  if(k == 0)
  {
    return 3.402823466e+38;  // FLOAT_MAX
  }

  // f_xi and f_yi is filled at "calc_slope_fit"

  G4double y = 0.0;
  for(G4int i = 0; i < G4int(f_yi.size()); ++i)
  {
    // if ( 1/a * ( 1 + k * f_xi [ i ] / a ) < 0 )
    if((1 + k * f_xi[i] / a) < 0)
    {
      y += 3.402823466e+38;  // FLOAT_MAX
    }
    else
    {
      y += (f_yi[i] - 1 / a * std::pow(1 + k * f_xi[i] / a, -1 / k - 1)) *
           (f_yi[i] - 1 / a * std::pow(1 + k * f_xi[i] / a, -1 / k - 1));
    }
  }

  return y;
}

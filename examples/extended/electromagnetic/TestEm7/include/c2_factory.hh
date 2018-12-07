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
/**
 *  \file electromagnetic/TestEm7/include/c2_factory.hh
 *  \brief Provides a factory class to avoid an infinite number of template
 *  \brief declarations.
 *
 *  \author Created by R. A. Weller and Marcus H. Mendenhall on 7/9/05.
 *  \author 2005 Vanderbilt University.
 *
 *         \version c2_factory.hh,v 1.13 2008/05/22 12:45:19 marcus Exp
 */
//

#ifndef __has_c2_factory_hh
#define __has_c2_factory_hh 1

#include "c2_function.hh"

/// \brief a factory of pre-templated c2_function generators
/// 
/// do \code
/// typedef c2_ptr<double> c2_p;
/// static c2_factory<double> c2;
/// c2_p f=c2.sin();
/// \endcode
/// \note The factory class doesn't contain any data.  
/// It can be statically instantiated at the top of a file,
/// and used everywhere inside, or even instantiated in your project's top-level
/// include file.
/// \see c2_math_factory
/// \ingroup factories
template <typename float_type> class c2_factory {
public:

/// make a *new object
static c2_classic_function_p<float_type> 
       &classic_function(float_type (*c_func)(float_type)) 
        { return *new c2_classic_function_p<float_type>(c_func); }
/// make a *new object
static c2_plugin_function_p<float_type> &plugin_function() 
        { return *new c2_plugin_function_p<float_type>(); }
/// make a *new object
static c2_plugin_function_p<float_type>
       &plugin_function(c2_function<float_type> &f) 
        { return *new c2_plugin_function_p<float_type>(f); }
/// make a *new object
static c2_const_plugin_function_p<float_type> &const_plugin_function() 
        { return *new c2_const_plugin_function_p<float_type>(); }
/// make a *new object
static c2_const_plugin_function_p<float_type>
       &const_plugin_function(const c2_function<float_type> &f) 
        { return *new c2_const_plugin_function_p<float_type>(f); }
/// make a *new object
static c2_scaled_function_p<float_type>
       &scaled_function(const c2_function<float_type> &outer, float_type scale) 
        { return *new c2_scaled_function_p<float_type>(outer, scale); }
/// make a *new object
static c2_cached_function_p<float_type>
       &cached_function(const c2_function<float_type> &func) 
        { return *new c2_cached_function_p<float_type>(func); }
/// make a *new object
static c2_constant_p<float_type>
       &constant(float_type x) { return *new c2_constant_p<float_type>(x); }
/// make a *new object
static interpolating_function_p<float_type> &
        interpolating_function() 
                {        return *new interpolating_function_p<float_type>(); }
/// make a *new object
static lin_log_interpolating_function_p<float_type> &
        lin_log_interpolating_function() 
                {return *new lin_log_interpolating_function_p<float_type>(); }
/// make a *new object
static log_lin_interpolating_function_p<float_type> &
        log_lin_interpolating_function() 
                {return *new log_lin_interpolating_function_p<float_type>(); }
/// make a *new object
static log_log_interpolating_function_p<float_type> &
        log_log_interpolating_function() 
                {return *new log_log_interpolating_function_p<float_type>(); }
/// make a *new object
static arrhenius_interpolating_function_p<float_type> &
        arrhenius_interpolating_function()
                {return *new arrhenius_interpolating_function_p<float_type>(); }
/// make a *new object
static c2_connector_function_p<float_type> &connector_function(
        float_type x0, const c2_function<float_type> &f0, float_type x2,
        const c2_function<float_type> &f2,  bool auto_center, float_type y1)
        {return *new c2_connector_function_p<float_type>(x0, f0, x2, f2,
                                                         auto_center, y1); }
/// make a *new object
static c2_connector_function_p<float_type> &connector_function(
        const c2_fblock<float_type> &fb0, const c2_fblock<float_type> &fb2, 
        bool auto_center, float_type y1)
        { return *new c2_connector_function_p<float_type>(fb0, fb2,
                                                          auto_center, y1); }
/// make a *new object
static c2_connector_function_p<float_type> &connector_function(
        float_type x0, float_type y0, float_type yp0, float_type ypp0,  
        float_type x2, float_type y2, float_type yp2, float_type ypp2,  
        bool auto_center, float_type y1)
        { return *new c2_connector_function_p<float_type>(
                x0, y0, yp0, ypp0, x2, y2, yp2, ypp2, auto_center, y1); }
/// make a *new object
static c2_piecewise_function_p<float_type> &piecewise_function() 
        { return *new c2_piecewise_function_p<float_type>(); }
/// make a *new object
static c2_sin_p<float_type> &sin() { return *new c2_sin_p<float_type>(); }
/// make a *new object
static c2_cos_p<float_type> &cos() { return *new c2_cos_p<float_type>(); }
/// make a *new object
static c2_tan_p<float_type> &tan() { return *new c2_tan_p<float_type>(); }
/// make a *new object
static c2_log_p<float_type> &log() { return *new c2_log_p<float_type>(); }
/// make a *new object
static c2_exp_p<float_type> &exp() { return *new c2_exp_p<float_type>(); }
/// make a *new object
static c2_sqrt_p<float_type> &sqrt() { return *new c2_sqrt_p<float_type>(); }
/// make a *new object
static c2_recip_p<float_type> &recip(float_type scale=1) 
        { return *new c2_recip_p<float_type>(scale); }
/// make a *new object
static c2_identity_p<float_type> &identity()
         { return *new c2_identity_p<float_type>(); }
/// make a *new object
static c2_linear_p<float_type>
       &linear(float_type x0, float_type y0, float_type slope) 
        { return *new c2_linear_p<float_type>(x0, y0, slope); }
/// make a *new object
static c2_quadratic_p<float_type> &quadratic(
             float_type x0, float_type y0, float_type xcoef, float_type x2coef) 
             { return *new c2_quadratic_p<float_type>(x0, y0, xcoef, x2coef); }
/// make a *new object
static c2_power_law_p<float_type> &power_law(float_type scale, float_type power)
        { return *new c2_power_law_p<float_type>(scale, power); }
/// make a *new object
static c2_inverse_function_p<float_type>
        &inverse_function(const c2_function<float_type> &source) 
        { return *new c2_inverse_function_p<float_type>(source); }

#if 0
/// \brief handle template for inverse_integrated_density_bins
/// \brief <float_type, Final<float_type> >(bincenters, binheights)
template <template <typename f_t> class Final > 
    static interpolating_function_p<float_type> 
    & inverse_integrated_density_bins(const std::vector<float_type> &bincenters,
                const std::vector<float_type> &binheights)
                throw(c2_exception) 
      {
        return ::inverse_integrated_density_bins<float_type, Final<float_type> >
                        (bincenters, binheights);
      }
/// \brief handle template for inverse_integrated_density_function
/// \brief <float_type, Final<float_type> >(bincenters, binheights)
template <template <typename f_t> class Final > 
    static interpolating_function_p<float_type>
    & inverse_integrated_density_function(
                const std::vector<float_type> &bincenters,
                const c2_function<float_type> &binheights)
                throw(c2_exception) 
      {
        return ::inverse_integrated_density_function<float_type,
                                                     Final<float_type> >
                        (bincenters, binheights);
      }
#endif

};

#endif

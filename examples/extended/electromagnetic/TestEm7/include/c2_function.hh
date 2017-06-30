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
 *  \file electromagnetic/TestEm7/include/c2_function.hh
 *  \brief Provides the headers for the general c2_function algebra which
 *  fast, flexible operations on piecewise-twice-differentiable functions
 *
 *  \author Created by R. A. Weller and Marcus H. Mendenhall on 7/9/05.
 *  \author Copyright 2005 __Vanderbilt University__. All rights reserved.
 *
 *         \version c2_function.hh 490 2012-04-10 19:05:40Z marcus 
 *  \see \ref c2_factory "Factory Functions" for information on constructing 
 */

//
// $Id: c2_function.hh 104041 2017-05-09 07:44:14Z gcosmo $

#ifndef __has_c2_function_hh
#define __has_c2_function_hh 1

// MSVC does not automatically define numerical constants such as M_PI 
// this came from the msdn website, so it should be right...
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#define c2_isnan _isnan
#define c2_isfinite _finite
#else
#define c2_isnan std::isnan
#define c2_isfinite std::isfinite
#endif

#include <cmath>
#include <vector>
#include <utility>
#include <string>
#include <stdexcept>
#include <typeinfo>
#include <sstream>
#include <limits> // fails under gcc-4.3 without this here, was ok 

/// \brief the exception class for c2_function operations.
class c2_exception : public std::exception {
public:
    /// \brief construct the exception with an error message
    /// \param msgcode the message
    c2_exception(const char msgcode[]) : info(msgcode) { }
    virtual ~c2_exception() throw() { }
    /** Returns a C-style character string describing the general cause
    *  of the current error.  */
    virtual const char* what() const throw() { return info.c_str(); }
private:
    std::string info;
};

// put these forward references here, and with a bogus typename to make swig
template <typename float_type> class c2_composed_function_p;
template <typename float_type> class c2_sum_p;
template <typename float_type> class c2_diff_p;
template <typename float_type> class c2_product_p;
template <typename float_type> class c2_ratio_p;
template <typename float_type> class c2_piecewise_function_p;
template <typename float_type> class c2_quadratic_p;
template <typename float_type> class c2_ptr;
/**
        \defgroup abstract_classes Abstract Classes
        \defgroup arithmetic_functions Arithmetic Functions
        \defgroup math_functions Mathemetical Functions
        \defgroup parametric_functions Parametric Families of Functions
        \defgroup interpolators Interpolating Functions
        \defgroup containers Functions which are containers for, or functions 
                  of, other functions
        \defgroup factories Factory classes which reduce silly template typing
        \defgroup transforms Classes which provide coordinate system 
                  transformations, 
        \defgroup with derivatives
*/

/// \brief structure used to hold evaluated function data at a point.  
///
/// Contains all the information for the function at one point. 
template <typename float_type> class c2_fblock 
{        
public:
        /// \brief the abscissa
        float_type x;
        /// \brief the value of the function at \a x
        float_type y;
        /// \brief the derivative at \a x
        float_type yp;
        /// \brief the second derivative at \a x
        float_type ypp; 
        /// flag, filled in by c2_function::fill_fblock(), 
        /// indicating the derivative is NaN of Inf
        bool ypbad;
        /// is NaN of Inf
        bool yppbad; 
};

/**
 \brief the parent class for all c2_functions.
 \ingroup abstract_classes
  c2_functions know their value, first, and second derivative at 
  almost every point.
  They can be efficiently combined with binary operators, 
  via c2_binary_function, 
  composed via c2_composed_function_,
  have their roots found via find_root(),
  and be adaptively integrated via partial_integrals() or integral().
  They also can carry information with them about how to find 'interesting' 
  points on the function.  
  This information is set with set_sampling_grid() and extracted with 
  get_sampling_grid().

  Particularly important subclasses are the interpolating functions classes,
  interpolating_function , lin_log_interpolating_function, 
  log_lin_interpolating_function, 
    log_log_interpolating_function, and arrhenius_interpolating_function,
    as well as the template functions
    inverse_integrated_density_function().
 
 For a discussion of memory management, see \ref memory_management
 
 */
template <typename float_type=double> class c2_function {
public:    
    /// \brief get versioning information for the header file
    /// \return the CVS Id string
        const std::string cvs_header_vers() const { return 
                "c2_function.hh 490 2012-04-10 19:05:40Z marcus ";
        }
        
    /// \brief get versioning information for the source file
    /// \return the CVS Id string
        const std::string cvs_file_vers() const ;
        
public:
    /// \brief destructor
  virtual ~c2_function() { 
    if(sampling_grid && !no_overwrite_grid) delete sampling_grid;         
    if(root_info) delete root_info; 
    if(owner_count) {
      std::ostringstream outstr;
      outstr << "attempt to delete an object with non-zero ownership in class";
      outstr << typeid(*this).name() << std::endl;
      //throw c2_exception(outstr.str().c_str());
    }
  }
        
  /// \brief get the value and derivatives. 
  ///
  /// There is required checking for null pointers on the derivatives, 
  /// and most implementations should operate faster if derivatives are not 
  /// \param[in] x the point at which to evaluate the function
  /// \param[out] yprime the first derivative (if pointer is non-null)
  /// \param[out] yprime2 the second derivative (if pointer is non-null)
  /// \return the value of the function
  virtual float_type value_with_derivatives(float_type x, float_type *yprime, 
    float_type *yprime2) const /* throw(c2_exception) */ =0 ; 
        
    /// \brief evaluate the function in the classic way, ignoring derivatives.
    /// \param x the point at which to evaluate
    /// \return the value of the function
        inline float_type operator () (float_type x) const /* throw(c2_exception) */ 
        { return value_with_derivatives(x, (float_type *)0, (float_type *)0); } 

        /// \brief get the value and derivatives. 
    ///
    /// \param[in] x the point at which to evaluate the function
    /// \param[out] yprime the first derivative (if pointer is non-null)
    /// \param[out] yprime2 the second derivative (if pointer is non-null)
    /// \return the value of the function
    inline float_type operator () (float_type x, float_type *yprime, 
                            float_type *yprime2) const /* throw(c2_exception) */ 
        { return value_with_derivatives(x, yprime, yprime2); } 
        
  /// \brief solve f(x)==value very efficiently, with explicit knowledge 
  /// of derivatives of the function
  ///
  /// find_root solves by iterated inverse quadratic extrapolation 
  /// for a solution 
  /// to f(x)=y.  It includes checks against bad convergence, so it should 
  /// never be 
  /// able to fail.  Unlike typical secant method or fancier Brent's 
  /// method finders, 
  /// this does not depend in any strong wasy on the brackets, 
  /// unless the finder has 
  /// to resort to successive approximations to close in on a root. 
  /// Often, it is possible 
  /// to make the brackets equal to the domain of the function, if there is
  /// any clue as to where the root lies, as given by the parameter \a start.  
  /// \param lower_bracket the lower bound for the search
  /// \param upper_bracket the upper bound for the search.  
  /// Function sign must be 
  /// opposite to that at \a lower_bracket
  /// \param start starting value for the search
  /// \param value the value of the function being sought 
  /// (solves f(x) = \a value)
  /// \param[out] error If pointer is zero, errors raise exception. 
  /// Otherwise, returns error here.
  /// \param[out] final_yprime If pointer is not zero, 
  /// return derivative of function 
  /// at root
  /// \param[out] final_yprime2 If pointer is not zero, 
  /// return second derivative of 
  /// function at root
  /// \return the position of the root.
  /// \see \ref rootfinder_subsec "Root finding sample" 
  float_type find_root(float_type lower_bracket, float_type upper_bracket, 
                       float_type start, 
                       float_type value, int *error=0, 
                       float_type *final_yprime=0, 
                       float_type *final_yprime2=0 ) const /* throw(c2_exception) */; 
  /// solve f(x)=value
  /// partial_integrals uses a method with an error O(dx**10) with 
  /// full information from 
  /// the derivatives, and falls back to lower order methods 
  /// if informed of incomplete 
  /// derivatives. It uses exact midpoint splitting of the intervals 
  /// for recursion, 
  /// resulting in no recomputation of the function during recursive 
  /// descent at previously 
  /// computed points.
  /// \param xgrid points between which to evaluate definite integrals.  
  /// \param partials if non-NULL, a vector in which to receive the 
  /// partial integrals.
  /// It will automatically be sized apprpropriately, 
  /// if provided, to contain \a n - 1 
  /// elements where \a n is the length of \a xgrid  
  /// \param abs_tol the absolute error bound for each segment
  /// \param rel_tol the fractional error bound for each segment.  
  /// If the error is smaller than either the relative or absolute tolerance, 
  /// the integration step is finished.
  /// \param derivs number of derivatives to trust, 
  /// which sets the order of the integrator.  
  /// The order is 3*\a derivs + 4. \a derivs can be 0, 1, or 2.
  /// \param adapt if true, use recursive adaptation, 
  /// otherwise do simple evaluation on 
  /// the grid provided with no error checking.
  /// \param extrapolate if true, use simple Richardson 
  /// extrapolation on the final 2 steps 
  /// to reduce the error. \return sum of partial integrals, 
  /// which is the definite integral 
  /// from the first value in \a xgrid to the last.
  float_type partial_integrals(std::vector<float_type> xgrid, 
                               std::vector<float_type> *partials = 0,
                               float_type abs_tol=1e-12, 
                               float_type rel_tol=1e-12, 
                               int derivs=2, bool adapt=true, 
                               bool extrapolate=true) 
    const /* throw(c2_exception) */;
        
  /// \brief a fully-automated integrator which uses the information 
  /// provided by the 
  /// get_sampling_grid() function to figure out what to do.
  ///
  /// It returns the integral of the function over the domain requested
  /// with error tolerances as specified.  It is just a front-end 
  /// to partial_integrals()
  /// 
  /// \param amin lower bound of the domain for integration
  /// \param amax upper bound of the domain for integration
  /// \param partials if non-NULL, a vector in which to receive 
  /// the partial integrals.
  /// It will automatically be sized appropriately, 
  /// if provided, to contain \a n - 1 
  /// elements where \a n is the length of \a xgrid  
  /// \param abs_tol the absolute error bound for each segment
  /// \param rel_tol the fractional error bound for each segment.  
  /// If the error is smaller than either the relative or absolute tolerance, 
  /// the integration 
  /// step is finished.
  /// \param derivs number of derivatives to trust, which sets the 
  /// order of the integrator.  
  /// The order is 3*\a derivs + 4. \a derivs can be 0, 1, or 2.
  /// \param adapt if true, use recursive adaptation, 
  /// otherwise do simple evaluation on 
  /// the grid provided with no error checking.
  /// \param extrapolate if true, use simple Richardson 
  /// extrapolation on the final 2 steps 
  /// to reduce the error. \return sum of partial integrals, 
  /// which is the definite integral 
  /// from the first value in \a xgrid to the last.
  float_type integral(float_type amin, float_type amax, 
                      std::vector<float_type> *partials = 0,
                      float_type abs_tol=1e-12, float_type rel_tol=1e-12, 
                      int derivs=2, bool adapt=true, bool extrapolate=true) 
    const /* throw(c2_exception) */;

  /// \brief create a c2_piecewise_function_p from 
  /// c2_connector_function_p segments which 
  /// is a representation of the parent function to the specified accuracy, 
  /// but maybe much cheaper to evaluate
  ///
  /// This method has three modes, depending on the \a derivs flag. 
  ///
  /// If \a derivs is 2,
  /// it computes a c2_piecewise_function_p representation of its 
  /// parent function, 
  /// which may be a much faster 
  /// function to use in codes if the parent function is expensive.  
  /// If \a xvals 
  /// and \a yvals are non-null,
  /// it will also fill them in with the function values at each grid point the 
  /// adaptive algorithm chooses.
  ///
  /// If \a derivs is 1, this does not create the connectors, 
  /// and returns an null pointer, but will fill in the \a xvals and \a yvals 
  /// vectors with values of the function at points such that the 
  /// linear interpolation 
  /// error between the points
  /// is bounded by the tolerance values given.  
  /// Because it uses derivative information 
  /// from the function to manage the 
  /// error control, it is almost completely free of issues with 
  /// missing periods of oscillatory functions,
  /// even with no information provided in the sampling grid.
  /// This is typically useful for sampling a function for plotting.
  ///
  /// If \a derivs is 0, this does something very like what it does 
  /// if \a derivs = 1, 
  /// but without derivatives.  
  /// Instead, to compute the intermediate value of the function 
  /// for error control, 
  /// it just uses 
  /// 3-point parabolic interpolation.  
  /// This is useful amost exclusively for converting 
  /// a non-c2_function,
  /// with no derivatives, but wrapped in a c2_classic_function wrapper, 
  /// into a table 
  /// of values to seed an interpolating_function_p.
  /// Note, however, that without derivatives, this is very 
  /// susceptible to missing 
  /// periods of oscillatory 
  /// functions, so it is important to set a sampling grid 
  /// which isn't too much coarser 
  /// than the typical oscillations.
  ///
  /// \note the \a sampling_grid of the returned function matches the 
  /// \a sampling_grid of its parent.
  /// \see \ref sample_function_for_plotting "Adaptive Sampling Examples"
  /// \param amin lower bound of the domain for sampling
  /// \param amax upper bound of the domain for sampling
  /// \param abs_tol the absolute error bound for each segment
  /// \param rel_tol the fractional error bound for each segment. 
  /// \param derivs if 0 or 1, return a useless function, 
  /// but fill in the \a xvals and 
  /// \a yvals vectors (if non-null).
  /// Also, if 0 or 1, tolerances refer to linear interpolation, not high-order 
  /// interpolation. 
  /// If 2, return a full piecewise collection of 
  /// c2_connector_function_p segments.  
  /// See discussion above.
  /// \param [in,out] xvals vector of abscissas at which the function 
  /// was actually 
  /// sampled (if non-null)
  /// \param [in,out] yvals vector of function values corresponding to \a xvals 
  /// (if non-null)
  /// \return a new, sampled representation, if \a derivs is 2.  
  /// A null pointer if \a derivs is 0 or 1.
  c2_piecewise_function_p<float_type> *adaptively_sample(
                 float_type amin, float_type amax,
                 float_type abs_tol=1e-12, float_type rel_tol=1e-12,
                 int derivs=2, std::vector<float_type> *xvals=0, 
                 std::vector<float_type> *yvals=0) const /* throw(c2_exception) */;
        
  inline float_type xmin() const { return fXMin; }
  inline float_type xmax() const { return fXMax; }
  void set_domain(float_type amin, float_type amax) { fXMin=amin; fXMax=amax; }
                
  /// and sampler do increment it.
  /// \return number of evaluations logged since last reset.
  size_t get_evaluations() const { return evaluations; }
  /// \brief reset the counter
  void reset_evaluations()  const { evaluations=0; } 
  /// \brief count evaluations
  inline void increment_evaluations() const { evaluations++; }
        
  /// \brief check that a vector is monotonic, throw an exception if not, 
  /// and return a flag if it is reversed
  ///
  /// \param data a vector of data points which are expected to be monotonic.
  /// \param message an informative string to include in an exception if 
  /// this throws c2_exception
  /// \return true if in decreasing order, false if increasing 
  bool check_monotonicity(const std::vector<float_type> &data, 
                          const char message[]) const /* throw(c2_exception) */;
        
  /// \brief establish a grid of 'interesting' points on the function.
  /// 
  /// The sampling grid describes a reasonable initial set of points 
  /// to look at the function.
  /// this should generally be set at a scale which is quite coarse, 
  /// and sufficient for initializing
  /// adaptive integration or possibly root bracketing. 
  /// For sampling a function to build a new 
  /// interpolating function, one may want to refine this for accuracy.  
  /// However, 
  /// interpolating_functions themselves return their original 
  /// X grid by default, so refining 
  /// the grid in this case might be a bad idea.
  /// \param grid a vector of abscissas.  
  /// The contents is copied into an internal vector, 
  /// so the \a grid can be discarded after passingin.
  virtual void set_sampling_grid(const std::vector<float_type> &grid) 
    /* throw(c2_exception) */; 
        
  /// \brief get the sampling grid, which may be a null pointer
  /// \return pointer to the sampling grid
  std::vector<float_type> *get_sampling_grid_pointer() const { 
    return sampling_grid; } 

  virtual void get_sampling_grid(float_type amin, float_type amax, 
                                 std::vector<float_type> &grid) const ;
        
  /// The grid is modified in place.
  void preen_sampling_grid(std::vector<float_type> *result) const;
  void refine_sampling_grid(std::vector<float_type> &grid, 
                            size_t refinement) const;
        
  /// \brief create a new c2_function from this one which 
  /// is normalized on the interval 
  /// \param amin lower bound of the domain for integration
  /// \param amax upper bound of the domain for integration
  /// \param norm the desired integral for the function over the region
  /// \return a new c2_function with the desired \a norm.
  c2_function<float_type> &normalized_function(
      float_type amin, 
      float_type amax, 
      float_type norm=1.0) const /* throw(c2_exception) */;
  c2_function<float_type> &square_normalized_function(
      float_type amin, float_type amax, 
      float_type norm=1.0) 
    const /* throw(c2_exception) */;
  /// \brief create a new c2_function from this one which is square-normalized 
  /// with the provided \a weight on the interval 
  /// \param amin lower bound of the domain for integration
  /// \param amax upper bound of the domain for integration
  /// \param weight a c2_function providing the weight
  /// \param norm the desired integral for the function over the region
  /// \return a new c2_function with the desired \a norm.
  c2_function<float_type> &square_normalized_function(
      float_type amin, float_type amax, const c2_function<float_type> &weight, 
      float_type norm=1.0)
    const /* throw(c2_exception) */;

  /// \brief factory function to create a c2_sum_p from a regular 
  /// algebraic expression.
  /// \param rhs the right-hand term of the sum
  /// \return a new c2_function
  c2_sum_p<float_type> &operator + (const c2_function<float_type> &rhs) const 
  { return *new c2_sum_p<float_type>(*this, rhs); }
  /// \brief factory function to create a c2_diff_p from a regular 
  /// algebraic expression.
  /// \param rhs the right-hand term of the difference
  /// \return a new c2_function
  c2_diff_p<float_type> &operator - (const c2_function<float_type> &rhs) const 
  { return *new c2_diff_p<float_type>(*this, rhs); }
  /// \brief factory function to create a c2_product_p from a 
  /// regular algebraic expression.
  /// \param rhs the right-hand term of the product
  /// \return a new c2_function
  c2_product_p<float_type> &operator * 
  (const c2_function<float_type> &rhs) const 
  { return *new c2_product_p<float_type>(*this, rhs); }
  c2_ratio_p<float_type> &operator / (const c2_function<float_type> &rhs) const 
  { return *new c2_ratio_p<float_type>(*this, rhs); }
  /// \brief compose this function outside another.
  /// \param inner the inner function
  /// \return the composed function
  /// \anchor compose_operator
  c2_composed_function_p<float_type> & operator ()
    (const c2_function<float_type> &inner) const 
  { return *new c2_composed_function_p<float_type>((*this), inner); }

  /// \brief Find out where a calculation ran into trouble, if it got a nan.
  /// If the most recent computation did not return a nan, this is undefined.
  /// \return \a x value of point at which something went wrong, if integrator 
  /// (or otherwise) returned a nan.
  float_type get_trouble_point() const { return bad_x_point; }
        
  /// \brief increment our reference count.  
  /// Destruction is only legal if the count is zero.
  void claim_ownership() const { owner_count++; }
  /// \brief decrement our reference count. Do not destroy at zero.
  /// \return final owner count, to check whether object should disappear.
  size_t release_ownership_for_return() const /* throw(c2_exception) */ { 
    if(!owner_count) {
      std::ostringstream outstr;
      outstr << "attempt to release ownership of an unowned function in class ";
      outstr << typeid(*this).name() << std::endl;
      throw c2_exception(outstr.str().c_str());
    }
    owner_count--; 
    return owner_count;
  }
  void release_ownership() const /* throw(c2_exception) */ { 
    if(!release_ownership_for_return()) delete this;
  }
  /// \brief get the reference count, mostly for debugging
  /// \return the count
  size_t count_owners() const { return owner_count; }

protected:
  c2_function(const c2_function<float_type>  &src) 
  : sampling_grid(0),
    no_overwrite_grid(false),
    fXMin(src.fXMin), fXMax(src.fXMax), root_info(0), owner_count(0)
  {} // copy constructor only copies domain, and is only for internal use
  c2_function() : 
    sampling_grid(0), no_overwrite_grid(0), 
    fXMin(-std::numeric_limits<float_type>::max()), 
    fXMax(std::numeric_limits<float_type>::max()), root_info(0), owner_count(0)
  {}

  virtual void set_sampling_grid_pointer(std::vector<float_type> &grid) 
  { 
    if (sampling_grid && !no_overwrite_grid) delete sampling_grid; 
    sampling_grid=&grid; no_overwrite_grid=1; 
  }
        
  std::vector<float_type> * sampling_grid;
  bool no_overwrite_grid;
        
  float_type fXMin, fXMax;
  mutable size_t evaluations;
  /// \brief this point may be used to record where a calculation 
  /// ran into trouble
  mutable float_type bad_x_point;
public: 
  /// \brief fill in a c2_fblock<float_type>... a 
  /// shortcut for the integrator & sampler
  /// \param [in,out] fb the block to fill in with information
  inline void fill_fblock(c2_fblock<float_type> &fb) const /* throw(c2_exception) */ 
  {
    fb.y=value_with_derivatives(fb.x, &fb.yp, &fb.ypp);
    fb.ypbad=c2_isnan(fb.yp) || !c2_isfinite(fb.yp);
    fb.yppbad=c2_isnan(fb.ypp) || !c2_isfinite(fb.ypp);
    increment_evaluations();
  }

private:
  /// \brief the data element for the internal recursion stack for 
  /// the sampler and integrator
  struct recur_item { 
    c2_fblock<float_type> f1; size_t depth; 
    float_type previous_estimate, abs_tol, step_sum; 
    bool done; 
    size_t f0index, f2index;
  };
        

  /// \brief structure used to pass information recursively in integrator.
  ///
  /// the \a abs_tol is scaled by a factor of two at each division.  
  /// Everything else is just passed down.
  struct c2_integrate_recur { 
    c2_fblock<float_type> *f0, *f1;
    float_type abs_tol, rel_tol, eps_scale, extrap_coef, extrap2, 
      dx_tolerance, abs_tol_min;
    std::vector< recur_item > *rb_stack;
    int  derivs;
    bool adapt, extrapolate, inited;
  };

  /// \brief structure used to pass information recursively in sampler.
  ///
  struct c2_sample_recur { 
    c2_fblock<float_type> *f0, *f1;
    float_type abs_tol, rel_tol, dx_tolerance, abs_tol_min;
    int derivs;
    c2_piecewise_function_p<float_type> *out;
    std::vector<float_type> *xvals, *yvals;
    std::vector< recur_item > *rb_stack;
    bool inited;
  };

  /// \brief structure used to hold root bracketing information
  ///
  struct c2_root_info {
    c2_fblock<float_type> lower, upper;
    bool inited;
  };
        
  /// \brief Carry out the recursive subdivision and integration.
  ///
  /// This passes information recursively through the \a recur block pointer
  /// to allow very efficient recursion.
  /// \param rb a pointer to the recur struct.
  float_type integrate_step(struct c2_integrate_recur &rb) 
    const /* throw(c2_exception) */;
    
  /// \brief Carry out the recursive subdivision for sampling.
  ///
  /// This passes information recursively through the \a recur block pointer
  /// to allow very efficient recursion.
  /// \param rb a pointer to the recur struct.
  void sample_step(struct c2_sample_recur &rb) const /* throw(c2_exception) */;

  /// this carry a memory of the last root bracketing,
  /// to avoid the necessity of evaluating the function on the 
  /// brackets every time
  /// if the brackets have not been changed. 
  /// it is declared as a pointer, since many c2_functions may 
  /// never need one allocated
  mutable struct  c2_root_info *root_info;
  
  mutable size_t owner_count;
};

/// \brief a container into which any conventional c-style 
/// function can be dropped,
/// to create a degenerate c2_function without derivatives. 
/// Mostly useful for sampling into interpolating functions.
/// construct a reference to this with c2_classic_function()
/// \ingroup containers
template <typename float_type=double> class c2_classic_function_p 
  : public c2_function<float_type> {
public:
  /// \brief construct the container 
  /// \param c_func a pointer to a conventional c-style function
  c2_classic_function_p(const float_type (*c_func)(float_type)) 
    : c2_function<float_type>(), func(c_func)  {}

  /// \copydoc c2_function::value_with_derivatives
  /// Uses the internal function pointer set by set_function().
  virtual float_type value_with_derivatives(
    float_type x, float_type *yprime, 
    float_type *yprime2) const /* throw(c2_exception) */ 
  {
    if(!func) 
      throw c2_exception("c2_classic_function called with null function");
    if(yprime) *yprime=0;
    if(yprime2) *yprime2=0;
    return func(x); 
  }
  virtual ~c2_classic_function_p() { }
        
protected:
  /// \brief pointer to our function
  const float_type (*func)(float_type);
};

/// \brief create a container for a c2_function which handles the 
/// reference counting. \ingroup containers
/// It is useful as a smart container to hold a c2_function and keep 
/// the reference count correct.  
/// The recommended way for a class to store a c2_function which is 
/// handed in from the outside 
/// is for it to have a c2_ptr member into which the passed-in 
/// function is stored.
/// This way, when the class instance is deleted, 
/// it will automatically dereference any function
/// which it was handed.
///
/// This class contains a copy constructor and operator=, 
/// to make it fairly easy to make 
/// a std::vector of these objects, and have it work as expected.
template <typename float_type> class c2_const_ptr {
public:
  /// \brief construct the container with no function
  c2_const_ptr() : func(0)  {}
  /// \brief construct the container with a pre-defined function
  /// \param f the function to store
  c2_const_ptr(const c2_function<float_type> &f) : func(0) 
  { this->set_function(&f); }
  /// \brief copy constructor
  /// \param src the container to copy
  c2_const_ptr(const c2_const_ptr<float_type> &src) : func(0)
  { this->set_function(src.get_ptr()); }
  void set_function(const c2_function<float_type> *f) 
  { 
    if(func) func->release_ownership();
    func=f; 
    if(func) func->claim_ownership();
  }
        
  /// \brief fill the container from another container
  /// \param f the container to copy
  const c2_const_ptr<float_type> & operator =
  (const c2_const_ptr<float_type> &f) 
  { this->set_function(f.get_ptr()); return f; }
  /// \brief fill the container with a function
  /// \param f the function
  const c2_function<float_type> & operator =
  (const c2_function<float_type> &f) 
  { this->set_function(&f); return f; }
  /// \brief release the function without destroying it, 
  /// so it can be returned from a function
  ///
  /// This is usually the very last line of a function 
  /// before the return statement, so that 
  /// any exceptions that happen during execution of the 
  /// function will cause proper cleanup.
  /// Once the function has been released from its container this way, 
  /// it is an orhpaned object
  /// until the caller claims it, so it could get lost if an exception happens.
  void release_for_return() /* throw(c2_exception) */
  {        
    if(func) func->release_ownership_for_return();
    func=0;
  }
  /// \brief clear the function
  ///
  /// Any attempt to use this c2_plugin_function_p throws an exception 
  /// if the saved 
  /// function is cleared.
  void unset_function(void) { this->set_function(0);  }
  /// \brief destructor
  ~c2_const_ptr() { this->set_function(0); }
        
  /// \brief get a reference to our owned function
  inline const c2_function<float_type> &get() const /* throw(c2_exception) */ 
  { 
    if(!func) throw c2_exception("c2_ptr accessed uninitialized");
    return *func;
  }
  /// \brief get an unchecked pointer to our owned function
  inline const c2_function<float_type> *get_ptr() const { return func; }
  /// \brief get a checked pointer to our owned function
  inline const c2_function<float_type> *operator -> () const 
  { return &get(); }
  /// \brief check if we have a valid function
  bool valid() const { return func != 0; }
        
  /// \brief type coercion operator which lets us use a pointer as if it were 
  /// a const c2_function
  operator const c2_function<float_type>& () const { return this->get(); }
        
  /// \brief convenience operator to make us look like a function
  /// \param x the value at which to evaluate the contained function
  /// \return the evaluated function

  float_type operator()(float_type x) const /* throw(c2_exception) */ 
  { return get()(x); }
  /// \brief convenience operator to make us look like a function
  /// \param x the value at which to evaluate the contained function
  /// \param yprime the derivative
  /// \param yprime2 the second derivative
  /// \return the evaluated function
  /// \note If you using this repeatedly, 
  /// do const c2_function<float_type> &func=ptr;
  /// and use func(x).  Calling this operator wastes some time, 
  /// since it checks the alidity of the pointer every time.
  float_type operator()(float_type x, float_type *yprime, 
                        float_type *yprime2) const /* throw(c2_exception) */ 
  { return get().value_with_derivatives(x, yprime, yprime2); }
  /// \brief factory function to create a c2_sum_p from a regular 
  /// algebraic expression.
  /// \param rhs the right-hand term of the sum
  /// \return a new c2_function
  c2_sum_p<float_type> &operator + (const c2_function<float_type> &rhs) 
    const /* throw(c2_exception) */
  { return *new c2_sum_p<float_type>(get(), rhs); }
  c2_diff_p<float_type> &operator - (const c2_function<float_type> &rhs)  
    const /* throw(c2_exception) */
  { return *new c2_diff_p<float_type>(get(), rhs); }
  c2_product_p<float_type> &operator * (const c2_function<float_type> &rhs) 
    const /* throw(c2_exception) */
  { return *new c2_product_p<float_type>(get(), rhs); }
  c2_ratio_p<float_type> &operator / (const c2_function<float_type> &rhs)  
    const /* throw(c2_exception) */
  { return *new c2_ratio_p<float_type>(get(), rhs); }
  /// \brief compose this function outside another.
  /// \param inner the inner function
  /// \return the composed function
  c2_composed_function_p<float_type> & 
  operator ()(const c2_function<float_type> &inner) 
    const /* throw(c2_exception) */
  { return *new c2_composed_function_p<float_type>(get(), inner); }
        
protected:
  const c2_function<float_type> * func;
};

template <typename float_type> class c2_ptr : public c2_const_ptr<float_type >
{
public:
  /// \brief construct the container with no function
  c2_ptr() : c2_const_ptr<float_type>()  {}
  /// \brief construct the container with a pre-defined function
  /// \param f the function to store
  c2_ptr(c2_function<float_type> &f) : 
    c2_const_ptr<float_type>() { this->set_function(&f); } 
  /// \brief copy constructor
  /// \param src the container to copy
  c2_ptr(const c2_ptr<float_type> &src) : 
    c2_const_ptr<float_type>() { this->set_function(src.get_ptr()); }
  /// \brief get a checked pointer to our owned function
  inline c2_function<float_type> &get() const /* throw(c2_exception) */ 
  { return *const_cast<c2_function<float_type>*>(
      &c2_const_ptr<float_type>::get()); }
  /// \brief get an unchecked pointer to our owned function
  inline c2_function<float_type> *get_ptr() const 
  { return const_cast<c2_function<float_type>*>(this->func); }
  /// \brief get a checked pointer to our owned function
  inline c2_function<float_type> *operator -> () const 
  { return &get(); }
  /// \brief fill the container from another container
  /// \param f the container to copy
  const c2_ptr<float_type> & operator =(const c2_ptr<float_type> &f) 
  { this->set_function(f.get_ptr()); return f; }
  /// \brief fill the container with a function
  /// \param f the function
  c2_function<float_type> & operator =(c2_function<float_type> &f) 
  { this->set_function(&f); return f; }
private:
  /// \brief hidden non-const-safe version of operator=
  void operator =(const c2_const_ptr<float_type> &) { }
  /// \brief hidden non-const-safe version of operator=
  void operator =(const c2_function<float_type> &) { }
};

template <typename float_type, template <typename> class c2_class > 
class c2_typed_ptr 
  : public c2_const_ptr<float_type> {
public:
  /// \brief construct the container with no function
  c2_typed_ptr() : c2_ptr<float_type>()  {}
  /// \brief construct the container with a pre-defined function
  /// \param f the function to store
  c2_typed_ptr(c2_class<float_type> &f)
    : c2_const_ptr<float_type>() { this->set_function(&f); }
  /// \brief copy constructor
  /// \param src the container to copy
  c2_typed_ptr(const c2_typed_ptr<float_type, c2_class> &src) 
    : c2_const_ptr<float_type>() { this->set_function(src.get_ptr()); }
        
  /// \brief get a reference to our owned function
  inline c2_class<float_type> &get() const /* throw(c2_exception) */ 
  { 
    return *static_cast<c2_class<float_type> *>
      (const_cast<c2_function<float_type>*>(&c2_const_ptr<float_type>::get()));
  }
  /// \brief get a checked pointer to our owned function
  inline c2_class<float_type> *operator -> () const 
                { return &get(); }
  /// \brief get an unchecked pointer to our owned function
  inline c2_class<float_type> *get_ptr() const 
  { return static_cast<c2_class<float_type> *>(
       const_cast<c2_function<float_type>*>(this->func)); }

  operator c2_class<float_type>&() const { return get(); }
  /// \brief fill the container from another container
  /// \param f the container to copy
  void operator =(const c2_typed_ptr<float_type, c2_class> &f) 
  { this->set_function(f.get_ptr()); }
  /// \brief fill the container with a function
  /// \param f the function
  void operator =(c2_class<float_type> &f) 
  { this->set_function(&f); }
private:
  void operator =(const c2_const_ptr<float_type> &) { }
  void operator =(const c2_function<float_type> &) { }
};

template <typename float_type=double> class c2_plugin_function_p : 
        public c2_function<float_type> {
public:
  /// \brief construct the container with no function
  c2_plugin_function_p() : c2_function<float_type>(), func()  {}
  /// \brief construct the container with a pre-defined function
  c2_plugin_function_p(c2_function<float_type> &f) : 
    c2_function<float_type>(),func(f)  { }

  void set_function(c2_function<float_type> *f) 
  {
    func.set_function(f);
    if(f) this->set_domain(f->xmin(), f->xmax());
  }
  /// \copydoc c2_function::value_with_derivatives
  /// Uses the internal function pointer set by set_function().
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */ 
  {
    if(!func.valid()) 
      throw c2_exception("c2_plugin_function_p called uninitialized");
    return func->value_with_derivatives(x, yprime, yprime2); 
  }
  /// \brief destructor
  virtual ~c2_plugin_function_p() { }
        
  /// \brief clear our function
  void unset_function() { func.unset_function(); }
        
  virtual void get_sampling_grid(float_type amin, float_type amax, 
                                 std::vector<float_type> &grid) const 
  {        
    if(!func.valid()) 
      throw c2_exception("c2_plugin_function_p called uninitialized");
    if(this->sampling_grid) 
      c2_function<float_type>::get_sampling_grid(amin, amax, grid);
    else  func->get_sampling_grid(amin, amax, grid);
  }
protected:
        c2_ptr<float_type> func;
};

template <typename float_type=double> class c2_const_plugin_function_p 
  : public c2_plugin_function_p<float_type> {
public:
  /// \brief construct the container with no function
  c2_const_plugin_function_p() : c2_plugin_function_p<float_type>()  {}
  /// \brief construct the container with a pre-defined function
  c2_const_plugin_function_p(const c2_function<float_type> &f) : 
    c2_plugin_function_p<float_type>()  { this->set_function(&f); }
  void set_function(const c2_function<float_type> *f) 
  { c2_plugin_function_p<float_type>::set_function(
       const_cast<c2_function<float_type>*>(f)); }
  /// \brief destructor
  virtual ~c2_const_plugin_function_p() { }
        
  /// \brief get a const reference to our owned function, for direct access
  const c2_function<float_type> &get() const /* throw(c2_exception) */ 
  { return this->func.get(); }
};

template <typename float_type=double> class c2_binary_function 
  : public c2_function<float_type> {
public:
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */ 
  { 
    if(stub) 
      throw c2_exception("attempt to evaluate a c2_binary_function stub");
    return this->combine(*Left.get_ptr(), *Right.get_ptr(), x, yprime, yprime2);
  }

  /// \brief destructor releases ownership of member functions
  ///
  virtual ~c2_binary_function() { }

protected:
  c2_binary_function( float_type (*combiner)(
     const c2_function<float_type> &left, 
     const c2_function<float_type> &right, 
     float_type x, float_type *yprime, 
     float_type *yprime2),
                      const c2_function<float_type> &left,  
                      const c2_function<float_type> &right) : 
    c2_function<float_type>(), combine(combiner), Left(left), 
    Right(right), stub(false)
  { 
    this->set_domain(
                     (left.xmin() > right.xmin()) ? left.xmin() : right.xmin(), 
                     (left.xmax() < right.xmax()) ? left.xmax() : right.xmax()
                     );
  } 
  c2_binary_function(float_type (*combiner)(
      const c2_function<float_type> &left, 
      const c2_function<float_type> &right, 
      float_type x, float_type *yprime, float_type *yprime2)
                     ) : c2_function<float_type>(), combine(combiner), 
    Left(), Right(), stub(true) { }
        
public:
  float_type (* const combine)(
      const c2_function<float_type> &left, 
      const c2_function<float_type> &right, 
      float_type x, float_type *yprime, float_type *yprime2);
        
protected:        
  const c2_const_ptr<float_type> Left,  Right;
  bool stub;
};

template <typename float_type=double> class c2_scaled_function_p 
  : public c2_function<float_type> {
public:
  /// \brief construct the function with its scale factor.
  ///
  /// \param outer the function to be scaled
  /// \param scale the multiplicative scale factor
  c2_scaled_function_p(const c2_function<float_type> &outer, 
                       float_type scale) : 
    c2_function<float_type>(), func(outer), yscale(scale) { }
        
  /// \brief set a new scale factor
  /// \param scale the new factor
  void reset(float_type scale) { yscale=scale; }
        
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw (c2_exception) */ 
  { 
    float_type y=this->func->value_with_derivatives(x, yprime, yprime2); 
    if(yprime) (*yprime)*=yscale; 
    if(yprime2) (*yprime2)*=yscale; 
    return y*yscale; 
  }

protected:
  c2_scaled_function_p<float_type>() : func() {} 
  /// \brief the scaling factor for the function
  const c2_const_ptr<float_type> func;
  float_type yscale;
};

/// \brief A container into which any other c2_function can be dropped.
/// \ingroup containers
/// It allows a function to be pre-evaluated at a point, 
/// and used at multiple places 
/// in an expression
/// efficiently. If it is re-evaluated at the previous point, 
/// it returns the remembered values;
/// otherwise, it re-evauates the function at the new point.
///
template <typename float_type=double> class c2_cached_function_p 
  : public c2_function<float_type> {
public:
  /// \brief construct the container
  ///
  /// \param f the function to be cached
  c2_cached_function_p(const c2_function<float_type> &f) 
  : c2_function<float_type>(),
    func(f), init(false)  {}
  /// \copydoc c2_function::value_with_derivatives
  ///
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */ 
  {
    if(!init || x != x0) {
      y=this->func->value_with_derivatives(x, &yp, &ypp);
      x0=x;
      init=true;
    }
    if(yprime) *yprime=yp;
    if(yprime2) *yprime2=ypp;
    return y; 
  }

protected:
  c2_cached_function_p() : func() {}
  const c2_const_ptr<float_type> func;
  mutable bool init;
  mutable float_type x0, y, yp, ypp;
        
};

template <typename float_type=double> class  c2_composed_function_p 
  : public c2_binary_function<float_type> {
public:
  c2_composed_function_p(const c2_function<float_type> &outer, 
                         const c2_function<float_type> &inner) : 
    c2_binary_function<float_type>(combine, outer, inner) { 
    this->set_domain(inner.xmin(), inner.xmax()); }
  /// \brief Create a stub just for the combiner to avoid statics. 
  c2_composed_function_p() : c2_binary_function<float_type>(combine) {} 

  /// \brief execute math necessary to do composition
  static float_type combine(const c2_function<float_type> &left, 
                            const c2_function<float_type> &right, 
                            float_type x, float_type *yprime, 
                            float_type *yprime2) /* throw(c2_exception) */
  {
    float_type y0, y1;
    if(yprime || yprime2) {
      float_type yp0, ypp0, yp1, ypp1;
      y0=right.value_with_derivatives(x, &yp0, &ypp0);
      y1=left.value_with_derivatives(y0, &yp1, &ypp1);
      if(yprime) *yprime=yp1*yp0;
      if(yprime2) *yprime2=ypp0*yp1+yp0*yp0*ypp1;
    } else {
      y0=right(x);
      y1=left(y0);
    }
    return y1;
  }        
};

template <typename float_type=double> class c2_sum_p 
  : public c2_binary_function<float_type> {
public:        
  /// \brief construct \a left + \a right
  /// \param left the left function 
  /// \param right the right function
  c2_sum_p(const c2_function<float_type> &left, 
           const c2_function<float_type> &right) 
    : c2_binary_function<float_type>(combine, left, right) {}
  /// \brief Create a stub just for the combiner to avoid statics. 
  c2_sum_p() : c2_binary_function<float_type>(combine) {} ; 

  /// \brief execute math necessary to do addition
  static float_type combine(const c2_function<float_type> &left, 
                            const c2_function<float_type> &right, 
                            float_type x, float_type *yprime, 
                            float_type *yprime2) /* throw(c2_exception) */
  {
    float_type y0, y1;
    if(yprime || yprime2) {
      float_type yp0, ypp0, yp1, ypp1;
      y0=left.value_with_derivatives(x, &yp0, &ypp0);
      y1=right.value_with_derivatives(x, &yp1, &ypp1);
      if(yprime) *yprime=yp0+yp1;
      if(yprime2) *yprime2=ypp0+ypp1;
    } else {
      y0=left(x);
      y1=right(x);
    }
    return y0+y1;
  }
};

template <typename float_type=double> class c2_diff_p 
  : public c2_binary_function<float_type> {
public:        
  /// \brief construct \a left - \a right
  /// \param left the left function 
  /// \param right the right function
  c2_diff_p(const c2_function<float_type> &left, 
            const c2_function<float_type> &right) 
    : c2_binary_function<float_type>(combine, left, right) {}
  /// \brief Create a stub just for the combiner to avoid statics. 
  c2_diff_p() : c2_binary_function<float_type>(combine) {} ; 

  /// \brief execute math necessary to do subtraction
  static float_type combine(const c2_function<float_type> &left, 
                            const c2_function<float_type> &right, 
                            float_type x, float_type *yprime, 
                            float_type *yprime2) /* throw(c2_exception) */
  {
    float_type y0, y1;
    if(yprime || yprime2) {
      float_type yp0, ypp0, yp1, ypp1;
      y0=left.value_with_derivatives(x, &yp0, &ypp0);
      y1=right.value_with_derivatives(x, &yp1, &ypp1);
      if(yprime) *yprime=yp0-yp1;
      if(yprime2) *yprime2=ypp0-ypp1;
    } else {
      y0=left(x);
      y1=right(x);
    }
    return y0-y1;
  }
};


/// \brief create a c2_function which is the product of two other c2_functions.
/// \ingroup arithmetic_functions
/// This should always be constructed using c2_function::operator*()
template <typename float_type=double> class c2_product_p 
  : public c2_binary_function<float_type> {
public:        
  /// \brief construct \a left * \a right
  /// \param left the left function 
  /// \param right the right function
  c2_product_p(const c2_function<float_type> &left, 
               const c2_function<float_type> &right) 
    : c2_binary_function<float_type>(combine, left, right) {}
  /// \brief Create a stub just for the combiner to avoid statics. 
  c2_product_p() : c2_binary_function<float_type>(combine) {} ; 
  
  /// \brief execute math necessary to do multiplication
  static float_type combine(const c2_function<float_type> &left, 
                            const c2_function<float_type> &right, 
                            float_type x, float_type *yprime, 
                            float_type *yprime2)  /* throw(c2_exception) */
    {
      float_type y0, y1;
      if(yprime || yprime2) {
        float_type yp0, ypp0, yp1, ypp1;
        y0=left.value_with_derivatives(x, &yp0, &ypp0);
        y1=right.value_with_derivatives(x, &yp1, &ypp1);
        if(yprime) *yprime=y1*yp0+y0*yp1;
        if(yprime2) *yprime2=ypp0*y1+2.0*yp0*yp1+ypp1*y0;
      } else {
        y0=left(x);
        y1=right(x);
      }
      return y0*y1;
    }
};


/// \brief create a c2_function which is the ratio of two other c2_functions.
/// \ingroup arithmetic_functions
/// This should always be constructed using c2_function::operator/()
template <typename float_type=double> class c2_ratio_p 
  : public c2_binary_function<float_type> {
public:        
  /// \brief construct \a left / \a right
  /// \param left the left function 
  /// \param right the right function
  c2_ratio_p(const c2_function<float_type> &left, 
             const c2_function<float_type> &right) 
    : c2_binary_function<float_type>(combine, left, right) {}
  /// \brief Create a stub just for the combiner to avoid statics. 
  c2_ratio_p() : c2_binary_function<float_type>(combine) {} ; 
        
  /// \brief execute math necessary to do division
  static float_type combine(const c2_function<float_type> &left, 
                            const c2_function<float_type> &right, 
                            float_type x, float_type *yprime, 
                            float_type *yprime2) /* throw(c2_exception) */
  {
    float_type y0, y1;
    if(yprime || yprime2) {
      float_type yp0, ypp0, yp1, ypp1;
      y0=left.value_with_derivatives(x, &yp0, &ypp0);
      y1=right.value_with_derivatives(x, &yp1, &ypp1);
      if(yprime) *yprime=(yp0*y1-y0*yp1)/(y1*y1); // first deriv of ratio
      if(yprime2) *yprime2=(y1*y1*ypp0+y0*(2*yp1*yp1-y1*ypp1)-2*y1*yp0*yp1)
                    /(y1*y1*y1); 
    } else {
      y0=left(x);
      y1=right(x);
    }
    return y0/y1;
  }
};

/// \brief a c2_function which is constant 
/// \ingroup parametric_functions
///
/// The factory function c2_factory::constant() creates *new c2_constant_p()
template <typename float_type> class c2_constant_p 
  : public c2_function<float_type> {
public:
  c2_constant_p(float_type x) : c2_function<float_type>(), value(x) {}
  void reset(float_type val) { value=val; }
  virtual float_type value_with_derivatives(
          float_type, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */ 
  { if(yprime) *yprime=0; if(yprime2) *yprime2=0; return value; }
        
private:
  float_type value;
};

/// \brief a transformation of a coordinate, including an inverse
/// \ingroup transforms
template <typename float_type> class c2_transformation {
public:
  /// \brief initialize all our function pointers
  /// \param transformed true if this function is not the identity
  /// \param xin input X transform
  /// \param xinp input X transform derivative
  /// \param xinpp input X transform second derivative
  /// \param xout output X transform, which MUST be the inverse of \a xin
  c2_transformation(bool transformed,
                    float_type (*xin)(float_type), 
                    float_type (*xinp)(float_type), 
                    float_type (*xinpp)(float_type), 
                    float_type (*xout)(float_type)
                    ) :
    fTransformed(transformed), fHasStaticTransforms(true),
    pIn(xin), pInPrime(xinp), pInDPrime(xinpp), pOut(xout) { }

  /// \brief initialize all our function pointers so that only the (overridden) 
  /// virtual functions can be called without an error
  /// \param transformed true if this function is nonlinear
  c2_transformation(bool transformed) :
    fTransformed(transformed), fHasStaticTransforms(false),
    pIn(report_error), pInPrime(report_error), pInDPrime(report_error), 
    pOut(report_error) { }
  /// \brief the destructor
  virtual ~c2_transformation() { }
  /// \brief flag to indicate if this transform is not the identity
  const bool fTransformed;
  /// \brief flag to indicate if the static function pointers can 
  /// be used for efficiency
  const bool fHasStaticTransforms;

  /// \note the pointers to functions allow highly optimized access when static 
  /// functions are available.
  /// They are only used inside value_with_derivatives(), 
  /// which is assumed to be the most critical routine.
  /// \brief non-virtual pointer to input X transform
  float_type (* const pIn)(float_type);
  /// \brief non-virtual pointer to input X transform derivative
  float_type (* const pInPrime)(float_type);
  /// \brief non-virtual pointer to input X transform second derivative
  float_type (* const pInDPrime)(float_type);
  /// \brief non-virtual pointer to output X transform 
  float_type (* const pOut)(float_type);
        
  /// \brief virtual input X transform
  virtual float_type fIn(float_type x) const { return pIn(x); }
  /// \brief virtual input X transform derivative
  virtual float_type fInPrime(float_type x) const { return pInPrime(x); }
  /// \brief virtual input X transform second derivative
  virtual float_type fInDPrime(float_type x) const { return pInDPrime(x); }
  /// \brief virtual output X transform 
  virtual float_type fOut(float_type x) const { return pOut(x); }

protected:
  /// \brief utility function for unimplemented conversion
  static float_type report_error(float_type x)  { 
    throw c2_exception("use of improperly constructed axis transform"); 
    return x; }
  /// \brief utility function f(x)=x useful in axis transforms
  static float_type ident(float_type x)  { return x; }
  /// \brief utility function f(x)=1 useful in axis transforms
  static float_type one(float_type)  { return 1; }
  /// \brief utility function f(x)=0 useful in axis transforms
  static float_type zero(float_type)  { return 0; }
  /// \brief utility function f(x)=1/x useful in axis transforms
  static float_type recip(float_type x)  { return 1.0/x; }
  /// \brief utility function f(x)=-1/x**2 useful in axis transforms
  static float_type recip_prime(float_type x)  { return -1/(x*x); }
  /// \brief utility function f(x)=2/x**3 useful in axis transforms
  static float_type recip_prime2(float_type x)  { return 2/(x*x*x); }
};

/// \brief the identity transform
/// \ingroup transforms
template <typename float_type> class c2_transformation_linear 
  : public c2_transformation<float_type> {
public:
  /// \brief constructor
  c2_transformation_linear() : c2_transformation<float_type>(
     false, this->ident, 
     this->one, this->zero, 
     this->ident) { }
  /// \brief destructor
  virtual ~c2_transformation_linear() { }
};
/// \brief log axis transform
/// \ingroup transforms
template <typename float_type> class c2_transformation_log 
  : public c2_transformation<float_type> {
public:
  /// \brief constructor
  c2_transformation_log() : c2_transformation<float_type>(
     true, std::log, this->recip, 
     this->recip_prime, std::exp) { }
  /// \brief destructor
  virtual ~c2_transformation_log() { }
};
/// \brief reciprocal axis transform
/// \ingroup transforms
template <typename float_type> class c2_transformation_recip 
  : public c2_transformation<float_type> {
public:
  /// \brief constructor
  c2_transformation_recip() : c2_transformation<float_type>(
     true, this->recip, 
     this->recip_prime, 
     this->recip_prime2, this->recip) { }
  /// \brief destructor
  virtual ~c2_transformation_recip() { }
};

/// \brief a transformation of a function in and out of a coordinate space, 
/// using 2 c2_transformations
///
/// This class is a container for two axis transforms, 
///  but also provides the critical evaluate()
/// function which converts a result in internal 
/// coordinates (with derivatives) into the 
/// external representation
/// \ingroup transforms
template <typename float_type> 
class c2_function_transformation {
public:
  /// \brief construct this from two c2_transformation instances
  /// \param xx the X axis transform
  /// \param yy the Y axis transform  
  c2_function_transformation(
     const c2_transformation<float_type> &xx, 
     const c2_transformation<float_type> &yy) :
    isIdentity(!(xx.fTransformed || yy.fTransformed)), X(xx), Y(yy) { }
  /// \brief destructor
  virtual ~c2_function_transformation() { delete &X; delete &Y; }
  virtual float_type evaluate(
          float_type xraw, 
          float_type y, float_type yp0, float_type ypp0,
          float_type *yprime, float_type *yprime2) const;
  const bool isIdentity;
  /// \brief the X axis transform
  const c2_transformation<float_type> &X;
  /// \brief the Y axis transform 
  const c2_transformation<float_type> &Y;
};

/// \brief a transformation of a function in and out of lin-lin space
///
/// \ingroup transforms
template <typename float_type> class c2_lin_lin_function_transformation : 
  public c2_function_transformation<float_type> {
public:
  c2_lin_lin_function_transformation() : 
    c2_function_transformation<float_type>(
       *new c2_transformation_linear<float_type>, 
       *new c2_transformation_linear<float_type>
                                           ) { }
  virtual ~c2_lin_lin_function_transformation() { }
};

/// \brief a transformation of a function in and out of log-log space
///
/// \ingroup transforms
template <typename float_type> class c2_log_log_function_transformation : 
        public c2_function_transformation<float_type> {
public:
        c2_log_log_function_transformation() : 
                c2_function_transformation<float_type>(
                        *new c2_transformation_log<float_type>, 
                        *new c2_transformation_log<float_type>
                ) { }
        virtual ~c2_log_log_function_transformation() { }
};

/// \brief a transformation of a function in and out of lin-log space
///
/// \ingroup transforms
template <typename float_type> class c2_lin_log_function_transformation : 
        public c2_function_transformation<float_type> {
public:
        c2_lin_log_function_transformation() : 
                c2_function_transformation<float_type>(
                        *new c2_transformation_linear<float_type>, 
                        *new c2_transformation_log<float_type>
                ) { }
        virtual ~c2_lin_log_function_transformation() { }
};

/// \brief a transformation of a function in and out of log-lin space
///
/// \ingroup transforms
template <typename float_type> class c2_log_lin_function_transformation : 
        public c2_function_transformation<float_type> {
public:
        c2_log_lin_function_transformation() : 
                c2_function_transformation<float_type>(
                        *new c2_transformation_log<float_type>, 
                        *new c2_transformation_linear<float_type>
                ) { }
        virtual ~c2_log_lin_function_transformation() { }
};

/// \brief a transformation of a function in and out of Arrhenius 
/// (1/x vs. log(y)) space
///
/// \ingroup transforms
template <typename float_type> class c2_arrhenius_function_transformation : 
        public c2_function_transformation<float_type> {
public:
        c2_arrhenius_function_transformation() : 
                c2_function_transformation<float_type>(
                        *new c2_transformation_recip<float_type>, 
                        *new c2_transformation_log<float_type>
                ) { }
        virtual ~c2_arrhenius_function_transformation() { }
};

/**
    \brief create a cubic spline interpolation of a set of (x,y) pairs
        \ingroup interpolators
    This is one of the main reasons for c2_function objects to exist.

    It provides support for cubic spline interpolation of data 
    provides from tables 
    of \a x, \a y pairs.
    It supports automatic, transparent linearization of the data 
    before storing in 
    its tables (through
    subclasses such as 
    log_lin_interpolating_function, lin_log_interpolating_function, and
    log_log_interpolating_function) to permit very high 
    accuracy representations of 
    data which have a suitable
    structure.  It provides utility functions 
    LinearInterpolatingGrid() and LogLogInterpolatingGrid()
    to create grids for mapping other functions onto a arithmetic 
    or geometric grid.

    In its simplest form, an untransformed cubic spline of a data set, 
    using natural boundary conditions
    (vanishing second derivative), is created as: \n
    \code 
        c2_ptr<double> c2p;
        c2_factory<double> c2;
    std::vector<double> xvals(10), yvals(10); 
    // < fill in xvals and yvals >
    c2p myfunc=c2.interpolating_function().load(xvals, yvals,true,0,true,0);
    // and it can be evaluated at a point for its value only by: 
    double y=myfunc(x); 
    // or it can be evaluated with its derivatives by
    double yprime, yprime2; 
    double y=myfunc(x,&yprime, &yprime2);
    \endcode
        
    The factory function c2_factory::interpolating_function() 
    creates *new interpolating_function_p()
*/

template <typename float_type=double> class interpolating_function_p  
  : public c2_function<float_type> {
public:
  /// \brief an empty linear-linear cubic-spline interpolating_function_p
  ///
  interpolating_function_p() : c2_function<float_type>(), 
    fTransform(*new  c2_lin_lin_function_transformation<float_type>) { } 

  /// \brief an empty cubic-spline interpolating_function_p with a 
  /// specific transform
  ///
  interpolating_function_p(const c2_function_transformation<float_type> &
                           transform) 
    : c2_function<float_type>(), 
      fTransform(transform) { } 

  /// \brief do the dirty work of constructing the spline from a function. 
  /// \param x the list of abscissas.  Must be either strictly 
  /// increasing or strictly decreasing.
  /// Strictly increasing is preferred, as less memory is used since 
  /// a copy is not 
  /// required for the sampling grid.
  /// \param f the list of function values.
  /// \param lowerSlopeNatural if true, set y''(first point)=0, 
  /// otherwise compute it 
  /// from \a lowerSope
  /// \param lowerSlope derivative of the function at the lower bound, 
  /// used only 
  /// if \a lowerSlopeNatural is false
  /// \param upperSlopeNatural if true, set y''(last point)=0, 
  /// otherwise compute 
  /// it from \a upperSope
  /// \param upperSlope derivative of the function at the upper bound, 
  /// used only 
  /// if \a upperSlopeNatural is false
  /// \param splined if true (default), use cubic spline, 
  /// if false, use linear interpolation.
  /// \return the same interpolating function, filled
  interpolating_function_p<float_type> & load(
      const std::vector<float_type> &x, 
      const std::vector<float_type> &f, 
      bool lowerSlopeNatural, float_type lowerSlope, 
      bool upperSlopeNatural, float_type upperSlope, bool splined=true
        ) /* throw(c2_exception) */;

  /// \brief do the dirty work of constructing the spline from a function. 
  /// \param data std::vector of std::pairs of x,y. 
  /// Will be sorted into x increasing order in place.
  /// \param lowerSlopeNatural if true, set y''(first point)=0, 
  /// otherwise compute it from \a lowerSope
  /// \param lowerSlope derivative of the function at the lower bound, 
  /// used only if \a lowerSlopeNatural is false
  /// \param upperSlopeNatural if true, set y''(last point)=0, 
  /// otherwise compute 
  /// it from \a upperSope
  /// \param upperSlope derivative of the function at the upper bound, 
  /// used only if \a upperSlopeNatural is false
  /// \param splined if true (default), use cubic spline, 
  /// if false, use linear interpolation.
  /// \return the same interpolating function, filled
  interpolating_function_p<float_type> & load_pairs(
     std::vector<std::pair<float_type, float_type> > &data,  
     bool lowerSlopeNatural, float_type lowerSlope, 
     bool upperSlopeNatural, float_type upperSlope, bool splined=true
                                                    ) /* throw(c2_exception) */;

  /// \brief do the dirty work of constructing the spline from a function.
  /// \param func a function without any requirement of valid derivatives 
  /// to sample 
  /// into an interpolating function.
  /// Very probably a c2_classic_function.
  /// \param amin the lower bound of the region to sample
  /// \param amax the upper bound of the region to sample
  /// \param abs_tol the maximum absolute error permitted when 
  /// linearly interpolating the points.
  /// the real error will be much smaller, 
  /// since this uses cubic splines at the end.
  /// \param rel_tol the maximum relative error 
  /// permitted when linearly interpolating the points.
  /// the real error will be much smaller, 
  /// since this uses cubic splines at the end.
  /// \param lowerSlopeNatural if true, set y'(first point) 
  /// from 3-point parabola, 
  /// otherwise compute it from \a lowerSope
  /// \param lowerSlope derivative of the function at the lower bound, 
  /// used only if \a lowerSlopeNatural is false
  /// \param upperSlopeNatural if true, 
  /// set y'(last point) from 3-point parabola, 
  /// otherwise compute it from \a upperSope
  /// \param upperSlope derivative of the function at the upper bound, 
  /// used only 
  /// if \a upperSlopeNatural is false
  /// \return the same interpolating function, filled
  /// \note If the interpolator being filled has a log vertical axis, 
  /// put the desired 
  /// relative error in 
  /// \a abs_tol, and 0 in \a rel_tol since the absolute error 
  /// on the log of a function 
  /// is the relative error
  /// on the function itself. 
  interpolating_function_p<float_type> &  
  sample_function(const c2_function<float_type> &func, 
                  float_type amin, float_type amax, 
                  float_type abs_tol, float_type rel_tol,
                  bool lowerSlopeNatural, float_type lowerSlope, 
                  bool upperSlopeNatural, float_type upperSlope
                  ) /* throw(c2_exception) */;
        
  /// \brief initialize from a grid of points and a 
  /// c2_function (un-normalized) to an 
  /// interpolator which, when evaluated with a 
  /// uniform random variate on [0,1] returns 
  /// random numbers
  /// distributed as the input function.
  /// \see  \ref random_subsec "Arbitrary random generation"
  /// inverse_integrated_density starts with a probability density  
  /// std::vector, 
  /// generates the integral, 
  /// and generates an interpolating_function_p  of the inverse function which, 
  /// when evaluated using a uniform random on [0,1] returns values
  /// with a density distribution equal to the input distribution
  /// If the data are passed in reverse order (large X first), 
  /// the integral is carried out 
  /// from the big end.
  /// \param bincenters the positions at which to sample the 
  /// function \a binheights
  /// \param binheights a function which describes the density 
  /// of the random number 
  /// distribution to be produced.
  /// \return an initialized interpolator, which 
  /// if evaluated randomly with a uniform variate on [0,1] produces numbers
  /// distributed according to \a binheights
  interpolating_function_p<float_type> & load_random_generator_function(
     const std::vector<float_type> &bincenters, 
     const c2_function<float_type> &binheights)
    /* throw(c2_exception) */;

  interpolating_function_p<float_type> & load_random_generator_bins(
    const std::vector<float_type> &bins, 
    const std::vector<float_type> &binheights, 
    bool splined=true)
    /* throw(c2_exception) */;
        
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */;
        
  /// \brief destructor
  virtual ~interpolating_function_p() { delete &fTransform; } 
        
  virtual interpolating_function_p<float_type> &clone() 
    const /* throw(c2_exception) */
  {        return *new interpolating_function_p<float_type>(); }
        
  void get_data(std::vector<float_type> &xvals, 
                std::vector<float_type> &yvals) const /* throw() */ ;
        
  void get_internal_data(
       std::vector<float_type> &xvals, 
       std::vector<float_type> &yvals, 
       std::vector<float_type> &y2vals) const 
  { xvals=X; yvals=F; y2vals=y2; }

  void set_lower_extrapolation(float_type bound);
  void set_upper_extrapolation(float_type bound);
        
  interpolating_function_p <float_type> & 
  unary_operator(const c2_function<float_type> &source) const;

  interpolating_function_p <float_type> & 
  binary_operator(const c2_function<float_type> &rhs,
                  const c2_binary_function<float_type> *combining_stub
                  ) const;
  interpolating_function_p <float_type> & 
  add_pointwise (const c2_function<float_type> &rhs) const { 
    return binary_operator(rhs, new c2_sum_p<float_type>()); }
  interpolating_function_p <float_type> & 
  subtract_pointwise (const c2_function<float_type> &rhs) const {
                return binary_operator(rhs, new c2_diff_p<float_type>()); }
  interpolating_function_p <float_type> & 
  multiply_pointwise (const c2_function<float_type> &rhs) const { 
                return binary_operator(rhs, new c2_product_p<float_type>()); }
  interpolating_function_p <float_type> & 
  divide_pointwise (const c2_function<float_type> &rhs) const { 
                return binary_operator(rhs, new c2_ratio_p<float_type>()); }
  void clone_data(const interpolating_function_p <float_type> &rhs) {
    Xraw=rhs.Xraw; X=rhs.X; F=rhs.F; y2=rhs.y2;
    set_sampling_grid_pointer(Xraw);
  }

  const c2_function_transformation<float_type> &fTransform;
        
protected:    
  /// \brief create the spline coefficients
  void spline(
              bool lowerSlopeNatural, float_type lowerSlope, 
              bool upperSlopeNatural, float_type upperSlope
              ) /* throw(c2_exception) */;

  static bool comp_pair(std::pair<float_type,float_type> const &i, 
                        std::pair<float_type,float_type> const &j) 
  {return i.first<j.first;}

  std::vector<float_type> Xraw, X, F, y2;
  c2_const_ptr<float_type> sampler_function;
  bool xInverted;
  mutable size_t lastKLow;
};

/// \brief A spline with X transformed into log space.  
/// \ingroup interpolators
///
template <typename float_type=double> class log_lin_interpolating_function_p 
  : public interpolating_function_p <float_type> {
public:
  /// \brief an empty log-linear cubic-spline interpolating_function_p
  ///
  log_lin_interpolating_function_p() : 
    interpolating_function_p<float_type>(
       *new c2_log_lin_function_transformation<float_type>)
  { }
  virtual interpolating_function_p<float_type> &clone() 
    const /* throw(c2_exception) */
  { return *new log_lin_interpolating_function_p<float_type>(); }
};


/// \brief A spline with Y transformed into log space.  
/// \ingroup interpolators
/// Most useful for functions looking like y=exp(x)
///
template <typename float_type=double> class lin_log_interpolating_function_p 
  : public interpolating_function_p <float_type> {
public:
    /// \brief an empty linear-log cubic-spline interpolating_function_p
    ///
  lin_log_interpolating_function_p() 
  : interpolating_function_p<float_type>(
      *new c2_lin_log_function_transformation<float_type>)
  { } 
  virtual interpolating_function_p<float_type> &clone() 
    const /* throw(c2_exception) */
  { return *new lin_log_interpolating_function_p<float_type>(); }
};


/// \brief A spline with X and Y transformed into log space.  
/// \ingroup interpolators
/// Most useful for functions looking like y=x^n or any other 
/// function with a huge X and Y dynamic range.
///
template <typename float_type=double> class log_log_interpolating_function_p 
  : public interpolating_function_p <float_type> {
public:        
    /// \brief an empty log-log cubic-spline interpolating_function_p
    ///
  log_log_interpolating_function_p() : 
    interpolating_function_p<float_type>(
      *new c2_log_log_function_transformation<float_type>)
  { } 
  virtual interpolating_function_p<float_type> &clone() 
    const /* throw(c2_exception) */
  { return *new log_log_interpolating_function_p<float_type>(); }
};


/// \brief A spline with X in reciprocal space and Y transformed in log space.  
/// \ingroup interpolators
/// Most useful for thermodynamic types of data where Y is roughly A*exp(-B/x). 
/// Typical examples are reaction rate data, and thermistor calibration data.
///
template <typename float_type=double> class arrhenius_interpolating_function_p 
  : public interpolating_function_p <float_type> {
public:
    /// \brief an empty arrhenius cubic-spline interpolating_function_p
    ///
  arrhenius_interpolating_function_p() 
  : interpolating_function_p<float_type>(
      *new c2_arrhenius_function_transformation<float_type>)
  { } 
  virtual interpolating_function_p<float_type> &clone() 
    const /* throw(c2_exception) */
  { return *new arrhenius_interpolating_function_p<float_type>(); }        
};

/// \brief compute sin(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::sin() creates *new c2_sin_p
template <typename float_type=double> class c2_sin_p 
  : public c2_function<float_type> {
public:
  /// \brief constructor.
  c2_sin_p() : c2_function<float_type>() {}
        
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { float_type q=std::sin(x); 
    if(yprime) *yprime=std::cos(x); 
    if(yprime2) *yprime2=-q; 
    return q; }
    
  virtual void get_sampling_grid(float_type amin, float_type amax,  
                                 std::vector<float_type> &grid) const; 
};

/// \brief compute cos(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::cos() creates *new c2_cos_p
template <typename float_type=double> class c2_cos_p 
  : public c2_sin_p<float_type> {
public:
  /// \brief constructor.
  c2_cos_p() : c2_sin_p<float_type>() {}
        
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { float_type q=std::cos(x); 
    if(yprime) *yprime=-std::sin(x); 
    if(yprime2) *yprime2=-q; 
    return q; }
};

/// \brief compute tan(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::tan() creates *new c2_tan_p
template <typename float_type=double> class c2_tan_p 
  : public c2_function<float_type> {          
public:
  /// \brief constructor.
  c2_tan_p() : c2_function<float_type>() {}
        
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  {
    float_type c=std::cos(x), ss=std::sin(x);
    float_type t=ss/c;
    float_type yp=1/(c*c);
    if(yprime) { *yprime=yp; } 
    if(yprime2){*yprime2=2*t*yp; } 
    return t; 
  }        
};

/// \brief compute log(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::log() creates *new c2_log_p
template <typename float_type=double> class c2_log_p 
  : public c2_function<float_type> {
public:
  /// \brief constructor.
  c2_log_p() : c2_function<float_type>() {}

  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { if(yprime) *yprime=1.0/x; 
    if(yprime2) *yprime2=-1.0/(x*x); 
    return std::log(x); }        
};

/// \brief compute exp(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::exp() creates *new c2_exp_p
template <typename float_type=double>  class c2_exp_p 
  : public c2_function<float_type> {
public:
  /// \brief constructor.
  c2_exp_p() : c2_function<float_type>() {}

  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { float_type q=std::exp(x); 
    if(yprime) *yprime=q; 
    if(yprime2) *yprime2=q; 
    return q; }
};

/// \brief compute sqrt(x) with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::sqrt() creates *new c2_sqrt_p()
template <typename float_type=double> class c2_sqrt_p 
  : public c2_function<float_type> {
public:
  /// \brief constructor.
  c2_sqrt_p() : c2_function<float_type>() {}
        
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { float_type q=std::sqrt(x); 
    if(yprime) *yprime=0.5/q; 
    if(yprime2) *yprime2=-0.25/(x*q); 
    return q; }
};

/// \brief compute scale/x with its derivatives.
/// \ingroup parametric_functions
///
/// The factory function c2_factory::recip() creates *new c2_recip_p
template <typename float_type=double> class c2_recip_p 
  : public c2_function<float_type> {
public:
  /// \brief constructor.
  c2_recip_p(float_type scale) : c2_function<float_type>(), rscale(scale) {}
        
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { 
    float_type q=1.0/x; 
    float_type y=rscale*q;
    if(yprime) *yprime=-y*q; 
    if(yprime2) *yprime2=2*y*q*q; 
    return y; 
  }
  /// \brief reset the scale factor
  /// \param scale the new numerator
  void reset(float_type scale) { rscale=scale; } 
private:
  float_type rscale;
};

/// \brief compute x with its derivatives.
/// \ingroup math_functions
///
/// The factory function c2_factory::identity() creates *new c2_identity_p
template <typename float_type=double> class c2_identity_p 
  : public c2_function<float_type> {
public:
  /// \brief constructor.
  c2_identity_p() : c2_function<float_type>() {}
        
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { if(yprime) *yprime=1.0; if(yprime2) *yprime2=0; return x; }
};

/**
 \brief create a linear mapping of another function 
 \ingroup parametric_functions
 for example, given a c2_function \a f 
 \code 
 c2_function<double> &F=c2_linear<double>(1.2, 2.0, 3.0)(f);
 \endcode
 produces a new c2_function F=2.0+3.0*(\a f - 1.2) 
 
 The factory function c2_factory::linear() creates *new c2_linear_p
*/
template <typename float_type=double> class c2_linear_p 
  : public c2_function<float_type> {
public:
  /// \brief Construct the operator f=y0 + slope * (x-x0)
  /// \param x0 the x offset
  /// \param y0 the y-intercept i.e. f(x0)
  /// \param slope the slope of the mapping
  c2_linear_p(float_type x0, float_type y0, float_type slope) : 
    c2_function<float_type>(), xint(x0), intercept(y0), m(slope) {}
  /// \brief Change the slope and intercepts after construction.
  /// \param x0 the x offset
  /// \param y0 the y-intercept 
  /// \param slope the slope of the mapping
  void reset(float_type x0, float_type y0, float_type slope) 
  { xint=x0; intercept=y0; m=slope; } 
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { if(yprime) *yprime=m; 
    if(yprime2) *yprime2=0; 
    return m*(x-xint)+intercept; }
        
private:
  float_type xint, intercept, m;
protected:
  c2_linear_p() {} 
};

/**
\brief create a quadratic mapping of another function 
 \ingroup parametric_functions
 for example, given a c2_function \a f 
 \code 
 c2_function<double> &F=c2_quadratic<double>(1.2, 2.0, 3.0, 4.0)(f);
 \endcode
 produces a new c2_function F=2.0 + 3.0*(f-1.2) + 4.0*(f-1.2)^2 

 note that the parameters are overdetermined, 
 but allows the flexibility of two different representations

 The factory function c2_factory::quadratic() creates *new c2_quadratic_p
 */
template <typename float_type=double> class c2_quadratic_p 
  : public c2_function<float_type> {
public:
    /// \brief Construct the operator
    /// \param x0 the center around which the powers are computed
    /// \param y0 the value of the function at \a x = \a x0
    /// \param xcoef the scale on the (\a x - \a x0) term
    /// \param x2coef the scale on the (\a x - \a x0)^2 term
  c2_quadratic_p(float_type x0, float_type y0, 
                 float_type xcoef, float_type x2coef) : 
    c2_function<float_type>(), intercept(y0), center(x0), a(x2coef), b(xcoef) {}
    /// \brief Modify the coefficients after construction
    /// \param x0 the new center around which the powers are computed
    /// \param y0 the new value of the function at \a x = \a x0
    /// \param xcoef the new scale on the (\a x - \a x0) term
    /// \param x2coef the new scale on the (\a x - \a x0)^2 term    
  void reset(float_type x0, float_type y0, float_type xcoef,
             float_type x2coef) { intercept=y0; center=x0; a=x2coef; b=xcoef; } 
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { float_type dx=x-center; 
    if(yprime) *yprime=2*a*dx+b; 
    if(yprime2) *yprime2=2*a; 
    return a*dx*dx+b*dx+intercept; }
        
private:
  float_type intercept, center, a, b;
protected:
  c2_quadratic_p() {} 
};

/**
\brief create a power law mapping of another function 
 \ingroup parametric_functions
 for example, given a c2_function \a f 
 \code 
 c2_power_law_p<double> PLaw(1.2, 2.5);
 c2_composed_function_p<double> &F=PLaw(f); 
 \endcode
 produces a new c2_function F=1.2 * f^2.5 
 
 The factory function c2_factory::power_law() creates *new c2_power_law_p
 */
template <typename float_type=double> class c2_power_law_p 
  : public c2_function<float_type> {
public:
    /// \brief Construct the operator
    /// \param scale the multipler
    /// \param power the exponent
  c2_power_law_p(float_type scale, float_type power) :
    c2_function<float_type>(), a(scale), b(power) {}
    /// \brief Modify the mapping after construction
    /// \param scale the new multipler
    /// \param power the new exponent
  void reset(float_type scale, float_type power) { a=scale; b=power; }
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */
  { float_type q=a*std::pow(x,b-2); 
    if(yprime) *yprime=b*q*x; 
    if(yprime2) *yprime2=b*(b-1)*q; 
    return q*x*x; }
        
private:
  float_type a, b;
protected:
  c2_power_law_p() {} 
};

/**
\brief create the formal inverse function of another function 
 \ingroup containers
 for example, given a c2_function \a f 
 \code 
 c2_inverse_function<double> inv(f);
 a=f(x);
 x1=inv(a);
 \endcode
 will return x1=x to machine precision.  The important part of this
 is that the resulting function is a first-class c2_function, so it knows its 
 derivatives, too, unlike the case of a simple root-finding inverse.  This means
 it can be integrated (for example) quite efficiently.

 \see \ref combined_inversion_hinting_sampling

*/
template <typename float_type=double> class c2_inverse_function_p 
  : public c2_function<float_type> {
public:
  /// \brief Construct the operator
  /// \param source the function to be inverted
  c2_inverse_function_p(const c2_function<float_type> &source);  
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw(c2_exception) */;
    
  /// \brief give the function a hint as to where to look for its inverse
  /// \param hint the likely value of the inverse, 
  /// which defaults to whatever the evaluation returned.
    void set_start_hint(float_type hint) const { start_hint=hint; }
    
    virtual float_type get_start_hint(float_type x) const 
  { return hinting_function.valid()? hinting_function(x) : start_hint; } 

  /// \brief set or unset the approximate function used to start the root finder
  /// \anchor set_hinting_function_discussion
  /// A hinting function is mostly useful if the evaluation of this inverse is
  /// going to be carried out in very non-local order, 
  /// so the root finder has to start over 
  /// for each step.  If most evaluations are going to be made 
  /// in fairly localized clusters 
  /// (scanning through the function, for example), the default mechanism used 
  /// (which just remembers the last point)
  /// is almost certainly faster.
  /// 
  /// Typically, the hinting function is likely to be set up by 
  /// creating the inverse function,
  /// and then adaptively sampling an interpolating function from it, 
  /// and then using the result
  /// to hint it.  Another way, if the parent function is already 
  /// an interpolating function, 
  /// is just to create a version of the parent with the x & y 
  /// coordinates reversed.
  /// 
  /// \see \ref combined_inversion_hinting_sampling
  ///
  /// \param hint_func the function that is an approximate inverse 
  /// of the parent of 
  /// this inverse_function
  void set_hinting_function(const c2_function<float_type> *hint_func) 
  { hinting_function.set_function(hint_func); }
  /// \brief set the hinting function from a pointer.
  /// 
  /// See \ref set_hinting_function_discussion "discussion"
  /// \param hint_func the container holding the function
  void set_hinting_function(const c2_const_ptr<float_type> hint_func) 
  { hinting_function=hint_func; }
        
protected:
  c2_inverse_function_p() {} 
  mutable float_type start_hint;
  const c2_const_ptr<float_type> func;
  c2_const_ptr<float_type> hinting_function;
};

/** 
  \brief
  An interpolating_function_p  which is the cumulative integral of a histogram.
  \ingroup interpolators
  Note than binedges should be one element longer than binheights, 
  since the lower & upper edges are specified. 
  Note that this is a malformed spline, 
  since the second derivatives are all zero, 
  so it has less continuity.
  Also, note that the bin edges can be given in backwards order to generate the 
  reversed accumulation (starting at the high end) 
*/

template <typename float_type=double>  class accumulated_histogram 
  : public interpolating_function_p <float_type> {
public:
  /// \brief Construct the integrated histogram
  /// \param binedges the edges of the bins in \a binheights.  
  /// It should have one more element than \a binheights
  /// \param binheights the number of counts in each bin.
  /// \param normalize if true, normalize integral to 1
  /// \param inverse_function if true, drop zero channels, 
  /// and return inverse function for random generation
  /// \param drop_zeros eliminate null bins before integrating, 
  /// so integral is strictly monotonic. 
  accumulated_histogram(const std::vector<float_type>binedges, 
                        const std::vector<float_type> binheights,
                        bool normalize=false, 
                        bool inverse_function=false, bool drop_zeros=true);
        
};

/// \brief create a c2_function which smoothly connects two other c2_functions.
/// \ingroup parametric_functions
/// This takes two points and generates a polynomial 
/// which matches two c2_function arguments 
/// at those two points, with two derivatives at each point, 
/// and an arbitrary value at the center of the 
/// region.  It is useful for splicing together functions 
/// over rough spots (0/0, for example).
/// 
/// If \a auto_center is true, the value at the midpoint is computed so 
/// that the resulting polynomial is
/// of order 5.  If \a auto_center is false, the 
/// value \a y1 is used at the midpoint, 
/// resulting in a 
/// polynomial of order 6.
///
/// This is usually used in conjunction 
/// with c2_piecewise_function_p to assemble an 
/// apparently seamless 
/// function from a series of segments. 
/// \see \ref piecewise_applications_subsec "Sample Applications" 
/// and \ref c2_function::adaptively_sample() "Adaptive sampling"
///
template <typename float_type=double> class c2_connector_function_p 
  : public c2_function<float_type> {
public:        
  /// \brief construct the container from two functions
  /// \param x0 the point at which to match \a f1 and its derivatives
  /// \param f0 the function on the left side to be connected 
  /// \param x2 the point at which to match \a f2 and its derivatives
  /// \param f2 the function on the right side to be connected
  /// \param auto_center if true, no midpoint value is specified.  
  /// If false, match the value \a y1 at the midpoint
  /// \param y1 the value to match at the midpoint, if \a auto_center is false
  /// \return a c2_function with domain (\a x0,\a x2) which smoothly 
  /// connects \a f0(x0) and \a f2(x2)
  c2_connector_function_p(float_type x0, const c2_function<float_type> &f0, 
                          float_type x2, const c2_function<float_type> &f2, 
                        bool auto_center, float_type y1);
  /// \brief construct the container from numerical values
  /// \param x0 the position of the left edge
  /// \param y0 the function derivative on the left boundary
  /// \param yp0 the function second derivative on the left boundary
  /// \param ypp0 the function value on the left boundary
  /// \param x2 the position of the right edge
  /// \param y2 the function derivative on the right boundary
  /// \param yp2 the function second derivative on the right boundary
  /// \param ypp2 the function value on the right boundary
  /// \param auto_center if true, no midpoint value is specified.  
  /// If false, match the value \a y1 at the midpoint
  /// \param y1 the value to match at the midpoint, if \a auto_center is false
  /// \return a c2_function with domain (\a x0,\a x2) 
  /// which smoothly connects the points described
  /// \anchor c2_connector_raw_init_docs
  c2_connector_function_p(
                float_type x0, float_type y0, float_type yp0, float_type ypp0, 
                float_type x2, float_type y2, float_type yp2, float_type ypp2, 
                bool auto_center, float_type y1);
  /// \brief construct the container from c2_fblock<float_type> objects
  /// \param fb0 the left edge
  /// \param fb2 the right edge
  /// \param auto_center if true, no midpoint value is specified.  
  /// If false, match the value \a y1 at the midpoint
  /// \param y1 the value to match at the midpoint, if \a auto_center is false
  /// \return a c2_function with domain (\a fb0.x,\a fb2.x) which smoothly 
  /// connects \a fb0 and \a fb2
  c2_connector_function_p(
                const c2_fblock<float_type> &fb0, 
                const c2_fblock<float_type> &fb2, 
                bool auto_center, float_type y1);

  /// \brief destructor
  virtual ~c2_connector_function_p();
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw (c2_exception) */;
protected:
        /// \brief fill container numerically
        void init(
                  const c2_fblock<float_type> &fb0, 
                  const c2_fblock<float_type> &fb2, 
                  bool auto_center, float_type y1);

  float_type fhinv, fy1, fa, fb, fc, fd, fe, ff;
};

/// \brief create a c2_function which is a piecewise assembly 
/// of other c2_functions.
/// \ingroup containers
/// The functions must have increasing, non-overlapping domains.  
/// Any empty space
/// between functions will be filled with a linear interpolation.
///
/// \note If you want a smooth connection, 
/// instead of the default linear interpolation,
/// create a c2_connector_function_p to bridge the gap.  
/// The linear interpolation is intended 
/// to be a barely intelligent bridge, and may never get used by anyone.
///
/// \note The creation of the container results in the 
/// creation of an explicit sampling grid.  
/// If this is used with functions with a large domain, 
/// or which generate very dense sampling grids,
/// it could eat a lot of memory.  Do not abuse this by using functions which 
/// can generate gigantic grids.
/// 
/// \see \ref piecewise_applications_subsec "Sample Applications" \n
/// c2_plugin_function_p page \n
/// c2_connector_function_p page \n
/// \ref c2_function::adaptively_sample() "Adaptive sampling"
///
template <typename float_type=double> class c2_piecewise_function_p 
  : public c2_function<float_type> {
public:        
  /// \brief construct the container
  c2_piecewise_function_p();
  /// \brief destructor
  virtual ~c2_piecewise_function_p();
  virtual float_type value_with_derivatives(
          float_type x, float_type *yprime, 
          float_type *yprime2) const /* throw (c2_exception) */;
  /// \brief append a new function to the sequence
  ///
  /// This takes a c2_function, and appends it onto the end of 
  /// the piecewise collection.
  /// The domain of the function (which MUST be set) 
  /// specifies the place it will be used in 
  /// the final function.  If the domain exactly abuts 
  /// the domain of the previous function, it
  /// will be directly attached.  If there is a gap, the gap will be filled 
  /// in by linear interpolation.
  /// \param func a c2_function with a defined domain to be 
  /// appended to the collection
  void append_function(const c2_function<float_type> &func) 
    /* throw (c2_exception) */;
protected:
  std::vector<c2_const_ptr<float_type> > functions;
  mutable int lastKLow;
};

#include "c2_function.icc"

#endif 

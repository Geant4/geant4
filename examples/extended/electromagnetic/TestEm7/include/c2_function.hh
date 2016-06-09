/**
 *  \file 
 *  \brief Provides the headers for the general c2_function algebra which supports 
 *  fast, flexible operations on piecewise-twice-differentiable functions
 *
 *  \author Created by R. A. Weller and Marcus H. Mendenhall on 7/9/05.
 *  \author Copyright 2005 __Vanderbilt University__. All rights reserved.
 *
 * 	\version c2_function.hh,v 1.53 2007/11/12 13:58:57 marcus Exp
 */

#ifndef __has_C2Functions_c2_h
#define __has_C2Functions_c2_h 1

#include <cmath>
#include <vector>
#include <string>

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

// put these forward references here, and with a bogus typename to make swig happy.
template <typename float_type> class c2_composed_function;
template <typename float_type> class c2_sum;
template <typename float_type> class c2_diff;
template <typename float_type> class c2_product;
template <typename float_type> class c2_ratio;

/**
 \brief the parent class for all c2_functions.

  c2_functions know their value, first, and second derivative at almost every point.
  They can be efficiently combined with binary operators, via c2_binary_function, 
  composed via c2_composed_function,
  have their roots found via find_root(),
  and be adaptively integrated via partial_integrals() or integral().
  They also can carry information with them about how to find 'interesting' points on the function.  
  This information is set with set_sampling_grid() and extracted with get_sampling_grid().

  Particularly important subclasses are the interpolating functions classes,
    interpolating_function , lin_log_interpolating_function, log_lin_interpolating_function, 
    log_log_interpolating_function, and arrhenius_interpolating_function,
    as well as the template functions
    inverse_integrated_density().
 
 \warning
 The composite flavors of c2_functions (c2_sum, c2_composed_function, c2_binary_function, e.g.) make no effort to manage 
 deletion of their component functions.
 These are just container classes, and the user (along with normal automatic variable semantics) 
 is responsible for the lifetime of components.
 Inappropriate attention to this can cause massive memory leaks.  
 However, in most cases these do exactly what is intended.
 The classes will be left this way since the only other option is to use copy constructors on everything, 
 which would make this all very slow.
 
 */
template <typename float_type=double> class c2_function {
public:    
    /// \brief get versioning information for the header file
    /// \return the CVS Id string
	const std::string cvs_header_vers() const { return 
		"c2_function.hh,v 1.53 2007/11/12 13:58:57 marcus Exp";
	}
	
    /// \brief get versioning information for the source file
    /// \return the CVS Id string
	const std::string cvs_file_vers() const ;
	
public:
    /// \brief destructor
	virtual ~c2_function() { if(sampling_grid && !no_overwrite_grid) delete sampling_grid; }
	
	/// \brief get the value and derivatives. 
    ///
    /// There is required checking for null pointers on the derivatives, 
    /// and most implementations should operate faster if derivatives are not needed.
    /// \param[in] x the point at which to evaluate the function
    /// \param[out] yprime the first derivative (if pointer is non-null)
    /// \param[out] yprime2 the second derivative (if pointer is non-null)
    /// \return the value of the function
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception) =0 ; // { return 0; };
	
    /// \brief evaluate the function in the classic way, ignoring derivatives.
    /// \param x the point at which to evaluate
    /// \return the value of the function
	inline float_type operator () (float_type x) const throw(c2_exception) 
        { return value_with_derivatives(x, (float_type *)0, (float_type *)0); } 

    /// \brief compose this function outside another.
    /// \param inner the inner function
    /// \return the composed function
	c2_composed_function<float_type> & operator ()(const c2_function<float_type> &inner) const 
		{ return *new c2_composed_function<float_type>((*this), inner); }

	/// \brief get the value and derivatives. 
    ///
    /// \param[in] x the point at which to evaluate the function
    /// \param[out] yprime the first derivative (if pointer is non-null)
    /// \param[out] yprime2 the second derivative (if pointer is non-null)
    /// \return the value of the function
    inline float_type operator () (float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception) 
        { return value_with_derivatives(x, yprime, yprime2); } 
	
	/// \brief solve f(x)==value very efficiently, with explicit knowledge of derivatives of the function
	///
    /// find_root solves by iterated inverse quadratic extrapolation for a solution to f(x)=y.  It
    /// includes checks against bad convergence, so it should never be able to fail.  Unlike typical
    /// secant method or fancier Brent's method finders, this does not depend in any strong wasy on the
    /// brackets, unless the finder has to resort to successive approximations to close in on a root.
    /// Often, it is possible to make the brackets equal to the domain of the function, if there is
    /// any clue as to where the root lies, as given by the parameter \a start.  
    /// \param lower_bracket the lower bound for the search
    /// \param upper_bracket the upper bound for the search.  Function sign must be 
    /// opposite to that at \a lower_bracket
    /// \param start starting value for the search
    /// \param value the value of the function being sought (solves f(x) = \a value)
    /// \param[out] error If pointer is zero, errors raise exception. Otherwise, returns error here.
    /// \param[out] final_yprime If pointer is not zero, return derivative of function at root
    /// \param[out] final_yprime2 If pointer is not zero, return second derivative of function at root
    /// \return the position of the root.
	float_type find_root(float_type lower_bracket, float_type upper_bracket, float_type start, 
        float_type value, int *error=0, 
        float_type *final_yprime=0, float_type *final_yprime2=0 ) const throw(c2_exception) ; // solve f(x)=value

	/// \brief for points in xgrid, adaptively return Integral[f(x),{x,xgrid[i],xgrid[i+1]}] and return in vector, along with sum
    ///
	/// partial_integrals uses a method with an error O(dx**10) with full information from the derivatives,
    /// and falls back to lower order methods if informed of incomplete derivatives.
    /// It uses exact midpoint splitting of the intervals for recursion, resulting in no recomputation of the function
    /// during recursive descent at previously computed points.
    /// \param xgrid points between which to evaluate definite integrals.  
    /// \param partials if non-NULL, a vector in which to receive the partial integrals.
	/// It will automatically be sized apprpropriately, if provided, to contain \a n - 1 elements where \a n is the length of \a xgrid  
    /// \param abs_tol the absolute error bound for each segment
    /// \param rel_tol the fractional error bound for each segment.  
    /// If the error is smaller than either the relative or absolute tolerance, the integration step is finished.
    /// \param derivs number of derivatives to trust, which sets the order of the integrator.  The order 
    /// is 3*\a derivs + 4. \a derivs can be 0, 1, or 2.
    /// \param adapt if true, use recursive adaptation, otherwise do simple evaluation on the grid provided
    /// with no error checking.
    /// \param extrapolate if true, use simple Richardson extrapolation on the final 2 steps to reduce the error.
    /// \return sum of partial integrals, whcih is the definite integral from the first value in \a xgrid to the last.
	float_type partial_integrals(std::vector<float_type> xgrid, std::vector<float_type> *partials = 0,
          float_type abs_tol=1e-12, float_type rel_tol=1e-12, int derivs=2, bool adapt=true, bool extrapolate=true) const;
	
    /// \brief a fully-automated integrator which uses the information provided by the get_sampling_grid() function
    /// to figure out what to do.
	///
	/// It returns the integral of the function over the domain requested
    /// with error tolerances as specified.  It is just a front-end to partial_integrals()
    /// 
    /// \param xmin lower bound of the domain for integration
	/// \param xmax upper bound of the domain for integration
    /// \param partials if non-NULL, a vector in which to receive the partial integrals.
	/// It will automatically be sized apprpropriately, if provided, to contain \a n - 1 elements where \a n is the length of \a xgrid  
    /// \param abs_tol the absolute error bound for each segment
    /// \param rel_tol the fractional error bound for each segment.  
    /// If the error is smaller than either the relative or absolute tolerance, the integration step is finished.
    /// \param derivs number of derivatives to trust, which sets the order of the integrator.  The order 
    /// is 3*\a derivs + 4. \a derivs can be 0, 1, or 2.
    /// \param adapt if true, use recursive adaptation, otherwise do simple evaluation on the grid provided
    /// with no error checking.
    /// \param extrapolate if true, use simple Richardson extrapolation on the final 2 steps to reduce the error.
    /// \return sum of partial integrals, whcih is the definite integral from the first value in \a xgrid to the last.
	float_type integral(float_type xmin, float_type xmax, std::vector<float_type> *partials = 0,
             float_type abs_tol=1e-12, float_type rel_tol=1e-12, int derivs=2, bool adapt=true, bool extrapolate=true) const;

	/// \brief return the lower bound of the domain for this function as set by set_domain()
	inline float_type xmin() const { return fXMin; }
	/// \brief return the upper bound of the domain for this function as set by set_domain()
	inline float_type xmax() const { return fXMax; }
	/// \brief set the domain for this function.
	void set_domain(float_type xmin, float_type xmax) { fXMin=xmin; fXMax=xmax; }
		
	/// \brief this is a counter owned by the function but which can be used to monitor efficiency of algorithms.
	///
	/// It is not maintained automatically in general!  
	/// The adaptive integrator and root finder do clear it at the start and update it for performance checking.
	/// \return number of evaluations logged since last reset.
	volatile int get_evaluations() const { return evaluations; }
	/// \brief reset the counter
	void reset_evaluations()  const { evaluations=0; } // evaluations are 'invisible' to constant
	/// \brief count evaluations
	inline void increment_evaluations() const { evaluations++; }
	
	/// \brief check that a vector is monotonic, throw an exception if not, and return a flag if it is reversed
	///
	/// \param data a vector of data points which are expected to be monotonic.
	/// \param message an informative string to include in an exception if this throws c2_exception
	/// \return true if in decreasing order, false if increasing 
	bool check_monotonicity(const std::vector<float_type> &data, const char message[]) throw(c2_exception);
	
	/// \brief establish a grid of 'interesting' points on the function.
	/// 
	/// The sampling grid describes a reasonable initial set of points to look at the function.
    /// this should generally be set at a scale which is quite coarse, and sufficient for initializing
    /// adaptive integration or possibly root bracketing. For sampling a function to build a new interpolating
    /// function, one may want to refine this for accuracy.  However, interpolating_functions themselves
    /// return their original X grid by default, so refining the grid in this case might be a bad idea.
	/// \param grid a vector of abscissas.  The contents is copied into an internal vector, so the \a grid can be discarded after passingin.
	virtual void set_sampling_grid(const std::vector<float_type> &grid) throw(c2_exception); 
	
    /// \brief return the grid of 'interesting' points along this function which lie in the region requested
	///
    /// if a sampling grid is defined, work from there, otherwise return vector of (xmin, xmax)
	/// \param xmin the lower bound for which the function is to be sampled
	/// \param xmax the upper bound for which the function is to be sampled
	/// \return a new vector containing the samplng grid.
	virtual std::vector<float_type> &get_sampling_grid(float_type xmin, float_type xmax) const ;
	
	/// \brief clean up endpoints on a grid of points
	///
	/// \param[in,out] result the sampling grid with excessively closely space endpoints removed.
	/// The grid is modified in place.
	void preen_sampling_grid(std::vector<float_type> *result) const;		
	/// \brief refine a grid by splitting each interval into more intervals
	/// \param grid the grid to refine
	/// \param refinement the number of new steps for each old step
	/// \return a new sampling grid with more points.
	std::vector<float_type> & refine_sampling_grid(const std::vector<float_type> &grid, size_t refinement) const;		
	
	/// \brief create a new c2_function from this one which is normalized on the interval 
    /// \param xmin lower bound of the domain for integration
	/// \param xmax upper bound of the domain for integration
	/// \param norm the desired integral for the function over the region
	/// \return a new c2_function with the desired \a norm.
	c2_function<float_type> &normalized_function(float_type xmin, float_type xmax, float_type norm=1.0);
	/// \brief create a new c2_function from this one which is square-normalized on the interval 
    /// \param xmin lower bound of the domain for integration
	/// \param xmax upper bound of the domain for integration
	/// \param norm the desired integral for the function over the region
	/// \return a new c2_function with the desired \a norm.
	c2_function<float_type> &square_normalized_function(float_type xmin, float_type xmax, float_type norm=1.0);
	/// \brief create a new c2_function from this one which is square-normalized with the provided \a weight on the interval 
    /// \param xmin lower bound of the domain for integration
	/// \param xmax upper bound of the domain for integration
	/// \param weight a c2_function providing the weight
	/// \param norm the desired integral for the function over the region
	/// \return a new c2_function with the desired \a norm.
	c2_function<float_type> &square_normalized_function(float_type xmin, float_type xmax, const c2_function<float_type> &weight, float_type norm=1.0);

	/// \brief factory function to create a c2_sum from an regular algebraic expression.
	/// \note
	/// be very wary of ownership issues if this is used in a complex expression.
	c2_sum<float_type> &operator + (const c2_function<float_type> &rhs)  { return *new c2_sum<float_type>(*this, rhs); }
	/// \brief factory function to create a c2_diff from an regular algebraic expression.
	/// \note
	/// be very wary of ownership issues if this is used in a complex expression.
	c2_diff<float_type> &operator - (const c2_function<float_type> &rhs)  { return *new c2_diff<float_type>(*this, rhs); }
	/// \brief factory function to create a c2_product from an regular algebraic expression.
	/// \note
	/// be very wary of ownership issues if this is used in a complex expression.
	c2_product<float_type> &operator * (const c2_function<float_type> &rhs)  { return *new c2_product<float_type>(*this, rhs); }
	/// \brief factory function to create a c2_ratio from an regular algebraic expression.
	/// \note
	/// be very wary of ownership issues if this is used in a complex expression.
	c2_ratio<float_type> &operator / (const c2_function<float_type> &rhs)  { return *new c2_ratio<float_type>(*this, rhs); }
	

	
	std::vector<float_type> *sampling_grid;
	bool no_overwrite_grid;
	    
protected:
	c2_function(const c2_function<float_type>  &src) : sampling_grid(0),
		no_overwrite_grid(false),
        fXMin(src.fXMin), fXMax(src.fXMax), rootInitialized(false)
        {} // copy constructor only copies domain, and is only for internal use
	c2_function() : 
			sampling_grid(0), no_overwrite_grid(0), 
        fXMin(-std::numeric_limits<float_type>::max()), 
        fXMax(std::numeric_limits<float_type>::max()),
        rootInitialized(false)
		{} // prevent accidental naked construction (impossible any since this has pure virtual methods)
	
	// this should only be called very early on, by a constructor, before anyone else
	// sets a sampling grid, or it will leak memory
	virtual void set_sampling_grid_pointer(std::vector<float_type> &grid) 
		{ 
			if (sampling_grid && !no_overwrite_grid) delete sampling_grid; // grid was ours, lose it. 
			sampling_grid=&grid; no_overwrite_grid=1; 
		}
	
	float_type fXMin, fXMax;
	mutable int evaluations;
	
private:
	/// \brief structure used for recursion in adaptive integrator.  
	///
	/// Contains all the information for the function at one point. 
	struct c2_integrate_fblock {	float_type x, y, yp, ypp; };
	/// \brief structure used to pass information recursively.
	///
	/// the \a abs_tol is scaled by a factor of two at each division.  
	/// Everything else is just passed down.
	struct c2_integrate_recur { struct c2_integrate_fblock *f0, *f1, *f2;
		float_type abs_tol, rel_tol, *lr, eps_scale, extrap_coef, extrap2;
		int depth, derivs;
		bool adapt, extrapolate;
	};
	
    /// \brief Carry out the recursive subdivision and integration.
    ///
    /// This passes information recursively through the \a recur block pointer
    /// to allow very efficient recursion.
    /// \param rb a pointer to the recur struct.
	float_type integrate_step(struct c2_integrate_recur &rb) const;
    
    /// these carry a memory of the last root bracketing,
    /// to avoid the necessity of evaluating the function on the brackets every time
    /// if the brackets have not been changed. 
    mutable float_type lastRootLowerX, lastRootUpperX, lastRootLowerY, lastRootUpperY;
    mutable int rootInitialized;
	
};

/// \brief a container into which any other c2_function can be dropped, to allow expressions
/// with replacable components.  
///
///It is useful for plugging different InterpolatingFunctions into a c2_function expression.
///It saves a lot of effort in other places with casting away const declarations.
///
/// It is also useful as a wrapper for a function if it is necessary to have a copy of a function
/// which has a different domain or sampling grid than the parent function.  This can be
/// be used, for example, to patch badly-behaved functions with c2_piecewise_function by 
/// taking the parent function, creating two plugins of it with domains on each side of the 
/// nasty bit, and then inserting a nice function in the hole.
template <typename float_type=double> class c2_plugin_function : public c2_function<float_type> {
public:
	/// \brief construct the container with no function
	c2_plugin_function() : c2_function<float_type>(), func(0), owns(false)  {}
	/// \brief construct the container with a pre-defined function
	c2_plugin_function(const c2_function<float_type> &f) : c2_function<float_type>(), func(0), owns(false) { set_function(f); }
	/// \brief fill the container with a new function
	void set_function(const c2_function<float_type> &f) 
		{ 
			if(owns && func) delete func;
			func=&f; set_domain(f.xmin(), f.xmax()); 
			set_sampling_grid_pointer(*f.sampling_grid); 
			owns=false;
		}
	/// \copydoc c2_function::value_with_derivatives
	/// Uses the internal function pointer set by set_function().
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception) 
		{
			if(!func) throw c2_exception("c2_plugin_function<float_type> called uninitialized");
			return this->func->value_with_derivatives(x, yprime, yprime2); 
		}
	/// \brief clear the function
	///
	/// Any attempt to use this c2_plugin_function throws an exception if the saved function is cleared.
	void unset_function(void) { if(owns && func) delete func; func=0; owns=false; }
	/// \brief destructor
	~c2_plugin_function() { if(owns && func) delete func; }
	/// \brief tell us we should delete the function at destruction.  NOT sticky when function is reset
	void set_ownership() { this->owns=true; }
	
protected:
	const c2_function<float_type> *func;
	bool owns;
	
};

/// \brief Provides support for c2_function objects which are constructed from two other c2_function
/// objects.  
///
/// It provides a very primitive ownership concept, so that the creator can tag a function
/// as being owned by this function, so that when this function is deleted, the owned function will be deleted, too.
/// This allows a piece of code to create various c2_function objects, combine them with binary operators,
/// appropriately mark wich ones have no other possible owners, and return the final function with
/// reasonable faith that everything will get cleaned up when the final function is deleted.  Note that
/// none of this marking is automatic, to keep this class very lightweight.
template <typename float_type=double> class c2_binary_function : public c2_function<float_type> {
public:
	
	
	///  \brief function to manage the binary operation, used by c2_binary_function::value_with_derivatives() 
    /// 
    /// Normally not used alone, but can be used to combine functions in other contexts. 
	/// See interpolating_function::binary_operator() for an example.
	/// \param left the function on the left of the binary operator or outside the composition
	/// \param right the function to the right of the operator or inside the composition
    /// \param[in] x the point at which to evaluate the function
    /// \param[out] yprime the first derivative (if pointer is non-null)
    /// \param[out] yprime2 the second derivative (if pointer is non-null)
    /// \return the value of the function
	virtual float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right, 
					   float_type x, float_type *yprime, float_type *yprime2) const =0;
	
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw (c2_exception) 
    { 
		return this->combine(this->Left, this->Right, x, yprime, yprime2);
    }

	/// \brief allow c2_binary_function to remember ownership of contained functions for automatic cleanup
	///
	/// upon destruction, this will cause disposal of the left member function
	void set_left_ownership(void) { leftown=true; }
	/// \brief allow c2_binary_function to remember ownership of contained functions for automatic cleanup
	///
	/// upon destruction, this will cause disposal of the right member function
	void set_right_ownership(void) { rightown=true; }
	/// \brief allow c2_binary_function to remember ownership of contained functions for automatic cleanup
	///
	/// upon destruction, this will cause disposal of both member functions
	void set_ownership(void) { leftown=rightown=true; }
	/// \brief destructor executes disposal of member functions if flagged
	///
	/// depends on judicious use of set_ownership(), set_right_ownership(), or set_left_ownership()
	virtual ~c2_binary_function() {
		if(leftown) delete &Left;
		if(rightown) delete &Right;
	}
		
protected:
	/// \brief construct the binary function
	/// \param left the c2_function to be used in the left side of the binary relation
	/// \param right the c2_function to be used in the right side of the binary relation
	c2_binary_function( const c2_function<float_type> &left,  const c2_function<float_type> &right) : 
	c2_function<float_type>(), 
	Left(left), Right(right), leftown(false), rightown(false) 
	{ 
			set_domain(
					   (left.xmin() > right.xmin()) ? left.xmin() : right.xmin(), 
					   (left.xmax() < right.xmax()) ? left.xmax() : right.xmax()
					   );
	} 
	
	/// \brief construct a 'stub' c2_binary_function, which provides access to the combine() function
	/// \note Do not evaluate a 'stub' ever.  It is only used so that combine() can be called
    c2_binary_function() : c2_function<float_type>(), 
		Left(*((c2_function<float_type> *)0)), Right(*((c2_function<float_type> *)0)) {}
	
	const c2_function<float_type> &Left,  &Right;
	bool leftown, rightown;
		
};

/// \brief Create a very lightweight method to return a scalar multiple of another function.
/// 
template <typename float_type=double> class c2_scaled_function : public c2_plugin_function<float_type> {
public:
	/// \brief construct the function with its scale factor.
	///
	/// \param outer the function to be scaled
	/// \param scale the multiplicative scale factor
	c2_scaled_function(const c2_function<float_type> &outer, float_type scale) : 
		c2_plugin_function<float_type>(outer), yscale(scale) { }
	
	/// \brief set a new scale factor
	/// \param scale the new factor
	void reset(float_type scale) { yscale=scale; }
	
	/// \copydoc c2_function::value_with_derivatives
	/// 
	/// provide our own value_with_derivatives which bypasses the combiner for quicker operation
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw (c2_exception) 
    { 
		float_type y=this->func->value_with_derivatives(x, yprime, yprime2); 
		if(yprime) (*yprime)*=yscale; 
		if(yprime2) (*yprime2)*=yscale; 
		return y*yscale; 
    }

protected:
    c2_scaled_function<float_type>() {} // hide default constructor, since its use is almost always an error.
	float_type yscale;
};

/// \brief A container into which any other c2_function can be dropped.
///
/// It allows a function to be pre-evaluated at a point, and used at multiple places in an expression
/// efficiently. If it is re-evaluated at the previous point, it returns the remembered values;
/// otherwise, it re-evauates the function at the new point.
///
template <typename float_type=double> class c2_cached_function : public c2_plugin_function<float_type> {
public:
	/// \brief construct the container
	///
	/// \param f the function to be cached
	c2_cached_function(const c2_function<float_type> &f) : c2_plugin_function<float_type>(f), init(false)  {}
	/// \copydoc c2_function::value_with_derivatives
	///
	/// Checks to see if the function is being re-evaluated at the previous point, and
	/// returns remembered values if so.
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception) 
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
    c2_cached_function() {} // hide default constructor, since its use is almost always an error.
	mutable bool init;
	mutable float_type x0, y, yp, ypp;
	
};

/// \brief Provides function composition (nesting)
///
/// This allows evaluation of \a f(g(x)) where \a f and \a g are c2_function objects.
///
/// \note See c2_binary_function for discussion of ownership.
template <typename float_type=double> class  c2_composed_function : public c2_binary_function<float_type> {
public:
	
	/// \brief construct \a outer( \a inner (x))
    /// \note See c2_binary_function for discussion of ownership.
	/// \param outer the outer function 
	/// \param inner the inner function
	c2_composed_function(const c2_function<float_type> &outer, const c2_function<float_type> &inner) : c2_binary_function<float_type>(outer, inner) {}
	/// \brief Create a stub just for the combiner to avoid statics. 
	c2_composed_function() : c2_binary_function<float_type>() {} ;

	virtual float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right, 
						  float_type x, float_type *yprime, float_type *yprime2) const
    {
		float_type y0, yp0, ypp0, y1, yp1, ypp1;
		float_type *yp0p, *ypp0p, *yp1p, *ypp1p;
		if(yprime || yprime2) {
			yp0p=&yp0; ypp0p=&ypp0; yp1p=&yp1; ypp1p=&ypp1;
		} else {
			yp0p=ypp0p=yp1p=ypp1p=0;
		}
		
		y0=right.value_with_derivatives(x, yp0p, ypp0p);
		y1=left.value_with_derivatives(y0, yp1p, ypp1p);
		if(yprime) *yprime=yp1*yp0;
		if(yprime2) *yprime2=ypp0*yp1+yp0*yp0*ypp1;
		return y1;
	}	
};

/// \brief create a c2_function which is the sum of two other c2_function objects.
///
/// \note See c2_binary_function for discussion of ownership.
template <typename float_type=double> class c2_sum : public c2_binary_function<float_type> {
public:	
	/// \brief construct \a left + \a right
    /// \note See c2_binary_function for discussion of ownership.
	/// \param left the left function 
	/// \param right the right function
	c2_sum(const c2_function<float_type> &left, const c2_function<float_type> &right) : c2_binary_function<float_type>(left, right) {}
	/// \brief Create a stub just for the combiner to avoid statics. 
	c2_sum() : c2_binary_function<float_type>() {} ; // create a stub just for the combiner to avoid statics

	// function to do derivative arithmetic for sums
	virtual float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right, 
						  float_type x, float_type *yprime, float_type *yprime2) const
    {
		float_type y0, yp0, ypp0, y1, yp1, ypp1;
		float_type *yp0p, *ypp0p, *yp1p, *ypp1p;
		if(yprime || yprime2) {
			yp0p=&yp0; ypp0p=& ypp0; yp1p=&yp1; ypp1p=&ypp1;
		} else {
			yp0p=ypp0p=yp1p=ypp1p=0;
		}
		y0=left.value_with_derivatives(x, yp0p, ypp0p);
		y1=right.value_with_derivatives(x, yp1p, ypp1p);
		if(yprime) *yprime=yp0+yp1;
		if(yprime2) *yprime2=ypp0+ypp1;
		return y0+y1;
	}
};


/// \brief create a c2_function which is the difference of two other c2_functions.
///
/// \note See c2_binary_function for discussion of ownership.
template <typename float_type=double> class c2_diff : public c2_binary_function<float_type> {
public:	
	/// \brief construct \a left - \a right
    /// \note See c2_binary_function for discussion of ownership.
	/// \param left the left function 
	/// \param right the right function
	c2_diff(const c2_function<float_type> &left, const c2_function<float_type> &right) : c2_binary_function<float_type>(left, right) {}
	/// \brief Create a stub just for the combiner to avoid statics. 
	c2_diff() : c2_binary_function<float_type>() {} ; // create a stub just for the combiner to avoid statics

	// function to do derivative arithmetic for diffs
	virtual float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right, 
							  float_type x, float_type *yprime, float_type *yprime2) const
    {
		float_type y0, yp0, ypp0, y1, yp1, ypp1;
		float_type *yp0p, *ypp0p, *yp1p, *ypp1p;
		if(yprime || yprime2) {
			yp0p=&yp0; ypp0p=&ypp0; yp1p=&yp1; ypp1p=&ypp1;
		} else {
			yp0p=ypp0p=yp1p=ypp1p=0;
		}
		y0=left.value_with_derivatives(x, yp0p, ypp0p);
		y1=right.value_with_derivatives(x, yp1p, ypp1p);
		if(yprime) *yprime=yp0-yp1;
		if(yprime2) *yprime2=ypp0-ypp1;
		return y0-y1;
	}
};


/// \brief create a c2_function which is the product of two other c2_functions.
///
/// \note See c2_binary_function for discussion of ownership.
template <typename float_type=double> class c2_product : public c2_binary_function<float_type> {
public:	
	/// \brief construct \a left * \a right
    /// \note See c2_binary_function for discussion of ownership.
	/// \param left the left function 
	/// \param right the right function
	c2_product(const c2_function<float_type> &left, const c2_function<float_type> &right) : c2_binary_function<float_type>(left, right) {}
	/// \brief Create a stub just for the combiner to avoid statics. 
	c2_product() : c2_binary_function<float_type>() {} ; // create a stub just for the combiner to avoid statics

	virtual float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right, 
							   float_type x, float_type *yprime, float_type *yprime2) const
    {
		float_type y0, yp0, ypp0, y1, yp1, ypp1;
		float_type *yp0p, *ypp0p, *yp1p, *ypp1p;
		if(yprime || yprime2) {
			yp0p=&yp0; ypp0p=&ypp0; yp1p=&yp1; ypp1p=&ypp1;
		} else {
			yp0p=ypp0p=yp1p=ypp1p=0;
		}
		y0=left.value_with_derivatives(x, yp0p, ypp0p);
		y1=right.value_with_derivatives(x, yp1p, ypp1p);
		if(yprime) *yprime=y1*yp0+y0*yp1;
		if(yprime2) *yprime2=ypp0*y1+2.0*yp0*yp1+ypp1*y0;
		return y0*y1;
	}
};


/// \brief create a c2_function which is the ratio of two other c2_functions.
///
/// \note See c2_binary_function for discussion of ownership.
template <typename float_type=double> class c2_ratio : public c2_binary_function<float_type> {
public:	
	/// \brief construct \a left / \a right
    /// \note See c2_binary_function for discussion of ownership.
	/// \param left the left function 
	/// \param right the right function
	c2_ratio(const c2_function<float_type> &left, const c2_function<float_type> &right) : c2_binary_function<float_type>(left, right) {}
	/// \brief Create a stub just for the combiner to avoid statics. 
	c2_ratio() : c2_binary_function<float_type>() {} ; // create a stub just for the combiner to avoid statics
	
	virtual float_type combine(const c2_function<float_type> &left, const c2_function<float_type> &right, 
							  float_type x, float_type *yprime, float_type *yprime2) const
    {
		float_type y0, yp0, ypp0, y1, yp1, ypp1;
		float_type *yp0p, *ypp0p, *yp1p, *ypp1p;
		if(yprime || yprime2) {
			yp0p=&yp0; ypp0p=&ypp0; yp1p=&yp1; ypp1p=&ypp1;
		} else {
			yp0p=ypp0p=yp1p=ypp1p=0;
		}
		y0=left.value_with_derivatives(x, yp0p, ypp0p);
		y1=right.value_with_derivatives(x, yp1p, ypp1p);
		if(yprime) *yprime=(yp0*y1-y0*yp1)/(y1*y1); // first deriv of ratio
		if(yprime2) *yprime2=(y1*y1*ypp0+y0*(2*yp1*yp1-y1*ypp1)-2*y1*yp0*yp1)/(y1*y1*y1); // second deriv of ratio 
		return y0/y1;
	}

};

/// \brief a c2_function which is constant : can do interpolating_function  f2=f1 + c2_constant(11.3)
template <typename float_type> class c2_constant : public c2_function<float_type> {
public:
	c2_constant(float_type x=0.0) : c2_function<float_type>(), value(x) {}
	void reset(float_type val) { value=val; }
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception) 
	{ if(yprime) *yprime=0; if(yprime2) *yprime2=0; return value; }
	
private:
		float_type value;
};

/**
    \brief create a cubic spline interpolation of a set of (x,y) pairs

    This is one of the main reasons for c2_function objects to exist.

    It provides support for cubic spline interpolation of data provides from tables of \a x, \a y pairs.
    It supports automatic, transparent linearization of the data before storing in its tables (through
    subclasses such as 
    log_lin_interpolating_function, lin_log_interpolating_function, and
    log_log_interpolating_function) to permit very high accuracy representations of data which have a suitable
    structure.  It provides utility functions LinearInterpolatingGrid() and LogLogInterpolatingGrid()
    to create grids for mapping other functions onto a arithmetic or geometric grid.

    In its simplest form, an untransformed cubic spline of a data set, using natural boundary conditions
    (vanishing second derivative), is created as: \n
    \code 
    std::vector<double> xvals(10), yvals(10); 
    // < fill in xvals and yvals >
    interpolating_function<double>  myfunc(xvals, yvals);
    // and it can be evaluated at a point for its value only by: 
    double y=myfunc(x); 
    // or it can be evaluated with its derivatives by
    double yprime, yprime2; 
    double y=myfunc(x,&yprime, &yprime2);
    \endcode
*/

template <typename float_type=double> class interpolating_function  : public c2_function<float_type> {
public:
    /// \brief create the most general interpolating_function which defaults to linear-linear space
    ///
    /// lots to say here, but see Numerical Recipes for a discussion of cubic splines.
    /// \param x the list of abscissas.  Must be either strictly increasing or strictly decreasing.
	/// Strictly increasing is preferred, as less memory is used since a copy is not required for the sampling grid.
    /// \param f the list of function values.
    /// \param lowerSlopeNatural if true, set y''(first point)=0, otherwise compute it from \a lowerSope
    /// \param lowerSlope derivative of the function at the lower bound, used only if \a lowerSlopeNatural is false
    /// \param upperSlopeNatural if true, set y''(last point)=0, otherwise compute it from \a upperSope
    /// \param upperSlope derivative of the function at the upper bound, used only if \a upperSlopeNatural is false
    /// \param inputXConversion a function (not a c2_function) which converts \a x into the internal space. \n
	/// If this is NULL, use linear space and ignore inputXConversionPrime, inputXConversionDPrime 
    /// \param inputYConversion a function (not a c2_function) which converts \a y into the internal space. \n
	/// If this is NULL, use linear space and ignore inputYConversionPrime, inputYConversionDPrime, outputYConversion 
    /// \param outputYConversion a function (not a c2_function) which converts \a y out of the internal space
    /// \param inputXConversionPrime the derivative of \a inputXConversion
    /// \param inputYConversionPrime the derivative of \a inputYConversion
    /// \param inputXConversionDPrime the second derivative of \a inputXConversion
    /// \param inputYConversionDPrime the second derivative of \a inputYConversion
    interpolating_function(const std::vector<float_type> &x, const std::vector<float_type> &f, 
						  bool lowerSlopeNatural=true, float_type lowerSlope=0.0, 
						  bool upperSlopeNatural=true, float_type upperSlope=0.0,
						  float_type (*inputXConversion)(float_type)=0, 
						  float_type (*inputYConversion)(float_type)=0, 
						  float_type (*outputYConversion)(float_type)=0, 
						  float_type (*inputXConversionPrime)(float_type)=0, 
						  float_type (*inputYConversionPrime)(float_type)=0, 
						  float_type (*inputXConversionDPrime)(float_type)=0, 
						  float_type (*inputYConversionDPrime)(float_type)=0 
						   ) throw(c2_exception) : c2_function<float_type>()
        { init(x, f, lowerSlopeNatural, lowerSlope, upperSlopeNatural, upperSlope, 
               inputXConversion, inputYConversion, outputYConversion,
               inputXConversionPrime, inputYConversionPrime, 
               inputXConversionDPrime, inputYConversionDPrime 
               );
        }
	
    /// \brief copy constructor
    /// \param rhs interpolating_function  to copy from
	interpolating_function(const interpolating_function <float_type> &rhs) : c2_function<float_type>(rhs), 
		Xraw(rhs.Xraw), X(rhs.X), F(rhs.F), y2(rhs.y2),
		fXin(rhs.fXin), fYin(rhs.fYin), fYout(rhs.fYout), 
		fXinPrime(rhs.fXinPrime), fYinPrime(rhs.fYinPrime), 
		fXinDPrime(rhs.fXinDPrime), fYinDPrime(rhs.fYinDPrime) ,
		xInverted(rhs.xInverted), lastKLow(-1)
	{	set_sampling_grid_pointer(Xraw); }
	
	virtual ~interpolating_function() { } // just to suppress warnings about no virtual destructor
	
    /// \brief retrieve copies of the x & y tables from which this was built
    ///
	/// This is often useful in the creation of new interpolating functions with transformed data.
	/// The vectorswill have their sizes set correctly on return.
	/// \param [in, out] xvals the abscissas 
	/// \param [in, out] yvals the ordinates
	void get_data(std::vector<float_type> &xvals, std::vector<float_type> &yvals) const throw() ;
	
    /// \brief enable extrapolation of the function below the tabulated data.
    ///
    /// This allows the interpolator to be extrapolated outside the bounds of the data,
    /// using whatever derivatives it already had at the lower bound.
    /// \param bound the abscissa to which the function should be extended.
	void set_lower_extrapolation(float_type bound);
    /// \brief enable extrapolation of the function above the tabulated data.
    ///
    /// This allows the interpolator to be extrapolated outside the bounds of the data,
    /// using whatever derivatives it already had at the upper bound.
    /// \param bound the abscissa to which the function should be extended.
	void set_upper_extrapolation(float_type bound);
	
	// these functions correctly combine the interpolating function with another interpolating function 
	// preserving the X bounds and mapping functions of the host (left hand) function.
	
    /// \brief create a new interpolating_function  which is the \a source 
    /// function applied to every point in the interpolating tables
    ///
    /// This carefully manages the derivative of the composed function at the two ends.
    /// \param source the function to apply
    /// \return a new interpolating_function  with the same mappings for x and y
	interpolating_function <float_type> & unary_operator(const c2_function<float_type> &source) const;

    /// \brief create a new interpolating_function  which is the parent interpolating_function  
    /// combined with \a rhs using \a combiner at every point in the interpolating tables
    ///
    /// This carefully manages the derivative of the composed function at the two ends.
    /// \param rhs the function to apply
    /// \param combining_stub a function which defines which binary operation to use.
    /// \return a new interpolating_function  with the same mappings for x and y
	interpolating_function <float_type> & binary_operator(const c2_function<float_type> &rhs,
           c2_binary_function<float_type> *combining_stub
           ) const;
	
	// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	// can be upcast back to a c2_function to produce unprocessed binaries.

	/// \brief produce a newly resampled interpolating_function  which is the specified sum.
    /// \param rhs the function to add, pointwise
    /// \return a new interpolating_function 
    /// \note
	/// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	/// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	/// can be upcast back to a c2_function to produce unprocessed binaries.
    interpolating_function <float_type> & operator + (const c2_function<float_type> &rhs) const { 
		return binary_operator(rhs, new c2_sum<float_type>()); }
	/// \brief produce a newly resampled interpolating_function  which is the specified difference.
    /// \param rhs the function to subtract, pointwise
    /// \return a new interpolating_function 
    /// \note
	/// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	/// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	/// can be upcast back to a c2_function to produce unprocessed binaries.
	interpolating_function <float_type> & operator - (const c2_function<float_type> &rhs) const {
		return binary_operator(rhs, new c2_diff<float_type>()); }
	/// \brief produce a newly resampled interpolating_function  which is the specified product.
    /// \param rhs the function to multiply, pointwise
    /// \return a new interpolating_function 
    /// \note
	/// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	/// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	/// can be upcast back to a c2_function to produce unprocessed binaries.
	interpolating_function <float_type> & operator * (const c2_function<float_type> &rhs) const { 
		return binary_operator(rhs, new c2_product<float_type>()); }
	/// \brief produce a newly resampled interpolating_function  which is the specified ratio.
    /// \param rhs the function to divide, pointwise
    /// \return a new interpolating_function 
    /// \note
	/// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	/// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	/// can be upcast back to a c2_function to produce unprocessed binaries.
	interpolating_function <float_type> & operator / (const c2_function<float_type> &rhs) const { 
		return binary_operator(rhs, new c2_ratio<float_type>()); }
	/// \brief produce a newly resampled interpolating_function  which is the specified sum.
    /// \param rhs a constant to add, pointwise
    /// \return a new interpolating_function 
    /// \note
	/// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	/// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	/// can be upcast back to a c2_function to produce unprocessed binaries.
	interpolating_function <float_type> & operator + (float_type rhs) const { return (*this)+c2_constant<float_type>(rhs); }
	/// \brief produce a newly resampled interpolating_function  which is the specified difference.
    /// \param rhs a constant to subtract, pointwise
    /// \return a new interpolating_function 
    /// \note
	/// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	/// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	/// can be upcast back to a c2_function to produce unprocessed binaries.
	interpolating_function <float_type> & operator - (float_type rhs) const { return (*this)-c2_constant<float_type>(rhs); }
	/// \brief produce a newly resampled interpolating_function  which is the specified product.
    /// \param rhs a constant to multiply, pointwise
    /// \return a new interpolating_function 
    /// \note
	/// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	/// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	/// can be upcast back to a c2_function to produce unprocessed binaries.
	interpolating_function <float_type> & operator * (float_type rhs) const { return (*this)*c2_constant<float_type>(rhs); }
	/// \brief produce a newly resampled interpolating_function  which is the specified ratio.
    /// \param rhs a constant to divide, pointwise
    /// \return a new interpolating_function 
    /// \note
	/// InterpolatingFunctions override the c2_function operators, since they explicitly re-generate the interpolation table
	/// when they are applied.  If this is not desired, these operators are not virtual, so the interpolating_function 
	/// can be upcast back to a c2_function to produce unprocessed binaries.
	interpolating_function <float_type> & operator / (float_type rhs) const { return (*this)/c2_constant<float_type>(rhs); }
        
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception);
	
	/// \brief move value & derivatives into our internal coordinates (use splint to go the other way!)
    /// \note why?
	void localize_derivatives(float_type xraw, float_type y, float_type yprime, float_type yprime2, float_type *y0, float_type *yp0, float_type *ypp0) const;
	
protected:
		
	interpolating_function() : c2_function<float_type>() { } // default constructor is never used, prevent accidents by protecting it.
    
    /// \brief do the dirty work of constructing the spline.  See interpolating_function  constructor for details.
	void init(const std::vector<float_type> &, const std::vector<float_type> &, 
			  bool lowerSlopeNatural, float_type lowerSlope, 
			  bool upperSlopeNatural, float_type upperSlope,
			  float_type (*inputXConversion)(float_type)=0, 
			  float_type (*inputXConversionPrime)(float_type)=0, 
			  float_type (*inputXConversionDPrime)(float_type)=0, 
			  float_type (*inputYConversion)(float_type)=0, 
			  float_type (*inputYConversionPrime)(float_type)=0, 
			  float_type (*inputYConversionDPrime)(float_type)=0, 
			  float_type (*outputYConversion)(float_type)=0
			  ) throw(c2_exception) ;

    std::vector<float_type> Xraw, X, F, y2;
	
	float_type (*fXin)(float_type), (*fYin)(float_type), (*fYout)(float_type);
	float_type (*fXinPrime)(float_type), (*fYinPrime)(float_type);
	float_type (*fXinDPrime)(float_type), (*fYinDPrime)(float_type);
	
	int xInverted;
    mutable int lastKLow;
};

/// \brief An interpolatingFunction with X transformed into log space.  
///
/// Most useful for functions looking like y=log(x) or any other function with a huge X dynamic range,
/// and a slowly varying Y.
template <typename float_type=double> class log_lin_interpolating_function : public interpolating_function <float_type> {
public:
    /// \brief Construct the function.
    /// \param x the list of abscissas.  Must be either strictly increasing or strictly decreasing.
	/// Strictly increasing is preferred, as less memory is used since a copy is not required for the sampling grid.
    /// \param f the list of function values.
    /// \param lowerSlopeNatural if true, set y''(first point)=0 in LogLin space, otherwise compute it from \a lowerSope
    /// \param lowerSlope derivative of the function at the lower bound, used only if \a lowerSlopeNatural is false
    /// \param upperSlopeNatural if true, set y''(last point)=0 in LogLin space, otherwise compute it from \a upperSope
    /// \param upperSlope derivative of the function at the upper bound, used only if \a upperSlopeNatural is false
    log_lin_interpolating_function(const std::vector<float_type> &x, const std::vector<float_type> &f, 
								bool lowerSlopeNatural=true, float_type lowerSlope=0.0, 
								bool upperSlopeNatural=true, float_type upperSlope=0.0);
protected:
	log_lin_interpolating_function() {} // do not allow naked construction... it is usually an accident.
};


/// \brief An interpolatingFunction with Y transformed into log space.  
///
/// Most useful for functions looking like y=exp(x)
template <typename float_type=double> class lin_log_interpolating_function : public interpolating_function <float_type> {
public:
    /// \brief Construct the function.
    /// \param x the list of abscissas.  Must be either strictly increasing or strictly decreasing.
	/// Strictly increasing is preferred, as less memory is used since a copy is not required for the sampling grid.
    /// \param f the list of function values.
    /// \param lowerSlopeNatural if true, set y''(first point)=0 in LinLog space, otherwise compute it from \a lowerSope
    /// \param lowerSlope derivative of the function at the lower bound, used only if \a lowerSlopeNatural is false
    /// \param upperSlopeNatural if true, set y''(last point)=0 in LinLog space, otherwise compute it from \a upperSope
    /// \param upperSlope derivative of the function at the upper bound, used only if \a upperSlopeNatural is false
    lin_log_interpolating_function(const std::vector<float_type> &x, const std::vector<float_type> &f, 
								bool lowerSlopeNatural=true, float_type lowerSlope=0.0, 
								bool upperSlopeNatural=true, float_type upperSlope=0.0);
protected:
	lin_log_interpolating_function() {} // do not allow naked construction... it is usually an accident.
};


/// \brief An interpolatingFunction with X and Y transformed into log space.  
///
/// Most useful for functions looking like y=x^n or any other function with a huge X and Y dynamic range.
template <typename float_type=double> class log_log_interpolating_function : public interpolating_function <float_type> {
public:
    /// \brief Construct the function.
    /// \param x the list of abscissas.  Must be either strictly increasing or strictly decreasing.
	/// Strictly increasing is preferred, as less memory is used since a copy is not required for the sampling grid.
    /// \param f the list of function values.
    /// \param lowerSlopeNatural if true, set y''(first point)=0 in LogLog space, otherwise compute it from \a lowerSope
    /// \param lowerSlope derivative of the function at the lower bound, used only if \a lowerSlopeNatural is false
    /// \param upperSlopeNatural if true, set y''(last point)=0 in LogLog space, otherwise compute it from \a upperSope
    /// \param upperSlope derivative of the function at the upper bound, used only if \a upperSlopeNatural is false
    log_log_interpolating_function(const std::vector<float_type> &x, const std::vector<float_type> &f, 
								bool lowerSlopeNatural=true, float_type lowerSlope=0.0, 
								bool upperSlopeNatural=true, float_type upperSlope=0.0);
protected:
	log_log_interpolating_function() {} // do not allow naked construction... it is usually an accident.
};


/// \brief An interpolating_function  with X in reciprocal space and Y transformed in log space.  
///
/// Most useful for thermodynamic types of data where Y is roughly A*exp(-B/x). 
/// Typical examples are reaction rate data, and thermistor calibration data.
template <typename float_type=double> class arrhenius_interpolating_function : public interpolating_function <float_type> {
public:
    /// \brief Construct the function.
    /// \param x the list of abscissas.  Must be either strictly increasing or strictly decreasing.
	/// Strictly increasing is preferred, as less memory is used since a copy is not required for the sampling grid.
    /// \param f the list of function values.
    /// \param lowerSlopeNatural if true, set y''(first point)=0 in Arrhenius space, otherwise compute it from \a lowerSope
    /// \param lowerSlope derivative of the function at the lower bound, used only if \a lowerSlopeNatural is false
    /// \param upperSlopeNatural if true, set y''(last point)=0 in Arrhenius space, otherwise compute it from \a upperSope
    /// \param upperSlope derivative of the function at the upper bound, used only if \a upperSlopeNatural is false
    arrhenius_interpolating_function(const std::vector<float_type> &x, const std::vector<float_type> &f, 
								bool lowerSlopeNatural=true, float_type lowerSlope=0.0, 
								bool upperSlopeNatural=true, float_type upperSlope=0.0);
protected:
	arrhenius_interpolating_function() {} // do not allow naked construction... it is usually an accident.
};

/**
 \brief create a linear-linear interpolating grid with both x & y set to 
 (xmin, xmin+dx, ... xmin + (count-1)*dx )

 very useful for transformaiton with other functions e.g. 
 \code
 f=c2_sin<double>::sin(LinearInterpolatingGrid(-0.1,0.1, 65)) 
 \endcode
 creates a spline table of sin(x) slightly beyond the first period
 \param xmin the starting point for the grid
 \param dx the step size for the grid
 \param count the number of points in the  grid
 \return an identity interpolating_function  with the requested grid 
 */
template <typename float_type> interpolating_function <float_type> &linear_interpolating_grid(float_type xmin, float_type dx, int count) {
	std::vector<float_type> x(count);
	for(int i=0; i<count; i++) x[i]=xmin + i * dx;
	return *new interpolating_function <float_type>(x,x);
}

/**
 \brief create a log-log interpolating grid with both x & y set to 
 (xmin, xmin*dx, ... xmin * dx^(count-1) )

 very useful for transformaiton with other functions e.g. 
 \code
 f=c2_log<double>::log(LogLogInterpolatingGrid(2, 1.1, 65)) 
 \endcode
 creates a spline table of log(x)
 \param xmin the starting point for the grid
 \param dx the ratio between points
 \param count the number of points in the  grid
 \return an identity log_log_interpolating_function with the requested grid 
*/
template <typename float_type> log_log_interpolating_function <float_type> &log_log_interpolating_grid(float_type xmin, float_type dx, int count) {
	std::vector<float_type> x(count);
	x[0]=xmin;
	for(int i=1; i<count; i++) x[i]=dx*x[i-1];
	return *new log_log_interpolating_function<float_type>(x,x);
}

/// \brief compute sin(x) with its derivatives.
///
/// Creates a singleton instance c2_sin::sin of itself for convenient access.
template <typename float_type=double> class c2_sin : public c2_function<float_type> {
public:
	/// \brief constructor.  There is alread a singleton c2_sin::sin, which usually suffices.
	c2_sin() {}
	
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ float_type q=std::sin(x); if(yprime) *yprime=std::cos(x); if(yprime2) *yprime2=-q; return q; }	
    
    /// \brief return a grid dynamically, suitable for use with trig functions with period 2*pi
    /// \param xmin the lower bound for the grid
    /// \param xmax upper bound for the grid
    /// \return a new sampling grid.
	virtual std::vector<float_type> &get_sampling_grid(float_type xmin, float_type xmax); 
    /// \brief the static singleton 
	static const c2_sin sin;
};
/// \brief compute cos(x) with its derivatives.
///
/// Creates a singleton instance c2_cos::cos of itself for convenient access.
template <typename float_type=double> class c2_cos : public c2_sin<float_type> {
public:
	/// \brief constructor.  There is already a singleton c2_cos::cos, which usually suffices.
	c2_cos() {}
	
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ float_type q=std::cos(x); if(yprime) *yprime=-std::sin(x); if(yprime2) *yprime2=-q; return q; }	
    /// \brief the static singleton 
	static const c2_cos cos;
};
/// \brief compute tan(x) with its derivatives.
///
/// Creates a singleton instance c2_tan::tan of itself for convenient access.
template <typename float_type=double> class c2_tan : public c2_function<float_type> {
public:
	/// \brief constructor.  There is already a singleton c2_tan::tan, which usually suffices.
	c2_tan() {}
	
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{
		float_type c=std::cos(x), s=std::sin(x);
		float_type t=s/c;
		float_type yp=1/(c*c);
		if(yprime) *yprime=yp; if(yprime2) *yprime2=2*t*yp; 
		return t; 
	}	
    /// \brief the static singleton 
	static const c2_tan tan;
};
/// \brief compute log(x) with its derivatives.
///
/// Creates a singleton instance c2_log::log of itself for convenient access.
template <typename float_type=double> class c2_log : public c2_function<float_type> {
public:
	/// \brief constructor.  There is already a singleton c2_log::log, which usually suffices.
	c2_log() {}

	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ if(yprime) *yprime=1.0/x; if(yprime2) *yprime2=-1.0/(x*x); return std::log(x); }	
    /// \brief the static singleton 
	static const c2_log log;
};
/// \brief compute exp(x) with its derivatives.
///
/// Creates a singleton instance c2_exp::exp of itself for convenient access.
template <typename float_type=double>  class c2_exp : public c2_function<float_type> {
public:
	/// \brief constructor.  There is already a singleton c2_exp::exp, which usually suffices.
	c2_exp() {}

	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ float_type q=std::exp(x); if(yprime) *yprime=q; if(yprime2) *yprime2=q; return q; }
    /// \ brief the static singleton 
	static const c2_exp exp;
};
/// \brief compute sqrt(x) with its derivatives.
///
/// Creates a singleton instance c2_sqrt::sqrt of itself for convenient access.
template <typename float_type=double> class c2_sqrt : public c2_function<float_type> {
public:
	/// \brief constructor.  There is already a singleton c2_sqrt::sqrt, which usually suffices.
	c2_sqrt() {}
	
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ float_type q=std::sqrt(x); if(yprime) *yprime=0.5/q; if(yprime2) *yprime2=-0.25/(x*q); return q; }
    /// \brief the static singleton 
	static const c2_sqrt sqrt;
};
/// \brief compute scale/x with its derivatives.
///
/// Creates a singleton instance c2_recip:recip of itself for convenient access.
template <typename float_type=double> class c2_recip : public c2_function<float_type> {
public:
	/// \brief constructor.  There is already a singleton c2_recip::recip, which usually suffices.
	c2_recip(float_type scale) : rscale(scale) {}
	
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
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
    /// \brief the static singleton 
	static const c2_recip recip;
private:
	float_type rscale;
};
/// \brief compute x with its derivatives.
///
/// Creates a singleton instance c2_identity::identity of itself for convenient access.
template <typename float_type=double> class c2_identity : public c2_function<float_type> {
public:
	/// \brief constructor.  There is already a singleton c2_identity::identity, which usually suffices.
	c2_identity() {}
	
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ if(yprime) *yprime=1.0; if(yprime2) *yprime2=0; return x; }
    /// \brief the static singleton 
	static const c2_identity identity;
};

/**
 \brief create a linear mapping of another function 
 
 for example, given a c2_function \a f 
 \code 
 c2_linear<double> L(1.2, 2.0, 3.0);
 c2_composed_function<double> &F=L(f); 
 \endcode
 produces a new c2_function F=2.0+3.0*(\a f - 1.2) 
*/
template <typename float_type=double> class c2_linear : public c2_function<float_type> {
public:
    /// \brief Construct the operator f=y0 + slope * (x-x0)
	/// \param x0 the x offset
    /// \param y0 the y-intercept i.e. f(x0)
    /// \param slope the slope of the mapping
	c2_linear(float_type x0, float_type y0, float_type slope) : xint(x0), intercept(y0), m(slope) {}
    /// \brief Change the slope and intercepts after construction.
	/// \param x0 the x offset
    /// \param y0 the y-intercept 
    /// \param slope the slope of the mapping
	void reset(float_type x0, float_type y0, float_type slope) { xint=x0; intercept=y0; m=slope; } 
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ if(yprime) *yprime=m; if(yprime2) *yprime2=0; return m*(x-xint)+intercept; }
	
private:
		float_type xint, intercept, m;
protected:
		c2_linear() {} // do not allow naked construction... it is usually an accident.
};

/**
\brief create a quadratic mapping of another function 
 
 for example, given a c2_function \a f 
 \code 
 c2_quadratic<double> Q(1.2, 2.0, 3.0, 4.0);
 c2_composed_function<double> &F=Q(f); 
 \endcode
 produces a new c2_function F=2.0 + 3.0*(f-1.2) + 4.0*(f-1.2)^2 

 note that the parameters are overdetermined, but allows the flexibility of two different representations

 */
template <typename float_type=double> class c2_quadratic : public c2_function<float_type> {
public:
    /// \brief Construct the operator
    /// \param x0 the center around which the powers are computed
    /// \param y0 the value of the function at \a x = \a x0
    /// \param xcoef the scale on the (\a x - \a x0) term
    /// \param x2coef the scale on the (\a x - \a x0)^2 term
	c2_quadratic(float_type x0, float_type y0, float_type xcoef, float_type x2coef) : intercept(y0), center(x0), a(x2coef), b(xcoef) {}
    /// Modify the mapping after construction
    /// \param x0 the new center around which the powers are computed
    /// \param y0 the new value of the function at \a x = \a x0
    /// \param xcoef the new scale on the (\a x - \a x0) term
    /// \param x2coef the new scale on the (\a x - \a x0)^2 term    
	void reset(float_type x0, float_type y0, float_type xcoef, float_type x2coef) { intercept=y0; center=x0; a=x2coef; b=xcoef; } 
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ float_type dx=x-center; if(yprime) *yprime=2*a*dx+b; if(yprime2) *yprime2=2*a; return a*dx*dx+b*dx+intercept; }
	
private:
		float_type intercept, center, a, b;
protected:
		c2_quadratic() {} // do not allow naked construction... it is usually an accident.
};

/**
\brief create a power law mapping of another function 
 
 for example, given a c2_function \a f 
 \code 
 c2_power_law<double> PLaw(1.2, 2.5);
 c2_composed_function<double> &F=PLaw(f); 
 \endcode
 produces a new c2_function F=1.2 * f^2.5 
 
 */
template <typename float_type=double> class c2_power_law : public c2_function<float_type> {
public:
    /// \brief Construct the operator
    /// \param scale the multipler
    /// \param power the exponent
	c2_power_law(float_type scale, float_type power) : a(scale), b(power) {}
    /// \brief Modify the mapping after construction
    /// \param scale the new multipler
    /// \param power the new exponent
	void reset(float_type scale, float_type power) { a=scale; b=power; }
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception)
	{ float_type q=a*std::pow(x,b-2); if(yprime) *yprime=b*q*x; if(yprime2) *yprime2=b*(b-1)*q; return q*x*x; }
	
private:
		float_type a, b;
protected:
		c2_power_law() {} // do not allow naked construction... it is usually an accident.
};

/**
\brief create the formal inverse function of another function 
 
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
 
 Note that it is a subclass of c2_scaled_function only to manage ownership of another c2_function.
 */
template <typename float_type=double> class c2_inverse_function : public c2_plugin_function<float_type> {
public:
    /// \brief Construct the operator
    /// \param source the function to be inverted
	c2_inverse_function(const c2_function<float_type> &source);  
    virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw(c2_exception);
    
    /// \brief give the function a hint as to where to look for its inverse
    /// \param hint the likely value of the inverse, which defaults to whatever the evaluation returned.
    void set_start_hint(float_type hint) const { start_hint=hint; }
    
    /// \brief get the starting hint.  
    /// 
    /// This is virtual so if there is a better way, this can be easily overridden.
    ///  It is used in value_with_derivatives() to guess where to start the root finder.
    /// \param x the abscissa for which an estimate is needed
    virtual float_type get_start_hint(float_type x) const { return start_hint; } 
    
protected:
    c2_inverse_function() {} // do not allow naked construction... it is usually an accident.
    mutable float_type start_hint;
};

/** 
  \brief
    An interpolating_function  which is the cumulative integral of a histogram.

    Note than binedges should be one element longer than binheights, since the lower & upper edges are specified. 
    Note that this is a malformed spline, since the second derivatives are all zero, so it has less continuity.
    Also, note that the bin edges can be given in backwards order to generate the 
    reversed accumulation (starting at the high end) 
*/

template <typename float_type=double>  class accumulated_histogram : public interpolating_function <float_type> {
public:
    /// \brief Construct the integrated histogram
    /// \param binedges the edges of the bins in \a binheights.  It should have one more element than \a binheights
    /// \param binheights the number of counts in each bin.
    /// \param normalize if true, normalize integral to 1
    /// \param inverse_function if true, drop zero channels, and return inverse function for random generation
    /// \param drop_zeros eliminate null bins before integrating, so integral is strictly monotonic. 
	accumulated_histogram(const std::vector<float_type>binedges, const std::vector<float_type> binheights,
						 bool normalize=false, bool inverse_function=false, bool drop_zeros=true);
	
};

/**
  \brief Construct a function useful for generation of random numbers from the given distribution
  
  inverse_integrated_density<InterpolatingFunctionFlavor>() starts with a probability density c2_function, generates the integral, 
  and generates an interpolating_function  which, when evaluated using a uniform random on [0,1] returns values
  with a density distribution equal to the input distribution
  If the data are passed in reverse order (large X first), the integral is carried out from the big end.
 
  \sa  template <typename Intermediate, typename Final> Final inverse_integrated_density(const std::vector, c2_function &)
 
  \param bincenters points at which to sample the c2_function \a binheights
  \param binheights a c2_function which describes the random number distribution to be produced.
  \return an interpolating_function  of the type requested in the template which,
     if evaluated randomly with a uniform variate on [0,1) produces numbers
     distributed according to \a binheights
*/

template <typename float_type, typename Final > 
	Final &inverse_integrated_density(const std::vector<float_type> &bincenters, c2_function<float_type> &binheights)
{	
	std::vector<float_type> integral;
	
	// integrate from first to last bin in original order, leaving results in integral
	// ask for relative error of 1e-6 on each bin, with absolute error set to 0 (since we don't know the data scale).
	float_type sum=binheights.partial_integrals(bincenters, &integral, 0.0, 1e-6); 
	// the integral vector now has partial integrals... it must be accumulated by summing
	integral.insert(integral.begin(), 0.0); // integral from start to start is 0
	float_type scale=1.0/sum;
	for(size_t i=1; i<integral.size(); i++) integral[i]=integral[i]*scale + integral[i-1];
	integral.back()=1.0; // force exact value on boundary
	
	return  *new Final(integral, bincenters, 
					   false, 1.0/(scale*binheights(bincenters.front() )), 
					   false, 1.0/(scale*binheights(bincenters.back() ))
					   ); // use integral as x axis in inverse function
}

/**
 \brief Construct a function useful for generation of random numbers from the given distribution

 \code
 template <typename Intermediate, typename Final> 
     Final & inverse_integrated_density(const std::vector &bincenters, const std::vector &binheights) 
 \endcode
 is a variant of \code
 template <typename Final> 
     Final & inverse_integrated_density(const std::vector &bincenters, c2_function &binheights)
 \endcode
 which takes two std::vectors and generates the intermediate interpolating_function  required for 
 inverse_integrated_density(), and then calls it.

 \param bincenters points at which \a binheights are defined
 \param binheights an std::vector which describes the random number distribution to be produced.
 \return an interpolating_function  of the type requested in the template which,
 if evaluated randomly with a uniform variate on [0,1) produces numbers
  distributed according to \a binheights
*/

template <typename float_type, typename Intermediate, typename Final> Final
    &inverse_integrated_density(const std::vector<float_type> &bincenters, const std::vector<float_type> &binheights) 
{	
	std::vector<float_type> be(bincenters), bh(binheights);
	
	if(be[1] < be[0]) { // reverse data for interpolator if x axis passed in backwards
		std::reverse(be.begin(), be.end());
		std::reverse(bh.begin(), bh.end());
	}
	
	Intermediate temp(be, bh); // create a temporary interpolating_function  to integrate
	Final &result=inverse_integrated_density<Final>(bincenters, temp);
	
	return result;
}

/// \brief create a c2_function which smoothly connects two other c2_functions.
///
/// This takes two points and generates a polynomial which matches two c2_function arguments 
/// at those two points, with two derivatives at each point, and an arbitrary value at the center of the 
/// region.  It is useful for splicing together functions over rough spots (0/0, for example).
/// 
/// If \a auto_center is true, the value at the midpoint is computed so that the resulting polynomial is
/// of order 5.  If \a auto_center is false, the value \a y1 is used at the midpoint, resulting in a 
/// polynomial of order 6.
template <typename float_type=double> class c2_connector_function : public c2_function<float_type> {
public:	
	/// \brief construct the container
	/// \param f1 the function on the left side to be connected 
	/// \param f2 the function on the right side to be connected
	/// \param x0 the point at which to match \a f1 and its derivatives
	/// \param x2 the point at which to match \a f2 and its derivatives
	/// \param auto_center if true, no midpoint value is specified.  If false, match the value \a y1 at the midpoint
	/// \param y1 the value to match at the midpoint, if \a auto_center is false
	/// \return a c2_function with domain (\a x0,\a x2) which smoothly connects \a f1 and \a f2
	c2_connector_function(const c2_function<float_type> &f1, const c2_function<float_type> &f2, float_type x0, float_type x2, 
			bool auto_center, float_type y1);
	/// \brief destructor
	virtual ~c2_connector_function();
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw (c2_exception);
protected:
	float_type fx1, fhinv, fdx, fy1, fa, fb, fc, fd, fe, ff;
};



/// \brief create a c2_function which is a piecewise assembly of other c2_functions.
///
/// The functions must have increasing, non-overlapping domains.  Any empty space
/// between functions will be filled with a linear interpolation.
/// \note The creation of the container results in the creation of an explicit sampling grid.  
/// If this is used with functions with a large domain, or which generate very dense sampling grids,
/// it could eat a lot of memory.  Do not abuse this by using functions which can generate gigantic grids.
/// 
/// See c2_plugin_function for a discussion of how this might be used.
template <typename float_type=double> class c2_piecewise_function : public c2_function<float_type> {
public:	
	/// \brief construct the container
	c2_piecewise_function();
	/// \brief destructor
	virtual ~c2_piecewise_function();
	virtual float_type value_with_derivatives(float_type x, float_type *yprime, float_type *yprime2) const throw (c2_exception);
	/// \brief append a new function to the sequence
	///
	/// This takes a c2_function, and appends it onto the end of the piecewise collection.
	/// The domain of the function (which MUST be set) specifies the place it will be used in 
	/// the final function.  If the domain exactly abuts the domain of the previous function, it
	/// will be directly attached.  If there is a gap, the gap will be filled in by linear interpolation.
	/// If the function being appended is to be deleted automatically when this container is deleted, set the pass_ownership flag.
	/// \param func a c2_function with a defined domain to be appended to the collection
	/// \param pass_ownership if set, \a func will be deleted when the container is destroyed
	void append_function(c2_function<float_type> &func, bool pass_ownership) throw (c2_exception);
protected:
	std::vector<c2_function<float_type> *> functions;
	std::vector<bool> owns;
	mutable int lastKLow;
};

#include "c2_function.icc"

#endif 

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
#ifndef G4SHARED_PTR_HH_
#define G4SHARED_PTR_HH_

//#if __cplusplus >= 201103L
#include <memory>

#define G4shared_ptr std::shared_ptr
#define G4weak_ptr std::weak_ptr
#define G4static_pointer_cast std::static_pointer_cast
#define G4const_pointer_cast std::const_pointer_cast
#define G4dynamic_pointer_cast std::dynamic_pointer_cast
#define G4enable_shared_from_this std::enable_shared_from_this
#define G4enable_shared_from_this2 std::enable_shared_from_this
/*
#else
#include "CLHEP/Utility/memory.h"

#define G4shared_ptr G4::shared_ptr
#define G4weak_ptr G4::weak_ptr
#define G4static_pointer_cast G4::static_pointer_cast
#define G4const_pointer_cast G4::const_pointer_cast
#define G4dynamic_pointer_cast G4::dynamic_pointer_cast
#define G4enable_shared_from_this G4::enable_shared_from_this
#define G4enable_shared_from_this2 G4::enable_shared_from_this2

namespace G4
{
using CLHEP::shared_ptr;
using CLHEP::weak_ptr;
using CLHEP::static_pointer_cast;
using CLHEP::const_pointer_cast;
using CLHEP::dynamic_pointer_cast;
}

namespace CLHEP
{

template< typename T >
class enable_shared_from_this
{
protected:
	enable_shared_from_this()
{ }

	~enable_shared_from_this()
	{ }

	enable_shared_from_this( enable_shared_from_this const & )
	{ }

	enable_shared_from_this &
	operator = ( enable_shared_from_this const & )
	{
		return *this;
	}

public:
	shared_ptr<T>
	shared_from_this()
	{
		shared_ptr<T>  p( weak_this_ );
		// assert( p.get() == this );
		return p;
	}

	shared_ptr<T const>
	shared_from_this() const
	{
		shared_ptr<T const>  p( weak_this_ );
		// assert( p.get() == this );
		return p;
	}

public: // actually private, but avoids compiler template friendship issues

	// Note: invoked automatically by shared_ptr; do not call
	template< typename X, typename Y >
	void
	_internal_accept_owner( shared_ptr<X> const * ppx, Y * py ) const
	{
		if( weak_this_.expired() )
			weak_this_ = shared_ptr<T>( *ppx, py );
	}

private:
	mutable weak_ptr<T>  weak_this_;
};  // enable_shared_from_this<>


//  enable_shared_from_this2.hpp

namespace detail
{

class esft2_deleter_wrapper
{
private:
	shared_ptr<void>  deleter_;

public:
	esft2_deleter_wrapper(){ }

	template< typename T >
	void
	set_deleter( shared_ptr<T> const & deleter )
	{
		deleter_ = deleter;
	}

	template< typename T >
	void
	operator () ( T* )
	{
		assert( deleter_.use_count() <= 1 );
		deleter_.reset();
	}
};

} // namespace detail

template< typename T >
class enable_shared_from_this2
{
protected:

	enable_shared_from_this2(){ }

	enable_shared_from_this2( enable_shared_from_this2 const & )
	{ }

	enable_shared_from_this2 & operator = ( enable_shared_from_this2 const & )
	{
		return *this;
	}

	~enable_shared_from_this2()
	{
		assert( shared_this_.use_count() <= 1 ); // ensure no dangling shared_ptrs
	}

private:
	mutable  weak_ptr<T>    weak_this_;
	mutable  shared_ptr<T>  shared_this_;

public:

	shared_ptr<T>
	shared_from_this()
	{
		init_weak_once();
		return shared_ptr<T>( weak_this_ );
	}

	shared_ptr<T const>
	shared_from_this() const
	{
		init_weak_once();
		return shared_ptr<T>( weak_this_ );
	}

private:

	void init_weak_once() const
	{
		if( weak_this_._empty() )
		{
			shared_this_.reset( static_cast< T* >( 0 )
					, detail::esft2_deleter_wrapper()
			);
			weak_this_ = shared_this_;
		}
	}

public:  // actually private, but avoids compiler template friendship issues

// Note: invoked automatically by shared_ptr; do not call
	template< typename X, typename Y >
	void
	_internal_accept_owner( shared_ptr<X> * ppx, Y * py ) const
	{
		assert( ppx != 0 );

		if( weak_this_.use_count() == 0 )
			weak_this_ = shared_ptr<T>( *ppx, py );
		else if( shared_this_.use_count() != 0 )
		{
			assert( ppx->unique() ); // no weak_ptrs should exist either, but there's no way to check that

			detail::esft2_deleter_wrapper *  pd
			= //boost::
					get_deleter<detail::esft2_deleter_wrapper>( shared_this_ );
			assert( pd != 0 );

			pd->set_deleter( *ppx );

			ppx->reset( shared_this_, ppx->get() );
			shared_this_.reset();
		}
	}
};  // enable_shared_from_this2<>
}

namespace G4
{
	using CLHEP::enable_shared_from_this;
	using CLHEP::enable_shared_from_this2;
}
#endif
*/

#endif /* G4SHARED_PTR_HH_ */

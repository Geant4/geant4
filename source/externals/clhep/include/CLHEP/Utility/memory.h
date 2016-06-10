#ifndef CLHEP_MEMORY_H
#define CLHEP_MEMORY_H

// ======================================================================
//
// memory - memory management utilities
//
// Note:  the following adaptation of the C++0X std::shared_ptr/weak_ptr
// interface and semantics has been customized for the specific internal
// needs of CLHEP/Random; it neither has nor needs the full generality
// of its namesake.
//
// Author:  W. E. Brown, 2010-03-19, adapted from the boost library's
// shared_ptr and related functionality whose internal attributions bear
// the following various notices:
//
//   (C) Copyright Greg Colvin and Beman Dawes 1998, 1999.
//   Copyright (c) 2001, 2002, 2003 Peter Dimov
//   Copyright (c) 2001, 2002, 2003 Peter Dimov and Multi Media Ltd.
//   Copyright (c) 2001-2008 Peter Dimov
//   Copyright (c) 2001-2009 Peter Dimov
//   Copyright 2002, 2009 Peter Dimov
//   Copyright 2004-2005 Peter Dimov
//   Copyright 2004-2008 Peter Dimov
//   Copyright 2005, 2006 Peter Dimov
//   Copyright 2008 Frank Mori Hess
//   Copyright 2008 Peter Dimov
//   Distributed under the Boost Software License, Version 1.0.
//   See http://www.boost.org/LICENSE_1_0.txt
//
// ======================================================================

// don't generate unnecessary warnings
#if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ > 6)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdelete-non-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#elif __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdelete-non-virtual-dtor"
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "CLHEP/Utility/defs.h"
#include "CLHEP/Utility/noncopyable.h"
#include "CLHEP/Utility/type_traits.h"

#include <algorithm>  // for swap
#include <cassert>    // for assert macro
#include <cstddef>    // for size_t
#include <exception>  // for exception
#include <functional> // for less
#include <iosfwd>     // for basic_ostream
#include <memory>     // for allocator, auto_ptr
#include <typeinfo>   // for bad_cast, type_info


namespace CLHEP {


// ----------------------------------------------------------------------
// forward declarations
// ----------------------------------------------------------------------

template< typename T > class shared_ptr;
template< typename T > class weak_ptr;
template< typename T > class enable_shared_from_this;
template< typename T > class enable_shared_from_this2;


// ----------------------------------------------------------------------
// bad_weak_ptr - exception thrown when a stale weak_ptr is encountered
// ----------------------------------------------------------------------

class bad_weak_ptr
  : public std::exception
{
public:
  inline  virtual  char const *  what() const throw();

};  // bad_weak_ptr

char const *
  bad_weak_ptr::what() const throw()
{
  return "bad_weak_ptr";
}


namespace sp {


// ----------------------------------------------------------------------
// abstract_ctrl_block - shared_ptr's counters and type-erased deleter
// ----------------------------------------------------------------------

class abstract_ctrl_block
  : public noncopyable
{
public:
  inline  void  class_invariant() const throw();
  // class class_invariant

  inline           abstract_ctrl_block();
  inline  virtual  ~abstract_ctrl_block() throw();
  // constructor and destructor

  inline           void    add_ref();
  inline           bool    add_ref_lock();
  inline           void    weak_add_ref() throw();
          virtual  void *  get_deleter( std::type_info const & ti ) = 0;
  inline           void    release() throw();
  inline           void    weak_release() throw();
          virtual  void    dispose() throw() = 0;
  inline  virtual  void    destroy() throw();
  // resource management functions

  inline  long use_count() const throw();
  // accessor

private:
  int  n_shared_ptrs;
  int  n_weak_ptrs;

};  // abstract_ctrl_block

void
  abstract_ctrl_block::class_invariant() const throw()
{
  assert( n_shared_ptrs == 0  || n_weak_ptrs >= 1 );
}

abstract_ctrl_block::abstract_ctrl_block()
  : n_shared_ptrs( 1 )
  , n_weak_ptrs  ( 1 )
{
  class_invariant();
}

abstract_ctrl_block::~abstract_ctrl_block() throw()
{
  class_invariant();
}

void
  abstract_ctrl_block::add_ref()
{
  class_invariant();
  ++n_shared_ptrs;
}

bool
  abstract_ctrl_block::add_ref_lock()  // true on success
{
  class_invariant();
  return n_shared_ptrs ? (++n_shared_ptrs != 0) : false;
}

void
  abstract_ctrl_block::weak_add_ref() throw()
{
  class_invariant();
  ++n_weak_ptrs;
}

void
  abstract_ctrl_block::release() throw()
{
  class_invariant();
  if( 0 == --n_shared_ptrs )
    dispose(), weak_release();
}

void
  abstract_ctrl_block::weak_release() throw()
{
  class_invariant();
  if( 0 == --n_weak_ptrs )
    destroy();
}

void
  abstract_ctrl_block::destroy() throw()
{
  assert( n_weak_ptrs == 0 );
  delete this;
}

long
  abstract_ctrl_block::use_count() const throw()
{
  class_invariant();
  return n_shared_ptrs;
}


// ----------------------------------------------------------------------
// concrete ctrl_block_* variations:
//   ctrl_block_p  : owned pointer only; no deleter, no allocator
//   ctrl_block_pd : owned pointer and deleter only; no allocator
//   ctrl_block_pda: owned pointer, deleter, and allocator
// ----------------------------------------------------------------------

template< typename P >  // P is pointee type
  class ctrl_block_p
    : public abstract_ctrl_block
{
  typedef ctrl_block_p<P> this_type;

public:
  inline  explicit ctrl_block_p( P * );
  inline           ~ctrl_block_p() throw();
  // constructor and destructor

  inline  void *  operator new ( std::size_t );
  inline  void    operator delete ( void * );
  // allocation functions

  inline  virtual  void *  get_deleter( std::type_info const & );
  inline  virtual  void    dispose() throw();
  // resource management functions

private:
  P *  owned_ptr;

};  // ctrl_block_p

template< typename P >
ctrl_block_p<P>::ctrl_block_p( P * p )
  : abstract_ctrl_block( )
  , owned_ptr( p )
{ }

template< typename P >
ctrl_block_p<P>::~ctrl_block_p() throw()
{ }

template< typename P >
void
  ctrl_block_p<P>::dispose() throw()
{
  delete owned_ptr;
}

template< typename P >
void *
  ctrl_block_p<P>::get_deleter( std::type_info const & )
{
  return 0;
}

template< typename P >
void *
  ctrl_block_p<P>::operator new ( std::size_t )
{
  return std::allocator<this_type>().allocate( 1 );
}

template< typename P >
void
  ctrl_block_p<P>::operator delete ( void * p )
{
  std::allocator<this_type>().deallocate( static_cast<this_type*>(p), 1 );
}

template< typename P  // pointee type
        , typename D  // deleter type
        >
  class ctrl_block_pd
    : public abstract_ctrl_block
{
  typedef ctrl_block_pd<P,D> this_type;

public:
  inline  ctrl_block_pd( P *, D );
  inline  ~ctrl_block_pd() throw();
  // constructor and destructor

  inline  void *  operator new ( std::size_t );
  inline  void    operator delete ( void * );
  // allocation functions

  inline  virtual  void *  get_deleter( std::type_info const & );
  inline  virtual  void    dispose() throw();
  // resource management functions

private:
  P *  owned_ptr;
  D    deleter;   // D's copy constructor must not throw, and
                  // call to deleter( owned_ptr ) must not throw

};  // ctrl_block_pd

template< typename P, typename D >
ctrl_block_pd<P,D>::ctrl_block_pd( P * p, D d )
  : abstract_ctrl_block( )
  , owned_ptr( p )
  , deleter  ( d )
{ }

template< typename P, typename D >
ctrl_block_pd<P,D>::~ctrl_block_pd() throw()
{ }

template< typename P, typename D >
void
  ctrl_block_pd<P,D>::dispose() throw()
{
  deleter( owned_ptr );
}

template< typename P, typename D >
void *
  ctrl_block_pd<P,D>::get_deleter( std::type_info const & ti )
{
  return ti == typeid(D) ? &reinterpret_cast<char&>( deleter ) : 0;
}

template< typename P, typename D >
void *
  ctrl_block_pd<P,D>::operator new ( std::size_t )
{
  return std::allocator<this_type>().allocate( 1 );
}

template< typename P, typename D >
void
  ctrl_block_pd<P,D>::operator delete ( void * p )
{
  std::allocator<this_type>().deallocate( static_cast<this_type*>(p), 1 );
}

template< typename P  // pointee type
        , typename D  // deleter type
        , typename A  // allocator type
        >
  class ctrl_block_pda
    : public abstract_ctrl_block
{
  typedef ctrl_block_pda<P,D,A> this_type;

public:
  inline  ctrl_block_pda( P *, D, A );
  inline  ~ctrl_block_pda() throw();
  // constructor and destructor

  inline  virtual  void *  get_deleter( std::type_info const & );
  inline  virtual  void    dispose() throw();
  inline  virtual  void    destroy() throw();
  // resource management functions

private:
  P *  owned_ptr;
  D    deleter;   // D's copy constructor must not throw, and
                  // call to deleter( owned_ptr ) must not throw
  A    allocator; // A's copy constructor must not throw

};  // ctrl_block_pda

template< typename P, typename D, typename A >
ctrl_block_pda<P,D,A>::ctrl_block_pda( P * p, D d, A a )
  : abstract_ctrl_block( )
  , owned_ptr( p )
  , deleter  ( d )
  , allocator( a )
{ }

template< typename P, typename D, typename A >
ctrl_block_pda<P,D,A>::~ctrl_block_pda() throw()
{ }

template< typename P, typename D, typename A >
void
  ctrl_block_pda<P,D,A>::dispose() throw()
{
  deleter( owned_ptr );
}

template< typename P, typename D, typename A >
void
  ctrl_block_pda<P,D,A>::destroy() throw()
{
  typename A::template rebind< this_type >::other  this_allocator( allocator );

  this_allocator.destroy( this );  // this->~this_type();
  this_allocator.deallocate( this, 1 );
}

template< typename P, typename D, typename A >
void *
  ctrl_block_pda<P,D,A>::get_deleter( std::type_info const & ti )
{
  return ti == typeid( D ) ? &reinterpret_cast<char&>( deleter ) : 0;
}


// ----------------------------------------------------------------------
// shared_ctrl_handle, weak_ctrl_handle - ctrl block handles
// ----------------------------------------------------------------------

class shared_ctrl_handle;
class weak_ctrl_handle;

struct sp_nothrow_tag { };

class shared_ctrl_handle
{
  friend  class weak_ctrl_handle;

public:
  inline  shared_ctrl_handle() throw();
  template< typename P >
    inline  explicit
    shared_ctrl_handle( P * );
  template< typename P, typename D >
    inline  shared_ctrl_handle( P *, D );
  template< typename P, typename D, typename A >
    inline  shared_ctrl_handle( P *, D, A );
  template< typename P >
    inline  explicit
    shared_ctrl_handle( std::auto_ptr<P> & );
  inline  ~shared_ctrl_handle() throw();
  // constructors and destructor

  inline  void  swap( shared_ctrl_handle & ) throw();
  inline  shared_ctrl_handle( shared_ctrl_handle const & ) throw();
  inline  shared_ctrl_handle &
    operator = ( shared_ctrl_handle const & ) throw();
  // copy functions

  inline  explicit
    shared_ctrl_handle( weak_ctrl_handle const & );
  inline  shared_ctrl_handle( weak_ctrl_handle const &, sp_nothrow_tag );
  // copy-like functions

  inline  void *  get_deleter( std::type_info const & ) const;
  inline  bool    unique() const throw();
  inline  bool    empty() const throw();
  inline  long    use_count() const throw();
  // accessors

  friend inline
  bool
    operator == ( shared_ctrl_handle const &, shared_ctrl_handle const & );
  friend inline
  bool
    operator < ( shared_ctrl_handle const &, shared_ctrl_handle const & );
  // comparisons

private:
  abstract_ctrl_block *  acb_ptr;

};  // shared_ctrl_handle

shared_ctrl_handle::shared_ctrl_handle() throw()
  : acb_ptr( 0 )
{ }

template< typename P >
  shared_ctrl_handle::shared_ctrl_handle( P * p )
  // a fctn-try block would be slightly more efficient here,
  // but some older compilers don't understand it
  : acb_ptr( 0 )
{
  try  {
    acb_ptr = new ctrl_block_p<P>(p);
  }
  catch(...)  {
    delete p;
    throw;
  }
}

template< typename P, typename D >
  shared_ctrl_handle::shared_ctrl_handle( P * p, D d )
  // a fctn-try block would be slightly more efficient here,
  // but some older compilers don't understand it
  : acb_ptr( 0 )
{
  try  {
    acb_ptr = new ctrl_block_pd<P,D>(p, d);
  }
  catch(...)  {
    d( p );
    throw;
  }
}

template< typename P, typename D, typename A >
  shared_ctrl_handle::shared_ctrl_handle( P * p, D d, A a )
  : acb_ptr( 0 )
{
  typedef  ctrl_block_pda<P,D,A>
           ctrl_block;
  typedef  typename A::template rebind<ctrl_block>::other
           ctrl_block_allocator;
  ctrl_block_allocator cba( a );

  try
  {
    acb_ptr = cba.allocate( 1 );
    new( static_cast<void*>(acb_ptr) ) ctrl_block(p, d, a);
  }
  catch(...)
  {
    d( p );
    if( acb_ptr != 0 )
      cba.deallocate( static_cast<ctrl_block*>( acb_ptr ), 1 );
    throw;
  }
}

template< typename P >
  shared_ctrl_handle::shared_ctrl_handle( std::auto_ptr<P> & p )
    : acb_ptr( new ctrl_block_p<P>( p.get() ) )
{
  p.release();
}

shared_ctrl_handle::~shared_ctrl_handle() throw()
{
  if( acb_ptr != 0 )
    acb_ptr->release();
}

void
  shared_ctrl_handle::swap( shared_ctrl_handle & other ) throw()
{
  abstract_ctrl_block * tmp = other.acb_ptr;
  other.acb_ptr = acb_ptr;
  acb_ptr = tmp;
}

shared_ctrl_handle::shared_ctrl_handle( shared_ctrl_handle const & other ) throw()
  : acb_ptr( other.acb_ptr )
{
  if( acb_ptr != 0 )
    acb_ptr->add_ref();
}

shared_ctrl_handle &
  shared_ctrl_handle::operator = ( shared_ctrl_handle const & other ) throw()
{
  abstract_ctrl_block * tmp = other.acb_ptr;

  if( tmp != acb_ptr )
  {
    if( tmp     != 0 ) tmp->add_ref();
    if( acb_ptr != 0 ) acb_ptr->release();
    acb_ptr = tmp;
  }

  return *this;
}

void *
  shared_ctrl_handle::get_deleter( std::type_info const & ti ) const
{
  return acb_ptr ? acb_ptr->get_deleter( ti ) : 0;
}

bool
  shared_ctrl_handle::unique() const throw()
{
  return 1L == use_count();
}

bool
  shared_ctrl_handle::empty() const throw()
{
  return acb_ptr == 0;
}

long
  shared_ctrl_handle::use_count() const throw()
{
  return acb_ptr == 0 ? 0L : acb_ptr->use_count();
}

bool
  operator == ( shared_ctrl_handle const & lhs, shared_ctrl_handle const & rhs )
{
  return lhs.acb_ptr == rhs.acb_ptr;
}

bool
  operator < ( shared_ctrl_handle const & lhs, shared_ctrl_handle const & rhs )
{
  return std::less<abstract_ctrl_block*>()( lhs.acb_ptr, rhs.acb_ptr );
}

class weak_ctrl_handle
{
  friend  class shared_ctrl_handle;

public:

  inline  weak_ctrl_handle() throw();
  inline  weak_ctrl_handle( shared_ctrl_handle const & ) throw();
  inline  ~weak_ctrl_handle() throw();
  // constructors and destructor

  inline  void  swap( weak_ctrl_handle & ) throw();
  inline  weak_ctrl_handle( weak_ctrl_handle const & ) throw();
  inline  weak_ctrl_handle & operator = ( shared_ctrl_handle const & ) throw();
  // copy functions

  inline  weak_ctrl_handle & operator = ( weak_ctrl_handle const & ) throw();
  // copy-like functions

  inline  bool  empty() const throw();
  inline  long  use_count() const throw();
  // accessors

  friend inline
  bool
    operator == ( weak_ctrl_handle const &, weak_ctrl_handle const & );
  friend inline
  bool
    operator < ( weak_ctrl_handle const &, weak_ctrl_handle const & );
  // comparisons

private:
  abstract_ctrl_block *  acb_ptr;

};  // weak_ctrl_handle

weak_ctrl_handle::weak_ctrl_handle() throw()
  : acb_ptr( 0 )
{ }

weak_ctrl_handle::weak_ctrl_handle( shared_ctrl_handle const & other ) throw()
  : acb_ptr( other.acb_ptr )
{
  if( acb_ptr != 0 )
    acb_ptr->weak_add_ref();
}

weak_ctrl_handle::~weak_ctrl_handle() throw()
{
  if( acb_ptr != 0 )
    acb_ptr->weak_release();
}

void
  weak_ctrl_handle::swap( weak_ctrl_handle & other ) throw()
{
  abstract_ctrl_block *  tmp = other.acb_ptr;
  other.acb_ptr = acb_ptr;
  acb_ptr = tmp;
}

weak_ctrl_handle::weak_ctrl_handle( weak_ctrl_handle const & other ) throw()
  : acb_ptr( other.acb_ptr )
{
  if( acb_ptr != 0 )
    acb_ptr->weak_add_ref();
}

weak_ctrl_handle &
  weak_ctrl_handle::operator = ( shared_ctrl_handle const & other ) throw()
{
  abstract_ctrl_block *  tmp = other.acb_ptr;

  if( tmp != acb_ptr )
  {
    if( tmp     != 0 ) tmp->weak_add_ref();
    if( acb_ptr != 0 ) acb_ptr->weak_release();
    acb_ptr = tmp;
}

  return *this;
}

weak_ctrl_handle &
  weak_ctrl_handle::operator = ( weak_ctrl_handle const & other ) throw()
{
  abstract_ctrl_block *  tmp = other.acb_ptr;

  if( tmp != acb_ptr )
{
    if( tmp     != 0 ) tmp->weak_add_ref();
    if( acb_ptr != 0 ) acb_ptr->weak_release();
    acb_ptr = tmp;
}

  return *this;
}

bool
  weak_ctrl_handle::empty() const throw()
{
  return acb_ptr == 0;
}

long
  weak_ctrl_handle::use_count() const throw()
{
  return acb_ptr == 0 ? 0L : acb_ptr->use_count();
}

bool
  operator == ( weak_ctrl_handle const & lhs, weak_ctrl_handle const & rhs )
{
  return lhs.acb_ptr == rhs.acb_ptr;
}

bool
  operator < ( weak_ctrl_handle const & lhs, weak_ctrl_handle const & rhs )
{
  return std::less<abstract_ctrl_block*>()( lhs.acb_ptr, rhs.acb_ptr );
}

shared_ctrl_handle::shared_ctrl_handle( weak_ctrl_handle const & other )
  : acb_ptr( other.acb_ptr )
{
  if( acb_ptr == 0  ||  ! acb_ptr->add_ref_lock() )
    throw bad_weak_ptr();
}

shared_ctrl_handle::shared_ctrl_handle( weak_ctrl_handle const & other
                                      , sp_nothrow_tag )
  : acb_ptr( other.acb_ptr )
{
  if( acb_ptr != 0  &&  ! acb_ptr->add_ref_lock() )
    acb_ptr = 0;
}


// ----------------------------------------------------------------------
// cast tags
// ----------------------------------------------------------------------

struct static_cast_tag      { };
struct const_cast_tag       { };
struct dynamic_cast_tag     { };
struct polymorphic_cast_tag { };


// ----------------------------------------------------------------------
// shared_ptr_traits - specify dependent types
// ----------------------------------------------------------------------

template< typename T >
  struct shared_ptr_traits
{
  typedef  T &  reference;
};

template<>
  struct shared_ptr_traits<void>
{
  typedef  void  reference;
};

template<>
  struct shared_ptr_traits<void const>
{
  typedef  void  reference;
};

template<>
  struct shared_ptr_traits<void volatile>
{
  typedef  void  reference;
};

template<>
  struct shared_ptr_traits<void const volatile>
{
  typedef  void  reference;
};


// ----------------------------------------------------------------------
// enable_shared_from_this support
// ----------------------------------------------------------------------

template< typename X, typename Y, typename T >
inline void
  sp_enable_shared_from_this( shared_ptr<X>              const * ppx
                            , Y                          const * py
                            , enable_shared_from_this<T> const * pe
                            )
{
  if( pe != 0 )
    pe->_internal_accept_owner( ppx, const_cast<Y*>( py ) );
}

template< typename X, typename Y, typename T >
inline void
  sp_enable_shared_from_this( shared_ptr<X>                     * ppx
                            , Y                           const * py
                            , enable_shared_from_this2<T> const * pe
                            )
{
  if( pe != 0 )
    pe->_internal_accept_owner( ppx, const_cast<Y*>( py ) );
}

inline void
  sp_enable_shared_from_this( ... )
{ }

}  // namespace sp


// ----------------------------------------------------------------------
// shared_ptr - "if you are the last person, please turn out the light"
// ----------------------------------------------------------------------

template< typename P >  // pointee type
  class shared_ptr
{
  typedef  shared_ptr<P>                                 this_type;
  typedef  typename sp::shared_ptr_traits<P>::reference  reference;

  template< typename >  friend  class shared_ptr;
  template< typename >  friend  class weak_ptr;

public:
  typedef  P  element_type;
  // pointee type

  shared_ptr() throw();
  template< typename P2 >
    inline  explicit
    shared_ptr( P2 * );
  template< typename P2, typename D >
    inline  shared_ptr( P2 *, D );
  template< typename P2, typename D, typename A >
    inline  shared_ptr( P2 *, D, A );
  // constructors

  inline  void          swap( shared_ptr<P> & ) throw();
  inline  shared_ptr &  operator = ( shared_ptr const & ) throw();
  // copy functions; generated copy constructor, destructor are fine

  template< typename P2 >
    inline  explicit
    shared_ptr( weak_ptr<P2> const & );
  template< typename P2 >
    inline  shared_ptr( weak_ptr<P2> const &, sp::sp_nothrow_tag ) throw();
  template< typename P2 >
    inline  shared_ptr( shared_ptr<P2> const &, P * ) throw();
  template< typename P2 >
    inline  shared_ptr( shared_ptr<P2> const &, sp::static_cast_tag );
  template< typename P2 >
    inline  shared_ptr( shared_ptr<P2> const &, sp::const_cast_tag );
  template< typename P2 >
    inline  shared_ptr( shared_ptr<P2> const &, sp::dynamic_cast_tag );
  template< typename P2 >
    inline  shared_ptr( shared_ptr<P2> const &, sp::polymorphic_cast_tag );
  template< typename P2 >
    inline  explicit
    shared_ptr( std::auto_ptr<P2> & );
  template< typename AP >
    inline  explicit
    shared_ptr( AP
              , typename enable_if_auto_ptr<AP,void*>::type = 0
              );
  template< typename P2 >
    inline
    shared_ptr( shared_ptr<P2> const &
              , typename enable_if_ptr_convertible<P2,P,void*>::type = 0
              ) throw();
  template< typename P2 >
    inline  shared_ptr &  operator = ( shared_ptr<P2> const & ) throw();
  template< typename P2 >
    inline  shared_ptr &  operator = ( std::auto_ptr<P2> & );
  template< typename AP >
    inline  typename enable_if_auto_ptr< AP, shared_ptr & >::type
    operator = ( AP );
  // copy-like functions

  inline  void reset() throw();
  template< typename P2 >
    inline  void  reset( P2 * );
  template< typename P2, typename D >
    inline  void  reset( P2 *, D );
  template< typename P2, typename D, typename A >
    inline  void  reset( P2 *, D, A );
  template< typename P2 >
    inline  void  reset( shared_ptr<P2> const &, P * );
  // reset functions

             inline  operator bool () const throw();
  inline  reference  operator *    () const throw();
  inline  P *        operator ->   () const throw();
  // pointer-like behavior

  inline  P *      get() const throw();
  inline  bool     unique() const throw();
  inline  long     use_count() const throw();
  // accessors

  template< typename P2 >
    inline  bool _internal_less( shared_ptr<P2> const & ) const;
  inline  void *  _internal_get_deleter( std::type_info const & ) const;
  inline  bool    _internal_equiv( shared_ptr const & ) const;
  // implementation helpers -- do not use

private:
  P *                    px;  // contained pointer
  sp::shared_ctrl_handle pn;  // control information

};  // shared_ptr

template< typename P, typename P2 >
  inline bool  operator == ( shared_ptr<P> const &, shared_ptr<P2> const & );
template< typename P, typename P2 >
  inline bool  operator != ( shared_ptr<P> const &, shared_ptr<P2> const & );
template< typename P, typename P2 >
  inline bool  operator <  ( shared_ptr<P> const &, shared_ptr<P2> const & );

template< typename P >
  inline void  swap( shared_ptr<P> &, shared_ptr<P> & );

template< typename P, typename P2 >
  inline shared_ptr<P>  static_pointer_cast( shared_ptr<P2> const & );
template< typename P, typename P2 >
  inline shared_ptr<P>  const_pointer_cast( shared_ptr<P2> const & );
template< typename P, typename P2 >
  inline shared_ptr<P>  dynamic_pointer_cast( shared_ptr<P2> const & );

template< typename P >
  inline P *  get_pointer( shared_ptr<P> const & );
template< typename D, typename P >
  inline D *  get_deleter( shared_ptr<P> const & );

template< typename C, typename T, typename P >
  inline std::basic_ostream<C,T> &  operator << ( std::basic_ostream<C,T> &
                                                , shared_ptr<P> const &
                                                );

template< typename P >
  shared_ptr<P>::shared_ptr() throw()
  : px( 0 )
  , pn(         )
{ }

template< typename P >
template< typename P2 >  // P2 must be a complete type
  shared_ptr<P>::shared_ptr( P2 * p )
  : px( p )
  , pn( p )
{
  sp::sp_enable_shared_from_this( this, p, p );
}

template< typename P >
template< typename P2, typename D >  // D's copy c'tor must not throw
  shared_ptr<P>::shared_ptr( P2 * p, D d )
  : px( p    )
  , pn( p, d )
{
  sp::sp_enable_shared_from_this( this, p, p );
}

template< typename P >
template< typename P2, typename D, typename A >  // D's, A's copy c'tors must not throw
  shared_ptr<P>::shared_ptr( P2 * p, D d, A a )
  : px( p       )
  , pn( p, d, a )
{
  sp::sp_enable_shared_from_this( this, p, p );
}

template< typename P >
  void
  shared_ptr<P>::swap( shared_ptr<P> & other ) throw()
{
  std::swap( px, other.px );
  pn.swap( other.pn );
}

template< typename P >
  shared_ptr<P> &
  shared_ptr<P>::operator = ( shared_ptr const & other ) throw()
{
  this_type( other ).swap( *this );
  return *this;
}

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( weak_ptr<P2> const & other )
  : px( 0  )  // temporarily
  , pn( other.pn )  // may throw
{
  px = other.px;  // safe to copy other.px, as pn(other.pn) did not throw
}

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( weak_ptr<P2> const & other
                           , sp::sp_nothrow_tag
                           ) throw()
  : px( 0                        )  // temporarily
  , pn( other.pn, sp::sp_nothrow_tag() )
{
  if( ! pn.empty() )
    px = other.px;
}

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( shared_ptr<P2> const & other
                           , P * p
                           ) throw()
  : px( p        )
  , pn( other.pn )
{ }

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( shared_ptr<P2> const & other
                           , sp::static_cast_tag
                           )
  : px( static_cast<element_type*>( other.px ) )
  , pn( other.pn                               )
{ }

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( shared_ptr<P2> const & other
                           , sp::const_cast_tag
                           )
  : px( const_cast<element_type*>( other.px ) )
  , pn( other.pn                              )
{ }

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( shared_ptr<P2> const & other
                           , sp::dynamic_cast_tag
                           )
  : px( dynamic_cast<element_type*>( other.px ) )
  , pn( other.pn                                )
{
  if( px == 0 )               // cast failed?
    pn = sp::shared_ctrl_handle();  // yes; need our own control information
}

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( shared_ptr<P2> const & other
                           , sp::polymorphic_cast_tag
                           )
  : px( dynamic_cast<element_type*>( other.px ) )
  , pn( other.pn                                )
{
  if( px == 0 )
    throw std::bad_cast();
}

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( std::auto_ptr<P2> & other )
  : px( other.get() )
  , pn(             )  // temporarily
{
  P2 *  tmp = other.get();
  pn = sp::shared_ctrl_handle( other );
  sp::sp_enable_shared_from_this( this, tmp, tmp );
}

template< typename P >
template< typename AP >
  shared_ptr<P>::shared_ptr( AP other
                           , typename enable_if_auto_ptr<AP,void*>::type
                           )
  : px( other.get() )
  , pn(             )  // temporarily
{
  typename AP::element_type *  tmp = other.get();
  pn = sp::shared_ctrl_handle( other );
  sp::sp_enable_shared_from_this( this, tmp, tmp );
}

template< typename P >
template< typename P2 >
  shared_ptr<P>::shared_ptr( shared_ptr<P2> const & other
                           , typename enable_if_ptr_convertible<P2,P,void*>::type
                           ) throw()
  : px( other.px )
  , pn( other.pn )
  { }

template< typename P >
template< typename P2 >
  shared_ptr<P> &
  shared_ptr<P>::operator = ( shared_ptr<P2> const & other ) throw()
{
  this_type( other ).swap( *this );
  return *this;
}

template< typename P >
template< typename P2 >
  shared_ptr<P> &
  shared_ptr<P>::operator = ( std::auto_ptr<P2> & other )
{
  this_type( other ).swap( *this );
  return *this;
}

template< typename P >
template< typename AP >
  typename enable_if_auto_ptr< AP, shared_ptr<P> & >::type
  shared_ptr<P>::operator = ( AP other )
{
  this_type( other ).swap( *this );
  return *this;
}

template< typename P >
  void
  shared_ptr<P>::reset() throw()
{
  this_type().swap( *this );
}

template< typename P >
template< typename P2 >
  void
  shared_ptr<P>::reset( P2 * p )  // P2 must be a complete type
{
  assert( p == 0 || p != px );  // oughtn't reset oneself
  this_type( p ).swap( *this );
}

template< typename P >
template< typename P2, typename D >
  void
  shared_ptr<P>::reset( P2 * p, D d )
{
  this_type( p, d ).swap( *this );
}

template< typename P >
template< typename P2, typename D, typename A >
  void
  shared_ptr<P>::reset( P2 * p, D d, A a )
{
  this_type( p, d, a ).swap( *this );
}

template< typename P >
template< typename P2 >
  void
  shared_ptr<P>::reset( shared_ptr<P2> const & other, P * p )
{
  this_type( other, p ).swap( *this );
}

template< typename P >
  shared_ptr<P>::operator bool () const throw()
{
  return px;
}

template< typename P >
  typename sp::shared_ptr_traits<P>::reference
  //typename shared_ptr<P>::reference
  shared_ptr<P>::operator * () const throw()
{
  assert( px != 0 );
  return *px;
}

template< typename P >
  P *
  shared_ptr<P>::operator -> () const throw()
{
  assert( px != 0 );
  return px;
}

template< typename P >
  P *
  shared_ptr<P>::get() const throw()
{
  return px;
}

template< typename P >
  bool
  shared_ptr<P>::unique() const throw()
{
  return pn.unique();
}

template< typename P >
  long
  shared_ptr<P>::use_count() const throw()
{
  return pn.use_count();
}

template< typename P >
template< typename P2 >
  bool
  shared_ptr<P>::_internal_less( shared_ptr<P2> const & rhs ) const
{
  return pn < rhs.pn;
}

template< typename P >
  void *
  shared_ptr<P>::_internal_get_deleter( std::type_info const & ti ) const
{
  return pn.get_deleter( ti );
}

template< typename P >
  bool
  shared_ptr<P>::_internal_equiv( shared_ptr const & other ) const
{
  return px == other.px && pn == other.pn;
}

template< typename P, typename P2 >
  bool
  operator == ( shared_ptr<P> const & a, shared_ptr<P2> const & b )
{
  return a.get() == b.get();
}

template< typename P, typename P2 >
  bool
  operator != ( shared_ptr<P> const & a, shared_ptr<P2> const & b )
{
  return a.get() != b.get();
}

template< typename P, typename P2 >
  bool
  operator < ( shared_ptr<P> const & a, shared_ptr<P2> const & b )
{
  return a._internal_less(b);
}

template< typename P >
  void
  swap( shared_ptr<P> & a, shared_ptr<P> & b )
{
  a.swap( b );
}

template< typename P, typename P2 >
  shared_ptr<P>
  static_pointer_cast( shared_ptr<P2> const & other )
{
  return shared_ptr<P>( other, sp::static_cast_tag() );
}

template< typename P, typename P2 >
  shared_ptr<P>
  const_pointer_cast( shared_ptr<P2> const & other )
{
  return shared_ptr<P>( other, sp::const_cast_tag() );
}

template< typename P, typename P2 >
  shared_ptr<P>
  dynamic_pointer_cast( shared_ptr<P2> const & other )
{
  return shared_ptr<P>( other, sp::dynamic_cast_tag() );
}

template< typename P >
  P *
  get_pointer( shared_ptr<P> const & p )
{
  return p.get();
}

template< typename D, typename P >
  D *
  get_deleter( shared_ptr<P> const & p )
{
  return static_cast<D*>( p._internal_get_deleter( typeid(D)) );
}

template< typename C, typename T, typename P >
  std::basic_ostream<C,T> &
  operator << ( std::basic_ostream<C,T> & os, shared_ptr<P> const & p )
{
  os << p.get();
  return os;
}


// ----------------------------------------------------------------------
// weak_ptr - non-owning handle from which a shared_ptr can be obtained
// ----------------------------------------------------------------------

template< typename P >
  class weak_ptr
{
  typedef  weak_ptr<P>  this_type;

  template< typename >  friend  class shared_ptr;
  template< typename >  friend  class weak_ptr;

public:
  typedef  P  element_type;

  inline  weak_ptr() throw();

  //  generated copy constructor, assignment, destructor are fine

  inline  void  swap( this_type & other ) throw();
  template< typename P2 >
  inline
  weak_ptr( weak_ptr<P2> const & r
          , typename enable_if_ptr_convertible<P2,P,void*>::type = 0
          ) throw();
  template< typename P2 >
  inline
  weak_ptr( shared_ptr<P2> const & r
          , typename enable_if_ptr_convertible<P2,P,void*>::type = 0
          ) throw();
  template< typename P2 >
  inline  weak_ptr &  operator = (weak_ptr<P2> const & r) throw();
  template< typename P2 >
  inline  weak_ptr &  operator = (shared_ptr<P2> const & r) throw();
  // copy-like functions

  inline  shared_ptr<P>  lock() const throw();
  inline  long           use_count() const throw();
  inline  bool           expired() const throw();
  inline  bool           _empty() const; // extension, not in std::weak_ptr
  inline  void           reset() throw();
  // accessors

  inline  void  _internal_assign( P * px2, sp::shared_ctrl_handle const & pn2 );
  template< typename P2 >
    inline  bool _internal_less( weak_ptr<P2> const & rhs ) const;

private:
  P *                   px;  // contained pointer
  sp::weak_ctrl_handle  pn;  // control information

};  // weak_ptr

template< typename P, typename P2 >
  inline  bool  operator < ( weak_ptr<P> const & a, weak_ptr<P2> const & b );

template< typename P >
  inline  void  swap( weak_ptr<P> & a, weak_ptr<P> & b );

template< typename P >
weak_ptr<P>::weak_ptr() throw()
  : px( 0 )
  , pn(         )
{ }

template< typename P >
template< typename P2 >
  weak_ptr<P>::weak_ptr( weak_ptr<P2> const & r
                       , typename enable_if_ptr_convertible<P2,P,void*>::type
                       ) throw()
  : px( r.lock().get() )  // same as r.px, but doesn't risk invalidation
  , pn( r.pn           )
{ }

template< typename P >
template< typename P2 >
  weak_ptr<P>::weak_ptr( shared_ptr<P2> const & r
                       , typename enable_if_ptr_convertible<P2,P,void*>::type
                       ) throw()
  : px( r.px )
  , pn( r.pn )
{ }

template< typename P >
template< typename P2 >
  weak_ptr<P> &
  weak_ptr<P>::operator = (weak_ptr<P2> const & r) throw()
{
  px = r.lock().get();
  pn = r.pn;
  return *this;
}

template< typename P >
template< typename P2 >
  weak_ptr<P> &
  weak_ptr<P>::operator = (shared_ptr<P2> const & r) throw()
{
  px = r.px;
  pn = r.pn;
  return *this;
}

template< typename P >
  shared_ptr<P>
  weak_ptr<P>::lock() const throw()
{
  return shared_ptr<element_type>( *this, sp::sp_nothrow_tag() );
}

template< typename P >
  long
  weak_ptr<P>::use_count() const throw()
{
  return pn.use_count();
}

template< typename P >
  bool
  weak_ptr<P>::expired() const throw()
{
  return pn.use_count() == 0;
}

template< typename P >
  bool
  weak_ptr<P>::_empty() const // extension, not in std::weak_ptr
{
  return pn.empty();
}

template< typename P >
  void
  weak_ptr<P>::reset() throw()
{
  this_type().swap(*this);
}

template< typename P >
  void
  weak_ptr<P>::swap( this_type & other ) throw()
{
  std::swap(px, other.px);
  pn.swap(other.pn);
}

template< typename P >
  void
  weak_ptr<P>::_internal_assign( P * px2, sp::shared_ctrl_handle const & pn2 )
{
  px = px2;
  pn = pn2;
}

template< typename P >
template< typename P2 >
  bool
  weak_ptr<P>::_internal_less( weak_ptr<P2> const & rhs ) const
{
  return pn < rhs.pn;
}

template< typename P, typename P2 >
  bool
  operator < ( weak_ptr<P> const & a, weak_ptr<P2> const & b )
{
  return a._internal_less(b);
}

template< typename P >
  void
  swap( weak_ptr<P> & a, weak_ptr<P> & b )
{
  a.swap(b);
}


// ----------------------------------------------------------------------
// do_nothing_deleter - for shared_ptrs not taking ownership
// ----------------------------------------------------------------------

struct do_nothing_deleter {
  inline  void  operator () ( void const * ) const;
};

void
do_nothing_deleter::operator () ( void const * ) const
{ }


}  // namespace CLHEP


#if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ > 6)
  #pragma GCC diagnostic pop
#elif __clang__
  #pragma clang diagnostic pop
#endif



#endif  // CLHEP_MEMORY_H
//
// ======================================================================


#if 0

//  enable_shared_from_this.hpp

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
    assert( p.get() == this );
    return p;
  }

  shared_ptr<T const>
    shared_from_this() const
  {
    shared_ptr<T const>  p( weak_this_ );
    assert( p.get() == this );
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
  esft2_deleter_wrapper()
  { }

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

  enable_shared_from_this2()
  { }

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
        = boost::get_deleter<detail::esft2_deleter_wrapper>( shared_this_ );
      assert( pd != 0 );

      pd->set_deleter( *ppx );

      ppx->reset( shared_this_, ppx->get() );
      shared_this_.reset();
    }
  }
};  // enable_shared_from_this2<>

#endif  // 0

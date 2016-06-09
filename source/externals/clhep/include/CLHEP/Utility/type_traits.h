#ifndef CLHEP_TYPE_TRAITS_H
#define CLHEP_TYPE_TRAITS_H

// ======================================================================
//
// type_traits - selected C++0X metaprogramming constructs
//
// Author:  W. E. Brown; 2010-03-05
//
// ======================================================================


#include "CLHEP/Utility/defs.h"

#include <memory>  // for auto_ptr

#if defined(__GXX_EXPERIMENTAL_CXX0X__)
#  define CLHEP_HAS_RVALUE_REFERENCE
#else
#  define CLHEP_NO_RVALUE_REFERENCE
#endif


namespace CLHEP  {

// ----------------------------------------------------------------------
// Contents:
// ----------------------------------------------------------------------

// helper class:
template< typename T, T v > struct integral_constant;
typedef  integral_constant<bool, true  >  true_type;
typedef  integral_constant<bool, false >  false_type;

// primary type categories:
template< typename T > struct is_void;
template< typename T > struct is_integral;
template< typename T > struct is_floating_point;
template< typename T > struct is_array;
template< typename T > struct is_pointer;
template< typename T > struct is_lvalue_reference;
template< typename T > struct is_rvalue_reference;
template< typename T > struct is_member_object_pointer;
template< typename T > struct is_member_function_pointer;
template< typename T > struct is_enum;
template< typename T > struct is_union;
template< typename T > struct is_class;
template< typename T > struct is_function;

// composite type categories:
template< typename T > struct is_reference;
template< typename T > struct is_arithmetic;
template< typename T > struct is_fundamental;
template< typename T > struct is_object;
template< typename T > struct is_scalar;
template< typename T > struct is_compound;
template< typename T > struct is_member_pointer;

// type properties:
template< typename T > struct is_const;
template< typename T > struct is_volatile;
#if 0
template< typename T > struct is_trivial;
template< typename T > struct is_trivially_copyable;
template< typename T > struct is_standard_layout;
template< typename T > struct is_pod;
template< typename T > struct is_literal_type;
template< typename T > struct is_empty;
template< typename T > struct is_polymorphic;
#endif  // 0
template< typename T > struct is_abstract;
#if 0
template< typename T, typename... Args > struct is_constructible;
template< typename T, typename... Args > struct is_nothrow_constructible;
template< typename T > struct has_default_constructor;
template< typename T > struct has_copy_constructor;
template< typename T > struct has_copy_assign;
template< typename T > struct has_move_constructor;
template< typename T > struct has_move_assign;
template< typename T > struct has_trivial_default_constructor;
template< typename T > struct has_trivial_copy_constructor;
template< typename T > struct has_trivial_move_constructor;
template< typename T > struct has_trivial_copy_assign;
template< typename T > struct has_trivial_move_assign;
template< typename T > struct has_trivial_destructor;
template< typename T > struct has_nothrow_default_constructor;
template< typename T > struct has_nothrow_copy_constructor;
template< typename T > struct has_nothrow_move_constructor;
template< typename T > struct has_nothrow_copy_assign;
template< typename T > struct has_nothrow_move_assign;
template< typename T > struct has_virtual_destructor;
#endif  // 0
template< typename T > struct is_signed;
template< typename T > struct is_unsigned;
#if 0
template< typename T > struct alignment_of;
#endif  // 0
template< typename T > struct rank;
template< typename T, unsigned I = 0 > struct extent;

// type relations:
template< typename T, typename U > struct is_same;
#if 0
template< typename Base, typename Derived > struct is_base_of;
#endif  // 0
template< typename From, typename To > struct is_convertible;
#if 0
template< typename From, typename To > struct is_explicitly_convertible;
#endif  // 0

// const-volatile modifications:
template< typename T > struct remove_const;
template< typename T > struct remove_volatile;
template< typename T > struct remove_cv;
template< typename T > struct add_const;
template< typename T > struct add_volatile;
template< typename T > struct add_cv;

// reference modifications:
template< typename T > struct remove_reference;
template< typename T > struct add_lvalue_reference;
template< typename T > struct add_rvalue_reference;

// sign modifications:
#if 0
template< typename T > struct make_signed;
template< typename T > struct make_unsigned;
#endif  // 0

// array modifications:
template< typename T > struct remove_extent;
template< typename T > struct remove_all_extents;

// pointer modifications:
template< typename T > struct remove_pointer;
template< typename T > struct add_pointer;

// other transformations:
#if 0
template< std::size_t Len, std::size_t Align > struct aligned_storage;
template< std::size_t Len, typename... Types > struct aligned_union;
template< typename T > struct decay;
#endif  // 0
template< bool, typename T = void > struct enable_if;
template< bool, typename T, typename F > struct conditional;
#if 0
template< typename... T > struct common_type;
template< typename T > struct underlying_type;
template< typename > typename result_of; // undefined
template< typename F, typename... ArgTypes > typename result_of<F(ArgTypes...)>;
#endif  // 0

// non-standard (but useful) extensions:
template< typename From, typename To > struct is_ptr_convertible;
template< typename From, typename To, typename R=void > struct enable_if_convertible;
template< typename From, typename To, typename R=void > struct enable_if_ptr_convertible;
template< typename P, typename R=void > struct enable_if_auto_ptr;


// ----------------------------------------------------------------------
// integral_constant - a helper class, useful in its own right
// ----------------------------------------------------------------------

template< typename T, T v >
  struct integral_constant
{
  typedef  T                       value_type;
  typedef  integral_constant<T,v>  type;

  static  value_type const  value  =  v;

  operator value_type()  { return value; }
};  // integral_constant<,>


// ----------------------------------------------------------------------
// yes_t, no_t - unimplemented types with distinct sizeof
// ----------------------------------------------------------------------

namespace tt  {

typedef  char (& yes_t);      // ref to char
typedef  char (& no_t ) [2];  // ref to 2-char array

}  // namespace tt


// ----------------------------------------------------------------------
// primary<,> - type classification helper
// ----------------------------------------------------------------------

namespace tt  {

enum primary_code
{ _unknown                 = 0u
, _void                    = 1u <<  0
, _integral                = 1u <<  1
, _floating_point          = 1u <<  2
, _array                   = 1u <<  3
, _pointer                 = 1u <<  4
, _lvalue_reference        = 1u <<  5
, _rvalue_reference        = 1u <<  6
, _member_object_pointer   = 1u <<  7
, _member_function_pointer = 1u <<  8
, _enum                    = 1u <<  9
, _union                   = 1u << 10  // Help, compiler!
, _class                   = 1u << 11
, _function                = 1u << 12
};  // primary_code

// Helpers to recognize classes:
template< typename U >  yes_t  isAclass( void(U::*)() );
template< typename U >  no_t   isAclass( ...          );

// Helpers to recognize functions:
template< typename U >
  no_t   isAfunction( U(*)[1] );  // arrays of non-{fctn/ref/void}s
template< typename U >
  yes_t  isAfunction( ...     );

// encode via helpers or by elimination:
//   enum
//   union  // need help, compiler!
//   class
//   function
template< typename T >
struct encode
{
  static  primary_code const  value
    =      ( sizeof(isAclass   <T>(0)) == sizeof(yes_t) )  ? _class
       : ( ( sizeof(isAfunction<T>(0)) == sizeof(yes_t) )  ? _function
       :   /* by elimination */                              _enum
         );
};  // encode<>

// encode cv-qualified type:
template< typename T >
  struct encode<T const> : public encode<T>  { };
template< typename T >
  struct encode<T volatile> : public encode<T>  { };
template< typename T >
  struct encode<T const volatile> : public encode<T>  { };

// encode array:
template< typename T >
  struct encode<T[]>
{ static  primary_code const  value  =  _array; };
template< typename T >
  struct encode<T const[]>
{ static  primary_code const  value  =  _array; };
template< typename T >
  struct encode<T volatile[]>
{ static  primary_code const  value  =  _array; };
template< typename T >
  struct encode<T const volatile[]>
{ static  primary_code const  value  =  _array; };
template< typename T, unsigned N >
  struct encode<T[N]>
{ static  primary_code const  value  =  _array; };
template< typename T, unsigned N >
  struct encode<T const[N]>
{ static  primary_code const  value  =  _array; };
template< typename T, unsigned N >
  struct encode<T volatile[N]>
{ static  primary_code const  value  =  _array; };
template< typename T, unsigned N >
  struct encode<T const volatile[N]>
{ static  primary_code const  value  =  _array; };

// encode floating_point:
template<>
  struct encode<float>
{ static  primary_code const  value  =  _floating_point; };
template<>
  struct encode<double>
{ static  primary_code const  value  =  _floating_point; };
template<>
  struct encode<long double>
{ static  primary_code const  value  =  _floating_point; };

// encode integral:
template<>
  struct encode<bool>
{ static  primary_code const  value  =  _integral; };
template<>
  struct encode<signed char>
{ static  primary_code const  value  =  _integral; };
template<>
  struct encode<char>
{ static  primary_code const  value  =  _integral; };
template<>
  struct encode<unsigned char>
{ static  primary_code const  value  =  _integral; };
#if 0
template<>
  struct encode<wchar_t>
{ static  primary_code const  value  =  _integral; };
#endif
template<>
  struct encode<short>
{ static  primary_code const  value  =  _integral; };
template<>
  struct encode<unsigned short>
{ static  primary_code const  value  =  _integral; };
template<>
  struct encode<int>
{ static  primary_code const  value  =  _integral; };
template<>
  struct encode<unsigned int>
{ static  primary_code const  value  =  _integral; };
template<>
  struct encode<long>
{ static  primary_code const  value  =  _integral; };
template<>
  struct encode<unsigned long>
{ static  primary_code const  value  =  _integral; };

// encode member_function_pointer:
template< typename T, typename C >
  struct encode<T (C::*)()>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C >
  struct encode<T (C::*)() const>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C >
  struct encode<T (C::*)() volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C >
  struct encode<T (C::*)() const volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1 >
  struct encode<T (C::*)(A1)>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1 >
  struct encode<T (C::*)(A1) const>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1 >
  struct encode<T (C::*)(A1) volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1 >
  struct encode<T (C::*)(A1) const volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2 >
  struct encode<T (C::*)(A1,A2)>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2 >
  struct encode<T (C::*)(A1,A2) const>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2 >
  struct encode<T (C::*)(A1,A2) volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2 >
  struct encode<T (C::*)(A1,A2) const volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3 >
  struct encode<T (C::*)(A1,A2,A3)>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3 >
  struct encode<T (C::*)(A1,A2,A3) const>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3 >
  struct encode<T (C::*)(A1,A2,A3) volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3 >
  struct encode<T (C::*)(A1,A2,A3) const volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3, typename A4 >
  struct encode<T (C::*)(A1,A2,A3,A4)>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3, typename A4 >
  struct encode<T (C::*)(A1,A2,A3,A4) const>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3, typename A4 >
  struct encode<T (C::*)(A1,A2,A3,A4) volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3, typename A4 >
  struct encode<T (C::*)(A1,A2,A3,A4) const volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3, typename A4, typename A5 >
  struct encode<T (C::*)(A1,A2,A3,A4,A5)>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3, typename A4, typename A5 >
  struct encode<T (C::*)(A1,A2,A3,A4,A5) const>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3, typename A4, typename A5 >
  struct encode<T (C::*)(A1,A2,A3,A4,A5) volatile>
{ static  primary_code const  value  =  _member_function_pointer; };
template< typename T, typename C
        , typename A1, typename A2, typename A3, typename A4, typename A5 >
  struct encode<T (C::*)(A1,A2,A3,A4,A5) const volatile>
{ static  primary_code const  value  =  _member_function_pointer; };

// encode member_object_pointer:
template< typename T, typename C >
  struct encode<T C::*>
{ static  primary_code const  value  =  _member_object_pointer; };

// encode pointer:
template< typename T >
  struct encode<T *>
{ static  primary_code const  value  =  _pointer; };

// encode lvalue_reference:
template< typename T >
  struct encode<T &>
{ static  primary_code const  value  =  _lvalue_reference; };

// encode rvalue_reference:
#if defined(CLHEP_HAS_RVALUE_REFERENCE)
template< typename T >
  struct encode<T&&>
{ static  primary_code const  value  =  _rvalue_reference; };
#endif  // CLHEP_HAS_RVALUE_REFERENCE

// encode void:
template<>
  struct encode<void>
{ static  primary_code const  value  =  _void; };

// apply encoding:
template< typename T, unsigned int p >
  struct primary : integral_constant<bool, bool(p & encode<T>::value)>  { };

}  // namespace tt


// ----------------------------------------------------------------------
// is_void - metaprogramming type trait detecting void types
// ----------------------------------------------------------------------

template< typename T >
  struct is_void
  : public tt::primary<T, tt::_void >  { };


// ----------------------------------------------------------------------
// is_integral - metaprogramming type trait detecting integer types
// ----------------------------------------------------------------------

template< typename T >
  struct is_integral
  : public tt::primary<T, tt::_integral >  { };


// ----------------------------------------------------------------------
// is_floating_point - metaprogramming type trait detecting real types
// ----------------------------------------------------------------------

template< typename T >
  struct is_floating_point
  : public tt::primary<T, tt::_floating_point >  { };

// ----------------------------------------------------------------------
// is_array - metaprogramming type trait detecting T[...] types
// ----------------------------------------------------------------------

template< typename T >
  struct is_array
  : public tt::primary<T, tt::_array >  { };


// ----------------------------------------------------------------------
// is_pointer - metaprogramming type trait detecting T* types
// ----------------------------------------------------------------------

template< typename T >
  struct is_pointer
  : public tt::primary<T, tt::_pointer >  { };


// ----------------------------------------------------------------------
// is_lvalue_reference - metaprogramming type trait detecting T& types
// ----------------------------------------------------------------------

template< typename T >
  struct is_lvalue_reference
  : public tt::primary<T, tt::_lvalue_reference >  { };


// ----------------------------------------------------------------------
// is_rvalue_reference - metaprogramming type trait detecting T&& types
// ----------------------------------------------------------------------

template< typename T >  struct is_rvalue_reference
  : public tt::primary<T, tt::_rvalue_reference >  { };


// ----------------------------------------------------------------------
// is_member_object_pointer - metaprogramming type trait
// ----------------------------------------------------------------------

template< typename T >  struct is_member_object_pointer
  : public conditional< is_member_function_pointer<T>::value
                      , false_type
                      , tt::primary<T, tt::_member_object_pointer>
                      >::type
{ };


// ----------------------------------------------------------------------
// is_member_function_pointer - metaprogramming type trait
// ----------------------------------------------------------------------

template< typename T >
  struct is_member_function_pointer
  : public tt::primary<T, tt::_member_function_pointer >  { };


// ----------------------------------------------------------------------
// is_enum - metaprogramming type trait detecting enumeration types
// ----------------------------------------------------------------------

template< typename T >
  struct is_enum
  : public tt::primary<T, tt::_enum >  { };


// ----------------------------------------------------------------------
// is_union - metaprogramming type trait detecting union types
// ----------------------------------------------------------------------

template< typename T >
  struct is_union
  : public tt::primary<T, tt::_union >  { };


// ----------------------------------------------------------------------
// is_class - metaprogramming type trait detecting class types
// ----------------------------------------------------------------------

template< typename T >
  struct is_class
  : public tt::primary<T, tt::_class >  { };


// ----------------------------------------------------------------------
// is_function - metaprogramming type trait detecting function types
// ----------------------------------------------------------------------

template< typename T >
  struct is_function
  : public tt::primary<T, tt::_function >  { };


// ----------------------------------------------------------------------
// is_reference - metaprogramming composite type trait
// ----------------------------------------------------------------------

template< typename T >
  struct is_reference
  : public tt::primary< T, tt::_lvalue_reference
                         | tt::_rvalue_reference
                         >
{ };


// ----------------------------------------------------------------------
// is_arithmetic - metaprogramming composite type trait
// ----------------------------------------------------------------------

template< typename T >
  struct is_arithmetic
  : public tt::primary< T, tt::_integral
                         | tt::_floating_point
                         >
{ };


// ----------------------------------------------------------------------
// is_fundamental - metaprogramming composite type trait
// ----------------------------------------------------------------------

template< typename T >
  struct is_fundamental
  : public tt::primary< T, tt::_integral
                         | tt::_floating_point
                         | tt::_void
                         >
{ };


// ----------------------------------------------------------------------
// is_object - metaprogramming composite type trait
// ----------------------------------------------------------------------

template< typename T >
  struct is_object
  : public tt::primary< T, tt::_array
                         | tt::_class
                         | tt::_enum
                         | tt::_floating_point
                         | tt::_integral
                         | tt::_member_object_pointer
                         | tt::_member_function_pointer
                         | tt::_pointer
                         | tt::_union
                         >
{ };


// ----------------------------------------------------------------------
// is_scalar - metaprogramming composite type trait
// ----------------------------------------------------------------------

template< typename T >
  struct is_scalar
  : public tt::primary< T, tt::_integral
                         | tt::_floating_point
                         | tt::_enum
                         | tt::_pointer
                         | tt::_member_object_pointer
                         | tt::_member_function_pointer
                         >
{ };


// ----------------------------------------------------------------------
// is_compound - metaprogramming composite type trait
// ----------------------------------------------------------------------

template< typename T >
  struct is_compound
  : public tt::primary< T, tt::_array
                         | tt::_pointer
                         | tt::_lvalue_reference
                         | tt::_rvalue_reference
                         | tt::_member_object_pointer
                         | tt::_member_function_pointer
                         | tt::_enum
                         | tt::_union
                         | tt::_class
                         | tt::_function
                         >
{ };


// ----------------------------------------------------------------------
// is_member_pointer - metaprogramming composite type trait
// ----------------------------------------------------------------------

template< typename T >
  struct is_member_pointer
  : public tt::primary< T, tt::_member_object_pointer
                         | tt::_member_function_pointer
                         >
{ };


// ----------------------------------------------------------------------
// cv<> - helper analyzing a type's cv-qualification(s)
// ----------------------------------------------------------------------

namespace tt  {

template< typename T >
  struct cv
{
  static   bool const  is_c  = false;
  static   bool const  is_v  = false;
  typedef  T const             add_c_type;
  typedef  T       volatile    add_v_type;
  typedef  T const volatile    add_cv_type;
  typedef  T                   rem_c_type;
  typedef  T                   rem_v_type;
  typedef  T                   rem_cv_type;
};

template< typename T >
  struct cv<T const>
{
  static   bool const  is_c  = true;
  static   bool const  is_v  = false;
  typedef  T const             add_c_type;
  typedef  T const volatile    add_v_type;
  typedef  T const volatile    add_cv_type;
  typedef  T                   rem_c_type;
  typedef  T const             rem_v_type;
  typedef  T                   rem_cv_type;
};

template< typename T >
  struct cv<T volatile>
{
  static   bool const  is_c  = false;
  static   bool const  is_v  = true;
  typedef  T const volatile    add_c_type;
  typedef  T       volatile    add_v_type;
  typedef  T const volatile    add_cv_type;
  typedef  T       volatile    rem_c_type;
  typedef  T                   rem_v_type;
  typedef  T                   rem_cv_type;
};

template< typename T >
  struct cv<T const volatile>
{
  static   bool const  is_c  = true;
  static   bool const  is_v  = true;
  typedef  T const volatile    add_c_type;
  typedef  T const volatile    add_v_type;
  typedef  T const volatile    add_cv_type;
  typedef  T       volatile    rem_c_type;
  typedef  T const             rem_v_type;
  typedef  T                   rem_cv_type;
};

template< typename T >
  struct cv<T &>
{
  static   bool const  is_c  = false;
  static   bool const  is_v  = false;
  typedef  T &                 add_c_type;
  typedef  T &                 add_v_type;
  typedef  T &                 add_cv_type;
  typedef  T &                 rem_c_type;
  typedef  T &                 rem_v_type;
  typedef  T &                 rem_cv_type;
};

}  // namespace tt


// ----------------------------------------------------------------------
// is_const - metaprogramming type trait detecting type constness
// ----------------------------------------------------------------------

template< typename T >
  struct is_const
  : public integral_constant<bool, tt::cv<T>::is_c >  { };


// ----------------------------------------------------------------------
// is_volatile - metaprogramming type trait detecting type volatility
// ----------------------------------------------------------------------

template< typename T >
  struct is_volatile
  : public integral_constant<bool, tt::cv<T>::is_v >  { };


// ----------------------------------------------------------------------
// is_abstract_class - helper detecting when a class is abstract
// ----------------------------------------------------------------------

namespace tt  {

template< typename, bool >
  struct is_abstract_class
  : public false_type  { };  // default: not a class, hence not abstract

template< typename C >
  struct is_abstract_class<C,true>  // C is known to be a class type
{
protected:
  template< typename T >
    static  no_t   take( T (*)[1] );  // can't form array of abstract T
  template< typename T >
    static  yes_t  take( ...      );

public:
  static  bool const  value  =  sizeof( take<C>(0) ) == sizeof(yes_t);
};  // is_abstract_class<,true>

}  // namespace tt


// ----------------------------------------------------------------------
// is_abstract - metaprogramming type trait detecting abstract classes
// ----------------------------------------------------------------------

template< typename T >
  struct is_abstract
  : public tt::is_abstract_class< T
                                , is_class<T>::value
                                >
  { };


// ----------------------------------------------------------------------
// is_signed - metaprogramming type trait detecting type signedness
// ----------------------------------------------------------------------

template< typename >
  struct is_signed
  : public false_type  { };

template<>
  struct is_signed<signed char>
  : public true_type  { };
template<>
  struct is_signed<short>
  : public true_type  { };
template<>
  struct is_signed<int>
  : public true_type  { };
template<>
  struct is_signed<long>
  : public true_type  { };

template< typename T >
  struct is_signed<T const>
  : public is_signed<T>  { };
template< typename T >
  struct is_signed<T volatile>
  : public is_signed<T>  { };
template< typename T >
  struct is_signed<T const volatile>
  : public is_signed<T>  { };


// ----------------------------------------------------------------------
// is_unsigned - metaprogramming type trait detecting type unsignedness
// ----------------------------------------------------------------------

template< typename >
  struct is_unsigned
  : public false_type  { };

template<>
  struct is_unsigned<unsigned char>
  : public true_type  { };
template<>
  struct is_unsigned<unsigned short>
  : public true_type  { };
template<>
  struct is_unsigned<unsigned int>
  : public true_type  { };
template<>
  struct is_unsigned<unsigned long>
  : public true_type  { };

template< typename T >
  struct is_unsigned<T const>
  : public is_unsigned<T>  { };
template< typename T >
  struct is_unsigned<T volatile>
  : public is_unsigned<T>  { };
template< typename T >
  struct is_unsigned<T const volatile>
  : public is_unsigned<T>  { };


// ----------------------------------------------------------------------
// arr<> - helper analyzing a type's array qualification(s)
// ----------------------------------------------------------------------

namespace tt  {

template< typename T >
  struct arr  // non-array
{
  typedef  T  rem_ext_type;
  typedef  T  rem_arr_type;

  static  int const  rank =  0;

  template< unsigned I >
    struct extent  { static  int const  value  =  0; };
};

template< typename T, unsigned N >
  struct arr<T[N]>
{
  typedef  T                                  rem_ext_type;
  typedef  typename tt::arr<T>::rem_arr_type  rem_arr_type;

  static  int const  rank =  1 + tt::arr<T>::rank;

  template< unsigned I >
  struct extent
  {
    static  int const  value  =  (I == rank)
                              ?  N
                              :  tt::arr<T>::template extent<I>::value;
  };
};

template< typename T >
  struct arr<T[]>
{
  typedef  T  rem_ext_type;
  typedef  T  rem_arr_type;

  static  int const  rank =  1;

  template< unsigned I >
    struct extent  { static  int const  value  =  0; };
};

}  // namespace tt


// ----------------------------------------------------------------------
// rank - metaprogramming type trait detecting array's rank
// ----------------------------------------------------------------------

template< typename T >
  struct rank
  : public integral_constant<int, tt::arr<T>::rank>  { };


// ----------------------------------------------------------------------
// extent - metaprogramming type trait detecting array's extent
// ----------------------------------------------------------------------

template< typename T, unsigned I >
  struct extent
  : public integral_constant<int, tt::arr<T>::template extent<I>::value>  { };


// ----------------------------------------------------------------------
// is_same - metaprogramming type trait detecting type identity
// ----------------------------------------------------------------------

template< typename T, typename U >
  struct is_same      : public false_type  { };
template< typename T >
  struct is_same<T,T> : public true_type  { };


// ----------------------------------------------------------------------
// any_conversion - helper to avoid passing a UDT through ... parameter
// ----------------------------------------------------------------------

namespace tt  {

struct any_conversion
{
  template< typename T >
    any_conversion( T const volatile & );
  template< typename T >
    any_conversion( T                & );  // no cv-qual on fctn-refs
};  // any_conversion

}  // namespace tt


// ----------------------------------------------------------------------
// converts_to - helper detecting convertability
// ----------------------------------------------------------------------

namespace tt  {

template< typename From, typename To, bool >
  struct converts
  : public false_type  { };  // default: can't convert to abstract To

template< typename From, typename To >
struct converts<From,To,false>  // To is non-abstract
{
protected:
  static  yes_t  take( To, int );
  static  no_t   take( any_conversion, ... );
  static  From  from;

public:
  static  bool const  value
    =  sizeof( take( from, 0 ) ) == sizeof(yes_t);
};  // converts<>

}  // namespace tt


// ----------------------------------------------------------------------
// is_convertible - metaprogramming type trait detecting convertability
// ----------------------------------------------------------------------

template< typename From, typename To >
  struct is_convertible
  : public tt::converts<From,To,is_abstract<To>::value>  { };

template< >  struct is_convertible<void,void>
  : public true_type  { };

template< typename T >
  struct is_convertible<T,void>
  : public true_type  { };

template< typename T >
  struct is_convertible<void,T>
  : public false_type  { };

template< >
  struct is_convertible<const void,const void>
  : public true_type  { };

template< typename T >
  struct is_convertible<T,const void>
  : public true_type  { };

template< typename T >
  struct is_convertible<const void,T>
  : public false_type  { };

template< >
  struct is_convertible<volatile void,volatile void>
  : public true_type  { };

template< typename T >
  struct is_convertible<T,volatile void>
  : public true_type  { };

template< typename T >
  struct is_convertible<volatile void,T>
  : public false_type  { };

template< >
  struct is_convertible<const volatile void,const volatile void>
  : public true_type  { };

template< typename T >
  struct is_convertible<T,const volatile void>
  : public true_type  { };

template< typename T >
  struct is_convertible<const volatile void,T>
  : public false_type  { };

template< typename From, int N, typename To >
  struct is_convertible<From[N],To>
  : public is_convertible<From*,To>  { };

template< typename From, typename To, int N >
  struct is_convertible<From,To[N]>
  : public false_type  { };


// ----------------------------------------------------------------------
// remove_const - metaprogramming type trait ensuring non-constness
// ----------------------------------------------------------------------

template< typename T >
  struct remove_const
{
  typedef  typename tt::cv<T>::rem_c_type  type;
};


// ----------------------------------------------------------------------
// remove_volatile - metaprogramming type trait ensuring non-volatility
// ----------------------------------------------------------------------

template< typename T >
  struct remove_volatile
{
  typedef  typename tt::cv<T>::rem_v_type  type;
};


// ----------------------------------------------------------------------
// remove_cv - metaprogramming type trait ensuring no cv-qualification
// ----------------------------------------------------------------------

template< typename T >
  struct remove_cv
{
  typedef  typename tt::cv<T>::rem_cv_type  type;
};


// ----------------------------------------------------------------------
// add_const - metaprogramming type trait ensuring constness
// ----------------------------------------------------------------------

template< typename T >
  struct add_const
{
  typedef  typename tt::cv<T>::add_c_type  type;
};


// ----------------------------------------------------------------------
// add_volatile - metaprogramming type trait ensuring volatility
// ----------------------------------------------------------------------

template< typename T >
  struct add_volatile
{
  typedef  typename tt::cv<T>::add_v_type  type;
};


// ----------------------------------------------------------------------
// add_cv - metaprogramming type trait ensuring constness & volatility
// ----------------------------------------------------------------------

template< typename T >
  struct add_cv
{
  typedef  typename tt::cv<T>::add_cv_type  type;
};


// ----------------------------------------------------------------------
// ref<> - helper analyzing a type's reference qualification
// ----------------------------------------------------------------------

namespace tt  {

template< typename T
        , primary_code = encode<T>::value
        >
  struct ref  // non-lref && non-rref && non-void
{
  typedef  T&   add_lref_type;
  #if defined(CLHEP_HAS_RVALUE_REFERENCE)
  typedef  T&&  add_rref_type;
  #endif  // CLHEP_HAS_RVALUE_REFERENCE
  typedef  T    rem_ref_type;
};

template< typename T >
  struct ref<T&,_lvalue_reference>
{
  typedef  T&  add_lref_type;
  typedef  T&  add_rref_type;
  typedef  T   rem_ref_type;
};

#if defined(CLHEP_HAS_RVALUE_REFERENCE)
template< typename T >
  struct ref<T&&,_rvalue_reference>
{
  typedef  T&   add_lref_type;
  typedef  T&&  add_rref_type;
  typedef  T    rem_ref_type;
};
#endif  // CLHEP_HAS_RVALUE_REFERENCE

template< typename T >
  struct ref<T,_void>
{
  typedef  T  add_lref_type;
  typedef  T  add_rref_type;
  typedef  T  rem_ref_type;
};

}  // namespace tt


// ----------------------------------------------------------------------
// remove_reference - metaprogramming type trait ensuring non-reference
// ----------------------------------------------------------------------

template< typename T >
  struct remove_reference
{
  typedef  typename tt::ref<T>::rem_ref_type  type;
};


// ----------------------------------------------------------------------
// add_lvalue_reference - metaprogramming type trait ensuring lvalue-ref
// ----------------------------------------------------------------------

template< typename T >
  struct add_lvalue_reference
{
  typedef  typename tt::ref<T>::add_lref_type  type;
};


// ----------------------------------------------------------------------
// add_rvalue_reference - metaprogramming type trait ensuring rvalue-ref
// ----------------------------------------------------------------------

template< typename T >
  struct add_rvalue_reference
{
  typedef  typename tt::ref<T>::add_rref_type  type;
};


// ----------------------------------------------------------------------
// ptr<> - helper analyzing a type's pointer qualification
// ----------------------------------------------------------------------

namespace tt  {

template< typename T >
  struct ptr
{
  typedef  typename tt::ref<T>::rem_ref_type *  add_ptr_type;
  typedef  T                                    rem_ptr_type;
};

template< typename T >
  struct ptr<T *>
{
  typedef  T * *  add_ptr_type;
  typedef  T      rem_ptr_type;
};

template< typename T >
  struct ptr<T * const>
{
  typedef  T * const *  add_ptr_type;
  typedef  T            rem_ptr_type;
};

template< typename T >
  struct ptr<T * volatile>
{
  typedef  T * volatile *  add_ptr_type;
  typedef  T               rem_ptr_type;
};

template< typename T >
  struct ptr<T * const volatile>
{
  typedef  T * const volatile *  add_ptr_type;
  typedef  T                     rem_ptr_type;
};

}  // namespace tt


// ----------------------------------------------------------------------
// remove_extent - metaprogramming type trait reducing an array's extent
// ----------------------------------------------------------------------

template< typename T >
  struct remove_extent
{
  typedef  typename tt::arr<T>::rem_ext_type  type;
};


// ----------------------------------------------------------------------
// remove_all_extents - metaprogramming type trait yielding a non-array
// ----------------------------------------------------------------------

template< typename T >
  struct remove_all_extents
{
  typedef  typename tt::arr<T>::rem_arr_type  type;
};


// ----------------------------------------------------------------------
// remove_pointer - metaprogramming type trait ensuring non-pointer
// ----------------------------------------------------------------------

template< typename T >
  struct remove_pointer
{
  typedef  typename tt::ptr<T>::rem_ptr_type  type;
};


// ----------------------------------------------------------------------
// add_pointer - metaprogramming type trait ensuring pointer
// ----------------------------------------------------------------------

template< typename T >
  struct add_pointer
{
  typedef  typename tt::ptr<T>::add_ptr_type  type;
};


// ----------------------------------------------------------------------
// enable_if - metaprogramming construct for applied SFINAE
// ----------------------------------------------------------------------

template< typename T > struct enable_if<true ,T>  { typedef T type; };
template< typename T > struct enable_if<false,T>  { };


// ----------------------------------------------------------------------
// conditional - metaprogramming construct for type selection
// ----------------------------------------------------------------------

template< typename T, typename F > struct conditional<true ,T,F>  { typedef T type; };
template< typename T, typename F > struct conditional<false,T,F>  { typedef F type; };


// ----------------------------------------------------------------------
// is_ptr_convertible - variant of is_convertible, based on ptrs to types
// ----------------------------------------------------------------------

template< typename From, typename To >
  struct is_ptr_convertible
{
protected:
  static  tt::yes_t  take( To * );
  static  tt::no_t   take( ...  );

public:
  static  bool const  value
    =  sizeof( take( static_cast<From*>(0) ) ) == sizeof(tt::yes_t);
};  // is_ptr_convertible<,>


// ----------------------------------------------------------------------
// enable_if_convertible - convenience metaprogramming type trait
// ----------------------------------------------------------------------

template< typename From  // type of conversion's source
        , typename To    // type of conversion's target
        , typename R     // result type if conversion is valid
        >
  struct enable_if_convertible
  : public enable_if< is_convertible<From,To>::value, R >  { };


// ----------------------------------------------------------------------
// enable_if_ptr_convertible - convenience metaprogramming type trait
// ----------------------------------------------------------------------

template< typename From  // type of conversion's source
        , typename To    // type of conversion's target
        , typename R     // result type if conversion is valid
        >
  struct enable_if_ptr_convertible
  : public enable_if< is_ptr_convertible<From,To>::value, R >  { };


// ----------------------------------------------------------------------
// enable_if_auto_ptr - convenience metaprogramming type trait
// ----------------------------------------------------------------------

template< typename P  // pointee type
        , typename R  // result type
        >
  struct enable_if_auto_ptr  { };

template< typename P, typename R >
  struct enable_if_auto_ptr< std::auto_ptr<P>, R >
{
  typedef  R  type;
};


// ----------------------------------------------------------------------

}  // namespace CLHEP


#endif  // HEP_TYPE_TRAITS_H

// ======================================================================

/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef MCGIDI_VECTOR_HPP
#define MCGIDI_VECTOR_HPP

#define CPU_MEM false
#define UVM_MEM true

#ifdef HAVE_OPENMP_TARGET
    #ifdef USE_OPENMP_NO_GPU
        #define VAR_MEM false
    #else
        #define VAR_MEM true
    #endif
#else
    #define VAR_MEM false
#endif

typedef int MCGIDI_VectorSizeType;

#define MCGIDI_SWAP(a,b,type) {type ttttttttt=a;a=b;b=ttttttttt;}

#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#endif

#if defined(__HIP__)
#include <hip/hip_version.h>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hip/hip_common.h>
#endif

#include <string.h>
#include <stdio.h>
#include "cassert"
#include <algorithm>
#include <LUPI_declareMacro.hpp>
#include <vector>

namespace MCGIDI {

template <class T>
class Vector 
{
 private:
   T* _data;
   std::size_t _capacity;
   std::size_t _size;
   bool _mem_type;

 public:
   typedef T* iterator;
   typedef T* const_iterator;

   LUPI_HOST_DEVICE Vector()        : _data(0), _capacity(0), _size(0), _mem_type(CPU_MEM) {};
   LUPI_HOST_DEVICE Vector( std::size_t s, bool mem_flag = CPU_MEM ) : _data(0), _capacity(s), _size(s), _mem_type(mem_flag)
   {
       
      if( s == 0 ){ _data = nullptr; return;}	
        switch ((int)_mem_type){
            case CPU_MEM:
                _data = new T [_capacity];
                break;
            case UVM_MEM:
            {
                void *ptr = nullptr;
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
                cudaMallocManaged(&ptr, _capacity*sizeof(T), cudaMemAttachGlobal);
#elif defined(__HIP__) && !defined(__HIP_DEVICE_COMPILE__)
                hipMallocManaged(&ptr, _capacity*sizeof(T), hipMemAttachGlobal);
#endif
                _data = new(ptr) T[_capacity]; 
                break;
            }
            default:
                _data = new T [_capacity];
                break;
        }
   }
   LUPI_HOST_DEVICE Vector( std::size_t s, const T& d, bool mem_flag = CPU_MEM ) : _data(0), _capacity(s), _size(s), _mem_type(mem_flag)
   { 
      if( s == 0 ){ _data = nullptr; return;}	
        switch ( (int) _mem_type){
            case CPU_MEM:
                _data = new T [_capacity];
                break;
            case UVM_MEM:
            {
                void *ptr = nullptr;
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
                cudaMallocManaged(&ptr, _capacity*sizeof(T), cudaMemAttachGlobal);
#elif defined(__HIP__) && !defined(__HIP_DEVICE_COMPILE__)
                hipMallocManaged(&ptr, _capacity*sizeof(T), hipMemAttachGlobal);
#endif
                _data = new(ptr) T[_capacity];
                break;
            }
            default:
                _data = new T [_capacity];
                break;
        }
      for (std::size_t ii = 0; ii < _capacity; ++ii)
         _data[ii] = d;
   }

   LUPI_HOST_DEVICE Vector(const Vector<T>& aa )
        : _data(0), _capacity(aa._capacity), _size(aa._size), _mem_type(aa._mem_type)
   {
      if( _capacity == 0 ){ _data = nullptr; return; }

        switch ( (int) _mem_type){
            case CPU_MEM:
                _data = new T [_capacity];
                break;
            case UVM_MEM:
            {
                void *ptr = nullptr;
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
                cudaMallocManaged(&ptr, _capacity*sizeof(T), cudaMemAttachGlobal);
#elif defined(__HIP__) && !defined(__HIP_DEVICE_COMPILE__)
                hipMallocManaged(&ptr, _capacity*sizeof(T), hipMemAttachGlobal);
#endif
                _data = new(ptr) T[_capacity];
                break;
            }
            default:
                _data = new T [_capacity];
                break;
        }
 
      for (std::size_t ii=0; ii<_size; ++ii)
         _data[ii] = aa._data[ii];
   }

   LUPI_HOST Vector(const std::vector<T>& aa )
        : _data(0), _capacity(aa.size()), _size(aa.size()), _mem_type(CPU_MEM)
   {
      if( _capacity == 0 ){ _data = nullptr; return;}	

        switch ( (int) _mem_type){
            case CPU_MEM:
                _data = new T [_capacity];
                break;
            case UVM_MEM:
            {
                void *ptr = nullptr;
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
                cudaMallocManaged(&ptr, _capacity*sizeof(T), cudaMemAttachGlobal);
#elif defined(__HIP__) && !defined(__HIP_DEVICE_COMPILE__)
                hipMallocManaged(&ptr, _capacity*sizeof(T), hipMemAttachGlobal);
#endif
                _data = new(ptr) T[_capacity];
                break;
            }
            default:
                _data = new T [_capacity];
                break;
        }
 
      for (std::size_t ii=0; ii<_size; ++ii)
         _data[ii] = aa[ii];
   }
   
   LUPI_HOST_DEVICE ~Vector() { 
        switch ( (int) _mem_type){
            case CPU_MEM:
                delete[] _data; 
                break;
            case UVM_MEM:
                 for (std::size_t i=0; i < _size; ++i)
                   _data[i].~T();
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
                cudaFree(_data);
#elif defined(__HIP__) && !defined(__HIP_DEVICE_COMPILE__)
                hipFree(_data);
#endif
                break;
            default:
                delete[] _data; 
                break;
        }
   }

   LUPI_HOST_DEVICE iterator begin() { return _data; }

   LUPI_HOST_DEVICE const_iterator begin() const { return _data; }

   LUPI_HOST_DEVICE iterator end() { return _data + _size; }

   LUPI_HOST_DEVICE const_iterator end() const { return _data + _size; }

   /// Needed for copy-swap idiom
   LUPI_HOST_DEVICE void swap(Vector<T>& other)
   {
      MCGIDI_SWAP(_data,     other._data,     T*);
      MCGIDI_SWAP(_capacity, other._capacity, std::size_t);
      MCGIDI_SWAP(_size,     other._size,     std::size_t);
      MCGIDI_SWAP(_mem_type, other._mem_type, bool);
   }
   
   /// Implement assignment using copy-swap idiom
   LUPI_HOST_DEVICE Vector<T>& operator=(const Vector<T>& aa)
   {
      if (&aa != this)
      {
         Vector<T> temp(aa);
         this->swap(temp);
      }
      return *this;
   }

   LUPI_HOST Vector<T>& operator=(const std::vector<T>& aa)
   {
      Vector<T> temp(aa);
      this->swap(temp);
      return *this;
   }
   
   LUPI_HOST_DEVICE int get_mem_type()
   {
	return _mem_type;
   }

   LUPI_HOST_DEVICE void push_back( const T& dataElem )
   {
      assert( _size < _capacity );
      _data[_size] = dataElem;
      _size++;
   }

   LUPI_HOST_DEVICE const T& operator[]( std::size_t index ) const
   {
      // assert( index < _capacity ); 
      // assert( index >= 0); comment out pointless assertion size_t type is >= 0 by definition
      return _data[index];
   }

   LUPI_HOST_DEVICE T& operator[]( std::size_t index )
   {
      // assert( index < _capacity );
      // assert( index >= 0); comment out pointless assertion size_t type is >= 0 by definition
      return _data[index];
   }
   
   LUPI_HOST_DEVICE std::size_t capacity() const
   {
      return _capacity;
   }

   LUPI_HOST_DEVICE std::size_t size() const
   {
      return _size;
   }
   
   LUPI_HOST_DEVICE T& back()
   {
      return _data[_size-1];
   }
   
   LUPI_HOST_DEVICE T& back() const
   {
      return _data[_size-1];
   }
   
   LUPI_HOST_DEVICE void reserve( std::size_t s, char ** address = nullptr, bool mem_flag = CPU_MEM )
   {
      if (s == _capacity) return;
      assert( _capacity == 0 );
      _capacity = s;
      _mem_type = mem_flag;
      if( s == 0 ){ _data = nullptr; return;}	
        switch ( (int) _mem_type){
            case CPU_MEM:
                if (address == nullptr || *address == nullptr) _data = new T [_capacity];
                else {
                    _data = new(*address) T [_capacity];
                    *address += sizeof(T) * _capacity;
                }
                break;
            case UVM_MEM:
            {
                void *ptr = nullptr;
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
                cudaMallocManaged(&ptr, _capacity*sizeof(T), cudaMemAttachGlobal);
#elif defined(__HIP__) && !defined(__HIP_DEVICE_COMPILE__)
                hipMallocManaged(&ptr, _capacity*sizeof(T), hipMemAttachGlobal);
#endif
                _data = new(ptr) T[_capacity];
                break;
            }
            default:
                if (address == nullptr || *address == nullptr) _data = new T [_capacity];
                else {
                    _data = new(*address) T [_capacity];
                    *address += sizeof(T) * _capacity;
                }
                break;
        }
   }

   LUPI_HOST_DEVICE void resize( std::size_t s, char ** address = nullptr, bool mem_flag = CPU_MEM )
   {
      if (_capacity != 0) { 
          assert( _capacity >= s);
          _size = s;
          return;
      }
      assert( _capacity == 0 );
      _capacity = s;
      _size = s;
      _mem_type = mem_flag;
      if( s == 0 ){ _data = nullptr; return;}	
        switch ( (int) _mem_type){
            case CPU_MEM:
                if (address == nullptr || *address == nullptr) {
                    _data = new T [_capacity];
                }
                else {
                    _data = new(*address) T [_capacity];
                    std::size_t delta = sizeof(T) * _capacity;
                    std::size_t sub = delta % 8;
                    if (sub != 0) delta += (8-sub);
                    *address += delta;
                }
                break;
            case UVM_MEM:
            {
                void *ptr = nullptr;
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
                cudaMallocManaged(&ptr, _capacity*sizeof(T), cudaMemAttachGlobal);
#elif defined(__HIP__) && !defined(__HIP_DEVICE_COMPILE__)
                hipMallocManaged(&ptr, _capacity*sizeof(T), hipMemAttachGlobal);
#endif
                _data = new(ptr) T[_capacity];
                break;
            }
            default:
                if (address == nullptr || *address == nullptr) _data = new T [_capacity];
                else {
                    _data = new(*address) T [_capacity];
                    std::size_t delta = sizeof(T) * _capacity;
                    std::size_t sub = delta % 8;
                    if (sub != 0) delta += (8-sub);
                    *address += delta;
                }
                break;
        }
   }

   LUPI_HOST_DEVICE void resize( std::size_t s, const T& d, char ** address = nullptr, bool mem_flag = CPU_MEM ) 
   { 
      assert( _capacity == 0 );
      _capacity = s;
      _size = s;
      _mem_type = mem_flag;
      if( s == 0 ){ _data = nullptr; return;}	
        switch ( (int) _mem_type){
            case CPU_MEM:
                if (address == nullptr || *address == nullptr) _data = new T [_capacity];
                else {
                    _data = new(*address) T [_capacity];
                    std::size_t delta = sizeof(T) * _capacity;
                    std::size_t sub = delta % 8;
                    if (sub != 0) delta += (8-sub);
                    *address += delta;
                }
                break;
            case UVM_MEM:
            {
                void *ptr = nullptr;
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)
                cudaMallocManaged(&ptr, _capacity*sizeof(T), cudaMemAttachGlobal);
#elif defined(__HIP__) && !defined(__HIP_DEVICE_COMPILE__)
                hipMallocManaged(&ptr, _capacity*sizeof(T), hipMemAttachGlobal);
#endif
                _data = new(ptr) T[_capacity];
                break;
            }
            default:
                if (address == nullptr || *address == nullptr) _data = new T [_capacity];
                else {
                    _data = new(*address) T [_capacity];
                    std::size_t delta = sizeof(T) * _capacity;
                    std::size_t sub = delta % 8;
                    if (sub != 0) delta += (8-sub);
                    *address += delta;
                    *address += sizeof(T) * _capacity;
                }
                break;
        }
      for (std::size_t ii = 0; ii < _capacity; ++ii)
         _data[ii] = d;
   }

   LUPI_HOST_DEVICE bool empty() const
   {
       return ( _size == 0 );
   }

   LUPI_HOST_DEVICE void eraseEnd( std::size_t NewEnd )
   {
       assert( NewEnd <= _size );
       _size = NewEnd;
   }

   LUPI_HOST_DEVICE  void pop_back()
   {
       assert(_size > 0);
       _size--;
   }

   LUPI_HOST_DEVICE void clear()
   {
       _size = 0;
   }

   LUPI_HOST_DEVICE void appendList( std::size_t listSize, T* list )
   {
       assert( _size + listSize < _capacity );

       for( std::size_t i = _size; i < _size + listSize; i++ )
       {
           _data[i] = list[ i-_size ];
       }

   }

   //Atomically retrieve an availible index then increment that index some amount
   LUPI_HOST_DEVICE std::size_t atomic_Index_Inc( std::size_t inc )
   {
       if (_size+inc > _capacity)
          {MCGIDI_PRINTF("inc too much (size %d, inc %d cap %d)\n", _size, inc, _capacity); abort(); }
       assert(_size+inc <= _capacity);
       std::size_t pos;

//       #include "mc_omp_atomic_capture.hh"
       {pos = _size; _size = _size + inc;}

       return pos;
   }

   // This will not work for a vector of base classes.
   LUPI_HOST_DEVICE std::size_t internalSize() const {
       std::size_t delta = sizeof(T) * _size;
       std::size_t sub = delta % 8;
       if (sub != 0) delta += (8-sub);
       return delta;
   }

   LUPI_HOST_DEVICE void forceCreate(std::size_t a_size, T* a_data) {
       _capacity = a_size;
       _size = a_size;
       _data = a_data;
   }
};

}
#endif

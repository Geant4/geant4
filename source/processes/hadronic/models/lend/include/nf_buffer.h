/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef nf_buffer_h_included
#define nf_buffer_h_included


#if defined __cplusplus

#include <iterator>


template<typename T>
class nf_Buffer {
  private:
    T *m_data;
    size_t m_length;

  public:

    using iterator = T*;
    using const_iterator = T const *;

    inline
    constexpr
    nf_Buffer() noexcept : m_data(nullptr), m_length(0) {}

    inline
    nf_Buffer(nf_Buffer const &c) :
      m_data(new T[c.m_length]),
      m_length(c.m_length)
    {
      for(size_t i = 0;i < m_length;++ i){
        m_data[i] = c.m_data[i];
      }
    }

    inline
    ~nf_Buffer() noexcept {
      deallocate();
    }

    inline
    constexpr
    size_t size() const noexcept { return m_length; }

    inline
    void clear(T value){
      for(size_t i = 0;i < m_length;++ i){
        m_data[i] = value;
      }
    }

    inline
    void allocate(size_t length){
      deallocate();
      m_length = length;
      m_data = new T[length];
    }

    inline
    void deallocate() noexcept {
      delete[] m_data;
      m_length = 0;
    }

    inline
    void resize(size_t length){
      allocate(length);
    }

    inline
    std::vector<T> vector() const {
      return std::vector<T>(cbegin(), cend());
    }

    inline
    T* data() noexcept {return m_data;}

    inline
    constexpr
    T const * data() const noexcept {return m_data;}

    template<typename I>
    inline
    T& operator[](I idx) noexcept {return m_data[idx];}

    template<typename I>
    inline
    constexpr
    T const & operator[](I idx) const noexcept {return m_data[idx];}


    inline
    iterator begin() noexcept { return m_data; }

    inline
    constexpr
    const_iterator begin() const noexcept { return m_data; }

    inline
    iterator end() noexcept { return m_data + m_length; }

    inline
    constexpr
    const_iterator end() const noexcept { return m_data + m_length; }

    inline
    constexpr
    const_iterator cbegin() const noexcept { return m_data; }

    inline
    constexpr
    const_iterator cend() const noexcept { return m_data + m_length; }
};




#endif

#endif          /* End of nf_buffer_h_included. */

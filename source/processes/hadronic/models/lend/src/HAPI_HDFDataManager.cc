/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "HAPI.hpp"
#include <vector>

#ifdef HAPI_USE_HDF5
namespace HAPI {

    // constructor
    HDFDataManager::HDFDataManager(std::string const &a_filename) :
        m_filename( a_filename ) {

#if defined (GIDIP_HAVE_COMPILER_FLOATING_POINT_EXCEPTIONS)
        LUPI_FPE_disable_and_clear( __FILE__, __LINE__ );   // disable sigfpe cores
#endif

        m_file_id = H5Fopen( a_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
        H5Eset_auto1( nullptr, nullptr );

        m_dataset_ints = H5Dopen2( m_file_id, "iData", H5P_DEFAULT );
        m_iDataPresent = m_dataset_ints != H5I_INVALID_HID;
        if( m_iDataPresent ) m_dataspace_ints = H5Dget_space( m_dataset_ints );

        m_dataset_doubles = H5Dopen2( m_file_id, "dData", H5P_DEFAULT );
        m_dDataPresent = m_dataset_doubles != H5I_INVALID_HID;
        if( m_dDataPresent ) m_dataspace_doubles = H5Dget_space( m_dataset_doubles );

#if defined (GIDIP_HAVE_COMPILER_FLOATING_POINT_EXCEPTIONS)
        // Re-enable floating point exception detection
        LUPI_FPE_test( __FILE__, __LINE__ );                // test sigfpe exception
        LUPI_FPE_enable( __FILE__, __LINE__ );              // reenable sigfpe cores
#endif

        m_stride[0] = 1;
        m_block[0] = 1;
    }

    HDFDataManager::~HDFDataManager()
    {
        if( m_iDataPresent ) {
            H5Dclose(m_dataset_ints);
            H5Sclose(m_dataspace_ints);
        }   
        if( m_iDataPresent ) {
            H5Dclose(m_dataset_doubles);
            H5Sclose(m_dataspace_doubles);
        }
      H5Fclose(m_file_id);
    }

    void HDFDataManager::getDoubles(nf_Buffer<double> &result, size_t startIndex, size_t endIndex)
    {
        if( !m_dDataPresent ) throw LUPI::Exception( "HDFDataManager::getDoubles: HDF5 file " + m_filename + " has no 'dData' dataset." );

        hid_t memspace;
        herr_t status;

        hsize_t size = endIndex - startIndex;

        hsize_t dims[] {size};
        hsize_t offset[] {startIndex};
        hsize_t count[] {size};

        result.resize(size);
        m_num_double_reads ++;
        m_num_double_elem += size;
        
        // now can we access the allocated array and read into that?

        memspace = H5Screate_simple(1, dims, nullptr);
        status = H5Sselect_hyperslab(m_dataspace_doubles, H5S_SELECT_SET, offset, m_stride, count, m_block);
        if( status != 0 ) throw "H5Sselect_hyperslab error in HDFDataManager::getDoubles.";

        status = H5Dread(m_dataset_doubles, H5T_NATIVE_DOUBLE, memspace, m_dataspace_doubles, H5P_DEFAULT, result.data());
        if( status != 0 ) throw "H5Dread error in HDFDataManager::getDoubles.";

        H5Sclose(memspace);

    }

    void HDFDataManager::getInts(nf_Buffer<int> &result, size_t startIndex, size_t endIndex)
    {
        if( !m_iDataPresent ) throw LUPI::Exception( "HDFDataManager::getInts: HDF5 file " + m_filename + " has no 'iData' dataset." );

        hid_t memspace;
        herr_t status;
        hsize_t size = endIndex - startIndex;

        hsize_t dims[] {size};
        hsize_t offset[] {startIndex};
        hsize_t count[] {size};

        result.resize(size);

        m_num_int_reads ++;
        m_num_int_elem += size;

        memspace = H5Screate_simple(1, dims, nullptr);
        status = H5Sselect_hyperslab(m_dataspace_ints, H5S_SELECT_SET, offset, m_stride, count, m_block);
        if( status != 0 ) throw "H5Sselect_hyperslab error in HDFDataManager::getDoubles.";

        status = H5Dread(m_dataset_ints, H5T_NATIVE_INT, memspace, m_dataspace_ints, H5P_DEFAULT, result.data());
        if( status != 0 ) throw "H5Dread error in HDFDataManager::getDoubles.";

        H5Sclose(memspace);

    }

}
#endif

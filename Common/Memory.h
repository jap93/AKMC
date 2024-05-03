#pragma once

const int OFFSET = 32;

#include <complex>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "Constants.h"

#include "mm_malloc.h"

/////////////////////////////////////////////////////////////////////////////////////////
// MATRIX functions start here
/////////////////////////////////////////////////////////////////////////////////////////

inline std::complex<double>* alloc_cvector(int n, std::string str)
{
    std::complex<double>* v = nullptr;

    v = (std::complex<double>*) malloc(n * sizeof(std::complex<double>));

    return v;
}

inline void free_cvector(std::complex<double>* v)
{

    if (v != nullptr)
    {

        free(v);

        v = nullptr;

    }

}

inline std::complex<double>* alloc_cvector_alligned(int n, std::string str)
{
    std::complex<double>* v = nullptr;

    v = (std::complex<double>*) _mm_malloc(n * sizeof(std::complex<double>), DOFFSET);

    return v;
}

inline void free_cvector_alligned(std::complex<double>* v)
{

    if (v != nullptr)
    {

        _mm_free(v);

        v = nullptr;

    }

}

/** make vector of double */
inline float* alloc_fvector(int n, std::string str, float iv = 0.0)
{
    int i;

    float* v = nullptr;

    try {

        #if defined NEW

            v = new float[n];

        #elif defined ALIGNINTEL

            v = (float*) _mm_malloc(n * sizeof(float), DOFFSET);

        #elif defined ALIGNPOSIX

            void *memptr;
            int error = posix_memalign(&memptr, DOFFSET, n * sizeof(float));
            v = (float*) memptr;

        #elif defined ALIGNMEM

            v = (float*)aligned_alloc(DOFFSET, n * sizeof(float));

        #else

            v = (float*)malloc(n * sizeof(float));

        #endif

    }
    catch (std::bad_alloc&)
    {
        std::cout << "could not allocate float vector" << std::endl;
        std::cout << "attempted to allocate " << str << std::endl;
        std::cout.flush();
        exit(EXIT_FAILURE);
    }

    // set initial value
    for (i = 0; i < n; i++)
        v[i] = iv;

    return v;

}

inline void free_fvector(float* v)
{

    if (v != nullptr)
    {
        #if defined NEW

            delete [] v;

        #elif defined ALIGNINTEL

            _mm_free(v);

        #elif defined ALIGNPOSIX

            free (v);

        #elif defined ALIGNMEM

            free(v);

        #else

            free(v);

        #endif

        v = nullptr;

    }

}

inline double* alloc_dvector(int n, std::string str, double iv = 0.0)
{
    int i;

    double* v = nullptr;

    try {

        #if defined NEW

            v = new double[n];

        #elif defined ALIGNINTEL

            v = (double*) _mm_malloc(n * sizeof(double), DOFFSET);

        #elif defined ALIGNPOSIX

            void *memptr;
            int error = posix_memalign(&memptr, DOFFSET, n * sizeof(double));
            v = (double*) memptr;

        #elif defined ALIGNMEM

            v = (double*)aligned_alloc(DOFFSET, n * sizeof(double));

        #else

            v = (double*)malloc(n * sizeof(double));

        #endif

    }
    catch (std::bad_alloc&)
    {
        std::cout << "could not allocate double3 vector" << std::endl;
        std::cout << "attempted to allocate " << str << std::endl;
        std::cout.flush();
        exit(EXIT_FAILURE);
    }

    // set initial value
    for (i = 0; i < n; i++)
        v[i] = iv;

    return v;

}

/** free memory of double vector */
inline void free_dvector(double* v)
{

    if (v != nullptr)
    {
        #if defined NEW

            delete [] v;

        #elif defined ALIGNINTEL

            _mm_free(v);

        #elif defined ALIGNPOSIX

            free (v);

        #elif defined ALIGNMEM

            free(v);

        #else

            free(v);

        #endif

        v = nullptr;

    }

}

/** make vector of double */
inline double* alloc_dvector_aligned(int n, std::string str, double iv = 0.0)
{
    int i;

    double* v = nullptr;

    try {

        v = (double*)_mm_malloc(n * sizeof(double), DOFFSET);

    }
    catch (std::bad_alloc&)
    {
        std::cout << "could not allocate double3 vector" << std::endl;
        std::cout << "attempted to allocate " << str << std::endl;
        std::cout.flush();
        exit(EXIT_FAILURE);
    }

    // set initial value
    for (i = 0; i < n; i++)
        v[i] = iv;

    return v;

}

/** free memory of double vector */
inline void free_dvector_aligned(double* v)
{

    if (v != nullptr)
    {

        _mm_free(v);

        v = nullptr;

    }

}


/** make an integer vector */
inline  int* alloc_ivector(int n, std::string str, int iv = 0)
{
    int* v = nullptr;

    try
    {

        v = (int*)malloc(n * sizeof(int));

    }
    catch (std::bad_alloc&)
    {
        std::cout << "could not allocate double3 vector" << std::endl;
        std::cout << "attempted to allocate " << str << std::endl;
        std::cout.flush();
        exit(EXIT_FAILURE);
    }


    for (int i = 0; i < n; i++)
        v[i] = iv;

    return v;

}

/** free memory for integer vector */
inline void free_ivector(int* v)
{

    if (v != nullptr)
    {

        free(v);

        v = nullptr;

    }

}

/** make an integer vector */
inline unsigned int* alloc_uivector(int n, std::string str, unsigned int iv = 0)
{
    unsigned int* v = nullptr;

    try
    {

        v = (unsigned int*)malloc(n * sizeof(int));

    }
    catch (std::bad_alloc&)
    {
        std::cout << "could not allocate double3 vector" << std::endl;
        std::cout << "attempted to allocate " << str << std::endl;
        std::cout.flush();
        exit(EXIT_FAILURE);
    }


    for (int i = 0; i < n; i++)
        v[i] = iv;

    return v;
}

/** free memory for integer vector */
inline void free_uivector(unsigned int* v)
{

    if (v != nullptr)
    {

        free(v);

        v = nullptr;

    }

}
/////////////////////////////////////////////////////////////////////////////////////////
// MATRIX functions start here
/////////////////////////////////////////////////////////////////////////////////////////

/** creates an array of pointers for matrix of double */
inline double** alloc_dmatrix(int nr, int nc, std::string str, double iv = 0.0)
{
    int
        i, j;

    double** m = nullptr;

    if (nr == 0 || nc == 0)
        return m;

    try
    {

        m = (double **)malloc(nr * sizeof(double *)); 
        for (i=0; i < nr; i++) 
            m[i] = (double *)malloc(nc * sizeof(double)); 

    }
    catch (std::bad_alloc&)
    {
        std::cout << "could not allocate double3 vector" << std::endl;
        std::cout << "attempted to allocate " << str << std::endl;
        std::cout.flush();
        exit(EXIT_FAILURE);
    }


    for (i = 0; i < nr; i++)
    {

        for (j = 0; j < nc; j++)
        {

            m[i][j] = iv;

        }

    }

    return m;

}


/** breaks down a matrix of type double */
inline void free_dmatrix(double** m)
{

    if (m != nullptr)
    {

        free(m);

        m = nullptr;

    }

}

/** makes an array of pointers for integer matrix */
inline int** alloc_imatrix(int nr, int nc, std::string str, int iv = 0)
{

    int
        i, j;

    int** m = nullptr;

    if (nr == 0 || nc == 0)
        return m;

    try
    {

        m = (int **)malloc(nr * sizeof(int *)); 
        for (i=0; i < nr; i++) 
            m[i] = (int *)malloc(nc * sizeof(int)); 

    }
    catch (std::bad_alloc&)
    {
        std::cout << "could not allocate int matrix" << std::endl;
        std::cout << "attempted to allocate " << str << std::endl;
        std::cout.flush();
        exit(EXIT_FAILURE);
    }


    for (i = 0; i < nr; i++)
    {

        for (j = 0; j < nc; j++)
        {

            m[i][j] = iv;

        }

    }

    return m;

}


/** breaks down integer matrix */
inline void free_imatrix(int** m)
{
    if (m != nullptr)
    {

        free(m);
        
        m = nullptr;

    }

}

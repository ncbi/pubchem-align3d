/*  $Id: shape_math_functions.hpp 685605 2024-07-26 12:29:33Z thiessen $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Authors:  Evan Bolton, Leonid Zaslavsky, Paul Thiessen
*
* ===========================================================================
*/

#ifndef ALIGN3D_SHAPE_MATH_FUNCTIONS__HPP
#define ALIGN3D_SHAPE_MATH_FUNCTIONS__HPP


namespace  Align3D {

  template< class T > T expApproxDouble( T x ){
    const double C = 0.0009765625; // 1.0/1024; 1024 = 2^10;
    double x_double = static_cast< double >(x);
    double y = static_cast< double>(1) + C*x_double;
    y = y*y; // 1
    y = y*y; // 2
    y = y*y; // 3
    y = y*y; // 4
    y = y*y; // 5 
    y = y*y; // 6
    y = y*y; // 7
    y = y*y; // 8
    y = y*y; // 9
    y = y*y; // 10    
    return static_cast< T >( y );
  }
  
  template< class T > T expLookup( T x ){
    x = -500.0 * x;
    if( x <= (T) 0 ){
      return (T) 1;
    }else{
      long lookup = std::round( (T) x );
      if( lookup < expMinusSize ){
	return (T) expMinus[ lookup ];
      }
    }
    return (T) 0;
  }

  template< class T > T expLookupTrunc( T x ){
    x = -500.0 * x;
    if( x <= (T) 0 ){
      return (T) 1;
    }else{
      long lookup = std::round( (T) x );
      if( lookup < expMinusSizeTrunc ){
	return (T) expMinusTrunc[ lookup ];
      }
    }
    return (T) 0;
  }

  template< class T > T exp( T x ){
    return expApproxDouble( x );
    //    return expLookup( x );
    //    return expLookupTrunc( x );
    //    return std::exp( x );
    //    return (T) std::exp( (float) x );
  }
}

#endif // ALIGN3D_SHAPE_MATH_FUNCTIONS__HPP
  

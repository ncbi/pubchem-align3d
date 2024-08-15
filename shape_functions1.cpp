/*  $Id: shape_functions1.cpp 685605 2024-07-26 12:29:33Z thiessen $
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

// C++ Includes

#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <iostream>
#include <algorithm>


#include "shape_functions.hpp"

namespace Align3D {

void VRotTransMatrixNew( Uint8& packed,
				   const float* m )
{
  packed = 0;

  double  trace = 1.0 + m[ 0 ] + m[ 4 ] + m[ 8 ];
  double  s, qw, qx, qy, qz; 
  if ( 0.1 < trace ) {  // Default
    s = 2.0 * sqrt( trace );
    qw = 0.25 * s;
    double  s_invert = 1.0 / s;
    qx = ( m[ 7 ] - m[ 5 ] ) * s_invert;
    qy = ( m[ 2 ] - m[ 6 ] ) * s_invert;
    qz = ( m[ 3 ] - m[ 1 ] ) * s_invert;
    ////// PUBCHEM_INFO( "0 trace: " << trace
    //              << "  (Q: " << qw << " " << qx << " " << qy << " " << qz << ")"
    //              << "  (M: " << m[0] << " " << m[4] << " " << m[8] << ")" );
  } else if ( 0.005 < ( m[ 0 ] - m[ 4 ] ) &&
	      0.005 < ( m[ 0 ] - m[ 8 ] ) ) {  // X is largest
    s = 2.0 * sqrt( 1.0 + m[ 0 ] - m[ 4 ] - m[ 8 ] );
    double  s_invert = 1.0 / s;
    qw = ( m[ 7 ] - m[ 5 ] ) * s_invert;
    qx = 0.25 * s;
    qy = ( m[ 3 ] + m[ 1 ] ) * s_invert;
    qz = ( m[ 2 ] + m[ 6 ] ) * s_invert;
    ////// PUBCHEM_INFO( "1 trace: " << trace
    //              << "  (Q: " << qw << " " << qx << " " << qy << " " << qz << ")"
    //              << "  (M: " << m[0] << " " << m[4] << " " << m[8] << ")" );
  } else if ( 0.005 < ( m[ 4 ] - m[ 8 ] ) ) {  // Y is largest
    s = 2.0 * sqrt( 1.0 + m[ 4 ] - m[ 0 ] - m[ 8 ] );
    double  s_invert = 1.0 / s;
    qw = ( m[ 2 ] - m[ 6 ] ) * s_invert;
    qx = ( m[ 3 ] + m[ 1 ] ) * s_invert;
    qy = 0.25 * s;
    qz = ( m[ 7 ] + m[ 5 ] ) * s_invert;
    ////// PUBCHEM_INFO( "2 trace: " << trace
    //              << "  (Q: " << qw << " " << qx << " " << qy << " " << qz << ")"
    //              << "  (M: " << m[0] << " " << m[4] << " " << m[8] << ")" );
  } else {    // Z is largest
    s = 2.0 * sqrt( 1.0 + m[ 8 ] - m[ 0 ] - m[ 4 ] );
    double  s_invert = 1.0 / s;
    qw = ( m[ 3 ] - m[ 1 ] ) * s_invert;
    qx = ( m[ 2 ] + m[ 6 ] ) * s_invert;
    qy = ( m[ 7 ] + m[ 5 ] ) * s_invert;
    qz = 0.25 * s;
    ////// PUBCHEM_INFO( "3 trace: " << trace
    //              << "  (Q: " << qw << " " << qx << " " << qy << " " << qz << ")"
    //              << "  (M: " << m[0] << " " << m[4] << " " << m[8] << ")" );
  }

  // Encode Quaternion
  static const unsigned  rbase  = 0xFF;  // Rotation base
  static const double  cr = static_cast< double >( rbase ) / 2.0;

  // Normalize
  s = cr / sqrt( qw * qw + qx * qx + qy * qy + qz * qz );

  Uint8  p1 = static_cast< Uint8 >( rint( cr + ( qw * s ) ) );
  Uint8  p2 = static_cast< Uint8 >( rint( cr + ( qx * s ) ) );
  Uint8  p3 = static_cast< Uint8 >( rint( cr + ( qy * s ) ) );
  Uint8  p4 = static_cast< Uint8 >( rint( cr + ( qz * s ) ) );
  if ( rbase < p1 ) {
    //// //// PUBCHEM_WARN( "Qw encoding is corrupt (0-255 range)!  " << p1 << " " << qw * s );
  }
  if ( rbase < p2 ) {
    //// //// PUBCHEM_WARN( "Qx encoding is corrupt (0-255 range)!  " << p2 << " " << qx * s );
  }
  if ( rbase < p3 ) {
    //// //// PUBCHEM_WARN( "Qy encoding is corrupt (0-255 range)!  " << p3 << " " << qy * s );
  }
  if ( rbase < p4 ) {
    //// //// PUBCHEM_WARN( "Qz encoding is corrupt (0-255 range)!  " << p4 << " " << qz * s );
  }

  // Encode Translation
  static const unsigned  tbase = static_cast< unsigned >( pow( static_cast< double >( 0x1FFFFFFF ), ( 1.0 / 3.0 ) ) );
  static const unsigned  tbase2 = tbase * tbase;
  static const double  itresolution = static_cast< double >( tbase ) / log( 101. );

  double  x = m[  9 ];  // X
  if ( x < 0.0 ) {
    packed |= 0x01;
    x = -x;
  }
  x = log( 1.0 + x );
  Uint8  p5 = static_cast< Uint8 >( rint( x * itresolution ) );
  if ( p5 >= tbase ) {  // Range check
    p5 = tbase - 1;
    //// //// PUBCHEM_WARN( "X axis translation abs(" << m[  9 ] << ") is beyond maximum range of " << PubChem::Data::expPackedNew2[ p5 ] );
  }

  double  y = m[ 10 ];  // Y
  if ( y < 0.0 ) {
    packed |= 0x02;
    y = -y;
  }
  y = log( 1.0 + y );
  Uint8  p6 = static_cast< Uint8 >( rint( y * itresolution ) );
  if ( p6 >= tbase ) {  // Range check
    p6 = tbase - 1;
    //// //// PUBCHEM_WARN( "Y axis translation abs(" << m[ 10 ] << ") is beyond maximum range of " << PubChem::Data::expPackedNew2[ p6 ] );
  }

  double  z = m[ 11 ];  // Z
  if ( z < 0.0 ) {
    packed |= 0x04;
    z = -z;
  }
  z = log( 1.0 + z );
  Uint8  p7 = static_cast< Uint8 >( rint( z * itresolution ) );
  if ( p7 >= tbase ) {  // Range check
    p7 = tbase - 1;
    //// //// PUBCHEM_WARN( "Z axis translation abs(" << m[ 11 ] << ") is beyond maximum range of " << PubChem::Data::expPackedNew2[ p7 ] );
  }

  Uint8  pt = p5 * tbase2 + p6 * tbase + p7;

  // Pack matrix
  packed |= ( ( p1 & rbase ) << 56 ); // Qw
  packed |= ( ( p2 & rbase ) << 48 ); // Qy
  packed |= ( ( p3 & rbase ) << 40 ); // Qx
  packed |= ( ( p4 & rbase ) << 32 ); // Qz
  packed |= ( ( pt & 0x1FFFFFFF ) << 3 );  // X,Y,Z
}



void VQuaternionTransArray2RotTransMatrix( double* m,
						     const double* q,
						     const bool normalize )
{
  double  q00 = q[ 0 ];
  double  q01 = q[ 1 ];
  double  q02 = q[ 2 ];
  double  q03 = q[ 3 ];

  // Normalize quaternion
  if ( true == normalize ) {
    double  s = 1.0 / sqrt( q00 * q00 + q01 * q01 + q02 * q02 + q03 * q03 );
    q00 *= s;
    q01 *= s;
    q02 *= s;
    q03 *= s;
  }

  double  sqw = q00 * q00;
  double  sqx = q01 * q01;
  double  sqy = q02 * q02;
  double  sqz = q03 * q03;

  double  qxqy = q01 * q02;
  double  qwqz = q00 * q03;

  double  qxqz = q01 * q03;
  double  qwqy = q00 * q02;

  double  qyqz = q02 * q03;
  double  qwqx = q00 * q01;

  // Left-handed rotation matrix
  m[  0 ] = sqw + sqx - sqy - sqz;  // X X
  m[  1 ] = 2.0 * ( qxqy - qwqz );  // X Y
  m[  2 ] = 2.0 * ( qxqz + qwqy );  // X Z
  m[  3 ] = 2.0 * ( qxqy + qwqz );  // Y X
  m[  4 ] = sqw - sqx + sqy - sqz;  // Y Y
  m[  5 ] = 2.0 * ( qyqz - qwqx );  // Y Z
  m[  6 ] = 2.0 * ( qxqz - qwqy );  // Z X
  m[  7 ] = 2.0 * ( qyqz + qwqx );  // Z Y
  m[  8 ] = sqw - sqx - sqy + sqz;  // Z Z
  m[  9 ] = q[ 4 ];
  m[ 10 ] = q[ 5 ];
  m[ 11 ] = q[ 6 ];
}



void VQuaternionTransArray2RotTransMatrix( float* m,
						     const double* q,
						     const bool normalize )
{
  double  q00 = q[ 0 ];
  double  q01 = q[ 1 ];
  double  q02 = q[ 2 ];
  double  q03 = q[ 3 ];

  // Normalize quaternion
  if ( true == normalize ) {
    double  s = 1.0 / sqrt( q00 * q00 + q01 * q01 + q02 * q02 + q03 * q03 );
    q00 *= s;
    q01 *= s;
    q02 *= s;
    q03 *= s;
  }

  double  sqw = q00 * q00;
  double  sqx = q01 * q01;
  double  sqy = q02 * q02;
  double  sqz = q03 * q03;

  double  qxqy = q01 * q02;
  double  qwqz = q00 * q03;

  double  qxqz = q01 * q03;
  double  qwqy = q00 * q02;

  double  qyqz = q02 * q03;
  double  qwqx = q00 * q01;

  // Left-handed rotation matrix
  m[  0 ] = sqw + sqx - sqy - sqz;  // X X
  m[  1 ] = 2.0 * ( qxqy - qwqz );  // X Y
  m[  2 ] = 2.0 * ( qxqz + qwqy );  // X Z
  m[  3 ] = 2.0 * ( qxqy + qwqz );  // Y X
  m[  4 ] = sqw - sqx + sqy - sqz;  // Y Y
  m[  5 ] = 2.0 * ( qyqz - qwqx );  // Y Z
  m[  6 ] = 2.0 * ( qxqz - qwqy );  // Z X
  m[  7 ] = 2.0 * ( qyqz + qwqx );  // Z Y
  m[  8 ] = sqw - sqx - sqy + sqz;  // Z Z
  m[  9 ] = q[ 4 ];
  m[ 10 ] = q[ 5 ];
  m[ 11 ] = q[ 6 ];
}




void  VCombineRotTransMatrixMatrix( double* t,
					      const double* t1,
					      const double* t2 )
{
  // Rotation Matrix
  t[ 0 ] = t1[ 0] * t2[ 0] + t1[ 1] * t2[ 3] + t1[ 2] * t2[ 6];
  t[ 1 ] = t1[ 0] * t2[ 1] + t1[ 1] * t2[ 4] + t1[ 2] * t2[ 7];
  t[ 2 ] = t1[ 0] * t2[ 2] + t1[ 1] * t2[ 5] + t1[ 2] * t2[ 8];

  t[ 3 ] = t1[ 3] * t2[ 0] + t1[ 4] * t2[ 3] + t1[ 5] * t2[ 6];
  t[ 4 ] = t1[ 3] * t2[ 1] + t1[ 4] * t2[ 4] + t1[ 5] * t2[ 7];
  t[ 5 ] = t1[ 3] * t2[ 2] + t1[ 4] * t2[ 5] + t1[ 5] * t2[ 8];

  t[ 6 ] = t1[ 6] * t2[ 0] + t1[ 7] * t2[ 3] + t1[ 8] * t2[ 6];
  t[ 7 ] = t1[ 6] * t2[ 1] + t1[ 7] * t2[ 4] + t1[ 8] * t2[ 7];
  t[ 8 ] = t1[ 6] * t2[ 2] + t1[ 7] * t2[ 5] + t1[ 8] * t2[ 8];

  // Translation Matrix
  t[  9 ] = t1[ 0] * t2[ 9] + t1[ 1] * t2[10] + t1[ 2] * t2[11] + t1[ 9];
  t[ 10 ] = t1[ 3] * t2[ 9] + t1[ 4] * t2[10] + t1[ 5] * t2[11] + t1[10];
  t[ 11 ] = t1[ 6] * t2[ 9] + t1[ 7] * t2[10] + t1[ 8] * t2[11] + t1[11];
}


void VRotTransMatrix2QuaternionTransArray( double* q,
						     const double* m )
{
  double  trace = 1.0 + m[ 0 ] + m[ 4 ] + m[ 8 ];
  double  s, qw, qx, qy, qz;
  if ( trace > 0.00001 ) {
    s = 2.0 * sqrt( trace );
    qw = 0.25 * s;
    double  s_invert = 1.0 / s;
    qx = ( m[ 7 ] - m[ 5 ] ) * s_invert;
    qy = ( m[ 2 ] - m[ 6 ] ) * s_invert;
    qz = ( m[ 3 ] - m[ 1 ] ) * s_invert;
    ////// PUBCHEM_INFO( "0 trace: " << trace << "  " << qw << " " << qx << " " << qy << " " << qz );
  } else if ( m[ 0 ] > m[ 4 ] &&
	      m[ 0 ] > m[ 8 ] ) {
    s = 2.0 * sqrt( 1.0 + m[ 0 ] - m[ 4 ] - m[ 8 ] );
    double  s_invert = 1.0 / s;
    qw = ( m[ 7 ] - m[ 5 ] ) * s_invert;
    qx = 0.25 * s;
    qy = ( m[ 3 ] + m[ 1 ] ) * s_invert;
    qz = ( m[ 2 ] + m[ 6 ] ) * s_invert;
    ////// PUBCHEM_INFO( "1 trace: " << trace << "  " << qw << " " << qx << " " << qy << " " << qz );
  } else if ( m[ 4 ] > m[ 8 ] ) {
    s = 2.0 * sqrt( 1.0 + m[ 4 ] - m[ 0 ] - m[ 8 ] );
    double  s_invert = 1.0 / s;
    qw = ( m[ 2 ] - m[ 6 ] ) * s_invert;
    qx = ( m[ 3 ] + m[ 1 ] ) * s_invert;
    qy = 0.25 * s;
    qz = ( m[ 7 ] + m[ 5 ] ) * s_invert;
    ////// PUBCHEM_INFO( "2 trace: " << trace << "  " << qw << " " << qx << " " << qy << " " << qz );
  } else {
    s = 2.0 * sqrt( 1.0 + m[ 8 ] - m[ 0 ] - m[ 4 ] );
    double  s_invert = 1.0 / s;
    qw = ( m[ 3 ] - m[ 1 ] ) * s_invert;
    qx = ( m[ 2 ] + m[ 6 ] ) * s_invert;
    qy = ( m[ 7 ] + m[ 5 ] ) * s_invert;
    qz = 0.25 * s;
    ////// PUBCHEM_INFO( "3 trace: " << trace << "  " << qw << " " << qx << " " << qy << " " << qz );
  }

  // Normalize
  s = 1.0 / sqrt( qw * qw + qx * qx + qy * qy + qz * qz );
  q[ 0 ] = qw * s;
  q[ 1 ] = qx * s;
  q[ 2 ] = qy * s;
  q[ 3 ] = qz * s;

  q[ 4 ] = m[  9 ];
  q[ 5 ] = m[ 10 ];
  q[ 6 ] = m[ 11 ];
}










void VStericFrameRotation( float* m,     // float[ 9 ]
				     const float* coord, // float[ n * 3 ]
				     const unsigned n,
				     float* ev )   // float[ 3 ]
{
  /*
    std::cerr << "*** COORDINATES: " << std::endl;

    for ( unsigned  i = 0;  i < n;  ++i ) {
      std::cerr << "\t" << coord[3*i]
		<< "\t" << coord[3*i+1]
		<< "\t" << coord[3*i+2]
		<< std::endl;
    }
  */
  
    // Compute initial inertial tensor
    double  it[ 3 ][ 3 ];
    double  sx = coord[ 0 ], sy = coord[ 1 ], sz = coord[ 2 ];
    double  xx = sx * sx, yy = sy * sy, zz = sz * sz;
    it[ 0 ][ 0 ] = yy + zz;
    it[ 1 ][ 1 ] = xx + zz;
    it[ 2 ][ 2 ] = xx + yy;
    it[ 0 ][ 1 ] = sx * sy;
    it[ 0 ][ 2 ] = sx * sz;
    it[ 1 ][ 2 ] = sy * sz;
    for ( unsigned  i = 1;  i < n;  ++i ) {
        sx = coord[     i * 3 ];
        sy = coord[ 1 + i * 3 ];
        sz = coord[ 2 + i * 3 ];

        xx = sx * sx;
        yy = sy * sy;
        zz = sz * sz;

        it[ 0 ][ 0 ] += yy + zz;
        it[ 1 ][ 1 ] += xx + zz;
        it[ 2 ][ 2 ] += xx + yy;
        it[ 0 ][ 1 ] += sx * sy;
        it[ 0 ][ 2 ] += sx * sz;
        it[ 1 ][ 2 ] += sy * sz;
    }

    /*
    std::cerr << "*** Begin MATRIX: " << std::endl;
    for ( unsigned  i = 0;  i < 3;  ++i ) {
          for ( unsigned  j = 0;  j < 3;  ++j ) {
	    std::cerr << "\t" << it[i][j];
	  }
	  std::cerr << std::endl;
	  }
    */
    

    // Diagonalize 3x3 steric inertial tensor matrix using Jacobi Diagonalization method

    // Initialize eigenvectors
    m[  0 ] = 1.0;
    m[  1 ] = 0.0;
    m[  2 ] = 0.0;
    m[  3 ] = 0.0;
    m[  4 ] = 1.0;
    m[  5 ] = 0.0;
    m[  6 ] = 0.0;
    m[  7 ] = 0.0;
    m[  8 ] = 1.0;
    m[  9 ] = 0.0;
    m[ 10 ] = 0.0;
    m[ 11 ] = 0.0;
  
    // Initialize eigenvalues
    double  e[ 3 ];
    e[ 0 ] = it[ 0 ][ 0 ];
    e[ 1 ] = it[ 1 ][ 1 ];
    e[ 2 ] = it[ 2 ][ 2 ];

    // Main iteration loop for numerical modification
    const unsigned  max_iter = 50;
    const double  convergence_threshold = 0.000001;
    double  sum;       // Sum of off-diagonal elements
    double  s, c, t;  // sin, cos, tan of rotation angle
    double  g, h, tz, tau, theta;  // More temporary storage
    double  thresh, thresh_val = 1.0 / ( 5.0 * 3.0 * 3.0 );
    if ( n > 1 ) {
        for ( unsigned  iter = 0;  iter++ < max_iter;  ) {
            // Test for convergence 
            sum = fabs( it[ 0 ][ 1 ] ) + fabs( it[ 0 ][ 2 ] ) + fabs( it[ 1 ][ 2 ] );
            if ( sum < convergence_threshold ) {
                break;  // Convergence reached
            }

            if ( 4 > iter ) {
                thresh = sum * thresh_val;
            } else {
                thresh = 0.0;
            }

            // Modify off diagonal elements in the attempt to make them (near) zero
            for ( unsigned  p = 0;  p < 2;  p++ ) {
                for ( unsigned  q = p + 1;  q < 3;  q++ ) {
                    g = 100.0 * fabs( it[ p ][ q ] );
                    if ( 4 < iter &&
                         g < convergence_threshold ) {
                        it[ p ][ q ] = 0.0;  // Almost or at zero... so make it zero
                    } else if ( thresh < fabs( it[ p ][ q ] ) ) {
                        // Calculate Jacobi transformation
                        if ( g < convergence_threshold ) {
                            t = it[ p ][ q ] / ( e[ q ] - e[ p ] );
                        } else {
                            theta = ( e[ q ] - e[ p ] ) / ( 2.0 * it[ p ][ q ] );
                            if ( 0.0 > theta ) {
                                t = -1.0 / ( sqrt( 1.0 + ( theta * theta ) ) - theta );
                            } else {
                                t = 1.0 / ( sqrt( 1.0 + ( theta * theta ) ) + theta );
                            }
                        }

                        c = 1.0 / sqrt( 1.0 + ( t * t ) );
                        s = t * c;
                        tau = s / ( 1.0 + c );
                        tz = t * it[ p ][ q ];

                        it[ p ][ q ] = 0.0;

                        // Apply Jacobi transformation
                        for ( unsigned  r = 0;  r < p;  r++ ) {
                            g = it[ r ][ p ];
                            h = it[ r ][ q ];
                            it[ r ][ p ] = g - s * ( h + g * tau );
                            it[ r ][ q ] = h + s * ( g - h * tau );
                        }
                        for ( unsigned  r = p + 1;  r < q;  r++ ) {
                            g = it[ p ][ r ];
                            h = it[ r ][ q ];
                            it[ p ][ r ] = g - s * ( h + g * tau );
                            it[ r ][ q ] = h + s * ( g - h * tau );
                        }
                        for ( unsigned  r = q + 1;  r < 3;  r++ ) {
                            g = it[ p ][ r ];
                            h = it[ q ][ r ];
                            it[ p ][ r ] = g - s * ( h + g * tau );
                            it[ q ][ r ] = h + s * ( g - h * tau );
                        }

                        // Update eigenvalues/eigenvectors
                        e[ p ] -= tz;
                        e[ q ] += tz;
                        for ( unsigned  r = 0;  r < 3;  r++ ) {
                            g = m[ r * 3 + p ];
                            h = m[ r * 3 + q ];
                            m[ r * 3 + p ] = g - s * ( h + g * tau );
                            m[ r * 3 + q ] = h + s * ( g - h * tau );
                        }
                    }
                }
            }
        }

        if ( sum > convergence_threshold ) {
	  /////            PUBCHEM_WARN( "Convergence (" << sum << ") failed for inertial tensor diagonalization" );
        }

        // Sort Eigenvalue/eigenvectors
        double  diff = e[ 0 ] - e[ 1 ];  // Be wary of degenerate axes
        if ( diff > 0.0001 ) {
            std::swap( e[ 0 ], e[ 1 ] );

            double  tmp;
            tmp = - m[ 0 ];
            m[ 0 ] = m[ 3 ];
            m[ 3 ] = tmp;

            tmp = - m[ 1 ];
            m[ 1 ] = m[ 4 ];
            m[ 4 ] = tmp;

            tmp = - m[ 2 ];
            m[ 2 ] = m[ 5 ];
            m[ 5 ] = tmp;
//PUBCHEM_INFO( "Swap x->y" );
        }
        diff = e[ 1 ] - e[ 2 ];
        if ( diff > 0.0001 ) {
            std::swap( e[ 1 ], e[ 2 ] );

            double  tmp;
            tmp = - m[ 3 ];
            m[ 3 ] = m[ 6 ];
            m[ 6 ] = tmp;

            tmp = - m[ 4 ];
            m[ 4 ] = m[ 7 ];
            m[ 7 ] = tmp;

            tmp = - m[ 5 ];
            m[ 5 ] = m[ 8 ];
            m[ 8 ] = tmp;
//PUBCHEM_INFO( "Swap y->z" );

            diff = e[ 0 ] - e[ 1 ];
            if ( diff > 0.0001 ) {
                std::swap( e[ 0 ], e[ 1 ] );

                tmp = - m[ 0 ];
                m[ 0 ] = m[ 3 ];
                m[ 3 ] = tmp;

                tmp = - m[ 1 ];
                m[ 1 ] = m[ 4 ];
                m[ 4 ] = tmp;

                tmp = - m[ 2 ];
                m[ 2 ] = m[ 5 ];
                m[ 5 ] = tmp;
//PUBCHEM_INFO( "Swap x->y" );
            }
        }
    }

    /*
    std::cerr << "*** End MATRIX: " << std::endl;

    for ( unsigned  i = 0;  i < 3;  ++i ) {
          for ( unsigned  j = 0;  j < 3;  ++j ) {
	    std::cerr << "\t" << it[i][j];
	  }
	  std::cerr << std::endl;
    }
    */

    if ( NULL != ev ) {
        ev[ 0 ] = e[ 0 ];
        ev[ 1 ] = e[ 1 ];
        ev[ 2 ] = e[ 2 ];
    }
}



void VTransform2StericFrame( float* m,     // float[ 12 ]
				       float* coord, // float[ n * 3 ]
				       const unsigned n,
				       float* ev )   // float[ 3 ]
{
  // Determine rotation necessary to put into inertial frame of reference
  VStericFrameRotation( m, coord, n, ev );

  // Rotate molecule into inertial frame of reference
  VApplyRotTransMatrix( coord, n, m );

  // Determine molecule steric center
  Get3dStericCenterOfGravity( &( m[ 9 ] ), coord, n );
  m[  9 ] = - m[  9 ];
  m[ 10 ] = - m[ 10 ];
  m[ 11 ] = - m[ 11 ];

  //    PUBCHEM_INFO( "PMI :: " << e[ 0 ] << " " << e[ 1 ] << " " << e[ 2 ] );
  //    PUBCHEM_INFO( "Rot :: " << m[ 0 ] << " " << m[ 1 ] << " " << m[ 2 ]
  //                     << " " << m[ 3 ] << " " << m[ 4 ] << " " << m[ 5 ]
  //                     << " " << m[ 6 ] << " " << m[ 7 ] << " " << m[ 8 ] );
  //    PUBCHEM_INFO( "Trans: " << m[ 9 ] << " " << m[ 10 ] << " " << m[ 11 ] );

  // Translate molecule to the steric center
  Translate3dCoords( coord, &( m[ 9 ] ), n );
}



void VApplyRotTransMatrix( float* c,
				     const unsigned n,
				     const float* m )
{
  return VApplyRotTransMatrix( c, c, n, m );
}
// Avoids need to copy array
void VApplyRotTransMatrix( float* c,
				     const float* d,
				     const unsigned n,
				     const float* m )
{
  double  m00 = m[  0 ];
  double  m01 = m[  1 ];
  double  m02 = m[  2 ];
  double  m03 = m[  3 ];
  double  m04 = m[  4 ];
  double  m05 = m[  5 ];
  double  m06 = m[  6 ];
  double  m07 = m[  7 ];
  double  m08 = m[  8 ];
  double  m09 = m[  9 ];
  double  m10 = m[ 10 ];
  double  m11 = m[ 11 ];

  double  x, y, z;
  unsigned  n3 = n * 3;
  for ( unsigned  i = 0;  i < n3;  i += 3 ) {
    x = d[     i ];
    y = d[ 1 + i ];
    z = d[ 2 + i ];

    c[     i ] = static_cast< float >( m09 + m00 * x + m01 * y + m02 * z );
    c[ 1 + i ] = static_cast< float >( m10 + m03 * x + m04 * y + m05 * z );
    c[ 2 + i ] = static_cast< float >( m11 + m06 * x + m07 * y + m08 * z );
  }
}




void TransformCoordToStericFrameAndCalculateEigenvalues(
								     const unsigned natom,
								     float* coord,        // float[ natom * 3 ]
								     float* ev            // float[ 3 ]
								     ){
  float m[12];
  VTransform2StericFrame( m, coord, natom, ev );
}

void TransformCoordToStericFrameAndCalculateQrat(
							      const unsigned natom,
							      const float threshold,
							      float* coord,        // float[ natom * 3 ]
							      unsigned& qrat
							      ){

  float ev[3];
  TransformCoordToStericFrameAndCalculateEigenvalues( natom, coord, ev );
  qrat = CalculateQrat( ev, threshold );
}

unsigned CalculateQrat( const float* ev_ncbi, const float threshold ){



  
  float double_ev_oe[3];
  double_ev_oe[0] = ev_ncbi[1] + ev_ncbi[2] - ev_ncbi[0];
  double_ev_oe[1] = ev_ncbi[0] + ev_ncbi[2] - ev_ncbi[1];
  double_ev_oe[2] = ev_ncbi[0] + ev_ncbi[1] - ev_ncbi[2];
  
  /////  std::cerr << "** 1 ** double_EV_oe: " << double_ev_oe[0] << " " << double_ev_oe[1] << " " << double_ev_oe[2] << std::endl;

  std::sort( double_ev_oe, double_ev_oe+3, std::greater<float>());

  ////  std::cerr << "** 2 ** double_EV_oe: " << double_ev_oe[0] << " " << double_ev_oe[1] << " " << double_ev_oe[2] << std::endl;


  Uint8  u_rqyx, u_rqzy;
  unsigned qrat = 1000;
  if( double_ev_oe[1] > 0 ){
    if ( threshold < ( double_ev_oe[1] / double_ev_oe[0] ) ) {
      u_rqyx = 1;
    } else {
      u_rqyx = 0;
    }
    if ( threshold < ( double_ev_oe[2] / double_ev_oe[1] ) ) {
      u_rqzy = 1;
    } else {
      u_rqzy = 0;
    }
    
    qrat = u_rqyx + u_rqzy;
  }
      
  return qrat;
}

} // namespace Align3D

!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Utilities for *Solver - routines
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 28 Sep 1998
! *
! *****************************************************************************/

!> Basic utilities for matrix stuff.
!------------------------------------------------------------------------------

!> \ingroup ElmerLib
!> \{


MODULE MatrixUtils

#include "../config.h"

   USE Types
   USE ListMatrix
   USE CRSMatrix
   USE BandMatrix
   
   IMPLICIT NONE


CONTAINS

!> Matrix vector multiplication of sparse matrices.
!------------------------------------------------------------------------------
   SUBROUTINE MatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_MatrixVectorMultiply( A,u,v )

     CASE( MATRIX_BAND,MATRIX_SBAND )
       CALL Band_MatrixVectorMultiply( A,u,v )

     CASE( MATRIX_LIST )
       CALL Warn('MatrixVectorMultiply','Not implemented for List matrix type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MatrixVectorMultiply


!------------------------------------------------------------------------------
!> Matrix vector multiplication of sparse matrices.
!------------------------------------------------------------------------------
   SUBROUTINE MaskedMatrixVectorMultiply( A,u,v,ActiveRow,ActiveCol )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
     LOGICAL, DIMENSION(:) :: ActiveRow
     LOGICAL, DIMENSION(:) :: ActiveCol
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_MaskedMatrixVectorMultiply( A,u,v,ActiveRow, ActiveCol )

     CASE DEFAULT
       CALL Fatal('MaskedMatrixVectorMultiply','Not implemented for List matrix type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MaskedMatrixVectorMultiply
!------------------------------------------------------------------------------


!> Matrix vector multiplication of sparse matrices.
!------------------------------------------------------------------------------
   SUBROUTINE TransposeMatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_TransposeMatrixVectorMultiply( A,u,v )

     CASE DEFAULT 
       CALL Fatal('TransposeMatrixVectorMultiply','Not implemented for other than CRS type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE TransposeMatrixVectorMultiply
!------------------------------------------------------------------------------

!> Zero matrix values.
!------------------------------------------------------------------------------
   SUBROUTINE ZeroMatrix( A )
     TYPE(Matrix_t) :: A
     
     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_ZeroMatrix( A )
       
     CASE( MATRIX_BAND,MATRIX_SBAND )
       CALL Band_ZeroMatrix( A )
       
     CASE( MATRIX_LIST )
       CALL List_ZeroMatrix( A % ListMatrix ) 
     END SELECT

   END SUBROUTINE ZeroMatrix
     

 END MODULE MatrixUtils
   

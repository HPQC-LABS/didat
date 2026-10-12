c***********************************************************************
	SUBROUTINE QROD(N,NR,NC,A,R,F,B,GC,GS)
C** Performs ORTHOGONAL DECOMPOSITION OF THE LINEAR LEAST-SQUARES    
C            EQUATION J * X = F TO A * X = B(TRANSPOSE) * F WHERE   
C            J IS THE JACOBIAN IN WHICH THE FIRST N ROWS AND COLUMNS
C            ARE TRANSFORMED TO THE UPPER TRIANGULAR MATRIX A      
C            (J = B * A), X IS THE INDEPENDENT VARIABLE VECTOR, AND
C            F IS THE DEPENDENT VARIABLE VECTOR. THE TRANSFORMATION
C            IS APPLIED TO ONE ROW OF THE JACOBIAN MATRIX AT A TIME.
C  PARAMETERS :                                                   
C      N   -  (INTEGER) DIMENSION OF A TO BE TRANSFORMED.        
C      NR  -  (INTEGER) ROW DIMENSION OF A DECLARED IN CALLING PROGRAM.
C      NC  -  (INTEGER) Column DIMENSION OF F DECLARED IN CALLING PROGRAM.
C      A   -  (REAL*8 ARRAY OF DIMENSIONS .GE. N*N) UPPER TRIANGULAR
C             TRANSFORMATION MATRIX.                               
C      R   -  (REAL*8 LINEAR ARRAY OF DIMENSION .GE. N) ROW OF    
C             JACOBIAN TO BE ADDED.                             
C      F   -  (REAL*8 LINEAR ARRAY .GE. TO THE ROW DIMENSION OF THE
C             JACOBIAN) TRANSFORMED DEPENDENT VARIABLE MATRIX.    
C      B   -  (REAL*8) VALUE OF F THAT CORRESPONDS TO THE ADDED  
C             JACOBIAN ROW.                                     
C     GC   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS COSINE TRANSFORMATIONS.
C     GS   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS SINE TRANSFORMATIONS. 
C--------------------------------------------------------------------
C  AUTHOR : MICHAEL DULICK, Department of Chemistry,
C           UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO N2L 3G1
C--------------------------------------------------------------------
      INTEGER  I,J,K,N,NC,NR
      REAL*8 A(NR,NC), R(N), F(NR), GC(N), GS(N), B, Z(2)
      DO 10 I = 1,N
          Z(1) = R(I)
          J = I - 1
          DO  K = 1,J
              Z(2) = GC(K) * A(K,I) + GS(K) * Z(1)
              Z(1) = GC(K) * Z(1) - GS(K) * A(K,I)
              A(K,I) = Z(2)
              ENDDO
          GC(I) = 1.D0
          GS(I) = 0.D0
          IF(Z(1) .EQ. 0.D0) GOTO 10
          IF(DABS(A(I,I)) .LT. DABS(Z(1))) THEN
              Z(2) = A(I,I) / Z(1)
              GS(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
              GC(I) = Z(2) * GS(I)
            ELSE
              Z(2) = Z(1) / A(I,I)
              GC(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
              GS(I) = Z(2) * GC(I)
            ENDIF
          A(I,I) = GC(I) * A(I,I) + GS(I) * Z(1)
          Z(2) = GC(I) * F(I) + GS(I) * B
          B = GC(I) * B - GS(I) * F(I)
          F(I) = Z(2)
   10     CONTINUE
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

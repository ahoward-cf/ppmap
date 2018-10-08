C**********************************************************************
C     SUBROUTINE INTRP2() PERFORMS BICUBIC INTERPOLATION ON A TWO
C     DIMENSIONAL EQUALLY SPACED GRID.  RETURNS INTERPOLATED VALUE
C     AND FIRST AND SECOND PARTIAL DERIVATIVES IN Z, ZX, ZXX, ZY, ZYY,
C     ZXY. INTERPOLATION INSURES CONTINUOUS Y AND YX ACROSS GRIDPOINTS.
C     USES CUBIC EXTRAPOLATION IF BEYOND GRID RANGE (DANGEROUS!).
C
C     USES BCUCOF(), BCUINT() FROM NUMERICAL RECIPES, P.99
C     *** NOTE: 2ND DERIVATIVES ADDED 5 JULY 1990 ***
C     *** NOTE: CORRECTED FOR XMIN, YMIN .NE. 0 23 SEP 1993 ***
C     AUTHOR: R. JONES  4 JUNE 1990
C**********************************************************************

      SUBROUTINE INTRP2(X,XMIN,XMAX,M,Y,YMIN,YMAX,N,GRID,
     &                  Z,ZX,ZY,ZXX,ZYY,ZXY)

      IMPLICIT DOUBLE PRECISION (A-H,K-L,O-Z)
      DIMENSION GRID( 0:M, 0:N )
      DIMENSION F(4),F1(4),F2(4),F12(4),C(4,4)

      H  = ( XMAX - XMIN ) / M
      XX = ( X    - XMIN ) / H
      MM = INT ( XX )
      K  = ( YMAX - YMIN ) / N
      YY = ( Y    - YMIN ) / K
      NN = INT ( YY )

C     Next line insures all elements of GRID used are in range 0-M,N
      MM = MAX ( 1 , MIN ( M - 2, MM ) )
      XR = XX - MM
      NN = MAX ( 1 , MIN ( N - 2, NN ) )
      YR = YY - NN

C     Assemble information needed by BCUINT()
      XL = MM * H + XMIN
      XU = XL + H   
      YL = NN * K + YMIN
      YU = YL + K
      DO 10 J = NN, NN+1
      DO 10 I = MM, MM+1
         II = 1 + (I-MM) + 3 * (J-NN)
         IF ( II .EQ. 5 ) II = 3
         F(II)   =  GRID(I,J)
         F1(II)  = (GRID(I+1,J) - GRID(I-1,J)) / (2D0 * H)
         F2(II)  = (GRID(I,J+1) - GRID(I,J-1)) / (2D0 * K)
         F12(II) = (GRID(I+1,J+1) - GRID(I+1,J-1) - GRID(I-1,J+1)
     &           +  GRID(I-1,J-1)) / ( 4D0 * H * K )
 10   CONTINUE

      CALL BCUINT( F, F1, F2, F12, XL, XU, YL, YU, X, Y, Z, ZX, ZY,
     &             ZXX, ZYY, ZXY )

      RETURN
      END

      SUBROUTINE BCUINT(Y,Y1,Y2,Y12,X1L,X1U,X2L,X2U,X1,X2,A,A1,A2,
     &                  A11,A22,A12)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(4),Y1(4),Y2(4),Y12(4),C(4,4)
      
      CALL BCUCOF(Y,Y1,Y2,Y12,X1U-X1L,X2U-X2L,C)
      T = (X1 - X1L) / (X1U - X1L)
      U = (X2 - X2L) / (X2U - X2L)
      A  = 0D0
      A1 = 0D0
      A2 = 0D0
      A11 = 0D0
      A22 = 0D0
      DO 10 I = 4, 1, -1
         A   = T * A  + ((C(I,4)*U + C(I,3))*U + C(I,2))*U + C(I,1)
         A2  = T * A2 + (3D0*C(I,4)*U + 2D0*C(I,3))*U + C(I,2)
         A22 = T * A22 + 6D0*C(I,4)*U + 2D0*C(I,3)
         A1  = U * A1 + (3D0*C(4,I)*T + 2D0*C(3,I))*T + C(2,I)
         A11 = U * A11 + 6D0*C(4,I)*T + 2D0*C(3,I)
 10   CONTINUE
      A12 =          C(2,2) + U*(2D0*C(2,3) + 3D0*C(2,4)*U) + 
     &      T * (2D0*C(3,2) + U*(4D0*C(3,3) + 6D0*C(3,4)*U) +
     &      T * (3D0*C(4,2) + U*(6D0*C(4,3) + 9D0*C(4,4)*U)))
      A1  = A1  /  (X1U - X1L)
      A11 = A11 / ((X1U - X1L) * (X1U - X1L))
      A2  = A2  /  (X2U - X2L)
      A22 = A22 / ((X2U - X2L) * (X2U - X2L))
      A12 = A12 / ((X1U - X1L) * (X2U - X2L))
      RETURN
      END
      
      SUBROUTINE BCUCOF(Y,Y1,Y2,Y12,D1,D2,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(4),Y1(4),Y2(4),Y12(4),C(4,4),CL(16),X(16),WT(16,16)
      
      DATA WT/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,
     &  10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,
     &  4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,
     &  10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,
     &  0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,
     &  10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,
     &  5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,
     &  10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
     
      D1D2 = D1 * D2
      DO 10 I = 1, 4
         X(I)    = Y(I)
         X(I+4)  = Y1(I)  * D1
         X(I+8)  = Y2(I)  * D2
         X(I+12) = Y12(I) * D1D2
 10   CONTINUE
      DO 20 I = 1, 16
         XX = 0D0
         DO 30 K = 1, 16
 30         XX = XX + WT(I,K) * X(K)
         CL(I) = XX
 20   CONTINUE
      L = 0
      DO 40 I = 1, 4
         DO 50 J = 1, 4
            L = L + 1
 50         C(I,J) = CL(L)
 40   CONTINUE
      RETURN
      END


!*************************************************************************
SUBROUTINE MALLA(NX,NY,NZ,LADO)
! AVISM MESH CREATION. THE CODE ASSUMES A CUBIC DOMAIN GIVEN BY LADO SUCH THAT
! THE FIRST CELL LEFT CORNER IS AT -LADO/2 AND THE LAST CELL RIGHT CORNER IS
! AT LADO/2. 
! RADi CONTAINS THE COORDINATE OF THE CENTRES OF THE CELLS IN THE i-AXIS
!*************************************************************************
      USE COMMONDATA
      IMPLICIT NONE

!input variables
      INTEGER NX,NY,NZ
      REAL*4 LADO
!local variables
      INTEGER I,J,K
      REAL*4 A,B,C
!ouptut
      REAL*4 RADXX(0:NX+1), RADYY(0:NY+1), RADZZ(0:NZ+1)
      REAL*4 DXX, DYY, DZZ

       A=-LADO/2.0
       B=LADO/2.0

!*     GRID

!*     X-AXIS
      C=(B-A)/(NX-1)
      RADXX(1)=A
      DO I=2,NX
        RADXX(I)=RADXX(1)+(I-1)*C
      END DO

!*     FICTICIUS CELLS
      RADXX(0)=RADXX(1)-C
      RADXX(NX+1)=RADXX(NX)+C

!*     Y-AXIS
      C=(B-A)/(NY-1)
      RADYY(1)=A
      DO J=2,NY
        RADYY(J)=RADYY(1)+(J-1)*C
      END DO

!*     FICTICIUS CELLS
      RADYY(0)=RADYY(1)-C
      RADYY(NY+1)=RADYY(NY)+C

!*     Z-AXIS
      C=(B-A)/(NZ-1)
      RADZZ(1)=A
      DO K=2,NZ
        RADZZ(K)=RADZZ(1)+(K-1)*C
      END DO

!*     FICTICIUS CELLS
      RADZZ(0)=RADZZ(1)-C
      RADZZ(NZ+1)=RADZZ(NZ)+C


      DXX=RADXX(2)-RADXX(1)
      DYY=RADYY(2)-RADYY(1)
      DZZ=RADZZ(2)-RADZZ(1)


      RADX(0:NX+1)=RADXX(0:NX+1)
      RADY(0:NY+1)=RADYY(0:NY+1)
      RADZ(0:NZ+1)=RADZZ(0:NZ+1)
      DX=DXX
      DY=DYY
      DZ=DZZ

      RETURN

!*************************************************************************
END SUBROUTINE MALLA
!*************************************************************************



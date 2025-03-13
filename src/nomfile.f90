
!**********************************************************************
       SUBROUTINE NOMFILE3(DIR, ITER,FILNOMO, FILNOMO1, FILNOMO2, FILNOMO3)
!**********************************************************************
       IMPLICIT NONE
       INTEGER ITER
       CHARACTER(LEN=*) DIR, FILNOMO, FILNOMO2, FILNOMO1, FILNOMO3
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO


       !FILNOMO='output_files/voids'//NOM !list of voids
       !FILNOMO1='output_files/map'//NOM
       !FILNOMO2='output_files/profiles'//NOM
       !FILNOMO3='output_files/profiles_ref'//NOM
       !FILNOMO4='output_files/subvoids'//NOM !list of subvoids
       !FILNOMO5='output_files/map_sub'//NOM

       FILNOMO=TRIM(ADJUSTL(DIR))//'/voids'//NOM !list of voids
       FILNOMO1=TRIM(ADJUSTL(DIR))//'/map'//NOM
       FILNOMO2=TRIM(ADJUSTL(DIR))//'/protovoids'//NOM
       FILNOMO3=TRIM(ADJUSTL(DIR))//'/void_pieces'//NOM

!**********************************************************************
     END SUBROUTINE NOMFILE3
!*************************************************************************


!**********************************************************************
       SUBROUTINE NOMFILEINERTIA(DIR,ITER,FILEINERTIA)
!**********************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER(LEN=*) DIR, FILEINERTIA
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILEINERTIA=TRIM(ADJUSTL(DIR))//'/inertia'//NOM

!*************************************************************************
END SUBROUTINE NOMFILEINERTIA
!*************************************************************************

Program ga

  USE Numtypes
  USE Genetics
  USE Evolution
  USE Statistics

  IMPLICIT NONE

  Integer, Parameter :: Ndim = 1 ! Numero de parametros
  Integer, Parameter :: Ng = 100 ! Numero de bichos por generacion

  Integer, Parameter :: Nt = 20
  Integer :: NVals(Nt)
  Real (kind=DP) :: v1, v2, h, X, RV(Ng)
  
  Type (Generation) :: G1, G2
  Character (len=100) :: fname

  Integer :: Ngene(4) = (/0, Ndim, 0, 0/), I, J, Idx(Ng), Npinta = Ng

  MutProb = 0.03_DP ! Probability that a gene mutes

  ! Inicio la generacio 1
  CALL Init_Generation(G1, Ng)

  ! Cada uno de los Ng organismos de la generacion son 
  ! iniciados con un valor 
  Do I = 1, Ng
     CALL Init_Organism(G1%Member(I), NGene, .true.)
!     G1%Member(I)%Genotype%Rgene = 20.0_DP
     CALL Feval(G1%Member(I))
  End Do

  ! Ahora empieza la evolucion
  Do I = 0, 100

     If (Mod(I, 1) == 0) Then
        CALL Elite(G1, idx)
        Do J = 1, Npinta
           Rv(J) = G1%Member(Idx(J))%Genotype%RGene(1)
        End Do
        Write(*,*)I, Mean(Rv(1:Npinta)), Stddev(Rv(1:Npinta))

        Write(fname, '(1A,1I5.5,1A)')'data/g=', I,'.dat'
        Open (Unit=44, File = Trim(fname))
        Do J = 1, Npinta
           Write(44,'(2ES33.25)')&
                & G1%Member(Idx(J))%Genotype%Rgene(1), &
                & G1%Member(Idx(J))%Fitness
        End Do
        Close(44)
     End If
     
     G2 = Next_Generation(G1)
     G1 = G2
  End Do


  Stop
End Program ga

Subroutine Feval(Or)
  USE NumTypes
  USE Genetics
  
  Type (Organism), Intent (inout) :: Or

  Or%Fitness = (Or%Genotype%RGene(1)**2) + &
       & 40.0_DP * Sin(Or%Genotype%RGene(1))

  Return
End Subroutine Feval
  

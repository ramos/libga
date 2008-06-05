!
! MODULE with routines to perform GA optimizations
!
! Copyright (C) 2005  Alberto Ramos <alberto@martin.ft.uam.es>
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA
! 

! $ v. 1.0; Released: 13/05/2008; $

! ***************************************************
! *
MODULE Evolution
! *
! ***************************************************

  USE NumTypes
  USE Error
  USE Constants
  USE NonNumeric
  USE Genetics

  IMPLICIT NONE

  ! The external function (i.e. problem dependent) Feval
  ! gives the fitness of an organism
  Interface Feval
     Subroutine Feval(Or)
       USE Genetics
       
       Type (Organism), Intent (inout) :: Or
     End Subroutine Feval
  End Interface

  Type Generation
     Type (Organism), Allocatable :: Member(:)
     Integer :: Nmembers
  End Type Generation

  Interface Assignment (=)
     Module Procedure Igual
  End Interface


  Real (kind=DP) :: CrssProb = 0.05_DP
  Integer ::  Nelite = 2


  Private Igual

CONTAINS

! ***********************************
! *
  Subroutine Igual(G1, G2)
! *
! ***********************************

    Type (Generation), intent (out) :: G1
    Type (Generation), intent (in) :: G2

    Integer :: I

    CALL Init_Generation(G1, G2%NMembers)
    Do I = 1, G2%NMembers
       G1%Member(I) = G2%Member(I)
    End Do
    
    Return
  End Subroutine Igual

! ***********************************
! *
  Subroutine Init_Generation(Gen, N)
! *
! ***********************************
    
    Type (Generation), Intent (out) :: Gen
    Integer, Intent (in) :: N

    Integer :: I

    If (Allocated(Gen%Member)) Then
       Do I = 1, Gen%Nmembers
          If (Allocated(Gen%Member(I)%Genotype%CGene)) &
               & Deallocate(Gen%Member(I)%Genotype%CGene)

          If (Allocated(Gen%Member(I)%Genotype%RGene)) &
               & Deallocate(Gen%Member(I)%Genotype%RGene)

          If (Allocated(Gen%Member(I)%Genotype%IGene)) &
               & Deallocate(Gen%Member(I)%Genotype%IGene)

          If (Allocated(Gen%Member(I)%Genotype%SGene)) &
               & Deallocate(Gen%Member(I)%Genotype%SGene)
       End Do
       Deallocate(Gen%Member)
    End If
    ALLOCATE(Gen%Member(N))
    Gen%Nmembers = N
    
    Return
  End Subroutine Init_Generation

! ***********************************
! *
  Function Next_Generation(Gen) Result (Next)
! *
! ***********************************

    Type (Generation), Intent (in) :: Gen
    Type (Generation) :: Next

    Type (Generation) :: Inter
    Real (kind=DP) :: RFitVals(Gen%Nmembers), Rmax, FSum, rnd, &
         & RAcum(Gen%Nmembers), pos
    Integer :: I, Indexx(Gen%Nmembers), Ipos

    CALL Init_Generation(Next, Gen%Nmembers)
    CALL Init_Generation(Inter, Gen%Nmembers)

    Do I = 1, Gen%Nmembers
       RFitVals(I) = Gen%Member(I)%Fitness
    End Do
    CALL QSort(RFitVals, Indexx)
    Rmax = RFitVals(Gen%NMembers)
    
    RFitVals(:) = Rmax - RFitVals

    FSum = Sum(RFitVals)
    CALL Random_Number(rnd)
    rnd = rnd * FSum
    ForAll (I=1:Gen%NMembers) &
         & RAcum(I) = Sum(RFitVals(1:I))

    Do I = 1, Gen%NMembers - Nelite
       pos = Mod(rnd + (I-1)*Fsum/Real(Gen%Nmembers,kind=DP), FSum)
       Ipos = &
            & Locate(RAcum, pos) + 1
       Inter%Member(I) = Gen%Member(Indexx(Ipos))
!       Write(*,*)'Pasa: ', Indexx(Ipos), Inter%Member(I)%Genotype&
!            &%RGene(1), RFitVals(Ipos)
    End Do
    Do I = Gen%NMembers - Nelite + 1, Gen%NMembers
       Inter%Member(I) = Gen%Member(Indexx(I-Gen%NMembers + Nelite))
    End Do

    CALL Permutation(Indexx)
    Do I = 1, Gen%Nmembers - Nelite
       Next%Member(I) = Inter%Member(Indexx(I))
    End Do
    Do I = Gen%NMembers - Nelite + 1, Gen%NMembers
       Next%Member(I) = Inter%Member(I)
    End Do

    Do I = 1, Int( (Gen%NMembers-Nelite)/2)
       CALL Random_Number(rnd)
       If (rnd < CrssProb) Then
          CALL Crossover(Next%Member(2*I-1), Next%Member(2*I))
       End If
    End Do

    Do I = 1, Gen%NMembers-Nelite
       CALL Mutate(Next%Member(I))
       CALL Feval(Next%Member(I))
    End Do

    Do I = Gen%NMembers-Nelite + 1, Gen%NMembers
       CALL Feval(Next%Member(I))
    End Do


    ! Free extra allocated space to avoid memory leaks
    Do I = 1, Gen%NMembers
       If (Allocated(Inter%Member(I)%Genotype%CGene)) Deallocate(Inter%Member(I)%Genotype%CGene)
       If (Allocated(Inter%Member(I)%Genotype%RGene)) Deallocate(Inter%Member(I)%Genotype%RGene)
       If (Allocated(Inter%Member(I)%Genotype%IGene)) Deallocate(Inter%Member(I)%Genotype%IGene)
       If (Allocated(Inter%Member(I)%Genotype%SGene)) Deallocate(Inter%Member(I)%Genotype%SGene)
    End Do
    Deallocate(Inter%Member)


    Return
  End Function Next_Generation

! ***********************************
! *
  Subroutine Elite(Gen, Indexx) 
! *
! ***********************************

    Type (Generation), Intent (in) :: Gen
    Integer, Intent (out) :: Indexx(:)

    Real (kind=DP) :: RFitVals(Gen%NMembers)
    Integer :: I

    Do I = 1, Gen%Nmembers
       RFitVals(I) = Gen%Member(I)%Fitness
    End Do
    CALL QSort(RFitVals, Indexx)
!    Write(*,*)(Gen%Member(Indexx(I))%Genotype%RGene(1), I=1, 10)


    Return
  End Subroutine Elite

End MODULE Evolution


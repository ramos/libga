!
! MODULE with routines to perform GA optimizations
!
! Copyright (C) 2008  Alberto Ramos <alberto@martin.ft.uam.es>
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
MODULE Genetics
! *
! ***************************************************

  USE NumTypes
  USE Error
  USE Constants
  USE Statistics
  USE Time

  Type Genotype
     Complex (kind=DPC), Allocatable :: CGene(:)
     Real (kind=DP), Allocatable :: RGene(:)
     Integer, Allocatable :: IGene(:)
     Character (len=1), Allocatable :: SGene(:)
     Integer :: NCGene, NRGene, NIGene, NSGene
  End Type Genotype

  Type Organism
     Type (Genotype) :: Genotype
     Real (kind=DP) :: Fitness
  End Type Organism

  Interface Assignment (=)
     Module Procedure Igual, Igual_GEN
  End Interface


  Interface RandomizeG
     Module Procedure RandomizeG_C, RandomizeG_R, &
          & RandomizeG_I, RandomizeG_S
  End Interface

  ! ************************
  ! * START of tuneable parameters
  ! * of the model
  ! ************************

  ! Probability that a single Gene mutates
  Real (kind=DP) :: MutProb = 0.1_DP

  ! Standard deviation of the parameter that mix two organism
  ! in the crossover operation. Must be greater that 0.5(?)
  Real (kind=DP) :: StdCrss = 0.5_DP

  ! When initializing organism, Use a normal dstribution with this 
  ! mean and standard deviations.
  Real (kind=DP) :: Rnd_Avg = 0.0_DP
  Real (kind=DP) :: Rnd_Stdd = 10.0_DP

  ! ************************
  ! * END of tuneable parameters
  ! * of the model
  ! ************************


  Private Igual, Igual_GEN, RandomizeG_C, RandomizeG_R, &
          & RandomizeG_I, RandomizeG_S

CONTAINS

! ********************************
! *
  Subroutine Igual_GEN(Gn1, Gn2)
! *
! ********************************

    Type (Genotype), Intent (in) :: Gn2
    Type (Genotype), Intent (out) :: Gn1

    Integer :: NG(4)
    
    Ng = (/Gn2%NCGene, Gn2%NRGene, Gn2%NIGene, Gn2%NSGene/)

    CALL Init_Genotype(Gn1, Ng)
    If (Gn2%NCGene > 0) Then
       Gn1%CGene = Gn2%CGene
    End If

    If (Gn2%NRGene > 0) Then
       Gn1%RGene = Gn2%RGene
    End If

    If (Gn2%NIGene > 0) Then
       Gn1%IGene = Gn2%IGene
    End If

    If (Gn2%NSGene > 0) Then
       Gn1%SGene = Gn2%SGene
    End If

    Return
  End Subroutine Igual_GEN

! ********************************
! *
  Subroutine Igual(Or1, Or2)
! *
! ********************************
    
    Type (Organism), Intent (in) :: Or2
    Type (Organism), Intent (out) :: Or1

    CALL Init_Organism(Or1, &
         & (/Or2%Genotype%NCgene, Or2%Genotype%NRgene, &
         &   Or2%Genotype%NIgene, Or2%Genotype%NSgene/))

    Or1%Genotype = Or2%Genotype
    Or1%Fitness  = Or2%Fitness

    Return
  End Subroutine Igual


! ********************************
! *
  Subroutine Init_Genotype(Gn, Ngene, Irnd)
! *
! ********************************
! * Ngene(1) #Gene of type Complex
! * Ngene(2) #Gene of type Real (DP)
! * Ngene(3) #Gene of type Integer
! * Ngene(4) #Gene of type Character
! ********************************

    Type (Genotype), Intent (out) :: Gn
    Integer, Intent (in) :: Ngene(:)
    Logical, Intent (in), Optional :: Irnd

    Gn%NCGene = 0
    Gn%NRGene = 0
    Gn%NIGene = 0
    Gn%NSGene = 0
    If (Ngene(1) > 0) Then
       If (Allocated(Gn%CGene)) Deallocate(Gn%CGene)
       Allocate (Gn%CGene(Ngene(1)))
       Gn%NCGene = Ngene(1)
       If (Present(Irnd)) Then
          If (Irnd) Then
             CALL RandomizeG(Gn%CGene)
          End If
       End If
    End If

    If (Ngene(2) > 0) Then
       If (Allocated(Gn%RGene)) Then
          Deallocate(Gn%RGene)
          Write(*,*)'Deallocating'
       End If
       Allocate (Gn%RGene(Ngene(2)))
       Gn%NRGene = Ngene(2)
       If (Present(Irnd)) Then
          If (Irnd) Then
             CALL RandomizeG(Gn%RGene)
          End If
       End If
    End If

    If (Ngene(3) > 0) Then
       If (Allocated(Gn%IGene)) Deallocate(Gn%IGene)
       Allocate (Gn%IGene(Ngene(3)))
       Gn%NIGene = Ngene(3)
       If (Present(Irnd)) Then
          If (Irnd) Then
             CALL RandomizeG(Gn%IGene)
          End If
       End If
    End If

    If (Ngene(4) > 0) Then
       If (Allocated(Gn%SGene)) Deallocate(Gn%SGene)
       Allocate (Gn%SGene(Ngene(4)))
       Gn%NSGene = Ngene(4)
       If (Present(Irnd)) Then
          If (Irnd) Then
             CALL RandomizeG(Gn%SGene)
          End If
       End If
    End If

    Return
  End Subroutine Init_Genotype

! ********************************
! *
  Subroutine Init_Organism(Or, Ng, Irnd)
! *
! ********************************

    Type (Organism), Intent (out) :: Or
    Integer, Intent (in) :: Ng(:)
    Logical, Intent (in), Optional :: Irnd

    CALL Init_Genotype(Or%Genotype, Ng, Irnd)
    Or%Fitness = -Huge(1.0_DP)

    Return
  End Subroutine Init_Organism

! ********************************
! *
  Subroutine RandomizeG_C(Arr)
! *
! ********************************

    Complex (kind=DPC), intent (out) :: Arr(:)

    Real (Kind=DP) :: RPart(Size(Arr)), IPart(Size(Arr))
    
    CALL Normal(RPart, 0.0_DP, 10.0_DP)
    CALL Normal(IPart, 0.0_DP, 10.0_DP)
    Arr = Cmplx(Rpart, Ipart, kind=DPC)
    
    Return
  End Subroutine RandomizeG_C

! ********************************
! *
  Subroutine RandomizeG_R(Arr)
! *
! ********************************

    Real (kind=DP), intent (out) :: Arr(:)

    CALL Normal(Arr, Rnd_Avg, Rnd_Stdd)
    
    Return
  End Subroutine RandomizeG_R

! ********************************
! *
  Subroutine RandomizeG_I(Arr)
! *
! ********************************

    Integer, Intent (out) :: Arr(:)

    Real (Kind=DP) :: RPart(Size(Arr))
    
    CALL Normal(RPart, 0.0_DP, 10.0_DP)
    Arr = Int(Rpart)
    
    Return
  End Subroutine RandomizeG_I

! ********************************
! *
  Subroutine RandomizeG_S(Arr)
! *
! ********************************

    Character (len=1), intent (out) :: Arr(:)

    Real (Kind=DP) :: RPart(Size(Arr))
    
    CALL Normal(RPart, 0.0_DP, 10.0_DP)
    Arr = Char(Int(Rpart))
    
    Return
  End Subroutine RandomizeG_S

! ********************************
! *
  Subroutine Mutate(Or)
! *
! ********************************

    Type (Organism), Intent (inout) :: Or

    Real (kind=DP) :: RRnd, CRnd
    Integer :: I

    Do I = 1, Or%Genotype%NCGene
       CALL Random_Number(RRnd)
       If (RRnd < MutProb) Then
          CALL Normal(RRnd)
          CALL Normal(CRnd)
          Or%Genotype%CGene(I) = Or%Genotype%CGene(I)*&
               & (Cmplx(1.0_DP,kind=DPC)+Cmplx(RRnd, CRnd, kind=DPC))
       End If
    End Do

    Do I = 1, Or%Genotype%NRGene
       CALL Random_Number(RRnd)
       If (RRnd < MutProb) Then
          CALL Normal(RRnd)
          Or%Genotype%RGene(I) = Or%Genotype%RGene(I)*&
               & (1.0_DP+RRnd)
       End If
    End Do

    Do I = 1, Or%Genotype%NIGene
       CALL Random_Number(RRnd)
       If (RRnd < MutProb) Then
          CALL Normal(RRnd)
          Or%Genotype%IGene(I) = Int(Or%Genotype%IGene(I)*&
               & (1.0_DP+RRnd))
       End If
    End Do

    Do I = 1, Or%Genotype%NSGene
       CALL Random_Number(RRnd)
       If (RRnd < MutProb) Then
          CALL Normal(RRnd)
          Or%Genotype%SGene(I) = Achar(Int(Or%Genotype%CGene(I)*&
               & (1.0_DP+RRnd)))
       End If
    End Do

    Return
  End Subroutine Mutate

! ********************************
! *
  Subroutine Crossover(Or1, Or2)
! *
! ********************************

    Type (Organism), Intent (inout) :: Or1, Or2

    Complex (kind=DPC) :: CCh1, CCh2
    Real (kind=DP) :: DeltaR, DeltaI, Rh1, Rh2
    Integer :: I, Ih1, Ih2
    Character (len=1) :: Sh1, Sh2

    Do I = 1, Min(Or1%Genotype%NCGene, Or1%Genotype%NCGene)
       CALL Normal(DeltaR, 0.5_DP, StdCrss)
       CALL Normal(DeltaI, 0.5_DP, StdCrss)
       CCh1 = Cmplx(DeltaR, DeltaI, kind=DPC) * &
            & Or1%Genotype%CGene(I) + &
            & Cmplx(1.0_DP-DeltaR, 1.0_DP-DeltaI, kind=DPC) * &
            & Or2%Genotype%CGene(I) 

       CALL Normal(DeltaR, 0.5_DP, StdCrss)
       CALL Normal(DeltaI, 0.5_DP, StdCrss)
       CCh2 = Cmplx(DeltaR, DeltaI, kind=DPC) * &
            & Or1%Genotype%CGene(I) + &
            & Cmplx(1.0_DP-DeltaR, 1.0_DP-DeltaI, kind=DPC) * &
            & Or2%Genotype%CGene(I) 

       Or1%Genotype%CGene(I) = CCh1
       Or2%Genotype%CGene(I) = CCh2
    End Do

    Do I = 1, Min(Or1%Genotype%NRGene, Or1%Genotype%NRGene)
       CALL Normal(DeltaR, 0.5_DP, StdCrss)
       Rh1 = DeltaR * Or1%Genotype%RGene(I) + &
             & (1.0_DP-DeltaR) * Or2%Genotype%RGene(I) 

       CALL Normal(DeltaR, 0.5_DP, StdCrss)
       Rh2 = DeltaR * Or1%Genotype%RGene(I) + &
             & (1.0_DP-DeltaR) * Or2%Genotype%RGene(I) 

       Or1%Genotype%RGene(I) = Rh1
       Or2%Genotype%RGene(I) = Rh2
    End Do

    Do I = 1, Min(Or1%Genotype%NIGene, Or1%Genotype%NIGene)
       CALL Normal(DeltaR, 0.5_DP, StdCrss)
       Ih1 = Int(DeltaR * Or1%Genotype%IGene(I) + &
             & (1.0_DP-DeltaR) * Or2%Genotype%IGene(I) )

       CALL Normal(DeltaR, 0.5_DP, StdCrss)
       Ih2 = Int(DeltaR * Or1%Genotype%IGene(I) + &
             & (1.0_DP-DeltaR) * Or2%Genotype%IGene(I) )

       Or1%Genotype%IGene(I) = Ih1
       Or2%Genotype%IGene(I) = Ih2
    End Do

    Do I = 1, Min(Or1%Genotype%NSGene, Or1%Genotype%NSGene)
       CALL Normal(DeltaR, 0.5_DP, StdCrss)
       Sh1 = Achar(Int(DeltaR * Or1%Genotype%IGene(I) + &
             & (1.0_DP-DeltaR) * Or2%Genotype%IGene(I) ))

       CALL Normal(DeltaR, 0.5_DP, StdCrss)
       Sh2 = Achar(Int(DeltaR * Or1%Genotype%IGene(I) + &
             & (1.0_DP-DeltaR) * Or2%Genotype%IGene(I) ))

       Or1%Genotype%SGene(I) = Sh1
       Or2%Genotype%SGene(I) = Sh2
    End Do

    Or1%Fitness = -Huge(1.0_DP)
    Or2%Fitness = -Huge(1.0_DP)


    Return
  End Subroutine Crossover

! ********************************
! *
  Subroutine Save_Organism(Or, fname, cmt)
! *
! ********************************

    Type (Organism), Intent (in) :: Or
    Character (len=*), Intent (in) :: fname
    Character (len=*), Intent (in), Optional :: cmt
    
    Integer :: I

    Open(Unit=69, File=Trim(fname), ACTION="WRITE")
    If (Present(cmt)) Then
       Write(69,*)Trim(cmt)
    Else
       Write(69,*)asctime(gettime())
    End If
    Write(69,25)Or%Genotype%NCgene, Or%Genotype%NRGene, &
         & Or%Genotype%NIGene, Or%Genotype%NSGene
    Write(69,15)Or%Fitness

    Write(69,15)(Or%Genotype%Cgene(I), I=1, Or%Genotype%NCgene)
    Write(69,15)(Or%Genotype%Rgene(I), I=1, Or%Genotype%NRgene)
    Write(69,25)(Or%Genotype%Igene(I), I=1, Or%Genotype%NIgene)
    Write(69,35)(Or%Genotype%Sgene(I), I=1, Or%Genotype%NSgene)


15  FORMAT((100ES33.25))
25  FORMAT((100I25))
35  FORMAT((100A1))
    Close(69)

    Return
  End Subroutine Save_Organism


! ********************************
! *
  Subroutine Read_Organism(Or, fname)
! *
! ********************************

    Type (Organism), Intent (out) :: Or
    Character (len=*), Intent (in) :: fname
    
    Integer :: I
    Integer :: Ng(4)

    Open(Unit=69, File=Trim(fname), ACTION="READ")
    Read(69,*)
    Read(69,25)(Ng(I), I=1, 4)
    CALL Init_Organism(Or, Ng, .false.)
    Read(69,15)Or%Fitness

    Read(69,15)(Or%Genotype%Cgene(I), I=1, Or%Genotype%NCgene)
    Read(69,15)(Or%Genotype%Rgene(I), I=1, Or%Genotype%NRgene)
    Read(69,25)(Or%Genotype%Igene(I), I=1, Or%Genotype%NIgene)
    Read(69,35)(Or%Genotype%Sgene(I), I=1, Or%Genotype%NSgene)


15  FORMAT((100ES33.25))
25  FORMAT((100I25))
35  FORMAT((100A1))
    Close(69)

    Return
  End Subroutine Read_Organism

End MODULE Genetics


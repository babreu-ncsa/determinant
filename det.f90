!!!!
!! File: det.f90
!! Description: Calculate determinant of a matrix using LAPACK's DGETRF
!! Author: Bruno R. de Abreu  |  babreu at illinois dot edu
!! National Center for Supercomputing Applications (NCSA)
!!  
!! Creation Date: Friday, 27th August 2021, 11:32:02 am
!! Last Modified: Friday, 27th August 2021, 12:17:32 pm
!!  
!! Copyright (c) 2021, Bruno R. de Abreu, National Center for Supercomputing Applications.
!! All rights reserved.
!! License: This program and the accompanying materials are made available to any individual
!!          under the citation condition that follows: On the event that the software is
!!          used to generate data that is used implicitly or explicitly for research
!!          purposes, proper acknowledgment must be provided in the citations section of
!!          publications. This includes both the author's name and the National Center
!!          for Supercomputing Applications. If you are uncertain about how to do
!!          so, please check this page: https://github.com/babreu-ncsa/cite-me.
!!          This software cannot be used for commercial purposes in any way whatsoever.
!!          Omitting this license when redistributing the code is strongly disencouraged.
!!          The software is provided without warranty of any kind. In no event shall the
!!          author or copyright holders be liable for any kind of claim in connection to
!!          the software and its usage.
!!!!

program determinant
    use, intrinsic :: iso_fortran_env
    implicit none
    integer, parameter :: dp = real64   !! double precision
    integer, parameter :: i32 = int32   !! 32 bits integer
    integer(i32), parameter :: ord = 100_i32    !! order of square matrix
    real(dp) :: startT, endT    !! measure calc time
    real(dp), dimension(ord,ord) :: m   !! matrix
    real(dp) :: det

    !! initialize m with random numbers
    call random_seed()
    call random_number(m)

    call cpu_time(startT)
    call getdet(ord,m,det)
    call cpu_time(endT)
    write(*,*) 'det M =', det
    write(*,*) 'Calculated in (s):', (endT-startT)

end program determinant

subroutine getdet(ord,m,det)
    implicit none
!    integer, parameter :: dp = real64   !! double precision
!    integer, parameter :: i32 = int32   !! 32 bits integer
    integer, intent(in) :: ord     !! order of matrix
    real*8, dimension(ord,ord), intent(in) :: m   !! matrix
    real*8, intent(out) :: det    !! determinant
    real*8 :: detp, detu
    integer, dimension(ord) :: ipiv    !! pivot indexes
    integer :: info    !! exit code 
    integer :: i       !! auxiliary
    real*8, dimension(ord,ord) :: m_copy  !! don't want to change m

    m_copy = m
    !! factorization m_copy = P*L*U
    call DGETRF(ord, ord, m_copy, ord, ipiv, info)

    if (info /= 0) then
        write(*,*) 'LU decomposition failed. INFO =', info
        RETURN
    endif
    
    !! calculate det P
    detp = 1.0
    do i = 1, ord
        if (i /= ipiv(i)) then
            detp = -detp
        endif
    enddo

    !! calculate det U
    detu = 1.0
    do i = 1, ord
        detu = detu * m_copy(i,i)
    enddo

    !! finally, det m
    det = detu * detp

end subroutine getdet





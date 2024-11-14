      
      module matsolve_serial
        implicit none

        ! External subroutines from LAPACK
        ! Driver routines
        !  Linear Solver General Matrix
        external SGESV ! Single precision
        external DGESV ! Double precision 
        !  Linear Solver for Band Matrix
        external SGBSV ! Single precision
        external DGBSV ! Double precision
        !  Linear Solve for Tri-diagonal Matrix
        external SGTSV ! Single precision
        external DGTSV ! Double precision

        ! Computational routines
        !  Linear Solver for General Matrix using Factorization
        external SGETRS ! Single precision
        external DGETRS ! Double precision
        !  Linear Solver for Band Matrix using Factorization
        external SGBTRS ! Single precision
        external DGBTRS ! Double precision
        !  Linear Solver for Tri-diagonal Matrix using Factorization
        external SGTTRS ! Single precision
        external DGTTRS ! Double precision
        !  Linear Solver for Triangular Matrix (Packed storage)
        external STRTRS ! Single precision
        external DTRTRS ! Double precision
        !  Inverse Calculator for General Matrix using Factorization
        external SGETRI ! Single precision
        external DGETRI ! Double precision
        !  Inverse Calculator for Triancular Matrix
        external STPTRI ! Single precision
        external DTPTRI ! Double precision


      contains
      ! ---------------------------------------------------------------|
      subroutine solve_linear_eqns_general_r8(workvec, mat, n)
      implicit none
      integer(4),intent(in) :: n  ! Size of the system
      real(8),intent(inout) :: workvec(1:n,1)  ! Input vector to be updated with the solution vector (2D array form is required in LAPACK routines)
      real(8),intent(inout) :: mat(1:n,1:n)  ! Coefficient matrix

      integer(4) :: ipiv(1:n)
      integer(4) :: info

      ! mat() becomes LU factorized matrix after DGESV()
      ! workvec becomes solution vector after DGESV()
      ! info will be zero if the solution is properly obtained.
      ! if info is non-zero, that means mat is singular.
      call DGESV(n,1,mat,n,ipiv,workvec,n, info)

      ! Sanitation check after DEGSV()
      if( info.gt.0 ) then
        write(*,*)'Error, solve_linear_eqns_general_r8(): ' 
        write(*,*)'  The diagonal element of the triangular factor,'
        write(*,*)'  U(',info,',',info,') is zero, so that'
        write(*,*)'  A is singular. The solution not obtained.'
        stop
      endif

      return
      end subroutine solve_linear_eqns_general_r8

      ! ---------------------------------------------------------------|

      subroutine solve_linear_eqns_band_r8(workvec, mat, n, nl, nu)
      implicit none
      ! Band matrix, for example, nl=1, nu=3:
      ! [a b c d 0 ...       0]
      ! [x a b c d 0 ...     0]
      ! [0 x a b c d 0 ...   0]
      ! [0 0 x a b c d 0 ... 0]
      ! [         ...         ]
      ! [      ... 0 x a b c d]
      ! [      ...   0 x a b c]
      ! [      ...     0 x a b]
      ! [      ...       0 x a]
      integer(4),intent(in) :: n  ! Size of the system
      integer(4),intent(in) :: nl ! the number of lower subdiagonals
      integer(4),intent(in) :: nu ! the number of upper subdiagonals
      real(8),intent(inout) :: workvec(1:n,1)  ! Input vector to be updated with the solution vector (2D array form is required in LAPACK routines)
      real(8),intent(inout) :: mat(1:n,1:n)  ! Coefficient matrix

      integer(4) :: ipiv(1:n)
      integer(4) :: info

      real(8) :: ab(1:2*nl+nu+1, 1:n)
      integer(4) :: ldab

      integer(4) :: i,j
      
      ldab = 2*nl+nu+1

      ! Band data storage scheme
      ! Reference: https://github.com/numericalalgorithmsgroup/LAPACK_Examples/blob/master/examples/source/dgbsv_example.f90
      do i=1,n
      do j=max(i-nl,1),max(i+2*nu+nl+1,n)
        ab(nl+nu+1 + i - j, j) = mat(i,j)
      enddo
      enddo

      ! mat() becomes LU factorized matrix after DGBSV()
      ! workvec becomes solution vector after DGBSV()
      ! info will be zero if the solution is properly obtained.
      ! if info is non-zero, that means mat is singular.
      call DGBSV(n,nl,nu,1,ab,ldab,ipiv,workvec,n, info)

      ! Sanitation check after DGBSV()
      if( info.gt.0 ) then
        write(*,*)'Error, solve_linear_eqns_band_r8(): ' 
        write(*,*)'  The diagonal element of the triangular factor,'
        write(*,*)'  U(',info,',',info,') is zero, so that'
        write(*,*)'  A is singular. The solution not obtained.'
        stop
      endif

      return
      end subroutine solve_linear_eqns_band_r8

      ! ---------------------------------------------------------------|

      subroutine solve_linear_eqns_tridiagonal_r8(workvec, dl, d, du, n)
      implicit none
      ! Tridiagonal matrix, n by n
      ! [d du 0 0    ... 0]
      ! [dl d du 0 0 ... 0]
      ! [0 dl d du 0 ... 0]
      ! [        ...      ]
      ! [0   ... 0 dl d du]
      ! [0      ... 0 dl d]
      integer(4),intent(in) :: n
      real(8),intent(inout) :: workvec(1:n,1)
      real(8),intent(in) :: dl(1:n-1)
      real(8),intent(in) :: d(1:n)
      real(8),intent(in) :: du(1:n-1)

      integer(4) :: info

      call DGTSV( n, 1, dl, d, du, workvec, n, info)

      ! Sanitation check after DGTSV()
      if( info.gt.0 ) then
        write(*,*)'Error, solve_linear_eqns_tridiagonal_r8(): ' 
        write(*,*)'  The diagonal element of the triangular factor,'
        write(*,*)'  U(',info,',',info,') is zero, so that'
        write(*,*)'  A is singular. The solution not obtained.'
        stop
      endif

      return
      end subroutine solve_linear_eqns_tridiagonal_r8

      ! ---------------------------------------------------------------|

      subroutine get_tridiagonal_from_matrix_r8(dl, d, du, mat, n)
      implicit none
      integer(4),intent(in) :: n
      real(8),intent(inout) :: dl(1:n-1)
      real(8),intent(inout) :: d(1:n)
      real(8),intent(inout) :: du(1:n-1)
      real(8),intent(in) :: mat(1:n,1:n)

      integer(4) :: i

      ! lower diagonal
      !dl(1) = mat(1+1,1)
      !dl(2) = mat(2+1,2)
      !dl(i) = mat(i+1,i)
      !dl(n-1) = mat(n-1 +1, n)
!
      ! diagonal
      !d(1) = mat(1,1)
      !d(2) = mat(2,2)
      !d(i) = mat(i,i)
      !d(n) = mat(n,n)

      ! upper diagonal
      !du(1) = mat(1,1+1)
      !du(2) = mat(2,2+1)
      !du(i) = mat(i,i+1)
      !du(n-1) = mat(n-1,n-1 + 1)

      do i=1,n-1
        dl(i) = mat(i+1,i)
        d(i) = mat(i,i)
        du(i) = mat(i,i+1) 
      enddo
      d(n) = mat(n,n)

      return
      end subroutine get_tridiagonal_from_matrix_r8

      ! ---------------------------------------------------------------|

      subroutine put_tridiagonal_to_matrix_r8(mat, dl,d,du,n)
      implicit none
      integer(4),intent(in) :: n
      real(8),intent(inout) :: mat(1:n,1:n)
      real(8),intent(in) :: dl(1:n-1)
      real(8),intent(in) :: d(1:n)
      real(8),intent(in) :: du(1:n-1)

      integer(4) :: i

      do i=1,n-1
        mat(i+1,i) = dl(i)
        mat(i,i) = d(i)
        mat(i,i+1) = du(i) 
      enddo
      mat(n,n) = d(n)      

      return
      end subroutine put_tridiagonal_to_matrix_r8

      ! ---------------------------------------------------------------|

      subroutine print_matrix_r8(mat, nrow, ncol, desc)
      implicit none
      integer(4) :: nrow
      integer(4) :: ncol
      real(8) :: mat(1:nrow,1:ncol)
      character(len=*) :: desc ! Description for the printed matrix

      integer(4) ::  i, j

      write(*,'(A)') desc
      do i=1,nrow
          write(*,'(11(:,1x,f15.8))') ( mat( i, j ), j = 1, ncol )
      enddo
      
      return
      end subroutine print_matrix_r8




      end module matsolve_serial

      program main

      use mpi
      use matsolve_serial

      implicit none

      integer(4),parameter :: n=3
      real(8) :: a(n,n)
      real(8) :: v(n,1)

      real(8) :: b(n)

      real(8) :: dl(n-1)
      real(8) :: d(n)
      real(8) :: du(n-1)

      a(1,1)=10.d0
      a(1,2)=2.d0
      a(1,3)=3.d0
      a(2,1)=4.d0
      a(2,2)=50.d0
      a(2,3)=6.d0
      a(3,1)=0.d0
      a(3,2)=8.d0
      a(3,3)=90.d0

      b(1)=0.1d0
      b(2)=0.2d0
      b(3)=0.3d0

      a(:,:) = 0.0d0

      dl(:) = 1.0d0
      d(:) = -2.0d0
      du(:) = 1.0d0

      call put_tridiagonal_to_matrix_r8(a, dl,d,du,n)
        
      v(:,1) = b(:)

      call print_matrix_r8(a, n,n, 'Coefficient matrix, a()')

      call print_matrix_r8(v, n,1, 'Work vector, v()')

      write(*,*) 'Matrix solver wrapper for LAPACK' 

      !call solve_linear_eqns_general_r8(v, a, n) ! Works
      !call solve_linear_eqns_band_r8(v, a, n, 1, 2)  ! Works

      call get_tridiagonal_from_matrix_r8(dl,d,du, a,n)
      call solve_linear_eqns_tridiagonal_r8(v, dl,d,du,n)

      call print_matrix_r8(a, n,n, 'After the solution, a()')

      call print_matrix_r8(v, n,1, 'Solution vector, v()')








      end program main
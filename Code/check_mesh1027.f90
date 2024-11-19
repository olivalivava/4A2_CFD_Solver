      
      subroutine check_mesh(g)

!     Check the cell area and facet length calculations before attempting to
!     solve the flow, make sure you do this for both the "bump" and "bend" test
!     cases

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_grid), intent(inout) :: g
      real :: area_min, dx_error, dy_error, tol
      integer :: ni, nj
!     *************
      integer :: i,j
      real :: sum_dx(g%ni-1,g%nj-1), sum_dy(g%ni-1,g%nj-1)
      real :: x1,x2,x3,x4,y1,y2,y3,y4
!     *************
!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Exact checking of floating point numbers never goes well, define a
!     small tolerance to use for a comparative operation instead
      tol = 1e-4 * g%l_min

!     Check that all of the cell areas are positive, either with the intrinsic
!     "minval" function or with nested do loops. Print the output to the screen
!     and flag negative numbers as an error with an if statement to "stop" the
!     program
!     INSERT
!     *******************
      area_min = minval(g%area)

!     To check if there is any negative cell area
      if (area_min < 0.0) then
         write(6,*) "ERROR: Negative cell area found!!!!!"
         stop
      else
         write(6,*) "All cell area are checked to be positive."
      end if

      write(6,*)
!     ******************* 
!     Next check that the sum of the edge vectors around every quadrilateral is 
!     very nearly zero in both the x and y-coordinate directions. You can
!     complete this with some elementwise addition of the arrays and use of the
!     "maxval" and "abs" intrinsic functions.
!     INSERT
!     ********************
      do i=1, ni-1
         do j=1, nj-1
            x1 = (g%x(i+1,j)-g%x(i,j))
            x2 = (g%x(i+1,j+1)-g%x(i+1,j))
            x3 = (g%x(i,j+1)-g%x(i+1,j+1))
            x4 = (g%x(i,j)-g%x(i,j+1))
            sum_dx(i,j) = x1 + x2 + x3 + x4
            y1 = (g%y(i+1,j)-g%y(i,j))
            y2 = (g%y(i+1,j+1)-g%y(i+1,j))
            y3 = (g%y(i,j+1)-g%y(i+1,j+1))
            y4 = (g%y(i,j)-g%y(i,j+1))
            sum_dy(i,j) = y1 + y2 + y3 + y4
         enddo
      enddo

      dx_error = maxval (abs(sum_dx))
      dy_error = maxval (abs(sum_dy))

      if (dx_error > tol) THEN
         print *, "dx_error = ", dx_error
         write(6,*) "ERROR: x-direction edge vectors do not sum up to zero."
         stop
      else if (dy_error > tol) THEN
         write(6,*) "ERROR: y-direction edge vectors do not sum up to zero."
         stop
      else
         write(6,*) "Edge vector sum checked."
      end if 
      
!     ********************      
!     It may be worthwhile to complete some other checks, the prevous call to
!     the "write_output" subroutine has written a file that you can read and
!     postprocess using the Python script plot_mesh.py. This program also has
!     access to all of the mesh parameters used within the solver that you could
!     inspect graphically.

!     Print a blank line
      write(6,*)

      end subroutine check_mesh

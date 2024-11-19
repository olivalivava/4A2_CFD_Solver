      
      subroutine generate_mesh(geom,g)

!     Create cells of the mesh to cover the domain defined by geometry curves,
!     the values of the node coordinates, x(i,j) and y(i,j) are stored in "g"

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_geometry), intent(in) :: geom
      type(t_grid), intent(inout) :: g
      real :: si_a(geom%ni_a), si_b(geom%ni_b), si(g%ni)
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERT
!     *************
      integer :: i,j
      real :: wt(g%nj),sj(g%nj)
!     (wt is the weight function in between 0 and 1 for nj entries.)
!     *************
      
!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Calculate the non-dimensional curve lengths of the geometry input and
!     generate linearly spaced points in i-direction at desired mesh resolution
      call dist(geom%x_a,geom%y_a,1,si_a)
      call dist(geom%x_b,geom%y_b,1,si_b)
      call linspace(0.0,1.0,si)

!     Interpolate the geometry curves to the required resolution in the 
!     i-direction, this allows the mesh to be refined without altering the 
!     geometry definition file, the data is stored at "j = 1" and "j = nj"
      call interp(si_a,geom%x_a,si,g%x(:,1))
      call interp(si_a,geom%y_a,si,g%y(:,1))
      call interp(si_b,geom%x_b,si,g%x(:,nj))
      call interp(si_b,geom%y_b,si,g%y(:,nj))

!     **************** 
!      write(6,*) 'x-coordinates'
!      write(6,*) g%x(1,1), g%x(ni,nj)
!     ****************
!     Calculate the coordinates of all the intermediate points within the mesh.
!     Create a new vector of non-dimensional spacings in the j-direction using 
!     "linspace", loop over the mesh in the i-direction and calculate the
!     intermediate coordinates from a weighted sum of the two boundaries
!     INSERT
!     *******************
!      call linspace(0.0, 1.0, wt)
!      do i=1, ni
!         do j=1, nj
!            g%x(i,j)= (1-wt(j))*(g%x(i,1))+ wt(j)*(g%x(i,nj))
!            g%y(i,j)= (1-wt(j))*(g%y(i,1))+ wt(j)*(g%y(i,nj))
!         enddo
!      enddo   
      call linspace(0.0, 1.0, sj)
      do i =1, ni
         g%x(i,:) = (1-sj)*g%x(i,1) + sj*g%x(i,nj)
         g%y(i,:) = (1-sj)*g%y(i,1) + sj*g%y(i,nj)
      enddo
      
!     ******************* 
!     In all of the test cases for the basic solver the the "j = 1" and "j = nj"
!     boundaries are walls, for the extensions you may need to return to this
!     and communicate the position of the walls to the solver in a more 
!     flexible way. The "wall" variable is an "ni * nj" logical array where 
!     "true" indicates the node is on a wall.
      g%wall = .false.
      g%wall(:,[1,g%nj]) = .true.

!     Print that the mesh has been created
      write(6,*) 'Interpolated mesh from the bounding geometry curves'
!      write(6,*) g%x(1,1), g%x(1,2),g%x(1,3),g%x(1,4),g%x(1,5)
!      write(6,*) g%x(2,1), g%x(2,2),g%x(2,3),g%x(2,4),g%x(2,5)
!      write(6,*) g%x(3,1), g%x(3,2),g%x(3,3),g%x(3,4),g%x(3,5)
!      write(6,*) g%x(ni,1), g%x(ni,2),g%x(ni,3),g%x(ni,4),g%x(ni,5) 
      
      end subroutine generate_mesh



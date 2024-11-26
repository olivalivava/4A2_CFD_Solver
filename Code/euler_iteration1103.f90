
      subroutine euler_iteration(av,g)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use flux_stencil
      use smooth_stencil
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      real, dimension(g%ni,g%nj-1) :: mass_i, flux_i, momx_i, momy_i
      real, dimension(g%ni-1,g%nj) :: mass_j, flux_j, momx_j, momy_j
      real, dimension(g%ni, g%nj) :: prop,dcell
      integer :: i, j, ni, nj

!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
!     INSERT
!     *************************
      do i=1, ni-1
         do j=1, nj-1
            mass_i(i,j) = g%ro(i,j)*(g%vx(i,j)*g%ly_i(i,j) + g%vy(i,j)*g%lx_i(i,j)) !ro*v*A, cell side not facet length
            mass_j(i,j) = g%ro(i,j)*(-g%vx(i,j)*g%ly_j(i,j) + g%vy(i,j)*g%lx_j(i,j))
         enddo
      enddo
!     *************************
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERT
!     *************************
      prop(:,:) = g%ro (:, :)
      call sum_fluxes(av,mass_i,mass_j,g%area,prop, dcell)
      g%dro = prop(1:ni-1,1:nj-1) - g%ro(1:ni-1, 1:nj-1)
      g%ro(:, :) = prop(:,:)
!      g%dro(:,:) = dcell(1:ni-1, 1:nj-1)
!     *************************      
!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERT
!     **************************
      flux_i(:,1:nj-1) = mass_i(:,1:nj-1)* g%hstag(:,1:nj-1)
      flux_j(1:ni-1,:) = mass_j(1:ni-1,:)* g%hstag(1:ni-1,:)
!     ************************** 
!     Update the internal energy with enthalpy fluxes
!     INSERT
!     **************************
      prop = g%roe(:,:)
      call sum_fluxes(av,flux_i,flux_j,g%area,prop, dcell)
      g%droe = prop(1:ni-1,1:nj-1) - g%roe(1:ni-1, 1:nj-1)
      g%roe(:,:) = prop
!     g%droe = dcell(1:ni-1, 1:nj-1)
!     **************************      
!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERT
!     *************************
      momx_i(:,1:nj-1) = mass_i(:,1:nj-1)* g%vx(:,1:nj-1)+ g%p(:,1:nj-1)*g%lx_i(:,1:nj-1)
      momx_j(1:ni-1,:) = mass_j(1:ni-1,:)* g%vx(1:ni-1,:)+ g%p(1:ni-1,:)*g%lx_j(1:ni-1,:)
!      write(6,*) momx_j(2,2), momx_i(2,17),flux_i(2,2), mass_i(2,2), flux_j(2,17), mass_j(2,17)
!     *************************      
!     Update the x-momentum with momentum flux
!     INSERT
!     *************************
      prop = g%rovx(:,:)
      call sum_fluxes(av,momx_i,momx_j,g%area,prop, dcell)
      g%drovx = prop(1:ni-1,1:nj-1) - g%rovx(1:ni-1, 1:nj-1)
      g%rovx(:,:) = prop
!      g%drovx = dcell(1:ni-1, 1:nj-1)
!     *************************      
!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERT
!     *************************
      momy_i(:,1:nj-1) = mass_i(:,1:nj-1)* g%vy(:,1:nj-1)+ g%p(:,1:nj-1)*g%ly_i(:,1:nj-1)
      momy_j(1:ni-1,:) = mass_j(1:ni-1,:)* g%vy(1:ni-1,:)+ g%p(1:ni-1,:)*g%ly_j(1:ni-1,:)
!     *************************      
!     Update the y-momentum with momentum flux
!     INSERT
!     *************************
      prop = g%rovy(:,:)
      call sum_fluxes(av,momy_i,momy_j,g%area,prop, dcell)
      g%drovy = prop(1:ni-1,1:nj-1) - g%rovy(1:ni-1, 1:nj-1)
      g%rovy(:,:) = prop
!      g%drovy = dcell(1:ni-1, 1:nj-1)
!     ************************* 
!     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g%ro)
      call smooth_array(av,g%roe)
      call smooth_array(av,g%rovx)
      call smooth_array(av,g%rovy)
      

      end subroutine euler_iteration



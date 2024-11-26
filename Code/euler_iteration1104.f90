
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
      real, dimension(g%ni,g%nj-1) :: mass_i, flux_i
      real, dimension(g%ni-1,g%nj) :: mass_j, flux_j
      real, dimension(g%ni, g%nj) :: prop, dnode
      integer :: i, j, ni, nj, smoothtype
      smoothtype = av%smoothtype

!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
!     INSERT
!     ************************
      mass_i = ((g%rovx(:,1:nj-1)+g%rovx(:,2:nj))*0.5)*g%lx_i(:,1:nj-1) &
           + ((g%rovy(:,1:nj-1)+g%rovy(:,2:nj))*0.5)*g%ly_i(:,1:nj-1)
      mass_j = ((g%rovx(1:ni-1,:)+g%rovx(2:ni,:))*0.5)*g%lx_j(1:ni-1,:) &
               + ((g%rovy(1:ni-1,:)+g%rovy(2:ni,:))*0.5)*g%ly_j(1:ni-1,:)
!     *************************
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERT
!     *************************
      prop = g%ro
      call sum_fluxes(av,mass_i,mass_j,g%area, g%ro, g%dro)
      g%ro = g%ro_start + (g%ro - prop)

!     *************************      
!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERT
!     **************************
      flux_i = mass_i* ((g%hstag(:,1:nj-1) + g%hstag(:,2:nj))/2.0)
      flux_j = mass_j* ((g%hstag(1:ni-1,:) + g%hstag(2:ni,:))/2.0)
!     ************************** 
!     Update the internal energy with enthalpy fluxes
!     INSERT
!     **************************
      prop = g%roe
      call sum_fluxes(av,flux_i,flux_j,g%area, g%roe, g%droe)
      g%roe = g%roe_start + (g%roe - prop)
!     **************************      
!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERT
!     *************************
      flux_i = mass_i* ((g%vx(:,1:nj-1)+g%vx(:, 2:nj))/2.0)+ ((g%p(:,1:nj-1)+g%p(:,2:nj))/2.0)*g%lx_i(:,1:nj-1)
      flux_j = mass_j* ((g%vx(1:ni-1,:)+g%vx(2:ni,:))/2.0)+ ((g%p(1:ni-1,:)+g%p(2:ni,:))/2.0)*g%lx_j(1:ni-1,:)
!     *************************      
!     Update the x-momentum with momentum flux
!     INSERT
!     *************************
      prop = g%rovx
      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovx, g%drovx)
      g%rovx = g%rovx_start + (g%rovx - prop)
!     *************************      
!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERT
!     *************************
      flux_i = mass_i* ((g%vy(:,1:nj-1)+g%vy(:, 2:nj))/2.0)+ ((g%p(:,1:nj-1)+g%p(:,2:nj))/2.0)*g%ly_i(:,1:nj-1)
      flux_j = mass_j* ((g%vy(1:ni-1,:)+g%vy(2:ni,:))/2.0)+ ((g%p(1:ni-1,:)+g%p(2:ni,:))/2.0)*g%ly_j(1:ni-1,:)
!     *************************      
!     Update the y-momentum with momentum flux
!     INSERT
!     *************************
      prop = g%rovy
      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovy_start, g%drovy)
      g%rovy = g%rovy_start + (g%rovy - prop)
!     ************************* 
!     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g%ro,g%corr_ro, smoothtype)
      call smooth_array(av,g%roe,g%corr_roe, smoothtype)
      call smooth_array(av,g%rovx,g%corr_rovx, smoothtype)
      call smooth_array(av,g%rovy,g%corr_rovy, smoothtype)

      end subroutine euler_iteration



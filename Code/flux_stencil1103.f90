
      module flux_stencil

!     Packaging a subroutine in a module allows it to recieve the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine sum_fluxes(av,flux_i,flux_j,area,prop,dcell)

!     This subroutine sums the fluxes into each cell, calculates the change in 
!     the cell property inside, distributes the change to the four nodes of the
!     cell and then adds it onto the flow property

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: flux_i(:,:), flux_j(:,:), area(:,:)
      real, intent(inout) :: prop(:,:)
      real, intent(out) :: dcell(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: dnode
      integer :: ni, nj, i, j

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)

!     Use the finite volume method to find the change in the variables "prop"
!     over the timestep "dt", save it in the array "dcell"
!     INSERT
!     ********************
!      do i=1, ni-1
!         do j=1, nj-1
!            dcell(i,j) = (flux_i(i,j) - flux_i(i+1,j) + flux_j(i,j) - flux_j(i,j+1))*av%dt/area(i,j)
!            // spatial derivative terms in Euler equation, using divergence, not time evolution            
!            print *, i, j, dcell(i,j)
!         enddo
!      enddo

      dcell(1:ni-1, 1:nj-1) = (flux_i(1:ni-1,1:nj-1) - flux_i(2:ni,1:nj-1) + &
           flux_j(1:ni-1,1:nj-1) - flux_j(1:ni-1,2:nj))*av%dt/area(1:ni-1,1:nj-1)
      
!     Top and bottom edges (j = 1&nj)
!      do i = 2, ni-1
!         dcell(i,1)  = dcell(i,2)   
!         dcell(i,nj) = dcell(i,nj-1) 
!      end do
      
!     Left and right edges (i = 1&ni)
!      do j = 2, nj-1
!         dcell(1,j)  = dcell(2,j)   
!         dcell(ni,j) = dcell(ni-1,j) 
!      end do

!     Corner cells
!      dcell(1,1)   = dcell(2,2)
!      dcell(1,nj)  = dcell(2,nj-1)
!      dcell(ni,1)  = dcell(ni-1,2)
!      dcell(ni,nj) = dcell(ni-1,nj-1)

!     [Used the neighbour values to boost the code efficiency without affect too much accuracy]     
!     ********************      
!     Now distribute the changes equally to the four corners of each cell. Each 
!     interior grid point receives one quarter of the change from each of the 
!     four cells adjacent to it.
!     INSERT
!     **********************
!      do i=2, ni-1
!         do j=2, nj-1
!            dnode(i,j) = 0.25 * (dcell(i,j) + dcell(i-1,j) + dcell(i,j-1) + dcell(i-1,j-1))
!         enddo
!      enddo

      dnode(2:ni-1,2:nj-1) = 0.25 * (dcell(2:ni-1,2:nj-1) + dcell(1:ni-2, 2:nj-1) +&
           dcell(2:ni-1, 1:nj-2) + dcell(1:ni-2, 1:nj-2))
!     ********************** 
!     Bounding edge nodes do not have four adjacent cells and so must be treated
!     differently, they only recieve half the change from each of the two
!     adjacent cells. Distribute the changes for the "i = 1 & ni" edges as well
!     as the "j = 1 & nj" edges. 
!     INSERT
!     **********************
!     Top and bottom edges (j = 1 & nj)
      do i = 2, ni-1
         dnode(i,1) = 0.5 * (dcell(i,1) + dcell(i-1,1))
         dnode(i,nj) = 0.5 * (dcell(i,nj-1) + dcell(i-1,nj-1))
      end do
      
!     Left and right edges (i = 1 & ni)
      do j = 2, nj-1
         dnode(1,j) = 0.5 * (dcell(1,j) + dcell(1,j-1))
         dnode(ni,j) = 0.5 * (dcell(ni-1,j) + dcell(ni-1,j-1))
      end do
      
!     **********************      
!     Finally distribute the changes to be to the four bounding corner points, 
!     these receive the full change from the single cell of which they form one 
!     corner.
!     INSERT
!     ***********************
      dnode(1,1)   = dcell(1,1)
      dnode(1,nj)  = dcell(1,nj-1)
      dnode(ni,1)  = dcell(ni-1,1)
      dnode(ni,nj) = dcell(ni-1,nj-1)
!     **********************      
!     Update the solution by adding the changes at the nodes "dnode" to the flow
!     property "prop"
!     INSERT
!     ***********************
      prop(:,:) = prop(:,:) + dnode(:,:)
!     ************************      
      end subroutine sum_fluxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module flux_stencil



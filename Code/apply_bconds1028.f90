      
      subroutine apply_bconds(av,g,bcs)

!     This subroutine applies both the inlet and outlet boundary conditions, as
!     it modifies both the primary and secondary flow variables they must be
!     calculated first

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(inout) :: bcs

!     Declare the other variables you need here
!     INSERT
!     ***************
      real :: rostag
      real :: t_in(g%nj),v_in(g%nj), E(g%nj)
      integer :: nj,ni
      real, dimension(g%ni, g%nj) :: t_static, v

      nj = g%nj; ni = g%ni
!     *************** 
!     At the inlet boundary the change in density is driven towards "rostag",
!     which is then used to obtain the other flow properties to match the
!     specified stagnation pressure, temperature and flow angle. 

!     To help prevent instabilities forming at the inlet boundary condition the 
!     changes in inlet density are relaxed by a factor "rfin" normally set to 
!     0.25 but it can be reduced further.

!     It is also worth checking if "ro" is greater than "rostag" and limiting 
!     the values to be slightly less than "rostag". This can prevent the solver 
!     crashing during severe transients.
      if(av%nstep == 1) then
          bcs%ro = g%ro(1,:)
      else
          bcs%ro = bcs%rfin * g%ro(1,:) + (1 - bcs%rfin) * bcs%ro
      endif
      bcs%ro = min(bcs%ro,0.9999 * bcs%rostag)

!     Calculate "p(1,:)", "rovx(1,:)", "rovy(1,:)" and "roe(1,:)" from the inlet 
!     "ro(:)", "pstag", "tstag" and "alpha". Also set "vx(1,:)", "vy(1,:)" and 
!     "hstag(1,:)"
!     INSERT
!     ****************
      !rostag = bcs%pstag/(av%rgas*bcs%tstag)
      !      write(6,*) bcs%pstag, g%ro(1,10), rostag, av%gam
      t_static(1,:) = bcs%tstag*(bcs%ro(:)/bcs%rostag)**(av%gam-1)
      v(1,:) = (2*av%cp*(bcs%tstag-t_static(1,:)))**0.5
      g%roe(1,:) = (av%cv *t_static(1,:) + 0.5*(v(1,:)**2))*bcs%ro(:)
      g%p(1,:) = bcs%ro(:) * av%rgas*t_static(1,:)
      g%rovx(1,:) = bcs%ro(:)*v(1,:)*cos(bcs%alpha)
      g%rovy(1,:) = bcs%ro(:)*v(1,:)*sin(bcs%alpha)
      
      !E(:) = av%cv*t_in(:) + 0.5*v_in(:)**2
      !g%roe(1,:) = bcs%ro*E(:)

      g%vx(1,:) = v(1,:)*cos(bcs%alpha)
      g%vy(1,:) = v(1,:)*sin(bcs%alpha)
      g%hstag(1,:) = av%cp * bcs%tstag
!      write(6,*) nj,g%p(1,nj),t_in(nj),v_in(nj),g%ro(1,nj),g%vx(1,nj),g%vy(1,nj),g%roe(1,nj),bcs%alpha,g%hstag(1,nj)

      
!     **************** 
!     For the outlet boundary condition set the value of "p(ni,:)" to the
!     specified value of static pressure "p_out" in "bcs"
!     INSERT
!     ****************
      g%p(ni,:) = bcs%p_out


!      write(6,*) "Boundary Condition set:"
!      write(6,*) "At inlet: p =", g%p(1,1), "ro = ", g%ro(1,1), "rovx =",g%rovx(1,1), "rovy = ", g%rovy(1,1), "roe = ", g%roe(1,1)
!      write(6,*) "At outlet: p = ", g%p(ni,1)
!     ****************      
      end subroutine apply_bconds



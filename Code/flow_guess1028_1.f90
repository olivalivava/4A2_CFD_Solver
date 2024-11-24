

      subroutine flow_guess(av,g,bcs,guesstype)

!     This calculates an initial guess of the primary flowfield variables and
!     stores them at the nodes within the mesh dataytype

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(in) :: bcs
      integer, intent(in) :: guesstype
      integer :: i, j, ni, nj, j_mid
      
!     Variables required for the crude guess
      real :: t_out, v_out, ro_out, lx, ly, l

!     Variables required for the improved guess, you will need to add to these
      real :: l_i(g%ni)
!     INSERT
!     ****************
      real :: l_i_lc(g%nj-1)
      real :: v_guess(g%ni), ro_guess(g%ni), t_guess(g%ni), M_guess(g%ni)
      real :: p_guess(g%ni)
      real :: R, gamma, M_out, p0, t0, mass_out, p_out, cp, cv
      real :: mach_lim, t_lim, t_new, tol, t_min, E_tot
      integer :: count, max_itr
!     **************** 
!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;
!     ***************
!     Store the values locally for convenience
      R = av%rgas; gamma = av%gam
      cp = av%cp; cv = av%cv
      p0 = bcs%pstag; t0 = bcs%tstag; p_out=bcs%p_out

!     Set tolerance for iteration:
      tol = 1.0e-4
!     **************
!     Assuming isentropic flow to the the exit plane calculate the static
!     temperature and the exit velocity
      t_out = bcs%tstag * (bcs%p_out / bcs%pstag)**av%fgam
      v_out = (2 * av%cp * (bcs%tstag - t_out))**0.5
      ro_out = bcs%p_out / (av%rgas * t_out)

!     Determine which guess calcation method to use by the value of "guesstype"
      if(guesstype == 1) then

!         Store the exit density and internal energy as if they were uniform
!         (upper two graphs have no change with displacement
          g%ro = ro_out 
          g%roe  = g%ro * (av%cv * t_out + 0.5 * v_out**2)

!         Calculate the gradient of the mesh lines in the centre of the domain
!         to determine the assumed direction of the flow
          j_mid = nj / 2
          do i = 1,ni-1
              lx = g%lx_j(i,j_mid); ly = g%ly_j(i,j_mid); 
              l = hypot(lx,ly)
              g%rovx(i,:) = g%ro(i,:) * v_out * ly / l
              g%rovy(i,:) = -g%ro(i,:) * v_out * lx / l
          end do

!         Copy the values to the "i = ni" nodes as an approximation
          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)

!         Print the guess that has been calculated
          write(6,*) 'Crude flow guess calculated'
          write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
               g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*)

      else if(guesstype == 2) then 

!         Calculate the length of each "i = const" line between the "j = 1" and 
!         "j = nj" boundaries of the domain and store it in the local variable
!         "l_i". You could calculate the length along each i-facet from the x 
!         and y projected lengths with "hypot" and then sum them up in the
!         second dimension with "sum". 
!         INSERT
!         ***************
         do i=1, ni
            do j=1, nj-1
               l_i_lc(j) = hypot(g%lx_i(i,j), g%ly_i(i,j))
            end do
            l_i(i) = sum (l_i_lc(:), dim=1)
!            write(6,*) "l_i(", i,") = ", l_i(i)
         enddo
       
!         write(6,*) "l_i = ", l_i(i)
!         ***************         
!         Use the exit temperature, density and velocity calculated for the 
!         crude guess with "l_i" to estimate the mass flow rate at the exit
!         INSERT
!         ***************
         M_out = v_out/ ((av%gam *av%rgas *t_out)**0.5)
         mass_out = ro_out*l_i(ni)*v_out
!         write(6,*) "ro*A*V value = ", mass_out
         
!         *************** 
!         Set a limit to the maximum allowable mach number in the initial
!         guess, call this "mach_lim", calculate the corresponding temperature,
!         called "t_lim"
!         INSERT
!         ******************
      !    mach_lim = 1.0
      !    M_guess = min(M_guess, 0.9999*mach_lim)
      !    t_lim = t0/(1+0.5*((gamma-1)*(mach_lim**2)))

      !    REMOVED THE SUBSONIC FLOW ONLY LIMITATION

!         ******************
!         Now estimate the velocity and density at every "i = const" line, call 
!         the velocity "v_guess(i)" and the density "ro_guess(i)":
!             1. Assume density is constant at the exit value
!             2. Use continuity and "l_i(i)" to estimate the velocity
!             3. Assume stagnation temperature is constant for static temp
!             4. Limit the static temperature, lookup intrinsic "max"
!             5. Calculate the density throughout "ro_guess(i)"
!             6. Update the estimate of the velocity "v_guess(i)" 
!         INSERT
!         *******************
        ! E_tot = cv*t_out + 0.5*((v_out)**2)
         do i=1, ni
            v_guess(i) = mass_out/(l_i(i)*ro_out)
            t_guess(i) = bcs%tstag - (v_guess(i)**2)/(2.0*av%cp)

!             if (t_guess(i) < t_lim) then
! !              write(6,*) 'Detected Mach Number > 1, and applied the Mach Limit Successfully.'
!                t_guess(i) = t_lim
!                v_guess(i) = sqrt(av%gam* av%rgas*t_guess(i))
!             endif
            
            M_guess(i) = v_guess(i)/sqrt(av%gam* av%rgas*t_guess(i))
            ro_guess(i) = bcs%rostag * ((t_guess(i)/bcs%tstag)**(1/(av%gam-1)))
            v_guess(i) = mass_out/(l_i(i)*ro_guess(i))
         end do

      !    t_min = minval(t_guess, dim=1)
      !    if (t_min < t_lim) then
      !       write(6,*) 'M=' ,maxval(M_guess), 'at position i =', maxloc(M_guess)
      !       write(6,*) "ERROR: Temperature exceeds the limit. Lower down the initial Mach number."
      !       stop
      !    else
      !       write(6,*) "Checked the intial mach number."
      !    end if
!        [Notes] Stagnation Temperature constant, so total energy kept constant         
         
!         *******************
!         Direct the calculated velocity to be parallel to the "j = const"
!         gridlines for all values of i and j. This can be achieved with a 
!         similar calculation to the "j = nj/2" one that was performed in the 
!         crude guess. Then set all of ro, roe, rovx and rovy, note that roe 
!         includes the kinetic energy component of the internal energy.
!         INSERT 
!         ********************
         do i=1, ni
            g%ro(i,:) = ro_guess(i)
            g%roe(i,:) = g%ro(i,:) * (av%cv*t_guess(i) + 0.5 * (v_guess(i)**2))
          end do  
                  
          do i = 1,ni-1
             do j=1, nj
              lx = g%lx_j(i,j); ly = g%ly_j(i,j); 
              l = hypot(lx,ly)
!              print *, i,j_mid,lx,ly,l, g%ro(i,18), p_guess(i), M_guess(i), v_guess(i), ro_guess(i)
              g%rovx(i,j) = g%ro(i,j) * v_guess(i) * ly / l
              g%rovy(i,j) = -g%ro(i,j) * v_guess(i) * lx / l
             end do
          end do

          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)
!         ********************         
!         Make sure the guess has been copied for the "i = ni" values too
!         INSERT
!         ********************
!          g%ro(ni,:) = g%ro(ni-1, :)
!          g%roe(ni,:) = g%roe(ni-1, :)
!          g%rovx(ni,:) = g%rovx(ni-1,:)
!          g%rovy(ni,:) = g%rovy(ni-1,:)
!          g%ro(:,nj) = g%ro(:, nj-1)
!          g%roe(:,nj) = g%roe(:, nj-1)
!          g%rovx(:,nj) = g%rovx(:,nj-1)
!          g%rovy(:,nj) = g%rovy(:,nj-1)
!         ********************
!         Print the first elements of the guess like for the crude guess
!         INSERT
!         ********************
          write(6,*)
          write(6,*) 'Option 2 flow guess calculated'
          write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
               g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*)
!         ********************         
      end if

!     The initial guess values derived from the boundary conditions are also
!     useful as a reference to non-dimensionalise the convergence
      av%ro_ref = sum(g%ro(1,:)) / (nj*1.0)
      av%roe_ref = sum(g%roe(1,:)) / (nj*1.0)
      av%rov_ref = max(sum(g%rovx(1,:)),sum(g%rovy(1,:))) / (nj*1.0)

      end subroutine flow_guess



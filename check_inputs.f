      subroutine check_inputs( ierr )

c  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
c    Author: Thomas S. Lund, lund@cora.nwra.com

c   *** Error codes:
c          0 - normal exit
c          6 - required input not found
c          7 - static dimensions too small
c          9 - inconsistent input parameter values

      include 'sam.h'

      ierr = 0

c   *** Make sure the static dimensions are sufficient

      if( Nx .gt. mLx .or. Ny .gt. mLy .or. Nz .gt. mLz ) then
         print *, 'ERROR: static dimensions are too small'
         print *, 'mLx, mLy, mLz = ', mLx, mLy, mLz
         ierr = 7
         return
      end if

c   *** Set default planes output positions

      i = index_param( 'n_skip_p', labels, n_inputs )

      if( found(i) .and. k_xy_plane(1) .eq. 0 ) then
         k_xy_plane(1) = Nz/2
         ii = index_param( 'k_xy_plane1', labels, n_inputs )
         values(ii) = k_xy_plane(1)
      end if

      if( found(i) .and. j_xz_plane(1) .eq. 0 ) then
         j_xz_plane(1) = Ny/2
         ii = index_param( 'j_xz_plane1', labels, n_inputs )
         values(ii) = j_xz_plane(1)
      end if

      if( found(i) .and. i_yz_plane(1) .eq. 0 ) then
         i_yz_plane(1) = Nx/2
         ii = index_param( 'i_yz_plane1', labels, n_inputs )
         values(ii) = i_yz_plane(1)
      end if

c   *** Determine the number of planes to be written

      do i=1, 9
         if( k_xy_plane(i) .eq. 0 ) goto 1
      end do
1     n_xy_planes = i-1

      do i=1, 9
         if( j_xz_plane(i) .eq. 0 ) goto 2
      end do
2     n_xz_planes = i-1

      do i=1, 9
         if( i_yz_plane(i) .eq. 0 ) goto 3
      end do
3     n_yz_planes = i-1

c   *** If this is a gravity wave problem check inputs for consistency
c   *** and determine which two defining parameters are specified.

      i_gw_type = 0

      if( abs(i_prob) .eq. 4 ) then
         i0 = index_param( 'lambda_x', labels, n_inputs )
         i1 = index_param( 'gam'     , labels, n_inputs )
         i2 = index_param( 'omega'   , labels, n_inputs )
         i3 = index_param( 'uo'      , labels, n_inputs )
         if( .not. found(i0) ) then
            print *, 'ERROR: lambda_x is a required input for ',
     &               'a GW problem'
            ierr = 9
            return
         end if
         if( found(i1) .and. found(i2) .and. found(i3) ) then
            print *, 'ERROR: You can only specify two of Gam, ',
     &               'omega, Uo for a GW problem'
            ierr = 9
            return
         end if
         if(      found(i1) .and. found(i2) ) then
            i_gw_type = 1
         else if( found(i1) .and. found(i3) ) then
            i_gw_type = 2
         else if( found(i2) .and. found(i3) ) then
            i_gw_type = 3
         else
            print *, 'ERROR: You must specify exactly two of Gam, ',
     &               'omega, Uo for a GW problem'
            ierr = 6
            return
         end if
         if( i_strat .ne. 1 ) then
            print *, 'ERROR: i_strat must equal 1 for a GW problem'
            ierr = 6
            return
         end if
         if( i_prob .eq. -4 .and. 
     &       (i_gw_type .eq. 1 .or. i_gw_type .eq. 2) .and.
     &       abs(xL/zL-Gam) .gt. 1.0e-8 ) then
            print *, 'ERROR: The computational box aspect ratio does ',
     &               'agree with the GW parameter Gam'
            ierr = 6
            return
         end if
      end if

c   *** Perform an initial check for missing required inputs

      do i=1, n_inputs
         if( required(i) .and. .not.found(i) ) goto 50
      end do

c   *** Promote selected optional inputs to required based on
c   *** other read-in parameter values.

      if( i_restart .eq. 1 ) then
         i = index_param( 'nt_restart', labels, n_inputs )
         required(i) = .true.
      end if

      if( i_les .eq. 1 ) then
         i = index_param( 'c_smag', labels, n_inputs )
         required(i) = .true.
      end if

      if( i_cfl .eq. 0 ) then
         i = index_param( 'dt0', labels, n_inputs )
         required(i) = .true.
      end if

      if( i_force .eq. 1 ) then
         i = index_param( 'k_force', labels, n_inputs )
         required(i) = .true.
      end if

      if( i_strat .eq. 1 ) then
         i = index_param( 'lapse0', labels, n_inputs )
         required(i) = .true.
         i = index_param( 'grav', labels, n_inputs )
         required(i) = .true.
         i = index_param( 'z0', labels, n_inputs )
         required(i) = .true.
         i = index_param( 'to', labels, n_inputs )
         required(i) = .true.
         i = index_param( 'pr', labels, n_inputs )
         required(i) = .true.
      end if

      if( i_prob .eq. 5 ) then
         i = index_param( 'shear', labels, n_inputs )
         required(i) = .true.
      end if

      if( i_fs .ne. 0 ) then
         i = index_param( 'lambda_z_fs', labels, n_inputs )
         required(i) = .true.
         i = index_param( 'Ri_fs', labels, n_inputs )
         required(i) = .true.
      end if

c   *** Check again for missing required inputs

      do i=1, n_inputs
         if( required(i) .and. .not.found(i) ) goto 50
      end do

c   *** Check for inconsistent parameter values

      if( i_restart .lt. 0 .or. i_restart .gt. 1 ) then
         i = index_param( 'i_restart', labels, n_inputs )
         goto 60
      end if

      if( i_les .lt. 0 .or. i_les .gt. 2 ) then
         i = index_param( 'i_les', labels, n_inputs )
         goto 60
      end if

      if( i_stat .lt. 0 .or. i_stat .gt. 2 ) then
         i = index_param( 'i_stat', labels, n_inputs )
         goto 60
      end if

      if( i_cfl .lt. 0 .or. i_cfl .gt. 2 ) then
         i = index_param( 'i_cfl', labels, n_inputs )
         goto 60
      end if

      if( nrk_max .lt. 0 .or. nrk_max .gt. 3 ) then
         i = index_param( 'nrk_max', labels, n_inputs )
         goto 60
      end if

      if( n_dealias .lt. 0 .or. n_dealias .gt. 3 ) then
         i = index_param( 'n_dealias', labels, n_inputs )
         goto 60
      end if

      if( i_prob .lt. -4 .or. i_prob .gt. 5 ) then
         i = index_param( 'i_prob', labels, n_inputs )
         goto 60
      end if

      if( i_force .lt. 0 .or. i_force .gt. 1 ) then
         i = index_param( 'i_force', labels, n_inputs )
         goto 60
      end if

      if( k_force .lt. 0 .or. k_force .gt. max(Nx/2,Ny/2,Nz/2) ) then
         i = index_param( 'k_force', labels, n_inputs )
         goto 60
      end if

      return

50    print *, 'ERROR: Required input ',trim(labels(i)),
     &         ' was not provided'
      ierr = 6
      return

60    print *, 'ERROR: Invalid value for parameter ',trim(labels(i)),
     &         ' =', values(i)
      ierr = 9
      return

      end


      function index_param( param, labels, n_params )

      character( *) param
      character(12) labels(n_params), label_lower
      character(len(param)) param_lower

      call to_lower( param, param_lower, len(param) )

      do i=1, n_params
         call to_lower( labels(i), label_lower, 12 )
         if( param_lower .eq. label_lower ) then
            index_param = i
            return
         end if
      end do

      index_param = 0

      return
      end

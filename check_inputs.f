      subroutine check_inputs( labels, values, found, required, ierr )

c  © 2002 – 2014 NorthWest Research Associates, Inc. All Rights Reserved
c    Author: Thomas S. Lund, lund@cora.nwra.com

c   *** Error codes:
c          0 - normal exit
c          6 - required input not found
c          7 - static dimensions too small
c          9 - inconsistent input parameter values

      include 'sam.h'

      real          values(*)
      logical        found(*), required(*)
      character(12) labels(*)

      ierr = 0

c   *** Make sure the static dimensions are sufficient

      if( Nx .gt. mLx .or. Ny .gt. mLy .or. Nz .gt. mLz ) then
         print *, 'ERROR: static dimensions are too small'
         print *, 'mLx, mLy, mLz = ', mLx, mLy, mLz
         ierr = 7
         return
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

      if( i_prob .lt. 0 .or. i_prob .gt. 2 ) then
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

      character(* ) param
      character(12) labels(n_params)

      do i=1, n_params
         if( trim(param) .eq. trim(labels(i)) ) then
            index_param = i
            return
         end if
      end do

      index_param = 0

      return
      end

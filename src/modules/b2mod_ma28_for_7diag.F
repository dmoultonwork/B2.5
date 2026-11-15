      module b2mod_ma28_for_7diag                                 !sv 11.01.00
      use b2mod_types
      use b2mod_ysmp_sdrv
      implicit none

      integer, allocatable, save :: iwrk_7p(:)
      real (kind=R8), allocatable, save ::
     1 aaf_7p(:), wrk0_7p(:), wrk1_7p(:), wrk2_7p(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE (iwrk_7p, aaf_7p, wrk0_7p, wrk1_7p, wrk2_7p)
#endif

      real (kind=R8), allocatable, save ::
     * r1_7p(:),
     * s1_7p(:),
     * a1h_7p(:),
     * w_7p(:),
     * res_7p(:),
     * a1h_save_7p(:)

#ifdef _OPENMP
!$OMP THREADPRIVATE(r1_7p, s1_7p, a1h_7p, w_7p, res_7p, a1h_save_7p)
#endif

      integer, allocatable, save ::
     * irn_7p(:),
     * icn_7p(:),
     * ish_7p(:),
     * jsh_7p(:),
     * ijsh_7p(:),
     * irnh_7p(:),
     * icnh_7p(:),
     * icnhkeep_7p(:),
     * ikeep_7p(:,:),
     * iw28_7p(:,:),
     * irnh_save_7p(:),
     * icnh_save_7p(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE(irn_7p, icn_7p, ish_7p, jsh_7p, ijsh_7p)
!$OMP THREADPRIVATE(irnh_7p, icnh_7p, icnhkeep_7p)
!$OMP THREADPRIVATE(ikeep_7p, iw28_7p, irnh_save_7p, icnh_save_7p)
#endif

      integer, save :: nc0_7p,nnd_7p
      integer, save :: n_nonzeros_7p,solvdim_7p,solvdim1_7p
      integer, save :: nltrsol_7p
      logical, save :: lfirst_ilu7 = .true.
#ifdef _OPENMP
!$OMP THREADPRIVATE(nc0_7p, nltrsol_7p, nnd_7p, n_nonzeros_7p)
!$OMP THREADPRIVATE(solvdim_7p, solvdim1_7p)
!$OMP THREADPRIVATE(lfirst_ilu7)
#endif

c     Per-species pivot state of the 7-diagonal solver. The thread
c     private state is overwritten by COPYIN at every parallel region
c     and lost at its end, and the master's state belongs to some other
c     species, so a threaded solve started from a poor factorisation
c     and fell back to a full analysis. Each species keeps the state
c     of its own previous solve here (pivot order and index arrays only;
c     the matrix values, work arrays and saved copies are rebuilt on
c     every call). The serial code reaches the same effect by solving
c     similar species one after the other.
      type ma28_7_species_state
        logical :: valid = .false.
        integer, allocatable :: ikeep(:,:)
        integer, allocatable :: icnhkeep(:), icnh(:), irnh(:)
      end type ma28_7_species_state
      integer, parameter :: ma28_7_nstate = 128
      type(ma28_7_species_state), save ::
     1 ma28_7_state(0:ma28_7_nstate-1)

      contains

      subroutine save_ma28_7_state(is)
c     Called by the thread that has just solved species is.
      implicit none
      integer, intent(in) :: is
      if (is.lt.0 .or. is.ge.ma28_7_nstate) return
      if (lfirst_ilu7 .or. .not.allocated(ikeep_7p)) return
      if (.not.allocated(ma28_7_state(is)%ikeep)) then
        allocate(ma28_7_state(is)%ikeep(size(ikeep_7p,1),
     1   size(ikeep_7p,2)))
        allocate(ma28_7_state(is)%icnhkeep(size(icnhkeep_7p)))
        allocate(ma28_7_state(is)%icnh(size(icnh_7p)))
        allocate(ma28_7_state(is)%irnh(size(irnh_7p)))
      end if
      ma28_7_state(is)%ikeep = ikeep_7p
      ma28_7_state(is)%icnhkeep = icnhkeep_7p
      ma28_7_state(is)%icnh = icnh_7p
      ma28_7_state(is)%irnh = irnh_7p
      ma28_7_state(is)%valid = .true.
      return
      end subroutine save_ma28_7_state

      subroutine load_ma28_7_state(is)
c     Called before the solve of species is. With no saved state the
c     thread keeps whatever state it has.
      implicit none
      integer, intent(in) :: is
      if (is.lt.0 .or. is.ge.ma28_7_nstate) return
      if (.not.ma28_7_state(is)%valid .or. .not.allocated(ikeep_7p))
     1 return
      ikeep_7p = ma28_7_state(is)%ikeep
      icnhkeep_7p = ma28_7_state(is)%icnhkeep
      icnh_7p = ma28_7_state(is)%icnh
      irnh_7p = ma28_7_state(is)%irnh
      lfirst_ilu7 = .false.
      return
      end subroutine load_ma28_7_state

      subroutine alloc_b2mod_ma28_for_7diag                           !sv 11.01.00
      implicit none

      if (.not.allocated(r1_7p))  then

        write(*,*) 'Allocating space for b2mod_ma28_for_7diag'        !sv 11.01.00
        write(*,*) 'nnd        = ',nnd_7p
        write(*,*) 'n_nonzeros = ',n_nonzeros_7p
        write(*,*) 'solvdim    = ',solvdim_7p
        write(*,*) 'solvdim1   = ',solvdim1_7p

        allocate(r1_7p(1:nnd_7p))
        allocate(s1_7p(1:nnd_7p))
        allocate(a1h_7p(1:solvdim_7p))
        allocate(w_7p(1:nnd_7p))
        allocate(res_7p(1:nnd_7p))
        allocate(a1h_save_7p(1:n_nonzeros_7p))
        allocate(irn_7p(1:n_nonzeros_7p))
        allocate(icn_7p(1:n_nonzeros_7p))
        allocate(ish_7p(1:n_nonzeros_7p))
        allocate(jsh_7p(1:n_nonzeros_7p))
        allocate(ijsh_7p(1:n_nonzeros_7p))
        allocate(irnh_7p(1:solvdim1_7p))
        allocate(icnh_7p(1:solvdim_7p))
        allocate(icnhkeep_7p(1:solvdim_7p))
        allocate(ikeep_7p(1:nnd_7p,5))
        allocate(iw28_7p(1:nnd_7p,8))
        allocate(irnh_save_7p(1:n_nonzeros_7p))
        allocate(icnh_save_7p(1:n_nonzeros_7p))

        write(*,*) 'Allocated space for b2mod_ma28_for_7diag'         !sv 11.01.00
      end if
      return
      end subroutine alloc_b2mod_ma28_for_7diag                       !sv 11.01.00

      subroutine dealloc_b2mod_ma28_for_7diag

#ifdef _OPENMP
!$OMP PARALLEL
#endif
      if (allocated(r1_7p)) then
        deallocate(r1_7p)
        deallocate(s1_7p)
        deallocate(a1h_7p)
        deallocate(w_7p)
        deallocate(res_7p)
        deallocate(a1h_save_7p)
        deallocate(irn_7p)
        deallocate(icn_7p)
        deallocate(ish_7p)
        deallocate(jsh_7p)
        deallocate(ijsh_7p)
        deallocate(irnh_7p)
        deallocate(icnh_7p)
        deallocate(icnhkeep_7p)
        deallocate(ikeep_7p)
        deallocate(iw28_7p)
        deallocate(irnh_save_7p)
        deallocate(icnh_save_7p)
      endif
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      return
      end subroutine dealloc_b2mod_ma28_for_7diag                 !sv 11.01.00

      subroutine restart_ma28_for_7                               !srv 11.01.00
      implicit none
      lfirst_ilu7=.true.
      write(*,*) 'Restarting ma28_for_7'
      return
      end subroutine restart_ma28_for_7

      subroutine solver_init_7p(nx, ny)
      use b2mod_indirect
      implicit none
      integer, intent(in) :: nx, ny
      integer, save :: style
      integer, save :: mult_nonzero, mult_solvdim, mult_solvdim1
      integer :: n2
      external ipgeti, xertst
      data style/2/

cxpb This routine to be called once to initialize the arrays used by
cxpb the matrix solver. Must be done once only and outside the iteration
cxpb loop to avoid errors when using OpenMP parallelization.

*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, itcnt
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..alocate arrays
#ifdef _OPENMP
!$OMP PARALLEL
!$OMP& SHARED(nx, ny, n2)
#endif
      call ipgeti ('b2ux7p_style', style)
      call xertst (0.le.style.and.style.lt.4,
     &  'faulty internal parameter style')
      write(*,*)
     .  'Allocating space for solver for 7-point stencil matrix'  !srv 11.01.00
      if(style.eq.2) then
        mult_nonzero=10
        mult_solvdim=15
        mult_solvdim1=10
        call ipgeti ('b2ux7p_mult_nonzero', mult_nonzero)
        call ipgeti ('b2ux7p_mult_solvdim', mult_solvdim)
        call ipgeti ('b2ux7p_mult_solvdim1', mult_solvdim1)
        nnd_7p=(nx+2)*(ny+2)
        n_nonzeros_7p=mult_nonzero*nnd_7p
        solvdim_7p=mult_solvdim*n_nonzeros_7p
        solvdim1_7p=mult_solvdim1*n_nonzeros_7p
        call alloc_b2mod_ma28_for_7diag
        call iniilu_7(nx+2,ny+2)
        nltrsol_7p=2
        call ipgeti ('b2ux7p_nltrsol', nltrsol_7p)
      else if (style.eq.3) then
        nnd_7p=(nx+2)*(ny+2)
        nnd_sdrv  = (nx+2)*(ny+2)
        n_nz_sdrv = 4*nnd_sdrv
        nsp = 3*nnd_7p + 16*n_nz_sdrv                             !srv 29.12.15
        call alloc_b2mod_ysmp_sdrv
        ia = 0                                                    !srv 15.05.07
        ja = 0                                                    !srv 15.05.07
      else if (style.eq.0) then
        call xerrab ('Style not applicable to 7-point stencil!')
*    ..obtain workspace for iluter
        allocate(aaf_7p(1:n2*5))
        allocate(wrk0_7p(1:n2))
        allocate(wrk1_7p(1:n2))
        allocate(wrk2_7p(1:n2))
        allocate(iwrk_7p(1:1))
      else if (style.eq.1) then
        call xerrab ('Style not applicable to 7-point stencil!')
*    ..obtain workspace for slv5pt
        allocate(aaf_7p(1:n2*(3*min(nx+2,ny+2)+1)))
        allocate(wrk0_7p(1:n2))
        allocate(wrk1_7p(1:1))
        allocate(wrk2_7p(1:1))
        allocate(iwrk_7p(1:n2))
      endif
      write(*,*)
     .  'Allocated space for solver for 7-point stencil matrix'   !srv 11.01.00
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      return
      end subroutine solver_init_7p

      end module b2mod_ma28_for_7diag                             !srv 11.01.00

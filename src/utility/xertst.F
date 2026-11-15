      subroutine xertst (lcond, messg)
      use b2mod_types
      implicit none
      logical lcond
      character messg*(*)
*     ------------------------------------------------------------------
*     XERTST tests a condition; if violated, calls xerrab to terminate
*     execution.
*     ------------------------------------------------------------------
      external xerrab
*     ------------------------------------------------------------------
      if (.not.lcond) then
       call xerrab (messg)
      endif
      return
*     ------------------------------------------------------------------
      end subroutine xertst

      subroutine b25_mpi_stop
#if ( defined(USE_MPI) && MPI_MOD > 0 )
      use mpi     ! IGNORE
#endif
      implicit none
#ifdef USE_MPI
#if MPI_MOD < 1
      include 'mpif.h'
#endif
      integer ierr_mpi
      logical init_mpi
#ifdef B25_EIRENE
      logical ltrigger
      integer req_trig

      ltrigger=.false.
#endif
      call mpi_initialized(init_mpi,ierr_mpi)
      if (init_mpi) then
#ifdef B25_EIRENE
c       non-blocking, to match the sleeping wait of the child ranks
        call mpi_ibcast(ltrigger,1,MPI_LOGICAL,0,
     .                   MPI_COMM_WORLD,req_trig,ierr_mpi)
        call mpi_wait(req_trig,MPI_STATUS_IGNORE,ierr_mpi)
#endif
        call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
        call mpi_finalize(ierr_mpi)
      end if
#endif
      return
      end subroutine b25_mpi_stop

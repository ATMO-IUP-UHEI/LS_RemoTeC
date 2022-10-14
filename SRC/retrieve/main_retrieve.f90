!------------------------------------------------------------------------------
!> \mainpage RemoTeC simulation package
!------------------------------------------------------------------------------

program main
   use header_module
   use wrapper_retrieve_module
   use read_miprep_module, only: open_miprep, close_miprep
   use omp_lib
   implicit none

   type(shared_data), pointer :: fixedData
   type(pixel_data), pointer :: varyingData
   type(output_data), pointer :: outputData
   character(stringlen), dimension(:), allocatable :: atm
   integer, dimension(:), allocatable :: pixelid                          ! Identifier for groundpixel
   character(stringlen) :: runpath, atm_file, arg, first_atm
   integer :: i, nfile, runid, ierr, io
   real(double) :: stoptime, starttime

   !*** Read in command line arguments runid, atm_file and runpath
   if (iargc() .eq. 2) then
      call getarg(1, arg)
      read (arg, *) runid
      if (runid .lt. 0 .or. runid .gt. 999999) call stopretrieval('RUNID must be between 0 and 999999')
      call getarg(2, atm_file)
      runpath = './'
   elseif (iargc() .eq. 3) then
      call getarg(1, arg)
      read (arg, *) runid
      call getarg(2, atm_file)
      call getarg(3, runpath)
   else
      print *, 'Give correct number of arguments'
      stop
   end if

   !*** Get nfile = total number of input spectra to be retrieved
   !*** read in first processing file to extract spectral grid
   open (newunit(io), file=trim(atm_file), action='read')
   nfile = 1
   read (io, *, iostat=ierr) first_atm
   do while (ierr .ne. -1)
      nfile = nfile + 1
      read (io, *, iostat=ierr)
   end do
   close (io)
   nfile = nfile - 1

   !*** Read shared data into memory
   call init_shared(fixedData, runpath, first_atm, runid, ierr)
   if (ierr .ne. 0) call stopretrieval("MAIN: error init fixed data")

   !*** Get filenames of input spectra to be retrieved
   allocate (atm(nfile), pixelid(nfile))
   open (newunit(io), file=trim(atm_file), action='read')
   if (fixedData%flag%atm == 3) then
      do i = 1, nfile
         read (io, *) atm(i), pixelid(i)
      end do
      call open_miprep(trim(fixedData%path%meteo)//atm(1), ierr)
      if (ierr .ne. 0) call stopretrieval("MAIN: error opening MIPrep files")
   else
      do i = 1, nfile
         read (io, '(A)') atm(i)
      end do
   end if
   close (io)

   !*** For measuring CPU time
   call cpu_time(starttime)

   !$OMP PARALLEL private(varyingData, outputData, ierr)
   !$OMP DO
   do i = 1, nfile
      call init_pixel(fixedData, varyingData, outputData, atm(i), pixelid(i), ierr)
      if (ierr == 0) then
         call retrieve_wrapper(fixedData, varyingData, outputData, ierr)
         !$OMP CRITICAL(writefile)
         call write_output(runID, fixedData, varyingData, fixedData%flag%atm, outputData)
         !$OMP END CRITICAL(writefile)
      end if
      call free_pixel(varyingData, outputData)
   end do
   !$OMP END DO
   !$OMP END PARALLEL

   if (fixedData%flag%atm == 3) then
      call close_miprep(atm(1), ierr)
      if (ierr .ne. 0) call stopretrieval("MAIN: error closing MIPrep files")
   end if
   call free_shared(fixedData)

   !*** Mean CPU time per retrieval
   call cpu_time(stoptime)
   write (6, '(a, f8.3, a )') ' Mean cpu time per retrieval = ', (stoptime - starttime)/nfile, ' s'
   !*** Number of retrievals
   write (6, '(a, i6 )') ' Number of retrievals = ', nfile

end program main

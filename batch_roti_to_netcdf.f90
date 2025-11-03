program batch_roti_to_netcdf
  use rinex_module
  use, intrinsic :: ieee_arithmetic
  use netcdf
  implicit none

  integer, parameter :: dp_local = kind(1.0d0)
  character(len=*), parameter :: path = "/mnt/c/Users/jijin/Desktop/JIF/test_ROTi/"
  character(len=512) :: file_list(1000)
  integer :: nfiles, i, ncid
  integer :: dim_file, dim_x, dim_y, dim_loc
  integer :: var_name, var_loc, var_roti
  integer :: ierr
  real(dp_local), allocatable :: tec_arr(:,:), roti_arr(:,:)
  real(dp_local) :: rec_loc(3)
  real(dp_local) :: dt
  character(len=512) :: fname_full, fname_only
  integer :: p

  ! ---------------------------------------------------------
  ! 1. Get file list
  ! ---------------------------------------------------------
  call get_file_list(path, "*.17_TEC", file_list, nfiles)
  if (nfiles == 0) then
     print *, "No matching files found in:", trim(path)
     stop
  end if
  print *, "Files found:", nfiles

  ! ---------------------------------------------------------
  ! 2. Create NetCDF file
  ! ---------------------------------------------------------
  ierr = nf90_create(trim(path)//"ROTI_output.nc", NF90_CLOBBER, ncid)
  call check_err(ierr, "create")

  ierr = nf90_def_dim(ncid, "file", nfiles, dim_file)
  ierr = nf90_def_dim(ncid, "x", 2880, dim_x)
  ierr = nf90_def_dim(ncid, "y", 32, dim_y)
  ierr = nf90_def_dim(ncid, "loc", 3, dim_loc)

  ierr = nf90_def_var(ncid, "filename", NF90_CHAR, (/dim_file/), var_name)
  ierr = nf90_def_var(ncid, "rec_loc", NF90_DOUBLE, (/dim_loc, dim_file/), var_loc)
  ierr = nf90_def_var(ncid, "roti_arr", NF90_DOUBLE, (/dim_x, dim_y, dim_file/), var_roti)

  ierr = nf90_enddef(ncid)

  ! ---------------------------------------------------------
  ! 3. Process each file
  ! ---------------------------------------------------------
  dt = 30.0_dp_local
  do i = 1, nfiles
     fname_full = trim(file_list(i))

     ! Extract only filename (no path) for NetCDF metadata
     p = index(fname_full, '/', back=.true.)
     if (p > 0) then
        fname_only = fname_full(p+1:)
     else
        fname_only = fname_full
     end if

     print *, "Processing:", trim(fname_full)

     ! Pass full path to the module
     call read_rinex2(trim(fname_full), rec_loc, tec_arr)
     call compute_roti(tec_arr, dt, roti_arr)

     ! Write metadata + data to NetCDF
     ierr = nf90_put_var(ncid, var_name, trim(fname_only), start=(/i/))
     ierr = nf90_put_var(ncid, var_loc, rec_loc, start=(/1, i/))
     ierr = nf90_put_var(ncid, var_roti, roti_arr, start=(/1, 1, i/))
  end do

  ! ---------------------------------------------------------
  ! 4. Close NetCDF
  ! ---------------------------------------------------------
  ierr = nf90_close(ncid)
  call check_err(ierr, "close")
  print *, "ROTI data saved to:", trim(path)//"ROTI_output.nc"

contains

  ! ---------------------------------------------------------
  ! Get file list using bash-compatible ls
  ! ---------------------------------------------------------
  subroutine get_file_list(path, pattern, files, n)
    character(len=*), intent(in) :: path, pattern
    character(len=*), intent(out) :: files(:)
    integer, intent(out) :: n
    integer :: ios, unit
    character(len=1024) :: cmd, line

    cmd = "ls -1 " // trim(path)//trim(pattern)//" > filelist.txt"
    call execute_command_line(cmd)

    n = 0
    open(newunit=unit, file="filelist.txt", status="old", action="read", iostat=ios)
    if (ios /= 0) then
       print *, "Error: could not open filelist.txt"
       return
    end if

    do
       read(unit,'(A)', iostat=ios) line
       if (ios /= 0) exit
       n = n + 1
       files(n) = adjustl(trim(line))
    end do
    close(unit)
    call execute_command_line("rm -f filelist.txt")
  end subroutine get_file_list

  ! ---------------------------------------------------------
  ! NetCDF error checker
  ! ---------------------------------------------------------
  subroutine check_err(status, msg)
    integer, intent(in) :: status
    character(len=*), intent(in) :: msg
    if (status /= NF90_NOERR) then
       print *, "NetCDF error in", trim(msg), ":", trim(nf90_strerror(status))
       stop
    end if
  end subroutine check_err

end program batch_roti_to_netcdf

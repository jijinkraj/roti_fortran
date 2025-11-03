module rinex_module
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
contains

  !-------------------[ getnum ]-------------------
  subroutine getnum(str, maxcount, nums, count)
    implicit none
    character(len=*), intent(in) :: str
    integer, intent(in) :: maxcount
    integer, allocatable, intent(out) :: nums(:)
    integer, intent(out) :: count
    integer :: i
    character(len=1) :: c
    character(len=256) :: temp

    allocate(nums(maxcount))
    count = 0
    temp = ''

    do i = 1, len_trim(str)
        c = str(i:i)
        if (c >= '0' .and. c <= '9') then
            temp = trim(temp)//c
        else
            if (len_trim(temp) > 0) then
                count = count + 1
                if (count > size(nums)) return
                read(temp,*) nums(count)
                temp = ''
            end if
        end if
    end do

    if (len_trim(temp) > 0) then
        count = count + 1
        if (count <= size(nums)) read(temp,*) nums(count)
    end if
  end subroutine getnum


  !-------------------[ read_rinex2 ]-------------------
  subroutine read_rinex2(filename, rec_loc, tec_arr)
    use, intrinsic :: ieee_arithmetic
    implicit none
    character(len=*), intent(in)        :: filename
    real(dp), intent(out)               :: rec_loc(3)
    real(dp), allocatable, intent(out)  :: tec_arr(:,:)

    integer, parameter :: maxepoch = 2880
    integer, parameter :: maxsat   = 32

    integer :: ios, pos, row
    integer :: num_count, nprn, i, prn_id, yy, mo, dd, hr, mn
    real(dp) :: ss, tec_value
    integer :: flag
    real(dp) :: NaN
    character(len=256) :: line, ssstr
    integer, allocatable :: prns(:)

    allocate(tec_arr(maxepoch, maxsat))
    NaN = ieee_value(0.0_dp, ieee_quiet_nan)
    tec_arr = NaN
    rec_loc = NaN

    open(unit=10, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error: cannot open file.'
        return
    end if

    do
        read(10,'(A)',iostat=ios) line
        if (ios /= 0) then
            print *, "Error: can't read or reached EOF"
            close(10);  return
        end if
        pos = index(line,'APPROX POSITION XYZ')
        if (pos > 0) then
            read(line,*) rec_loc(1), rec_loc(2), rec_loc(3)
        end if
        if (index(line,'END OF HEADER') > 0) exit
    end do

    do
        read(10,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle

        if (count([(verify(line(i:i),'0123456789')==0, i=1,len_trim(line))]) < 3) cycle

        read(line,*,err=100) yy,mo,dd,hr,mn,ss
        row = (hr*120 + mn*2 + int(ss)/30)+1
        if (row < 1 .or. row > maxepoch) cycle

        if (len_trim(line) > 31) then
            ssstr = line(31:)
        else
            ssstr = ''
        end if

        call getnum(ssstr, maxsat, prns, num_count)

        if (num_count > 0) then
        else
          print *, 'No numbers found in ssstr.'
        end if
        
        nprn = prns(1)
        do i = 1, nprn
            read(10,*,iostat=ios) tec_value, flag
            if (ios /= 0) exit
            prn_id = prns(i+1)
            if (prn_id < 1 .or. prn_id > maxsat) cycle
            if (flag > 0) then
                tec_arr(row, prn_id) = NaN
            else
                tec_arr(row, prn_id) = tec_value
            end if
        end do

        if (allocated(prns)) deallocate(prns)
        100 continue
    end do

    close(10)
  end subroutine read_rinex2


  !-------------------[ compute_roti ]-------------------
  subroutine compute_roti(tec_arr, dt, roti_arr)
    use, intrinsic :: ieee_arithmetic
    implicit none
    real(dp), intent(in)  :: tec_arr(:,:)
    real(dp), intent(in)  :: dt
    real(dp), allocatable, intent(out) :: roti_arr(:,:)
    integer :: nrow, ncol, i, j, win
    real(dp), allocatable :: rot(:,:)
    real(dp) :: mean1, mean2, NaN

    NaN = ieee_value(0.0_dp, ieee_quiet_nan)
    nrow = size(tec_arr,1); ncol = size(tec_arr,2)
    allocate(rot(nrow-1,ncol), roti_arr(nrow,ncol))
    roti_arr = NaN

    ! Compute ROT
    do j = 1, ncol
      do i = 1, nrow - 1
        if (.not.ieee_is_nan(tec_arr(i,j)) .and. .not.ieee_is_nan(tec_arr(i+1,j))) then
          rot(i,j) = (tec_arr(i+1,j) - tec_arr(i,j)) / dt
        else
          rot(i,j) = NaN
        end if
      end do
    end do

    ! Compute ROTI using 10-point window
    win = 10
    do j = 1, ncol
      do i = win, nrow - 1
        if (count(.not.ieee_is_nan(rot(i-win+1:i,j))) == win) then
          mean1 = sum(rot(i-win+1:i,j)) / win
          mean2 = sum(rot(i-win+1:i,j)**2) / win
          roti_arr(i,j) = sqrt(max(0.0_dp, mean2 - mean1**2))
        end if
      end do
    end do
  end subroutine compute_roti

end module rinex_module

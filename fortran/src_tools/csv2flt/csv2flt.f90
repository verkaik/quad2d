program create_hdr_tiles
  ! modules
  use utilsmod, only: MXSLEN, I4B, R8B, tCsv, errmsg, logmsg, ta
  use hdrModule, only: writeflt
  !
  implicit none
  !
! -- locals
  type(tCSV), pointer :: csv => null()
  character(len=MXSLEN) :: f_csv, f_flt, s, ic_s, ir_s, iv_s
  integer(I4B) :: na, nc, nr, ic, ir, i4mv, i, iv, nmax
  integer(I4B), dimension(:), allocatable :: ic_arr, ir_arr, iv_arr, cnt
  integer(I4B), dimension(:,:), allocatable :: xi4
  real(R8B) :: xll, yll, cs
! ------------------------------------------------------------------------------
  !
  na = nargs() - 1
  if (na < 11) then
    call errmsg('Number of arguments must be at least 11.')
  end if
  call getarg(1,f_csv)
  call getarg(2,ic_s)
  call getarg(3,ir_s)
  call getarg(4,iv_s)
  call getarg(5,f_flt)
  call getarg(6,s); read(s,*) nc
  call getarg(7,s); read(s,*) nr
  call getarg(8,s); read(s,*) xll
  call getarg(9,s); read(s,*) yll
  call getarg(10,s); read(s,*) cs
  call getarg(11,s); read(s,*) i4mv
  !
  allocate(csv)
  call csv%read(f_csv)
  !
  allocate(xi4(nc,nr)); xi4 = i4mv
  call csv%get_column(key=ic_s, i4a=ic_arr)
  call csv%get_column(key=ir_s, i4a=ir_arr)
  call csv%get_column(key=iv_s, i4a=iv_arr)
  ! 
  nmax = maxval(iv_arr)
  allocate(cnt(nmax)); cnt = 0
  !
  do i = 1, size(ic_arr)
    ic = ic_arr(i); ir = ir_arr(i); iv = iv_arr(i)
    xi4(ic,ir) = iv
    cnt(iv) = 1
  end do
  !
  call logmsg('# unique values: '//ta([sum(cnt)]))
  
  call writeflt(f_flt, xi4, nc, nr, xll, yll, cs, i4mv)
  !
end program

program create_hdr_tiles
  ! modules
  use utilsmod, only: MXSLEN, I4B, R4B, logmsg, errmsg, ta, open_file, tBB, tBBX
  use hdrModule, only: tHdr, tHdrHdr, i_flt, i_env, writeflt
  use ieee_arithmetic, only: ieee_is_nan
  !
  implicit none
  !
! -- locals
  type(tBb) :: bbi
  type(tBbx) :: bbx
  type(tHdr), pointer :: hdrg => null()
  type(tHdrHdr), pointer :: hdr => null()
  logical :: write_envi, write_tile_index, ok, write_tilenumber
  integer(I4B) :: na, ntile, nr_blk, nc_blk, nr, nc, gnc, gnr, itile
  integer(I4B) :: iu, irb, icb, ir, ic, ir0, ir1, ic0, ic1, i_file_type, mv_envi
  integer(I4B), dimension(:,:), allocatable :: tilenumber
  character(len=MXSLEN) :: s, fp_in, fp_out, f
! ------------------------------------------------------------------------------
  !
  ! read the input file
  na = nargs() - 1
  if (na < 5) then
    call errmsg('Number of arguments must be at least 3.')
  end if
  call getarg(1,fp_in)
  call getarg(2,fp_out)
  call getarg(3,s); read(s,*) nr_blk
  call getarg(4,s); read(s,*) nc_blk
  if ((nr_blk < 0).or.(nc_blk < 0)) then
    write_tilenumber = .true.
  else
    write_tilenumber = .false.
  end if
  nr_blk = abs(nr_blk); nc_blk = abs(nc_blk)
  !
  call getarg(5,s)
  if (trim(s) == 'envi') then
    write_envi = .true.
    if (na == 6) then
      call getarg(6,s); read(s,*) mv_envi
    else
      call errmsg('ENVI missing value not set.')
    end if
  else
    write_envi = .false.
  end if
  !
  if (na == 6) then
    write_tile_index = .true.
  else
    write_tile_index = .false.
  end if
  !
  ! read the flt header
  allocate(hdr); call hdr%read(fp_in)
  gnc = hdr%ncol; gnr = hdr%nrow
  call logmsg('Processing headered file '//trim(fp_in)//' of size: '// &
    ta((/gnc/))//' cols x '//ta((/gnr/))//' rows')
  nc = gnc/nc_blk; nr = gnr/nr_blk
  !
  ntile = nr_blk*nc_blk
  call logmsg('Writing '//ta((/ntile/))//' tiles of size: '//ta((/nc/))// &
    ' cols x '// ta((/nr/))//' rows')
  !
  if (write_tilenumber) then
    allocate(tilenumber(gnc,gnr)); tilenumber = 0
    call hdr%get_bbx(bbx)
  end if
  !
  itile = 0
  do irb = 1, nr_blk ! loop over blocks
    do icb = 1, nc_blk ! loop over blocks
      allocate(hdrg)
      call hdrg%init()
      call logmsg('... Reading block ('//ta((/icb,irb/))//')')
      ir0 = (irb-1)*nr + 1; ir1 = ir0 + nr -1
      ic0 = (icb-1)*nc + 1; ic1 = ic0 + nc -1
      bbi%ic0 = ic0; bbi%ic1 = ic1; bbi%ir0 = ir0; bbi%ir1 = ir1 
      call hdrg%read_clip_grid(fp_in, bbi) ! read block
      !
      call logmsg('... Replacing NaN values with nodata.')
      ok = .false.
      if (allocated(hdrg%dat%xi4)) then
        ok = .true.
        do ir = 1, size(hdrg%dat%xi4,2); do ic = 1, size(hdrg%dat%xi4,1)
          if (ieee_is_nan(real(hdrg%dat%xi4(ic,ir),R4B))) then
            hdrg%dat%xi4(ic,ir) = hdrg%hdr%mvi4
          end if
        end do; end do
      end if
      if (allocated(hdrg%dat%xr4)) then
        ok = .true.
        do ir = 1, size(hdrg%dat%xr4,2); do ic = 1, size(hdrg%dat%xr4,1)
          if (ieee_is_nan(hdrg%dat%xr4(ic,ir))) then
            hdrg%dat%xr4(ic,ir) = hdrg%hdr%mvr4
          end if
        end do; end do
      end if
      if (allocated(hdrg%dat%xr8)) then
        ok = .true.
        do ir = 1, size(hdrg%dat%xr8,2); do ic = 1, size(hdrg%dat%xr8,1)
          if (ieee_is_nan(hdrg%dat%xr8(ic,ir))) then
            hdrg%dat%xr8(ic,ir) = hdrg%hdr%mvr8
          end if
        end do; end do
      end if
      if (.not.ok) then
        call errmsg('Unsupported file type.')
      end if
      !
      call hdrg%hdr%clip(ic0, ic1, ir0, ir1)
      i_file_type = hdrg%get_file_type()
      if (i_file_type == i_env) then
        call hdrg%set_mv(mvi4=mv_envi)
      end if
      call logmsg('... Checking for empty blocks.')
      if (.not.hdrg%check_empty()) then
        itile = itile + 1
        if (write_tile_index) then
          f = trim(fp_out)//'_itile_r'//ta((/irb/),'(i3.3)')//'_c'//ta((/icb/),'(i3.3)')
          hdrg%dat%xi4 = itile
        else
          f = trim(fp_out)//'_r'//ta((/irb/),'(i3.3)')//'_c'//ta((/icb/),'(i3.3)')
        end if
        !
        if (write_envi) then
          call hdrg%write(f, file_type='envi')
        else
          call hdrg%write(f)
        end if
        !
        if (write_tilenumber) then
          do ir = ir0, ir1; do ic = ic0, ic1
            tilenumber(ic,ir) = itile
          end do; end do
        end if
      else
        call logmsg('*** block is empty! ****')
      end if
      call hdrg%clean()
      deallocate(hdrg)
    end do
  end do
  !
  if (write_tilenumber) then
    f = trim(fp_out)//'_tilenumber'
    call writeflt(f, tilenumber, gnc, gnr, bbx%xll, bbx%yll, bbx%cs, 0)
  end if
  !
  close(iu)
  !
end program

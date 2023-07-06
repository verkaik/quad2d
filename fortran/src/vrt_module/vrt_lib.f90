module vrt_module
  ! modules
  use utilsmod, only: I1B, I2B, I4B, I8B, R4B, R8B, &
    i_i1, i_i2, i_i4, i_r4, i_r8, I4ZERO, R4ZERO, R8ZERO, R8HALF, &
    logmsg, errmsg, tBb, tBbX, tBbObj, MXSLEN, open_file, read_line, &
    renumber, bbi_intersect, bbx_intersect, base_name, get_xy, get_icr,&
    ta, swap_slash, point_in_bb, strip_ext, get_ext, change_case
  use hdrModule, only: tHdrHdr, tHdr, i_uscl_nodata, i_dscl_nodata

  implicit none
  
  private
  
  ! begin part
  integer(I4B), parameter :: i_VRTDataset_beg    =  1
  integer(I4B), parameter :: i_GeoTransform      =  2
  integer(I4B), parameter :: i_VRTRasterBand_beg =  3
  integer(I4B), parameter :: i_NoDataValue       =  4
  integer(I4B), parameter :: n_raw_beg = i_NoDataValue
  !
  ! middle part for each tile
  integer(I4B), parameter :: i_ComplexSource_beg = 1
  integer(I4B), parameter :: i_SourceFilename    = 2
  integer(I4B), parameter :: i_SourceBand        = 3
  integer(I4B), parameter :: i_SourceProperties  = 4
  integer(I4B), parameter :: i_SrcRect           = 5
  integer(I4B), parameter :: i_DstRect           = 6
  integer(I4B), parameter :: i_NODATA            = 7
  integer(I4B), parameter :: i_ComplexSource_end = 8
  integer(I4B), parameter :: n_raw_mid = i_ComplexSource_end
  
  ! end part
  integer(I4B), parameter :: i_VRTRasterBand_end =  1
  integer(I4B), parameter :: i_VRTDataset_end    =  2
  integer(I4B), parameter :: n_raw_end = i_VRTDataset_end
  
  type :: tVrtRawLine
    integer(I4B) :: ni = 0 ! number of indent
    integer(I4B) :: ns = 0 ! number of string items
    integer(I4B) :: nv = 0 ! number of values to set
    integer(I4B), dimension(:), pointer :: iv     => null()
    logical,      dimension(:), pointer :: iv_set => null()
    character(MXSLEN), dimension(:), pointer :: s => null()
  contains
    procedure :: init       => tVrtRawLine_init
    procedure :: tVrtRawLine_set_by_i
    procedure :: tVrtRawLine_set_by_s
    generic   :: set        => tVrtRawLine_set_by_i, &
                               tVrtRawLine_set_by_s
    procedure :: all_set    => tVrtRawLine_all_set
    procedure :: get        => tVrtRawLine_get_by_i
    procedure :: get_string => tVrtRawLine_get_string
    procedure :: clean      => tVrtRawLine_clean
  end type tVrtRawLine
  !
  type :: tVrtRawLineSet
    type(tVrtRawLine), dimension(:), pointer :: raw => null()
  end type tVrtRawLineSet
  !
  type :: tVrtRaw
    type(tVrtRawLine),    dimension(:), pointer :: raw_beg => null()
    type(tVrtRawLineSet), dimension(:), pointer :: raw_mid => null()
    type(tVrtRawLine),    dimension(:), pointer :: raw_end => null()
  contains
    procedure :: init  => tVrtRaw_init
    procedure :: clean => tVrtRaw_clean
    procedure :: read  => tVrtRaw_read
    procedure :: write => tVrtRaw_write
  end type tVrtRaw
  !
  type :: tVrtTile
    type(tHdr), pointer :: hdrg => null()
    type(tBbX) :: src_bbx
    type(tBb)  :: src_bbi
    type(tBbX) :: dst_bbx
    type(tBb)  :: dst_bbi
    character(len=MXSLEN) :: file_ext
    character(len=MXSLEN) :: file_name
    logical :: active
  contains
    procedure :: init  => tVrtTile_init
    procedure :: clean => tVrtTile_clean
    procedure :: set   => tVrtTile_set
    procedure :: get   => tVrtTile_get
  end type tVrtTile
  !
  type, public :: tVrt
    integer(I4B) :: iu = 0
    type(tVrtRaw), pointer :: raw => null()
    character(len=MXSLEN) :: f
    type(tBbX) :: dst_bbx
    type(tBb)  :: dst_bbi
    integer(I4B) :: ntiles
    logical ::overlap
    logical :: full_data_read
    type(tVrtTile), dimension(:), pointer :: tiles => null()
  contains
    procedure :: init                  => tVrt_init
    procedure :: clean                 => tVrt_clean
    procedure :: clean_x               => tVrt_clean_x
    procedure :: overlapping_tiles     => tVrt_overlapping_tiles
    procedure :: read_ntiles           => tVrt_read_ntiles
    procedure :: get_ntiles            => tVrt_get_ntiles
    procedure :: get_bb                => tVrt_get_bb
    procedure :: get_tile_file         => tVrt_get_tile_file
    procedure :: init_tile_raster      => tVrt_init_tile_raster
    procedure :: read_full_tile        => tVrt_read_full_tile
    procedure :: read_extent           => tVrt_read_extent
    procedure :: read_extent_native_cs => tVrt_read_extent_native_cs
    procedure :: get_act_tile          => tVrt_get_act_tile
    procedure :: read_xy               => tVrt_read_xy
    procedure :: init_write            => tVrt_init_write
    procedure :: set_nodata            => tVrt_set_nodata
    procedure :: write                 => tVrt_write
    procedure :: buffer_data           => tVrt_buffer_data
  end type tVrt
  !
  type, public :: tVrtArray
    character(len=MXSLEN) :: f
    integer(I4B) :: n = 0
    integer(I4B), dimension(:), allocatable :: iconst
    real(R4B), dimension(:), allocatable :: r4const
    type(tVrt), dimension(:), pointer :: vrta => null()
  contains
    procedure :: init           => tVrtArray_init
    procedure :: clean          => tVrtArray_clean
    procedure :: check_constant => tVrtArray_check_constant
    procedure :: get_vrt        => tVrtArray_get_vrt
  end type tVrtArray
  !
  contains
  !
! ==============================================================================
! ==============================================================================
! Type: tVrtArray
! ==============================================================================
! ==============================================================================
  
  subroutine tVrtArray_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtArray) :: this
    ! -- local
    type(tVrt), pointer :: vrt => null()
    integer(I4B), parameter :: MAXLINE = 1000
    character(len=MXSLEN), dimension(MAXLINE) :: sa
    character(len=MXSLEN) :: s, f
    logical :: flag
    integer(I4B) :: i, iu, ios
! ------------------------------------------------------------------------------
    !
    f = this%f
    if (len_trim(f) == 0) then
      call errmsg('tVrtArray_init: file not present.')
    end if
    !
    call this%clean()
    !
    this%n = 0
    !
    call open_file(f, iu, 'r')
    do while(.true.)
      read(unit=iu,iostat=ios,fmt='(a)') s
      if (ios /= 0) then
        exit
      end if
      if (len_trim(s) == 0) cycle
      !
      this%n = this%n + 1
      if (this%n > MAXLINE) then
        call errmsg('tVrtArray_init: increase MAXLINE.')
      else
        sa(this%n) = trim(s)
      end if
    end do
    close(iu)
    !
    if (this%n == 0) then
      call logmsg('tVrtArray_init: nothing to do, returning...')
      return
    end if
    !
    allocate(this%vrta(this%n))
    allocate(this%iconst(this%n)); this%iconst = 0
    allocate(this%r4const(this%n)); this%r4const = R4ZERO
    
    do i = 1, this%n
      ! check for constant
      s = adjustl(change_case(sa(i), 'l'))
      if (s(1:8) == 'constant') then
        this%iconst(i) = 1
        read(sa(i),*) s, this%r4const(i)
      else
        vrt => this%vrta(i)
        call vrt%init(sa(i))
      end if
    end do
    !
    return
  end subroutine tVrtArray_init
  
  subroutine tVrtArray_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtArray) :: this
    ! -- local
    type(tVrt), pointer :: vrt => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    if (associated(this%vrta)) then
      do i = 1, this%n
        vrt => this%vrta(i)
        call vrt%clean()
      end do
      deallocate(this%vrta); this%vrta => null()
    end if
    if (allocated(this%iconst)) deallocate(this%iconst)
    if (allocated(this%r4const)) deallocate(this%r4const)
    !
    return
  end subroutine tVrtArray_clean
  
  subroutine tVrtArray_check_constant(this, i, const, r4const)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtArray) :: this
    integer(I4B), intent(in) :: i
    logical, intent(out) :: const
    real(R4B), intent(out) :: r4const
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (this%iconst(i) == 1) then
      const   = .true.
      r4const = this%r4const(i)
    else
      const   = .false.
      r4const = R4ZERO
    end if
    !
    return
  end subroutine tVrtArray_check_constant
  
  function tVrtArray_get_vrt(this, i) result(vrt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtArray) :: this
    integer(I4B), intent(in) :: i
    type(tVrt), pointer :: vrt
    ! -- local
! ------------------------------------------------------------------------------
    !
    vrt => null()
    if ((i <= 0).or.(i > this%n)) then
      call errmsg('tVrtArray_get_vrt: index out of range.')
    end if
    vrt => this%vrta(i)
    !
    return
  end function tVrtArray_get_vrt
  
! ==============================================================================
! ==============================================================================
! Type: tVrtRawLine
! ==============================================================================
! ==============================================================================
  
  subroutine tVrtRawLine_init(this, ni, sa)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRawLine) :: this
    integer(I4B), intent(in) :: ni
    character(len=*), dimension(:), intent(in) :: sa
    ! -- local
    integer(I4B) :: i, ns, nv
! ------------------------------------------------------------------------------
    !
    this%ni = ni
    
    ns = size(sa)
    if (ns == 0) then
      call errmsg('tVrtRawLine_init: ns = 0')
    end if
    !
    this%ns = ns
    !
    ! first, count the values to set
    nv = 0
    do i = 1, ns
      if (len_trim(sa(i)) == 0) then
        nv = nv + 1
      end if
    end do
    this%nv = nv
    !
    ! allocate and store
    allocate(this%s(this%ns))
    if (this%nv > 0) then
      allocate(this%iv(this%nv), this%iv_set(this%nv))
      this%iv = 0; this%iv_set = .false.
    end if
    !
    nv = 0
    do i = 1, ns
      if (len_trim(sa(i)) == 0) then
        nv = nv + 1
        this%iv(nv) = i
      end if
      this%s(i) = trim(adjustl(sa(i)))
    end do
    !
    return
  end subroutine tVrtRawLine_init
  
  subroutine tVrtRawLine_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRawLine) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (associated(this%iv))     deallocate(this%iv)
    if (associated(this%iv_set)) deallocate(this%iv_set)
    if (associated(this%s))      deallocate(this%s)
    !
    this%ni = 0
    this%ns = 0
    this%nv = 0
    this%iv     => null()
    this%iv_set => null()
    this%s      => null()
    !
    return
  end subroutine tVrtRawLine_clean
    
  subroutine tVrtRawLine_set_by_i(this, i, sv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRawLine) :: this
    integer(I4B), intent(in) :: i
    character(len=*), intent(in) :: sv
    ! -- local
    integer(I4B) :: iv
! ------------------------------------------------------------------------------
    !
    iv = this%iv(i)
    this%s(iv) = trim(sv)
    this%iv_set(i) = .true.
    !
    return
  end subroutine tVrtRawLine_set_by_i

  subroutine tVrtRawLine_get_by_i(this, iv, i1v, i2v, i4v, r4v, r8v, r8a, sv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRawLine) :: this
    integer(I4B), intent(in) :: iv
    integer(I1B),                        intent(out),   optional :: i1v
    integer(I2B),                        intent(out),   optional :: i2v
    integer(I4B),                        intent(out),   optional :: i4v
    real(R4B),                           intent(out),   optional :: r4v
    real(R8B),                           intent(out),   optional :: r8v
    real(R8B),             dimension(:), intent(inout), optional :: r8a
    character(len=MXSLEN),               intent(out),   optional :: sv
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    i = this%iv(iv)
    if (present(i1v)) read(this%s(i),*) i1v
    if (present(i2v)) read(this%s(i),*) i2v
    if (present(i4v)) read(this%s(i),*) i4v
    if (present(r4v)) read(this%s(i),*) r4v
    if (present(r8v)) read(this%s(i),*) r8v
    if (present(r8a)) read(this%s(i),*) r8a
    if (present(sv)) sv = this%s(i)
    !
    return
  end subroutine tVrtRawLine_get_by_i
  
  subroutine tVrtRawLine_set_by_s(this, s)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRawLine) :: this
    character(len=*), intent(in) :: s
    ! -- local
    integer(I4B) :: i, iv, i1, i2
! ------------------------------------------------------------------------------
    !
    do iv = 1, this%nv
      i = this%iv(iv)
      i1 = index(s,trim(this%s(i-1)))
      if (i1 == 0) then
        call logmsg('"'//trim(s)//'"')
        call errmsg('tVrtRawLine_set_by_s: could not parse string.')
      else
        i1 = i1 + len_trim(this%s(i-1))
      end if
      i2 = index(s,trim(this%s(i+1)))
      if (i2 == 0) then
        call logmsg('"'//trim(s)//'"')
        call errmsg('tVrtRawLine_set_by_s: could not parse string.')
      else
        i2 = i2 - 1
      end if
      this%s(i) = s(i1:i2)
      this%iv_set(iv) = .true.
    end do
    !
    return
  end subroutine tVrtRawLine_set_by_s
  
  function tVrtRawLine_all_set(this) result(all_set)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRawLine) :: this
    logical :: all_set
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    all_set = .true.
    do i = 1, this%nv
      if (this%iv_set(i) == .false.) then
        all_set = .false.
        return
      end if
    end do
    !
    return
  end function tVrtRawLine_all_set
  
  function tVrtRawLine_get_string(this) result(s)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRawLine) :: this
    character(len=MXSLEN) :: s
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    s = ''
    do i = 1, this%ni
      s = trim(s)//'$$'
    end do
    do i = 1, this%ns
      s = trim(s)//trim(this%s(i))
    end do
    !
    do i = 1, len_trim(s)
      if (s(i:i) == '$') then
        s(i:i) = ' '
      end if
    end do
    !
    return
  end function tVrtRawLine_get_string
    
! ==============================================================================
! ==============================================================================
! Type: tVrtRaw
! ==============================================================================
! ==============================================================================
  
  subroutine tVrtRaw_init(this, ntiles, source)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRaw) :: this
    integer(I4B), intent(in) :: ntiles
    character(len=*) :: source
    ! -- local
    type(tVrtRawLine), pointer :: vl => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    ! allocate 
    allocate(this%raw_beg(n_raw_beg))
    allocate(this%raw_mid(ntiles))
    do i = 1, ntiles
      allocate(this%raw_mid(i)%raw(n_raw_mid))
    end do
    allocate(this%raw_end(n_raw_end))
    !
    vl => this%raw_beg(i_VRTDataset_beg)
    call vl%init(0, (/'<VRTDataset rasterXSize="','', '" rasterYSize="','', &
      '">'/))

    vl => this%raw_beg(i_GeoTransform)
    call vl%init(1,(/'<GeoTransform> ','','</GeoTransform>'/))
    
    vl => this%raw_beg(i_VRTRasterBand_beg)
    call vl%init(1,(/'<VRTRasterBand dataType="', '','" band="','','">'/))
    call vl%set(2,'1') ! DEFAULT
    
    vl => this%raw_beg(i_NoDataValue)
    call vl%init(2,(/'<NoDataValue>','','</NoDataValue>'/))
    
    do i = 1, ntiles
      vl => this%raw_mid(i)%raw(i_ComplexSource_beg)
      !call vl%init(2,(/'<ComplexSource>'/))
      call vl%init(2,(/'<'//trim(source)//'>'/))

      vl => this%raw_mid(i)%raw(i_SourceFilename)
      call vl%init(3,(/'<SourceFilename relativeToVRT="','','">','',&
        '</SourceFilename>'/))
      call vl%set(1,'0') ! DEFAULT

      vl => this%raw_mid(i)%raw(i_SourceBand)
      call vl%init(3,(/'<SourceBand>','','</SourceBand>'/))
      call vl%set(1,'1') ! DEFAULT
      
      vl => this%raw_mid(i)%raw(i_SourceProperties)
      call vl%init(3,(/'<SourceProperties RasterXSize="','','" RasterYSize="',&
        '','" DataType="','','" BlockXSize="','','" BlockYSize="','','" />'/))
      
      vl => this%raw_mid(i)%raw(i_SrcRect)
      call vl%init(3,(/'<SrcRect xOff="','','" yOff="','','" xSize="','',&
        '" ySize="','','" />'/))
      
      vl => this%raw_mid(i)%raw(i_DstRect)
      call vl%init(3,(/'<DstRect xOff="','','" yOff="','','" xSize="','',&
        '" ySize="','','" />'/))
      
      vl => this%raw_mid(i)%raw(i_NODATA)
      call vl%init(3,(/'<NODATA>','','</NODATA>'/))
      
      vl => this%raw_mid(i)%raw(i_ComplexSource_end)
      !call vl%init(2,(/'</ComplexSource>'/))
      call vl%init(2,(/'</'//trim(source)//'>'/))
    end do
    
    vl => this%raw_end(i_VRTRasterBand_end)
    call vl%init(1,(/'</VRTRasterBand>'/))

    vl => this%raw_end(i_VRTDataset_end)
    call vl%init(1,(/'</VRTDataset>'/))
    !
    return
  end subroutine tVrtRaw_init
  
  subroutine tVrtRaw_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRaw) :: this
    ! -- local
    type(tVrtRawLine), pointer :: vl => null()
    integer(I4B) :: i, j
! ------------------------------------------------------------------------------
    if (associated(this%raw_beg)) then
      do i = 1, n_raw_beg
        vl => this%raw_beg(i)
        call vl%clean()
      end do
      deallocate(this%raw_beg)
    end if
    if (associated(this%raw_mid)) then
      do j = 1, size(this%raw_mid)
        do i = 1, n_raw_mid
          vl => this%raw_mid(j)%raw(i)
          call vl%clean()
        end do
        deallocate(this%raw_mid(j)%raw)
      end do
      deallocate(this%raw_mid)
    end if
    if (associated(this%raw_end)) then
      do i = 1, n_raw_end
        vl => this%raw_end(i)
        call vl%clean()
      end do
      deallocate(this%raw_end)
    end if
    !
    this%raw_beg => null()
    this%raw_mid => null()
    this%raw_end => null()
    !
    return
  end subroutine tVrtRaw_clean
    
  subroutine tVrtRaw_read(this, iu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRaw) :: this
    integer(I4B), intent(in) :: iu
    ! -- local
    type(tVrtRawLine), pointer :: vl => null()
    logical :: found, first
    integer(I4B) :: i, j, ios
    character(len=MXSLEN) :: s
! ------------------------------------------------------------------------------
    !
    ! rewind
    rewind(iu)
    
    do i = 1, n_raw_beg
      vl => this%raw_beg(i)
      found =.false.
      do while(.not.found)
        call read_line(iu, s, ios)
        if (ios /= 0) exit
        if (index(s,trim(vl%s(1))) > 0) found = .true.
      end do
      if (found) then
        call vl%set(s)
        if (.not.vl%all_set()) then
          call logmsg('"'//trim(s)//'"')
          call errmsg('Could not set all values.')
        end if
      else
        call logmsg('Could not find: '//trim(vl%s(1)))
      end if
    end do
    first = .true.
    do j = 1, size(this%raw_mid)
      do i = 1, n_raw_mid
        vl => this%raw_mid(j)%raw(i)
        if (first) then
          found =.false.
          do while(.not.found)
            call read_line(iu, s, ios)
            if (ios /= 0) exit
            if (index(s,trim(vl%s(1))) > 0) found = .true.
          end do
          first = .false.
        else
          call read_line(iu, s, ios)
        end if
        if (found) then
          call vl%set(s)
          if (.not.vl%all_set()) then
            call logmsg('"'//trim(s)//'"')
            call errmsg('Could not set all values.')
          end if
        else
          call logmsg('Could not find: '//trim(vl%s(1)))
        end if
      end do
    end do
    do i = 1, n_raw_end
      vl => this%raw_end(i)
      found =.false.
      do while(.not.found)
        call read_line(iu, s, ios)
        if (ios /= 0) exit
        if (index(s,trim(vl%s(1))) > 0) found = .true.
      end do
      call vl%set(s)
      if (.not.vl%all_set()) then
        call logmsg('"'//trim(s)//'"')
        call errmsg('Could not set all values.')
      end if
    end do
    !
    return
  end subroutine tVrtRaw_read

  subroutine tVrtRaw_write(this, iu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtRaw) :: this
    integer(I4B), intent(in) :: iu
    ! -- local
    type(tVrtRawLine), pointer :: vl => null()
    integer(I4B) :: i, j
    character(len=MXSLEN) :: s
! ------------------------------------------------------------------------------
    !
    ! rewind
    rewind(iu)
    
    do i = 1, n_raw_beg
      vl => this%raw_beg(i)
      s = vl%get_string()
      write(iu,'(a)') trim(s)
    end do
    do j = 1, size(this%raw_mid)
      do i = 1, n_raw_mid
        vl => this%raw_mid(j)%raw(i)
        s = vl%get_string()
        write(iu,'(a)') trim(s)
      end do
    end do
    do i = 1, n_raw_end
      vl => this%raw_end(i)
      s = vl%get_string()
      write(iu,'(a)') trim(s)
    end do
    !
    return
  end subroutine tVrtRaw_write
  
! ==============================================================================
! ==============================================================================
! Type: tVrtTile
! ==============================================================================
! ==============================================================================
  
  subroutine tVrtTile_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtTile) :: this
 ! ------------------------------------------------------------------------------
    !
    if (.not.associated(this%hdrg)) then
      allocate(this%hdrg)
    end if
    call this%src_bbx%init()
    call this%src_bbi%init()
    call this%dst_bbx%init()
    call this%dst_bbx%init()
    this%file_ext = ''
    this%file_name = ''
    this%active = .true.
    !
    return
  end subroutine tVrtTile_init
  
  subroutine tVrtTile_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtTile) :: this
 ! ------------------------------------------------------------------------------
    !
    if (associated(this%hdrg)) then
      call this%hdrg%clean()
    end if
    call this%init()
    !
    return
  end subroutine tVrtTile_clean
  
  subroutine tVrtTile_set(this, src_bbx, src_bbi, dst_bbx, dst_bbi, &
    file_ext, file_name, active)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtTile) :: this
    type(tBbX), intent(in), optional :: src_bbx
    type(tBb),  intent(in), optional :: src_bbi
    type(tBbX), intent(in), optional :: dst_bbx
    type(tBb),  intent(in), optional :: dst_bbi
    character(len=*), intent(in), optional :: file_ext
    character(len=*), intent(in), optional :: file_name
    logical, intent(in), optional :: active
 ! ------------------------------------------------------------------------------
    !
    if (present(src_bbi))   this%src_bbi   = src_bbi
    if (present(src_bbx))   this%src_bbx   = src_bbx
    if (present(dst_bbi))   this%dst_bbi   = dst_bbi
    if (present(dst_bbx))   this%dst_bbx   = dst_bbx
    if (present(file_ext))  this%file_ext  = file_ext
    if (present(file_name)) this%file_name = file_name
    if (present(active))    this%active    = active
    !
    return
  end subroutine tVrtTile_set
  
  subroutine tVrtTile_get(this, src_bbx, src_bbi, dst_bbx, dst_bbi, &
    file_ext, file_name, active, hdrg, full_file_name)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrtTile) :: this
    type(tBbX), intent(out), optional :: src_bbx
    type(tBb),  intent(out), optional :: src_bbi
    type(tBbX), intent(out), optional :: dst_bbx
    type(tBb),  intent(out), optional :: dst_bbi
    character(len=*), intent(out), optional :: file_ext
    character(len=*), intent(out), optional :: file_name
    logical, intent(out), optional :: active
    type(tHdr), pointer, intent(out), optional :: hdrg
    character(len=*), intent(out), optional :: full_file_name
 ! ------------------------------------------------------------------------------
    !
    if (present(src_bbi))   src_bbi   = this%src_bbi
    if (present(src_bbx))   src_bbx   = this%src_bbx
    if (present(dst_bbi))   dst_bbi   = this%dst_bbi
    if (present(dst_bbx))   dst_bbx   = this%dst_bbx
    if (present(file_ext))  file_ext  = this%file_ext
    if (present(file_name)) file_name = this%file_name
    if (present(active))    active    = this%active
    if (present(hdrg))      hdrg     => this%hdrg
    if (present(full_file_name)) then
      full_file_name = trim(this%file_name)//trim(this%file_ext)
    end if
    !
    return
    end subroutine tVrtTile_get
    
! ==============================================================================
! ==============================================================================
! Type: tVrt
! ==============================================================================
! ==============================================================================
  
  subroutine tVrt_init_write(this, fp, data_type, f_tile, mv_tile, &
    dst_bbi, dst_bbx, tile_topol)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    character(len=*), intent(in) :: fp
    character(len=*), intent(in) :: data_type
    character(len=*), dimension(:), intent(in) :: f_tile
    character(len=*), dimension(:), intent(in) :: mv_tile
    type(tBb),  dimension(:), pointer, intent(in) :: dst_bbi
    type(tBbx), dimension(:), pointer, intent(in) :: dst_bbx
    integer(I4B), dimension(:,:), intent(in), optional :: tile_topol
    ! -- local
    type(tBB), pointer :: bbi_f => null(), bbi => null()
    type(tBBX), pointer :: bbx_f => null(), bbx => null()
    type(tVrtRawLine), pointer :: raw_line => null()
    type(tVrtRawLineSet), pointer :: raw_line_set => null()
    character(len=MXSLEN) :: mv, s
    integer(I4B), dimension(:), allocatable :: nr_tile, nc_tile
    integer(I4B), dimension(:), allocatable :: ir_off,  ic_off
    integer(I4B) :: ntiles, itile, nc, nr, bnc, bnr, ic, ir
    integer(I4B) :: xoff, yoff, xsize, ysize
    real(R8B), dimension(6) :: gt
    real(R8B), dimension(:), allocatable :: r8wk
    !
    logical :: vert_stacked
! ------------------------------------------------------------------------------
    !
    this%f = trim(fp)//'.vrt'
    ntiles = size(dst_bbi)
    allocate(this%raw)
    call this%raw%init(ntiles, source='ComplexSource')
    !
    vert_stacked = .true.
    allocate(r8wk(ntiles))
    do itile = 1, ntiles
      bbx => dst_bbx(ntiles)
      r8wk(itile) = bbx%cs
    end do
    if (minval(r8wk) == maxval(r8wk)) then
      vert_stacked = .false.
      if (.not.present(tile_topol)) then
        call errmsg('tVrt_init_write: tile_topol not found.')
      end if
    else
      vert_stacked = .true.
    end if
    deallocate(r8wk)
    !
    if (vert_stacked) then
      bbi_f => dst_bbi(ntiles)
      bbx_f => dst_bbx(ntiles)
    else
      allocate(bbi_f, bbx_f)
      call bbi_f%init(); call bbx_f%init()
      nr = maxval(tile_topol(1,:)); nc = maxval(tile_topol(2,:))
      allocate(nr_tile(nr), nc_tile(nc))
      do itile = 1, ntiles
        bbi => dst_bbi(itile); bbx => dst_bbx(itile)
        ir = tile_topol(1,itile); ic = tile_topol(2,itile)
        nr_tile(ir) = bbi%nrow; nc_tile(ic) = bbi%ncol;
        !
        bbi_f%ic0 = min(bbi_f%ic0, bbi%ic0); bbi_f%ic1 = max(bbi_f%ic1, bbi%ic1)
        bbi_f%ir0 = min(bbi_f%ir0, bbi%ir0); bbi_f%ir1 = max(bbi_f%ir1, bbi%ir1)
        bbi_f%ncol = bbi_f%ic1 - bbi_f%ic0 + 1
        bbi_f%nrow = bbi_f%ir1 - bbi_f%ir0 + 1
        !
        bbx_f%xll = min(bbx_f%xll, bbx%xll); bbx_f%xur = max(bbx_f%xur, bbx%xur)
        bbx_f%yll = min(bbx_f%yll, bbx%yll); bbx_f%yur = max(bbx_f%yur, bbx%yur)
        bbx_f%cs = bbx%cs
      end do
      !
      ! set the tile offsets
      allocate(ir_off(nr), ic_off(nc))
      ir_off = 0; ic_off = 0
      do ir = 2, nr
        ir_off(ir) = ir_off(ir-1) + nr_tile(ir-1)
      end do
      do ic = 2, nc
        ic_off(ic) = ic_off(ic-1) + nc_tile(ic-1)
      end do
      deallocate(nr_tile, nc_tile)
    end if
    !
    raw_line => this%raw%raw_beg(1)
    call raw_line%set(1, ta([bbi_f%ncol]))
    call raw_line%set(2, ta([bbi_f%nrow]))
    !
    raw_line => this%raw%raw_beg(2)
    !GT(1) x-coordinate of the upper-left corner of the upper-left pixel.
    !GT(2) w-e pixel resolution / pixel width.
    !GT(3) row rotation (typically zero).
    !GT(4) y-coordinate of the upper-left corner of the upper-left pixel.
    !GT(5) column rotation (typically zero).
    !GT(6) n-s pixel resolution / pixel height (negative value for a north-up image).
    gt(1) = bbx_f%xll
    gt(2) = bbx_f%cs
    gt(3) = R8ZERO
    gt(4) = bbx_f%yur
    gt(5) = R8ZERO
    gt(6) = -bbx_f%cs
    call raw_line%set(1, ta(arr=gt,sep_in=','))
    !
    raw_line => this%raw%raw_beg(3)
    call raw_line%set(1, trim(data_type))
    !
    raw_line => this%raw%raw_beg(4)
    call raw_line%set(1, trim(mv_tile(1))) ! using the missing value of the first tile !
    !
    do itile = 1, ntiles
      bbi => dst_bbi(itile)
      !
      raw_line_set => this%raw%raw_mid(itile)
      raw_line => raw_line_set%raw(2)
      call raw_line%set(1, ta([0]))
      call raw_line%set(2, trim(f_tile(itile)))
      !
      raw_line => raw_line_set%raw(4)
      call raw_line%set(1, ta([bbi%ncol]))
      call raw_line%set(2, ta([bbi%nrow]))
      call raw_line%set(3, trim(data_type))
      call raw_line%set(4, ta([bbi%ncol]))
      call raw_line%set(5, ta([1]))
      !
      xoff = 0; yoff = 0
      xsize = bbi%ncol; ysize = bbi%nrow
      raw_line => raw_line_set%raw(5)
      call raw_line%set(1, ta([xoff]))
      call raw_line%set(2, ta([yoff]))
      call raw_line%set(3, ta([xsize]))
      call raw_line%set(4, ta([ysize]))
      !
      if (vert_stacked) then
        xoff = 0; yoff = 0
        xsize = bbi_f%ncol; ysize = bbi_f%nrow
      else
        ir = tile_topol(1,itile); ic = tile_topol(2,itile)
        xoff = ic_off(ic); yoff = ir_off(ir)
        xsize = bbi%ncol; ysize = bbi%nrow
      end if
      raw_line => raw_line_set%raw(6)
      call raw_line%set(1, ta([xoff]))
      call raw_line%set(2, ta([yoff]))
      call raw_line%set(3, ta([xsize]))
      call raw_line%set(4, ta([ysize]))
      !
      raw_line => raw_line_set%raw(7)
      call raw_line%set(1, trim(mv_tile(itile)))
    end do
    !
    if (allocated(ir_off)) deallocate(ir_off)
    if (allocated(ic_off)) deallocate(ic_off)
    !
    return
  end subroutine tVrt_init_write
    
  subroutine tVrt_init(this, f)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    character(len=*), intent(in) :: f
    ! -- local
    type(tVrtRawLine), pointer :: vl => null()
    type(tVrtRaw), pointer :: raw => null()
    type(tVrtTile), pointer :: tile => null()
    type(tBb)  :: gbbi, src_lbbi, dst_lbbi
    type(tBbX) :: gbbx, src_lbbx, dst_lbbx
    !
    character(len=MXSLEN) :: s, tile_source, file_name, file_ext
    integer(I4B) :: itile, i, n, ir0, ir1
    real(R8B), dimension(6) :: gt
! ------------------------------------------------------------------------------
    !
    this%f = f
    call open_file(this%f, this%iu, 'r')
    call this%read_ntiles(tile_source)
    
    allocate(this%raw)
    raw => this%raw
    !
    call raw%init(this%ntiles, tile_source)
    call raw%read(this%iu)
    close(this%iu); this%iu = -1
    !
    vl => raw%raw_beg(i_VRTDataset_beg)
    call vl%get(1, i4v=gbbi%ncol); call vl%get(2, i4v=gbbi%nrow)
    gbbi%ic0 = 1; gbbi%ic1 = gbbi%ncol
    gbbi%ir0 = 1; gbbi%ir1 = gbbi%nrow
    this%dst_bbi = gbbi
    !
    vl => raw%raw_beg(i_GeoTransform); call vl%get(1, r8a=gt)
    !GT(1) x-coordinate of the upper-left corner of the upper-left pixel.
    !GT(2) w-e pixel resolution / pixel width.
    !GT(3) row rotation (typically zero).
    !GT(4) y-coordinate of the upper-left corner of the upper-left pixel.
    !GT(5) column rotation (typically zero).
    !GT(6) n-s pixel resolution / pixel height (negative value for a north-up image).
    gbbx%xll = gt(1); gbbx%yur = gt(4); gbbx%cs = gt(2)
    gbbx%xur = gbbx%xll + gbbi%ncol*gbbx%cs
    gbbx%yll = gbbx%yur - gbbi%nrow*gbbx%cs
    this%dst_bbx = gbbx
    
    allocate(this%tiles(this%ntiles))
    do itile = 1, this%ntiles
      !
      tile => this%tiles(itile)
      call tile%init()
      !
      vl => raw%raw_mid(itile)%raw(i_SourceFilename); call vl%get(2, sv=s)
      file_name = strip_ext(s); file_ext = get_ext(s)
      
      !i = index(s,'.',back=.true.); n = len_trim(s)
      !tile%file_name = s(1:i-1) !!!
      !tile%file_ext = s(i:n) !!!
      !
      ! destination (global) bounding box
      vl => raw%raw_mid(itile)%raw(i_DstRect)
      call vl%get(1, i4v=dst_lbbi%ic0)
      call vl%get(2, i4v=dst_lbbi%ir0)
      call vl%get(3, i4v=dst_lbbi%ncol)
      call vl%get(4, i4v=dst_lbbi%nrow)
      dst_lbbi%ic0 = dst_lbbi%ic0 + 1
      dst_lbbi%ir0 = dst_lbbi%ir0 + 1
      dst_lbbi%ic1 = dst_lbbi%ic0 + dst_lbbi%ncol - 1
      dst_lbbi%ir1 = dst_lbbi%ir0 + dst_lbbi%nrow - 1
      !tile%dst_bbi = dst_lbbi !!!
      !
      dst_lbbx%cs = gbbx%cs
      dst_lbbx%xll = gbbx%xll + (dst_lbbi%ic0-1)*dst_lbbx%cs
      dst_lbbx%yur = gbbx%yur - (dst_lbbi%ir0-1)*dst_lbbx%cs
      dst_lbbx%xur = dst_lbbx%xll + dst_lbbi%ncol*dst_lbbx%cs
      dst_lbbx%yll = dst_lbbx%yur - dst_lbbi%nrow*dst_lbbx%cs
      !tile%dst_bbx = dst_lbbx !!!
      !
      ! source (local) bounding box
      vl => raw%raw_mid(itile)%raw(i_SrcRect)
      call vl%get(3, i4v=src_lbbi%ncol)
      call vl%get(4, i4v=src_lbbi%nrow)
      src_lbbi%ic0 = 1; src_lbbi%ic1 = src_lbbi%ncol
      src_lbbi%ir0 = 1; src_lbbi%ir1 = src_lbbi%nrow
      !tile%src_bbi = src_lbbi !!!
      !
      src_lbbx%cs = gbbx%cs *(dst_lbbi%ncol/src_lbbi%ncol) ! calculate tile cell size
      src_lbbx%xll = dst_lbbx%xll
      src_lbbx%xur = dst_lbbx%xur
      src_lbbx%yll = dst_lbbx%yll
      src_lbbx%yur = dst_lbbx%yur
      !tile%src_bbx = src_lbbx !!!
      !
      call tile%set(file_name=file_name, file_ext=file_ext, &
        dst_bbi=dst_lbbi, dst_bbx=dst_lbbx, &
        src_bbi=src_lbbi, src_bbx=src_lbbx)
    end do
    !
    this%overlap = this%overlapping_tiles()
    this%full_data_read = .false.
    if (this%overlap) then
      call logmsg('*** The VRT contains are overlapping rasters ***')
    end if
    !
    return
  end subroutine tVrt_init
!
  subroutine tVrt_write(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    call open_file(this%f, this%iu, 'w')
    !
    call this%raw%write(this%iu)
    !
    close(this%iu)
    !
    return
  end subroutine tVrt_write
  !
  subroutine tVrt_read_full_tile(this, itile, hdrg, xi1, xi2, xi4, &
    mvi1, mvi2, mvi4, renum, nid, f_csv, clean_hdrg)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    integer(I4B), intent(in) :: itile
    !
    type(tHdr), intent(out), optional, pointer :: hdrg
    !
    integer(I1B), dimension(:,:), allocatable, intent(inout), optional :: xi1
    integer(I2B), dimension(:,:), allocatable, intent(inout), optional :: xi2
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: xi4
    !
    integer(I1B), intent(out), optional :: mvi1
    integer(I2B), intent(out), optional :: mvi2
    integer(I4B), intent(out), optional :: mvi4
    !
    logical, intent(in), optional :: renum
    integer(I4B), intent(inout), optional :: nid
    character(len=*), intent(in), optional :: f_csv
    !
    logical, intent(in), optional :: clean_hdrg
    
    ! -- local
    logical :: clean
    character(len=MXSLEN) :: f
    type(tVrtTile), pointer :: tile => null()
    logical :: renum_used
    integer(I4B) :: i
    !
! ------------------------------------------------------------------------------
    if (present(renum)) then
      renum_used = renum
    else
      renum_used = .false.
    end if
    !
    if (present(clean_hdrg)) then
      clean = clean_hdrg
    else
      clean = .true.
    end if
    !
    if (renum_used.and.(.not.present(nid))) then
      call errmsg('tVrt_read_full_tile: nid not specified.')
    end if
    !
    tile => this%tiles(itile)
    !
    select case(tile%file_ext)
      case('.flt','.FLT')
        !if (associated(tile%hdrg)) then
        !  call errmsg('tVrt_read_full_tile: program error')
        !end if
        !allocate(tile%hdrg)
        if (.not.associated(tile%hdrg)) then
          allocate(tile%hdrg)
        end if
        call tile%hdrg%read_full_grid(tile%file_name)
        !
        if (renum_used) then
          i = tile%hdrg%get_data_type()
          if (i /= i_i4) then
            call errmsg('tVrt_read_full_tile: integer*4 is only supported.')
          end if
          if (present(f_csv)) then
            call renumber(tile%hdrg%dat%xi4, tile%hdrg%hdr%mvi4, nid, f_csv)
          else
            call renumber(tile%hdrg%dat%xi4, tile%hdrg%hdr%mvi4, nid)
          end if
          if (.false.) then !debug
            f = trim(tile%file_name)//'.debug'
            call tile%hdrg%write(f)
          end if
        end if
        !
        if (present(xi4)) then
          if (allocated(xi4)) deallocate(xi4)
          if (.not.present(mvi4)) call errmsg('tVrt_read_full_tile: mvi4.')
          call tile%hdrg%get_grid(xi4=xi4, mvi4=mvi4)
        end if
        if (present(hdrg)) then
          hdrg => tile%hdrg
        end if
        if (clean) then
          call tile%hdrg%clean()
          deallocate(tile%hdrg)
        end if
        !
      case default
        call errmsg('tVrt_read_full_tile: not supported raster type.')
    end select
    !
    return
  end subroutine tVrt_read_full_tile
    
  subroutine tVrt_read_extent_native_cs(this, bbx)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    !
    type(tBbx), intent(in) :: bbx
    ! -- local
    type(tVrtTile), pointer :: tile => null()
    type(tHdr), pointer :: hdrg => null()
    type(tBbx) :: bbx_tile, bbx_read
    
    real(R8B), dimension(:), allocatable :: act_tile_cs
    integer(I4B), dimension(:), allocatable :: act_tile
    integer(I4B) :: i_act, n_act, itile
! ------------------------------------------------------------------------------
    
    ! get the active tiles
    call this%get_act_tile(bbx, act_tile, act_tile_cs)
    n_act = size(act_tile)
    
    ! read the tiles
    do i_act = 1, n_act
      itile = act_tile(i_act) 
      tile => this%tiles(itile); hdrg => tile%hdrg;
      !
      ! set the extent and native cell size
      call tile%get(src_bbx=bbx_tile)
      bbx_read = bbx; bbx_read%cs = bbx_tile%cs
      
      call hdrg%read_extent(tile%file_name, bbx_read, &
        i_uscl_nodata, i_dscl_nodata)
      !
      tile%dst_bbx = bbx_read
    end do
    !
    ! clean up
    deallocate(act_tile, act_tile_cs)
    !
    return
  end subroutine tVrt_read_extent_native_cs
  
  subroutine tVrt_read_extent(this, bbx, xi1, xi2, xi4, xi8, xr4, xr8, &
    mvi1, mvi2, mvi4, mvi8, mvr4, mvr8, i_uscl, i_dscl, xtile, &
    f_bin, clean_tile)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    !
    type(tBbx), intent(in) :: bbx
    !
    integer(I1B), dimension(:,:), allocatable, intent(inout), optional :: xi1
    integer(I2B), dimension(:,:), allocatable, intent(inout), optional :: xi2
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: xi4
    integer(I8B), dimension(:,:), allocatable, intent(inout), optional :: xi8
    real(R4B),    dimension(:,:), allocatable, intent(inout), optional :: xr4
    real(R8B),    dimension(:,:), allocatable, intent(inout), optional :: xr8
    !
    integer(I1B), intent(out), optional :: mvi1
    integer(I2B), intent(out), optional :: mvi2
    integer(I4B), intent(out), optional :: mvi4
    integer(I8B), intent(out), optional :: mvi8
    real(R4B),    intent(out), optional :: mvr4
    real(R8B),    intent(out), optional :: mvr8
    !
    integer(I4B), intent(in), optional :: i_uscl
    integer(I4B), intent(in), optional :: i_dscl
    !
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: xtile
    character(len=MXSLEN), dimension(:), allocatable, intent(inout), optional :: f_bin
    !
    logical, intent(in), optional :: clean_tile
    !
    ! -- local
    type(tVrtTile), pointer :: tile => null()
    type(tHdr), pointer :: hdrg => null()
    type(tHdrHdr), pointer :: hdr => null()
    type(tBbx) :: tbbx
    character(len=MXSLEN) :: f, fp
    logical :: loverlap, clean_tile_loc
    integer(I4B), dimension(:), allocatable :: act_tile
    integer(I4B) :: itile, n_act, i_act, nc, nr, ir, ic, jr, jc, n
    integer(I4B), dimension(:,:), allocatable :: xtile_loc
    real(R8B), dimension(:), allocatable :: act_tile_cs
    real(R8B) :: x, y
    !
    integer(I4B) :: i4v
    real(R4B)    :: r4v
    real(R8B)    :: r8v
 ! ------------------------------------------------------------------------------
    !
    if (present(clean_tile)) then
      clean_tile_loc = clean_tile
    else
      clean_tile_loc = .true.
    end if
    !
    ! get the active tiles
    call this%get_act_tile(bbx, act_tile, act_tile_cs)
    n_act = size(act_tile)
    !
    ! read the tiles
    do i_act = 1, n_act
      itile = act_tile(i_act) 
      tile => this%tiles(itile); hdrg => tile%hdrg;
      call hdrg%read_extent(tile%file_name, bbx, i_uscl, i_dscl)
      hdr => hdrg%hdr; call hdr%get_bbx(tbbx)
      tile%dst_bbx = tbbx
      if (.false.) then ! debug
        fp = trim(hdrg%fp)//'_scaled'
        call hdrg%write(fp)
      end if
    end do
    !
    ! append the grids tiles
    nc = (bbx%xur - bbx%xll)/bbx%cs
    nr = (bbx%yur - bbx%yll)/bbx%cs
    !
    allocate(xtile_loc(nc,nr)); xtile_loc = 0
    !
    if (present(xi1)) then
      call errmsg('tVrt_read_extent: xi1 not yet supported')
    end if
    if (present(xi2)) then
      call errmsg('tVrt_read_extent: xi2 not yet supported')
    end if
    if (present(xi4)) then
      if (allocated(xi4)) deallocate(xi4)
      allocate(xi4(nc,nr))
      itile = act_tile(1); hdrg => tile%hdrg; hdr => hdrg%hdr
      mvi4 = hdr%mvi4; xi4 = mvi4 
      do i_act = 1, n_act
        itile = act_tile(i_act); tile => this%tiles(itile)
        hdrg => tile%hdrg; hdr => hdrg%hdr
        call tile%get(dst_bbx=tbbx)
        do ir = 1, hdr%nrow; do ic = 1, hdr%ncol
          i4v = hdrg%dat%xi4(ic,ir)
          if (i4v /= hdr%mvi4) then
            call get_xy(x, y, ic, ir, tbbx%xll, tbbx%yur, tbbx%cs)
            call get_icr(jc, jr, x, y, bbx%xll, bbx%yur, bbx%cs)
            if (xi4(jc,jr) == mvi4) then
              xi4(jc,jr) = i4v
              xtile_loc(jc,jr) = itile
            end if
          end if
        end do; end do
      end do
    end if
    if (present(xi8)) then
      call errmsg('tVrt_read_extent: xi8 not yet supported')
    end if
    if (present(xr4)) then
      if (allocated(xr4)) deallocate(xr4)
      allocate(xr4(nc,nr))
      itile = act_tile(1); hdrg => tile%hdrg; hdr => hdrg%hdr
      mvr4 = hdr%mvr4; xr4 = mvr4 
      do i_act = 1, n_act
        itile = act_tile(i_act); tile => this%tiles(itile)
        hdrg => tile%hdrg; hdr => hdrg%hdr
        call tile%get(dst_bbx=tbbx)
        do ir = 1, hdr%nrow; do ic = 1, hdr%ncol
          r4v = hdrg%dat%xr4(ic,ir)
          if (r4v /= hdr%mvr4) then
            call get_xy(x, y, ic, ir, tbbx%xll, tbbx%yur, tbbx%cs)
            call get_icr(jc, jr, x, y, bbx%xll, bbx%yur, bbx%cs)
            if (xr4(jc,jr) == mvr4) then
              xr4(jc,jr) = r4v
              xtile_loc(jc,jr) = itile
            end if
          end if
        end do; end do
      end do
    end if
    if (present(xr8)) then
      if (allocated(xr8)) deallocate(xr8)
      allocate(xr8(nc,nr))
      itile = act_tile(1); hdrg => tile%hdrg; hdr => hdrg%hdr
      mvr8 = hdr%mvr8; xr8 = mvr8 
      do i_act = 1, n_act
        itile = act_tile(i_act); tile => this%tiles(itile)
        hdrg => tile%hdrg; hdr => hdrg%hdr
        call tile%get(dst_bbx=tbbx)
        do ir = 1, hdr%nrow; do ic = 1, hdr%ncol
          r8v = hdrg%dat%xr8(ic,ir)
          if (r8v /= hdr%mvr8) then
            call get_xy(x, y, ic, ir, tbbx%xll, tbbx%yur, tbbx%cs)
            call get_icr(jc, jr, x, y, bbx%xll, bbx%yur, bbx%cs)
            if (xr8(jc,jr) == mvr8) then
              xr8(jc,jr) = r8v
              xtile_loc(jc,jr) = itile
            end if
          end if
        end do; end do
      end do
    end if
    !
    if (present(xtile)) then
      if (allocated(xtile)) deallocate(xtile)
      allocate(xtile, source=xtile_loc)
    end if
    if (present(f_bin)) then
      n = maxval(xtile_loc)
      if (allocated(f_bin)) deallocate(f_bin)
      allocate(f_bin(n)); f_bin = ''
      do ir = 1, nr; do ic = 1, nc
        itile = xtile_loc(ic,ir)
        if (itile > 0) then
          if (len_trim(f_bin(itile)) == 0) then
            tile => this%tiles(itile)
            hdrg => tile%hdrg; hdr => hdrg%hdr
            f_bin(itile) = hdr%f_bin
            if (len_trim(f_bin(itile)) == 0) then
              call errmsg('Could not read BIN_DATA in envi file.')
            end if
          end if
        end if
      end do; end do
    end if
    !
    ! clean up
    deallocate(xtile_loc)
    if (clean_tile_loc) then
      do i_act = 1, n_act
        itile = act_tile(i_act) 
        tile => this%tiles(itile); hdrg => tile%hdrg
        call hdrg%clean()
        tile%dst_bbx = tile%src_bbx
      end do
    end if
    deallocate(act_tile, act_tile_cs)
    !
    return
  end subroutine tVrt_read_extent
  
  subroutine tVrt_get_act_tile(this, bbx, act_tile, act_tile_cs)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    type(tBbx), intent(in) :: bbx
    integer(I4B), intent(inout), dimension(:), allocatable :: act_tile
    real(R8B), intent(inout), dimension(:), allocatable :: act_tile_cs
    ! -- local
    type(tVrtTile), pointer :: tile => null()
    type(tBbx) :: tbbx
    character(len=MXSLEN) :: f
    integer(I4B) :: itile, iact, n_act
! ------------------------------------------------------------------------------
    
     ! first, check the extents
    do itile = 1, this%ntiles
      tile => this%tiles(itile)
      call tile%get(src_bbx=tbbx)
      if(.not.tbbx%match(bbx)) then
        call errmsg('tVrt_get_act_tile: non-matching grids.')
      end if
    end do
    !
    ! second, intersect the tiles
    do iact = 1, 2
      n_act = 0
      do itile = this%ntiles, 1, -1
        tile => this%tiles(itile)
        call tile%get(full_file_name=f)
        !
        call tile%get(src_bbx=tbbx)
        if (.not.(bbx_intersect(tbbx, bbx))) then
          call tile%set(active=.false.)
        else
          n_act = n_act + 1
          if (iact == 2) then
            act_tile(n_act) = itile
            act_tile_cs(n_act) = tbbx%cs
            call tile%get(full_file_name=f) 
            call logmsg('tVrt_get_act_tile: intersection found for '//trim(base_name(f)))
          end if
        end if
      end do
      if (iact == 1) then
        if (n_act > 0) then
          if (allocated(act_tile)) deallocate(act_tile)
          if (allocated(act_tile_cs)) deallocate(act_tile_cs)
          allocate(act_tile(n_act), act_tile_cs(n_act))
        else
          call errmsg('tVrt_get_act_tile: no intersected tile(s) found.')
        end if
      end if
    end do
    !
    return
  end subroutine tVrt_get_act_tile
    
  subroutine tVrt_buffer_data(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    !
    ! -- local
    type(tVrtTile), pointer :: tile => null()
    type(tHdr), pointer :: hdrg => null()
    logical :: active
    integer(I4B) :: itile
 ! ------------------------------------------------------------------------------
    !
    do itile = 1, this%ntiles
      tile => this%tiles(itile)
      call tile%get(active=active)
      if (active) then
        hdrg => tile%hdrg
        if (associated(hdrg%dat_buf)) then
          call hdrg%dat_buf%clean()
          deallocate(hdrg%dat_buf); hdrg%dat_buf => null()
        end if
        allocate(hdrg%dat_buf)
        call hdrg%dat%copy(hdrg%dat_buf)
      end if
    end do
    !
    return
  end subroutine tVrt_buffer_data
    
  subroutine tVrt_read_xy(this, x, y, cs, i1a, i2a, i4a, i8a, r4a, r8a, &
    mvi1, mvi2, mvi4, mvi8, mvr4, mvr8, i_uscl, i_dscl, xtile, f_bin)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    !
    real(R8B), dimension(:), intent(in) :: x
    real(R8B), dimension(:), intent(in) :: y
    real(R8B), dimension(:), intent(in) :: cs
    !
    integer(I1B), dimension(:), allocatable, intent(inout), optional :: i1a
    integer(I2B), dimension(:), allocatable, intent(inout), optional :: i2a
    integer(I4B), dimension(:), allocatable, intent(inout), optional :: i4a
    integer(I8B), dimension(:), allocatable, intent(inout), optional :: i8a
    real(R4B),    dimension(:), allocatable, intent(inout), optional :: r4a
    real(R8B),    dimension(:), allocatable, intent(inout), optional :: r8a
    !
    integer(I1B), intent(out), optional :: mvi1
    integer(I2B), intent(out), optional :: mvi2
    integer(I4B), intent(out), optional :: mvi4
    integer(I8B), intent(out), optional :: mvi8
    real(R4B),    intent(out), optional :: mvr4
    real(R8B),    intent(out), optional :: mvr8
    !
    integer(I4B), intent(in), optional :: i_uscl
    integer(I4B), intent(in), optional :: i_dscl
    !
    integer(I4B), dimension(:), allocatable, intent(inout), optional :: xtile
    character(len=MXSLEN), dimension(:), allocatable, intent(inout), optional :: f_bin
    !
    ! -- local
    type(tVrtTile), pointer :: tile => null()
    type(tHdr), pointer :: hdrg => null()
    type(tHdrHdr), pointer :: hdr => null()
    type(tBbx) :: tbbx
    !
    logical :: compressed_loc
    real(R4B) :: r4v, mvr4_tile
    real(R8B) :: r8v, mvr8_tile
    integer(I4B), dimension(:), allocatable :: xtile_loc
    integer(I4B) :: np, itile, i, n, ntile_act
    integer(I1B), dimension(:), allocatable :: flg1, flg3
    integer(I1B), dimension(:,:), allocatable :: flg2
 ! ------------------------------------------------------------------------------
    !
    if (present(f_bin)) then
      compressed_loc = .true.
    else
      compressed_loc = .false.
    end if
    !
    np = size(x)
    if ((np /= size(y)).or.(np /= size(cs))) then
      call errmsg('tLayerModel_read_xy: invalid dimensions.')
    end if
    !
    allocate(xtile_loc(np)); xtile_loc = 0
    !
    if (present(r4a)) then
      if (allocated(r4a)) deallocate(r4a)
      allocate(r4a(np))
    elseif (present(r8a)) then
      if (allocated(r8a)) deallocate(r8a)
      allocate(r8a(np))
    else
      call errmsg('tVrt_read_xy: ony r4 is yet supported.')
    end if
    !
    ! first, check the extents
    allocate(flg1(this%ntiles), flg2(np,this%ntiles), flg3(np))
    flg1 = 0; flg3 = 0
    !
    ntile_act = 0
    do itile = 1, this%ntiles
      tile => this%tiles(itile)
      call tile%get(src_bbx=tbbx)
      !
      n = 0
      do i = 1, np
        if (point_in_bb(reshape([x(i),y(i)],[2,1]), tbbx)) then
          flg2(i,itile) = 1; n = n + 1
          xtile_loc(i) = itile
        else
          flg2(i,itile) = 0
        end if
      end do
      if (n > 0) then
        ntile_act = ntile_act + 1
        flg1(itile) = 1
      end if
    end do
    !
    if (.not.this%full_data_read) then
      do itile = 1, this%ntiles
        tile => this%tiles(itile)
        if (flg1(itile) == 1) then
          hdrg => tile%hdrg
          if (compressed_loc) then
            call hdrg%init_read_xy(tile%file_name, mvr8=R8ZERO)
          else
            call hdrg%init_read_xy(tile%file_name)
          end if
        end if
      end do
    end if
    !
    do itile = this%ntiles, 1, -1
      if (flg1(itile) == 0) cycle
      hdrg => tile%hdrg; hdr => hdrg%hdr
      if (present(r4a)) then
        do i = 1, np
          if (flg2(i,itile) == 1) then
            mvr4_tile = hdr%mvr4
            call hdrg%read_val_xy(x(i), y(i), cs(i), r4v=r4v, &
              i_uscl=i_uscl, i_dscl=i_dscl)
            if (flg3(i) == 0) then
              mvr4 = mvr4_tile
              r4a(i) = mvr4
            end if
            if (r4v == mvr4_tile) r4v = mvr4
            if (r4v /= mvr4) then
              r4a(i) = r4v
            end if
            flg3(i) = 1
          end if
        end do
      end if
      if (present(r8a)) then
        do i = 1, np
          if (flg2(i,itile) == 1) then
            mvr8_tile = hdr%mvr8
            call hdrg%read_val_xy(x(i), y(i), cs(i), r8v=r8v, &
              i_uscl=i_uscl, i_dscl=i_dscl)
            if (flg3(i) == 0) then
              mvr8 = mvr8_tile
              r8a(i) = mvr8
            end if
            if (r8v == mvr8_tile) r8v = mvr8
            if (r8v /= mvr8) then
              r8a(i) = r8v
            end if
            flg3(i) = 1
          end if
        end do
      end if
    end do
    !
    if (present(xtile)) then
      if (allocated(xtile)) deallocate(xtile)
      allocate(xtile, source=xtile_loc)
    end if
    if (present(f_bin)) then
      if (allocated(f_bin)) deallocate(f_bin)
      allocate(f_bin(ntile_act)); f_bin = ''
      do itile = 1, this%ntiles
        if (flg1(itile) == 1) then
          if (len_trim(f_bin(itile)) == 0) then
            tile => this%tiles(itile)
            hdrg => tile%hdrg; hdr => hdrg%hdr
            f_bin(itile) = hdr%f_bin
            if (len_trim(f_bin(itile)) == 0) then
              call errmsg('Could not read BIN_DATA in envi file.')
            end if
          end if
        end if
      end do
    end if
    
    !do itile = 1, this%ntiles
    !  tile => this%tiles(itile)
    !  if (flg1(itile) == 1) then
    !    hdrg => tile%hdrg
    !    call hdrg%clean_x()
    !  end if
    !end do
    !
    if (allocated(flg1)) deallocate(flg1)
    if (allocated(flg2)) deallocate(flg2)
    if (allocated(flg3)) deallocate(flg3)
    if (allocated(xtile_loc)) deallocate(xtile_loc)
    !
    return
  end subroutine tVrt_read_xy
    
  subroutine tVrt_set_nodata(this, mvr8)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    !
    real(R8B), intent(in), optional :: mvr8
    ! -- local
    type(tHdr), pointer :: hdrg => null()
    type(tHdrHdr), pointer :: hdr => null()
    type(tVrtTile), pointer :: tile => null()
    integer(I4B) :: itile
! ------------------------------------------------------------------------------
    !
    do itile = 1, this%ntiles
      tile => this%tiles(itile); hdrg => tile%hdrg; hdr => hdrg%hdr
      
    end do
    !
    return
  end subroutine tVrt_set_nodata
  
  subroutine tVrt_init_tile_raster(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    ! -- local
    type(tHdr), pointer :: hdrg => null()
    type(tVrtTile), pointer :: tile => null()
    !
    integer(I4B) :: itile
! ------------------------------------------------------------------------------
    !
    do itile = 1, this%ntiles
      tile => this%tiles(itile)
      select case(tile%file_ext)
      case('.flt','.FLT')
        allocate(tile%hdrg)
        hdrg => tile%hdrg
        call hdrg%hdr_get_val_init(tile%file_name)
      case default
        call errmsg('tVrt_init_tile_raster: not supported raster type.')
      end select
    end do
    !
    return
  end subroutine tVrt_init_tile_raster
  
  subroutine tvrt_read_ntiles(this, tile_source)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    character(len=MXSLEN), intent(out) :: tile_source
    ! -- local    
    integer(I4B) :: ios
    character(len=MXSLEN) :: s
    integer(I4B) :: n
! ------------------------------------------------------------------------------
    this%ntiles = 0
    do while(.true.)
      read(unit=this%iu,iostat=ios,fmt='(a)') s
      if (ios /= 0) exit
      if (len_trim(s) == 0) cycle
      if (index(s,'Source>') > 0) then
        s = adjustl(s)
        if ((s(1:1) == '<').and.(s(2:2) /= '/')) then
          this%ntiles = this%ntiles + 1
          n = len_trim(s)
          tile_source = s(2:n-1)
        end if
      end if
    end do
    rewind(this%iu)
    !
    return
  end subroutine tvrt_read_ntiles 
!
  function tVrt_get_ntiles(this) result(ntiles)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    integer(I4B) :: ntiles
    ! -- local    
! ------------------------------------------------------------------------------
    ntiles = this%ntiles
    !
    return
  end function tVrt_get_ntiles 

  subroutine tVrt_get_bb(this, itile, src_bbi, src_bbx, dst_bbi, dst_bbx)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    integer(I4B), intent(in), optional :: itile
    type(tBb),  intent(out), optional :: src_bbi
    type(tBbX), intent(out), optional :: src_bbx
    type(tBb),  intent(out), optional :: dst_bbi
    type(tBbX), intent(out), optional :: dst_bbx
    ! -- local
    type(tVrtTile), pointer :: tile => null()
! ------------------------------------------------------------------------------
    if (present(itile)) then
      tile => this%tiles(itile)
      if (present(src_bbi)) then
        src_bbi = tile%src_bbi
      end if
      if (present(src_bbx)) then
        src_bbx = tile%src_bbx
      end if
      if (present(dst_bbi)) then
        dst_bbi = tile%dst_bbi
      end if
      if (present(dst_bbx)) then
        dst_bbx = tile%dst_bbx
      end if
    else
      if (present(dst_bbi)) then
        dst_bbi = this%dst_bbi
      end if
      if (present(dst_bbx)) then
        dst_bbx = this%dst_bbx
      end if
    end if
    !
    return
  end subroutine tVrt_get_bb
  
  function tVrt_get_tile_file(this, itile) result(f)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    integer(I4B), intent(in) :: itile
    character(len=MXSLEN) :: f
    ! -- local
    type(tVrtTile), pointer :: tile => null()
! ------------------------------------------------------------------------------
    tile => this%tiles(itile)
    f = tile%file_name
    !
    return
  end function tVrt_get_tile_file
  
  function tVrt_overlapping_tiles(this) result(overlap)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    logical :: overlap
    ! -- local
    type(tVrtTile), pointer :: t1 => null(), t2 => null()
    integer(I4B) :: i, j
! ------------------------------------------------------------------------------
    overlap = .false.
    do i = 1, this%ntiles
      t1 => this%tiles(i)
      do j = i + 1, this%ntiles
        t2 => this%tiles(j)
        overlap = bbi_intersect(t1%dst_bbi, t2%dst_bbi)
        if (overlap) exit
      end do
      if (overlap) exit
    end do
    !
    return
  end function tVrt_overlapping_tiles
  
  subroutine tVrt_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    ! -- local
    logical :: lop
    type(tVrtTile), pointer :: tile => null()
    !
    integer(I4B) :: itile
! ------------------------------------------------------------------------------
    !
    if (this%iu == 0) then
      lop = .false.
    else
      inquire(this%iu, opened=lop)
    end if
    if (lop) then
      close(this%iu)
    end if
    !
    if (associated(this%raw)) then
      call this%raw%clean()
      deallocate(this%raw)
    end if
    !
    if (associated(this%tiles)) then
      do itile = 1, this%ntiles
        tile => this%tiles(itile)
        if (associated(tile%hdrg)) then
          call tile%hdrg%clean()
          deallocate(tile%hdrg)
        end if
      end do
      deallocate(this%tiles)
    end if
    !
    this%iu       = 0
    this%raw      => null()
    this%f        = ''
    !this%dst_bbi = 
    !this%dst_bbx =
    this%ntiles   = 0
    !
    return
  end subroutine tVrt_clean

  subroutine tVrt_clean_x(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVrt) :: this
    ! -- local
    type(tVrtTile), pointer :: tile => null()
    !
    integer(I4B) :: itile
! ------------------------------------------------------------------------------
    if (associated(this%tiles)) then
      do itile = 1, this%ntiles
        tile => this%tiles(itile)
        if (associated(tile%hdrg)) then
          call tile%hdrg%clean_dat()
          call tile%hdrg%clean_dat_buf()
        end if
      end do
    end if
    !
    return
  end subroutine tVrt_clean_x
  
end module vrt_module

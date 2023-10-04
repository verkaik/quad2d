module hdrModule
  ! modules
  use utilsmod, only: I1B, I2B, I4B, I8B, R4B, R8B, i_i1, i_i2, i_i4, i_i8, i_r4, i_r8, &
    MXSLEN, change_case, logmsg, errmsg, tNum, open_file, parse_line, ta, &
    get_bb_extent, R8HALF, R8ONE, tBb, tBbX, get_slash, get_xy, cast_2darray, strip_ext, &
    i_p, i_n, i_s, i_w, i_e, i_nw, i_ne, i_sw, i_se, bilinear_interpolation, get_stencil, &
    check_nan
  !
  implicit none
  !
  private
  !
  integer(I1B), parameter :: I1ZERO = 0
  integer(I2B), parameter :: I2ZERO = 0
  integer(I4B), parameter :: I4ZERO = 0
  integer(I8B), parameter :: I8ZERO = 0
  real(R4B),    parameter :: R4ZERO = 0.0
  real(R8B),    parameter :: R8ZERO = 0.d0
  !
  ! upscaling
  integer(I4B), parameter, public :: i_uscl_nodata    = 0 ! no data
  integer(I4B), parameter, public :: i_uscl_spec      = 1 ! special/ibound
  integer(I4B), parameter, public :: i_uscl_arith     = 2 ! arithmetic <====
  integer(I4B), parameter, public :: i_uscl_geom      = 3 ! geometric <====
  integer(I4B), parameter, public :: i_uscl_sumq      = 4 ! sum(Q)
  integer(I4B), parameter, public :: i_uscl_sumcdr    = 5 ! sum(cond)*ratio <====
  integer(I4B), parameter, public :: i_uscl_invc      = 6 ! inverse (c)
  integer(I4B), parameter, public :: i_uscl_mostfr    = 7 ! most freq. occ
  integer(I4B), parameter, public :: i_uscl_suminvcvr = 8 ! sum(1/c)*ratio
  integer(I4B), parameter, public :: i_uscl_perc      = 9 ! percentile
  integer(I4B), parameter, public :: i_uscl_arithnd   = 10 ! average including locations with nodata (rch/evt)  <====
  integer(I4B), parameter, public :: n_uscl = i_uscl_arithnd
  
  character(len=MXSLEN), dimension(n_uscl), public, parameter :: &
    uscl_names = ['spec', 'arith', 'geom', 'sumq', 'sumcdr', &
      'invc', 'mostfr', 'suminvcvr', 'perc', 'arithnd']
  
  ! downscaling
  integer(I4B), parameter, public :: i_dscl_nodata       = 0 ! no data
  integer(I4B), parameter, public :: i_dscl_nointp       = 1 ! no interpolation <====
  integer(I4B), parameter, public :: i_dscl_sample_scale = 2 ! sample with cell size scaling <====
  integer(I4B), parameter, public :: i_dscl_intp         = 3 ! interpolation <====
  integer(I4B), parameter, public :: n_dscl = i_dscl_intp
  !
  character(len=MXSLEN), dimension(n_dscl), public, parameter :: &
    dscl_names = ['nointp','sample_scale','intp']
  !
  integer(I4B), parameter, public :: i_flt = 1
  integer(I4B), parameter, public :: i_env = 2
  
  type, public :: tHdrHdr
    integer(I4B) :: i_file_type
    integer(I4B) :: i_data_type
    !
    integer(I4B) :: ncol
    integer(I4B) :: nrow
    integer(I4B) :: nbits
    character(len=MXSLEN) :: pixeltype
    integer(I4B) :: i_uscl_type
    integer(I4B) :: i_dscl_type
    !
    real(R4B) :: xllr4
    real(R4B) :: xurr4
    real(R4B) :: yllr4 
    real(R4B) :: yurr4 
    real(R4B) :: csr4
    !
    real(R8B) :: xllr8
    real(R8B) :: xurr8
    real(R8B) :: yllr8
    real(R8B) :: yurr8
    real(R8B) :: csr8
    !
    logical :: lmvi1
    logical :: lmvi2
    logical :: lmvi4
    logical :: lmvi8
    logical :: lmvr4
    logical :: lmvr8
    !
    integer(I1B) :: mvi1
    integer(I2B) :: mvi2
    integer(I4B) :: mvi4
    integer(I8B) :: mvi8
    real(R4B)    :: mvr4
    real(R8B)    :: mvr8
    !
    logical :: l_read_full_grid = .false.
    !
    integer(I4B) :: ic0
    integer(I4B) :: ic1
    integer(I4B) :: ir0
    integer(I4B) :: ir1
    !
    character(len=MXSLEN) :: f_bin
  contains
    procedure :: init  => hdrhdr_init
    procedure :: clean => hdrhdr_clean
    procedure :: read  => hdrhdr_read
    procedure :: hdrhdr_read_ehdr, hdrhdr_read_envi
    procedure :: clip    => hdrhdr_clip
    procedure :: get_bbx => hdrhdr_get_bbx
    !procedure :: write => hdrhdr_write
  end type thdrHdr
  
  type :: tHdrExt
    integer(I1B), dimension(:), allocatable :: xi1
    integer(I2B), dimension(:), allocatable :: xi2
    integer(I4B), dimension(:), allocatable :: xi4
    integer(I8B), dimension(:), allocatable :: xi8
    real(R4B),    dimension(:), allocatable :: xr4
    real(R8B),    dimension(:), allocatable :: xr8
  contains
    procedure :: clean => tHdrExt_clean
    procedure :: copy  => tHdrExt_copy
  end type tHdrExt
  !
  type tHdrData
    integer(I1B), dimension(:,:), allocatable :: xi1
    integer(I2B), dimension(:,:), allocatable :: xi2
    integer(I4B), dimension(:,:), allocatable :: xi4
    integer(I8B), dimension(:,:), allocatable :: xi8
    real(R4B),    dimension(:,:), allocatable :: xr4
    real(R8B),    dimension(:,:), allocatable :: xr8
    !
    logical :: ext_active = .false.
    type(thdrExt), pointer :: ext_n  => null()
    type(thdrExt), pointer :: ext_s  => null()
    type(thdrExt), pointer :: ext_w  => null()
    type(thdrExt), pointer :: ext_e  => null()
    type(thdrExt), pointer :: ext_nw => null()
    type(thdrExt), pointer :: ext_ne => null()
    type(thdrExt), pointer :: ext_sw => null()
    type(thdrExt), pointer :: ext_se => null()
  contains
    procedure :: init  => tHdrData_init
    procedure :: clean => tHdrData_clean
    procedure :: copy  => tHdrData_copy
  end type tHdrData
  
  type, public :: tHdr
    character(len=MXSLEN)      :: fp = ''
    type(thdrHdr), pointer     :: hdr        => null()
    type(thdrHdr), pointer     :: hdr_src => null()
    type(thdrHdr), pointer     :: hdr_no_ext => null()
    integer(I4B) :: iu_bin = 0
    integer(I4B) :: i_data_type = 0
    !
    type(tHdrData), pointer :: dat     => null()
    type(tHdrData), pointer :: dat_buf => null()
    
    !integer(I1B), dimension(:,:), allocatable :: xi1
    !integer(I2B), dimension(:,:), allocatable :: xi2
    !integer(I4B), dimension(:,:), allocatable :: xi4
    !integer(I8B), dimension(:,:), allocatable :: xi8
    !real(R4B),    dimension(:,:), allocatable :: xr4
    !real(R8B),    dimension(:,:), allocatable :: xr8
    !!
    !logical :: ext_active = .false.
    !type(thdrExt), pointer :: ext_n  => null()
    !type(thdrExt), pointer :: ext_s  => null()
    !type(thdrExt), pointer :: ext_w  => null()
    !type(thdrExt), pointer :: ext_e  => null()
    !type(thdrExt), pointer :: ext_nw => null()
    !type(thdrExt), pointer :: ext_ne => null()
    !type(thdrExt), pointer :: ext_sw => null()
    !type(thdrExt), pointer :: ext_se => null()
    !
    logical :: map_active = .false.
    logical :: mapped_all = .false.
  contains
 !   generic   :: read_hdr => hdr_read_grid_hdr
    procedure :: init            => hdr_init
    procedure :: init_map        => hdr_init_map
    procedure :: init_read_xy    => hdr_init_read_xy
    procedure :: clean           => hdr_clean
    procedure :: clean_dat       => hdr_clean_dat
    procedure :: clean_dat_buf   => hdr_clean_dat_buf
    procedure :: set_mv          => hdr_set_mv
    procedure :: get_file_type   => hdr_get_file_type
    procedure :: get_data_type   => hdr_get_data_type
    procedure :: set_i_data_type => hdr_set_i_data_type
    procedure :: check_empty     => hdr_check_empty
    procedure :: read_grid       => hdr_read_grid
    procedure :: map_grid        => hdr_map_grid
    !
    !generic :: read_extent => hdr_read_extent_i4_r8
    !procedure :: hdr_read_extent_i4_r8
    !
    !generic   :: read_block => hdr_read_block_i4_r8
    !procedure :: hdr_read_block_i4_r8
    !
    procedure :: read_extent       => hdr_read_extent
    procedure :: read_clip_grid    => hdr_read_clip_grid
    procedure :: read_full_grid    => hdr_read_full_grid
    procedure :: up_scale          => hdr_up_scale
    procedure :: down_scale_nointp => hdr_down_scale_nointp
    procedure :: down_scale_intp   => hdr_down_scale_intp
    !
    procedure :: get_global_bb  => hdr_get_global_bb
    procedure :: get_grid       => hdr_get_grid
    !
    procedure :: replace_grid => hdr_replace_grid
    !
    generic :: set_grid => hdr_set_grid_i4,  hdr_set_grid_r4, hdr_set_grid_i8
    procedure :: hdr_set_grid_i4, hdr_set_grid_r4, hdr_set_grid_i8
    !
    procedure :: hdr_uscl_read_check_r4, hdr_uscl_read_check_r8
    !
    procedure :: hdr_get_val_init!, hdr_get_val_clean
    !generic :: hdr_get_val => hdr_get_val_r8
    procedure :: hdr_get_val_r8
    procedure :: read_val    => hdr_read_val
    procedure :: read_val_xy => hdr_read_val_xy
    procedure :: read_arr    => hdr_read_arr

    procedure :: write => hdr_write
    !
    procedure :: include_ext => hdr_include_ext
    procedure :: exclude_ext => hdr_exclude_ext
    procedure :: get_ext_val => hdr_get_ext_val
    procedure :: get_nbr_val => hdr_get_nbr_val
    !
    procedure :: check_val => hdr_check_val
    !
    procedure :: get_base_name => hdr_get_base_name
  end type tHdr
  !
  interface writeflt
    module procedure :: writeflt_i1_r4
    module procedure :: writeflt_i2_r4
    module procedure :: writeflt_i4_r4
    module procedure :: writeflt_r4_r4
    module procedure :: writeflt_i1_r8
    module procedure :: writeflt_i2_r8
    module procedure :: writeflt_i4_r8
    module procedure :: writeflt_r4_r8
  end interface
  private :: writeflt_i1_r4, writeflt_i2_r4, writeflt_i4_r4, writeflt_r4_r4
  private :: writeflt_i1_r8, writeflt_i2_r8, writeflt_i4_r8, writeflt_r4_r8
  !  
  ! work arrays
  integer(I1B), dimension(:,:), pointer :: i1wk2d => null()
  integer(I2B), dimension(:,:), pointer :: i2wk2d => null()
  integer(I4B), dimension(:,:), pointer :: i4wk2d => null()
  integer(R4B), dimension(:,:), pointer :: r4wk2d => null()
  !
  ! set public
  public :: writeflt, write_env
  
contains
  
! ==============================================================================
! ==============================================================================
! tHdrExt
! ==============================================================================
  subroutine tHdrExt_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdrExt) :: this
    !
    ! -- local
! ------------------------------------------------------------------------------
    if (allocated(this%xi1)) deallocate(this%xi1)
    if (allocated(this%xi2)) deallocate(this%xi2)
    if (allocated(this%xi4)) deallocate(this%xi4)
    if (allocated(this%xi8)) deallocate(this%xi8)
    if (allocated(this%xr4)) deallocate(this%xr4)
    if (allocated(this%xr8)) deallocate(this%xr8)
    !
    return
  end subroutine tHdrExt_clean
  
  subroutine tHdrExt_copy(this, tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdrExt) :: this
    type(tHdrExt), intent(inout) :: tgt
    !
    ! -- local
! ------------------------------------------------------------------------------
    if (allocated(this%xi1)) allocate(tgt%xi1, source=this%xi1)
    if (allocated(this%xi2)) allocate(tgt%xi2, source=this%xi2)
    if (allocated(this%xi4)) allocate(tgt%xi4, source=this%xi4)
    if (allocated(this%xi8)) allocate(tgt%xi8, source=this%xi8)
    if (allocated(this%xr4)) allocate(tgt%xr4, source=this%xr4)
    if (allocated(this%xr8)) allocate(tgt%xr8, source=this%xr8)
    !
    return
  end subroutine tHdrExt_copy
  
! ==============================================================================
! ==============================================================================
! tHdrData
! ==============================================================================
    
  subroutine tHdrData_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdrData) :: this
    !
    ! -- local
! ------------------------------------------------------------------------------
    call this%clean()
    !
    this%ext_active = .false.
    this%ext_n  => null()
    this%ext_s  => null()
    this%ext_w  => null()
    this%ext_e  => null()
    this%ext_nw => null()
    this%ext_ne => null()
    this%ext_sw => null()
    this%ext_se => null()
    !
    return
  end subroutine tHdrData_init
  
  subroutine tHdrData_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdrData) :: this
    !
    ! -- local
! ------------------------------------------------------------------------------
    if (allocated(this%xi1)) deallocate(this%xi1)
    if (allocated(this%xi2)) deallocate(this%xi2)
    if (allocated(this%xi4)) deallocate(this%xi4)
    if (allocated(this%xi8)) deallocate(this%xi8)
    if (allocated(this%xr4)) deallocate(this%xr4)
    if (allocated(this%xr8)) deallocate(this%xr8)
    !
    if (associated(this%ext_n))  call this%ext_n%clean()
    if (associated(this%ext_s))  call this%ext_s%clean()
    if (associated(this%ext_w))  call this%ext_w%clean()
    if (associated(this%ext_e))  call this%ext_e%clean()
    if (associated(this%ext_nw)) call this%ext_nw%clean()
    if (associated(this%ext_ne)) call this%ext_ne%clean()
    if (associated(this%ext_sw)) call this%ext_sw%clean()
    if (associated(this%ext_se)) call this%ext_se%clean()
    !
    if (associated(this%ext_n))  deallocate(this%ext_n)
    if (associated(this%ext_s))  deallocate(this%ext_s)
    if (associated(this%ext_w))  deallocate(this%ext_w)
    if (associated(this%ext_e))  deallocate(this%ext_e)
    if (associated(this%ext_nw)) deallocate(this%ext_nw)
    if (associated(this%ext_ne)) deallocate(this%ext_ne)
    if (associated(this%ext_sw)) deallocate(this%ext_sw)
    if (associated(this%ext_se)) deallocate(this%ext_se)
    !
    this%ext_n  => null()
    this%ext_s  => null()
    this%ext_w  => null()
    this%ext_e  => null()
    this%ext_nw => null()
    this%ext_ne => null()
    this%ext_sw => null()
    this%ext_se => null()
    !
    return
  end subroutine tHdrData_clean
  
  subroutine tHdrData_copy(this, tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdrData) :: this
    type(tHdrData) :: tgt
! ------------------------------------------------------------------------------
    if (allocated(this%xi1)) allocate(tgt%xi1, source=this%xi1)
    if (allocated(this%xi2)) allocate(tgt%xi2, source=this%xi2)
    if (allocated(this%xi4)) allocate(tgt%xi4, source=this%xi4)
    if (allocated(this%xi8)) allocate(tgt%xi8, source=this%xi8)
    if (allocated(this%xr4)) allocate(tgt%xr4, source=this%xr4)
    if (allocated(this%xr8)) allocate(tgt%xr8, source=this%xr8)
    !
    tgt%ext_active = this%ext_active
    if (associated(this%ext_n)) then
      allocate(tgt%ext_n); call this%ext_n%copy(tgt%ext_n)
    end if
    if (associated(this%ext_s)) then
      allocate(tgt%ext_s); call this%ext_s%copy(tgt%ext_s)
    end if
    if (associated(this%ext_w)) then
      allocate(tgt%ext_w); call this%ext_w%copy(tgt%ext_w)
    end if
    if (associated(this%ext_e)) then
      allocate(tgt%ext_e); call this%ext_e%copy(tgt%ext_e)
    end if
    if (associated(this%ext_nw)) then
      allocate(tgt%ext_nw); call this%ext_nw%copy(tgt%ext_nw)
    end if
    if (associated(this%ext_ne)) then
      allocate(tgt%ext_ne); call this%ext_ne%copy(tgt%ext_ne)
    end if
    if (associated(this%ext_sw)) then
      allocate(tgt%ext_sw); call this%ext_sw%copy(tgt%ext_sw)
    end if
    if (associated(this%ext_se)) then
      allocate(tgt%ext_se); call this%ext_se%copy(tgt%ext_se)
    end if
    !
    return
  end subroutine tHdrData_copy
  
! ==============================================================================
! ==============================================================================
! tHdrHdr
! ==============================================================================
! ==============================================================================
  subroutine hdrhdr_init(this, ncol, nrow, nbits, pixeltype, i_uscl_type, &
    i_dscl_type, xllr4, xurr4, yllr4, yurr4, csr4, xllr8, xurr8, yllr8, yurr8, &
    csr8, mvi1, mvi2, mvi4, mvi8, mvr4, mvr8, l_read_full_grid, ic0, ic1, ir0, ir1)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdrHdr) :: this
    !
    integer(I4B), intent(in), optional :: ncol
    integer(I4B), intent(in), optional :: nrow
    integer(I4B), intent(in), optional :: nbits
    character(len=*), intent(in), optional :: pixeltype
    integer(I4B), intent(in), optional :: i_uscl_type
    integer(I4B), intent(in), optional :: i_dscl_type
    !
    real(R4B), intent(in), optional :: xllr4
    real(R4B), intent(in), optional :: xurr4
    real(R4B), intent(in), optional :: yllr4 
    real(R4B), intent(in), optional :: yurr4 
    real(R4B), intent(in), optional :: csr4
    !
    real(R8B), intent(in), optional :: xllr8
    real(R8B), intent(in), optional :: xurr8
    real(R8B), intent(in), optional :: yllr8
    real(R8B), intent(in), optional :: yurr8
    real(R8B), intent(in), optional :: csr8
    !
    integer(I1B), intent(in), optional :: mvi1
    integer(I2B), intent(in), optional :: mvi2
    integer(I4B), intent(in), optional :: mvi4
    integer(I8B), intent(in), optional :: mvi8
    real(R4B), intent(in), optional    :: mvr4
    real(R8B), intent(in), optional    :: mvr8
    !
    logical, intent(in), optional :: l_read_full_grid
    !
    integer(I4B), intent(in), optional :: ic0
    integer(I4B), intent(in), optional :: ic1
    integer(I4B), intent(in), optional :: ir0
    integer(I4B), intent(in), optional :: ir1
    !
    ! -- local
! ------------------------------------------------------------------------------
    !
    this%i_file_type = I4ZERO
    !
    this%ncol = I4ZERO
    this%nrow = I4ZERO
    this%nbits = I4ZERO
    this%pixeltype = ''
    this%i_uscl_type = i_uscl_nodata
    this%i_dscl_type = i_dscl_nodata
    !
    this%xllr4 = huge(R4ZERO)
    this%xurr4 = huge(R4ZERO)
    this%yllr4 = huge(R4ZERO)
    this%yurr4 = huge(R4ZERO)
    this%csr4  = R4ZERO
    !
    this%xllr8 = huge(R8ZERO)
    this%xurr8 = huge(R8ZERO)
    this%yllr8 = huge(R8ZERO)
    this%yurr8 = huge(R8ZERO)
    this%csr8  = R8ZERO
    !
    this%lmvi1 = .false.
    this%lmvi2 = .false.
    this%lmvi4 = .false.
    this%lmvi8 = .false.
    this%lmvr4 = .false.
    this%lmvr8 = .false.
    !
    this%mvi1 = huge(I1ZERO)
    this%mvi2 = huge(I2ZERO)
    this%mvi4 = huge(I4ZERO)
    this%mvi8 = huge(I8ZERO)
    this%mvr4 = huge(R4ZERO)
    this%mvr8 = huge(R8ZERO)
    !
    this%l_read_full_grid = .false.
    !
    this%ic0 = huge(I4ZERO)
    this%ic1 = I4ZERO
    this%ir0 = huge(I4ZERO)
    this%ir1 = I4ZERO
    !
    ! Overwrite default for optional arguments.
    if(present(ncol))        this%ncol        = ncol
    if(present(nrow))        this%nrow        = nrow
    if(present(nbits))       this%nbits       = nbits
    if(present(pixeltype))   this%pixeltype   = pixeltype
    if(present(i_uscl_type)) this%i_uscl_type = i_uscl_type
    if(present(i_dscl_type)) this%i_dscl_type = i_dscl_type
    !
    if(present(xllr4)) this%xllr4 = xllr4
    if(present(xurr4)) this%xurr4 = xurr4
    if(present(yllr4)) this%yllr4 = yllr4
    if(present(yurr4)) this%yurr4 = yurr4
    if(present(csr4))  this%csr4  = csr4
    !
    if(present(xllr8)) this%xllr8 = xllr8
    if(present(xurr8)) this%xurr8 = xurr8
    if(present(yllr8)) this%yllr8 = yllr8
    if(present(yurr8)) this%yurr8 = yurr8
    if(present(csr8))  this%csr8  = csr8
    !
    if(present(mvi1)) then
      this%mvi1 = mvi1; this%lmvi1 = .true.
    end if
    if(present(mvi2)) then
      this%mvi2 = mvi2; this%lmvi2 = .true.
    end if
    if(present(mvi4)) then
      this%mvi4 = mvi4; this%lmvi4 = .true.
    end if
    if(present(mvi8)) then
      this%mvi8 = mvi8; this%lmvi8 = .true.
    end if
    if(present(mvr4)) then
      this%mvr4 = mvr4; this%lmvr4 = .true.
    end if
    if(present(mvr8)) then
      this%mvr8 = mvr8; this%lmvr8 = .true.
    end if
    !
    if(present(l_read_full_grid)) this%l_read_full_grid = l_read_full_grid
    !
    if(present(ic0)) this%ic0 = ic0
    if(present(ic1)) this%ic1 = ic1
    if(present(ir0)) this%ir0 = ir0
    if(present(ir1)) this%ir1 = ir1
    !
    this%f_bin = ''
    !
    return
 end subroutine hdrhdr_init
  
subroutine hdrhdr_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdrHdr) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    call this%init()
    !
    return
  end subroutine hdrhdr_clean

  subroutine hdrhdr_read(this, fp)
! ******************************************************************************
    ! -- arguments
    class(thdrHdr) :: this
    character(len=*), intent(in) :: fp
    ! -- locals
    integer(I4B) :: iu
    character(len=MXSLEN) :: f, s_r, s
! ------------------------------------------------------------------------------
    !
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'R')
    !
    ! read the first line
    read(iu,'(a)') s_r; rewind(iu)
    s = change_case(s_r, 'l')
    !
    if (s == 'envi') then
      this%i_file_type = i_env
      call this%hdrhdr_read_envi(iu)
    else
      this%i_file_type = i_flt
      call this%hdrhdr_read_ehdr(iu)
    end if
    close(iu)
    !
    return
  end subroutine hdrhdr_read

  subroutine hdrhdr_read_envi(this, iu)
! ******************************************************************************
    ! -- arguments
    class(thdrHdr) :: this
    integer(I4B), intent(in) :: iu
    ! -- locals
    !
    ! parameters
    integer(I4B), parameter :: i_samples                 = 1
    integer(I4B), parameter :: i_lines                   = 2
    integer(I4B), parameter :: i_data_type               = 3
    integer(I4B), parameter :: i_map_info_pixel_easting  = 4
    integer(I4B), parameter :: i_map_info_pixel_northing = 5
    integer(I4B), parameter :: i_map_info_x_pixel_size   = 6
    integer(I4B), parameter :: i_map_info_y_pixel_size   = 7
    integer(I4B), parameter :: i_bin_data                = 8 ! NOT OFFICIAL
    integer(I4B), parameter :: nkey = i_map_info_y_pixel_size
    !
    type(tNum) :: cs, nodata, xul, yul
    !
    logical :: lerr
    !
    character(len=MXSLEN) :: s, k
    character(len=MXSLEN), dimension(:), allocatable :: sa1, sa2
    !
    integer(I4B) :: ios, nrewind, nfound, iu_bin, i, i0, i1
    integer(I4B), dimension(nkey) :: flag
    real(R8B) :: csx, csy
! ------------------------------------------------------------------------------
    !
    call cs%init()
    call nodata%init()
    call xul%init()
    call yul%init()
    !
    flag = 0
    do while(.true.)
      read(unit=iu,iostat=ios,fmt='(a)') s
      if (ios == 0) then
        if (index(s,'=') <= 0) cycle
        call parse_line(s, sa1, token_in='=')
        k = change_case(sa1(1), 'l')
      end if
      !
      select case(k)
      case('samples')
        read(sa1(2),*) this%ncol; flag(i_samples) = 1
      case('lines')
        read(sa1(2),*) this%nrow; flag(i_lines) = 1
      case('data type')
        read(sa1(2),*) i; flag(i_data_type) = 1
        ! see: https://www.l3harrisgeospatial.com/docs/enviheaderfiles.html
        select case(i)
        case(1) ! 1 = Byte: 8-bit unsigned integer
          call errmsg('Not supported: 8-bit unsigned integer')
        case(2) ! 2 = Integer: 16-bit signed integer
          this%nbits = 16; this%i_data_type = i_i2
        case(3) ! 3 = Long: 32-bit signed integer
          this%nbits = 32; this%i_data_type = i_i4
        case(4) ! 4 = Floating-point: 32-bit single-precision
          this%nbits = 32; this%i_data_type = i_r4
        case(5) ! 5 = Double-precision: 64-bit double-precision floating-point
          this%nbits = 64; this%i_data_type = i_r8
        case(6) ! 6 = Complex: Real-imaginary pair of single-precision floating-point
          call errmsg('Not supported: Complex: Real-imaginary pair of single-precision floating-point')
        case(9) ! 9 = Double-precision complex: Real-imaginary pair of double precision floating-point
          call errmsg('Not supported: Double-precision complex: Real-imaginary pair of double precision floating-point')
        case(12) ! 12 = Unsigned integer: 16-bit
          call errmsg('Not supported: Unsigned integer: 16-bit')
        case(13) ! 13 = Unsigned long integer: 32-bit
          call errmsg('Not supported: Unsigned long integer: 32-bit')
        case(14) ! 14 = 64-bit long integer (signed)
          this%nbits = 64; this%i_file_type = i_i8
        case(15) ! 15 = 64-bit unsigned long integer (unsigned)
          call errmsg('Not supported: 64-bit unsigned long integer (unsigned)')
        end select
      case('map info')
        i0 = index(s,'{'); i1 = index(s,'}')
        call parse_line(s(i0+1:i1-1), sa2, token_in=',')
        flag(i_map_info_pixel_easting) = 1
        read(sa2(4),*) this%xllr8 !XLL
        flag(i_map_info_pixel_northing) = 1
        read(sa2(5),*) this%yurr8 !YUR
        flag(i_map_info_x_pixel_size) = 1
        read(sa2(6),*) csx
        flag(i_map_info_y_pixel_size) = 1
        read(sa2(7),*) csy
        if (csx /= csy) then
          call errmsg('x pixel size differs from y pixel size.')
        else
          this%csr8 = csx
        end if
      end select
      !
      nfound = sum(flag)
      if (nfound == nkey) exit
      if (ios /= 0) then
        if (nrewind < nfound) then
          nrewind = nrewind + 1
          rewind(iu)
        else
          if (nfound == nrewind) then
            exit
          end if
        end if
      end if
    end do
    !
    ! optional data
    rewind(iu)
    do while(.true.)
      read(unit=iu,iostat=ios,fmt='(a)') s
      if (ios == 0) then
        if (index(s,'=') <= 0) cycle
        call parse_line(s, sa1, token_in='=')
        k = change_case(sa1(1), 'l')
      end if
      !
      select case(k)
      case('bin_data')
        this%f_bin = sa1(2)
      end select
      if (ios /= 0) then
        exit
      end if
    end do
    !
    ! set the coordinates
    this%xurr8 = this%xllr8 + this%csr8*this%ncol
    this%yllr8 = this%yurr8 - this%csr8*this%nrow
    !
    this%xllr4 = real(this%xllr8,R4B)
    this%xurr4 = real(this%xurr8,R4B)
    this%yllr4 = real(this%yllr8,R4B)
    this%yurr4 = real(this%yurr8,R4B)
    this%csr4  = real(this%csr8,R4B)
    !
    return
  end subroutine hdrhdr_read_envi
  
  subroutine hdrhdr_read_ehdr(this, iu)
! ******************************************************************************
    ! -- arguments
    class(thdrHdr) :: this
    integer(I4B), intent(in) :: iu
    ! -- locals
    !
    ! parameters
    integer(I4B), parameter :: i_byteorder     =  1
    integer(I4B), parameter :: i_layout        =  2
    integer(I4B), parameter :: i_nrows         =  3 
    integer(I4B), parameter :: i_ncols         =  4
    integer(I4B), parameter :: i_nbands        =  5
    integer(I4B), parameter :: i_nbits         =  6
    integer(I4B), parameter :: i_bandrowbytes  =  7
    integer(I4B), parameter :: i_totalrowbytes =  8
    integer(I4B), parameter :: i_pixeltype     =  9
    integer(I4B), parameter :: i_ulxmap        = 10
    integer(I4B), parameter :: i_ulymap        = 11
    integer(I4B), parameter :: i_xllcorner     = 12
    integer(I4B), parameter :: i_yllcorner     = 13
    integer(I4B), parameter :: i_xdim          = 14
    integer(I4B), parameter :: i_ydim          = 15
    integer(I4B), parameter :: i_cellsize      = 16
    integer(I4B), parameter :: i_nodata        = 17
    integer(I4B), parameter :: i_nodata_value  = 18
    integer(I4B), parameter :: nkey = i_nodata_value
    !
    logical :: lerr
    character(len=MXSLEN) :: s, k, v
    character(len=MXSLEN), dimension(:), allocatable :: sa
    integer(I4B) :: ios, nrewind, nfound
    integer(I4B), dimension(nkey) :: flag
    !
    type(tNum) :: cs, nodata, xll, yll, xcul, ycul
! ------------------------------------------------------------------------------
    !
    call cs%init()
    call nodata%init()
    call xll%init()
    call yll%init()
    call xcul%init()
    call ycul%init()
    !
    flag = 0
    nrewind = 0
    do while(.true.)
      read(unit=iu,iostat=ios,fmt='(a)') s
      if (ios == 0) then
        call parse_line(s, sa)
        k = change_case(sa(1), 'l')
        v = sa(2)
        select case(k)
        case('ncols')
          read(v,*) this%ncol; flag(i_ncols) = 1
        case('nrows')
          read(v,*) this%nrow; flag(i_nrows) = 1
        case('xllcorner')
          call xll%set_val_by_s(v); flag(i_xllcorner) = 1; flag(i_ulxmap) = 1
        case('ulxmap') ! center
          call xcul%set_val_by_s(v); flag(i_xllcorner) = 1; flag(i_ulxmap) = 1
        case('yllcorner')
          call yll%set_val_by_s(v); flag(i_yllcorner) = 1; flag(i_ulymap) = 1
        case('ulymap') ! center
          call ycul%set_val_by_s(v); flag(i_yllcorner) = 1; flag(i_ulymap) = 1
        case('nodata','nodata_value')
          call nodata%set_val_by_s(v); flag(i_nodata) = 1; flag(i_nodata_value) = 1
        case('cellsize','xdim', 'ydim')
          call cs%set_val_by_s(v); flag(i_cellsize) = 1; flag(i_xdim) = 1; flag(i_ydim) = 1
        case('nbits')
          read(v,*) this%nbits; flag(i_nbits) = 1
        case('pixeltype')
          read(v,*) this%pixeltype; flag(i_pixeltype) = 1
          this%pixeltype = change_case(this%pixeltype, 'l')
        case('dlt_uscltype')
          v = change_case(v, 'l')
          select case(v)
          case('arith')
            this%i_uscl_type = i_uscl_arith
          case('geom')
            this%i_uscl_type = i_uscl_geom
          case('sumcdr')
            this%i_uscl_type = i_uscl_sumcdr
          case default
            call errmsg('Invalid up-scaling method')
          end select
        case('dlt_dscltype')
          v = change_case(v, 'l')
          select case(v)
          case('nointp')
            this%i_dscl_type = i_dscl_nointp
          case default
            call errmsg('Invalid down-scaling method')
          end select
          
        end select
      end if
      !
      nfound = sum(flag)
      if (nfound == nkey) exit
      if (ios /= 0) then
        if (nrewind < nfound) then
          nrewind = nrewind + 1
          rewind(iu)
        else
          if (nfound == nrewind) then
            exit
          end if
        end if
      end if
    end do
    !
    ! check
    lerr = .false.
    if (flag(i_ncols) == 0)        lerr = .true.
    if (flag(i_nrows) == 0)        lerr = .true.
    if ((flag(i_cellsize) == 0).and.(flag(i_xdim) == 0).and.(flag(i_ydim) == 0)) lerr = .true.
    if ((flag(i_xllcorner) == 0).and.(flag(i_ulxmap) == 0)) lerr = .true.
    if ((flag(i_yllcorner) == 0).and.(flag(i_ulymap) == 0)) lerr = .true.
    if (flag(i_nodata_value) == 0) lerr = .true.
    if (lerr) then
      call logmsg('hdrhdr_read: incomplete header.')
    end if
    !
    ! r4 coordinates
    if (cs%flg(i_r4)) then
      this%csr4 = cs%r4v
    else
      call logmsg('Warning, could not set cs in r4 precision.')
    end if
    if (xll%flg(i_r4)) then
      this%xllr4 = xll%r4v
    end if
    if (xcul%flg(i_r4).and.(this%csr4 /= R4ZERO)) then
!      xll = xcul - cs/2.d0
      this%xllr4 = xcul%r4v - this%csr4 / 2.
    end if
    if (yll%flg(i_r4)) then
      this%yllr4 = yll%r4v
    end if
    if (ycul%flg(i_r4).and.(this%csr4 /= R4ZERO)) then
!      yll = ycul - cs*this%nrow + cs/2.d0
      this%yllr4 = ycul%r4v - this%csr4*this%nrow + this%csr4/2.
    end if
    !
    ! r8 coordinates
    if (cs%flg(i_r8)) then
      this%csr8 = cs%r8v
    else
      call logmsg('Warning, could not set cs in r8 precision.')
    end if
    if (xll%flg(i_r8)) then
      this%xllr8 = xll%r8v
    end if
    if (xcul%flg(i_r8).and.(this%csr8 /= R8ZERO)) then
      this%xllr8 = xcul%r8v - this%csr8 / 2.d0
    end if
    if (yll%flg(i_r8)) then
      this%yllr8 = yll%r8v
    end if
    if (ycul%flg(i_r8).and.(this%csr8 /= R8ZERO)) then
      this%yllr8 = ycul%r8v - this%csr8*this%nrow + this%csr8/2.d0
    end if
    !
    if (nodata%flg(i_i1)) then
      this%mvi1 = nodata%i1v; this%lmvi1 = .true.
    end if
    if (nodata%flg(i_i2)) then
      this%mvi2 = nodata%i2v; this%lmvi2 = .true.
    end if
    if (nodata%flg(i_i4)) then
      this%mvi4 = nodata%i4v; this%lmvi4 = .true.
    end if
    if (nodata%flg(i_r4)) then
      this%mvr4 = nodata%r4v; this%lmvr4 = .true.
    end if
    if (nodata%flg(i_r8)) then
      this%mvr8 = nodata%r8v; this%lmvr8 = .true.
    end if
    !
    this%xurr4 = this%xllr4 + this%ncol*this%csr4
    this%yurr4 = this%yllr4 + this%nrow*this%csr4
    this%xurr8 = this%xllr8 + this%ncol*this%csr8
    this%yurr8 = this%yllr8 + this%nrow*this%csr8
    !
    select case(this%pixeltype)
    case('signedint')
       select case(this%nbits)
       case(8)
         this%i_data_type = i_i1
       case(16)
         this%i_data_type = i_i2
       case(32)
         this%i_data_type = i_i4
       case default
         call errmsg('Error, unsupported nbits for signedint.')
       end select
    case('float')
       select case(this%nbits)
       case(32)
         this%i_data_type = i_r4
       case default
         call errmsg('Error, unsupported nbits for float.')
       end select
    case default
      call errmsg('Error, not supported pixeltype.')
    end select
    !
    return
  end subroutine hdrhdr_read_ehdr
  
  subroutine hdrhdr_get_bbx(this, bbx)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdrHdr) :: this
    type(tBBX), intent(out) :: bbx
    ! -- local
! ------------------------------------------------------------------------------
    bbx%xll = this%xllr8
    bbx%xur = this%xurr8
    bbx%yll = this%yllr8
    bbx%yur = this%yurr8
    bbx%cs  = this%csr8
    !
    return
  end subroutine hdrhdr_get_bbx
  
  subroutine hdrhdr_clip(this, ic0, ic1, ir0, ir1)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdrHdr) :: this
    integer(I4B), intent(in) :: ic0, ic1, ir0, ir1
    ! -- local
! ------------------------------------------------------------------------------
    !
    !call logmsg('Setting header to clipping area')
    this%xllr4 = this%xllr4 + (ic0-1)*this%csr4
    this%xllr8 = this%xllr8 + (ic0-1)*this%csr8
    this%yllr4 = this%yllr4 + (this%nrow-ir1)*this%csr4
    this%yllr8 = this%yllr8 + (this%nrow-ir1)*this%csr8
    this%ncol = ic1 - ic0 + 1
    this%nrow = ir1 - ir0 + 1
    !
    this%xurr4 = this%xllr4 + this%ncol*this%csr4
    this%yurr4 = this%yllr4 + this%nrow*this%csr4
    this%xurr8 = this%xllr8 + this%ncol*this%csr8
    this%yurr8 = this%yllr8 + this%nrow*this%csr8
    !
    return
  end subroutine hdrhdr_clip
  
! ==============================================================================
! ==============================================================================
! tHdr
! ==============================================================================
! ==============================================================================
 
  subroutine hdr_init(this, fp, iu_bin, i_data_type, ext_active, map_active, &
    mapped_all)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    character(len=*), intent(in), optional :: fp
    integer(I4B), intent(in), optional :: iu_bin
    integer(I4B), intent(in), optional :: i_data_type
    logical, intent(in), optional :: ext_active
    logical, intent(in), optional :: map_active
    logical, intent(in), optional :: mapped_all
    ! -- local
! ------------------------------------------------------------------------------
    !
    this%fp           = ''
    this%iu_bin       = 0
    this%i_data_type  = 0
    this%map_active   = .false.
    this%mapped_all   = .false.
    !
    allocate(this%dat)
    call this%dat%init()
    !
    ! overrule the defaults
    if(present(fp))          this%fp          = fp
    if(present(iu_bin))      this%iu_bin      = iu_bin
    if(present(i_data_type)) this%i_data_type = i_data_type
    if(present(ext_active)) then
      this%dat%ext_active  = ext_active
    end if
    if(present(map_active))  this%map_active  = map_active
    if(present(mapped_all))  this%mapped_all  = mapped_all
    !
    return
  end subroutine hdr_init

  subroutine hdr_init_map(this)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    call this%init()
    
    
    write(*,*) 'TODO'
    !
    return
  end subroutine hdr_init_map
  
  subroutine hdr_set_mv(this, mvi1, mvi2, mvi4, mvi8, mvr4, mvr8)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I1B), intent(in), optional :: mvi1
    integer(I2B), intent(in), optional :: mvi2
    integer(I4B), intent(in), optional :: mvi4
    integer(I8B), intent(in), optional :: mvi8
    real(R4B), intent(in), optional    :: mvr4
    real(R8B), intent(in), optional    :: mvr8
    
    ! -- local
    type(thdrHdr), pointer :: hdr
    !
    logical :: ldone
    character(len=MXSLEN) :: s
    integer(I4B) :: ios
    !
    integer(I1B) :: i1v
    integer(I2B) :: i2v
    integer(I4B) :: i4v
    integer(I8B) :: i8v
    real(R4B)    :: r4v
    real(R8B)    :: r8v
! ------------------------------------------------------------------------------
    hdr => this%hdr
    !
    ldone = .true.
    if (present(mvi1)) then
      write(s,*) mvi1; ldone = .false.
    end if
    if (present(mvi2)) then
      write(s,*) mvi2; ldone = .false.
    end if
    if (present(mvi4)) then
      write(s,*) mvi4; ldone = .false.
    end if
    if (present(mvi8)) then
      write(s,*) mvi8; ldone = .false.
    end if
    if (present(mvr4)) then
      write(s,*) mvr4; ldone = .false.
    end if
    if (present(mvr8)) then
      write(s,*) mvr8; ldone = .false.
    end if
    !
    if (.not.ldone) then
      s = adjustl(s)
      read(s,*,iostat=ios) i1v
      if (ios == 0) then 
        hdr%mvi1 = i1v; hdr%lmvi1 = .true.
      end if
      read(s,*,iostat=ios) i2v
      if (ios == 0) then
        hdr%mvi2 = i2v; hdr%lmvi2 = .true.
      end if
      read(s,*,iostat=ios) i4v
      if (ios == 0) then
        hdr%mvi4 = i4v; hdr%lmvi4 = .true.
      end if
      read(s,*,iostat=ios) i8v
      if (ios == 0) then
        hdr%mvi8 = i8v; hdr%lmvi8 = .true.
      end if
      read(s,*,iostat=ios) r4v
      if (ios == 0) then
        hdr%mvr4 = r4v; hdr%lmvr4 = .true.
      end if
      read(s,*,iostat=ios) r8v
      if (ios == 0) then
        hdr%mvr8 = r8v; hdr%lmvr8 = .true.
      end if
    end if
    !
    return
  end subroutine hdr_set_mv

  function hdr_get_file_type(this) result(i_file_type)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I4B) :: i_file_type
    ! -- local
    type(thdrHdr), pointer :: hdr
! ------------------------------------------------------------------------------
    hdr => this%hdr
    i_file_type = hdr%i_file_type
    !
    return
  end function hdr_get_file_type
  ! 
  function hdr_get_data_type(this) result(i_data_type)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I4B) :: i_data_type
    ! -- local
! ------------------------------------------------------------------------------
    i_data_type = this%i_data_type
    !
    return
  end function hdr_get_data_type
  
  function hdr_check_empty(this) result(empty)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    logical :: empty
    ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    integer(I4B) :: irow, icol, nrow, ncol
! ------------------------------------------------------------------------------
    hdr => this%hdr; dat => this%dat
    !
    empty = .true.
    !
    select case(this%i_data_type)
    case(i_i1)
      nrow = size(dat%xi1,2); ncol = size(dat%xi1,1)
      if (.not.hdr%lmvi1) call errmsg('Missing value mvi1.')
      do irow = 1, nrow
        do icol = 1, ncol
          if (dat%xi1(icol,irow) /= hdr%mvi1) then
            empty = .false.
            exit
          end if
        end do
      end do
    case(i_i2)
      nrow = size(dat%xi2,2); ncol = size(dat%xi2,1)
      if (.not.hdr%lmvi2) call errmsg('Missing value mvi2.')
      do irow = 1, nrow
        do icol = 1, ncol
          if (dat%xi2(icol,irow) /= hdr%mvi2) then
            empty = .false.
            exit
          end if
        end do
      end do
    case(i_i4)
      nrow = size(dat%xi4,2); ncol = size(dat%xi4,1)
      if (.not.hdr%lmvi4) call errmsg('Missing value mvi4.')
      do irow = 1, nrow
        do icol = 1, ncol
          if (dat%xi4(icol,irow) /= hdr%mvi4) then
            empty = .false.
            exit
          end if
        end do
      end do
    case(i_i8)
      nrow = size(dat%xi8,2); ncol = size(dat%xi8,1)
      if (.not.hdr%lmvi8) call errmsg('Missing value mvi8.')
      do irow = 1, nrow
        do icol = 1, ncol
          if (dat%xi8(icol,irow) /= hdr%mvi8) then
            empty = .false.
            exit
          end if
        end do
      end do
    case(i_r4)
      nrow = size(dat%xr4,2); ncol = size(dat%xr4,1)
      if (.not.hdr%lmvr4) call errmsg('Missing value mvr4.')
      do irow = 1, nrow
        do icol = 1, ncol
          if (dat%xr4(icol,irow) /= hdr%mvr4) then
            empty = .false.
            exit
          end if
        end do
      end do
    case(i_r8)
      nrow = size(dat%xr8,2); ncol = size(dat%xr8,1)
      if (.not.hdr%lmvr8) call errmsg('Missing value mvr8.')
      do irow = 1, nrow
        do icol = 1, ncol
          if (dat%xr8(icol,irow) /= hdr%mvr8) then
            empty = .false.
            exit
          end if
        end do
      end do
    end select
    !
    return
  end function hdr_check_empty
  
  subroutine hdr_set_i_data_type(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    ! -- local
    type(thdrHdr), pointer :: hdr
! ------------------------------------------------------------------------------
    hdr => this%hdr
    this%i_data_type = hdr%i_data_type
    !
    return
  end subroutine hdr_set_i_data_type
  
  subroutine hdr_clean_dat(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdr) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (associated(this%dat)) then
      call this%dat%clean()
      deallocate(this%dat); this%dat => null()
    end if
    !
    return
  end subroutine hdr_clean_dat
  
  subroutine hdr_clean_dat_buf(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdr) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (associated(this%dat_buf)) then
      call this%dat_buf%clean()
      deallocate(this%dat_buf); this%dat_buf => null()
    end if
    !
    return
  end subroutine hdr_clean_dat_buf
  
  subroutine hdr_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (associated(this%hdr)) then
      call this%hdr%clean()
      deallocate(this%hdr)
    end if
    if (associated(this%hdr_src)) then
      call this%hdr_src%clean()
      deallocate(this%hdr_src)
    end if
    if (associated(this%hdr_no_ext)) deallocate(this%hdr_no_ext)
    !
    this%hdr => null()
    this%hdr_src => null()
    this%hdr_no_ext => null()
    !
    call this%clean_dat()
    call this%clean_dat_buf()
    !
    this%fp           = ''
    this%iu_bin       = 0
    this%i_data_type  = 0
    this%map_active   = .false.
    this%mapped_all   = .false.
    !
    return
  end subroutine hdr_clean
  
  subroutine hdr_read_full_grid(this, fp)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdr) :: this
    character(len=*), intent(inout) :: fp
    ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    character(len=MXSLEN) :: f
    integer(I4B) :: ncol, nrow, iu, icol, irow
    logical :: lop, lalloc
! ------------------------------------------------------------------------------
    !
    call this%init()
    !
    this%fp = strip_ext(fp)
    if (.not.associated(this%hdr)) then
      allocate(this%hdr)
      call this%hdr%init()
      call this%hdr%read(this%fp)
    end if
    hdr => this%hdr; dat => this%dat
    !
    if (hdr%l_read_full_grid) then
      return
    end if
    !
    ncol = hdr%ncol
    nrow = hdr%nrow
    !
    if (this%iu_bin == 0) then
      lop = .false.
    else
      inquire(this%iu_bin, opened=lop)
    end if
    !
    if (.not.lop) then
      if (hdr%i_file_type == i_flt) then
        f = trim(this%fp)//'.flt'
      elseif (hdr%i_file_type == i_env) then
        f = trim(this%fp)
      end if
      call open_file(f, this%iu_bin, 'r', .true.)
    end if
    !
    lalloc = .false.
    call this%set_i_data_type()
    select case(this%i_data_type)
    case(i_i1)
      if (allocated(dat%xi1)) lalloc = .true.
    case(i_i2)
      if (allocated(dat%xi2)) lalloc = .true.
    case(i_i4)
      if (allocated(dat%xi4)) lalloc = .true.
    case(i_i8)
      if (allocated(dat%xi8)) lalloc = .true.
    case(i_r4)
      if (allocated(dat%xr4)) lalloc = .true.
    case(i_r8)
      if (allocated(dat%xr8)) lalloc = .true.
    end select
    !
    if (lalloc) then
      call logmsg('Error: hdr should be cleaned.')
    end if
    !
    select case(this%i_data_type)
    case(i_i1)
      allocate(dat%xi1(ncol,nrow))
      read(this%iu_bin)((dat%xi1(icol,irow),icol=1,ncol),irow=1,nrow)
    case(i_i2)
      allocate(dat%xi2(ncol,nrow))
      read(this%iu_bin)((dat%xi2(icol,irow),icol=1,ncol),irow=1,nrow)
    case(i_i4)
      allocate(dat%xi4(ncol,nrow))
      read(this%iu_bin)((dat%xi4(icol,irow),icol=1,ncol),irow=1,nrow)
    case(i_i8)
      allocate(dat%xi8(ncol,nrow))
      read(this%iu_bin)((dat%xi8(icol,irow),icol=1,ncol),irow=1,nrow)
    case(i_r4)
      allocate(dat%xr4(ncol,nrow))
      read(this%iu_bin)((dat%xr4(icol,irow),icol=1,ncol),irow=1,nrow)
    case(i_r8)
      allocate(dat%xr8(ncol,nrow))
      read(this%iu_bin)((dat%xr8(icol,irow),icol=1,ncol),irow=1,nrow)
    end select
    !
    if (.not.lop) then
      close(this%iu_bin)
      this%iu_bin = -1
    else
      rewind(this%iu_bin)
    end if
    !
    call logmsg('Done reading grid ('//ta(arr=(/ncol,nrow/),sep_in=',')//')...')
    !
    this%hdr%l_read_full_grid = .true.
    !
    return
  end subroutine hdr_read_full_grid

  subroutine hdr_include_ext(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdr) :: this
    ! -- local
    type(tHdrData), pointer :: dat => null()
    logical :: l_n, l_s, l_w, l_e
    integer(I4B) :: nc, nr, mc, mr, ic, ir, r_off, c_off
! ------------------------------------------------------------------------------
    dat => this%dat
    dat%ext_active = .true.
    if (associated(this%hdr_no_ext)) deallocate(this%hdr_no_ext)
    allocate(this%hdr_no_ext)
    this%hdr_no_ext = this%hdr ! copy the values only!
    !
    nc = this%hdr%ncol; nr = this%hdr%nrow
    !
    l_n = .false.; l_s = .false.; l_w = .false.; l_e = .false.
    if (associated(dat%ext_n)) l_n = .true.
    if (associated(dat%ext_s)) l_s = .true.
    if (associated(dat%ext_w)) l_w = .true.
    if (associated(dat%ext_e)) l_e = .true.
    !
    mc = nc; mr = nr
    if (l_n) mr = mr + 1
    if (l_s) mr = mr + 1
    if (l_w) mc = mc + 1
    if (l_e) mc = mc + 1
    !
    ! offset
    r_off = 0; c_off = 0
    if (l_n) r_off = 1; if (l_w) c_off = 1
    !
    ! extend the matrix
    select case(this%i_data_type)
      case(i_i1)
        ! TODO
        call errmsg('i1 not yet supported')
      case(i_i2)
        ! TODO
        call errmsg('i2 not yet supported')
      case(i_i4)
        ! copy interior cells to work array
        nc = size(dat%xi4,1); nr = size(dat%xi4,2)
        allocate(i4wk2d(nc,nr))
        do ir = 1, nr
          do ic = 1, nc
            i4wk2d(ic,ir) = dat%xi4(ic,ir)
          end do
        end do
        deallocate(dat%xi4); allocate(dat%xi4(mc,mr))
        !
        ! set for interior cells
        do ir = 1, nr
          do ic = 1, nc
            dat%xi4(ic+c_off,ir+r_off) = i4wk2d(ic,ir)
          end do
        end do
        deallocate(i4wk2d)
        !
        ! set for exterior cells
        if (l_n) then
          do ic = 1, nc  
            dat%xi4(ic+c_off,1) = dat%ext_n%xi4(ic)
          end do
        end if
        if (l_s) then
          do ic = 1, nc  
            dat%xi4(ic+c_off,mr) = dat%ext_s%xi4(ic)
          end do
        end if
        if (l_w) then
          do ir = 1, nr  
            dat%xi4(1,ir+r_off) = dat%ext_w%xi4(ir)
          end do
        end if
        if (l_e) then
          do ir = 1, nr  
            dat%xi4(mc,ir+r_off) = dat%ext_e%xi4(ir)
          end do
        end if
        if (l_n.and.l_w) then
          dat%xi4(1,1) = dat%ext_nw%xi4(1)
        end if
        if (l_n.and.l_e) then
          dat%xi4(mc,1) = dat%ext_ne%xi4(1)
        end if
        if (l_s.and.l_w) then
          dat%xi4(1,mr) = dat%ext_sw%xi4(1)
        end if
        if (l_s.and.l_e) then
          dat%xi4(mc,mr) = dat%ext_se%xi4(1)
        end if
      case(i_R4)
        ! TODO
        call errmsg('i1 not yet supported')
    end select
    !
    this%hdr%ncol = mc
    this%hdr%nrow = mr
    if (l_w) then
      this%hdr%xllr4 = this%hdr%xllr4 - this%hdr%csr4
      this%hdr%xllr8 = this%hdr%xllr8 - this%hdr%csr8
      this%hdr%ic0 = this%hdr%ic0 - 1
    end if
    if (l_e) then
      this%hdr%ic1 = this%hdr%ic1 + 1
    end if
    if (l_s) then
      this%hdr%yllr4 = this%hdr%yllr4 - this%hdr%csr4
      this%hdr%yllr8 = this%hdr%yllr8 - this%hdr%csr8
      this%hdr%ir1 = this%hdr%ir1 + 1
    end if
    if (l_n) then
      this%hdr%ir0 = this%hdr%ir0 - 1
    end if
    !
    return
  end subroutine hdr_include_ext
  !
  subroutine hdr_exclude_ext(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    ! -- local
    type(tHdrData), pointer :: dat => null()
    logical :: l_n, l_s, l_w, l_e
    integer(I4B) :: nc, nr, mc, mr, ic, ir, r_off, c_off, ic0, ic1, ir0, ir1
! ------------------------------------------------------------------------------
    dat => this%dat
    dat%ext_active = .false.
    if (.not.associated(this%hdr_no_ext)) then
      call errmsg('hdr_exclude_ext: could not find original header')
    end if
    !
    ! restore the original header
    this%hdr = this%hdr_no_ext
    !
    mc = this%hdr%ncol; mr = this%hdr%nrow
    !
    ! offsets
    r_off = 0; c_off = 0
    if (associated(dat%ext_n)) r_off = 1
    if (associated(dat%ext_w)) c_off = 1
    ! extend the matrix
    select case(this%i_data_type)
      case(i_i1)
        nc = size(dat%xi1,1); nr = size(dat%xi1,2)
        allocate(i1wk2d(nc,nr))
        i1wk2d = dat%xi1
        !
        deallocate(dat%xi1); allocate(dat%xi1(mc,mr))
        ic0 = 1 + c_off; ic1 = ic0 + mc - 1
        ir0 = 1 + r_off; ir1 = ir0 + mr - 1
        dat%xi1 = i1wk2d(ic0:ic1,ir0:ir1)
        deallocate(i1wk2d)
      case(i_i2)
        nc = size(dat%xi2,1); nr = size(dat%xi2,2)
        allocate(i2wk2d(nc,nr))
        i2wk2d = dat%xi2
        !
        deallocate(dat%xi2); allocate(dat%xi2(mc,mr))
        ic0 = 1 + c_off; ic1 = ic0 + mc - 1
        ir0 = 1 + r_off; ir1 = ir0 + mr - 1
        dat%xi2 = i2wk2d(ic0:ic1,ir0:ir1)
        deallocate(i2wk2d)
      case(i_i4)
        nc = size(dat%xi4,1); nr = size(dat%xi4,2)
        allocate(i4wk2d(nc,nr))
        i4wk2d = dat%xi4
        !
        deallocate(dat%xi4); allocate(dat%xi4(mc,mr))
        ic0 = 1 + c_off; ic1 = ic0 + mc - 1
        ir0 = 1 + r_off; ir1 = ir0 + mr - 1
        dat%xi4 = i4wk2d(ic0:ic1,ir0:ir1)
        deallocate(i4wk2d)
      case(i_r4)
        nc = size(dat%xr4,1); nr = size(dat%xr4,2)
        allocate(r4wk2d(nc,nr))
        r4wk2d = dat%xr4
        !
        deallocate(dat%xr4); allocate(dat%xr4(mc,mr))
        ic0 = 1 + c_off; ic1 = ic0 + mc - 1
        ir0 = 1 + r_off; ir1 = ir0 + mr - 1
        dat%xr4 = r4wk2d(ic0:ic1,ir0:ir1)
        deallocate(r4wk2d)
    end select
    !
    return
  end subroutine hdr_exclude_ext
    
  subroutine hdr_get_ext_val(this, i_dir, ic, ir, &
    i1v, i2v, i4v, i8v, r4v, r8v)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    !
    integer(I4B), intent(in), optional :: ic
    integer(I4B), intent(in), optional :: ir
    integer(I4B), intent(in) :: i_dir
    !
    integer(I1B), intent(out), optional :: i1v
    integer(I2B), intent(out), optional :: i2v
    integer(I4B), intent(out), optional :: i4v
    integer(I8B), intent(out), optional :: i8v
    real(R4B),    intent(out), optional :: r4v
    real(R8B),    intent(out), optional :: r8v
    !
    ! -- local
    logical :: done
    type(thdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    integer(I4B) :: dt, jc, jr
! ------------------------------------------------------------------------------
    hdr => this%hdr; dat => this%dat
    !
    if ((i_dir == i_n).or.(i_dir == i_s)) then
      if (.not.present(ic)) then
        call errmsg('hdr_get_ext_val: ic not present.')
      end if
    end if
    if ((i_dir == i_w).or.(i_dir == i_e)) then
      if (.not.present(ir)) then
        call errmsg('hdr_get_ext_val: ir not present.')
      end if
    end if
    if ((i_dir == i_nw).or.(i_dir == i_ne).or.&
        (i_dir == i_sw).or.(i_dir == i_se)) then
      if ((.not.present(ic)).or.(.not.present(ir))) then
        call errmsg('hdr_get_ext_val: ic/ir not present.')
      end if
    end if
    !
    dt = this%get_data_type()
    if (present(i1v)) then
      call errmsg('hdr_get_ext_val: i1 not yet supported.')
    end if
    if (present(i2v)) then
      call errmsg('hdr_get_ext_val: i2 not yet supported.')
    end if
    if (present(i4v)) then
      if (dt /= i_i4) call errmsg('hdr_get_ext_val: data is not i4.')
      i4v = this%hdr%mvi4
      if ((i_dir == i_n)  .and. associated(dat%ext_n))  i4v = dat%ext_n%xi4(ic)
      if ((i_dir == i_s)  .and. associated(dat%ext_s))  i4v = dat%ext_s%xi4(ic)
      if ((i_dir == i_w)  .and. associated(dat%ext_w))  i4v = dat%ext_w%xi4(ir)
      if ((i_dir == i_e)  .and. associated(dat%ext_e))  i4v = dat%ext_e%xi4(ir)
      !
      ! corners
      if ((i_dir == i_nw)) then
        if ((ic == 1).and.(ir == 1)) then
          if (associated(dat%ext_nw)) i4v = dat%ext_nw%xi4(1)
        else
          if ((ic == 1       ).and.(associated(dat%ext_w))) i4v = dat%ext_w%xi4(ir-1)
          if ((ir == 1       ).and.(associated(dat%ext_n))) i4v = dat%ext_n%xi4(ic-1)
        end if
      end if
      if ((i_dir == i_ne)) then
        if ((ic == hdr%ncol).and.(ir == 1)) then
          if (associated(dat%ext_ne)) i4v = dat%ext_ne%xi4(1)
        else
          if ((ic == hdr%ncol).and.(associated(dat%ext_e))) i4v = dat%ext_e%xi4(ir-1)
          if ((ir == 1       ).and.(associated(dat%ext_n))) i4v = dat%ext_n%xi4(ic+1)
        end if
      end if
      if ((i_dir == i_sw)) then
        if ((ic == 1).and.(ir == hdr%nrow)) then
          if (associated(dat%ext_sw)) i4v = dat%ext_sw%xi4(1)
        else
          if ((ic == 1)       .and.(associated(dat%ext_w))) i4v = dat%ext_w%xi4(ir+1)
          if ((ir == hdr%nrow).and.(associated(dat%ext_s))) i4v = dat%ext_s%xi4(ic-1)
        end if
      end if
      if ((i_dir == i_se)) then
        if ((ic == hdr%ncol).and.(ir == hdr%nrow)) then
          if (associated(dat%ext_se)) i4v = dat%ext_se%xi4(1)
        else
          if ((ic == hdr%ncol).and.(associated(dat%ext_e))) i4v = dat%ext_e%xi4(ir+1)
          if ((ir == hdr%nrow).and.(associated(dat%ext_s))) i4v = dat%ext_s%xi4(ic+1)
        end if
      end if
    end if
    if (present(i8v)) then
      call errmsg('hdr_get_ext_val: i8 not yet supported.')
    end if
    if (present(r4v)) then
      if (dt /= i_r4) call errmsg('hdr_get_ext_val: data is not r4.')
      r4v = this%hdr%mvr4
      if ((i_dir == i_n)  .and. associated(dat%ext_n))  r4v = dat%ext_n%xr4(ic)
      if ((i_dir == i_s)  .and. associated(dat%ext_s))  r4v = dat%ext_s%xr4(ic)
      if ((i_dir == i_w)  .and. associated(dat%ext_w))  r4v = dat%ext_w%xr4(ir)
      if ((i_dir == i_e)  .and. associated(dat%ext_e))  r4v = dat%ext_e%xr4(ir)
    end if
    if (present(r8v)) then
      call errmsg('hdr_get_ext_val: r8 not yet supported.')
    end if
    !
    return
  end subroutine hdr_get_ext_val
  
  subroutine hdr_get_nbr_val(this, ic, ir, i_dir, &
    i1v, i2v, i4v, i8v, r4v, r8v)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    !
    integer(I4B), intent(in) :: ic
    integer(I4B), intent(in) :: ir
    integer(I4B), intent(in) :: i_dir
    !
    integer(I1B), intent(out), optional :: i1v
    integer(I2B), intent(out), optional :: i2v
    integer(I4B), intent(out), optional :: i4v
    integer(I8B), intent(out), optional :: i8v
    real(R4B),    intent(out), optional :: r4v
    real(R8B),    intent(out), optional :: r8v
     ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    logical :: ok
    integer(I4B) :: jc, jr
    integer(I4B), dimension(:,:), allocatable ::sten
! ------------------------------------------------------------------------------
    hdr => this%hdr; dat => this%dat
    !
    ! checks
    ok = .true.
    if (present(i1v).and.(.not.hdr%lmvi1)) ok = .false.
    if (present(i2v).and.(.not.hdr%lmvi2)) ok = .false.
    if (present(i4v).and.(.not.hdr%lmvi4)) ok = .false.
    if (present(i8v).and.(.not.hdr%lmvi8)) ok = .false.
    if (present(r4v).and.(.not.hdr%lmvr4)) ok = .false.
    if (present(r8v).and.(.not.hdr%lmvr8)) ok = .false.
    if (.not.ok) call errmsg('hdr_get_nbr_val: not matching data types.')
    !
    call get_stencil(ic, ir, hdr%ncol, hdr%nrow, 9, sten)
    jc = sten(1,i_dir); jr = sten(2,i_dir)
    !
    if (jc > 0) then
      if (present(i1v)) i1v = dat%xi1(jc,jr)
      if (present(i2v)) i2v = dat%xi2(jc,jr)
      if (present(i4v)) i4v = dat%xi4(jc,jr)
      if (present(i8v)) i8v = dat%xi8(jc,jr)
      if (present(r4v)) r4v = dat%xr4(jc,jr)
      if (present(r8v)) r8v = dat%xr8(jc,jr)
    else
      if (present(i1v)) call this%get_ext_val(ic=ic, ir=ir, i_dir=i_dir, i1v=i1v)
      if (present(i2v)) call this%get_ext_val(ic=ic, ir=ir, i_dir=i_dir, i2v=i2v)
      if (present(i4v)) call this%get_ext_val(ic=ic, ir=ir, i_dir=i_dir, i4v=i4v)
      if (present(i8v)) call this%get_ext_val(ic=ic, ir=ir, i_dir=i_dir, i8v=i8v)
      if (present(r4v)) call this%get_ext_val(ic=ic, ir=ir, i_dir=i_dir, r4v=r4v)
      if (present(r8v)) call this%get_ext_val(ic=ic, ir=ir, i_dir=i_dir, r8v=r8v)
    end if
    !
    return
  end subroutine hdr_get_nbr_val
    
  function hdr_check_val(this, vi1, vi2, vi4, vr4, vr8) result(found)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    logical :: found
    integer(I1B), intent(in), optional :: vi1
    integer(I2B), intent(in), optional :: vi2
    integer(I4B), intent(in), optional :: vi4
    integer(R4B), intent(in), optional :: vr4
    integer(R8B), intent(in), optional :: vr8
    ! -- local
    type(tHdrData), pointer :: dat => null()
    integer(I4B) :: ir, ic, nr, nc
! ------------------------------------------------------------------------------
    found = .false.
    nr = this%hdr%nrow; nc = this%hdr%ncol
    dat => this%dat
    !
    if(present(vi4)) then
      if(allocated(dat%xi4)) then
        do ir = 1, nr
          do ic = 1, nc
            if (dat%xi4(ic,ir) == vi4) then
              found = .true.
              exit
            end if
          end do
        end do
      end if
    end if
    !
    return
  end function hdr_check_val

  function hdr_get_base_name(this) result(s)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    character(len=MXSLEN) :: s
    ! -- local
    integer(I4B) :: i
    character(len=1) :: slash
! ------------------------------------------------------------------------------
    !
    slash = get_slash()
    i = index(this%fp, slash, back=.true.)
    if (i > 0) then
      s = this%fp(i+1:)
    else
      s = this%fp
    end if
    !
    return
  end function hdr_get_base_name
  
  subroutine hdr_get_global_bb(this, ic0, ic1, ir0, ir1)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I4B), intent(out) :: ic0
    integer(I4B), intent(out) :: ic1
    integer(I4B), intent(out) :: ir0
    integer(I4B), intent(out) :: ir1
 ! ------------------------------------------------------------------------------
    ic0 = this%hdr%ic0
    ic1 = this%hdr%ic1
    ir0 = this%hdr%ir0
    ir1 = this%hdr%ir1
    !
    return
  end subroutine hdr_get_global_bb
  
!  subroutine hdr_read_extent_i4_r8(this, fp, xmin, xmax, ymin, ymax)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    real(R8B), intent(in) :: xmin, xmax, ymin, ymax
!    ! -- local
!    type(tBB) :: bbi
!    integer(I4B) :: ic0, ic1, ir0, ir1
!    type(thdrHdr), pointer :: chdr, mhdr
!    integer(I4B) :: nr, nc
!    real(R8B) :: yul
!! ------------------------------------------------------------------------------
!    !
!    allocate(this%hdr, this%hdr_src)
!    chdr => this%hdr; mhdr => this%hdr_src
!    call chdr%init(); call mhdr%init()
!    call chdr%read(fp); call mhdr%read(fp)
!    !
!    call get_bb_extent(chdr%xllr8, chdr%yllr8, chdr%csr8, &
!      chdr%ncol, chdr%nrow, xmin, xmax, ymin, ymax, &
!      ic0, ic1, ir0, ir1)
!    nc = ic1 - ic0 + 1; nr = ir1 - ir0 + 1
!    !
!    bbi%ic0 = ic0; bbi%ic1 = ic1; bbi%ir0 = ir0; bbi%ir1 = ir1
!    call this%read_clip_grid(fp, bbi)
!    
!    ! set new header data
!    chdr%ncol = nc; chdr%nrow = nr
!    chdr%xllr8 = mhdr%xllr8 + (ic0-1)*chdr%csr8
!    chdr%xllr4 = real(chdr%xllr8,R4B)
!    
!    chdr%yllr8 = mhdr%yllr8 + (mhdr%nrow-ir1)*chdr%csr8
!    chdr%yllr4 = real(chdr%yllr8,R4B)
!    !
!    return
!  end subroutine hdr_read_extent_i4_r8
  
  subroutine hdr_init_read_xy(this, fp, mvr8)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    character(len=*), intent(inout) :: fp
    real(R8B), intent(in), optional :: mvr8
    ! -- local
    type(thdrHdr), pointer :: hdr => null()
    logical :: lop
    character(len=MXSLEN) :: f
! ------------------------------------------------------------------------------
    !
    this%fp = strip_ext(fp)
    if (.not.associated(this%hdr)) then
      allocate(this%hdr)
      call this%hdr%init()
      call this%hdr%read(this%fp)
    end if
    hdr => this%hdr
    !
    allocate(this%dat)
    call this%dat%init()
    !
    if (present(mvr8)) then
      if (hdr%i_file_type /= i_env) then
        call errmsg('hdr_init_read_xy: setting mv only allowed for ENVI.')
      end if
      hdr%mvr8 = mvr8
      hdr%lmvr8 = .true.
    end if
    !
    if (this%iu_bin == 0) then
      lop = .false.
    else
      inquire(this%iu_bin, opened=lop)
    end if
    !
    if (.not.lop) then
      if (hdr%i_file_type == i_flt) then
        f = trim(this%fp)//'.flt'
      elseif (hdr%i_file_type == i_env) then
        f = trim(this%fp)
      end if
      call open_file(f, this%iu_bin, 'r', .true.)
    end if
    !
    return
  end subroutine hdr_init_read_xy
  
  subroutine hdr_read_extent(this, fp, dst_bbx, i_uscl, i_dscl)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    character(len=*), intent(inout) :: fp
    type(tBBX), intent(in) :: dst_bbx
    integer(I4B), intent(in), optional :: i_uscl
    integer(I4B), intent(in), optional :: i_dscl
    ! -- local
    type(tBB) :: bbi
    integer(I4B) :: ic0, ic1, ir0, ir1, ius, ids
    type(thdrHdr), pointer :: dst_hdr, src_hdr
    integer(I4B) :: nr, nc
    real(R8B) :: yul
! ------------------------------------------------------------------------------
    !
    call this%init(fp=fp)
    !
    if (associated(this%hdr)) then
      call this%hdr%clean()
      deallocate(this%hdr)
    end if
    if (associated(this%hdr_src)) then
      call this%hdr_src%clean()
      deallocate(this%hdr_src)
    end if
    allocate(this%hdr, this%hdr_src)
    dst_hdr => this%hdr; src_hdr => this%hdr_src
    call src_hdr%init(); call dst_hdr%init()
    call src_hdr%read(fp); call dst_hdr%read(fp)
    !
    call get_bb_extent(src_hdr%xllr8, src_hdr%yllr8, src_hdr%csr8, &
      src_hdr%ncol, src_hdr%nrow, &
      dst_bbx%xll, dst_bbx%xur, dst_bbx%yll, dst_bbx%yur, &
      ic0, ic1, ir0, ir1)
    nc = ic1 - ic0 + 1; nr = ir1 - ir0 + 1
    !
    bbi%ic0 = ic0; bbi%ic1 = ic1; bbi%ir0 = ir0; bbi%ir1 = ir1
    if (associated(this%dat_buf)) then
      call logmsg('Re-using data...')
      call this%clean_dat()
      allocate(this%dat)
      call this%dat_buf%copy(this%dat)
      call this%set_i_data_type()
    else
      call this%read_clip_grid(fp, bbi)
    end if
    call this%hdr%clip(ic0, ic1, ir0, ir1)
    !
    if (present(i_uscl)) then
      ius = i_uscl
    else
      ius = i_uscl_arith
    end if
    !
    if (present(i_dscl)) then
      ids = i_dscl
    else
      ids = i_dscl_nointp
    end if
    if (dst_bbx%cs == src_hdr%csr8) then ! no scaling
      call logmsg('---> No scaling applied <---')
    else if (dst_bbx%cs > src_hdr%csr8) then ! upscale
      call this%up_scale(ius, dst_bbx%cs)
    else ! downscaling
      select case(ids)
      case(i_dscl_nointp)
        call this%down_scale_nointp(ids, dst_bbx%cs)
      case(i_dscl_sample_scale)
        call this%down_scale_nointp(ids, dst_bbx%cs, correct_for_area=.true.)
      case(i_dscl_intp)
        call this%down_scale_intp(ids, dst_bbx%cs)
      case(i_dscl_nodata)
        call errmsg('hdr_read_extent: invalid downscaling option.')
      end select
    end if
    !
    return
  end subroutine hdr_read_extent
  
  subroutine hdr_up_scale(this, i_uscl, dst_cs)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I4B), intent(in) :: i_uscl
    real(R8B), intent(in) :: dst_cs
    ! -- local
    type(tHdrHdr), pointer :: src_hdr => null()
    type(tHdrData), pointer :: dat => null()
    !
    integer(I4B), dimension(:,:), allocatable :: xi4
    real(R4B),    dimension(:,:), allocatable :: xr4
    !
    integer(I4B) :: i4v, i4s
    real(R4B)    :: r4v, r4s
    !
    real(R8B) :: es
    !
    integer(I4B) :: ns, mc, mr, ic, ir, jc, jr
    integer(I4B) :: ic0, ic1, ir0, ir1, n
! ------------------------------------------------------------------------------
    !
    select case(i_uscl)
    case(i_uscl_arith, i_uscl_nodata)
      call logmsg('---> Upscaling: artihmetic <---')
    case(i_uscl_geom)
      call logmsg('---> Upscaling: geometric <---')
    case(i_uscl_sumcdr)
      call logmsg('---> Upscaling: sum conductance <---')
    case default
      call errmsg('Not supported upscaling method.')
    end select
    !
    src_hdr => this%hdr; dat => this%dat
    ns = dst_cs/src_hdr%csr8
    mc = src_hdr%ncol/ns; mr = src_hdr%nrow/ns
    !
    select case(this%i_data_type)
    case(i_i1)
      call errmsg('hdr_up_scale: i1 not yet supported')
    case(i_i2)
      call errmsg('hdr_up_scale: i2 not yet supported')
    case(i_i4)
      allocate(xi4(mc,mr)); xi4 = src_hdr%mvi4
      do jr = 1, mr; do jc = 1, mc
        ir0 = (jr-1)*ns + 1; ir1 = ir0 + ns - 1
        ic0 = (jc-1)*ns + 1; ic1 = ic0 + ns - 1
        n = 0; i4s = 0; es = R8ZERO
        do ir = ir0, ir1; do ic = ic0, ic1
          i4v = dat%xi4(ic,ir)
          if (i4v /= src_hdr%mvi4) then
            n =  n + 1; i4s = i4s + i4v; es = es + log(real(i4v,R8B))
          end if
        end do; end do
        !
        select case(i_uscl)
        case(i_uscl_arith, i_uscl_nodata)
          if (n > 0) then
            xi4(jc,jr) = i4s/n
          end if
        case(i_uscl_geom)
         if (n > 0) then
           xi4(jc,jr) = int(exp(es/n),i4b)
         end if
        case(i_uscl_sumcdr)
          xi4(jc,jr) = i4s
        end select
      end do; end do
      !
      call this%replace_grid(xi4=xi4, mvi4=src_hdr%mvi4)
    case(i_r4)
      allocate(xr4(mc,mr)); xr4 = src_hdr%mvr4
      do jr = 1, mr; do jc = 1, mc
        ir0 = (jr-1)*ns + 1; ir1 = ir0 + ns - 1
        ic0 = (jc-1)*ns + 1; ic1 = ic0 + ns - 1
        n = 0; r4s = 0; es = R8ZERO
        do ir = ir0, ir1; do ic = ic0, ic1
          r4v = dat%xr4(ic,ir)
          if (r4v /= src_hdr%mvr4) then
            n =  n + 1; r4s = r4s + r4v; es = es + log(real(r4v,R8B))
          end if
        end do; end do
        !
        select case(i_uscl)
        case(i_uscl_arith, i_uscl_nodata)
          if (n > 0) then
            xr4(jc,jr) = r4s/n
          end if
        case(i_uscl_geom)
          if (n > 0) then
            xr4(jc,jr) = real(exp(es/n),r4b)
          end if  
        case(i_uscl_sumcdr)
          xr4(jc,jr) = r4s
        end select
      end do; end do
      !
      call this%replace_grid(xr4=xr4, mvr4=src_hdr%mvr4)
    case(i_r8)
      call errmsg('hdr_up_scale: r8 not yet supported')
    end select
    !
    ! update the header
    src_hdr%ncol = mc; src_hdr%nrow = mr
    src_hdr%csr8 = dst_cs
    src_hdr%csr4 = real(src_hdr%csr8,R4B)
    !
    return
  end subroutine hdr_up_scale
  
  subroutine hdr_down_scale_nointp(this, i_dscl, dst_cs, correct_for_area)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdr) :: this
    integer(I4B), intent(in) :: i_dscl
    real(R8B), intent(in) :: dst_cs
    logical, intent(in), optional :: correct_for_area
    ! -- local
    type(thdrHdr), pointer :: src_hdr => null()
    type(tHdrData), pointer :: dat => null()
    !
    logical :: cfa
    integer(I4B), dimension(:,:), allocatable :: xi4
    real(R4B),    dimension(:,:), allocatable :: xr4
    real(R8B),    dimension(:,:), allocatable :: xr8
    !
    integer(I4B) :: ns, mc, mr, jr, jc, ir, ic
    real(R8B) :: cfr8
! ------------------------------------------------------------------------------
    if (present(correct_for_area)) then
      cfa = correct_for_area
    else
      cfa = .false.
    end if
    !
    src_hdr => this%hdr; dat => this%dat
    ns = src_hdr%csr8/dst_cs
    mc = src_hdr%ncol*ns; mr = src_hdr%nrow*ns
    !
    if (cfa) then
      call logmsg('---> Down-scaling: sampling with area correction <---')
      cfr8 = R8ONE/(ns**2)
    else
      call logmsg('---> Down-scaling: sampling <---')
      cfr8 = R8ONE
    end if
    !
    select case(this%i_data_type)
    case(i_i1)
      call errmsg('hdr_down_scale_nointp: i1 not yet supported')
    case(i_i2)
      call errmsg('hdr_down_scale_nointp: i2 not yet supported')
    case(i_i4)
      allocate(xi4(mc,mr)); xi4 = src_hdr%mvi4
      do jr = 1, mr; do jc = 1, mc
        ir = ceiling(real(jr,R8B)/ns); ic = ceiling(real(jc,R8B)/ns)
        if (dat%xi4(ic,ir) /= src_hdr%mvi4) then
          xi4(jc,jr) = dat%xi4(ic,ir) * int(cfr8,I4B)
        end if
      end do; end do
      call this%replace_grid(xi4=xi4, mvi4=src_hdr%mvi4)
    case(i_r4)
      allocate(xr4(mc,mr)); xr4 = src_hdr%mvr4
      do jr = 1, mr; do jc = 1, mc
        ir = ceiling(real(jr,R8B)/ns); ic = ceiling(real(jc,R8B)/ns)
        if (dat%xr4(ic,ir) /= src_hdr%mvr4) then
          xr4(jc,jr) = dat%xr4(ic,ir) * real(cfr8,r4B)
        end if
      end do; end do
      call this%replace_grid(xr4=xr4, mvr4=src_hdr%mvr4)
    case(i_r8)
      allocate(xr8(mc,mr)); xr8 = src_hdr%mvr8
      do jr = 1, mr; do jc = 1, mc
        ir = ceiling(real(jr,R8B)/ns); ic = ceiling(real(jc,R8B)/ns)
        if (dat%xr8(ic,ir) /= src_hdr%mvr8) then
          xr8(jc,jr) = dat%xr8(ic,ir) * cfr8
        end if
      end do; end do
      call this%replace_grid(xr8=xr8, mvr8=src_hdr%mvr8)
    end select
    !
    ! update the header
    src_hdr%ncol = mc; src_hdr%nrow = mr
    src_hdr%csr8 = dst_cs
    src_hdr%csr4 = real(src_hdr%csr8,R4B)
    !
    return
  end subroutine hdr_down_scale_nointp
  
  subroutine hdr_down_scale_intp(this, i_dscl, dst_cs)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I4B), intent(in) :: i_dscl
    real(R8B), intent(in) :: dst_cs
    ! -- local
    type(tHdrHdr), pointer :: src_hdr
    type(tHdrData), pointer :: dat => null()
    !
    integer(I4B), dimension(:,:), allocatable :: xi4
    real(R4B),    dimension(:,:), allocatable :: xr4
    !
    integer(I4B) :: i4v
    real(R4B) :: r4v
    !
    real(R8B), dimension(2,1) :: x_src, x_dst
    real(R8B), dimension(2,1) :: xll
    real(R8B), dimension(2,1) :: xur
    real(R8B), dimension(4)   :: va
    real(R8B)                 :: mv, v, xs, ys, xd, yd, src_cs
    !
    integer(I4B) :: ns, mc, mr, jr, jc, ir, ic
    !
    integer(I4B) :: ic1, ir1, ic2, ir2, ic3, ir3, ic4, ir4
! ------------------------------------------------------------------------------
    !
    call logmsg('---> Down-scaling: bi-linear interpolation <---')
    !
    src_hdr => this%hdr; dat => this%dat
    ns = src_hdr%csr8/dst_cs
    mc = src_hdr%ncol*ns; mr = src_hdr%nrow*ns
    src_cs = src_hdr%csr8
    !
    select case(this%i_data_type)
    case(i_i1)
      call errmsg('hdr_down_scale_intp: i1 not yet supported')
    case(i_i2)
      call errmsg('hdr_down_scale_intp: i2 not yet supported')
    case(i_i4)
      call errmsg('hdr_down_scale_intp: i4 not yet supported')
      !allocate(xi4(mc,mr)); xi4 = src_hdr%mvi4
      !mv = real(src_hdr%mvi4,R8B)
      !do jr = 1, mr; do jc = 1, mc
      !  ir = ceiling(real(jr,R8B)/ns); ic = ceiling(real(jc,R8B)/ns)
      !  !
      !  ! linear interpolation
      !  call get_xy(xs, ys, ic, ir, src_hdr%xllr8, src_hdr%yurr8, src_cs)
      !  call get_xy(xd, yd, jc, jr, src_hdr%xllr8, src_hdr%yurr8, dst_cs)
      !  x_dst = reshape([xd, yd], shape(x_dst))
      !  !
      !  ! check
      !  if ((xd == xs).and.(yd == ys)) then
      !    call errmsg('hdr_down_scale_intp: program error.')
      !  end if
      !  !
      !  va = mv
      !  if (xd < xs) then ! west
      !    if (yd > ys) then ! north
      !      xll = reshape([xs - src_cs, ys], shape(xll))
      !      xur = reshape([xs, ys + src_cs], shape(xur))
      !      call this%get_nbr_val(ic, ir, i_w,  i4v=i4v); va(1) = real(i4v,R8B)
      !      call this%get_nbr_val(ic, ir, i_nw, i4v=i4v); va(2) = real(i4v,R8B)
      !      call this%get_nbr_val(ic, ir, i_n,  i4v=i4v); va(3) = real(i4v,R8B)
      !      va(4) = real(this%xi4(ic,ir),R8B)
      !    else ! south
      !      xll = reshape([xs - src_cs, ys - src_cs], shape(xll))
      !      xur = reshape([xs, ys], shape(xur))
      !      call this%get_nbr_val(ic, ir, i_sw, i4v=i4v); va(1) = real(i4v,R8B)
      !      call this%get_nbr_val(ic, ir, i_w,  i4v=i4v); va(2) = real(i4v,R8B)
      !      va(3) = real(this%xi4(ic,ir),R8B)
      !      call this%get_nbr_val(ic, ir, i_s,  i4v=i4v); va(4) = real(i4v,R8B)
      !    end if
      !  else ! east
      !    if (yd > ys) then ! north
      !      xll = reshape([xs, ys], shape(xll))
      !      xur = reshape([xs + src_cs, ys + src_cs], shape(xur))
      !      va(1) = real(this%xi4(ic,ir),R8B)
      !      call this%get_nbr_val(ic, ir, i_n,  i4v=i4v); va(2) = real(i4v,R8B)
      !      call this%get_nbr_val(ic, ir, i_ne, i4v=i4v); va(3) = real(i4v,R8B)
      !      call this%get_nbr_val(ic, ir, i_e,  i4v=i4v); va(4) = real(i4v,R8B)
      !    else ! south
      !      xll = reshape([xs, ys - src_cs], shape(xll))
      !      xur = reshape([xs + src_cs, ys], shape(xur))
      !      call this%get_nbr_val(ic, ir, i_s,  i4v=i4v); va(1) = real(i4v,R8B)
      !      va(2) = real(this%xi4(ic,ir),R8B)
      !      call this%get_nbr_val(ic, ir, i_ne, i4v=i4v); va(3) = real(i4v,R8B)
      !      call this%get_nbr_val(ic, ir, i_se, i4v=i4v); va(4) = real(i4v,R8B)
      !    end if
      !  end if
      !  !
      !  v = bilinear_interpolation(x_dst, xll, xur, va, mv)
      !  xi4(jc,jr) = int(v,I4B)
      !end do; end do
      !call this%replace_grid(xi4=xi4, mvi4=src_hdr%mvi4)
    case(i_i8)
      call errmsg('hdr_down_scale_intp: i8 not yet supported')
    case(i_r4)
      allocate(xr4(mc,mr)); xr4 = src_hdr%mvr4
      mv = real(src_hdr%mvr4,R8B)
      do jr = 1, mr; do jc = 1, mc
        ir = ceiling(real(jr,R8B)/ns); ic = ceiling(real(jc,R8B)/ns)
        !
        ! linear interpolation
        call get_xy(xs, ys, ic, ir, src_hdr%xllr8, src_hdr%yurr8, src_cs)
        call get_xy(xd, yd, jc, jr, src_hdr%xllr8, src_hdr%yurr8, dst_cs)
        x_dst = reshape([xd, yd], shape(x_dst))
        !
        ! check
        if ((xd == xs).and.(yd == ys)) then
          call errmsg('hdr_down_scale_intp: program error.')
        end if
        !
        va = mv
        if (xd < xs) then ! west
          if (yd > ys) then ! north
            xll = reshape([xs - src_cs, ys], shape(xll))
            xur = reshape([xs, ys + src_cs], shape(xur))
            call this%get_nbr_val(ic, ir, i_w,  r4v=r4v); va(1) = real(r4v,R8B)
            call this%get_nbr_val(ic, ir, i_nw, r4v=r4v); va(2) = real(r4v,R8B)
            call this%get_nbr_val(ic, ir, i_n,  r4v=r4v); va(3) = real(r4v,R8B)
            va(4) = real(dat%xr4(ic,ir),R8B)
          else ! south
            xll = reshape([xs - src_cs, ys - src_cs], shape(xll))
            xur = reshape([xs, ys], shape(xur))
            call this%get_nbr_val(ic, ir, i_sw, r4v=r4v); va(1) = real(r4v,R8B)
            call this%get_nbr_val(ic, ir, i_w,  r4v=r4v); va(2) = real(r4v,R8B)
            va(3) = real(dat%xr4(ic,ir),R8B)
            call this%get_nbr_val(ic, ir, i_s,  r4v=r4v); va(4) = real(r4v,R8B)
          end if
        else ! east
          if (yd > ys) then ! north
            xll = reshape([xs, ys], shape(xll))
            xur = reshape([xs + src_cs, ys + src_cs], shape(xur))
            va(1) = real(dat%xr4(ic,ir),R8B)
            call this%get_nbr_val(ic, ir, i_n,  r4v=r4v); va(2) = real(r4v,R8B)
            call this%get_nbr_val(ic, ir, i_ne, r4v=r4v); va(3) = real(r4v,R8B)
            call this%get_nbr_val(ic, ir, i_e,  r4v=r4v); va(4) = real(r4v,R8B)
          else ! south
            xll = reshape([xs, ys - src_cs], shape(xll))
            xur = reshape([xs + src_cs, ys], shape(xur))
            call this%get_nbr_val(ic, ir, i_s,  r4v=r4v); va(1) = real(r4v,R8B)
            va(2) = real(dat%xr4(ic,ir),R8B)
            call this%get_nbr_val(ic, ir, i_ne, r4v=r4v); va(3) = real(r4v,R8B)
            call this%get_nbr_val(ic, ir, i_se, r4v=r4v); va(4) = real(r4v,R8B)
          end if
        end if
        !
        v = bilinear_interpolation(x_dst, xll, xur, va, mv)
        xr4(jc,jr) = real(v,R4B)
      end do; end do
      call this%replace_grid(xr4=xr4, mvr4=src_hdr%mvr4)
    case(i_r8)
      call errmsg('hdr_down_scale_intp: r8 not yet supported')
    end select
    !
    ! update the header
    src_hdr%ncol = mc; src_hdr%nrow = mr
    src_hdr%csr8 = dst_cs
    src_hdr%csr4 = real(src_hdr%csr8,R4B)
    !
    return
  end subroutine hdr_down_scale_intp
    
  subroutine hdr_read_clip_grid(this, fp, bbi)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    character(len=*), intent(inout) :: fp
    type(tBB), intent(in) :: bbi
    ! -- local
    type(thdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    logical :: cliphdr, l_n, l_s, l_w, l_e
    character(len=MXSLEN) :: f
    integer(I4B) :: ic0, ic1, ir0, ir1
    integer(I4B) :: nc, nr, mc, mr, ic, ir, jc, jr, n
! ------------------------------------------------------------------------------
    !
    if (.not.associated(this%hdr)) then
      allocate(this%hdr)
      call this%hdr%init()
      call this%hdr%read(fp)
    end if
    hdr => this%hdr; dat => this%dat
    !
    ic0 = bbi%ic0
    ic1 = bbi%ic1
    ir0 = bbi%ir0
    ir1 = bbi%ir1
    !
    hdr%ic0 = ic0
    hdr%ic1 = ic1
    hdr%ir0 = ir0
    hdr%ir1 = ir1
    !
    nc = ic1 - ic0 + 1; nr = ir1 - ir0 + 1
    mc = hdr%ncol; mr = hdr%nrow
    !
    l_n = .false.; l_s = .false.; l_w = .false.; l_e = .false.
    if (ir0 > 1)  l_n = .true.
    if (ir1 < mr) l_s = .true.
    if (ic0 > 1)  l_w = .true.
    if (ic1 < mc) l_e = .true.
    !
    if (hdr%i_file_type == i_flt) then
      f = trim(fp)//'.flt'
    else
      f = trim(fp)
    end if
    call open_file(f, this%iu_bin, 'R', .true.)
    !
    call this%set_i_data_type()
    select case(this%i_data_type)
    case(i_i1)
      call errmsg('hdr_read_clip_grid: i1 not yet supported')
    case(i_i2)
      call errmsg('hdr_read_clip_grid: i2 not yet supported')
    case(i_i4)
      !
      ! read the boundaries
      if (l_n) then
        allocate(dat%ext_n); allocate(dat%ext_n%xi4(nc))
        call this%read_arr(icol=ic0, irow=ir0-1, i4a=dat%ext_n%xi4)
      end if
      ! read the boundaries
      if (l_s) then
        allocate(dat%ext_s); allocate(dat%ext_s%xi4(nc))
        call this%read_arr(icol=ic0, irow=ir1+1, i4a=dat%ext_s%xi4)
      end if
      if (l_w) then
        allocate(dat%ext_w); allocate(dat%ext_w%xi4(nr))
        do ir = ir0, ir1
          call this%read_val(icol=ic0-1, irow=ir, i4v=dat%ext_w%xi4(ir-ir0+1))
        end do
      end if
      if (l_e) then
        allocate(dat%ext_e); allocate(dat%ext_e%xi4(nr))
        do ir = ir0, ir1
          call this%read_val(icol=ic1+1, irow=ir, i4v=dat%ext_e%xi4(ir-ir0+1))
        end do
      end if
      if (l_n.and.l_w) then
        allocate(dat%ext_nw); allocate(dat%ext_nw%xi4(1))
        call this%read_val(icol=ic0-1, irow=ir0-1, i4v=dat%ext_nw%xi4(1))
      end if
      if (l_n.and.l_e) then
        allocate(dat%ext_ne); allocate(dat%ext_ne%xi4(1))
        call this%read_val(icol=ic1+1, irow=ir0-1, i4v=dat%ext_ne%xi4(1))
      end if
      if (l_s.and.l_w) then
        allocate(dat%ext_sw); allocate(dat%ext_sw%xi4(1))
        call this%read_val(icol=ic0-1, irow=ir1+1, i4v=dat%ext_sw%xi4(1))
      end if
      if (l_s.and.l_e) then
        allocate(dat%ext_se); allocate(dat%ext_se%xi4(1))
        call this%read_val(icol=ic1+1, irow=ir1+1, i4v=dat%ext_se%xi4(1))
      end if
      !
      ! read the block
      allocate(dat%xi4(nc,nr))
      do ir = 1, nr
        jr = ir0 + ir - 1; jc = ic0
        call this%read_arr(icol=jc, irow=jr, i4a=dat%xi4(:,ir))
      end do
    case(i_i8)
      call errmsg('hdr_read_clip_grid: i8 not yet supported')
    case(i_r4)
      ! read the boundaries
      if (l_n) then
        allocate(dat%ext_n); allocate(dat%ext_n%xr4(nc))
        call this%read_arr(icol=ic0, irow=ir0-1, r4a=dat%ext_n%xr4)
      end if
      ! read the boundaries
      if (l_s) then
        allocate(dat%ext_s); allocate(dat%ext_s%xr4(nc))
        call this%read_arr(icol=ic0, irow=ir1+1, r4a=dat%ext_s%xr4)
      end if
      if (l_w) then
        allocate(dat%ext_w); allocate(dat%ext_w%xr4(nr))
        do ir = ir0, ir1
          call this%read_val(icol=ic0-1, irow=ir, r4v=dat%ext_w%xr4(ir-ir0+1))
        end do
      end if
      if (l_e) then
        allocate(dat%ext_e); allocate(dat%ext_e%xr4(nr))
        do ir = ir0, ir1
          call this%read_val(icol=ic1+1, irow=ir, r4v=dat%ext_e%xr4(ir-ir0+1))
        end do
      end if
      if (l_n.and.l_w) then
        allocate(dat%ext_nw); allocate(dat%ext_nw%xr4(1))
        call this%read_val(icol=ic0-1, irow=ir0-1, r4v=dat%ext_nw%xr4(1))
      end if
      if (l_n.and.l_e) then
        allocate(dat%ext_ne); allocate(dat%ext_ne%xr4(1))
        call this%read_val(icol=ic1+1, irow=ir0-1, r4v=dat%ext_ne%xr4(1))
      end if
      if (l_s.and.l_w) then
        allocate(dat%ext_sw); allocate(dat%ext_sw%xr4(1))
        call this%read_val(icol=ic0-1, irow=ir1+1, r4v=dat%ext_sw%xr4(1))
      end if
      if (l_s.and.l_e) then
        allocate(dat%ext_se); allocate(dat%ext_se%xr4(1))
        call this%read_val(icol=ic1+1, irow=ir1+1, r4v=dat%ext_se%xr4(1))
      end if
      !
      ! read the block
      allocate(dat%xr4(nc,nr))
      do ir = 1, nr
        jr = ir0 + ir - 1; jc = ic0
        call this%read_arr(icol=jc, irow=jr, r4a=dat%xr4(:,ir))
      end do
    case(i_r8)
      ! read the boundaries
      if (l_n) then
        allocate(dat%ext_n); allocate(dat%ext_n%xr8(nc))
        call this%read_arr(icol=ic0, irow=ir0-1, r8a=dat%ext_n%xr8)
      end if
      ! read the boundaries
      if (l_s) then
        allocate(dat%ext_s); allocate(dat%ext_s%xr8(nc))
        call this%read_arr(icol=ic0, irow=ir1+1, r8a=dat%ext_s%xr8)
      end if
      if (l_w) then
        allocate(dat%ext_w); allocate(dat%ext_w%xr8(nr))
        do ir = ir0, ir1
          call this%read_val(icol=ic0-1, irow=ir, r8v=dat%ext_w%xr8(ir-ir0+1))
        end do
      end if
      if (l_e) then
        allocate(dat%ext_e); allocate(dat%ext_e%xr8(nr))
        do ir = ir0, ir1
          call this%read_val(icol=ic1+1, irow=ir, r8v=dat%ext_e%xr8(ir-ir0+1))
        end do
      end if
      if (l_n.and.l_w) then
        allocate(dat%ext_nw); allocate(dat%ext_nw%xr8(1))
        call this%read_val(icol=ic0-1, irow=ir0-1, r8v=dat%ext_nw%xr8(1))
      end if
      if (l_n.and.l_e) then
        allocate(dat%ext_ne); allocate(dat%ext_ne%xr8(1))
        call this%read_val(icol=ic1+1, irow=ir0-1, r8v=dat%ext_ne%xr8(1))
      end if
      if (l_s.and.l_w) then
        allocate(dat%ext_sw); allocate(dat%ext_sw%xr8(1))
        call this%read_val(icol=ic0-1, irow=ir1+1, r8v=dat%ext_sw%xr8(1))
      end if
      if (l_s.and.l_e) then
        allocate(dat%ext_se); allocate(dat%ext_se%xr8(1))
        call this%read_val(icol=ic1+1, irow=ir1+1, r8v=dat%ext_se%xr8(1))
      end if
      !
      ! read the block
      allocate(dat%xr8(nc,nr))
      do ir = 1, nr
        jr = ir0 + ir - 1; jc = ic0
        call this%read_arr(icol=jc, irow=jr, r8a=dat%xr8(:,ir))
      end do
    end select
    close(this%iu_bin)
    !
    return
  end subroutine hdr_read_clip_grid
  
  subroutine hdr_map_grid(child, prent, mvi4)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: child
    type(thdr), intent(in) :: prent
    integer(I4B), intent(in), optional :: mvi4
    ! -- local
    logical :: l_in, l_prent_data, l_child_data, l_set
    type(thdrHdr), pointer :: chdr => null()
    type(thdrHdr), pointer :: phdr => null()
    integer(I4B) :: ic, ir, jc, jr, ic0, ic1, ir0, ir1
    integer(I4B) :: i4v_prent, i4v_child
    real(R8B) :: x, y !relative x and y
! ------------------------------------------------------------------------------
    !
    chdr => child%hdr; phdr => prent%hdr
    !
    ! check for downscaling
    if (chdr%csr8 < phdr%csr8) then
      call errmsg('hdr_map_grid: downscaling is not (yet)supported.')
    end if
    !
    ! check for data types
    if (child%i_data_type /= prent%i_data_type) then
      call errmsg('hdr_map_grid: different types are not (yet) supported.')
    end if
    !
    select case(child%i_data_type)
    case(i_i4)
      !if (.not.phdr%lmvi4) then
      !  call errmsg('hdr_map_grid: mvi4 of parent not set.')
      !end if
      if (.not.present(mvi4)) then
        call errmsg('hdr_map_grid: mvi4 not set.')
      end if
      do ir = 1, chdr%nrow
        do ic = 1, chdr%ncol
          call get_xy(x, y, ic, ir, chdr%xllr8, chdr%yurr8, chdr%csr8)
          call get_icir(jc, jr, x, y, phdr%xllr8, phdr%xurr8, &
            phdr%yllr8, phdr%yurr8, phdr%csr8, &
            phdr%ncol, phdr%nrow, l_in)
          if (l_in) then
            i4v_prent = prent%dat%xi4(jc,jr)
            i4v_child = child%dat%xi4(ic,ir)
            if (i4v_prent /= phdr%mvi4) then ! data value
              l_prent_data = .true.
            else
              l_prent_data = .false.
            end if
            if ((i4v_child /= chdr%mvi4).and.(i4v_child /= mvi4)) then ! data value set
              l_child_data = .true.
              if (l_prent_data) then
                if (i4v_prent /= i4v_child) then
                  call logmsg('hdr_map_grid: WARNING; cell data is already set.')
                  call logmsg('(x,y) = '//ta(arr=[real(x,R4B),real(y,R4B)],sep_in=','))
                  call logmsg('values: '//ta(arr=[i4v_prent,i4v_child],sep_in=','))
                end if
              end if 
            else
              l_child_data = .false.
            end if
            !
            l_set = .false.
            if (i4v_child == mvi4) then
              l_set = .true.
            else
              if (l_prent_data.and.(.not.l_child_data)) then
                l_set = .true.
              end if
            end if
            !
            if (l_set) then
              child%dat%xi4(ic,ir) = i4v_prent
            end if
          end if
        end do
      end do
    case default
      call errmsg('hdr_map_grid: type is not (yet) supported.')
    end select
    !
    return
  end subroutine hdr_map_grid
  
  subroutine hdr_read_grid(this, fp, xi1, xi2, xi4, xi8, xr4, xr8, &
    mvi1, mvi2, mvi4, mvi8, mvr4, mvr8, &
    xllr4, yllr4, csr4, xllr8, yllr8, csr8)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdr) :: this
    !
    character(len=*), intent(inout) :: fp
    integer(I1B), dimension(:,:), intent(out), allocatable, optional :: xi1
    integer(I2B), dimension(:,:), intent(out), allocatable, optional :: xi2
    integer(I4B), dimension(:,:), intent(out), allocatable, optional :: xi4
    integer(I8B), dimension(:,:), intent(out), allocatable, optional :: xi8
    real(R4B),    dimension(:,:), intent(out), allocatable, optional :: xr4
    real(R8B),    dimension(:,:), intent(out), allocatable, optional :: xr8
    !
    integer(I1B), intent(out), optional :: mvi1
    integer(I2B), intent(out), optional :: mvi2
    integer(I4B), intent(out), optional :: mvi4
    integer(I8B), intent(out), optional :: mvi8
    real(R4B),    intent(out), optional :: mvr4
    real(R8B),    intent(out), optional :: mvr8
    !
    real(R4B), intent(out), optional :: xllr4
    real(R4B), intent(out), optional :: yllr4
    real(R4B), intent(out), optional :: csr4
    !
    real(R8B), intent(out), optional :: xllr8
    real(R8B), intent(out), optional :: yllr8
    real(R8B), intent(out), optional :: csr8
    !
    ! -- local
    type(tHdrHdr), pointer :: hdr
! ------------------------------------------------------------------------------
    !
    call this%read_full_grid(fp)
    hdr => this%hdr
    !
    if (present(xi1)) then
      if (.not.present(mvi1)) call errmsg('hdr_read_grid: mvi1 missing.')
      call this%get_grid(xi1=xi1, mvi1=mvi1)
    end if
    if (present(xi2)) then
      if (.not.present(mvi2)) call errmsg('hdr_read_grid: mvi2 missing.')
      call this%get_grid(xi2=xi2, mvi2=mvi2)
    end if
    if (present(xi4)) then
      if (.not.present(mvi4)) call errmsg('hdr_read_grid: mvi4 missing.')
      call this%get_grid(xi4=xi4, mvi4=mvi4)
    end if
    if (present(xi8)) then
      if (.not.present(mvi8)) call errmsg('hdr_read_grid: mvi8 missing.')
      call this%get_grid(xi8=xi8, mvi8=mvi8)
    end if
    if (present(xr4)) then
      if (.not.present(mvr4)) call errmsg('hdr_read_grid: mvr4 missing.')
      call this%get_grid(xr4=xr4, mvr4=mvr4)
    end if
    if (present(xr8)) then
      if (.not.present(mvr8)) call errmsg('hdr_read_grid: mvr8 missing.')
      call this%get_grid(xr8=xr8, mvr8=mvr8)
    end if
    !
    if (present(xllr4)) xllr4 = hdr%xllr4
    if (present(yllr4)) yllr4 = hdr%yllr4
    if (present(csr4))  csr4  = hdr%csr4
    !
    if (present(xllr8)) xllr8 = hdr%xllr8
    if (present(yllr8)) yllr8 = hdr%yllr8
    if (present(csr8))  csr8  = hdr%csr8
    !
    call this%clean()
    !
    return
  end subroutine hdr_read_grid
  
  subroutine hdr_get_grid(this, xi1, xi2, xi4, xi8, xr4, xr8, &
    mvi1, mvi2, mvi4, mvi8, mvr4, mvr8)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdr) :: this
    integer(I1B), dimension(:,:), intent(out), allocatable, optional :: xi1
    integer(I2B), dimension(:,:), intent(out), allocatable, optional :: xi2
    integer(I4B), dimension(:,:), intent(out), allocatable, optional :: xi4
    integer(I8B), dimension(:,:), intent(out), allocatable, optional :: xi8
    real(R4B),    dimension(:,:), intent(out), allocatable, optional :: xr4
    real(R8B),    dimension(:,:), intent(out), allocatable, optional :: xr8
    !
    integer(I1B), intent(out), optional :: mvi1
    integer(I2B), intent(out), optional :: mvi2
    integer(I4B), intent(out), optional :: mvi4
    integer(I8B), intent(out), optional :: mvi8
    real(R4B),    intent(out), optional :: mvr4
    real(R8B),    intent(out), optional :: mvr8
    !
    ! -- local
    type(tHdrHdr), pointer :: hdr
    type(tHdrData), pointer :: dat
    logical :: flg
    integer(I4B) :: ncol, nrow, icol, irow
! ------------------------------------------------------------------------------
    !
    hdr => this%hdr; dat => this%dat
    ncol = hdr%ncol; nrow = hdr%nrow
    !
    flg = .false.
    if (present(xi1)) then
      if (allocated(xi1)) deallocate(xi1); allocate(xi1(ncol,nrow))
      if (this%i_data_type /= i_i1) flg = .true.
    end if
    if (present(xi2)) then
      if (allocated(xi2)) deallocate(xi2); allocate(xi2(ncol,nrow))
      if (this%i_data_type /= i_i2) flg = .true.
    end if
    if (present(xi4)) then
      if (allocated(xi4)) deallocate(xi4); allocate(xi4(ncol,nrow))
      if (this%i_data_type /= i_i4) flg = .true.
    end if
    if (present(xi8)) then
      if (allocated(xi8)) deallocate(xi8); allocate(xi8(ncol,nrow))
      if (this%i_data_type /= i_i8) flg = .true.
    end if
    if (present(xr4)) then
      if (allocated(xr4)) deallocate(xr4); allocate(xr4(ncol,nrow))
      if (this%i_data_type /= i_r4) flg = .true.
    end if
    if (present(xr8)) then
      if (allocated(xr8)) deallocate(xr8); allocate(xr8(ncol,nrow))
      if (this%i_data_type /= i_r8) flg = .true.
    end if
    
    if (flg) then
      call logmsg('Warning, original headered file has different precision.')
    end if
    !
    select case(this%i_data_type)
    case(i_i1)
      call errmsg('Not yet supported.')
    case(i_i4)
      if (present(xi1)) call cast_2darray(xi4_in=dat%xi4, mvi4_in=hdr%mvi4, xi1_out=xi1, mvi1_out=mvi1)
      if (present(xi2)) call cast_2darray(xi4_in=dat%xi4, mvi4_in=hdr%mvi4, xi2_out=xi2, mvi2_out=mvi2)
      if (present(xi4)) call cast_2darray(xi4_in=dat%xi4, mvi4_in=hdr%mvi4, xi4_out=xi4, mvi4_out=mvi4)
      if (present(xi8)) call cast_2darray(xi4_in=dat%xi4, mvi4_in=hdr%mvi4, xi8_out=xi8, mvi8_out=mvi8)
      if (present(xr4)) call cast_2darray(xi4_in=dat%xi4, mvi4_in=hdr%mvi4, xr4_out=xr4, mvr4_out=mvr4)
      if (present(xr8)) call cast_2darray(xi4_in=dat%xi4, mvi4_in=hdr%mvi4, xr8_out=xr8, mvr8_out=mvr8)
    case(i_r4)
      if (present(xi1)) call cast_2darray(xr4_in=dat%xr4, mvr4_in=hdr%mvr4, xi1_out=xi1, mvi1_out=mvi1)
      if (present(xi2)) call cast_2darray(xr4_in=dat%xr4, mvr4_in=hdr%mvr4, xi2_out=xi2, mvi2_out=mvi2)
      if (present(xi4)) call cast_2darray(xr4_in=dat%xr4, mvr4_in=hdr%mvr4, xi4_out=xi4, mvi4_out=mvi4)
      if (present(xi8)) call cast_2darray(xr4_in=dat%xr4, mvr4_in=hdr%mvr4, xi8_out=xi8, mvi8_out=mvi8)
      if (present(xr4)) call cast_2darray(xr4_in=dat%xr4, mvr4_in=hdr%mvr4, xr4_out=xr4, mvr4_out=mvr4)
      if (present(xr8)) call cast_2darray(xr4_in=dat%xr4, mvr4_in=hdr%mvr4, xr8_out=xr8, mvr8_out=mvr8)
    !  pass
    !case(i_r8)
    !  pass
    end select
    !
    return
    end subroutine hdr_get_grid
  
!  subroutine hdr_read_grid_i1_r4(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    integer(I1B), dimension(:,:), intent(out), allocatable :: x
!    integer(I1B), intent(out) :: mv
!    real(R4B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    mv = hdr%mvi1; xll = hdr%xllr4; yll = hdr%yllr4; cs = hdr%csr4
!    call this%hdr_get_grid_i1(x, mv)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_i1_r4
!
!  subroutine hdr_read_grid_i1_r8(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    integer(I1B), dimension(:,:), intent(out), allocatable :: x
!    integer(I1B), intent(out) :: mv
!    real(R8B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    xll = hdr%xllr8; yll = hdr%yllr8; cs = hdr%csr8
!    call this%hdr_get_grid_i1(x, mv)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_i1_r8
!  
!  subroutine hdr_read_grid_i2_r4(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    integer(I2B), dimension(:,:), intent(out), allocatable :: x
!    integer(I2B), intent(out) :: mv
!    real(R4B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    xll = hdr%xllr4; yll = hdr%yllr4; cs = hdr%csr4
!    call this%hdr_get_grid_i2(x, mv)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_i2_r4
!
!  subroutine hdr_read_grid_i2_r8(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    integer(I2B), dimension(:,:), intent(out), allocatable :: x
!    integer(I2B), intent(out) :: mv
!    real(R8B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    xll = hdr%xllr8; yll = hdr%yllr8; cs = hdr%csr8
!    call this%hdr_get_grid_i2(x, mv)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_i2_r8
!  
!  subroutine hdr_read_grid_i4_r4(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    integer(I4B), dimension(:,:), intent(out), allocatable :: x
!    integer(I4B), intent(out) :: mv
!    real(R4B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    xll = hdr%xllr4; yll = hdr%yllr4; cs = hdr%csr4
!    call this%hdr_get_grid_i4(x, mv)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_i4_r4
!
!  subroutine hdr_read_grid_i4_r8(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    integer(I4B), dimension(:,:), intent(out), allocatable :: x
!    integer(I4B), intent(out) :: mv
!    real(R8B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    xll = hdr%xllr8; yll = hdr%yllr8; cs = hdr%csr8
!    call this%hdr_get_grid_i4(x, mv)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_i4_r8
!  
!  subroutine hdr_read_grid_r4_r4(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    real(R4B), dimension(:,:), intent(out), allocatable :: x
!    real(R4B), intent(out) :: mv
!    real(R4B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    xll = hdr%xllr4; yll = hdr%yllr4; cs = hdr%csr4
!    call this%hdr_get_grid_r4(x, mv=mv, bs=1)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_r4_r4
!
!  subroutine hdr_read_grid_r4_r8(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    real(R4B), dimension(:,:), intent(inout), allocatable :: x
!    real(R4B), intent(inout) :: mv
!    real(R8B), intent(inout) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: ic0, ir0, ic1, ir1, bs
!! ------------------------------------------------------------------------------
!    !
!    allocate(this%hdr)
!    call this%hdr%init()
!    call this%hdr%read(fp)
!    !
!    if (allocated(x)) then
!      call logmsg('*** Clipping and scaling ***')
!      call this%hdr_uscl_read_check_r8(xll, yll, cs, size(x,1), size(x,2), &
!        ic0, ic1, ir0, ir1, bs)
!      call this%read_clip_grid(fp, ic0, ic1, ir0, ir1)
!    else
!      call this%read_full_grid(fp)
!      hdr => this%hdr
!      xll = hdr%xllr8; yll = hdr%yllr8; cs = hdr%csr8
!      bs = 1
!    end if
!    call this%hdr_get_grid_r4(x, mv=mv, bs=bs)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_r4_r8
!  
!  subroutine hdr_read_grid_r8_r4(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    real(R8B), dimension(:,:), intent(out), allocatable :: x
!    real(R8B), intent(out) :: mv
!    real(R4B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    xll = hdr%xllr4; yll = hdr%yllr4; cs = hdr%csr4
!    call this%hdr_get_grid_r8(x, mv=mv, bs=1)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_r8_r4
!
!  subroutine hdr_read_grid_r8_r8(this, fp, x, mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    real(R8B), dimension(:,:), intent(out), allocatable :: x
!    real(R8B), intent(out) :: mv
!    real(R8B), intent(out) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!! ------------------------------------------------------------------------------
!    !
!    call this%read_full_grid(fp)
!    hdr => this%hdr
!    !
!    xll = hdr%xllr8; yll = hdr%yllr8; cs = hdr%csr8
!    call this%hdr_get_grid_r8(x, mv=mv, bs=1)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_grid_r8_r8
!  
!  subroutine hdr_read_block_i4_r8(this, fp, ir0, ir1, ic0, ic1, x, &
!    mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    integer(I4B), intent(in) :: ir0
!    integer(I4B), intent(in) :: ir1
!    integer(I4B), intent(in) :: ic0
!    integer(I4B), intent(in) :: ic1
!    real(R4B), dimension(:,:), intent(inout), allocatable :: x
!    real(R4B), intent(inout) :: mv
!    real(R8B), intent(inout) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: nr, nc
!    
!! ------------------------------------------------------------------------------
!    !
!    allocate(this%hdr)
!    call this%hdr%init()
!    call this%hdr%read(fp)
!    !
!    if (allocated(x)) deallocate(x)
!    nr = ir1 - ir0 + 1; nc = ic1 - ic0 + 1
!    allocate(x(nc,nr))
!    !
!    call this%hdr_read_clip_grid(fp, ic0, ic1, ir0, ir1)
!    mv = this%hdr%mvr8; xll = this%hdr%xllr8; yll = this%hdr%yllr8; cs = this%hdr%csr8
!    call this%hdr_get_grid_r4(x, 1)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_block_i4_r8

!  subroutine hdr_read_block_i4_r8(this, fp, ir0, ir1, ic0, ic1, x, &
!    mv, xll, yll, cs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    character(len=*), intent(inout) :: fp
!    integer(I4B), intent(in) :: ir0
!    integer(I4B), intent(in) :: ir1
!    integer(I4B), intent(in) :: ic0
!    integer(I4B), intent(in) :: ic1
!    real(R4B), dimension(:,:), intent(inout), allocatable :: x
!    real(R4B), intent(inout) :: mv
!    real(R8B), intent(inout) :: xll, yll, cs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: nr, nc
!    
!! ------------------------------------------------------------------------------
!    !
!    allocate(this%hdr)
!    call this%hdr%init()
!    call this%hdr%read(fp)
!    !
!    if (allocated(x)) deallocate(x)
!    nr = ir1 - ir0 + 1; nc = ic1 - ic0 + 1
!    allocate(x(nc,nr))
!    !
!    call this%read_clip_grid(fp, ic0, ic1, ir0, ir1)
!    xll = this%hdr%xllr8; yll = this%hdr%yllr8; cs = this%hdr%csr8
!    call this%hdr_get_grid_r4(x, bs=1, mv=mv)
!    call this%clean()
!    !
!    return
!  end subroutine hdr_read_block_i4_r8
    
    
!  subroutine hdr_get_grid_i1(this, x, mv)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    integer(I1B), dimension(:,:), intent(out), allocatable :: x
!    integer(I1B), intent(out) :: mv
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: ncol, nrow, icol, irow
!! ------------------------------------------------------------------------------
!    !
!    hdr => this%hdr
!    ncol = hdr%ncol; nrow = hdr%nrow
!    if (allocated(x)) deallocate(x)
!    allocate(x(ncol,nrow))
!    if (this%i_data_type /= i_i1) then
!      call logmsg('Warning, original hdr file has different precision.')
!    end if
!    !
!    if (.not.hdr%lmvi1) call errmsg('Missing value mvi1.')
!    mv = hdr%mvi1
!    select case(this%i_data_type)
!    case(i_i1)
!      if (.not.hdr%lmvi1) call errmsg('Missing value mvi1.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi1(icol,irow) /= hdr%mvi1) then
!            x(icol,irow) = int(this%xi1(icol,irow),I1B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i2)
!    if (.not.hdr%lmvi2) call errmsg('Missing value mvi2.')      
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi2(icol,irow) /= hdr%mvi2) then
!            x(icol,irow) = int(this%xi2(icol,irow),I1B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i4)
!    if (.not.hdr%lmvi4) call errmsg('Missing value mvi4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi4(icol,irow) /= hdr%mvi4) then
!            x(icol,irow) = int(this%xi4(icol,irow),I1B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i8)
!    if (.not.hdr%lmvi8) call errmsg('Missing value mvi8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi8(icol,irow) /= hdr%mvi8) then
!            x(icol,irow) = int(this%xi8(icol,irow),I1B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r4)
!      if (.not.hdr%lmvr4) call errmsg('Missing value mvr4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xr4(icol,irow) /= hdr%mvr4) then
!            x(icol,irow) = int(this%xr4(icol,irow),I1B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r8)
!      if (.not.hdr%lmvr8) call errmsg('Missing value mvr8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xr8(icol,irow) /= hdr%mvr8) then
!            x(icol,irow) = int(this%xr8(icol,irow),I1B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    end select
!    !
!    return
!  end subroutine hdr_get_grid_i1
!  
!  subroutine hdr_get_grid_i2(this, x, mv)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    integer(I2B), dimension(:,:), intent(out), allocatable :: x
!    integer(I2B), intent(out) :: mv
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: ncol, nrow, icol, irow
!! ------------------------------------------------------------------------------
!    !
!    hdr => this%hdr
!    ncol = hdr%ncol; nrow = hdr%nrow
!    if (allocated(x)) deallocate(x)
!    allocate(x(ncol,nrow))
!    if (this%i_data_type /= i_i2) then
!      call logmsg('Warning, original hdr file has different precision.')
!    end if
!    !
!    if (.not.hdr%lmvi2) call errmsg('Missing value mvi2.')
!    mv = hdr%mvi2
!    select case(this%i_data_type)
!    case(i_i1)
!      if (.not.hdr%lmvi1) call errmsg('Missing value mvi1.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi1(icol,irow) /= hdr%mvi1) then
!            x(icol,irow) = int(this%xi1(icol,irow),I2B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i2)
!      if (.not.hdr%lmvi2) call errmsg('Missing value mvi2.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi2(icol,irow) /= hdr%mvi2) then
!            x(icol,irow) = int(this%xi2(icol,irow),I2B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i4)
!      if (.not.hdr%lmvi4) call errmsg('Missing value mvi4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi4(icol,irow) /= hdr%mvi4) then
!            x(icol,irow) = int(this%xi4(icol,irow),I2B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i8)
!      if (.not.hdr%lmvi8) call errmsg('Missing value mvi8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi8(icol,irow) /= hdr%mvi8) then
!            x(icol,irow) = int(this%xi8(icol,irow),I2B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r4)
!      if (.not.hdr%lmvr4) call errmsg('Missing value mvr4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xr4(icol,irow) /= hdr%mvr4) then
!            x(icol,irow) = int(this%xr4(icol,irow),I2B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r8)
!      if (.not.hdr%lmvr8) call errmsg('Missing value mvr8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xr8(icol,irow) /= hdr%mvr8) then
!            x(icol,irow) = int(this%xr8(icol,irow),I2B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    end select
!    !
!    return
!  end subroutine hdr_get_grid_i2
!  
!  subroutine hdr_get_grid_i4(this, x, mv)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    integer(I4B), dimension(:,:), intent(out), allocatable :: x
!    integer(I4B), intent(out) :: mv
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: ncol, nrow, icol, irow
!! ------------------------------------------------------------------------------
!    !
!    hdr => this%hdr
!    ncol = hdr%ncol; nrow = hdr%nrow
!    if (allocated(x)) deallocate(x)
!    allocate(x(ncol,nrow))
!    if (this%i_data_type /= i_i4) then
!      call logmsg('Warning, original hdr file has different precision.')
!    end if
!    !
!    if (.not.hdr%lmvi4) call errmsg('Missing value mvi4.')
!    mv = hdr%mvi4
!    select case(this%i_data_type)
!    case(i_i1)
!      if (.not.hdr%lmvi1) call errmsg('Missing value mvi1.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi1(icol,irow) /= hdr%mvi1) then
!            x(icol,irow) = int(this%xi1(icol,irow),I4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i2)
!      if (.not.hdr%lmvi2) call errmsg('Missing value mvi2.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi2(icol,irow) /= hdr%mvi2) then
!            x(icol,irow) = int(this%xi2(icol,irow),I4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i4)
!      if (.not.hdr%lmvi4) call errmsg('Missing value mvi4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi4(icol,irow) /= hdr%mvi4) then
!            x(icol,irow) = int(this%xi4(icol,irow),I4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i8)
!      if (.not.hdr%lmvi8) call errmsg('Missing value mvi8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi8(icol,irow) /= hdr%mvi8) then
!            x(icol,irow) = int(this%xi8(icol,irow),I4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r4)
!      if (.not.hdr%lmvr4) call errmsg('Missing value mvr4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xr4(icol,irow) /= hdr%mvr4) then
!            x(icol,irow) = int(this%xr4(icol,irow),I4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r8)
!      if (.not.hdr%lmvr8) call errmsg('Missing value mvr8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xr8(icol,irow) /= hdr%mvr8) then
!            x(icol,irow) = int(this%xr8(icol,irow),I4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    end select
!    !
!    return
!  end subroutine hdr_get_grid_i4
!  
!  subroutine hdr_get_grid_i8(this, x, mv)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    integer(I8B), dimension(:,:), intent(out), allocatable :: x
!    integer(I8B), intent(out) :: mv
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: ncol, nrow, icol, irow
!! ------------------------------------------------------------------------------
!    !
!    hdr => this%hdr
!    ncol = hdr%ncol; nrow = hdr%nrow
!    if (allocated(x)) deallocate(x)
!    allocate(x(ncol,nrow))
!    if (this%i_data_type /= i_i8) then
!      call logmsg('Warning, original hdr file has different precision.')
!    end if
!    !
!    if (.not.hdr%lmvi8) call errmsg('Missing value mvi8.')
!    mv = hdr%mvi8
!    select case(this%i_data_type)
!    case(i_i1)
!      if (.not.hdr%lmvi1) call errmsg('Missing value mvi1.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi1(icol,irow) /= hdr%mvi1) then
!            x(icol,irow) = int(this%xi1(icol,irow),I8B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i2)
!      if (.not.hdr%lmvi2) call errmsg('Missing value mvi2.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi2(icol,irow) /= hdr%mvi2) then
!            x(icol,irow) = int(this%xi2(icol,irow),I8B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i4)
!      if (.not.hdr%lmvi4) call errmsg('Missing value mvi4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi4(icol,irow) /= hdr%mvi4) then
!            x(icol,irow) = int(this%xi4(icol,irow),I8B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i8)
!      if (.not.hdr%lmvi8) call errmsg('Missing value mvi8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi8(icol,irow) /= hdr%mvi8) then
!            x(icol,irow) = int(this%xi8(icol,irow),I8B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r4)
!      if (.not.hdr%lmvr4) call errmsg('Missing value mvr4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xr4(icol,irow) /= hdr%mvr4) then
!            x(icol,irow) = int(this%xr4(icol,irow),I8B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r8)
!      if (.not.hdr%lmvr8) call errmsg('Missing value mvr8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xr8(icol,irow) /= hdr%mvr8) then
!            x(icol,irow) = int(this%xr8(icol,irow),I8B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    end select
!    !
!    return
!  end subroutine hdr_get_grid_i8
!  
!  subroutine hdr_get_grid_r4(this, x, mv, bs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    real(R4B), dimension(:,:), intent(inout), allocatable :: x
!    real(R4B), intent(out) :: mv
!    integer(I4B), intent(in) :: bs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: ncol, nrow, icol, irow, jcol, jrow, ir0, ir1, ic0, ic1, n
!    real(R4B) :: r4v
!    real(R8B) :: r8v
!    logical :: lscale
!! ------------------------------------------------------------------------------
!    !
!    hdr => this%hdr
!    if (.not.hdr%lmvr4) call errmsg('Missing value mvr4.')
!    mv = hdr%mvr4
!    !
!    if (allocated(x)) then
!      lscale = .true.
!      ncol = size(x,1); nrow = size(x,2)
!    else
!      lscale = .false.
!      ncol = hdr%ncol; nrow = hdr%nrow
!      allocate(x(ncol,nrow))
!    end if
!    mv = hdr%mvr4
!    do irow = 1, nrow
!      do icol = 1, ncol
!        x(icol,irow) = mv
!      end do
!    end do
!    if (this%i_data_type /= i_r4) then
!      call logmsg('Warning, original hdr file has different precision.')
!    end if
!    !
!    select case(this%i_data_type)
!    case(i_i1)
!      if (.not.hdr%lmvi1) call errmsg('Missing value mvi1.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi1(icol,irow) /= hdr%mvi1) then
!            x(icol,irow) = real(this%xi1(icol,irow),R4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i2)
!      if (.not.hdr%lmvi2) call errmsg('Missing value mvi2.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi2(icol,irow) /= hdr%mvi2) then
!            x(icol,irow) = real(this%xi2(icol,irow),R4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i4)
!      if (.not.hdr%lmvi4) call errmsg('Missing value mvi4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi4(icol,irow) /= hdr%mvi4) then
!            x(icol,irow) = real(this%xi4(icol,irow),R4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i8)
!      if (.not.hdr%lmvi8) call errmsg('Missing value mvi8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi8(icol,irow) /= hdr%mvi8) then
!            x(icol,irow) = real(this%xi4(icol,irow),R4B)
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r4)
!      if (.not.hdr%lmvr4) call errmsg('Missing value mvr4.')
!      do irow = 1, nrow ! target grid
!        do icol = 1, ncol ! target grid
!          ir0 = (irow-1)*bs + 1; ir1 = ir0 + bs - 1
!          ic0 = (icol-1)*bs + 1; ic1 = ic0 + bs - 1
!          !
!          n = 0; r4v = R4ZERO; x(icol,irow) = mv
!          do jrow = ir0, ir1
!            do jcol = ic0, ic1
!              if (this%xr4(jcol,jrow) /= hdr%mvr4) then
!                if (hdr%i_uscl_type == i_uscl_geom) then
!                  r4v = r4v + log(this%xr4(icol,irow))
!                else
!                  r4v = r4v + this%xr4(icol,irow)
!                end if
!                n = n + 1
!              end if
!            end do
!          end do
!          select case(hdr%i_uscl_type)
!            case(i_uscl_arith, i_uscl_nodata)
!              if (n > 0) then
!                r4v = r4v / n
!              end if
!            case(i_uscl_geom)
!              if (n > 0) then
!                r4v = exp(r4v/n)
!              end if
!            case(i_uscl_sumcdr)
!              ! nothing
!          end select
!          x(icol,irow) = real(r4v,R4B)
!        end do
!      end do
!    case(i_r8)
!      if (.not.hdr%lmvr8) call errmsg('Missing value mvr8.')
!      do irow = 1, nrow ! target grid
!        do icol = 1, ncol ! target grid
!          ir0 = (irow-1)*bs + 1; ir1 = ir0 + bs - 1
!          ic0 = (icol-1)*bs + 1; ic1 = ic0 + bs - 1
!          !
!          n = 0; r8v = R8ZERO; x(icol,irow) = mv
!          do jrow = ir0, ir1
!            do jcol = ic0, ic1
!              if (this%xr8(jcol,jrow) /= hdr%mvr8) then
!                if (hdr%i_uscl_type == i_uscl_geom) then
!                  r8v = r8v + log(this%xr8(icol,irow))
!                else
!                  r8v = r8v + this%xr8(icol,irow)
!                end if
!                n = n + 1
!              end if
!            end do
!          end do
!          select case(hdr%i_uscl_type)
!            case(i_uscl_arith, i_uscl_nodata)
!              if (n > 0) then
!                r8v = r8v / n
!              end if
!            case(i_uscl_geom)
!              if (n > 0) then
!                r8v = exp(r8v/n)
!              end if
!            case(i_uscl_sumcdr)
!              ! nothing
!          end select
!          x(icol,irow) = real(r8v,R4B)
!        end do
!      end do
!    end select
!    !
!    return
!  end subroutine hdr_get_grid_r4
!  
!  subroutine hdr_get_grid_r8(this, x, mv, bs)
!! ******************************************************************************
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(thdr) :: this
!    real(R8B), dimension(:,:), intent(out), allocatable :: x
!    real(R8B), intent(out) :: mv
!    integer(I4B), intent(in) :: bs
!    ! -- local
!    type(thdrHdr), pointer :: hdr
!    integer(I4B) :: ncol, nrow, icol, irow, jcol, jrow, ir0, ir1, ic0, ic1, n
!    real(R4B) :: r4v
!    real(R8B) :: r8v
!! ------------------------------------------------------------------------------
!    !
!    hdr => this%hdr
!    ncol = hdr%ncol; nrow = hdr%nrow
!    if (allocated(x)) deallocate(x)
!    allocate(x(ncol,nrow))
!    if (this%i_data_type /= i_r8) then
!      call logmsg('Warning, original hdr file has different precision.')
!    end if
!    !
!    if (.not.hdr%lmvr8) call errmsg('Missing value mvr8.')
!    mv = hdr%mvr8 !!!
!    select case(this%i_data_type)
!    case(i_i1)
!      if (.not.hdr%lmvi1) call errmsg('Missing value mvi1.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi1(icol,irow) /= hdr%mvi1) then
!            x(icol,irow) = real(this%xi1(icol,irow),R8B) !!!
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i2)
!      if (.not.hdr%lmvi2) call errmsg('Missing value mvi2.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi2(icol,irow) /= hdr%mvi2) then
!            x(icol,irow) = real(this%xi2(icol,irow),R8B) !!!
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i4)
!      if (.not.hdr%lmvi4) call errmsg('Missing value mvi4.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi4(icol,irow) /= hdr%mvi4) then
!            x(icol,irow) = real(this%xi4(icol,irow),R8B) !!!
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_i8)
!      if (.not.hdr%lmvi8) call errmsg('Missing value mvi8.')
!      do irow = 1, nrow
!        do icol = 1, ncol
!          if (this%xi8(icol,irow) /= hdr%mvi8) then
!            x(icol,irow) = real(this%xi8(icol,irow),R8B) !!!
!          else
!            x(icol,irow) = mv
!          end if
!        end do
!      end do
!    case(i_r4)
!      if (.not.hdr%lmvr4) call errmsg('Missing value mvr4.')
!      do irow = 1, nrow ! target grid
!        do icol = 1, ncol ! target grid
!          ir0 = (irow-1)*bs + 1; ir1 = ir0 + bs - 1
!          ic0 = (icol-1)*bs + 1; ic1 = ic0 + bs - 1
!          !
!          n = 0; r4v = R4ZERO; x(icol,irow) = mv
!          do jrow = ir0, ir1
!            do jcol = ic0, ic1
!              if (this%xr4(jcol,jrow) /= hdr%mvr4) then
!                if (hdr%i_uscl_type == i_uscl_geom) then
!                  r4v = r4v + log(this%xr4(icol,irow))
!                else
!                  r4v = r4v + this%xr4(icol,irow)
!                end if
!                n = n + 1
!              end if
!            end do
!          end do
!          select case(hdr%i_uscl_type)
!            case(i_uscl_arith, i_uscl_nodata)
!              if (n > 0) then
!                r4v = r4v / n
!              end if
!            case(i_uscl_geom)
!              if (n > 0) then
!                r4v = exp(r4v) / n
!              end if
!            case(i_uscl_sumcdr)
!              ! nothing
!          end select
!          x(icol,irow) = real(r4v, R8B)
!        end do
!      end do
!    case(i_r8)
!      if (.not.hdr%lmvr8) call errmsg('Missing value mvr8.')
!      do irow = 1, nrow ! target grid
!        do icol = 1, ncol ! target grid
!          ir0 = (irow-1)*bs + 1; ir1 = ir0 + bs - 1
!          ic0 = (icol-1)*bs + 1; ic1 = ic0 + bs - 1
!          !
!          n = 0; r8v = r8ZERO; x(icol,irow) = mv
!          do jrow = ir0, ir1
!            do jcol = ic0, ic1
!              if (this%xr8(jcol,jrow) /= hdr%mvr8) then
!                if (hdr%i_uscl_type == i_uscl_geom) then
!                  r8v = r8v + log(this%xr8(icol,irow))
!                else
!                  r8v = r8v + this%xr8(icol,irow)
!                end if
!                n = n + 1
!              end if
!            end do
!          end do
!          select case(hdr%i_uscl_type)
!            case(i_uscl_arith, i_uscl_nodata)
!              if (n > 0) then
!                r8v = r8v / n
!              end if
!            case(i_uscl_geom)
!              if (n > 0) then
!                r8v = exp(r8v) / n
!              end if
!            case(i_uscl_sumcdr)
!              ! nothing
!          end select
!          x(icol,irow) = real(r8v, R8B)
!        end do
!      end do
!    end select
!    !
!    return
!  end subroutine hdr_get_grid_r8
  
  subroutine hdr_replace_grid(this, xi1, xi2, xi4, xi8, xr4, xr8, &
    mvi1, mvi2, mvi4, mvi8, mvr4, mvr8)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tHdr) :: this
    !
    integer(I1B), dimension(:,:), intent(in), optional :: xi1
    integer(I2B), dimension(:,:), intent(in), optional :: xi2
    integer(I4B), dimension(:,:), intent(in), optional :: xi4
    integer(I8B), dimension(:,:), intent(in), optional :: xi8
    real(R4B),    dimension(:,:), intent(in), optional :: xr4
    real(R8B),    dimension(:,:), intent(in), optional :: xr8
    !
    integer(I1B), intent(in), optional :: mvi1
    integer(I2B), intent(in), optional :: mvi2
    integer(I4B), intent(in), optional :: mvi4
    integer(I8B), intent(in), optional :: mvi8
    real(R4B),    intent(in), optional :: mvr4
    real(R8B),    intent(in), optional :: mvr8
    ! -- local
    type(thdrData), pointer :: dat => null()
! ------------------------------------------------------------------------------
    dat => this%dat
    !
    select case(this%i_data_type)
    case(i_i1)
      if ((.not.present(xi1)).or.(.not.present(mvi1))) then
        call errmsg('hdr_replace_grid: xi1 or mvi1 not found.')
      end if
      if (allocated(dat%xi1)) deallocate(dat%xi1)
      allocate(dat%xi1, source=xi1)
    case(i_i2)
      if ((.not.present(xi2)).or.(.not.present(mvi2))) then
        call errmsg('hdr_replace_grid: xi2 or mvi2 not found.')
      end if
      if (allocated(dat%xi2)) deallocate(dat%xi2)
      allocate(dat%xi2, source=xi2)
    case(i_i4)
      if ((.not.present(xi4)).or.(.not.present(mvi4))) then
        call errmsg('hdr_replace_grid: xi4 or mvi4 not found.')
      end if
      if (allocated(dat%xi4)) deallocate(dat%xi4)
      allocate(dat%xi4, source=xi4)
    case(i_i8)
      if ((.not.present(xi8)).or.(.not.present(mvi8))) then
        call errmsg('hdr_replace_grid: xi8 or mvi8 not found.')
      end if
      if (allocated(dat%xi8)) deallocate(dat%xi8)
      allocate(dat%xi8, source=xi8)
    case(i_r4)
      if ((.not.present(xr4)).or.(.not.present(mvr4))) then
        call errmsg('hdr_replace_grid: xr4 or mvr4 not found.')
      end if
      if (allocated(dat%xr4)) deallocate(dat%xr4)
      allocate(dat%xr4, source=xr4)
    case(i_r8)
      if ((.not.present(xr8)).or.(.not.present(mvr8))) then
        call errmsg('hdr_replace_grid: xr8 or mvr8 not found.')
      end if
      if (allocated(dat%xr8)) deallocate(dat%xr8)
      allocate(dat%xr8, source=xr8)
    end select
    !
    return
  end subroutine hdr_replace_grid
    
  subroutine hdr_set_grid_i4(this, x, mv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I4B), dimension(:,:), intent(in) :: x !change type!
    integer(I4B), intent(in) :: mv !change type!
    ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    integer(I4B) :: ncol, nrow, icol, irow
! ------------------------------------------------------------------------------
    !
    hdr => this%hdr; dat => this%dat
    ncol = hdr%ncol; nrow = hdr%nrow
    !
    select case(this%i_data_type)
    case(i_i1)
      hdr%mvi1 = int(mv,I1B); hdr%lmvi1 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi1(icol,irow) = int(x(icol,irow),I1B)
      end do; end do
    case(i_i2)
      hdr%mvi2 = int(mv,I2B); hdr%lmvi2 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi2(icol,irow) = int(x(icol,irow),I2B)
      end do; end do
    case(i_i4)
      hdr%mvi4 = int(mv,I4B); hdr%lmvi4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi4(icol,irow) = int(x(icol,irow),I4B)
      end do; end do
    case(i_i8)
      hdr%mvi8 = int(mv,I8B); hdr%lmvi8 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi8(icol,irow) = int(x(icol,irow),I8B)
      end do; end do
    case(i_r4)
      hdr%mvr4 = int(mv,R4B); hdr%lmvr4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xr4(icol,irow) = real(x(icol,irow),R4B)
      end do; end do
    case(i_r8)
      hdr%mvr8 = int(mv,R8B); hdr%lmvr4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xr8(icol,irow) = real(x(icol,irow),R8B)
      end do; end do
    end select
    !
    return
  end subroutine hdr_set_grid_i4
 
  subroutine hdr_set_grid_r4(this, x, mv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    real(R4B), dimension(:,:), intent(in) :: x !change type!
    real(R4B), intent(in) :: mv !change type!
    ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    integer(I4B) :: ncol, nrow, icol, irow
! ------------------------------------------------------------------------------
    !
    hdr => this%hdr; dat => this%dat
    ncol = hdr%ncol; nrow = hdr%nrow
    !
    select case(this%i_data_type)
    case(i_i1)
      hdr%mvi1 = int(mv,I1B); hdr%lmvi1 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi1(icol,irow) = int(x(icol,irow),I1B)
      end do; end do
    case(i_i2)
      hdr%mvi2 = int(mv,I2B); hdr%lmvi2 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi2(icol,irow) = int(x(icol,irow),I2B)
      end do; end do
    case(i_i4)
      hdr%mvi4 = int(mv,I4B); hdr%lmvi4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi4(icol,irow) = int(x(icol,irow),I4B)
      end do; end do
    case(i_i8)
      hdr%mvi8 = int(mv,I8B); hdr%lmvi8 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi8(icol,irow) = int(x(icol,irow),I8B)
      end do; end do
    case(i_r4)
      hdr%mvr4 = int(mv,R4B); hdr%lmvr4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xr4(icol,irow) = real(x(icol,irow),R4B)
      end do; end do
    case(i_r8)
      hdr%mvr8 = int(mv,R8B); hdr%lmvr4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xr8(icol,irow) = real(x(icol,irow),R8B)
      end do; end do
    end select
    !
    return
  end subroutine hdr_set_grid_r4
  
  subroutine hdr_set_grid_i8(this, x, mv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I8B), dimension(:,:), intent(in) :: x !change type!
    integer(I8B), intent(in) :: mv !change type!
    ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    integer(I4B) :: ncol, nrow, icol, irow
! ------------------------------------------------------------------------------
    !
    hdr => this%hdr; dat => this%dat
    ncol = hdr%ncol; nrow = hdr%nrow
    !
    select case(this%i_data_type)
    case(i_i1)
      hdr%mvi1 = int(mv,I1B); hdr%lmvi1 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi1(icol,irow) = int(x(icol,irow),I1B)
      end do; end do
    case(i_i2)
      hdr%mvi2 = int(mv,I2B); hdr%lmvi2 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi2(icol,irow) = int(x(icol,irow),I2B)
      end do; end do
    case(i_i4)
      hdr%mvi4 = int(mv,I4B); hdr%lmvi4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi4(icol,irow) = int(x(icol,irow),I4B)
      end do; end do
    case(i_i8)
      hdr%mvi8 = int(mv,I8B); hdr%lmvi8 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xi8(icol,irow) = int(x(icol,irow),I8B)
      end do; end do
    case(i_r4)
      hdr%mvr4 = int(mv,R4B); hdr%lmvr4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xr4(icol,irow) = real(x(icol,irow),R4B)
      end do; end do
    case(i_r8)
      hdr%mvr8 = int(mv,R8B); hdr%lmvr4 = .true.
      do irow = 1, nrow; do icol = 1, ncol
        dat%xr8(icol,irow) = real(x(icol,irow),R8B)
      end do; end do
    end select
    !
    return
  end subroutine hdr_set_grid_i8
  
  subroutine hdr_uscl_read_check_r4(this, xll_t, yll_t, cs_t, nc, nr, &
    ic0, ic1, ir0, ir1, bs)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    real(R4B), intent(in) :: xll_t
    real(R4B), intent(in) :: yll_t
    real(R4B), intent(in) :: cs_t
    integer(I4B), intent(in) :: nc
    integer(I4B), intent(in) :: nr
    integer(I4B), intent(out) :: ic0
    integer(I4B), intent(out) :: ic1
    integer(I4B), intent(out) :: ir0
    integer(I4B), intent(out) :: ir1
    integer(I4B), intent(out) :: bs
    ! -- local
    type(thdrHdr), pointer :: hdr
    logical :: lxll, lyll, lcs 
    real(R4B), pointer :: xll_s, yll_s, cs_s
    integer(I4B) :: n, mc, mr
! ------------------------------------------------------------------------------
    !
    hdr => this%hdr
    xll_s => hdr%xllr4; yll_s => hdr%yllr4; cs_s => hdr%csr4
    !
    if (cs_t > cs_s) then
      lcs  = (mod(cs_t,cs_s).eq.R4ZERO)
    else
      lcs  = (mod(cs_s,cs_t).eq.R4ZERO)
    end if
    lxll = (mod(xll_s-xll_t,cs_s).eq.R4ZERO)
    lyll = (mod(yll_s-yll_t,cs_s).eq.R4ZERO)
    !
    if (.not.lcs)  call errmsg('Not coinciding cell sizes.')
    if (.not.lxll) call errmsg('Not coinciding lower-left x-coordinates.')
    if (.not.lyll) call errmsg('Not coinciding lower-left y-coordinates.')
    !
    bs = int(cs_t/cs_s)
    !
    ic0 = int((xll_t-xll_s)/cs_s)+1
    ic1 = ic0 + bs*nc - 1
    ir1 = hdr%nrow - int((yll_t-yll_s)/cs_s)
    ir0 = ir1 - bs*nr + 1
    !
    ! check
    mc = ic1 - ic0 + 1; mr = ir1 - ir0 + 1
    if (mod(mc,bs).ne.0) call errmsg('Program error')
    if (mod(mr,bs).ne.0) call errmsg('Program error')
    !
    return
  end subroutine hdr_uscl_read_check_r4
  
  subroutine hdr_uscl_read_check_r8(this, xll_t, yll_t, cs_t, nc, nr, &
    ic0, ic1, ir0, ir1, bs)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    real(R8B), intent(in) :: xll_t
    real(R8B), intent(in) :: yll_t
    real(R8B), intent(in) :: cs_t
    integer(I4B), intent(in) :: nc
    integer(I4B), intent(in) :: nr
    integer(I4B), intent(out) :: ic0
    integer(I4B), intent(out) :: ic1
    integer(I4B), intent(out) :: ir0
    integer(I4B), intent(out) :: ir1
    integer(I4B), intent(out) :: bs
    ! -- local
    type(thdrHdr), pointer :: hdr
    logical :: lxll, lyll, lcs 
    real(R8B), pointer :: xll_s, yll_s, cs_s
    integer(I4B) :: n, mc, mr
! ------------------------------------------------------------------------------
    !
    hdr => this%hdr
    xll_s => hdr%xllr8; yll_s => hdr%yllr8; cs_s => hdr%csr8
    !
    if (cs_t > cs_s) then
      lcs  = (mod(cs_t,cs_s).eq.R8ZERO)
    else
      call errmsg('Error: downscaling not yet implemented')
    end if
    lxll = (mod(xll_s-xll_t,cs_s).eq.R8ZERO)
    lyll = (mod(yll_s-yll_t,cs_s).eq.R8ZERO)
    !
    if (.not.lcs)  call errmsg('Not coinciding cell sizes.')
    if (.not.lxll) call errmsg('Not coinciding lower-left x-coordinates.')
    if (.not.lyll) call errmsg('Not coinciding lower-left y-coordinates.')
    !
    bs = int(cs_t/cs_s)
    !
    ic0 = int((xll_t-xll_s)/cs_s)+1
    ic1 = ic0 + bs*nc - 1
    ir1 = hdr%nrow - int((yll_t-yll_s)/cs_s)
    ir0 = ir1 - bs*nr + 1
    !
    ! check
    mc = ic1 - ic0 + 1; mr = ir1 - ir0 + 1
    if (mod(mc,bs).ne.0) call errmsg('Program error')
    if (mod(mr,bs).ne.0) call errmsg('Program error')
    !
    return
  end subroutine hdr_uscl_read_check_r8
  
  subroutine hdr_get_val_init(this, fp)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    character(len=*), intent(inout) :: fp
    ! -- local
    character(len=MXSLEN) :: f
! ------------------------------------------------------------------------------
    call this%init()
    this%fp = fp
    !
    allocate(this%hdr)
    call this%hdr%init()
    call this%hdr%read(fp)
    call this%set_i_data_type()
    !
    f = trim(fp)//'.flt'
    call open_file(f, this%iu_bin, 'r', .true.)
    !
    return
  end subroutine hdr_get_val_init
  
  subroutine hdr_get_val_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    ! -- local
    logical :: lop
! ------------------------------------------------------------------------------
    call this%hdr%clean()
    deallocate(this%hdr)
    !
    inquire(this%iu_bin, opened=lop)
    if (lop) then
    
      close(this%iu_bin)
    end if
    call this%clean()
    !
    return
  end subroutine hdr_get_val_clean
  
 function hdr_get_val_r8(this, icol, irow, cs_t, lmv) result(x)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I4B), intent(in) :: icol, irow ! coarse index
    real(R8B), intent(in) :: cs_t
    logical, intent(out) :: lmv
    real(R8B) :: x
    ! -- local
    logical :: lcs
    integer(I4B) :: ir0, ir1, ic0, ic1, bs, n, ir, ic, nr, nc
    integer(I1B) :: i1v
    integer(I2B) :: i2v
    integer(I4B) :: i4v
    real(R4B) :: r4v
    real(R4B), dimension(:), allocatable :: r4a
    real(R8B) :: r8v, xexp
! ------------------------------------------------------------------------------
    !
    bs = 1; ir0 = irow; ir1 = irow; ic0 = icol; ic1 = icol
    if (cs_t >= this%hdr%csr8) then
      lcs = (mod(cs_t,this%hdr%csr8).eq.R8ZERO)
      if (.not.lcs) call errmsg('Not coinciding cell sizes.')
      bs = int(cs_t/this%hdr%csr8)
      ir0 = (irow-1)*bs + 1; ir1 = ir0 + bs - 1
      ic0 = (icol-1)*bs + 1; ic1 = ic0 + bs - 1
      nr = ir1 - ir0 + 1; nc = ic1 - ic0 + 1
    else
      call errmsg('Error: downscaling not yet supported')
    end if
    !
    x = R8ZERO; xexp = R8ZERO; n = 0
    !
    select case(this%i_data_type)
    case(i_i1)
      do ir = ir0, ir1
        do ic = ic0, ic1
          call this%read_val(icol=ic, irow=ir, i1v=i1v)
          if (i1v /= this%hdr%mvi1) then
            n = n + 1
            x = x + real(i1v,R8B); xexp = xexp + log(real(i1v,R8B))
          end if
        end do
      end do
    case(i_i2)
      do ir = ir0, ir1
        do ic = ic0, ic1
          call this%read_val(icol=ic, irow=ir, i2v=i2v)
          if (i2v /= this%hdr%mvi2) then
            n = n + 1
            x = x + real(i2v,R8B); xexp = xexp + log(real(i2v,R8B))
          end if
        end do
      end do
    case(i_i4)
      do ir = ir0, ir1
        do ic = ic0, ic1
          call this%read_val(icol=ic, irow=ir, i4v=i4v)
          if (i4v /= this%hdr%mvi4) then
            n = n + 1
            x = x + real(i4v,R8B); xexp = xexp + log(real(i4v,R8B))
          end if
        end do
      end do
    case(i_r4)
      do ir = ir0, ir1
        do ic = ic0, ic1
          call this%read_val(icol=ic, irow=ir, r4v=r4v)
          if (r4v /= this%hdr%mvr4) then
            n = n + 1
            x = x + real(r4v,R8B); xexp = xexp + log(real(r4v,R8B))
          end if
        end do
      end do
    end select
    !
    if (n == 0) then
      lmv = .true.
    else
      lmv = .false.
    end if
    !
    select case(this%hdr%i_uscl_type)
      case(i_uscl_arith, i_uscl_nodata)
        if (n > 0) then
          x = x / n
         end if
       case(i_uscl_geom)
         if (n > 0) then
           x = exp(xexp/n)
         end if
       case(i_uscl_sumcdr)
         ! nothing
    end select
    return
  end function hdr_get_val_r8
  !
  subroutine hdr_read_val_xy(this, xd, yd, csd, &
    i_uscl, i_dscl, i1v, i2v, i4v, i8v, r4v, r8v)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    !
    real(R8B), intent(in) :: xd
    real(R8B), intent(in) :: yd
    real(R8B), intent(in) :: csd
    integer(I4B), intent(in), optional :: i_uscl
    integer(I4B), intent(in), optional :: i_dscl
    !
    integer(I1B), intent(out), optional :: i1v
    integer(I2B), intent(out), optional :: i2v
    integer(I4B), intent(out), optional :: i4v
    integer(I8B), intent(out), optional :: i8v
    real(R4B),    intent(out), optional :: r4v
    real(R8B),    intent(out), optional :: r8v
    !
    ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    logical :: l_in
    real(R8B), dimension(2,1) :: x_src, x_dst
    real(R8B), dimension(2,1) :: xll
    real(R8B), dimension(2,1) :: xur
    real(R8B), dimension(4)   :: va
    real(R8B) :: xs, ys, v, mv, src_cs
    integer(I4B), dimension(:,:), allocatable ::sten
    integer(I4B) :: bs, ic, ir
    integer(I4B) :: ic1, ic2, ic3, ic4, ir1, ir2, ir3, ir4
 ! ------------------------------------------------------------------------------
    !
    hdr => this%hdr
    if (csd > hdr%csr8) then ! upscale
      if (mod(csd, hdr%csr8) /= 2) then
        call errmsg('hdr_read_val_xy: inconsistent cell sizes.')
      end if
      bs =  csd/hdr%csr8
    else
      call get_icir(ic, ir, xd, yd, hdr%xllr8, hdr%xurr8, &
        hdr%yllr8, hdr%yurr8, hdr%csr8, &
        hdr%ncol, hdr%nrow, l_in)
      
      if (i_dscl == i_dscl_nointp) then
        if (present(r4v)) then
          r4v = hdr%mvr4
          if (l_in) call this%read_val(icol=ic, irow=ir, r4v=r4v, &
            try_reuse=.true.)
        end if
        if (present(r8v)) then
          r8v = hdr%mvr8
          if (l_in) call this%read_val(icol=ic, irow=ir, r8v=r8v, &
            try_reuse=.true.)
        end if
      else
        ! linear interpolation
        src_cs = hdr%csr8
        call get_xy(xs, ys, ic, ir, hdr%xllr8, hdr%yurr8, hdr%csr8)
        x_dst = reshape([xd, yd], shape(x_dst))
        !
        ! check
        if ((xd == xs).and.(yd == ys)) then
          call errmsg('hdr_down_scale_intp: program error.')
        end if
        !
        mv = hdr%mvr8
        va = mv
        call get_stencil(ic, ir, hdr%ncol, hdr%nrow, 5, sten)
        !
        if (xd < xs) then ! west
          if (yd > ys) then ! north
            xll = reshape([xs - src_cs, ys], shape(xll))
            xur = reshape([xs, ys + src_cs], shape(xur))
            ic1 = sten(1,i_w);  ir1 = sten(2,i_w)
            ic2 = sten(1,i_nw); ir2 = sten(2,i_nw)
            ic3 = sten(1,i_n);  ir3 = sten(2,i_n)
            ic4 = sten(1,i_p);  ir4 = sten(2,i_p)
          else ! south
            xll = reshape([xs - src_cs, ys - src_cs], shape(xll))
            xur = reshape([xs, ys], shape(xur))
            ic1 = sten(1,i_sw); ir1 = sten(2,i_sw)
            ic2 = sten(1,i_w);  ir2 = sten(2,i_w)
            ic3 = sten(1,i_p);  ir3 = sten(2,i_p)
            ic4 = sten(1,i_s);  ir4 = sten(2,i_s)
          end if
        else ! east
          if (yd > ys) then ! north
            xll = reshape([xs, ys], shape(xll))
            xur = reshape([xs + src_cs, ys + src_cs], shape(xur))
            ic1 = sten(1,i_p);  ir1 = sten(2,i_p)
            ic2 = sten(1,i_n);  ir2 = sten(2,i_n)
            ic3 = sten(1,i_ne); ir3 = sten(2,i_ne)
            ic4 = sten(1,i_e);  ir4 = sten(2,i_e)
          else ! south
            xll = reshape([xs, ys - src_cs], shape(xll))
            xur = reshape([xs + src_cs, ys], shape(xur))
            ic1 = sten(1,i_s);  ir1 = sten(2,i_s)
            ic2 = sten(1,i_p);  ir2 = sten(2,i_p)
            ic3 = sten(1,i_ne); ir3 = sten(2,i_ne)
            ic4 = sten(1,i_se); ir4 = sten(2,i_se)
          end if
        end if
        !
        if (present(r4v)) then
          call this%read_val(icol=ic1, irow=ir1, r4v=r4v, try_reuse=.true.); va(1) = real(r4v,R8B)
          call this%read_val(icol=ic2, irow=ir2, r4v=r4v, try_reuse=.true.); va(2) = real(r4v,R8B)
          call this%read_val(icol=ic3, irow=ir3, r4v=r4v, try_reuse=.true.); va(3) = real(r4v,R8B)
          call this%read_val(icol=ic4, irow=ir4, r4v=r4v, try_reuse=.true.); va(4) = real(r4v,R8B)
        end if
        !
        v = bilinear_interpolation(x_dst, xll, xur, va, mv)
        !
        if (present(r4v)) then
          r4v = real(v,R4B)
        end if
        if (present(r8v)) then
          r8v = v
        end if
      end if
    end if
    !
    return
  end subroutine hdr_read_val_xy
 
  subroutine hdr_read_val(this, icol, irow, &
    i1v, i2v, i4v, i8v, r4v, r8v, try_reuse)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    !
    integer(I4B), intent(in) :: icol
    integer(I4B), intent(in) :: irow
    !
    integer(I1B), intent(out), optional :: i1v
    integer(I2B), intent(out), optional :: i2v
    integer(I4B), intent(out), optional :: i4v
    integer(I8B), intent(out), optional :: i8v
    real(R4B),    intent(out), optional :: r4v
    real(R8B),    intent(out), optional :: r8v
    !
    logical, intent(in), optional :: try_reuse
    ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
    logical :: tr, dar
    integer(I8B) :: p
! ------------------------------------------------------------------------------
    !
    if ((icol < 1).or.(icol > this%hdr%ncol).or. &
        (irow < 1).or.(irow > this%hdr%nrow)) then
      hdr => this%hdr
      if (present(i1v)) i1v = hdr%mvi1
      if (present(i2v)) i2v = hdr%mvi2
      if (present(i4v)) i4v = hdr%mvi4
      if (present(i8v)) i8v = hdr%mvi8
      if (present(r4v)) r4v = hdr%mvr4
      if (present(r8v)) r8v = hdr%mvr8
      return
    end if
    !
    tr = .false.
    if (present(try_reuse)) then
      tr = try_reuse
    end if
    !
    dat => this%dat
    if (present(i1v)) then
      dar = .true.
      if (tr) then
        if (allocated(dat%xi1)) then
          i1v = dat%xi1(icol,irow); dar = .false.
        end if
      end if
      if (dar) then
        p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*I1B
        read(unit=this%iu_bin,pos=p) i1v
      end if
    end if
    !
    if (present(i2v)) then
      dar = .true.
      if (tr) then
        if (allocated(dat%xi2)) then
          i2v = dat%xi2(icol,irow); dar = .false.
        end if
      end if
      if (dar) then
        p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*I2B
        read(unit=this%iu_bin,pos=p) i2v
      end if
    end if
    !
    if (present(i4v)) then
      dar = .true.
      if (tr) then
        if (allocated(dat%xi4)) then
          i4v = dat%xi4(icol,irow); dar = .false.
        end if
      end if
      if (dar) then
        p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*I4B
        read(unit=this%iu_bin,pos=p) i4v
      end if
    end if
    !
    if (present(i8v)) then
      dar = .true.
      if (tr) then
        if (allocated(dat%xi8)) then
          i8v = dat%xi8(icol,irow); dar = .false.
        end if
      end if
      if (dar) then
        p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*I8B
       read(unit=this%iu_bin,pos=p) i8v
      end if
    end if
    !
    if (present(r4v)) then
      dar = .true.
      if (tr) then
        if (allocated(dat%xr4)) then
          r4v = dat%xr4(icol,irow); dar = .false.
        end if
      end if
      if (dar) then
        p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*R4B
        read(unit=this%iu_bin,pos=p) r4v
      end if
    end if
    !
    if (present(r8v)) then
      dar = .true.
      if (tr) then
        if (allocated(dat%xr8)) then
          r8v = dat%xr8(icol,irow); dar = .false.
        end if
      end if
      if (dar) then
        p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*R8B
        read(unit=this%iu_bin,pos=p) r8v
      end if
    end if
    !
    return
  end subroutine hdr_read_val
  
  subroutine hdr_read_arr(this, icol, irow, i1a, i2a, i4a, i8a, r4a, r8a)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    integer(I4B), intent(in) :: icol
    integer(I4B), intent(in) :: irow
    !
    integer(I1B), dimension(:), intent(out), optional :: i1a
    integer(I2B), dimension(:), intent(out), optional :: i2a
    integer(I4B), dimension(:), intent(out), optional :: i4a
    integer(I8B), dimension(:), intent(out), optional :: i8a
    real(R4B),    dimension(:), intent(out), optional :: r4a
    real(R8B),    dimension(:), intent(out), optional :: r8a
    ! -- local
    type(tHdrHdr), pointer :: hdr => null()
    logical :: has_nan
    integer(I4B) :: i, n
    integer(I8B) :: p
! ------------------------------------------------------------------------------
    !
    hdr => this%hdr
    !
    if (present(i1a)) then
      p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*I1B
      read(unit=this%iu_bin,pos=p) i1a
    end if
    !
    if (present(i2a)) then
      p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*I2B
      read(unit=this%iu_bin,pos=p) i2a
    end if
    !
    if (present(i4a)) then
      p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*I4B
      read(unit=this%iu_bin,pos=p) i4a
    end if
    !
    if (present(i8a)) then
      p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*I8B
      read(unit=this%iu_bin,pos=p) i8a
    end if
    !
    if (present(r4a)) then
      p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*R4B
      read(unit=this%iu_bin,pos=p) r4a
      !
      ! check for nans
      has_nan = .false.; n = 0
      do i = 1, size(r4a)
        if (check_nan(r4a(i))) then
          has_nan = .true.
          n = n + 1
        end if
      end do
      if (has_nan) then
        if (.not.hdr%lmvr4) then
          call errmsg('hdr_read_arr: cannot set missing value for r4a')
        end if
        do i = 1, size(r4a)
          if (check_nan(r4a(i))) r4a(i) = hdr%mvr4
        end do
      end if
    end if
    !
    if (present(r8a)) then
      p = 1+((int(irow,I8B)-1)*int(this%hdr%ncol,I8B)+int(icol,I8B)-1)*R8B
      read(unit=this%iu_bin,pos=p) r8a
    end if
    !
    return
  end subroutine hdr_read_arr
  
  subroutine hdr_write(this, fp, file_type)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(thdr) :: this
    character(len=*), intent(in) :: fp
    character(len=*), intent(in), optional :: file_type
    ! -- local
    integer(I4B) :: i_file_type
    type(tHdrHdr), pointer :: hdr => null()
    type(tHdrData), pointer :: dat => null()
! ------------------------------------------------------------------------------
    !
    if (present(file_type)) then
      if ((file_type == 'envi').or.(file_type == 'ENVI')) then
        i_file_type = i_env
      end if
    else
      i_file_type = i_flt
    end if
    !
    hdr => this%hdr; dat => this%dat
    !
    if (i_file_type == i_flt) then
      select case(this%i_data_type)
      case(i_i1)
        if (.not.allocated(dat%xi4)) then
          call errmsg('Could not write '//trim(fp))
        end if
        call writeflt_i1_r8(fp, dat%xi1, hdr%ncol, hdr%nrow, &
          hdr%xllr8, hdr%yllr8, hdr%csr8, hdr%mvi1)
      case(i_i2)
        if (.not.allocated(dat%xi4)) then
          call errmsg('Could not write '//trim(fp))
        end if
        call writeflt_i2_r8(fp, dat%xi2, hdr%ncol, hdr%nrow, &
          hdr%xllr8, hdr%yllr8, hdr%csr8, hdr%mvi2)
      case(i_i4)
        if (.not.allocated(dat%xi4)) then
          call errmsg('Could not write '//trim(fp))
        end if
        call writeflt_i4_r8(fp, dat%xi4, hdr%ncol, hdr%nrow, &
          hdr%xllr8, hdr%yllr8, hdr%csr8, hdr%mvi4)
      case(i_r4)
        if (.not.allocated(dat%xr4)) then
          call errmsg('Could not write '//trim(fp))
        end if
        call writeflt_r4_r8(fp, dat%xr4, hdr%ncol, hdr%nrow, &
          hdr%xllr8, hdr%yllr8, hdr%csr8, hdr%mvr4)
      end select
    else
      select case(this%i_data_type)
      case(i_i1)
        call errmsg('Integer 1 not supported: could not write '//trim(fp))
      case(i_i2)
        call write_env(fp, hdr%xllr8, hdr%yllr8, hdr%csr8, xi2=dat%xi2)
      case(i_i4)
        call write_env(fp, hdr%xllr8, hdr%yllr8, hdr%csr8, xi4=dat%xi4)
      case(i_i8)
        call write_env(fp, hdr%xllr8, hdr%yllr8, hdr%csr8, xi8=dat%xi8)
      case(i_r4)
        call write_env(fp, hdr%xllr8, hdr%yllr8, hdr%csr8, xr4=dat%xr4)
      case(i_r8)
        call write_env(fp, hdr%xllr8, hdr%yllr8, hdr%csr8, xr8=dat%xr8)
      end select
      !
    end if
    !
    return
  end subroutine hdr_write
  
! ==============================================================================
! ==============================================================================
! general subroutine
! ==============================================================================
! ==============================================================================

  subroutine get_icir(ic, ir, x, y, xll, xur, yll, yur, cs, ncol, nrow, l_in)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), intent(out) :: ic
    integer(I4B), intent(out) :: ir
    real(R8B), intent(in) :: x
    real(R8B), intent(in) :: y
    real(R8B), intent(in) :: xll
    real(R8B), intent(in) :: xur
    real(R8B), intent(in) :: yll
    real(R8B), intent(in) :: yur
    real(R8B), intent(in) :: cs
    integer(I4B), intent(in) :: ncol
    integer(I4B), intent(in) :: nrow
    logical, intent(inout) :: l_in
    ! -- local
! ------------------------------------------------------------------------------
    !
    l_in = .true.
    !
    ! check 1: xy
    if (x < xll) then
      l_in = .false.; ic = 0; ir = 0; return
    end if
    if (x > xur) then
      l_in = .false.; ic = 0; ir = 0; return
    end if
    if (y < yll) then
      l_in = .false.; ic = 0; ir = 0; return
    end if
    if (y > yur) then
      l_in = .false.; ic = 0; ir = 0; return
    end if
    !
    ic = int((x - xll)/cs) + 1
    ir = int((yur - y)/cs) + 1
    !
    if (ic < 1) then
      l_in = .false.; ic = 0; ir = 0; return
    end if
    if (ic > ncol) then
      l_in = .false.; ic = 0; ir = 0; return
    end if
    if (ir < 1) then
      l_in = .false.; ic = 0; ir = 0; return
    end if
    if (ir > nrow) then
      l_in = .false.; ic = 0; ir = 0; return
    end if
    !
    return
 end subroutine get_icir
  
 subroutine writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, nodata, &
  nbits, pixeltype, hdrKeys)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: iu, ncol, nrow
    character(len=*), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    integer(I4B), intent(in) :: nbits
    character(len=*), intent(in) :: pixeltype
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    write(iu,'(a)') 'ncols '//ta((/ncol/))
    write(iu,'(a)') 'nrows '//ta((/nrow/))
    write(iu,'(a)') 'xllcorner '//ta((/xll/))
    write(iu,'(a)') 'yllcorner '//ta((/yll/))
    write(iu,'(a)') 'cellsize '//ta((/cs/))
    write(iu,'(a)') 'nodata_value '//trim(nodata)
    write(iu,'(a)') 'nbits '//ta((/nbits/))
    write(iu,'(a)') 'pixeltype '//trim(pixeltype)
    write(iu,'(a)') 'byteorder lsbfirst'
    if (present(hdrkeys)) then
      do i = 1, size(hdrkeys)
        write(iu,'(a)') trim(hdrkeys(i))
      end do
    end if
    close(iu)
    !
    return
  end subroutine writeflt_header_r4
  
  subroutine writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, nodata, &
  nbits, pixeltype, hdrKeys)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: iu, ncol, nrow
    character(len=*), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    integer(I4B), intent(in) :: nbits
    character(len=*), intent(in) :: pixeltype
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    write(iu,'(a)') 'ncols '//ta((/ncol/))
    write(iu,'(a)') 'nrows '//ta((/nrow/))
    write(iu,'(a)') 'xllcorner '//ta((/xll/))
    write(iu,'(a)') 'yllcorner '//ta((/yll/))
    write(iu,'(a)') 'cellsize '//ta((/cs/))
    write(iu,'(a)') 'nodata_value '//trim(nodata)
    write(iu,'(a)') 'nbits '//ta((/nbits/))
    write(iu,'(a)') 'pixeltype '//trim(pixeltype)
    write(iu,'(a)') 'byteorder lsbfirst'
    if (present(hdrkeys)) then
      do i = 1, size(hdrkeys)
        write(iu,'(a)') trim(hdrkeys(i))
      end do
    end if
    close(iu)
    !
    return
  end subroutine writeflt_header_r8
!
  subroutine writeflt_i1_r4(fp, x, ncol, nrow, xll, yll, cs, nodata, hdrKeys)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/hdr.html#raster-hdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I1B), dimension(ncol,nrow), intent(in) :: x
    integer(I1B), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    if (present(hdrKeys)) then
      call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        8, 'signedint', hdrKeys)
    else
      call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        8, 'signedint')
    end if
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    !
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i1_r4
  
  subroutine writeflt_i2_r4(fp, x, ncol, nrow, xll, yll, cs, nodata, hdrKeys)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/hdr.html#raster-hdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I2B), dimension(ncol,nrow), intent(in) :: x
    integer(I2B), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    if (present(hdrKeys)) then
      call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        8, 'signedint', hdrKeys)
    else
      call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        8, 'signedint')
    end if
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    !
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i2_r4
    
  subroutine writeflt_i4_r4(fp, x, ncol, nrow, xll, yll, cs, nodata, hdrKeys)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/hdr.html#raster-hdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I4B), dimension(ncol,nrow), intent(in) :: x
    integer(I4B), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    if (present(hdrKeys)) then
      call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        32, 'signedint', hdrKeys)
    else
      call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        32, 'signedint')
    end if
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    !
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i4_r4
 
  subroutine writeflt_r4_r4(fp, x, ncol, nrow, xll, yll, cs, nodata, hdrKeys)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/hdr.html#raster-hdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    real(R4B), dimension(ncol,nrow), intent(in) :: x
    real(R4B), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    if (present(hdrKeys)) then
      call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        32, 'float', hdrKeys)
    else
      call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        32, 'float')
    end if
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    if (.false.)then
      f = trim(fp)//'.flt.asc'
      call open_file(f, iu, 'w')
      write(iu,*)((x(icol,irow),icol=1,ncol),irow=1,nrow)
      close(iu)
    end if
    return
  end subroutine writeflt_r4_r4
  
  subroutine writeflt_i1_r8(fp, x, ncol, nrow, xll, yll, cs, nodata, hdrKeys)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/hdr.html#raster-hdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I1B), dimension(ncol,nrow), intent(in) :: x
    integer(I1B), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    if (present(hdrKeys)) then
      call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        8, 'signedint', hdrKeys)
    else
      call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        8, 'signedint')
    end if
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i1_r8
  
  subroutine writeflt_i2_r8(fp, x, ncol, nrow, xll, yll, cs, nodata, hdrKeys)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/hdr.html#raster-hdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I2B), dimension(ncol,nrow), intent(in) :: x
    integer(I2B), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    if (present(hdrKeys)) then
      call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        8, 'signedint', hdrKeys)
    else
      call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        8, 'signedint')
    end if
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i2_r8
  
  subroutine writeflt_i4_r8(fp, x, ncol, nrow, xll, yll, cs, nodata, hdrKeys)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/hdr.html#raster-hdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I4B), dimension(ncol,nrow), intent(in) :: x
    integer(I4B), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    if (present(hdrKeys)) then
      call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        32, 'signedint', hdrKeys)
    else
      call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
        32, 'signedint')
    end if
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i4_r8
 
  subroutine writeflt_r4_r8(fp, x, ncol, nrow, xll, yll, cs, nodata, hdrKeys)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/hdr.html#raster-hdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    real(R4B), dimension(ncol,nrow), intent(in) :: x
    real(R4B), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    character(len=*), dimension(:), intent(in), optional :: hdrKeys
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    if (present(hdrKeys)) then
      call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
      32, 'float', hdrKeys)
    else
      call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), &
      32, 'float')
    end if
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_r4_r8
!
  subroutine write_env(fp,  xllr8, yllr8, csr8, xi2, xi4, xi8, xr4, xr8)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: fp
    real(R8B), intent(in) :: xllr8, yllr8, csr8
    !
    integer(I2B), dimension(:,:), intent(in), optional :: xi2
    integer(I4B), dimension(:,:), intent(in), optional :: xi4
    integer(I8B), dimension(:,:), intent(in), optional :: xi8 ! NOT SUPPORTED IN QGIS
    real(R4B),    dimension(:,:), intent(in), optional :: xr4
    real(R8B),    dimension(:,:), intent(in), optional :: xr8
    !
    ! -- locals
    logical :: ldone
    character(len=MXSLEN) :: f, s
    integer(I4B) :: iu, ic, ir, nc, nr, env_data_type
    real(R8B) :: yurr8
! ------------------------------------------------------------------------------
    !
    ldone = .true.
    if (present(xi2)) then
      nc = size(xi2,1); nr = size(xi2,2); env_data_type = 2; ldone = .false.
    end if
    if (present(xi4)) then
      nc = size(xi4,1); nr = size(xi4,2); env_data_type = 3; ldone = .false.
    end if
    if (present(xi8)) then
      nc = size(xi8,1); nr = size(xi8,2); env_data_type = 14; ldone = .false.
    end if
    if (present(xr4)) then
      nc = size(xr4,1); nr = size(xr4,2); env_data_type = 4; ldone = .false.
    end if
    if (present(xr8)) then
      nc = size(xr8,1); nr = size(xr8,2); env_data_type = 5; ldone = .false.
    end if
    if (ldone) then
      call logmsg('Nothing to do...'); return
    end if
    !
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    write(iu,'(a)') 'ENVI'
    write(iu,'(a)') 'samples = '//ta((/nc/))
    write(iu,'(a)') 'lines   = '//ta((/nr/))
    write(iu,'(a)') 'bands   = 1'
    write(iu,'(a)') 'header offset = 0'
    write(iu,'(a)') 'file type = ENVI Standard'
    write(iu,'(a)') 'data type = '//ta((/env_data_type/))
    write(iu,'(a)') 'interleave = bsq'
    write(iu,'(a)') 'byte order = 0'
    yurr8 = yllr8 + nr*csr8
    s = '{, 1, 1, '//ta((/xllr8/))//', '//ta((/yurr8/))//', '//ta((/csr8/))//', '//ta((/csr8/))//',}'
    write(iu,'(a)') 'map info = '//trim(s)
    close(iu)
    !
    f = trim(fp)
    call open_file(f, iu, 'w', .true.)
    !
    if (present(xi2)) write(iu)((xi2(ic,ir),ic=1,nc),ir=1,nr)
    if (present(xi4)) write(iu)((xi4(ic,ir),ic=1,nc),ir=1,nr)
    if (present(xi8)) write(iu)((xi8(ic,ir),ic=1,nc),ir=1,nr)
    if (present(xr4)) write(iu)((xr4(ic,ir),ic=1,nc),ir=1,nr)
    if (present(xr8)) write(iu)((xr8(ic,ir),ic=1,nc),ir=1,nr)
    !
    close(iu)
    !
    return
  end subroutine write_env
  
end module hdrModule

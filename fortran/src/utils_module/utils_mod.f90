module utilsmod
  use, intrinsic :: iso_fortran_env , only: error_unit, output_unit, &
    I1B => int8, I2B => int16, I4B => int32, I8B => int64, &
    R4B => real32, R8B => real64

  implicit none

  integer(I4B), parameter :: MXSLEN = 1024
  integer(I4B), parameter :: MXSLENLONG = 10240
  
  integer(I4B), parameter :: i_c  = 0 ! character
  integer(I4B), parameter :: i_i1 = 1 ! integer 1 byte
  integer(I4B), parameter :: i_i2 = 2 ! integer 2 bytes
  integer(I4B), parameter :: i_i4 = 3 ! integer 4 bytes
  integer(I4B), parameter :: i_i8 = 4 ! integer 8 bytes
  integer(I4B), parameter :: i_r4 = 5 ! real 4 bytes (float/single)
  integer(I4B), parameter :: i_r8 = 6 ! real 8 bytes (double)
  integer(I4B), parameter :: i_l4 = 7 ! logical 4 bytes
  integer(I4B), parameter :: i_num = i_l4
  !
  character(len=MXSLEN), dimension(i_num), parameter :: num_names = &
    ['integer_1b', 'integer_2b', 'integer_4b', &
     'integer_8b', 'real_4b', 'real_8b', 'logical_4b']
  !
  integer(I4B), parameter :: i_p  = 1
  integer(I4B), parameter :: i_n  = 2
  integer(I4B), parameter :: i_s  = 3
  integer(I4B), parameter :: i_w  = 4
  integer(I4B), parameter :: i_e  = 5
  integer(I4B), parameter :: i_nw = 6 
  integer(I4B), parameter :: i_ne = 7
  integer(I4B), parameter :: i_sw = 8
  integer(I4B), parameter :: i_se = 9
  integer(I4B), parameter :: i_sten = i_se
  !
  integer(I4B), dimension(2,i_sten) :: sten_dir
  data sten_dir/0, 0, 0, -1, 0, 1, -1, 0, 1, 0, -1, -1, 1, -1, -1, 1, 1, 1/
  !
  type tNum
    logical, dimension(i_num) :: flg
    integer(I1B) :: i1v
    integer(I2B) :: i2v
    integer(I4B) :: i4v
    integer(I8B) :: i8v
    real(R4B)    :: r4v
    real(R8B)    :: r8v
    logical      :: l4v
  contains
    procedure :: init          => tNum_init
    procedure :: set_val_by_s  => tNum_set_val_by_s
    procedure :: set_val       => tNum_set_val
    procedure :: get_val       => tNum_get_val
    procedure :: copy          => tNum_copy
  end type tNum
  
  integer(I4B), dimension(:,:), allocatable :: i4wk2d

#ifdef LINUX
  integer(I4B), parameter :: os = 2
#else
  integer(I4B), parameter :: os = 1
#endif

  character(len=1), parameter :: win_slash = '\'
  character(len=1), parameter :: lin_slash = '/'

  character(len=MXSLEN), dimension(100) :: sa
  character(len=1), parameter :: comment = '#'

  integer(I4B), parameter :: IZERO = 0
  real(R4B), parameter    :: RZERO = 0.0
  real(R4B), parameter    :: RONE = 1.0
  real(R8B), parameter    :: DZERO = 0.d0
  real(R8B), parameter    :: DONE  = 1.d0
  real(R8B), parameter    :: DHALF  = 0.5d0
  !
  integer(I1B), parameter :: I1ZERO = 0
  integer(I2B), parameter :: I2ZERO = 0
  integer(I4B), parameter :: I4ZERO = 0
  integer(I8B), parameter :: I8ZERO = 0
  real(R4B),    parameter :: R4ZERO = 0.0
  real(R8B),    parameter :: R8ZERO = 0.d0
  !
  integer(I1B), parameter :: I1ONE = 1
  integer(I2B), parameter :: I2ONE = 1
  integer(I4B), parameter :: I4ONE = 1
  integer(I8B), parameter :: I8ONE = 1
  real(R4B),    parameter :: R4ONE = 1.0
  real(R8B),    parameter :: R8ONE = 1.d0
  !
  integer(I1B), parameter :: I1MINONE = -1
  integer(I2B), parameter :: I2MINONE = -1
  integer(I4B), parameter :: I4MINONE = -1
  integer(I8B), parameter :: I8MINONE = -1
  real(R4B),    parameter :: R4MINONE = -1.0
  real(R8B),    parameter :: R8MINONE = -1.d0
  !
  real(R4B),    parameter :: R4HALF = 0.5
  real(R8B),    parameter :: R8HALF = 0.5d0
  !
  real(R4B),    parameter :: R4TINY = TINY(R4ZERO)
  real(R8B),    parameter :: R8TINY = TINY(R8ZERO)
  !
  type tVal
    logical :: lnum = .true.
    type(tNum), pointer           :: x => null()
    character(len=:), allocatable :: s
  contains
    procedure :: init        => tVal_init
    procedure :: clean       => tVal_clean
    procedure :: set_val     => tVal_set_val
    procedure :: set_val_s   => tVal_set_val_s
    procedure :: get_val     => tVal_get_val
    procedure :: get_val_arr => tVal_get_val_arr
    procedure :: get_val_s   => tVal_get_val_s
    procedure :: copy        => tVal_copy
  end type tVal
  !
  type tIniSect
    character(len=MXSLEN)                            :: name
    integer(I4B)                                     :: n = 0
    character(len=MXSLEN), dimension(:), allocatable :: keys
    type(tVal), dimension(:), pointer                :: v => null()
  contains
    procedure :: init  => tIniSect_init
    procedure :: clean => tIniSect_clean
    procedure :: set   => tIniSect_set
  end type tIniSect
  !
  type tIni
    character(len=MXSLEN) :: f = ''
    integer(I4B) :: n_sect
    character(len=MXSLEN), dimension(:), allocatable :: sect_names
    type(tIniSect), dimension(:), pointer :: sect => null()
  contains
    procedure :: init    => tIni_init
    procedure :: clean   => tIni_clean
    procedure :: read    => tIni_read
    procedure :: get_val => tIni_get_val
  end type tIni
  !
  type tCSV_hdr
    character(len=MXSLEN) :: key = ''
    integer(I4B)          :: i_type = 0
  contains
    procedure :: init => tCSV_hdr_init
    procedure :: set  => tCSV_hdr_set
    procedure :: get  => tCSV_hdr_get
    procedure :: copy => tCSV_hdr_copy
  end type tCSV_hdr
  type tCSV
    character(len=MXSLEN) :: file = ''
    integer(I4B) :: nc     = 0
    integer(I4B) :: nr     = 0
    integer(I4B) :: nc_max = 0
    integer(I4B) :: nr_max = 0
    type(tCSV_hdr), dimension(:),   pointer :: hdr => null()
    type(tVal),     dimension(:,:), pointer :: val => null()
  contains
    procedure :: init          => tCSV_init
    procedure :: read          => tCSV_read
    procedure :: write         => tCSV_write
    procedure :: get_nc        => tCSV_get_nc
    procedure :: get_nr        => tCSV_get_nr
    procedure :: get_col       => tCSV_get_col
    procedure :: get_row       => tCSV_get_row
    procedure :: exist_col     => tCSV_exist_col
    procedure :: add_key       => tCSV_add_key
    procedure :: get_key       => tCSV_get_key
    procedure :: read_hdr      => tCSV_read_hdr
    procedure :: set_hdr       => tCSV_set_hdr
    procedure :: add_hdr       => tCSV_add_hdr
    procedure :: set_val       => tCSV_set_val
    procedure :: get_val       => tCSV_get_val
    procedure :: get_column    => tCSV_get_column
    procedure :: get_matrix    => tCSV_get_matrix
    procedure :: get_selection => tCSV_get_selection
    procedure :: clean      => tCSV_clean
    procedure :: clean_and_init_row => tCSV_clean_and_init_row
  end type tCSV
  !
  interface get_unique
    module procedure :: get_unique_i4
  end interface get_unique
  private :: get_unique_i4
  
  interface fillgap
    module procedure :: fillgap_r4
  end interface fillgap
  private :: fillgap_r4

  interface fill_with_nearest
    module procedure :: fill_with_nearest_i4
    module procedure :: fill_with_nearest_r4
  end interface fill_with_nearest
  private :: fill_with_nearest_r4
  
  interface readidf_block
    module procedure :: readidf_block_i4
    module procedure :: readidf_block_r4
    module procedure :: readidf_block_r8
  end interface readidf_block
  private :: readidf_block_i4, readidf_block_r4, readidf_block_r8

  interface readidf
    module procedure :: readidf_i_r
    !module procedure :: readidf_r_r
    module procedure :: readidf_r4_r4
    module procedure :: readidf_r_d
  end interface
  private :: readidf_i_r, readidf_r_d, readidf_r_r, readidf_r4_r4

  interface writebin
    module procedure :: writebin_i
  end interface writebin
  private writebin_i

!  interface writeidf
!    module procedure :: writeidf_i_r
!    module procedure :: writeidf_r_r
!    module procedure :: writeidf_r_d
!  end interface
!  private :: writeidf_i_r, writeidf_r_r, writeidf_r_d
  interface writeidf
     module procedure :: writeidf_i1_r8
     module procedure :: writeidf_i2_r8
     module procedure :: writeidf_i4_r8
     module procedure :: writeidf_i8_r8
     module procedure :: writeidf_r4_r8
     module procedure :: writeidf_r8_r8
  end interface
  private :: writeidf_i1_r8, writeidf_i2_r8, writeidf_i4_r8, writeidf_i8_r8, writeidf_r4_r8, writeidf_r8_r8

  interface readasc
    module procedure :: readasc_r_r
    module procedure :: readasc_r_d
  end interface
  private :: readasc_r_r, readasc_r_d

  interface writeasc
    module procedure :: writeasc_i4_r4
    module procedure :: writeasc_i4_r8
    module procedure :: writeasc_r4_r8
    module procedure :: writeasc_r8_r8
  end interface
  private :: writeasc_i4_r4, writeasc_i4_r8, writeasc_r4_r8, writeasc_r8_r8

  interface readflt
    module procedure :: readflt_i1
    module procedure :: readflt_i4
    module procedure :: readflt_r4
  end interface
  private :: readflt_i1, readflt_i4, readflt_r4
  
  interface writeflt
    module procedure :: writeflt_i1_r4
    module procedure :: writeflt_i4_r4
    module procedure :: writeflt_r4_r4
    module procedure :: writeflt_i1_r8
    module procedure :: writeflt_i4_r8
    module procedure :: writeflt_r4_r8
    module procedure :: writeflt_r8_r8
  end interface
  private :: writeflt_i1_r4, writeflt_i4_r4, writeflt_r4_r4
  private :: writeflt_i1_r8, writeflt_i4_r8, writeflt_r4_r8, writeflt_r8_r8
  
  interface addboundary
    module procedure :: addboundary_i
    module procedure :: addboundary_r
    module procedure :: addboundary_d
    module procedure :: addboundary_i_list
  end interface
  private :: addboundary_i, addboundary_r, addboundary_d, addboundary_i_list

  interface calc_unique
    module procedure :: calc_unique_i
    module procedure :: calc_unique_r
  end interface calc_unique
  private :: calc_unique_i, calc_unique_r

  interface get_grid_bb
    module procedure :: get_r4grid_bb
  end interface get_grid_bb
  private :: get_r4grid_bb
  
  interface get_bb_extent
    module procedure :: get_bb_extent_r8
  end interface get_bb_extent
  private :: get_bb_extent_r8
  
  interface writetofile
    module procedure :: writetofile_i4
    module procedure :: writetofile_r4
    module procedure :: writetofile_r8
  end interface writetofile
  private :: writetofile_i4, writetofile_r4, writetofile_r8

  interface renumber
    module procedure :: renumber_i4
  end interface renumber
  private :: renumber_i4
  
  interface ta
    module procedure :: ta_i1
    module procedure :: ta_i2
    module procedure :: ta_i4
    module procedure :: ta_i8
    module procedure :: ta_r4
    module procedure :: ta_r8
    module procedure :: ta_c
  end interface
  private :: ta_i1, ta_i2, ta_i4, ta_i8, ta_r4, ta_r8

  type tPol
    integer(I4B) :: id = 0
    integer(I4B) :: n = 0
    real(R4B) :: xmin =  huge(0.), ymin =  huge(0.)
    real(R4B) :: xmax = -huge(0.), ymax = -huge(0.)
    real(R8B), dimension(:,:), allocatable :: xy
  end type

  type tBbX
    real(R8B) :: xll = huge(R8ZERO)
    real(R8B) :: xur = R8ZERO
    real(R8B) :: yur = huge(R8ZERO)
    real(R8B) :: yll = R8ZERO
    real(R8B) :: cs = R8ZERO
  contains
    procedure :: init  => tBbX_init
    procedure :: read  => tBbX_read
    procedure :: write => tBbX_write
    procedure :: match => tBbX_match
  end type tBbX
  
  type tBb
    integer(I4B) :: ic0  = huge(I4ZERO)
    integer(I4B) :: ic1  = I4ZERO
    integer(I4B) :: ir0  = huge(I4ZERO)
    integer(I4B) :: ir1  = I4ZERO
    integer(I4B) :: ncol = I4ZERO
    integer(I4B) :: nrow = I4ZERO
  contains
    procedure :: init  => tBb_init
    procedure :: read  => tBb_read
    procedure :: write => tBb_write
  end type tBb
  
  type tBbObj
    type(tBB)  :: prent_bbi
    type(tBB)  :: child_bbi
    type(tBBX) :: prent_bbx
    type(tBBX) :: child_bbx
  contains
    procedure :: init => tBbObj_init
    procedure :: set  => tBbObj_set
  end type tBbObj

  type tI4grid
     integer(I4B), dimension(:,:), allocatable :: x
     real(R4B) :: xll, yll, cs, nodata
     type(tBb), pointer :: bb => null()
  end type tI4grid
  public :: tI4grid
  !
  integer(I4B), parameter :: i_no_compress    = 0
  integer(I4B), parameter :: i_line_compress  = 1
  integer(I4B), parameter :: i_cgrid_compress = 2
  !
  type tGrid
    integer(I4B) :: i_data_type = 0
    !
    integer(I4B) :: nnod_dat = 0
    integer(I4B), dimension(:), allocatable :: nod_dat
    !
    integer(I4B) :: nc = 0
    integer(I4B) :: nr = 0
    !
    type(tBBx) :: bbx
    !
    integer(I1B), dimension(:,:), allocatable :: xi1
    integer(I2B), dimension(:,:), allocatable :: xi2
    integer(I4B), dimension(:,:), allocatable :: xi4
    integer(I8B), dimension(:,:), allocatable :: xi8
    real(R4B),    dimension(:,:), allocatable :: xr4
    real(R8B),    dimension(:,:), allocatable :: xr8
    !
    integer(I1B), allocatable :: mvi1
    integer(I2B), allocatable :: mvi2
    integer(I4B), allocatable :: mvi4
    integer(I8B), allocatable :: mvi8
    real(R4B),    allocatable :: mvr4
    real(R8B),    allocatable :: mvr8
    !
    integer(I4B)                              :: nmask
    integer(I4B), dimension(:,:), allocatable :: mask
  contains
    procedure :: set_const     => tGrid_set_const
    procedure :: set_mv        => tGrid_set_mv
    procedure :: set_arr       => tGrid_set_arr
    procedure :: set_val       => tGrid_set_val
    procedure :: set_nod_dat   => tGrid_set_nod_dat
    procedure :: get           => tGrid_get
    procedure :: get_nod_dat   => tGrid_get_nod_dat
    procedure :: init          => tGrid_init
    procedure :: clean         => tGrid_clean
    procedure :: clean_xi      => tGrid_clean_xi
    procedure :: clean_nod_dat => tGrid_clean_nod_dat
    procedure :: any_neg       => tGrid_any_neg
    procedure :: any_pos       => tGrid_any_pos
    procedure :: count_data    => tGrid_count_data
    ! mask:
    procedure :: set_mask      => tGrid_set_mask
    procedure :: get_mask      => tGrid_get_mask
    procedure :: write_mask    => tGrid_write_mask
    procedure :: read_mask     => tGrid_read_mask
  end type tGrid
  !
  type tGridArr
    integer(I4B) :: n = 0
    type(tGrid), dimension(:), pointer :: x => null()
  contains
    procedure :: init  => tGridArr_init
    procedure :: clean => tGridArr_clean
  end type
  !
  type tUnp
    integer :: n = 0

    integer :: ic0 = huge(0)
    integer :: ic1 = 0
    integer :: ir0 = huge(0)
    integer :: ir1 = 0
    !
    integer :: ncol = 0
    integer :: nrow = 0
    !
    integer :: nbnd = 0
    integer, dimension(:,:), allocatable :: bnd

    integer, dimension(:), allocatable :: out_itnn
    integer, dimension(:,:), allocatable :: out_bndnn
    double precision, dimension(:), allocatable :: out_d
    !
    ! ingoing links
    integer :: nin = 0
    integer, dimension(:), allocatable :: in_map
    integer, dimension(:,:), allocatable :: in_bndnn
    double precision, dimension(:), allocatable :: in_d
    integer, dimension(:), allocatable :: in_flag
    integer, dimension(:), allocatable :: in_srcp
    integer, dimension(:), allocatable :: in_srci

    ! global all-2-all connection for the first node
    integer, dimension(:,:), allocatable :: all_bndnn
    double precision, dimension(:), allocatable :: all_d

    integer, dimension(:,:), allocatable :: bndmap
  end type tUnp

  type tTimeSeries
    logical, pointer                    :: act    => null()
    logical, pointer                    :: read  => null()
    character(len=MXSLEN), pointer      :: rawhdr => null()
    character(len=MXSLEN), pointer      :: raw    => null()
    character(len=MXSLEN), pointer      :: id     => null()
    real(R8B), pointer                  :: x      => null()
    real(R8B), pointer                  :: y      => null()
    integer(I4B), pointer               :: ic     => null()
    integer(I4B), pointer               :: ir     => null()
    integer(I4B), pointer               :: im     => null()
    integer(I4B), dimension(:), pointer :: nod    => null()
    real(R8B), pointer                  :: glev   => null()
    real(R8B), dimension(:,:), pointer  :: val    => null()
    integer(I4B), pointer               :: sm_corr=> null()
    integer(I4B), pointer               :: nlay   => null()
  contains
    procedure :: clean => timeseries_clean
  end type tTimeSeries
  
  save

  contains

  subroutine grid_load_imbalance(xi4, mvi4, wi4, imbal, nparts, ximbal)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), dimension(:,:), intent(in) :: xi4
    integer(I4B)                , intent(in) :: mvi4
    integer(I4B), dimension(:,:), intent(in), optional :: wi4
    real(R8B),                    intent(out) :: imbal
    integer(I4B),                 intent(out) :: nparts
    real(R4B), dimension(:,:),    allocatable, intent(inout), optional :: ximbal
    ! -- local
    integer(I4B) :: nc, nr, ic, ir, minid, maxid, id, ip, ngid, n, i
    integer(I4B), dimension(:), allocatable :: ids
    integer(I4B), dimension(:,:), allocatable :: wi4_loc
    real(R8B), dimension(:), allocatable :: load, loadimbal
    real(R8B) :: totload, tpwgts
! ------------------------------------------------------------------------------
    nc = size(xi4,1); nr = size(xi4,2)
    !
    minid = huge(minid); maxid = -huge(maxid)
    do ir = 1, nr; do ic = 1, nc
      id = xi4(ic,ir)
      if (id /= mvi4) then
        minid = min(minid,id); maxid = max(maxid,id)
      end if
    end do; end do
    !
    ngid = maxid - minid + 1
    !
    allocate(ids(maxid)); ids = 0
    do ir = 1, nr; do ic = 1, nc
      id = xi4(ic,ir)
      if (id /= mvi4) then
        ids(id - minid + 1) = 1
      end if
    end do; end do
    !
    nparts = sum(ids); tpwgts = R8ONE/real(nparts,R8B)
    !
    n = 0
    do i = 1, ngid
      if (ids(i) == 1) then
        n = n + 1
        ids(i) = n
      end if
    end do
    !
    allocate(load(nparts)); load = R8ZERO
    if (present(wi4)) then
      allocate(wi4_loc,source=wi4)
    else
      allocate(wi4_loc(nc,nr))
      wi4_loc = 1
    end if
    !
    do ir = 1, nr; do ic = 1, nc
      id = xi4(ic,ir)
      if (id /= mvi4) then
        i = ids(id - minid + 1)
        load(i) = load(i) + real(wi4_loc(ic,ir),R8B)
      end if
    end do; end do
    !
    totload = sum(load)
    allocate(loadimbal, source=load)
    do ip = 1, nparts
      loadimbal(ip) = loadimbal(ip)/(totload*tpwgts)
    end do
    imbal = real(maxval(loadimbal),R4B)
    !
    if (present(ximbal)) then
      if (allocated(ximbal)) deallocate(ximbal)
      allocate(ximbal(nc,nr)); ximbal = R4ZERO
      do ir = 1, nr; do ic = 1, nc
        id = xi4(ic,ir)
        if (id /= mvi4) then
          i = ids(id - minid + 1)
          ximbal(ic,ir) = real(loadimbal(i),R4B)
        end if
      end do; end do
    end if
    ! clean up
    deallocate(ids, wi4_loc, load, loadimbal)
    !
    return
  end subroutine grid_load_imbalance
  
  subroutine extrapolate_mfo(xr4, mvr4)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    real(R4B), dimension(:,:), intent(inout) :: xr4
    real(R4B),                    intent(in) :: mvr4
    ! -- local
    logical :: lprocess
    integer(I4B), parameter :: nsten = 9
    integer(I4B), dimension(nsten) :: idx
    integer(I4B), dimension(:,:), allocatable ::sten, i4wk
    integer(I4B) :: ir, ic, jr, jc, nr, nc, i, n, iter, cnt
    real(R4B) :: r4v
    real(R4B), dimension(nsten) :: r4a
    real(R4B), dimension(:,:), allocatable :: r4wk
! ------------------------------------------------------------------------------
    !
    call logmsg('BEGIN focal statistics...')
    !
    nr = size(xr4,2); nc = size(xr4,1)
    allocate(i4wk(nc,nr), r4wk(nc,nr)); i4wk = 0
    !
    iter = 0
    do while(.true.)
      iter = iter + 1
      !call logmsg('Iteration: '//ta([iter],'(i4.4)'))
      !
      r4wk = R4ZERO
      do ir = 1, nr; do ic = 1, nc
        lprocess = .true.
        if (iter > 1) then
          if (i4wk(ic,ir) == 0) then
            lprocess = .false.
          end if
        else
          if (xr4(ic,ir) /= mvr4) then
            lprocess = .false.
          end if
        end if
        !
        if (lprocess) then
          call get_stencil(ic, ir, nc, nr, nsten, sten)
          n = 0
          do i = 2, nsten
            jc = sten(1,i); jr = sten(2,i)
            if ((jc == 0).or.(jr == 0)) cycle
            r4v = xr4(jc,jr)
            if (r4v /= mvr4) then
              n = n + 1
              idx(n) = n
              r4a(n) = r4v
            end if
          end do
          !
          if (n > 0) then
            i4wk(ic,ir) = 1
            select case(n)
            case(1,2)
              r4wk(ic,ir) = r4a(1)
            case default
              call quicksort_r(r4a, idx, n)
              !r4v = get_most_freq([1., 2., 2., 3., 4., 4., 4., 5.], 8)
              r4v = get_most_freq(r4a, n)
              r4wk(ic,ir) = r4v
            end select
          end if
        end if
      end do; end do
      !
      ! set the values
      do ir = 1, nr; do ic = 1, nc
        if (i4wk(ic,ir) == 1) then
          xr4(ic,ir) = r4wk(ic,ir)
        end if
      end do; end do
      !
      ! set the target missing values
      do ir = 1, nr; do ic = 1, nc
        if (i4wk(ic,ir) == 1) then
          call get_stencil(ic, ir, nc, nr, nsten, sten)
          do i = 2, nsten
            jc = sten(1,i); jr = sten(2,i)
            if ((jc == 0).or.(jr == 0)) cycle
            r4v = xr4(jc,jr)
            if (r4v == mvr4) then
              if (i4wk(jc,jr) == 0) then
                i4wk(jc,jr) = 2
              end if
            end if
          end do
        end if
      end do; end do
      !
      cnt = 0
      do ir = 1, nr; do ic = 1, nc
        if (i4wk(ic,ir) == 1) then
          i4wk(ic,ir) = 0
        end if
        if (i4wk(ic,ir) == 2) then
          i4wk(ic,ir) = 1
          cnt = cnt + 1
        end if
      end do; end do
      !
      if (cnt == 0) exit
    end do
    !
    call logmsg('END focal statistics...')
    !
    return
  end subroutine extrapolate_mfo
  
  function get_most_freq(freq, nfreq) result(most_freq)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    real(R4B), dimension(:), intent(in) :: freq
    integer(I4B), intent(in) :: nfreq
    real(R4B) :: most_freq
    ! -- local
    integer(I4B) :: i, mi, ni
! ------------------------------------------------------------------------------
   
     ni = 1  !number of unique
     mi = ni !max. number of unique
     most_freq = freq(ni)
     
     do i = 2, nfreq
       if (freq(i) /= freq(i-1)) then
         if (ni > mi) then
           most_freq = freq(i-1)
           mi = ni
         end if
         ni = 1
       else
         ni = ni + 1
        end if
     end do
     !
     !test final
     if(ni > mi) most_freq = freq(nfreq)
     !
    return
  end function get_most_freq
  
  subroutine get_compiler(cdate, cversion)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use iso_fortran_env, only: compiler_options, compiler_version
    ! -- dummy
    character(len=MXSLEN), intent(out) :: cdate
    character(len=MXSLEN), intent(out) :: cversion 
    ! -- local
    logical :: cfound
! ------------------------------------------------------------------------------
    cfound = .false.
#ifdef __GFORTRAN__
    cfound =.true.
#endif
#ifdef __INTEL_COMPILER
    cfound =.true.
#endif
#ifdef _CRAYFTN
    cfound =.true.
#endif
    !
    ! -- set compiler strings
    if (.not.cfound) then
      cdate = '??? ?? ???? ??:??:??'
      cversion = 'UNKNOWN COMPILER'
    else
      cdate = trim(adjustl(__DATE__//' '//__TIME__))
      cversion = trim(adjustl(compiler_version()))
    end if
    !
    return
  end subroutine get_compiler
  !
  subroutine coarse_to_fine_grid(bbx_c, bbx_f, xi4c, xr4c, xi4f, xr4f)
! ******************************************************************************
  
  
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    type(tBBX), intent(in) :: bbx_c
    type(tBBX), intent(in) :: bbx_f
    !
    integer(I4B), dimension(:,:), intent(in), optional :: xi4c
    real(R4B),    dimension(:,:), intent(in), optional :: xr4c
    !
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: xi4f
    real(R4B),    dimension(:,:), allocatable, intent(inout), optional :: xr4f
    !
    ! -- local
    integer(I4B) :: i, ic, ir, jc, jr, nc, nr, mc, mr, n
    real(R8B) :: x, y
! ------------------------------------------------------------------------------
    if (bbx_c%cs < bbx_f%cs) then
      call errmsg('coarse_to_fine_grid: invalid cell sizes.')
    end if
    if (.not.bbx_c%match(bbx_f)) then
      call errmsg('coarse_to_fine_grid: non matching grids.')
    end if
    n = 0
    if (present(xi4c)) then
      if (.not.present(xi4f)) call errmsg('coarse_to_fine_grid: xi4f not found.')
      n = n + 1
    end if
    if (present(xr4c)) then
      if (.not.present(xr4f)) call errmsg('coarse_to_fine_grid: xr4f not found.')
      n = n + 1
    end if
    if (n /= 1) then
      call errmsg('coarse_to_fine_grid: no input found.')
    end if
    !
    nc = (bbx_c%xur - bbx_c%xll)/bbx_c%cs
    nr = (bbx_c%yur - bbx_c%yll)/bbx_c%cs
    mc = (bbx_f%xur - bbx_f%xll)/bbx_f%cs
    mr = (bbx_f%yur - bbx_f%yll)/bbx_f%cs
    !
    if (present(xi4c)) then
      if ((size(xi4c,1) /= nc).or.(size(xi4c,2) /= nr)) then
        call errmsg('coarse_to_fine_grid: invalid dimensions.')
      end if
      if (allocated(xi4f)) deallocate(xi4f)
      allocate(xi4f(mc,mr))
    end if
    if (present(xr4c)) then
      if ((size(xr4c,1) /= nc).or.(size(xr4c,2) /= nr)) then
        call errmsg('coarse_to_fine_grid: invalid dimensions.')
      end if
      if (allocated(xr4f)) deallocate(xr4f)
      allocate(xr4f(mc,mr))
    end if
    !
    if (present(xi4c)) then
      do ir = 1, mr; do ic = 1, mc
        call get_xy(x, y, ic, ir,  bbx_f%xll, bbx_f%yur, bbx_f%cs)
        call get_icr(jc, jr, x, y, bbx_c%xll, bbx_c%yur, bbx_c%cs)
        if ((jc >= 1).and.(jc <= nc).and.&
            (jr >= 1).and.(jr <= nr)) then
          xi4f(ic,ir) = xi4c(jc,jr)
        else
          call errmsg('coarse_to_fine_grid')
        end if
      end do; end do
    end if
    
    return
  end subroutine coarse_to_fine_grid
  
  subroutine get_stencil(ic, ir, nc, nr, np, sten)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), intent(in) :: ic
    integer(I4B), intent(in) :: ir
    integer(I4B), intent(in) :: nc
    integer(I4B), intent(in) :: nr
    integer(I4B), intent(in) :: np
    integer(I4B), dimension(:,:), intent(inout), allocatable :: sten
    ! -- local
    integer(I4B) :: i, jc, jr
! ------------------------------------------------------------------------------
    if ((np /= 5).and.(np /= 9)) then
      call errmsg('get_stencil: only 5-point of 9-point stencil allowed.')
    end if
    !
    if (allocated(sten)) deallocate(sten)
    allocate(sten(2,np))
    !
    do i = 1, np
      jc = ic + sten_dir(1,i)
      jr = ir + sten_dir(2,i)
      !
      if ((jc < 1).or.(jc > nc).or.(jr < 1).or.(jr > nr)) then
        jc = 0; jr = 0
      end if
      sten(1,i) = jc; sten(2,i) = jr
    end do
    !
  end subroutine get_stencil
  
  function linear_interpolation(x, x1, x2, v1, v2) result(v)
! ******************************************************************************
! See: https://en.wikipedia.org/wiki/Linear_interpolation
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    real(R8B), intent(in) :: x
    real(R8B), intent(in) :: x1
    real(R8B), intent(in) :: x2
    real(R8B), intent(in) :: v1
    real(R8B), intent(in) :: v2
    real(R8B) :: v ! result
    
    ! -- local
    real(R8B) :: w1, w2, dx2x1
! ------------------------------------------------------------------------------
    dx2x1 = x2 - x1
    !
    w1 = (x2 - x )/dx2x1
    w2 = (x  - x1)/dx2x1
    !
    v = v1*w1 + v2*w2
    !
  end function linear_interpolation
  !
  function bilinear_interpolation(xp, xll, xur, va, mv) result(v)
! ******************************************************************************
! See: https://en.wikipedia.org/wiki/Bilinear_interpolation
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    !  y2--va(2) ----- va(3) 
    !         |           |
    !         |           |
    !         |           |
    !  y1--va(1) ----- va(4)
    !         |           |
    !        x1          x2
    !
    ! -- dummy
    real(R8B), dimension(2,1), intent(in) :: xp
    real(R8B), dimension(2,1), intent(in) :: xll
    real(R8B), dimension(2,1), intent(in) :: xur
    real(R8B), dimension(4),   intent(in) :: va
    real(R8B),                 intent(in) :: mv
    real(R8B) :: v ! result
    ! -- local
    integer(I4B), dimension(1) :: mloc
    integer(I4B) :: i, nmv, i_nearest, i0, i1
    integer(I4B), dimension(4) :: flg
    real(R8B), dimension(4) :: va_used
    real(R8B), dimension(2) :: a, b, s
    real(R8B) :: x, y, x1, x2, y1, y2, vv, xproj, yproj
    real(R8B) :: dx, dy, dxy, d
    real(R8B) :: dx2x, dxx1, dx2x1
    real(R8B) :: dy2y, dyy1, dy2y1
    real(R8B) :: f, w11, w12, w21, w22
! ------------------------------------------------------------------------------
    x  =  xp(1,1);  y =  xp(2,1)
    x1 = xll(1,1); y1 = xll(2,1)
    x2 = xur(1,1); y2 = xur(2,1)
    !
    ! compute the weights
    dx2x = x2 - x; dxx1 = x - x1; dx2x1 = x2 - x1
    dy2y = y2 - y; dyy1 = y - y1; dy2y1 = y2 - y1
    !
    if (dx2x1 /= dy2y1) then
      call errmsg('bilinear_interpolation: no equidistant.')
    end if
    !
    ! checks
    if ((x < x1).or.(x > x2).or.(y < y1).or.(y > y2)) then
      call errmsg('bilinear_interpolation: coordinate out of box.')
    end if
    va_used = va
    !
    ! treat missing values
    nmv = 0; flg = 0
    do i = 1, 4
      if (va(i) == mv) then
        nmv = nmv + 1
        flg(i) = 1
      end if
    end do
    !
    select case(nmv)
    case(0)
      !
      f = dx2x1*dy2y1 + R8TINY
      !
      w11 = dx2x*dy2y/f
      w12 = dx2x*dyy1/f
      w21 = dxx1*dy2y/f
      w22 = dxx1*dyy1/f
      !
      ! evaluate the solution
      v = w11*va(1) + w12*va(2) + w21*va(3) + w22*va(4)
    case(1)
      if (dxx1 < dx2x) then
        if (dyy1 < dy2y) then
          i_nearest = 1
        else
          i_nearest = 2
        end if
      else
        if (dyy1 < dy2y) then
          i_nearest = 4
        else
          i_nearest = 3
        end if
      end if
      mloc = maxloc(flg); i = mloc(1)
      v = va(i)
      if (i == i_nearest) then
        v = mv
      else
        if (flg(1) == 1) then
          if ((i_nearest == 2).or.(i_nearest == 3)) then
            vv = linear_interpolation(y, y1, y2, va(4), va(3))
            v  = linear_interpolation(x, x1, x2, va(2), vv)
          else
            v = linear_interpolation(y, y1, y2, va(4), va(3))
          end if
        end if
        if (flg(4) == 1) then
          if ((i_nearest == 2).or.(i_nearest == 3)) then
            vv = linear_interpolation(y, y1, y2, va(1), va(2))
            v  = linear_interpolation(x, x1, x2, vv,    va(2))
          else
            v = linear_interpolation(y, y1, y2, va(1), va(2))
          end if
        end if
        if (flg(2) == 1) then
          if ((i_nearest == 1).or.(i_nearest == 4)) then
            vv = linear_interpolation(y, y1, y2, va(4), va(3))
            v  = linear_interpolation(x, x1, x2, va(1), vv)
          else
            v = linear_interpolation(y, y1, y2, va(4), va(3))
          end if
        end if
        if (flg(3) == 1) then
          if ((i_nearest == 1).or.(i_nearest == 4)) then
            vv = linear_interpolation(y, y1, y2, va(1), va(2))
            v  = linear_interpolation(x, x1, x2, vv, va(4))
          else
            v = linear_interpolation(y, y1, y2, va(1), va(2))
          end if
        end if
      end if
    case(2) ! linear interpolation
      if ((flg(2) == 1).and.(flg(3) == 1)) then
        v = linear_interpolation(x, x1, x2, va(1), va(4))
      end if
      if ((flg(1) == 1).and.(flg(3) == 4)) then
        v = linear_interpolation(x, x1, x2, va(2), va(3))
      end if
      if ((flg(1) == 1).and.(flg(2) == 1)) then
        v = linear_interpolation(y, y1, y2, va(4), va(3))
      end if
      if ((flg(4) == 1).and.(flg(3) == 1)) then
        v = linear_interpolation(y, y1, y2, va(1), va(2))
      end if
      if (((flg(2) == 1).and.(flg(4) == 1)) .or. &
          ((flg(1) == 1).and.(flg(3) == 1))) then
        !
        if ((flg(1) == 1).and.(flg(3) == 1)) then
          i0 = 1; i1 = 3
          a(1) = dx2x1; a(2) = dy2y1
          b(1) = dxx1;  b(2) = dyy1
        else
          i0 = 2; i1 = 4
          a(1) = dx2x1; a(2) = -dy2y1
          b(1) = dxx1;  b(2) = -dy2y
        end if
        dxy = sqrt(a(1)**2 + a(2)**2)
        f = (a(1)*b(1) + a(2)*b(2))/(dxy**2)
        s = a*f; d = sqrt(s(1)**2 + s(2)**2)
        v = linear_interpolation(d, R8ZERO, dxy, va(i0), va(i1))
      end if
    case(3)
      mloc = minloc(flg); i = mloc(1)
      v = va(i)
    case(4)
      v = mv
    end select
    !
  end function bilinear_interpolation
    
  subroutine tBb_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tBb) :: this
! ------------------------------------------------------------------------------
    this%ic0  = huge(I4ZERO)
    this%ic1  = I4ZERO
    this%ir0  = huge(I4ZERO)
    this%ir1  = I4ZERO
    this%ncol = I4ZERO
    this%nrow = I4ZERO
    !
    return
  end subroutine tBb_init
    
  subroutine tBb_read(this, iu, pos)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tBb) :: this
    !
    integer(I4B), intent(in) :: iu
    integer(I8B), intent(in), optional :: pos
    ! -- local
! ------------------------------------------------------------------------------
    if (.not.present(pos)) then
      read(iu) this%ic0, this%ic1, this%ir0, this%ir1, this%ncol, this%nrow
    else
      read(iu,pos=pos) this%ic0, this%ic1, this%ir0, this%ir1, this%ncol, this%nrow
    end if
    !
    return
  end subroutine tBb_read
  
  subroutine tBb_write(this, iu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tBb) :: this
    !
    integer(I4B), intent(in) :: iu
! ------------------------------------------------------------------------------
    write(iu) this%ic0, this%ic1, this%ir0, this%ir1, this%ncol, this%nrow
    !
    return
  end subroutine tBb_write
  
  subroutine tBbX_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tBbX) :: this
! ------------------------------------------------------------------------------
    this%xll = huge(R8ZERO)
    this%xur = R8ZERO
!    this%yur = huge(R8ZERO)
!    this%yll = R8ZERO
    this%yur = R8ZERO
    this%yll = huge(R8ZERO)
    this%cs  = R8ZERO
    return
  end subroutine tBbX_init
  
  subroutine tBbX_read(this, iu, pos)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tBbX) :: this
    !
    integer(I4B), intent(in) :: iu
    integer(I8B), intent(in), optional :: pos
! ------------------------------------------------------------------------------
    if (.not.present(pos)) then
      read(iu) this%xll, this%xur, this%yur, this%yll, this%cs
    else
      read(iu,pos=pos) this%xll, this%xur, this%yur, this%yll, this%cs
    end if
    !
    return
  end subroutine tBbX_read
  
  subroutine tBbX_write(this, iu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tBbX) :: this
    !
    integer(I4B), intent(in) :: iu
! ------------------------------------------------------------------------------
    write(iu) this%xll
    write(iu) this%xur
    write(iu) this%yur
    write(iu) this%yll
    write(iu) this%cs
    !
    return
  end subroutine tBbX_write
  
  function tBbX_match(this, bbx) result(match)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tBbX) :: this
    type(tBbX) :: bbx
    logical :: match
! ------------------------------------------------------------------------------
    match = .true.
    if (bbx%cs > this%cs) then
      if (mod(bbx%cs, this%cs) > R8TINY) then
        match = .false.; return
      end if
    else
      if (mod(this%cs, bbx%cs) > R8TINY) then
        match = .false.; return
      end if
    end if
    !
    if (abs(mod(this%xll - bbx%xll, this%cs)) > R8TINY) then
      match = .false.; return
    end if
    if (abs(mod(this%xur - bbx%xur, this%cs)) > R8TINY) then
      match = .false.; return
    end if
    if (abs(mod(this%yur - bbx%yur, this%cs)) > R8TINY) then
      match = .false.; return
    end if
    if (abs(mod(this%yll - bbx%yll, this%cs)) > R8TINY) then
      match = .false.; return
    end if
    !
    return
  end function tBbX_match
  
  subroutine renumber_i4(x, mv, nid, f_csv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), dimension(:,:), intent(inout) :: x
    integer(I4B), intent(in) :: mv
    integer(I4B), intent(inout) :: nid
    character(len=*), intent(in), optional :: f_csv
    ! -- local
    character(len=MXSLEN) :: f
    integer(I4B), dimension(:), allocatable :: i4wk
    integer(I4B) :: ir, ic, nc, nr, xmin, xmax, n_max, n, id, jd, i, iu
! ------------------------------------------------------------------------------
    nc = size(x,1); nr = size(x,2)
    !
    xmin = huge(I4ZERO); xmax = -huge(I4ZERO)
    do ir = 1, nr
      do ic = 1, nc
        id = x(ic,ir)
        if (id /= mv) then
          xmin = min(xmin, id); xmax = max(xmax, id)
        end if
      end do
    end do
    !
    n_max = xmax - xmin + 1
    if (n_max == 0) then
      call errmsg('renumber_i4: nothing to renumber.')
    end if
    !
    allocate(i4wk(n_max)); i4wk = 0
    do ir = 1, nr
      do ic = 1, nc
        id = x(ic,ir)
        if (id /= mv) then
          jd = id - xmin + 1; i4wk(jd) = 1
        end if
      end do
    end do
    !
    ! set the new ids
    n = 0
    do jd = 1, n_max
      if (i4wk(jd) == 1) then
        n = n + 1
        nid = nid + 1; i4wk(jd) = nid
      end if
    end do
    !
    do ir = 1, nr
      do ic = 1, nc
        id = x(ic,ir)
        if (id /= mv) then
          jd = id - xmin + 1; x(ic,ir)= i4wk(jd)
        end if
      end do
    end do
    !
    !call logmsg('# max. IDs: '//ta((/n_max/))//'; # IDs set: '//ta((/n/)))
    !
    if (present(f_csv)) then
      f = f_csv; call open_file(f, iu, 'w')
      write(iu,'(a)') 'id_new, id_old'
      do jd = 1, n_max
        if (i4wk(jd) /= 0) then
          write(iu,'(a)') ta(arr=(/i4wk(jd),jd/),sep_in=',')
        end if
      end do
      close(iu)
    end if
    !
    deallocate(i4wk)
    !
    return
  end subroutine renumber_i4
  
  subroutine tGrid_set_const(this, i1v, i2v, i4v, i8v, r4v, r8v)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I1B), intent(in), optional :: i1v
    integer(I2B), intent(in), optional :: i2v
    integer(I4B), intent(in), optional :: i4v
    integer(I8B), intent(in), optional :: i8v
    real(R4B),    intent(in), optional :: r4v
    real(R8B),    intent(in), optional :: r8v
    
    ! -- local
    integer(I4B) :: ir, ic
! ------------------------------------------------------------------------------
    if (present(i1v).and.allocated(this%xi1)) then
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xi1(ic,ir) = i1v
      end do; end do
    end if
    if (present(i2v).and.allocated(this%xi2)) then
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xi2(ic,ir) = i2v
      end do; end do
    end if
    if (present(i4v).and.allocated(this%xi4)) then
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xi4(ic,ir) = i4v
      end do; end do
    end if
    if (present(i8v).and.allocated(this%xi8)) then
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xi8(ic,ir) = i8v
      end do; end do
    end if
    if (present(r4v).and.allocated(this%xr4)) then
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xr4(ic,ir) = r4v
      end do; end do
    end if
    if (present(r8v).and.allocated(this%xr8)) then
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xr8(ic,ir) = r8v
      end do; end do
    end if
    !
    return
    end subroutine tGrid_set_const
    
  subroutine tGrid_set_mask(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    ! -- local
    logical :: l_store_first, l_store_last
    integer(I1B), dimension(:,:), allocatable :: wrk
    integer(I4B) :: ic, ir, n, ncd, nrd
    real(R4B) :: f
! ------------------------------------------------------------------------------
    allocate(wrk(this%nc, this%nr))
    wrk = 0
    !
    select case(this%i_data_type)
    case(i_i1)
      do ir = 1, this%nr; do ic = 1, this%nc
        if (this%xi1(ic,ir) /= this%mvi1) wrk(ic,ir) = 1
      end do; end do
    case(i_i2)
      do ir = 1, this%nr; do ic = 1, this%nc
        if (this%xi2(ic,ir) /= this%mvi2) wrk(ic,ir) = 1
      end do; end do
    case(i_i4)
      do ir = 1, this%nr; do ic = 1, this%nc
        if (this%xi4(ic,ir) /= this%mvi4) wrk(ic,ir) = 1
      end do; end do
    case(i_i8)
      do ir = 1, this%nr; do ic = 1, this%nc
        if (this%xi8(ic,ir) /= this%mvi8) wrk(ic,ir) = 1
      end do; end do
    case(i_r4)
      do ir = 1, this%nr; do ic = 1, this%nc
        if (this%xr4(ic,ir) /= this%mvr4) wrk(ic,ir) = 1
      end do; end do
    case(i_r8)
      do ir = 1, this%nr; do ic = 1, this%nc
        if (this%xr8(ic,ir) /= this%mvr8) wrk(ic,ir) = 1
      end do; end do
    end select
    !
    ! count
    n = 0
    do ir = 1, this%nr
      do ic = 1, this%nc
        if (wrk(ic,ir) == 1) n = n + 1
      end do
    end do
    !
    ! scan along the columns
    ncd = 0
    do ir = 1, this%nr; do ic = 1, this%nc
      if (wrk(ic,ir) == 1) then
        if (ic == 1) then
          ncd = ncd + 1
        else
          if (wrk(ic-1,ir) == 0) then
            ncd = ncd + 1
          end if
        end if
      end if
    end do; end do
      !
    ! scan along the rows
    nrd = 0
    do ic = 1, this%nc; do ir = 1, this%nr
      if (wrk(ic,ir) == 1) then
        if (ir == 1) then
          nrd = nrd + 1
        else
          if (wrk(ic,ir-1) == 0) then
            nrd = nrd + 1
          end if
        end if
      end if
    end do; end do
    !
    if (ncd <= nrd) then
      if (ncd > 0) then
        allocate(this%mask(2,ncd))
      end if
      !
      this%nmask = 0
      do ir = 1, this%nr
        do ic = 1, this%nc
          l_store_first = .false.
          l_store_last  = .false.
          !
          if (wrk(ic,ir) == 1) then
            if (ic == 1) then
              this%nmask = this%nmask + 1; l_store_first = .true.
            else
              if (wrk(ic-1,ir) == 0) then
                this%nmask = this%nmask + 1; l_store_first = .true.
              end if
            end if
            if (ic == this%nc) then
              l_store_last = .true.
            else
              if (wrk(ic+1,ir) == 0) then
                l_store_last = .true.
              end if
            end if
          end if
          !
          if (l_store_first.or.l_store_last) then
            call icrl_to_node(n, ic, ir, 1, this%nc, this%nr)
          end if
          if (l_store_first) this%mask(1,this%nmask) = n
          if (l_store_last)  this%mask(2,this%nmask) = n
        end do
      end do
    else
      if (nrd > 0) then
        allocate(this%mask(2,nrd))
      end if
      !
      this%nmask = 0
      do ic = 1, this%nc
        do ir = 1, this%nr
          l_store_first = .false.
          l_store_last  = .false.
          !
          if (wrk(ic,ir) == 1) then
            if (ir == 1) then
              this%nmask = this%nmask + 1; l_store_first = .true.
            else
              if (wrk(ic,ir-1) == 0) then
                this%nmask = this%nmask + 1; l_store_first = .true.
              end if
            end if
            if (ir == this%nr) then
              l_store_last = .true.
            else
              if (wrk(ic,ir+1) == 0) then
                l_store_last = .true.
              end if
            end if
          end if
          !
          if (l_store_first.or.l_store_last) then
            call icrl_to_node(n, ic, ir, 1, this%nc, this%nr)
          end if
          if (l_store_first) this%mask(1,this%nmask) = n
          if (l_store_last)  this%mask(2,this%nmask) = n
        end do
      end do
    end if
    !
    !call logmsg('Line compression: stored '//ta([n])//' cells as 2 x '// &
    !  ta([this%nmask])//' = '//ta([2*this%nmask])//' cells...')
    !
    deallocate(wrk)
    !
    return
  end subroutine tGrid_set_mask

  subroutine tGrid_get_mask(this, mask)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I4B), dimension(:,:), intent(inout), allocatable :: mask
    !
    ! -- local
    integer(I4B) :: i, n, ic, ir, ic0, ir0, ic1, ir1, idummy
! ------------------------------------------------------------------------------
    if (this%nmask == 0) then
      call logmsg('tGrid_get_mask: nothing to do, returning...')
      return
    end if
    if (.not.allocated(this%mask)) then
      call logmsg('tGrid_get_mask: program error 1.')
    end if
    !
    if (allocated(mask)) deallocate(mask)
    allocate(mask(this%nc,this%nr))
    mask = 0
    do i = 1, this%nmask
       call node_to_icrl(this%mask(1,i), ic0, ir0, idummy, this%nc, this%nr)
       call node_to_icrl(this%mask(2,i), ic1, ir1, idummy, this%nc, this%nr)
       if ((ir0 /= ir1).and.(ic0 /= ic1)) then
         call errmsg('tGrid_get_mask: program error.')
       end if
       !
       if (ir0 /= ir1) then ! along the rows
         ic = ic0
         do ir = ir0, ir1
           mask(ic,ir) = 1
         end do
       else ! along the columns
         ir = ir0
         do ic = ic0, ic1
           mask(ic,ir) = 1
         end do
       end if
    end do
    !
    return
  end subroutine tGrid_get_mask
  
  subroutine tGrid_write_mask(this, iu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I4B), intent(in) :: iu
! ------------------------------------------------------------------------------
    !
    if ((this%nmask == 0).or.(.not.allocated(this%mask))) then
      call errmsg('tGrid_write_mask')
    end if
    !
    write(iu) this%nmask !I4B
    write(iu) this%mask  !I4B*2*mask%nmask
    !
    return
  end subroutine tGrid_write_mask
  
  subroutine tGrid_read_mask(this, iu, pos)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I4B), intent(in) :: iu
    integer(I8B), intent(in), optional :: pos
    !
    ! -- local
    integer(I8B) :: p
! ------------------------------------------------------------------------------
    if(.not.present(pos)) then
      read(iu) this%nmask !I4B
    else
      p = pos
      read(iu,pos=p) this%nmask !I4B
    end if
    !
    if (this%nmask <= 0) call errmsg('tGrid_read_mask')
    if (allocated(this%mask)) deallocate(this%mask)
    allocate(this%mask(2,this%nmask))
    !
    if(.not.present(pos)) then
      read(iu) this%mask 
    else
      p = p + sizeof(this%nmask)
      read(iu,pos=p) this%mask
    end if
    !
    return
  end subroutine tGrid_read_mask
  
  subroutine tGrid_set_mv(this, mvi1, mvi2, mvi4, mvi8, mvr4, mvr8)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I1B), intent(in), optional :: mvi1
    integer(I2B), intent(in), optional :: mvi2
    integer(I4B), intent(in), optional :: mvi4
    integer(I8B), intent(in), optional :: mvi8
    real(R4B),    intent(in), optional :: mvr4
    real(R8B),    intent(in), optional :: mvr8
    !
    ! -- local
    integer(I4B) :: n
! ------------------------------------------------------------------------------
    !
    n = 0
    if (present(mvi1)) n = n + 1
    if (present(mvi2)) n = n + 1
    if (present(mvi4)) n = n + 1
    if (present(mvi8)) n = n + 1
    if (present(mvr4)) n = n + 1
    if (present(mvr8)) n = n + 1
    if (n > 1) then
      call errmsg('tGrid_set_mv: only one array is allowed to set.')
    end if
    !
    if (present(mvi1)) then
      allocate(this%mvi1); this%mvi1 = mvi1; this%i_data_type = i_i1
    end if
    if (present(mvi2)) then
      allocate(this%mvi2); this%mvi2 = mvi2; this%i_data_type = i_i2
    end if
    if (present(mvi4)) then
      allocate(this%mvi4); this%mvi4 = mvi4; this%i_data_type = i_i4
    end if
    if (present(mvi8)) then
      allocate(this%mvi8); this%mvi8 = mvi8; this%i_data_type = i_i8
    end if
    if (present(mvr4)) then
      allocate(this%mvr4); this%mvr4 = mvr4; this%i_data_type = i_r4
    end if
    if (present(mvr8)) then
      allocate(this%mvr8); this%mvr8 = mvr8; this%i_data_type = i_r8
    end if
    !
    return
  end subroutine tGrid_set_mv
    
  subroutine tGrid_set_arr(this, xi1, xi2, xi4, xi8, xr4, xr8)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I1B), dimension(:,:), intent(in), optional :: xi1
    integer(I2B), dimension(:,:), intent(in), optional :: xi2
    integer(I4B), dimension(:,:), intent(in), optional :: xi4
    integer(I8B), dimension(:,:), intent(in), optional :: xi8
    real(R4B),    dimension(:,:), intent(in), optional :: xr4
    real(R8B),    dimension(:,:), intent(in), optional :: xr8
    !
    ! -- local
    integer(I4B) :: ir, ic, n
! ------------------------------------------------------------------------------
    !
    ! check
    n = 0
    if (present(xi1)) n = n + 1
    if (present(xi2)) n = n + 1
    if (present(xi4)) n = n + 1
    if (present(xi8)) n = n + 1
    if (present(xr4)) n = n + 1
    if (present(xr8)) n = n + 1
    if (n /= 1) then
      call errmsg('tGrid_set_arr: only one array is allowed to set.')
    end if
    !
    if (present(xi1)) then
      this%nc = size(xi1,1); this%nr = size(xi1,2)
      allocate(this%xi1(this%nc,this%nr))
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xi1(ic,ir) = xi1(ic,ir)
      end do; end do
    end if
    if (present(xi2)) then
      this%nc = size(xi2,1); this%nr = size(xi2,2)
      allocate(this%xi2(this%nc,this%nr))
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xi2(ic,ir) = xi2(ic,ir)
      end do; end do
    end if
    if (present(xi4)) then
      this%nc = size(xi4,1); this%nr = size(xi4,2)
      allocate(this%xi4(this%nc,this%nr))
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xi4(ic,ir) = xi4(ic,ir)
      end do; end do
    end if
    if (present(xi8)) then
      this%nc = size(xi8,1); this%nr = size(xi8,2)
      allocate(this%xi8(this%nc,this%nr))
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xi8(ic,ir) = xi8(ic,ir)
      end do; end do
    end if
    if (present(xr4)) then
      this%nc = size(xr4,1); this%nr = size(xr4,2)
      allocate(this%xr4(this%nc,this%nr))
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xr4(ic,ir) = xr4(ic,ir)
      end do; end do
    end if
    if (present(xr8)) then
      this%nc = size(xr8,1); this%nr = size(xr8,2)
      allocate(this%xr8(this%nc,this%nr))
      do ir = 1, this%nr; do ic = 1, this%nc
        this%xr8(ic,ir) = xr8(ic,ir)
      end do; end do
    end if
    !
    return
  end subroutine tGrid_set_arr
    
  subroutine tGrid_set_val(this, icir, i1v, i2v, i4v, i8v, r4v, r8v)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I4B), dimension(:,:), intent(in) :: icir
    integer(I1B), intent(in), optional :: i1v
    integer(I2B), intent(in), optional :: i2v
    integer(I4B), intent(in), optional :: i4v
    integer(I8B), intent(in), optional :: i8v
    real(R4B),    intent(in), optional :: r4v
    real(R8B),    intent(in), optional :: r8v
    
    ! -- local
    integer(I4B) :: i, n, ir, ic
! ------------------------------------------------------------------------------
    n = size(icir,2)
    if (n == 0) return
    !
    if (present(i1v).and.allocated(this%xi1)) then
      do i = 1, n
        ic = icir(1,i); ir = icir(2,i); this%xi1(ic,ir) = i1v
      end do
    end if
    if (present(i2v).and.allocated(this%xi2)) then
      do i = 1, n
        ic = icir(1,i); ir = icir(2,i); this%xi2(ic,ir) = i2v
      end do
    end if
    if (present(i4v).and.allocated(this%xi4)) then
      do i = 1, n
        ic = icir(1,i); ir = icir(2,i); this%xi4(ic,ir) = i4v
      end do
    end if
    if (present(i8v).and.allocated(this%xi8)) then
      do i = 1, n
        ic = icir(1,i); ir = icir(2,i); this%xi8(ic,ir) = i8v
      end do
    end if
    if (present(r4v).and.allocated(this%xr4)) then
      do i = 1, n
        ic = icir(1,i); ir = icir(2,i); this%xr4(ic,ir) = r4v
      end do
    end if
    if (present(r8v).and.allocated(this%xr8)) then
      do i = 1, n
        ic = icir(1,i); ir = icir(2,i); this%xr8(ic,ir) = r8v
      end do
    end if
    !
    return
  end subroutine tGrid_set_val
    
  subroutine tGrid_set_nod_dat(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    ! -- local
    integer(I4B), dimension(:), allocatable :: wk
    integer(I4B) :: il, ir, ic, nc, nr, n, i
! ------------------------------------------------------------------------------
    !
    il = 1; nc = this%nc; nr = this%nr
    allocate(wk(nc*nr))
    !
    this%nnod_dat = 0
    !
    if (allocated(this%xi1)) then
      if(.not.allocated(this%mvi1)) call errmsg('tGrid_set_nod_dat: i1')
      do ir = 1, nr; do ic = 1, nc
        if (this%xi1(ic,ir) /= this%mvi1) then
          this%nnod_dat = this%nnod_dat + 1; n = this%nnod_dat
          call icrl_to_node(wk(n), ic, ir, il, nc, nr)
        end if
      end do; end do
    end if
    if (allocated(this%xi2)) then
      if(.not.allocated(this%mvi2)) call errmsg('tGrid_set_nod_dat: i2')
      do ir = 1, nr; do ic = 1, nc
        if (this%xi2(ic,ir) /= this%mvi2) then
          this%nnod_dat = this%nnod_dat + 1; n = this%nnod_dat
          call icrl_to_node(wk(n), ic, ir, il, nc, nr)
        end if
      end do; end do
    end if
    if (allocated(this%xi4)) then
      if(.not.allocated(this%mvi4)) call errmsg('tGrid_set_nod_dat: i4')
      do ir = 1, nr; do ic = 1, nc
        if (this%xi4(ic,ir) /= this%mvi4) then
          this%nnod_dat = this%nnod_dat + 1; n = this%nnod_dat
          call icrl_to_node(wk(n), ic, ir, il, nc, nr)
        end if
      end do; end do
    end if
    if (allocated(this%xi8)) then
      if(.not.allocated(this%mvi8)) call errmsg('tGrid_set_nod_dat: i8')
      do ir = 1, nr; do ic = 1, nc
        if (this%xi8(ic,ir) /= this%mvi8) then
          this%nnod_dat = this%nnod_dat + 1; n = this%nnod_dat
          call icrl_to_node(wk(n), ic, ir, il, nc, nr)
        end if
      end do; end do
    end if
    if (allocated(this%xr4)) then
      if(.not.allocated(this%mvr4)) call errmsg('tGrid_set_nod_dat: r4')
      do ir = 1, nr; do ic = 1, nc
        if (this%xr4(ic,ir) /= this%mvr4) then
          this%nnod_dat = this%nnod_dat + 1; n = this%nnod_dat
          call icrl_to_node(wk(n), ic, ir, il, nc, nr)
        end if
      end do; end do
    end if
    if (allocated(this%xr8)) then
      if(.not.allocated(this%mvr8)) call errmsg('tGrid_set_nod_dat: r8')
      do ir = 1, nr; do ic = 1, nc
        if (this%xr8(ic,ir) /= this%mvr8) then
          this%nnod_dat = this%nnod_dat + 1; n = this%nnod_dat
          call icrl_to_node(wk(n), ic, ir, il, nc, nr)
        end if
      end do; end do
    end if
    !
    if (this%nnod_dat == 0) then
      deallocate(wk); return
    else
      if (allocated(this%nod_dat)) deallocate(this%nod_dat)
      allocate(this%nod_dat(this%nnod_dat))
      do i = 1, this%nnod_dat
        this%nod_dat(i) = wk(i)
      end do
      deallocate(wk)
    end if
    !
    return
  end subroutine tGrid_set_nod_dat
  
  subroutine tGrid_get(this, i_data_type)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I4B), intent(out), optional :: i_data_type
! ------------------------------------------------------------------------------
    if (present(i_data_type)) then
      if (this%i_data_type == 0) then
        call errmsg('tGrid_get: i_data_type not set.')
      end if
      i_data_type = this%i_data_type
    end if
    !
    return
  end subroutine tGrid_get
    
  subroutine tGrid_get_nod_dat(this, icir, nc, nr, n)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I4B), dimension(:,:), allocatable, intent(inout) :: icir
    integer(I4B), intent(in) :: nc
    integer(I4B), intent(in) :: nr
    integer(I4B), intent(out) :: n
    
    ! -- local
    integer(I4B) :: i, ic, ir, il
! ------------------------------------------------------------------------------
    !
    if (allocated(icir)) deallocate(icir)
    n = this%nnod_dat
    if ((n == 0).or.(.not.allocated(this%nod_dat))) return
    allocate(icir(2,n))
    !
    do i = 1, n
      call node_to_icrl(this%nod_dat(i), ic, ir, il, nc, nr)
      if ((il /= 1).or.(ic < 1).or.(ic > nc).or.(ir < 1).or.(ir > nr)) then
        call errmsg('tGrid_set_nod_dat: ic/ir/il out of range.')
      else
        icir(1,i) = ic; icir(2,i) = ir
      end if
    end do
    !
    return
  end subroutine tGrid_get_nod_dat
    
  subroutine tGrid_init(this, nc, nr, bbx, mvi1, mvi2, mvi4, mvi8, mvr4, mvr8)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I4B), intent(in), optional :: nc
    integer(I4B), intent(in), optional :: nr
    type(tBBx),   intent(in), optional :: bbx
    integer(I1B), intent(in), optional :: mvi1
    integer(I2B), intent(in), optional :: mvi2
    integer(I4B), intent(in), optional :: mvi4
    integer(I8B), intent(in), optional :: mvi8
    real(R4B),    intent(in), optional :: mvr4
    real(R8B),    intent(in), optional :: mvr8
    ! -- local
! ------------------------------------------------------------------------------
    call this%clean()
    !
    this%i_data_type = 0
    this%nnod_dat    = 0
    this%nc          = 0
    this%nr          = 0
    call this%bbx%init()
    !
    call this%set_mv(mvi1, mvi2, mvi4, mvi8, mvr4, mvr8)
    if (present(nc)) then
      if (.not.present(nr)) call errmsg('tGrid_init: nr should be set.')
      this%nc = nc
    end if
    if (present(nr)) then
      if (.not.present(nc)) call errmsg('tGrid_init: nr should be set.')
      this%nr = nr
    end if
    if (present(bbx)) then
      this%bbx = bbx
      this%nc = (bbx%xur - bbx%xll)/bbx%cs
      this%nr = (bbx%yur - bbx%yll)/bbx%cs
    end if
    !
    select case(this%i_data_type)
    case(i_i1)
      allocate(this%xi1(this%nc,this%nr)); this%xi1 = mvi1
    case(i_i2)
      allocate(this%xi2(this%nc,this%nr)); this%xi2 = mvi2
    case(i_i4)
      allocate(this%xi4(this%nc,this%nr)); this%xi4 = mvi4
    case(i_i8)
      allocate(this%xi8(this%nc,this%nr)); this%xi8 = mvi8
    case(i_r4)
      allocate(this%xr4(this%nc,this%nr)); this%xr4 = mvr4
    case(i_r8)
      allocate(this%xr8(this%nc,this%nr)); this%xr8 = mvr8
    end select
    !
    return
  end subroutine tGrid_init
  !
  subroutine tGrid_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    ! -- local
! ------------------------------------------------------------------------------
    call this%clean_nod_dat()
    call this%clean_xi()
    !
    if (allocated(this%mask)) deallocate(this%mask)
    this%nmask = 0
    !
    return
    end subroutine tGrid_clean
    
  subroutine tGrid_clean_xi(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (allocated(this%xi1)) deallocate(this%xi1)
    if (allocated(this%xi2)) deallocate(this%xi2)
    if (allocated(this%xi4)) deallocate(this%xi4)
    if (allocated(this%xi8)) deallocate(this%xi8)
    if (allocated(this%xr4)) deallocate(this%xr4)
    if (allocated(this%xr8)) deallocate(this%xr8)
    !
    if (allocated(this%mvi1)) deallocate(this%mvi1)
    if (allocated(this%mvi2)) deallocate(this%mvi2)
    if (allocated(this%mvi4)) deallocate(this%mvi4)
    if (allocated(this%mvi8)) deallocate(this%mvi8)
    if (allocated(this%mvr4)) deallocate(this%mvr4)
    if (allocated(this%mvr8)) deallocate(this%mvr8)
    !
    return
  end subroutine tGrid_clean_xi
  
  subroutine tGrid_clean_nod_dat(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (allocated(this%nod_dat)) deallocate(this%nod_dat)
    !
    this%nnod_dat = 0
    !
    return
  end subroutine tGrid_clean_nod_dat
    
  function tGrid_any_pos(this) result(found)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    logical :: found
    ! -- local
    integer(I4B) :: ir, ic
! ------------------------------------------------------------------------------
    found = .false.
    if (allocated(this%xi4)) then
      do ir =  1, this%nr
        do ic = 1, this%nc
          if (this%xi4(ic,ir) > 0) then
            found = .true.; exit
          end if
        end do
      end do
    end if
    !
  end function tGrid_any_pos
    
  function tGrid_any_neg(this) result(found)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    logical :: found
    ! -- local
    integer(I4B) :: ir, ic
! ------------------------------------------------------------------------------
    found = .false.
    if (allocated(this%xi4)) then
      do ir =  1, this%nr
        do ic = 1, this%nc
          if (this%xi4(ic,ir) < 0) then
            found = .true.; exit
          end if
        end do
      end do
    end if
    !
  end function tGrid_any_neg
    
  function tGrid_count_data(this) result(n)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGrid) :: this
    integer(I4B) :: n
    ! -- local
    integer(I4B) :: ir, ic
! ------------------------------------------------------------------------------
    if (this%i_data_type == 0) then
      call errmsg('tGrid_data_set: data_type not set.')
    end if
    !
    n = 0
    select case(this%i_data_type)
    case(i_i1)
      call errmsg('tGrid_data_set: i1 not yet supported.')
    case(i_i2)
      call errmsg('tGrid_data_set: i4 not yet supported.')
    case(i_i4)
      if ((.not.allocated(this%xi4)).or.(.not.allocated(this%mvi4))) then
        call errmsg('tGrid_data_set: xi4/mvi4 not found.')
      end if
      do ir = 1, this%nr; do ic = 1, this%nc
        if (this%xi4(ic,ir) /= this%mvi4) then
          n = n + 1
        end if
      end do; end do
    case(i_r4)
      if ((.not.allocated(this%xr4)).or.(.not.allocated(this%mvr4))) then
        call errmsg('tGrid_data_set: xr4/mvr4 not found.')
      end if
      do ir = 1, this%nr; do ic = 1, this%nc
        if (this%xr4(ic,ir) /= this%mvr4) then
          n = n + 1
        end if
      end do; end do
    case(i_r8)
      call errmsg('tGrid_data_set: r8 not yet supported.')
    end select
    !
  end function tGrid_count_data
  
! ==============================================================================
! ==============================================================================
! tGridArr
! ==============================================================================
! ==============================================================================
    
  subroutine tGridArr_init(this, nlev)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGridArr) :: this
    integer(I4B), intent(in), optional :: nlev
    ! -- local
! ------------------------------------------------------------------------------
    !
    call this%clean()
    if (present(nlev)) then
      this%n = nlev
    else
      this%n = 1
    end if
    allocate(this%x(this%n))
    !
    return
  end subroutine tGridArr_init
    
  subroutine tGridArr_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGridArr) :: this
    ! -- local
    type(tGrid), pointer :: g => null()
    integer(I4B) :: i
    
! ------------------------------------------------------------------------------
    !
    if (associated(this%x)) then
      do i = 1, this%n
        g => this%x(i)
        call g%clean()
      end do
      deallocate(this%x)
    end if
    this%n = 0
    !
    return
  end subroutine tGridArr_clean
    
! ==============================================================================
! ==============================================================================
! tCSV_hdr
! ==============================================================================
! ==============================================================================
  subroutine tCSV_hdr_init(this)
! ******************************************************************************  
    ! -- arguments
    class(tCSV_hdr) :: this
    ! -- locals
! ------------------------------------------------------------------------------
    this%key = ''
    this%i_type = I4ZERO
    !
    return
  end subroutine tCSV_hdr_init
  
  subroutine tCSV_hdr_set(this, key, i_type)
! ******************************************************************************  
    ! -- arguments
    class(tCSV_hdr) :: this
    character(len=*), intent(in), optional :: key
    integer(I4B), intent(in), optional :: i_type
    ! -- locals
! ------------------------------------------------------------------------------
    !
    if (present(key)) then
      this%key = key
      this%key = change_case(this%key, 'l') ! convert to lower case
    end if
    if (present(i_type)) then
      this%i_type = i_type
    end if
    !
    return
  end subroutine tCSV_hdr_set
    
  subroutine tCSV_hdr_get(this, key, i_type)
! ******************************************************************************  
    ! -- arguments
    class(tCSV_hdr) :: this
    character(len=*), intent(out), optional :: key
    integer(I4B), intent(out), optional :: i_type
    ! -- locals
! ------------------------------------------------------------------------------
    !
    if (present(key)) then
      key = this%key 
    end if
    if (present(i_type)) then
      i_type = this%i_type
    end if
    !
    return
  end subroutine tCSV_hdr_get
  
  subroutine tCSV_hdr_copy(this, tgt)
! ******************************************************************************  
    ! -- arguments
    class(tCSV_hdr) :: this
    type(tCSV_hdr), intent(inout), pointer :: tgt
    ! -- locals
! ------------------------------------------------------------------------------
    tgt%key = this%key
    tgt%i_type = tgt%i_type
    !
    return
  end subroutine tCSV_hdr_copy
  
! ==============================================================================
! ==============================================================================
! tVal
! ==============================================================================
! ==============================================================================
  
  subroutine tVal_init(this)
! ******************************************************************************  
    ! -- arguments
    class(tVal) :: this
! ------------------------------------------------------------------------------
    !
    this%lnum = .true.
    allocate(this%x)
    !
    return
  end subroutine tVal_init

  subroutine tVal_clean(this)
! ******************************************************************************  
    ! -- arguments
    class(tVal) :: this
! ------------------------------------------------------------------------------
    !
    if (associated(this%x)) deallocate(this%x)
    if (allocated(this%s)) deallocate(this%s)
    !
    this%x => null()
    !
    return
  end subroutine tVal_clean
  
  subroutine tVal_set_val(this, i1v, i2v, i4v, i8v, r4v, r8v, cv, create_s)
! ******************************************************************************  
    ! -- arguments
    class(tVal) :: this
    !
    integer(I1B),     intent(in), optional :: i1v
    integer(I2B),     intent(in), optional :: i2v
    integer(I4B),     intent(in), optional :: i4v
    integer(I8B),     intent(in), optional :: i8v
    real(R4B),        intent(in), optional :: r4v
    real(R8B),        intent(in), optional :: r8v
    character(len=*), intent(in), optional :: cv
    logical,          intent(in), optional :: create_s
    ! -- locals
    type(tNum), pointer :: num => null()
    logical :: create
    character(len=MXSLEN) :: s
    integer(I4B) :: n, ns
! ------------------------------------------------------------------------------
    !
    if (present(create_s)) then
      create = create_s
    else
      create = .false.
    end if
    !
    num => this%x
    !
    n = 0
    if (present(i1v)) then
      call num%set_val(i1v=i1v); n = n + 1
      if (create) write(s,*) i1v
    end if
    if (present(i2v)) then
      call num%set_val(i2v=i2v); n = n + 1
      if (create) write(s,*) i2v
    end if
    if (present(i4v)) then
      call num%set_val(i4v=i4v); n = n + 1
      if (create) write(s,*) i4v
    end if
    if (present(i8v)) then
      call num%set_val(i8v=i8v); n = n + 1
      if (create) write(s,*) i8v
    end if
    if (present(r4v)) then
      call num%set_val(r4v=r4v); n = n + 1
      if (create) write(s,*) r4v
    end if
    if (present(r8v)) then
      call num%set_val(r8v=r8v); n = n + 1
      if (create) write(s,*) r8v
    end if
    if (present(cv)) then
      this%lnum = .false.
      if (allocated(this%s)) deallocate(this%s)
      allocate(this%s, source=trim(cv))
    else
      if (create) then
        s = adjustl(s)
        ns = len_trim(s)
        if (allocated(this%s)) deallocate(this%s)
        allocate(character(ns) :: this%s)
        this%s = trim(s)
      end if
      this%lnum = .true.
    end if
    !
    if (n > 1) then
      call errmsg('tVal_set_val: too many values set.')
    end if
    !
    return
  end subroutine tVal_set_val
  
  subroutine tVal_set_val_s(this, s)
! ******************************************************************************  
    ! -- arguments
    class(tVal) :: this
    character(len=*), intent(in) :: s
    ! -- locals
    type(tNum), pointer :: num => null()
! ------------------------------------------------------------------------------
    num => this%x
    call num%set_val_by_s(s)
    allocate(this%s, source=trim(s))
    !
    return
  end subroutine tVal_set_val_s
  
  subroutine tVal_get_val(this, i1v, i2v, i4v, i8v, r4v, r8v, l4v, cv)
! ******************************************************************************  
    ! -- arguments
    class(tVal) :: this
    !
    integer(I1B),     intent(out), optional :: i1v
    integer(I2B),     intent(out), optional :: i2v
    integer(I4B),     intent(out), optional :: i4v
    integer(I8B),     intent(out), optional :: i8v
    real(R4B),        intent(out), optional :: r4v
    real(R8B),        intent(out), optional :: r8v
    logical,          intent(out), optional :: l4v
    character(len=*), intent(out), optional :: cv
    ! -- locals
    type(tNum), pointer :: num => null()
    character(len=MXSLEN) :: s
    integer(I4B) :: n
    logical :: lset
! ------------------------------------------------------------------------------
    !
    num => this%x
    !
    n = 0
    if (present(i1v)) then
      call num%get_val(i1v=i1v, lset=lset); n = n + 1
    end if
    if (present(i2v)) then
      call num%get_val(i2v=i2v, lset=lset); n = n + 1
    end if
    if (present(i4v)) then
      call num%get_val(i4v=i4v, lset=lset); n = n + 1
    end if
    if (present(i8v)) then
      call num%get_val(i8v=i8v, lset=lset); n = n + 1
    end if
    if (present(r4v)) then
      call num%get_val(r4v=r4v, lset=lset); n = n + 1
    end if
    if (present(r8v)) then
      call num%get_val(r8v=r8v, lset=lset); n = n + 1
    end if
    if (present(l4v)) then
      call num%get_val(l4v=l4v, lset=lset); n = n + 1
    end if
    if (present(cv)) then
      if (.not.allocated(this%s)) then
        call logmsg('tVal_get_val: character data not present.')
      else
        cv = this%s
      end if
    end if
    !
    if (n > 1) then
      call logmsg('tVal_get_val: get multiple values.')
    end if
    !
    return
  end subroutine tVal_get_val
  
  subroutine tVal_get_val_arr(this, ir, i1a, i2a, i4a, i8a, r4a, r8a, ca)
! ******************************************************************************  
    ! -- arguments
    class(tVal) :: this
    !
    integer(I4B), intent(in) :: ir
    integer(I1B),     dimension(:), allocatable, intent(inout), optional :: i1a
    integer(I2B),     dimension(:), allocatable, intent(inout), optional :: i2a
    integer(I4B),     dimension(:), allocatable, intent(inout), optional :: i4a
    integer(I8B),     dimension(:), allocatable, intent(inout), optional :: i8a
    real(R4B),        dimension(:), allocatable, intent(inout), optional :: r4a
    real(R8B),        dimension(:), allocatable, intent(inout), optional :: r8a
    character(len=*), dimension(:), allocatable, intent(inout), optional :: ca
    ! -- locals
    type(tNum), pointer :: num => null()
    character(len=MXSLEN) :: s
    integer(I4B) :: n
    logical :: lset
! ------------------------------------------------------------------------------
    !
    num => this%x
    !
    n = 0
    if (present(i1a)) then
      call num%get_val(i1v=i1a(ir), lset=lset); n = n + 1
    end if
    if (present(i2a)) then
      call num%get_val(i2v=i2a(ir), lset=lset); n = n + 1
    end if
    if (present(i4a)) then
      call num%get_val(i4v=i4a(ir), lset=lset); n = n + 1
    end if
    if (present(i8a)) then
      call num%get_val(i8v=i8a(ir), lset=lset); n = n + 1
    end if
    if (present(r4a)) then
      call num%get_val(r4v=r4a(ir), lset=lset); n = n + 1
    end if
    if (present(r8a)) then
      call num%get_val(r8v=r8a(ir), lset=lset); n = n + 1
    end if
    if (present(ca)) then
      if (.not.allocated(this%s)) then
        call logmsg('tVal_get_val_arr: character data not present.')
      else
        ca(ir) = this%s
      end if
    end if
    !
    if (n > 1) then
      call logmsg('tVal_get_val_arr: get multiple values.')
    end if
    !
    return
  end subroutine tVal_get_val_arr
  
  function tVal_get_val_s(this, i_type, &
    i1mv, i2mv, i4mv, i8mv, r4mv, r8mv, cmv) result(s)
! ******************************************************************************  
    ! -- arguments
    class(tVal) :: this
    character(len=MXSLEN) :: s
    !
    integer(I4B),     intent(in) :: i_type
    integer(I1B),     intent(in), optional :: i1mv
    integer(I2B),     intent(in), optional :: i2mv
    integer(I4B),     intent(in), optional :: i4mv
    integer(I8B),     intent(in), optional :: i8mv
    real(R4B),        intent(in), optional :: r4mv
    real(R8B),        intent(in), optional :: r8mv
    character(len=*), intent(in), optional :: cmv
    ! -- locals
    type(tNum), pointer :: num => null()
    !
    logical :: apply_mv
    integer(I4B) :: n_mv
    integer(I1B)          :: i1v
    integer(I2B)          :: i2v
    integer(I4B)          :: i4v
    integer(I8B)          :: i8v
    real(R4B)             :: r4v
    real(R8B)             :: r8v
    character(len=MXSLEN) :: cv
! ------------------------------------------------------------------------------
    !
    ! set the missing values
    n_mv = 0
    if (present(i1mv)) then
      n_mv = n_mv + 1
      call cast_number(i1v, i2v, i4v, i8v, r4v, r8v, cv, i1v_in=i1mv)
    end if
    if (present(i2mv)) then
      n_mv = n_mv + 1
      call cast_number(i1v, i2v, i4v, i8v, r4v, r8v, cv, i2v_in=i2mv)
    end if
    if (present(i4mv)) then
      n_mv = n_mv + 1
      call cast_number(i1v, i2v, i4v, i8v, r4v, r8v, cv, i4v_in=i4mv)
    end if
    if (present(i8mv)) then
      n_mv = n_mv + 1
      call cast_number(i1v, i2v, i4v, i8v, r4v, r8v, cv, i8v_in=i8mv)
    end if
    if (present(r4mv)) then
      n_mv = n_mv + 1
      call cast_number(i1v, i2v, i4v, i8v, r4v, r8v, cv, r4v_in=r4mv)
    end if
    if (present(r8mv)) then
      n_mv = n_mv + 1
      call cast_number(i1v, i2v, i4v, i8v, r4v, r8v, cv, r8v_in=r8mv)
    end if
    if (present(cmv)) then
      n_mv = n_mv + 1
      call cast_number(i1v, i2v, i4v, i8v, r4v, r8v, cv, cv_in=cmv)
    end if
    if (n_mv > 1) then
      call errmsg('tVal_get_val_s: too many missing values specified.')
    else if (n_mv == 1) then
      apply_mv = .true.
    else
      apply_mv = .false.
    end if
    !
    s = ''
    ! 
    if (i_type == i_c) then
      if (.not.apply_mv) then
        if(.not.allocated(this%s)) then
          call errmsg('tVal_get_val_s: could not retrieve value.')
        else
          s = this%s
        end if
      end if
      return
    end if
    
    num => this%x
    !
    if (.not.num%flg(i_type)) then
      call errmsg('tVal_get_val_s: could not retrieve value.')
    end if
    !
    select case(i_type)
    case(i_i1)
      if (apply_mv) then
        if (num%i1v /= i1v) write(s,*) num%i1v
      else
        write(s,*) num%i1v
      end if
    case(i_i2)
      if (apply_mv) then
        if (num%i2v /= i2v) write(s,*) num%i2v
      else
        write(s,*) num%i2v
      end if
    case(i_i4)
      if (apply_mv) then
        if (num%i4v /= i4v) write(s,*) num%i4v
      else
        write(s,*) num%i4v
      end if
    case(i_i8)
      if (apply_mv) then
        if (num%i8v /= i8v) write(s,*) num%i8v
      else
        write(s,*) num%i8v
      end if
    case(i_r4)
      if (apply_mv) then
        if (num%r4v /= r4v) write(s,*) num%r4v
      else
        write(s,*) num%r4v
      end if
    case(i_r8)
      if (apply_mv) then
        if (num%r8v /= r8v) write(s,*) num%r8v
      else
        write(s,*) num%r8v
      end if
    end select
    !
    if (this%lnum) then
      s = adjustl(s); s = change_case(s,'l')
      s = remove_trailing_zeros(s)
    end if
    !
    return
  end function tVal_get_val_s
 
  subroutine tVal_copy(this, tgt)
! ******************************************************************************  
    ! -- arguments
    class(tVal) :: this
    type(tVal), pointer, intent(inout) :: tgt
    ! -- locals
    type(tNum), pointer :: x, x_tgt
! ------------------------------------------------------------------------------
    !
    call tgt%clean()
    !
    tgt%lnum = this%lnum
    allocate(tgt%x); x_tgt => tgt%x; x => this%x
    call x%copy(x_tgt)
    allocate(tgt%s, source=this%s)
    !
    return
  end subroutine tVal_copy
    
! ==============================================================================
! ==============================================================================
! tIni
! ==============================================================================
! ==============================================================================

  subroutine tIni_init(this, f)
! ******************************************************************************  
    ! -- arguments
    class(tIni) :: this
    character(len=*), intent(in) :: f
    ! -- locals
! ------------------------------------------------------------------------------
    call this%clean()
    !
    this%f = f
    call this%read()
    !
    return
  end subroutine tIni_init
  
  subroutine tIni_clean(this)
! ******************************************************************************  
    ! -- arguments
    class(tIni) :: this
    ! -- locals
    type(tIniSect), pointer :: sect => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    if (allocated(this%sect_names)) then
      deallocate(this%sect_names)
    end if
    !
    if (associated(this%sect)) then
      do i = 1, this%n_sect
        sect => this%sect(i)
        call sect%clean()
      end do
      deallocate(this%sect); this%sect => null()
    end if
    !
    this%f = ''
    this%n_sect = 0
    !
    return
  end subroutine tIni_clean
    
  subroutine tIni_read(this)
! ******************************************************************************  
    ! -- arguments
    class(tIni) :: this
    ! -- locals
    integer(I4B), parameter :: n_raw_max = 10000
    !
    type(tIniSect), pointer :: sect => null()
    character(len=MXSLEN), dimension(:), allocatable :: raw, sa
    character(len=MXSLEN) :: s
    integer(I4B), dimension(1) :: loc
    integer(I4B), dimension(:), allocatable :: sect_idx
    integer(I4B) :: iu, i_raw, n_raw, i_sect, i0, i1, n, ios
! ------------------------------------------------------------------------------
    !
    allocate(raw(n_raw_max))
    !
    call open_file(this%f, iu, 'r')
    !
    n_raw = 0
    do while(.true.)
      n_raw = n_raw + 1
      if (n_raw > n_raw_max) call errmsg('tIni_read: please increase n_raw_max.')
      call read_line(iu, s, ios, [';','#'])
      if (len_trim(s) == 0) then
        n_raw = n_raw - 1
        exit
      end if
      raw(n_raw) = s
    end do
    close(iu)
    !
    if (n_raw == 0) then
      call logmsg('tIni_read: no data found, returning.')
      return
    end if
    !
    allocate(sect_idx(n_raw)); sect_idx = 0
    do i_raw = 1, n_raw
       s = raw(i_raw)
       if (len_trim(get_section_name(s)) > 0) then
          this%n_sect = this%n_sect + 1
       end if
       sect_idx(i_raw) = this%n_sect
    end do
    !
    if (this%n_sect == 0) then
      call logmsg('tIni_read: no sections found, returning.')
      return
    end if
    !
    allocate(this%sect(this%n_sect))
    allocate(this%sect_names(this%n_sect))
    !
    do i_sect = 1, this%n_sect
      loc = findloc(sect_idx, i_sect); i0 = loc(1)
      if (i_sect < this%n_sect) then
        loc = findloc(sect_idx, i_sect+1); i1 = loc(1) - 1
      else
        loc = findloc(sect_idx, i_sect, back=.true.); i1 = loc(1)
      end if
      !
      sect => this%sect(i_sect)
      this%sect_names(i_sect) = get_section_name(raw(i0))
      call sect%init(get_section_name(raw(i0)), i1 - i0)
      n = 0
      do i_raw = i0 + 1, i1
        n = n + 1
        call split_str(raw(i_raw), '=', sa)
        call sect%set(n, sa(1), sa(2))
      end do
    end do
    !
    ! obselete:
    if (allocated(raw)) deallocate(raw)
    if (allocated(sect_idx)) deallocate(sect_idx)
    !
    return
  end subroutine tIni_read
    
  subroutine tIni_get_val(this, sect_name, key, &
    i1v, i2v, i4v, i8v, r4v, r8v, l4v, cv, &
    i1v_def, i2v_def, i4v_def, i8v_def, r4v_def, r8v_def, l4v_def, cv_def)
! ******************************************************************************
    ! -- arguments
    class(tIni) :: this
    ! -- locals
    character(len=*), intent(in) :: sect_name
    character(len=*), intent(in) :: key
    !
    integer(I1B), intent(out), optional :: i1v
    integer(I2B), intent(out), optional :: i2v
    integer(I4B), intent(out), optional :: i4v
    integer(I8B), intent(out), optional :: i8v
    real(R4B),    intent(out), optional :: r4v
    real(R8B),    intent(out), optional :: r8v
    logical,      intent(out), optional :: l4v
    character(len=MXSLEN), intent(out), optional :: cv
    !
    integer(I1B), intent(in), optional :: i1v_def
    integer(I2B), intent(in), optional :: i2v_def
    integer(I4B), intent(in), optional :: i4v_def
    integer(I8B), intent(in), optional :: i8v_def
    real(R4B),    intent(in), optional :: r4v_def
    real(R8B),    intent(in), optional :: r8v_def
    logical,      intent(in), optional :: l4v_def
    character(len=*), intent(in), optional :: cv_def
    !
    ! -- locals
    type(tIniSect), pointer :: sect => null()
    type(tVal), pointer :: v => null()
    logical :: found, use_def
    character(len=MXSLEN) :: s
    integer(I4B), dimension(1) :: loc
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    found = .true.
    loc = findloc(this%sect_names, trim(sect_name)); i = loc(1)
    if (i <= 0) then
      found = .false.
    end if
    !
    sect => this%sect(i)
    loc = findloc(sect%keys, trim(key)); i = loc(1)
    if (i <= 0) then
      found = .false.
    end if
    !
    if (.not.found) then
      use_def = .false.
      if (present(i1v).and.present(i1v_def)) then
        use_def = .true.; i1v = i1v_def; s = ta([i1v])
      end if
      if (present(i2v).and.present(i2v_def)) then
        use_def = .true.; i2v = i2v_def; s = ta([i2v])
      end if
      if (present(i4v).and.present(i4v_def)) then
        use_def = .true.; i4v = i4v_def; s = ta([i4v])
      end if
      if (present(i8v).and.present(i8v_def)) then
        use_def = .true.; i8v = i8v_def; s = ta([i8v])
      end if
      if (present(r4v).and.present(r4v_def)) then
        use_def = .true.; r4v = r4v_def; s = ta([r4v])
      end if
      if (present(r8v).and.present(r8v_def)) then
        use_def = .true.; r8v = r8v_def; s = ta([r8v])
      end if
      if (present(l4v).and.present(l4v_def)) then
        use_def = .true.; l4v = l4v_def; write(s,*) l4v
      end if
      if (present(cv).and.present(cv_def)) then
        use_def = .true.; cv = trim(cv_def); s = cv
      end if
      if (.not. use_def) then
        call errmsg('tIni_get_val: could set value for key '//trim(key)//&
          ' in ['//trim(sect%name)//'].')
      else
        call logmsg('Using default value for key '//trim(key)//&
          ' in ['//trim(sect%name)//']: '//trim(adjustl(s)))
      end if
    else
      v => sect%v(i)
      call v%get_val(i1v, i2v, i4v, i8v, r4v, r8v, l4v, cv)
    end if
    !
    return
  end subroutine tIni_get_val
  
  function get_section_name(s) result(name)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: s
    character(len=MXSLEN) :: name
    ! -- locals
    integer(I4B) :: i0, i1
! ------------------------------------------------------------------------------
    i0 = index(s,'['); i1 = index(s,']')
    if ((i0 > 0).and.(i1 > 0)) then
      name = trim(adjustl(s(i0+1:i1-1)))
    else
      name = ''
    end if
    !
    return
  end function get_section_name
  
  subroutine tIniSect_init(this, name, n)
! ******************************************************************************  
    ! -- arguments
    class(tIniSect) :: this
    character(len=*), intent(in) :: name
    integer(I4B), intent(in) :: n
    ! -- locals
! ------------------------------------------------------------------------------
    call this%clean()
    !
    this%name = name
    this%n = n
    if (this%n > 0) then
      allocate(this%keys(this%n))
      allocate(this%v(this%n))
    end if
    !
    return
  end subroutine tIniSect_init
  
  subroutine tIniSect_clean(this)
! ******************************************************************************  
    ! -- arguments
    class(tIniSect) :: this
    ! -- locals
    type(tVal), pointer :: v => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    this%name = ''
    this%n = 0
    if (allocated(this%keys)) then
      deallocate(this%keys)
    end if
    !
    if (associated(this%v)) then
      do i = 1, size(this%v)
        v => this%v(i)
        call v%clean()
      end do
      deallocate(this%v); this%v => null()
    end if
    !
    return
  end subroutine tIniSect_clean
  
  subroutine tIniSect_set(this, i, key, cv)
! ******************************************************************************  
    ! -- arguments
    class(tIniSect) :: this
    !
    integer(I4B), intent(in) :: i
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: cv
    ! -- locals
    type(tVal), pointer :: v => null()
! ------------------------------------------------------------------------------
    this%keys(i) = key
    !
    v => this%v(i)
    call v%init()
    call v%set_val_s(cv)
    !
    return
  end subroutine tIniSect_set
  
! ==============================================================================
! ==============================================================================
! tCSV
! ==============================================================================
! ==============================================================================
    
  subroutine tCSV_init(this, file, hdr_keys, nr, nc, nr_max, nc_max, hdr_i_type)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*), intent(in) :: file
    character(len=*), dimension(:), intent(in) :: hdr_keys
    integer(I4B), intent(in), optional :: nr
    integer(I4B), intent(in), optional :: nc
    integer(I4B), intent(in), optional :: nr_max
    integer(I4B), intent(in), optional :: nc_max
    integer(I4B), dimension(:), intent(in), optional :: hdr_i_type
    ! -- locals
    type(tVal), pointer :: v => null()
    integer(I4B) :: ic, ir
! ------------------------------------------------------------------------------
    !
    this%file = trim(file)
    this%nc     = 0
    this%nr     = 0
    this%nc_max = 0
    this%nr_max = 0
    !
    if (present(nr)) then
      this%nr = nr
    end if
    if (present(nr_max)) then
      this%nr_max = nr_max
    else
      this%nr_max = this%nr
    end if
    !
    if (present(nc)) then
      this%nc = nc
    else
      this%nc = size(hdr_keys)
    end if
    if (present(nc_max)) then
      this%nc_max = nc_max
    else
      this%nc_max = this%nc
    end if
    !
    allocate(this%hdr(this%nc_max))
    !
    if (present(hdr_i_type)) then
      call this%set_hdr(keys=hdr_keys, i_type=hdr_i_type)
    else
      call this%set_hdr(keys=hdr_keys)
    end if
    !
    !
    if ((this%nr_max > 0).and.(this%nc_max > 0)) then
      allocate(this%val(this%nc_max, this%nr_max))
      do ir = 1, this%nr_max
        do ic = 1, this%nc_max
          v => this%val(ic,ir)
          call v%init()
        end do
      end do
    end if
    !
    return
  end subroutine tCSV_init

  subroutine tCSV_clean(this)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    ! -- locals
    type(tCSV_hdr), pointer :: hdr => null()
    type(tVal), pointer :: v => null()
    integer(I4B) :: i, ir, ic
! ------------------------------------------------------------------------------
    !
    if (associated(this%hdr)) then
      deallocate(this%hdr); this%hdr => null()
    end if
    !
    if (associated(this%val)) then
      do ir = 1, size(this%val,2)
        do ic = 1, size(this%val,1)
          v => this%val(ic,ir)
          call v%clean()
        end do
      end do
      deallocate(this%val); this%val => null()
    end if
    !
    this%file = ''
    this%nc     = I4ZERO
    this%nr     = I4ZERO
    this%nr_max = I4ZERO
    !
    return
  end subroutine tCSV_clean
  
  subroutine tCSV_clean_and_init_row(this, ir)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    integer(I4B), intent(in) :: ir
    ! -- locals
    type(tVal), pointer :: v => null()
    integer(I4B) :: i, ic
! ------------------------------------------------------------------------------
    if (associated(this%val)) then
      do ic = 1, size(this%val,1)
        v => this%val(ic,ir)
        call v%clean()
        call v%init()
      end do
    end if
    !
    return
  end subroutine tCSV_clean_and_init_row
  
  subroutine tCSV_set_hdr(this, keys, i_type)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*), dimension(:), intent(in) :: keys
    integer(I4B), dimension(:), intent(in), optional :: i_type
    ! -- locals
    type(tCSV_hdr), pointer :: hdr => null()
    integer(I4B) :: i, i_type_set
! ------------------------------------------------------------------------------
    !
    do i = 1, size(keys)
      hdr => this%hdr(i)
      !
      if (present(i_type)) then
        i_type_set = i_type(i)
      else
        i_type_set = 0
      end if
      call hdr%set(keys(i), i_type_set)
    end do
    !
    return
  end subroutine tCSV_set_hdr
  
  subroutine tCSV_add_hdr(this, keys, i_type)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*), dimension(:), intent(in) :: keys
    integer(I4B), dimension(:), intent(in), optional :: i_type
    ! -- locals
    type(tCSV_hdr), pointer :: hdr => null()
    integer(I4B) :: i, j, i_type_set
! ------------------------------------------------------------------------------
    !
    do j = 1, size(keys)
      i = this%nc + j
      if (i > size(this%hdr)) then
        call errmsg('tCSV_add_hdr: column out of range.')
      end if
      !
      hdr => this%hdr(i)
      !
      if (present(i_type)) then
        i_type_set = i_type(j)
      else
        i_type_set = 0
      end if
      call hdr%set(keys(j), i_type_set)
    end do
    !
    this%nc = this%nc + size(keys)
    return
  end subroutine tCSV_add_hdr
  
  subroutine tCSV_read_hdr(this, file, hdr_keys)
    ! -- arguments
    class(tCSV) :: this
    character(len=*) :: file
    character(len=MXSLEN), dimension(:), allocatable, intent(inout) :: hdr_keys
    ! -- locals
    character(len=MXSLENLONG) :: s
    integer(I4B) :: iu
! ------------------------------------------------------------------------------
    !
    call open_file(file, iu, 'r')
    read(iu,'(a)') s; close(iu)
    if (allocated(hdr_keys)) deallocate(hdr_keys)
    call parse_line(s, hdr_keys, token_in=',')
    !
    return
  end subroutine tCSV_read_hdr
 
  subroutine tCSV_read(this, file, nr_max, nc_max, ir_offset)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*) :: file
    integer(I4B), intent(in), optional :: nr_max
    integer(I4B), intent(in), optional :: nc_max
    integer(I4B), intent(in), optional :: ir_offset
    ! -- locals
    type(tVal), pointer :: v => null()
    character(len=MXSLENLONG) :: s
    character(len=MXSLEN), dimension(:), allocatable :: hdr_keys, sa
    integer(I4B) :: iu, ios, nr, nc, ir, ic, nr_max_loc, nc_max_loc
! ------------------------------------------------------------------------------
    !
    call open_file(file, iu, 'r')
    read(iu,'(a)') s
    call parse_line(s, hdr_keys, token_in=',')
    nc = size(hdr_keys)
    !
    nr = 0
    do while(.true.)
      read(unit=iu,iostat=ios,fmt='(a)') s
      if (ios == 0) then
        if (len_trim(s) > 0) then
          nr = nr + 1
        end if
      else
        exit
      end if
    end do
    rewind(iu); read(iu,'(a)') s
    !
    if (present(nr_max)) then
      nr_max_loc = max(nr, nr_max)
    else
      nr_max_loc = nr
    end if
    if (present(nc_max)) then
      nc_max_loc = max(nc, nc_max)
    else
      nc_max_loc = nc
    end if
    !
    if (present(ir_offset)) then
      ir = ir_offset
    else
      call this%init(file=file, hdr_keys=hdr_keys, nr=nr, nc=nc, nr_max=nr_max_loc, nc_max=nc_max)
      ir = 0
    end if
    do while(.true.)
      read(unit=iu,iostat=ios,fmt='(a)') s
      if (ios == 0) then
        if (len_trim(s) > 0) then
          ir = ir + 1
          call parse_line(s, sa, token_in=',')
          if (size(sa) /= this%nc) then
            call errmsg('tCSV_read: could not read line.')
          end if
          do ic = 1, this%nc
            v => this%val(ic,ir)
            v%lnum = .false. !24-4-23
            call v%set_val_s(sa(ic))
          end do
        end if
      else
        exit
      end if
    end do
    !
    close(iu)
    !
    ! cleanup
    deallocate(hdr_keys, sa)
    !
    return
  end subroutine tCSV_read
  
  subroutine tCSV_write(this, i1mv, i2mv, i4mv, i8mv, r4mv, r8mv, cmv, ir0, ir1)
! *****************************************************************************
    ! -- arguments
    class(tCSV) :: this
    integer(I1B), intent(in), optional :: i1mv
    integer(I2B), intent(in), optional :: i2mv
    integer(I4B), intent(in), optional :: i4mv
    integer(I8B), intent(in), optional :: i8mv
    real(R4B),    intent(in), optional :: r4mv
    real(R8B),    intent(in), optional :: r8mv
    character(len=*), intent(in), optional :: cmv
    integer(I4B), intent(in), optional :: ir0
    integer(I4B), intent(in), optional :: ir1
    ! -- locals
    type(tCSV_hdr), pointer :: hdr => null()
    type(tVal), pointer :: v => null()
    character(len=MXSLENLONG) :: s
    character(len=MXSLEN), dimension(:), allocatable :: sa
    integer(I4B) :: iu, ir, ic, i_type, n_mv, ir0_loc, ir1_loc
! ------------------------------------------------------------------------------
    !
    ! check
    n_mv = 0
    if (present(i1mv)) n_mv = n_mv + 1
    if (present(i2mv)) n_mv = n_mv + 1
    if (present(i4mv)) n_mv = n_mv + 1
    if (present(i8mv)) n_mv = n_mv + 1
    if (present(r4mv)) n_mv = n_mv + 1
    if (present(r8mv)) n_mv = n_mv + 1
    if (present(cmv))  n_mv = n_mv + 1
    if (n_mv > 1) then
      call errmsg('tCSV_write: only 1 missing value allowed.')
    end if
    if (present(ir0)) then
      ir0_loc = ir0
    else
      ir0_loc = 1
    end if
    if (present(ir1)) then
      ir1_loc = ir1
    else
      ir1_loc = this%nr
    end if
    ir0_loc = max(1,ir0_loc); ir0_loc = min(ir0_loc,this%nr)
    ir1_loc = max(1,ir1_loc); ir1_loc = min(ir1_loc,this%nr)
    !
    call open_file(this%file, iu, 'w')
    !
    allocate(sa(this%nc))
    !
    ! write header
    do ic = 1, this%nc
      hdr => this%hdr(ic)
      call hdr%get(key=sa(ic))
    end do
    s = ta(sa, sep_in=',')
    write(iu,'(a)') trim(s)
    !
    ! write data
    do ir = ir0_loc, ir1_loc
      do ic = 1, this%nc
        v => this%val(ic,ir); hdr => this%hdr(ic)
        call hdr%get(i_type=i_type)
        if ((n_mv == 1).and.v%lnum) then
          if (present(i1mv)) sa(ic) = v%get_val_s(i_type, i1mv=i1mv)
          if (present(i2mv)) sa(ic) = v%get_val_s(i_type, i2mv=i2mv)
          if (present(i4mv)) sa(ic) = v%get_val_s(i_type, i4mv=i4mv)
          if (present(i8mv)) sa(ic) = v%get_val_s(i_type, i8mv=i8mv)
          if (present(r4mv)) sa(ic) = v%get_val_s(i_type, r4mv=r4mv)
          if (present(r8mv)) sa(ic) = v%get_val_s(i_type, r8mv=r8mv)
          if (present(cmv))  sa(ic) = v%get_val_s(i_type, cmv=cmv)
        else
          sa(ic) = v%get_val_s(i_type)
        end if
      end do
      s = ta(sa, sep_in=',')
      write(iu,'(a)') trim(s)
    end do
    !
    close(iu)
    !
    ! cleanup
    deallocate(sa)
    !
    return
  end subroutine tCSV_write
  
  function tCSV_get_nc(this) result(nc)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    integer(I4B) :: nc
 ! ------------------------------------------------------------------------------
    nc = this%nc
    !
    return
  end function tCSV_get_nc
   
  function tCSV_get_nr(this) result(nr)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    integer(I4B) :: nr
 ! ------------------------------------------------------------------------------
    nr = this%nr
    !
    return
  end function tCSV_get_nr
  
  function tCSV_get_col(this, hdr_key) result(ic)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*), intent(in) :: hdr_key
    integer(I4B) :: ic ! result
    ! -- locals
    character(len=MXSLEN) :: key
    logical :: found
    type(tCSV_hdr), pointer :: hdr => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    found = .false.; ic = 0
    do i = 1, this%nc
      hdr => this%hdr(i)
      call hdr%get(key=key)
      if (trim(key) == trim(hdr_key)) then
        found = .true.
        ic = i
        exit
      end if
    end do
    !
    if (.not.found) then
      call errmsg('tCSV_get_col: '//trim(hdr_key)//' not found.')
    end if
    !
    return
  end function tCSV_get_col
  
  function tCSV_get_row(this, key_val, hdr_key) result(ir)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*), intent(in) :: key_val
    character(len=*), intent(in) :: hdr_key
    integer(I4B) :: ir ! result
    ! -- locals
    character(len=MXSLEN) :: cv
    integer(I4B) :: ic, jr
! ------------------------------------------------------------------------------
    ir = 0
    ic = this%get_col(hdr_key)
    do jr = 1, this%nr
      call this%get_val(ir=jr, ic=ic, cv=cv)
      if (trim(key_val) == trim(cv)) then
        ir = jr; exit
      end if
    end do
  
    return
  end function tCSV_get_row
  
  function tCSV_exist_col(this, hdr_key) result(exist)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*), intent(in) :: hdr_key
    logical :: exist ! result
    ! -- locals
    character(len=MXSLEN) :: key
    logical :: found
    type(tCSV_hdr), pointer :: hdr => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    exist = .false.
    do i = 1, this%nc
      hdr => this%hdr(i)
      call hdr%get(key=key)
      if (trim(key) == trim(hdr_key)) then
        exist = .true.
        exit
      end if
    end do
    !
    return
  end function tCSV_exist_col
  
  function tCSV_get_key(this, ic) result(hdr_key)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    integer(I4B), intent(in) :: ic
    character(len=MXSLEN) :: hdr_key
    ! -- locals
    type(tCSV_hdr), pointer :: hdr => null()
! ------------------------------------------------------------------------------
    hdr => this%hdr(ic)
    call hdr%get(key=hdr_key)
    !
    return
  end function tCSV_get_key
  
  subroutine tCSV_add_key(this, key)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*), intent(in) :: key
    ! -- locals
    type(tCSV_hdr), pointer :: hdr_src => null(), hdr_tgt => null()
    type(tVal),     pointer :: val_src => null(), val_tgt => null()
    !
    type(tCSV_hdr), dimension(:),   pointer :: hdra => null()
    type(tVal),     dimension(:,:), pointer :: vala => null()
    !
    integer(I4B) :: ir, ic, nr, nc
! ------------------------------------------------------------------------------
    !
    nr = this%nr; nc = this%nc
    !
    allocate(hdra(nc+1))
    !
    ! copy the header and initialize
    do ic = 1, nc
      hdr_src => this%hdr(ic); hdr_tgt => hdra(ic)
      call hdr_src%copy(hdr_tgt)
    end do
    hdr_tgt => hdra(nc+1)
    call hdr_tgt%set(key)
    !
    ! set the new header
    deallocate(this%hdr)
    allocate(this%hdr(nc+1))
    do ic = 1, nc + 1
      hdr_src => hdra(ic); hdr_tgt => this%hdr(ic)
      call hdr_src%copy(hdr_tgt)
    end do
    !
    ! copy the values and initialize
    allocate(vala(nc+1,nr))
    !
    do ir = 1, nr; do ic = 1, nc
      val_src => this%val(ic,ir); val_tgt => vala(ic,ir)
      call val_src%copy(val_tgt)
      call val_src%clean()
    end do; end do
    ic = nc + 1
    do ir = 1, nr
      val_tgt => vala(ic,ir)
      call val_tgt%init()
    end do
    !
    ! set the new values
    deallocate(this%val)
    allocate(this%val(nc+1,nr))
    do ir = 1, nr; do ic = 1, nc + 1
      val_src => vala(ic,ir); val_tgt => this%val(ic,ir)
      call val_src%copy(val_tgt)
      call val_src%clean()
    end do; end do
    !
    ! set the number of columns
    this%nc = this%nc + 1
    !
    ! clean up
    deallocate(hdra)
    deallocate(vala)
    !
    return
  end subroutine tCSV_add_key
  
  subroutine tCSV_set_val(this, ic, ir, &
    i1v, i2v, i4v, i8v, r4v, r8v, cv, create_s)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    !
    integer(I4B),     intent(in)           :: ic
    integer(I4B),     intent(in)           :: ir
    integer(I1B),     intent(in), optional :: i1v
    integer(I2B),     intent(in), optional :: i2v
    integer(I4B),     intent(in), optional :: i4v
    integer(I8B),     intent(in), optional :: i8v
    real(R4B),        intent(in), optional :: r4v
    real(R8B),        intent(in), optional :: r8v
    character(len=*), intent(in), optional :: cv
    logical,          intent(in), optional :: create_s
    ! -- locals
    type(tVal), pointer :: v => null()
    type(tNum), pointer :: num => null()
    character(len=MXSLEN) :: s
    integer(I4B) :: n
! ------------------------------------------------------------------------------
    !
    v => this%val(ic,ir)
    !
    call v%set_val(i1v, i2v, i4v, i8v, r4v, r8v, cv, create_s)
    !
    return
  end subroutine tCSV_set_val
  
  subroutine tCSV_get_column(this, ic, key, i1a, i2a, i4a, i8a, r4a, r8a, ca)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    !
    integer(I4B),  intent(in), optional :: ic
    character(len=*), intent(in), optional :: key
    integer(I1B),     dimension(:), allocatable, intent(inout), optional :: i1a
    integer(I2B),     dimension(:), allocatable, intent(inout), optional :: i2a
    integer(I4B),     dimension(:), allocatable, intent(inout), optional :: i4a
    integer(I8B),     dimension(:), allocatable, intent(inout), optional :: i8a
    real(R4B),        dimension(:), allocatable, intent(inout), optional :: r4a
    real(R8B),        dimension(:), allocatable, intent(inout), optional :: r8a
    character(len=*), dimension(:), allocatable, intent(inout), optional :: ca
    !
    ! -- locals
    type(tVal), pointer :: v => null()
    character(len=MXSLEN) :: s
    integer(I4B) :: jc, nr, nc, ir
! ------------------------------------------------------------------------------
    !
    if (present(ic)) then
      jc = ic
    else
      jc = this%get_col(key)
    end if
    nr = this%get_nr(); nc = this%get_nc()
    !
    ! checks
    if ((jc < 1).or.(jc > nc)) then
      call errmsg('tCSV_get_column: invalid column number.')
    end if
    !
    if (present(i1a)) then
      if (allocated(i1a)) deallocate(i1a); allocate(i1a(nr))
    end if
    if (present(i2a)) then
      if (allocated(i2a)) deallocate(i2a); allocate(i2a(nr))
    end if
    if (present(i4a)) then
      if (allocated(i4a)) deallocate(i4a); allocate(i4a(nr))
    end if
    if (present(i8a)) then
      if (allocated(i8a)) deallocate(i8a); allocate(i8a(nr))
    end if
    if (present(r4a)) then
      if (allocated(r4a)) deallocate(r4a); allocate(r4a(nr))
    end if
    if (present(r8a)) then
      if (allocated(r8a)) deallocate(r8a); allocate(r8a(nr))
    end if
    if (present(ca)) then
      if (allocated(ca)) deallocate(ca); allocate(ca(nr))
    end if
    !
    do ir = 1, nr
      v => this%val(jc,ir)
      call v%get_val_arr(ir, i1a, i2a, i4a, i8a, r4a, r8a, ca)
    end do
    !
    return
  end subroutine tCSV_get_column

  subroutine tCSV_get_matrix(this, i1x, i2x, i4x, i8x, r4x, r8x, cx)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    !
    integer(I1B),     dimension(:,:), allocatable, intent(inout), optional :: i1x
    integer(I2B),     dimension(:,:), allocatable, intent(inout), optional :: i2x
    integer(I4B),     dimension(:,:), allocatable, intent(inout), optional :: i4x
    integer(I8B),     dimension(:,:), allocatable, intent(inout), optional :: i8x
    real(R4B),        dimension(:,:), allocatable, intent(inout), optional :: r4x
    real(R8B),        dimension(:,:), allocatable, intent(inout), optional :: r8x
    character(len=*), dimension(:,:), allocatable, intent(inout), optional :: cx
    ! -- locals
    integer(I1B),          dimension(:), allocatable :: i1a
    integer(I2B),          dimension(:), allocatable :: i2a
    integer(I4B),          dimension(:), allocatable :: i4a
    integer(I8B),          dimension(:), allocatable :: i8a
    real(R4B),             dimension(:), allocatable :: r4a
    real(R8B),             dimension(:), allocatable :: r8a
    character(len=MXSLEN), dimension(:), allocatable :: ca
    type(tVal), pointer :: v => null() 
    character(len=MXSLEN) :: s
    integer(I4B) ::nr, nc, ic
! ------------------------------------------------------------------------------
    nr = this%get_nr(); nc = this%get_nc()
    !
    if (present(i1x)) then
      if (allocated(i1x)) deallocate(i1x); allocate(i1x(nc,nr))
      do ic = 1, nc
        call this%get_column(ic=ic, i1a=i1a); i1x(ic,:) = i1a
      end do
    end if
    if (present(i2x)) then
      if (allocated(i2x)) deallocate(i2x); allocate(i2x(nc,nr))
      do ic = 1, nc
        call this%get_column(ic=ic, i2a=i2a); i2x(ic,:) = i2a
      end do
    end if
    if (present(i4x)) then
      if (allocated(i4x)) deallocate(i4x); allocate(i4x(nc,nr))
      do ic = 1, nc
        call this%get_column(ic=ic, i4a=i4a); i4x(ic,:) = i4a
      end do
    end if
    if (present(i8x)) then
      if (allocated(i8x)) deallocate(i8x); allocate(i8x(nc,nr))
      do ic = 1, nc
        call this%get_column(ic=ic, i8a=i8a); i8x(ic,:) = i8a
      end do
    end if
    if (present(r4x)) then
      if (allocated(r4x)) deallocate(r4x); allocate(r4x(nc,nr))
      do ic = 1, nc
        call this%get_column(ic=ic, r4a=r4a); r4x(ic,:) = r4a
      end do
    end if
    if (present(r8x)) then
      if (allocated(r8x)) deallocate(r8x); allocate(r8x(nc,nr))
      do ic = 1, nc
        call this%get_column(ic=ic, r8a=r8a); r8x(ic,:) = r8a
      end do
    end if
    if (present(cx)) then
      if (allocated(cx)) deallocate(cx); allocate(cx(nc,nr))
      do ic = 1, nc
        call this%get_column(ic=ic, ca=ca); cx(ic,:) = ca
      end do
    end if
    !
    if (allocated(i1a)) deallocate(i1a)
    if (allocated(i2a)) deallocate(i2a)
    if (allocated(i4a)) deallocate(i4a)
    if (allocated(i8a)) deallocate(i8a)
    if (allocated(r4a)) deallocate(r4a)
    if (allocated(r4a)) deallocate(r4a)
    if (allocated(ca))  deallocate(ca)
    !
    return
  end subroutine tCSV_get_matrix
  
  subroutine tCSV_get_val(this, ir, ic, key, &
    i1v, i2v, i4v, i8v, r4v, r8v, l4v, cv)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    !
    integer(I4B), intent(in) :: ir
    integer(I4B), intent(in), optional :: ic
    character(len=*), intent(in), optional :: key
    integer(I1B), intent(out), optional :: i1v
    integer(I2B), intent(out), optional :: i2v
    integer(I4B), intent(out), optional :: i4v
    integer(I8B), intent(out), optional :: i8v
    real(R4B), intent(out), optional :: r4v
    real(R8B), intent(out), optional :: r8v
    logical, intent(out), optional :: l4v
    character(len=*), intent(out), optional :: cv
    !
    ! -- locals
    type(tVal), pointer :: v => null()
    character(len=MXSLEN) :: s
    integer(I4B) :: jc
! ------------------------------------------------------------------------------
    !
    if (present(ic)) then
      jc = ic
    else
      jc = this%get_col(key)
    end if
    !
    ! checks
    if ((ir < 1).or.(ir > this%nr)) then
      call errmsg('tCSV_get_val: invalid row number.')
    end if
    if ((jc < 1).or.(jc > this%nc)) then
      call errmsg('tCSV_get_val: invalid column number.')
    end if
    !
    v => this%val(jc,ir)
    call v%get_val(i1v, i2v, i4v, i8v, r4v, r8v, l4v, cv)
    !
    return
  end subroutine tCSV_get_val

  subroutine tCSV_get_selection(this, key, cv, csv_sel, sel_rows)
! ******************************************************************************  
    ! -- arguments
    class(tCSV) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: cv
    type(tCSV), intent(inout), pointer :: csv_sel
    integer(I4B), dimension(:), allocatable, intent(inout) :: sel_rows
    !
    ! -- local
    type(tCSV_hdr), pointer :: hdr => null(), hdr_tgt => null()
    type(tVal), pointer :: val => null(), val_tgt => null()
    character(len=MXSLEN), dimension(:), allocatable :: sa
    integer(I4B) :: ir, nr, ic, nhdr, i
! ------------------------------------------------------------------------------
    ! 
    if (associated(csv_sel)) then
      call csv_sel%clean(); deallocate(csv_sel); csv_sel => null()
    end if
    if (allocated(sel_rows)) then
      deallocate(sel_rows)
    end if
    !
    call this%get_column(key=key, ca=sa)
    nr = 0
    do ir = 1, this%nr
      if (sa(ir) == cv) then
        nr = nr + 1
      end if
    end do
    !
    if (nr == 0) then
      return
    end if
    allocate(csv_sel, sel_rows(nr))
    !
    csv_sel%file = this%file
    csv_sel%nr = nr; csv_sel%nr_max = nr; csv_sel%nc = this%nc
    nhdr = size(this%hdr)
    allocate(csv_sel%hdr(nhdr))
    do i = 1, nhdr
      hdr => this%hdr(i); hdr_tgt => csv_sel%hdr(i)
      call hdr%copy(hdr_tgt)
    end do
    allocate(csv_sel%val(this%nc,nr))
    !
    nr = 0
    do ir = 1, this%nr
      if (sa(ir) == cv) then
        nr = nr + 1
        sel_rows(nr) = ir
        do ic = 1, this%nc
          val => this%val(ic,ir); val_tgt => csv_sel%val(ic,nr)
          call val%copy(val_tgt)
        end do
      end if
    end do
    !
    return
  end subroutine tCSV_get_selection
    
! ==============================================================================
! ==============================================================================
! generic functions
! ==============================================================================
! ==============================================================================
  
  subroutine cast_number_from_string(cv, i1v, i2v, i4v, i8v, r4v, r8v, l4v, &
    i1mv, i2mv, i4mv, i8mv, r4mv, r8mv, l4mv)
! ******************************************************************************  
    ! -- arguments
    character(len=*), intent(in), optional :: cv
    !
    integer(I1B),     intent(in), optional :: i1mv
    integer(I2B),     intent(in), optional :: i2mv
    integer(I4B),     intent(in), optional :: i4mv
    integer(I8B),     intent(in), optional :: i8mv
    real(R4B),        intent(in), optional :: r4mv
    real(R8B),        intent(in), optional :: r8mv
    logical,          intent(in), optional :: l4mv
    ! 
    integer(I1B),     intent(out), optional :: i1v
    integer(I2B),     intent(out), optional :: i2v
    integer(I4B),     intent(out), optional :: i4v
    integer(I8B),     intent(out), optional :: i8v
    real(R4B),        intent(out), optional :: r4v
    real(R8B),        intent(out), optional :: r8v
    logical,          intent(out), optional :: l4v
    ! -- locals
    integer(I4B) :: ios, n, ierr
! ------------------------------------------------------------------------------
    n = 0
    if (present(i1v)) then
      read(cv,*,iostat=ios) i1v; n = n + 1
    end if
    if (present(i2v)) then
      read(cv,*,iostat=ios) i2v; n = n + 1
    end if
    if (present(i4v)) then
      read(cv,*,iostat=ios) i4v; n = n + 1
    end if
    if (present(i8v)) then
      read(cv,*,iostat=ios) i8v; n = n + 1
    end if
    if (present(r4v)) then
      read(cv,*,iostat=ios) r4v; n = n + 1
    end if
    if (present(r8v)) then
      read(cv,*,iostat=ios) r8v; n = n + 1
    end if
    if (present(l4v)) then
      read(cv,*,iostat=ios) l4v; n = n + 1
    end if
    if (n /= 1) then
      call errmsg('cast_number_from_string: too many input arguments.')
    end if
    !
    ierr = 0
    if (ios /= 0) then
      if (present(i1v)) then
        if (present(i1mv)) then
          i1v = i1mv
        else
          ierr = 1
        end if
      end if
      if (present(i2v)) then
        if (present(i2mv)) then
          i2v = i2mv
        else
          ierr = 1
        end if
      end if
      if (present(i4v)) then
        if (present(i4mv)) then
          i4v = i4mv
        else
          ierr = 1
        end if
      end if
      if (present(i8v)) then
        if (present(i8mv)) then
          i8v = i8mv
        else
          ierr = 1
        end if
      end if
      if (present(r4v)) then
        if (present(r4mv)) then
          r4v = r4mv
        else
          ierr = 1
        end if
      end if
      if (present(r8v)) then
        if (present(r8mv)) then
          r8v = r8mv
        else
          ierr = 1
        end if
      end if
      if (present(l4v)) then
        if (present(l4mv)) then
          l4v = l4mv
        else
          ierr = 1
        end if
      end if
    end if
    if (ierr == 1) then
      call errmsg('cast_number_from_string: missing data not found.')
    end if
    !
    return
  end subroutine cast_number_from_string
    
  subroutine cast_number(i1v_out, i2v_out, i4v_out, i8v_out, &
    r4v_out, r8v_out, cv_out, &
    i1v_in, i2v_in, i4v_in, i8v_in, r4v_in, r8v_in, cv_in)
! ******************************************************************************  
    ! -- arguments
    integer(I1B),     intent(out) :: i1v_out
    integer(I2B),     intent(out) :: i2v_out
    integer(I4B),     intent(out) :: i4v_out
    integer(I8B),     intent(out) :: i8v_out
    real(R4B),        intent(out) :: r4v_out
    real(R8B),        intent(out) :: r8v_out
    character(len=*), intent(out) :: cv_out
    ! 
    integer(I1B),     intent(in), optional :: i1v_in
    integer(I2B),     intent(in), optional :: i2v_in
    integer(I4B),     intent(in), optional :: i4v_in
    integer(I8B),     intent(in), optional :: i8v_in
    real(R4B),        intent(in), optional :: r4v_in
    real(R8B),        intent(in), optional :: r8v_in
    character(len=*), intent(in), optional :: cv_in
    ! -- locals
    integer(I4B) :: ios
! ------------------------------------------------------------------------------
    if (present(i1v_in)) then
      i1v_out = int(i1v_in,I1B); i2v_out = int(i1v_in,I2B)
      i4v_out = int(i1v_in,I4B); i8v_out = int(i1v_in,I8B)
      r4v_out = int(i1v_in,R4B); r8v_out = int(i1v_in,R8B)
      write(cv_out,*) i1v_in; cv_out = adjustl(cv_out)
      return
    end if
    if (present(i2v_in)) then
      i1v_out = int(i2v_in,I1B); i2v_out = int(i2v_in,I2B)
      i4v_out = int(i2v_in,I4B); i8v_out = int(i2v_in,I8B)
      r4v_out = int(i2v_in,R4B); r8v_out = int(i2v_in,R8B)
      write(cv_out,*) i2v_in; cv_out = adjustl(cv_out)
      return
    end if
    if (present(i4v_in)) then
      i1v_out = int(i4v_in,I1B); i2v_out = int(i4v_in,I2B)
      i4v_out = int(i4v_in,I4B); i8v_out = int(i4v_in,I8B)
      r4v_out = int(i4v_in,R4B); r8v_out = int(i4v_in,R8B)
      write(cv_out,*) i4v_in; cv_out = adjustl(cv_out)
      return
    end if
    if (present(i8v_in)) then
      i1v_out = int(i8v_in,I1B); i2v_out = int(i8v_in,I2B)
      i4v_out = int(i8v_in,I4B); i8v_out = int(i8v_in,I8B)
      r4v_out = int(i8v_in,R4B); r8v_out = int(i8v_in,R8B)
      write(cv_out,*) i8v_in; cv_out = adjustl(cv_out)
      return
    end if
    if (present(r4v_in)) then
      i1v_out = int(r4v_in,I1B); i2v_out = int(r4v_in,I2B)
      i4v_out = int(r4v_in,I4B); i8v_out = int(r4v_in,I8B)
      r4v_out = int(r4v_in,R4B); r8v_out = int(r4v_in,R8B)
      write(cv_out,*) r4v_in; cv_out = adjustl(cv_out)
      return
    end if
    if (present(r8v_in)) then
      i1v_out = int(r8v_in,I1B); i2v_out = int(r8v_in,I2B)
      i4v_out = int(r8v_in,I4B); i8v_out = int(r8v_in,I8B)
      r4v_out = int(r8v_in,R4B); r8v_out = int(r8v_in,R8B)
      write(cv_out,*) r8v_in; cv_out = adjustl(cv_out)
      return
    end if
    if (present(cv_in)) then
      read(cv_in,*,iostat=ios) i1v_out
      if (ios /= 0) call logmsg('cast_number: could not convert to i1')
      read(cv_in,*,iostat=ios) i2v_out
      if (ios /= 0) call logmsg('cast_number: could not convert to i2')
      read(cv_in,*,iostat=ios) i4v_out
      if (ios /= 0) call logmsg('cast_number: could not convert to i4')
      read(cv_in,*,iostat=ios) i8v_out
      if (ios /= 0) call logmsg('cast_number: could not convert to i8')
      read(cv_in,*,iostat=ios) r4v_out
      if (ios /= 0) call logmsg('cast_number: could not convert to r4')
      read(cv_in,*,iostat=ios) r8v_out
      if (ios /= 0) call logmsg('cast_number: could not convert to r8')
      cv_out = cv_in
    end if
    !
    return
  end subroutine cast_number
    
  subroutine cast_2darray(xi1_in, xi2_in, xi4_in, xi8_in, xr4_in, xr8_in, &
    xi1_out, xi2_out, xi4_out, xi8_out, xr4_out, xr8_out, &
    mvi1_in, mvi2_in, mvi4_in, mvi8_in, mvr4_in, mvr8_in, &
    mvi1_out, mvi2_out, mvi4_out, mvi8_out, mvr4_out, mvr8_out)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
 ! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I1B), dimension(:,:), intent(in), allocatable, optional :: xi1_in
    integer(I2B), dimension(:,:), intent(in), allocatable, optional :: xi2_in
    integer(I4B), dimension(:,:), intent(in), allocatable, optional :: xi4_in
    integer(I8B), dimension(:,:), intent(in), allocatable, optional :: xi8_in
    real(R4B),    dimension(:,:), intent(in), allocatable, optional :: xr4_in
    real(R8B),    dimension(:,:), intent(in), allocatable, optional :: xr8_in
    !
    integer(I1B), dimension(:,:), intent(inout), allocatable, optional :: xi1_out
    integer(I2B), dimension(:,:), intent(inout), allocatable, optional :: xi2_out
    integer(I4B), dimension(:,:), intent(inout), allocatable, optional :: xi4_out
    integer(I8B), dimension(:,:), intent(inout), allocatable, optional :: xi8_out
    real(R4B),    dimension(:,:), intent(inout), allocatable, optional :: xr4_out
    real(R8B),    dimension(:,:), intent(inout), allocatable, optional :: xr8_out
    !
    integer(I1B), intent(in), optional :: mvi1_in
    integer(I2B), intent(in), optional :: mvi2_in
    integer(I4B), intent(in), optional :: mvi4_in
    integer(I8B), intent(in), optional :: mvi8_in
    real(R4B),    intent(in), optional :: mvr4_in
    real(R8B),    intent(in), optional :: mvr8_in
    !
    integer(I1B), intent(out), optional :: mvi1_out
    integer(I2B), intent(out), optional :: mvi2_out
    integer(I4B), intent(out), optional :: mvi4_out
    integer(I8B), intent(out), optional :: mvi8_out
    real(R4B),    intent(out), optional :: mvr4_out
    real(R8B),    intent(out), optional :: mvr8_out
     ! -- locals
    logical :: lmv
    character(len=MXSLEN) :: s
    integer(I1B) :: mvi1, i1v
    integer(I2B) :: mvi2, i2v
    integer(I4B) :: mvi4, i4v
    integer(I8B) :: mvi8, i8v
    real(R4B)    :: mvr4, r4v
    real(R8B)    :: mvr8, r8v
    !
    integer(I4B) :: n, nc, nr, ic, ir
    integer(I4B), dimension(6) :: flg_in, flg_out
! ------------------------------------------------------------------------------
    !
    ! checks
    flg_in = 0
    if (present(xi1_in)) then
      flg_in(1) = 1
      if (.not.present(mvi1_in)) call errmsg('cast_2darray: missing mvi1.')
    end if
    if (present(xi2_in)) then
      flg_in(2) = 1
      if (.not.present(mvi2_in)) call errmsg('cast_2darray: missing mvi2.')
    end if
    if (present(xi4_in)) then
      flg_in(3) = 1
      if (.not.present(mvi4_in)) call errmsg('cast_2darray: missing mvi4.')
    end if
    if (present(xi8_in)) then
      flg_in(4) = 1
      if (.not.present(mvi8_in)) call errmsg('cast_2darray: missing mvi8.')
    end if
    if (present(xr4_in)) then
      flg_in(5) = 1
      if (.not.present(mvr4_in)) call errmsg('cast_2darray: missing mvr4.')
    end if
    if (present(xr8_in)) then
      flg_in(6) = 1
      if (.not.present(mvr8_in)) call errmsg('cast_2darray: missing mvr8.')
    end if
    if (sum(flg_in) == 0) then
      call logmsg('cast_2darray: no input; nothing to do...')
      return
    end if
    if (sum(flg_in) > 1) then
      call logmsg('cast_2darray: too many input arrays.')
    end if
    
    flg_out = 0
    if (present(xi1_out)) flg_out(1) = 1
    if (present(xi2_out)) flg_out(2) = 1 
    if (present(xi4_out)) flg_out(3) = 1
    if (present(xi8_out)) flg_out(4) = 1
    if (present(xr4_out)) flg_out(5) = 1
    if (present(xr8_out)) flg_out(6) = 1
    if (sum(flg_out) == 0) then
      call logmsg('cast_2darray: no output; nothing to do...')
      return
    end if
    if (sum(flg_out) > 1) then
      call logmsg('cast_2darray: too many output arrays.')
    end if
    !
    call cast_number(mvi1, mvi2, mvi4, mvi8, mvr4, mvr8, s, &
      mvi1_in, mvi2_in, mvi4_in, mvi8_in, mvr4_in, mvr8_in)
    !
    nc = 0; nr = 0
    if (present(xi1_in)) then; nc = size(xi1_in,1); nr = size(xi1_in,2); end if
    if (present(xi2_in)) then; nc = size(xi2_in,1); nr = size(xi2_in,2); end if
    if (present(xi4_in)) then; nc = size(xi4_in,1); nr = size(xi4_in,2); end if
    if (present(xi8_in)) then; nc = size(xi8_in,1); nr = size(xi8_in,2); end if
    if (present(xr4_in)) then; nc = size(xr4_in,1); nr = size(xr4_in,2); end if
    if (present(xr8_in)) then; nc = size(xr8_in,1); nr = size(xr8_in,2); end if
    if ((nc == 0).or.(nr == 0)) then
      call errmsg('cast_2darray: invalid nr/nc.')
    end if
    !
    if (present(xi1_out)) then; if (allocated(xi1_out)) deallocate(xi1_out); end if
    if (present(xi2_out)) then; if (allocated(xi2_out)) deallocate(xi2_out); end if
    if (present(xi4_out)) then; if (allocated(xi4_out)) deallocate(xi4_out); end if
    if (present(xi8_out)) then; if (allocated(xi8_out)) deallocate(xi8_out); end if
    if (present(xr4_out)) then; if (allocated(xr4_out)) deallocate(xr4_out); end if
    if (present(xr8_out)) then; if (allocated(xr8_out)) deallocate(xr8_out); end if
    !
    if (present(xi1_out)) allocate(xi1_out(nc,nr))
    if (present(xi2_out)) allocate(xi2_out(nc,nr))
    if (present(xi4_out)) allocate(xi4_out(nc,nr))
    if (present(xi8_out)) allocate(xi8_out(nc,nr))
    if (present(xr4_out)) allocate(xr4_out(nc,nr))
    if (present(xr8_out)) allocate(xr8_out(nc,nr))
    !
    if (present(xi1_out)) then
      call errmsg('Not yet supported')
    end if
    if (present(xi2_out)) then
      call errmsg('Not yet supported')
    end if
    if (present(xi4_out)) then
      mvi4_out = mvi4
      do ir = 1, nr; do ic = 1, nc
        if (flg_in(1)==1) i4v = int(xi1_in(ic,ir),I4B)
        if (flg_in(2)==1) i4v = int(xi2_in(ic,ir),I4B)
        if (flg_in(3)==1) i4v = int(xi4_in(ic,ir),I4B)
        if (flg_in(4)==1) i4v = int(xi8_in(ic,ir),I4B)
        if (flg_in(5)==1) i4v = int(xr4_in(ic,ir),I4B)
        if (flg_in(6)==1) i4v = int(xr8_in(ic,ir),I4B)
        !
        if (flg_in(1)==1) lmv = (xi1_in(ic,ir) == mvi1_in)
        if (flg_in(2)==1) lmv = (xi2_in(ic,ir) == mvi2_in)
        if (flg_in(3)==1) lmv = (xi4_in(ic,ir) == mvi4_in)
        if (flg_in(4)==1) lmv = (xi8_in(ic,ir) == mvi8_in)
        if (flg_in(5)==1) lmv = (xr4_in(ic,ir) == mvr4_in)
        if (flg_in(6)==1) lmv = (xr8_in(ic,ir) == mvr8_in)
        !
        if (lmv) i4v = mvi4_out; xi4_out(ic,ir) = i4v
      end do; end do
    end if
    if (present(xi8_out)) then
      call errmsg('Not yet supported')
    end if
    if (present(xr4_out)) then
      mvr4_out = mvr4
      do ir = 1, nr; do ic = 1, nc
        if (flg_in(1)==1) r4v = real(xi1_in(ic,ir),R4B)
        if (flg_in(2)==1) r4v = real(xi2_in(ic,ir),R4B)
        if (flg_in(3)==1) r4v = real(xi4_in(ic,ir),R4B)
        if (flg_in(4)==1) r4v = real(xi8_in(ic,ir),R4B)
        if (flg_in(5)==1) r4v = real(xr4_in(ic,ir),R4B)
        if (flg_in(6)==1) r4v = real(xr8_in(ic,ir),R4B)
        !
        if (lmv) r4v = mvr4_out; xr4_out(ic,ir) = r4v
        !if (rep_mv.and.(r4v == mvr4_in)) r4v = mvr4; xr4_out(ic,ir) = r4v
      end do; end do
    end if
    if (present(xr8_out)) then
      call errmsg('Not yet supported')
    end if
    !
    return
  end subroutine cast_2darray
    
  function remove_trailing_zeros(s_in) result(s_out)
! ******************************************************************************  
    ! -- arguments
    character(len=*), intent(in) :: s_in
    character(len=MXSLEN) :: s_out
    ! -- locals
    integer(I4B) :: i, j, n
! ------------------------------------------------------------------------------
    !
    n = len_trim(s_in)
    !
    if (n == 0) then
      s_out = ''
      return
    end if
    !
    i = index(s_in, '.')
    if (i == 0) then
      s_out = s_in
      return
    end if
    !
    i = index(s_in, 'e+000')
    if (i == 0) then
      i = n
    else
      i = i - 1
    end if
    if (s_in(i:i) == '0') then
      do while(.true.)
        j = i - 1
        if (j > 0) then
          if (s_in(j:j) /= '0') exit
        else
          exit
        end if
        i = i - 1
      end do
    end if
    s_out = s_in(1:i)
    !
    return
  end function remove_trailing_zeros
  
  subroutine read_line(iu, s, ios, comment)
! ******************************************************************************  
    ! -- arguments
    integer(I4B), intent(in) :: iu
    character(len=MXSLEN), intent(out) :: s
    integer(I4B), intent(out) :: ios
    character(len=*), intent(in), dimension(:), optional :: comment
    ! -- locals
    logical :: lcomment, chk_comment
    integer(I4B) :: i, n
! ------------------------------------------------------------------------------
    chk_comment = .false.
    if (present(comment)) then
      chk_comment = .true.
    end if
    !
    s = ''
    do while(.true.)
      read(unit=iu,iostat=ios,fmt='(a)') s
      if (ios /= 0) then
        exit
      end if
      s = adjustl(s)
      !
      lcomment = .false.
      if (chk_comment) then
        do i = 1, size(comment)
          n = len_trim(comment(i))
          if (s(1:n) == trim(comment(i))) then
            lcomment = .true.
          end if
        end do
      end if
      !
      if ((.not.lcomment).and.(len_trim(s) /= 0)) then
        exit
      end if
    end do
  end subroutine read_line
  
  subroutine get_xy(x, y, ic, ir, xll, yur, cs)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    real(R8B), intent(out) :: x
    real(R8B), intent(out) :: y
    integer(I4B), intent(in) :: ic
    integer(I4B), intent(in) :: ir
    real(R8B), intent(in) :: xll
    real(R8B), intent(in) :: yur
    real(R8B), intent(in) :: cs
    ! -- local
! ------------------------------------------------------------------------------
    x = xll + ic*cs - R8HALF*cs
    y = yur - ir*cs + R8HALF*cs
    !
    return
  end subroutine get_xy
  
  subroutine get_icr(ic, ir, x, y, xll, yur, cs)
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
    real(R8B), intent(in) :: yur
    real(R8B), intent(in) :: cs
    ! -- local
    real(R8B) :: dx, dy
! ------------------------------------------------------------------------------
    ic = 0; ir = 0
    dx = x - xll; dy = yur - y
    if ((dx < R8ZERO).or.(dy < DZERO)) then
      return
    end if
    !
    ic = int(dx/cs)+1
    ir = int(dy/cs)+1
    !
    if ((ic == 0).or.(ir == 0)) then
      call logmsg('get_icr: ic=0 or ir=0 detected.')
    end if
    
    return
  end subroutine get_icr
  
  function valid_icr(ic, ir, nc, nr) result(valid)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), intent(in) :: ic
    integer(I4B), intent(in) :: ir
    integer(I4B), intent(in) :: nc
    integer(I4B), intent(in) :: nr
    logical :: valid
! ------------------------------------------------------------------------------
    if ((ic >= 1).and.(ic <= nc).and.(ir >= 1).and.(ir <= nr)) then
      valid = .true.
    else
      valid = .false.
    end if
    !
    return
  end function valid_icr
  
  subroutine get_mapped_icr(bbx_src, icir_src, bbx_tgt, icir_tgt)
! ******************************************************************************
    ! -- arguments
    type(tBbx), intent(in) :: bbx_src
    integer(I4B), dimension(2), intent(in) :: icir_src
    type(tBbx), intent(in) :: bbx_tgt
    integer(I4B), dimension(:,:), allocatable, intent(inout) :: icir_tgt
  
    ! -- locals
    integer(I4B) :: ic, ir, ic0, ic1, ir0, ir1, nc, nr, n
    real(R8B) :: x, y, xmin, xmax, ymin, ymax, halfcs, dxy
! ------------------------------------------------------------------------------
    !
    ! get the source mid-point coordinate
    call get_xy(x, y, icir_src(1), icir_src(2), &
      bbx_src%xll, bbx_src%yur, bbx_src%cs)
    !
    ! clean up
    if (allocated(icir_tgt)) deallocate(icir_tgt)
    !
    if (bbx_src%cs < bbx_tgt%cs) then ! source fine, target coarse
      allocate(icir_tgt(2,1))
      call get_icr(icir_tgt(1,1), icir_tgt(2,1), x, y, &
        bbx_tgt%xll, bbx_tgt%yur, bbx_tgt%cs)
    else ! source coarse, target fine
      halfcs = R8HALF*bbx_src%cs
      xmin = max(x-halfcs, bbx_tgt%xll); xmax = min(x+halfcs, bbx_tgt%xur)
      ymin = max(y-halfcs, bbx_tgt%yll); ymax = min(y+halfcs, bbx_tgt%yur)
      !
      ! slightly disturb
      dxy = 0.001d0*bbx_tgt%cs ! 1%% of the cell-sze
      xmin = xmin + dxy; xmax = xmax - dxy
      ymin = ymin + dxy; ymax = ymax - dxy
      !
      call get_icr(ic0, ir1, xmin, ymin, bbx_tgt%xll, bbx_tgt%yur, bbx_tgt%cs)
      call get_icr(ic1, ir0, xmax, ymax, bbx_tgt%xll, bbx_tgt%yur, bbx_tgt%cs)
      !
      nc = ic1 - ic0 + 1; nr = ir1 - ir0 + 1
      !if ((nc > 2).or.(nr > 2)) then
      !  call logmsg('get_mapped_icr: nc > 2 or nr > 2')
      !end if
      n = nc*nr; allocate(icir_tgt(2,n)); n = 0
      do ir = ir0, ir1
        do ic = ic0, ic1
          n = n + 1
          icir_tgt(1,n) = ic; icir_tgt(2,n) = ir
        end do
      end do
    end if
    !
    return
  end subroutine get_mapped_icr
  
  subroutine icrl_to_node(n, ic, ir, il, nc, nr)
! ******************************************************************************  
    ! -- arguments
    integer(I4B), intent(out) :: n
    integer(I4B), intent(in) :: ic
    integer(I4B), intent(in) :: ir
    integer(I4B), intent(in) :: il
    integer(I4B), intent(in) :: nr
    integer(I4B), intent(in) :: nc
    ! -- locals
    integer(I4B) :: nrc
! ------------------------------------------------------------------------------
    !
    nrc = nr*nc
    n = ic + (ir-1)*nc + (il-1)*nrc
    !
    return
  end subroutine icrl_to_node
  
  subroutine node_to_icrl(n, ic, ir, il, nc, nr)
! ******************************************************************************  
    ! -- arguments
    integer(I4B), intent(in) :: n
    integer(I4B), intent(out) :: ic
    integer(I4B), intent(out) :: ir
    integer(I4B), intent(out) :: il
    integer(I4B), intent(in) :: nr
    integer(I4B), intent(in) :: nc
    ! -- locals
    integer(I4B) :: nrc
! ------------------------------------------------------------------------------
    !
    nrc = nr*nc
    !
    il = int((n - 1)/nrc) + 1
    ir = int((n - (il-1)*nrc - 1)/nc) + 1
    ic = n - (ir-1)*nc - (il-1)*nrc
    !
    return
  end subroutine node_to_icrl
  
  subroutine tBbObj_init(this)
! ******************************************************************************  
    ! -- arguments
    class(tBbObj) :: this
! ------------------------------------------------------------------------------
    call this%prent_bbi%init()
    call this%prent_bbx%init()
    call this%child_bbi%init()
    call this%child_bbx%init()
    !
    return
  end subroutine tBbObj_init
    
  subroutine tBbObj_set(this, prent_bbi, child_bbi, prent_bbx, child_bbx)
! ******************************************************************************  
    ! -- arguments
    class(tBbObj) :: this
    ! -- locals
    type(tBb), intent(in), optional :: prent_bbi
    type(tBb), intent(in), optional :: child_bbi
    type(tBbX), intent(in), optional :: prent_bbx
    type(tBbX), intent(in), optional :: child_bbx
! ------------------------------------------------------------------------------
    !
    if (present(prent_bbi)) then
      this%prent_bbi = prent_bbi
    end if
    if (present(child_bbi)) then
      this%child_bbi = child_bbi
    end if
    if (present(prent_bbx)) then
      this%prent_bbx = prent_bbx
    end if
    if (present(child_bbx)) then
      this%child_bbx = child_bbx
    end if
    !
    return
  end subroutine tBbObj_set
  
! ==============================================================================
! ==============================================================================
! tNum
! ==============================================================================
! ==============================================================================
  subroutine tNum_init(this)
! ******************************************************************************  
    ! -- arguments
    class(tNum) :: this
    ! -- locals
! ------------------------------------------------------------------------------
    this%flg = .false.
    
    this%i1v = I1ZERO
    this%i2v = I2ZERO
    this%i4v = I4ZERO
    this%i8v = I8ZERO
    this%r4v = R4ZERO
    this%r8v = R8ZERO
    this%l4v = .false.
    !
    return
  end subroutine tNum_init
  
  subroutine tNum_set_val(this, i1v, i2v, i4v, i8v, r4v, r8v, l4v)
! ******************************************************************************  
    ! -- arguments
    class(tNum) :: this
    !
    integer(I1B), intent(in), optional :: i1v
    integer(I2B), intent(in), optional :: i2v
    integer(I4B), intent(in), optional :: i4v
    integer(I8B), intent(in), optional :: i8v
    real(R4B),    intent(in), optional :: r4v
    real(R8B),    intent(in), optional :: r8v
    logical,      intent(in), optional :: l4v
    ! -- locals
! ------------------------------------------------------------------------------
    !
    this%flg = .false.
    !
    if (present(i1v)) then
      this%flg(i_i1) = .true.; this%i1v = i1v
    end if
    if (present(i2v)) then
      this%flg(i_i2) = .true.; this%i2v = i2v
    end if
    if (present(i4v)) then
      this%flg(i_i4) = .true.; this%i4v = i4v
    end if
    if (present(i8v)) then
      this%flg(i_i8) = .true.; this%i8v = i8v
    end if
    if (present(r4v)) then
      this%flg(i_r4) = .true.; this%r4v = r4v
    end if
    if (present(r8v)) then
      this%flg(i_r8) = .true.; this%r8v = r8v
    end if
    if (present(l4v)) then
      this%flg(i_l4) = .true.; this%l4v = l4v
    end if
    !
    return
  end subroutine tNum_set_val
  
  subroutine tNum_get_val(this, i1v, i2v, i4v, i8v, r4v, r8v, l4v, lset)
! ******************************************************************************  
    ! -- arguments
    class(tNum) :: this
    integer(I1B), intent(out), optional :: i1v
    integer(I2B), intent(out), optional :: i2v
    integer(I4B), intent(out), optional :: i4v
    integer(I8B), intent(out), optional :: i8v
    real(R4B),    intent(out), optional :: r4v
    real(R8B),    intent(out), optional :: r8v
    logical,      intent(out), optional :: l4v
    logical, intent(out) :: lset
    ! -- locals
    logical :: flg
! ------------------------------------------------------------------------------
    !
    ! check the flag
    lset = .false.
    if (present(i1v)) lset = this%flg(i_i1)
    if (present(i2v)) lset = this%flg(i_i2)
    if (present(i4v)) lset = this%flg(i_i4)
    if (present(i8v)) lset = this%flg(i_i8)
    if (present(r4v)) lset = this%flg(i_r4)
    if (present(r8v)) lset = this%flg(i_r8)
    if (present(l4v)) lset = this%flg(i_l4)
    !
    if (.not.lset) then
      !return
      call errmsg('tNum_get_val: data not set.')
    end if
    !
    if (present(i1v)) i1v = this%i1v
    if (present(i2v)) i2v = this%i2v
    if (present(i4v)) i4v = this%i4v
    if (present(i8v)) i8v = this%i8v
    if (present(r4v)) r4v = this%r4v
    if (present(r8v)) r8v = this%r8v
    if (present(l4v)) l4v = this%l4v
    !
    return
  end subroutine tNum_get_val
  
  subroutine tNum_set_val_by_s(this, s)
! ******************************************************************************  
    ! -- arguments
    class(tNum) :: this
    character(len=*), intent(in) :: s
    !
    integer(I1B) :: i1v
    integer(I2B) :: i2v
    integer(I4B) :: i4v
    integer(I8B) :: i8v
    real(R4B) :: r4v
    real(R8B) :: r8v
    logical :: l4v
    !
    integer(I4B) :: ios
! ------------------------------------------------------------------------------
    !
    this%flg = .false.
    !
    read(s,*,iostat=ios) i1v
    if (ios == 0) then
      this%i1v = i1v; this%flg(i_i1) = .true.
    end if
    read(s,*,iostat=ios) i2v
    if (ios == 0) then
      this%i2v = i2v; this%flg(i_i2) = .true.
    end if
    read(s,*,iostat=ios) i4v
    if (ios == 0) then
      this%i4v = i4v; this%flg(i_i4) = .true.
    end if
    read(s,*,iostat=ios) i8v
    if (ios == 0) then
      this%i8v = i8v; this%flg(i_i8) = .true.
    end if
    read(s,*,iostat=ios) r4v
    if (ios == 0) then
      this%r4v = r4v; this%flg(i_r4) = .true.
    end if
    read(s,*,iostat=ios) r8v
    if (ios == 0) then
      this%r8v = r8v; this%flg(i_r8) = .true.
    end if
    read(s,*,iostat=ios) l4v
    if (ios == 0) then
      this%l4v = l4v; this%flg(i_l4) = .true.
    end if
    !
    return
  end subroutine tNum_set_val_by_s
  
  subroutine tNum_copy(this, tgt)
! ******************************************************************************  
    ! -- arguments
    class(tNum) :: this
    type(tNum), pointer, intent(inout) :: tgt
! ------------------------------------------------------------------------------
    
    tgt%flg = this%flg
    tgt%i1v = this%i1v
    tgt%i2v = this%i2v
    tgt%i4v = this%i4v
    tgt%i8v = this%i8v
    tgt%r4v = this%r4v
    tgt%r8v = this%r8v
    tgt%l4v = this%l4v
    !
    return
  end subroutine tNum_copy
  
! ==============================================================================
! ==============================================================================
!
! ==============================================================================
! ==============================================================================
  
  subroutine timeseries_clean(this)
! ******************************************************************************  
    ! -- arguments
    class(tTimeSeries) :: this
! ------------------------------------------------------------------------------
    if (associated(this%act))    deallocate(this%act)
    if (associated(this%read))   deallocate(this%read)
    if (associated(this%rawhdr)) deallocate(this%rawhdr)
    if (associated(this%raw))    deallocate(this%raw)
    if (associated(this%id))     deallocate(this%id)
    if (associated(this%x))      deallocate(this%x)
    if (associated(this%y))      deallocate(this%y)
    if (associated(this%ic))     deallocate(this%ic)
    if (associated(this%ir) )    deallocate(this%ir)
    if (associated(this%im) )    deallocate(this%im)
    if (associated(this%nod))    deallocate(this%nod)
    if (associated(this%glev))   deallocate(this%glev)
    if (associated(this%val))    deallocate(this%val)
    if (associated(this%sm_corr))deallocate(this%sm_corr)
    if (associated(this%nlay))   deallocate(this%nlay)
    !
    this%act    => null()
    this%rawhdr => null()
    this%raw    => null()
    this%id     => null()
    this%x      => null()
    this%y      => null()
    this%ic     => null()
    this%ir     => null()
    this%im     => null()
    this%nod    => null()
    this%glev   => null()
    this%val    => null()
    this%sm_corr=> null()
    this%nlay   => null()
    !
    return
  end subroutine timeseries_clean
  
  function point_in_bb(xy, bbx) result(l_in)
! ******************************************************************************
    ! -- arguments
    real(R8B), dimension(:,:), intent(in) :: xy
    type(tBbx), intent(in) :: bbx
    logical :: l_in
    ! -- locals
    integer(I4B) :: i
    real(R8B) :: x, y
! ------------------------------------------------------------------------------
    l_in = .true.
    do i = 1, size(xy,2)
      x = xy(1,i); y = xy(2,i)
      if (x < bbx%xll) then
        l_in = .false.; exit
      end if
      if (x > bbx%xur) then
        l_in = .false.; exit
      end if
      if (y < bbx%yll) then
        l_in = .false.; exit
      end if
      if (y > bbx%yur) then
        l_in = .false.; exit
      end if
    end do
    !
    return
  end function point_in_bb
  
  function bbi_intersect(bbi1, bbi2) result(l_is)
! ******************************************************************************
    ! -- arguments
    type(tBB), intent(in) :: bbi1
    type(tBB), intent(in) :: bbi2
    logical :: l_is
    ! -- locals
! ------------------------------------------------------------------------------
    l_is = .false.
    if (bbi1%ic0 > bbi2%ic1) return
    if (bbi1%ic1 < bbi2%ic0) return
    if (bbi1%ir0 > bbi2%ir1) return
    if (bbi1%ir1 < bbi2%ir0) return
    l_is = .true.
    !
    return
  end function bbi_intersect
  
  function bbx_intersect(bbx1, bbx2) result(l_is)
! ******************************************************************************
    ! -- arguments
    type(tBBX), intent(in) :: bbx1
    type(tBBX), intent(in) :: bbx2
    logical :: l_is
    ! -- locals
! ------------------------------------------------------------------------------
    l_is = .false.
    !
    if (bbx1%xll > bbx2%xur) return
    if (bbx1%xur < bbx2%xll) return
    if (bbx1%yur < bbx2%yll) return
    if (bbx1%yll > bbx2%yur) return
    l_is = .true.
    !
    return
  end function bbx_intersect
  !
  subroutine apply_mask(mask, bbx_m, bbx_x, xi4, xr4, mvi4, mvr4)
! ******************************************************************************
    ! -- arguments
    integer(I4B), dimension(:,:), intent(in) :: mask
    type(tBBX), intent(in) :: bbx_m
    type(tBBX), intent(in) :: bbx_x
    !
    integer(I4B), dimension(:,:), intent(inout), optional :: xi4
    real(R4B),    dimension(:,:), intent(inout), optional :: xr4
    integer(I4B), intent(in), optional :: mvi4
    real(R4B), intent(in), optional :: mvr4
    ! -- locals
    real(R8B) :: x, y
    integer(I4B) :: ir, ic, jr, jc, ncm, nrm, ncx, nrx
! ------------------------------------------------------------------------------
    !
    ncm = size(mask,1); nrm = size(mask,2)
    if (present(xi4)) then
      if (.not.present(mvi4)) call errmsg('apply_mask: mvi4 not found.')
      ncx = size(xi4,1); nrx = size(xi4,2)
    end if
    if (present(xr4)) then
      if (.not.present(mvr4)) call errmsg('apply_mask: mvr4 not found.')
      ncx = size(xr4,1); nrx = size(xr4,2)
    end if
    
    if (bbx_m%cs <= bbx_x%cs) then
      if (present(xi4)) then
        do ir = 1, nrm; do ic = 1, ncm
          if (mask(ic,ir) == 0) then
            call get_xy(x, y, ic, ir, bbx_m%xll, bbx_m%yur, bbx_m%cs)
            call get_icr(jc, jr, x, y, bbx_x%xll, bbx_x%yur, bbx_x%cs)
            if ((jc >= 1).and.(jc <= ncx).and.(jr >= 1).and.(jr <= nrx)) then
              xi4(jc,jr) = mvi4
            end if
          end if
        end do; end do
      end if
      if (present(xr4)) then
        do ir = 1, nrm; do ic = 1, ncm
          if (mask(ic,ir) == 0) then
            call get_xy(x, y, ic, ir, bbx_m%xll, bbx_m%yur, bbx_m%cs)
            call get_icr(jc, jr, x, y, bbx_x%xll, bbx_x%yur, bbx_x%cs)
            if ((jc >= 1).and.(jc <= ncx).and.(jr >= 1).and.(jr <= nrx)) then
              xr4(jc,jr) = mvr4
            end if
          end if
        end do; end do
      end if
    else
      if (present(xi4)) then
        do ir = 1, nrx; do ic = 1, ncx
          if (xi4(ic,ir) /= mvi4) then
            call get_xy(x, y, ic, ir, bbx_x%xll, bbx_x%yur, bbx_x%cs)
            call get_icr(jc, jr, x, y, bbx_m%xll, bbx_m%yur, bbx_m%cs)
            if ((jc >= 1).and.(jc <= ncm).and.(jr >= 1).and.(jr <= nrm)) then
              if (mask(jc,jr) == 0) then
                xi4(ic,ir) = mvi4
              end if
            end if
          end if
        end do; end do
      end if
      if (present(xr4)) then
        do ir = 1, nrx; do ic = 1, ncx
          if (xr4(ic,ir) /= mvr4) then
            call get_xy(x, y, ic, ir, bbx_x%xll, bbx_x%yur, bbx_x%cs)
            call get_icr(jc, jr, x, y, bbx_m%xll, bbx_m%yur, bbx_m%cs)
            if ((jc >= 1).and.(jc <= ncm).and.(jr >= 1).and.(jr <= nrm)) then
              if (mask(jc,jr) == 0) then
                xr4(ic,ir) = mvr4
              end if
            end if
          end if
        end do; end do
      end if
    end if
    !
    return
  end subroutine apply_mask
  
  subroutine remove_tab(s)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: s
    
    ! -- locals
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    do i = 1, len_trim(s)
      if (s(i:i) == achar(9)) then
        s(i:i) = ' '
      end if
    end do
    !
    return
  end subroutine remove_tab

  subroutine insert_tab(s)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: s
    
    ! -- locals
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    do i = 1, len_trim(s)
      if (s(i:i) == ' ') then
        s(i:i) = achar(9)
      end if
    end do
    !
    return
  end subroutine insert_tab
  
  subroutine parse_line(s, sa, i4a, token_in)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: s
    character(len=MXSLEN), dimension(:), allocatable, intent(inout), optional :: sa
    integer(I4B), dimension(:), allocatable, intent(inout), optional :: i4a
    character(len=1), intent(in), optional :: token_in
    ! -- locals
    character(len=1) :: token
    character(len=MXSLENLONG) :: st
    logical :: ltab
    integer(I4B) :: n, m, i, iact
! ------------------------------------------------------------------------------
    if (present(token_in)) then
      token = token_in
    else
      token = ' '
    end if
    !
    ! find tabs and overule token
    ltab = .false.
    do i = 1, len_trim(s)
      if (s(i:i) == achar(9)) then
        ltab = .true.
      end if
    end do
    if (ltab) then
      token = achar(9)
    end if
    !
    do iact = 1, 2
      n = 0; st = adjustl(s)
      m = len_trim(st)
      if (st(m+1:m+1) /= token) then
        st(m+1:m+1) = token
      end if
      do while(.true.)
        m = len_trim(st)
        if (m == 0) then
          exit
        else
          n = n + 1
        end if
        i = index(st,token)
        if ((i <= 0)) then
          exit
        end if
        if (iact == 2) then
          if (present(sa)) then
            sa(n) = st(1:i-1)
          end if
          if (present(i4a)) then
            call cast_number_from_string(cv=st(1:i-1), i4v=i4a(n))
          end if
        end if
        st = adjustl(st(i:))
        if (token /= ' ') then
          if (st(1:1) == token) then
            st = adjustl(st(2:))
          end if
        end if
      end do
      if (iact == 1) then
        if (n > 0) then
          if (present(sa)) then
            if (allocated(sa)) deallocate(sa)
            allocate(sa(n))
          end if
          if (present(i4a)) then
            if (allocated(i4a)) deallocate(i4a)
            allocate(i4a(n))
          end if
        end if
      end if
    end do
    !
    if (n == 0) then
      call errmsg('parse_line: empty string')
    end if
    !
    return
  end subroutine parse_line
  
  function get_args() result(args)
! ******************************************************************************
    ! -- arguments
    character(len=MXSLEN), dimension(:), allocatable :: args
    ! -- locals
    integer(I4B) :: na, i
! ------------------------------------------------------------------------------
    na = nargs()-1
    if (allocated(args)) deallocate(args)
    allocate(args(na))
    do i = 1, na
      call getarg(i, args(i))
    end do
    !
    return
  end function get_args
  
  subroutine linear_regression(x, y, slope, yint, corr)
! ******************************************************************************
    ! -- arguments
    real(R8B), dimension(:), intent(in) :: x
    real(R8B), dimension(:), intent(in) :: y
    real(R8B), intent(out) :: slope
    real(R8B), intent(out) :: yint
    real(R8B), intent(out) :: corr
    !
    ! -- locals
    integer(I4B) :: i, n
    real(R8B) :: r8n, sumx, sumx2, sumxy, sumy, sumy2
! ------------------------------------------------------------------------------
    sumx  = DZERO  ! sum of x
    sumx2 = DZERO  ! sum of x**2
    sumxy = DZERO  ! sum of x * y
    sumy  = DZERO  ! sum of y
    sumy2 = DZERO  ! sum of y**2
    
    n = size(x); r8n = real(n,R8B)
    if (n /= size(y)) then
      call errmsg('linear_regression: size of x and y do not match.')
    end if
    !
    do i = 1, n
      sumx  = sumx  + x(i)
      sumx2 = sumx2 + x(i)*x(i)
      sumxy = sumxy + x(i)*y(i)
      sumy  = sumy  + y(i)
      sumy2 = sumy2 + y(i)*y(i)
    end do
    
    slope = (r8n  * sumxy -  sumx * sumy)  / (r8n * sumx2 - sumx**2) ! slope
    yint  = (sumy * sumx2 -  sumx * sumxy) / (r8n * sumx2 - sumx**2) ! y-intercept
    corr  = (sumxy - sumx * sumy / r8n) / &                          ! correlation coefficient
            sqrt((sumx2 - sumx**2/r8n) * (sumy2 - sumy**2/r8n))
    !
    return
  end subroutine linear_regression
  
  recursive subroutine label_node(ia, ja, id1, i4wk1d, ireg)
! ******************************************************************************
    ! -- arguments
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: ja
    integer(I4B), intent(in) :: id1
    integer(I4B), dimension(:), intent(inout) :: i4wk1d
    integer(I4B), intent(in) :: ireg
    !
    ! -- locals
    integer(I4B) :: i, id2
! ------------------------------------------------------------------------------
    !
    i4wk1d(id1) = ireg
    !
    do i = ia(id1)+1, ia(id1+1)-1
      id2 = ja(i)
      if (i4wk1d(id2) == 0) then
        call label_node(ia, ja, id2, i4wk1d, ireg)
      end if
    end do
    !
    return
  end subroutine label_node

  subroutine get_r4grid_bb(a, mv, id, bb)
! ******************************************************************************
    ! -- arguments
    real(R4B), dimension(:,:), intent(in) :: a
    real(R4B), intent(in) :: mv
    integer(I4B), dimension(:), allocatable, intent(out) :: id
    type(tBb), dimension(:), allocatable, intent(out) :: bb
    !
    ! -- locals
    integer(I4B) :: nr, nc, ir, ic, i, j, n, mxid, nid, i4v
    !
    integer(I4B), dimension(:), allocatable :: i4wk1d
! ------------------------------------------------------------------------------
    nc = size(a,1); nr = size(a,2)
    !
    mxid = 0
    do ir = 1, nc
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          i4v = int(a(ic,ir),I4B)
          mxid = max(mxid,i4v)
        end if
      end do
    end do
    !
    if (mxid <= 0) call errmsg('get_r4grid_bb: error 1')
    allocate(i4wk1d(mxid))
    do i = 1, mxid
      i4wk1d(i) = 0
    end do
    !
    do ir = 1, nc
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          i4v = int(a(ic,ir),I4B)
          i4wk1d(i4v) = 1
        end if
      end do
    end do
    !
    nid = 0
    do i = 1, mxid
      if (i4wk1d(i) == 1) then
        nid = nid + 1
        i4wk1d(i) = nid
      end if
    end do
    !
    if (nid == 0) call errmsg('get_r4grid_bb: error 2')
    allocate(id(nid), bb(nid))
    !
    do i = 1, mxid
      j = i4wk1d(i)
      if (j > 0) then
        id(j) = i
      end if
    end do
    !
    do ir = 1, nc
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          i4v = int(a(ic,ir),I4B)
          i = i4wk1d(i4v)
          bb(i)%ic0 = min(bb(i)%ic0,ic); bb(i)%ic1 = max(bb(i)%ic1,ic)
          bb(i)%ir0 = min(bb(i)%ir0,ir); bb(i)%ir1 = max(bb(i)%ir1,ir)
        end if
      end do
    end do
    do i = 1, nid
      bb(i)%ncol = bb(i)%ic1-bb(i)%ic0+1
      bb(i)%nrow = bb(i)%ir1-bb(i)%ir0+1
    end do
    !
    return
  end subroutine get_r4grid_bb
  
  subroutine get_bb_extent_r8(xll, yll, cs, ncol, nrow, xmin, xmax, ymin, ymax, &
    ic0, ic1, ir0, ir1)
! ******************************************************************************
    ! -- arguments
    real(R8B), intent(in) :: xll, yll, cs
    integer(I4B), intent(in) :: ncol, nrow
    real(R8B), intent(in) :: xmin, xmax, ymin, ymax
    integer(I4B), intent(out) :: ic0, ic1, ir0, ir1
    ! -- locals
    real(R8B) :: yul
! ------------------------------------------------------------------------------
    yul = yll + nrow*cs
    ic0 = int((xmin-xll)/cs)+1; ic1 = int((xmax-xll)/cs)
    ir0 = int((yul-ymax)/cs)+1; ir1 = int((yul-ymin)/cs)
    ic0 = max(ic0,1); ic0 = min(ic0,ncol)
    ic1 = max(ic1,1); ic1 = min(ic1,ncol)
    ir0 = max(ir0,1); ir0 = min(ir0,nrow)
    ir1 = max(ir1,1); ir1 = min(ir1,nrow)
    !
    return
  end subroutine get_bb_extent_r8
  
  subroutine grid_to_graph(a, mv, ia, ja, id, bb, label_bnd_id)
! ******************************************************************************
    ! -- arguments
    real(R4B),    dimension(:,:),              intent(in)  :: a
    real(R4B),                                 intent(in)  :: mv
    integer(I4B), dimension(:),   allocatable, intent(out) :: ia
    integer(I4B), dimension(:),   allocatable, intent(out) :: ja
    integer(I4B), dimension(:),   allocatable, intent(out) :: id
    type(tBb),    dimension(:),   allocatable, intent(out) :: bb
    logical, optional                                      :: label_bnd_id
    !
    ! -- locals
    integer(I4B), parameter :: jp = 1
    integer(I4B), parameter :: jn = 2
    integer(I4B), parameter :: js = 3
    integer(I4B), parameter :: je = 4
    integer(I4B), parameter :: jw = 5
    integer(I4B), parameter :: ns = jw
    !
    integer(I4B), dimension(2,ns) :: st
    data st/ 0,  0, &
             0, -1, &
             0,  1, &
             1,  0, &
            -1,  0/
    !
    logical :: lidbnd, lbnd
    !
    integer(I4B) :: i, j, nc, nr, ic, ir, jc, jr, mxid, i4v, nid, nbr, nja
    integer(I4B) :: ir0, ir1, ic0, ic1, iact
    !
    integer(I4B), dimension(:),   allocatable :: i4wk1d
    integer(I4B), dimension(:),   allocatable :: i4wk1d2
    integer(I4B), dimension(:,:), allocatable :: i4wk2d
! ------------------------------------------------------------------------------
    !
    ! clean
    if (allocated(ia)) deallocate(ia)
    if (allocated(ja)) deallocate(ja)
    if (allocated(id)) deallocate(id)
    if (allocated(bb)) deallocate(bb)
    if (present(label_bnd_id)) then
      lidbnd = label_bnd_id
    else
      lidbnd = .false.
    end if
    !
    nc = size(a,1); nr = size(a,2)
    !
    ! determine the maximum
    mxid = 0
    do ir = 1, nr
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          i4v = int(a(ic,ir),I4B)
          mxid = max(mxid,i4v)
        end if
      end do
    end do
    !
    ! label the ids
    allocate(i4wk1d(mxid))
    do i = 1, mxid
      i4wk1d(i) = 0
    end do
    do ir = 1, nr
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          i4v = int(a(ic,ir),I4B)
          i4wk1d(i4v) = 1
        end if
      end do
    end do
    !
    nid = 0
    do i = 1, mxid
      if (i4wk1d(i) == 1) then
        nid = nid + 1
        i4wk1d(i) = nid
       end if
    end do
    !
    allocate(id(nid))
    do i = 1, mxid
      j = i4wk1d(i)
      if (j > 0) then
        id(j) = i
      end if
    end do
    !
    ! create local matrix; set bb
    allocate(i4wk2d(nc,nr))
    allocate(bb(nid))
    !
    do ir = 1, nr
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          i4v = int(a(ic,ir),I4B)
          i = i4wk1d(i4v) ! local id
          i4wk2d(ic,ir) = i
          bb(i)%ic0 = min(bb(i)%ic0,ic); bb(i)%ic1 = max(bb(i)%ic1,ic)
          bb(i)%ir0 = min(bb(i)%ir0,ir); bb(i)%ir1 = max(bb(i)%ir1,ir)
        else
          i4wk2d(ic,ir) = 0
        end if
      end do
    end do
    do i = 1, nid
      bb(i)%ncol = bb(i)%ic1-bb(i)%ic0+1
      bb(i)%nrow = bb(i)%ir1-bb(i)%ir0+1
    end do
    !
    ! determine the number of neighbors
    allocate(i4wk1d2(nid))
    do iact = 1, 2
      nja = 0
      do i = 1, nid
        do j = 1, nid
          i4wk1d2(j) = 0
        end do
        lbnd = .false.
        do ir = max(1,bb(i)%ir0-1), min(nr,bb(i)%ir1+1)
          do ic = max(1,bb(i)%ic0-1), min(nc,bb(i)%ic1+1)
           if (i4wk2d(ic,ir) == i) then
             do j = 2, ns
               jc = ic+st(1,j); jc = max(jc,1); jc = min(jc,nc)
               jr = ir+st(2,j); jr = max(jr,1); jr = min(jr,nr)
               i4v = i4wk2d(jc,jr)
               if (i4v > 0) then
                 if (i4v /= i) then
                   i4wk1d2(i4v) = 1 ! neighbor found
                end if
               else
                 lbnd = .true.
               end if
             end do
           end if
          end do
        end do
        !
        ! label catchment
        if (lbnd.and.lidbnd) then
          id(i) = -abs(id(i))
        end if
        !
        nja = nja + 1
        if (iact == 2) then ! center id
          ja(nja) = i
        end if
        nbr = 0
        do j = 1, nid
          if (i4wk1d2(j) == 1) then
            nbr = nbr + 1; nja = nja + 1
            if (iact == 2) then ! neighbor id
              ja(nja) = j
            end if
          end if
        end do
        if (iact == 2) then
          ia(i+1) = ia(i) + nbr + 1
        end if
      end do
      if (iact == 1) then
        allocate(ia(nid+1), ja(nja))
        ia(1) = 1
      end if
    end do !act
    !
    ! clean up
    deallocate(i4wk1d, i4wk1d2, i4wk2d)
    !
    return
  end subroutine grid_to_graph
  
  recursive subroutine balance_graph(ia, ja, id1, lev1, lev, n)
! ******************************************************************************
    ! -- arguments
    integer(I4B), dimension(:), intent(in)    :: ia
    integer(I4B), dimension(:), intent(in)    :: ja
    integer(I4B),               intent(in)    :: id1
    integer(I4B),               intent(in)    :: lev1
    integer(I4B), dimension(:), intent(inout) :: lev
    integer(I4B),               intent(inout) :: n
    !
    ! -- locals
    integer(I4B) :: i, id2, lev2, tlev, dlev
! ------------------------------------------------------------------------------
    !
    lev(id1) = lev1
    !
    do i = ia(id1)+1, ia(id1+1)-1
      id2 = ja(i); lev2 = lev(id2); dlev = lev1-lev2
      if (dlev > 1) then
        tlev = lev2 + 1
        call balance_graph(ia, ja, id2, tlev, lev, n)
      end if
      if (abs(dlev) > 1) then
        n = n + 1
      end if
    end do
    !
    return
  end subroutine balance_graph
  
  function get_jd(y, m, d) result(jd)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: y
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: d
    real(R8B) :: jd
    ! -- locals
    integer(I4B) :: date
    ! -- functions

! ------------------------------------------------------------------------------
    date = y*10000+m*100+d
    call cfn_datehms2mjd(date,0,0,0,jd)
    !
    return
  end function get_jd

 subroutine get_ymd_from_jd(jd, date, y, m, d)
! ******************************************************************************
    ! -- arguments
    real(R8B), intent(in) :: jd
    integer(I4B), intent(out) :: date
    integer(I4B), intent(out) :: y
    integer(I4B), intent(out) :: m
    integer(I4B), intent(out) :: d
    ! -- locals
    integer(I4B) :: idum
    character(len=100) :: s
! ------------------------------------------------------------------------------
    call cfn_mjd2datehms(jd,date,idum,idum,idum)
    write(s,*) date
    s = adjustl(s)
    read(s(1:4),*) y
    read(s(5:6),*) m
    read(s(7:8),*) d
    !
    return
  end subroutine get_ymd_from_jd

  subroutine jd_next_month(jd)
! ******************************************************************************
    ! -- arguments
    real(R8B), intent(inout) :: jd
    ! -- locals
    integer(I4B) :: date, y, m, d
! ------------------------------------------------------------------------------
    call get_ymd_from_jd(jd, date, y, m, d)

    if (m == 12) then
      y = y + 1
      m = 1
    else
      m = m + 1
    end if
    jd = get_jd(y, m, d)
    !
    return
  end subroutine jd_next_month
!  
  function get_month_days_s(date) result(nd)
! ******************************************************************************
    ! -- arguments
    character(len=MXSLEN), intent(in) :: date
    integer(I4B) :: nd
    ! -- locals
    real(R8B) :: jd
    ! -- functions
    integer(I4B) :: y, m, d, ios
! ------------------------------------------------------------------------------
    !
    read(date(1:4),*,iostat=ios) y
    read(date(5:6),*,iostat=ios) m
    read(date(7:8),*,iostat=ios) d
    if (ios /= 0) then
      d = 1
    end if
    jd = get_jd(y, m, d)
    nd = get_month_days(jd)
    !
    return
  end function get_month_days_s
  !
  function get_month_days(jd) result(nd)
! ******************************************************************************
    ! -- arguments
    real(R8B), intent(in) :: jd
    integer(I4B) :: nd
    ! -- locals
    integer(I4B) :: date1, date2
    real(R8B) :: jd1, jd2, djd
    ! -- functions
    integer(I4B) :: y, m, d
    real(R8B) :: cfn_jd_delta
! ------------------------------------------------------------------------------
    call get_ymd_from_jd(jd, date1, y, m, d)
    if (m == 12) then
      date2 = (y+1)*10000+1*100+1
    else
      date2 = y*10000+(m+1)*100+1
    end if
    !
    call cfn_datehms2mjd(date1,0,0,0,jd1)
    call cfn_datehms2mjd(date2,0,0,0,jd2)
    djd = cfn_jd_delta(jd2,jd1)
    !
    nd = nint(djd)
    !
    return
  end function get_month_days

  function readgen(f) result(p)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: f
    type(tPol), dimension(:), allocatable :: p
    ! -- locals
    integer(I4B) :: iu, id, n, m, ios, iact, i
    real(R4B) :: x, y
    character(len=MXSLEN) :: s
! ------------------------------------------------------------------------------
    !
    call open_file(f, iu, 'r')
    !
    do iact = 1, 3
      n = 0
      do while(.true.)
        read(iu,'(a)',iostat=ios) s
        if (ios /= 0) exit
        if (s == 'END') exit
        n = n + 1
        read(s,*) id
        if (iact == 2) p(n)%id = id
        m = 0
        do while(.true.)
          read(iu,'(a)',iostat=ios) s
          if (ios /= 0) exit
          if (s == 'END') exit
          m = m + 1
          if (iact == 2) then
            p(n)%n = p(n)%n + 1
          end if
          if (iact == 3) then
            read(s,*) x, y
            p(n)%xmin = min(x, p(n)%xmin); p(n)%ymin = min(y, p(n)%ymin)
            p(n)%xmax = max(x, p(n)%xmax); p(n)%ymax = max(y, p(n)%ymax)
            p(n)%xy(1,m) = x; p(n)%xy(2,m) = y
          end if
        end do
      end do
      if (iact == 1) then
        rewind iu
        allocate(p(max(1,n)))
      end if
      if (iact == 2) then
        rewind iu
        do i = 1, n
          m = max(1,p(i)%n)
          allocate(p(i)%xy(2,m))
        end do
      end if
    end do
    !
    close(iu)
    !
    return
  end function readgen

  subroutine filtergen_i1(p, i1a, xmin, ymin, cs, i1mv)
! ******************************************************************************
    ! -- arguments
    type(tPol), dimension(:), intent(in) :: p
    integer(I1B), dimension(:,:), intent(inout) :: i1a
    real(R8B), intent(in) :: xmin, ymin, cs
    integer(I1B) :: i1mv
    ! -- locals
    integer(I4B) :: nc, nr, i, ic0, ic1, ir0, ir1, ir, ic
    real(R4B) :: xmax, ymax, x0, x1, y0, y1
    real(R8B), dimension(2) :: point
    logical :: lin
! ------------------------------------------------------------------------------
    !
    write(*,*) 'Filtering for gen-file...'
    !
    nc = size(i1a,1); nr = size(i1a,2)
    xmax = xmin + nc*cs; ymax = ymin + nr*cs
    !
    do i = 1, size(p)
      ! determine bounding boxes
      x0 = p(i)%xmin; x1 = p(i)%xmax
      y0 = p(i)%ymin; y1 = p(i)%ymax
      x0 = max(xmin, x0); y0 = max(ymin, y0)
      x1 = min(xmax, x1); y1 = min(ymax, y1)
      ic0 = (x0 - xmin)/cs; ic0 = max(ic0-1,1)
      ic1 = (x1 - xmin)/cs; ic1 = min(ic1+1,nc)
      ir0 = (ymax - y1)/cs; ir0 = max(ir0-1,1)
      ir1 = (ymax - y0)/cs; ir1 = min(ir1+1,nr)
      !
      do ir = ir0, ir1
        do ic = ic0, ic1
          point(1) = xmin + ic*cs - cs/2
          point(2) = ymax - ir*cs + cs/2
          call polygon_contains_point_2d (p(i)%n, p(i)%xy, point, lin)
          if (lin) then
            i1a(ic,ir) = i1mv
          end if
        end do
      end do
    end do
    !
    return
  end subroutine filtergen_i1

  subroutine writebin_i(f, x, nodata)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: f
    integer(I4B), dimension(:,:), intent(in) :: x
    integer(I4B), intent(in) :: nodata
    ! -- locals
    integer(I4B) :: nc, nr, n, ic, ir, iu
! ------------------------------------------------------------------------------

    call open_file(f, iu, 'w', .true.)
    !
    ! count
    nc = size(x,1); nr = size(x,2); n = 0
    do ir = 1, nr
      do ic = 1, nc
        if (x(ic,ir) /= nodata) then
          n = n + 1
        end if
      end do
    end do
    !
    ! write
    write(iu) n
    do ir = 1, nr
      do ic = 1, nc
        if (x(ic,ir) /= 0) then
          write(iu) ic, ir, x(ic,ir)
        end if
      end do
    end do
    close(iu)
    !
    return
  end subroutine writebin_i

  function tas(s_in) result(s_out)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: s_in
    character(len=:), allocatable :: s_out
    ! -- locals
! ------------------------------------------------------------------------------
    s_out = trim(adjustl(s_in))
    !
    return
  end function tas

  function ta_i1(arr, fmt_in, sep_in) result(s)
! ******************************************************************************
    ! -- arguments
    integer(I1B), dimension(:), intent(in) :: arr
    character(len=:), allocatable :: s
    character(len=*), intent(in), optional :: fmt_in
    character(len=1), intent(in), optional :: sep_in
    ! -- locals
    logical :: lfmt
    integer(I4B) :: i, n
    character(len=MXSLENLONG) :: s_tmp
    character(len=MXSLEN) :: w, fmt
    character(len=1) :: sep
! ------------------------------------------------------------------------------
    lfmt = .false.
    if (present(fmt_in)) then
      fmt = fmt_in
      lfmt = .true.
    end if
    if (lfmt) then
      write(w,fmt) arr(1)
      s_tmp = trim(w)
    else
      write(w,*) arr(1)
      s_tmp = trim(adjustl(w))
    end if
    if (present(sep_in)) then
      sep = sep_in
    else
      sep = ' '
    end if
    do i = 2, size(arr)
      if (lfmt) then
        write(w,fmt) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(w)
      else
        write(w,*) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(adjustl(w))
      end if
    end do
    !
    n = len_trim(s_tmp)
    allocate(character(len=n) :: s)
    s = s_tmp(1:n)
    !
    return
  end function ta_i1
  
  function ta_i2(arr, fmt_in, sep_in) result(s)
! ******************************************************************************
    ! -- arguments
    integer(I2B), dimension(:), intent(in) :: arr
    character(len=:), allocatable :: s
    character(len=*), intent(in), optional :: fmt_in
    character(len=1), intent(in), optional :: sep_in
    ! -- locals
    logical :: lfmt
    integer(I4B) :: i, n
    character(len=MXSLENLONG) :: s_tmp
    character(len=MXSLEN) :: w, fmt
    character(len=1) :: sep
! ------------------------------------------------------------------------------
    lfmt = .false.
    if (present(fmt_in)) then
      fmt = fmt_in
      lfmt = .true.
    end if
    if (lfmt) then
      write(w,fmt) arr(1)
      s_tmp = trim(w)
    else
      write(w,*) arr(1)
      s_tmp = trim(adjustl(w))
    end if
    if (present(sep_in)) then
      sep = sep_in
    else
      sep = ' '
    end if
    do i = 2, size(arr)
      if (lfmt) then
        write(w,fmt) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(w)
      else
        write(w,*) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(adjustl(w))
      end if
    end do
    !
    n = len_trim(s_tmp)
    allocate(character(len=n) :: s)
    s = s_tmp(1:n)
    !
    return
  end function ta_i2
  
  function ta_i4(arr, fmt_in, sep_in) result(s)
! ******************************************************************************
    ! -- arguments
    integer(I4B), dimension(:), intent(in) :: arr
    character(len=:), allocatable :: s
    character(len=*), intent(in), optional :: fmt_in
    character(len=1), intent(in), optional :: sep_in
    ! -- locals
    logical :: lfmt
    integer(I4B) :: i, n
    character(len=MXSLENLONG) :: s_tmp
    character(len=MXSLEN) :: w, fmt
    character(len=1) :: sep
! ------------------------------------------------------------------------------
    lfmt = .false.
    if (present(fmt_in)) then
      fmt = fmt_in
      lfmt = .true.
    end if
    if (lfmt) then
      write(w,fmt) arr(1)
      s_tmp = trim(w)
    else
      write(w,*) arr(1)
      s_tmp = trim(adjustl(w))
    end if
    if (present(sep_in)) then
      sep = sep_in
    else
      sep = ' '
    end if
    do i = 2, size(arr)
      if (lfmt) then
        write(w,fmt) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(w)
      else
        write(w,*) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(adjustl(w))
      end if
    end do
    !
    n = len_trim(s_tmp)
    allocate(character(len=n) :: s)
    s = s_tmp(1:n)
    !
    return
  end function ta_i4
  
  function ta_i8(arr, fmt_in, sep_in) result(s)
! ******************************************************************************
    ! -- arguments
    integer(I8B), dimension(:), intent(in) :: arr
    character(len=:), allocatable :: s
    character(len=*), intent(in), optional :: fmt_in
    character(len=1), intent(in), optional :: sep_in
    ! -- locals
    logical :: lfmt
    integer(I4B) :: i, n
    character(len=MXSLENLONG) :: s_tmp
    character(len=MXSLEN) :: w, fmt
    character(len=1) :: sep
! ------------------------------------------------------------------------------
    lfmt = .false.
    if (present(fmt_in)) then
      fmt = fmt_in
      lfmt = .true.
    end if
    if (lfmt) then
      write(w,fmt) arr(1)
      s_tmp = trim(w)
    else
      write(w,*) arr(1)
      s_tmp = trim(adjustl(w))
    end if
    if (present(sep_in)) then
      sep = sep_in
    else
      sep = ' '
    end if
    do i = 2, size(arr)
      if (lfmt) then
        write(w,fmt) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(w)
      else
        write(w,*) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(adjustl(w))
      end if
    end do
    !
    n = len_trim(s_tmp)
    allocate(character(len=n) :: s)
    s = s_tmp(1:n)
    !
    return
  end function ta_i8
  
  function ta_r4(arr, fmt_in, sep_in) result(s)
! ******************************************************************************
    ! -- arguments
    real(R4B), dimension(:), intent(in) :: arr
    character(len=:), allocatable :: s
    character(len=*), intent(in), optional :: fmt_in
    character(len=1), intent(in), optional :: sep_in
    ! -- locals
    logical :: lfmt
    integer(I4B) :: i, n
    character(len=MXSLENLONG) :: s_tmp
    character(len=MXSLEN) :: w, fmt
    character(len=1) :: sep
! ------------------------------------------------------------------------------
    lfmt = .false.
    if (present(fmt_in)) then
      fmt = fmt_in
      lfmt = .true.
    end if
    if (lfmt) then
      write(w,fmt) arr(1)
      s_tmp = trim(w)
    else
      write(w,*) arr(1)
      s_tmp = trim(adjustl(w))
    end if
    if (present(sep_in)) then
      sep = sep_in
    else
      sep = ' '
    end if
    do i = 2, size(arr)
      if (lfmt) then
        write(w,fmt) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(w)
      else
        write(w,*) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(adjustl(w))
      end if
    end do
    !
    n = len_trim(s_tmp)
    allocate(character(len=n) :: s)
    s = s_tmp(1:n)
    !
    return
  end function ta_r4
  
  function ta_r8(arr, fmt_in, sep_in) result(s)
! ******************************************************************************
    ! -- arguments
    real(R8B), dimension(:), intent(in) :: arr
    character(len=:), allocatable :: s
    character(len=*), intent(in), optional :: fmt_in
    character(len=1), intent(in), optional :: sep_in
    ! -- locals
    logical :: lfmt
    integer(I4B) :: i, n
    character(len=MXSLENLONG) :: s_tmp
    character(len=MXSLEN) :: w, fmt
    character(len=1) :: sep
! ------------------------------------------------------------------------------
    lfmt = .false.
    if (present(fmt_in)) then
      fmt = fmt_in
      lfmt = .true.
    end if
    if (lfmt) then
      write(w,fmt) arr(1)
      s_tmp = trim(w)
    else
      write(w,*) arr(1)
      s_tmp = trim(adjustl(w))
    end if
    if (present(sep_in)) then
      sep = sep_in
    else
      sep = ' '
    end if
    do i = 2, size(arr)
      if (lfmt) then
        write(w,fmt) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(w)
      else
        write(w,*) arr(i)
        s_tmp = trim(s_tmp)//sep//trim(adjustl(w))
      end if
    end do
    !
    n = len_trim(s_tmp)
    allocate(character(len=n) :: s)
    s = s_tmp(1:n)
    !
    return
  end function ta_r8
  
  function ta_c(arr, add_quotes, sep_in) result(s)
! ******************************************************************************
    ! -- arguments
    character(len=*), dimension(:), intent(in) :: arr
    logical, intent(in), optional :: add_quotes
    character(len=1), intent(in), optional :: sep_in
    character(len=:), allocatable :: s
    ! -- locals
    integer(I4B) :: i
    character(len=1) :: q, sep
    character(len=MXSLEN) :: w
! ------------------------------------------------------------------------------
    if (present(add_quotes)) then
      q = '"'
    else
      q =''
    end if
    if (present(sep_in)) then
      sep = sep_in
    else
      sep = ' '
    end if
    !
    w = arr(1)
    s = trim(q)//trim(adjustl(w))//trim(q)
    do i = 2, size(arr)
      w = arr(i)
      s = s//sep//trim(q)//trim(adjustl(w))//trim(q)
    end do
    !
    return
  end function ta_c

  function bb_overlap(bb1, bb2) result(loverlap)
! ******************************************************************************
    ! -- arguments
    type(tBB), pointer, intent(in) :: bb1
    type(tBB), pointer, intent(in) :: bb2
    logical :: loverlap
    ! -- locals
! ------------------------------------------------------------------------------
    loverlap = .true.
    if ((bb1%ic0 >= bb2%ic1).or.(bb2%ic0 >= bb1%ic1)) then
      loverlap = .false.
    end if
    if ((bb1%ir0 >= bb2%ir1).or.(bb2%ir0 >= bb1%ir1)) then
      loverlap = .false.
    end if
    !
    return
  end function bb_overlap
  
  subroutine swap_slash(s)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: s
    ! -- locals
    character(len=1) :: src_slash, tgt_slash
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    if (os == 1) then ! windows
      src_slash = lin_slash
      tgt_slash = win_slash
    else ! linux
      src_slash = win_slash
      tgt_slash = lin_slash
    end if
    !
    do i = 1, len_trim(s)
      if (s(i:i) == src_slash) then
        s(i:i) = tgt_slash
      end if
    end do
    !
    return
  end subroutine swap_slash
  !
  function get_slash() result(slash)
! ******************************************************************************
    ! -- arguments
    character(len=1) :: slash
    ! -- locals
! ------------------------------------------------------------------------------
    if (os == 1) then ! windows
      slash = win_slash
    else ! linux
      slash = lin_slash
    end if
    !
    return
  end function get_slash

  function count_dir_files(d) result(n)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: d
    integer(I4B) :: n
    ! -- locals
    character(len=MXSLEN) :: tf, s
    integer(I4B) :: iu, ios
! ------------------------------------------------------------------------------
    tf = 'tmp.txt'
    call swap_slash(d)
    if (os == 1) then ! windows
      call system('dir /b '//trim(d)//' > '//trim(tf))
    else
      call system('ls -1 '//trim(d)//' > '//trim(tf))
    end if
    !
    n = 0
    call open_file(tf, iu, 'r')
    do while(.true.)
      read(unit=iu,fmt='(a)',iostat=ios) s
      if (ios /= 0) exit
      n = n + 1
    end do
    close(iu,status='delete')
    !
    return
  end function count_dir_files
  
  function get_dir_files(d) result(f)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: d
    character(len=MXSLEN), dimension(:), allocatable :: f ! result
    ! -- locals
    character(len=MXSLEN) :: tf, s
    integer(I4B) :: iu, ios, n, iact
! ------------------------------------------------------------------------------
    tf = 'tmp.txt'
    call swap_slash(d)
    if (os == 1) then ! windows
      call system('dir /b/s '//trim(d)//' > '//trim(tf))
    else
      call errmsg('get_dir_files: linux not yet supported.')
    end if
    !
    call open_file(tf, iu, 'r')
    do iact = 1, 2
      n = 0
      do while(.true.)
        read(unit=iu,fmt='(a)',iostat=ios) s
        if (ios /= 0) exit
        n = n + 1
        if (iact == 2) then
          f(n) = trim(adjustl(s))
        end if
      end do
      if (iact == 1) then
        allocate(f(n))
        rewind(iu)
      end if
    end do
    !
    close(iu,status='delete')
    !
    return
  end function get_dir_files
  
  subroutine get_rel_up(f, n)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: f
    integer(I4B), intent(in) :: n
    ! -- locals
    integer(I4B) :: i, j
    character(len=1) :: slash
! ------------------------------------------------------------------------------
    call swap_slash(f)
    slash = get_slash()
    do i = 1, n
      j = index(f, slash)
      if (j < 0) then
        call errmsg('Could not determine relative path')
      end if
      f = f(j+1:)
    end do
    f = '.'//slash//trim(f)
    !
    return
  end subroutine get_rel_up

  subroutine create_dir(d, lverb_in)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: d
    logical, intent(in), optional :: lverb_in
    ! -- locals
    logical :: ldirexist, lverb
    integer(I4B) :: ios
! ------------------------------------------------------------------------------
    if (present(lverb_in)) then
      lverb = lverb_in
    else
      lverb = .false.
    end if
    !
    call swap_slash(d)
    inquire(directory=d, exist=ldirexist, iostat=ios)
    if (ios.ne.0) ldirexist=.false.
    if (ldirexist) then
      if (.not.lverb) then
        call logmsg('Directory '//trim(d)//' already already exists.')
      end if
      return
    end if
    !
    if (.not.lverb) then
      call logmsg('Creating directory '//trim(d)//'.')
    end if
    if (os == 1) then !windows
      call system('mkdir '//trim(d))
    end if
    if (os == 2) then !linux
      call system('mkdir -p '//trim(d))
    end if
    !
    return
  end subroutine create_dir
  !
  subroutine change_work_dir(d, lverb_in)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: d
    logical, intent(in), optional :: lverb_in
    ! -- locals
    character(len=MXSLEN) :: cd
    logical :: lverb, ldirexist
    integer(I4B) :: ios
! ------------------------------------------------------------------------------
    if (present(lverb_in)) then
      lverb = lverb_in
    else
      lverb = .false.
    end if
    !
    call swap_slash(d)
    inquire(directory=d, exist=ldirexist, iostat=ios)
    if (ios.ne.0) ldirexist=.false.
    if (.not.ldirexist) then
      call errmsg('Directory '//trim(d)//' does not exist.')
    end if
    !
    if (.not.lverb) then
      call get_work_dir(cd)
      call logmsg('Changing working directory '//trim(cd)//' -> '//trim(d)//'.')
    end if
    call chdir(d)
    !
    return
  end subroutine change_work_dir
  !
  subroutine get_work_dir(d)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: d
! ------------------------------------------------------------------------------
    call getcwd(d)
    !
    return
  end subroutine get_work_dir
  !
  subroutine get_abs_path(f)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: f
    ! -- locals
    character(len=1) :: slash
    character(len=MXSLEN) :: d, s
    integer(I4B) :: i, j, n
! ------------------------------------------------------------------------------
    !
    call swap_slash(f)
    if (f(1:1) /= '.') return
    !
    call get_work_dir(d)
    !
    s = f
    !
    slash = get_slash()
    !
    n = 0
    if (s(1:3) == '..'//slash) then
      n = 1; s = s(4:)
      do while (s(1:3) == '..'//slash)
        n = n + 1; s = s(4:)
      end do
    end if
    !
    if (n == 0) then 
      f = trim(d)//trim(s(3:))
    else
      do i = 1, n
        j = index(d, slash, back=.true.)
        if (j <= 0) then
          call errmsg('get_abs_path: '//trim(f)//'.')
        end if
        d = d(1:j-1)
      end do
      f = trim(d)//slash//trim(s)
    end if
    !
    return
  end subroutine
  !
  subroutine open_file(f, iu, act_in, lbin_in, pos_in)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: f
    integer(I4B), intent(inout) :: iu
    character(len=1), optional :: act_in
    logical, intent(in), optional :: lbin_in
    character(len=*), intent(in), optional :: pos_in
    ! -- locals
    character(len=1) :: act
    character(len=MXSLEN) :: pos, stat
    logical :: lbin, lex, lop, lverb
! ------------------------------------------------------------------------------
    pos = 'asis'
    if (present(pos_in)) then
      pos = pos_in
    end if
    if ((trim(pos) /= 'asis').and. &
        (trim(pos) /= 'rewind').and. &
        (trim(pos) /= 'append')) then
      call errmsg('open_file: invalid position specifier.')
    end if
    lverb = .false.
    if (present(act_in)) then
      if ((act_in == 'R').or.(act_in == 'W')) then
        lverb = .true.
      end if
      act = change_case(act_in, 'l')
    else
      act = 'r'
    end if
    if (present(lbin_in)) then
      lbin = lbin_in
    else
      lbin= .false.
    end if
    !
    call swap_slash(f)
    !
    inquire(file=f, exist=lex)
    if (lex) then
      inquire(file=f, opened=lop, number=iu)
      if (lop) then
        call logmsg('Warning: file '//trim(f)//' is already opened.')
        return
      end if
    end if
    !
    if (act == 'r') then
      inquire(file=f, exist=lex)
      if (.not.lex) then
        call errmsg('File '//trim(f)//' does not exist.')
      end if
    end if

    iu = getlun()
    if ((act == 'r') .and.(.not.lbin)) then
      if (.not.lverb) then
        call logmsg('Reading ascii file '//trim(f)//'...')
      end if
      open(unit=iu, file=f, form='formatted', access='sequential', &
        action='read', status='old',share='denynone', position=pos)
    else if ((act == 'w') .and.(.not.lbin)) then
      if (.not.lverb) then
        call logmsg('Writing ascii file '//trim(f)//'...')
      end if
      open(unit=iu, file=f, form='formatted', access='sequential', &
        action='write', status='replace', position=pos)
    else if ((act == 'r') .and.(lbin)) then
      if (.not.lverb) then
        call logmsg('Reading binary file '//trim(f)//'...')
      end if
      open(unit=iu, file=f, form='unformatted', access='stream', &
        action='readwrite', status='old',share='denynone', position=pos)
    else if ((act == 'w') .and.(lbin)) then
      if (.not.lverb) then        
        call logmsg('Writing binary file '//trim(f)//'...')
      end if
      if (pos == 'append') then
        stat = 'old'
      else
        stat = 'replace'
      end if
      open(unit=iu, file=f, form='unformatted', access='stream', &
        action='write', status=stat, position=pos)
    else
      call errmsg('Subroutine open_file called with invalid option')
    end if
    !
    return
  end subroutine open_file

  subroutine checkdim(nc1, nr1, nc2, nr2)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: nc1, nr1, nc2, nr2
! ------------------------------------------------------------------------------
    if (nc1 /= nc2) then
      call errmsg('Inconsistent number of columns.')
    end if
    if (nr1 /= nr2) then
      call errmsg('Inconsistent number of rows.')
    end if
    !
    return
  end subroutine checkdim

  function readline(iu, so) result(ios)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: iu
    character(len=*), intent(out), optional :: so
    integer(I4B) :: ios
    ! -- locals
    character(len=MXSLEN) :: s
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    do while(.true.)
      read(unit=iu, iostat=ios, fmt='(a)') s
      if (ios /= 0) exit
      so = trim(adjustl(s))
      if ((so(1:1) /= comment) .and. (len_trim(so) > 0)) then
        i = index(so, comment, back=.true.)
        if (i > 0) then
          so = so(1:i-1)
        end if
        exit
      end if
    end do
    !
    return
  end function readline

  function createtoken() result(t)
! ******************************************************************************
    ! -- arguments
    character(len=1) :: t
! ------------------------------------------------------------------------------
    t = '?'
    !
    return
  end function createtoken

  function counttoken(s, t) result(n)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: s
    character(len=1), optional, intent(in) :: t
    integer(I4B) :: n
    ! -- locals
    character(len=1) :: tloc
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    tloc = createtoken()
    if (present(t)) then
      tloc = t
    end if
    !
    n = 0
    do i = 1, len(s)
      if (s(i:i) == tloc) then
        n = n + 1
      endif
    enddo
    !
    return
  end function counttoken

  subroutine replacetoken(s, t, i)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(inout) :: s
    character(len=1), intent(in) :: t
    integer(I4B) :: i
    ! -- locals
    integer(I4B) :: i1, i2, n
    character(len=MXSLEN) :: ns, is, fmt
! ------------------------------------------------------------------------------
    i1 = index(s, t)
    i2 = index(s, t, back=.true.)
    n = i2 - i1 + 1
    write(ns,*) n
    fmt = '(i'//trim(ns)//'.'//trim(ns)//')'
    write(is,fmt) i
    is = adjustl(is)
    s(i1:i2) = is(1:n)
    !
    return
  end subroutine

  subroutine getminmax(key, sep, token, imin, imax)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: key
    character(len=1), intent(in) :: sep
    character(len=*), intent(in) :: token
    integer(I4B), intent(out) :: imin
    integer(I4B), intent(out) :: imax
    ! -- locals
    character(len=MXSLEN) :: s
    integer(I4B) :: i, j, n, ios1, ios2, ival1, ival2
    character(len=MXSLEN), dimension(:), allocatable :: words
! ------------------------------------------------------------------------------
    words = getwords(key, sep)
    n = size(words)
    if (n == 0) return
    !
    imin = 0
    imax = 0
    do i = 1, n
      s = words(i)
      if (s(1:1) == token) then
        j = index(s,':')
        if (j == 0) then
          read(s(2:),*,iostat=ios1) ival1
          if (ios1 == 0) then
            imin = ival1
            imax = imin
          else
            !write(*,*) '@@@ debug'
          end if
        else
          read(s(2:j-1),*,iostat=ios1) ival1
          read(s(j+1:),*,iostat=ios2) ival2
          if ((ios1 == 0).and.(ios2 == 0)) then
            imin = ival1
            imax = ival2
          else
            call errmsg("Could not read "//trim(s))
          end if
        end if
      end if
    end do
    !
    return
  end subroutine getminmax

  function getwords(s_in, token) result(words)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in), optional :: s_in
    character(len=MXSLEN), dimension(:), allocatable :: words
    character(len=1), optional, intent(in) :: token
    ! -- locals
    integer(I4B) :: i, j, n, i0, i1
    character(len=1) :: tokenLocal
    character(len=1), parameter :: quote = '"'
    character(len=MXSLEN) :: s, s1, s2
    integer(I4B), parameter :: maxwords = 100
    integer(I4B), dimension(maxwords) :: ind
    logical :: lquote
! ------------------------------------------------------------------------------
    !
    s = s_in
    !
    tokenLocal = ' '
    if (present(token)) then
      tokenLocal = token
    endif
    !
    ! first check for quotes
    lquote = .false.
    i0 = index(s_in,quote)
    if (i0 > 0) then
      i1 = index(s,quote,back=.true.)
      if ((i1 > 0).and.(i0 /= i1)) lquote = .true.
    end if
    if (lquote) then
      do i = i0, i1
        if (s(i:i) == tokenLocal) s(i:i) = quote
      end do
    end if
    !
    i = index(s_in,comment)
    if (i > 0) then
      !s1 = trim(s_in(1:i))
      s1 = trim(s(1:i))
    else
      !s1 = trim(s_in)
      s1 = trim(s)
    endif
    s1 = adjustl(s1)
    if (len_trim(s1) == n) then
       call errmsg('getwords.')
    end if
    if (len_trim(s1) == 1) then
      n = 1
      allocate(words(1))
      words(1) = s1
      return
    endif
    !
    ! count
    ind(1) = 1
    n = 1
    do i = 2, len_trim(s1)
      if ((s1(i:i) == tokenLocal).and. (s1(i-1:i-1) /= tokenLocal)) then
        n = n + 1
        ind(n) = i
      end if
    end do
    ind(n+1) = len_trim(s1)+1
    if (n > 0) then
      if (allocated(words)) deallocate(words)
      allocate(words(n))
    endif
    do i = 1, n
      read(s1(ind(i):ind(i+1)-1),'(a)') s2
      j = index(trim(s2), tokenLocal, back=.true.)+1
!      words(i) = s2(j:)
      s2 = s2(j:)
      if (lquote) then
        do j = 1, len_trim(s2)
          if (s2(j:j) == quote) s2(j:j) = ' '
          s2 = adjustl(s2)
        end do
      end if
      words(i) = s2
    end do
    !
    return
  end function getwords

  subroutine writetofile_i4(lun, arr, pre, post)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: lun
    integer(I4B), dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: pre
    character(len=*), intent(in), optional :: post
    ! -- locals
    integer(I4B) :: i, n
    character(len=MXSLEN) :: fmt
! ------------------------------------------------------------------------------
    n = size(arr)
    do i = 1, n
      write(sa(i),*) arr(i)
    end do
    !
    if ((.not.present(pre)).and.(.not.present(post))) then
      write(fmt,'(a,i,a)') '(', n-1, '(a,1x),a)'
      write(lun,fmt) (tas(sa(i)),i=1,n)
    end if
    if ((present(pre)).and.(.not.present(post))) then
      write(fmt,'(a,i,a)') '(a,1x', n-1, '(a,1x),a)'
      write(lun,fmt) tas(pre), (tas(sa(i)),i=1,n)
    end if
    if ((.not.present(pre)).and.(present(post))) then
      write(fmt,'(a,i,a)') '(', n-1, '(a,1x),a,1x,a)'
      write(lun,fmt) (tas(sa(i)),i=1,n), tas(post)
    end if
    if ((present(pre)).and.(present(post))) then
      write(fmt,'(a,i,a)') '(a,1x,', n-1, '(a,1x),a,1x,a)'
      write(lun,fmt) tas(pre), (tas(sa(i)),i=1,n), &
        trim(adjustl(post))
    end if
    flush(lun)
    !
    return
  end subroutine writetofile_i4

   subroutine writetofile_r4(lun, arr, pre, post)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: lun
    real(R4B), dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: pre
    character(len=*), intent(in), optional :: post
    ! -- locals
    integer(I4B) :: i, n
    character(len=MXSLEN) :: fmt
! ------------------------------------------------------------------------------
    n = size(arr)
    do i = 1, n
      write(sa(i),*) arr(i)
    end do
    !
    if ((.not.present(pre)).and.(.not.present(post))) then
      write(fmt,'(a,i,a)') '(', n-1, '(a,1x),a)'
      write(lun,fmt) (trim(adjustl(sa(i))),i=1,n)
    end if
    if ((present(pre)).and.(.not.present(post))) then
      write(fmt,'(a,i,a)') '(a,1x', n-1, '(a,1x),a)'
      write(lun,fmt) trim(adjustl(pre)), (trim(adjustl(sa(i))),i=1,n)
    end if
    if ((.not.present(pre)).and.(present(post))) then
      write(fmt,'(a,i,a)') '(', n-1, '(a,1x),a,1x,a)'
      write(lun,fmt) (trim(adjustl(sa(i))),i=1,n),  trim(adjustl(post))
    end if
    if ((present(pre)).and.(present(post))) then
      write(fmt,'(a,i,a)') '(a,1x,', n-1, '(a,1x),a,1x,a)'
      write(lun,fmt) trim(adjustl(pre)), (trim(adjustl(sa(i))),i=1,n), &
        trim(adjustl(post))
    end if
    flush(lun)
    !
    return
  end subroutine writetofile_r4

   subroutine writetofile_r8(lun, arr, pre, post)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: lun
    real(R8B), dimension(:), intent(in) :: arr
    character(len=*), intent(in), optional :: pre
    character(len=*), intent(in), optional :: post
    ! -- locals
    integer(I4B) :: i, n
    character(len=MXSLEN) :: fmt
! ------------------------------------------------------------------------------
    n = size(arr)
    do i = 1, n
      write(sa(i),*) arr(i)
    end do
    !
    if ((.not.present(pre)).and.(.not.present(post))) then
      write(fmt,'(a,i,a)') '(', n-1, '(a,1x),a)'
      write(lun,fmt) (trim(adjustl(sa(i))),i=1,n)
    end if
    if ((present(pre)).and.(.not.present(post))) then
      write(fmt,'(a,i,a)') '(a,1x', n-1, '(a,1x),a)'
      write(lun,fmt) trim(adjustl(pre)), (trim(adjustl(sa(i))),i=1,n)
    end if
    if ((.not.present(pre)).and.(present(post))) then
      write(fmt,'(a,i,a)') '(', n-1, '(a,1x),a,1x,a)'
      write(lun,fmt) (trim(adjustl(sa(i))),i=1,n),  trim(adjustl(post))
    end if
    if ((present(pre)).and.(present(post))) then
      write(fmt,'(a,i,a)') '(a,1x,', n-1, '(a,1x),a,1x,a)'
      write(lun,fmt) trim(adjustl(pre)), (trim(adjustl(sa(i))),i=1,n), &
        trim(adjustl(post))
    end if
    flush(lun)
    !
    return
  end subroutine writetofile_r8

  subroutine setuniweight(ncol, nrow, rwrk, nodata, iwrk)
! ******************************************************************************
    ! -- arguments
    integer, intent(in) :: ncol, nrow
    real, dimension(ncol,nrow), intent(in) :: rwrk
    real, intent(in) :: nodata
    integer, dimension(:,:), allocatable, intent(out) :: iwrk
    ! -- locals
    integer :: irow, icol
! ------------------------------------------------------------------------------
    if (.not.allocated(iwrk)) then
      allocate(iwrk(ncol,nrow))
    end if

    do irow = 1, nrow
      do icol = 1, ncol
        if (rwrk(icol,irow) /= nodata) then
          iwrk(icol,irow) = 1
        else
          iwrk(icol,irow) = 0
        end if
      end do
    end do
    !
    return
  end subroutine setuniweight

  function getlun() result(lun)
! ******************************************************************************
    ! -- arguments
    integer :: lun
    ! -- locals
    logical :: lex
! ------------------------------------------------------------------------------
    do lun=20, 5000
      inquire(unit=lun,opened=lex)
      if(.not.lex)exit
    end do
    !
    return
  end function getlun

  function fileexist(fname) result(lex)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: fname
    logical :: lex
! ------------------------------------------------------------------------------
    inquire(file=fname,exist=lex)
    !
    return
  end function fileexist
  
  subroutine chkexist(fname)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: fname
    ! -- locals
    logical :: lex
! ------------------------------------------------------------------------------
    inquire(file=fname,exist=lex)
    if (.not.lex) then
      call errmsg('Cannot find '//trim(fname))
    end if
    !
    return
  end subroutine chkexist

  function getdigits(n) result(ndig)
! ******************************************************************************
    ! -- arguments
    integer(kind=8), intent(in) :: n
    integer :: ndig
! ------------------------------------------------------------------------------
    ndig = 1
    if (abs(n) > 10) ndig = 2
    if (abs(n) > 100) ndig = 3
    if (abs(n) > 1000) ndig = 4
    if (n < 0) ndig = ndig + 1
    !
    return
  end function

  subroutine writeascheader(lun, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- arguments
    integer, intent(in) :: lun, ncol, nrow
    double precision, intent(in) :: xll, yll, cs, nodata
    ! -- locals
    character(len=MXSLEN) :: s
! ------------------------------------------------------------------------------
    !
    write(s,*) ncol; write(lun,'(a)') 'ncols '//trim(adjustl(s))
    write(s,*) nrow; write(lun,'(a)') 'nrows '//trim(adjustl(s))
    write(s,*) xll; write(lun,'(a)') 'xllcorner '//trim(adjustl(s))
    write(s,*) yll; write(lun,'(a)') 'yllcorner '//trim(adjustl(s))
    write(s,*) cs; write(lun,'(a)') 'cellsize '//trim(adjustl(s))
    write(s,*) nodata; write(lun,'(a)') 'nodata_value '//trim(adjustl(s))
    !
    return
  end subroutine writeascheader

  subroutine writeasc_i4_r4(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    integer(I4B), dimension(ncol,nrow), intent(in) :: x
    real(R4B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    integer(I4B) :: lun, icol, irow
! ------------------------------------------------------------------------------

    write(*,'(1x,a,1x,2a)') 'Writing',trim(f),'...'

    lun = getlun(); open(unit=lun,file=f,status='replace')
    call writeascheader(lun, ncol, nrow, real(xll,R8B), real(yll,R8B), real(cs,R8B),        &
      dble(nodata))
    write(lun,*)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(lun)
    !
    return
  end subroutine writeasc_i4_r4

  subroutine writeasc_r4_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    real(R4B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    integer(I4B) :: lun, icol, irow
! ------------------------------------------------------------------------------

    write(*,'(1x,a,1x,2a)') 'Writing',trim(f),'...'

    lun = getlun(); open(unit=lun,file=f,status='replace')
    call writeascheader(lun, ncol, nrow, xll, yll, cs, nodata)
    write(lun,*)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(lun)
    !
    return
  end subroutine writeasc_r4_r8

  subroutine readflt_header(iu, ncol, nrow, xll, yll, cs, nodata, &
    nbits, pixeltype)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(out) :: ncol, nrow
    character(len=*), intent(inout) :: nodata
    real(R8B), intent(out) :: xll, yll, cs
    integer(I4B), intent(out) :: nbits
    character(len=*), intent(inout) :: pixeltype
    ! -- locals
    !
    integer(I4B), parameter :: i_nrows         =  1 
    integer(I4B), parameter :: i_ncols         =  2
    integer(I4B), parameter :: i_nbits         =  3
    integer(I4B), parameter :: i_pixeltype     =  4
    integer(I4B), parameter :: i_ulxmap        =  5
    integer(I4B), parameter :: i_ulymap        =  6
    integer(I4B), parameter :: i_xllcorner     =  7
    integer(I4B), parameter :: i_yllcorner     =  8
    integer(I4B), parameter :: i_xdim          =  9
    integer(I4B), parameter :: i_ydim          = 10
    integer(I4B), parameter :: i_cellsize      = 11
    integer(I4B), parameter :: i_nodata        = 12
    integer(I4B), parameter :: i_nodata_value  = 13
    integer(I4B), parameter :: nkey = i_nodata_value
    
    logical :: lx, ly
    character(len=1) :: cdum
    character(len=MXSLEN) :: s, k
    integer(I4B) :: ios, nrewind, nfound
    integer(I4B), dimension(nkey) :: flag
    real(R8B) :: xcul, ycul
! ------------------------------------------------------------------------------
    flag = 0
    nrewind = 0
    lx = .false.; ly = .false.
    do while(.true.)
      read(unit=iu,iostat=ios,fmt='(a)') s
      if (ios == 0) then
        read(s,*) k
        k = change_case(k, 'l')
        select case(k)
        case('ncols')
          read(s,*) cdum, ncol; flag(i_ncols) = 1
        case('nrows')
          read(s,*) cdum, nrow; flag(i_nrows) = 1
        case('xllcorner')
          read(s,*) cdum, xll; flag(i_xllcorner) = 1; flag(i_ulxmap) = 1
        case('ulxmap') ! center
          read(s,*) cdum, xcul; flag(i_xllcorner) = 1; flag(i_ulxmap) = 1
          lx = .true.
        case('yllcorner')
          read(s,*) cdum, yll; flag(i_yllcorner) = 1; flag(i_ulymap) = 1
        case('ulymap') ! center
          read(s,*) cdum, ycul; flag(i_yllcorner) = 1; flag(i_ulymap) = 1
          ly = .true.
        case('nodata','nodata_value')
          read(s,*) cdum, nodata; flag(i_nodata) = 1; flag(i_nodata_value) = 1
        case('cellsize','xdim', 'ydim')
          read(s,*) cdum, cs; flag(i_cellsize) = 1; flag(i_xdim) = 1; flag(i_ydim) = 1
        case('nbits')
          read(s,*) cdum, nbits; flag(i_nbits) = 1
        case('pixeltype')
          read(s,*) cdum, pixeltype; flag(i_pixeltype) = 1
          pixeltype = change_case(pixeltype, 'l')
        end select
      end if
      !
      nfound = sum(flag)
      if (nfound == nkey) exit
      if (ios /= 0) then
        if (nrewind < nfound) then
          nrewind = nrewind + 1; rewind(iu)
        else
          if ((flag(i_nodata) == 0).or.(flag(i_nodata_value) == 0)) then
            flag(i_nodata) = 1; flag(i_nodata_value) = 1
            nodata = '0'
            nrewind = nrewind + 1; rewind(iu)
          else
            call errmsg('Could not parse Ehdr header file.')
          end if
        end if
      end if
    end do
    close(iu)
    !
    if (lx) then
      xll = xcul - cs/2.d0
    end if
    if (ly) then
      yll = ycul - cs*nrow + cs/2.d0
    end if
    !
    return
    end subroutine readflt_header
  
  subroutine readflt_i1(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(out) :: ncol, nrow
    integer(I1B), dimension(:,:), allocatable, intent(inout) :: x
    integer(I1B), intent(out) :: nodata
    real(R8B), intent(out) :: xll, yll, cs
    ! -- locals
    character(len=1) :: cdum
    character(len=MXSLEN) :: f, nodata_s, pixeltype
    integer(I4B) :: iu, nbits, icol, irow
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu)
    call readflt_header(iu, ncol, nrow, xll, yll, cs, nodata_s, nbits, pixeltype)
    if ((nbits /= 8).or.(pixeltype /= 'signedint')) then
      call errmsg('Could not read '//trim(f))
    end if
    read(nodata_s,*) nodata
    !
    if (allocated(x)) deallocate(x)
    allocate(x(ncol,nrow))
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'r', .true.)
    read(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine readflt_i1
  
  subroutine readflt_i4(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(out) :: ncol, nrow
    integer(I4B), dimension(:,:), allocatable, intent(inout) :: x
    integer(I4B), intent(out) :: nodata
    real(R8B), intent(out) :: xll, yll, cs
    ! -- locals
    character(len=1) :: cdum
    character(len=MXSLEN) :: f, nodata_s, pixeltype
    integer(I4B) :: iu, nbits, icol, irow
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu)
    call readflt_header(iu, ncol, nrow, xll, yll, cs, nodata_s, nbits, pixeltype)
    if ((nbits /= 32).or.(pixeltype /= 'signedint')) then
      call errmsg('Could not read '//trim(f))
    end if
    read(nodata_s,*) nodata
    !
    if (allocated(x)) deallocate(x)
    allocate(x(ncol,nrow))
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'r', .true.)
    read(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine readflt_i4
  
  subroutine readflt_r4(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(out) :: ncol, nrow
    real(R4B), dimension(:,:), allocatable, intent(inout) :: x
    real(R4B), intent(out) :: nodata
    real(R8B), intent(out) :: xll, yll, cs
    ! -- locals
    character(len=1) :: cdum
    character(len=MXSLEN) :: f, nodata_s, pixeltype
    integer(I4B) :: iu, nbits, icol, irow
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu)
    call readflt_header(iu, ncol, nrow, xll, yll, cs, nodata_s, nbits, pixeltype)
    if ((nbits /= 32).or.(pixeltype /= 'float')) then
      call errmsg('Could not read '//trim(f))
    end if
    read(nodata_s,*) nodata
    !
    if (allocated(x)) deallocate(x)
    allocate(x(ncol,nrow))
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'r', .true.)
    read(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine readflt_r4
  
  subroutine writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, nodata, &
  nbits, pixeltype)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: iu, ncol, nrow
    character(len=*), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    integer(I4B), intent(in) :: nbits
    character(len=*), intent(in) :: pixeltype
    ! -- locals
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
    close(iu)
    !
    return
  end subroutine writeflt_header_r4
  
  subroutine writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, nodata, &
  nbits, pixeltype)
! ******************************************************************************
    ! -- arguments
    integer(I4B), intent(in) :: iu, ncol, nrow
    character(len=*), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    integer(I4B), intent(in) :: nbits
    character(len=*), intent(in) :: pixeltype
    ! -- locals
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
    close(iu)
    !
    return
  end subroutine writeflt_header_r8
!
  subroutine writeflt_i1_r4(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/ehdr.html#raster-ehdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I1B), dimension(ncol,nrow), intent(in) :: x
    integer(I1B), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), 8, 'signedint')
    close(iu)
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
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i1_r4
  
 subroutine writeflt_i4_r4(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/ehdr.html#raster-ehdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I4B), dimension(ncol,nrow), intent(in) :: x
    integer(I4B), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), 32, 'signedint')
    close(iu)
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
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i4_r4
 
  subroutine writeflt_r4_r4(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/ehdr.html#raster-ehdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    real(R4B), dimension(ncol,nrow), intent(in) :: x
    real(R4B), intent(in) :: nodata
    real(R4B), intent(in) :: xll, yll, cs
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    call writeflt_header_r4(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), 32, 'float')
    close(iu)
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
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_r4_r4
  
  subroutine writeflt_i1_r8(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/ehdr.html#raster-ehdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I1B), dimension(ncol,nrow), intent(in) :: x
    integer(I1B), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), 8, 'signedint')
    close(iu)
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
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i1_r8
  
 subroutine writeflt_i4_r8(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/ehdr.html#raster-ehdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    integer(I4B), dimension(ncol,nrow), intent(in) :: x
    integer(I4B), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), 32, 'signedint')
    close(iu)
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
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_i4_r8
 
  subroutine writeflt_r4_r8(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/ehdr.html#raster-ehdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    real(R4B), dimension(ncol,nrow), intent(in) :: x
    real(R4B), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/nodata/)), 32, 'float')
    close(iu)
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
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    write(iu)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    !
    return
  end subroutine writeflt_r4_r8
!
  subroutine writeflt_r8_r8(fp, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
! see: https://gdal.org/drivers/raster/ehdr.html#raster-ehdr
    ! -- arguments
    character(len=*), intent(in) :: fp
    integer(I4B), intent(in) :: ncol, nrow
    real(R8B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: nodata
    real(R8B), intent(in) :: xll, yll, cs
    ! -- locals
    character(len=MXSLEN) :: f
    integer(I4B) :: iu, icol, irow, n
    real(R4B), dimension(:,:), allocatable :: r4x
    real(R4B) :: r4nodata
! ------------------------------------------------------------------------------
    f = trim(fp)//'.hdr'
    call open_file(f, iu, 'w')
    r4nodata = real(nodata,R4B)
    call writeflt_header_r8(iu, ncol, nrow, xll, yll, cs, ta((/r4nodata/)), 32, 'float')
    close(iu)
    !
    f = trim(fp)//'.flt'
    call open_file(f, iu, 'w', .true.)
    allocate(r4x(ncol,nrow))
    do irow = 1, nrow
      do icol = 1, ncol
        if (x(icol,irow) /= nodata) then
          r4x(icol,irow) = real(x(icol,irow),R4B)
        else
          r4x(icol,irow) = r4nodata
        end if
      end do
    end do
    !
    ! count
    n = 0
    do irow = 1,nrow
      do icol = 1, ncol
        if (r4x(icol,irow) /= r4nodata) n = n + 1
      end do
    end do
    call logmsg('# data cells: '//ta((/n/)))
    !
    write(iu)((r4x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(iu)
    deallocate(r4x)
    !
    return
  end subroutine writeflt_r8_r8
!
  subroutine writeidf_i1_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    integer(I1B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    real(R8B), dimension(:,:), allocatable :: wrk
    integer(I4B) :: ic, ir
! ------------------------------------------------------------------------------
    allocate(wrk(ncol,nrow))
    do ir = 1, nrow
      do ic = 1, ncol
        wrk(ic,ir) = real(x(ic,ir),R8B)
      end do
    end do
    !
    call imod_utl_printtext('Writing '//trim(f),0)
    if (.not.idfwrite_wrapper(ncol, nrow, wrk, (/cs/), (/cs/), xll, yll, nodata, '', f, 4)) then
      call imod_utl_printtext('Could not write '//trim(f),2)
    end if
    deallocate(wrk)
    !
    return
  end subroutine writeidf_i1_r8

  subroutine writeidf_i2_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    integer(I2B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    real(R8B), dimension(:,:), allocatable :: wrk
    integer(I4B) :: ic, ir
! ------------------------------------------------------------------------------
    allocate(wrk(ncol,nrow))
    do ir = 1, nrow
      do ic = 1, ncol
        wrk(ic,ir) = real(x(ic,ir),R8B)
      end do
    end do
    !
    call imod_utl_printtext('Writing '//trim(f),0)
    if (.not.idfwrite_wrapper(ncol, nrow, wrk, (/cs/), (/cs/), xll, yll, nodata, '', f, 4)) then
      call imod_utl_printtext('Could not write '//trim(f),2)
    end if
    deallocate(wrk)
    !
    return
  end subroutine writeidf_i2_r8

  subroutine writeidf_i4_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    integer(I4B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    real(R8B), dimension(:,:), allocatable :: wrk
    integer(I4B) :: ic, ir
! ------------------------------------------------------------------------------
    allocate(wrk(ncol,nrow))
    do ir = 1, nrow
      do ic = 1, ncol
        wrk(ic,ir) = real(x(ic,ir),R8B)
      end do
    end do
    !
    call imod_utl_printtext('Writing '//trim(f),0)
    if (.not.idfwrite_wrapper(ncol, nrow, wrk, (/cs/), (/cs/), xll, yll, nodata, '', f, 4)) then
      call imod_utl_printtext('Could not write '//trim(f),2)
    end if
    deallocate(wrk)
    !
    return
  end subroutine writeidf_i4_r8

  subroutine writeidf_i8_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    integer(I8B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    real(R8B), dimension(:,:), allocatable :: wrk
    integer(I4B) :: ic, ir
! ------------------------------------------------------------------------------
    allocate(wrk(ncol,nrow))
    do ir = 1, nrow
      do ic = 1, ncol
        wrk(ic,ir) = real(x(ic,ir),R8B)
      end do
    end do
    !
    call imod_utl_printtext('Writing '//trim(f),0)
    if (.not.idfwrite_wrapper(ncol, nrow, wrk, (/cs/), (/cs/), xll, yll, nodata, '', f, 8)) then
      call imod_utl_printtext('Could not write '//trim(f),2)
    end if
    deallocate(wrk)
    !
    return
  end subroutine writeidf_i8_r8

  subroutine writeidf_r4_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    real(R4B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    real(R8B), dimension(:,:), allocatable :: wrk
    integer(I4B) :: ic, ir
! ------------------------------------------------------------------------------
    allocate(wrk(ncol,nrow))
    do ir = 1, nrow
      do ic = 1, ncol
        wrk(ic,ir) = real(x(ic,ir),R8B)
      end do
    end do
    !
    call imod_utl_printtext('Writing '//trim(f),0)
    if (.not.idfwrite_wrapper(ncol, nrow, wrk, (/cs/), (/cs/), xll, yll, nodata, '', f, 4)) then
      call imod_utl_printtext('Could not write '//trim(f),2)
    end if
    deallocate(wrk)
    !
    return
  end subroutine writeidf_r4_r8

  subroutine writeidf_r8_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    real(R8B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
! ------------------------------------------------------------------------------
    call imod_utl_printtext('Writing '//trim(f),0)
    if (.not.idfwrite_wrapper(ncol, nrow, x, (/cs/), (/cs/), xll, yll, nodata, '', f, 8)) then
      call imod_utl_printtext('Could not write '//trim(f),2)
    end if
    !
    return
  end subroutine writeidf_r8_r8

  subroutine writeasc_i4_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    integer(I4B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    integer(I4B) :: lun, icol, irow
! ------------------------------------------------------------------------------

    write(*,'(1x,a,1x,2a)') 'Writing',trim(f),'...'

    lun = getlun(); open(unit=lun,file=f,status='replace')
    call writeascheader(lun, ncol, nrow, xll, yll, cs, nodata)
    write(lun,*)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(lun)
    !
    return
  end subroutine writeasc_i4_r8

  subroutine writeasc_r8_r8(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ncol, nrow
    real(R8B), dimension(ncol,nrow), intent(in) :: x
    real(R8B), intent(in) :: xll, yll, cs, nodata
    ! -- locals
    integer(I4B) :: lun, icol, irow
! ------------------------------------------------------------------------------

    write(*,'(1x,a,1x,2a)') 'Writing',trim(f),'...'

    lun = getlun(); open(unit=lun,file=f,status='replace')
    call writeascheader(lun, ncol, nrow, xll, yll, cs, nodata)
    write(lun,*)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(lun)
    !
    return
  end subroutine writeasc_r8_r8

  subroutine readasc_r_d(f, x, ncol, nrow, xll, yll, cs, nodata, idebug)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: f
    real, dimension(:,:), allocatable :: x
    integer, intent(out) :: ncol, nrow
    double precision, intent(out) :: xll, yll, cs, nodata
    integer, intent(in), optional :: idebug
    ! --- local
    character(len=1) :: cdum
    integer :: lun, icol, irow, ideb
! ------------------------------------------------------------------------------
    ideb = 0
    if (present(idebug)) then
      ideb = idebug
    end if

    if (allocated(x)) then
      deallocate(x)
    end if

    write(*,*) 'Reading ',trim(f), '...'
    open(newunit=lun,file=f,action='read',status='old')
    read(lun,*) cdum, ncol
    read(lun,*) cdum, nrow
    read(lun,*) cdum, xll
    read(lun,*) cdum, yll
    read(lun,*) cdum, cs
    read(lun,*) cdum, nodata
    allocate(x(ncol,nrow))
    if(ideb == 1) then
      !ncol = ncol/10; nrow = nrow/10; cs=cs*10
      !deallocate(x); allocate(x(ncol,nrow))
      do irow = 1, nrow
        do icol = 1, ncol
          x(icol,irow) = 1
        end do
      end do
      write(*,*) '@@@@ debug mode readasc, setting grid to 1!'
      return
    end if
    read(lun,*)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(lun)
    !
    return
  end subroutine readasc_r_d

  subroutine readasc_r_r(f, x, ncol, nrow, xll, yll, cs, nodata, idebug)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: f
    real, dimension(:,:), allocatable :: x
    integer, intent(out) :: ncol, nrow
    real, intent(out) :: xll, yll, cs, nodata
    integer, intent(in), optional :: idebug
    ! --- local
    character(len=1) :: cdum
    integer :: lun, icol, irow, ideb
! ------------------------------------------------------------------------------
    ideb = 0
    if (present(idebug)) then
      ideb = idebug
    end if

    if (allocated(x)) then
      deallocate(x)
    end if

    write(*,*) 'Reading ',trim(f), '...'
    open(newunit=lun,file=f,action='read',status='old')
    read(lun,*) cdum, ncol
    read(lun,*) cdum, nrow
    read(lun,*) cdum, xll
    read(lun,*) cdum, yll
    read(lun,*) cdum, cs
    read(lun,*) cdum, nodata
    allocate(x(ncol,nrow))
    read(lun,*)((x(icol,irow),icol=1,ncol),irow=1,nrow)
    close(lun)
    !
    return
  end subroutine readasc_r_r

  function readidf_r8val(idf, icol, irow) result (r8val)
! ******************************************************************************
    ! -- modules
    use imod_idf, only: idfobj
    ! -- arguments

    type(idfobj), intent(in) :: idf
    integer(I4B), intent(in) :: icol, irow
    real(R8B) :: r8val
    ! --- local
    integer(I4B), parameter :: icf = 1
    integer(I4B) :: irec
    real(R4B) :: r4val
! ------------------------------------------------------------------------------
    if ((icol < 1).or.(icol > idf%ncol)) then
      call errmsg('readidf_val: invalid column: '//ta((/icol/)))
    end if
    if ((irow < 1).or.(irow > idf%nrow)) then
      call errmsg('readidf_val: invalid row: '//ta((/irow/)))
    end if
    !
    irec=icf +10  +abs(idf%ieq-1) *2    +idf%ieq*(idf%nrow+idf%ncol) +idf%itb*2
    irec=irec+  ((irow-1)*idf%ncol)+icol
    !
    if (idf%itype == 4) then
      read(idf%iu,rec=irec) r4val
      r8val = real(r4val,R8B)
    else
      read(idf%iu,rec=irec) r8val
    end if
    !
    return
  end function readidf_r8val
  
  subroutine readidf_block_i4(idf, ir0, ir1, ic0, ic1, arr, nodata_in)
! ******************************************************************************
    ! -- modules
    use imod_idf, only: idfobj
    ! -- arguments
    type(idfobj), intent(in) :: idf
    integer(I4B), intent(in) :: ir0, ir1, ic0, ic1
    integer(I4B), dimension(:,:), pointer, intent(inout) :: arr
    integer(I4B), intent(in), optional :: nodata_in
    ! --- local
    integer(I4B) :: nc, nr, ic, ir, jc, jr
    integer(I4B) :: nodata
    real(R4B) :: r4val
    real(R8B) :: r8val
! ------------------------------------------------------------------------------
    if (present(nodata_in)) then
      nodata = nodata_in
    else
      nodata = int(idf%nodata)
    end if
    !
    nc = ic1-ic0+1; nr = ir1-ir0+1
    write(sa(1),*) nc; write(sa(2),*) nr;
    call logmsg('Reading i4 block for '//trim(idf%fname)//' ('//trim(adjustl(sa(1)))//&
      ','//trim(adjustl(sa(2)))//')...')
    if (associated(arr)) deallocate(arr)
    allocate(arr(nc,nr))
    !
    do ir = ir0, ir1
      do ic = ic0, ic1
        jr = ir-ir0+1; jc = ic-ic0+1
        r8val = readidf_r8val(idf, ic, ir)
        if (r8val /= idf%nodata) then
          arr(jc,jr) = int(r8val)
        else
          arr(jc,jr) = nodata
        end if
      end do
    end do
    call logmsg('Done reading i4 block')
    !
    return
  end subroutine readidf_block_i4

  subroutine readidf_block_r4(idf, ir0, ir1, ic0, ic1, arr, nodata_in)
! ******************************************************************************
    ! -- modules
    use imod_idf, only: idfobj
    ! -- arguments
    type(idfobj), intent(in) :: idf
    integer(I4B), intent(in) :: ir0, ir1, ic0, ic1
    real(R4B), dimension(:,:), pointer, intent(inout) :: arr
    real(R4B), intent(in), optional :: nodata_in
    ! --- local
    integer(I4B) :: nc, nr, ic, ir, jc, jr
    real(R4B) :: nodata
    real(R4B) :: r4val
    real(R8B) :: r8val
! ------------------------------------------------------------------------------
    if (present(nodata_in)) then
      nodata = nodata_in
    else
      nodata = real(idf%nodata,R4B)
    end if
    !
    nc = ic1-ic0+1; nr = ir1-ir0+1
    write(sa(1),*) nc; write(sa(2),*) nr;
    call logmsg('Reading r4 block for '//trim(idf%fname)//' ('//trim(adjustl(sa(1)))//&
      ','//trim(adjustl(sa(2)))//')...')
    if (associated(arr)) deallocate(arr)
    allocate(arr(nc,nr))
    !
    do ir = ir0, ir1
      do ic = ic0, ic1
        jr = ir-ir0+1; jc = ic-ic0+1
        r8val = readidf_r8val(idf, ic, ir)
        if (r8val /= idf%nodata) then
          arr(jc,jr) = real(r8val,R4B)
        else
          arr(jc,jr) = nodata
        end if
      end do
    end do
    call logmsg('Done reading r4 block')
    !
    return
  end subroutine readidf_block_r4

  subroutine readidf_block_r8(idf, ir0, ir1, ic0, ic1, arr, nodata_in)
! ******************************************************************************
    ! -- modules
    use imod_idf, only: idfobj
    ! -- arguments
    type(idfobj), intent(in) :: idf
    integer(I4B), intent(in) :: ir0, ir1, ic0, ic1
    real(R8B), dimension(:,:), pointer, intent(inout) :: arr
    real(R8B), intent(in), optional :: nodata_in
    ! --- local
    integer(I4B) :: nc, nr, ic, ir, jc, jr
    real(R8B) :: nodata
    real(R4B) :: r4val
    real(R8B) :: r8val
! ------------------------------------------------------------------------------
    if (present(nodata_in)) then
      nodata = nodata_in
    else
      nodata = int(idf%nodata)
    end if
    !
    nc = ic1-ic0+1; nr = ir1-ir0+1
    write(sa(1),*) nc; write(sa(2),*) nr;
    call logmsg('Reading r8 block for '//trim(idf%fname)//' ('//trim(adjustl(sa(1)))//&
      ','//trim(adjustl(sa(2)))//')...')
    if (associated(arr)) deallocate(arr)
    allocate(arr(nc,nr))
    !
    do ir = ir0, ir1
      do ic = ic0, ic1
        jr = ir-ir0+1; jc = ic-ic0+1
        r8val = readidf_r8val(idf, ic, ir)
        if (r8val /= idf%nodata) then
          arr(jc,jr) = r8val
        else
          arr(jc,jr) = nodata
        end if
      end do
    end do
    call logmsg('Done reading r8 block')
    !
    return
  end subroutine readidf_block_r8

  subroutine readidf_i_r(f, x, ncol, nrow, xll, yll, cs, nodata, &
     in_ir0, in_ir1, in_ic0, in_ic1, in_idebug)
! ******************************************************************************
    ! -- modules
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    integer, dimension(:,:), allocatable :: x
    integer, intent(out) :: ncol, nrow
    real, intent(out) :: xll, yll, cs, nodata
    integer, intent(in), optional :: in_ir0, in_ir1, in_ic0, in_ic1, in_idebug
    ! --- local
    integer :: idebug, ir0, ir1, ic0, ic1
    integer :: rddata, icol, irow, jcol, jrow
    type(idfobj) :: idf
! ------------------------------------------------------------------------------
    idebug = 0
    if (present(in_idebug)) then
      idebug = in_idebug
    end if
    !
    rddata = 1
    if (present(in_ir0).or.present(in_ir1).or. &
        present(in_ic0).or.present(in_ic1)) then
      rddata = 0
    end if
    call imod_utl_printtext('Reading '//trim(f),0)
    if (.not.idfread(idf,f,rddata)) then
      call imod_utl_printtext('Could not read '//trim(f),2)
    end if
    if (rddata == 1) close(idf%iu)
    !
    ir0 = 1; if (present(in_ir0)) ir0 = max(1, in_ir0)
    ir1 = 1; if (present(in_ir1)) ir1 = min(idf%nrow, in_ir1)
    ic0 = 1; if (present(in_ic0)) ic0 = max(1, in_ic0)
    ic1 = 1; if (present(in_ic1)) ic1 = min(idf%ncol, in_ic1)
    !
    if (allocated(x)) then
      deallocate(x)
    end if
    !
    ncol   = ic1 - ic0 + 1
    nrow   = ir1 - ir0 + 1
    cs     = idf%dx
    xll    = idf%xmin + (ic0-1)*cs
    yll    = idf%ymin + (idf%nrow - ir1)*cs
    nodata = idf%nodata
    allocate(x(ncol,nrow))
    if(idebug == 1) then
      do irow = 1, nrow
        do icol = 1, ncol
          x(icol,irow) = 1
        end do
      end do
      write(*,*) '@@@@ debug mode readasc, setting grid to 1!'
      call idfdeallocatex(idf)
      return
    end if
    !
    if (rddata == 0) then
      allocate(idf%x(ncol,nrow))
      do irow = ir0, ir1
        do icol = ic0, ic1
          jrow = irow - ir0 + 1; jcol = icol - ic0 + 1
          idf%x(jcol,jrow) = readidf_r8val(idf, icol, irow)
        end do
      end do
      close(idf%iu)
    end if
    !
    do irow = 1, nrow
      do icol = 1, ncol
        if (idf%x(icol,irow) /= idf%nodata) then
          x(icol,irow) = int(idf%x(icol,irow))
        else
           x(icol,irow) = 0
        end if
      end do
    end do
    nodata = 0.
    call idfdeallocatex(idf)
    !
    return
  end subroutine readidf_i_r

  subroutine readidf_r_r(f, x, ncol, nrow, xll, yll, cs, nodata, idebug)
! ******************************************************************************
    ! -- modules
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    real, dimension(:,:), allocatable :: x
    integer, intent(out) :: ncol, nrow
    real, intent(out) :: xll, yll, cs, nodata
    integer, intent(in), optional :: idebug
    ! --- local
    integer :: icol, irow, ideb
    type(idfobj) :: idf
! ------------------------------------------------------------------------------
    ideb = 0
    if (present(idebug)) then
      ideb = idebug
    end if

    if (allocated(x)) then
      deallocate(x)
    end if

    call imod_utl_printtext('Reading '//trim(f),0)
    if (.not.idfread(idf,f,1)) then
      call imod_utl_printtext('Could not read '//trim(f),2)
    end if
    close(idf%iu)

    ncol   = idf%ncol
    nrow   = idf%nrow
    xll    = idf%xmin
    yll    = idf%ymin
    cs     = idf%dx
    nodata = idf%nodata
    allocate(x(ncol,nrow))
    if(ideb == 1) then
      do irow = 1, nrow
        do icol = 1, ncol
          x(icol,irow) = 1
        end do
      end do
      write(*,*) '@@@@ debug mode readasc, setting grid to 1!'
    else
      do irow = 1, nrow
        do icol = 1, ncol
          if (idf%x(icol,irow) /= idf%nodata) then
            x(icol,irow) = idf%x(icol,irow)
          else
            x(icol,irow) = 0.
          end if
        end do
      end do
      nodata = 0.d0
    end if
    call idfdeallocatex(idf)
    !
    return
  end subroutine readidf_r_r

  subroutine readidf_r4_r4(f, x, ncol, nrow, xll, yll, cs, nodata)
! ******************************************************************************
    ! -- modules
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    real(R4B), dimension(:,:), allocatable :: x
    integer(I4B), intent(out) :: ncol, nrow
    real(R4B), intent(out) :: xll, yll, cs
    real(R4B), intent(in) :: nodata
    ! --- local
    integer :: icol, irow
    type(idfobj) :: idf
! ------------------------------------------------------------------------------

    if (allocated(x)) then
      deallocate(x)
    end if

    call imod_utl_printtext('Reading '//trim(f),0)
    if (.not.idfread(idf,f,1)) then
      call imod_utl_printtext('Could not read '//trim(f),2)
    end if
    close(idf%iu)

    ncol   = idf%ncol
    nrow   = idf%nrow
    xll    = idf%xmin
    yll    = idf%ymin
    cs     = idf%dx
    allocate(x(ncol,nrow))
    do irow = 1, nrow
      do icol = 1, ncol
        if (idf%x(icol,irow) /= idf%nodata) then
          x(icol,irow) = idf%x(icol,irow)
        else
          x(icol,irow) = nodata
        end if
      end do
    end do
    call idfdeallocatex(idf)
    !
    return
  end subroutine readidf_r4_r4
  
  subroutine readidf_r_d(f, x, ncol, nrow, xll, yll, cs, nodata, idebug)
! ******************************************************************************
    ! -- modules
    use imod_idf
    ! -- arguments
    character(len=*), intent(in) :: f
    real, dimension(:,:), allocatable :: x
    integer, intent(out) :: ncol, nrow
    double precision, intent(out) :: xll, yll, cs, nodata
    integer, intent(in), optional :: idebug
    ! --- local
    integer :: icol, irow, ideb
    type(idfobj) :: idf
! ------------------------------------------------------------------------------

    ideb = 0
    if (present(idebug)) then
      ideb = idebug
    end if

    if (allocated(x)) then
      deallocate(x)
    end if

    call imod_utl_printtext('Reading '//trim(f),0)
    if (.not.idfread(idf,f,1)) then
      call imod_utl_printtext('Could not read '//trim(f),2)
    end if
    close(idf%iu)

    ncol   = idf%ncol
    nrow   = idf%nrow
    xll    = idf%xmin
    yll    = idf%ymin
    cs     = idf%dx
    nodata = idf%nodata
    allocate(x(ncol,nrow))
    if(ideb == 1) then
      do irow = 1, nrow
        do icol = 1, ncol
          x(icol,irow) = 1
        end do
      end do
      write(*,*) '@@@@ debug mode readasc, setting grid to 1!'
    else
      do irow = 1, nrow
        do icol = 1, ncol
          if (idf%x(icol,irow) /= idf%nodata) then
            x(icol,irow) = idf%x(icol,irow)
          else
            x(icol,irow) = 0.
          end if
        end do
      end do
      nodata = 0.d0
    end if
    call idfdeallocatex(idf)
    !
    return
  end subroutine readidf_r_d

  subroutine errmsg(msg)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: msg
! ------------------------------------------------------------------------------
    write(*,'(a)') 'Error: '//trim(msg)
    stop 1
  end subroutine errmsg

  subroutine logmsg(msg)
! ******************************************************************************
    ! -- arguments
    character(len=*), intent(in) :: msg
! ------------------------------------------------------------------------------
    write(*,'(a)') trim(msg)
    !
    return
  end subroutine logmsg

! $Id: quicksort.f90 558 2015-03-25 13:44:47Z larsnerger $

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Collection of subroutines to sort and return a one-dimensional array
!!!    as well as corresponding sorted index of the array a. Original code
!!!    (distributed under GNU Free licence 1.2) was taken from
!!!    http://rosettacode.org/wiki/Quicksort#Fortran and modified to
!!!    also return sorted index of the original array a.
!!!    Copyright (C) 2015  Sanita Vetra-Carvalho
!!!
!!!    This program is distributed under the Lesser General Public License (LGPL) version 3,
!!!    for more details see <https://www.gnu.org/licenses/lgpl.html>.
!!!
!!!    Email: s.vetra-carvalho @ reading.ac.uk
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> subroutine to sort using the quicksort algorithm
!! @param[in,out] a, an array of doubles to be sorted
!! @param[out] idx_a, an array of sorted indecies of original array a
!! @param[in] na, dimension of the array a

recursive subroutine quicksort_d(a,idx_a,na)

! DUMMY ARGUMENTS
integer(I4B), intent(in) :: na ! nr or items to sort
real(R8B), dimension(nA), intent(inout) :: a ! vector to be sorted
integer(I4B), dimension(nA), intent(inout) :: idx_a ! sorted indecies of a

! LOCAL VARIABLES
integer(I4B) :: left, right, mid
real(R8B) :: pivot, temp
integer(I4B) :: marker, idx_temp

if (nA > 1) then
! insertion sort limit of 47 seems best for sorting 10 million
! integers on Intel i7-980X CPU.  Derived data types that use
! more memory are optimized with smaller values - around 20 for a 16
! -byte type.
  if (nA > 47) then
  ! Do quicksort for large groups
  ! Get median of 1st, mid, & last points for pivot (helps reduce
  ! long execution time on some data sets, such as already
  ! sorted data, over simple 1st point pivot)
    mid = (nA+1)/2
    if (a(mid) >= a(1)) then
      if (a(mid) <= a(nA)) then
        pivot = a(mid)
      else if (a(nA) > a(1)) then
        pivot = a(nA)
      else
        pivot = a(1)
      end if
    else if (a(1) <= a(nA)) then
      pivot = a(1)
    else if (a(nA) > a(mid)) then
      pivot = a(nA)
    else
      pivot = a(mid)
    end if

    left = 0
    right = nA + 1

    do while (left < right)
      right = right - 1
      do while (A(right) > pivot)
        right = right - 1
      end do
      left = left + 1
      do while (A(left) < pivot)
        left = left + 1
      end do
      if (left < right) then
        temp = A(left)
        idx_temp = idx_a(left)
        A(left) = A(right)
        idx_a(left) = idx_a(right)
        A(right) = temp
        idx_a(right) = idx_temp
      end if
    end do

    if (left == right) then
      marker = left + 1
    else
      marker = left
    end if

    call quicksort_d(A(:marker-1),idx_A(:marker-1),marker-1)
    call quicksort_d(A(marker:),idx_A(marker:),nA-marker+1)

  else
      call InsertionSort_d(A,idx_a,nA)    ! Insertion sort for small groups is
      !  faster than Quicksort
  end if
end if

return
end subroutine quicksort_d

!> subroutine to sort using the insertionsort algorithm and return indecies
!! @param[in,out] a, an array of doubles to be sorted
!! @param[in,out] idx_a, an array of integers of sorted indecies
!! @param[in] na, dimension of the array a
subroutine InsertionSort_d(a,idx_a,na)

  ! DUMMY ARGUMENTS
  integer(I4B),intent(in) :: na
  real(R8B), dimension(nA), intent(inout) :: a
  integer(I4B),dimension(nA), intent(inout) :: idx_a

! LOCAL VARIABLES
  real(R8B) :: temp
  integer(I4B):: i, j
  integer(I4B):: idx_tmp

  do i = 2, nA
     j = i - 1
     temp = A(i)
     idx_tmp = idx_a(i)
     do
        if (j == 0) exit
        if (a(j) <= temp) exit
        A(j+1) = A(j)
        idx_a(j+1) = idx_a(j)
        j = j - 1
     end do
     a(j+1) = temp
     idx_a(j+1) = idx_tmp
  end do
  return
end subroutine InsertionSort_d

recursive subroutine quicksort_r(a,idx_a,na)

! DUMMY ARGUMENTS
integer(I4B), intent(in) :: na ! nr or items to sort
real(R4B), dimension(nA), intent(inout) :: a ! vector to be sorted
integer(I4B), dimension(nA), intent(inout) :: idx_a ! sorted indecies of a

! LOCAL VARIABLES
integer(I4B) :: left, right, mid
real(R4B) :: pivot, temp
integer(I4B) :: marker, idx_temp

if (nA > 1) then
! insertion sort limit of 47 seems best for sorting 10 million
! integers on Intel i7-980X CPU.  Derived data types that use
! more memory are optimized with smaller values - around 20 for a 16
! -byte type.
  if (nA > 47) then
  ! Do quicksort for large groups
  ! Get median of 1st, mid, & last points for pivot (helps reduce
  ! long execution time on some data sets, such as already
  ! sorted data, over simple 1st point pivot)
    mid = (nA+1)/2
    if (a(mid) >= a(1)) then
      if (a(mid) <= a(nA)) then
        pivot = a(mid)
      else if (a(nA) > a(1)) then
        pivot = a(nA)
      else
        pivot = a(1)
      end if
    else if (a(1) <= a(nA)) then
      pivot = a(1)
    else if (a(nA) > a(mid)) then
      pivot = a(nA)
    else
      pivot = a(mid)
    end if

    left = 0
    right = nA + 1

    do while (left < right)
      right = right - 1
      do while (A(right) > pivot)
        right = right - 1
      end do
      left = left + 1
      do while (A(left) < pivot)
        left = left + 1
      end do
      if (left < right) then
        temp = A(left)
        idx_temp = idx_a(left)
        A(left) = A(right)
        idx_a(left) = idx_a(right)
        A(right) = temp
        idx_a(right) = idx_temp
      end if
    end do

    if (left == right) then
      marker = left + 1
    else
      marker = left
    end if

    call quicksort_r(A(:marker-1),idx_A(:marker-1),marker-1)
    call quicksort_r(A(marker:),idx_A(marker:),nA-marker+1)

  else
      call InsertionSort_r(A,idx_a,nA)    ! Insertion sort for small groups is
      !  faster than Quicksort
  end if
end if

return
end subroutine quicksort_r

subroutine InsertionSort_r(a,idx_a,na)

  ! DUMMY ARGUMENTS
  integer(I4B),intent(in) :: na
  real(R4B), dimension(nA), intent(inout) :: a
  integer(I4B),dimension(nA), intent(inout) :: idx_a

! LOCAL VARIABLES
  real(R4B) :: temp
  integer(I4B):: i, j
  integer(I4B):: idx_tmp

  do i = 2, nA
     j = i - 1
     temp = A(i)
     idx_tmp = idx_a(i)
     do
        if (j == 0) exit
        if (a(j) <= temp) exit
        A(j+1) = A(j)
        idx_a(j+1) = idx_a(j)
        j = j - 1
     end do
     a(j+1) = temp
     idx_a(j+1) = idx_tmp
  end do
  return
end subroutine InsertionSort_r

subroutine get_unique_i4(x, xu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- dummy
  integer(I4B), dimension(:), intent(in) :: x
  integer(I4B), dimension(:), allocatable, intent(inout) :: xu
  ! -- local
  real(R4B), dimension(:), allocatable :: xr4
  integer(I4B), dimension(:), allocatable :: idx
  integer(I4B) :: n, i, nu, iact
  real(R4B) :: r4v
! ------------------------------------------------------------------------------
  !
  if (allocated(xu)) deallocate(xu)
  n = size(x)
  !
  if (n == 1) then
    allocate(xu(1))
    xu(1) = x(1)
    return
  end if
  !
  allocate(idx(n), xr4(n))
  do i = 1, n
    idx(i) = i; xr4(i) = real(x(i),R4B)
  end do
  !
  call quicksort_r(xr4, idx, n)
  !
  ! count and fill
  do iact = 1, 2
    r4v = xr4(1); nu = 1
    if (iact == 2) then
      xu(nu) = int(r4v,I4B)
    end if
    do i = 1, n
      if (r4v /= xr4(i)) then
        nu = nu + 1
        r4v = xr4(i)
        if (iact == 2) then
          xu(nu) = int(r4v,I4B)
        end if
      end if
    end do
    if (iact == 1) then
      if (nu > 0) then
        allocate(xu(nu))
      end if
    end if
  end do
  !
  deallocate(xr4, idx)
  !
  return
end subroutine get_unique_i4

subroutine addboundary_i(wrk, ncol, nrow)
! ******************************************************************************
  integer, intent(in) :: ncol, nrow
  integer, dimension(ncol,nrow), intent(inout) :: wrk

  integer :: icol, irow, jp

  !write(*,*) 'ADDING CONSTANT HEAD BOUNDARY!'
  do irow = 1, nrow
    do icol = 1, ncol
      jp = wrk(icol,irow)
      if (jp > 0) then
        ! N
        if (irow > 1) then
          if (wrk(icol,irow-1) == 0) then
            wrk(icol,irow-1) = -jp
          end if
        end if
        ! S
        if (irow < nrow) then
          if (wrk(icol,irow+1) == 0) then
            wrk(icol,irow+1) = -jp
          end if
        end if
        ! W
        if (icol > 1) then
          if (wrk(icol-1,irow) == 0) then
            wrk(icol-1,irow) = -jp
          end if
        end if
        ! E
        if (icol < ncol) then
          if (wrk(icol+1,irow) == 0) then
            wrk(icol+1,irow) = -jp
          end if
        end if
      end if
    end do
  end do

  !do irow = 1, nrow
  !  do icol = 1, ncol
  !    jp = wrk(icol,irow)
  !    wrk(icol,irow) = abs(jp)
  !  end do
  !end do
  !
  return
end subroutine addboundary_i

subroutine addboundary_i_list(wrk, ncol, nrow, icir, nicir)
! ******************************************************************************
  integer, intent(in) :: ncol, nrow
  integer, dimension(ncol,nrow), intent(inout) :: wrk
  integer, intent(in) :: nicir
  integer, dimension(:,:), intent(in) :: icir
  
  integer :: i, icol, irow, jp

  !write(*,*) 'ADDING CONSTANT HEAD BOUNDARY!'
  do i = 1, nicir
    icol = icir(1,i); irow = icir(2,i)
    jp = wrk(icol,irow)
    if (jp > 0) then
      ! N
      if (irow > 1) then
        if (wrk(icol,irow-1) == 0) then
          wrk(icol,irow-1) = -jp
        end if
      end if
      ! S
      if (irow < nrow) then
        if (wrk(icol,irow+1) == 0) then
          wrk(icol,irow+1) = -jp
        end if
      end if
      ! W
      if (icol > 1) then
        if (wrk(icol-1,irow) == 0) then
          wrk(icol-1,irow) = -jp
        end if
      end if
      ! E
      if (icol < ncol) then
        if (wrk(icol+1,irow) == 0) then
          wrk(icol+1,irow) = -jp
        end if
      end if
    end if
  end do

  !do irow = 1, nrow
  !  do icol = 1, ncol
  !    jp = wrk(icol,irow)
  !    wrk(icol,irow) = abs(jp)
  !  end do
  !end do
  !
  return
end subroutine addboundary_i_list

subroutine addboundary_d(wrk, ncol, nrow, nodata)
! ******************************************************************************
  integer, intent(in) :: ncol, nrow
  double precision, dimension(ncol,nrow), intent(inout) :: wrk
  double precision, intent(in) :: nodata

  integer :: icol, irow, jp

  !write(*,*) 'ADDING CONSTANT HEAD BOUNDARY!'
  do irow = 1, nrow
    do icol = 1, ncol
      jp = int(wrk(icol,irow))
      if (jp > 0) then
        ! N
        if (irow > 1) then
          if (wrk(icol,irow-1) == nodata) then
            wrk(icol,irow-1) = -dble(jp)
          end if
        end if
        ! S
        if (irow < nrow) then
          if (wrk(icol,irow+1) == nodata) then
            wrk(icol,irow+1) = -dble(jp)
          end if
        end if
        ! W
        if (icol > 1) then
          if (wrk(icol-1,irow) == nodata) then
            wrk(icol-1,irow) = -dble(jp)
          end if
        end if
        ! E
        if (icol < ncol) then
          if (wrk(icol+1,irow) == nodata) then
            wrk(icol+1,irow) = -dble(jp)
          end if
        end if
      end if
    end do
  end do

  !do irow = 1, nrow
  !  do icol = 1, ncol
  !    jp = wrk(icol,irow)
  !    wrk(icol,irow) = abs(jp)
  !  end do
  !end do
  !
  return
end subroutine addboundary_d

subroutine addboundary_r(wrk, nodata)
! ******************************************************************************
  real, dimension(:,:), intent(inout) :: wrk
  real, intent(in) :: nodata

  integer :: ncol, nrow, icol, irow, jp

  ncol = size(wrk,1); nrow = size(wrk,2)

  !write(*,*) 'ADDING CONSTANT HEAD BOUNDARY!'
  do irow = 1, nrow
    do icol = 1, ncol
      if (wrk(icol,irow) /= nodata) then
        jp = int(wrk(icol,irow))
        if (jp > 0) then
          ! N
          if (irow > 1) then
            if (wrk(icol,irow-1) == nodata) then
              wrk(icol,irow-1) = -real(jp)
            end if
          end if
          ! S
          if (irow < nrow) then
            if (wrk(icol,irow+1) == nodata) then
              wrk(icol,irow+1) = -real(jp)
            end if
          end if
          ! W
          if (icol > 1) then
            if (wrk(icol-1,irow) == nodata) then
              wrk(icol-1,irow) = -real(jp)
            end if
          end if
          ! E
          if (icol < ncol) then
            if (wrk(icol+1,irow) == nodata) then
              wrk(icol+1,irow) = -real(jp)
            end if
          end if
        end if
      end if
    end do
  end do

  !do irow = 1, nrow
  !  do icol = 1, ncol
  !    jp = wrk(icol,irow)
  !    wrk(icol,irow) = abs(jp)
  !  end do
  !end do
  !
  return
end subroutine addboundary_r

  subroutine calc_unique_i(p, ps, pu, unp, id, nbnd, xll, yll, cs, lbnd_in)
! ******************************************************************************
    ! -- arguments
    integer, dimension(:,:), intent(in) :: p
    integer, intent(in) :: ps
    integer, dimension(:,:), allocatable, intent(inout) :: pu
    type(tUnp), dimension(:), allocatable, intent(inout) :: unp
    integer, intent(out) :: id
    integer, intent(out) :: nbnd
    real, intent(in) :: xll, yll, cs
    logical, intent(in), optional :: lbnd_in
    ! --- local
    integer, parameter ::jp = 1, jn = 2, js = 3, jw = 4, je = 5
    integer, parameter ::jnw = 6, jne = 7, jsw = 8, jse = 9
    integer, parameter :: nsten = jse
    integer, dimension(2,nsten) :: s1

    integer, parameter :: nnmax = 100
    integer :: nst, ic, ir, jc, jr, ncol, nrow, n1, n2, i, j, k, iact, n, m
    integer :: ics, irs, ict, irt, is, it, jt, ic0, ic1, ir0, ir1, itmin, imin
    integer :: ictmin, irtmin, nc, nr, nlst
    integer, dimension(:,:), allocatable :: lst1, lst2, wrk
    double precision, dimension(:), allocatable :: ds
    integer, dimension(:), allocatable :: dsi
    double precision :: d, dmin
    logical :: ldone, lbnd
! ------------------------------------------------------------------------------
    if (ps == 5) then
      nst = je
    else
      nst = nsten
    end if

    lbnd = .false.
    if (present(lbnd_in)) then
      lbnd = lbnd_in
    end if

    ncol = size(p,1); nrow = size(p,2)
    if (allocated(pu)) then
      deallocate(pu)
    end if
    allocate(pu(ncol,nrow))
    
    nlst = max(nst,ncol*nrow)
    allocate(lst1(2,nlst), lst2(2,nlst), wrk(ncol,nrow))

    do ir = 1, nrow
      do ic = 1, ncol
        pu(ic,ir) = 0
        wrk(ic,ir) = 0
      end do
    end do

    !write(*,*) 'Computing unique parts...'
    !write(*,*) 'Min/max=',minval(p), maxval(p)

    id = 0
    do ir = 1, nrow
      do ic = 1, ncol
        if ((p(ic,ir) /= 0) .and. (pu(ic,ir) == 0)) then

          ! set stencil
          s1(1,jp) = ic;             s1(2,jp) = ir
          s1(1,jn) = ic;             s1(2,jn) = max(1,   ir-1)
          s1(1,js) = ic;             s1(2,js) = min(nrow,ir+1)
          s1(1,jw) = max(1,   ic-1); s1(2,jw) = ir
          s1(1,je) = min(ncol,ic+1); s1(2,je) = ir
          s1(1,jnw) = s1(1,jw); s1(2,jnw) = s1(2,jn)
          s1(1,jne) = s1(1,je); s1(2,jne) = s1(2,jn)
          s1(1,jsw) = s1(1,jw); s1(2,jsw) = s1(2,js)
          s1(1,jse) = s1(1,je); s1(2,jse) = s1(2,js)
          !
          id = id + 1

          n1 = 0
          jc = s1(1,1); jr = s1(2,1)
          pu(jc,jr) = id
          do i = 2, nst
            jc = s1(1,i); jr = s1(2,i)
            if ((abs(p(jc,jr)) > 0) .and. (pu(jc,jr) == 0)) then
              n1 = n1 + 1
              lst1(1,n1) = jc; lst1(2,n1) = jr
            end if
          end do

          ldone = .false.
          do while (.not.ldone)
            n2 = 0
            do i = 1, n1
              jc = lst1(1,i); jr = lst1(2,i)
              s1(1,jp) = jc;             s1(2,jp) = jr
              s1(1,jn) = jc;             s1(2,jn) = max(1,   jr-1)
              s1(1,js) = jc;             s1(2,js) = min(nrow,jr+1)
              s1(1,jw) = max(1,   jc-1); s1(2,jw) = jr
              s1(1,je) = min(ncol,jc+1); s1(2,je) = jr
              s1(1,jnw) = s1(1,jw); s1(2,jnw) = s1(2,jn)
              s1(1,jne) = s1(1,je); s1(2,jne) = s1(2,jn)
              s1(1,jsw) = s1(1,jw); s1(2,jsw) = s1(2,js)
              s1(1,jse) = s1(1,je); s1(2,jse) = s1(2,js)
              pu(jc,jr) = id
              do j = 2, nst
                jc = s1(1,j); jr = s1(2,j)
                if ((abs(p(jc,jr)) > 0) .and. (pu(jc,jr) == 0) .and. (wrk(jc,jr) == 0)) then
                  n2 = n2 + 1
                  lst2(1,n2) = jc; lst2(2,n2) = jr
                  wrk(jc,jr) = 1
                end if
              end do
            end do
            !
            if (n2 == 0) then
              !write(*,*) 'id =',id
              ldone = .true.
              exit
            end if
            !
            ! set list 1
            do i = 1, n1
              jc = lst1(1,i); jr = lst1(2,i)
              pu(jc,jr) = id
            end do
            !
            ! copy list, set work
            do i = 1, n2
              jc = lst2(1,i); jr = lst2(2,i)
              lst1(1,i) = jc; lst1(2,i) = jr
              wrk(jc,jr) = 0
            end do
            n1 = n2
          end do
        end if
      end do
    end do
    !
    !write(*,*) '# unique parts found:',id
    !
    ! cleanup
    deallocate(lst1, lst2)
    !
    ! determine bounding box
    if (allocated(unp)) then
      deallocate(unp)
    end if
    allocate(unp(id))
    !
    do i = 1, id
      unp(i)%ic0 = ncol
      unp(i)%ic1 = 0
      unp(i)%ir0 = nrow
      unp(i)%ir1 = 0
    end do

    do ir = 1, nrow
      do ic = 1, ncol
        i = pu(ic,ir)
        if (i /= 0) then
          unp(i)%n = unp(i)%n + 1
          if (p(ic,ir) < 0) then
            pu(ic,ir) = -pu(ic,ir)
          end if
          unp(i)%ic0 = min(unp(i)%ic0, ic)
          unp(i)%ic1 = max(unp(i)%ic1, ic)
          unp(i)%ir0 = min(unp(i)%ir0, ir)
          unp(i)%ir1 = max(unp(i)%ir1, ir)
        end if
      end do
    end do
    !
    !n = 0
    !do ir = 1, nrow
    !  do ic = 1, ncol
    !    i = p(ic,ir)
    !    if (i /= 0) then
    !      n = n + 1
    !    end if
    !  end do
    !end do

    do i = 1, id
      unp(i)%ncol = unp(i)%ic1 - unp(i)%ic0 + 1
      unp(i)%nrow = unp(i)%ir1 - unp(i)%ir0 + 1
      unp(i)%ncol = max(1, unp(i)%ncol)
      unp(i)%nrow = max(1, unp(i)%nrow)
    end do
    !
    ! count the bounds
    do is = 1, id
      ir0 = unp(is)%ir0; ir1 = unp(is)%ir1
      ic0 = unp(is)%ic0; ic1 = unp(is)%ic1
      do ir = max(ir0-1,1), min(ir1+1,nrow)
        do ic = max(ic0-1,1), min(ic1+1,ncol)
          wrk(ic,ir) = 0
        end do
      end do
      do ir = ir0, ir1
        do ic = ic0, ic1
          j = pu(ic,ir)
          if ((j > 0).and.(is == j)) then
            ! set stencil
            s1(1,jp) = ic;             s1(2,jp) = ir
            s1(1,jn) = ic;             s1(2,jn) = max(1,   ir-1)
            s1(1,js) = ic;             s1(2,js) = min(nrow,ir+1)
            s1(1,jw) = max(1,   ic-1); s1(2,jw) = ir
            s1(1,je) = min(ncol,ic+1); s1(2,je) = ir
            s1(1,jnw) = s1(1,jw); s1(2,jnw) = s1(2,jn)
            s1(1,jne) = s1(1,je); s1(2,jne) = s1(2,jn)
            s1(1,jsw) = s1(1,jw); s1(2,jsw) = s1(2,js)
            s1(1,jse) = s1(1,je); s1(2,jse) = s1(2,js)
            do i = 2, nst
              jc = s1(1,i); jr = s1(2,i)
              if ((pu(jc,jr) == 0).and.(wrk(jc,jr) == 0)) then
                unp(is)%nbnd = unp(is)%nbnd + 1
                wrk(jc,jr) = 1
              end if
            end do
          end if
        end do
      end do
    end do
    !
    deallocate(wrk)

    !call writeidf('unp.idf',pu,ncol,nrow,xll,yll,cs,0.)
    !stop
    !
    ! set boundary
    nbnd = 0

    if (.not.lbnd) return

    do i = 1, id
      do iact = 1, 2
        unp(i)%nbnd = 0
        do ir = unp(i)%ir0, unp(i)%ir1
          do ic = unp(i)%ic0, unp(i)%ic1
            if (pu(ic,ir) == -i) then
              unp(i)%nbnd = unp(i)%nbnd + 1
              if (iact == 2) then
                j = unp(i)%nbnd
                unp(i)%bnd(1,j) = ic
                unp(i)%bnd(2,j) = ir
              end if
            end if
          end do
        end do
        if (iact == 1) then
          n = max(unp(i)%nbnd, 1)
          allocate(unp(i)%bnd(2,n))
          allocate(unp(i)%out_bndnn(2,n))
          allocate(unp(i)%out_itnn(n))
          allocate(unp(i)%out_d(n))
          do j = 1, n
            unp(i)%out_itnn(j) = 0
            unp(i)%out_d(j) = 0
          end do
          nbnd = nbnd + n
        end if
      end do !iact
    end do !i
    !
    ! all to all
    do is = 1, id
      allocate(unp(is)%all_bndnn(2,id))
      allocate(unp(is)%all_d(id))
      ics = unp(is)%bnd(1,1); irs = unp(is)%bnd(2,1)
      do it = 1, id
        ict = unp(it)%bnd(1,1); irt = unp(it)%bnd(2,1)
        unp(is)%all_bndnn(1,it) = ict
        unp(is)%all_bndnn(2,it) = irt
        unp(is)%all_d(it) = sqrt( real((ics-ict)**2 + (irs-irt)**2 )) !P
      end do
    end do

    allocate(ds(id), dsi(id))
    !
    do is = 1, id
      do j = 1, unp(is)%nbnd
        do i = 1, id
          ds(i) = huge(0)
          dsi(i) = i
        end do
        !
        ics = unp(is)%bnd(1,j); irs = unp(is)%bnd(2,j)
        !
        ! distances to bounding box of neighboring parts
        do it = 1, id
          if (it == is) then
            ds(is) = huge(0.d0)
            cycle
          end if
          ic0 = unp(it)%ic0; ic1 =  unp(it)%ic1
          ir0 = unp(it)%ir0; ir1 =  unp(it)%ir1
          nc = unp(it)%ncol; nr = unp(it)%nrow

          ds(it) = min(ds(it), sqrt( real((ics-(ic0+nc/2))**2 + (irs-(ir0+nr/2))**2 ))) !P
          ds(it) = min(ds(it), sqrt( real((ics-ic0)**2 + (irs-ir0)**2 ))) !NW
          ds(it) = min(ds(it), sqrt( real((ics-ic1)**2 + (irs-ir0)**2 ))) !NE
          ds(it) = min(ds(it), sqrt( real((ics-ic0)**2 + (irs-ir1)**2 ))) !SW
          ds(it) = min(ds(it), sqrt( real((ics-ic1)**2 + (irs-ir1)**2 ))) !SE
        end do
        !
        call quicksort_d(ds, dsi, id)
        !
        dmin = huge(0)
!        do k = 1, nnmax
        do k = 1, id
          it = dsi(k)
          do i = 1, unp(it)%nbnd
            ict = unp(it)%bnd(1,i); irt = unp(it)%bnd(2,i)
            d = sqrt( real((ics-ict)**2 + (irs-irt)**2 ))
            if (d <= dmin) then
              itmin = it; imin = i; dmin = d
              ictmin = ict; irtmin = irt
              unp(is)%out_bndnn(1,j) = ict
              unp(is)%out_bndnn(2,j) = irt
              unp(is)%out_d(j)       = dmin
              unp(is)%out_itnn(j)    = itmin
            end if
          end do
        end do
      end do
    end do

    ! allocate for boundary map and fill
    do i = 1, id
      !unp(i)%ncol = unp(i)%ic1 - unp(i)%ic0 + 1
      !unp(i)%nrow = unp(i)%ir1 - unp(i)%ir0 + 1
      !unp(i)%ncol = max(1, unp(i)%ncol)
      !unp(i)%nrow = max(1, unp(i)%nrow)
      !
      allocate(unp(i)%bndmap(unp(i)%ncol, unp(i)%nrow))
      do ir = 1, unp(i)%nrow
        do ic = 1, unp(i)%ncol
          unp(i)%bndmap(ic,ir) = 0
        end do
      end do
      do j = 1, unp(i)%nbnd
        ic = unp(i)%bnd(1,j)
        ir = unp(i)%bnd(2,j)
        ic = ic - unp(i)%ic0 + 1
        ir = ir - unp(i)%ir0 + 1
        if (unp(i)%bndmap(ic,ir) /= 0) then
          write(*,*) 'Program error'
          stop 1
        end if
        unp(i)%bndmap(ic,ir) = j
      end do
    end do

    do iact = 1, 2
      m = 0
      do is = 1, id
        do j = 1, unp(is)%nbnd
          ic = unp(is)%bnd(1,j);       ir = unp(is)%bnd(2,j)
          jc = unp(is)%out_bndnn(1,j); jr = unp(is)%out_bndnn(2,j)
          it = unp(is)%out_itnn(j)
          !
          jc = jc - unp(it)%ic0 + 1; jr = jr - unp(it)%ir0 + 1
          if ((jc < 0).or.(jc > unp(it)%ncol)) then
            write(*,*) 'Program error, column out of range.'
            stop 1
          end if
          if ((jr < 0).or.(jr > unp(it)%nrow)) then
            write(*,*) 'Program error, row out of range.'
            stop 1
          end if
          k = unp(it)%bndmap(jc,jr)
          if (k < 0) then
            write(*,*) 'Program error, mapping out of range.'
            stop 1
          end if
          jc = unp(it)%out_bndnn(1,k); jr = unp(it)%out_bndnn(2,k)
          jt = unp(it)%out_itnn(k)
          !
          if ((is == jt).and.(ic == jc).and.(ir == jr)) then
            ! do nothing
          else
            unp(jt)%nin = unp(jt)%nin + 1
            if (iact == 2) then
              m = m + 1
              n = unp(jt)%nin
              unp(jt)%in_map(n) = k
              unp(jt)%in_srcp(n) = it
              unp(jt)%in_srci(n) = j
              unp(jt)%in_bndnn(1,n) = ic
              unp(jt)%in_bndnn(2,n) = ir
              unp(jt)%in_d(n) = unp(is)%out_d(j)
            end if
          end if
        end do
      end do
      if (iact == 1) then
        do i = 1, id
          n = max(1, unp(i)%nin)
          allocate(unp(i)%in_map(n))
          allocate(unp(i)%in_flag(n))
          allocate(unp(i)%in_srcp(n))
          allocate(unp(i)%in_srci(n))
          allocate(unp(i)%in_bndnn(2,n))
          allocate(unp(i)%in_d(n))
          unp(i)%nin = 0
          do j = 1, n
            unp(i)%in_flag(j) = 0
          end do
        end do
      end if
    end do
    !
    return
  end subroutine calc_unique_i

  subroutine calc_unique_r(p, ps, pu, unp, id, nbnd, xll, yll, cs, lbnd_in)
! ******************************************************************************
    ! -- arguments
    real, dimension(:,:), intent(in) :: p
    integer, intent(in) :: ps
    integer, dimension(:,:), allocatable, intent(inout) :: pu
    type(tUnp), dimension(:), allocatable, intent(inout) :: unp
    integer, intent(out) :: id
    integer, intent(out) :: nbnd
    real, intent(in) :: xll, yll, cs
    logical, intent(in), optional :: lbnd_in
    ! --- local
    integer, parameter ::jp = 1, jn = 2, js = 3, jw = 4, je = 5
    integer, parameter ::jnw = 6, jne = 7, jsw = 8, jse = 9
    integer, parameter :: nsten = jse
    integer, dimension(2,nsten) :: s1

    integer, parameter :: nnmax = 100
    integer :: nst, ic, ir, jc, jr, ncol, nrow, n1, n2, i, j, k, iact, n, m
    integer :: ics, irs, ict, irt, is, it, jt, ic0, ic1, ir0, ir1, itmin, imin
    integer :: ictmin, irtmin, nc, nr
    integer, dimension(:,:), allocatable :: lst1, lst2, wrk
    double precision, dimension(:), allocatable :: ds
    integer, dimension(:), allocatable :: dsi
    double precision :: d, dmin
    logical :: ldone, lbnd
! ------------------------------------------------------------------------------
    if (ps == 5) then
      nst = je
    else
      nst = nsten
    end if

    lbnd = .false.
    if (present(lbnd_in)) then
      lbnd = lbnd_in
    end if

    ncol = size(p,1); nrow = size(p,2)
    if (allocated(pu)) then
      deallocate(pu)
    end if
    allocate(pu(ncol,nrow))
    allocate(lst1(2,ncol*nrow), lst2(2,ncol*nrow), wrk(ncol,nrow))

    do ir = 1, nrow
      do ic = 1, ncol
        pu(ic,ir) = 0
        wrk(ic,ir) = 0
      end do
    end do

    write(*,*) 'Computing unique parts...'
    !write(*,*) 'Min/max=',minval(p), maxval(p)

    id = 0
    do ir = 1, nrow
      do ic = 1, ncol
        if ((int(p(ic,ir)) /= 0) .and. (pu(ic,ir) == 0)) then

          ! set stencil
          s1(1,jp) = ic;             s1(2,jp) = ir
          s1(1,jn) = ic;             s1(2,jn) = max(1,   ir-1)
          s1(1,js) = ic;             s1(2,js) = min(nrow,ir+1)
          s1(1,jw) = max(1,   ic-1); s1(2,jw) = ir
          s1(1,je) = min(ncol,ic+1); s1(2,je) = ir
          s1(1,jnw) = s1(1,jw); s1(2,jnw) = s1(2,jn)
          s1(1,jne) = s1(1,je); s1(2,jne) = s1(2,jn)
          s1(1,jsw) = s1(1,jw); s1(2,jsw) = s1(2,js)
          s1(1,jse) = s1(1,je); s1(2,jse) = s1(2,js)
          !
          id = id + 1

          n1 = 0
          jc = s1(1,1); jr = s1(2,1)
          pu(jc,jr) = id
          do i = 2, nst
            jc = s1(1,i); jr = s1(2,i)
            if ((abs(int(p(jc,jr))) > 0) .and. (pu(jc,jr) == 0)) then
              n1 = n1 + 1
              lst1(1,n1) = jc; lst1(2,n1) = jr
            end if
          end do

          ldone = .false.
          do while (.not.ldone)
            n2 = 0
            do i = 1, n1
              jc = lst1(1,i); jr = lst1(2,i)
              s1(1,jp) = jc;             s1(2,jp) = jr
              s1(1,jn) = jc;             s1(2,jn) = max(1,   jr-1)
              s1(1,js) = jc;             s1(2,js) = min(nrow,jr+1)
              s1(1,jw) = max(1,   jc-1); s1(2,jw) = jr
              s1(1,je) = min(ncol,jc+1); s1(2,je) = jr
              s1(1,jnw) = s1(1,jw); s1(2,jnw) = s1(2,jn)
              s1(1,jne) = s1(1,je); s1(2,jne) = s1(2,jn)
              s1(1,jsw) = s1(1,jw); s1(2,jsw) = s1(2,js)
              s1(1,jse) = s1(1,je); s1(2,jse) = s1(2,js)
              pu(jc,jr) = id
              do j = 2, nst
                jc = s1(1,j); jr = s1(2,j)
                if ((abs(p(jc,jr)) > 0) .and. (pu(jc,jr) == 0) .and. (wrk(jc,jr) == 0)) then
                  n2 = n2 + 1
                  lst2(1,n2) = jc; lst2(2,n2) = jr
                  wrk(jc,jr) = 1
                end if
              end do
            end do
            !
            if (n2 == 0) then
              !write(*,*) 'id =',id
              ldone = .true.
              exit
            end if
            !
            ! set list 1
            do i = 1, n1
              jc = lst1(1,i); jr = lst1(2,i)
              pu(jc,jr) = id
            end do
            !
            ! copy list, set work
            do i = 1, n2
              jc = lst2(1,i); jr = lst2(2,i)
              lst1(1,i) = jc; lst1(2,i) = jr
              wrk(jc,jr) = 0
            end do
            n1 = n2
          end do
        end if
      end do
    end do
    !
    write(*,*) '# unique parts found:',id
    !
    ! cleanup
    deallocate(lst1, lst2)
    !
    ! determine bounding box
    if (allocated(unp)) then
      deallocate(unp)
    end if
    allocate(unp(id))
    !
    do i = 1, id
      unp(i)%ic0 = ncol
      unp(i)%ic1 = 0
      unp(i)%ir0 = nrow
      unp(i)%ir1 = 0
    end do

    do ir = 1, nrow
      do ic = 1, ncol
        i = pu(ic,ir)
        if (i /= 0) then
          unp(i)%n = unp(i)%n + 1
          if (p(ic,ir) < 0) then
            pu(ic,ir) = -pu(ic,ir)
          end if
          unp(i)%ic0 = min(unp(i)%ic0, ic)
          unp(i)%ic1 = max(unp(i)%ic1, ic)
          unp(i)%ir0 = min(unp(i)%ir0, ir)
          unp(i)%ir1 = max(unp(i)%ir1, ir)
        end if
      end do
    end do
    !
    !n = 0
    !do ir = 1, nrow
    !  do ic = 1, ncol
    !    i = p(ic,ir)
    !    if (i /= 0) then
    !      n = n + 1
    !    end if
    !  end do
    !end do

    do i = 1, id
      unp(i)%ncol = unp(i)%ic1 - unp(i)%ic0 + 1
      unp(i)%nrow = unp(i)%ir1 - unp(i)%ir0 + 1
      unp(i)%ncol = max(1, unp(i)%ncol)
      unp(i)%nrow = max(1, unp(i)%nrow)
    end do
    !
    ! count the bounds
    do is = 1, id
      ir0 = unp(is)%ir0; ir1 = unp(is)%ir1
      ic0 = unp(is)%ic0; ic1 = unp(is)%ic1
      do ir = max(ir0-1,1), min(ir1+1,nrow)
        do ic = max(ic0-1,1), min(ic1+1,ncol)
          wrk(ic,ir) = 0
        end do
      end do
      do ir = ir0, ir1
        do ic = ic0, ic1
          j = pu(ic,ir)
          if ((j > 0).and.(is == j)) then
            ! set stencil
            s1(1,jp) = ic;             s1(2,jp) = ir
            s1(1,jn) = ic;             s1(2,jn) = max(1,   ir-1)
            s1(1,js) = ic;             s1(2,js) = min(nrow,ir+1)
            s1(1,jw) = max(1,   ic-1); s1(2,jw) = ir
            s1(1,je) = min(ncol,ic+1); s1(2,je) = ir
            s1(1,jnw) = s1(1,jw); s1(2,jnw) = s1(2,jn)
            s1(1,jne) = s1(1,je); s1(2,jne) = s1(2,jn)
            s1(1,jsw) = s1(1,jw); s1(2,jsw) = s1(2,js)
            s1(1,jse) = s1(1,je); s1(2,jse) = s1(2,js)
            do i = 2, nst
              jc = s1(1,i); jr = s1(2,i)
              if ((pu(jc,jr) == 0).and.(wrk(jc,jr) == 0)) then
                unp(is)%nbnd = unp(is)%nbnd + 1
                wrk(jc,jr) = 1
              end if
            end do
          end if
        end do
      end do
    end do
    !
    deallocate(wrk)

    !call writeidf('unp.idf',pu,ncol,nrow,xll,yll,cs,0.)
    !stop
    !
    ! set boundary
    nbnd = 0

    if (.not.lbnd) return

    do i = 1, id
      do iact = 1, 2
        unp(i)%nbnd = 0
        do ir = unp(i)%ir0, unp(i)%ir1
          do ic = unp(i)%ic0, unp(i)%ic1
            if (pu(ic,ir) == -i) then
              unp(i)%nbnd = unp(i)%nbnd + 1
              if (iact == 2) then
                j = unp(i)%nbnd
                unp(i)%bnd(1,j) = ic
                unp(i)%bnd(2,j) = ir
              end if
            end if
          end do
        end do
        if (iact == 1) then
          n = max(unp(i)%nbnd, 1)
          allocate(unp(i)%bnd(2,n))
          allocate(unp(i)%out_bndnn(2,n))
          allocate(unp(i)%out_itnn(n))
          allocate(unp(i)%out_d(n))
          do j = 1, n
            unp(i)%out_itnn(j) = 0
            unp(i)%out_d(j) = 0
          end do
          nbnd = nbnd + n
        end if
      end do !iact
    end do !i
    !
    ! all to all
    do is = 1, id
      allocate(unp(is)%all_bndnn(2,id))
      allocate(unp(is)%all_d(id))
      ics = unp(is)%bnd(1,1); irs = unp(is)%bnd(2,1)
      do it = 1, id
        ict = unp(it)%bnd(1,1); irt = unp(it)%bnd(2,1)
        unp(is)%all_bndnn(1,it) = ict
        unp(is)%all_bndnn(2,it) = irt
        unp(is)%all_d(it) = sqrt( real((ics-ict)**2 + (irs-irt)**2 )) !P
      end do
    end do

    allocate(ds(id), dsi(id))
    !
    do is = 1, id
      do j = 1, unp(is)%nbnd
        do i = 1, id
          ds(i) = huge(0)
          dsi(i) = i
        end do
        !
        ics = unp(is)%bnd(1,j); irs = unp(is)%bnd(2,j)
        !
        ! distances to bounding box of neighboring parts
        do it = 1, id
          if (it == is) then
            ds(is) = huge(0.d0)
            cycle
          end if
          ic0 = unp(it)%ic0; ic1 =  unp(it)%ic1
          ir0 = unp(it)%ir0; ir1 =  unp(it)%ir1
          nc = unp(it)%ncol; nr = unp(it)%nrow

          ds(it) = min(ds(it), sqrt( real((ics-(ic0+nc/2))**2 + (irs-(ir0+nr/2))**2 ))) !P
          ds(it) = min(ds(it), sqrt( real((ics-ic0)**2 + (irs-ir0)**2 ))) !NW
          ds(it) = min(ds(it), sqrt( real((ics-ic1)**2 + (irs-ir0)**2 ))) !NE
          ds(it) = min(ds(it), sqrt( real((ics-ic0)**2 + (irs-ir1)**2 ))) !SW
          ds(it) = min(ds(it), sqrt( real((ics-ic1)**2 + (irs-ir1)**2 ))) !SE
        end do
        !
        call quicksort_d(ds, dsi, id)
        !
        dmin = huge(0)
!        do k = 1, nnmax
        do k = 1, id
          it = dsi(k)
          do i = 1, unp(it)%nbnd
            ict = unp(it)%bnd(1,i); irt = unp(it)%bnd(2,i)
            d = sqrt( real((ics-ict)**2 + (irs-irt)**2 ))
            if (d <= dmin) then
              itmin = it; imin = i; dmin = d
              ictmin = ict; irtmin = irt
              unp(is)%out_bndnn(1,j) = ict
              unp(is)%out_bndnn(2,j) = irt
              unp(is)%out_d(j)       = dmin
              unp(is)%out_itnn(j)    = itmin
            end if
          end do
        end do
      end do
    end do

    ! allocate for boundary map and fill
    do i = 1, id
      !unp(i)%ncol = unp(i)%ic1 - unp(i)%ic0 + 1
      !unp(i)%nrow = unp(i)%ir1 - unp(i)%ir0 + 1
      !unp(i)%ncol = max(1, unp(i)%ncol)
      !unp(i)%nrow = max(1, unp(i)%nrow)
      !
      allocate(unp(i)%bndmap(unp(i)%ncol, unp(i)%nrow))
      do ir = 1, unp(i)%nrow
        do ic = 1, unp(i)%ncol
          unp(i)%bndmap(ic,ir) = 0
        end do
      end do
      do j = 1, unp(i)%nbnd
        ic = unp(i)%bnd(1,j)
        ir = unp(i)%bnd(2,j)
        ic = ic - unp(i)%ic0 + 1
        ir = ir - unp(i)%ir0 + 1
        if (unp(i)%bndmap(ic,ir) /= 0) then
          write(*,*) 'Program error'
          stop 1
        end if
        unp(i)%bndmap(ic,ir) = j
      end do
    end do

    do iact = 1, 2
      m = 0
      do is = 1, id
        do j = 1, unp(is)%nbnd
          ic = unp(is)%bnd(1,j);       ir = unp(is)%bnd(2,j)
          jc = unp(is)%out_bndnn(1,j); jr = unp(is)%out_bndnn(2,j)
          it = unp(is)%out_itnn(j)
          !
          jc = jc - unp(it)%ic0 + 1; jr = jr - unp(it)%ir0 + 1
          if ((jc < 0).or.(jc > unp(it)%ncol)) then
            write(*,*) 'Program error, column out of range.'
            stop 1
          end if
          if ((jr < 0).or.(jr > unp(it)%nrow)) then
            write(*,*) 'Program error, row out of range.'
            stop 1
          end if
          k = unp(it)%bndmap(jc,jr)
          if (k < 0) then
            write(*,*) 'Program error, mapping out of range.'
            stop 1
          end if
          jc = unp(it)%out_bndnn(1,k); jr = unp(it)%out_bndnn(2,k)
          jt = unp(it)%out_itnn(k)
          !
          if ((is == jt).and.(ic == jc).and.(ir == jr)) then
            ! do nothing
          else
            unp(jt)%nin = unp(jt)%nin + 1
            if (iact == 2) then
              m = m + 1
              n = unp(jt)%nin
              unp(jt)%in_map(n) = k
              unp(jt)%in_srcp(n) = it
              unp(jt)%in_srci(n) = j
              unp(jt)%in_bndnn(1,n) = ic
              unp(jt)%in_bndnn(2,n) = ir
              unp(jt)%in_d(n) = unp(is)%out_d(j)
            end if
          end if
        end do
      end do
      if (iact == 1) then
        do i = 1, id
          n = max(1, unp(i)%nin)
          allocate(unp(i)%in_map(n))
          allocate(unp(i)%in_flag(n))
          allocate(unp(i)%in_srcp(n))
          allocate(unp(i)%in_srci(n))
          allocate(unp(i)%in_bndnn(2,n))
          allocate(unp(i)%in_d(n))
          unp(i)%nin = 0
          do j = 1, n
            unp(i)%in_flag(j) = 0
          end do
        end do
      end if
    end do
    !
    return
  end subroutine calc_unique_r

  subroutine calc_unique_grid_r4(a, mv, id, min_id, max_id, bb_a, gnc)
! ******************************************************************************
    ! -- arguments
    real(R4B), dimension(:,:), intent(in) :: a
    real(R4B), intent(in) :: mv
    integer(I8B), dimension(:,:), allocatable, intent(out) :: id
    integer(I8B), intent(out) :: min_id
    integer(I8B), intent(out) :: max_id
    type(tBb), intent(in) :: bb_a
    integer(I4B), intent(in) :: gnc
    ! --- local
    logical :: lfnd, lfull
    !
    type(tUnp), dimension(:), allocatable :: unp
    type(tBb), dimension(:), allocatable :: bb
    integer(I4B) :: nc, nr, ir, ic, mxid, nid, i, j, i4v, mc, mr, jr, jc, i4dum
    integer(I4B) :: newid
    integer(I4B) :: n
    integer(I8B) :: uid
    !
    integer(I4B), dimension(:), allocatable :: i4wk1d, i4wk1d2, i4wk1d3
    integer(I4B), dimension(:,:), allocatable :: i4wk2d, pu
    !
    real(R4B) :: r4dum
! ------------------------------------------------------------------------------
    !
    nc = size(a,1); nr = size(a,2)
    !
    if (allocated(id)) deallocate(id)
    allocate(id(nc,nr))
    do ir = 1, nr
      do ic = 1, nc
        id(ic,ir) = 0
      end do
    end do
    !
    mxid = 0
    do ir = 1, nr
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          mxid = max(mxid,int(a(ic,ir),I4B))
        end if
      end do
    end do
    newid = mxid
    !
    allocate(i4wk1d(mxid))
    do i = 1, mxid
      i4wk1d(i) = 0
    end do
    do ir = 1, nr
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          i4v = int(a(ic,ir),I4B)
          i4wk1d(i4v) = 1
        end if
      end do
    end do
    !
    nid = 0
    do i = 1, mxid
      if (i4wk1d(i) == 1) then
        nid = nid + 1; i4wk1d(i) = nid
      end if
    end do
    !
    allocate(bb(nid), i4wk1d3(nid))
    do ir = 1, nr
      do ic = 1, nc
        if (a(ic,ir) /= mv) then
          i4v = int(a(ic,ir),I4B)
          i = i4wk1d(i4v)
          i4wk1d3(i) = i4v ! local id --> global id
          bb(i)%ir0 = min(ir,bb(i)%ir0); bb(i)%ir1 = max(ir,bb(i)%ir1)
          bb(i)%ic0 = min(ic,bb(i)%ic0); bb(i)%ic1 = max(ic,bb(i)%ic1)
        end if
      end do
    end do
    do i = 1, nid
      call logmsg('Processing '//ta((/i/))//'/'//ta((/nid/))//'...')
      bb(i)%ncol = bb(i)%ic1-bb(i)%ic0+1; bb(i)%nrow = bb(i)%ir1-bb(i)%ir0+1
      mc = bb(i)%ncol; mr = bb(i)%nrow
      allocate(i4wk2d(mc,mr))
      lfull = .true.
      do ir = bb(i)%ir0, bb(i)%ir1
        do ic = bb(i)%ic0, bb(i)%ic1
          jr = ir-bb(i)%ir0+1; jc = ic-bb(i)%ic0+1
          if (a(ic,ir) == i4wk1d3(i)) then
            i4wk2d(jc,jr) = 1
          else
            i4wk2d(jc,jr) = 0
            lfull = .false.
          end if
        end do
      end do
      !
      if (lfull) then
        jr = 1; jc = 1
        ir = jr+bb(i)%ir0-1; ic = jc+bb(i)%ic0-1 ! index a
        jr = ir+bb_a%ir0-1; jc = ic+bb_a%ic0-1 ! global extent
        !uid = int(jc,I8B)+ int((jr-1),I8B)*int(gnc,I8B)
        uid = int(i4wk1d3(i),I8B)
        do ir = bb(i)%ir0, bb(i)%ir1
          do ic = bb(i)%ic0, bb(i)%ic1
            if (a(ic,ir) == i4wk1d3(i)) then
              id(ic,ir) = uid
            end if
          end do
        end do
      else
        n = 0
        call calc_unique_i(i4wk2d, 9, pu, unp, n, i4dum, r4dum, r4dum, r4dum)
        jc = -1; jr = -1
        do j = 1, n ! loop over separate parts
          lfnd = .false.
          do ir = unp(j)%ir0, unp(j)%ir1
            do ic = unp(j)%ic0, unp(j)%ic1
              if (pu(ic,ir) == j) then
                jc = ic; jr = ir
                lfnd = .true.
              end if
              if (lfnd) exit
            end do
            if (lfnd) exit
          end do
          if ((jc < 0).or.(jr < 0)) then
            call errmsg('Program error')
          end if
          ir = jr+bb(i)%ir0-1; ic = jc+bb(i)%ic0-1 ! index a
          jr = ir+bb_a%ir0-1; jc = ic+bb_a%ic0-1 ! global extent
          !uid = int(jc,I8B)+ int((jr-1),I8B)*int(gnc,I8B)
          if (j > 1) then
            newid = newid + 1
            uid = int(newid,I8B)
          else
            uid = int(i4wk1d3(i),I8B)
          end if
          !
          do ir = unp(j)%ir0, unp(j)%ir1
            do ic = unp(j)%ic0, unp(j)%ic1
              if (pu(ic,ir) == j) then
                jr = ir+bb(i)%ir0-1; jc = ic+bb(i)%ic0-1
                id(jc,jr) = uid
              end if
            end do
          end do
        end do
      end if
      deallocate(i4wk2d)
    end do
    !
    min_id = huge(min_id); max_id = -huge(max_id)
    do ir = 1, nr
      do ic = 1, nc
        if (id(ic,ir) /= 0) then
          min_id = min(min_id,id(ic,ir))
          max_id = max(max_id,id(ic,ir))
        end if
      end do
    end do
    !
    call logmsg('-->'//ta((/(max_id-min_id)/1000000/))//' M')
    ! 
    return
  end subroutine calc_unique_grid_r4
  
  subroutine getidmap(wrk, ir0, ir1, ic0, ic1, maxid, idmap, ncat, idmapinv, catarea, idbb)
! ******************************************************************************
    ! -- arguments
    integer, dimension(:,:), intent(in) :: wrk
    integer, intent(out) :: maxid
    integer, intent(in) :: ir0, ir1, ic0, ic1
    integer, dimension(:), allocatable, intent(inout) :: idmap
    integer, intent(out) :: ncat
    integer, dimension(:), allocatable, intent(inout) :: idmapinv
    double precision, dimension(:), allocatable, intent(inout) :: catarea
    integer, dimension(:,:), allocatable, intent(inout) :: idbb
    ! --- local
    integer :: nc, nr, ic, ir, j, k, jc, jr, id
! ------------------------------------------------------------------------------
    nc = size(wrk,1); nr = size(wrk,2)
    !
    ! determine the mappings from/to the catchment IDs
    maxid = maxval(wrk)
    if (allocated(idmap)) deallocate(idmap)
    allocate(idmap(maxid))
    do j = 1, maxid
      idmap(j) = 0
    end do
    do ir = 1, nr
      do ic = 1, nc
        id = abs(wrk(ic,ir))
        if (id /= 0) then
          idmap(id) = idmap(id) + 1
        end if
      end do
    end do
    ncat = 0
    do j = 1, maxid
      if (idmap(j) > 0) then
        ncat = ncat + 1
      end if
    end do
    if (allocated(catarea)) deallocate(catarea)
    allocate(catarea(ncat))
    if (allocated(idmapinv)) deallocate(idmapinv)
    allocate(idmapinv(ncat))
    ncat = 0
    do j = 1, maxid
      if (idmap(j) > 0) then
        ncat = ncat + 1
        idmapinv(ncat) = j
        catarea(ncat) = idmap(j)
        idmap(j) = ncat
      end if
    end do
    !
    ! determine bounding box
    if (allocated(idbb)) deallocate(idbb)
    allocate(idbb(4,ncat))
    do j = 1, ncat
      idbb(1,j) = nr !ir0
      idbb(2,j) = 0  !ir1
      idbb(3,j) = nc !ic0
      idbb(4,j) = 0  !ic1
    end do
    !
    do ir = ir0, ir1
      do ic = ic0, ic1
        jr = ir-ir0+1; jc = ic-ic0+1
        j = abs(wrk(jc,jr))
        if (j /= 0) then
          k = idmap(j)
          idbb(1,k) = min(idbb(1,k), ir)
          idbb(2,k) = max(idbb(2,k), ir)
          idbb(3,k) = min(idbb(3,k), ic)
          idbb(4,k) = max(idbb(4,k), ic)
        end if
      end do
    end do
    !
    return
  end subroutine getidmap
  
  subroutine split_str(s, token, sa)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    character(*), intent(in) :: s
    character(len=1), intent(in) :: token
    character(len=MXSLEN), dimension(:), allocatable, intent(inout) :: sa
    ! -- local
    logical :: set
    integer(I4B) :: n, i, i0, i1
! ------------------------------------------------------------------------------
    if (allocated(sa)) deallocate(sa)
    !
    n = 1
    do i = 1, len_trim(s)
      if (s(i:i) == token) then
        n = n + 1
      end if
    end do
    !
    allocate(sa(n))
    !
    set = .false.; i0 = 1; n = 0
    do i = 1, len_trim(s)
      if (s(i:i) == token) then
        set = .true.
        i1 = i - 1
      end if
      if (i == len_trim(s)) then
        set = .true.
        i1 = i
      end if
      if (set) then
        n = n + 1
        sa(n) = trim(adjustl(s(i0:i1)))
        i0 = i + 1
        set = .false.
      end if
    end do
    !
    return
  end subroutine split_str

  function strip_ext(s_in) result(s_out)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    character(*), intent(in) :: s_in
    character(len=1) :: opt
    character(len(s_in)) :: s_out
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    i = index(s_in,'.',back=.true.)
    if (i > 0) then
      s_out = s_in(1:i-1)
    else
      s_out = s_in
    end if
    !
    return
  end function strip_ext

  function get_ext(s_in) result(s_out)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    character(*), intent(in) :: s_in
    character(len=1) :: opt
    character(len(s_in)) :: s_out
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    i = index(s_in,'.',back=.true.)
    if (i > 0) then
      s_out = change_case(trim(s_in(i:)),'l')
    else
      s_out = ''
    end if
    !
    return
  end function get_ext
  
  function get_dir(f) result(d)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    character(*), intent(in) :: f
    character(len(f))        :: d
    ! -- local
    character(len=1) :: slash
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    slash = get_slash()
    i = index(f, slash, back=.true.)
    if (i > 0) then
      d = f(1:i)
    else
      d = f
    end if
    !
    return
  end function get_dir
  
  function base_name(f) result (bn)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    character(*), intent(in) :: f
    character(len(f))        :: bn
    ! -- local
    character(len=1) :: slash
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    slash = get_slash()
    i = index(f, slash, back=.true.)
    if (i > 0) then
      bn = f(i+1:)
    else
      bn = f
    end if
    !
    return
  end function base_name
  
  !###====================================================================
  function change_case(str, opt) result (string)
  !###====================================================================
    character(*), intent(in) :: str
    character(len=1) :: opt
    character(len(str))      :: string

    integer :: ic, i

    character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

  !   Capitalize each letter if it is lowecase
    string = str
    if ((opt == 'U') .or. (opt == 'u')) then
      do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
      end do
    end if
    if ((opt == 'L') .or. (opt == 'l')) then
      do i = 1, LEN_TRIM(str)
        ic = INDEX(cap, str(i:i))
        if (ic > 0) string(i:i) = low(ic:ic)
      end do
    end if
    !
    return
  end function change_case
  
  function count_i1a(a, v) result(n)
! ******************************************************************************
    ! -- arguments
    integer(I1B), dimension(:), intent(in) :: a
    integer(I1B), intent(in) :: v 
    integer(I4B) :: n
    ! -- locals
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    n = 0
    do i = 1, size(a)
      if (a(i) == v) n = n + 1
    end do
    !
    return
  end function count_i1a
  
  subroutine fillgap_r4(x, nodata, xtgt, nine_point)
! ******************************************************************************
    ! -- arguments
    real(R4B), dimension(:,:), intent(inout) :: x
    real(R4B), intent(in) :: nodata
    real(R4B), intent(in) :: xtgt
    logical, intent(in), optional :: nine_point
    ! -- locals
    integer(I4B), parameter :: maxiter = 1000
    !
    integer(I4B), parameter :: i_e = 1
    integer(I4B), parameter :: i_w = 2
    integer(I4B), parameter :: i_s = 3
    integer(I4B), parameter :: i_n = 4
    integer(I4B), parameter :: i_sw = 5
    integer(I4B), parameter :: i_se = 6
    integer(I4B), parameter :: i_nw = 7
    integer(I4B), parameter :: i_ne = 8
    integer(I4B), parameter :: nsten = i_ne
    integer(I4B), dimension(2,nsten) :: sicir
    logical, dimension(nsten) :: lsten
    !
    logical :: tgt_flag, np
    integer(I1B), dimension(:,:), allocatable :: i1wrk
    integer(I4B) :: nc, nr, ic, ir, n, m, ic0, ic1, ir0, ir1, nbr, i, j, maxcnt
    integer(I4B) :: ntgt, iter, nnodata, jc, jr
    integer(I4B) :: bbic0, bbic1, bbir0, bbir1, bbjc0, bbjc1, bbjr0, bbjr1
    integer(I4B), dimension(:,:), allocatable :: i4wrk
    integer(I4B), dimension(8) :: i4idx, ucnt
    real(R4B), dimension(:), allocatable :: r4wrk
    real(R4B), dimension(8) :: r4nbr, r4ucnt
    real(R4B) :: r4huge, rval, rvalp, rval_min, rval_max
    real(R4B), parameter :: my_nodata  = -12345.
! ------------------------------------------------------------------------------
    !
    if (present(nine_point)) then
      np = nine_point
    else
      np = .true.
    end if
    !
    r4huge = huge(r4huge)
    !
    nc = size(x,1); nr = size(x,2)
    allocate(i1wrk(nc,nr),i4wrk(2,nc*nr), r4wrk(nc*nr))
    bbir0 = 1; bbir1 = nr; bbic0 = 1; bbic1 = nc
    !
    iter = 0; n = 1
    do while(.true.)
      iter = iter + 1
      !
      do ir = 1, nr
        do ic = 1, nc
          i1wrk(ic,ir) = 0
        end do
      end do
      !
      if (iter == 1) then
        rvalp = xtgt
       else
        rvalp = my_nodata
      end if
      !
      n = 0
      do ir = bbir0, bbir1
        do ic = bbic0, bbic1
          if (x(ic,ir) == rvalp) then
            ic0 = ic - 1; ic1 = ic + 1
            ir0 = ir - 1; ir1 = ir + 1
            nbr = 0; lsten = .false.
            !
            if (ic1 <= nc) then ! E
              jc = ic1; jr = ir; rval = x(jc,jr); sicir(1,i_e) = jc; sicir(2,i_e) = jr; lsten(i_e) = .true.
              if ((rval /= xtgt).and.(rval /= nodata).and.(rval /= my_nodata)) then
                nbr = nbr + 1; r4nbr(nbr) = rval
              end if
              if (np.and.(ir1 <= nr)) then !S
                jc = ic1; jr = ir1; rval = x(jc,jr); sicir(1,i_se) = jc; sicir(2,i_se) = jr; lsten(i_se) = .true.
                if ((rval /= xtgt).and.(rval /= nodata).and.(rval /= my_nodata)) then
                  nbr = nbr + 1; r4nbr(nbr) = rval
                end if
              end if
              if (np.and.(ir0 >= 1)) then !N
                jc = ic1; jr = ir0; rval = x(jc,jr); sicir(1,i_ne) = jc; sicir(2,i_ne) = jr; lsten(i_ne) = .true.
                if ((rval /= xtgt).and.(rval /= nodata).and.(rval /= my_nodata)) then
                  nbr = nbr + 1; r4nbr(nbr) = rval
                end if
              end if
            end if
            if (ic0 >= 1) then ! W
              jc = ic0; jr = ir; rval = x(jc,jr); sicir(1,i_w) = jc; sicir(2,i_w) = jr; lsten(i_w) = .true.
              if ((rval /= xtgt).and.(rval /= nodata).and.(rval /= my_nodata)) then
                nbr = nbr + 1; r4nbr(nbr) = rval
              end if
              if (np.and.(ir1 <= nr)) then !S
                jc = ic0; jr = ir1; rval = x(jc,jr); sicir(1,i_sw) = jc; sicir(2,i_sw) = jr; lsten(i_sw) = .true.
                if ((rval /= xtgt).and.(rval /= nodata).and.(rval /= my_nodata)) then
                  nbr = nbr + 1; r4nbr(nbr) = rval
                end if
              end if
              if (np.and.(ir0 >= 1)) then !N
                jc = ic0; jr = ir0; rval = x(jc,jr); sicir(1,i_nw) = jc; sicir(2,i_nw) = jr; lsten(i_nw) = .true.
                if ((rval /= xtgt).and.(rval /= nodata).and.(rval /= my_nodata)) then
                  nbr = nbr + 1; r4nbr(nbr) = rval
                end if
              end if
            end if
            if (ir1 <= nr) then !S
              jc = ic; jr = ir1; rval = x(jc,jr); sicir(1,i_s) = jc; sicir(2,i_s) = jr; lsten(i_s) = .true.
              if ((rval /= xtgt).and.(rval /= nodata).and.(rval /= my_nodata)) then
                nbr = nbr + 1; r4nbr(nbr) = rval
              end if
            end if
            if (ir0 >= 1) then !N
              jc = ic; jr = ir0; rval = x(jc,jr); sicir(1,i_n) = jc; sicir(2,i_n) = jr; lsten(i_n) = .true.
              rval = x(ic,ir0)
              if ((rval /= xtgt).and.(rval /= nodata).and.(rval /= my_nodata)) then
                nbr = nbr + 1; r4nbr(nbr) = rval
              end if
            end if
            !
            ! neighbors found
            if (nbr > 0) then
              do i = 1, nsten
                if (lsten(i)) then
                  jc = sicir(1,i); jr = sicir(2,i)
                  if (x(jc,jr) == xtgt) then
                    i1wrk(jc,jr) = 1
                  end if
                end if
              end do
              !
              n = n + 1
              i4wrk(1,n) = ic; i4wrk(2,n) = ir
              if (nbr <= 2) then
                r4wrk(n) = r4nbr(1)
              else
                ! check if all are the same
                rval_min = r4huge; rval_max = -r4huge
                do i = 1, nbr
                  rval_min = min(rval_min,r4nbr(i))
                  rval_max = max(rval_max,r4nbr(i))
                end do
                if (rval_min == rval_max) then
                  r4wrk(n) = rval_min
                else
                  do i = 1, nbr
                    i4idx(i) = i
                  end do
                  call quicksort_r(r4nbr, i4idx, nbr)
                  !
                  ucnt = 0; r4ucnt(1) = r4nbr(1); m = 1
                  do i = 1, nbr
                    if (r4nbr(i) /= r4ucnt(m)) then
                      m = m + 1
                      r4ucnt(m) = r4nbr(i)
                    end if
                    ucnt(m) = ucnt(m) + 1
                  end do
                  if (m == 1) then
                    r4wrk(n) = r4ucnt(1)
                  else
                    maxcnt = 0
                    do i = 1, m
                      maxcnt = max(maxcnt,ucnt(i))
                    end do
                    do i = 1, m
                      if (ucnt(i) == maxcnt) then
                        r4wrk(n) = r4ucnt(i)
                        exit
                      end if
                    end do
                  end if
                end if
              end if
            end if
          end if
        end do
      end do
      !
      ! set the target value
      !call logmsg('# cells labeled: '//ta((/100.*n/(6000*6000)/)))
      do i = 1, n
        ic = i4wrk(1,i); ir = i4wrk(2,i)
        x(ic,ir) = r4wrk(i)
        i1wrk(ic,ir) = 0
      end do
      !
      bbjr0 = nr+1; bbjr1 = 0; bbjc0 = nc+1; bbjc1 = 0
      do ir = 1, nr
        do ic = 1, nc
          if (i1wrk(ic,ir) == 1) then
            x(ic,ir) = my_nodata
            bbjr0 = min(bbjr0,ir); bbjr1 = max(bbjr1,ir)
            bbjc0 = min(bbjc0,ic); bbjc1 = max(bbjc1,ic)
          end if
        end do
      end do
      !
      ! set loop bounding box
      bbir0 = max(bbjr0-1,1); bbir1 = min(bbjr1+1,nr)
      bbic0 = max(bbjc0-1,1); bbic1 = min(bbjc1+1,nc)
      !
      ! count the remaining target values
      ntgt = 0; nnodata = 0
      do ir = 1, nr
        do ic = 1, nc
          if (x(ic,ir) == xtgt) then
            ntgt = ntgt + 1
          end if
          if (x(ic,ir) == my_nodata) then
            nnodata = nnodata + 1
          end if
        end do
      end do
      call logmsg('Iteration '//ta((/iter/))//'; # added: '//ta((/n/))//'; # remaining: '//ta((/ntgt,nnodata/)))
      ntgt = ntgt + nnodata
      !
      !if (iter == 1) exit
      if (iter == maxiter) then
        call errmsg('fillgap_r4: maximum iterations of '//ta((/iter/))//' reached.')
        exit
      end if
      if (n == 0) exit
    end do
    call logmsg('Total iterations: '//ta((/iter/))//'; # not filled: '//ta((/ntgt/)))
    !
    deallocate(i1wrk,i4wrk,r4wrk)
    return
  end subroutine fillgap_r4
  
  subroutine fill_with_nearest_i4(x, nodata, xtgt)
! ******************************************************************************
    ! -- arguments
    integer(I4B), dimension(:,:), intent(inout) :: x
    integer(I4B), intent(in) :: nodata
    integer(I4B), intent(in) :: xtgt
    ! -- locals
    logical :: lfound
    integer(I4B), dimension(:), allocatable :: cnt
    integer(I4B), dimension(:,:), allocatable :: icir
    integer(I4B), dimension(1) :: mloc
    integer(I4B) :: n, m, nc, nr, mc, mr, ic, ir, jc, jr, kc, kr, ic0, ic1, ir0, ir1, &
      jc0, jc1, jr0, jr1,  ntgt, iact, i, j, k, nb, id0, id1, id
    integer(I4B) :: i4v, i4vmin, i4vmax
    integer(I4B), dimension(:), allocatable :: i4vi, i4vb
! ------------------------------------------------------------------------------
    !
    nc = size(x,1); nr = size(x,2)
    allocate(i4vb(2*nc + 2*nr))
    !
    ! store the location to interpolate
    do iact = 1, 2
      ntgt = 0
      do ir = 1, nr
        do ic = 1, nc
          i4v = x(ic,ir)
          if (i4v /= nodata) then
            if (i4v == xtgt) then
              ntgt = ntgt + 1
              if (iact == 2) then
                icir(1,ntgt) = ic
                icir(2,ntgt) = ir
              end if
            end if
          end if
        end do
      end do
      if (iact == 1) then
        if (ntgt > 0) then
          allocate(icir(2,ntgt), i4vi(ntgt))
          do i = 1, ntgt
            i4vi(i) = nodata
          end do
        end if
      end if
    end do
    !
    if (ntgt == 0) then
      return
    else
      call logmsg('# interpolation cells: '//ta((/ntgt/)))
    end if
    !
    do i = 1, ntgt
      jc = icir(1,i); jr = icir(2,i)
      !
      lfound = .false.; n = 0
      do while(.not.lfound)
        n = n + 1
        ir0 = jr - n; ir1 = jr + n; ic0 = jc - n; ic1 = jc + n
        ir0 = max(1,ir0); ir1 = min(nr,ir1); ic0 = max(1,ic0); ic1 = min(nc,ic1)
        nb = 0; i4vmin = huge(i4vmin); i4vmax = -huge(i4vmax)
        !
        do j = 1, 4
          select case(j)
          case(1) !N
            jr0 = ir0; jr1 = ir0; jc0 = ic0; jc1 = ic1
          case(2) !S
            jr0 = ir1; jr1 = ir1; jc0 = ic0; jc1 = ic1
          case(3) !W
            jr0 = ir0 + 1; jr1 = ir1 - 1; jc0 = ic0; jc1 = ic0
          case(4) !E
            jr0 = ir0 + 1; jr1 = ir1 - 1; jc0 = ic1; jc1 = ic1
           end select
          !
          mr = ir1 - ir0 + 1; mc = ic1 - ic0 + 1
          if (.not.lfound) then
            do ir = jr0, jr0
              do ic = jc0, jc1
                i4v = x(ic,ir)
                if ((i4v /= nodata).and.(i4v /= xtgt)) then
                  lfound = .true.
                  nb = nb + 1
                  i4vb(nb) = i4v
                  i4vmin = min(i4vmin, i4v); i4vmax = max(i4vmax, i4v)
                end if
              end do
            end do
          end if
        end do
        !
        if (lfound) then
          if (nb == 1) then
            i4vi(i) = i4vb(1)
          else
            if (i4vmin == i4vmax) then
              i4vi(i) = i4vb(1)
            else
              id0 = int(i4vmin,I4B); id1 = int(i4vmax,I4B)
              m = id1 - id0 + 1
              allocate(cnt(m))
              do k = 1, m
                cnt(k) = 0
              end do
              do k = 1, nb
                id = int(i4vb(k),I4B) - id0 + 1
                cnt(id) = cnt(id) + 1
              end do
              mloc = maxloc(cnt); id = mloc(1) + id0 - 1
              i4vi(i) = id
              deallocate(cnt)
            end if
          end if
        end if
      end do
    end do
    !
    do i = 1, ntgt
      jc = icir(1,i); jr = icir(2,i)
      i4v = i4vi(i)
      if (i4v == nodata) then
        call errmsg('Invalid interpolated value')
      end if
      x(jc,jr) = i4v
    end do
    !
    if (allocated(icir)) deallocate(icir)
    if (allocated(i4vi)) deallocate(i4vi)
    !
    return
  end subroutine fill_with_nearest_i4
  
  subroutine fill_with_nearest_r4(x, nodata, xtgt)
! ******************************************************************************
    ! -- arguments
    real(R4B), dimension(:,:), intent(inout) :: x
    real(R4B), intent(in) :: nodata
    real(R4B), intent(in) :: xtgt
    ! -- locals
    logical :: lfound
    integer(I4B), dimension(:), allocatable :: cnt
    integer(I4B), dimension(:,:), allocatable :: icir
    integer(I4B), dimension(1) :: mloc
    integer(I4B) :: n, m, nc, nr, mc, mr, ic, ir, jc, jr, kc, kr, ic0, ic1, ir0, ir1, &
      jc0, jc1, jr0, jr1,  ntgt, iact, i, j, k, nb, id0, id1, id
    real(R4B) :: r4v, r4vmin, r4vmax
    real(R4B), dimension(:), allocatable :: r4vi, r4vb
! ------------------------------------------------------------------------------
    !
    nc = size(x,1); nr = size(x,2)
    allocate(r4vb(2*nc + 2*nr))
    !
    ! store the location to interpolate
    do iact = 1, 2
      ntgt = 0
      do ir = 1, nr
        do ic = 1, nc
          r4v = x(ic,ir)
          if (r4v /= nodata) then
            if (r4v == xtgt) then
              ntgt = ntgt + 1
              if (iact == 2) then
                icir(1,ntgt) = ic
                icir(2,ntgt) = ir
              end if
            end if
          end if
        end do
      end do
      if (iact == 1) then
        if (ntgt > 0) then
          allocate(icir(2,ntgt), r4vi(ntgt))
          do i = 1, ntgt
            r4vi(i) = nodata
          end do
        end if
      end if
    end do
    !
    if (ntgt == 0) then
      return
    else
      call logmsg('# interpolation cells: '//ta((/ntgt/)))
    end if
    !
    do i = 1, ntgt
      jc = icir(1,i); jr = icir(2,i)
      !
      lfound = .false.; n = 0
      do while(.not.lfound)
        n = n + 1
        ir0 = jr - n; ir1 = jr + n; ic0 = jc - n; ic1 = jc + n
        ir0 = max(1,ir0); ir1 = min(nr,ir1); ic0 = max(1,ic0); ic1 = min(nc,ic1); 
        nb = 0; r4vmin = huge(r4vmin); r4vmax = -huge(r4vmax)
        !
        do j = 1, 4
          select case(j)
          case(1) !N
            jr0 = ir0; jr1 = ir0; jc0 = ic0; jc1 = ic1
          case(2) !S
            jr0 = ir1; jr1 = ir1; jc0 = ic0; jc1 = ic1
          case(3) !W
            jr0 = ir0 + 1; jr1 = ir1 - 1; jc0 = ic0; jc1 = ic0
          case(4) !E
            jr0 = ir0 + 1; jr1 = ir1 - 1; jc0 = ic1; jc1 = ic1
           end select
          !
          mr = ir1 - ir0 + 1; mc = ic1 - ic0 + 1
          if (.not.lfound) then
            do ir = jr0, jr0
              do ic = jc0, jc1
                r4v = x(ic,ir)
                if ((r4v /= nodata).and.(r4v /= xtgt)) then
                  lfound = .true.
                  nb = nb + 1
                  r4vb(nb) = r4v
                  r4vmin = min(r4vmin, r4v); r4vmax = max(r4vmax, r4v)
                end if
              end do
            end do
          end if
        end do
        !
        if (lfound) then
          if (nb == 1) then
            r4vi(i) = r4vb(1)
          else
            if (r4vmin == r4vmax) then
              r4vi(i) = r4vb(1)
            else
              id0 = int(r4vmin,I4B); id1 = int(r4vmax,I4B)
              m = id1 - id0 + 1
              allocate(cnt(m))
              do k = 1, m
                cnt(k) = 0
              end do
              do k = 1, nb
                id = int(r4vb(k),I4B) - id0 + 1
                cnt(id) = cnt(id) + 1
              end do
              mloc = maxloc(cnt); id = mloc(1) + id0 - 1
              r4vi(i) = real(id,R4B)
              deallocate(cnt)
            end if
          end if
        end if
      end do
    end do
    !
    do i = 1, ntgt
      jc = icir(1,i); jr = icir(2,i)
      r4v = r4vi(i)
      if (r4v == nodata) then
        call errmsg('Invalid interpolated value')
      end if
      x(jc,jr) = r4v
    end do
    !
    if (allocated(icir)) deallocate(icir)
    if (allocated(r4vi)) deallocate(r4vi)
    !
    return
  end subroutine fill_with_nearest_r4
    
end module utilsmod

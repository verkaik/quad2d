module multigrid_module
  
  use utilsmod, only: I1B, I2B, I4B, I8B, R4B, R8B, tGrid, tGridArr, tBbx, tBb, ta, &
    I_I1, I_I2, I_I4, I_I8, I_R4, I_R8, base_name, MXSLEN, errmsg, logmsg, swap_slash
  use vrt_module, only: tVrt
  use hdrModule, only: writeflt
  
  implicit none
  
  ! modules
  !
  type, public :: tMultiGrid
    integer(I4B), dimension(:), allocatable :: grid_count
    integer(I4B) :: ngrid
    type(tGrid), dimension(:), pointer :: grid => null()
  contains
    procedure :: init           => tMultiGrid_init
    procedure :: clean          => tMultiGrid_clean
    procedure :: write_vrt      => tMultiGrid_write_vrt
    procedure :: set_grid_count => tMultiGrid_set_grid_count
  end type tMultiGrid
  !
  type, public :: tMultiGridArray
    integer(I4B) :: nmgrid
    type(tMultiGrid), dimension(:), pointer :: mga => null()
  contains
    procedure :: init      => tMultiGridArray_init
    procedure :: clean     => tMultiGridArray_clean
    procedure :: get_mg    => tMultiGridArray_get_mg
    procedure :: write_vrt => tMultiGridArray_write_vrt
  end type tMultiGridArray
  !
  save
  
  contains
  
! ==============================================================================
! ==============================================================================
! tMultiGrid
! ==============================================================================
! ==============================================================================

  subroutine tMultiGrid_init(this, n, xll, xur, yll, yur, csa, &
     mvi1, mvi2, mvi4, mvi8, mvr4, mvr8)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMultiGrid) :: this
    !
    integer(I4B), intent(in) :: n
    real(R8B), intent(in), optional :: xll
    real(R8B), intent(in), optional :: xur
    real(R8B), intent(in), optional :: yll
    real(R8B), intent(in), optional :: yur
    real(R8B), dimension(:), intent(in), optional :: csa
    !
    integer(I1B), intent(in), optional :: mvi1
    integer(I2B), intent(in), optional :: mvi2
    integer(I4B), intent(in), optional :: mvi4
    integer(I8B), intent(in), optional :: mvi8
    real(R4B),    intent(in), optional :: mvr4
    real(R8B),    intent(in), optional :: mvr8
    
    ! -- local
    type(tGrid), pointer :: g => null()
    type(tBbx) :: bbx
    integer(I4B) :: i, m
    logical :: usebbx
! ------------------------------------------------------------------------------
    !
    m = 0
    if (present(xll)) m = m + 1
    if (present(xur)) m = m + 1
    if (present(yll)) m = m + 1
    if (present(xur)) m = m + 1
    if (present(csa)) m = m + 1
    if (m == 5) then
      usebbx = .true.
    else
      usebbx = .false.
    end if
    !
    call this%clean()
    !
    this%ngrid = n
    !
    allocate(this%grid_count(this%ngrid))
    this%grid_count = 0
    allocate(this%grid(this%ngrid))
    !
    if (usebbx) then
      call bbx%init()
      bbx%xll = xll; bbx%xur = xur; bbx%yll = yll; bbx%yur = yur
      do i = 1, this%ngrid
        g => this%grid(i)
        bbx%cs = csa(i)
        call g%init(bbx=bbx, mvi1=mvi1, mvi2=mvi2, mvi4=mvi4, mvi8=mvi8, &
          mvr4=mvr4, mvr8=mvr8)
      end do
    else
      do i = 1, this%ngrid
        g => this%grid(i)
        call g%init(mvi1=mvi1, mvi2=mvi2, mvi4=mvi4, mvi8=mvi8, &
          mvr4=mvr4, mvr8=mvr8)
      end do
    end if
    !
    return
   end subroutine tMultiGrid_init

   subroutine tMultiGrid_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMultiGrid) :: this
    !
    type(tGrid), pointer :: g => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    if (allocated(this%grid_count)) deallocate(this%grid_count)
    if (associated(this%grid)) then
      do i = 1, this%ngrid
        g => this%grid(i)
        call g%clean()
      end do
      deallocate(this%grid)
    end if
    !
    this%ngrid = 0
    this%grid => null()
    !
    return
  end subroutine tMultiGrid_clean
   
   subroutine tMultiGrid_write_vrt(this, fp)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMultiGrid) :: this
    character(len=*), intent(in) :: fp
    !
    ! -- local
    type(tVrt), pointer :: vrt => null()
    type(tGridArr), pointer :: ga => null()
    type(tGrid),     pointer :: g  => null()
    !
    type(tBb),  dimension(:), pointer :: dst_bbi => null()
    type(tBbx), dimension(:), pointer :: dst_bbx => null()
    type(tBb), pointer :: bbi => null()
    type(tBbX), pointer :: bbx => null()
    character(len=MXSLEN), dimension(:), allocatable :: f_tile, mv_tile
    character(len=MXSLEN) :: f, data_type
    integer(I4B), dimension(:), allocatable :: wrk
    integer(I4B) :: ig, nc, nr, ng
! ------------------------------------------------------------------------------
    !
    ! count the number of grids having data
    ng = 0
    do ig = 1, this%ngrid
      g => this%grid(ig)
      if (g%count_data() > 0) then
        ng = ng + 1
      end if
    end do
    if (ng == 0) then
      call logmsg('tMultiGrid_write_vrt: no data found for '//base_name(fp))
      return
    end if
    !
    ! first, check if all data_type are the same
    allocate(wrk(ng)); ng = 0
    do ig = 1, this%ngrid
      g => this%grid(ig)
      if (g%count_data() > 0) then
        ng = ng + 1
        call g%get(i_data_type=wrk(ng))
      end if
    end do
    if (minval(wrk) /= maxval(wrk)) then
      call errmsg('tMultiGrid_write_vrt: inconsistent data types.')
    end if
    deallocate(wrk)
    !
    allocate(dst_bbi(ng))
    allocate(dst_bbx(ng))
    allocate(f_tile(ng))
    allocate(mv_tile(ng))
    !
    ng = 0
    do ig = 1, this%ngrid
      g => this%grid(ig)
      if (g%count_data() > 0) then
        ng = ng + 1
        bbi => dst_bbi(ng); bbx => dst_bbx(ng)
        nc = g%nc; nr = g%nr
        bbi%ic0 = 1; bbi%ic1 = nc; bbi%ir0 = 1; bbi%ir1 = nr
        bbi%ncol = nc; bbi%nrow = nr
        bbx = g%bbx
        !
        f = trim(fp)//'_g'//ta([ig],'(i2.2)')
        call swap_slash(f)
        !
        select case(g%i_data_type)
        case(i_i1)
          call errmsg('tMultiGrid_write_vrt: i1 not yet supported.')
        case(i_i2)
          call errmsg('tMultiGrid_write_vrt: i2 not yet supported.')
        case(i_i4)
          call writeflt(f, g%xi4, bbi%ncol, bbi%nrow, &
            bbx%xll, bbx%yll, bbx%cs, g%mvi4)
          mv_tile(ng) = ta([g%mvi4]); data_type = 'Int32'
        case(i_i8)
          call errmsg('tMultiGrid_write_vrt: i8 not yet supported.')
        case(i_r4)
          call writeflt(f, g%xr4, bbi%ncol, bbi%nrow, &
            bbx%xll, bbx%yll, bbx%cs, g%mvr4)
          mv_tile(ng) = ta([g%mvr4]); data_type = 'Float32'
        case(i_r8)
          call errmsg('tMultiGrid_write_vrt: r8 not yet supported.')
        end select
        !
        f_tile(ng) = trim(f)//'.flt'
        !
      end if
    end do
    !
    ! write the vrt-file
    allocate(vrt)
    call vrt%init_write(fp, data_type, f_tile, mv_tile, dst_bbi, dst_bbx)
    call vrt%write()
    call vrt%clean()
    deallocate(vrt)
    deallocate(f_tile, mv_tile)
    !
    return
  end subroutine tMultiGrid_write_vrt
    
  subroutine tMultiGrid_set_grid_count(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMultiGrid) :: this
    !
    ! -- local
    type(tGrid), pointer :: g  => null()
    integer(I4B) :: ig
! ------------------------------------------------------------------------------
    !
    do ig = 1, this%ngrid
      g => this%grid(ig)
      this%grid_count(ig) = g%count_data()
    end do
    !
    return
  end subroutine tMultiGrid_set_grid_count
   
! ==============================================================================
! ==============================================================================
! tMultiGridArray
! ==============================================================================
! ==============================================================================
  
  subroutine tMultiGridArray_init(this, nmgrid)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMultiGridArray) :: this
    integer(I4B), intent(in) :: nmgrid
    ! -- local
! ------------------------------------------------------------------------------
    call this%clean()
    this%nmgrid = nmgrid
    allocate(this%mga(this%nmgrid))
    !
    return
  end subroutine tMultiGridArray_init
   
  subroutine tMultiGridArray_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMultiGridArray) :: this
    ! -- local
    integer(I4B) :: i 
! ------------------------------------------------------------------------------
    if (associated(this%mga)) then
      do i = 1, this%nmgrid
        call this%mga(i)%clean()
      end do
      deallocate(this%mga)
    end if
    this%nmgrid = 0
    !
    return
  end subroutine tMultiGridArray_clean
   
  subroutine tMultiGridArray_write_vrt(this, fp)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMultiGridArray) :: this
    character(len=MXSLEN), intent(in) :: fp
    ! -- local
    type(tMultiGrid), pointer :: mg => null()
    character(len=MXSLEN) :: f
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    do i = 1, this%nmgrid
      mg => this%get_mg(i)
      f = trim(fp)//ta([i],'(i2.2)')
      call mg%write_vrt(f)
    end do
    !
    return
  end subroutine tMultiGridArray_write_vrt
  
  function tMultiGridArray_get_mg(this, i) result(mg)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMultiGridArray) :: this
    integer(I4B), intent(in) :: i
    type(tMultiGrid), pointer :: mg
    ! -- local
! ------------------------------------------------------------------------------
    mg => this%mga(i)
    return
  end function tMultiGridArray_get_mg
  
end module multigrid_module

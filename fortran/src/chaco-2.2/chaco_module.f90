module chaco_module
  ! -- modules
  use utilsmod, only: I4B, I8B, MXSLEN, logmsg, errmsg, open_file, ta
  use metis_module, only: tMetis
  !
  implicit none
  !
  private
  !
  integer(I4B), parameter         :: i_nodata              = 0
  !
  integer(I4B), parameter, public :: i_glob_part_linear    = 4
  integer(I4B), parameter, public :: i_glob_part_random    = 5
  integer(I4B), parameter, public :: i_glob_part_scattered = 6
  !
  integer(I4B), parameter, public :: i_loc_ref_kl          = 1
  integer(I4B), parameter, public :: i_loc_ref_none        = 2
  !
  integer(I4B), parameter, public :: i_part_dim_bisection      = 1
  integer(I4B), parameter, public :: i_part_dim_quadrisection  = 2
  integer(I4B), parameter, public :: i_part_dim_octasection    = 3
  !
  type, public :: tChaco
    character(len=MXSLEN) :: graph_file
    character(len=MXSLEN) :: output_file
    integer(I4B)          :: i_glob_part
    integer(I4B)          :: i_loc_ref
    integer(I4B)          :: i_part_dim
    integer(I4B)          :: nparts ! size of 1-D mesh
    integer(I4B)          :: nvtxs
    integer(I8B), dimension(:), allocatable :: part
  contains
    procedure :: init            => tChaco_init
    procedure :: init_from_metis => tChaco_init_from_metis
    procedure :: write_graph     => tChaco_write_graph
    procedure :: set_opt         => tChaco_set_opt
    procedure :: run             => tChaco_run
    procedure :: read_output     => tChaco_read_output
    procedure :: set_metis       => tChaco_set_metis
    procedure :: clean           => tChaco_clean
  end type tChaco

  save

contains

  subroutine tChaco_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
!
    ! -- dummy
    class(tChaco) :: this
    ! -- local
! ------------------------------------------------------------------------------
    this%graph_file  = 'chaco_graph.txt'
    this%output_file = 'chaco_output.txt'
    !call this%set_opt(i_glob_part_linear, i_loc_ref_kl, i_part_dim_bisection)
    call this%set_opt(i_glob_part_random, i_loc_ref_kl, i_part_dim_bisection)
    !call this%set_opt(i_glob_part_scattered, i_loc_ref_kl, i_part_dim_bisection)
    !
    return
  end subroutine
    
  subroutine tChaco_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
!
    ! -- dummy
    class(tChaco) :: this
    ! -- local
! ------------------------------------------------------------------------------
    this%graph_file = ''
    this%output_file = ''
    this%i_glob_part = i_nodata
    this%i_loc_ref   = i_nodata
    this%i_part_dim  = i_nodata
    this%nparts = 0
    this%nvtxs  = 0
    if (allocated(this%part)) deallocate(this%part)
    !
    return
  end subroutine tChaco_clean
    
  subroutine tChaco_init_from_metis(this, met)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
!
    ! -- dummy
    class(tChaco) :: this
    type(tMetis), pointer, intent(in) :: met
    ! -- local
! ------------------------------------------------------------------------------
    call this%clean()
    call this%init()
    this%nparts = met%nparts
    call this%write_graph(met%xadj, met%adjncy, met%vwgt, met%adjwgt)
    return
  end subroutine tChaco_init_from_metis

  subroutine tChaco_write_graph(this, xadj, adjncy, vwgt, adjwgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
!
    ! -- dummy
    class(tChaco) :: this
    integer(I8B), dimension(:), intent(in) :: xadj
    integer(I8B), dimension(:), intent(in) :: adjncy
    integer(I8B), dimension(:), intent(in) :: vwgt
    integer(I8B), dimension(:), intent(in) :: adjwgt
    ! -- local
    integer(I4B) :: nja, iu, i, j, k, n, n_max
    integer(I4B), dimension(:), allocatable :: i4a
! ------------------------------------------------------------------------------
    this%nvtxs = size(vwgt); nja = size(adjwgt)
    !
    call open_file(this%graph_file, iu, 'w')
    write(iu,'(a)') ta([this%nvtxs, nja])//' 011'
    !
    ! determine the maximum number of neighbors
    n_max = 0
    do i = 1, this%nvtxs
      n = xadj(i+1) - xadj(i)+1 + 1
      n_max = max(n_max, n)
    end do
    !
    n = 1 + 2*n_max
    allocate(i4a(n))
    !
    do i = 1, this%nvtxs
      n = 1; i4a(n) = vwgt(i)
      do j = xadj(i)+1, xadj(i+1)
        k = adjncy(j) + 1
        n = n + 1; i4a(n) = k
        n = n + 1; i4a(n) = adjwgt(j)
      end do
      !
      write(iu,'(a)') ta(i4a(1:n))
    end do
    !
    close(iu)
    deallocate(i4a)
    !
    return
  end subroutine tChaco_write_graph
  
  subroutine tChaco_set_opt(this, i_glob_part, i_loc_ref, i_part_dim)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
!
    ! -- dummy
    class(tChaco) :: this
    integer(I4B), intent(in), optional :: i_glob_part
    integer(I4B), intent(in), optional :: i_loc_ref
    integer(I4B), intent(in), optional :: i_part_dim
    ! -- local
! ------------------------------------------------------------------------------
    if (present(i_glob_part)) then
      this%i_glob_part = i_glob_part
    end if
    if (present(i_loc_ref)) then
      this%i_loc_ref   = i_loc_ref
    end if
    if (present(i_part_dim)) then
      this%i_part_dim  = i_part_dim
    end if
    !
    return
  end subroutine tChaco_set_opt
  
  subroutine tChaco_run(this, f_exe)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
!
    ! -- dummy
    class(tChaco) :: this
    character(len=*) :: f_exe
    ! -- local
    character(len=MXSLEN) :: f, s
    integer(I4B) :: iu
! ------------------------------------------------------------------------------
    !
    f = 'User_Params'
    call open_file(f, iu, 'w')
    write(iu,'(a)') 'OUTPUT_ASSIGN = TRUE'
    write(iu,'(a)') 'OUT_ASSIGN_INV = TRUE'
    write(iu,'(a)') 'ARCHITECTURE = 1'
    close(iu)
    !
    f = 'chaco.inp'
    call open_file(f, iu, 'w')
    write(iu,'(a)') trim(this%graph_file)
    write(iu,'(a)') trim(this%output_file)
    write(iu,'(a)') ta([this%i_glob_part])
    write(iu,'(a)') ta([this%i_loc_ref])
    write(iu,'(a)') ta([this%nparts])
    write(iu,'(a)') ta([this%i_part_dim])
    write(iu,'(a)') 'n'
    close(iu)
    
    s = trim(f_exe)//' < '//trim(f)
    call system(s)
    
    return
  end subroutine tChaco_run
  
  subroutine tChaco_read_output(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
!
    ! -- dummy
    class(tChaco) :: this
    ! -- local
    character(len=MXSLEN) :: s
    integer(I4B) :: iu, ipart, i, j, n
! ------------------------------------------------------------------------------
    
    call open_file(this%output_file, iu, 'r')
    !
    allocate(this%part(this%nvtxs))
    this%part = 0
    
    do ipart = 1, this%nparts
      read(iu,*) s; read(s,*) n
      do i = 1, n
        read(iu,*) s; read(s,*) j
        this%part(j) = ipart
      end do
    end do
    !
    close(iu)
    !
    if (minval(this%part) == 0) then
      call errmsg('tChaco_read_output')
    end if
    !
    return
  end subroutine tChaco_read_output
  
  subroutine tChaco_set_metis(this, met)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
!
    ! -- dummy
    class(tChaco) :: this
    type(tMetis), pointer, intent(inout) :: met
    
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    call this%read_output()
    !
    do i = 1, met%nvtxs
      met%part(i) = this%part(i) - 1
    end do
    !
    call met%calc_imbal()
    !
    return
  end subroutine tChaco_set_metis
end module chaco_module

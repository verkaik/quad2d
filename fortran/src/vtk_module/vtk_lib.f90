module vtk_module
  use utilsmod, only: I1B, I4B, R4B, R8B, R4ZERO, mxslen, open_file, errmsg, ta, base_name
  use kdtree2_module
  
  implicit none

  private
  
  integer(I4B), parameter, public :: i_vtk_voxel = 11
  
  type :: tMesh
    integer(I4B) :: nv = 0 ! number of vertices
    real(R8B), dimension(:,:), allocatable :: x  ! vertex coordinates
    integer(I4B) :: nc = 0 ! number of cells
    integer(I4B), dimension(:), allocatable :: ia ! cell --> vertices ia (nc+1)
    integer(I4B), dimension(:), allocatable :: ja ! cell --> vertices ja (nja)
    integer(I4B) :: nja
    integer(I4B), dimension(:), allocatable :: celltype ! VTK cell type 
  contains
    procedure :: init         => tMesh_init
    procedure :: clean        => tMesh_clean
    procedure :: set_uniform  => tMesh_set_uniform
  end type tMesh
  
  type :: tVtuVar
    character(len=MXSLEN) :: name = ''
    real(R4B), dimension(:), allocatable :: v
  contains
    procedure :: init  => tVtuVar_init
    procedure :: clean => tVtuVar_clean
    procedure :: set   => tVtuVar_set
  end type tVtuVar
  !
  type, public :: tVtu
    character(len=MXSLEN) :: f
    integer(I4B)          :: nvar = 0
    type(tMesh), pointer  :: mesh => null()
    type(tVtuVar), dimension(:), pointer ::  var => null()
  contains
    procedure :: init          => tVtu_init
    procedure :: clean         => tVtu_clean
    procedure :: set_var       => tVtu_set_var
    procedure :: write_asc     => tVtu_write_asc
    !procedure :: write_bin_b64 => tVtu_write_bin_b64
    procedure :: write_bin_raw => tVtu_write_bin_raw
  end type tVtu
  
  type :: tPvdDataPart
    integer(I4B) :: part = 0
    character(len=MXSLEN) :: partname = ''
    character(len=MXSLEN) :: f_vtu = ''
  end type tPvdDataPart
  
  type :: tPvdDataTime
    real(R4B) :: timestep = R4ZERO
    integer(I4B) :: npart = 0
    type(tPvdDataPart), dimension(:), pointer :: dspart => null() ! length npart
  contains
    procedure :: init  => tPvdDataTime_init
    procedure :: clean => tPvdDataTime_clean
    procedure :: set   => tPvdDataTime_set
    procedure :: get   => tPvdDataTime_get
  end type tPvdDataTime
  !
  type, public :: tPvd
    character(len=MXSLEN) :: f = ''
    integer(I4B) :: ndstime = 0
    type(tPvdDataTime), dimension(:), pointer :: dstime => null() ! length dstime
  contains
    procedure :: init  => tPvd_init
    procedure :: set   => tPvd_set
    procedure :: write => tPvd_write
    procedure :: clean => tPvd_clean
  end type tPvd
  !
  public :: kdtree_remove_doubles
  !
  contains
  
! ==============================================================================
! ==============================================================================
! tMesh
! ==============================================================================
! ==============================================================================

  subroutine tMesh_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMesh) :: this
! ------------------------------------------------------------------------------
    !
    this%nv = 0
    this%nc = 0
    this%nja = 0
    !
    return
  end subroutine tMesh_init

  subroutine tMesh_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMesh) :: this
! ------------------------------------------------------------------------------
    this%nv = 0
    this%nc = 0
    this%nja = 0
    !
    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%ia)) deallocate(this%ia)
    if (allocated(this%ja)) deallocate(this%ja)
    if (allocated(this%celltype)) deallocate(this%celltype)
    !
    return
  end subroutine tMesh_clean
  
  subroutine tMesh_set_uniform(this, celltype, x, nc, ja)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMesh) :: this
    integer(I4B), intent(in) :: celltype
    real(R8B), dimension(:,:), intent(in) :: x
    integer(I4B), intent(in) :: nc
    integer(I4B), dimension(:), intent(in) :: ja
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    this%nv = size(x,2)
    this%nc = nc
    allocate(this%x, source=x)
    !
    ! checks
    if (celltype == i_vtk_voxel) then
      this%nja = size(ja)
      if (8*nc /= this%nja) then
        call errmsg('vtk_voxel: invalid size of ja.')
      end if
    else
      call errmsg('Only vtk_voxel is supported.')
    end if
    !
    allocate(this%ja, source=ja)
    allocate(this%celltype(this%nc)); this%celltype = celltype
    !
    allocate(this%ia(this%nc+1))
    this%ia(1) = 1
    do i = 2, this%nc + 1
      this%ia(i) = this%ia(i-1) + 8
    end do
    !
    return
  end subroutine tMesh_set_uniform
  
! ==============================================================================
! ==============================================================================
! tVtuVar
! ==============================================================================
! ==============================================================================
  
  subroutine tVtuVar_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtuVar) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tVtuVar_init

  subroutine tVtuVar_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtuVar) :: this
! ------------------------------------------------------------------------------
    if (allocated(this%v)) then
      deallocate(this%v)
    end if
    this%name = ''
    !
    return
  end subroutine tVtuVar_clean
  
  subroutine tVtuVar_set(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtuVar) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tVtuVar_set
  
! ==============================================================================
! ==============================================================================
! tVtu
! ==============================================================================
! ==============================================================================
  
  subroutine tVtu_init(this, f)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
    character(len=*), intent(in), optional :: f
! ------------------------------------------------------------------------------
    if (present(f)) then
      this%f = f
    end if
    !
    allocate(this%mesh)
    call this%mesh%init()
    !
    return
  end subroutine tVtu_init

  subroutine tVtu_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
    ! -- local
    integer(I4B) :: ivar
! ------------------------------------------------------------------------------
    this%f = ''
    this%nvar = 0
    call this%mesh%clean()
    deallocate(this%mesh); this%mesh => null()
    do ivar = 1, this%nvar
      call this%var(ivar)%clean()
    end do
    deallocate(this%var); this%var => null()
    !
    return
  end subroutine tVtu_clean
  
  subroutine tVtu_set_var(this, names, var)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
    character(len=*), dimension(:), intent(in) :: names
    real(R4B), dimension(:,:), intent(in) :: var
    ! -- local
    type(tVtuVar), pointer :: v => null()
    integer(I4B) :: ivar
! ------------------------------------------------------------------------------
    !
    this%nvar = size(var,1)
    allocate(this%var(this%nvar))
    !
    do ivar = 1, this%nvar
      v => this%var(ivar)
      v%name = names(ivar)
      allocate(v%v, source=var(ivar,:))
    end do
    !
    return
  end subroutine tVtu_set_var
  
  subroutine tVtu_write_asc(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
    !
    ! -- local
    type(tMesh), pointer :: mesh => null()
    type(tVtuVar), pointer :: var => null()
    integer(I4B) :: iu, i, j
    character(len=MXSLEN) :: s, s1, s2
! ------------------------------------------------------------------------------
    mesh => this%mesh
    !
    call open_file(this%f, iu, 'w')
    !
    write(iu,'(a)') '<?xml version="1.0"?>'
    write(iu,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(iu,'(a)') '<UnstructuredGrid>'
    write(s1,*) mesh%nv
    write(s2,*) mesh%nc
    write(iu,'(5a)') '<Piece NumberOfPoints="',trim(adjustl(s1)),'" NumberOfCells="',trim(adjustl(s2)),'">'
    write(iu,'(a)') '<Points>'
    write(iu,'(a)') '<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="ascii">'
    !-------------------
    ! points
    !-------------------
    do i = 1, mesh%nv
      write(iu,'(a,1x,a,1x,a)') ta(mesh%x(:,i))
    end do
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '</Points>'
    write(iu,'(a)') '<Cells>'
    write(iu,'(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
    !-------------------
    ! connectivity
    !-------------------
    write(iu,*) (mesh%ja(i)-1,i=1,mesh%nja)
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
    !-------------------
    ! offsets
    !-------------------
    write(iu,*) (mesh%ia(i)-1,i=2,mesh%nc+1)
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '<DataArray type="Int32" Name="types" format="ascii">'
    !-------------------
    ! types
    !-------------------
    write(iu,*) (mesh%celltype(i),i=1,mesh%nc)
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '</Cells>'
    write(iu,'(a)') '<CellData>'
    !-------------------
    ! variables
    !-------------------
    do i = 1, this%nvar
      var => this%var(i)
      write(s,'(3a)')  '<DataArray type="Float32" Name="',trim(var%name),'" NumberOfComponents="1" format="ascii">'
      write(iu,'(a)') trim(s)
      ! variable
      write(iu,*) var%v
      write(iu,'(a)') '</DataArray>'
    end do
    !
    write(iu,'(a)') '</CellData>'
    write(iu,'(a)') '</Piece>'
    write(iu,'(a)') '</UnstructuredGrid>'
    write(iu,'(a)') '</VTKFile>'
    close(iu)
    !
    return
  end subroutine tVtu_write_asc
  
!  subroutine tVtu_write_bin_b64(this)
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- modules
!    use penf
!    use befor64
!    ! -- dummy
!    class(tVtu) :: this
!    !
!    ! -- local
!    type(tMesh), pointer :: mesh => null()
!    type(tVtuVar), pointer :: var => null()
!    integer(I4B) :: iu, i, j
!    character(len=MXSLEN) :: s, s1, s2
!    !
!    integer(I4B) :: n
!    character(1), parameter:: endrec = char(10) !< End-character for binary-record finalize. 
!    character(len=:), allocatable:: code64 ! base64 encoded string
!    integer(I8P), dimension(:), allocatable :: i8wrk
!    real(R8P), dimension(:), allocatable :: r8wrk
!    integer(I1P), dimension(:), allocatable :: varp
!! ------------------------------------------------------------------------------
!    mesh => this%mesh
!    !
!    call b64_init()
!    call open_file(this%f, iu, 'w', .true.)
!    !
!    write(iu) '<?xml version="1.0"?>'//endrec
!    write(iu) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//endrec
!    write(iu) '<UnstructuredGrid>'//endrec
!    write(s1,*) mesh%nv
!    write(s2,*) mesh%nc
!    write(s,'(5a)') '<Piece NumberOfPoints="',trim(adjustl(s1)),'" NumberOfCells="',trim(adjustl(s2)),'">'
!    write(iu) trim(s)//endrec
!    write(iu) '<Points>'//endrec
!    write(iu) '<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="binary">'//endrec
!    !-------------------
!    ! points
!    !-------------------
!    allocate(r8wrk(mesh%nv*3))
!    n = 0
!    do i = 1, mesh%nv
!      n = n + 1; r8wrk(n) = mesh%x(1,i)
!      n = n + 1; r8wrk(n) = mesh%x(2,i)
!      n = n + 1; r8wrk(n) = mesh%x(3,i)
!    end do
!    call pack_data(a1=[int(3*mesh%nv*BYR8P,I4P)],a2=r8wrk,packed=varp)
!    call b64_encode(n=varp,code=code64)
!    write(iu) trim(code64)//endrec
!    deallocate(r8wrk,varp,code64)
!    !-------------------
!    write(iu) '</DataArray>'//endrec
!    write(iu) '</Points>'//endrec
!    write(iu) '<Cells>'//endrec
!    write(iu) '<DataArray type="Int64" Name="connectivity" format="binary">'//endrec
!    !-------------------
!    ! connectivity
!    !-------------------
!    allocate(i8wrk(mesh%nja))
!    do i = 1, mesh%nja
!      i8wrk(i) = mesh%ja(i) - 1 
!    end do
!    call pack_data(a1=[int(mesh%nja*BYI8P,I4P)],a2=i8wrk,packed=varp)
!    call b64_encode(n=varp,code=code64)
!    write(iu) code64//endrec
!    deallocate(i8wrk,varp,code64)
!    !-------------------
!    write(iu) '</DataArray>'//endrec
!    write(iu) '<DataArray type="Int64" Name="offsets" format="binary">'//endrec
!    !-------------------
!    ! offsets
!    !-------------------
!    allocate(i8wrk(mesh%nc))
!    do i = 2, mesh%nc+1
!      i8wrk(i-1) = mesh%ia(i) - 1 
!    end do
!    call pack_data(a1=[int(mesh%nc*BYI8P,I4P)],a2=i8wrk,packed=varp)
!    call b64_encode(n=varp,code=code64)
!    write(iu) code64//endrec
!    deallocate(i8wrk,varp,code64)
!    !-------------------
!    write(iu) '</DataArray>'//endrec
!    write(iu) '<DataArray type="Int64" Name="types" format="binary">'//endrec
!    !-------------------
!    ! types
!    !-------------------
!    allocate(i8wrk(mesh%nc))
!    do i = 1, mesh%nc
!      i8wrk(i) = mesh%celltype(i)
!    end do
!    call pack_data(a1=[int(mesh%nc*BYI8P,I4P)],a2=i8wrk,packed=varp)
!    call b64_encode(n=varp,code=code64)
!    write(iu) code64//endrec
!    deallocate(i8wrk,varp,code64)
!    !-------------------
!    write(iu) '</DataArray>'//endrec
!    write(iu) '</Cells>'//endrec
!    write(iu) '<CellData>'//endrec
!    !-------------------
!    ! variables
!    !-------------------
!    allocate(r8wrk(mesh%nc))
!    do i = 1, this%nvar
!      var => this%var(i)
!      write(s,'(3a)')  '<DataArray type="Float32" Name="',trim(var%name),'" NumberOfComponents="1" format="binary">'
!      write(iu) trim(s)//endrec
!      !-------------------
!      do j = 1, mesh%nc
!        r8wrk(j) = var%v(j)
!      end do
!      call pack_data(a1=[int(mesh%nc*BYR8P,I4P)],a2=r8wrk,packed=varp)
!      call b64_encode(n=varp,code=code64)
!      write(iu) code64//endrec
!      deallocate(varp, code64)
!      write(iu) '</DataArray>'//endrec
!    end do
!    deallocate(r8wrk)
!    !
!    write(iu) '</CellData>'//endrec
!    write(iu) '</Piece>'//endrec
!    write(iu) '</UnstructuredGrid>'//endrec
!    write(iu) '</VTKFile>'//endrec
!    close(iu)
!    !
!    return
!  end subroutine tVtu_write_bin_b64
  
  subroutine tVtu_write_bin_raw(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
    !
    ! -- local
    type(tMesh), pointer :: mesh => null()
    type(tVtuVar), pointer :: var => null()
    integer(I4B) :: iu, i, j, offset
    character(len=MXSLEN) :: s, s1, s2
    character(1), parameter:: endrec = char(10) !< End-character for binary-record finalize. 
! ------------------------------------------------------------------------------
    mesh => this%mesh
    !
    call open_file(this%f, iu, 'w', .true.)
    !
    write(iu) '<?xml version="1.0"?>'//endrec
    write(iu) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//endrec
    write(iu) '<UnstructuredGrid>'//endrec
    write(s1,*) mesh%nv; write(s2,*) mesh%nc
    write(iu) '<Piece NumberOfPoints="'//trim(adjustl(s1))//'" NumberOfCells="',trim(adjustl(s2))//'">'//endrec
    write(iu) '<Points>'//endrec
    offset = 0; write(s,*) offset
    write(iu) '<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="append" offset="'//trim(adjustl(s))//'"/>'//endrec
    offset = offset + 3*mesh%nv*R8B + I4B; write(s,*) offset
    write(iu) '</Points>'//endrec
    write(iu) '<Cells>'//endrec
    write(iu) '<DataArray type="Int32" Name="connectivity" format="append" offset="'//trim(adjustl(s))//'"/>'//endrec
    offset = offset + mesh%nja*I4B + I4B; write(s,*) offset
    write(iu) '<DataArray type="Int32" Name="offsets" format="append" offset="'//trim(adjustl(s))//'"/>'//endrec
    offset = offset + mesh%nc*I4B + I4B; write(s,*) offset
    write(iu) '<DataArray type="Int32" Name="types" format="append" offset="'//trim(adjustl(s))//'"/>'//endrec
    offset = offset + mesh%nc*I4B + I4B; write(s,*) offset
    write(iu) '</Cells>'//endrec
    write(iu) '<CellData>'//endrec
    !
    do i = 1, this%nvar
      var => this%var(i)
      write(iu) '<DataArray type="Float32" Name="'//trim(var%name)// &
        '" NumberOfComponents="1" format="append" offset="'//trim(adjustl(s))//'"/>'//endrec
      offset = offset + mesh%nc*R4B + I4B; write(s,*) offset
    end do
    !
    write(iu) '</CellData>'//endrec
    write(iu) '</Piece>'//endrec
    write(iu) '</UnstructuredGrid>'//endrec
    write(iu) '<AppendedData encoding="raw">'//endrec
    !-------------------
    write(iu) '_'
    write(iu) int(3*mesh%nv*R8B,I4B), ((mesh%x(i,j),i=1,3),j=1,mesh%nv)
    write(iu) int(mesh%nja*I4B,I4B), (mesh%ja(i)-1,i=1,mesh%nja)
    write(iu) int(mesh%nc*I4B,I4B), (mesh%ia(i)-1,i=2,mesh%nc+1)
    write(iu) int(mesh%nc*I4B,I4B),(mesh%celltype(i),i=1,mesh%nc)
    do i = 1, this%nvar
      var => this%var(i)
      write(iu) int(mesh%nc*R4B,I4B),(var%v(j),j=1,mesh%nc)
    end do
    !-------------------
    write(iu) endrec//'</AppendedData>'//endrec
    write(iu) '</VTKFile>'
    close(iu)
    !
    return
  end subroutine tVtu_write_bin_raw
  
! ==============================================================================
! ==============================================================================
! tPvdDataTime
! ==============================================================================
! ==============================================================================
  
  subroutine tPvdDataTime_init(this, npart)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvdDataTime) :: this
    integer(I4B), intent(in) :: npart
! ------------------------------------------------------------------------------
    !
    this%timestep = 0
    this%npart = npart
    allocate(this%dspart(this%npart))
    !
    return
  end subroutine tPvdDataTime_init

  subroutine tPvdDataTime_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvdDataTime) :: this
    ! 
    type(tPvdDataPart), pointer :: dsp => null()
    integer(I4B) :: ip
! ------------------------------------------------------------------------------
    !
    deallocate(this%dspart)
    this%npart = 0
    this%timestep = R4ZERO
    !
    return
  end subroutine tPvdDataTime_clean
  
  subroutine tPvdDataTime_set(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvdDataTime) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tPvdDataTime_set
  
  subroutine tPvdDataTime_get(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvdDataTime) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tPvdDataTime_get
  
! ==============================================================================
! ==============================================================================
! tPvd
! ==============================================================================
! ==============================================================================
  
  subroutine tPvd_init(this, f, ndstime, npart)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvd) :: this
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: ndstime
    integer(I4B), intent(in) :: npart
    ! -- local
    type(tPvdDataTime), pointer :: dst => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    this%f = f
    this%ndstime = ndstime
    allocate(this%dstime(ndstime))
    !
    do i = 1, ndstime
      dst => this%dstime(i)
      call dst%init(npart=npart)
    end do
    !
    return
  end subroutine tPvd_init

  subroutine tPvd_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvd) :: this
    ! -- local
    type(tPvdDataTime), pointer :: dst => null()
    integer(I4B) :: it
! ------------------------------------------------------------------------------
    do it = 1, this%ndstime
      dst => this%dstime(it)
      call dst%clean()
    end do
    deallocate(this%dstime); this%dstime => null()
    !
    this%f = ''
    this%ndstime = 0
    !
    return
  end subroutine tPvd_clean
  
  subroutine tPvd_set(this, itimestep, ipart, f_vtu, timestep, partname)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvd) :: this
    integer(I4B), intent(in) :: itimestep
    integer(I4B), intent(in) :: ipart
    character(len=*), intent(in) :: f_vtu
    real(R4B), intent(in), optional :: timestep
    character(len=*), intent(in), optional :: partname
    ! -- local
    type(tPvdDataPart), pointer :: dsp => null()
    type(tPvdDataTime), pointer :: dst => null()
    character(len=MXSLEN) :: name
    real(R4B) :: ts
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    if (present(timestep)) then
      ts = timestep
    else
      ts = real(itimestep,R4B)
    end if
    if (present(partname)) then
      name = partname
    else
      name = ta([ipart])
    end if
    !
    dst => this%dstime(itimestep)
    dst%timestep = ts
    dsp => dst%dspart(ipart)
    dsp%part = ipart
    dsp%partname = name
    dsp%f_vtu = base_name(f_vtu)
    !
    return
  end subroutine tPvd_set
  
  subroutine tPvd_write(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvd) :: this
    ! -- local
    type(tPvdDataPart), pointer :: dsp => null()
    type(tPvdDataTime), pointer :: dst => null()
    character(len=MXSLEN) :: s
    character(len=MXSLEN), dimension(9) :: sa
    integer(I4B) :: iu, it, ip
! ------------------------------------------------------------------------------
    
    call open_file(this%f, iu, 'w')
    !
    write(iu,'(a)') '<?xml version="1.0"?>'
    write(iu,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
    write(iu,'(a)') '<Collection>'
    sa(1) = '<DataSet timestep="' 
    sa(3) = '" part="'
    sa(5) = '" name="'
    sa(7) = '" group="" file="'
    sa(9) = '"/>'
    do it = 1, this%ndstime
      dst => this%dstime(it)
      sa(2) = ta([dst%timestep])
      do ip = 1, dst%npart
        dsp => dst%dspart(ip)
        sa(4) = ta([ip])
        sa(6) = trim(dsp%partname)
        sa(8) = trim(dsp%f_vtu)
        s = ta(sa,trim_sep=.true.)
        write(iu,'(a)') trim(s)
      end do
    end do
    write(iu,'(a)') '</Collection>'
    write(iu,'(a)') '</VTKFile> '
    !
    close(iu)
    !
    return
  end subroutine tPvd_write
  
! ==============================================================================
! ==============================================================================
! general
! ==============================================================================
! ==============================================================================
  
  subroutine kdtree_remove_doubles(xy, nxy, xyi)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer, intent(inout) :: nxy
    real(R8B), dimension(3, nxy), intent(inout) :: xy
    integer(I4B), dimension(:), intent(out) :: xyi
    ! -- locals
    integer(I4B), parameter :: ns = 4
    real(R4B), parameter :: eps = 1e-20
    type(kdtree2), pointer :: wrktree
    real(kdkind), dimension(:,:), allocatable :: rwrk
    integer(I4B), dimension(:), allocatable :: iwrk
    type(kdtree2_result), allocatable :: results(:)
    integer(I4B), dimension(:), allocatable :: ifnd
    integer(I4B) :: i, j, idx, idx_first, idx_new, mxy, n, j_first
! ------------------------------------------------------------------------------
    
    allocate(rwrk(3,nxy))
    do i = 1, nxy
      rwrk(1,i) = xy(1,i)
      rwrk(2,i) = xy(2,i)
      rwrk(3,i) = xy(3,i)
    end do
    wrktree => kdtree2_create(rwrk,sort=.true.,rearrange=.false.)  ! this is how you create a tree. 
    write(*,*) 'Done building kD-tree...'
    
    ! remove duplicates
    allocate(iwrk(nxy)); iwrk = 1
    allocate(results(ns),ifnd(ns))
    xyi = 0
    do i = 1, nxy
      call kdtree2_n_nearest(tp=wrktree,qv=rwrk(:,i),nn=ns,results=results)
      ! double
      ifnd = 0; j_first = 0; idx_first = 0
      do j = 1, ns
        if (results(j)%dis < eps) then
           idx = results(j)%idx 
           if (j_first == 0) then
             j_first = j
             idx_first = idx
           end if
           ifnd(j) = 1
        end if
      end do
      if (idx_first == 0) then
        call errmsg('kdtree_remove_doubles: program error 1.')
      end if
      if (sum(ifnd) > 1) then ! handle double
        do j = 1, ns
          if (ifnd(j) == 1) then 
            idx = results(j)%idx
            if ((j /= j_first).and.(iwrk(idx_first) == 1)) then
              iwrk(idx) = 0
            end if
            xyi(idx) = idx_first
          end if
        end do
      else
        if (idx_first /= i) then
          call errmsg('kdtree_remove_doubles: program error 2.')
        end if
        xyi(i) = i
      end if
    end do
    ! number of unique coordinates
    mxy = sum(iwrk)
    n = 0
    do i = 1, nxy
      if (iwrk(i) == 1) then
        n = n + 1
        iwrk(i) = n
        xy(1,n) = rwrk(1,i)
        xy(2,n) = rwrk(2,i)
        xy(3,n) = rwrk(3,i)
      end if
    end do
    !
    do i = 1, nxy
      idx = xyi(i)
      idx_new = iwrk(idx)
      if (idx_new == 0) then
        call errmsg('kdtree_remove_doubles: program error 3.')
      end if
      xyi(i) = idx_new
    end do
    !
    nxy = mxy
    deallocate(iwrk, rwrk)
    call kdtree2_destroy(wrktree) ! detroy tree  
    !
    return
  end subroutine kdtree_remove_doubles
    
end module vtk_module
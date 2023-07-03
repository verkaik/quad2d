module vtk_lib
  use utilsmod, only: I4B, R4B, R8B, mxslen, open_file, errmsg, ta
  use kdtree2_module
  
  implicit none

  private
  
  integer(I4B), parameter :: vtk_voxel = 11
  
  type :: tMesh
    integer(I4B) :: part = 0
    !
    integer(I4B) :: nv = 0 ! number of vertices
    real(R8B), dimension(:,:), allocatable :: x  ! vertex coordinates
    integer(I4B) :: nc = 0 ! number of cells
    real(R8B), dimension(:,:), allocatable :: xc ! centroid coordinates
    integer(I4B), dimension(:), allocatable :: ia ! cell --> vertices ia (nc+1)
    integer(I4B), dimension(:), allocatable :: ja ! cell --> vertices ja (nja)
    integer(I4B) :: nja
    integer(I4B), dimension(:), allocatable :: celltype ! VTK cell type 
    integer(I4B), dimension(:), allocatable :: cellmap ! mapping VTK cell(nc) --> MODFLOW cell
  contains
    procedure :: init  => tMesh_init
    procedure :: clean => tMesh_clean
    procedure :: set   => tMesh_set
  end type tMesh
  
  type :: tVtuVar
    character(len=MXSLEN) :: name
    integer(I4B), dimension(:), allocatable :: mask
    real(R8B), dimension(:), allocatable :: v
  contains
    procedure :: init  => tVtuVar_init
    procedure :: clean => tVtuVar_clean
    procedure :: set   => tVtuVar_set
  end type tVtuVar
  !
  type :: tVtu
    character(len=MXSLEN) :: f
    integer(I4B) :: nvar
    type(tVtuVar), dimension(:), pointer ::  var => null()
  contains
    procedure :: init  => tVtu_init
    procedure :: clean => tVtu_clean
    procedure :: set   => tVtu_set
    procedure :: write => tVtu_write
  end type tVtu
    
  type :: tPvdData
    real(R4B)             :: timestep = 1
    real(R4B)             :: part = 1
    character(len=MXSLEN) :: f_vtu
    type(tVtu), pointer   :: vtu => null()
  contains
    procedure :: init  => tPvdData_init
    procedure :: clean => tPvdData_clean
    procedure :: set   => tPvdData_set
  end type tPvdData
  !
  type, public :: tPvd
    character(len=MXSLEN) :: f
    integer(I4B) :: npartmesh
    integer(I4B) :: ndataset
    type(tMesh), dimension(:), pointer    :: partmesh => null() ! length npart
    type(tPvdData), dimension(:), pointer :: dataset  => null() ! length ndataset
  contains
    procedure :: init  => tPvd_init
    procedure :: clean => tPvd_clean
  end type tPvd
  
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
    return
  end subroutine tMesh_clean
  
  subroutine tMesh_set(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMesh) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tMesh_set
  
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
  
  subroutine tVtu_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tVtu_init

  subroutine tVtu_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tVtu_clean
  
  subroutine tVtu_set(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tVtu_set
  
  subroutine tVtu_write(this, mesh)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVtu) :: this
    type(tMesh), intent(in) :: mesh
    !
    ! -- local
    type(tVtuVar), pointer :: var => null()
    integer(I4B) :: iu, i, j
    character(len=MXSLEN) :: s, s1, s2
! ------------------------------------------------------------------------------
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
    ! points
    do i = 1, mesh%nv
      write(iu,'(a,1x,a,1x,a)') ta(mesh%x(:,i))
    end do
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '</Points>'
    write(iu,'(a)') '<Cells>'
    write(iu,'(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
    ! connectivity
    write(iu,*) (mesh%ja(i)-1,i=1,mesh%nja)
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
    ! offsets
    write(iu,*) (mesh%ia(i)-1,i=2,mesh%nc+1)
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '<DataArray type="Int32" Name="types" format="ascii">'
    ! types
    write(iu,*) (mesh%celltype(i),i=1,mesh%nc)
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '</Cells>'
    write(iu,'(a)') '<CellData>'
    ! variables
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
  end subroutine tVtu_write
  
! ==============================================================================
! ==============================================================================
! tPvdData
! ==============================================================================
! ==============================================================================
  
  subroutine tPvdData_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvdData) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tPvdData_init

  subroutine tPvdData_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvdData) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tPvdData_clean
  
  subroutine tPvdData_set(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvdData) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tPvdData_set
 
  
! ==============================================================================
! ==============================================================================
! tPvd
! ==============================================================================
! ==============================================================================
  
  subroutine tPvd_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvd) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tPvd_init

  subroutine tPvd_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvd) :: this
! ------------------------------------------------------------------------------
    return
  end subroutine tPvd_clean
  
  subroutine tPvd_write(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tPvd) :: this
  
    return
  end subroutine tPvd_write
  
! ==============================================================================
! ==============================================================================
! general
! ==============================================================================
! ==============================================================================
  
  subroutine kdtree_remove_doubles(xy, nxy)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer, intent(inout) :: nxy
    real(r8b), dimension(3, nxy), intent(inout) :: xy
    ! -- locals
    integer, parameter :: ns = 4
    real(R4B), parameter :: eps = 1e-20
    type(kdtree2), pointer :: wrktree
    real(kdkind), dimension(:,:), allocatable :: rwrk
    integer, dimension(:), allocatable :: iwrk
    type(kdtree2_result), allocatable :: results(:)
    integer(I4B), dimension(:), allocatable :: ifnd
    integer(I4B) :: i, j, idx, mxy, n, jfirst
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
    do i = 1, nxy
      call kdtree2_n_nearest(tp=wrktree,qv=rwrk(:,i),nn=ns,results=results)
      ! double
      ifnd = 0; jfirst = 0
      do j = 1, ns
        if (results(j)%dis < eps) then
           if (jfirst.eq.0) jfirst = j
           idx = results(j)%idx 
           ifnd(j) = 1         
        end if   
      end do
      if (sum(ifnd) > 1) then ! handle double
         idx = results(jfirst)%idx ! first found
         if (iwrk(idx).eq.1) then
           do j = 1, ns
             if (ifnd(j) == 1) then 
               if (j /= jfirst) then
                 idx = results(j)%idx
                 iwrk(idx) = 0
               end if
             end if
           end do
         end if
      end if
    end do
    ! number of unique coordinates
    mxy = sum(iwrk)
      
    n = 0
    do i = 1, nxy
      if (iwrk(i) == 1) then
         n = n + 1
         xy(1,n) = rwrk(1,i)
         xy(2,n) = rwrk(2,i)
         xy(3,n) = rwrk(3,i)
      end if
    end do
    nxy = mxy
    deallocate(iwrk, rwrk)
    call kdtree2_destroy(wrktree) ! detroy tree  
    
  end subroutine kdtree_remove_doubles
    
end module vtk_lib
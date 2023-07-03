module utils
  use, intrinsic :: iso_fortran_env , only: error_unit, output_unit, &
    i4b => int32, r4b => real32, r8b => real64
  
  implicit none
  
  integer(i4b), parameter :: mxslen = 1000
  integer(i4b), parameter :: idis = 1, idisv = 2, idisu = 3

  integer, parameter :: vtk_voxel = 11
  integer, parameter :: vtk_wedge = 13
  
  type tMesh
    integer(i4b) :: nv = 0 ! number of vertices
    real(r8b), dimension(:,:), allocatable :: x  ! vertex coordinates
    integer(i4b) :: nc = 0 ! number of cells
    real(r8b), dimension(:,:), allocatable :: xc ! centroid coordinates
    integer(i4b), dimension(:), allocatable :: ia ! cell --> vertices ia (nc+1)
    integer(i4b), dimension(:), allocatable :: ja ! cell --> vertices ja (nja)
    integer(i4b) :: nja
    integer(i4b), dimension(:), allocatable :: celltype ! VTK cell type 
    integer(i4b), dimension(:), allocatable :: cellmap ! mapping VTK cell(nc) --> MODFLOW cell
  end type tMesh
  
  type tDis
    type(tMesh), pointer :: mesh => null()
  contains
    procedure :: read => grb_read_dis
  end type tDis

  type tDisv
    type(tMesh), pointer :: mesh => null()
  contains
    procedure :: read => grb_read_disv
  end type tDisv
  
  type tGrb
    integer(i4b) :: iu
    integer(i4b) :: itype = 0
    integer(i4b) :: version = 0
    integer(i4b) :: ntxt = 0
    integer(i4b) :: lentxt = 0
    type(tDis),  pointer :: dis => null()
    type(tDisv), pointer :: disv => null()
!    type(tDisu), pointer :: disu => null()
  contains
    procedure :: init => grb_init
    procedure :: read => grb_read
    procedure :: get  => grb_get
  end type tGrb
  
  type tMf6
    type(tGrb), pointer :: grb => null()
  contains
    procedure :: init => mf6_init
    procedure :: write => mf6_write_vtu
  end type tMf6
  
contains

! ******************************************************************************
  subroutine mf6_init(this, f)
! ******************************************************************************
    ! -- dummy
    class(tMf6) :: this
    character(len=*), intent(in) :: f
! ------------------------------------------------------------------------------
    
    allocate(this%grb)
    call this%grb%init(f)
    call this%grb%read()
  
  end subroutine mf6_init
  
! ******************************************************************************
  subroutine mf6_write_vtu(this, vtk_f, hds_f, var, kper, kstp, vtk_type)
! ******************************************************************************
    ! -- dummy
    class(tMf6) :: this
    character(len=*), intent(in) :: hds_f
    character(len=*), intent(in) :: vtk_f
    character(len=*), intent(in) :: var
    integer(i4b), intent(in) :: kper, kstp
    character(len=*), intent(in) :: vtk_type
    ! -- local
    type(tMesh), pointer :: mesh
    real(r8b), dimension(:), allocatable :: v
! ------------------------------------------------------------------------------
  
    mesh => this%grb%get()
    call getvar_ulasav(v, mesh%cellmap, hds_f, var, kper, kstp)
    select case(vtk_type)
    case('asc')
      call write_vtu_asc(vtk_f, mesh, var, v)
    case('bin_b64')
      call write_vtu_bin_b64(vtk_f, mesh, var, v)
    case('bin_raw')
      call write_vtu_bin_raw(vtk_f, mesh, var, v)
    end select
    deallocate(v)
    
  end subroutine mf6_write_vtu
  
! ******************************************************************************
  subroutine getvar_ulasav(v, cellmap, f, var, kper, kstp)
! ******************************************************************************
    ! -- dummy
    real(r8b), dimension(:), allocatable, intent(inout) :: v
    integer(i4b), dimension(:), intent(in) :: cellmap
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: f
    integer(i4b), intent(in) :: kper, kstp
    ! -- local
    integer(i4b) :: n, iu, p, pp, ios, i, icol, irow
    !
    ! ulasav
    character(len=16) :: rtext
    integer(i4b) ::  rkstp, rkper, rncol, rnrow, rilay
    real(r8b) ::  rpertim, rtotim
    real(r8b), dimension(:,:), allocatable :: wrk
    
! ------------------------------------------------------------------------------
    n = size(cellmap)
    allocate(v(n))
    v = 0.d0
    
    open(newunit=iu, file=f, form='unformatted', access='stream', status='old')
    p = 1
    do while(1)
      !        4      4      8        8       16     4      4      4 = 52 bytes 
      read(iu, pos=p,iostat=ios) rkstp, rkper, rpertim, rtotim, rtext, rncol, rnrow, rilay
      p = p + 52
      !allocate(wrk(rncol,rnrow))
      !read(iu)((wrk(icol,irow),icol=1,rncol),irow=1,rnrow)
      if (ios /= 0) exit
      if ((rkstp == kstp).and.(rkper == kper).and.(trim(rtext) == trim(var))) then
        do i = 1, n
          pp = p + 8*(cellmap(i)-1)
          read(iu,pos=pp) v(i)
        end do
      end if
      p = p + 8*rncol*rnrow
    end do
    close(iu)
  
  end subroutine getvar_ulasav
  
! ******************************************************************************
  subroutine write_vtu_asc(vtk_f, mesh, var, v)
! ******************************************************************************
    ! -- dummy
    character(len=*), intent(in) :: vtk_f
    type(tMesh), pointer :: mesh
    character(len=*) :: var
    real(r8b), dimension(:), intent(in) :: v
    ! -- local
    integer(i4b) :: iu, i, j
    character(len=mxslen) :: s
    character(len=mxslen), dimension(3) :: sa
! ------------------------------------------------------------------------------
    
    open(newunit=iu, file=vtk_f, status='replace')
    
    write(iu,'(a)') '<?xml version="1.0"?>'
    write(iu,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(iu,'(a)') '<UnstructuredGrid>'
    write(sa(1),*) mesh%nv
    write(sa(2),*) mesh%nc
    write(iu,'(5a)') '<Piece NumberOfPoints="',trim(adjustl(sa(1))),'" NumberOfCells="',trim(adjustl(sa(2))),'">'
    write(iu,'(a)') '<Points>'
    write(iu,'(a)') '<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="ascii">'
    ! points
    do i = 1, mesh%nv
      write(sa(1),*) mesh%x(1,i)
      write(sa(2),*) mesh%x(2,i)
      write(sa(3),*) mesh%x(3,i)
      write(iu,'(a,1x,a,1x,a)')(trim(adjustl(sa(j))),j=1,3)
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
    write(s,'(3a)')  '<DataArray type="Float32" Name="',trim(var),'" NumberOfComponents="1" format="ascii">'
    write(iu,'(a)') trim(s)
    ! variable
    write(iu,*) v
    write(iu,'(a)') '</DataArray>'
    write(iu,'(a)') '</CellData>'
    write(iu,'(a)') '</Piece>'
    write(iu,'(a)') '</UnstructuredGrid>'
    write(iu,'(a)') '</VTKFile>'
    close(iu)
  
  end subroutine write_vtu_asc

! ******************************************************************************
  subroutine write_vtu_bin_b64(vtk_f, mesh, var, v)
! ******************************************************************************
    ! -- modules
    use penf
    use befor64
    ! -- dummy
    character(len=*), intent(in) :: vtk_f
    type(tMesh), pointer :: mesh
    character(len=*) :: var
    real(r8b), dimension(:), intent(in) :: v
    ! -- local
    integer(i4b) :: iu, i, j, n, i1, i2
    character(len=mxslen) :: s
    character(len=mxslen), dimension(3) :: sa
    character(1), parameter:: endrec = char(10) !< End-character for binary-record finalize. 
    !
    character(len=:), allocatable:: code64 ! base64 encoded string
    integer(i8p), dimension(:), allocatable :: iwrk
    real(r8p), dimension(:), allocatable :: dwrk
    integer(I1P), dimension(:), allocatable :: varp
! ------------------------------------------------------------------------------
    
    call b64_init()
    open(newunit=iu, file=vtk_f, action='write', status='replace', form='unformatted', access='stream')
    
    write(iu) '<?xml version="1.0"?>'//endrec
    write(iu) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//endrec
    write(iu) '<UnstructuredGrid>'//endrec
    write(sa(1),*) mesh%nv
    write(sa(2),*) mesh%nc
    write(s,'(5a)') '<Piece NumberOfPoints="',trim(adjustl(sa(1))),'" NumberOfCells="',trim(adjustl(sa(2))),'">'
    write(iu) trim(s)//endrec
    write(iu) '<Points>'//endrec
    write(iu) '<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="binary">'//endrec
    !-------------------
    ! points
    !-------------------
    allocate(dwrk(mesh%nv*3))
    n = 0
    do i = 1, mesh%nv
      n = n + 1; dwrk(n) = mesh%x(1,i)
      n = n + 1; dwrk(n) = mesh%x(2,i)
      n = n + 1; dwrk(n) = mesh%x(3,i)
    end do
    call pack_data(a1=[int(3*mesh%nv*BYR8P,I4P)],a2=dwrk,packed=varp)
    call b64_encode(n=varp,code=code64)
    write(iu) trim(code64)//endrec
    deallocate(dwrk,varp,code64)
    !-------------------
    write(iu) '</DataArray>'//endrec
    write(iu) '</Points>'//endrec
    write(iu) '<Cells>'//endrec
    write(iu) '<DataArray type="Int64" Name="connectivity" format="binary">'//endrec
    !-------------------
    ! connectivity
    !-------------------
    allocate(iwrk(mesh%nja))
    do i = 1, mesh%nja
      iwrk(i) = mesh%ja(i) - 1 
    end do
    call pack_data(a1=[int(mesh%nja*BYI8P,I4P)],a2=iwrk,packed=varp)
    call b64_encode(n=varp,code=code64)
    write(iu) code64//endrec
    deallocate(iwrk,varp,code64)
    !-------------------
    write(iu) '</DataArray>'//endrec
    write(iu) '<DataArray type="Int64" Name="offsets" format="binary">'//endrec
    !-------------------
    ! offsets
    !-------------------
    allocate(iwrk(mesh%nc))
    do i = 2, mesh%nc+1
      iwrk(i-1) = mesh%ia(i) - 1 
    end do
    call pack_data(a1=[int(mesh%nc*BYI8P,I4P)],a2=iwrk,packed=varp)
    call b64_encode(n=varp,code=code64)
    write(iu) code64//endrec
    deallocate(iwrk,varp,code64)
    !-------------------
    write(iu) '</DataArray>'//endrec
    write(iu) '<DataArray type="Int64" Name="types" format="binary">'//endrec
    !-------------------
    ! types
    !-------------------
    allocate(iwrk(mesh%nc))
    do i = 1, mesh%nc
      iwrk(i) = mesh%celltype(i)
    end do
    call pack_data(a1=[int(mesh%nc*BYI8P,I4P)],a2=iwrk,packed=varp)
    call b64_encode(n=varp,code=code64)
    write(iu) code64//endrec
    deallocate(iwrk,varp,code64)
    !-------------------
    write(iu) '</DataArray>'//endrec
    write(iu) '</Cells>'//endrec
    write(iu) '<CellData>'//endrec
    write(s,'(3a)') '<DataArray type="Float64" Name="',trim(var),'" NumberOfComponents="1" format="binary">'
    write(iu) trim(s)//endrec
    !-------------------
    ! variable
    !-------------------
    allocate(dwrk(mesh%nc))
    do i = 1, mesh%nc
      dwrk(i) = v(i)
    end do
    call pack_data(a1=[int(mesh%nc*BYR8P,I4P)],a2=dwrk,packed=varp)
    call b64_encode(n=varp,code=code64)
    write(iu) code64//endrec
    deallocate(dwrk,varp,code64)
    !-------------------
    write(iu) '</DataArray>'//endrec
    write(iu) '</CellData>'//endrec
    write(iu) '</Piece>'//endrec
    write(iu) '</UnstructuredGrid>'//endrec
    write(iu) '</VTKFile>'
    close(iu)
  
  end subroutine write_vtu_bin_b64
  
! ******************************************************************************
  subroutine write_vtu_bin_raw(vtk_f, mesh, var, v)
! ******************************************************************************
    ! -- modules
    use penf
    use befor64
    ! -- dummy
    character(len=*), intent(in) :: vtk_f
    type(tMesh), pointer :: mesh
    character(len=*) :: var
    real(r8b), dimension(:), intent(in) :: v
    ! -- local
    integer, parameter :: nstp = 1000
    integer(i4b) :: iu, i, j, n, i1, i2, offset
    character(len=mxslen) :: s
    character(len=mxslen), dimension(3) :: sa
    character(1), parameter:: endrec = char(10) !< End-character for binary-record finalize. 
    !
    character(len=:), allocatable:: code64 ! base64 encoded string
    integer(i8p), dimension(:), allocatable :: iwrk
    real(r8p), dimension(:), allocatable :: dwrk
    integer(I1P), dimension(:), allocatable :: varp
! ------------------------------------------------------------------------------
    
    call b64_init()    
    open(newunit=iu, file=vtk_f, action='write', status='replace', form='unformatted', access='stream')
    
    write(iu) '<?xml version="1.0"?>'//endrec
    write(iu) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//endrec
    write(iu) '<UnstructuredGrid>'//endrec
    write(sa(1),*) mesh%nv; write(sa(2),*) mesh%nc
    write(iu) '<Piece NumberOfPoints="'//trim(adjustl(sa(1)))//'" NumberOfCells="',trim(adjustl(sa(2)))//'">'//endrec
    write(iu) '<Points>'//endrec
    offset = 0; write(sa(1),*) offset
    write(iu) '<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="append" offset="'//trim(adjustl(sa(1)))//'"/>'//endrec
    offset = offset + 3*mesh%nv*byr8p + byi4p; write(sa(1),*) offset
    write(iu) '</Points>'//endrec
    write(iu) '<Cells>'//endrec
    write(iu) '<DataArray type="Int32" Name="connectivity" format="append" offset="'//trim(adjustl(sa(1)))//'"/>'//endrec
    offset = offset + mesh%nja*byi4p + byi4p; write(sa(1),*) offset
    write(iu) '<DataArray type="Int32" Name="offsets" format="append" offset="'//trim(adjustl(sa(1)))//'"/>'//endrec
    offset = offset + mesh%nc*byi4p + byi4p; write(sa(1),*) offset
    write(iu) '<DataArray type="Int32" Name="types" format="append" offset="'//trim(adjustl(sa(1)))//'"/>'//endrec
    offset = offset + mesh%nc*byi4p + byi4p; write(sa(1),*) offset
    write(iu) '</Cells>'//endrec
    write(iu) '<CellData>'//endrec
    
    do i=1, nstp
      write(sa(2),'(a,i4.4)') '_',i
      write(iu) '<DataArray type="Float64" Name="'//trim(var)//trim(adjustl(sa(2)))//'" NumberOfComponents="1" format="append" offset="'//trim(adjustl(sa(1)))//'"/>'//endrec
      offset = offset + mesh%nc*r8p + byi4p; write(sa(1),*) offset    
    end do
     
    write(iu) '</CellData>'//endrec
    write(iu) '</Piece>'//endrec
    write(iu) '</UnstructuredGrid>'//endrec
    write(iu) '<AppendedData encoding="raw">'//endrec
    !-------------------
    write(iu) '_'
    write(iu) int(3*mesh%nv*byr8p,i4p), ((mesh%x(i,j),i=1,3),j=1,mesh%nv)
    write(iu) int(mesh%nja*byi4p,i4p), (mesh%ja(i)-1,i=1,mesh%nja)
    write(iu) int(mesh%nc*byi4p,i4p), (mesh%ia(i)-1,i=2,mesh%nc+1)
    write(iu) int(mesh%nc*byi4p,i4p),(mesh%celltype(i),i=1,mesh%nc)
    
    do j=1, nstp
      write(iu) int(mesh%nc*byr8p,i4p),(j*v(i),i=1,mesh%nc)
    enddo
    !-------------------
    write(iu) endrec//'</AppendedData>'//endrec
    write(iu) '</VTKFile>'
    close(iu)
  
  end subroutine write_vtu_bin_raw
  
! ******************************************************************************
  subroutine grb_init(this, f)
! ******************************************************************************
    ! -- dummy
    class(tGrb) :: this
    character(len=*), intent(in) :: f
    ! --local
    integer(i4b) :: iu, i
    character(len=50), dimension(4) :: hdr
    character(len=4) :: s
! ------------------------------------------------------------------------------
    open(newunit=this%iu, file=f, form='unformatted', access='stream', status='old')
    do i = 1, 4
      read(this%iu) hdr(i)
    end do
    read(hdr(1)(5:),*) s
    select case(s)
      case('DIS')
        this%itype = idis
      case('DISV')
        this%itype = idisv
      case('DISU')
        this%itype = idisu
    end select
    
    read(hdr(2)(8:),*) this%version
    read(hdr(3)(5:),*) this%ntxt
    read(hdr(4)(7:),*) this%lentxt
      
    rewind(this%iu)
      
  end subroutine grb_init
  
! ******************************************************************************
  subroutine grb_read(this)
! ******************************************************************************
    ! -- dummy
    class(tGrb) :: this
! ------------------------------------------------------------------------------
    select case(this%itype)
      case(idis)
        allocate(this%dis)
        call this%dis%read(this%iu, this%ntxt, this%lentxt)
      case(idisv)
        allocate(this%disv)
        call this%disv%read(this%iu, this%ntxt, this%lentxt)
      case(idisu)
        write(*,*) '*** DISU is not yet supported ***'
        stop 1
    end select
  
  end subroutine grb_read
  
! ******************************************************************************
  function grb_get(this) result(mesh)
! ******************************************************************************
    ! -- dummy
    class(tGrb) :: this
    type(tMesh), pointer :: mesh
! ------------------------------------------------------------------------------
    select case(this%itype)
      case(idis)
         mesh => this%dis%mesh  
      case(idisv)
        write(*,*) '*** DISV is not yet supported ***'
        stop 1
      case(idisu)
        write(*,*) '*** DISV is not yet supported ***'
        stop 1
    end select
  
  end function grb_get
  
! ******************************************************************************
  subroutine grb_read_dis(this, iu, ntxt, lentxt)
! ******************************************************************************
    use kdtree2_module
    ! -- dummy
    class(tDis) :: this
    integer(i4b), intent(in) :: iu, ntxt, lentxt
    ! -- local
    integer(i4b) :: is, i, icol, irow, ilay, n, iact, nmf
    real(r8b), dimension(:,:), allocatable :: wrk
    real(r8b), dimension(:), allocatable :: cdelr, cdelc
    real(r8b) :: yupper, x1, x2, y1, y2, z1, z2
    type(tMesh), pointer :: mesh
    !
    integer(i4b) :: ncells, nlay, nrow, ncol, nja
    real(r8b) :: xorigin, yorigin
    real(r8b), dimension(:), allocatable :: delr, delc
    real(r8b), dimension(:,:), allocatable :: top
    real(r8b), dimension(:,:,:), allocatable :: botm
    integer(i4b), dimension(:,:,:), allocatable :: idomain
    
    integer, parameter :: ns = 1
    real, parameter :: eps = 1e-20
    type(kdtree2), pointer :: wrktree
    type(kdtree2_result), allocatable :: results(:)
! ------------------------------------------------------------------------------
  
    ! read data
    is = 1 + 200 + ntxt*lentxt
    i = is
    read(unit=iu,pos=i) ncells
    i = is + 1*i4b
    read(unit=iu,pos=i) nlay
    i = is + 2*i4b
    read(unit=iu,pos=i) nrow
    i = is + 3*i4b
    read(unit=iu,pos=i) ncol
    i = is + 4*i4b
    read(unit=iu,pos=i) nja
    i = is + 5*i4b
    read(unit=iu,pos=i) xorigin
    i = is + 5*i4b + 1*r8b
    read(unit=iu,pos=i) yorigin
    allocate(delr(ncol))
    i = is + 5*i4b + 3*r8b
    read(unit=iu,pos=i) (delr(icol),icol=1,ncol)
    allocate(delc(nrow))
    i = is + 5*i4b + 3*r8b + nrow*r8b
    read(unit=iu,pos=i) (delc(irow),irow=1,nrow)
    allocate(top(ncol,nrow))
    i = is + 5*i4b + 3*r8b + nrow*r8b + ncol*r8b
    read(unit=iu,pos=i) ((top(icol,irow),icol=1,ncol),irow=1,nrow)
    allocate(botm(ncol,nrow,nlay))
    i = is + 5*i4b + 3*r8b + nrow*r8b + ncol*r8b + ncol*nrow*r8b
    read(unit=iu,pos=i) (((botm(icol,irow,ilay),icol=1,ncol),irow=1,nrow),ilay=1,nlay)
    allocate(idomain(ncol,nrow,nlay))
    i = is + 5*i4b + 3*r8b + nrow*r8b + ncol*r8b + ncol*nrow*r8b + ncells*r8b + i4b*(ncells+1) + i4b*nja
    read(unit=iu,pos=i) (((idomain(icol,irow,ilay),icol=1,ncol),irow=1,nrow),ilay=1,nlay)
    close(iu)
    !
    ! ---- DEBUG
    !nlay = 1
    !nrow = 2
    !ncol = 2
    !ncells = nlay*nrow*ncol
    ! ---- DEBUG
    
    allocate(cdelr(ncol+1))
    cdelr = 0.
    do icol = 1, ncol
      cdelr(icol+1) = cdelr(icol)+delr(icol)
    end do
    allocate(cdelc(nrow+1))
    cdelc = 0.
    do irow = 1, nrow
      cdelc(irow+1) = cdelc(irow)+delc(irow)
    end do
    yupper = yorigin + cdelc(nrow+1)
    
    allocate(wrk(3,8*ncells))
    allocate(this%mesh)
    mesh => this%mesh
    do iact = 1, 2
      n = 0
      nmf = 0
      mesh%nc = 0
      do ilay = 1, nlay
        do irow = 1, nrow
          do icol = 1, ncol
            nmf = nmf + 1
            if (idomain(icol,irow,ilay) /= 0) then
              mesh%nc = mesh%nc + 1
              x1 = xorigin + cdelr(icol)
              x2 = x1 + delr(icol)
              y1 = yupper - cdelc(irow)
              y2 = y1 - delc(irow)
              if (ilay == 1) then
                z1 = top(icol,irow)
              else
                z1 = botm(icol,irow,ilay-1)
              end if
              z2 = botm(icol,irow,ilay)
              if (iact == 1) then
                n = n + 1; wrk(1,n) = x1; wrk(2,n) = y2; wrk(3,n) = z2 !0
                n = n + 1; wrk(1,n) = x2; wrk(2,n) = y2; wrk(3,n) = z2 !1
                n = n + 1; wrk(1,n) = x1; wrk(2,n) = y1; wrk(3,n) = z2 !2
                n = n + 1; wrk(1,n) = x2; wrk(2,n) = y1; wrk(3,n) = z2 !3
                n = n + 1; wrk(1,n) = x1; wrk(2,n) = y2; wrk(3,n) = z1 !4
                n = n + 1; wrk(1,n) = x2; wrk(2,n) = y2; wrk(3,n) = z1 !5
                n = n + 1; wrk(1,n) = x1; wrk(2,n) = y1; wrk(3,n) = z1 !6
                n = n + 1; wrk(1,n) = x2; wrk(2,n) = y1; wrk(3,n) = z1 !7
              else
                mesh%cellmap(mesh%nc) = nmf
                mesh%ia(mesh%nc+1) = mesh%ia(mesh%nc) + 8
                call kdtree2_n_nearest(tp=wrktree,qv=(/x1,y2,z2/),nn=ns,results=results); n = n + 1; mesh%ja(n) = results(1)%idx
                call kdtree2_n_nearest(tp=wrktree,qv=(/x2,y2,z2/),nn=ns,results=results); n = n + 1; mesh%ja(n) = results(1)%idx
                call kdtree2_n_nearest(tp=wrktree,qv=(/x1,y1,z2/),nn=ns,results=results); n = n + 1; mesh%ja(n) = results(1)%idx
                call kdtree2_n_nearest(tp=wrktree,qv=(/x2,y1,z2/),nn=ns,results=results); n = n + 1; mesh%ja(n) = results(1)%idx
                call kdtree2_n_nearest(tp=wrktree,qv=(/x1,y2,z1/),nn=ns,results=results); n = n + 1; mesh%ja(n) = results(1)%idx
                call kdtree2_n_nearest(tp=wrktree,qv=(/x2,y2,z1/),nn=ns,results=results); n = n + 1; mesh%ja(n) = results(1)%idx
                call kdtree2_n_nearest(tp=wrktree,qv=(/x1,y1,z1/),nn=ns,results=results); n = n + 1; mesh%ja(n) = results(1)%idx
                call kdtree2_n_nearest(tp=wrktree,qv=(/x2,y1,z1/),nn=ns,results=results); n = n + 1; mesh%ja(n) = results(1)%idx
              end if 
            end if
          end do
        end do
      end do
      if (iact == 1) then
        call kdtree_remove_doubles(wrk, n)
        mesh%nv = n
        allocate(mesh%x(3,mesh%nv))
        do i = 1, mesh%nv
          mesh%x(1,i) = wrk(1,i)
          mesh%x(2,i) = wrk(2,i)
          mesh%x(3,i) = wrk(3,i)
        end do
        deallocate(wrk)
        wrktree => kdtree2_create(mesh%x,sort=.true.,rearrange=.false.)
        allocate(results(ns)) 
        
        allocate(mesh%celltype(mesh%nc))
        mesh%celltype = vtk_voxel
        allocate(mesh%cellmap(mesh%nc))
        allocate(mesh%ia(mesh%nc+1))
        mesh%ia(1) = 1
        mesh%nja = 8*mesh%nc
        allocate(mesh%ja(mesh%nja))
      end if
    end do
    ! determine mesh
  end subroutine grb_read_dis

! ******************************************************************************
  subroutine grb_read_disv(this, iu, ntxt, lentxt)
! ******************************************************************************
    use kdtree2_module
    ! -- dummy
    class(tDisv) :: this
    integer(i4b), intent(in) :: iu, ntxt, lentxt
    ! -- local
    integer(i4b) :: is, i, j, k, ipos, n, ivert, nvcell, iact, nmf, nv, mwrk
    real(r8b) :: x, y, z1, z2
    real(r8b), dimension(:,:), allocatable :: wrk
    type(tMesh), pointer :: mesh
    !
    integer(i4b) :: ncells, nlay, ncpl, nvert, njavert, nja
    real(r8b) :: xorigin, yorigin
    real(r8b), dimension(:), allocatable :: top
    real(r8b), dimension(:), allocatable :: botm
    real(r8b), dimension(:,:), allocatable :: vertices
    real(r8b), dimension(:), allocatable :: cellx
    real(r8b), dimension(:), allocatable :: celly
    integer(i4b), dimension(:), allocatable :: iavert
    integer(i4b), dimension(:), allocatable :: javert
    integer(i4b), dimension(:), allocatable :: idomain
    
    integer, parameter :: ns = 1
    real, parameter :: eps = 1e-20
    type(kdtree2), pointer :: wrktree
    type(kdtree2_result), allocatable :: results(:)
! ------------------------------------------------------------------------------
  
    write(*,*) '@@@1'
    
    ! read data
    is = 1 + 200 + ntxt*lentxt
    i = is
    read(unit=iu,pos=i) ncells
    i = is + 1*i4b
    read(unit=iu,pos=i) nlay
    i = is + 2*i4b
    read(unit=iu,pos=i) ncpl
    i = is + 3*i4b
    read(unit=iu,pos=i) nvert
    i = is + 4*i4b
    read(unit=iu,pos=i) njavert
    i = is + 5*i4b
    read(unit=iu,pos=i) nja
    i = is + 6*i4b
    read(unit=iu,pos=i) xorigin
    i = is + 6*i4b + 1*r8b
    read(unit=iu,pos=i) yorigin
    
    write(*,*) '@@@2'
    
    allocate(top(ncpl))
    i = is + 6*i4b + 3*r8b
    read(unit=iu,pos=i) (top(i),i=1,ncpl)
    
    allocate(botm(ncells))
    i = is + 6*i4b + 3*r8b + 8*ncpl
    read(unit=iu,pos=i) (botm(i),i=1,ncells)
    
    allocate(vertices(2,nvert))
    i = is + 6*i4b + 3*r8b + 8*ncpl + 8*ncells
    read(unit=iu,pos=i) ((vertices(j,k),j=1,2),k=1,nvert)
    
    allocate(cellx(ncells))
    i = is + 6*i4b + 3*r8b + 8*ncpl + 8*ncells + 2*8*nvert
    read(unit=iu,pos=i) (cellx(i),i=1,ncells)

    write(*,*) '@@@3'
    
    allocate(celly(ncells))
    i = is + 6*i4b + 3*r8b + 8*ncpl + 8*ncells + 2*8*nvert + 8*ncpl
    read(unit=iu,pos=i) (celly(i),i=1,ncells)
    
    allocate(iavert(ncpl+1))
    i = is + 6*i4b + 3*r8b + 8*ncpl + 8*ncells + 2*8*nvert + 8*ncpl + 8*ncpl
    read(unit=iu,pos=i) (iavert(i),i=1,ncpl+1)
    
    write(*,*) '@@@4'
    allocate(javert(njavert))
    i = is + 6*i4b + 3*r8b + 8*ncpl + 8*ncells + 2*8*nvert + 8*ncpl + 8*ncpl + 4*(ncpl+1)
    write(*,*) '@@@-->',njavert
    
    !read(unit=iu,pos=i) (javert(i),i=1,njavert)
    write(*,*) '@@@5'

    allocate(idomain(ncells))
!    i = is + 6*i4b + 3*r8b + 8r8b*ncpl + r8b*ncells + 2*r8b*nvert + r8b*ncpl + r8b*ncpl + i4b*(ncpl+1) + i4b*njavert + i4b*(ncells+1) + i4b*nja
    read(unit=iu,pos=i) (idomain(i),i=1,ncells)

    write(*,*) '!!!!!!!!!!!!!!'
    
    allocate(this%mesh)
    mesh => this%mesh
    mwrk = 4*nlay*ncpl*nvert*2
    allocate(wrk(3,mwrk))
    do iact = 1, 2
      nmf = 0
      mesh%nc = 0
      nv = 0
      do k = 1, nlay
        do n = 1, ncpl
          nmf = nmf + 1
          if (idomain(nmf) /= 0) then
            mesh%nc = mesh%nc + 1
            nvcell = iavert(n+1) - iavert(n) - 1
            if (iact == 2) then
              select case(nvcell)
                case(3)
                  mesh%celltype(mesh%nc) = vtk_wedge
                case(4)
                  mesh%celltype(mesh%nc) = vtk_voxel
                case default
                   write(*,*) 'Error, Not supported cell type.'
                   stop 1
              end select
            end if
            
            do ipos =  iavert(n + 1) - 1, iavert(n) + 1, -1
              ivert = javert(ipos)
              x = vertices(1,ivert)
              y = vertices(2,ivert)
              if (k == 1) then 
                z1 = top(nmf)
                z2 = botm(nmf)
              else
                z1 = botm(nmf)
                z2 = botm(nmf-ncpl+1)
              end if
              if (iact == 1) then
                nv = nv + 1
                if (nv > mwrk) then
                   write(*,*) 'Error, out of bound.'
                   stop 1
                end if
                wrk(1,nv) = x; wrk(2,nv) = y; wrk(3,nv) = z2
                nv = nv + 1
                if (nv > mwrk) then
                   write(*,*) 'Error, out of bound.'
                   stop 1
                end if
                wrk(1,nv) = x; wrk(2,nv) = y; wrk(3,nv) = z2
              else
                
              end if
            enddo
          end if
        end do
      end do
      if (iact == 1) then
        allocate(mesh%celltype(mesh%nc))
        write(*,*) '@@@@',nv
        call kdtree_remove_doubles(wrk(:,1:nv), nv)
        write(*,*) '@@@@',nv
        mesh%nv = nv
      end if
    
    end do
   
    
    
    
    
    

    ! determine mesh
  end subroutine grb_read_disv
  
! ******************************************************************************
  subroutine kdtree_remove_doubles(xy, nxy)
! ******************************************************************************
    use kdtree2_module
    ! arguments
    integer, intent(inout) :: nxy
    real(r8b), dimension(3, nxy), intent(inout) :: xy
    ! locals
    integer, parameter :: ns = 4
    real, parameter :: eps = 1e-20
    type(kdtree2), pointer :: wrktree
    real(kdkind), dimension(:,:), allocatable :: rwrk
    integer, dimension(:), allocatable :: iwrk
    type(kdtree2_result), allocatable :: results(:)
    integer, dimension(:), allocatable :: ifnd
    integer :: i, j, idx, mxy, n, jfirst
    ! body
    
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
  
end module utils
  
! ******************************************************************************
program mf6tovtk
use utils
implicit none

! -- local
character(len=mxslen) :: grb_f, hds_f, vtk_f
integer(i4b) :: kper, kstp
type(tMf6), pointer :: m
! ------------------------------------------------------------------------------

write(*,*) '@@@Start!'

grb_f = 'd:\codes\develop\modflow6-examples\mf6\test031_many_gwf\m1_1\model.dis.grb'
hds_f = 'd:\codes\develop\modflow6-examples\mf6\test031_many_gwf\m1_1.hds'
vtk_f = 'd:\codes\develop\modflow6-examples\mf6\test031_many_gwf\test.vtu'

!grb_f = 'd:\codes\develop\modflow6\trunk\examples\ex11-disvmesh\ci.disv.grb'
!hds_f = 'd:\codes\develop\modflow6\trunk\examples\ex11-disvmesh\ci.output.hds'

allocate(m)
call m%init(grb_f)

kper = 1; kstp = 1
vtk_f = 'd:\codes\develop\modflow6-examples\mf6\test031_many_gwf\test_asc.vtu'
!call m%write( vtk_f, hds_f, 'HEAD', kper, kstp, 'asc')
vtk_f = 'd:\codes\develop\modflow6-examples\mf6\test031_many_gwf\test_bin_b64.vtu'
!call m%write( vtk_f, hds_f, 'HEAD', kper, kstp, 'bin_b64')
vtk_f = 'd:\codes\develop\modflow6-examples\mf6\test031_many_gwf\test_bin_raw.vtu'
call m%write( vtk_f, hds_f, 'HEAD', kper, kstp, 'bin_raw')

end program mf6tovtk
  
  
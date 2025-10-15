module mf6_post_module
  ! -- modules
  use utilsmod, only: I2B, I4B, I8B, R4B, R8B, R8ZERO, MXSLEN, &
    errmsg, logmsg, ta, open_file, node_to_icrl
  !
  implicit none
  !
  private
  !
  type, public :: tPostMod
    character(len=MXSLEN)              :: f
    integer(I4B)                       :: kper_beg
    integer(I4B)                       :: kper_end
    integer(I4B), dimension(:), allocatable :: kper_map
    integer(I4B)                       :: nper
    integer(I4B)                       :: nodes
    real(R8B), dimension(:,:), pointer :: xr8_read => null()
  contains
    procedure :: init                   => mf6_post_mod_init
    procedure :: clean                  => mf6_post_mod_clean
    procedure :: set_kper_map           => mf6_post_set_kper_map
    procedure :: read_ulasav            => mf6_post_mod_read_ulasav
    procedure :: read_ulasav_selection  => mf6_post_mod_read_ulasav_selection
    procedure :: read_budget            => mf6_post_mod_read_budget
    generic   :: get_grid               => mf6_get_r4grid
    procedure :: mf6_get_r4grid
  end type tPostMod
  !
  contains
  !
  subroutine mf6_post_mod_init(this, f, kper_beg, kper_end, nodes)
! ******************************************************************************
    ! -- arguments
    class(tPostMod) :: this
    character(len=*), intent(in) :: f
    integer(I4B), intent(in) :: kper_beg
    integer(I4B), intent(in) :: kper_end
    integer(I4B), intent(in) :: nodes
! ------------------------------------------------------------------------------
    call this%clean()
    !
    this%f = f
    this%kper_beg = kper_beg
    this%kper_end = kper_end
    this%nodes    = nodes
    this%nper = kper_end - kper_beg + 1
    if (this%nper <= 0) then
      call errmsg('mf6_post_init: invalid nper.')
    end if
    !
    allocate(this%xr8_read(this%nper, this%nodes))
    this%xr8_read = R8ZERO
  end subroutine mf6_post_mod_init
  
  subroutine mf6_post_mod_clean(this)
! ******************************************************************************
    ! -- arguments
    class(tPostMod) :: this
! ------------------------------------------------------------------------------
    this%f = ''
    this%kper_beg  = 0
    this%kper_end  = 0
    this%nper      = 0
    this%nodes     = 0
    !
    if (allocated(this%kper_map)) then
      deallocate(this%kper_map)
    end if
    if (associated(this%xr8_read)) then
      deallocate(this%xr8_read)
    end if
    this%xr8_read => null()
  end subroutine mf6_post_mod_clean
  !
  subroutine mf6_post_set_kper_map(this)
! ******************************************************************************
    ! -- arguments
    class(tPostMod) :: this
    ! --- local
    integer(I4B) :: i, iper
! ------------------------------------------------------------------------------
    if (allocated(this%kper_map)) deallocate(this%kper_map)
    allocate(this%kper_map(this%kper_end))
    this%kper_map = 0
    iper = 0
    do i = this%kper_beg, this%kper_end
      iper = iper + 1
      this%kper_map(i) = iper
    end do
    !
    return
  end subroutine mf6_post_set_kper_map
  
  subroutine mf6_post_mod_read_ulasav(this)
! ******************************************************************************
    ! -- arguments
    class(tPostMod) :: this
    ! --- local
    ! mf6 header:
    character(len=16) :: text_in ! ulasav
    integer(I4B) :: kstp_in, kper_in, ncol_in, nrow_in, ilay_in
    real(R8B) :: pertim_in, totim_in
    !
    integer(i8b), parameter :: NBHDR = I4B + I4B + R8B + R8B + 16 + I4B + I4B + I4B
    logical :: lop, lread_all
    integer(I4B) :: iu, ios, i, iper, kper_beg, kper_end
    integer(I8B) :: p
! ------------------------------------------------------------------------------
    !
    call this%set_kper_map()
    !
    call open_file(this%f, iu, 'r', .true.)
    p = 1
    do while(.true.) 
      read(unit=iu,iostat=ios, pos=p) kstp_in, kper_in, pertim_in, totim_in, &
        text_in, ncol_in, nrow_in, ilay_in
      p = p + NBHDR
      !
      if (ios /= 0) then
        exit
      end if
      !
      if (ncol_in /= this%nodes) then
        call errmsg('mf6_post_read_ulasav: invalid number of nodes.')
      end if
      !
      if (kper_in <= this%kper_end) then
        iper = this%kper_map(kper_in)
        if (iper > 0) then
          read(unit=iu, iostat=ios, pos=p)(this%xr8_read(iper,i),i=1,this%nodes)
          if (ios /= 0) then
            call errmsg('mf6_post_read_ulasav: could not read data.')
          end if
          this%kper_map(kper_in) = -abs(this%kper_map(kper_in))
        end if
      end if
      p = p + R8B*this%nodes
    end do
    close(iu)
    !
    ! check if all stress periods are read
    if (maxval(this%kper_map) > 0) then
      call logmsg('mf6_post_read_ulasav: data not read for all stress periods.')
    end if
    this%kper_map = abs(this%kper_map)
    !
    return
  end subroutine mf6_post_mod_read_ulasav

  subroutine mf6_post_mod_read_budget(this, ilay, map_ilay,ia, ja, ihc)
! ******************************************************************************
    ! -- arguments
    class(tPostMod) :: this
    integer(I4B), intent(in) :: ilay
    integer(I2B), dimension(:), intent(in) :: map_ilay
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: ja
    integer(I4B), dimension(:), intent(in) :: ihc
    !
    ! --- local
    ! mf6 header:
    character(len=16) :: text_in, txt1id1_in, txt2id1_in, txt1id2_in, txt2id2_in
    integer(I4B) :: kstp_in, kper_in, ndim1_in, ndim2_in, ndim3_in
    integer(I4B) :: imeth_in, ndat_in, nlist_in
    real(R8B) :: delt_in, pertim_in, totim_in
    !
    real(R8B), dimension(:), allocatable :: flowja
    integer(i8b), parameter :: NBHDR1 = I4B + I4B + 16  + I4B + I4B + I4B
    integer(i8b), parameter :: NBHDR2 = I4B + R8B + R8B + R8B
    !
    integer(I4B), dimension(:), allocatable :: nod_lay
    logical :: lop, lread_all
    integer(I4B) :: iu, ios, i, j, iper, kper_beg, kper_end, nja, jlay
    integer(I4B) :: nod_n, nod_m, nodes, nodes_lay, iact, n, m
    integer(I8B) :: p
! ------------------------------------------------------------------------------
    !
    nodes = size(map_ilay)
    nja = size(ja)
    
    do iact = 1, 2
      nodes_lay = 0
      do n = 1, nodes
        if (map_ilay(n) == ilay) then
          nodes_lay = nodes_lay + 1
          if (iact == 2) then
            nod_lay(nodes_lay) = n
          end if
        end if
      end do
      if ((iact == 1).and.(nodes_lay > 0))  then
        allocate(nod_lay(nodes_lay))
      end if
    end do
    !
    nja = size(ja)
    allocate(flowja(nja))
    !
    call this%set_kper_map()
    !
    call open_file(this%f, iu, 'r', .true.)
    p = 1
    do while(.true.) 
      ! record 1:
      read(unit=iu,iostat=ios, pos=p) kstp_in, kper_in, &
        text_in, ndim1_in, ndim2_in, ndim3_in
      
      if (ios /= 0) then
        exit
      end if
      
      p = p + NBHDR1
      ! record 2:
      read(unit=iu,iostat=ios, pos=p) imeth_in, delt_in, pertim_in, totim_in
      p = p + NBHDR2
      !
      if (adjustl(text_in) /= 'FLOW-JA-FACE') then
        select case(imeth_in)
        case(1)
          p = p + R8B*nja
        case(6)
          read(unit=iu,iostat=ios, pos=p) txt1id1_in; p = p + 16
          read(unit=iu,iostat=ios, pos=p) txt2id1_in; p = p + 16
          read(unit=iu,iostat=ios, pos=p) txt1id2_in; p = p + 16
          read(unit=iu,iostat=ios, pos=p) txt2id2_in; p = p + 16
          read(unit=iu,iostat=ios, pos=p) ndat_in; p = p + I4B
          p = p + 16 * (ndat_in - 1) ! skip: (AUXTXT(N),N=1,NDAT-1)
          read(unit=iu,iostat=ios, pos=p) nlist_in; p = p + I4B
          p = p + (I4B + I4B + R8B * ndat_in) * nlist_in ! skip: ((ID1(N),ID2(N),(DATA2D(I,N),I=1,NDAT)),N=1,NLIST)
        case default
          call errmsg('mf6_post_mod_read_budget: non supported IMETH')
        end select
        cycle
      end if
      
      ! additional checks
      if (imeth_in /= 1) then
        call errmsg('mf6_post_mod_read_budget: non supported IMETH.')
      end if
      if (ndim1_in /= nja) then
        call errmsg('mf6_post_mod_read_budget: inconsistent NDIM1.')
      end if
      !
      if (ios /= 0) then
        exit
      end if
      !
      if (kper_in <= this%kper_end) then
        iper = this%kper_map(kper_in)
        if (iper > 0) then
          read(unit=iu, iostat=ios, pos=p) flowja
          if (ios /= 0) then
            call errmsg('mf6_post_mod_read_budget: could not read data.')
          end if
          this%kper_map(kper_in) = -abs(this%kper_map(kper_in))
          !
          ! fill the budgets
          do i = 1, nodes_lay
            n = nod_lay(i)
            do j = ia(n)+1, ia(n+1)-1
              m = ja(j)
              if (ihc(j) == 0) then ! vertial connection
                jlay = map_ilay(m)
                if (jlay == ilay) then
                  call errmsg('mf6_post_mod_read_budget: program error layer.')
                end if
                if (jlay > ilay) then
                  this%xr8_read(iper,n) = this%xr8_read(iper,n) + flowja(j)
                end if
              end if
            end do
          end do
        end if
      end if
      
      p = p + R8B*nja
    end do
    close(iu)
    !
    ! check if all stress periods are read
    if (maxval(this%kper_map) > 0) then
      call errmsg('mf6_post_mod_read_budget: data not read for all stress periods.')
    end if
    this%kper_map = abs(this%kper_map)
    !
    ! clean up
    if (allocated(nod_lay)) deallocate(nod_lay)
    if (allocated(flowja)) deallocate(flowja)
    !
    return
  end subroutine mf6_post_mod_read_budget
  
  subroutine mf6_post_mod_read_ulasav_selection(this, nodes_read, heads_read)
! ******************************************************************************
    ! -- arguments
    class(tPostMod) :: this
    integer(I4B), dimension(:), intent(in) :: nodes_read
    real(R8B), dimension(:,:), intent(out) :: heads_read
    ! --- local
    ! mf6 header:
    character(len=16) :: text_in ! ulasav
    integer(I4B) :: kstp_in, kper_in, ncol_in, nrow_in, ilay_in
    real(R8B) :: pertim_in, totim_in
    !
    integer(I8B), parameter :: NBHDR = I4B + I4B + R8B + R8B + 16 + I4B + I4B + I4B
    integer(I4B) :: iu, ios, i, iper, kper_beg, kper_end, nod
    integer(I8B) :: p, p_nod
    !
    real(R8B) :: r8v
! ------------------------------------------------------------------------------
    !
    call this%set_kper_map()
    !
    call open_file(this%f, iu, 'r', .true.)
    !
    p = 1
    do while(.true.) 
      read(unit=iu,iostat=ios, pos=p) kstp_in, kper_in, pertim_in, totim_in, &
        text_in, ncol_in, nrow_in, ilay_in
      !
      if (ios /= 0) then
        exit
      end if
      if (ncol_in /= this%nodes) then
        call errmsg('mf6_post_read_ulasav_selection: invalid number of nodes.')
      end if
      p = p + NBHDR
      !
      if (kper_in <= this%kper_end) then
        iper = this%kper_map(kper_in)
        if (iper > 0) then
          do i = 1, size(nodes_read)
            nod = nodes_read(i)
            p_nod = p + (nod-1)*R8B
            read(unit=iu, iostat=ios, pos=p_nod) r8v
            heads_read(iper,i) = r8v
            if (ios /= 0) then
              call errmsg('mf6_post_read_ulasav_selection: could nod read head.')
            end if
          end do
          this%kper_map(iper) = -abs(this%kper_map(iper))
        end if
      end if
      p = p + R8B*this%nodes
    end do
    !
    close(iu)
    !
    ! check if all stress periods are read
    if (maxval(this%kper_map) > 0) then
      call logmsg('WARNING mf6_post_mod_read_ulasav_selection: '// &
        'data not read for all stress periods.')
    end if
    this%kper_map = abs(this%kper_map)
    
    return
  end subroutine mf6_post_mod_read_ulasav_selection
  
  subroutine mf6_get_r4grid(this, kper, nod_map, mvr4, xr4, nodes_offset, &
    node_area)
! ******************************************************************************
    ! -- arguments
    class(tPostMod) :: this
    integer(I4B), intent(in) :: kper
    integer(I4B), dimension(:,:), intent(in) :: nod_map
    real(R4B), intent(in) :: mvr4
    real(R4B), dimension(:,:), allocatable, intent(inout) :: xr4
    integer(I4B), intent(in), optional :: nodes_offset
    real(R8B), dimension(:), intent(in), optional :: node_area
    ! --- local
    integer(I4B) :: iper, nc, nr, ic, ir, nod, offset
    real(R8B) :: r8v, area
    logical :: larea
! ------------------------------------------------------------------------------
    if (present(nodes_offset)) then
      offset = nodes_offset
    else
      offset = 0
    end if
    !
    if (present(node_area)) then
      larea = .true.
    else
      larea = .false.
    end if
    !
    iper = kper - this%kper_beg + 1
    !
    ! check
    if ((iper <= 0).or.(iper > this%nper)) then
      call errmsg('mf6_get_r4grid: invalid kper.')
    end if
    !
    nc = size(nod_map,1); nr = size(nod_map,2)
    if (allocated(xr4)) deallocate(xr4)
    allocate(xr4(nc,nr)); xr4 = mvr4
    do ir = 1, nr; do ic = 1, nc
      nod = nod_map(ic,ir)
      if (nod > 0) then
        nod = nod + offset
        if (nod <= this%nodes) then
          if (larea) then
            area = node_area(nod)
            if (area == R8ZERO) then
              call errmsg('mf6_get_r4grid: program error, zero area.')
            end if
            r8v = this%xr8_read(iper,nod)
            r8v = r8v / area
          else
            r8v = this%xr8_read(iper,nod)
          end if
          xr4(ic,ir) = real(r8v, R4B)
        end if
      end if
    end do; end do
    !
    return
  end subroutine mf6_get_r4grid
  
end module mf6_post_module
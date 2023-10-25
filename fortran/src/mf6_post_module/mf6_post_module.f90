module mf6_post_module
  ! -- modules
  use utilsmod, only: I4B, I8B, R4B, R8B, R8ZERO, MXSLEN, &
    errmsg, logmsg, ta, open_file
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
      call errmsg('mf6_post_read_ulasav: data not read for all stress periods.')
    end if
    this%kper_map = abs(this%kper_map)
    !
    return
  end subroutine mf6_post_mod_read_ulasav
  
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
    integer(i8b), parameter :: NBHDR = I4B + I4B + R8B + R8B + 16 + I4B + I4B + I4B
    integer(I4B) :: iu, ios, i, iper, kper_beg, kper_end, nod
    integer(I8B) :: p, p_nod
    
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
            read(unit=iu,iostat=ios, pos=p_nod) heads_read(iper,i)
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
      call errmsg('mf6_post_mod_read_ulasav_selection: '// &
        'data not read for all stress periods.')
    end if
    this%kper_map = abs(this%kper_map)
    
    return
  end subroutine mf6_post_mod_read_ulasav_selection
  
  subroutine mf6_get_r4grid(this, kper, nod_map, mvr4, xr4, nodes_offset)
! ******************************************************************************
    ! -- arguments
    class(tPostMod) :: this
    integer(I4B), intent(in) :: kper
    integer(I4B), dimension(:,:), intent(in) :: nod_map
    real(R4B), intent(in) :: mvr4
    real(R4B), dimension(:,:), allocatable, intent(inout) :: xr4
    integer(I4B), intent(in), optional :: nodes_offset
    ! --- local
    integer(I4B) :: iper, nc, nr, ic, ir, nod, offset
! ------------------------------------------------------------------------------
    if (present(nodes_offset)) then
      offset = nodes_offset
    else
      offset = 0
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
          xr4(ic,ir) = this%xr8_read(iper,nod)
        end if
      end if
    end do; end do
    !
    return
  end subroutine mf6_get_r4grid
  
end module mf6_post_module
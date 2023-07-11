module mf6_wbd_mod
  ! modules
  use utilsmod, only: I1B, I2B, I4B, I8B, R4B, R8B, &
    I_I1, I_I2, I_I4, I_I8, I_R4, I_R8, I_C, &
    MXSLEN, I4ZERO, I4ONE, I8ZERO, R8ZERO, num_names, &
    errmsg, logmsg, ta, open_file, tCSV
  !
  implicit none
  !
  private
  !
  integer(I4B), parameter, public :: i_asc    = 1
  integer(I4B), parameter, public :: i_bin    = 2
  integer(I4B), parameter, public :: i_binpos = 3
  integer(I4B), parameter :: max_nr_csv = 1000
  integer(I8B), parameter :: mv = I8ZERO

  !
  ! work arrays
  integer(I1B), dimension(:), allocatable   :: i1a_wk
  integer(I2B), dimension(:), allocatable   :: i2a_wk
  integer(I4B), dimension(:), allocatable   :: i4a_wk
  integer(I8B), dimension(:), allocatable   :: i8a_wk
  real(R4B),    dimension(:), allocatable   :: r4a_wk
  real(R8B),    dimension(:), allocatable   :: r8a_wk
  !
  integer(I1B), dimension(:,:), allocatable :: i1x_wk
  integer(I2B), dimension(:,:), allocatable :: i2x_wk
  integer(I4B), dimension(:,:), allocatable :: i4x_wk
  integer(I8B), dimension(:,:), allocatable :: i8x_wk
  real(R4B),    dimension(:,:), allocatable :: r4x_wk
  real(R8B),    dimension(:,:), allocatable :: r8x_wk
  !
  type, public :: tMf6Wbd
    type(tCsv), pointer :: csv => null()
    !
    logical :: use_binpos
    character(len=MXSLEN) :: f_binpos
    integer(I4B) :: iu_binpos
  contains
    procedure :: init        => tMf6Wbd_init
    procedure :: clean       => tMf6Wbd_clean
    procedure :: write_array => tMf6Wbd_write_array
    procedure :: read_array  => tMf6Wbd_read_array
    procedure :: write_list  => tMf6Wbd_write_list
    procedure :: write_csv   => tMf6Wbd_write_csv
    procedure :: read_csv    => tMf6Wbd_read_csv
  end type tMf6Wbd
  
contains
  
! ==============================================================================
! ==============================================================================
! tMf6Wbd
! ==============================================================================
! ==============================================================================
  
  subroutine tMf6Wbd_init(this, f_csv, f_binpos)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMf6Wbd) :: this
    character(len=*), intent(in) :: f_csv
    character(len=*), intent(in), optional :: f_binpos
    ! -- local
    logical :: lex
    integer(I4B) :: iu
    type(tCSV), pointer :: csv => null()
! ------------------------------------------------------------------------------
    call this%clean()
    !
    allocate(this%csv); csv => this%csv
    !
    call csv%init(file=f_csv, &
      hdr_keys=['id', 'nr', 'file', 'file_type', 'data_type', &
       'binpos_beg', 'binpos_end'],&
      nr=max_nr_csv, &
        hdr_i_type=[I_C, I_I8, I_C, I_C, I_C, I_I8, I_I8])
    csv%nr = 0
    !
    if (present(f_binpos)) then
      this%use_binpos = .true.
      this%f_binpos = f_binpos
      !
      ! remove the old file when existing
      inquire(file=this%f_binpos, exist=lex)
      if (lex) then
        call open_file(this%f_binpos, this%iu_binpos, 'r', .true.)
        close(this%iu_binpos, status='delete'); this%iu_binpos = -1
      end if
    else
      this%use_binpos = .false.
    end if
    !
    return
  end subroutine tMf6Wbd_init
  
  subroutine tMf6Wbd_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMf6Wbd) :: this
! ------------------------------------------------------------------------------
    if (associated(this%csv)) then
      call this%csv%clean()
      deallocate(this%csv)
    end if
    
    this%csv => null()
    this%use_binpos = .false.
    this%f_binpos = ''
    this%iu_binpos = -1
    !
    return
  end subroutine tMf6Wbd_clean
  !
  subroutine tMf6Wbd_write_array(this, id, nod_ptr, &
    i1a, i2a, i4a, i8a, r4a, r8a, &
    f_asc, f_bin)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMf6Wbd) :: this
    !
    character(len=*), intent(in) :: id
    integer(I1B), dimension(:), intent(in), optional :: nod_ptr
    !
    integer(I1B), dimension(:), intent(in), optional :: i1a
    integer(I2B), dimension(:), intent(in), optional :: i2a
    integer(I4B), dimension(:), intent(in), optional :: i4a
    integer(I8B), dimension(:), intent(in), optional :: i8a
    real(R4B),    dimension(:), intent(in), optional :: r4a
    real(R8B),    dimension(:), intent(in), optional :: r8a
    !
    character(len=*), intent(in), optional :: f_asc
    character(len=*), intent(in), optional :: f_bin
    ! -- local
    ! dummy header:
    integer(I4B), parameter :: I4DUM = I4ONE
    real(R8B), parameter :: R8DUM = R8ZERO
    character(len=16) :: C16DUM = ''
    !
    type(tCsv), pointer :: csv => null()
    character(len=MXSLEN) :: f
    logical :: lex, lop
    integer(I4B) :: iu, n, i, i_dat, i_data_type, ir_csv
    integer(I8B) :: binpos_beg, binpos_end
    !
! ------------------------------------------------------------------------------
    !
    ! checks
    n = 0
    if (present(i1a)) n = n + 1
    if (present(i2a)) n = n + 1
    if (present(i4a)) n = n + 1
    if (present(i8a)) n = n + 1
    if (present(r4a)) n = n + 1
    if (present(r8a)) n = n + 1
    if (n /= 1) call errmsg('tMf6Wbd_write_array: error arguments.')
    n = 0
    if (present(f_asc)) n = n + 1
    if (present(f_bin)) n = n + 1
    if (n > 1) call errmsg('tMf6Wbd_write_array: error arguments.')
    !
    i_data_type = -1
    if (present(i1a)) i_data_type = I_I1
    if (present(i2a)) i_data_type = I_I2
    if (present(i4a)) i_data_type = I_I4
    if (present(i8a)) i_data_type = I_I8
    if (present(r4a)) i_data_type = I_R4
    if (present(r8a)) i_data_type = I_R8
    !
    ! count and set the output array
    n = 0
    if (present(nod_ptr)) then
      do i = 1, size(nod_ptr)
        if (nod_ptr(i) == 1) n = n + 1
      end do
      if (n == 0) then
        call errmsg('tMf6Wbd_write_array: no active cells found in nod_ptr.')
      end if
      if (present(i1a)) then
        allocate(i1a_wk(n)); n = 0
        do i = 1, size(nod_ptr)
          if (nod_ptr(i) == 1) then
            n = n + 1; i1a_wk(n) = i1a(i)
          end if
        end do
      end if
      if (present(i2a)) then
        allocate(i2a_wk(n)); n = 0
        do i = 1, size(nod_ptr)
          if (nod_ptr(i) == 1) then
            n = n + 1; i2a_wk(n) = i2a(i)
          end if
        end do
      end if
      if (present(i4a)) then
        allocate(i4a_wk(n)); n = 0
        do i = 1, size(nod_ptr)
          if (nod_ptr(i) == 1) then
            n = n + 1; i4a_wk(n) = i4a(i)
          end if
        end do
      end if
      if (present(i8a)) then
        allocate(i8a_wk(n)); n = 0
        do i = 1, size(nod_ptr)
          if (nod_ptr(i) == 1) then
            n = n + 1; i8a_wk(n) = i8a(i)
          end if
        end do
      end if
      if (present(r4a)) then
        allocate(r4a_wk(n)); n = 0
        do i = 1, size(nod_ptr)
          if (nod_ptr(i) == 1) then
            n = n + 1; r4a_wk(n) = r4a(i)
          end if
        end do
      end if
      if (present(r8a)) then
        allocate(r8a_wk(n)); n = 0
        do i = 1, size(nod_ptr)
          if (nod_ptr(i) == 1) then
            n = n + 1; r8a_wk(n) = r8a(i)
          end if
        end do
      end if
    else
      if (present(i1a)) then
        n = size(i1a); allocate(i1a_wk, source=i1a)
      end if
      if (present(i2a)) then
        n = size(i2a); allocate(i2a_wk, source=i2a)
      end if
      if (present(i4a)) then
        n = size(i4a); allocate(i4a_wk, source=i4a)
      end if
      if (present(i8a)) then
        n = size(i8a); allocate(i8a_wk, source=i8a)
      end if
      if (present(r4a)) then
        n = size(r4a); allocate(r4a_wk, source=r4a)
      end if
      if (present(r8a)) then
        n = size(r8a); allocate(r8a_wk, source=r8a)
      end if
    end if
    !
    ! determine the output type
    i_dat = i_binpos
    if (present(f_asc)) i_dat = i_asc
    if (present(f_bin)) i_dat = i_bin
    !
    ! save the metadata
    csv => this%csv
    ir_csv = csv%get_row(key_val=trim(id), hdr_key=csv%hdr(1)%key)
    if (ir_csv > 0) then
      call csv%clean_and_init_row(ir_csv)
    else
      csv%nr = csv%nr + 1
      ir_csv = csv%nr
    end if
    !
    call csv%set_val(ic=1, ir=ir_csv, cv=trim(id)) ! id
    call csv%set_val(ic=2, ir=ir_csv, i8v=int(n,I8B)) ! nr
    call csv%set_val(ic=5, ir=ir_csv, cv=trim(num_names(i_data_type))) ! data_type
    !
    select case(i_dat)
    case(i_binpos)
      if (len_trim(this%f_binpos) == 0) then
        call errmsg('tMf6Wbd_write_array: not binpos file specified.')
      end if
      inquire(file=this%f_binpos, exist=lex)
      if (lex) then
        call open_file(this%f_binpos, this%iu_binpos, 'w', .true., 'append')
      else
        call open_file(this%f_binpos, this%iu_binpos, 'w', .true.)
      end if
      !
      ! write the dummy header
      inquire(this%iu_binpos, pos=binpos_beg)
      write(this%iu_binpos) I4DUM, I4DUM, R8DUM, R8DUM, C16DUM, n, I4DUM, I4DUM
      if (present(i1a)) write(this%iu_binpos) i1a_wk
      if (present(i2a)) write(this%iu_binpos) i2a_wk
      if (present(i4a)) write(this%iu_binpos) i4a_wk
      if (present(i8a)) write(this%iu_binpos) i8a_wk
      if (present(r4a)) write(this%iu_binpos) r4a_wk
      if (present(r8a)) write(this%iu_binpos) r8a_wk
      close(this%iu_binpos)
      call open_file(this%f_binpos, this%iu_binpos, 'w', .true., 'append')
      inquire(this%iu_binpos, pos=binpos_end)
      close(this%iu_binpos)
      !
      ! save the metadata
      call csv%set_val(ic= 3, ir=ir_csv, cv=trim(this%f_binpos)) ! file
      call csv%set_val(ic= 4, ir=ir_csv, cv='binpos') ! file_type
      call csv%set_val(ic= 6, ir=ir_csv, i8v=binpos_beg) ! binpos_beg
      call csv%set_val(ic= 7, ir=ir_csv, i8v=binpos_end) ! binpos_end
    case(i_bin)
      f =  f_bin; call open_file(f, iu, 'w', .true.)
      write(iu) I4DUM, I4DUM, R8DUM, R8DUM, C16DUM, n, I4DUM, I4DUM
      if (present(i1a)) write(iu) i1a_wk
      if (present(i2a)) write(iu) i2a_wk
      if (present(i4a)) write(iu) i4a_wk
      if (present(i8a)) write(iu) i8a_wk
      if (present(r4a)) write(iu) r4a_wk
      if (present(r8a)) write(iu) r8a_wk
      close(iu)
      !
      ! save the metadata
      call csv%set_val(ic= 3, ir=ir_csv, cv=trim(f_bin)) ! file
      call csv%set_val(ic= 4, ir=ir_csv, cv='bin') ! file_type
      call csv%set_val(ic= 6, ir=ir_csv, i8v=mv) ! binpos_beg
      call csv%set_val(ic= 7, ir=ir_csv, i8v=mv) ! binpos_end
    case(i_asc)
      f = f_asc; call open_file(f, iu, 'w')
      if (present(i1a)) then
        do i = 1, size(i1a_wk)
          write(iu,'(a)') ta([i1a_wk(i)])
        end do
      end if
      if (present(i2a)) then
        do i = 1, size(i2a_wk)
          write(iu,'(a)') ta([i2a_wk(i)])
        end do
      end if
      if (present(i4a)) then
        do i = 1, size(i4a_wk)
          write(iu,'(a)') ta([i4a_wk(i)])
        end do
      end if
      if (present(i8a)) then
        do i = 1, size(i8a_wk)
          write(iu,'(a)') ta([i8a_wk(i)])
        end do
      end if
      if (present(r4a)) then
        do i = 1, size(r4a_wk)
          write(iu,'(a)') ta([r4a_wk(i)])
        end do
      end if
      if (present(r8a)) then
        do i = 1, size(r8a_wk)
          write(iu,'(a)') ta([r8a_wk(i)])
        end do
      end if
      close(iu)
      !
      call csv%set_val(ic= 3, ir=ir_csv, cv=trim(f_asc)) ! file
      call csv%set_val(ic= 4, ir=ir_csv, cv='asc') ! file_type
      call csv%set_val(ic= 6, ir=ir_csv, i8v=mv) ! binpos_beg
      call csv%set_val(ic= 7, ir=ir_csv, i8v=mv) ! binpos_end
    end select
    !
    ! cleanup
    if (allocated(i1a_wk)) deallocate(i1a_wk)
    if (allocated(i2a_wk)) deallocate(i2a_wk)
    if (allocated(i4a_wk)) deallocate(i4a_wk)
    if (allocated(i8a_wk)) deallocate(i8a_wk)
    if (allocated(r4a_wk)) deallocate(r4a_wk)
    if (allocated(r8a_wk)) deallocate(r8a_wk)
    !
    return
  end subroutine tMf6Wbd_write_array
  !
  subroutine tMf6Wbd_read_array(this, id, i1a, i2a, i4a, r8a)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMf6Wbd) :: this
    !
    character(len=*), intent(in) :: id
    integer(I1B), dimension(:), intent(inout), pointer, optional :: i1a
    integer(I2B), dimension(:), intent(inout), pointer, optional :: i2a
    integer(I4B), dimension(:), intent(inout), pointer, optional :: i4a
    real(R8B),    dimension(:), intent(inout), pointer, optional :: r8a
    !
    ! -- local
    type(tCsv), pointer :: csv => null()
    character(len=MXSLEN), dimension(:), allocatable :: ids
    character(len=MXSLEN) :: file, file_type, data_type
    integer(I4B) :: ir, ic, i, i_dat, n, n_read, iu
    integer(I8B) :: binpos_beg
    !
    ! header
    integer(I4B) :: I4DUM
    real(R8B) :: R8DUM
    character(len=16) :: C16DUM
! ------------------------------------------------------------------------------
    !
    csv => this%csv
    
    ! first, find the row corresponding to the id
    call csv%get_column(key='id', ca=ids)
    ir = 0
    do i = 1, size(ids)
      if (trim(ids(i)) == trim(id)) then
        ir = i; exit
      end if
    end do
    if (ir == 0) then
      call errmsg('tMf6Wbd_read_array: id '//trim(id)//' not found.')
    end if
    !
    ! determine the data type
    call csv%get_val(ir=ir, ic=csv%get_col('nr'), i4v=n)
    if (n <= 0) then
      call errmsg('tMf6Wbd_read_array: array dimension <= 0.')
    end if
    !
    call csv%get_val(ir=ir, ic=csv%get_col('data_type'), cv=data_type)
    !
    if (present(i1a)) then
      if (data_type /= 'integer_1b') then
        call errmsg('tMf6Wbd_read_array: inconsistent argument.')
      end if
      if (associated(i1a)) deallocate(i1a)
      allocate(i1a(n))
    end if
    if (present(i2a)) then
      if (data_type /= 'integer_2b') then
        call errmsg('tMf6Wbd_read_array: inconsistent argument.')
      end if
      if (associated(i2a)) deallocate(i2a)
      allocate(i2a(n))
    end if
    if (present(i4a)) then
      if (data_type /= 'integer_4b') then
        call errmsg('tMf6Wbd_read_array: inconsistent argument.')
      end if
      if (associated(i4a)) deallocate(i4a)
      allocate(i4a(n))
    end if
    if (present(r8a)) then
      if (data_type /= 'real_8b') then
        call errmsg('tMf6Wbd_read_array: inconsistent argument.')
      end if
      if (associated(r8a)) deallocate(r8a)
      allocate(r8a(n))
    end if
    !
    call csv%get_val(ir=ir, ic=csv%get_col('file'), cv=file)
    call csv%get_val(ir=ir, ic=csv%get_col('file_type'), cv=file_type)
    !
    select case(file_type)
    case('binpos')
      call csv%get_val(ir=ir, ic=csv%get_col('binpos_beg'), i8v=binpos_beg)
      call open_file(file, iu, 'r', .true.)
      write(unit=iu, pos=binpos_beg)
      inquire(unit=iu, pos=binpos_beg)
      read(iu) I4DUM, I4DUM, R8DUM, R8DUM, C16DUM, n_read, I4DUM, I4DUM
      if (n /= n_read) then
        call errmsg('tMf6Wbd_read_array: inconsistent dimensions.')
      end if
      if (present(i1a)) read(iu) i1a
      if (present(i2a)) read(iu) i2a
      if (present(i4a)) read(iu) i4a
      if (present(r8a)) read(iu) r8a
      close(iu)
    case('asc')
      call open_file(file, iu, 'r')
      if (present(i1a)) read(iu,*) i1a
      if (present(i2a)) read(iu,*) i2a
      if (present(i4a)) read(iu,*) i4a
      if (present(r8a)) read(iu,*) r8a
      close(iu)
    case default
      call errmsg('tMf6Wbd_read_array: not yet supported.')
    end select
    !
    return
  end subroutine tMf6Wbd_read_array
  
  subroutine tMf6Wbd_write_list(this, id, nod_ptr, &
    i4a1, i4a2, i4a3, &
    r8a1, r8a2, r8a3, r8x, f_asc, f_bin)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMf6Wbd) :: this
    !
    character(len=*), intent(in) :: id
    integer(I1B), dimension(:), intent(in), optional :: nod_ptr
    integer(I4B), dimension(:), intent(in), optional :: i4a1
    integer(I4B), dimension(:), intent(in), optional :: i4a2
    integer(I4B), dimension(:), intent(in), optional :: i4a3
    real(R8B),    dimension(:), intent(in), optional :: r8a1
    real(R8B),    dimension(:), intent(in), optional :: r8a2
    real(R8B),    dimension(:), intent(in), optional :: r8a3
    real(R8B),    dimension(:,:), intent(in), optional :: r8x
    character(len=*), intent(in), optional :: f_asc
    character(len=*), intent(in), optional :: f_bin
    ! -- local
    type(tCsv), pointer :: csv => null()
    character(len=MXSLEN) :: f, fmt
    logical :: lex, lop
    integer(I4B) :: ni4a, nr8a, ncr8x
    integer(I4B) :: siz, n, ir, ic, i_dat, iu, i, j, ir_csv
    integer(I4B), dimension(3) :: i4siz, r8siz
    integer(I8B) :: binpos_beg, binpos_end
! ------------------------------------------------------------------------------
    !
    ! checks
    ni4a = 0; nr8a = 0
    if (present(i4a1)) then
      ni4a = ni4a + 1; i4siz(ni4a) = size(i4a1)
    end if
    if (present(i4a2)) then
      ni4a = ni4a + 1; i4siz(ni4a) = size(i4a2)
    end if
    if (present(i4a3)) then
      ni4a = ni4a + 1; i4siz(ni4a) = size(i4a3)
    end if
    if (present(r8a1)) then
      nr8a = nr8a + 1; r8siz(nr8a) = size(r8a1)
    end if
    if (present(r8a2)) then
      nr8a = nr8a + 1; r8siz(nr8a) = size(r8a2)
    end if
    if (present(r8a3)) then
      nr8a = nr8a + 1; r8siz(nr8a) = size(r8a3)
    end if
    if (present(r8x)) then
      if (present(r8a1).or.(present(r8a2).or.present(r8a3))) then
        call errmsg('tMf6Wbd_write_list: invalid arguments.')
      end if
      ncr8x = size(r8x,1)
      do i = 1, ncr8x
        nr8a = nr8a + 1; r8siz(nr8a) = size(r8x,2)
      end do
    end if
    !
    ! check the sizes
    !TODO
    !
    if (ni4a == 3) then ! exchanges
      if (nr8a /= 3) then
        call logmsg('tMf6Wbd_write_list: invalid call.')
      end if
      if (present(nod_ptr)) then
        call logmsg('tMf6Wbd_write_list: invalid call.')
      end if
    end if
    if (ni4a == 2) then ! HFB
      if (nr8a /= 1) then
        call logmsg('tMf6Wbd_write_list: invalid call.')
      end if
    end if
    !
    if (present(nod_ptr)) then
      n = 0
      do ir = 1, size(nod_ptr)
        if (nod_ptr(ir) == 1) n = n + 1
      end do
      if (n == 0) then
        call errmsg('tMf6Wbd_write_list: no active cells found in nod_ptr.')
      end if
      allocate(i4x_wk(ni4a,n), r8x_wk(nr8a,n))
      ic = 0
      if (present(i4a1)) then
        ic = ic + 1; n = 0
        do ir = 1, size(nod_ptr)
          if (nod_ptr(ir) == 1) then
            n = n + 1; i4x_wk(ic,n) = i4a1(ir)
          end if
        end do
      end if
      if (present(i4a2)) then
        ic = ic + 1; n = 0
        do ir = 1, size(nod_ptr)
          if (nod_ptr(ir) == 1) then
            n = n + 1; i4x_wk(ic,n) = i4a2(ir)
          end if
        end do
      end if
      if (present(i4a3)) then
        ic = ic + 1; n = 0
        do ir = 1, size(nod_ptr)
          if (nod_ptr(ir) == 1) then
            n = n + 1; i4x_wk(ic,n) = i4a3(ir)
          end if
        end do
      end if
      ic = 0
      if (present(r8a1)) then
        ic = ic + 1; n = 0
        do ir = 1, size(nod_ptr)
          if (nod_ptr(ir) == 1) then
            n = n + 1; r8x_wk(ic,n) = r8a1(ir)
          end if
        end do
      end if
      if (present(r8a2)) then
        ic = ic + 1; n = 0
        do ir = 1, size(nod_ptr)
          if (nod_ptr(ir) == 1) then
            n = n + 1; r8x_wk(ic,n) = r8a2(ir)
          end if
        end do
      end if
      if (present(r8a3)) then
        ic = ic + 1; n = 0
        do ir = 1, size(nod_ptr)
          if (nod_ptr(ir) == 1) then
            n = n + 1; r8x_wk(ic,n) = r8a3(ir)
          end if
        end do
      end if
      if (present(r8x)) then
        do ic = 1, ncr8x
          do ir = 1, size(nod_ptr)
            if (nod_ptr(ir) == 1) then
              n = n + 1; r8x_wk(ic,n) = r8x(ic,ir)
            end if
          end do
        end do
      end if
    else
      n = size(i4a1)
      allocate(i4x_wk(ni4a,n), r8x_wk(nr8a,n))
      ic = 0
      if (present(i4a1)) then
        ic = ic + 1
        i4x_wk(ic,:) = i4a1
      end if
      if (present(i4a2)) then
        ic = ic + 1
        i4x_wk(ic,:) = i4a2
      end if
      if (present(i4a3)) then
        ic = ic + 1
        i4x_wk(ic,:) = i4a3
      end if
      ic = 0
      if (present(r8a1)) then
        ic = ic + 1
        r8x_wk(ic,:) = r8a1
      end if
      if (present(r8a2)) then
        ic = ic + 1
        r8x_wk(ic,:) = r8a2
      end if
      if (present(r8a3)) then
        ic = ic + 1
        r8x_wk(ic,:) = r8a3
      end if
      if (present(r8x)) then
        do ic = 1, ncr8x
          r8x_wk(ic,:) = r8x(ic,:)
        end do
      end if
    end if
    !
    ! determine the output type
    i_dat = i_binpos
    if (present(f_asc)) i_dat = i_asc
    if (present(f_bin)) i_dat = i_bin
    !
    ! save the metadata
    csv => this%csv
    ir_csv = csv%get_row(key_val=trim(id), hdr_key=csv%hdr(1)%key)
    if (ir_csv > 0) then
      call csv%clean_and_init_row(ir_csv)
    else
      csv%nr = csv%nr + 1
      ir_csv = csv%nr
    end if
    !
    call csv%set_val(ic= 1, ir=ir_csv, cv=trim(id)) ! id
    call csv%set_val(ic= 2, ir=ir_csv, i8v=int(n,I8B)) ! nr
    !
    select case(i_dat)
    case(i_binpos)
      if (len_trim(this%f_binpos) == 0) then
        call errmsg('tMf6Wbd_write_array: not binpos file specified.')
      end if
      inquire(file=this%f_binpos, exist=lex)
      if (lex) then
        call open_file(this%f_binpos, this%iu_binpos, 'w', .true., 'append')
      else
        call open_file(this%f_binpos, this%iu_binpos, 'w', .true.)
      end if
      !
      inquire(this%iu_binpos, pos=binpos_beg)
      write(this%iu_binpos) (((i4x_wk(j,i), j = 1, ni4a), &
        (r8x_wk(j,i), j = 1, nr8a)), i = 1, n)
      close(this%iu_binpos)
      call open_file(this%f_binpos, this%iu_binpos, 'w', .true., 'append')
      inquire(this%iu_binpos, pos=binpos_end)
      close(this%iu_binpos)
      !
      ! save the metadata
      call csv%set_val(ic= 3, ir=ir_csv, cv=trim(this%f_binpos)) ! file
      call csv%set_val(ic= 4, ir=ir_csv, cv='binpos') ! file_type
      call csv%set_val(ic= 5, ir=ir_csv, cv='mixed') ! data_type
      call csv%set_val(ic= 6, ir=ir_csv, i8v=binpos_beg) ! binpos_beg
      call csv%set_val(ic= 7, ir=ir_csv, i8v=binpos_end) ! binpos_end
    
    case(i_bin)
      f =  f_bin; call open_file(f, iu, 'w', .true.)
      !call open_file(f, iu, 'w')
      !fmt = '('//trim(ta([ni4a]))//'i4,'//trim(ta([nr8a]))//'f4.0)'
      !write(iu,fmt) (((i4x_wk(j,i), j = 1, ni4a), &
      ! (r8x_wk(j,i), j = 1, nr8a)), i = 1, n)
      write(iu) (((i4x_wk(j,i), j = 1, ni4a), &
        (r8x_wk(j,i), j = 1, nr8a)), i = 1, n)
      close(iu)
      !
      ! save the metadata
      call csv%set_val(ic= 3, ir=ir_csv, cv=trim(f_bin)) ! file
      call csv%set_val(ic= 4, ir=ir_csv, cv='bin') ! file_type
      call csv%set_val(ic= 5, ir=ir_csv, cv='mixed') ! data_type
      call csv%set_val(ic= 6, ir=ir_csv, i8v=mv) ! binpos_beg
      call csv%set_val(ic= 7, ir=ir_csv, i8v=mv) ! binpos_end
    case(i_asc)
      f = f_asc; call open_file(f, iu, 'w')
      do i = 1, n
        write(iu,'(2a)') ta([i4x_wk(:,i)], sep_in=' ')//' '// &
                         ta([r8x_wk(:,i)], sep_in=' ')
      end do
      close(iu)
      !
      call csv%set_val(ic= 3, ir=ir_csv, cv=trim(f_asc)) ! file
      call csv%set_val(ic= 4, ir=ir_csv, cv='asc') ! file_type
      call csv%set_val(ic= 5, ir=ir_csv, cv='mixed') ! data_type
      call csv%set_val(ic= 6, ir=ir_csv, i8v=mv) ! binpos_beg
      call csv%set_val(ic= 7, ir=ir_csv, i8v=mv) ! binpos_end
    end select
    !
    ! clean up
    if (allocated(i4x_wk)) deallocate(i4x_wk)
    if (allocated(r8x_wk)) deallocate(r8x_wk)
    !
    return
  end subroutine tMf6Wbd_write_list
  
  subroutine tMf6Wbd_write_csv(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMf6Wbd) :: this
! ------------------------------------------------------------------------------
    !
    call this%csv%write(i8mv=mv)
    !call this%csv%write()
    !
    return
  end subroutine tMf6Wbd_write_csv
  
  subroutine tMf6Wbd_read_csv(this, f_csv, nr_max)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMf6Wbd) :: this
    character(len=*), intent(in) :: f_csv
    integer(I4B), intent(in), optional :: nr_max
    ! -- local
    type(tCSV), pointer :: csv
    integer(I4B) :: ic
! ------------------------------------------------------------------------------
    !
    call this%clean()
    allocate(this%csv); csv => this%csv
    if (present(nr_max)) then
      call csv%read(f_csv, nr_max=nr_max)
    else
      call csv%read(f_csv)
    end if
    !
    ! set the data types
    csv%hdr(1)%i_type = I_C
    csv%hdr(2)%i_type = I_I8
    csv%hdr(3)%i_type = I_C
    csv%hdr(4)%i_type = I_C
    csv%hdr(5)%i_type = I_C
    csv%hdr(6)%i_type = I_I8
    csv%hdr(7)%i_type = I_I8
    !
    return
  end subroutine tMf6Wbd_read_csv
  
end module
module quad2dModule
  ! modules
  use utilsmod, only: I1B, I2B, I4B, I8B, R4B, R8B, &
    I1ZERO, I2ZERO, I4ZERO, R8ZERO, I1ONE, I4ONE, R4ZERO, R4ONE, R8ONE, I4MINONE, &
    R8HALF, R4TINY, MXSLEN, MXSLENSHORT, MXSLENLONG, &
    i_c, i_i1, i_i2, i_i4, i_i8, i_r4, i_r8,  &
    logmsg, errmsg, tBb, tBbX, tBbObj, ta, icrl_to_node, node_to_icrl, &
    tGridArr, tGrid, get_xy, get_icr, valid_icr, get_mapped_icr, point_in_bb, &
    tCSV, quicksort_d, swap_slash, open_file, bbi_intersect, bbx_intersect, &
    strip_ext, apply_mask, &
    coarse_to_fine_grid, base_name, parse_line, get_slash, create_dir, get_unique, &
    cast_number_from_string, get_ext, change_case, fileexist, get_dir, &
    quicksort_r
  use hdrModule, only: tHdr, tHdrHdr, writeflt, &
    i_uscl_arith, i_dscl_nointp, &
    uscl_names, dscl_names, n_uscl, n_dscl, i_uscl_nodata, i_dscl_nodata
  use kdtree2_module
  use vrt_module, only: tVrt, tVrtArray
  use mf6_wbd_mod, only: tMf6Wbd, i_binpos, i_asc, i_bin, MAX_NR_CSV
  use multigrid_module, only: tMultiGrid, tMultiGridArray
  use mf6_post_module, only: tPostMod
  !
  implicit none
  !
  private
  !
  logical :: LDUM
  !
  real(R8b), parameter :: R8MV = -99999999.d0
  !
  integer(I4B), parameter :: i_p  = 1   ! ic,   ir 
  integer(I4B), parameter :: i_e  = 2   ! ic+1, ir
  integer(I4B), parameter :: i_w  = 3   ! ic-1, ir
  integer(I4B), parameter :: i_s  = 4   ! ic,   ir+1
  integer(I4B), parameter :: i_n  = 5   ! ic,   ir-1
  integer(I4B), parameter :: i_ne = 6   ! ic+1, ir-1
  integer(I4B), parameter :: i_nw = 7   ! ic-1, ir-1
  integer(I4B), parameter :: i_se = 8   ! ic+1, ir+1
  integer(I4B), parameter :: i_sw = 9   ! ic-1, ir+1
  integer(I4B), parameter :: nst_fp = i_n
  integer(I4B), parameter :: nst_np = i_sw
  logical, dimension(nst_fp) :: stl_fp
  logical, dimension(nst_np) :: stl_np
  integer(I4B), dimension(2,nst_fp) :: sticir_fp
  integer(I4B), dimension(2,nst_np) :: sticir_np
  !
  type, public :: tMF6Disu
    integer(I4B) :: nodes_intf
    integer(I4B) :: nodes
    integer(I4B) :: nja
    !
    real(R8B), dimension(:), pointer    :: top     => null() ! nodes
    real(R8B), dimension(:), pointer    :: bot     => null() ! nodes
    real(R8B), dimension(:), pointer    :: area    => null() ! nodes
    integer(I4B), dimension(:), pointer :: idomain => null() ! nodes
    !
    integer(I4B), dimension(:), pointer :: ia   => null() ! nodes + 1
    integer(I4B), dimension(:), pointer :: iac  => null() ! nodes
    integer(I4B), dimension(:), pointer :: ja   => null() ! nja
    integer(I4B), dimension(:), pointer :: ihc  => null() ! nja
    real(R8B), dimension(:), pointer    :: cl12 => null() ! nja
    real(R8B), dimension(:), pointer    :: hwva => null() ! nja
    !
    type(tMf6Wbd), pointer :: wbd => null()
    !
    integer(I4B) :: ngrid
    integer(I4B) :: nlay
    integer(I4B) :: nlay_act
    integer(I4B), dimension(:), allocatable :: lay_act
    real(R8B) :: cs_min
    real(R8B) :: cs_min_rea
    real(R8B) :: cs_max_rea
    !
    ! grid mapping
    integer(I2B), dimension(:), pointer :: map_ilay  => null() ! nodes
    integer(I1B), dimension(:), pointer :: map_igrid => null() ! nodes
    integer(I4B), dimension(:), pointer :: map_nod   => null() ! nodes
    !
    integer(I4B), dimension(:,:,:), allocatable :: grid_x_nod
    real(R4B),    dimension(:,:,:), allocatable :: grid_x_top
    real(R4B),    dimension(:,:,:), allocatable :: grid_x_bot
    type(tMultiGridArray),          pointer     :: grid_mga_nod => null()
    type(tMultiGridArray),          pointer     :: grid_mga_top => null()
    type(tMultiGridArray),          pointer     :: grid_mga_bot => null()
  contains
    procedure :: init                => tMF6Disu_init
    procedure :: init_csv            => tMF6Disu_init_csv
    procedure :: clean               => tMF6Disu_clean
    procedure :: set_cs_rea_ngrid    => tMF6Disu_set_cs_rea_ngrid
    procedure :: set                 => tMF6Disu_set
    procedure :: write               => tMF6Disu_write
    procedure :: read                => tMF6Disu_read
    procedure :: x_to_mga            => tMF6Disu_x_to_mga
    procedure :: x_to_top_nodes      => tMF6Disu_x_to_top_nodes
    procedure :: x_to_top_nodes_grid => tMF6Disu_x_to_top_nodes_grid
    procedure :: mga_to_x            => tMF6Disu_mga_to_x
    procedure :: create_mga          => tMF6Disu_create_mga
    procedure :: x_get_bb            => tMF6Disu_x_get_bb
    procedure :: get_nodes           => tMF6Disu_get_nodes
  end type tMF6Disu
  
  type, public :: tVertIntf
    integer(I4B) :: lm
    integer(I4B) :: nlay_act
    integer(I4B) :: n_cells
    real(R4B)    :: mv
    !
    integer(I4B), dimension(:), allocatable :: lay_act ! length nlay_act
    !
    integer(I4B), dimension(:,:), allocatable :: nod
    real(R4B), dimension(:,:), allocatable :: zp
  contains
    procedure :: init  => tVertIntf_init
    procedure :: clean => tVertIntf_clean
    procedure :: copy  => tVertIntf_copy
  end type tVertIntf
  
  type, public :: tMF6Exchange
    integer(I4B)                            :: nexg
    integer(I4B), dimension(:), allocatable :: cellidm1
    integer(I4B), dimension(:), allocatable :: cellidm2
    integer(I4B)                            :: ihc
    real(R8B)                               :: cl1
    real(R8B)                               :: cl2
    real(R8B)                               :: hwva
  contains
    procedure :: init  => tMF6Exchange_init
    procedure :: clean => tMF6Exchange_clean
    procedure :: set   => tMF6Exchange_set
    procedure :: copy  => tMF6Exchange_copy
    procedure :: write => tMF6Exchange_write
  end type tMF6Exchange
  
  type, public :: tNbrIntf
    logical :: active
    integer(I4B) :: nbr_gid
    integer(I4B) :: nbr_lid
    integer(I4B) :: n_cells
    integer(I4B), dimension(:,:), allocatable :: lcell_nr ! from --> to (2D)
    !
    integer(I4B) :: n_exg
    integer(I4B), dimension(:,:), allocatable :: cellid ! from --> to
    !
    integer(I4B) :: vintf_active
    type(tVertIntf), pointer :: vintf_from => null()
    type(tVertIntf), pointer :: vintf_to   => null()
    !
    integer(I4B) :: xch_active
    type(tMF6Exchange), pointer :: xch => null() ! this only exists for nbr_lid > lid !
    !
  contains
    procedure :: copy                => tNbrIntf_copy
    procedure :: append              => tNbrIntf_append
    procedure :: init                => tNbrIntf_init
    procedure :: clean               => tNbrIntf_clean
    procedure :: includes_fp_stencil => tNbrIntf_includes_fp_stencil
  end type tNbrIntf
  !
  type, public :: tIntf
    integer(I4B) :: n_nbr
    integer(I4B) :: mo_nbr
    integer(I4B) :: n_mv
    integer(I4B) :: n_fill
    integer(I4B) :: ncol
    integer(I4B) :: nrow
    type(tNbrIntf), dimension(:), pointer   :: nbr_intf => null()
    integer(I4B), dimension(:,:), allocatable :: mv_lcell_nr ! missing value cells
    integer(I4B), dimension(:),   allocatable :: fill_lcell_nr ! fill cells
  contains
    procedure :: init          => tIntf_init
    procedure :: get_nbr_intf  => tIntf_get_nbr_intf
    procedure :: get_nbytes    => tIntf_get_nbytes
    procedure :: replace       => tIntf_replace
    procedure :: remove        => tIntf_remove
    procedure :: copy          => tIntf_copy
    procedure :: clean         => tIntf_clean
    procedure :: read          => tIntf_read
    procedure :: write         => tIntf_write
    procedure :: get_mask      => tIntf_get_mask
    procedure :: print         => tIntf_print
  end type tIntf
  !
  type, public :: tGraph
    integer(I8B)                            :: nvtxs
    integer(I8B)                            :: nja  ! not used for METIS !
    integer(I8B), dimension(:), allocatable :: xadj
    integer(I8B), dimension(:), allocatable :: adjncy
    integer(I8B), dimension(:), allocatable :: vwgt
    integer(I8B), dimension(:), allocatable :: adjwgt
    !
    integer(I8B)                            :: gnod_min
    integer(I8B)                            :: gnod_max
    integer(I8B), dimension(:), allocatable :: g2lnod
    integer(I8B), dimension(:), allocatable :: l2gnod
  contains
    procedure :: clean         => tGraph_clean
    procedure :: setval_g2lnod => tGraph_setval_g2lnod
    procedure :: getval_g2lnod => tGraph_getval_g2lnod
    procedure :: disconnect    => tGraph_disconnect
    procedure :: read          => tGraph_read
    procedure :: write         => tGraph_write
    procedure :: get_nbytes    => tGraph_get_nbytes
    procedure :: balance       => tGraph_balance
  end type tGraph
  !
  type, public :: tLookupTable
    character(len=MXSLEN), dimension(:), allocatable :: keys
    integer(I4B), dimension(:,:), allocatable :: table
  contains
    procedure :: init  => tLookupTable_init
    procedure :: clean => tLookupTable_clean
  end type tLookupTable
  !
  type, public :: tLayerModel
    logical :: compressed
    character(MXSLEN) :: name
    integer(I4B) :: nlay
    type(tVrt),  pointer                :: top   => null()
    type(tVrt),  dimension(:),  pointer :: bots  => null()
    type(tVrt), pointer :: comp_dat_ptr => null()
  contains
    procedure :: init           => tLayerModel_init
    procedure :: clean          => tLayerModel_clean
    procedure :: read_extent    => tLayerModel_read_extent
    procedure :: read_xy        => tLayerModel_read_xy
    procedure :: read_full_grid => tLayerModel_read_full_grid
  end type tLayerModel
  !
  type, public :: tLayerModels
    integer(I4B) :: n_inp_mod
    type(tLookupTable), pointer :: lookup => null()
    type(tLayerModel), dimension(:), pointer :: lay_mods => null()
  contains
    procedure :: init  => tLayerModels_init
    procedure :: clean => tLayerModels_clean
    procedure :: get   => tLayerModels_get
  end type tLayerModels
  !
  integer(I4B), parameter :: i_not_assigned       = 0
  integer(I4B), parameter :: i_assign_exact_layer = 1
  integer(I4B), parameter :: i_assign_first_layer = 2
  integer(I4B), parameter :: i_assign_intersect   = 3
  !
  integer(I4B), parameter :: i_type_list  = 1
  integer(I4B), parameter :: i_type_array = 2
  !
  integer(I4B), parameter :: i_vrt       = 1
  integer(I4B), parameter :: i_vrt_array = 2
  integer(I4B), parameter :: i_csv       = 3
  !
  type :: tData
    integer(I4B) :: i_type
    integer(I4B) :: i_uscl
    integer(I4B) :: i_dscl
    type(tVrt),      pointer :: vrt  => null()
    type(tVrtArray), pointer :: vrta => null()
    type(tCSV),      pointer :: csv  => null()
  end type tData
  !
  type, public :: tDataModelData
    character(len=MXSLEN) :: id
    !
    integer(I4B) :: i_type
    integer(I4B) :: i_in_file_type
    integer(I4B) :: i_out_file_type
    integer(I4B) :: i_assign
    integer(I4B) :: ilay
    integer(I4B) :: itop
    integer(I4B) :: ibot
    integer(I4B) :: idist
    !
    logical(I4B) :: lconst
    real(R8B) :: r8const
    !
    integer(I4B) :: ndat
    type(tData), dimension(:), pointer :: dat  => null()
    !
    logical :: ilay0_sfac
  contains
    procedure :: init  => tDataModelData_init
    procedure :: clean => tDataModelData_clean
    procedure :: mf6_get_data_vrt_list        => tDataModelData_mf6_get_data_vrt_list
    procedure :: mf6_get_data_vrt_array       => tDataModelData_mf6_get_data_vrt_array
    procedure :: mf6_get_data_vrt_array_array => tDataModelData_mf6_get_data_vrt_array_array
    procedure :: mf6_get_data_csv_list        => tDataModelData_mf6_get_data_csv_list
  end type tDataModelData
  !
  type, public :: tDataModel
    character(MXSLEN) :: name
    type(tCSV), pointer :: param_list => null()
    type(tCSV), pointer :: param_map  => null()
  contains
    procedure :: init     => tDataModel_init
    procedure :: clean    => tDataModel_clean
    procedure :: get_ndat => tDataModel_get_ndat
    procedure :: get_dat  => tDataModel_get_dat
  end type tDataModel
  !
  integer(I4B), parameter :: NMAX_UNI_MAP = 100000
  type, public :: tDataModels
    integer(I4B) :: n_inp_mod
    type(tDataModel), dimension(:), pointer :: dat_mods => null()
    integer(I4B) :: n_uni_map
    character(len=MXSLENSHORT), dimension(:), allocatable :: uni_map_id
    integer(I4B), dimension(:,:), allocatable :: uni_map_ir
  contains
    procedure :: init           => tDataModels_init
    procedure :: clean          => tDataModels_clean
    procedure :: create_uni_map => tDataModels_create_uni_map
  end type tDataModels
  !
  integer(I4B), parameter, public :: i_lid                = 1
  integer(I4B), parameter, public :: i_gid                = 2
  integer(I4B), parameter, public :: i_i_graph            = 3
  integer(I4B), parameter, public :: i_lay_mod            = 4
  integer(I4B), parameter, public :: i_tgt_cs_min         = 5
  integer(I4B), parameter, public :: i_tgt_cs_min_lay_beg = 6
  integer(I4B), parameter, public :: i_tgt_cs_min_lay_end = 7
  integer(I4B), parameter, public :: i_tgt_cs_min_dz_top  = 8
  integer(I4B), parameter, public :: i_tgt_cs_max         = 9
  integer(I4B), parameter, public :: i_head               = 10
  integer(I4B), parameter, public :: i_merge              = 11
  integer(I4B), parameter, public :: n_prop_field         = i_merge
  !
  type, public :: tProps
    type(tCSV), pointer                            :: csv => null() 
    character(len=MXSLEN), dimension(n_prop_field) :: fields
  contains
    procedure :: init      => tProps_init
    procedure :: clean     => tProps_clean
    procedure :: get_field => tProps_get_field
  end type tProps
  !
  type, public :: tQuads
    character(len=MXSLEN) :: dir
    character(len=MXSLEN) :: uuid
    !
    type(tProps), pointer               :: props => null() 
    integer(I4B) :: n
    integer(I4B) :: n_act
    type(tQuad), dimension(:), pointer  :: xq => null()
    type(tGraph), pointer               :: qg  => null()
    type(tGraph), dimension(:), pointer :: qgd => null() ! disconnected graphs
    type(tLayerModels), pointer         :: lay_mods => null()
    type(tDataModels),  pointer         :: dat_mods => null()
  contains
    procedure :: init              => tQuads_init
    !procedure :: init_select       => tQuads_init_select
    procedure :: init_select_bb    => tQuads_init_select_bb
    procedure :: init_select_intf  => tQuads_init_select_intf
    procedure :: init_select_graph => tQuads_init_select_graph
    procedure :: generate_uuid     => tQuads_generate_uuid
    procedure :: clean             => tQuads_clean
    procedure :: get_quad          => tQuads_get_quad
    procedure :: get               => tQuads_get
    procedure :: get_file          => tQuads_get_file
    procedure :: get_number_active => tQuads_get_number_active
    !
    procedure :: construct_graph  => tQuads_construct_graph
    procedure :: disconnect_graph => tQuads_disconnect_graph
    !
    procedure :: join             => tQuads_join
    procedure :: split            => tQuads_split
    !
    !procedure :: read_mask        => tQuads_read_mask
    !procedure :: write_mask       => tQuads_write_mask
    procedure :: read_bb          => tQuads_read_bb
    procedure :: write_bb         => tQuads_write_bb
    procedure :: reorder_intf     => tQuads_reorder_intf
    procedure :: read_intf        => tQuads_read_intf
    procedure :: write_intf       => tQuads_write_intf
    procedure :: write_props      => tQuads_write_props
    procedure :: read_props       => tQuads_read_props
    procedure :: get_props_ptr    => tQuads_get_props_ptr
    procedure :: write_graphs     => tQuads_write_graphs
    procedure :: read_graph       => tQuads_read_graph
    procedure :: balance_graphs   => tQuads_balance_graphs
    procedure :: open_file        => tQuads_open_file
    procedure :: add_lm_intf      => tQuads_add_lm_intf
    procedure :: write_mf6_xch_intf   => tQuads_write_mf6_xch_intf
    procedure :: write_mf6_heads      => tQuads_write_mf6_heads
    procedure :: write_mf6_bnd_heads  => tQuads_write_mf6_bnd_heads
    procedure :: set_mod_dir          => tQuads_set_mod_dir
    procedure :: merge_disu           => tQuads_merge_disu
    !
  end type tQuads
  !
  type, public :: tQuad
    logical :: mapped_gids
    logical :: active
    logical :: work_done
    logical :: generated
    !
    integer(I4B) :: gid
    integer(I4B) :: lid
    !
    integer(I4B)                            :: gid_prent
    integer(I4B)                            :: lid_prent
    integer(I4B)                            :: n_child
    integer(I4B), dimension(:), allocatable :: gid_child
    integer(I4B), dimension(:), allocatable :: lid_child 
    !
    type(tProps), pointer :: props => null() 
    integer(I4B) :: area ! cell count
    real(R8B)    :: weight
    real(R8B)    :: bb_ratio
    real(R8B)    :: bb_xm
    real(R8B)    :: bb_ym
    integer(i4B) :: ntile
    integer(I4B), dimension(:), allocatable :: itile
    integer(I4B) :: itile_bb_xm
    integer(I4B) :: i_graph
    integer(I4B) :: i_prop
    integer(I4B), dimension(:), allocatable :: hlev
    !
    type(tBbObj) :: bbo
    !
    type(tHdr), pointer :: hdr_gids  => null()
    integer(I4B) :: gids_mv_map
    real(R8B) :: cs_max ! maximum coarsening resolution
    !
    type(tIntf), pointer :: intf => null()
    !
    type(tMultiGrid), pointer :: c2f_grid => null()
    type(tMultiGrid), pointer :: f2c_grid => null()
    !
    type(tLayerModels), pointer :: lay_mods => null()
    type(tDataModel),  pointer  :: dat_mod  => null()
    type(tMF6Disu), pointer     :: disu => null()
    !
    character(len=MXSLEN) :: mod_dir
  contains
    procedure :: init            => tQuad_init
    procedure :: init_map_gids   => tQuad_init_map_gids
    procedure :: clean           => tQuad_clean
    !
    procedure :: map_gids        => tQuad_map_gids
    procedure :: clean_gids      => tQuad_clean_gids
    procedure :: calc_prop       => tQuad_calc_prop
    procedure :: get_prop        => tQuad_get_prop
    procedure :: get_area        => tQuad_get_area
    procedure :: get_hlev        => tQuad_get_hlev
    procedure :: get_bb_ratio    => tQuad_get_bb_ratio
    procedure :: calc_interface  => tQuad_calc_interface
    !
    procedure :: set             => tQuad_set
    procedure :: set_cgrid       => tQuad_set_cgrid
    procedure :: clean_cgrid     => tQuad_clean_cgrid
    procedure :: set_flag        => tQuad_set_flag
    procedure :: grid_gen        => tQuad_grid_gen
    procedure :: grid_init       => tQuad_grid_init
    procedure :: mf6_write_data  => tQuad_mf6_write_data
    !
    procedure :: get               => tQuad_get
    procedure :: get_grid          => tQuad_get_grid
    procedure :: get_grid_dim      => tQuad_get_grid_dim
    procedure :: get_grid_bbx      => tQuad_get_grid_bbx
    procedure :: get_mo_nbr        => tQuad_get_mo_nbr
    procedure :: get_bb            => tQuad_get_bb
    !procedure :: get_grid_count  => tQuad_get_grid_count
    procedure :: get_flag        => tQuad_get_flag
    procedure :: get_intf        => tQuad_get_intf_ptr
    procedure :: get_prop_csv    => tQuad_get_prop_csv
    procedure :: set_prop_csv    => tQuad_set_prop_csv
    
    !procedure :: store_mask      => tQuad_store_mask
    !procedure :: store_gid_mask  => tQuad_store_gid_mask
    procedure :: check_gid_mask  => tQuad_check_gid_mask
    procedure :: get_nod_map     => tQuad_get_nod_map
    procedure :: mf6_get_data    => tQuad_mf6_get_data
  end type tQuad

  ! expose
  !public :: tQuads, tQuad, tNbrIntf, tIntf
  !public :: tLayerModels, tLayerModel, tLookupTable
  public :: get_number_of_levels, get_refinement_level
  public :: mf6_data_write
  public :: valid_icir
  
  save
  
  contains
  
! ==============================================================================
! ==============================================================================
! tMF6Disu
! ==============================================================================
! ==============================================================================
  
   subroutine tMF6Disu_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    ! -- local
! ------------------------------------------------------------------------------
    call this%clean
    !
    this%nodes_intf = I4ZERO
    this%nodes      = I4ZERO
    this%nja        = I4ZERO
    !
    this%ngrid      = I4ZERO
    this%nlay       = I4ZERO
    this%nlay_act   = I4ZERO
    this%cs_min     = R8ZERO
    this%cs_min_rea = R8ZERO
    this%cs_max_rea = R8ZERO
    !
    return
   end subroutine tMF6Disu_init
   
   subroutine tMF6Disu_init_csv(this, f_csv, f_bin, reuse)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    character(len=*), intent(in) :: f_csv
    character(len=*), intent(in), optional :: f_bin
    logical, intent(in), optional :: reuse
    ! -- local
    logical :: reuse_loc
! ------------------------------------------------------------------------------
    !
    if (present(reuse)) then
      reuse_loc = reuse
    else
      reuse_loc = .false.
    end if
    !
    allocate(this%wbd)
    if (present(f_bin)) then
      call this%wbd%init(f_csv=f_csv, f_binpos=f_bin, reuse=reuse_loc)
    else
      call this%wbd%init(f_csv=f_csv, reuse=reuse_loc)
    end if
    !
    return
   end subroutine tMF6Disu_init_csv
  
   subroutine tMF6Disu_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (associated(this%top))       deallocate(this%top)
    if (associated(this%bot))       deallocate(this%bot)
    if (associated(this%area))      deallocate(this%area)
    if (associated(this%idomain))   deallocate(this%idomain)
    if (associated(this%ia))        deallocate(this%ia)
    if (associated(this%iac))       deallocate(this%iac)
    if (associated(this%ja))        deallocate(this%ja)
    if (associated(this%ihc))       deallocate(this%ihc)
    if (associated(this%cl12))      deallocate(this%cl12)
    if (associated(this%hwva))      deallocate(this%hwva)
    if (associated(this%map_ilay))  deallocate(this%map_ilay)
    if (associated(this%map_igrid)) deallocate(this%map_igrid)
    if (associated(this%map_nod))   deallocate(this%map_nod)
    !
    this%top       => null()
    this%bot       => null()
    this%area      => null()
    this%idomain   => null()
    this%ia        => null()
    this%iac       => null()
    this%ja        => null()
    this%ihc       => null()
    this%cl12      => null()
    this%hwva      => null()
    this%map_ilay  => null()
    this%map_igrid => null()
    this%map_nod   => null()
    !
    if (associated(this%wbd)) then
      call this%wbd%clean()
      deallocate(this%wbd); this%wbd => null()
    end if
    if (allocated(this%grid_x_nod)) then
      deallocate(this%grid_x_nod)
    end if
    if (associated(this%grid_mga_nod)) then
      call this%grid_mga_nod%clean()
      deallocate(this%grid_mga_nod)
      this%grid_mga_nod => null()
    end if
    if (associated(this%grid_mga_top)) then
      call this%grid_mga_top%clean()
      deallocate(this%grid_mga_top)
      this%grid_mga_top => null()
    end if
    if (associated(this%grid_mga_bot)) then
      call this%grid_mga_bot%clean()
      deallocate(this%grid_mga_bot)
      this%grid_mga_bot => null()
    end if
    if (allocated(this%lay_act)) then
      deallocate(this%lay_act)
    end if
    !
    return
  end subroutine tMF6Disu_clean
  
  function tMF6Disu_x_get_bb(this, x) result(bba)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    class(tMF6Disu) :: this
    integer(I4B), dimension(:,:,:), intent(in) :: x
    type(tBB), dimension(:), pointer :: bba
    ! -- dummy
    ! -- local
    type(tBB), pointer :: bb => null()
    integer(I4B) :: nl, nr, nc, nodes, il, ir, ic, nod, i, nc_max, nr_max
! ------------------------------------------------------------------------------
    nc = size(x,1); nr = size(x,2); nl = size(x,3)
    nodes = maxval(x)

    ! determine the bounding box
    allocate(bba(nodes))
    !
    do i = 1, nodes
      call bba(i)%init()
    end do
    nc_max = 0; nr_max = 0
    do il = 1, nl; do ir = 1, nr; do ic = 1, nc
      nod = x(ic,ir,il)
      if (nod > 0) then
        bb => bba(nod)
        bb%ic0 = min(ic, bb%ic0); bb%ic1 = max(ic, bb%ic1)
        bb%ir0 = min(ir, bb%ir0); bb%ir1 = max(ir, bb%ir1)
        bb%ncol = bb%ic1 - bb%ic0 + 1; bb%nrow = bb%ir1 - bb%ir0 + 1
        nc_max = max(nc_max, bb%ncol)
        nr_max = max(nr_max, bb%nrow)
      end if
    end do; end do; end do
    !
    ! checks
    do i = 1, nodes
      bb => bba(i)
      if (bb%ncol /= bb%nrow) then
        call errmsg('tMF6Disu_x_get_bb: program error 1.')
      end if
      if (bb%ncol > 1) then
        if (mod(bb%ncol,2) /= 0) then
          call errmsg('tMF6Disu_x_get_bb: program error 2.')
        end if
      end if
    end do
    !
    return
  end function tMF6Disu_x_get_bb
  
  subroutine tMF6Disu_x_to_mga(this, x, ngrid, bbx)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    integer(I4B), dimension(:,:,:), intent(in) :: x
    integer(I4B), intent(in) :: ngrid
    type(tBbx), intent(in) :: bbx
    ! -- local
    type(tMultiGridArray), pointer :: mga => null()
    type(tMultiGrid), pointer      :: mg  => null()
    type(tGrid),      pointer      :: g   => null()
    type(tBB), dimension(:), pointer :: bba => null()
    type(tBB), pointer :: bb => null()
    type(tBBx) :: bbx_set
    integer(I4B), dimension(:), allocatable :: flg
    integer(I4B) :: nl, nr, nc, mr, mc, nodes, il, ir, ic, jr, jc
    integer(I4B) :: ig, i, nod, bs, ig_max_I1, il_max_I2, ig_max
    real(R8B) :: xc, yc
    real(R8B), dimension(:), allocatable :: csa
! ------------------------------------------------------------------------------
    !
    if (allocated(this%grid_x_nod)) then
      call errmsg('tMF6Disu_x_to_mga: grid_x is already being used.')
    end if
    allocate(this%grid_x_nod, source=x)
    !
    nc = size(x,1); nr = size(x,2); nl = size(x,3)
    nodes = maxval(x)
    !
    ! determine the levels and layers
    allocate(this%map_ilay(nodes), this%map_igrid(nodes), this%map_nod(nodes))
    il_max_I2 = huge(this%map_ilay(1))
    ig_max_I1 = huge(this%map_igrid(1))
    !
    this%map_ilay  = 0_I2B
    this%map_igrid = 0_I1B
    this%map_nod   = 0_I4B
    !
    ! get the bounding box
    bba => this%x_get_bb(x)
    !
    ! determine map_ilay and map_igrid
    ig_max = 0
    do il = 1, nl; do ir = 1, nr; do ic = 1, nc
      nod = x(ic,ir,il)
      if (nod > 0) then
        bb => bba(nod)
        if (il <= il_max_I2) then
          this%map_ilay(nod) = int(il, I2B)
        else
          call errmsg('tMF6Disu_x_to_mga: ilay out of range: '//ta([il]))
        end if
        ig = ngrid - int(log(real(bb%ncol, R8B))/log(2.d0), I4B)
        ig_max = max(ig, ig_max)
        if (ig <= ig_max_I1) then
          if (ig <= 0) then
            call errmsg('tMF6Disu_x_to_mga: program error.')
          end if
          this%map_igrid(nod) = int(ig, I1B)
        else
          call errmsg('tMF6Disu_x_to_mga: igrid out of range.')
        end if
      end if
    end do; end do; end do
    !
    ! check
    if (ig_max > this%ngrid) then
      call errmsg('tMF6Disu_x_to_mga: program error, ig_max > ngrid.')
    end if
    !
    if (associated(bba)) deallocate(bba)
    bba => null()
    !
    allocate(this%grid_mga_nod); mga => this%grid_mga_nod
    call mga%init(nl)
    !
    allocate(csa(this%ngrid))
    do i = 1, ngrid
      ig = ngrid - i + 1
      bs = 2**(ngrid-ig)
      if (ig <= this%ngrid) then
        csa(ig) = bbx%cs * bs
      end if
    end do
    do il = 1, nl
      mg => mga%get_mg(il)
      call mg%init(this%ngrid, xll=bbx%xll, xur=bbx%xur, &
        yll=bbx%yll, yur=bbx%yur, csa=csa, mvi4=0)
    end do
    deallocate(csa)
    !
    !do il = 1, nl
    !  mg => mga%get_mg(il)
    !  call mg%init(ngrid, xll=bbx%xll, xur=bbx%xur, yll=bbx%yll, yur=bbx%yur)
    !  ! coarse --> fine
    !  do i = 1, ngrid
    !    mc = nc/(2**(i-1)); mr = nr/(2**(i-1))
    !    ig = ngrid - i + 1
    !    bs = 2**(ngrid-ig)
    !    g => mg%grid(ig)
    !    bbx_set = bbx; bbx_set%cs = bbx%cs * bs
    !    call g%init(nc=mc, nr=mr, mvi4=0, bbx=bbx_set)
    !  end do
    !end do
    
    !
    if (allocated(flg)) deallocate(flg)
    allocate(flg(nodes)); flg = 0
    do il = 1, nl
      mg => mga%get_mg(il)
      do ig = 1, this%ngrid
        g => mg%grid(ig)
        call g%set_const(i4v=0)
      end do
      ! loop over the finest level
      do ir = 1, nr; do ic = 1, nc
        nod = x(ic,ir,il)
        if (nod > 0) then
          if (flg(nod) == 0) then
            ig = int(this%map_igrid(nod), I4B)
            bs = 2**(ngrid-ig)
            !
            ! get the coarser icol/irow
            call get_xy( xc, yc, ic, ir, bbx%xll, bbx%yur, bbx%cs)
            call get_icr(jc, jr, xc, yc, bbx%xll, bbx%yur, bbx%cs*bs)
            !
            ! store the node number
            g => mg%grid(ig); g%xi4(jc,jr) = nod
            call icrl_to_node(this%map_nod(nod), jc, jr, 1, g%nc, g%nr)
            !
            ! set the flag
            flg(nod) = 1
          end if
        end if
      end do; end do
    end do
    !
    if (allocated(flg)) deallocate(flg)
    !
    return
  end subroutine tMF6Disu_x_to_mga
   
  subroutine tMF6Disu_x_to_top_nodes(this, top_nodes)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    integer(I1B), dimension(:), allocatable, intent(inout) :: top_nodes
    ! -- local
    integer(I4B) :: il, ir, ic, nl, nr, nc, n
! ------------------------------------------------------------------------------
    if (.not.allocated(this%grid_x_nod)) then
      call errmsg('tMF6Disu_x_to_top_nodes: program error, '// &
        'could not find grid_x_nod.')
    end if
    !
    if (allocated(top_nodes)) deallocate(top_nodes)
    allocate(top_nodes(this%nodes))
    top_nodes = 0
    !
    nc = size(this%grid_x_nod,1)
    nr = size(this%grid_x_nod,2)
    nl = size(this%grid_x_nod,3)
    !
    do ir = 1, nr; do ic = 1, nc
      do il = 1, nl
        n = this%grid_x_nod(ic,ir,il)
        if (n > 0) then
          top_nodes(n) = 1
          exit
        end if
      end do
    end do; end do
    !
    return
  end subroutine tMF6Disu_x_to_top_nodes
  
  subroutine tMF6Disu_x_to_top_nodes_grid(this, top_nodes_grid)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    integer(I4B), dimension(:,:), allocatable, intent(inout) :: top_nodes_grid
    ! -- local
    integer(I4B) :: il, ir, ic, nl, nr, nc, n
! ------------------------------------------------------------------------------
    if (.not.allocated(this%grid_x_nod)) then
      call errmsg('tMF6Disu_x_to_top_nodes_grid: program error, '// &
        'could not find grid_x_nod.')
    end if
    !
    nc = size(this%grid_x_nod,1)
    nr = size(this%grid_x_nod,2)
    nl = size(this%grid_x_nod,3)
    !
    if (allocated(top_nodes_grid)) deallocate(top_nodes_grid)
    allocate(top_nodes_grid(nc,nr))
    top_nodes_grid = 0
    !
    do ir = 1, nr; do ic = 1, nc
      do il = 1, nl
        n = this%grid_x_nod(ic,ir,il)
        if (n > 0) then
          top_nodes_grid(ic,ir) = n
          exit
        end if
      end do
    end do; end do
    !
    return
  end subroutine tMF6Disu_x_to_top_nodes_grid
  
  subroutine tMF6Disu_create_mga(this, bbi, bbx, cs_min_tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    type(tBb), intent(in) :: bbi
    type(tBbx), intent(in) :: bbx
    real(R8B), intent(in) :: cs_min_tgt
    ! -- local
    type(tMultiGridArray), pointer :: mga_nod => null()
    type(tMultiGridArray), pointer :: mga_top => null()
    type(tMultiGridArray), pointer :: mga_bot => null()
    type(tMultiGrid), pointer      :: mg  => null()
    type(tGrid), pointer           :: g   => null()
    integer(I4B) :: ierr, ig, il, n, i, nc, nr, nl, mr, mc, ir, ic
    integer(I4B) :: idummy, f, bs
    real(R8B) :: cs
    real(R8B), dimension(:), allocatable :: csa
    real(R4B) :: top, bot
! ------------------------------------------------------------------------------
    !
    nc = bbi%ncol; nr = bbi%nrow
    if (cs_min_tgt < bbx%cs) then
      f = int(bbx%cs/cs_min_tgt, I4B)
      nc = nc*f; nr = nr*f
    end if
    
    ! check is everything that is used was set
    ierr = 0
    if (this%nodes == I4ZERO)    ierr = 1
    if (this%ngrid == I4ZERO)    ierr = 1
    if (this%nlay_act == I4ZERO) ierr = 1
    if (.not.associated(this%top))       ierr = 1
    if (.not.associated(this%bot))       ierr = 1
    if (.not.associated(this%map_ilay))  ierr = 1
    if (.not.associated(this%map_igrid)) ierr = 1
    if (.not.associated(this%map_nod))   ierr = 1
    if (ierr == 1) then
      call errmsg('tMF6Disu_create_mga: missing data.')
    end if
    !
    if (associated(this%grid_mga_nod)) then
      call this%grid_mga_nod%clean()
      deallocate(this%grid_mga_nod); this%grid_mga_nod => null()
    end if
    if (associated(this%grid_mga_top)) then
      call this%grid_mga_top%clean()
      deallocate(this%grid_mga_top); this%grid_mga_top => null()
    end if
    if (associated(this%grid_mga_top)) then
      call this%grid_mga_bot%clean()
      deallocate(this%grid_mga_bot); this%grid_mga_bot => null()
    end if
    !
    ! initialize the grid
    allocate(this%grid_mga_nod); mga_nod => this%grid_mga_nod
    allocate(this%grid_mga_top); mga_top => this%grid_mga_top
    allocate(this%grid_mga_bot); mga_bot => this%grid_mga_bot
    !
    call mga_nod%init(this%nlay_act)
    call mga_top%init(this%nlay_act)
    call mga_bot%init(this%nlay_act)
    
    allocate(csa(this%ngrid))
    do i = 1, this%ngrid
      ig = this%ngrid - i + 1
      bs = 2**(this%ngrid-ig)
      cs = cs_min_tgt * bs
      csa(ig) = cs
    end do
    !
    do il = 1, this%nlay_act
      mg => mga_nod%get_mg(il)
      call mg%init(this%ngrid, xll=bbx%xll, xur=bbx%xur, &
        yll=bbx%yll, yur=bbx%yur, csa=csa, mvi4=0)
      !
      mg => mga_top%get_mg(il)
      call mg%init(this%ngrid, xll=bbx%xll, xur=bbx%xur, &
        yll=bbx%yll, yur=bbx%yur, csa=csa, mvr4=0.)
      !
      mg => mga_bot%get_mg(il)
      call mg%init(this%ngrid, xll=bbx%xll, xur=bbx%xur, &
        yll=bbx%yll, yur=bbx%yur, csa=csa, mvr4=0.)
      
    end do
    deallocate(csa)
    !
    !do il = 1, this%nlay_act
    !  mg => mga%get_mg(il)
    !  call mg%init(this%ngrid, xll=bbx%xll, xur=bbx%xur, yll=bbx%yll, yur=bbx%yur)
    !  ! coarse --> fine
    !  do i = 1, this%ngrid
    !    mc = nc/(2**(i-1)); mr = nr/(2**(i-1))
    !    ig = this%ngrid - i + 1
    !    g => mg%grid(ig)
    !    call g%init(nc=mc, nr=mr, mvi4=0)
    !  end do
    !end do
    
    ! fill the grid
    do i = 1, this%nodes
      il = this%map_ilay(i); ig = this%map_igrid(i)
      n = this%map_nod(i); top = this%top(i); bot = this%bot(i)
      if ((il <= 0).or.(il > this%nlay_act)) then
        call errmsg('tMF6Disu_create_mga: program error map_ilay.')
      end if
      if ((ig <= 0).or.(ig > this%ngrid)) then
        call errmsg('tMF6Disu_create_mga: program error map_igrid.')
      end if
      if (n <= 0) then
        call errmsg('tMF6Disu_create_mga: program error map_nod.')
      end if
      !
      mg => mga_nod%get_mg(il); g => mg%grid(ig)
      call node_to_icrl(n, ic, ir, idummy, g%nc, g%nr)
      !
      if ((ic <= 0).or.(ic > g%nc).or. &
          (ir <= 0).or.(ir > g%nr)) then
        call errmsg('tMF6Disu_create_mga: program error.')
      end if
      g%xi4(ic,ir) = i
      !
      mg => mga_top%get_mg(il); g => mg%grid(ig); g%xr4(ic,ir) = top
      mg => mga_bot%get_mg(il); g => mg%grid(ig); g%xr4(ic,ir) = bot
    end do
    !
    ! set the counts
    do il = 1, this%nlay_act
      mg => mga_nod%get_mg(il); call mg%set_grid_count()
      mg => mga_top%get_mg(il); call mg%set_grid_count()
      mg => mga_bot%get_mg(il); call mg%set_grid_count()
    end do
    !
    return
  end subroutine tMF6Disu_create_mga
  
  subroutine tMF6Disu_mga_to_x(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    ! -- local
    type(tMultiGridArray), pointer :: mga => null()
    type(tMultiGrid), pointer      :: mg  => null()
    type(tGrid), pointer           :: g   => null()
    type(tBBx) :: bbx
    integer(I4B) :: nl, nr, nc, ngrid, bs, nod, ic0
    integer(I4B) :: ic1, ir0, ir1, il, ir, ic, ig, jr, jc, mc, mr
    real(R8B) :: xc, yc, xul, yul, xlr, ylr, cs_min, cs, quar_cs_min
! ------------------------------------------------------------------------------
    mga => this%grid_mga_nod
    if (.not.associated(mga)) then
      call errmsg('tMF6Disu_mga_to_x: grid_mga not found.')
    end if
    !
    mg => mga%get_mg(1)
    ngrid = mg%ngrid; g => mg%grid(ngrid); bbx = g%bbx
    cs_min = bbx%cs
    quar_cs_min = 0.25d0*cs_min
    nl = mga%nmgrid; nr = g%nr; nc = g%nc
    !
    if (allocated(this%grid_x_nod)) deallocate(this%grid_x_nod)
    allocate(this%grid_x_nod(nc,nr,nl))
    this%grid_x_nod = 0
    !
    do il = 1, nl
      mg => mga%get_mg(il)
      bs = 1
      do ig = 1, ngrid
        bs = 2**(ngrid-ig)
        g => mg%grid(ig); bbx = g%bbx; cs = bbx%cs
        do ir = 1, g%nr; do ic = 1, g%nc
          nod = g%xi4(ic,ir)
          if (nod /= g%mvi4) then
            call get_xy(xc, yc, ic, ir, bbx%xll, bbx%yur, cs)
            !
            xul = xc - R8HALF*cs + quar_cs_min
            yul = yc + R8HALF*cs - quar_cs_min
            xlr = xc + R8HALF*cs - quar_cs_min
            ylr = yc - R8HALF*cs + quar_cs_min
            !
            call get_icr(ic0, ir0, xul, yul, bbx%xll, bbx%yur, cs_min)
            call get_icr(ic1, ir1, xlr, ylr, bbx%xll, bbx%yur, cs_min)
            !
            mc = ic1 - ic0 + 1; mr = ir1 - ir0 + 1
            if ((mc /= bs).or.(mr /= bs)) then
              call errmsg('tMF6Disu_mga_to_x: program error.')
            end if
            do jr = ir0, ir1; do jc = ic0, ic1
              this%grid_x_nod(jc,jr,il) = nod
            end do; end do;
          end if
        end do; end do
      end do
    end do
    !
    return
  end subroutine tMF6Disu_mga_to_x
  
  subroutine tMF6Disu_get_nodes(this, ilay, xp, yp, nod, x, y, cs, top, bot)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    !
    integer(I4B), intent(in), optional :: ilay
    real(R8B), dimension(:), intent(in), optional :: xp
    real(R8B), dimension(:), intent(in), optional :: yp
    !
    integer(I4B), dimension(:), allocatable, intent(inout) :: nod
    real(R8B), dimension(:), allocatable, intent(inout) :: x
    real(R8B), dimension(:), allocatable, intent(inout) :: y
    real(R8B), dimension(:), allocatable, intent(inout) :: cs
    real(R4B), dimension(:), allocatable, intent(inout) :: top
    real(R4B), dimension(:), allocatable, intent(inout) :: bot
! ------------------------------------------------------------------------------
    !
    ! clean
    if (allocated(nod)) deallocate(nod)
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(cs)) deallocate(cs)
    if (allocated(top)) deallocate(top)
    if (allocated(bot)) deallocate(bot)
    !
    return
  end subroutine tMF6Disu_get_nodes
  
  subroutine tMF6Disu_set_cs_rea_ngrid(this, x, cs_min)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    integer(I4B), dimension(:,:,:), intent(in) :: x
    real(R8B), intent(in) :: cs_min
    ! -- local
    type(tBB), pointer :: bb => null()
    type(tBB), dimension(:), pointer :: bba => null()
    integer(I4B) :: i, nodes
    real(R8B) :: cs
! ------------------------------------------------------------------------------

    ! get the bounding box
    bba => this%x_get_bb(x); nodes = maxval(x)
    !
    this%cs_min_rea = huge(R8ONE); this%cs_max_rea = R8ZERO
    !
    do i = 1, nodes
      bb => bba(i)
      cs = bb%ncol*cs_min
      !
      this%cs_min_rea = min(this%cs_min_rea,cs)
      this%cs_max_rea = max(this%cs_max_rea,cs)
    end do
    !
    ! realised ngrid:
    this%ngrid = int(log(this%cs_max_rea/this%cs_min_rea)/log(2.d0),I4B) + 1
    !
    ! clean up
    if (associated(bba)) deallocate(bba)
    !
    return
  end subroutine tMF6Disu_set_cs_rea_ngrid
  
  subroutine tMF6Disu_set(this, zp, cs_min)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    real(R4B), dimension(:,:,:), intent(in) :: zp
    real(R8B), intent(in) :: cs_min
    !
    ! -- local
    type(tBB), pointer :: bb => null(), bb_nbr => null()
    type(tBB), dimension(:), pointer :: bba => null()
    real(R8B) :: cl12, hwva, cl12_hor, cl12_ver, hwva_hor, hwva_ver
    real(R4B) :: top, bot, thk, cs, cs_nbr, cs_min_rea, cs_max_rea
    logical :: lfound
    integer(I4B), dimension(:), allocatable :: i4wk, i4wk2
    integer(I4B) :: nl, nr, nc, il, jl, ir, ic, i, j, nod, ig_min, ig_max, nbr, inbr
    integer(I4B) :: ic0, ic1, ir0, ir1, il0, il1, n, iact
    integer(I4B) :: ihc, ihc_hor, ihc_ver
    !
    integer(I4B), dimension(:,:,:), allocatable :: x ! work
! ------------------------------------------------------------------------------
    !
    if (.not.allocated(this%grid_x_nod)) then
      call errmsg('tMF6Disu_set: grid_x not found.')
    end if
    allocate(x, source=this%grid_x_nod)
    !
    this%nodes = maxval(x)
    !
    ! get the bounding box
    bba => this%x_get_bb(x)
    !
    nc = size(x,1); nr = size(x,2); nl = size(x,3)
    if ((size(zp,1) /= nc).or.(size(zp,2) /= nr).or.(size(zp,3) /= nl+1)) then
      call errmsg('tMF6Disu_set: inconsistent dimensions.')
    end if
    !
    ! determine nja
    allocate(i4wk(nc*nr*nl))
    !
    do iact = 1, 2
      this%nja = 0
      do i = 1, this%nodes
        bb => bba(i)
        cs = bb%ncol*cs_min
        !
        il = int(this%map_ilay(i), I4B)
        !
        ! get the neighbors
        il0 = il - 1;     il1 = il + 1
        ic0 = bb%ic0 - 1; ic1 = bb%ic1 + 1
        ir0 = bb%ir0 - 1; ir1 = bb%ir1 + 1
        !
        n = 0
        if (ic0 >= 1) then !W
          do ir = bb%ir0, bb%ir1
            if (x(ic0,ir,il) > 0) then
              n = n + 1; i4wk(n) = x(ic0,ir,il)
            end if
          end do
        end if
        if (ic1 <= nc) then !E
          do ir = bb%ir0, bb%ir1
            if (x(ic1,ir,il) > 0) then
              n = n + 1; i4wk(n) = x(ic1,ir,il)
            end if
          end do
        end if
        if (ir0 >= 1) then !N
          do ic = bb%ic0, bb%ic1
            if (x(ic,ir0,il) > 0) then
              n = n + 1; i4wk(n) = x(ic,ir0,il)
            end if
          end do
        end if
        if (ir1 <= nr) then !S
          do ic = bb%ic0, bb%ic1
            if (x(ic,ir1,il) > 0) then
              n = n + 1; i4wk(n) = x(ic,ir1,il)
            end if
          end do
        end if
        !if (il0 >= 1) then !T
        !  do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        !    if (x(ic,ir,il0) > 0) then
        !      n = n + 1; i4wk(n) = x(ic,ir,il0)
        !    end if
        !  end do; end do
        !end if
        !if (il1 <= nl) then !B
        !  do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        !    if (x(ic,ir,il1) > 0) then
        !      n = n + 1; i4wk(n) = x(ic,ir,il1)
        !    end if
        !  end do; end do
        !end if
        !
        if (il0 >= 1) then !T
          do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
            do jl = il0, 1, -1
              if (x(ic,ir,jl) > 0) then
                n = n + 1; i4wk(n) = x(ic,ir,jl); exit
              end if
            end do
          end do; end do
        end if
        if (il1 <= nl) then !B
          do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
            do jl = il1, nl
              if (x(ic,ir,jl) > 0) then
                n = n + 1; i4wk(n) = x(ic,ir,jl); exit
              end if
            end do
          end do; end do
        end if
        !
        if (n > 0) then
          call get_unique(i4wk(1:n), i4wk2)
          nbr = size(i4wk2)
        else
          nbr = 0
        end if
        !
        if (iact == 2) then
          !
          this%iac(i) = 1 + nbr
          this%nja = this%nja + 1
          !
          ! center point
          this%ja(this%nja)   = i
          this%ihc(this%nja)  = 0 ! dummy
          this%cl12(this%nja) = R8ZERO ! dummy
          this%hwva(this%nja) = R8ZERO ! dummy
          !
          ! determine parameters
          ihc_hor = 1; ihc_ver = 0
          top = zp(bb%ic0, bb%ir0, il)
          bot = zp(bb%ic0, bb%ir0, il+1)
          thk = top - bot
          if (abs(thk) < R4TINY) then
            call errmsg('tMF6Disu_set: program error.')
          end if
          cl12_hor = bb%ncol*cs_min*R8HALF
          cl12_ver = real(thk,R8B)*R8HALF
          
          ! set nodal data
          this%top(i)  = top
          this%bot(i)  = bot
          this%area(i) = cs**2
          this%iac(i)  = 1 + nbr
          this%ia(i+1) = this%ia(i) + nbr + 1
          if (nbr == 0) then
            call logmsg('=====> setting idomain=0 for node '//ta([i])//' <=====')
            this%idomain(i) = 0
          end if
          !
          ! set graph data
          do inbr = 1, nbr
            nod = i4wk2(inbr)
            if (nod == 0) then
              call errmsg('tMF6Disu_set: program error.')
            end if
            bb_nbr => bba(nod)
            cs_nbr = bb_nbr%ncol*cs_min
            this%nja = this%nja + 1
            !
            lfound = .false.
            if (ic0 == bb_nbr%ic1) then !W
              ihc = 1; cl12 = cl12_hor
              hwva = min(cs, cs_nbr)
              lfound = .true.
            end if
            if (ic1 == bb_nbr%ic0) then !E
              ihc = 1; cl12 = cl12_hor
              hwva = min(cs, cs_nbr)
              lfound = .true.
            end if
            if (ir0 == bb_nbr%ir1) then !N
              ihc = 1; cl12 = cl12_hor
              hwva = min(cs, cs_nbr)
              lfound = .true.
            end if
            if (ir1 == bb_nbr%ir0) then !S
              ihc = 1; cl12 = cl12_hor
              hwva = min(cs, cs_nbr)
              lfound = .true.
            end if
            if (il > int(this%map_ilay(nod),I4B)) then !T
            !if (il0 == int(this%map_ilay(nod),I4B)) then !T
              ihc = 0; cl12 = cl12_ver
              hwva = min(cs, cs_nbr)**2
              lfound = .true.
            end if
            if (il < int(this%map_ilay(nod),I4B)) then !B
            !if (il1 == int(this%map_ilay(nod),I4B)) then !B
              ihc = 0; cl12 = cl12_ver
              hwva = min(cs, cs_nbr)**2
              lfound = .true.
            end if
            if (.not.lfound) then
              call errmsg('tMF6Disu_set: program error.')
            end if
            !
            ! set the data
            this%ja(this%nja)   = nod
            this%ihc(this%nja)  = ihc
            this%cl12(this%nja) = cl12
            this%hwva(this%nja) = hwva
          end do
        else ! iact = 1
          this%nja = this%nja + 1 + nbr
        end if
      end do
      !
      if (iact == 1) then
        allocate(this%top(this%nodes))
        allocate(this%bot(this%nodes))
        allocate(this%area(this%nodes))
        allocate(this%idomain(this%nodes))
        this%idomain = 1
        !
        allocate(this%ia(this%nodes + 1))
        this%ia(1) = 1
        allocate(this%iac(this%nodes))
        allocate(this%ja(this%nja))
        allocate(this%ihc(this%nja))
        allocate(this%cl12(this%nja))
        allocate(this%hwva(this%nja))
      end if
    end do !iact
    !
    ! clean up
    if (associated(bba))  deallocate(bba)
    if (allocated(i4wk))  deallocate(i4wk)
    if (allocated(i4wk2)) deallocate(i4wk2)
    if (allocated(x))     deallocate(x)
    !
    return
  end subroutine tMF6Disu_set
   
  subroutine tMF6Disu_write(this, write_opt, write_bin, write_asc, d_out, id_pref)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    integer(I4B), intent(in) :: write_opt
    logical, intent(in), optional :: write_bin
    logical, intent(in), optional :: write_asc
    character(len=*), intent(in), optional :: d_out
    character(len=*), intent(in), optional :: id_pref
    ! -- local
    type(tMf6Wbd), pointer :: wbd => null()
    character(len=MXSLEN) :: d, f_bin, f_asc, id_pref_loc, id_str
    character(len=MXSLEN), dimension(:), allocatable :: id
    logical :: wbin, wasc, wbnp
    integer(I4B) :: i, nid
    !
    integer(I1B), dimension(:), pointer :: i1a
    integer(I2B), dimension(:), pointer :: i2a
    integer(I4B), dimension(:), pointer :: i4a
    real(R8B), dimension(:), pointer :: r8a
! ------------------------------------------------------------------------------
    if (present(id_pref)) then
      id_pref_loc = id_pref
    else
      id_pref_loc = ''
    end if
    !
    wbnp = .true.
    wbin = .false.
    wasc = .false.
    !
    if (present(write_bin)) then
      wbin = write_bin; wbnp = .false.
    end if
    if (present(write_asc)) then
      wasc = write_asc; wbnp = .false.
    end if
    !
    if (wbin.or.wasc) then
      if (.not.present(d_out)) then
        call errmsg('tMF6Disu_write: argument d_out not specified.')
      end if
      d = d_out
    end if
    !
    wbd => this%wbd
    !
    ! write all R8 data
    if (write_opt == 1) then
      allocate(id(5)); id = ['top','bot','area','cl12','hwva']; nid = 5
    else
      nid = 0
    end if
    do i = 1, nid
      select case(id(i))
      case('top')
        r8a => this%top
      case('bot')
        r8a => this%bot
      case('area')
        r8a => this%area
      case('cl12')
        r8a => this%cl12
      case('hwva')
        r8a => this%hwva
      end select
      if (.not.associated(r8a)) cycle
      id_str = trim(id(i))//trim(id_pref)
      if (wbnp) call wbd%write_array(id=id_str, r8a=r8a)
      if (wbin) call wbd%write_array(id=id_str, r8a=r8a, f_bin=trim(d)//trim(id_str)//'.bin')
      if (wasc) call wbd%write_array(id=id_str, r8a=r8a, f_asc=trim(d)//trim(id_str)//'.asc')
    end do
    if (allocated(id)) deallocate(id)
    !
    ! write all I1 data
    if (write_opt == 1) then
      nid = 0
    else
      allocate(id(1)); id = ['map_igrid']; nid = 1
    end if
    do i = 1, nid
      select case(id(i))
      case('map_igrid')
        i1a => this%map_igrid
      end select
      if (.not.associated(i1a)) cycle
      id_str = trim(id(i))//trim(id_pref)
      if (wbnp) call wbd%write_array(id=id_str, i1a=i1a)
      if (wbin) call wbd%write_array(id=id_str, i1a=i1a, f_bin=trim(d)//trim(id_str)//'.bin')
      if (wasc) call wbd%write_array(id=id_str, i1a=i1a, f_asc=trim(d)//trim(id_str)//'.asc')
    end do
    if (allocated(id)) deallocate(id)
    
    ! write all I2 data
    if (write_opt == 1) then
      nid = 0
    else
      allocate(id(1)); id = ['map_ilay']; nid = 1
    end if
    do i = 1, nid
      select case(id(i))
      case('map_ilay')
        i2a => this%map_ilay
      end select
      if (.not.associated(i2a)) cycle
      id_str = trim(id(i))//trim(id_pref)
      if (wbnp) call wbd%write_array(id=id_str, i2a=i2a)
      if (wbin) call wbd%write_array(id=id_str, i2a=i2a, f_bin=trim(d)//trim(id_str)//'.bin')
      if (wasc) call wbd%write_array(id=id_str, i2a=i2a, f_asc=trim(d)//trim(id_str)//'.asc')
    end do
    if (allocated(id)) deallocate(id)
    
    ! write all I4 data
    if (write_opt == 1) then
      allocate(id(4)); id = ['idomain','iac','ja','ihc']; nid = 4
    else
      allocate(id(2)); id = ['ia','map_nod']; nid = 2
    end if
    do i = 1, nid
      select case(id(i))
      case('idomain')
        i4a => this%idomain
      case('ia')
        i4a => this%ia
      case('iac')
        i4a => this%iac
      case('ja')
        i4a => this%ja
      case('ihc')
        i4a => this%ihc
      case('map_nod')
        i4a => this%map_nod
      end select
      if (.not.associated(i4a)) cycle
      id_str = trim(id(i))//trim(id_pref)
      if (wbnp) call wbd%write_array(id=id_str, i4a=i4a)
      if (wbin) call wbd%write_array(id=id_str, i4a=i4a, f_bin=trim(d)//trim(id_str)//'.bin')
      if (wasc) call wbd%write_array(id=id_str, i4a=i4a, f_asc=trim(d)//trim(id_str)//'.asc')
    end do
    if (allocated(id)) deallocate(id)
    !
    ! write the csv-file
    call wbd%write_csv()
    !
    return
  end subroutine tMF6Disu_write
   
  subroutine tMF6Disu_read(this, id_pref)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Disu) :: this
    character(len=*), intent(in), optional :: id_pref
    ! -- local
    type(tMf6Wbd), pointer :: wbd => null()
    character(len=MXSLEN) :: id_pref_loc
    character(len=MXSLEN), dimension(:), allocatable :: id
    integer(I4B) :: i
    !
    integer(I1B), dimension(:), pointer :: i1a
    integer(I2B), dimension(:), pointer :: i2a
    integer(I4B), dimension(:), pointer :: i4a
    real(R8B),    dimension(:), pointer :: r8a
! ------------------------------------------------------------------------------
    !
    if (present(id_pref)) then
      id_pref_loc = id_pref
    else
      id_pref_loc = ''
    end if
    !
    wbd => this%wbd
    !
    allocate(id(5))
    id = ['top','bot','map_igrid','map_ilay','map_nod']
    do i = 1, size(id)
      select case(id(i))
      case('top')
        call wbd%read_array(trim(id(i))//trim(id_pref_loc), r8a=r8a)
        allocate(this%top, source=r8a)
      case('bot')
        call wbd%read_array(trim(id(i))//trim(id_pref_loc), r8a=r8a)
        allocate(this%bot, source=r8a)
      case('map_igrid')
        call wbd%read_array(trim(id(i))//trim(id_pref_loc), i1a=i1a)
        allocate(this%map_igrid, source=i1a)
      case('map_ilay')
        call wbd%read_array(trim(id(i))//trim(id_pref_loc), i2a=i2a)
        allocate(this%map_ilay, source=i2a)
      case('map_nod')
        call wbd%read_array(trim(id(i))//trim(id_pref_loc), i4a=i4a)
        allocate(this%map_nod, source=i4a)
      end select
    end do
    !
    ! clean up
    deallocate(id)
    !
    return
  end subroutine tMF6Disu_read
  
! ==============================================================================
! ==============================================================================
! tProps
! ==============================================================================
! ==============================================================================
 
   subroutine tProps_init(this, f, fields)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tProps) :: this
    character(len=*), intent(in) :: f
    character(len=*), dimension(:), intent(in) :: fields
    !
    ! -- local
    type(tCSV), pointer :: csv => null()
! ------------------------------------------------------------------------------
    !
    allocate(this%csv); csv => this%csv
    call csv%read(f)
    this%fields = fields
    !
    return
  end subroutine tProps_init
  
   subroutine tProps_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tProps) :: this
    !
! ------------------------------------------------------------------------------
    !
    if (associated(this%csv)) then
      call this%csv%clean(); deallocate(this%csv); this%csv => null()
    end if
    !
    return
  end subroutine tProps_clean
  
  subroutine tProps_get_field(this, i, field, empty)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tProps) :: this
    integer(I4B), intent(in) :: i
    character(len=MXSLEN), intent(out) :: field ! result
    logical, optional, intent(out) :: empty
    !
    ! -- local
    logical :: empty_loc
! ------------------------------------------------------------------------------
    empty_loc = .false.
    !
    if (i > size(this%fields)) then
      call errmsg('tProps_get_field: index out of range.')
    end if
    field = this%fields(i)
    if (len_trim(field) == 0) then
      if (present(empty)) then
        empty_loc = .true.
      else
        call errmsg('tProps_get_field: field not found.')
      end if
    end if
    !
    if (present(empty)) then
      empty = empty_loc 
    end if
    !
    return
  end subroutine tProps_get_field
  
! ==============================================================================
! ==============================================================================
! tDataModels
! ==============================================================================
! ==============================================================================
   subroutine tDataModels_init(this, n_inp_mod)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModels) :: this
    integer(I4B), intent(in) :: n_inp_mod
    ! -- local
! ------------------------------------------------------------------------------
    call this%clean()
    this%n_inp_mod = n_inp_mod
    allocate(this%dat_mods(n_inp_mod))
    allocate(this%uni_map_id(NMAX_UNI_MAP))
    allocate(this%uni_map_ir(n_inp_mod,NMAX_UNI_MAP))
    this%n_uni_map = 0
    this%uni_map_id = ''
    this%uni_map_ir = 0
    !
    return
  end subroutine tDataModels_init
  !
   subroutine tDataModels_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModels) :: this
    ! -- local
    type(tdataModel), pointer :: dat_mod
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    if (associated(this%dat_mods)) then
      do i = 1, this%n_inp_mod
        dat_mod => this%dat_mods(i)
        if (associated(dat_mod)) then
          call dat_mod%clean()
        end if 
      end do
      deallocate(this%dat_mods); this%dat_mods => null()
    end if
    !
    this%n_inp_mod = 0
    !
    if (allocated(this%uni_map_id)) deallocate(this%uni_map_id)
    if (allocated(this%uni_map_ir)) deallocate(this%uni_map_ir)
    !
    return
  end subroutine tDataModels_clean
  
  subroutine tDataModels_create_uni_map(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModels) :: this
    ! -- local
    type(tdataModel), pointer :: dat_mod => null()
    type(tCSV), pointer :: param_map => null()
    !
    logical :: found
    character(len=MXSLEN) :: id
    character(len=MXSLEN), dimension(:), allocatable :: id_arr
    integer(I4B) :: i, j, k, p
! ------------------------------------------------------------------------------
    do i = 1, this%n_inp_mod
      dat_mod => this%dat_mods(i)
      param_map => dat_mod%param_map
      !
      call param_map%get_column(key='id', ca=id_arr)
      !
      do j = 1, size(id_arr)
        id = id_arr(j)
        found = .false.
        do k = 1, this%n_uni_map
          if (id == this%uni_map_id(k)) then
            found = .true.; p = k
            exit
          end if
        end do
        !
        if (.not.found) then
          this%n_uni_map = this%n_uni_map + 1
          p = this%n_uni_map
          this%uni_map_id(p) = id
        end if
        this%uni_map_ir(i,p) = j
      end do
    end do
    !
    return
  end subroutine tDataModels_create_uni_map
   
! ==============================================================================
! ==============================================================================
! tDataModel
! ==============================================================================
! ==============================================================================
   
   subroutine tDataModel_init(this, name, f_list, f_map, n_inp_mod)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModel) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: f_list
    character(len=*), intent(in) :: f_map
    integer(I4B), intent(in) :: n_inp_mod
    ! -- local
    type(tCSV), pointer :: csv
    integer(I4B) :: ic
! ------------------------------------------------------------------------------
    call this%clean()
    !
    allocate(this%param_list)
    csv => this%param_list
    call csv%read(f_list)
    !
    allocate(this%param_map)
    csv => this%param_map
    call csv%read(f_map)
!    do ic = 1, n_inp_mod
!      call csv%add_key('inp_mod_'//ta([ic]))
!    end do
    !
    this%name = trim(name)
    !
    return
  end subroutine tDataModel_init
   
   subroutine tDataModel_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModel) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (associated(this%param_list)) then
      call this%param_list%clean()
      deallocate(this%param_list); this%param_list => null()
    end if
    if (associated(this%param_map)) then
      call this%param_map%clean()
      deallocate(this%param_map); this%param_map => null()
    end if
    this%name = ''
    !
    return
  end subroutine tDataModel_clean
   
  function tDataModel_get_ndat(this) result(ndat)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModel) :: this
    integer(I4B) :: ndat ! result
! ------------------------------------------------------------------------------
    !
    ndat = this%param_map%get_nr()
    !
    return
  end function tDataModel_get_ndat
  !
  subroutine tDataModel_get_dat(this, idat, dmdat, init_read)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModel) :: this
    integer(I4B), intent(in) :: idat
    type(tDataModelData), pointer, intent(inout) :: dmdat
    logical, intent(in), optional :: init_read
    ! -- local
    integer(I4B), parameter :: max_ndat = 3
    !
    type(tCSV), pointer :: parl => null(), parm => null(), csv => null()
    type(tData), pointer :: dat => null()
    type(tVrt), pointer :: vrt => null()
    type(tVrtArray), pointer :: vrta => null()
    !
    character(len=MXSLEN), dimension(max_ndat) :: param
    character(len=MXSLEN), dimension(:), allocatable :: ids
    character(len=MXSLEN) :: id, output_file_type, scltype_up, scltype_down, ext
    character(len=MXSLEN) :: s, s1, s2, s3, file, file_type, assign_method
    !
    logical :: write_nod, init_read_loc
    integer(I4B) :: i, j, ir, i_uscl, i_dscl, i_min, i_max
! ------------------------------------------------------------------------------
    if (present(init_read)) then
      init_read_loc = init_read
    else
      init_read_loc = .true.
    end if
    !
    if (associated(dmdat)) then
      call dmdat%clean(); deallocate(dmdat); dmdat => null()
    end if
    allocate(dmdat)
    !
    parl => this%param_list
    parm => this%param_map
    !
    ! read the field of the mapping csv-file
    call parm%get_val(ir=idat, ic=parm%get_col('id'), cv=dmdat%id)
    !
    call parm%get_val(ir=idat, ic=parm%get_col('out_file_type'), cv=s)
    output_file_type = change_case(s, 'l')
    select case(output_file_type)
    case('binpos')
      dmdat%i_out_file_type = i_binpos
    case('bin')
      dmdat%i_out_file_type = i_bin
    case('asc')
      dmdat%i_out_file_type = i_asc
    case default
      call errmsg('Output_type not recognized.')
    end select
    !
    call parm%get_val(ir=idat, ic=parm%get_col('write_nod'), l4v=write_nod)
    if (write_nod) then
      dmdat%i_type = i_type_list
    else
      dmdat%i_type = i_type_array
    end if
    !
    do i = 1, max_ndat
      call parm%get_val(ir=idat, ic=parm%get_col('param_'//ta([i])), cv=s)
      if (len_trim(s) > 0) then
        dmdat%ndat = dmdat%ndat + 1
        param(dmdat%ndat) = s
      end if
    end do
    if (dmdat%ndat == 0) then
      call errmsg('No parameters are set.')
    end if
    !
    call parl%get_column(key='id', ca=ids)
    allocate(dmdat%dat(dmdat%ndat))
    do i = 1, dmdat%ndat
      dat => dmdat%dat(i)
      !
      ir = 0
      do j = 1, size(ids)
        if (trim(ids(j)) == trim(param(i))) then
          ir = j; exit
        end if
      end do
      if (ir == 0) then
        call errmsg('tDataModel_get_dat: id '//trim(param(i))//' not found.')
      end if
      !
      call parl%get_val(ir=ir, ic=parl%get_col('file'),         cv=file)
      ext = change_case(get_ext(file), 'l')
      select case(ext)
      case('.vrt')
        dmdat%i_in_file_type = i_vrt
        dat%i_type = i_vrt
        allocate(dat%vrt); vrt => dat%vrt
        if (init_read_loc) call dat%vrt%init(file)
      case('.txt')
        dmdat%i_in_file_type = i_vrt_array
        dat%i_type = i_vrt_array
        allocate(dat%vrta); vrta => dat%vrta
        vrta%f = file
        if (init_read_loc) call vrta%init()
      case('.csv')
        dmdat%i_in_file_type = i_csv
        dat%i_type = i_csv
        allocate(dat%csv); csv => dat%csv
        if (init_read_loc) call csv%read(file)
      case default
        call errmsg('File not recognized: '//trim(file))
      end select
      !
      ! up- and downscaling
      call parl%get_val(ir=ir, ic=parl%get_col('scltype_up'),   cv=s)
      scltype_up = change_case(s, 'l')
      call parl%get_val(ir=ir, ic=parl%get_col('scltype_down'), cv=s)
      scltype_down = change_case(s, 'l')
      !
      i_uscl = i_uscl_nodata; i_dscl = i_dscl_nodata
      do j = 1, n_uscl
        if (trim(scltype_up) == trim(uscl_names(j))) then
          i_uscl = j; exit
        end if
      end do
      do j = 1, n_dscl
        if (trim(scltype_down) == trim(dscl_names(j))) then
          i_dscl = j; exit
        end if
      end do
      if ((dat%i_type == i_vrt).or.(dat%i_type == i_vrt_array)) then
        if ((i_uscl == 0).or.(i_dscl == 0)) then
          call errmsg('Invalid scaling option found.')
         end if
      else
        i_uscl = i_uscl_nodata; i_dscl = i_dscl_nodata
      end if
      dat%i_uscl = i_uscl; dat%i_dscl = i_dscl
    end do
    !
    call parm%get_val(ir=idat, ic=parm%get_col('ilay'), cv=s)
    call cast_number_from_string(cv=s, i4v=dmdat%ilay, i4mv=0)
    if (dmdat%ilay > 0) then
      dmdat%i_assign = i_assign_exact_layer
    end if
    if (dmdat%ilay == 0) then
      dmdat%i_assign = i_assign_intersect
    end if
    if (dmdat%ilay < 0) then
      dmdat%i_assign = i_assign_first_layer
    end if
    !
    dmdat%lconst = .false.
    if ((dmdat%i_type == i_type_array).and.&
        (dmdat%i_assign == i_assign_first_layer)) then
      if (.not.parm%exist_col('fill_value')) then
        call errmsg('Field fill_value not found.')
      end if
      call parm%get_val(ir=idat, ic=parm%get_col('fill_value'), cv=s)
      if (len_trim(s) == 0) then
        call errmsg('No value for field fill_value  found.')
      end if
      dmdat%lconst = .true.
      read(s,*) dmdat%r8const
    end if
    !
    dmdat%ilay0_sfac = .false.
    if ((dmdat%i_in_file_type == i_csv).and.(dmdat%ilay == 0).and. &
        (parm%exist_col('ilay0_sfac'))) then
      call parm%get_val(ir=idat, ic=parm%get_col('ilay0_sfac'), cv=s)
      if (len_trim(s) > 0) then
        dmdat%ilay0_sfac = .true.
        ! find the ID
        ir = 0
        do j = 1, size(ids)
          if (trim(ids(j)) == trim(s)) then
            ir = j; exit
          end if
        end do
        if (ir == 0) then
          call errmsg('tDataModel_get_dat: id '//trim(s)//' not found.')
        end if
        call parl%get_val(ir=ir, ic=parl%get_col('file'),         cv=file)
        ext = change_case(get_ext(file), 'l')
        if (ext == '.txt') then
          dat => dmdat%dat(1)
          ! up- and downscaling
          call parl%get_val(ir=ir, ic=parl%get_col('scltype_up'),   cv=s)
          scltype_up = change_case(s, 'l')
          call parl%get_val(ir=ir, ic=parl%get_col('scltype_down'), cv=s)
          scltype_down = change_case(s, 'l')
          i_uscl = i_uscl_nodata; i_dscl = i_dscl_nodata
          do j = 1, n_uscl
            if (trim(scltype_up) == trim(uscl_names(j))) then
              i_uscl = j; exit
            end if
          end do
          do j = 1, n_dscl
            if (trim(scltype_down) == trim(dscl_names(j))) then
              i_dscl = j; exit
            end if
          end do
          if ((i_uscl == 0).or.(i_dscl == 0)) then
            call errmsg('Invalid scaling option found.')
          end if
          dat%i_uscl = i_uscl; dat%i_dscl = i_dscl
          dat%i_type = i_vrt_array
          allocate(dat%vrta); vrta => dat%vrta
          vrta%f = file
        else
          call errmsg('tDataModel_get_dat: only VRT files allowed for ilay0_sfac.')
        end if
      end if
    end if
    !
    s1 = ''; s2 = ''; s3 = ''
    if (parm%exist_col('itop')) then
      call parm%get_val(ir=idat, ic=parm%get_col('itop'), cv=s1)
    end if
    if (parm%exist_col('ibot')) then
      call parm%get_val(ir=idat, ic=parm%get_col('ibot'), cv=s2)
    end if
    if (parm%exist_col('idist')) then
      call parm%get_val(ir=idat, ic=parm%get_col('idist'), cv=s3)
    end if
    if (dmdat%i_in_file_type == i_csv) then
      call cast_number_from_string(cv=s1, i4v=dmdat%itop, i4mv=4)
      call cast_number_from_string(cv=s2, i4v=dmdat%ibot, i4mv=5)
      call cast_number_from_string(cv=s3, i4v=dmdat%idist, i4mv=3)
    else
      select case(dmdat%ndat)
      case(1)
        call cast_number_from_string(cv=s1, i4v=dmdat%itop,  i4mv=1)
        call cast_number_from_string(cv=s2, i4v=dmdat%ibot,  i4mv=1)
        call cast_number_from_string(cv=s3, i4v=dmdat%idist, i4mv=1)
      case(2)
        call cast_number_from_string(cv=s1, i4v=dmdat%itop,  i4mv=1)
        call cast_number_from_string(cv=s2, i4v=dmdat%ibot,  i4mv=1)
        call cast_number_from_string(cv=s3, i4v=dmdat%idist, i4mv=2)
      case(3)
        call cast_number_from_string(cv=s1, i4v=dmdat%itop,  i4mv=1)
        call cast_number_from_string(cv=s2, i4v=dmdat%ibot,  i4mv=3)
        call cast_number_from_string(cv=s3, i4v=dmdat%idist, i4mv=2)
      end select
    end if
    !
    ! checks: TODO
    i_min = huge(I4ZERO); i_max = -huge(I4ZERO)
    do i = 1, dmdat%ndat
      dat => dmdat%dat(i)
      i_min = min(i_min, dat%i_type)
      i_max = max(i_max, dat%i_type)
    end do
    if (i_min /= i_max) then
      call errmsg('Mixed file types are not allowed.')
    end if
    !
    ! clean up
    if (allocated(ids)) deallocate(ids)
    !
    return
  end subroutine tDataModel_get_dat
  
! ==============================================================================
! ==============================================================================
! tDataModels
! ==============================================================================
! ==============================================================================
   subroutine tDataModelData_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModelData) :: this
! ------------------------------------------------------------------------------

    
    return
  end subroutine tDataModelData_init
   
   subroutine tDataModelData_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModelData) :: this
! ------------------------------------------------------------------------------

    
    return
  end subroutine tDataModelData_clean
   
! ==============================================================================
! ==============================================================================
! tLayerModels
! ==============================================================================
! ==============================================================================
  
   subroutine tLayerModels_init(this, f_csv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLayerModels) :: this
    character(len=*), intent(in) :: f_csv
    !
    ! -- local
    type(tCSV) :: csv
    type(tLookupTable), pointer :: lookup => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    call csv%read(f_csv)
    this%n_inp_mod = csv%get_nc()
    !
    if (this%n_inp_mod <= 0) then
      call errmsg('tLayerModels_init: invalid number of layer models.')
    end if
    !
    allocate(this%lookup)
    lookup => this%lookup
    call csv%get_matrix(i4x=lookup%table)
    !
    allocate(lookup%keys(this%n_inp_mod))
    do i = 1, this%n_inp_mod
      lookup%keys(i) = csv%get_key(i)
    end do
    !
    allocate(this%lay_mods(this%n_inp_mod))
    !
    return
  end subroutine tLayerModels_init
 
   subroutine tLayerModels_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLayerModels) :: this
    !
    ! -- local
    type(tLayerModel), pointer :: lay_mod => null()
    integer(I4B) :: il
! ------------------------------------------------------------------------------
    !
    if (associated(this%lookup)) then
      call this%lookup%clean(); deallocate(this%lookup)
    end if
    if (associated(this%lay_mods)) then
      do il = 1, this%n_inp_mod
        lay_mod => this%lay_mods(il)
        call lay_mod%clean()
      end do
      deallocate(this%lay_mods)
    end if
    !
    this%n_inp_mod  = 0
    this%lookup     => null()
    this%lay_mods   => null()
    !
    return
  end subroutine tLayerModels_clean
  
  subroutine tLayerModels_get(this, n_inp_mod, lookup_keys, nlay)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLayerModels) :: this
    !
    integer(I4B), intent(out), optional :: n_inp_mod
    character(len=MXSLEN), dimension(:), allocatable, &
      intent(inout), optional :: lookup_keys
    integer(I4B), dimension(:), allocatable, intent(inout), optional :: nlay
    ! -- local
    type(tLookupTable), pointer :: lookup => null()
    integer(I4B), dimension(:,:), pointer :: x => null()
    integer(I4B) :: ic, nc, nr
! ------------------------------------------------------------------------------
    if (present(n_inp_mod)) then
      n_inp_mod = this%n_inp_mod
    end if
    !
    if (present(lookup_keys)) then
      if (.not.associated(this%lookup)) then
        call errmsg('tLayerModels_get: lookup not set.')
      end if
      lookup => this%lookup
      if (allocated(lookup_keys)) deallocate(lookup_keys)
      allocate(lookup_keys, source=lookup%keys)
    end if
    !
    if (present(nlay)) then
      if (.not.associated(this%lookup)) then
        call errmsg('tLayerModels_get: lookup not set.')
      end if
      lookup => this%lookup
      if (allocated(nlay)) deallocate(nlay)
      x => lookup%table; nc = size(x,1); nr = size(x,2)
      allocate(nlay(nc))
      do ic = 1, nc
        nlay(ic) = maxval(x(ic,:))
      end do
    end if
    !
    return
  end subroutine tLayerModels_get
   
! ==============================================================================
! ==============================================================================
! tLayerModels
! ==============================================================================
! ==============================================================================
  
   subroutine tLayerModel_init(this, name, nlay, f_vrt, compressed)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLayerModel) :: this
    character(len=*), intent(in) :: name
    integer(I4B), intent(in) :: nlay
    character(len=*), dimension(:), intent(in) :: f_vrt
    logical, intent(in), optional :: compressed
    !
    ! -- local
    type(tVrt), pointer :: vrt => null()
    integer(I4B) :: il
! ------------------------------------------------------------------------------
    !
    this%compressed = .false.
    if (present(compressed)) then
      this%compressed = compressed
    end if
    !
    this%name = name
    this%nlay = nlay
    !
    if (size(f_vrt) /= (nlay + 1)) then
      call errmsg('tLayerModel_init: invalid data.')
    end if
    !
    if (this%compressed) then
      allocate(this%comp_dat_ptr)
      vrt => this%comp_dat_ptr
      call vrt%init(f_vrt(1))
    else
      allocate(this%top)
      vrt => this%top
      allocate(this%bots(nlay))
      call vrt%init(f_vrt(1))
      do il = 1, nlay
        vrt => this%bots(il)
        call vrt%init(f_vrt(il+1))
      end do
    end if
    !
    return
  end subroutine tLayerModel_init
 
   subroutine tLayerModel_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLayerModel) :: this
    !
    ! -- local
    type(tVrt), pointer :: vrt => null()
    integer(I4B) :: il
! ------------------------------------------------------------------------------
    !
    if (associated(this%top)) then
      call this%top%clean()
      deallocate(this%top)
    end if
    !
    if (associated(this%bots)) then
      do il = 1, this%nlay
        vrt => this%bots(il)
        call vrt%clean()
      end do
      deallocate(this%bots)
    end if
    !
    if (associated(this%comp_dat_ptr)) then
      call this%comp_dat_ptr%clean()
      deallocate(this%comp_dat_ptr)
    end if
    !
    this%compressed = .false.
    this%name = ''
    this%nlay = 0
    this%top => null()
    this%bots => null()
    !
    return
  end subroutine tLayerModel_clean
   
  subroutine tLayerModel_read_full_grid(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLayerModel) :: this
    !
    ! -- local
    type(tVrt), pointer :: vrt => null()
    integer(I4B) :: il, it
! ------------------------------------------------------------------------------
    !
    do il = 1, this%nlay + 1
      if (il == 1) then
        vrt => this%top
      else
        vrt => this%bots(il-1)
      end if
      do it = 1, vrt%ntiles
        call vrt%read_full_tile(itile=it, clean_hdrg=.false.)
      end do
      vrt%full_data_read = .true.
    end do
    !
    return
  end subroutine tLayerModel_read_full_grid
   
  subroutine tLayerModel_read_xy(this, x, y, cs, zp, mv, lay_act, nl_act)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLayerModel) :: this
    !
    real(R8B), dimension(:), intent(in) :: x
    real(R8B), dimension(:), intent(in) :: y
    real(R8B), dimension(:), intent(in) :: cs
    !
    real(R4B), dimension(:,:), allocatable, intent(inout) :: zp
    real(R4B), intent(out) :: mv
    integer(I4B), dimension(:), allocatable, intent(inout) :: lay_act
    integer(I4B), intent(out) :: nl_act
    !
    ! -- local
    type(tVrt), pointer :: vrt => null()
    character(len=MXSLEN) :: f
    character(len=MXSLEN), dimension(:), allocatable :: f_bin
    logical :: lop
    real(R4B), dimension(:,:), allocatable :: zp_all
    real(R4B), dimension(:), allocatable :: r4wk, z_read
    real(R4B) :: mvr4_read, top, bot, top_next, thk
    real(R8B), dimension(:), allocatable :: r8wk
    real(R8B) :: mvr8_read
    integer(I4B), dimension(:), allocatable :: iu, lay_read, xtile, flag
    integer(I4B) ::  np, i_uscl, i_dscl, il, jl, i, ip, ntile, it, nl
    integer(I8B) :: p
! ------------------------------------------------------------------------------
    ! up- and downscaling parameters
    i_uscl = i_uscl_arith
    i_dscl = i_dscl_nointp
    !
    np = size(x)
    if ((np /= size(y)).or.(np /= size(cs))) then
      call errmsg('tLayerModel_read_xy: invalid dimensions.')
    end if
    !
    if (allocated(zp)) deallocate(zp)
    !
    if (this%compressed) then
      vrt => this%comp_dat_ptr
      !
      call vrt%read_xy(x=x, y=y, cs=cs, r8a=r8wk, mvr8=mvr8_read, &
        i_uscl=i_uscl_nodata, i_dscl=i_dscl_nointp, xtile=xtile, f_bin=f_bin)
      !
      ntile = size(f_bin); allocate(iu(ntile))
      do it = 1, ntile
        f = f_bin(it)
        if (len_trim(f) > 0) then
          call open_file(f, iu(it), 'R', .true.)
        end if
      end do
      !
      ! first, determine the number of layers
      if (allocated(lay_act)) deallocate(lay_act)
      allocate(lay_act(this%nlay)); lay_act = 0
      allocate(lay_read(this%nlay), z_read(this%nlay))
      !
      do ip = 1, np
        p = int(r8wk(ip),I8B)
        if (p > 0) then
          it = xtile(ip)
          if (it <= 0) then
            call errmsg('tLayerModel_read_xy: program error.')
          end if
          read(iu(it), pos=p) nl; p = p + I4B
          lay_read = 0
          read(iu(it), pos=p) lay_read(1:nl)
          do i = 1, nl
            il = lay_read(i)
            if ((il < 0).or.(il > this%nlay)) then
              call errmsg('tLayerModel_read_xy: program error.')
            end if
            lay_act(il) = 1
          end do
        end if
      end do
      !
      nl_act = sum(lay_act)
      nl_act = 0
      do i = 1, this%nlay
        if (lay_act(i) == 1) then
          nl_act = nl_act + 1
          lay_act(i) = nl_act
        end if
      end do
      !
      if (nl_act == 0) then
        call errmsg('tLayerModel_read_xy: no layers found.')
      end if
      !
      mv = -99999.
      allocate(zp(np,nl_act+1)); zp = mv
      !
      do ip = 1, np
        p = int(r8wk(ip),I8B)
        if (p > 0) then
          it = xtile(ip)
          read(iu(it), pos=p) nl; p = p + I4B
          lay_read = 0; z_read = mv
          read(iu(it), pos=p) lay_read(1:nl); p = p + nl*I4B
          read(iu(it), pos=p) z_read(1:nl+1)
          do i = 1, nl
            il = lay_read(i); jl = lay_act(il)
            top = z_read(i); bot = z_read(i+1)
            zp(ip,jl)   = top
            zp(ip,jl+1) = bot
          end do
        end if
      end do
      
      ! remove missing values
      do ip = 1, np
        do il = 1, nl_act
          top = zp(ip,il)
          bot = zp(ip,il+1)
          if ((top == mv).and.(bot == mv)) then
            top_next = mv
            do jl = il, nl_act
              top_next = zp(ip,jl)
              if (top_next /= mv) exit
            end do
            if (top_next == mv) then
              !call errmsg('No layer model found.')
            else
              zp(ip,il) = top_next
              zp(ip,il+1) = top_next
            end if
          end if
          if ((top == mv).and.(bot /= mv)) then
            zp(ip,il) = bot
          end if
          if ((bot == mv).and.(top /= mv)) then
            zp(ip,il+1) = top
          end if
        end do
      end do
      !
      ! clean up
      if (allocated(iu)) then
        do it = 1, ntile
          inquire(unit=iu(it), opened=lop)
          if (lop) then
            close(iu(it))
          end if
        end do
        deallocate(iu)
      end if
      if (allocated(lay_read)) deallocate(lay_read)
      if (allocated(z_read)) deallocate(z_read)
    else
      allocate(zp_all(np,this%nlay+1))
      if (allocated(lay_act)) deallocate(lay_act)
      allocate(lay_act(this%nlay)); lay_act = 0
      !
      do il = 1, this%nlay + 1
        if (il == 1) then
          vrt => this%top
        else
          vrt => this%bots(il-1)
        end if
        !
        call vrt%read_xy(x=x, y=y, cs=cs, r4a=r4wk, mvr4=mvr4_read, &
          i_uscl=i_uscl_arith, i_dscl=i_dscl_nointp)
        !
        zp_all(:,il) = r4wk
        !
        if (il == 1) then
          mv = mvr4_read
        else
          ! replace missing value
          do i = 1, np
            if (zp_all(i,il) == mvr4_read) then
              zp_all(i,il) = mv
            end if
          end do
        end if
        !
        if (il > 1) then
          do i = 1, np
            top = zp_all(i,il-1); bot = zp_all(i,il)
            if ((top /= mv).and.(bot /= mv)) then
              thk = top - bot
              if (abs(thk) > R4TINY) then
                lay_act(il-1) = 1
                exit
              end if
            end if
          end do
        end if
      end do
      !
      nl_act = 0
      do il = 1, this%nlay
        if (lay_act(il) > 0) then
          nl_act = nl_act + 1
          lay_act(il) = nl_act
        end if
      end do
      if (nl_act == 0) then
        call errmsg('tLayerModel_read_xy: no layers found.')
      end if
      !
      ! fill for active layers only
      allocate(zp(np,nl_act+1)); zp = mv
      do il = 1, this%nlay
        jl = lay_act(il)
        if (jl > 0) then
          zp(:,jl)   = zp_all(:,il)
          zp(:,jl+1) = zp_all(:,il+1)
        end if
      end do
      !
      deallocate(zp_all)
    end if
    !
    return
  end subroutine tLayerModel_read_xy
    
  subroutine tLayerModel_read_extent(this, bbx, zp, mv, lay_act, nl_act)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLayerModel) :: this
    type(tBBX), intent(in) :: bbx
    real(R4B), dimension(:,:,:), allocatable, intent(inout) :: zp
    real(R4B), intent(out) :: mv
    integer(I4B), dimension(:), allocatable, intent(inout) :: lay_act
    integer(I4B), intent(out) :: nl_act
    !
    ! -- local
    type(tVrt), pointer :: vrt => null()
    character(len=MXSLEN) :: f
    character(len=MXSLEN), dimension(:), allocatable :: f_bin
    logical :: lop, lfound
    real(R4B) :: mvr4_read, thk, top, bot, top_next
    real(R4B), dimension(:), allocatable :: z_read
    real(R4B), dimension(:,:), allocatable :: r4wk
    real(R4B), dimension(:,:,:), allocatable :: zp_read
    real(R8B) :: mvr8_read
    real(R8B), dimension(:,:), allocatable :: r8wk
    integer(I4B) :: il, jl, i_uscl, i_dscl, nc, nr, mc, mr, ir, ic
    integer(I4B) :: ntile, it, nl, i
    integer(I4B), dimension(:), allocatable :: iu, lay_read
    integer(I4B), dimension(:,:), allocatable :: xtile
    integer(I8B) :: p
! ------------------------------------------------------------------------------
    ! up- and downscaling parameters
    i_uscl = i_uscl_arith
    i_dscl = i_dscl_nointp
    !
    nc = (bbx%xur - bbx%xll)/bbx%cs
    nr = (bbx%yur - bbx%yll)/bbx%cs
    !
    if (allocated(zp)) deallocate(zp)
    if (allocated(lay_act)) deallocate(lay_act)
    !
    if (this%compressed) then
      vrt => this%comp_dat_ptr
      call vrt%read_extent(xr8=r8wk, mvr8=mvr8_read, bbx=bbx, &
        i_uscl=i_uscl_nodata, i_dscl=i_dscl_nointp, xtile=xtile, f_bin=f_bin)
      ntile = size(f_bin); allocate(iu(ntile))
      do it = 1, ntile
        f = f_bin(it)
        if (len_trim(f) > 0) then
          call open_file(f, iu(it), 'r', .true.)
        end if
      end do
      !
      mc = size(r8wk,1); mr = size(r8wk,2)
      if ((mc /= nc).or.(mr /= nr)) then
        call errmsg('tLayerModel_read: inconsistent dimensions.')
      end if
      
      ! first, determine the number of layers
      allocate(lay_act(this%nlay)); lay_act = 0
      allocate(lay_read(this%nlay), z_read(this%nlay))
      !
      do ir = 1, nr; do ic = 1, nc
        p = int(r8wk(ic,ir),I8B)
        if (p > 0) then
          it = xtile(ic,ir)
          read(iu(it), pos=p) nl; p = p + I4B
          lay_read = 0
          read(iu(it), pos=p) lay_read(1:nl)
          do i = 1, nl
            il = lay_read(i)
            if ((il < 0).or.(il > this%nlay)) then
              call errmsg('tLayerModel_read_extent: program error.')
            end if
            lay_act(il) = 1
          end do
        end if
      end do; end do
      !
      ! count the number of layers; create the mapping
      nl_act = sum(lay_act)
      nl_act = 0
      do i = 1, this%nlay
        if (lay_act(i) == 1) then
          nl_act = nl_act + 1
          lay_act(i) = nl_act
        end if
      end do
      !
      if (nl_act == 0) then
        call errmsg('tLayerModel_read_extent: no layers found.')
      end if
      !
      mv = -99999.
      allocate(zp(nc,nr,nl_act+1)); zp = mv
      !
      do ir = 1, nr; do ic = 1, nc
        p = int(r8wk(ic,ir),I8B)
        if (p > 0) then
          it = xtile(ic,ir)
          read(iu(it), pos=p) nl; p = p + I4B
          lay_read = 0; z_read = mv
          read(iu(it), pos=p) lay_read(1:nl); p = p + nl*I4B
          read(iu(it), pos=p) z_read(1:nl+1)
          do i = 1, nl
            il = lay_read(i); jl = lay_act(il)
            top = z_read(i); bot = z_read(i+1)
            zp(ic,ir,jl)   = top
            zp(ic,ir,jl+1) = bot
          end do
        end if
      end do; end do
      !
      ! remove missing values
      do ir = 1, nr; do ic = 1, nc
        do il = 1, nl_act
          top = zp(ic,ir,il)
          bot = zp(ic,ir,il+1)
          if ((top == mv).and.(bot == mv)) then
            top_next = mv
            do jl = il, nl_act
              top_next = zp(ic,ir,jl)
              if (top_next /= mv) exit
            end do
            if (top_next == mv) then
              !call errmsg('No layer model found.')
            else
              zp(ic,ir,il) = top_next
              zp(ic,ir,il+1) = top_next
            end if
          end if
          if ((top == mv).and.(bot /= mv)) then
            zp(ic,ir,il) = bot
          end if
          if ((bot == mv).and.(top /= mv)) then
            zp(ic,ir,il+1) = top
          end if
        end do
      end do; end do
      !
      ! clean up
      if (allocated(iu)) then
        do it = 1, ntile
          inquire(unit=iu(it), opened=lop)
          if (lop) then
            close(iu(it))
          end if
        end do
        deallocate(iu)
      end if
      if (allocated(lay_read)) deallocate(lay_read)
      if (allocated(z_read)) deallocate(z_read)
    else
      allocate(zp_read(nc,nr,this%nlay+1))
      !
      ! read for all layers first
      do il = 1, this%nlay + 1
        if (il == 1) then
          vrt => this%top
        else
          vrt => this%bots(il-1)
        end if
        !
        if (allocated(r4wk)) deallocate(r4wk)
        !
        call vrt%read_extent(xr4=r4wk, mvr4=mvr4_read, bbx=bbx, &
          i_uscl=i_uscl_arith, i_dscl=i_dscl_nointp)
        call vrt%clean_x()
        !
        mc = size(r4wk,1); mr = size(r4wk,2)
        if ((mc /= nc).or.(mr /= nr)) then
          call errmsg('tLayerModel_read: inconsistent dimensions.')
        end if
        !
        if (il == 1) then
          mv = mvr4_read
        else
          ! replace missing value
          do ir = 1, nr; do ic = 1, nc
            if (r4wk(ic,ir) == mvr4_read) then
              r4wk(ic,ir) = mv
            end if
          end do; end do
        end if
        !
        zp_read(:,:,il) = r4wk
      end do
      !
      ! determine the active number of layers
      allocate(lay_act(this%nlay)); lay_act = 0
      nl_act = 0
      do il = 1, this%nlay
        lfound = .false.
        do ir = 1, nr
          do ic = 1, nc
            top = zp_read(ic,ir,il); bot = zp_read(ic,ir,il+1)
            if ((top /= mv).and.(bot /= mv)) then
              thk = top - bot
              if (abs(thk) > R4TINY) then
                nl_act = nl_act + 1
                lay_act(il) = nl_act
                lfound = .true.
              end if
            end if
            if (lfound) exit
          end do
          if (lfound) exit
        end do
      end do
      !
      if (nl_act == 0) then
        call errmsg('tLayerModel_read_extent: no layers found.')
      end if
      !
      allocate(zp(nc,nr,nl_act+1))
      do il = 1, this%nlay
        jl = lay_act(il)
        if (jl > 0) then
          zp(:,:,jl)   = zp_read(:,:,il) !top
          zp(:,:,jl+1) = zp_read(:,:,il+1) !bot
        end if
      end do
      !
      ! clean up
      if (allocated(zp_read)) deallocate(zp_read)
      if (allocated(r4wk)) deallocate(r4wk)
    end if
    !
    ! check
    do ir = 1, nr; do ic = 1, nc
      top = zp(ic,ir,1)
      do il = 2, nl_act + 1
        bot = zp(ic,ir,il)
        if ((top /= mv).and.(bot /= mv)) then
          if (top < bot) then
            call errmsg('tLayerModel_read: inconsistent layer model')
          end if
        end if
        top = bot
      end do
    end do; end do
    !
    return
  end subroutine tLayerModel_read_extent
    
! ==============================================================================
! ==============================================================================
! tLookupTable
! ==============================================================================
! ==============================================================================
   
   subroutine tLookupTable_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLookupTable) :: this
    !
    ! -- local
! ------------------------------------------------------------------------------
    return
  end subroutine tLookupTable_init
 
   subroutine tLookupTable_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tLookupTable) :: this
    !
    ! -- local
! ------------------------------------------------------------------------------
    if (allocated(this%keys)) then
      deallocate(this%keys)
    end if
    if (allocated(this%table)) then
      deallocate(this%table)
    end if
    !
    return
  end subroutine tLookupTable_clean
  
! ==============================================================================
! ==============================================================================
! tGraph
! ==============================================================================
! ==============================================================================
  
  subroutine tGraph_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGraph) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (allocated(this%xadj))     deallocate(this%xadj)
    if (allocated(this%adjncy))   deallocate(this%adjncy)
    if (allocated(this%vwgt))     deallocate(this%vwgt)
    if (allocated(this%adjwgt))   deallocate(this%adjwgt)
    if (allocated(this%g2lnod))   deallocate(this%g2lnod)
    if (allocated(this%l2gnod))   deallocate(this%l2gnod)
    !
    return
  end subroutine tGraph_clean
  
  subroutine tGraph_setval_g2lnod(this, gnod, v)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGraph) :: this
    integer(I4B), intent(in) :: gnod
    integer(I8B), intent(in) :: v
    ! -- local
    integer(I8B) :: n
! ------------------------------------------------------------------------------
    n = int(gnod,I8B) - this%gnod_min + 1
    if (n <= 0) then
      call errmsg('tGraph_set_mapped_nod')
    else
      this%g2lnod(n) = v
    end if
    !
    return
  end subroutine tGraph_setval_g2lnod
  
  function tGraph_getval_g2lnod(this, gnod) result(v)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGraph) :: this
    integer(I4B), intent(in) :: gnod
    integer(I4B) :: v
    ! -- local
    integer(I8B) :: n
! ------------------------------------------------------------------------------
    n = int(gnod,I8B) - this%gnod_min + 1
    if (n <= 0) then
      call errmsg('tGraph_get_mapped_nod')
    else
      v = int(this%g2lnod(n),I4B)
    end if
    !
    return
  end function tGraph_getval_g2lnod

  function tGraph_get_nbytes(this) result(nbytes)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGraph) :: this
    integer(I8B) :: nbytes
! ------------------------------------------------------------------------------
    nbytes = sizeof(this%nvtxs) + sizeof(this%nja) + &
             sizeof(this%xadj) + sizeof(this%adjncy) + sizeof(this%l2gnod)
    !
    return
  end function tGraph_get_nbytes
    
  subroutine tGraph_write(this, iu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGraph) :: this
    integer(I4B), intent(in) :: iu
! ------------------------------------------------------------------------------
    write(iu) this%nvtxs, this%nja
    write(iu) this%xadj
    if (this%nja > 0) then
      write(iu) this%adjncy
    end if
    write(iu) this%l2gnod
    !
    return
  end subroutine tGraph_write

  subroutine tGraph_read(this, iu, pos)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGraph) :: this
    integer(I4B), intent(in) :: iu
    integer(I8B), intent(in), optional :: pos
    ! -- local
    integer(I8B) :: p
! ------------------------------------------------------------------------------
    call this%clean()
    !
    if (.not.present(pos)) then
      read(iu) this%nvtxs, this%nja
    else
      p = pos
      read(iu,pos=p) this%nvtxs, this%nja
      p = p + sizeof(this%nvtxs) + sizeof(this%nja)
    end if
    !
    if (this%nvtxs <= 0) then
      call errmsg('tGraph_read')
    end if
    allocate(this%xadj(this%nvtxs+1))
    if (this%nja > 0) then
      allocate(this%adjncy(this%nja))
    end if
    allocate(this%l2gnod(this%nvtxs))
    !
    if (.not.present(pos)) then
      read(iu) this%xadj
      if (this%nja > 0) then
        read(iu) this%adjncy
      end if
      read(iu) this%l2gnod
    else
      read(iu,pos=p) this%xadj;   p = p + sizeof(this%xadj)
      read(iu,pos=p) this%adjncy; p = p + sizeof(this%adjncy)
      read(iu,pos=p) this%l2gnod; p = p + sizeof(this%l2gnod)
    end if
    !
    return
  end subroutine tGraph_read
  
  recursive subroutine tGraph_disconnect(this, i_start, i_disc, xd)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGraph) :: this
    integer(I4B), intent(in) :: i_start
    integer(I4B), intent(in) :: i_disc
    integer(I4B), dimension(:), intent(inout) :: xd
    ! -- local
    integer(I4B) :: i, j, i_nbr
! ------------------------------------------------------------------------------
    !
    i = i_start
    xd(i) = i_disc
    do j = this%xadj(i)+1, this%xadj(i+1)
      i_nbr = this%adjncy(j)
      if (xd(i_nbr) == 0) then
        call this%disconnect(i_nbr, i_disc, xd)
      end if
    end do
    !
    return
  end subroutine tGraph_disconnect
  
  recursive subroutine tGraph_balance(this, prent, lev_prent, lev_max, n)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tGraph) :: this
    integer(I8B),               intent(in)    :: prent
    integer(I8B),               intent(in)    :: lev_prent
    integer(I8B),               intent(in)    :: lev_max
    integer(I4B),               intent(inout) :: n
    ! -- local
    integer(I4B) :: j
    integer(I8B) :: child, lev_child, tlev, dlev
! ------------------------------------------------------------------------------
    !
    if (lev_prent > lev_max) then
      return
    end if
    !
    this%vwgt(prent) = lev_prent
    !
    do j = this%xadj(prent)+1, this%xadj(prent+1)
      child = this%adjncy(j)
      lev_child = this%vwgt(child)
      dlev = lev_prent - lev_child
      if (dlev > 1) then
        tlev = lev_child + 1
        call this%balance(child, tlev, lev_max, n)
      end if
      if (abs(dlev) > 1) then
        n = n + 1
      end if
    end do
    !
    return
  end subroutine tGraph_balance
  
! ==============================================================================
! ==============================================================================
! tQuads
! ==============================================================================
! ==============================================================================
  !
  subroutine tQuads_init(this, n, n_max, uuid, dir)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    integer(I4B), intent(in), optional     :: n
    integer(I4B), intent(in), optional     :: n_max
    character(len=*), intent(in), optional :: uuid
    character(len=*), intent(in), optional :: dir
    ! -- local
! ------------------------------------------------------------------------------
    !
    this%dir  = ''
    this%uuid = ''
    this%props => null()
    this%n     = 0
    this%n_act = 0
    this%xq    => null()
    this%qg    => null()
    this%qgd   => null()
    !
    ! overwrite with the optional arguments
    if (present(n)) then
      this%n = n
    end if
    !
    if (present(n_max)) then
      allocate(this%xq(n_max))
    end if
    !
    if (present(uuid)) then
      this%uuid = uuid
    end if
    !
    if (present(dir)) then
      this%dir = dir
      call swap_slash(this%dir)
    end if
    !
    return
  end subroutine tQuads_init
  !
!  subroutine tQuads_init_select(this, n_max, f_prop, uuid, dir, fields, &
!    fp_intf, read_vintf)
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(tQuads) :: this
!    integer(I4B), intent(in)     :: n_max
!    character(len=*), intent(in) :: f_prop
!    character(len=*), intent(in) :: uuid
!    character(len=*), intent(in) :: dir
!    character(len=*), dimension(:), intent(in) :: fields
!    character(len=*), intent(in) :: fp_intf
!    logical, intent(in), optional :: read_vintf
!    ! -- local
!    type(tProps), pointer :: props => null()
!    type(tQuad), pointer :: q => null(), q_nbr => null()
!    type(tBbObj) :: bbo
!    type(tIntf), pointer :: intf_read => null()
!    type(tIntf), pointer :: intf => null()
!    type(tNbrIntf), pointer :: nintf => null()
!    type(tGraph), pointer :: g_read => null(), g => null()
!    !
!    character(len=MXSLEN) :: key
!    !
!    integer(I4B), dimension(:), allocatable :: lid_read, gid_read, i_graph_read
!    integer(I4B), dimension(:), allocatable :: laymod_read
!    integer(I4B), dimension(:), allocatable :: i4wk, g2lnod, lev
!    integer(I4B) :: ic, lid, lid_nbr, lid_max, gid, i, j, k, nlid_read, n_disc, i_disc
!    integer(I4B) :: iu_bb, iu_intf, iu_graph, iact, xadj, inbr, n
!! ------------------------------------------------------------------------------
!    call this%init(n_max=n_max, uuid=uuid, dir=dir)
!    call this%read_props(f_prop, fields)
!    props => this%get_props_ptr()
!    nlid_read = props%csv%get_nr()
!    !
!    call props%get_field(i_lid, key)
!    call props%csv%get_column(key=key,i4a=lid_read)
!    !
!    call props%get_field(i_gid, key)
!    call props%csv%get_column(key=key,i4a=gid_read)
!    
!    call props%get_field(i_lay_mod, key)
!    call props%csv%get_column(key=key,i4a=laymod_read)
!    !
!    this%n = maxval(lid_read)
!    !
!    ! set the active flag
!    do lid = 1, this%n
!      q => this%get_quad(lid)
!      call q%init()
!      call q%set_flag(active=.false.)
!    end do
!    do i = 1, nlid_read
!      lid = lid_read(i); q => this%get_quad(lid)
!      call q%set_flag(active=.true.)
!      this%n_act = this%n_act + 1
!    end do
!    !
!    iu_bb   = this%open_file('.bb.bin',  'r',.true.)
!    iu_intf = this%open_file(trim(fp_intf),'r',.true.)
!    !
!    do i = 1, nlid_read
!      lid = lid_read(i); gid = gid_read(i)
!      q => this%get_quad(lid)
!      !
!      ! read the bounding box
!      bbo = this%read_bb(iu_bb, lid)
!      !
!      call q%init(gid=gid, lid=lid, bbo=bbo); q%props => props
!      call q%set(i_prop=i)
!      !
!      ! read the interface
!      allocate(intf_read)
!      call this%read_intf(iu_intf, lid, intf_read, read_vintf)
!      !
!      ! check for inactive interfaces
!      n = 0
!      do inbr = 1, intf_read%n_nbr
!        nintf => intf_read%nbr_intf(inbr)
!        lid_nbr = nintf%nbr_lid
!        if (lid_nbr > this%n) then
!          nintf%active = .false.; n = n + 1
!        else
!          q_nbr => this%get_quad(lid_nbr)
!          if (.not.q_nbr%get_flag(active=LDUM)) then
!            nintf%active = .false.; n = n + 1
!          end if
!        end if
!      end do
!      !
!      allocate(q%intf)
!      intf => q%intf
!      call intf_read%copy(intf)
!      call intf_read%clean(); deallocate(intf_read); intf_read => null()
!    end do
!    close(iu_bb); close(iu_intf)
!    !
!    ! read the graphs
!    call props%get_field(i_i_graph, key)
!    call props%csv%get_column(key='i_graph',i4a=i_graph_read)
!    n_disc = maxval(i_graph_read)
!    allocate(i4wk(n_disc)); i4wk = 0
!    do i = 1, nlid_read
!      i4wk(i_graph_read(i)) = 1
!    end do
!    allocate(this%qgd(n_disc))
!    !
!    iu_graph = this%open_file('.graphs.bin','r',.true.)
!    g_read => null()
!    do i_disc = 1, n_disc
!      if (i4wk(i_disc) == 1) then
!        g => this%qgd(i_disc)
!        !
!        if (associated(g_read)) then
!          call g_read%clean(); deallocate(g_read); g_read => null()
!        end if
!        allocate(g_read)
!        call this%read_graph(iu_graph, i_disc, g_read)
!        !
!        lid_max = maxval(g_read%l2gnod)
!        if (allocated(g2lnod)) deallocate(g2lnod)
!        allocate(g2lnod(lid_max)); g2lnod = 0
!        !
!        g%nvtxs = 0
!        do i = 1, g_read%nvtxs
!          lid = g_read%l2gnod(i)
!          if (lid > this%n) then
!            g_read%l2gnod(i) = -abs(g_read%l2gnod(i))
!          else
!            q => this%get_quad(lid)
!            if (.not.q%get_flag(active=LDUM)) then
!              g_read%l2gnod(i) = -abs(g_read%l2gnod(i))
!            else
!              g%nvtxs = g%nvtxs + 1
!              g2lnod(lid) = g%nvtxs ! new mapping
!            end if
!          end if
!        end do
!        !
!        ! determine the nja
!        do iact = 1, 2
!          g%nvtxs = 0
!          do i = 1, g_read%nvtxs
!            if (g_read%l2gnod(i) > 0) then
!              g%nvtxs = g%nvtxs + 1
!              if (iact == 2) then
!                g%xadj(g%nvtxs+1) = g%xadj(g%nvtxs)
!                g%l2gnod(g%nvtxs) = g_read%l2gnod(i)
!              end if
!              do j = g_read%xadj(i)+1, g_read%xadj(i+1)
!                k = g_read%adjncy(j)
!                if (g_read%l2gnod(k) > 0) then
!                  g%nja = g%nja + 1
!                  if (iact == 2) then
!                    g%xadj(g%nvtxs+1) = g%xadj(g%nvtxs+1) + 1
!                    xadj = g%xadj(g%nvtxs+1)
!                    lid = g_read%l2gnod(k)
!                    if (g2lnod(lid) == 0) then
!                      call errmsg('tQuads_init_select: program error.')
!                    end if
!                    g%adjncy(xadj) = g2lnod(lid)
!                  end if
!                end if
!              end do
!            end if
!          end do
!          if (iact == 1) then
!            allocate(g%xadj(g%nvtxs+1))
!            if (g%nja > 0) then
!              allocate(g%adjncy(g%nja))
!            end if
!            allocate(g%l2gnod(g%nvtxs))
!            g%xadj(1) = 0
!          end if
!        end do
!      end if
!    end do
!    close(iu_graph)
!    !
!    ! clean up
!    if (allocated(lid_read))     deallocate(lid_read)
!    if (allocated(gid_read))     deallocate(gid_read)
!    if (allocated(i_graph_read)) deallocate(i_graph_read)
!    if (allocated(laymod_read))   deallocate(laymod_read)
!    if (allocated(g2lnod))       deallocate(g2lnod)
!    if (associated(g_read)) then
!      call g_read%clean(); deallocate(g_read); g_read => null()
!    end if
!    !
!    return
!  end subroutine tQuads_init_select
  
  subroutine tQuads_init_select_bb(this, n_max, f_prop, uuid, dir, fields)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    integer(I4B), intent(in)     :: n_max
    character(len=*), intent(in) :: f_prop
    character(len=*), intent(in) :: uuid
    character(len=*), intent(in) :: dir
    character(len=*), dimension(:), intent(in) :: fields
    ! -- local
    type(tProps), pointer :: props => null()
    type(tQuad), pointer :: q => null()
    type(tBbObj) :: bbo
    !
    character(len=MXSLEN) :: key
    !
    integer(I4B), dimension(:), allocatable :: lid_read, gid_read
    integer(I4B) :: lid, gid, i, nlid_read
    integer(I4B) :: iu_bb
! ------------------------------------------------------------------------------
    call this%init(n_max=n_max, uuid=uuid, dir=dir)
    call this%read_props(f_prop, fields)
    props => this%get_props_ptr()
    nlid_read = props%csv%get_nr()
    !
    call props%get_field(i_lid, key)
    call props%csv%get_column(key=key,i4a=lid_read)
    !
    call props%get_field(i_gid, key)
    call props%csv%get_column(key=key,i4a=gid_read)
    !
    this%n = maxval(lid_read)
    !
    ! set the active flag
    do lid = 1, this%n
      q => this%get_quad(lid)
      call q%init()
      call q%set_flag(active=.false.)
    end do
    do i = 1, nlid_read
      lid = lid_read(i); q => this%get_quad(lid)
      call q%set_flag(active=.true.)
      this%n_act = this%n_act + 1
    end do
    !
    iu_bb   = this%open_file('.bb.bin',  'r',.true.)
    do i = 1, nlid_read
      lid = lid_read(i); gid = gid_read(i)
      q => this%get_quad(lid)
      !
      ! read the bounding box
      bbo = this%read_bb(iu_bb, lid)
      !
      call q%init(gid=gid, lid=lid, bbo=bbo); q%props => props
      call q%set(i_prop=i)
    end do
    close(iu_bb)
    !
    ! clean up
    if (allocated(lid_read)) deallocate(lid_read)
    if (allocated(gid_read)) deallocate(gid_read)
    !
    return
  end subroutine tQuads_init_select_bb
    
  subroutine tQuads_init_select_intf(this, fp_intf, read_vintf)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    character(len=*), intent(in) :: fp_intf
    logical, intent(in), optional :: read_vintf
    ! -- local
    type(tQuad), pointer :: q => null(), q_nbr => null()
    type(tIntf), pointer :: intf_read => null()
    type(tIntf), pointer :: intf => null()
    type(tNbrIntf), pointer :: nintf => null()
    !
    integer(I4B) :: iu_intf, i, lid, inbr, lid_nbr
! ------------------------------------------------------------------------------
    !
    iu_intf = this%open_file(trim(fp_intf),'r',.true.)
    !
    do i = 1, this%n
      q => this%get_quad(i)
      if (q%get_flag(active=LDUM)) then
        lid = q%lid
        !
        ! read the interface
        allocate(intf_read)
        call this%read_intf(iu_intf, lid, intf_read, read_vintf)
        !
        ! check for inactive interfaces
        do inbr = 1, intf_read%n_nbr
          nintf => intf_read%nbr_intf(inbr)
          lid_nbr = nintf%nbr_lid
          if (lid_nbr > this%n) then
            nintf%active = .false.
          else
            q_nbr => this%get_quad(lid_nbr)
            if (.not.q_nbr%get_flag(active=LDUM)) then
              nintf%active = .false.
            end if
          end if
        end do
        !
        allocate(q%intf)
        intf => q%intf
        call intf_read%copy(intf)
        call intf_read%clean(); deallocate(intf_read); intf_read => null()
      end if
    end do
    !
    close(iu_intf)
    !
    return
  end subroutine tQuads_init_select_intf
  
  subroutine tQuads_init_select_graph(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    ! -- local
    type(tProps), pointer :: props => null()
    type(tQuad), pointer :: q => null()
    type(tGraph), pointer :: g_read => null(), g => null()
    !
    character(len=MXSLEN) :: key
    !
    integer(I4B), dimension(:), allocatable :: i4wk, g2lnod, i_graph_read
    integer(I4B) :: n_disc, i_disc, nlid_read, lid, lid_max, iact, i, j, k
    integer(I4B) :: xadj, iu_graph
! ------------------------------------------------------------------------------
    !
    ! read the graphs
    props => this%get_props_ptr()
    call props%get_field(i_i_graph, key)
    call props%csv%get_column(key='i_graph',i4a=i_graph_read)
    n_disc = maxval(i_graph_read)
    allocate(i4wk(n_disc)); i4wk = 0
    do i = 1, nlid_read
      i4wk(i_graph_read(i)) = 1
    end do
    allocate(this%qgd(n_disc))
    !
    iu_graph = this%open_file('.graphs.bin','r',.true.)
    g_read => null()
    do i_disc = 1, n_disc
      if (i4wk(i_disc) == 1) then
        g => this%qgd(i_disc)
        !
        if (associated(g_read)) then
          call g_read%clean(); deallocate(g_read); g_read => null()
        end if
        allocate(g_read)
        call this%read_graph(iu_graph, i_disc, g_read)
        !
        lid_max = maxval(g_read%l2gnod)
        if (allocated(g2lnod)) deallocate(g2lnod)
        allocate(g2lnod(lid_max)); g2lnod = 0
        !
        g%nvtxs = 0
        do i = 1, g_read%nvtxs
          lid = g_read%l2gnod(i)
          if (lid > this%n) then
            g_read%l2gnod(i) = -abs(g_read%l2gnod(i))
          else
            q => this%get_quad(lid)
            if (.not.q%get_flag(active=LDUM)) then
              g_read%l2gnod(i) = -abs(g_read%l2gnod(i))
            else
              g%nvtxs = g%nvtxs + 1
              g2lnod(lid) = g%nvtxs ! new mapping
            end if
          end if
        end do
        !
        ! determine the nja
        do iact = 1, 2
          g%nvtxs = 0
          do i = 1, g_read%nvtxs
            if (g_read%l2gnod(i) > 0) then
              g%nvtxs = g%nvtxs + 1
              if (iact == 2) then
                g%xadj(g%nvtxs+1) = g%xadj(g%nvtxs)
                g%l2gnod(g%nvtxs) = g_read%l2gnod(i)
              end if
              do j = g_read%xadj(i)+1, g_read%xadj(i+1)
                k = g_read%adjncy(j)
                if (g_read%l2gnod(k) > 0) then
                  g%nja = g%nja + 1
                  if (iact == 2) then
                    g%xadj(g%nvtxs+1) = g%xadj(g%nvtxs+1) + 1
                    xadj = g%xadj(g%nvtxs+1)
                    lid = g_read%l2gnod(k)
                    if (g2lnod(lid) == 0) then
                      call errmsg('tQuads_init_select: program error.')
                    end if
                    g%adjncy(xadj) = g2lnod(lid)
                  end if
                end if
              end do
            end if
          end do
          if (iact == 1) then
            allocate(g%xadj(g%nvtxs+1))
            if (g%nja > 0) then
              allocate(g%adjncy(g%nja))
            end if
            allocate(g%l2gnod(g%nvtxs))
            g%xadj(1) = 0
          end if
        end do
      end if
    end do
    close(iu_graph)
    !
    ! clean up
    if (allocated(i_graph_read)) deallocate(i_graph_read)
    if (allocated(g2lnod))       deallocate(g2lnod)
    if (associated(g_read)) then
      call g_read%clean(); deallocate(g_read); g_read => null()
    end if
    !
    return
  end subroutine tQuads_init_select_graph
  
  function tQuads_get_props_ptr(this) result(props)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    type(tProps), pointer :: props
    ! -- local
! ------------------------------------------------------------------------------
    props => this%props
    !
    return
  end function tQuads_get_props_ptr
  
  subroutine tQuads_generate_uuid(this, luse_uuid, uuid_out)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! - modules
    use uuid_module
    ! -- dummy
    class(tQuads) :: this
    logical, intent(in) :: luse_uuid
    character(len=*), intent(in) :: uuid_out
    ! -- local
! ------------------------------------------------------------------------------
    this%uuid = trim(uuid_out)
    !
    ! overwrite if case a unique uuid is being generated
    if (luse_uuid) then
      this%uuid = generate_uuid(version=1)
    end if
    !
    return
  end subroutine tQuads_generate_uuid
  
  subroutine tQuads_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    ! -- local
    type(tProps), pointer :: p => null()
    type(tQuad), pointer :: q => null()
    type(tGraph), pointer :: g => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    if (associated(this%props)) then
      p => this%props; call p%clean()
      deallocate(this%props); this%props => null()
    end if
    !
    if (associated(this%xq)) then
      do i = 1, this%n
        q => this%xq(i)
        call q%clean()
      end do
      deallocate(this%xq)
    end if
    !
    if (associated(this%qg)) then
      g => this%qg; call g%clean()
      deallocate(this%qg); this%qg => null()
    end if
    !
    if (associated(this%qgd)) then
      do i = 1, size(this%qgd)
        g => this%qgd(i); call g%clean()
      end do
      deallocate(this%qgd); this%qgd => null()
    end if
    !
    if (associated(this%lay_mods)) then
      call this%lay_mods%clean(); this%lay_mods => null()
    end if
    if (associated(this%dat_mods)) then
      call this%dat_mods%clean(); this%dat_mods => null()
    end if
    !
    call this%init()
    !
    return
  end subroutine tQuads_clean
  
  function tQuads_get_quad(this, i) result(q)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    integer(I4B), intent(in) :: i
    type(tQuad), pointer :: q
    ! -- local
! ------------------------------------------------------------------------------
    !
    q => this%xq(i)
    !
    return
  end function tQuads_get_quad
  
  subroutine tQuads_construct_graph(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    ! -- local
    type(tGraph), pointer :: g => null()
    type(tQuad), pointer :: q => null()
    type(tQuad), pointer :: q_nbr => null()
    type(tIntf), pointer :: intf => null()
    type(tNbrIntf), pointer :: nintf => null()
    !
    logical, parameter :: LDUM = .true.
    integer(I4B) :: id, lid, iq, i, n, pos, xadj
! ------------------------------------------------------------------------------
    !
    ! !!! FOR MAPPED DATA ONLY !!!
    !
    ! init the graph
    if (associated(this%qg)) then
      call this%qg%clean()
      deallocate(this%qg)
    end if
    allocate(this%qg); g => this%qg
    !
    ! first, determine the maximum local id, contruct the graph mapping
    g%gnod_min = huge(1); g%gnod_max = 0
    do iq = 1, this%n
      q => this%xq(iq)
      g%gnod_min = min(g%gnod_min, q%lid)
      g%gnod_max = max(g%gnod_max, q%lid)
    end do
    allocate(g%g2lnod(g%gnod_max-g%gnod_min+1)); g%g2lnod = 0
    g%nvtxs = 0
    do iq = 1, this%n
      q => this%xq(iq)
      if (q%get_flag(active=LDUM)) then
        if (q%get_flag(mapped_gids=LDUM)) then
          g%nvtxs = g%nvtxs + 1
          call g%setval_g2lnod(q%lid, g%nvtxs) ! label node
        end if
      else
        call logmsg('Skipping gid '//ta([q%gid])//', area '// &
          ta([q%get_area()])//'...')
      end if
    end do
    allocate(g%l2gnod(g%nvtxs)); g%l2gnod = 0
    !
    ! determine nja
    g%nja = 0
    do iq = 1, this%n
      q => this%xq(iq)
      lid = g%getval_g2lnod(q%lid)
      if (lid > 0) then
        g%l2gnod(lid) = q%lid
        intf => q%intf
        do i = 1, intf%n_nbr
          nintf => intf%nbr_intf(i)
          if (nintf%includes_fp_stencil(intf%ncol, intf%nrow)) then
            lid = nintf%nbr_lid
            if (g%getval_g2lnod(lid) > 0) then
              g%nja = g%nja + 1
            end if
          end if
        end do
      end if
    end do
    !
    ! allocate the arrays
    allocate(g%xadj(g%nvtxs+1), g%adjncy(g%nja))
    g%xadj(1) = 0
    !
    ! construct the graph
    g%nvtxs = 0
    do iq = 1, this%n
      q => this%xq(iq)
      if (g%getval_g2lnod(q%lid) > 0) then
        intf => q%intf
        g%nvtxs = g%nvtxs + 1
        g%xadj(g%nvtxs+1) = g%xadj(g%nvtxs)
        do i = 1, intf%n_nbr
          nintf => intf%nbr_intf(i)
          if (nintf%includes_fp_stencil(intf%ncol, intf%nrow)) then
            id = g%getval_g2lnod(nintf%nbr_lid)
            if (id > 0) then
              g%xadj(g%nvtxs+1) = g%xadj(g%nvtxs+1) + 1
              xadj = g%xadj(g%nvtxs+1)
              g%adjncy(xadj) = id ! - 1 !!! minus 1 ???
            end if
          end if
        end do
      end if
    end do
    !
    return
  end subroutine tQuads_construct_graph
    
  subroutine tQuads_disconnect_graph(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    ! -- local
    type(tGraph), pointer :: g => null(), gd => null()
    type(tQuad), pointer :: q => null()
    !
    integer(I4B), dimension(:), allocatable :: x_disc, idx_area_disc
    real(R8B), dimension(:), allocatable :: area_disc
    integer(I4B), dimension(1) :: idum
    integer(I4B) :: lid, n, i, j, k, xadj, i_disc, j_disc, n_disc, iact
! ------------------------------------------------------------------------------
    !
    call logmsg('Disconnecting the graph...')
    !
    g => this%qg
    allocate(x_disc(g%nvtxs)); x_disc = 0
    ! NB. x_disc used GRAPH node numbers!
    !
    n_disc = 0
    do while(.true.)
      n_disc = n_disc + 1
      idum = minloc(x_disc); i = idum(1)
      call g%disconnect(i, n_disc, x_disc) ! calling recursive subroutine
      if (minval(x_disc) > 0) exit
    end do
    call logmsg('# Disconnected graphs: '//ta([n_disc]))
    !
    ! determine the areas and sort
    allocate(area_disc(n_disc), idx_area_disc(n_disc))
    area_disc = 0
    do i_disc = 1, n_disc
      idx_area_disc(i_disc) = i_disc
      do i = 1, g%nvtxs
        if (x_disc(i) == i_disc) then
          lid = g%l2gnod(i); q => this%xq(lid)
          call q%set(i_graph=i_disc) !set the graph number
          area_disc(i_disc) = area_disc(i_disc) - real(q%get_area(),R8B)
        end if
      end do
    end do
    call quicksort_d(area_disc, idx_area_disc, n_disc)
    deallocate(area_disc)
    !
    if (n_disc > 0) then
      ! allocate the disconnected graphs
      if (associated(this%qgd)) then
        do i = 1, size(this%qgd)
          call this%qgd(i)%clean()
        end do
        deallocate(this%qgd)
      end if
      allocate(this%qgd(n_disc))
      !
      do j_disc = 1, n_disc
        i_disc = idx_area_disc(j_disc)
        do iact = 1, 2
          gd => this%qgd(i_disc)
          gd%nvtxs = 0; gd%nja = 0
          do i = 1, g%nvtxs ! loop over the graph node numbers
            if (x_disc(i) == i_disc) then
              gd%nvtxs = gd%nvtxs + 1
              if (iact == 2) then
                gd%xadj(gd%nvtxs+1) = gd%xadj(gd%nvtxs)
                gd%l2gnod(gd%nvtxs) = g%l2gnod(i)
              end if
              do j = g%xadj(i)+1, g%xadj(i+1)
                gd%nja = gd%nja + 1
                if (iact == 2) then
                  gd%xadj(gd%nvtxs+1) = gd%xadj(gd%nvtxs+1) + 1
                  xadj = gd%xadj(gd%nvtxs+1)
                  n = g%adjncy(j)
                  gd%adjncy(xadj) = g%l2gnod(n) ! global ID !
                end if
              end do
            end if
          end do
          !
          if (iact == 1) then
            allocate(gd%xadj(gd%nvtxs+1))
            gd%xadj(1) = 0
            if (gd%nja > 0) then
              allocate(gd%adjncy(gd%nja))
            end if
            allocate(gd%l2gnod(gd%nvtxs)); gd%l2gnod = 0
          else
            ! correct the node numbers
            n = maxval(gd%l2gnod)
            allocate(gd%g2lnod(n)); gd%g2lnod = 0
            do i = 1, gd%nvtxs
              j = gd%l2gnod(i)
              gd%g2lnod(j) = i
            end do
            do i = 1, gd%nja
              j = gd%adjncy(i)
              k = gd%g2lnod(j)
              if (k == 0) then
                call errmsg('tQuads_disconnect_graph: program error')
              else
                gd%adjncy(i) = k
              end if
            end do
            deallocate(gd%g2lnod)
          end if
        end do !iact
      end do ! i_disc
    else
      call logmsg('Graph is fully connected.')
    end if
    !
    ! clean up
    deallocate(idx_area_disc)
    !
    return
  end subroutine tQuads_disconnect_graph
  
  subroutine tQuads_join(this, area_min)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(tQuads) :: this
    integer(I4B), intent(in) :: area_min
    ! -- local
    type(tQuad), pointer :: q => null(), q_nbr => null()
    type(tGraph), pointer :: gd => null()
    type(tIntf), pointer    :: intf => null(), intf_nbr => null(), intf_nbr_tmp => null()
    type(tNbrIntf), pointer :: nintf => null(), nintf_nbr => null(), nintf_nbr_tmp => null()
    type(tNbrIntf), allocatable :: nintf_test
    type(tBbObj) :: bbo, bbo_nbr
    type(tBb) :: bbi, bbi_nbr, bbi_nbr_org
    type(tBbx) :: bbx, bbx_nbr
    !
    logical :: LDUM, flag
    integer(I4B), dimension(:,:), pointer :: lcell_nr => null()
    integer(I4B) :: i_disc, i, j, k, i_nbr, n, m, mm, lid, lid_nbr, prent_lid, prent_gid, n0, n1
    integer(I4B) :: idummy, icell, ncell
    integer(I4B) :: ic, ir, jc, jr
    integer(I4B) :: ic0, ic1, ir0, ir1
    integer(I4B) :: jc0, jc1, jr0, jr1
    integer(I4B) :: n_join
    real(R8B) :: area, prent_area
! ------------------------------------------------------------------------------
    !
    do i_disc = 1, size(this%qgd)
      !
      gd => this%qgd(i_disc)
      if (gd%nvtxs == 1) then
        call logmsg('Nothing to join for disconnected 1-noded graph '// &
          ta([i_disc])//'...')
        cycle
      end if
      !
      ! count the number of parts
      n = 0
      do i = 1, gd%nvtxs
        lid = gd%l2gnod(i)
        q => this%get_quad(lid)
        if (q%get_area() < area_min) then
          n = n + 1
        end if
      end do
      if (n > 0) then
        call logmsg('Joining '//ta([n])//' quad(s)...')
        n = 0
        do while(.true.)
          n_join = 0
          do i = 1, gd%nvtxs
            lid = gd%l2gnod(i)
            q => this%get_quad(lid)
            if (q%get_flag(active=LDUM)) then
              if (q%get_area() < area_min) then
                n_join = n_join + 1
                ! loop over the neighbors
                prent_area = R8ZERO; prent_lid = -1
                do j = gd%xadj(i)+1, gd%xadj(i+1)
                  k = gd%adjncy(j); lid_nbr = gd%l2gnod(k)
                  q_nbr => this%get_quad(lid_nbr)
                  !
                  area = q_nbr%get_area()
                  if ((area >= area_min).and.(area > prent_area)) then
                    prent_lid = lid_nbr
                    prent_area = area
                  end if
                end do
                !
                ! modify the parent neighbor
                if (prent_lid > 0) then
                  q_nbr => this%get_quad(prent_lid)
                  prent_gid = q_nbr%gid
                  intf => q_nbr%intf
                  if (.not.allocated(q_nbr%gid_child)) then
                    allocate(q_nbr%gid_child(intf%n_nbr))
                    allocate(q_nbr%lid_child(intf%n_nbr))
                    q_nbr%gid_child = 0; q_nbr%lid_child = 0
                  end if
                  do j = 1, size(q_nbr%gid_child)
                    if (q_nbr%gid_child(j) == 0) then
                      q_nbr%gid_child(j) = q%gid
                      q_nbr%lid_child(j) = q%lid
                      q_nbr%n_child = q_nbr%n_child + 1
                      exit
                    end if
                  end do
                  !
                  intf     => q%intf
                  intf_nbr => q_nbr%intf
                  !
                  allocate(intf_nbr_tmp)
                  intf_nbr_tmp%n_nbr = (intf_nbr%n_nbr - 1) + (intf%n_nbr - 1)
                  allocate(intf_nbr_tmp%nbr_intf(intf_nbr_tmp%n_nbr))
                  !
                  ! initialize bounding box
                  bbo = q%bbo; bbo_nbr = q_nbr%bbo
                  bbi = q%bbo%child_bbi; bbi_nbr_org = q_nbr%bbo%child_bbi; bbi_nbr = bbi_nbr_org
                  bbx = q%bbo%child_bbx; bbx_nbr = q_nbr%bbo%child_bbx
                  !
                  bbi_nbr%ic0 = min(bbi_nbr%ic0, bbi%ic0)
                  bbi_nbr%ic1 = max(bbi_nbr%ic1, bbi%ic1)
                  bbi_nbr%ir0 = min(bbi_nbr%ir0, bbi%ir0)
                  bbi_nbr%ir1 = max(bbi_nbr%ir1, bbi%ir1)
                  bbi_nbr%ncol = bbi_nbr%ic1 - bbi_nbr%ic0 + 1
                  bbi_nbr%nrow = bbi_nbr%ir1 - bbi_nbr%ir0 + 1
                  !
                  bbx_nbr%xll = min(bbx_nbr%xll, bbx%xll)
                  bbx_nbr%xur = max(bbx_nbr%xur, bbx%xur)
                  bbx_nbr%yur = max(bbx_nbr%yur, bbx%yur)
                  bbx_nbr%yll = min(bbx_nbr%yll, bbx%yll)
                  !
                  ! copy original interfaces of the parent quad
                  m = 0
                  do k = 1, intf_nbr%n_nbr
                    nintf_nbr => intf_nbr%nbr_intf(k)
                    if (nintf_nbr%nbr_gid == q%gid) cycle
                    m = m + 1
                    nintf_nbr_tmp => intf_nbr_tmp%nbr_intf(m)
                    call nintf_nbr%copy(nintf_nbr_tmp)
                    !
                    ! set new node numbers
                    lcell_nr => nintf_nbr_tmp%lcell_nr
                    do icell = 1, size(lcell_nr,2)
                      n0 = lcell_nr(1,icell); n1 = lcell_nr(2,icell)
                      if (n0 > 0) then
                        flag = .true.
                      else
                        flag = .false.
                      end if
                      n0 = abs(n0)
                      !
                      ! first, determine the global index
                      call node_to_icrl(n0, jc0, jr0, idummy, bbi_nbr_org%ncol, bbi_nbr_org%nrow)
                      call node_to_icrl(n1, jc1, jr1, idummy, bbi_nbr_org%ncol, bbi_nbr_org%nrow)
                      ic0 = bbi_nbr_org%ic0 + jc0 - 1; ir0 = bbi_nbr_org%ir0 + jr0 - 1
                      ic1 = bbi_nbr_org%ic0 + jc1 - 1; ir1 = bbi_nbr_org%ir0 + jr1 - 1
                      !
                      ! second, convert the global index to the new local index of the parent quad
                      ic0 = ic0 - bbi_nbr%ic0 + 1; ir0 = ir0 - bbi_nbr%ir0 + 1
                      ic1 = ic1 - bbi_nbr%ic0 + 1; ir1 = ir1 - bbi_nbr%ir0 + 1
                      call icrl_to_node(lcell_nr(1,icell), ic0, ir0, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                      call icrl_to_node(lcell_nr(2,icell), ic1, ir1, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                      !
                      if (.not.flag) then
                        lcell_nr(1,icell) = -lcell_nr(1,icell)
                      end if
                    end do
                  end do
                  !
                  ! copy interfaces of quad subject to removal
                  do k = 1, intf%n_nbr
                    nintf => intf%nbr_intf(k)
                    if (nintf%nbr_gid == q_nbr%gid) cycle
                    m = m + 1
                    nintf_nbr_tmp => intf_nbr_tmp%nbr_intf(m)
                    call nintf_nbr%copy(nintf_nbr_tmp)
                    !
                    ! set new node numbers
                    lcell_nr => nintf_nbr_tmp%lcell_nr
                    do icell = 1, size(lcell_nr,2)
                      n0 = lcell_nr(1,icell); n1 = lcell_nr(2,icell)
                      if (n0 > 0) then
                        flag = .true.
                      else
                        flag = .false.
                      end if
                      n0 = abs(n0)
                      !
                      ! first, determine the global index
                      call node_to_icrl(n0, jc0, jr0, idummy, bbi%ncol, bbi%nrow)
                      call node_to_icrl(n1, jc1, jr1, idummy, bbi%ncol, bbi%nrow)
                      ic0 = bbi%ic0 + jc0 - 1; ir0 = bbi%ir0 + jr0 - 1
                      ic1 = bbi%ic0 + jc1 - 1; ir1 = bbi%ir0 + jr1 - 1
                      !
                      ! second, convert the global index to the new local index of the parent quad
                      ic0 = ic0 - bbi_nbr%ic0 + 1; ir0 = ir0 - bbi_nbr%ir0 + 1
                      ic1 = ic1 - bbi_nbr%ic0 + 1; ir1 = ir1 - bbi_nbr%ir0 + 1
                      call icrl_to_node(lcell_nr(1,icell), ic0, ir0, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                      call icrl_to_node(lcell_nr(2,icell), ic1, ir1, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                      !
                      if (.not.flag) then
                        lcell_nr(1,icell) = -lcell_nr(1,icell)
                      end if
                    end do
                  end do
                  !
                  ! determine most occuring neighbor
                  m = 0
                  do k = 1, intf_nbr_tmp%n_nbr
                    nintf => intf_nbr_tmp%nbr_intf(k)
                    if (nintf%n_cells > m) then
                      m = nintf%n_cells
                      intf_nbr_tmp%mo_nbr = nintf%nbr_gid
                    end if
                  end do
                  !
                  ! missing values
                  intf_nbr_tmp%n_mv = 0
                  if (intf_nbr%n_mv > 0) then
                    intf_nbr_tmp%n_mv = intf_nbr%n_mv
                  end if
                  if (intf%n_mv > 0) then
                    intf_nbr_tmp%n_mv = intf_nbr_tmp%n_mv + intf%n_mv
                  end if
                  if (intf_nbr_tmp%n_mv > 0) then
                    allocate(intf_nbr_tmp%mv_lcell_nr(2,intf_nbr_tmp%n_mv))
                    m = 0
                    !
                    ! add missing values from parent quad
                    lcell_nr => intf_nbr%mv_lcell_nr
                    do icell = 1, intf_nbr%n_mv
                      m = m + 1
                      !
                      n0 = lcell_nr(1,icell); n1 = lcell_nr(2,icell)
                      if (n0 > 0) then
                        flag = .true.
                      else
                        flag = .false.
                      end if
                      n0 = abs(n0)
                      !
                      call node_to_icrl(n0, jc0, jr0, idummy, bbi_nbr_org%ncol, bbi_nbr_org%nrow)
                      call node_to_icrl(n1, jc1, jr1, idummy, bbi_nbr_org%ncol, bbi_nbr_org%nrow)
                      ic0 = bbi%ic0 + jc0 - 1; ir0 = bbi%ir0 + jr0 - 1
                      ic1 = bbi%ic0 + jc1 - 1; ir1 = bbi%ir0 + jr1 - 1
                      
                      ic0 = ic0 - bbi_nbr%ic0 + 1; ir0 = ir0 - bbi_nbr%ir0 + 1
                      ic1 = ic1 - bbi_nbr%ic0 + 1; ir1 = ir1 - bbi_nbr%ir0 + 1
                      call icrl_to_node(lcell_nr(1,icell), ic0, ir0, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                      call icrl_to_node(lcell_nr(2,icell), ic1, ir1, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                      !
                      if (.not.flag) then
                        lcell_nr(1,icell) = -lcell_nr(1,icell)
                      end if
                    end do
                    !
                    ! add missing values from quad to be removed
                    lcell_nr => intf%mv_lcell_nr
                    do icell = 1, intf%n_mv
                      m = m + 1
                      !
                      n0 = lcell_nr(1,icell); n1 = lcell_nr(2,icell)
                      if (n0 > 0) then
                        flag = .true.
                      else
                        flag = .false.
                      end if
                      n0 = abs(n0)
                      !
                      call node_to_icrl(n0, jc0, jr0, idummy, bbi%ncol, bbi%nrow)
                      call node_to_icrl(n1, jc1, jr1, idummy, bbi%ncol, bbi%nrow)
                      ic0 = bbi%ic0 + jc0 - 1; ir0 = bbi%ir0 + jr0 - 1
                      ic1 = bbi%ic0 + jc1 - 1; ir1 = bbi%ir0 + jr1 - 1
                      ic0 = ic0 - bbi_nbr%ic0 + 1; ir0 = ir0 - bbi_nbr%ir0 + 1
                      ic1 = ic1 - bbi_nbr%ic0 + 1; ir1 = ir0 - bbi_nbr%ir0 + 1
                      call icrl_to_node(lcell_nr(1,icell), ic0, ir0, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                      call icrl_to_node(lcell_nr(2,icell), ic1, ir1, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                      !
                      if (.not.flag) then
                        lcell_nr(1,icell) = -lcell_nr(1,icell)
                      end if
                    end do
                  end if
                  !
                  call intf_nbr%clean()
                  call intf_nbr_tmp%copy(intf_nbr)
                  nintf_nbr_tmp => null()
                  call intf_nbr_tmp%clean(); deallocate(intf_nbr_tmp); intf_nbr_tmp => null()
                  !
                  ! update all other interfaces
                  do j = gd%xadj(i)+1, gd%xadj(i+1)
                    k = gd%adjncy(j); lid_nbr = gd%l2gnod(k)
                    if (lid_nbr == prent_lid) cycle ! skip the parent
                    q_nbr => this%get_quad(lid_nbr)
                    intf_nbr => q_nbr%intf
                    !
                    ! check if there already exist and interface to the parent
                    nintf_nbr => null(); nintf => null()
                    call intf_nbr%get_nbr_intf(nintf_nbr, prent_gid)
                    call intf_nbr%get_nbr_intf(nintf, q%gid)
                    !
                    if (associated(nintf)) then
                      if (.not.associated(nintf_nbr)) then
                        nintf%nbr_gid = prent_gid
                        nintf%nbr_lid = prent_lid
                      else
                        allocate(intf_nbr_tmp)
                        call intf_nbr%copy(intf_nbr_tmp)
                        call intf_nbr_tmp%remove(q%gid)
                        !
                        ! create new interface array
                        allocate(nintf_nbr_tmp)
                        call nintf_nbr%copy(nintf_nbr_tmp)
                        call nintf_nbr_tmp%append(nintf%lcell_nr)
                        ! replace
                        call intf_nbr_tmp%replace(prent_gid, nintf_nbr_tmp)
                        call intf_nbr%clean()
                        call intf_nbr_tmp%copy(intf_nbr)
                        ! cleanup
                        call nintf_nbr_tmp%clean(); deallocate(nintf_nbr_tmp); nintf_nbr_tmp => null()
                      end if
                    end if
                  
                  end do
                  !
                  ! deactivate the removed quad
                  call q%set_flag(active=.false.)
                  !
                  ! check the interfaces
                  q_nbr => this%get_quad(prent_lid)
                  !if (.not.q_nbr%check_gid_mask()) then
                  ! call errmsg('tQuads_join: inconsistent mask.')
                  !end if
                end if
              end if
            end if
          end do
          call logmsg('# quad joined = '//ta([n_join])//'...')
          if (n_join == 0) then
            exit
          end if
        end do
      else
        call logmsg('Nothing to join for disconnected graph '// &
          ta([i_disc])//'...')
      end if
    end do
    !
    ! cleanup
    !
    return
  end subroutine tQuads_join
    
  subroutine tQuads_split(this, area_max, g2lid)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use metis_module, only: tMetis
    !
    ! -- dummy
    class(tQuads) :: this
    !
    integer(I4B), intent(in) :: area_max
    integer(I4B), dimension(:), intent(inout) :: g2lid
    ! -- local
    type(tQuad), pointer :: qp => null(), qc => null(), qn => null()
    type(tMetis), pointer :: m => null()
    type(tIntf), pointer :: intfp => null(), intfn => null()
    type(tIntf) :: intf_new
    type(tNbrIntf), pointer :: nbr => null(), nnbr => null()
    !
    type(tBb) :: bbi, cbbi, nbbi
    type(tBbx) :: bbx
    !
    integer(I4B) :: nlid, lid, gid, gid_max, idummy
    integer(I4B) :: inbr, jnbr, iq, jq, itile
    integer(I4B) :: nc, nr, jc, jr, ic, ir, gic, gir, ic0, ic1, ir0, ir1
    integer(I4B) :: area_tot, area, n_mapped, nparts, icell, nlid_start
    integer(I4B), dimension(:,:), allocatable :: wrk2d, gids
    integer(I4B), dimension(:,:), pointer :: lcell_nr => null()
    real(R8B) :: area_avg
! ------------------------------------------------------------------------------
    !
    call logmsg('***** BEGIN splitting quads *****')
    !
    nlid = this%n
    !
    do iq = 1, this%n
      qp => this%xq(iq)
      if (qp%get_flag(mapped_gids=LDUM)) then
        call qp%get_prop(area=area)
        !
        ! split the quad
        if (area > area_max) then
          call this%get(gid_max=gid_max)
          nparts = max(area/area_max,2)
          !
          intfp => qp%intf
          call intfp%get_mask(wrk2d)
          !
          nc = size(wrk2d,1); nr = size(wrk2d,2)
          allocate(m)
          call m%init(load=wrk2d, nparts=nparts)
          call m%recur()
          call logmsg('***** Done METIS *****')
          call m%set_ids(ids=wrk2d)
          !
          call qp%get_prop(gid=gid)
          !
          allocate(gids(nc,nr))
          gids = 0
          !
          ! --- STEP 1: create the interfaces for the new quads
          call logmsg('***** Step 1 *****')
          !
          ! set the neighbors
          do inbr = 1, intfp%n_nbr
            nbr => intfp%nbr_intf(inbr)
            lcell_nr => nbr%lcell_nr
            do icell = 1, size(lcell_nr,2)
              call node_to_icrl(lcell_nr(2,icell), ic, ir, idummy, nc, nr)
              gids(ic,ir) = nbr%nbr_gid
            end do
          end do
          !
          nlid_start = nlid
          do jq = 1, nparts
            lid = nlid + jq; gid = gid_max + jq
            g2lid(gid) = lid
            !
            qc => this%xq(lid)
            call qc%init(gid=gid, lid=lid, bbo=qp%bbo, &
              gid_prent=qp%gid, lid_prent=qp%lid)
            !
            do ir = 1, nr
              do ic = 1, nc
                if (wrk2d(ic,ir) == jq) then
                  gids(ic,ir) = gid
                end if
              end do
            end do
          end do
          !
          do jq = 1, nparts
            lid = nlid + jq
            qc => this%xq(lid)
            !call qc%store_gid_mask(gids_in=gids)
            call qc%set_flag(mapped_gids=.true.)
            call qc%calc_interface(g2lid, gids, 0, qc%gid)
            call qc%calc_prop()
            !
            if (.not.qc%check_gid_mask(gids_in=gids)) then
              call errmsg('Inconsistent mask.')
            end if
            !
          end do
          !
          ! clean up
          if (allocated(wrk2d)) deallocate(wrk2d)
          if (allocated(gids))  deallocate(gids)
          !
          ! --- STEP 2: update the interfaces for the neighboring quads
          call logmsg('***** Step 2 *****')
          
          do inbr = 1, intfp%n_nbr
            nbr => intfp%nbr_intf(inbr)
            qn => this%xq(nbr%nbr_lid)
            if (.not.qn%get_flag(mapped_gids=LDUM)) cycle
            !
            nbbi = qn%bbo%child_bbi

            intfn => qn%intf
            !
            if (allocated(gids)) deallocate(gids)
            call intfn%get_mask(gids)
            ! set the global ID of the quad
            nc = size(gids,1); nr = size(gids,2)
            do ir = 1, nr; do ic = 1, nc
              if (gids(ic,ir) == 1) gids(ic,ir) = qn%gid
            end do; end do
            !
            ! set the neighbors different than the splitted quad
            do jnbr = 1, intfn%n_nbr
              nnbr => intfn%nbr_intf(jnbr)
              if (nnbr%nbr_gid /= qp%gid) then
                lcell_nr => nnbr%lcell_nr
                do icell = 1, size(lcell_nr,2)
                  call node_to_icrl(lcell_nr(2,icell), ic, ir, idummy, nc, nr)
                  gids(ic,ir) = nnbr%nbr_gid
                end do
              end if
            end do
            !
            ! add the new quads
            do jq = 1, nparts
              lid = nlid + jq; qc => this%xq(lid)
              if (allocated(wrk2d)) deallocate(wrk2d)
              call qc%intf%get_mask(wrk2d)
              cbbi = qc%bbo%child_bbi
              do ir = 1, size(wrk2d,2); do ic = 1, size(wrk2d,1)
                if (wrk2d(ic,ir) == 1) then
                  gic = cbbi%ic0 + ic - 1; gir = cbbi%ir0 + ir - 1
                  jc = gic - nbbi%ic0 + 1; jr = gir - nbbi%ir0 + 1
                  if ((jc >= 1).and.(jc <= nbbi%ncol).and. &
                      (jr >= 1).and.(jr <= nbbi%nrow)) then
                    gids(jc,jr) = qc%gid
                  end if
                end if
              end do; end do
            end do
            !
            ! replace the old
            call intfn%clean()
            call qn%calc_interface(g2lid, gids, 0, qn%gid)
            if (.not.qn%check_gid_mask(gids_in=gids)) then
              call errmsg('Inconsistent mask.')
            end if
          end do
          !
          ! increase the local IDs
          nlid = nlid + nparts
          !
          call qp%set_flag(active=.false.)
          !
          if (allocated(wrk2d)) deallocate(wrk2d)
          if (allocated(gids))  deallocate(gids)
        end if
      end if
    end do
    !
    this%n = nlid
    !
    call logmsg('***** END splitting quads *****')
    !
    return
  end subroutine tQuads_split
  
  !if (.false.) then
  !do it = it0, it1
  !  call vrt%get_bb(itile=it, dst_bbi=bbi_t)
  !  call vrt%read_full_tile(itile=it, hdrg=hdrg1)
  !  call hdrg1%get_grid(xi4=xid, mvi4=gid_mv)
  !  !
  !  allocate(xif(bbi_t%ncol,bbi_t%nrow))
  !  do ir = 1, bbi_t%nrow
  !    do ic = 1, bbi_t%ncol
  !      xif(ic,ir) = 0
  !    end do
  !  end do
  !  !
  !  do lid = 1, nlid
  !    q => xq%get_quad(lid); gid = l2gid(lid)
  !    if (q%get_flag(active=LDUM)) then
  !      call q%get_bb(child_bbi=bbi)
  !      call q%get_interface(intf)
  !      do i = 1, intf%n_nbr ! loop over the neightbors
  !        nintf => intf%nbr_intf(i)
  !        do j = 1, nintf%n_cells
  !          call node_to_icrl(nintf%lcell_nr(1,j), jc0, jr0, idummy, bbi%ncol, bbi%nrow)
  !          call node_to_icrl(nintf%lcell_nr(2,j), jc1, jr1, idummy, bbi%ncol, bbi%nrow)
  !          ic0 = bbi%ic0 + jc0 - 1; ir0 = bbi%ir0 + jr0 - 1 ! global index
  !          ic1 = bbi%ic0 + jc1 - 1; ir1 = bbi%ir0 + jr1 - 1 ! global index
  !          ic0 = ic0 - bbi_t%ic0 + 1; ir0 = ir0 - bbi_t%ir0 + 1 ! tile index
  !          ic1 = ic1 - bbi_t%ic0 + 1; ir1 = ir1 - bbi_t%ir0 + 1 ! tile index
  !          !
  !          if (xif(ic0,ir0) == 0) then
  !            xif(ic0,ir0) = gid
  !          else
  !            if (xif(ic0,ir0) /= gid) then
  !              call errmsg('Program error')
  !            end if
  !          end if
  !          if (xif(ic1,ir1) == 0) then
  !            xif(ic1,ir1) = nintf%nbr_gid
  !          else
  !            if (xif(ic1,ir1) /= nintf%nbr_gid) then
  !              call errmsg('Program error')
  !            end if
  !          end if
  !        end do
  !      end do
  !    end if
  !  end do
  !  call hdrg1%set_grid(xif, 0)
  !  fb = hdrg1%get_base_name()
  !  fp = 'f:\models\lhm\LHM-server\BRK\brk_1250_tiles_filt\'//trim(fb)//'_intf'
  !  call hdrg1%write(fp)
  !end do
  !end if
  !
  function tQuads_open_file(this, fb, act_in, lbin_in) result(iu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    character(len=*), intent(in) :: fb
    character(len=1), intent(in), optional :: act_in
    logical, intent(in), optional :: lbin_in
    integer(I4B) :: iu
    ! -- local
    character(len=MXSLEN) :: f
! ------------------------------------------------------------------------------
    f = this%get_file(trim(fb))
    call open_file(f, iu, act_in, lbin_in)
    !
    return
  end function tQuads_open_file
  
  subroutine tQuads_get(this, area_tot, n_mapped, gid_max)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    !
    integer(I4B), intent(out), optional :: area_tot
    integer(I4B), intent(out), optional :: n_mapped
    integer(I4B), intent(out), optional :: gid_max
    ! -- local
    type(tQuad), pointer :: q => null()
    integer(I4B) :: i, gid
! ------------------------------------------------------------------------------
    !
    if (present(area_tot)) then
      area_tot = I4ZERO
      do i = 1, this%n
        q => this%xq(i)
        if (q%get_flag(mapped_gids=LDUM)) then
          area_tot = area_tot + q%get_area()
        end if
      end do
    end if
    !
    if (present(n_mapped)) then
      n_mapped = I4ZERO
      do i = 1, this%n
        q => this%xq(i)
        if (q%get_flag(mapped_gids=LDUM)) then
          n_mapped = n_mapped + 1
        end if
      end do
    end if
    !
    if (present(gid_max)) then
      gid_max = I4ZERO
      do i = 1, this%n
        q => this%xq(i)
        call q%get_prop(gid=gid)
        gid_max = max(gid_max, gid)
      end do
    end if
    !
    return
  end subroutine tQuads_get
  
  function tQuads_get_file(this, s) result(f)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    character(len=*), intent(in) :: s
    !
    character(len=MXSLEN) :: f
 ! ------------------------------------------------------------------------------
    !
    f = trim(this%dir)//trim(this%uuid)//trim(s)
    !
    return
  end function tQuads_get_file

  function tQuads_get_number_active(this) result(n)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    integer(I4B) :: n
    ! -- local
    type(tQuad), pointer :: q => null()
    integer(I4B) :: i
 ! ------------------------------------------------------------------------------
    !
    n = 0
    do i = 1, this%n
      q => this%get_quad(i)
      if (q%get_flag(active=LDUM)) then
        n = n + 1
      end if
    end do
    !
    return
  end function tQuads_get_number_active
  
  subroutine tQuads_write_graphs(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    ! -- local
    type(tGraph), pointer :: g => null()
    character(len=MXSLEN) :: f
    integer(I8B), dimension(:), allocatable :: xadj
    integer(I8B) :: bytes
    integer(I4B) :: iu, ig, i, ng, lid
    integer(I4B), dimension(:), allocatable :: wrk
! ------------------------------------------------------------------------------
    !
    ! set pointer data
    ng = size(this%qgd)
    allocate(xadj(ng+1))
    bytes = sizeof(ng) + sizeof(xadj)
    !
    xadj(1) = bytes
    do ig = 1, ng
      g => this%qgd(ig)
      bytes = g%get_nbytes()
      xadj(ig+1) = xadj(ig) + bytes
    end do
    !
    f = this%get_file('.graphs.bin')
    call open_file(f, iu, 'w', .true.)
    !
    ! write pointer data
    write(iu) ng !I4B
    write(iu) xadj !I8B*size(xadj)
    !
    allocate(wrk(this%n)); wrk = 0
    do ig = 1, ng
       g => this%qgd(ig)
       call g%write(iu)
       !
       do i = 1, size(g%l2gnod)
         lid = g%l2gnod(i)
         if ((lid == 0).or.(lid > this%n)) then
           call errmsg('tQuads_write_graphs: program error.')
         else
           wrk(lid) = ig
         end if
       end do
    end do
    !
    close(iu)
    !
    !f = this%get_file('.lid2graph.bin')
    !call open_file(f, iu, 'w', .true.)
    !write(iu) this%n
    !write(iu) wrk
    !close(iu)
    !
    deallocate(xadj, wrk)
    !
    !iu = this%open_file('.graphs.bin','r',.true.)
    !ig = 2
    !call this%read_graph(iu, ig, g)
    
    return
  end subroutine tQuads_write_graphs
  
  subroutine tQuads_write_intf(this, fp)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    character(len=*), intent(in) :: fp
    ! -- local
    type(tQuad), pointer :: q => null()
    type(tIntf), pointer :: intf => null()
    !
    character(len=MXSLEN) :: f
    integer(I8B), dimension(:), allocatable :: xadj
    integer(I8B) :: bytes
    integer(I4B) :: lid, iu
! ------------------------------------------------------------------------------
    !
    ! set pointer data
    allocate(xadj(this%n+1))
    bytes = sizeof(this%n) + sizeof(xadj)
    !
    xadj(1) = bytes
    do lid = 1, this%n
      q => this%get_quad(lid)
      if (.not.q%get_flag(active=LDUM)) then
        bytes = 0
      else
        intf => q%intf
        bytes = intf%get_nbytes()
      end if
      xadj(lid+1) = xadj(lid) + bytes
    end do
    !
    f = this%get_file(trim(fp))
    call open_file(f, iu, 'w', .true.)
    !
    ! write pointer data
    write(iu) this%n !I4B
    write(iu) xadj !I8B*size(xadj)
    !
    do lid = 1, this%n
      q => this%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        intf => q%intf
        call intf%write(iu)
      end if
    end do
    !
    deallocate(xadj) 
    close(iu)
    !
    return
  end subroutine tQuads_write_intf
  
  subroutine tQuads_reorder_intf(this)
! ******************************************************************************
!
! lcell_nr <--> lcell_nr should be 1-one-1
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    !
    ! -- local
    type(tQuad), pointer    :: q     => null(), q_nbr     => null()
    type(tIntf), pointer    :: intf  => null(), intf_nbr  => null()
    type(tNbrIntf), pointer :: nintf => null(), nintf_nbr => null()
    type(tBB) :: bbi, bbi_nbr
    type(tBBx) :: bbx, bbx_nbr
    !
    character(len=MXSLEN) :: fp
    logical :: found
    integer(I4B), dimension(:,:,:), allocatable :: has_ups_nod
    integer(I4B) :: lid, lid_nbr, inbr, jnbr, icell, idummy
    integer(I4B) :: n1, n2, m1, an1, an2, am1, am2
    integer(I4B) :: ic1, ir1, gic1, gir1, ic2, ir2, gic2, gir2
    integer(I4B) :: jc1, jr1, jc2, jr2, i_dir
! ------------------------------------------------------------------------------
    do lid = 1, this%n
      q => this%get_quad(lid)
      !
      if (q%get_flag(active=LDUM)) then
        intf => q%intf
        call q%get_bb(child_bbi=bbi, child_bbx=bbx)
        do inbr = 1, intf%n_nbr
          nintf => intf%nbr_intf(inbr)
          lid_nbr = nintf%nbr_lid
          q_nbr => this%get_quad(lid_nbr)
          if (q_nbr%get_flag(active=LDUM)) then
            !
            if (q%lid < q_nbr%lid) then
              intf_nbr => q_nbr%intf
              call q_nbr%get_bb(child_bbi=bbi_nbr, child_bbx=bbx_nbr)
              !
              ! find the interface
              found = .false.
              do jnbr = 1,intf_nbr%n_nbr
                nintf_nbr => intf_nbr%nbr_intf(jnbr)
                if (nintf_nbr%nbr_gid == q%gid) then
                  found = .true.
                  exit
                end if
              end do
              if (.not.found) then
                call errmsg('tQuads_reorder_intf: program error 1.')
              end if
              !
              if (nintf_nbr%n_cells /= nintf%n_cells) then
                call errmsg('tQuads_reorder_intf: program error 2.')
              end if
              !
              ! modify the interface cells of the neighbor
              if (allocated(has_ups_nod)) deallocate(has_ups_nod)
              allocate(has_ups_nod(nst_np, intf_nbr%ncol, intf_nbr%nrow)); has_ups_nod = 0
              !
              ! first, flag the "from" cell of the neighbor having upstream nodes
              do icell = 1, nintf%n_cells
                n1 = nintf_nbr%lcell_nr(1,icell); an1 = abs(n1)
                n2 = nintf_nbr%lcell_nr(2,icell); an2 = abs(n2)
                call node_to_icrl(an1, ic1, ir1, idummy, intf_nbr%ncol, intf_nbr%nrow)
                call node_to_icrl(an2, ic2, ir2, idummy, intf_nbr%ncol, intf_nbr%nrow)
                i_dir = get_direction(ic1, ir1, ic2, ir2)
                !
                if (has_ups_nod(i_dir,ic1,ir1) /= 0) then
                  call errmsg('tQuads_reorder_intf: program error 3.')
                end if
                if (n1 > 0) then
                  has_ups_nod(i_dir,ic1,ir1) = 1
                else
                  has_ups_nod(i_dir,ic1,ir1) = -1
                end if
              end do
              nintf_nbr%lcell_nr = 0
              !
              !if (.false.) then
              !  fp = 'f:\models\lhm\LHM-Flex\pre-processing\debug\'//&
              !    ta([q_nbr%gid])//'has_ups_nod'
              !  call writeflt(fp, has_ups_nod, size(has_ups_nod,1), size(has_ups_nod,2), &
              !    bbx_nbr%xll, bbx_nbr%yll, bbx_nbr%cs, 0)
              !end if
              !
              ! second, set the cells
              do icell = 1, nintf%n_cells ! local node number for q!
                n1 = nintf%lcell_nr(1,icell); n2 = nintf%lcell_nr(2,icell) ! from/to
                !
                ! get the local icol/irow
                an1 = abs(n1); an2 = abs(n2) ! absolute values
                call node_to_icrl(an1, ic1, ir1, idummy, intf%ncol, intf%nrow)
                call node_to_icrl(an2, ic2, ir2, idummy, intf%ncol, intf%nrow)
                !
                ! convert to icol/irow of the global extent
                gic1 = bbi%ic0 + ic1 - 1; gir1 = bbi%ir0 + ir1 - 1
                gic2 = bbi%ic0 + ic2 - 1; gir2 = bbi%ir0 + ir2 - 1
                !
                ! convert to icol/irow of the neighbor; invert nodes
                jc2 = gic1 - bbi_nbr%ic0 + 1; jr2 = gir1 - bbi_nbr%ir0 + 1
                jc1 = gic2 - bbi_nbr%ic0 + 1; jr1 = gir2 - bbi_nbr%ir0 + 1
                i_dir = get_direction(jc1, jr1, jc2, jr2)
                !
                ! checks
                if ((.not.valid_icr(jc1, jr1, bbi_nbr%ncol, bbi_nbr%nrow)).or. &
                    (.not.valid_icr(jc2, jr2, bbi_nbr%ncol, bbi_nbr%nrow))) then
                  call errmsg('tQuads_reorder_intf: program error 4.')
                end if
                if (has_ups_nod(i_dir,jc1,jr1) == 0) then
                  call errmsg('tQuads_reorder_intf: program error 5.')
                end if
                !
                ! convert to the node number
                call icrl_to_node(am1, jc1, jr1, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                call icrl_to_node(am2, jc2, jr2, 1, bbi_nbr%ncol, bbi_nbr%nrow)
                !
                ! check for upstream present and store
                if (has_ups_nod(i_dir,jc1,jr1) < 0) then
                  m1 = -am1
                else
                  m1 = am1
                end if
                nintf_nbr%lcell_nr(1,icell) = m1
                nintf_nbr%lcell_nr(2,icell) = am2
              end do
            end if
          else
            call errmsg('tQuads_reorder_intf: program error 5.')
          end if
        end do
      end if
    end do
    !
    ! clean up
    if (allocated(has_ups_nod)) deallocate(has_ups_nod)
    !
    return
  end subroutine tQuads_reorder_intf
  
  subroutine tQuads_read_intf(this, iu, lid, intf, read_vintf)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: lid
    type(tIntf), intent(inout), pointer :: intf
    logical, intent(in), optional :: read_vintf
    !
    ! -- local
    integer(I4B) :: nlid
    integer(I8B) :: p, p0, p1
! ------------------------------------------------------------------------------
    !
    read(iu,pos=1) nlid
    if (lid > nlid) then
      call errmsg('tQuads_read_intf: lid > nlid')
    end if
    !
    read(iu,pos=1+I4B+(lid-1)*I8B) p0; p0 = p0 + 1
    read(iu,pos=1+I4B+lid*I8B)     p1; p1 = p1 - 1
    !
    if (p0 /= p1) then
      p = p0
      call intf%read(iu, pos=p0, read_vintf=read_vintf)
      if (intf%n_nbr == 0) then
        call logmsg('No neighbors found for lid = '//ta([lid]))
      end if
    end if
    !
    return
  end subroutine tQuads_read_intf
  
  subroutine tQuads_write_bb(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    ! -- local
    type(tQuad), pointer :: q => null()
    type(tBb) :: bbi
    type(tBbx):: bbx
    !
    character(len=MXSLEN) :: f
    integer(I8B), dimension(:), allocatable :: xadj
    integer(I8B) :: pos, bytes
    integer(I4B) :: lid, iu
! ------------------------------------------------------------------------------
    !
    ! set pointer data
    allocate(xadj(this%n+1))
    bytes = sizeof(this%n) + sizeof(xadj)
    !
    xadj(1) = bytes
    do lid = 1, this%n
      q => this%get_quad(lid)
      if (.not.q%get_flag(active=LDUM)) then
        bytes = 0
      else
        bytes = 0
        call q%get_bb(prent_bbi=bbi); bytes = bytes + sizeof(bbi)
        call q%get_bb(prent_bbx=bbx); bytes = bytes + sizeof(bbx)
        call q%get_bb(child_bbi=bbi); bytes = bytes + sizeof(bbi)
        call q%get_bb(child_bbx=bbx); bytes = bytes + sizeof(bbx)
      end if
      xadj(lid+1) = xadj(lid) + bytes
    end do
    !
    f = this%get_file('.bb.bin')
    call open_file(f, iu, 'w', .true.)
    !
    ! write pointer data
    write(iu) this%n !I4B
    write(iu) xadj !I8B*size(xadj)
    deallocate(xadj) 
    !
    do lid = 1, this%n
      q => this%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        call q%get_bb(prent_bbi=bbi); call bbi%write(iu)
        call q%get_bb(prent_bbx=bbx); call bbx%write(iu)
        call q%get_bb(child_bbi=bbi); call bbi%write(iu)
        call q%get_bb(child_bbx=bbx); call bbx%write(iu)
      end if
    end do
    !
    close(iu)
    !
    return
  end subroutine tQuads_write_bb
  
  function tQuads_read_bb(this, iu, lid) result(bbo)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: lid
    type(tBbObj) :: bbo
    ! -- local
    integer(I4B) :: nlid
    integer(I8B) :: p, p0, p1
! ------------------------------------------------------------------------------
    !
    read(iu,pos=1) nlid
    if (lid > nlid) then
      call errmsg('tQuads_read_bb: lid > nlid')
    end if
    !
    read(iu,pos=1+I4B+(lid-1)*I8B) p0; p0 = p0 + 1
    read(iu,pos=1+I4B+lid*I8B)     p1; p1 = p1 - 1
    !
    if (p0 /= p1) then
      p = p0
      call bbo%prent_bbi%read(iu, pos=p); p = p + sizeof(bbo%prent_bbi)
      call bbo%prent_bbx%read(iu, pos=p); p = p + sizeof(bbo%prent_bbx)
      call bbo%child_bbi%read(iu, pos=p); p = p + sizeof(bbo%child_bbi)
      call bbo%child_bbx%read(iu, pos=p)
    end if
    !
    return
  end function tQuads_read_bb
  
  subroutine tQuads_read_graph(this, iu, i_graph, g)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: i_graph
    type(tGraph), pointer    :: g
    ! -- local
    integer(I4B) :: n_disc
    integer(I8B) :: p, p0, p1
! ------------------------------------------------------------------------------
    read(iu,pos=1) n_disc
    if (i_graph > n_disc) then
      call errmsg('tQuads_read_graph: i_graph > n_disc')
    end if
    !
    read(iu,pos=1+I4B+(i_graph-1)*I8B) p0; p0 = p0 + 1
    read(iu,pos=1+I4B+i_graph*I8B)     p1; p1 = p1 - 1
    !
    if (p0 /= p1) then
      p = p0
      call g%read(iu, pos=p0)
    end if
    !
    return
  end subroutine tQuads_read_graph
  
  subroutine tQuads_balance_graphs(this, lwrite_props, f_bal_csv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    logical, intent(in), optional :: lwrite_props
    character(len=MXSLEN), intent(in), optional :: f_bal_csv
    ! -- local
    type(tProps), pointer :: props => null()
    type(tQuad), pointer :: q => null()
    type(tGraph), pointer :: g => null()
    integer(I4B) :: ig, i, lid, cnt, lay_mod, ic, ir
    integer(I8B) :: n, max_lay_mod
! ------------------------------------------------------------------------------
    !
    props => this%props
    !
    do ig = 1, size(this%qgd)
      g => this%qgd(ig)
      if (allocated(g%vwgt)) deallocate(g%vwgt)
      allocate(g%vwgt(g%nvtxs))
      !
      do i = 1, g%nvtxs
        lid = g%l2gnod(i)
        q => this%get_quad(lid)
        if (.not.q%get_flag(active=LDUM)) then
          call errmsg('tQuads_balance_graphs: program error.')
        end if
        call q%get_prop_csv(ikey=i_lay_mod, i4v=lay_mod)
        g%vwgt(i) = max(1,lay_mod)
      end do
      max_lay_mod = maxval(g%vwgt)
      !
      do while(.true.)
        cnt = 0
        do n = 1, g%nvtxs
          if (g%vwgt(n) > 0) then
            call g%balance(n, g%vwgt(n), max_lay_mod, cnt)
          end if
         end do
        call logmsg('Vertices left to be smoothed: '//ta((/cnt/)))
        if (cnt == 0) exit
      end do
      !
      ! map to properties
      do i = 1, g%nvtxs
        lid = g%l2gnod(i)
        q => this%get_quad(lid)
        lay_mod = g%vwgt(i)
        call q%set_prop_csv(ikey=i_lay_mod, i4v=lay_mod)
      end do
      !
      deallocate(g%vwgt)
      !
    end do
    !
    if (present(lwrite_props)) then
      if (.not.present(f_bal_csv)) then
        call errmsg('Could not write balanced file.')
      end if
      this%props%csv%file = trim(f_bal_csv)
      call this%props%csv%write()
    end if
    !
    return
  end subroutine tQuads_balance_graphs
    
  subroutine tQuads_read_props(this, f, fields)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    character(len=*), intent(in) :: f
    character(len=*), dimension(:), intent(in) :: fields
    ! -- local
    type(tQuad), pointer :: q => null()
    type(tProps), pointer :: props => null()
    type(tCSV), pointer :: csv => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    if (associated(this%props)) then
      props => this%props; call props%clean()
      deallocate(this%props); this%props => null()
    end if
    allocate(this%props); props => this%props
    call props%init(f, fields)
    !
    return
  end subroutine tQuads_read_props
      
  subroutine tQuads_write_props(this, f_csv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this
    character(len=*), intent(in), optional :: f_csv
    ! -- local
    type(tQuad), pointer :: q => null()
    type(tCSV), pointer :: csv => null()
    type(tProps), pointer :: props => null() 
    logical :: LDUM, active, mapped_gids, generated
    character(len=MXSLEN) :: s, f
    integer(I4B) :: i, lid, gid, lid_prent, gid_prent, area, i_graph, i_bnd
    real(R8B) :: xm, ym, cs_max
! ------------------------------------------------------------------------------
    !
    if (associated(this%props)) then
      call this%props%clean()
      deallocate(this%props); this%props => null()
    end if
    allocate(this%props); props => this%props
    allocate(props%csv); csv => props%csv
    !
    if (present(f_csv)) then
      f = f_csv
    else
      f = trim(this%dir)//trim(this%uuid)//'.props.csv'
    end if
    !
    call csv%init(file=f, &
      hdr_keys=['lid', 'gid', 'lid_prent', 'gid_prent', &
        'xm', 'ym', 'cs_max', 'area', 'i_graph', 'has_bnd'],&
      nr=this%n, hdr_i_type=[i_i4, i_i4, i_i4, i_i4, &
        i_r8, i_r8, i_r8, i_i4, i_i4, i_i4])
    !
    do i = 1, this%n
      q => this%xq(i)
      active      = q%get_flag(active=LDUM)
      mapped_gids = q%get_flag(mapped_gids=LDUM)
      generated   = q%get_flag(generated=LDUM)
      !
      lid_prent = int(R8MV,I4B)
      gid_prent = int(R8MV,I4B)
      area      = int(R8MV,I4B)
      xm        = R8MV
      ym        = R8MV
      cs_max    = R8MV
      i_graph   = int(R8MV,I4B)
      i_bnd     = int(R8MV,I4B)
      !
      if (mapped_gids) then
        if (active) then
          call q%get_prop(lid=lid, gid=gid, area=area, &
            xm=xm, ym=ym, cs_max=cs_max, i_graph=i_graph)
          if (q%get_flag(has_bnd=LDUM)) then
             i_bnd = 1
          else
             i_bnd = 0
          end if
          if (generated) then
            call q%get(lid_prent=lid_prent)
            call q%get(gid_prent=gid_prent)
          end if
        else
          ! inactive (merged to other quads)
          call q%get_prop(lid=lid, gid=gid)
        end if
      else
        ! not mapped
        call q%get_prop(lid=lid, gid=gid)
      end if
      !
      call q%get_prop(lid=lid, gid=gid)
      call csv%set_val(ic= 1, ir=i, i4v=lid)
      call csv%set_val(ic= 2, ir=i, i4v=gid)
      call csv%set_val(ic= 3, ir=i, i4v=lid_prent)
      call csv%set_val(ic= 4, ir=i, i4v=gid_prent)
      call csv%set_val(ic= 5, ir=i, r8v=xm)
      call csv%set_val(ic= 6, ir=i, r8v=ym)
      call csv%set_val(ic= 7, ir=i, r8v=cs_max)
      call csv%set_val(ic= 8, ir=i, i4v=area)
      call csv%set_val(ic= 9, ir=i, i4v=i_graph)
      call csv%set_val(ic=10, ir=i, i4v=i_bnd)
      !
    end do
    call csv%write(r8mv=R8MV)
    !call csv%clean(); deallocate(csv)
    !
    return
  end subroutine tQuads_write_props
  
subroutine tQuads_add_lm_intf(this, f_out_csv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this 
    character(len=*), optional :: f_out_csv
    ! -- local
    logical :: read_all_lm = .false.
    integer(I4B), parameter :: perc_intv = 10
    !
    type(tBb) :: bbi
    type(tBbX) :: bbx
    type(tQuad), pointer :: q => null(), q_nbr => null()
    type(tIntf), pointer :: intf => null(), intf_nbr => null()
    type(tNbrIntf), pointer :: nintf => null(), nintf_nbr => null()
    type(tLayerModel), pointer :: lay_mod  => null()
    type(tVertIntf), pointer :: vintf => null(), vintf_nbr => null()
    type(tMF6Exchange), pointer :: xch => null()
    type(tCsv), pointer :: csv => null()
    type(tLookupTable), pointer :: lookup => null()
    !
    character(len=:), allocatable :: s
    character(len=MXSLEN) :: f_flt
    logical :: found, fp, lwrite_props
    real(R8B), dimension(:), allocatable :: x, y, cs
    real(R4B), dimension(:,:), allocatable :: zp
    real(R4B) :: mv, top, bot, thk, perc
    integer(I4B), dimension(:,:), allocatable :: i4wk, nod
    integer(I4B), dimension(:), allocatable :: nl_act, lay_act, nbr_lid, nbr_nexg
    integer(I4B) :: lid, inbr, jnbr, lid_nbr, lm, icell, jcell, gicell
    integer(I4B) :: q_nact, q_nact_delta
    integer(I4B) :: i, n, ic, ir, il, idummy, nl, nl_intf, nr, nc, ncell, iact
! ------------------------------------------------------------------------------
    !
    lwrite_props = .false.
    if (present(f_out_csv)) then
      lwrite_props = .true.
    end if
    !
    if (lwrite_props) then
      if (len_trim(f_out_csv) == 0) then
        call errmsg('f_out_csv not set.')
      end if
    end if
    !
    ! read all the layers models
    if (read_all_lm) then
      do lm = 1, this%lay_mods%n_inp_mod
        lay_mod => this%lay_mods%lay_mods(lm)
        call lay_mod%read_full_grid()
      end do
    end if
    !
    ! determine the number of active quads
    q_nact_delta = max(1,int(real(this%n_act,R4B)/perc_intv,I4B))
    q_nact = 0
    !
    ! first, process all "from" interface nodes
    do lid = 1, this%n
      q => this%get_quad(lid)
      !
      if (q%get_flag(active=LDUM)) then
        q_nact = q_nact + 1
        if ((q_nact == 1).or.(mod(q_nact,q_nact_delta)==0).or.(q_nact==this%n_act)) then
          perc = 100.*real(q_nact,R4B)/this%n_act
          call logmsg('***** Calculate "from" interface nodes ('//ta([perc],fmt_in='(f6.2)')//'%) *****')
        end if
        !call logmsg('***** Processing quad '//ta([q%gid])//' *****')
        call q%get_bb(child_bbx=bbx, child_bbi=bbi)
        if (allocated(i4wk)) deallocate(i4wk)
        nc = bbi%ncol; nr = bbi%nrow
        allocate(i4wk(nc, nr)); i4wk = 0
        !
        ! first, label the entire interface (5-point stencil)
        intf => q%intf
        ncell = 0
        do inbr = 1, intf%n_nbr ! loop over interfaces
          nintf => intf%nbr_intf(inbr)
          if (nintf%active) then
            do icell = 1, nintf%n_cells
              fp = fp_dir(nintf%lcell_nr(:,icell), intf%ncol, intf%nrow) ! 5-point stencil only !
              n = nintf%lcell_nr(1,icell)
              call node_to_icrl(abs(n), ic, ir, idummy, intf%ncol, intf%nrow)
              if (i4wk(ic,ir) == 0) then
                ncell = ncell + 1
                if (fp) then
                  i4wk(ic,ir) = ncell
                else
                  i4wk(ic,ir) = -ncell
                end if
              end if
              if ((i4wk(ic,ir) < 0).and.(fp)) then
                i4wk(ic,ir) = -i4wk(ic,ir)
              end if
            end do
          end if
        end do
        if (.false.) then
          f_flt = 'f:\models\lhm\LHM-Flex\pre-processing\debug\'//ta([q%gid])//'_number'
          call writeflt(f_flt, i4wk, size(i4wk,1), size(i4wk,2), &
            bbx%xll, bbx%yll, bbx%cs, 0)
        end if
        
        !do icell = 1, intf%n_mv ! loop over missing values
        !  fp = fp_dir(intf%mv_lcell_nr(:,icell), intf%ncol, intf%nrow) ! 5-point stencil only !
        !  n = intf%mv_lcell_nr(1,icell)
        !  call node_to_icrl(abs(n), ic, ir, idummy, intf%ncol, intf%nrow)
        !  if (i4wk(ic,ir) == 0) then
        !    ncell = ncell + 1
        !    if (fp) then
        !      i4wk(ic,ir) = ncell
        !    else
        !      i4wk(ic,ir) = -ncell
        !    end if
        !    if ((i4wk(ic,ir) < 0).and.(fp)) then
        !      i4wk(ic,ir) = -i4wk(ic,ir)
        !    end if
        !  end if
        !end do
        !
        ! do nothing in case no interface cells are found
        if (ncell == 0) then
          call logmsg('Skipping, no interface cells for lid = '//ta([q%lid])//'...')
          cycle
        end if
        !
        ! determine the layer model data for ALL "from" nodes
        if (allocated(x)) deallocate(x)
        if (allocated(y)) deallocate(y)
        if (allocated(cs)) deallocate(cs)
        allocate(x(ncell), y(ncell), cs(ncell))
        !
        do inbr = 1, intf%n_nbr
          nintf => intf%nbr_intf(inbr)
          if (nintf%active) then
            do icell = 1, nintf%n_cells
              n = nintf%lcell_nr(1,icell)
              call node_to_icrl(abs(n), ic, ir, idummy, intf%ncol, intf%nrow)!
              jcell = abs(i4wk(ic,ir))
              call get_xy(x(jcell), y(jcell), ic, ir, bbx%xll, bbx%yur, bbx%cs)
              cs(jcell) = bbx%cs
            end do
          end if
        end do
        !do icell = 1, intf%n_mv
        !  n = intf%mv_lcell_nr(1,icell)
        !  jcell = abs(i4wk(ic,ir))
        !  call node_to_icrl(abs(n), ic, ir, idummy, intf%ncol, intf%nrow)
        !  call get_xy(x(jcell), y(jcell), ic, ir, bbx%xll, bbx%yur, bbx%cs)
        !  cs(jcell) = bbx%cs
        !end do
        !
        call q%get_prop_csv(ikey=i_lay_mod, i4v=lm)
        lay_mod => q%lay_mods%lay_mods(lm)
        call lay_mod%read_xy(x, y, cs, zp, mv, nl_act, nl)
        if (allocated(lay_act)) deallocate(lay_act)
        allocate(lay_act(nl)) 
        do il = 1, lay_mod%nlay
          i = nl_act(il)
          if (i > 0) lay_act(i) = il
        end do
        !
        if (allocated(nod)) deallocate(nod)
        allocate(nod(ncell,nl))
        nod = 0
        !
        ! set a node number for a valid cell (5-p stencil)
        n = 0
        do il = 1, nl
          do ir = 1, nr; do ic = 1, nc
            icell = i4wk(ic,ir)
            if (icell > 0) then
              top = zp(icell,il); bot = zp(icell,il+1)
              if ((top /= mv).and.(bot /= mv)) then
                thk = top - bot
                if (abs(thk) > R4TINY) then
                  n = n + 1
                  nod(icell,il) = n
                end if
              end if
            end if
          end do; end do
        end do
        i4wk = abs(i4wk)
        !
        do inbr = 1, intf%n_nbr
          nintf => intf%nbr_intf(inbr)
          if (nintf%active) then
            nintf%vintf_active = 1
            allocate(nintf%vintf_from)
            vintf => nintf%vintf_from
            call vintf%init(lm=lm, nlay_act=nl, n_cells=nintf%n_cells)
            !
            nl_intf = 0
            do il = 1, lay_mod%nlay
              i = nl_act(il)
              if (i > 0) then
                nl_intf = nl_intf + 1
                vintf%lay_act(nl_intf) = il
              end if
            end do
            !
            do icell = 1, nintf%n_cells ! loop over all interface cells
              fp = fp_dir(nintf%lcell_nr(:,icell), intf%ncol, intf%nrow) ! 5-point stencil only !
              n = nintf%lcell_nr(1,icell)
              call node_to_icrl(abs(n), ic, ir, idummy, intf%ncol, intf%nrow)!
              jcell = i4wk(ic,ir)
              if (jcell == 0) then
                call errmsg('tQuads_add_lm_intf: program error 1.')
              end if
              !
              if (fp) then
                do i = 1, nl
                  vintf%nod(icell,i) = nod(jcell,i)
                end do
              end if
              !
              do i = 1, nl
                vintf%zp(icell,i) = zp(jcell,i)
              end do
              vintf%zp(icell,nl+1) = zp(jcell,nl+1)
            end do
            !
          end if
          !
        end do
      end if
    end do
    !
    ! second, process all "to" interface nodes by copying
    q_nact = 0
    do lid = 1, this%n
      q => this%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        q_nact = q_nact + 1
        if ((q_nact == 1).or.(mod(q_nact,q_nact_delta)==0).or.(q_nact==this%n_act)) then
          perc = 100.*real(q_nact,R4B)/this%n_act
          call logmsg('***** Copy "from" interface nodes ('//ta([perc],fmt_in='(f6.2)')//'%) *****')
        end if
        intf => q%intf
        do inbr = 1, intf%n_nbr ! loop over interfaces
          nintf => intf%nbr_intf(inbr)
          if (nintf%active) then
            vintf => nintf%vintf_from
            if (.not.associated(vintf)) then
              call errmsg('tQuads_add_lm_intf: program error 1.')
            end if
            !
            ! find the interface
            lid_nbr = nintf%nbr_lid
            q_nbr => this%get_quad(lid_nbr)
            !
            if (.not.q_nbr%get_flag(active=LDUM)) then
              call errmsg('tQuads_add_lm_intf: program error 2.')
            end if
            intf_nbr => q_nbr%intf
            found = .false.
            do jnbr = 1, intf_nbr%n_nbr ! loop over interfaces
              nintf_nbr => intf_nbr%nbr_intf(jnbr)
              if (nintf_nbr%active) then
                if (nintf_nbr%nbr_gid == q%gid) then
                  found = .true.
                  exit
                end if
              end if
            end do
            if (.not.found) then
              call errmsg('tQuads_add_lm_intf: program error 3.')
            end if
            !
            ! copy the data
            allocate(nintf%vintf_to, source=nintf_nbr%vintf_from)
            !
          end if
        end do
      end if
    end do
    !
    ! determine the exchanges
    q_nact = 0
    lookup => this%lay_mods%lookup
    do lid = 1, this%n
      q => this%get_quad(lid)
      !if (q%gid /= 260076) cycle
      if (q%get_flag(active=LDUM)) then
        q_nact = q_nact + 1
        if ((q_nact == 1).or.(mod(q_nact,q_nact_delta)==0).or.(q_nact==this%n_act)) then
          perc = 100.*real(q_nact,R4B)/this%n_act
          call logmsg('***** Calculating exchanges ('//ta([perc],fmt_in='(f6.2)')//'%) *****')
        end if
        !call logmsg('***** Calculating exchanges for '//ta([q%gid])//' *****')
        !
        call q%get_bb(child_bbx=bbx)
        intf => q%intf
        do inbr = 1, intf%n_nbr ! loop over interfaces
          nintf => intf%nbr_intf(inbr)
          if (nintf%active) then
            ! only set the exchanges for the symmetric part
            if ((q%lid < nintf%nbr_lid).and.(nintf%vintf_active == 1)) then
              nintf%xch_active = 1
              allocate(nintf%xch); xch => nintf%xch
              call xch%init(ihc=1, cl1=bbx%cs*R8HALF, cl2=bbx%cs*R8HALF, hwva=bbx%cs)
              call xch%set(nintf%vintf_from, nintf%vintf_to, lookup)
              !
              ! deactivate in case no exchanges are found
              if (xch%nexg == 0) then
                call nintf%xch%clean(); deallocate(nintf%xch); nintf%xch => null()
                nintf%xch_active = 0
              end if
            end if
          end if
        end do
      end if
    end do
    !
    ! write
    call this%write_intf('.intf.lm.bin')
    if (lwrite_props) then
      csv => this%props%csv
      call csv%add_key('nr_of_neighbors')
      call csv%add_key('neighbor_lid')
      call csv%add_key('neighbor_nr_cells')
      !
      do lid = 1, this%n
        q => this%get_quad(lid)
        if (q%get_flag(active=LDUM)) then
          intf => q%intf
          do iact = 1, 2
            n = 0
            do inbr = 1, intf%n_nbr ! loop over interfaces
              nintf => intf%nbr_intf(inbr)
              if (nintf%active) then
                xch => null()
                if (nintf%xch_active == 1) then
                  n = n + 1
                  xch => nintf%xch
                else
                  ! TODO: the following search should be refactored:
                  q_nbr => this%get_quad(nintf%nbr_lid); intf_nbr => q_nbr%intf
                  found = .false.
                  do jnbr = 1, intf_nbr%n_nbr ! loop over interfaces
                    nintf_nbr => intf_nbr%nbr_intf(jnbr)
                    if (nintf_nbr%active) then
                      if (nintf_nbr%nbr_gid == q%gid) then
                        found = .true.
                        exit
                      end if
                    end if
                  end do
                  if (.not.found) then
                    call errmsg('tQuads_add_lm_intf: program error 3.')
                  end if
                  if (nintf_nbr%xch_active == 1) then
                    n = n + 1
                    xch => nintf_nbr%xch
                  end if
                end if
                if ((iact == 2).and.associated(xch)) then
                  nbr_lid(n)  = nintf%nbr_lid
                  nbr_nexg(n) = xch%nexg
                end if
              end if
            end do
            if (iact == 1) then
              if (allocated(nbr_lid)) deallocate(nbr_lid)
              if (allocated(nbr_nexg)) deallocate(nbr_nexg)
              if (n > 0) then
                allocate(nbr_lid(n), nbr_nexg(n))
              end if
            else
              call q%set_prop_csv(key='nr_of_neighbors', i4v=n)
              if (n > 0) then
                call q%set_prop_csv(key='neighbor_lid', cv=ta(nbr_lid, sep_in=';'))
                call q%set_prop_csv(key='neighbor_nr_cells', cv=ta(nbr_nexg, sep_in=';'))
                deallocate(nbr_lid, nbr_nexg)
              end if
            end if
          end do
        end if
      end do
      csv%file = trim(f_out_csv); call csv%write()
    end if
    !
    return
  end subroutine tQuads_add_lm_intf

  subroutine tQuads_write_mf6_xch_intf(this, root_dir, xch_id_field, f_out_csv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this 
    character(len=*), intent(in) :: root_dir
    character(len=*), intent(in) :: xch_id_field
    character(len=*), optional :: f_out_csv
    ! -- local
    type(tQuad), pointer :: q_m1 => null(), q_m2 => null()
    type(tIntf), pointer :: intf => null()
    type(tNbrIntf), pointer :: nintf => null()
    type(tMf6Wbd), pointer :: wbd => null()
    !
    character(len=1) :: slash
    character(len=MXSLEN) :: d, f, m1_id, m2_id, id
    integer(I4B) :: lid, inbr
! ------------------------------------------------------------------------------
    !
    ! create the directory for output
    d = root_dir
    call create_dir(d, .true.)
    !
    allocate(wbd)
    call wbd%init(f_out_csv)
    !
    slash = get_slash()
    do lid = 1, this%n
      q_m1 => this%get_quad(lid)
      if (q_m1%get_flag(active=LDUM)) then
        intf => q_m1%intf
        do inbr = 1, intf%n_nbr ! loop over interfaces
          nintf => intf%nbr_intf(inbr)
          if (nintf%active.and.(nintf%xch_active == 1)) then
            q_m2 => this%get_quad(nintf%nbr_lid)
            call q_m1%get_prop_csv(key=xch_id_field, cv=m1_id)
            call q_m2%get_prop_csv(key=xch_id_field, cv=m2_id)
            id = trim(m1_id)//'-'//trim(m2_id)
            f = trim(root_dir)//slash//'exchangedata_'//trim(id)//'.asc'
            call nintf%xch%write(f, id, wbd)
          end if
        end do
      end if
    end do
    !
    ! write the csv-file
    call wbd%write_csv()
    call wbd%clean()
    deallocate(wbd); wbd => null()
    !
    return
  end subroutine tQuads_write_mf6_xch_intf
  !
  subroutine tQuads_write_mf6_heads(this, lid0, lid1, heads_layer, &
    kper_beg, kper_end, tile_nc, tile_nr, f_vrt_pref, write_nod_map, &
    f_in_csv_merge)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this 
    integer(I4B), intent(in) :: lid0
    integer(I4B), intent(in) :: lid1
    character(len=*) :: heads_layer
    integer(I4B), intent(in) :: kper_beg
    integer(I4B), intent(in) :: kper_end
    integer(I4B), intent(in) :: tile_nc
    integer(I4B), intent(in) :: tile_nr
    character(len=*), intent(in) :: f_vrt_pref
    logical, intent(in) :: write_nod_map
    character(len=*), intent(in), optional :: f_in_csv_merge
    ! -- local
    real(R4B), parameter :: mvr4 = -99999.
    !
    type(tMF6Disu), pointer :: disu => null()
    !
    type(tBB),  dimension(:), pointer :: bbi_tile => null()
    type(tBBx), dimension(:), pointer :: bbx_tile => null()
    type(tBB) :: bbi, gbbi
    type(tBB), pointer :: bbip => null()
    type(tBBx) :: bbx, gbbx, cbbx
    type(tBBx), pointer :: bbxp => null()
    !
    type(tQuad), pointer :: q => null(), q2 => null()
    !
    type(tPostMod), pointer :: post => null()
    !
    type(tVrt), pointer :: vrt => null()
    !
    type(tCsv), pointer :: csv => null()
    !
    character(len=MXSLEN), dimension(:), allocatable :: pfix_tile, f_tile, mv_tile
    character(len=MXSLEN), dimension(:), allocatable :: f_tile_nodmap, mv_tile_nodmap
    character(len=MXSLEN) :: f, d, f_csv_dat
    character(len=MXSLENLONG) :: s_long
    integer(I4B) :: lid, nmod, nact, ntile, itile, bnc, bnr
    integer(I4B) :: ic0, ir0, ic1, ir1, ic, ir, jc, jr, kc, kr, ilm, il, jl
    integer(I4B) :: i, nodes, nodes_nodmap, kper
    integer(I4B) :: nod, nod_min, nod_max, nod_off, n, nc, nr, bs
    integer(I4B) :: im, nim, nlid, nodes_merge, lid2, nodes_offset
    integer(I4B), dimension(:), allocatable :: output_layer, i4wk
    integer(I4B), dimension(:), allocatable :: im_arr, lid2im_arr, lid_arr
    integer(I4B), dimension(:,:), allocatable :: nodmap_q, nodmap_t, tile_topol
    real(R4B), dimension(:,:), allocatable :: xr4_q, xr4_t
    real(R4B) :: r4v
    real(R8B) :: xll, yur, cs, cs_min, cs_min_rea, cxll, cxur, cyll, cyur
    !
    logical :: read_data, merge
! ------------------------------------------------------------------------------
    !
    if (present(f_in_csv_merge)) then
      merge = .true.
      allocate(csv)
      call csv%read(f_in_csv_merge)
      !
      call csv%get_column(key='lid', i4a=im_arr)
      nim = size(im_arr)
      !
      ! first, flag the quads and set im
      allocate(lid2im_arr(this%n)); lid2im_arr = 0
      do im = 1, nim
        ! get the local lids
        call csv%get_val(ir=im, ic=csv%get_col('lid_merged'), cv=s_long)
        call parse_line(s=s_long, i4a=lid_arr, token_in=';'); nlid = size(lid_arr)
        do i = 1, nlid
          lid = lid_arr(i)
          lid2im_arr(lid) = im
        end do
      end do
    else
      merge = .false.
    end if
    !
    call parse_line(s=heads_layer, i4a=output_layer, token_in=',')
    call this%lay_mods%get(n_inp_mod=nmod)
    if (size(output_layer) /= nmod) then
      call errmsg('Invalid value for heads_layer in [mf6_post]')
    end if
    !
    ! check if the output file is present
    do lid = lid0, lid1
      q => this%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        if (merge) then
          im = lid2im_arr(lid)
          if (im == 0) then
            call errmsg('tQuads_write_mf6_heads: program error 1, im=0.')
          end if
          call csv%get_val(ir=im, ic=csv%get_col(this%props%fields(i_head)), cv=f)
        else
          call q%get_prop_csv(ikey=i_head, cv=f)
        end if
        call swap_slash(f)
        if (.not.fileexist(f)) then
          call q%set_flag(active=.false.)
          this%n_act = this%n_act - 1
        end if
      end if
    end do
    !
    ! first, determine the minimum output grid size
    cs_min = HUGE(cs_min)
    do lid = lid0, lid1
      q => this%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        call q%get_prop_csv(key='cs_min_rea', r8v=cs_min_rea)
        cs_min = min(cs_min,cs_min_rea)
      end if
    end do
    !
    ! determing the bounding box and tile layout
    call gbbi%init(); call gbbx%init(); gbbx%cs = cs_min
    do lid = lid0, lid1
      q => this%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        call q%get_bb(child_bbx=bbx)
        gbbx%xll = min(gbbx%xll, bbx%xll); gbbx%xur = max(gbbx%xur, bbx%xur)
        gbbx%yll = min(gbbx%yll, bbx%yll); gbbx%yur = max(gbbx%yur, bbx%yur)
      end if
    end do
    gbbi%ncol = (gbbx%xur-gbbx%xll)/gbbx%cs
    gbbi%nrow = (gbbx%yur-gbbx%yll)/gbbx%cs
    gbbi%ic0 = 1; gbbi%ic1 = gbbi%ncol
    gbbi%ir0 = 1; gbbi%ir1 = gbbi%nrow
    !
    ! determine the tiles
    if ((tile_nc < gbbi%ncol).or.(tile_nr < gbbi%nrow)) then
      ntile = 0
      bnr = ceiling(real(gbbi%nrow,R4B)/real(tile_nr,R4B))
      bnc = ceiling(real(gbbi%ncol,R4B)/real(tile_nc,R4B))
      ntile = bnr*bnc
      allocate(bbi_tile(ntile), bbx_tile(ntile), pfix_tile(ntile), f_tile(ntile))
      allocate(tile_topol(2,ntile))
      ir0 = gbbi%ir0; ic0 = gbbi%ic0 ! upper-left
      xll = gbbx%xll; yur = gbbx%yur; cs = gbbx%cs
      !
      itile = 0
      do ir = 1, bnr; do ic = 1, bnc
        itile = itile + 1
        tile_topol(1,itile) = ir; tile_topol(2,itile) = ic
        !
        bbi_tile(itile)%ic0 = ic0 + (ic-1)*tile_nc
        bbi_tile(itile)%ic1 = bbi_tile(itile)%ic0 + tile_nc - 1
        bbi_tile(itile)%ir0 = ir0 + (ir-1)*tile_nr
        bbi_tile(itile)%ir1 = bbi_tile(itile)%ir0 + tile_nr - 1
        bbi_tile(itile)%ncol = bbi_tile(itile)%ic1 - bbi_tile(itile)%ic0 + 1
        bbi_tile(itile)%nrow = bbi_tile(itile)%ir1 - bbi_tile(itile)%ir0 + 1
        !
        bbx_tile(itile)%xll = xll + (ic-1)*tile_nc*cs
        bbx_tile(itile)%xur = bbx_tile(itile)%xll + tile_nc*cs
        bbx_tile(itile)%yur = yur - (ir-1)*tile_nr*cs
        bbx_tile(itile)%yll = bbx_tile(itile)%yur - tile_nr*cs
        bbx_tile(itile)%cs = cs
        !
        pfix_tile(itile) = '_r'//ta([ir],'(i3.3)')// &
                           '_c'//ta([ic],'(i3.3)') 
      end do; end do
    else
      ntile = 1
      allocate(bbi_tile(ntile), bbx_tile(ntile), pfix_tile(ntile), f_tile(ntile))
      allocate(tile_topol(2,ntile))
      bbi_tile = gbbi; bbx_tile = gbbx; pfix_tile = ''
      tile_topol = 1
    end if
    if (write_nod_map) then
      allocate(f_tile_nodmap(ntile))
    end if
    if (ntile > 1) then
      call logmsg('Number of tiles being used: '//ta([ntile])// &
        ' ('//ta([bnc])//', '//ta([bnr])//')')
    end if
    !
    if (.false.) then
      do lid = lid0, lid1
        q => this%get_quad(lid)
        if (q%get_flag(active=LDUM)) then
          call q%get_bb(child_bbi=bbi, child_bbx=bbx)
          call q%grid_init()
          do il = 1, q%disu%nlay_act
            jl = q%disu%lay_act(il)
            f = trim(f_vrt_pref)//'nod_map_lid'//ta([lid])//'_l'//ta([jl])
            call writeflt(f, q%disu%grid_x_nod(:,:,il), bbi%ncol, bbi%nrow, &
              bbx%xll, bbx%yll, bbx%cs, 0)
          end do
          call q%disu%clean(); deallocate(q%disu); q%disu => null()
        end if
      end do
      stop
    end if
    !
    ! loop over the tiles, read, and write
    do kper = kper_beg, kper_end
      nod_off = 0
      do itile = 1, ntile
        bbip => bbi_tile(itile)
        bbxp => bbx_tile(itile)
        !
        if (allocated(xr4_t)) deallocate(xr4_t)
        allocate(xr4_t(bbip%ncol,bbip%nrow)); xr4_t = mvr4
        if (write_nod_map) then
          if (allocated(nodmap_t)) deallocate(nodmap_t)
          allocate(nodmap_t(bbip%ncol,bbip%nrow)); nodmap_t = I4ZERO
        end if
        !
        do lid = lid0, lid1
          q => this%get_quad(lid)
          if (q%get_flag(active=LDUM)) then
            call q%get_bb(child_bbx=bbx)
            call q%get_prop_csv(key='cs_min_rea', r8v=cs_min_rea)
            if (bbx_intersect(bbx, bbxp)) then
              call q%get_prop_csv(ikey=i_lay_mod, i4v=ilm) ! layer model
              il = output_layer(ilm) ! layer for output
              !
              if (merge) then
                im = lid2im_arr(lid)
                call csv%get_val(ir=im, ic=csv%get_col('lid_merged'), cv=s_long)
                call parse_line(s=s_long, i4a=lid_arr, token_in=';'); nlid = size(lid_arr)
                nodes_offset = 0
                do i = 1, nlid
                   lid2 = lid_arr(i)
                   if (lid == lid2) exit
                   q2 => this%get_quad(lid2)
                   call q2%get_prop_csv(key='nodes', i4v=nodes)
                   nodes_offset = nodes_offset + nodes
                end do
                call csv%get_val(ir=im, ic=csv%get_col('csv_dat'), cv=f)
                call q%get_nod_map(il, nodmap_q, nodes, f, nodes_offset)
              else
                call q%get_nod_map(il, nodmap_q, nodes)
              end if
              !
              ! check dimensions
              nc = (bbx%xur-bbx%xll)/cs_min_rea; nr = (bbx%yur-bbx%yll)/cs_min_rea
              if ((size(nodmap_q,1) /= nc).or.(size(nodmap_q,2) /= nr)) then
                call errmsg('tQuads_write_mf6_heads: program error 2.')
              end if
              !
              if (allocated(nodmap_q)) then
                allocate(post)
                call logmsg('Reading heads for lid '//trim(ta([lid]))//'...')
                if (merge) then
                  im = lid2im_arr(lid) 
                  call csv%get_val(ir=im, ic=csv%get_col(this%props%fields(i_head)), cv=f)
                  call csv%get_val(ir=im, ic=csv%get_col('nodes'), i4v=nodes_merge)
                  call post%init(f, kper, kper, nodes_merge)
                  call post%read_ulasav()
                  call post%get_grid(kper, nodmap_q, mvr4, xr4_q)
                else
                  call q%get_prop_csv(ikey=i_head, cv=f); call swap_slash(f)
                  call post%init(f, kper, kper, nodes)
                  call post%read_ulasav()
                  call post%get_grid(kper, nodmap_q, mvr4, xr4_q)
                end if
                !
                if (.false.) then
                  if (write_nod_map .and. (kper == kper_beg)) then
                    f = trim(f_vrt_pref)//'nod_map_lid'//ta([lid])
                    call writeflt(f, nodmap_q, bbi%ncol, bbi%nrow, &
                    bbx%xll, bbx%yll, bbx%cs, 0)
                  end if
                end if
                !
                bs = cs_min_rea/bbxp%cs
                !
                ! set the data
                do ir = 1, nr; do ic = 1, nc ! loop over cells with resolution cs_min_rea
                  ! get the cell bounds
                  cxll = bbx%xll + (ic-1)*cs_min_rea; cxur = cxll + cs_min_rea
                  cyur = bbx%yur - (ir-1)*cs_min_rea; cyll = cyur - cs_min_rea
                  cbbx%cs = cs_min_rea
                  cbbx%xll = cxll; cbbx%xur = cxur; cbbx%yll = cyll; cbbx%yur = cyur
                  !
                  if (bbx_intersect(cbbx, bbxp)) then
                    call get_icr(ic0, ir0, cxll + bbxp%cs*R8HALF, &
                                           cyur - bbxp%cs*R8HALF, &
                      bbxp%xll, bbxp%yur, bbxp%cs)
                    call get_icr(ic1, ir1, cxur - bbxp%cs*R8HALF, &
                                           cyll + bbxp%cs*R8HALF, &
                      bbxp%xll, bbxp%yur, bbxp%cs)
                    !
                    if ((valid_icir(ic0, ir0, bbip%ncol, bbip%nrow)).and. &
                        (valid_icir(ic1, ir1, bbip%ncol, bbip%nrow))) then
                      ! check
                      if (((ic1-ic0+1)/=bs).or.((ir1-ir0+1)/=bs)) then
                        call errmsg('tQuads_write_mf6_heads: program error 3.')
                      end if
                      !
                      r4v = xr4_q(ic,ir)
                      if (r4v /= mvr4) then
                        do jr = ir0, ir1; do jc = ic0, ic1
                          xr4_t(jc,jr) = r4v
                        end do; end do
                      end if
                    end if
                  else
                    if (write_nod_map) nodmap_q(ic,ir) = 0
                  end if
                end do; end do
                !
                if (write_nod_map .and. (kper == kper_beg)) then
                  nod_min = huge(I4ZERO); nod_max = 0
                  do ir = 1, nr; do ic = 1, nc
                    nod = nodmap_q(ic,ir)
                    if (nod /= 0) then
                      nod_min = min(nod,nod_min)
                      nod_max = max(nod,nod_max)
                    end if
                  end do; end do
                  if (nod_min /= nod_max) then
                    n = nod_max - nod_min + 1
                    allocate(i4wk(n)); i4wk = 0
                    do ir = 1, nr; do ic = 1, nc
                      nod = nodmap_q(ic,ir)
                      if (nod /= 0) then
                        i = nod - nod_min + 1
                        i4wk(i) = 1
                      end if
                    end do; end do
                    nod = nod_off
                    nodes_nodmap = 0
                    do i = 1, n
                      if (i4wk(i) == 1) then
                        nod = nod + 1; nodes_nodmap = nodes_nodmap + 1
                        i4wk(i) = nod
                      end if
                    end do
                    nod_off = nod_off + nodes_nodmap
                    !
                    do ir = 1, nr; do ic = 1, nc ! loop over cells with resolution cs_min_rea
                      ! get the cell bounds
                      cxll = bbx%xll + (ic-1)*cs_min_rea; cxur = cxll + cs_min_rea
                      cyur = bbx%yur - (ir-1)*cs_min_rea; cyll = cyur - cs_min_rea
                      cbbx%cs = cs_min_rea
                      cbbx%xll = cxll; cbbx%xur = cxur; cbbx%yll = cyll; cbbx%yur = cyur
                      !
                      if (bbx_intersect(cbbx, bbxp)) then
                        call get_icr(ic0, ir0, cxll + bbxp%cs*R8HALF, &
                                               cyur - bbxp%cs*R8HALF, &
                          bbxp%xll, bbxp%yur, bbxp%cs)
                        call get_icr(ic1, ir1, cxur - bbxp%cs*R8HALF, &
                                               cyll + bbxp%cs*R8HALF, &
                          bbxp%xll, bbxp%yur, bbxp%cs)
                        !
                        if ((valid_icir(ic0, ir0, bbip%ncol, bbip%nrow)).and. &
                            (valid_icir(ic1, ir1, bbip%ncol, bbip%nrow))) then
                          n = nodmap_q(ic,ir)
                          if (n /= 0) then
                            nod = i4wk(n - nod_min + 1)
                            if (nod == 0) then
                              call errmsg('tQuads_write_mf6_heads: program error 4.')
                            end if
                            do jr = ir0, ir1; do jc = ic0, ic1
                              nodmap_t(jc,jr) = nod
                            end do; end do
                          end if
                        end if
                      end if
                    end do; end do
                    if (allocated(i4wk)) deallocate(i4wk)
                  end if
                end if
                !
                ! clean-up
                call post%clean(); deallocate(post); post => null()
              end if
            end if
          end if
        end do
        !
        ! write the tile
        f = trim(f_vrt_pref)//trim(pfix_tile(itile))//'_p'//ta([kper],'(i3.3)')
        d = get_dir(f)
        call create_dir(d, .true.)
        call writeflt(f, xr4_t, bbip%ncol, bbip%nrow, &
            bbxp%xll, bbxp%yll, bbxp%cs, mvr4)
        f_tile(itile) = trim(f)//'.flt'
        !
        if (write_nod_map .and. (kper == kper_beg)) then
          f = trim(f_vrt_pref)//trim(pfix_tile(itile))//'_nod_map'
          call writeflt(f, nodmap_t, bbip%ncol, bbip%nrow, &
              bbxp%xll, bbxp%yll, bbxp%cs, 0)
          f_tile_nodmap(itile) = trim(f)//'.flt'
        end if
      end do
      !
      ! write the VRT
      allocate(vrt)
      allocate(mv_tile(ntile)); mv_tile = ta([mvr4])
      f = trim(f_vrt_pref)//'_p'//ta([kper],'(i3.3)')
      if (allocated(tile_topol)) then
        call vrt%init_write(f, 'Float32', f_tile, mv_tile, bbi_tile, bbx_tile, tile_topol)
      else
        call vrt%init_write(f, 'Float32', f_tile, mv_tile, bbi_tile, bbx_tile)
      end if
      call vrt%write()
      deallocate(mv_tile)
      call vrt%clean(); deallocate(vrt); vrt => null()
      !
      if (write_nod_map .and. (kper == kper_beg)) then
        allocate(vrt)
        allocate(mv_tile(ntile)); mv_tile = ta([0])
        f = trim(f_vrt_pref)//'_nod_map'
        if (allocated(tile_topol)) then
          call vrt%init_write(f, 'Int32', f_tile_nodmap, mv_tile, bbi_tile, bbx_tile, tile_topol)
        else
          call vrt%init_write(f, 'Int32', f_tile_nodmap, mv_tile, bbi_tile, bbx_tile)
        end if
        call vrt%write()
        deallocate(mv_tile)
        call vrt%clean(); deallocate(vrt); vrt => null()
      end if
    end do ! kper
    !
    ! clean-up
    if (allocated(xr4_q))         deallocate(xr4_q)
    if (allocated(xr4_t))         deallocate(xr4_t)
    if (allocated(nodmap_q))      deallocate(nodmap_q)
    if (allocated(nodmap_t))      deallocate(nodmap_t)
    if (allocated(pfix_tile))     deallocate(pfix_tile)
    if (allocated(f_tile))        deallocate(f_tile)
    if (allocated(f_tile_nodmap)) deallocate(f_tile_nodmap)
    if (associated(bbi_tile))     deallocate(bbi_tile)
    if (associated(bbx_tile))     deallocate(bbx_tile)
    if (allocated(tile_topol))    deallocate(tile_topol)
    !
    return
  end subroutine tQuads_write_mf6_heads
  
  subroutine tQuads_write_mf6_bnd_heads(this, kper_beg, kper_end)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this 
    integer(I4B), intent(in) :: kper_beg
    integer(I4B), intent(in) :: kper_end
    ! -- local
    type(tQuad), pointer :: q => null(), q_nbr => null()
    type(tIntf), pointer :: intf => null(), intf_nbr => null()
    type(tNbrIntf), pointer :: nintf => null(), nintf_nbr => null()
    type(tMF6Exchange), pointer :: xch => null()
    type(tMf6Wbd), pointer :: wbd => null()
    type(tMF6Disu), pointer :: disu => null()
    type(tPostMod), pointer :: post => null()
    character(len=MXSLEN) :: f, f_bin, id
    logical :: LDUM, found
    integer(I4B), dimension(:), pointer :: m1_nodes => null()
    integer(I4B), dimension(:), allocatable :: nodes_flag, nodes_read
    integer(I4B) :: lid, nper, iper, kper, inbr, jnbr, i, nod, n
    real(R8B), dimension(:,:), allocatable :: heads_read
! ------------------------------------------------------------------------------
    nper = kper_end - kper_beg + 1
    !
    ! check if the output file is present
    do lid = 1, this%n
      q => this%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        call q%get_prop_csv(ikey=i_head, cv=f); call swap_slash(f)
        if (.not.fileexist(f)) then
          call q%set_flag(active=.false.)
          this%n_act = this%n_act - 1
        end if
      end if
    end do
    !
    do lid = 1, this%n
      q => this%get_quad(lid); disu => null(); post => null()
      if (q%get_flag(active=LDUM)) then
        ! get the grid
        call q%grid_init(); disu => q%disu; wbd => disu%wbd
        !
        call q%get_prop_csv(ikey=i_head, cv=f); call swap_slash(f)
        allocate(post)
        call post%init(f, kper_beg, kper_end, disu%nodes)
        !
        if (allocated(nodes_flag)) deallocate(nodes_flag)
        allocate(nodes_flag(disu%nodes))
        nodes_flag = 0
        !
        ! flag the nodes for output
        intf => q%intf
        do inbr = 1, intf%n_nbr
          nintf => intf%nbr_intf(inbr)
          q_nbr => this%get_quad(nintf%nbr_lid)
          if (.not.q_nbr%get_flag(active=LDUM)) then ! neighbor is inactive
            if (nintf%xch_active == 1) then
              m1_nodes => nintf_nbr%xch%cellidm1
            else
              ! find the exchange
              intf_nbr => q_nbr%intf
              found = .false.
              do jnbr = 1, intf_nbr%n_nbr
                nintf_nbr => intf_nbr%nbr_intf(jnbr)
                if (nintf_nbr%nbr_lid == q%lid) then
                  found = .true.; exit
                end if
              end do
              if (.not.found) then
                call errmsg('tQuads_write_mf6_bnd_heads: program error 1')
              end if
              if (nintf_nbr%xch_active == 0) then
                call errmsg('tQuads_write_mf6_bnd_heads: program error 2')
              end if
              m1_nodes => nintf_nbr%xch%cellidm2
            end if
            !
            do i = 1, size(m1_nodes)
              nod = m1_nodes(i)
              nodes_flag(nod) = 1
            end do
            m1_nodes => null()
          else
            !call logmsg('Doing nothing.')
          end if
        end do
        !
        n = sum(nodes_flag)
        if (n > 0) then
          if (allocated(nodes_read)) deallocate(nodes_read)
          allocate(nodes_read(n), heads_read(nper,n))
          n = 0
          do i = 1, disu%nodes
            if (nodes_flag(i) == 1) then
              n = n + 1
              nodes_read(n) = i
            end if
          end do
          !
          ! read for the boundary nodes
          call post%read_ulasav_selection(nodes_read, heads_read)
          !
          ! write the heads and the csv
          do iper = 1, nper
            kper = post%kper_map(iper)
            id = 'head_bnd_sp'//ta([kper],'(i3.3)')
            f_bin = trim(q%mod_dir)//'\'//trim(id)//'.bin'
            call wbd%write_list(id=id, i4a1=nodes_read, r8a1=heads_read(iper,:), f_bin=f_bin)
          end do
          call wbd%write_csv()
          !
        end if
        !
      end if
      !
      ! clean-up
      if (allocated(nodes_read)) deallocate(nodes_read)
      if (allocated(heads_read)) deallocate(heads_read)
      if (associated(disu)) then
        call disu%clean(); deallocate(q%disu); q%disu => null()
      end if
      if (associated(post)) then
        call post%clean(); deallocate(post); post => null()
      end if
    end do
    !
    return
  end subroutine tQuads_write_mf6_bnd_heads
    
  subroutine tQuads_set_mod_dir(this, root_dir, sub_dir_fields)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this 
    character(len=*), intent(in) :: root_dir
    character(len=*), intent(in) :: sub_dir_fields
    ! -- local
    type(tQuad), pointer :: q => null()
    character(len=MXSLEN), dimension(:), allocatable :: fields, da
    character(len=MXSLEN) :: d, sd, key
    character(len=1) :: slash
    integer(I4B) :: lid, i, j
! ------------------------------------------------------------------------------
    !
    call parse_line(sub_dir_fields, fields, token_in=',')
    allocate(da(size(fields)+1)); da = ''
    slash = get_slash()
    !
    do lid = 1, this%n
      q => this%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        da(1) = root_dir; j = 1
        do i = 1, size(fields)
          key = fields(i)
          if (len_trim(key) > 0) then
            call q%get_prop_csv(key=key, cv=sd)
          else
            sd = ''
          end if
          if (len_trim(sd) > 0) then
            j = j + 1
            da(j) = sd
          end if
        end do
        d = ta(da, sep_in=slash)
        call swap_slash(d)
        !
        ! store and create the directory
        q%mod_dir = d
        call create_dir(d, .true.)
      end if
    end do
    !
    if (allocated(fields)) deallocate(fields)
    if (allocated(da)) deallocate(da)
    !
    return
  end subroutine tQuads_set_mod_dir

  subroutine tQuads_merge_disu(this, lids, disu_arr, disu_merge)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuads) :: this 
    integer(I4B), dimension(:), intent(in) :: lids
    type(tMF6Disu), dimension(:), intent(in), pointer :: disu_arr
    type(tMF6Disu), intent(inout), pointer :: disu_merge
    ! -- local
    type(tQuad), pointer :: q_m1 => null(), q_m2 => null()
    type(tIntf), pointer :: intf => null()
    type(tNbrIntf), pointer :: nintf => null()
    type(tMF6Exchange), pointer :: xch => null()
    type(tMF6Disu), pointer :: disu => null()
    integer(I4B) :: nm, nodes, i, j, j0, k, offset, lid, inbr, lid_nbr, nexg
    integer(I4B) :: n, m, n0, n1, m0, m1, im0, im1, p, iac_max, nbr
    integer(I4B), dimension(:), allocatable :: nodes_arr, nodes_offset, lids2im
    integer(I4B), dimension(:), allocatable :: pos, iac_chk
    !
    ! work arrays
    integer(I4B), dimension(:), allocatable :: wrk_ia, wrk_ja, wrk_ihc
    integer(I4B), dimension(:), allocatable :: i4w
    real(R4B), dimension(:), allocatable :: r4w
    real(R8B), dimension(:), allocatable :: wrk_cl12, wrk_hwva
! ------------------------------------------------------------------------------
    !
    nm = size(lids)
    allocate(nodes_arr(nm), nodes_offset(nm))
    !
    do i = 1, nm
      disu => disu_arr(i)
      nodes_arr(i) = disu%nodes
    end do
    disu_merge%nodes = sum(nodes_arr)
    !
    allocate(disu_merge%top    (disu_merge%nodes))
    allocate(disu_merge%bot    (disu_merge%nodes))
    allocate(disu_merge%area   (disu_merge%nodes))
    allocate(disu_merge%idomain(disu_merge%nodes))
    allocate(disu_merge%iac    (disu_merge%nodes))
    !
    offset = 0
    do i = 1, nm
      nodes_offset(i) = offset
      disu => disu_arr(i)
      do j = 1, disu%nodes
        k = j + offset
        disu_merge%top(k)     = disu%top(j)
        disu_merge%bot(k)     = disu%bot(j)
        disu_merge%area(k)    = disu%area(j)
        disu_merge%idomain(k) = disu%idomain(j)
        disu_merge%iac(k)     = disu%iac(j)
      end do
      offset = offset + disu%nodes
    end do
    !
    allocate(lids2im(this%n)); lids2im = 0
    do i = 1, size(lids)
      lid = lids(i)
      lids2im(lid) = i
    end do
    !
    do im0 = 1, size(lids)
      lid = lids(im0)
      q_m1 => this%get_quad(lid)
      if (q_m1%get_flag(active=LDUM)) then
        intf => q_m1%intf
        do inbr = 1, intf%n_nbr ! loop over interfaces
          nintf => intf%nbr_intf(inbr)
          lid_nbr = nintf%nbr_lid
          im1 = lids2im(lid_nbr)
          if (nintf%active.and.(im1 > 0)) then ! internal interface
            if (nintf%xch_active == 1) then
              xch => nintf%xch
              do j = 1, xch%nexg
                n0 = xch%cellidm1(j); n1 = xch%cellidm2(j)
                m0 = n0 + nodes_offset(im0); m1 = n1 + nodes_offset(im1)
                disu_merge%iac(m0) = abs(disu_merge%iac(m0)) + 1
                disu_merge%iac(m1) = abs(disu_merge%iac(m1)) + 1
                !
                ! label the nodes that having exchange neigbbors
                disu_merge%iac(m0) = -disu_merge%iac(m0)
                disu_merge%iac(m1) = -disu_merge%iac(m1)
              end do
            end if
          end if
        end do
      end if
    end do
    disu_merge%nja = sum(abs(disu_merge%iac))
    allocate(iac_chk(disu_merge%nodes)); iac_chk = 0
    !
    allocate(wrk_ja  (disu_merge%nja)); wrk_ja   = 0
    allocate(wrk_ihc (disu_merge%nja)); wrk_ihc  = -1
    allocate(wrk_cl12(disu_merge%nja)); wrk_cl12 = -R8ONE
    allocate(wrk_hwva(disu_merge%nja)); wrk_hwva = -R8ONE
    !
    allocate(pos(disu_merge%nodes))
    pos(1) = 1
    do i = 2, disu_merge%nodes
      pos(i) = pos(i-1) + abs(disu_merge%iac(i-1))
    end do
    do i = 1, nm
      disu => disu_arr(i)
      do j = 1, disu%nodes
        k = disu%ia(j); n0 = disu%ja(k); m0 = n0 + nodes_offset(i)
        p = pos(m0)
        do k = disu%ia(j), disu%ia(j+1)-1
          n = disu%ja(k); m = n + nodes_offset(i)
          wrk_ja(p)   = m
          wrk_ihc(p)  = disu%ihc(k)
          wrk_cl12(p) = disu%cl12(k)
          wrk_hwva(p) = disu%hwva(k)
          iac_chk(m) = iac_chk(m) + 1
          p = p + 1
        end do
        pos(m0) = p
      end do
    end do
    !
    do im0 = 1, size(lids)
      lid = lids(im0)
      q_m1 => this%get_quad(lid)
      if (q_m1%get_flag(active=LDUM)) then
        intf => q_m1%intf
        do inbr = 1, intf%n_nbr ! loop over interfaces
          nintf => intf%nbr_intf(inbr)
          lid_nbr = nintf%nbr_lid
          im1 = lids2im(lid_nbr)
          if (nintf%active.and.(im1 > 0)) then ! internal interface
            if (nintf%xch_active == 1) then
              xch => nintf%xch
              do j = 1, xch%nexg
                n0 = xch%cellidm1(j); n1 = xch%cellidm2(j)
                m0 = n0 + nodes_offset(im0); m1 = n1 + nodes_offset(im1)
                p = pos(m0)
                wrk_ja(p)   = m1
                wrk_ihc(p)  = xch%ihc
                wrk_cl12(p) = xch%cl1
                wrk_hwva(p) = xch%hwva
                pos(m0) = p + 1
                iac_chk(m0) = iac_chk(m0) + 1
                p = pos(m1)
                wrk_ja(p)   = m0
                wrk_ihc(p)  = xch%ihc
                wrk_cl12(p) = xch%cl1
                wrk_hwva(p) = xch%hwva
                iac_chk(m1) = iac_chk(m1) + 1
                pos(m1) = p + 1
              end do
            end if
          end if
        end do
      end if
    end do
    !
    ! checks
    do i = 1, disu_merge%nodes
      if (abs(disu_merge%iac(i)) /= iac_chk(i)) then
        call errmsg('tQuads_merge_disu: program error 1.')
      end if
    end do
    !
    !
    allocate(disu_merge%ja  (disu_merge%nja))
    allocate(disu_merge%ihc (disu_merge%nja))
    allocate(disu_merge%cl12(disu_merge%nja))
    allocate(disu_merge%hwva(disu_merge%nja))
    
    ! create the ia array
    allocate(wrk_ia(disu_merge%nodes+1))
    wrk_ia(1) = 1
    do i = 2, disu_merge%nodes + 1
      wrk_ia(i) = wrk_ia(i-1) + abs(disu_merge%iac(i-1))
    end do
    !
    iac_max = maxval(abs(disu_merge%iac))
    allocate(r4w(iac_max), i4w(iac_max))
    !
    ! set the nodal adjacency data
    do i = 1, disu_merge%nodes
      ! set the node data
      j0 = wrk_ia(i)
      !
      ! set the noda data
      disu_merge%ja(j0)   = wrk_ja(j0)
      disu_merge%ihc(j0)  = wrk_ihc(j0)
      disu_merge%cl12(j0) = wrk_cl12(j0)
      disu_merge%hwva(j0) = wrk_hwva(j0)
      !
      ! set the neighbor data
      if (disu_merge%iac(i) > 0) then
        do j = j0 + 1, wrk_ia(i+1)-1
          disu_merge%ja(j)   = wrk_ja(j)
          disu_merge%ihc(j)  = wrk_ihc(j)
          disu_merge%cl12(j) = wrk_cl12(j)
          disu_merge%hwva(j) = wrk_hwva(j)
        end do
      else ! sort
        nbr = 0; r4w = R4ZERO; i4w = I4ZERO
        do j = j0 + 1, wrk_ia(i+1)-1
          nbr = nbr + 1
          r4w(nbr) = real(wrk_ja(j),R4B)
          i4w(nbr) = nbr
        end do
        !
        ! sort the neighbor node numbers
        call quicksort_r(r4w, i4w, nbr)
        !
        ! set the neighbor data
        n = 0
        do j = j0 + 1, wrk_ia(i+1)-1
          n = n + 1; k = i4w(n)
          disu_merge%ja(j)   = wrk_ja(j0+k)
          disu_merge%ihc(j)  = wrk_ihc(j0+k)
          disu_merge%cl12(j) = wrk_cl12(j0+k)
          disu_merge%hwva(j) = wrk_hwva(j0+k)
        end do
        !
      end if
      disu_merge%iac(i) = abs(disu_merge%iac(i))
    end do
    !
    ! checks
    do i = 1, disu_merge%nja
      if ((disu_merge%ja(i) == 0).or.(disu_merge%ihc(i) == -1).or.&
        (disu_merge%cl12(i) < R8ZERO).or.(disu_merge%hwva(i) < R8ZERO)) then
        call errmsg('tQuads_merge_disu: program error 2.')
      end if
    end do
    !
    ! clean-up
    deallocate(wrk_ja, wrk_ihc, wrk_cl12, wrk_hwva, r4w, i4w)
    deallocate(nodes_arr, nodes_offset, lids2im, pos, iac_chk)
    !
    return
  end subroutine tQuads_merge_disu
    
! ==============================================================================
! ==============================================================================
! tMF6Exchange
! ==============================================================================
! ==============================================================================
  
  subroutine tMF6Exchange_init(this, ihc, cl1, cl2, hwva)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Exchange) :: this
    !
    integer(I4B), intent(in), optional :: ihc
    real(R8B), intent(in), optional :: cl1
    real(R8B), intent(in), optional :: cl2
    real(R8B), intent(in), optional :: hwva
    ! -- local
! ------------------------------------------------------------------------------
    call this%clean()
    !
    this%nexg = I4ZERO
    !
    this%ihc  = I4ZERO
    this%cl1  = R8ZERO
    this%cl2  = R8ZERO
    this%hwva = R8ZERO
    !
    if (present(ihc))  this%ihc  = ihc
    if (present(cl1))  this%cl1  = cl1
    if (present(cl2))  this%cl2  = cl2
    if (present(hwva)) this%hwva = hwva
    !
    return
  end subroutine tMF6Exchange_init
  
  subroutine tMF6Exchange_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Exchange) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (allocated(this%cellidm1)) deallocate(this%cellidm1)
    if (allocated(this%cellidm2)) deallocate(this%cellidm2)
    !
    return
  end subroutine tMF6Exchange_clean
  
  subroutine tMF6Exchange_set(this, vintf_from, vintf_to, lookup)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Exchange) :: this
    type(tVertIntf), pointer, intent(in) :: vintf_from
    type(tVertIntf), pointer, intent(in) :: vintf_to
    type(tLookupTable), pointer, intent(in) :: lookup
    ! -- local
    logical, parameter :: check_duplicates = .false.
    integer(I1B), dimension(:,:), allocatable :: i1wk
    integer(I4B) :: max_lay_act, ic, ir, il, jl, il_fr, il_to, n_cells
    integer(I4B) :: ip, i, j, iact, nconn, ios
    integer(I4B) :: n_fr, n_to, nr, nc, nl_fr, nl_to, stat
    integer(I4B), dimension(:), allocatable :: lay_fr, lay_to, lay_ptr
    integer(I4B), dimension(:,:), allocatable :: lay_connect, nod_fr, nod_to
    integer(I4B), dimension(:,:,:), allocatable :: i4wk
! ------------------------------------------------------------------------------
    !
    ! check
    if (vintf_from%n_cells /= vintf_from%n_cells) then
      call errmsg('tMF6Exchange_set: program error.')
    end if
    n_cells = vintf_from%n_cells
    !
    if (vintf_from%lm == vintf_to%lm) then
      max_lay_act = max(maxval(vintf_from%lay_act), maxval(vintf_to%lay_act))
      !
      allocate(i4wk(2, n_cells, max_lay_act))
      i4wk = 0
      !
      ! fill the from/M1 cells
      do i = 1, vintf_from%nlay_act
        il = vintf_from%lay_act(i)
        do j = 1, n_cells
          i4wk(1,j,il) = vintf_from%nod(j,i)
        end do
      end do
      ! fill the to/M2 cells
      do i = 1, vintf_to%nlay_act
        il = vintf_to%lay_act(i)
        do j = 1, n_cells
          i4wk(2,j,il) = vintf_to%nod(j,i)
        end do
      end do
      !
      ! set the connections
      do iact = 1, 2
        this%nexg = 0
        do i = 1, max_lay_act
          do j = 1, n_cells ! all interface cells including 9-point
            n_fr = i4wk(1,j,i); n_to = i4wk(2,j,i)
            if ((n_fr > 0).and.(n_to > 0)) then
              this%nexg = this%nexg + 1
              if (iact == 2) then
                this%cellidm1(this%nexg) = n_fr
                this%cellidm2(this%nexg) = n_to
              end if
            end if
          end do
        end do
        if (iact == 1) then
          if (this%nexg > 0) then
            allocate(this%cellidm1(this%nexg))
            allocate(this%cellidm2(this%nexg))
          end if
        end if
      end do
    else
      if (allocated(lay_fr)) deallocate(lay_fr) !from
      if (allocated(lay_to)) deallocate(lay_to) !to
      nconn = size(lookup%table,2)
      allocate(lay_fr(nconn), lay_to(nconn))
      !
      lay_fr = lookup%table(vintf_from%lm,:) !from
      lay_to = lookup%table(vintf_to%lm,:) !to
      nl_fr = maxval(lay_fr); nl_to = maxval(lay_to)
      !
      ! blank the inactive layers
      allocate(lay_ptr(nconn)) ! allocate a little larger
      !
      lay_ptr = 0
      do jl = 1, vintf_from%nlay_act
        il_fr = vintf_from%lay_act(jl)
        lay_ptr(il_fr) = 1
      end do
      do i = 1, nconn
        il_fr = lay_fr(i)
        if (lay_ptr(il_fr) == 0) then
          lay_fr(i) = 0
        end if
      end do
      !
      lay_ptr = 0
      do jl = 1, vintf_to%nlay_act
        il_to = vintf_to%lay_act(jl)
        lay_ptr(il_to) = 1
      end do
      do i = 1, nconn
        il_to = lay_to(i)
        if (lay_ptr(il_to) == 0) then
          lay_to(i) = 0
        end if
      end do
      !
      if (allocated(lay_connect)) deallocate(lay_connect)
      allocate(lay_connect(nl_to,nl_fr)); lay_connect = 0
      !
      do i = 1, nconn
        ic = lay_to(i); ir = lay_fr(i)
        if ((ic > 0).and.(ir > 0)) then
          lay_connect(ic,ir) = 1
        end if
      end do
      !
      allocate(nod_fr(n_cells, nl_fr)); nod_fr = 0
      allocate(nod_to(n_cells, nl_to)); nod_to = 0
      !
      ! fill the from cells
      do jl = 1, vintf_from%nlay_act
        il_fr = vintf_from%lay_act(jl)
        do ip = 1, n_cells
          nod_fr(ip,il_fr) = vintf_from%nod(ip,jl)
        end do
      end do
      !
      ! fill the to cells
      do jl = 1, vintf_to%nlay_act
        il_to = vintf_to%lay_act(jl)
        do ip = 1, n_cells
          nod_to(ip,il_to) = vintf_to%nod(ip,jl)
        end do
      end do
      !
      do iact = 1, 2
        this%nexg = 0
        do il_fr = 1, nl_fr
          do ip = 1, n_cells
            do il_to = 1, nl_to
              !
              ! layers can be connected
              if (lay_connect(il_to,il_fr) == 1) then
                ! get the node numbers
                n_fr = nod_fr(ip,il_fr); n_to = nod_to(ip,il_to)
                ! positive number numbers mean connected
                if ((n_fr > 0).and.(n_to > 0)) then
                  this%nexg = this%nexg + 1
                  if (iact == 2) then
                    this%cellidm1(this%nexg) = n_fr
                    this%cellidm2(this%nexg) = n_to
                  end if
                end if
              end if
            end do
          end do
        end do
        if (iact == 1) then
          if (this%nexg > 0) then
            allocate(this%cellidm1(this%nexg))
            allocate(this%cellidm2(this%nexg))
          end if
        end if
      end do
      !
      ! clean up
      if (allocated(lay_fr)) deallocate(lay_fr)
      if (allocated(lay_to)) deallocate(lay_to)
      if (allocated(lay_ptr)) deallocate(lay_ptr)
      if (allocated(nod_fr)) deallocate(nod_fr)
      if (allocated(nod_to)) deallocate(nod_to)
      if (allocated(lay_connect)) deallocate(lay_connect)
    end if
    !
    ! check for duplicates
    if ((this%nexg > 0).and.check_duplicates) then
      nr = maxval(this%cellidm1); nc = maxval(this%cellidm2)
      allocate(i1wk(nc,nr), stat=ios)
      if (ios == 0) then
        i1wk = 0
        do i = 1, this%nexg
          ir = this%cellidm1(i); ic = this%cellidm2(i)
          if (i1wk(ic,ir) == 0) then
            i1wk(ic,ir) = 1
          else
            call errmsg('tMF6Exchange_set: program error, duplicates found.')
          end if
        end do
        deallocate(i1wk)
      else
        call logmsg('tMF6Exchange_set: skipping check.')
      end if
    end if
    !
    return
  end subroutine tMF6Exchange_set
  
  subroutine tMF6Exchange_copy(this, tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Exchange) :: this
    type(tMF6Exchange), pointer, intent(inout) :: tgt
    ! -- local
! ------------------------------------------------------------------------------
    call tgt%clean()
    !
    tgt%nexg = this%nexg
    allocate(tgt%cellidm1, source=this%cellidm1)
    allocate(tgt%cellidm2, source=this%cellidm2)
    tgt%ihc  = this%ihc
    tgt%cl1  = this%cl1
    tgt%cl2  = this%cl2
    tgt%hwva =this%hwva
    !
    return
  end subroutine tMF6Exchange_copy
  
  subroutine tMF6Exchange_write(this, f_asc, id, wbd)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tMF6Exchange) :: this
    character(len=*), intent(in) :: f_asc
    character(len=*), intent(in) :: id
    type(tMf6Wbd), pointer, intent(inout) :: wbd
    ! -- local
    character(len=MXSLEN) :: s
    integer(I4B) :: i
    integer(I4B), dimension(:), allocatable :: ihc
    real(R8B), dimension(:), allocatable :: cl1, cl2, hwva
! ------------------------------------------------------------------------------
    allocate(ihc(this%nexg));   ihc = this%ihc
    allocate(cl1(this%nexg));   cl1 = this%cl1
    allocate(cl2(this%nexg));   cl2 = this%cl2
    allocate(hwva(this%nexg)); hwva = this%hwva
    !
    call wbd%write_list(id=id, &
      i4a1=this%cellidm1, i4a2=this%cellidm2, i4a3=ihc, &
      r8a1=cl1, r8a2=cl2, r8a3=hwva, f_asc=f_asc)
    !
    deallocate(ihc,cl1,cl2,hwva)
    !
    return
  end subroutine tMF6Exchange_write
  
! ==============================================================================
! ==============================================================================
! tNbrIntf 
! ==============================================================================
! ==============================================================================
  
  subroutine tNbrIntf_copy(this, tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tNbrIntf) :: this
    type(tNbrIntf), intent(inout), pointer :: tgt
    ! -- local
! ------------------------------------------------------------------------------
    !
    tgt%active        = this%active
    tgt%nbr_gid       = this%nbr_gid
    tgt%nbr_lid       = this%nbr_lid
    tgt%n_cells       = this%n_cells
    tgt%vintf_active  = this%vintf_active
    tgt%xch_active    = this%xch_active
    !
    if (allocated(tgt%lcell_nr)) deallocate(tgt%lcell_nr)
    allocate(tgt%lcell_nr, source=this%lcell_nr)
    !
    if (associated(this%vintf_from)) then
      allocate(tgt%vintf_from)
      call this%vintf_from%copy(tgt%vintf_from)
    end if
    if (associated(this%vintf_to)) then
      allocate(tgt%vintf_to)
      call this%vintf_to%copy(tgt%vintf_to)
    end if
    if (associated(this%xch)) then
      allocate(tgt%xch)
      call this%xch%copy(tgt%xch)
    end if
    !
    return
  end subroutine tNbrIntf_copy
  
  function tNbrIntf_includes_fp_stencil(this, nc, nr) result(included)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tNbrIntf) :: this
    integer(I4B), intent(in) :: nc
    integer(I4B), intent(in) :: nr
    logical :: included
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    included = .false.
    do i = 1, this%n_cells
      if (fp_dir(this%lcell_nr(:,i), nc, nr)) then
        included = .true.
        exit
      end if
    end do
    !
    return
  end function tNbrIntf_includes_fp_stencil
    
  function fp_dir(n, nc, nr) result(res)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), dimension(2), intent(in) :: n
    integer(I4B), intent(in) :: nc, nr
    logical :: res
    ! -- local
    integer(I4B) :: n0, n1, ic0, ir0, ic1, ir1, i_dir, idummy
! ------------------------------------------------------------------------------
    n0 = abs(n(1)); n1 = abs(n(2))
    !
    call node_to_icrl(n0, ic0, ir0, idummy, nc, nr)
    call node_to_icrl(n1, ic1, ir1, idummy, nc, nr)
    !
    i_dir = get_direction(ic0, ir0, ic1, ir1)
    if ((i_dir == i_e).or.(i_dir == i_w).or. &
        (i_dir == i_s).or.(i_dir == i_n)) then
      res = .true.
    else
      res = .false.
    end if 
    !
    return
  end function fp_dir
  
  subroutine tNbrIntf_append(this, lcell_nr)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tNbrIntf) :: this
    integer(I4B), dimension(:,:), intent(in) :: lcell_nr
    ! -- local
    integer(I4B) :: n, m, i
    integer(I4B), dimension(:,:), allocatable :: i4wk2d
! ------------------------------------------------------------------------------
    n = size(lcell_nr,2)
    allocate(i4wk2d(2,this%n_cells+n))
    !
    m = 0
    do i = 1, this%n_cells
      m = m + 1
      i4wk2d(1,m) = this%lcell_nr(1,i)
      i4wk2d(2,m) = this%lcell_nr(2,i)
    end do
    do i = 1, n
      m = m + 1
      i4wk2d(1,m) = lcell_nr(1,i)
      i4wk2d(2,m) = lcell_nr(2,i)
    end do
    !
    this%n_cells = this%n_cells + n
    deallocate(this%lcell_nr)
    allocate(this%lcell_nr(2,this%n_cells))
    !
    do i = 1, this%n_cells
      this%lcell_nr(1,i) = i4wk2d(1,i)
      this%lcell_nr(2,i) = i4wk2d(2,i)
    end do
    deallocate(i4wk2d)
    !
    return
  end subroutine tNbrIntf_append
  
  subroutine tNbrIntf_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tNbrIntf) :: this
    ! -- local
! ------------------------------------------------------------------------------
    this%nbr_gid      = I4ZERO
    this%nbr_lid      = I4ZERO
    this%n_cells      = I4ZERO
    this%active       = .true.
    this%vintf_active = I4ZERO
    this%xch_active   = I4ZERO
    !
    return
  end subroutine tNbrIntf_init
    
  subroutine tNbrIntf_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tNbrIntf) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (allocated(this%lcell_nr)) deallocate(this%lcell_nr)
    !
    if (associated(this%vintf_from)) then
      call this%vintf_from%clean()
      deallocate(this%vintf_from); this%vintf_from => null()
    end if
    if (associated(this%vintf_to)) then
      call this%vintf_to%clean()
      deallocate(this%vintf_to); this%vintf_to => null()
    end if
    if (associated(this%xch)) then
      call this%xch%clean()
      deallocate(this%xch); this%xch => null()
    end if
    !
    return
  end subroutine tNbrIntf_clean
  
! ==============================================================================
! ==============================================================================
! tVertIntf 
! ==============================================================================
! ==============================================================================
  
  subroutine tVertIntf_init(this, lm, nlay_act, n_cells)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVertIntf) :: this
    !
    integer(I4B), intent(in), optional :: lm
    integer(I4B), intent(in), optional :: nlay_act
    integer(I4B), intent(in), optional :: n_cells
    ! -- local
    integer(I4B) :: nl, nc
! ------------------------------------------------------------------------------
    !
    call this%clean()
    !
    this%lm = 0
    this%nlay_act = 0
    this%n_cells = 0
    this%mv = -99999.
    !
    if (present(lm)) then
      this%lm = lm
    end if
    if (present(nlay_act)) then
      this%nlay_act = nlay_act
    end if
    if (present(n_cells)) then
      this%n_cells = n_cells
    end if
    !
    if (this%nlay_act > 0) then
      nl = this%nlay_act
      allocate(this%lay_act(nl)); this%lay_act = 0
    end if
    if ((this%n_cells > 0).and.(this%nlay_act > 0)) then
      nl = this%nlay_act; nc = this%n_cells
      allocate(this%nod(nc,nl)); this%nod = 0
      allocate(this%zp(nc,nl+1)); this%zp = this%mv
    end if
    !
    return
  end subroutine tVertIntf_init
  
  subroutine tVertIntf_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVertIntf) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (allocated(this%lay_act)) deallocate(this%lay_act)
    if (allocated(this%nod))     deallocate(this%nod)
    if (allocated(this%zp))      deallocate(this%zp)
    !
    return
  end subroutine tVertIntf_clean
  
  subroutine tVertIntf_copy(this, tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tVertIntf) :: this
    type(tVertIntf), intent(inout), pointer :: tgt
    ! -- local
! ------------------------------------------------------------------------------
    tgt%lm       = this%lm
    tgt%nlay_act = this%nlay_act
    tgt%n_cells  = this%n_cells
    tgt%mv       = this%mv
    allocate(tgt%lay_act, source=this%lay_act)
    allocate(tgt%nod, source=this%nod)
    allocate(tgt%zp, source=this%zp)
    !
    return
  end subroutine tVertIntf_copy
  
! ==============================================================================
! ==============================================================================
! tIntf 
! ==============================================================================
! ==============================================================================
  
  subroutine tIntf_init(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    this%n_nbr  = I4ZERO
    this%mo_nbr = I4ZERO
    this%n_mv   = I4ZERO
    this%n_fill = I4ZERO
    this%ncol   = I4ZERO
    this%nrow   = I4ZERO
    this%nbr_intf => null()
    !
    return
  end subroutine tIntf_init
  
  subroutine tIntf_print(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    ! -- local
    type(tNbrIntf), pointer :: nbr => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    call logmsg('======================================================')
    call logmsg('n_nbr  : '//ta([this%n_nbr]))
    call logmsg('mo_nbr : '//ta([this%mo_nbr]))
    call logmsg('n_mv   : '//ta([this%n_mv]))
    call logmsg('ncol   : '//ta([this%ncol]))
    call logmsg('nrow   : '//ta([this%nrow]))
    do i = 1, this%n_nbr
      nbr => this%nbr_intf(i)
      call logmsg('------------------------------------------------------')
      call logmsg('neighbor '//ta([i])//':')
      call logmsg(' nbr_gid '//ta([nbr%nbr_gid]))
      call logmsg(' nbr_lid '//ta([nbr%nbr_lid]))
      call logmsg(' n_cells '//ta([nbr%n_cells]))
    end do
    call logmsg('======================================================')
    !
    return
  end subroutine tIntf_print
  
  subroutine tIntf_copy(this, tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    type(tIntf), intent(inout), pointer :: tgt
    ! -- local
    type(tNbrIntf), pointer :: nbr => null(), nbr_tgt => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    tgt%n_nbr  = this%n_nbr
    tgt%mo_nbr = this%mo_nbr
    tgt%n_mv   = this%n_mv
    tgt%n_fill = this%n_fill
    tgt%ncol   = this%ncol
    tgt%nrow   = this%nrow
    !
    if (associated(tgt%nbr_intf)) deallocate(tgt%nbr_intf)
    allocate(tgt%nbr_intf(tgt%n_nbr))
    !
    do i = 1, tgt%n_nbr
      nbr => this%nbr_intf(i)
      nbr_tgt => tgt%nbr_intf(i)
      call nbr%copy(nbr_tgt)
    end do
    !
    if (this%n_mv > 0) then
      allocate(tgt%mv_lcell_nr, source=this%mv_lcell_nr)
    end if
    !
    if (this%n_fill > 0) then
      allocate(tgt%fill_lcell_nr, source=this%fill_lcell_nr)
    end if
    !
    return
  end subroutine tIntf_copy
  
  subroutine tIntf_get_mask(this, mask, bbx, gid)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    integer(I4B), dimension(:,:), intent(inout), allocatable :: mask
    type(tBbx), intent(in), optional :: bbx
    integer(I4B), intent(in), optional :: gid
    ! -- local
    logical, parameter :: debug = .true.
    type(tNbrIntf), pointer :: nbr => null()
    logical :: fill
    character(len=MXSLEN) :: fp
    integer(I4B), dimension(:,:), pointer :: lcell_nr => null()
    integer(I4B) :: i, inbr, jr, jc, ir, ic, ir0, ic0, ir1, ic1
    integer(I4B) :: n, n0, n1, icell, idummy, ist
    integer(I4B), dimension(:,:), allocatable :: wrk
! ------------------------------------------------------------------------------
    !
    if (allocated(mask)) then
      if (debug) then
        allocate(wrk(this%ncol,this%nrow))
        wrk = mask
      end if
      deallocate(mask)
    end if
    allocate(mask(this%ncol,this%nrow))
    mask = 0
    !
    ! neighbor contribution
    do inbr = 1, this%n_nbr
      nbr => this%nbr_intf(inbr)
      lcell_nr => nbr%lcell_nr
      do icell = 1, nbr%n_cells
        n0 = lcell_nr(1,icell); n1 = lcell_nr(2,icell)
        call node_to_icrl(abs(n0), ic0, ir0, idummy, this%ncol, this%nrow)
        call node_to_icrl(n1,      ic1, ir1, idummy, this%ncol, this%nrow)
        mask(ic0,ir0) = 1
      end do
    end do
    ! missing value contribution
    if (this%n_mv > 0) then
      do icell = 1, this%n_mv
        n0 = this%mv_lcell_nr(1,icell)
        n1 = this%mv_lcell_nr(2,icell)
        call node_to_icrl(abs(n0), ic0, ir0, idummy, this%ncol, this%nrow)
        call node_to_icrl(n1, ic1, ir1, idummy, this%ncol, this%nrow)
        mask(ic0,ir0) = 1
      end do
    end if
    !
    ! neighbor contribution
    do inbr = 1, this%n_nbr
      nbr => this%nbr_intf(inbr)
      lcell_nr => nbr%lcell_nr
      do icell = 1, nbr%n_cells
        n0 = lcell_nr(1,icell); n1 = lcell_nr(2,icell)
        call node_to_icrl(abs(n0), ic0, ir0, idummy, this%ncol, this%nrow)
        call node_to_icrl(n1,      ic1, ir1, idummy, this%ncol, this%nrow)
        !
        if (n0 > 0) then
          call get_upstream_icir(ic0, ir0, ic1, ir1, ic, ir)
          if (valid_icir(ic, ir, this%ncol, this%nrow)) then
            if (mask(ic,ir) /= 1) then
              mask(ic,ir) = 2
            end if
          end if
        end if
      end do
    end do
    !
    ! missing value contribution
    if (this%n_mv > 0) then
      do icell = 1, this%n_mv
        n0 = this%mv_lcell_nr(1,icell); n1 = this%mv_lcell_nr(2,icell)
        call node_to_icrl(abs(n0), ic0, ir0, idummy, this%ncol, this%nrow)
        call node_to_icrl(n1, ic1, ir1, idummy, this%ncol, this%nrow)
        !
        if (n0 > 0) then
          call get_upstream_icir(ic0, ir0, ic1, ir1, ic, ir)
          if (valid_icir(ic, ir, this%ncol, this%nrow)) then
            if (mask(ic,ir) /= 1) mask(ic,ir) = 2
          end if
        end if
      end do
    end if
    !
    if (present(gid)) then
      fp = 'e:\LHM\selection\simulations\run_output\mask'
      call writeflt(fp, mask, this%ncol, this%nrow, bbx%xll, bbx%yll, bbx%cs, 0)
    end if
    !
    do while(.true.)
      n = 0
      do ir = 1, this%nrow
        do ic = 1, this%ncol
          if (mask(ic,ir) == 2) then
            ! get 9-point stencil
            sticir_np = reshape([ic, ir, ic+1, ir, ic-1, ir, ic, ir+1, ic, ir-1, &
              ic+1, ir-1, ic-1, ir-1, ic+1, ir+1, ic-1, ir+1], shape(sticir_np))
            !
            ! loop over the neighbors
            do ist = 2, nst_np
              jc = sticir_np(1,ist); jr = sticir_np(2,ist)
              if (valid_icir(jc, jr, this%ncol, this%nrow)) then
                if (mask(jc,jr) == 0) then
                  mask(jc,jr) = 2; n = n + 1
                end if
              end if
            end do
            !
            mask(ic,ir) = 1
          end if
        end do
      end do
      
      !fp = 'e:\LHM\selection\simulations\run_output\mask_1'
      !call writeflt(fp, mask, this%ncol, this%nrow, bbx%xll, bbx%yll, bbx%cs, 0)
      !stop
      
      if (n == 0) then
        exit
      end if
    end do
    !
    ! add repair values
    do i = 1, this%n_fill
      n = this%fill_lcell_nr(i)
      call node_to_icrl(n, ic, ir, idummy, this%ncol, this%nrow)
      mask(ic,ir) = 1
    end do
    !
    if (present(gid)) then
      fp = 'e:\LHM\selection\simulations\run_output\mask_final'
      call writeflt(fp, mask, this%ncol, this%nrow, bbx%xll, bbx%yll, bbx%cs, 0)
      stop
    end if
    !
    ! check
    if (allocated(wrk)) then
      n = 0
      do ir = 1,  this%nrow
        do ic = 1, this%ncol
          if (mask(ic,ir) /= wrk(ic,ir)) then
            n = n + 1
          end if
        end do
      end do
      if (n > 0) then
        call errmsg('tIntf_get_mask program error.')
      end if
      deallocate(wrk)
    end if
    !
    return
  end subroutine tIntf_get_mask
    
  subroutine tIntf_write(this, iu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    integer(I4B), intent(in) :: iu
    ! -- local
    type(tNbrIntf), pointer :: nbr => null()
    type(tVertIntf), pointer :: vintf => null()
    type(tMF6Exchange), pointer :: xch => null()
    integer(I4B) :: i, j
! ------------------------------------------------------------------------------
    !
    write(iu) this%n_nbr, this%mo_nbr, this%n_mv, this%n_fill, &
      this%ncol, this%nrow
    do i = 1, this%n_nbr
      nbr => this%nbr_intf(i)
      write(iu) nbr%nbr_gid, nbr%nbr_lid, nbr%n_cells, nbr%vintf_active, &
        nbr%xch_active
      write(iu) nbr%lcell_nr
      if (nbr%vintf_active == 1) then
        do j = 1, 2
          if (j == 1) then
            vintf => nbr%vintf_from
          else
            vintf => nbr%vintf_to
          end if
          write(iu) vintf%lm
          write(iu) vintf%nlay_act
          write(iu) vintf%n_cells
          write(iu) vintf%mv
          write(iu) vintf%lay_act
          write(iu) vintf%nod
          write(iu) vintf%zp
        end do
      end if
      if (nbr%xch_active == 1) then
        xch => nbr%xch
        write(iu) xch%nexg
        write(iu) xch%cellidm1
        write(iu) xch%cellidm2
        write(iu) xch%ihc
        write(iu) xch%cl1
        write(iu) xch%cl2
        write(iu) xch%hwva
      end if
    end do
    !
    if (this%n_mv > 0) then
      write(iu) this%mv_lcell_nr
    end if
    !
    if (this%n_fill > 0) then
      write(iu) this%fill_lcell_nr
    end if
    !
    return
  end subroutine tIntf_write
  
  subroutine tIntf_read(this, iu, pos, read_vintf)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    integer(I4B), intent(in) :: iu
    integer(I8B), intent(in), optional :: pos
    logical, intent(in), optional :: read_vintf
    ! -- local
    type(tNbrIntf), pointer :: nbr => null()
    type(tVertIntf), pointer :: vintf => null()
    type(tMF6Exchange), pointer :: xch => null()
    logical :: read_vintf_loc
    integer(I8B) :: p
    integer(I4B) :: i, j
! ------------------------------------------------------------------------------
    if(.not.present(pos)) then
      read(iu) this%n_nbr, this%mo_nbr, this%n_mv, this%n_fill, &
        this%ncol, this%nrow
    else
      p = pos
      read(iu, pos=p) this%n_nbr, this%mo_nbr, this%n_mv, this%n_fill, &
        this%ncol, this%nrow
      p = p + sizeof(this%n_nbr) + sizeof(this%mo_nbr) + sizeof(this%n_mv) + &
              sizeof(this%n_fill) + sizeof(this%ncol)  + sizeof(this%nrow)
    end if
    if (present(read_vintf)) then
      read_vintf_loc = read_vintf
    else
      read_vintf_loc = .true.
    end if
    !
    if (this%n_nbr <= 0) then
      !call logmsg('tIntf_read: n_nbr <= 0.')
    else
      if (associated(this%nbr_intf)) deallocate(this%nbr_intf)
      allocate(this%nbr_intf(this%n_nbr))
    end if
    !
    do i = 1, this%n_nbr
      nbr => this%nbr_intf(i)
      call nbr%init()
      if(.not.present(pos)) then
        read(iu) nbr%nbr_gid, nbr%nbr_lid, nbr%n_cells, nbr%vintf_active, &
          nbr%xch_active
      else
        read(iu, pos=p) nbr%nbr_gid, nbr%nbr_lid, nbr%n_cells, &
          nbr%vintf_active, nbr%xch_active
        p = p + sizeof(nbr%nbr_gid) + sizeof(nbr%nbr_lid) + &
          sizeof(nbr%n_cells) + sizeof(nbr%vintf_active) + &
          sizeof(nbr%vintf_active)
      end if
      !
      if (nbr%n_cells <= 0) then
        call errmsg('tIntf_read: n_cells <= 0.')
      else
        if (allocated(nbr%lcell_nr)) deallocate(nbr%lcell_nr)
        allocate(nbr%lcell_nr(2,nbr%n_cells))
      end if
      !
      if(.not.present(pos)) then
        read(iu) nbr%lcell_nr
      else
        read(iu,pos=p) nbr%lcell_nr
        p = p + sizeof(nbr%lcell_nr)
      end if
      !
      ! check for zero values
      if (minval(abs(nbr%lcell_nr)) == 0) then
        call errmsg('tIntf_read: zero cell numbers found.')
      end if
      !
      if (nbr%vintf_active == 1) then
        do j = 1, 2
          if (j == 1) then
            if (associated(nbr%vintf_from)) then
              call nbr%vintf_from%clean()
              deallocate(nbr%vintf_from)
            end if
            allocate(nbr%vintf_from); vintf => nbr%vintf_from
          else
            if (associated(nbr%vintf_to)) then
              call nbr%vintf_to%clean()
              deallocate(nbr%vintf_to)
            end if
            allocate(nbr%vintf_to); vintf => nbr%vintf_to
          end if
          !
          if(.not.present(pos)) then
            read(iu) vintf%lm
            read(iu) vintf%nlay_act
            read(iu) vintf%n_cells
            read(iu) vintf%mv
          else
            read(iu,pos=p) vintf%lm
            p = p + sizeof(vintf%lm)
            read(iu,pos=p) vintf%nlay_act
            p = p + sizeof(vintf%nlay_act)
            read(iu,pos=p) vintf%n_cells
            p = p + sizeof(vintf%n_cells)
            read(iu,pos=p) vintf%mv
            p = p + sizeof(vintf%mv)
          end if
          !
          if (.not.read_vintf) then
            if (vintf%nlay_act > 0) p = p + I4B * vintf%nlay_act !vintf%lay_act
            if (vintf%n_cells > 0) then
              p = p + I4B * vintf%n_cells*vintf%nlay_act !vintf%nod
              p = p + R4B * vintf%n_cells*(vintf%nlay_act+1) !vintf%zp
            end if
          else
             if (vintf%nlay_act > 0) then
               allocate(vintf%lay_act(vintf%nlay_act))
             end if
             if (vintf%n_cells > 0) then
               allocate(vintf%nod(vintf%n_cells,vintf%nlay_act))
               allocate(vintf%zp(vintf%n_cells,vintf%nlay_act+1))
             end if
             if(.not.present(pos)) then
               read(iu) vintf%lay_act
               read(iu) vintf%nod
               read(iu) vintf%zp
             else
               read(iu,pos=p) vintf%lay_act
               p = p + sizeof(vintf%lay_act)
               read(iu,pos=p) vintf%nod
               p = p + sizeof(vintf%nod)
               read(iu,pos=p) vintf%zp
               p = p + sizeof(vintf%zp)
             end if
          end if
        end do
        !
        if (.not.read_vintf) then
          call nbr%vintf_from%clean()
          call nbr%vintf_to%clean()
          deallocate(nbr%vintf_from, nbr%vintf_to)
          nbr%vintf_from => null()
          nbr%vintf_to   => null()
          nbr%vintf_active = 0
        end if
      end if
      if (nbr%xch_active == 1) then
        if (associated(nbr%xch)) then
          call nbr%xch%clean(); deallocate(nbr%xch)
        end if
        allocate(nbr%xch); xch => nbr%xch
        !
        if(.not.present(pos)) then
          read(iu) xch%nexg
        else
          read(iu,pos=p) xch%nexg
          p = p + sizeof(xch%nexg)
        end if
        !
        allocate(xch%cellidm1(xch%nexg))
        allocate(xch%cellidm2(xch%nexg))
        !
        if(.not.present(pos)) then
          read(iu) xch%cellidm1
          read(iu) xch%cellidm2
          read(iu) xch%ihc
          read(iu) xch%cl1
          read(iu) xch%cl2
          read(iu) xch%hwva
        else
          read(iu,pos=p) xch%cellidm1
          p = p + sizeof(xch%cellidm1)
          read(iu,pos=p) xch%cellidm2
          p = p + sizeof(xch%cellidm2)
          read(iu,pos=p) xch%ihc
          p = p + sizeof(xch%ihc)
          read(iu,pos=p) xch%cl1
          p = p + sizeof(xch%cl1)
          read(iu,pos=p) xch%cl2
          p = p + sizeof(xch%cl2)
          read(iu,pos=p) xch%hwva
          p = p + sizeof(xch%hwva)
        end if
      end if
    end do
    !
    if (this%n_mv < 0) then
      call errmsg('tIntf_read: n_mv < 0.')
    end if
    if (this%n_mv > 0) then
      if (allocated(this%mv_lcell_nr)) deallocate(this%mv_lcell_nr)
      allocate(this%mv_lcell_nr(2,this%n_mv))
    end if
    if (this%n_mv > 0) then
      if(.not.present(pos)) then
        read(iu) this%mv_lcell_nr
      else
        read(iu,pos=p) this%mv_lcell_nr
        p = p + sizeof(this%mv_lcell_nr)
      end if
    end if
    !
    if (this%n_fill < 0) then
      call errmsg('tIntf_read: n_fill < 0.')
    end if
    if (this%n_fill > 0) then
      if (allocated(this%fill_lcell_nr)) deallocate(this%fill_lcell_nr)
      allocate(this%fill_lcell_nr(this%n_fill))
    end if
    if (this%n_fill > 0) then
      if(.not.present(pos)) then
        read(iu) this%fill_lcell_nr
      else
        read(iu,pos=p) this%fill_lcell_nr
        p = p + sizeof(this%fill_lcell_nr)
      end if
    end if
    !
    return
  end subroutine tIntf_read

  subroutine tIntf_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    ! -- local
    type(tNbrIntf), pointer :: nbr => null()
    integer(I4B) :: inbr
! ------------------------------------------------------------------------------
    !
    if (associated(this%nbr_intf)) then
      do inbr = 1, this%n_nbr
        nbr => this%nbr_intf(inbr)
        call nbr%clean()
      end do
      deallocate(this%nbr_intf); this%nbr_intf => null()
    end if
    !
    if (allocated(this%mv_lcell_nr)) then
      deallocate(this%mv_lcell_nr)
    end if
    if (allocated(this%fill_lcell_nr)) then
      deallocate(this%fill_lcell_nr)
    end if
    !
    call this%init()
    !
    return
  end subroutine tIntf_clean
  
  subroutine tIntf_replace(this, nbr_gid, nbr_tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    integer(I4B), intent(in) :: nbr_gid
    type(tNbrIntf), pointer, intent(in) :: nbr_tgt
    ! -- local
    type(tNbrIntf), pointer :: nbr => null()
    type(tNbrIntf), pointer :: nbr_src => null()
    logical :: found
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    found = .false.
    do i = 1, this%n_nbr
      nbr => this%nbr_intf(i)
      if (nbr%nbr_gid == nbr_gid) then
        found = .true.
        nbr_src => this%nbr_intf(i)
        exit
      end if
    end do
    !
    if (.not.found) then
      call logmsg('tIntf_replace:: nothing to be removed (id not found).')
      return
    end if
    !
    call nbr_src%clean(); call nbr_src%init()
    call nbr_tgt%copy(nbr_src)
    !
    return
  end subroutine tIntf_replace
  
  subroutine tIntf_remove(this, nbr_gid)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    integer(I4B), intent(in) :: nbr_gid
    ! -- local
    logical :: found
    integer(I4B) :: i, n
    type(tNbrIntf), pointer :: nbr => null()
    type(tNbrIntf), pointer :: nbr_src => null(), nbr_tgt => null()
    type(tNbrIntf), dimension(:), pointer :: wrk => null()
! ------------------------------------------------------------------------------
    !
    if (this%n_nbr <= 1) then
      call logmsg('tIntf_remove: nothing to be removed (n_nbr <= 1)')
      return
    end if
    !
    found = .false.
    do i = 1, this%n_nbr
      nbr => this%nbr_intf(i)
      if (nbr%nbr_gid == nbr_gid) then
        found = .true.; exit
      end if
    end do
    !
    if (.not.found) then
      call logmsg('tIntf_remove:: nothing to be removed (id not found).')
      return
    end if
    !
    allocate(wrk(this%n_nbr-1))
    n = 0
    do i = 1, this%n_nbr
      nbr_src => this%nbr_intf(i)
      if (nbr_src%nbr_gid /= nbr_gid) then
        n = n + 1
        nbr_tgt => wrk(n)
        call nbr_src%copy(nbr_tgt)
      end if
    end do
    !
    do i = 1, this%n_nbr
      nbr => this%nbr_intf(i)
      call nbr%clean()
    end do
    deallocate(this%nbr_intf)
    !
    ! copy the work array
    this%n_nbr = this%n_nbr - 1
    n = 0
    allocate(this%nbr_intf(this%n_nbr))
    do i = 1, this%n_nbr
      nbr_src => wrk(i)
      nbr_tgt => this%nbr_intf(i)
      call nbr_src%copy(nbr_tgt)
      !
      if (nbr_tgt%n_cells > n) then
        n = nbr_tgt%n_cells
        this%mo_nbr = nbr_tgt%nbr_gid
      end if
    end do
    !
    ! clean up
    do i = 1, this%n_nbr
      nbr => wrk(i)
      call nbr%clean()
    end do
    deallocate(wrk)
    !
    return
  end subroutine tIntf_remove
  
  function tIntf_get_nbytes(this) result(nbytes)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    integer(I8B) :: nbytes
    ! -- local
    type(tNbrIntf), pointer :: nbr => null()
    type(tVertIntf), pointer :: vintf => null()
    type(tMF6Exchange), pointer :: xch => null()
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    nbytes = 0
    nbytes = nbytes + sizeof(this%n_nbr)
    nbytes = nbytes + sizeof(this%mo_nbr)
    nbytes = nbytes + sizeof(this%n_mv)
    nbytes = nbytes + sizeof(this%n_fill)
    nbytes = nbytes + sizeof(this%ncol)
    nbytes = nbytes + sizeof(this%nrow)
    !
    do i = 1, this%n_nbr
      nbr => this%nbr_intf(i)
      
      nbytes = nbytes + sizeof(nbr%nbr_gid)
      nbytes = nbytes + sizeof(nbr%nbr_lid)
      nbytes = nbytes + sizeof(nbr%n_cells)
      nbytes = nbytes + sizeof(nbr%vintf_active)
      nbytes = nbytes + sizeof(nbr%xch_active)
      !
      nbytes = nbytes + sizeof(nbr%lcell_nr)
      !
      if (nbr%vintf_active == 1) then
        vintf => nbr%vintf_from
        nbytes = nbytes + sizeof(vintf%lm)
        nbytes = nbytes + sizeof(vintf%nlay_act)
        nbytes = nbytes + sizeof(vintf%n_cells)
        nbytes = nbytes + sizeof(vintf%mv)
        nbytes = nbytes + sizeof(vintf%lay_act)
        nbytes = nbytes + sizeof(vintf%nod)
        nbytes = nbytes + sizeof(vintf%zp)
        !
        vintf => nbr%vintf_to
        nbytes = nbytes + sizeof(vintf%lm)
        nbytes = nbytes + sizeof(vintf%nlay_act)
        nbytes = nbytes + sizeof(vintf%n_cells)
        nbytes = nbytes + sizeof(vintf%mv)
        nbytes = nbytes + sizeof(vintf%lay_act)
        nbytes = nbytes + sizeof(vintf%nod)
        nbytes = nbytes + sizeof(vintf%zp)
      end if
      !
      if (nbr%xch_active == 1) then
        xch => nbr%xch
        nbytes = nbytes + sizeof(xch%nexg)
        nbytes = nbytes + sizeof(xch%cellidm1)
        nbytes = nbytes + sizeof(xch%cellidm2)
        nbytes = nbytes + sizeof(xch%ihc)
        nbytes = nbytes + sizeof(xch%cl1)
        nbytes = nbytes + sizeof(xch%cl2)
        nbytes = nbytes + sizeof(xch%hwva)
      end if
    end do
    !
    if (this%n_mv > 0) then
      nbytes = nbytes + sizeof(this%mv_lcell_nr)
    end if
    !
    if (this%n_fill > 0) then
      nbytes = nbytes + sizeof(this%fill_lcell_nr)
    end if
    !
    !call logmsg('Wrong: '//ta([sizeof(this)])//' correct: '//ta([nbytes]))
    !
    return
  end function tIntf_get_nbytes
  
  subroutine tIntf_get_nbr_intf(this, nbr_res, nbr_gid)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tIntf) :: this
    integer(I4B), intent(in) :: nbr_gid
    type(tNbrIntf), pointer, intent(inout) :: nbr_res
    ! -- local
    integer(I4B) :: i
    type(tNbrIntf), pointer :: nbr => null()
! ------------------------------------------------------------------------------
    ! pointer should NOT be allocated
    if (associated(nbr_res)) then
      call nbr_res%clean(); 
      deallocate(nbr_res)
      nbr_res => null()
    end if
    !
    do i = 1, this%n_nbr
      nbr => this%nbr_intf(i)
      if (nbr%nbr_gid == nbr_gid) then
        allocate(nbr_res)
        nbr_res => nbr
        exit
      end if
    end do
    !
    return
  end subroutine tIntf_get_nbr_intf
  
  function get_direction(ic0, ir0, ic1, ir1) result(i_dir)
! ******************************************************************************
!
!    SPECIFICATIONS:
!
!   From --> to direction
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), intent(in) :: ic0
    integer(I4B), intent(in) :: ir0
    integer(I4B), intent(in) :: ic1
    integer(I4B), intent(in) :: ir1
    integer(I4B) :: i_dir
! ------------------------------------------------------------------------------
    i_dir = 0
    if (ir0 == ir1) then ! west or east
      if (ic0 < ic1) then
        i_dir = i_e; return
      else
        i_dir = i_w; return
      end if
    end if
    if (ic0 == ic1) then ! south or north
      if (ir0 < ir1) then
        i_dir = i_s; return
      else
        i_dir = i_n; return
      end if
    end if
    if ((ic0 < ic1).and.(ir0 > ir1)) then ! north-east
       i_dir = i_ne; return
    end if
    if ((ic0 > ic1).and.(ir0 > ir1)) then ! north-west
       i_dir = i_nw; return
    end if
    if ((ic0 < ic1).and.(ir0 < ir1)) then ! south-east
       i_dir = i_se; return
    end if
    if ((ic0 > ic1).and.(ir0 < ir1)) then ! south-west
       i_dir = i_sw; return
    end if
    !
    call errmsg('get_direction: program error.')
    !
    return
  end function get_direction
  
  function get_opposite_direction(i_dir_in) result(i_dir)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B) :: i_dir_in
    integer(I4B) :: i_dir
! ------------------------------------------------------------------------------
    select case(i_dir_in)
    case(i_e)
      i_dir = i_w
    case(i_w)
      i_dir = i_e
    case(i_s)
      i_dir = i_n
    case(i_n)
      i_dir = i_s
    case(i_ne)
      i_dir = i_sw
    case(i_nw)
      i_dir = i_se
    case(i_se)
      i_dir = i_nw
    case(i_sw)
      i_dir = i_ne
    case default
      call errmsg('get_opposite_direction: invalid direction.')
    end select
    !
    return
  end function get_opposite_direction
  
  subroutine get_upstream_icir(ic0, ir0, ic1, ir1, ic, ir)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), intent(in) :: ic0
    integer(I4B), intent(in) :: ir0
    integer(I4B), intent(in) :: ic1
    integer(I4B), intent(in) :: ir1
    integer(I4B), intent(out) :: ic
    integer(I4B), intent(out) :: ir
    ! -- local
    integer(I4B) :: i_dir
! ------------------------------------------------------------------------------
    ! get the direction
    i_dir = get_direction(ic0, ir0, ic1, ir1)
    i_dir =  get_opposite_direction(i_dir)
    !
    ! get 9-point stencil
    sticir_np = reshape([ic0, ir0, ic0+1, ir0, ic0-1, ir0, ic0, ir0+1, ic0, ir0-1, &
            ic0+1, ir0-1, ic0-1, ir0-1, ic0+1, ir0+1, ic0-1, ir0+1], shape(sticir_np))
     !
    ic = sticir_np(1,i_dir); ir = sticir_np(2,i_dir)
    !
    return
  end subroutine get_upstream_icir
 
  function valid_icir(ic, ir, nc, nr) result(valid)
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
  end function valid_icir
  
  function has_upstream_val(x, x_src, ic0, ir0, ic1, ir1) result(found)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), dimension(:,:), intent(in) :: x
    integer(I4B), intent(in) :: x_src
    integer(I4B), intent(in) :: ic0
    integer(I4B), intent(in) :: ir0
    integer(I4B), intent(in) :: ic1
    integer(I4B), intent(in) :: ir1
    logical :: found
    ! -- local
    integer(I4B) :: nc, nr, ic, ir, i_dir
! ------------------------------------------------------------------------------
    !
    found = .false.
    nc = size(x,1); nr = size(x,2)
    call get_upstream_icir(ic0, ir0, ic1, ir1, ic, ir)
    !
    if (valid_icir(ic,ir,nc,nr)) then
      if (x(ic,ir) == x_src) found = .true.
    end if
    !
    return
  end function has_upstream_val
  
  subroutine get_bnd_i4_5p(x, x_src, x_mv, int_icir, ext_icir, ext_v, &
    int_mv_icir, ext_mv_icir)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), dimension(:,:), intent(in) :: x
    integer(I4B), intent(in) :: x_src
    integer(I4B), intent(in) :: x_mv
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: int_icir
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: ext_icir
    integer(I4B), dimension(:),   allocatable, intent(inout), optional :: ext_v
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: int_mv_icir
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: ext_mv_icir
    ! -- local
    integer(I4B), dimension(:,:), allocatable :: int_icir_wk
    integer(I4B), dimension(:,:), allocatable :: ext_icir_wk
    integer(I4B), dimension(:,:), allocatable :: int_mv_icir_wk
    integer(I4B), dimension(:,:), allocatable :: ext_mv_icir_wk
    integer(I4B), dimension(:),   allocatable :: ext_v_wk
    logical :: skip
    integer(I4B) :: v
    integer(I4B) :: i, n, m, ist, ic, ir, jc, jr, kc, kr, nc, nr, n_max
! ------------------------------------------------------------------------------
    !
    nc = size(x,1); nr = size(x,2)
    n_max = 5*nc*nr
    allocate(int_icir_wk(2,n_max), ext_icir_wk(2,n_max), ext_v_wk(n_max))
    allocate(int_mv_icir_wk(2,n_max), ext_mv_icir_wk(2,n_max))
    !
    n = 0; m = 0
    do ir = 1, nr
      do ic = 1, nc
        if (x(ic,ir) == x_src) then
          !
          ! get 5-point stencil
          sticir_fp = reshape([ic, ir, ic+1, ir, ic-1, ir, ic, ir+1, ic, ir-1], shape(sticir_fp))
          !
          ! loop over the neighbors
          do ist = 2, nst_fp
            jc = sticir_fp(1,ist); jr = sticir_fp(2,ist)
            if ((jc >= 1).and.(jc <= nc).and.(jr >= 1).and.(jr <= nr)) then
              v = x(jc,jr)
              if (v /= x_src) then
                if (v /= x_mv) then
                  n = n + 1
                  if (n > n_max) then
                    call errmsg('tIntf_get_bnd_i4: program error.')
                  end if
                  int_icir_wk(1,n) = ic; int_icir_wk(2,n) = ir
                  ext_icir_wk(1,n) = jc; ext_icir_wk(2,n) = jr
                  ext_v_wk(n) = v
                else
                  m = m + 1
                  int_mv_icir_wk(1,m) = ic; int_mv_icir_wk(2,m) = ir
                  ext_mv_icir_wk(1,m) = jc; ext_mv_icir_wk(2,m) = jr
                end if
              end if
            end if
          end do
        end if
      end do
    end do
    !
    if (present(int_icir)) then
      if (allocated(int_icir)) deallocate(int_icir)
      if (n > 0) then
        allocate(int_icir(2,n))
        do i = 1, n
          int_icir(1,i) = int_icir_wk(1,i)
          int_icir(2,i) = int_icir_wk(2,i)
        end do
      end if
    end if
    if (present(ext_icir)) then
      if (allocated(ext_icir)) deallocate(ext_icir)
      if (n > 0) then
        allocate(ext_icir(2,n))
        do i = 1, n
          ext_icir(1,i) = ext_icir_wk(1,i)
          ext_icir(2,i) = ext_icir_wk(2,i)
        end do
      end if
    end if
    if (present(ext_v)) then
      if (allocated(ext_v)) deallocate(ext_v)
      if (n > 0) then
        allocate(ext_v(n))
        do i = 1, n
          ext_v(i) = ext_v_wk(i)
        end do
      end if
    end if
    if (present(int_mv_icir)) then
      if (allocated(int_mv_icir)) deallocate(int_mv_icir)
      if (m > 0) then
        allocate(int_mv_icir(2,m))
        do i = 1, m
          int_mv_icir(1,i) = int_mv_icir_wk(1,i)
          int_mv_icir(2,i) = int_mv_icir_wk(2,i)
        end do
      end if
    end if
    if (present(ext_mv_icir)) then
      if (allocated(ext_mv_icir)) deallocate(ext_mv_icir)
      if (m > 0) then
        allocate(ext_mv_icir(2,m))
        do i = 1, m
          ext_mv_icir(1,i) = ext_mv_icir_wk(1,i)
          ext_mv_icir(2,i) = ext_mv_icir_wk(2,i)
        end do
      end if
    end if
    !
    ! cleanup
    deallocate(int_icir_wk, ext_icir_wk, ext_v_wk)
    deallocate(int_mv_icir_wk, ext_mv_icir_wk)
    !
    return
  end subroutine get_bnd_i4_5p
    
  subroutine get_bnd_i4_9p(x, x_src, x_mv, int_icir, ext_icir, ext_v, &
    int_mv_icir, ext_mv_icir)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), dimension(:,:), intent(in) :: x
    integer(I4B), intent(in) :: x_src
    integer(I4B), intent(in) :: x_mv
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: int_icir
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: ext_icir
    integer(I4B), dimension(:),   allocatable, intent(inout), optional :: ext_v
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: int_mv_icir
    integer(I4B), dimension(:,:), allocatable, intent(inout), optional :: ext_mv_icir
    ! -- local
    integer(I4B), dimension(:,:), allocatable :: int_icir_wk
    integer(I4B), dimension(:,:), allocatable :: ext_icir_wk
    integer(I4B), dimension(:,:), allocatable :: int_mv_icir_wk
    integer(I4B), dimension(:,:), allocatable :: ext_mv_icir_wk
    integer(I4B), dimension(:),   allocatable :: ext_v_wk
    logical :: skip
    integer(I4B) :: v
    integer(I4B) :: i, n, m, ist, ic, ir, jc, jr, kc, kr, nc, nr, n_max
! ------------------------------------------------------------------------------
    !
    nc = size(x,1); nr = size(x,2)
    n_max = 9*nc*nr
    allocate(int_icir_wk(2,n_max), ext_icir_wk(2,n_max), ext_v_wk(n_max))
    allocate(int_mv_icir_wk(2,n_max), ext_mv_icir_wk(2,n_max))
    !
    n = 0; m = 0
    do ir = 1, nr
      do ic = 1, nc
        if (x(ic,ir) == x_src) then
          !
          ! get 9-point stencil
          sticir_np = reshape([ic,   ir,   & !i_p
                               ic+1, ir,   & !i_e
                               ic-1, ir,   & !i_w
                               ic,   ir+1, & !i_s
                               ic,   ir-1, & !i_n
                               ic+1, ir-1, & !i_ne
                               ic-1, ir-1, & !i_nw
                               ic+1, ir+1, & !i_sw
                               ic-1, ir+1], &
            shape(sticir_np))
          !
          ! loop over the neighbors
          do ist = 2, nst_np
            jc = sticir_np(1,ist); jr = sticir_np(2,ist)
            if ((jc >= 1).and.(jc <= nc).and.(jr >= 1).and.(jr <= nr)) then
              v = x(jc,jr)
              if (v /= x_src) then
                if (v /= x_mv) then
                  n = n + 1
                  if (n > n_max) then
                    call errmsg('tIntf_get_bnd_i4: program error.')
                  end if
                  int_icir_wk(1,n) = ic; int_icir_wk(2,n) = ir
                  ext_icir_wk(1,n) = jc; ext_icir_wk(2,n) = jr
                  ext_v_wk(n) = v
                else
                  m = m + 1
                  int_mv_icir_wk(1,m) = ic; int_mv_icir_wk(2,m) = ir
                  ext_mv_icir_wk(1,m) = jc; ext_mv_icir_wk(2,m) = jr
                end if
              end if
            end if
          end do
        end if
      end do
    end do
    !
    if (present(int_icir)) then
      if (allocated(int_icir)) deallocate(int_icir)
      if (n > 0) then
        allocate(int_icir(2,n))
        do i = 1, n
          int_icir(1,i) = int_icir_wk(1,i)
          int_icir(2,i) = int_icir_wk(2,i)
        end do
      end if
    end if
    if (present(ext_icir)) then
      if (allocated(ext_icir)) deallocate(ext_icir)
      if (n > 0) then
        allocate(ext_icir(2,n))
        do i = 1, n
          ext_icir(1,i) = ext_icir_wk(1,i)
          ext_icir(2,i) = ext_icir_wk(2,i)
        end do
      end if
    end if
    if (present(ext_v)) then
      if (allocated(ext_v)) deallocate(ext_v)
      if (n > 0) then
        allocate(ext_v(n))
        do i = 1, n
          ext_v(i) = ext_v_wk(i)
        end do
      end if
    end if
    if (present(int_mv_icir)) then
      if (allocated(int_mv_icir)) deallocate(int_mv_icir)
      if (m > 0) then
        allocate(int_mv_icir(2,m))
        do i = 1, m
          int_mv_icir(1,i) = int_mv_icir_wk(1,i)
          int_mv_icir(2,i) = int_mv_icir_wk(2,i)
        end do
      end if
    end if
    if (present(ext_mv_icir)) then
      if (allocated(ext_mv_icir)) deallocate(ext_mv_icir)
      if (m > 0) then
        allocate(ext_mv_icir(2,m))
        do i = 1, m
          ext_mv_icir(1,i) = ext_mv_icir_wk(1,i)
          ext_mv_icir(2,i) = ext_mv_icir_wk(2,i)
        end do
      end if
    end if
    !
    ! cleanup
    deallocate(int_icir_wk, ext_icir_wk, ext_v_wk)
    deallocate(int_mv_icir_wk, ext_mv_icir_wk)
    !
    return
  end subroutine get_bnd_i4_9p
  
  subroutine tQuad_init(this, gid, lid, bbo, lid_prent, gid_prent, ntile, cs_max)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in), optional :: gid
    integer(I4B), intent(in), optional :: lid
    type(tBbObj), intent(in), optional :: bbo
    integer(I4B), intent(in), optional :: gid_prent
    integer(I4B), intent(in), optional :: lid_prent
    integer(I4B), intent(in), optional :: ntile
    real(R8B),    intent(in), optional :: cs_max
    ! -- local
! ------------------------------------------------------------------------------
    !
    this%mapped_gids = .false.
    this%active      = .true.
    this%work_done   = .false.
    this%generated   = .false.
    !
    this%n_child     = I4ZERO
    this%props       => null()
    this%area        = I4ZERO
    this%weight      = R8ZERO
    this%bb_ratio    = R8ZERO
    this%bb_xm       = R8ZERO
    this%bb_ym       = R8ZERO
    this%itile_bb_xm = I4ZERO
    this%ntile       = I4ZERO
    if (present(ntile)) then
      this%ntile = ntile
      allocate(this%itile(ntile))
      this%itile = I4ZERO
    end if
    this%i_graph     = I4ZERO
    this%i_prop      = I4ZERO
    this%gids_mv_map = I4ZERO
    this%cs_max      = R8ZERO
    this%mod_dir     = ''
    this%hdr_gids    => null()
    call this%bbo%init()
    !
    ! optional arguments
    if (present(gid)) then
      this%gid = gid
    end if
    if (present(lid)) then
      this%lid = lid
    end if
    if (present(gid_prent)) then
      this%gid_prent = gid_prent
    end if
    if (present(lid_prent)) then
      this%lid_prent = lid_prent
    end if
    if (present(gid_prent).or.present(lid_prent)) then
      this%generated = .true.
    end if
    if (present(bbo)) then
      this%bbo = bbo
    end if
    if (present(cs_max)) then
      this%cs_max = cs_max
    end if
    !
    this%intf      => null()
    this%c2f_grid  => null()
    this%f2c_grid  => null()
    this%lay_mods  => null()
    this%dat_mod   => null()
    this%disu      => null()
    !
    return
  end subroutine tQuad_init
  
  subroutine tQuad_init_map_gids(this, mv, mv_map)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in) :: mv
    integer(I4B), intent(in) :: mv_map
    ! -- local
    type(tHdr), pointer :: gids => null()
    type(tHdrHdr), pointer :: hdr => null()
    type(tBb) :: bbi
    type(tBbX) :: bbx 
! ------------------------------------------------------------------------------
    if (associated(this%hdr_gids)) deallocate(this%hdr_gids)
    allocate(this%hdr_gids)
    !
    gids => this%hdr_gids
    allocate(gids%hdr); hdr => gids%hdr
    !
    bbi = this%bbo%child_bbi; bbx = this%bbo%child_bbx
    call hdr%init(ncol=bbi%ncol, nrow=bbi%nrow, nbits=32, pixeltype='signedint', &
      xllr8=bbx%xll, yllr8=bbx%yll, xurr8=bbx%xur, yurr8=bbx%yur, csr8=bbx%cs, mvi4=mv)
    call gids%init(map_active=.true., mapped_all=.false., i_data_type=i_i4)
    !
    allocate(gids%dat%xi4(hdr%ncol,hdr%nrow))
    this%gids_mv_map = mv_map
    gids%dat%xi4 = mv_map
    !
    this%mapped_gids = gids%mapped_all
    return
  end subroutine tQuad_init_map_gids
  
  subroutine tQuad_map_gids(this, prent_gids, itile)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    type(tHdr), intent(in) :: prent_gids
    integer(I4B), intent(in) :: itile
    ! -- local
    type(tHdr), pointer :: child_gids => null()
! ------------------------------------------------------------------------------
    !
    if (this%itile(itile) == 0) then
      return
    end if
    this%itile(itile) = -1
    !
    child_gids => this%hdr_gids
    call child_gids%map_grid(prent_gids, mvi4=this%gids_mv_map)
    if (.not.(any(this%itile == 1))) then
      this%mapped_gids = .true.; this%itile = abs(this%itile)
      child_gids%mapped_all = (.not.child_gids%check_val(vi4=this%gids_mv_map))
      if (.not.child_gids%mapped_all) then
        !call child_gids%write('f:\models\lhm\LHM-server\MOZART\HKV\'//ta([this%gid]))
        call logmsg('tQuad_map_gids: WARNING some data was not read for gid = ' &
          //ta([this%gid])//'!')
      end if
    else
      this%mapped_gids = .false.
    end if
    !
    return
  end subroutine tQuad_map_gids

  function tQuad_check_gid_mask(this, gids_in) result(ok)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), dimension(:,:), intent(in), target, optional :: gids_in
    logical :: ok
    !
    ! -- local
    type(tIntf), pointer :: intf => null()
    type(tBbX) :: bbx
    type(tBb) :: bbi
    !
    character(len=MXSLEN) :: fp
    integer(I4B), dimension(:,:), allocatable :: mask
    integer(I4B), dimension(:,:), pointer :: gids => null()
    integer(I4B) :: nc, nr, ic, ir, n
! ------------------------------------------------------------------------------
    !
    if (present(gids_in)) then
      gids => gids_in
    else
      if(.not.associated(this%hdr_gids)) then
        call errmsg('tQuad_check_gid_mask: read gids not found.')
      end if
      gids => this%hdr_gids%dat%xi4
    end if
    if(.not.associated(this%intf)) then
      call errmsg('tQuad_check_gid_mask: interface not found.')
    end if
    !
    intf => this%intf
!    if (this%gid == 5529638) then
!      call intf%get_mask(mask, bbx=this%bbo%child_bbx, gid=this%gid)
!    else
       call intf%get_mask(mask, bbx=this%bbo%child_bbx)
!    end if
    !
    ok = .true.
    if (size(gids,1) /= size(mask,1)) then
      ok = .false.
    end if
    if (size(gids,2) /= size(mask,2)) then
      ok = .false.
    end if
    if (.not.ok) then
      deallocate(mask)
      return
    end if
    !
    nc = size(gids,1); nr = size(gids,2)
    !
    n = 0
    do ir = 1, nr
      do ic = 1, nc
        if (gids(ic,ir) == this%gid) then
          if (mask(ic,ir) /= 1) then
            mask(ic,ir) = -1
            n = n + 1
          end if
        else
          if (mask(ic,ir) /= 0) then
            mask(ic,ir) = -1
            n = n + 1
          end if
        end if
      end do
    end do
    !
    if (n > 0) then
      ok = .false.
      call logmsg('Adding '//ta([n])//' fill value(s) to interface for gid = '//ta([this%gid]))
      intf%n_fill = n
      allocate(intf%fill_lcell_nr(intf%n_fill))
      call this%get_bb(child_bbi=bbi)
      n = 0
      do ir = 1, nr
        do ic = 1, nc
          if (mask(ic,ir) == -1) then
            n = n + 1
            call icrl_to_node(intf%fill_lcell_nr(n), ic, ir, 1, bbi%ncol, bbi%nrow)
          end if
        end do
      end do
      !fp = 'f:\models\pcr-globwb-flex\data\hblev12_tiles\mask_error'
      !bbx = this%bbo%child_bbx
      !call writeflt(fp, mask, nc, nr, bbx%xll, bbx%yll, bbx%cs, 0)
    end if
    !
    deallocate(mask)
    !
    return
  end function tQuad_check_gid_mask
!  
  subroutine tQuad_get_nod_map(this, il, nod_map, nodes, f_csv_dat, nodes_offset)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in) :: il
    integer(I4B), dimension(:,:), allocatable, intent(inout) :: nod_map
    integer(I4B), intent(out) :: nodes
    character(len=*), intent(in), optional :: f_csv_dat
    integer(I4B), optional :: nodes_offset
    ! -- local
    type(tMF6Disu), pointer :: disu => null()
    type(tBB) :: bbi
    integer(I4B) :: jl, ir, ic, n
! ------------------------------------------------------------------------------
    ! read the grid data
    if (present(f_csv_dat)) then
      call this%grid_init(f_csv_dat, '_'//trim(ta([this%lid])))
    else
      call this%grid_init()
    end if
    !
    disu => this%disu
    nodes = disu%nodes
    call this%get_bb(child_bbi=bbi)
    !
    if (allocated(nod_map)) deallocate(nod_map)
    !
    if (il == 0) then
      call disu%x_to_top_nodes_grid(nod_map)
    else
      do jl = 1, disu%nlay_act
        if (il == disu%lay_act(jl)) then
          allocate(nod_map, source=disu%grid_x_nod(:,:,jl))
          exit
        end if
      end do
    end if
    !
    if (present(nodes_offset)) then
      do ir = 1, size(nod_map,2)
        do ic = 1, size(nod_map,1)
          n = nod_map(ic,ir)
          if (n > 0) then
            nod_map(ic,ir) = n + nodes_offset
          end if
        end do
      end do
    end if
    !
    !clean-up grid data
    call disu%clean(); deallocate(this%disu); this%disu => null()
    !
    return
  end subroutine tQuad_get_nod_map
  
  subroutine tQuad_clean_gids(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
  
    ! -- local
    type(tHdr), pointer :: gids => null()
! ------------------------------------------------------------------------------
    gids => this%hdr_gids
    if (associated(gids)) then
      call gids%clean()
      deallocate(this%hdr_gids)
      this%hdr_gids => null()
    end if
    !
    return
  end subroutine tQuad_clean_gids
      
  function tQuad_get_area(this) result(area)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B) :: area
    ! -- local
! ------------------------------------------------------------------------------
    !
    area = this%area
    !
    return
  end function tQuad_get_area

  function tQuad_get_hlev(this, i) result(hlev)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in) :: i
    integer(I4B) :: hlev ! result
    ! -- local
! ------------------------------------------------------------------------------
    !
    hlev = 0
    if (allocated(this%hlev)) then
      if ((i >= 1).and.(i <= size(this%hlev))) then
        hlev = this%hlev(i)
      end if
    end if
    !
    return
  end function tQuad_get_hlev
  
  subroutine tQuad_calc_prop(this, itile, tile_bbx, mask, weight)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in), optional :: itile
    type(tBbx), intent(in), optional :: tile_bbx
    integer(I4B), dimension(:,:), intent(in), optional :: mask
    integer(I4B), dimension(:,:), intent(in), optional :: weight
    ! -- local
    type(kdtree2), pointer :: tree
    type(kdtree2_result), dimension(:), allocatable :: res
    type(tBb) :: bbi, bbit
    type(tBbX) :: bbx
    character(len=MXSLEN) :: fp
    integer(I4B) :: nc, nr, ir, ic, irm, icm,  mc, mr, irm0, irm1, icm0, icm1
    integer(I4B) :: i, n, n_gid, n0, n1
    integer(I4B), dimension(:,:), allocatable :: xi4
    real(R8B) :: xm, ym, w
    real(R8B), dimension(:,:), allocatable :: xy
! ------------------------------------------------------------------------------
    !
    bbi = this%bbo%child_bbi
    bbx = this%bbo%child_bbx
    if (present(mask)) then
      allocate(xi4,source=mask)
    else
      if (.not.associated(this%intf)) then
        call errmsg('tQuad_calc_prop: program error, interface not found.')
      end if
      call this%intf%get_mask(xi4, bbx=bbx)
    end if
    if ((size(xi4,1) /= bbi%ncol).or.(size(xi4,2) /= bbi%nrow)) then
      call errmsg('tQuad_calc_prop: program error, invalid dimensions.')
    end if
    !
    ! calculate area
    n_gid = 0
    do ir = 1, bbi%nrow
      do ic = 1, bbi%ncol
        if (xi4(ic,ir) == 1) n_gid = n_gid + 1
      end do
    end do
    !
    this%weight = R8ZERO
    if (present(weight)) then
      do ir = 1, bbi%nrow
        do ic = 1, bbi%ncol
          if (xi4(ic,ir) == 1) then
            this%weight = this%weight + real(weight(ic,ir),R8B)
          end if
        end do
      end do
    end if
    !
    if (n_gid == 0) then
      fp = 'f:\models\pcr-globwb-30arcsec\hydrobasins\lev12_tight\gid_'//ta([this%gid])
      call writeflt(fp, xi4, bbi%ncol, bbi%nrow, bbx%xll, bbx%yll, bbx%cs, 0)
      !call this%hdr_gids%write(f)
      call errmsg('tQuad_calc_area: no cells found for id = '//ta([this%gid]))
    end if
    this%area = n_gid
    !
    allocate(xy(2,n_gid))
    !
    ! compute tight bounding box and fill coordinates
    n = 0; xm = R8ZERO; ym = R8ZERO
    do ir = 1, bbi%nrow
      do ic = 1, bbi%ncol
        if (xi4(ic,ir) == 1) then
          bbit%ic0 = min(bbit%ic0, ic); bbit%ic1 = max(bbit%ic1, ic)
          bbit%ir0 = min(bbit%ir0, ir); bbit%ir1 = max(bbit%ir1, ir)
          n = n + 1
          call get_xy(xy(1,n), xy(2,n), ic, ir, bbx%xll, bbx%yur, bbx%cs)
          xm = xm + xy(1,n); ym = ym + xy(2,n)
        end if
      end do
    end do
    bbit%ncol = bbit%ic1 - bbit%ic0 + 1; bbit%nrow = bbit%ir1 - bbit%ir0 + 1
    !
    ! compute the ratio
    if (bbit%ncol > bbit%nrow) then
      this%bb_ratio = real(bbit%ncol,R8B)/real(bbit%nrow,R8B)
    else
      this%bb_ratio = real(bbit%nrow,R8B)/real(bbit%ncol,R8B)
    end if
    !
    ! compute the center
    xm = xm / n_gid; ym = ym / n_gid
    !
    if (n == 1) then
      this%bb_xm = xy(1,1); this%bb_ym = xy(2,1)
    else
      allocate(res(1))
      tree => kdtree2_create(xy,sort=.false.,rearrange=.false.)
      call kdtree2_n_nearest(tp=tree,qv=([xm, ym]),nn=1,results=res)
      i = res(1)%idx; this%bb_xm = xy(1,i); this%bb_ym = xy(2,i)
      call kdtree2_destroy(tree)
      if (allocated(res)) deallocate(res)
    end if
    if (allocated(xy)) deallocate(xy)
    !
    ! determine if (xm,ym) is within the tile
    if (present(itile)) then
      if (.not.present(tile_bbx)) then
        call errmsg('tQuad_calc_prop: tile_bbx is not present.')
      end if
      if (point_in_bb(xy=reshape([this%bb_xm, this%bb_ym], [2, 1]), bbx=tile_bbx)) then
        this%itile_bb_xm = itile
      end if
    end if
    !
    !clean up
    if (allocated(xi4)) deallocate(xi4)

    !mc = nint((real(nc,R8B)-1.d0)/2.d0)
    !mr = nint((real(nr,R8B)-1.d0)/2.d0)
    !icm = bb%ic0 + mc; irm = bb%ir0 + mr
    !icm = max(icm, bb%ic0); icm = min(icm, bb%ic1)
    !irm = max(irm, bb%ir0); irm = min(irm, bb%ir1)
    !!
    !if (gids%xi4(icm,irm) /= this%gid) then
    !  n0 = 0; n1 = 0
    !  if (nc > nr) then ! loop over the rows
    !    do ir = irm, bb%ir0, -1
    !      if (gids%xi4(icm,ir) == this%gid) then
    !        n0 = n0 + 1; irm0 = ir; exit
    !      end if
    !    end do
    !    do ir = irm, bb%ir1
    !      if (gids%xi4(icm,ir) == this%gid) then
    !        n1 = n1 + 1; irm1 = ir; exit
    !      end if
    !    end do
    !    if (n0 > n1) then
    !      irm = irm0
    !    else
    !      irm = irm1
    !    end if
    !  else
    !    do ic = icm, bb%ic0, -1
    !      if (gids%xi4(ic,irm) == this%gid) then
    !        n0 = n0 + 1; icm0 = ic; exit
    !      end if
    !    end do
    !    do ic = icm, bb%ic1
    !      if (gids%xi4(ic,irm) == this%gid) then
    !        n1 = n1 + 1; icm1 = ic; exit
    !      end if
    !    end do
    !    if (n0 > n1) then
    !      icm = icm0
    !    else
    !      icm = icm1
    !    end if
    !    !
    !  end if
    !  ! check
    !  if (gids%xi4(icm,irm) /= this%gid) then
    !    call errmsg('tQuad_calc_bb_prop: could not find center.')
    !  end if
    !end if
    !
    !
    return
  end subroutine tQuad_calc_prop
  
  subroutine tQuad_get_prop_csv(this, ic, key, ikey, i1v, i2v, i4v, i8v, &
    r4v, r8v, cv, empty)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    !
    integer(I4B),     intent(in), optional  :: ic
    character(len=*), intent(in), optional  :: key
    integer(I4B),     intent(in), optional  :: ikey
    integer(I1B),     intent(out), optional :: i1v
    integer(I2B),     intent(out), optional :: i2v
    integer(I4B),     intent(out), optional :: i4v
    integer(I8B),     intent(out), optional :: i8v
    real(R4B),        intent(out), optional :: r4v
    real(R8B),        intent(out), optional :: r8v
    character(len=*), intent(out), optional :: cv
    logical, intent(out), optional :: empty
    !
    ! -- local
    type(tProps), pointer :: props => null()
    character(len=MXSLEN) :: key_loc
    logical :: empty_loc
    integer(I4B) :: jc, n
 ! ------------------------------------------------------------------------------
    !
    props => this%props
    if (.not.associated(props)) then
      call errmsg('tQuad_get_prop_csv: properties not found.')
    end if
    !
    if (this%i_prop == 0) then
      call errmsg('tQuad_get_prop_csv: index to properties not set.')
    end if
    !
    n = 0; empty_loc = .true.
    if (present(ic)) then
      jc = ic; n = n + 1
      empty_loc = .false.
    end if
    if (present(key)) then
      jc = props%csv%get_col(key); n = n + 1
      empty_loc = .false.
    end if
    if (present(ikey)) then
      call props%get_field(ikey, key_loc, empty_loc)
      if (.not.empty_loc) then
        jc = props%csv%get_col(key_loc)
        n = n + 1
      end if
    end if
    !
    if (present(empty)) then
      empty = empty_loc
      if (empty_loc) return
    end if
    !
    if (empty_loc) then
      call errmsg('tQuad_get_prop_csv: key not found.')
    end if
    if (n /= 1) then
      call errmsg('tQuad_get_prop_csv: too many input defined.')
    end if
    !
    call props%csv%get_val(ic=jc, ir=this%i_prop, i1v=i1v, i2v=i2v, i4v=i4v, i8v=i8v, &
      r4v=r4v, r8v=r8v, cv=cv)
    !
    return
  end subroutine tQuad_get_prop_csv
  
  subroutine tQuad_set_prop_csv(this, ic, key, ikey, i1v, i2v, i4v, i8v, &
    r4v, r8v, cv)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    !
    integer(I4B),     intent(in), optional  :: ic
    character(len=*), intent(in), optional  :: key
    integer(I4B),     intent(in), optional  :: ikey
    integer(I1B),     intent(in), optional :: i1v
    integer(I2B),     intent(in), optional :: i2v
    integer(I4B),     intent(in), optional :: i4v
    integer(I8B),     intent(in), optional :: i8v
    real(R4B),        intent(in), optional :: r4v
    real(R8B),        intent(in), optional :: r8v
    character(len=*), intent(in), optional :: cv
    !
    ! -- local
    type(tProps), pointer :: props => null()
    character(len=MXSLEN) :: key_loc
    integer(I4B) :: jc, n
 ! ------------------------------------------------------------------------------
    !
    props => this%props
    if (.not.associated(props)) then
      call errmsg('tQuad_set_prop_csv: properties not found.')
    end if
    !
    if (this%i_prop == 0) then
      call errmsg('tQuad_set_prop_csv: index to properties not set.')
    end if
    !
    n = 0
    if (present(ic)) then
      jc = ic; n = n + 1
    end if
    if (present(key)) then
      jc = props%csv%get_col(key); n = n + 1
    end if
    if (present(ikey)) then
      call props%get_field(ikey, key_loc)
      jc = props%csv%get_col(key_loc)
      n = n + 1
    end if
    if (n /= 1) then
      call errmsg('tQuad_get_prop_csv: too many input defined.')
    end if
    !
    call props%csv%set_val(ic=jc, ir=this%i_prop, i1v=i1v, i2v=i2v, i4v=i4v, i8v=i8v, &
      r4v=r4v, r8v=r8v, cv=cv, create_s=.true.)
    !
    return
  end subroutine tQuad_set_prop_csv
    
  subroutine tQuad_get_prop(this, lid, gid, lid_prent, gid_prent, &
    area, weight, bb_ratio, xm, ym, cs_max, itile_bb_xm, i_graph, i_prop)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    !
    integer(I4B), intent(out), optional :: lid
    integer(I4B), intent(out), optional :: gid
    integer(I4B), intent(out), optional :: lid_prent
    integer(I4B), intent(out), optional :: gid_prent
    integer(I4B), intent(out), optional :: area
    real(R8B),    intent(out), optional :: weight
    real(R8B),    intent(out), optional :: bb_ratio
    real(R8B),    intent(out), optional :: xm
    real(R8B),    intent(out), optional :: ym
    real(R8B),    intent(out), optional :: cs_max
    integer(I4B), intent(out), optional :: itile_bb_xm
    integer(I4B), intent(out), optional :: i_graph
    integer(I4B), intent(out), optional :: i_prop
    ! -- local
! ------------------------------------------------------------------------------
    if (present(lid))         lid         = this%lid
    if (present(gid))         gid         = this%gid
    if (present(lid_prent))   lid_prent   = this%lid_prent
    if (present(gid_prent))   gid_prent   = this%gid_prent
    if (present(area))        area        = this%area
    if (present(weight))      weight      = this%weight
    if (present(bb_ratio))    bb_ratio    = this%bb_ratio
    if (present(xm))          xm          = this%bb_xm
    if (present(ym))          ym          = this%bb_ym
    if (present(cs_max))      cs_max      = this%cs_max
    if (present(itile_bb_xm)) itile_bb_xm = this%itile_bb_xm
    if (present(i_graph))     i_graph     = this%i_graph
    if (present(i_prop))      i_prop      = this%i_prop
    !
    return
  end subroutine tQuad_get_prop
  
  subroutine tQuad_calc_bb_ratio(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    ! -- local
    type(tBb) :: bb
    type(tHdr), pointer :: gids => null()
    type(tHdrHdr), pointer :: hdr => null()
    integer(I4B) :: nc, nr, ir, ic
! ------------------------------------------------------------------------------
    !
    gids => this%hdr_gids; hdr => gids%hdr
    !
    ! compute tight bounding box
    do ir = 1, size(gids%dat%xi4,2)
      do ic = 1, size(gids%dat%xi4,1)
        if (gids%dat%xi4(ic,ir) == this%gid) then
          bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
          bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
        end if
      end do
    end do
    nc = bb%ic1 - bb%ic0 + 1; nr = bb%ir1 - bb%ir0 + 1
    if (nc > nr) then
      this%bb_ratio = real(nc,R8B)/real(nr,R8B)
    else
      this%bb_ratio = real(nr,R8B)/real(nc,R8B)
    end if
    !
    return
  end subroutine tQuad_calc_bb_ratio
  
  function tQuad_get_bb_ratio(this) result(bb_ratio)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    real(R8B) :: bb_ratio
    ! -- local
! ------------------------------------------------------------------------------
    !
    bb_ratio = this%bb_ratio
    !
    return
  end function tQuad_get_bb_ratio
  
  subroutine tQuad_calc_interface(this, g2lid, &
    xi4_in, mvi4_in, prent_gid_in)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    !
    integer(I4B), dimension(:), intent(in) :: g2lid
    integer(I4B), dimension(:,:), intent(in), target, optional :: xi4_in
    integer(I4B), intent(in), optional :: mvi4_in
    integer(I4B), intent(in), optional :: prent_gid_in
    ! -- local
    type(tHdr), pointer :: gids => null()
    type(tHdrHdr), pointer :: hdr => null()
    type(tIntf), pointer :: intf => null()
    type(tNbrIntf), pointer :: nbr => null()
    !
    character(len=MXSLEN) :: f, fp
    logical :: flag
    integer(I4B), dimension(1) :: mloc
    integer(I4B), dimension(:), allocatable :: i4wk
    integer(I4B), dimension(:), allocatable :: ext_gid
    integer(I4B), dimension(:,:), pointer :: xi4
    integer(I4B), dimension(:,:), allocatable :: int_icir, ext_icir
    integer(I4B), dimension(:,:), allocatable :: int_mv_icir, ext_mv_icir
    integer(I4B) :: min_gid, max_gid, prent_gid, gid, lid, i, n, n_check, j, inbr
    integer(I4B) :: mvi4, ncol, nrow, ic, ir, ic0, ir0, ic1, ir1, n0, n1
! ------------------------------------------------------------------------------
    if (present(xi4_in)) then
      if ((.not.present(mvi4_in)).or.(.not.present(prent_gid_in))) then
        call errmsg('tQuad_calc_interface: mvi4_in or prent_gid not found.')
      end if
      xi4 => xi4_in
      prent_gid = prent_gid_in; mvi4 = mvi4_in
      ncol = size(xi4,1); nrow = size(xi4,2)
    else
      if(.not.associated(this%hdr_gids)) then
        call errmsg('tQuad_calc_interface: gids nog found.')
      end if
      gids => this%hdr_gids; hdr => gids%hdr; xi4 => gids%dat%xi4
      prent_gid = this%gid; mvi4 = hdr%mvi4; ncol = hdr%ncol; nrow = hdr%nrow
      if (.false.) then
        fp = 'f:\models\lhm\LHM-Flex\pre-processing\debug\'//&
          ta([this%gid])//'_gid'
        call writeflt(fp, xi4, ncol, nrow, hdr%xllr8, hdr%yllr8, hdr%csr8, mvi4)
      end if
    end if
    !
    if (associated(this%intf)) then
      intf => this%intf
      call intf%clean()
      deallocate(this%intf); this%intf => null()
    end if
    !
    allocate(this%intf)
    intf => this%intf
    call intf%init()
    intf%ncol = ncol; intf%nrow = nrow
    !
    ! get interface cells
    call get_bnd_i4_9p(xi4, prent_gid, mvi4, &
      int_icir=int_icir, ext_icir=ext_icir, ext_v=ext_gid, &
      int_mv_icir=int_mv_icir, ext_mv_icir=ext_mv_icir)
    !
    ! process the missing values neigboring cells and store (e.g. for sea-CHD)
    if (allocated(ext_mv_icir)) then
      intf%n_mv = size(ext_mv_icir,2)
      if (allocated(intf%mv_lcell_nr)) deallocate(intf%mv_lcell_nr)
      allocate(intf%mv_lcell_nr(2,intf%n_mv))
      do i = 1, intf%n_mv
        ic0 = int_mv_icir(1,i); ir0 = int_mv_icir(2,i)
        ic1 = ext_mv_icir(1,i); ir1 = ext_mv_icir(2,i)
        !
        call icrl_to_node(n0, ic0, ir0, 1, ncol, nrow)
        call icrl_to_node(n1, ic1, ir1, 1, ncol, nrow)
        !
        flag = has_upstream_val(xi4, prent_gid, ic0, ir0, ic1, ir1)
        if (flag) then
          intf%mv_lcell_nr(1,i) = n0
        else
          intf%mv_lcell_nr(1,i) = -n0
        end if
        intf%mv_lcell_nr(2,i) = n1
      end do
    end if
    !
    ! check exchange interfaces; return if nothing has to be done
    if (.not.allocated(int_icir)) then
      n_check = 0
      !if (this%gid == 485029) then
      !  write(*,*) '@@@ not allocated',xi4
      !  f = 'f:\models\pcr-globwb-30arcsec\hydrobasins\lev12_tight\gid_'//ta([this%gid])
      !  call this%hdr_gids%write(f)
      !  call errmsg('No interface cells found.')
      !end if
    else
      n_check = size(int_icir,2)
      !if (this%gid == 485029) then
      !  write(*,*) '@@@ n_check', n_check
      !end if
    end if
    if (n_check == 0) then
      ! cleanup
      if (allocated(i4wk))        deallocate(i4wk)
      if (allocated(int_icir))    deallocate(int_icir)
      if (allocated(ext_icir))    deallocate(ext_icir)
      if (allocated(int_mv_icir)) deallocate(int_mv_icir)
      if (allocated(ext_mv_icir)) deallocate(ext_mv_icir)
      if (allocated(ext_gid))     deallocate(ext_gid)
      return
    end if
    !if (this%gid == 485029) then
    !  f = 'f:\models\pcr-globwb-30arcsec\hydrobasins\lev12_tight\gid_correct_'//ta([this%gid])
    !  call this%hdr_gids%write(f)
    !end if
    !
    ! determine most occuring neighbor ID
    min_gid = minval(ext_gid); max_gid = maxval(ext_gid)
    if (min_gid == max_gid) then
      intf%mo_nbr = min_gid
      allocate(i4wk(1)); i4wk = size(ext_gid)
    else
      n = max_gid - min_gid + 1
      allocate(i4wk(n)); i4wk = 0
      do i = 1, size(ext_gid)
        lid = ext_gid(i) - min_gid + 1
        i4wk(lid) = i4wk(lid) + 1
      end do
      mloc = maxloc(i4wk); i = mloc(1)
      intf%mo_nbr = min_gid + i - 1
    end if
    !
    ! allocate the interfaces data structure
    intf%n_nbr = 0
    do lid = 1, size(i4wk)
      if (i4wk(lid) > 0) then
        intf%n_nbr = intf%n_nbr + 1
        i4wk(lid) = intf%n_nbr
      end if
    end do
    allocate(intf%nbr_intf(intf%n_nbr))
    !
    ! determine cells for each neighbor
    do inbr = 1, intf%n_nbr
      nbr => intf%nbr_intf(inbr)
      nbr%n_cells = 0
    end do
    do i = 1, size(ext_gid)
      gid = ext_gid(i); lid = gid - min_gid + 1
      inbr = i4wk(lid)
      nbr => intf%nbr_intf(inbr)
      nbr%nbr_gid = gid
      nbr%nbr_lid = g2lid(gid)
      if (nbr%nbr_lid == 0) then
        call errmsg('tQuad_calc_interface: program error.')
      end if
      nbr%n_cells = nbr%n_cells + 1
    end do
    !
    ! allocate
    do inbr = 1, intf%n_nbr
      nbr => intf%nbr_intf(inbr)
      !
      allocate(nbr%lcell_nr(2,nbr%n_cells))
      nbr%lcell_nr = 0
      nbr%n_cells = 0
    end do
    !
    ! fill
    do i = 1, size(ext_gid)
      gid = ext_gid(i); lid = gid - min_gid + 1
      inbr = i4wk(lid)
      nbr => intf%nbr_intf(inbr)
      !
      nbr%n_cells = nbr%n_cells + 1
      !
      ic0 = int_icir(1,i); ir0 = int_icir(2,i)
      ic1 = ext_icir(1,i); ir1 = ext_icir(2,i)
      !
      ! interior local node number
      call icrl_to_node(n0, ic0, ir0, 1, ncol, nrow)
      ! exterior local node number
      call icrl_to_node(n1, ic1, ir1, 1, ncol, nrow)
      !
      flag = has_upstream_val(xi4, prent_gid, ic0, ir0, ic1, ir1)
      if (flag) then
        nbr%lcell_nr(1,nbr%n_cells) = n0
      else
        nbr%lcell_nr(1,nbr%n_cells) = -n0
      end if
      !  
      nbr%lcell_nr(2,nbr%n_cells) = n1
    end do
    !
    ! cleanup
    if (allocated(i4wk))        deallocate(i4wk)
    if (allocated(int_icir))    deallocate(int_icir)
    if (allocated(ext_icir))    deallocate(ext_icir)
    if (allocated(int_mv_icir)) deallocate(int_mv_icir)
    if (allocated(ext_mv_icir)) deallocate(ext_mv_icir)
    if (allocated(ext_gid))     deallocate(ext_gid)
    !
    return
  end subroutine tQuad_calc_interface
  !
  subroutine tQuad_get_intf_ptr(this, intf)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    type(tIntf), intent(out), pointer :: intf
    ! -- local
! ------------------------------------------------------------------------------
    intf => this%intf
    !
    return
  end subroutine tQuad_get_intf_ptr
  
  subroutine tQuad_clean(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (associated(this%hdr_gids)) then
      call this%hdr_gids%clean()
      deallocate(this%hdr_gids)
    end if
    !
    if (associated(this%intf)) then
      call this%intf%clean(); deallocate(this%intf); this%intf => null()
    end if
    ! 
    if (allocated(this%itile)) then
      deallocate(this%itile)
    end if
    if (allocated(this%hlev)) then
      deallocate(this%hlev)
    end if
    !
    if (associated(this%disu)) then
      call this%disu%clean()
      deallocate(this%disu); this%disu => null()
    end if
    !
    call this%init()
    !
    return
  end subroutine tQuad_clean

  subroutine tQuad_set(this, lid_prent, gid_prent, itile, tile_bbi, i_graph, &
    i_prop, hlev, weight, child_bbi)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in), optional :: lid_prent
    integer(I4B), intent(in), optional :: gid_prent
    integer(I4B), intent(in), optional :: itile
    type(tBB), intent(in), optional    :: tile_bbi
    integer(I4B), intent(in), optional :: i_graph
    integer(I4B), intent(in), optional :: i_prop
    integer(I4B), dimension(:), intent(in), optional :: hlev
    real(R8B), intent(in), optional :: weight
    type(tBb), intent(in), optional :: child_bbi
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (present(lid_prent)) then
      this%lid_prent = lid_prent
    end if
    if (present(gid_prent)) then
      this%gid_prent = gid_prent
    end if
    if (present(itile)) then
      if (.not.present(tile_bbi)) then
        call errmsg('tQuad_set: itile was set, but tile_bbi not.')
      end if
      if (bbi_intersect(tile_bbi, this%bbo%child_bbi)) then
        this%itile(itile) = 1
      end if
    end if
    if (present(i_graph)) then
      this%i_graph = i_graph
    end if
    if (present(i_prop)) then
      this%i_prop = i_prop
    end if
    if (present(hlev)) then
      if (allocated(this%hlev)) deallocate(this%hlev)
      allocate(this%hlev,source=hlev)
    end if
    if (present(weight)) then
      this%weight = weight
    end if
    if (present(child_bbi)) then
      this%bbo%child_bbi = child_bbi
    end if
    !
    return
  end subroutine tQuad_set
  !
  subroutine tQuad_get(this, lid_prent, gid_prent)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(out), optional :: lid_prent
    integer(I4B), intent(out), optional :: gid_prent
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (present(lid_prent)) then
      if (.not.this%generated) then
        call errmsg('tQuad_get: quad was not generated.')
      end if 
      lid_prent = this%lid_prent
    end if
    if (present(gid_prent)) then
      if (.not.this%generated) then
        call errmsg('tQuad_get: quad was not generated.')
      end if 
      gid_prent = this%gid_prent
    end if
    !
    return
  end subroutine tQuad_get
  !
  subroutine tQuad_set_cgrid(this, nlev_in, nlev_out)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in) :: nlev_in
    integer(I4B), intent(out) :: nlev_out
    ! -- local
    type(tGrid), pointer :: g0 => null(), g1 => null(), g => null()
    type(tGridArr), pointer :: xlev => null(), xlev_wk => null()
    type(tBbx) :: bbx_mask
    type(tMultiGrid), pointer :: mg => null()
    !
    character(len=MXSLEN) :: f
    logical :: f2c
    integer(I4B), dimension(:,:), allocatable :: icir0, icir1, icir00
    integer(I4B) :: n, nc0, nr0, nc, nr,  nlev, mlev, ilev, jlev
    integer(I4B) :: ilevr, ilev0, ilev1, ilevi, ilevd, ilev1_tmp
    integer(I4B) :: jlev0, jlev1, jlevi, ngrid
    integer(I4B) :: i, j, k, ir, ic, jr, jc, kr, kc, igrid, jgrid, n_check, fi4
    real(R8B) :: cs, cs0, fr8
! ------------------------------------------------------------------------------
    !
    !if (this%gid == 26658) then
    !  write(*,*) '@@@@'
    !  f = 'f:\models\pcr-globwb-flex\data\hblev12_tiles\debug\gid_'//ta([this%gid])
    !  call this%hdr_gids%write(f)
    !end if
    
    if (nlev_in < 0) then
      f2c = .true.
      if (.not.associated(this%f2c_grid)) then
        allocate(this%f2c_grid)
      end if
      mg => this%f2c_grid
    else
      f2c = .false.
      if (.not.associated(this%c2f_grid)) then
        allocate(this%c2f_grid)
      end if
      mg => this%c2f_grid
    end if
    !
    ! initialize: first grid is at the finest level
    nlev = abs(nlev_in)
    allocate(xlev); call xlev%init(nlev)
    !
    ! set the reference level
    if (f2c) then
      ilevr = 1;    ilev0 = 2;      ilev1 = nlev; ilevi =  1; ilevd = -1
    else
      ilevr = nlev; ilev0 = nlev-1; ilev1 = 1;    ilevi = -1; ilevd =  1
    end if
    !
    g => xlev%x(ilevr)
    g%bbx = this%bbo%child_bbx
    !
    call this%get_bb(child_bbx=bbx_mask)
    call this%intf%get_mask(g%xi4, bbx=bbx_mask)
    g%nc = size(g%xi4,1); g%nr = size(g%xi4,2)
    nc0 = g%nc; nr0 = g%nr; cs0 = g%bbx%cs
    !
    ! first, set all values -1
    g%xi4 = -g%xi4
    !
    nc = nc0; nr = nr0; cs = cs0
    do ilev = ilev0, ilev1, ilevi
      g => xlev%x(ilev)
      if (f2c) then
        if ((mod(nc,2) /= 0).or.(mod(nr,2) /= 0)) then
          call logmsg('Maximum coarseness reached.')
          exit
        end if
        ilev1_tmp = ilev
        nc = nc/2; nr = nr/2; cs = cs*2.d0
      else
        nc = nc*2; nr = nr*2; cs = cs/2.d0
      end if
      !
      g%bbx = this%bbo%child_bbx
      g%bbx%cs = cs; g%nc = 0; g%nr = 0
      if ((nc > 0).and.(nr > 0)) then
        allocate(g%xi4(nc,nr)); g%nc = nc; g%nr = nr
        call g%set_const(i4v=I4ZERO)
      end if
      !
    end do
    !
    ! label the cells using the finest grid
    do ilev = ilev0, ilev1, ilevi
      !
      g0 => xlev%x(ilevr); g1 => xlev%x(ilev)
      !
      do ir = 1, g0%nr
        do ic = 1, g0%nc
          if (g0%xi4(ic,ir) /= 0) then
            call get_mapped_icr(g0%bbx, [ic, ir], g1%bbx, icir1)
            n = size(icir1,2)
            do i = 1, n
              jc = icir1(1,i); jr = icir1(2,i)
              if ((jc >= 1).and.(jc <= g1%nc).and.(jr >= 1).and.(jr <= g1%nr)) then
                g1%xi4(jc,jr) = -1
              else
                call errmsg('tQuad_set_cgrid: program error.')
              end if
            end do
          end if
        end do
      end do
    end do
    !
    mlev = 1
    do ilev = ilev0, ilev1, ilevi
      g0 => xlev%x(ilev+ilevd)
      !
      ! get the boundary of ilev-1
      call get_bnd_i4_5p(g0%xi4, -1, 1, int_icir=icir0)
      n_check = size(icir0,2)
      !
      ! check
      if (.not.g0%any_neg()) then
        !mlev = ilev+ilevd
        mlev = ilev + ilevd + 1 !25-08-23
        exit
      end if
      !
      call g0%set_val(icir0, i4v=1)
      !
      ! set to zero from ilev, ..., nlev
      ! label the finer cells
      if (f2c) then
        jlev0 = ilev; jlev1 = nlev; jlevi =  1
      else
        jlev0 = ilev; jlev1 = 1;   jlevi = -1
      end if
      !
      do jlev = jlev0, jlev1, jlevi
        g1 => xlev%x(jlev)
        do i = 1, size(icir0,2) ! cell of ilev-1
          ic = icir0(1,i); ir = icir0(2,i) ! boundary cell of ilev-1
          !call get_mapped_icr(bbx(ilev+ilevd), [ic, ir], bbx(jlev), icir1) ! fine -> coarse
          call get_mapped_icr(g0%bbx, [ic, ir], g1%bbx, icir1) ! fine -> coarse
          do j = 1, size(icir1,2)
            jc = icir1(1,j); jr = icir1(2,j) ! coarse index
            g1%xi4(jc,jr) = 0
            !
            if (jlev == ilev) then
              ! label the cells of ilev-1
              !call get_mapped_icr(bbx(ilev), [jc, jr], bbx(ilev+ilevd), icir00) ! coarse -> fine
              call get_mapped_icr(g1%bbx, [jc, jr], g0%bbx, icir00) ! coarse -> fine
              do k = 1, size(icir00,2)
                kc = icir00(1,k); kr = icir00(2,k) ! fine index
                if (g0%xi4(kc,kr) == -1) then
                  g0%xi4(kc,kr) = 1
                end if
              end do
            end if
          end do
        end do
      end do
      !
      mlev = ilev
    end do
    !
    if (f2c) then
      ilev0 = 1; ilev1 = mlev
      ngrid = mlev
      nlev_out = -ngrid
    else
      !ilev0 = nlev; ilev1 = mlev
      !ngrid = nlev-mlev+1
      ! nlev_out = nlev_in
      !nlev_out = ngrid !14-07-23
      !
      ngrid = nlev - mlev + 1 !25-08-23
      ilev0 = nlev; ilev1 = ilev0 - ngrid + 1 !25-08-23
      nlev_out = ngrid !25-08-23
    end if
    !
    ! contruct the grids, write and count
    allocate(xlev_wk); call xlev_wk%init(ngrid)
    call mg%init(ngrid)
    !allocate(mg%grid_count(mg%ngrid))
    !call mg%grid%init(mg%ngrid)
    !
    ! generate and store the coarse grid
    do igrid = 1, ngrid
      n = 0
      ! copy the levels
      jlev = 0
      do ilev = ilev0, ilev1, ilevi
        jlev = jlev + 1
        call mg%grid(jlev)%set_mv(mvi4=0)
        call mg%grid(jlev)%set_arr(xi4=xlev%x(ilev)%xi4)
        mg%grid(jlev)%bbx = xlev%x(ilev)%bbx
      end do
      !
      do jgrid = 1, igrid
        g => mg%grid(jgrid)
        if (jgrid == igrid) then
          g%xi4 = abs(g%xi4)
        else
          g%xi4 = max(g%xi4,0)
        end if
        !
        do ir = 1, g%nr
          do ic = 1, g%nc
            if (g%xi4(ic,ir) == 1) then
              n = n + 1; g%xi4(ic,ir) = n
            end if
          end do
        end do
        !
        mg%grid_count(igrid) = n
      end do
      !
      ! store the coarse grid
      do jgrid = 1, ngrid
        g => mg%grid(jgrid)
        call g%set_nod_dat()
        call g%clean_xi()
      end do
    end do
    !
    ! clean up
    if (associated(xlev)) then
      call xlev%clean(); deallocate(xlev)
    end if
    if (associated(xlev_wk)) then
      call xlev_wk%clean(); deallocate(xlev_wk)
    end if
    if (allocated(icir0)) deallocate(icir0)
    if (allocated(icir00)) deallocate(icir00)
    if (allocated(icir00)) deallocate(icir00)
    !
    return
  end subroutine tQuad_set_cgrid
  
  subroutine tQuad_clean_cgrid(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (associated(this%f2c_grid)) then
      call this%f2c_grid%clean()
      deallocate(this%f2c_grid); this%f2c_grid => null()
    end if
    if (associated(this%c2f_grid)) then
      call this%c2f_grid%clean()
      deallocate(this%c2f_grid); this%c2f_grid => null()
    end if
    !
    return
  end subroutine tQuad_clean_cgrid
  
  subroutine tQuad_get_grid_dim(this, ig, nlev, nc, nr)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in) :: ig
    integer(I4B), intent(in) :: nlev
    integer(I4B), intent(out) :: nc
    integer(I4B), intent(out) :: nr
    ! -- local
    integer(I4B) :: fac
! ------------------------------------------------------------------------------
    nc = this%bbo%child_bbi%ncol
    nr = this%bbo%child_bbi%nrow
    !
    fac = 2**(ig-1)
    !
    if (nlev < 0) then
      nc = nc/fac; nr = nr/fac
    else
      nc = nc*fac; nr = nr*fac
    end if
    !
    return
  end subroutine tQuad_get_grid_dim
  
  subroutine tQuad_get_grid_bbx(this, ig, nlev, bbx)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in) :: ig
    integer(I4B), intent(in) :: nlev
    type(tBbx), intent(out) :: bbx
    ! -- local
    integer(I4B) :: fac
! ------------------------------------------------------------------------------
    bbx = this%bbo%child_bbx
    !
    fac = 2**(ig-1)
    !
    if (nlev < 0) then
      bbx%cs = bbx%cs * real(fac,R8B)
    else
      bbx%cs = bbx%cs / real(fac,R8B)
    end if
    !
    return
  end subroutine tQuad_get_grid_bbx
  
  subroutine tQuad_get_grid(this, ig_in, renumber, xi4, ncell, bbx_f, n_offset)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(in) :: ig_in
    logical, intent(in) :: renumber
    integer(I4B), dimension(:,:), allocatable, intent(inout) :: xi4
    integer(I4B), intent(out) :: ncell
    type(tBbx), intent(out) :: bbx_f
    integer(I4B), intent(in), optional :: n_offset
    ! -- local
    type(tGrid), pointer :: g => null()
    type(tBbx) ::bbx_ig, bbx_jg
    type(tMultiGrid), pointer :: mg => null()
    !
    logical :: flag
    integer(I4B), dimension(:), allocatable :: i4wk
    integer(I4B), dimension(:,:), allocatable :: icir_ig, icir_jg, icir_f, xi4ig, xi4jg
    integer(I4B) :: nc_f, nr_f, nc_ig, nr_ig, nc_jg, nr_jg, jg, kg, ic, ir, jc, jr, kc, kr
    integer(I4B) :: i, j, k, n, m, m_f, nn, i4v_src, i4v_tgt, jg0, jg1, jgi, ig, nlev
    integer(I4B) :: nos, n_jg
! ------------------------------------------------------------------------------
    !
    if (present(n_offset)) then
      nos = n_offset
    else
      nos = 0
    end if
    !
    if (ig_in == 0) then
      if (associated(this%f2c_grid)) then
        mg => this%f2c_grid
      else
        mg => this%c2f_grid
      end if
      ig = 1
    end if
    if (ig_in < 0) then !f2c
      !mg => this%f2c_grid
      !nlev = -mg%ngrid
      !ig = abs(ig_in) + 1
      !kg = 1 ! identify the finest grid
      mg => this%f2c_grid
      nlev = -mg%ngrid
      !ig = abs(ig_in)
      ig = abs(ig_in) + 1 ! changed 07-06-23
      kg = 1 ! identify the finest grid
    else !c2f
      mg => this%c2f_grid
      nlev = mg%ngrid - 1
      ig = ig_in + 1
      kg = mg%ngrid
    end if 
    !
    call this%get_grid_dim(kg, nlev, nc_f, nr_f)
    call this%get_grid_bbx(kg, nlev, bbx_f)
    if (allocated(xi4)) deallocate(xi4)
    allocate(xi4(nc_f,nr_f)); xi4 = 0
    !
    ! allocate for ig and set
    call this%get_grid_dim(ig, nlev, nc_ig, nr_ig)
    call this%get_grid_bbx(ig, nlev, bbx_ig)
    if (allocated(xi4ig)) deallocate(xi4ig)
    allocate(xi4ig(nc_ig,nr_ig))
    xi4ig = 0
    !
    ! map all the coarser grids
    n = 0
    do jg = 1, mg%ngrid
      !
      ! grid jg:
      g => mg%grid(jg)
      nc_jg = g%nc; nr_jg = g%nr; bbx_jg = g%bbx
      !call this%get_grid_dim(jg, nlev, nc_jg, nr_jg)
      !call this%get_grid_bbx(jg, nlev, bbx_jg)
      call g%get_nod_dat(icir_jg, nc_jg, nr_jg, n_jg)
      !
      if (n_jg == 0) then
        cycle
      end if
      !
      ! loop over stored cells of jg:
      do i = 1, size(icir_jg,2) ! jg cells
        ic = icir_jg(1,i); ir = icir_jg(2,i)
        call get_mapped_icr(bbx_jg, [ic, ir], bbx_ig, icir_ig) ! jg --> ig
        call get_mapped_icr(bbx_jg, [ic, ir], bbx_f,  icir_f)  ! jg --> f
        !
        flag = (bbx_jg%cs <= bbx_ig%cs)
        if (nlev > 0) flag = (.not.flag)
        !
        if (flag) then ! jg is finer than ig: take index of jg
          n = n + 1
          do j = 1, size(icir_f,2)
            jc = icir_f(1,j); jr = icir_f(2,j)
            xi4(jc,jr) = n
            if (jg == 1) xi4(jc,jr) = -abs(xi4(jc,jr))
          end do
        else ! jg is coarser than ig: take index of ig
          ! loop over the fine indices
          do j = 1, size(icir_f,2)
            jc = icir_f(1,j); jr = icir_f(2,j)
            call get_mapped_icr(bbx_f, [jc, jr], bbx_ig, icir_ig) ! jg --> ig
            if (size(icir_ig,2) /= 1) then
              call errmsg('tQuad_get_grid: program error 1.')
            else
              kc = icir_ig(1,1); kr = icir_ig(2,1)
            end if
            if (xi4ig(kc,kr) == 0) then
              n = n + 1
              xi4ig(kc,kr) = n
            end if 
            m = xi4ig(kc,kr)
            if (m == 0) then
              call errmsg('tQuad_get_grid: program error 2.')
            end if
            xi4(jc,jr) = m
            if (jg == 1) xi4(jc,jr) = -abs(xi4(jc,jr))
          end do
        end if
      end do ! jg cells
      !
      ! check
      m = maxval(abs(xi4))
      if (m /= n) then
        call errmsg('tQuad_get_grid: program error 3.')
      end if
    end do
    !
    ! some checks
    ncell = n
    if (ncell /= mg%grid_count(ig)) then
      call errmsg('tQuad_get_grid: program error 4.')
    end if
    !
    ! renumber
    if (renumber) then
      allocate(i4wk(ncell)); i4wk = 0
      do ir = 1, nr_f
        do ic = 1, nc_f
          i4v_src = xi4(ic,ir)
          if (i4v_src < 0) then
            i4wk(abs(i4v_src)) = 1
          end if
        end do
      end do
      !
      ! count
      m_f = 0
      do i = 1, ncell
        if (i4wk(i) == 1) then
          m_f = m_f + 1
          i4wk(i)= m_f
        end if
      end do
      !
      m = m_f
      do ir = 1, nr_f
        do ic = 1, nc_f
          i4v_src = xi4(ic,ir)
          if (i4v_src > 0) then
            if (i4wk(i4v_src) == 0) then
              m = m + 1; i4wk(i4v_src) = m
              i4v_tgt = m
            else
              i4v_tgt = i4wk(i4v_src)
            end if
            xi4(ic,ir) = i4v_tgt
          end if
          xi4(ic,ir) = abs(i4v_src)
        end do
      end do
    end if
    !
    ! apply offset
    if (nos > 0) then
      do ir = 1, nr_f
        do ic = 1, nc_f
          if (xi4(ic,ir) > 0) then
            xi4(ic,ir) = xi4(ic,ir) + nos
          end if
        end do
      end do
    end if
    ! check
    if (minval(i4wk) == 0) then
      call errmsg('tQuad_get_grid: program error')
    end if
    !
    ! clean-up
    if (allocated(icir_ig)) deallocate(icir_ig)
    if (allocated(icir_jg)) deallocate(icir_jg)
    if (allocated(icir_f))  deallocate(icir_f)
    if (allocated(i4wk)) deallocate(i4wk)
    iF (allocated(xi4ig)) deallocate(xi4ig)
    if (allocated(xi4jg)) deallocate(xi4jg)
    !
    return
  end subroutine tQuad_get_grid

  subroutine tQuad_grid_gen(this, nl_act, lay_out, lskip, disu)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B), intent(out) :: nl_act
    integer(I4B), dimension(:), allocatable, intent(inout) :: lay_out
    logical, intent(out) :: lskip
    type(tMF6Disu), intent(inout), pointer :: disu
    ! -- local
    integer(I4B), parameter :: i_w  = 1
    integer(I4B), parameter :: i_e  = 2
    integer(I4B), parameter :: i_n  = 3
    integer(I4B), parameter :: i_s  = 4
    integer(I4B), parameter :: i_t  = 5
    integer(I4B), parameter :: i_b  = 6
    integer(I4B), parameter :: n_sten = i_b
    !
    logical, parameter :: debug = .true.
    character(len=MXSLEN), parameter :: debug_d = 'f:\models\lhm\LHM-Flex\pre-processing\all_partition_test\full_metis\np02048_weight\debug\'
    integer(I4B), parameter :: n_try_max = 10
    real(R4B), parameter :: mv_mask = -12345. 
    !
    type(tBBX) :: bbx_gid, bbx, bbx_read, bbx_wk
    type(tIntf), pointer :: intf => null()
    type(tNbrIntf), pointer :: nintf => null()
    type(tVertIntf), pointer :: vintf => null()
    type(tLayerModel), pointer :: lay_mod => null()
    type(tProps), pointer :: props => null()
    !
    logical :: use_found, found, use_lay_intv, lc2f, lf2c, lok
    !
    character(len=MXSLEN) :: fp, f_binpos
    !
    real(r4B),    dimension(:,:,:), allocatable :: zp
    integer(I4B), dimension(:),     allocatable :: igrid_tgt
    integer(I4B), dimension(:),     allocatable :: lay
    integer(I4B), dimension(:),     allocatable :: rlev
    integer(I4B), dimension(:),     allocatable :: rlev_new
    integer(I4B), dimension(:),     allocatable :: flg
    integer(I4B), dimension(:),     allocatable :: nod
    integer(I4B), dimension(:),     allocatable :: lay_act
    integer(I4B), dimension(:),     allocatable :: nod2cell
    integer(I4B), dimension(:,:),   allocatable :: mask
    integer(I4B), dimension(:,:),   allocatable :: xg
    integer(I4B), dimension(:,:),   allocatable :: xgf
    integer(I4B), dimension(:,:),   allocatable :: bb
    integer(I4B), dimension(:,:,:), allocatable :: x
    !
    real(R8B) :: cs_gid, cs_min_tgt, cs_max_tgt, cs, xc, yc
    real(R4B) :: mv, top, bot, thk, thk_cs_min, z_min, z_max
    !
    integer(I4B) :: ilm, il, jl, il0, il1, ir, ic, nl, nr, nc, n, m, ng, i, j, k
    integer(I4B) :: icell, jcell, ncell, mr, mc, inbr, idummy
    integer(I4B) :: ic0, ic1, ir0, ir1, my_rlev, nbr_rlev, old_rlev, bs
    integer(I4B) :: jc, jr, jc0, jc1, jr0, jr1, jl0, jl1, dl
    integer(I4B) :: ig, ig_min, ig_max, i_sten, n_try
    integer(I4B) :: nlev_c2f, nlev_f2c, max_nlev_c2f, max_nlev_f2c
    integer(I4B) :: nodes, ngrid
! ------------------------------------------------------------------------------
    !
    props => this%props
    call this%get_prop_csv(ikey=i_lay_mod,           i4v=ilm)
    call this%get_prop_csv(ikey=i_tgt_cs_min,        r8v=cs_min_tgt)
    call this%get_prop_csv(ikey=i_tgt_cs_min_dz_top, r4v=thk_cs_min, empty=use_lay_intv)
    if (use_lay_intv) then ! layer interval
      call this%get_prop_csv(ikey=i_tgt_cs_min_lay_beg, i4v=il0)
      call this%get_prop_csv(ikey=i_tgt_cs_min_lay_end, i4v=il1)
    end if
    call this%get_prop_csv(ikey=i_tgt_cs_max,         r8v=cs_max_tgt)
    !
    ! debug
    !il0 = 1
    !il1 = 1
    !cs_min_tgt = 6.25d0
    !cs_max_tgt = 50.d0
    !
    call this%get_bb(child_bbx=bbx_gid); cs_gid = bbx_gid%cs
    !
    ! read the mask
    intf => this%intf
    call intf%get_mask(mask, bbx=bbx_gid)
    if (.false.) then
      fp = 'e:\LHM\selection\simulations\run_output\'//&
        ta([this%gid])//'_mask'
      call writeflt(fp, mask, size(mask,1), size(mask,2), &
        bbx_gid%xll, bbx_gid%yll, bbx_gid%cs, 0)
      stop
    end if
    !
    ! coarse --> fine
    if (cs_gid >= cs_min_tgt) then ! changed 14-07-23
      lc2f = .true.
      max_nlev_c2f = get_number_of_levels(cs_gid, cs_min_tgt)
      call this%set_cgrid(max_nlev_c2f, nlev_c2f)
      cs_min_tgt = cs_gid*2**(real(-nlev_c2f+1,R8B))
    else
      lc2f = .false.
      nlev_c2f = 0
    end if
    !
    ! read the layer models at the FINEST resolution cs_min_tgt or cs_gid
    lay_mod => this%lay_mods%lay_mods(ilm)
    bbx_read = bbx_gid
    if (cs_min_tgt < bbx_gid%cs) then
      bbx_read%cs = cs_min_tgt
    end if
    call lay_mod%read_extent(bbx_read, zp, mv, lay_act, nl_act)
    !
    if (allocated(lay_out)) deallocate(lay_out)
    allocate(lay_out(nl_act))
    do il = 1, size(lay_act)
      i = lay_act(il)
      if (i > 0) lay_out(i) = il
    end do
    !
    ! get the finest number of columns/rows, and the number of active layers
    nc = size(zp,1); nr = size(zp,2); nl = nl_act
    !
    ! apply the mask
    do il = 1, nl + 1
      call apply_mask(mask=mask, bbx_m=bbx_gid, bbx_x=bbx_read, &
        xr4=zp(:,:,il), mvr4=mv_mask)
    end do
    !
    ! fine --> coarse
!    if (cs_gid <= cs_max_tgt) then ! changed 6-6-23
    if (cs_gid < cs_max_tgt) then
      lf2c = .true.
      max_nlev_f2c = get_number_of_levels(cs_gid, cs_max_tgt)
      call this%set_cgrid(max_nlev_f2c, nlev_f2c)
    else
      lf2c = .false.
      nlev_f2c = 0
    end if
    !
    if (lc2f.and.lf2c) then
      ig_min = nlev_f2c + 1 ! coarsest grid
      ig_max = nlev_c2f - 1 ! finest grid
      ngrid = abs(ig_max) + abs(ig_min) + 1
    else
      if (lc2f) then
        ig_max  = nlev_c2f - 1 ! finest grid
        !ig_min = ig_max
        !do ig = 1, ig_max 
        !  cs = cs_gid*2**(real(-ig,R8B))
        !  if (cs_max_tgt == cs) then
        !    ig_min = ig; exit
        !  end if
        !end do
        do ig = ig_max, 0, -1 ! changed 14-7-23
          cs = cs_gid/2**(real(ig,R8B))
          if (cs_max_tgt == cs) then
            ig_min = ig; exit
          end if
        end do
        !
        ngrid = abs(ig_max) + 1
      else
        ig_min = nlev_f2c + 1 ! coarsest grid
        ig_max = ig_min
        do ig = ig_min, -1
          cs = cs_gid*2**(real(abs(ig),R8B))
          if (cs_min_tgt == cs) then
            ig_max = ig; exit
          end if
        end do
        ngrid = abs(ig_min) + 1
      end if
    end if
    !
    ! first, dertermine the target grids
    allocate(igrid_tgt(nl))
    igrid_tgt = 0
    !
    if (use_lay_intv) then ! layer interval
      do il = il0, il1
        igrid_tgt(il) = ig_max
      end do
      do il = il0-1, 1, -1
        igrid_tgt(il) = max(igrid_tgt(il+1)-1, ig_min)
      end do
      do il = il1+1, nl
        igrid_tgt(il) = max(igrid_tgt(il-1)-1, ig_min)
      end do
    else
      igrid_tgt = ig_max
      !
      lok = .true.
      il1 = 1
      do ir = 1, nr; do ic = 1, nc
        thk = 0.
        do il = 1, nl
          top = zp(ic,ir,il); bot = zp(ic,ir,il+1)
          if ((top /= mv_mask).and.(bot /= mv_mask)) then
            if ((top == mv).or.(bot == mv)) then
              lok = .false.
            else
              thk = thk + top - bot
              if (thk <= thk_cs_min) then
                igrid_tgt(il) = ig_max
                il1 = max(il1,il)
              else
                exit
              end if
            end if
          end if
        end do
      end do; end do
      !
      if (.not.lok) then
        call logmsg('WARNING: layer model might not be present.')
      end if
      !
      z_max = -huge(R4ZERO); z_min = huge(R4ZERO)
      do il = 1, nl; do ir = 1, nr; do ic = 1, nc
        top = zp(ic,ir,il); bot = zp(ic,ir,il+1)
        if ((top /= mv_mask).and.(bot /= mv_mask)) then
          if (igrid_tgt(il) == ig_max) then
            z_min = min(z_min, bot)
            z_max = max(z_max, top)
          end if
        end if
      end do; end do; end do
      !call logmsg('Refinement thickness: '//ta([z_max-z_min]))
      !
      do il = il1+1, nl
        igrid_tgt(il) = max(igrid_tgt(il-1)-1, ig_min)
      end do
    end if
    !
    allocate(x(nc,nr,nl))
    x = 0; ncell = 0
    do il = 1, nl
      ig = igrid_tgt(il) 
      call this%get_grid(ig, .true., xg, ng, bbx, n_offset=ncell) ! get grid at the finest level
      ncell = ncell + ng
      call coarse_to_fine_grid(bbx, bbx_read, xi4c=xg, xi4f=xgf)
      if ((size(xgf,1) /= nc).and.(size(xgf,2) /= nr)) then
        call errmsg('tQuad_grid_gen: invalid grid dimensions.')
      end if
      x(:,:,il) = xgf
    end do
    !
    ! determine the bounding boxes
    if ((maxval(x) /= ncell)) then
      call logmsg('maxval(x)='//ta([maxval(x)]))
      call errmsg('tQuad_grid_gen: program error 1.')
    end if
    !
    ! bounding box
    allocate(bb(4,ncell), lay(ncell), rlev(ncell), rlev_new(ncell))
    do i = 1, ncell
      lay(i) = 0; rlev(i) = 0
      bb(1,i) = huge(I4ZERO); bb(2,i) = I4ZERO
      bb(3,i) = huge(I4ZERO); bb(4,i) = I4ZERO
    end do
    !
    do il = 1, nl; do ir = 1, nr; do ic = 1, nc
      i = x(ic,ir,il)
      if (i > 0) then
        lay(i) = il
        bb(1,i) = min(bb(1,i), ic); bb(2,i) = max(bb(2,i), ic)
        bb(3,i) = min(bb(3,i), ir); bb(4,i) = max(bb(4,i), ir)
      end if
    end do; end do; end do
    do i = 1, ncell
      n = bb(2,i) - bb(1,i) + 1
      rlev(i) = int(log(real(n,R8B))/log(2.d0),I4B)
    end do
    !
    ! check for zero layer thickness
    lskip = .false.
    do i = 1, ncell
      il = lay(i); ic0 = bb(1,i); ic1 = bb(2,i); ir0 = bb(3,i); ir1 = bb(4,i)
      do ir = ir0, ir1
        do ic = ic0, ic1
          top = zp(ic,ir,il); bot = zp(ic,ir,il+1)
          if ((top == mv).or.(bot == mv)) then
            !call errmsg('tQuad_grid_gen: program error 2.')
            lskip = .true. ! no layer model is found !
          end if
          thk = top - bot
          if (abs(thk) < R4TINY) then
            lay(i) = -lay(i)
            exit
          end if
        end do 
        if (lay(i) < 0) exit
      end do
      if (lay(i) < 0) then
        do ir = ir0, ir1
          do ic = ic0, ic1
            x(ic,ir,il) = 0
          end do
        end do
      end if
    end do
    !
    if (lskip) then
      !if (associated(mga)) then
        !do il = 1, nl
        !  mg => mga(il)
        !  call mg%clean()
        !end do
      !  call mga%clean()
      !  deallocate(mga); mga => null()
      !end if
      if (allocated(zp       )) deallocate(zp       )
      if (allocated(igrid_tgt)) deallocate(igrid_tgt)
      if (allocated(lay      )) deallocate(lay      )
      if (allocated(rlev     )) deallocate(rlev     )
      if (allocated(rlev_new )) deallocate(rlev_new )
      if (allocated(flg      )) deallocate(flg      )
      if (allocated(nod      )) deallocate(nod      )
      if (allocated(lay_act  )) deallocate(lay_act  )
      if (allocated(nod2cell )) deallocate(nod2cell )
      if (allocated(mask     )) deallocate(mask     )
      if (allocated(xg       )) deallocate(xg       )
      if (allocated(xgf      )) deallocate(xgf      )
      if (allocated(bb       )) deallocate(bb       )
      if (allocated(x        )) deallocate(x        )
      !
      return
    end if
    !
    n_try = 0
    do while(.true.)
      n_try = n_try + 1;
      if (n_try == n_try_max) exit
      n = 0; rlev_new = -1
      do i = 1, ncell
        il = lay(i)
        if (il > 0) then
          ic0 = bb(1,i); ic1 = bb(2,i); ir0 = bb(3,i); ir1 = bb(4,i)
          my_rlev = rlev(i)
          do i_sten = 1, n_sten
            select case(i_sten)
            case(i_w)
              jr0 = ir0; jr1 = ir1; jc0 = ic0-1; jc1 = ic0-1
              jl0 = il; jl1 = il; dl = 1; use_found = .false.
            case(i_e)
              jr0 = ir0; jr1 = ir1; jc0 = ic1+1; jc1 = ic1+1
              jl0 = il; jl1 = il; dl = 1; use_found = .false.
            case(i_n)
              jr0 = ir0-1; jr1 = ir0-1; jc0 = ic0; jc1 = ic1
              jl0 = il; jl1 = il; dl = 1; use_found = .false.
            case(i_s)
              jr0 = ir1+1; jr1 = ir1+1; jc0 = ic0; jc1 = ic1
              jl0 = il; jl1 = il; dl = 1; use_found = .false.
            case(i_t)
              jr0 = ir0; jr1 = ir1; jc0 = ic0; jc1 = ic1
              jl0 = il-1; jl1 = 1; dl = -1; use_found = .true.
            case(i_b)
              jr0 = ir0; jr1 = ir1; jc0 = ic0; jc1 = ic1
              jl0 = il+1; jl1 = nl; dl = 1; use_found = .true.
            end select
            !
            if ((jc0 < 1).or.(jc0 > nc).or.(jc1 < 1).or.(jc1 > nc).or. &
                (jr0 < 1).or.(jr0 > nr).or.(jr1 < 1).or.(jr1 > nr).or. &
                (jl0 < 1).or.(jl0 > nl).or.(jl1 < 1).or.(jl1 > nl)) then
              cycle
            end if
            !
            found = .false.
            do jl = jl0, jl1, dl
              do jr = jr0, jr1; do jc = jc0, jc1
                j = x(jc,jr,jl)
                if (j > 0) then
                  found = .true.
                  nbr_rlev = rlev(j)
                  if (my_rlev /= nbr_rlev) then
                    if (my_rlev > (nbr_rlev + 1)) then ! neighbor is finer
                      if (rlev_new(j) < 0) then
                        rlev_new(j) = my_rlev - 1
                      else
                        rlev_new(j) = max(rlev_new(j), my_rlev - 1)
                      end if
                    end if
                    if (nbr_rlev > (my_rlev + 1)) then ! neighbor is coarse
                      !if (rlev_new(j) < 0) then
                      !  rlev_new(j) = my_rlev + 1
                      !else
                      !  rlev_new(j) = min(rlev_new(j), my_rlev + 1)
                      !end if
                    end if
                  end if
                  !if ((my_rlev - nbr_rlev) > 1) then ! neighbor is finer
                  !  rlev_new(j) = max(rlev_new(j), my_rlev - 1)
                  !end if
                end if
              end do; end do
              if (use_found .and. found) exit
            end do
          end do
        end if
      end do
      !
      n = 0
      do j = 1, ncell
        if ((lay(j) > 0).and.(rlev_new(j) > 0)) then
          n = n + 1
        end if
      end do
      !
      if (n == 0) then
        call logmsg('Grid smoothing done...')
        exit
      else
        call logmsg('Grid smoothing: n = '//ta([n]))
      end if
      !
      do j = 1, ncell
        my_rlev = rlev_new(j)
        if ((lay(j) > 0).and.(my_rlev > 0)) then
          old_rlev = rlev(j)
          if (old_rlev > my_rlev) then
            call errmsg('tQuad_grid_gen: program error 3.')
          end if
          ic0 = bb(1,j); ic1 = bb(2,j); ir0 = bb(3,j); ir1 = bb(4,j)
          bs = 2**my_rlev
          call get_xy(xc, yc, ic0, ir0, bbx_read%xll,  bbx_read%yur, bbx_read%cs)
          call get_icr(jc, jr, xc, yc, bbx_read%xll, bbx_read%yur, bbx_read%cs*bs)
          jc0 = (jc-1)*bs + 1; jc1 = jc0 + bs - 1
          jr0 = (jr-1)*bs + 1; jr1 = jr0 + bs - 1
          !
          ! cancel other cells
          il = lay(j)
          do ir = jr0, jr1; do ic = jc0, jc1
            k = x(ic,ir,il)
            if (k > 0) then
              lay(k) = -abs(lay(k)) ! deactivate
              x(ic,ir,il) = j
            end if
          end do; end do
          !
          ! update my cell
          rlev(j) = rlev_new(j)
          bb(1,j) = jc0; bb(2,j) = jc1; bb(3,j) = jr0; bb(4,j) = jr1
          !
        end if
      end do
      !
      ! check for zero layer thickness
      do i = 1, ncell
        il = lay(i); ic0 = bb(1,i); ic1 = bb(2,i); ir0 = bb(3,i); ir1 = bb(4,i)
        if (il < 0) cycle
        do ir = ir0, ir1
          do ic = ic0, ic1
            top = zp(ic,ir,il); bot = zp(ic,ir,il+1)
            if ((top == mv).or.(bot == mv)) then
              call errmsg('tQuad_grid_gen: program error 2.')
            end if
            thk = top - bot
            if (abs(thk) < R4TINY) then
              lay(i) = -lay(i)
              exit
            end if
          end do 
          if (lay(i) < 0) exit
        end do
        if (lay(i) < 0) then
          do ir = ir0, ir1
            do ic = ic0, ic1
              x(ic,ir,il) = 0
            end do
          end do
        end if
      end do
    end do
    !
    if (n_try == n_try_max) then
      call logmsg('tQuad_grid_gen: maximum number of iterations reached '// &
      'for balancing grid '//ta([this%gid]))
    end if
    !
    ! count the number of cells
    !ncell_fin = 0
    !do icell = 1, ncell
    !  if (lay(icell) > 0) then
    !    ncell_fin = ncell_fin + 1
    !  end if
    !end do
    !call logmsg('# cells for gid '//ta([this%gid])//': '//ta([ncell_fin]))
    !
    ! ============================
    ! generate the node numbers
    ! ============================
    ! label the interior interface node
    allocate(nod(ncell))
    nod = 0
    !
    do inbr = 1, intf%n_nbr
      nintf => intf%nbr_intf(inbr)
      if (nintf%vintf_active == 1) then
        vintf => nintf%vintf_from
        !
        do i = 1, vintf%nlay_act
          jl = vintf%lay_act(i); il = lay_act(jl)
          if (il == 0) then
            call errmsg('tQuad_grid_gen: program error 1.')
          end if
          do jcell = 1, nintf%n_cells
            n = nintf%lcell_nr(1,jcell)
            call node_to_icrl(abs(n), ic, ir, idummy, intf%ncol, intf%nrow)
            !mask(ic,ir) = -1
            m = vintf%nod(jcell,i)
            if (m > 0) then
              call get_xy(xc, yc, ic, ir, bbx_gid%xll, bbx_gid%yur, bbx_gid%cs)
              call get_icr(ic0, ir0, xc-bbx_gid%cs/4.d0, yc+bbx_gid%cs/4.d0, &
                bbx_read%xll, bbx_read%yur, bbx_read%cs)
              call get_icr(ic1, ir1, xc+bbx_gid%cs/4.d0, yc-bbx_gid%cs/4.d0, &
                bbx_read%xll, bbx_read%yur, bbx_read%cs)
              !
              if (minval(x(ic0:ic1,ir0:ir1,il)) /= maxval(x(ic0:ic1,ir0:ir1,il))) then
                call errmsg('tQuad_grid_gen: program error 2.')
              end if
              icell = x(ic0,ir0,il)
              if (lay(icell) /= il) then
                call errmsg('tQuad_grid_gen: program error 3.')
              end if
              nod(icell) = m
            end if
          end do
        end do
      end if
    end do
    !
    nodes = maxval(abs(nod)); disu%nodes_intf = nodes
    do il = 1, nl; do ir = 1, nr; do ic = 1, nc
      icell = x(ic,ir,il)
      if (icell > 0) then
        if ((lay(icell) > 0).and.(nod(icell) == 0)) then
          nodes = nodes + 1
          nod(icell) = nodes
        end if
      end if
    end do; end do; end do
    !
    allocate(nod2cell(nodes)); nod2cell = 0
    do icell = 1, ncell
      n = nod(icell)
      if (n > 0) then
        nod2cell(n) = icell
      end if
    end do
    !
    if (minval(nod2cell) == 0) then
      call errmsg('tQuad_grid_gen: program error.')
    end if
    !
    ! set the node numbers
    allocate(flg(ncell)); flg = 0
    do il = 1, nl; do ir = 1, nr; do ic = 1, nc
      icell = x(ic,ir,il)
      if (icell > 0) then
        if (lay(icell) > 0) then
          x(ic,ir,il) = nod(icell)
          flg(nod(icell)) = flg(nod(icell)) + 1
        else
          x(ic,ir,il)  = 0
        end if
      end if
    end do; end do; end do
    !
    if (.false.) then
      do ig = ig_min, ig_max
        call this%get_grid(ig, .true., xg, ng, bbx) ! get grid at the finest level
        fp = trim(debug_d)//ta([this%gid])//'_f2c_grid_'//ta([ig])
        call writeflt(fp, xg, size(xg,1), size(xg,2), bbx%xll, bbx%yll, bbx%cs, 0)
      end do
      do il = 1, nl
        fp = trim(debug_d)//ta([this%gid])//'_l'//ta([il],'(i2.2)')
        call writeflt(fp, x(:,:,il), nc, nr, &
          bbx_read%xll, bbx_read%yll, bbx_read%cs, 0)
      end do
      fp = trim(debug_d)//ta([this%gid])//'_mask'
      call writeflt(fp, mask, size(mask,1), size(mask,2), &
        bbx%xll, bbx%yll, bbx%cs, 0)
      fp = trim(debug_d)//ta([this%gid])//'_top'
      call writeflt(fp, zp(:,:,1), nc, nr, &
        bbx_read%xll, bbx_read%yll, bbx_read%cs, mv)
      stop
    end if
    !
    ! clean up part 1
    if (allocated(igrid_tgt)) deallocate(igrid_tgt)
    if (allocated(lay      )) deallocate(lay      )
    if (allocated(rlev     )) deallocate(rlev     )
    if (allocated(rlev_new )) deallocate(rlev_new )
    if (allocated(flg      )) deallocate(flg      )
    if (allocated(nod      )) deallocate(nod      )
    if (allocated(lay_act  )) deallocate(lay_act  )
    if (allocated(nod2cell )) deallocate(nod2cell )
    if (allocated(mask     )) deallocate(mask     )
    if (allocated(xg       )) deallocate(xg       )
    if (allocated(xgf      )) deallocate(xgf      )
    if (allocated(bb       )) deallocate(bb       )
    !
    ! create the grid and store the arrays
    !allocate(disu)
    !if(.not.lwrite_asc) then
    !  f_binpos = trim(strip_ext(f_csv_dat))//'.binpos'
    !  call disu%init(f_csv_dat, f_binpos)
    !else
    !  call disu%init(f_csv_dat)
    !end if
    call disu%set_cs_rea_ngrid(x, bbx_read%cs)
    call disu%x_to_mga(x, ngrid, bbx_read)
    !ngrid = disu%ngrid
    !if (.false.) then
    !  fp = trim(debug_d)//ta([this%gid])//'_l'
    !  call disu%grid_mga_nod%write_vrt(fp)
    !end if
!   ! call disu%set(zp, cs_min_tgt)
    call disu%set(zp, bbx_read%cs) ! changed 07-06-23
    !cs_min_rea = disu%cs_min_rea; cs_max_rea = disu%cs_max_rea
    !
    !if (lwrite_disu) then
    !  if (lwrite_asc) then
    !    call disu%write(write_asc=.true., d_out=trim(this%mod_dir)//'\')
    !  else
    !    call disu%write()
    !  end if
    !end if
    !ncell_fin = disu%nodes
    !nja_fin = disu%nja
    !call disu%clean(); deallocate(disu)
    !
    ! clean up part 2
    if (allocated(zp       )) deallocate(zp       )
    if (allocated(x        )) deallocate(x        )
    if (allocated(x)) deallocate(x)
    !
    ! clean the coarse grid
    call this%clean_cgrid()
    !
    return
  end subroutine tQuad_grid_gen
  
  subroutine tQuad_grid_init(this, f_csv_dat, id_pref, nr_max)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    character(len=*), intent(in), optional :: f_csv_dat
    character(len=*), intent(in), optional :: id_pref
    integer(I4B), intent(in), optional :: nr_max
    ! -- local
    character(len=MXSLEN), parameter :: debug_d = &
      'f:\models\lhm\LHM-Flex\pre-processing\debug_read\'
    !
    type(tBb) :: bbi
    type(tBbx) :: bbx
    type(tMF6Disu), pointer :: disu
    type(tLayerModel),  pointer :: lay_mod  => null()
    type(tMf6Wbd), pointer :: wbd => null()
    !
    character(len=MXSLEN) :: f_csv_dat_loc, id_pref_loc, f_binpos, s, fp
    character(len=MXSLEN), dimension(:), allocatable :: sa
    integer(I4B) :: i, lm, ic
    real(R8B) :: cs_min, cs_min_tgt
! ------------------------------------------------------------------------------
    !
    if (present(id_pref)) then
      id_pref_loc = id_pref
    else
      id_pref_loc = ''
    end if
    !
    call this%get_bb(child_bbi=bbi, child_bbx=bbx)
    call this%get_prop_csv(key='cs_min_rea', r8v=cs_min_tgt) !25-08
    cs_min = min(cs_min_tgt, bbx%cs) !25-08
    !
    allocate(this%disu); disu => this%disu
    !
    ! Get the data csv-file, and expand
    if (present(f_csv_dat)) then
      f_csv_dat_loc = f_csv_dat
    else
      call this%get_prop_csv(key='csv_dat', cv=f_csv_dat_loc)
    end if
    !
    f_binpos = trim(strip_ext(f_csv_dat_loc))//'.binpos'
    allocate(disu%wbd); wbd => disu%wbd
    if (present(nr_max)) then
      call wbd%read_csv(f_csv_dat_loc, nr_max=nr_max)
    else
      call wbd%read_csv(f_csv_dat_loc, nr_max=MAX_NR_CSV)
    end if
    wbd%f_binpos = f_binpos
    !
    ! read the grid data
    call disu%read(id_pref_loc)
    !
    ! Set the metadata
    call this%get_prop_csv(key='nodes', i4v=disu%nodes)
    call this%get_prop_csv(key='nlay_act', i4v=disu%nlay_act)
    call this%get_prop_csv(key='lay_act', cv=s)
    disu%cs_min = cs_min !25-08
    call this%get_prop_csv(ikey=i_lay_mod, i4v=lm)
    lay_mod => this%lay_mods%lay_mods(lm)
    disu%nlay = lay_mod%nlay
    call this%get_prop_csv(key='ngrid_lev', i4v=disu%ngrid)
    !
    call parse_line(s, sa, token_in=';')
    if (size(sa) /= disu%nlay_act) then
      call errmsg('tQuad_grid_read: inconsistent data.')
    end if
    allocate(disu%lay_act(disu%nlay_act))
    do i = 1, disu%nlay_act
      read(sa(i),*) disu%lay_act(i)
    end do
    !
    ! determine the mga
    !call this%get_prop_csv(ikey=i_tgt_cs_min, r8v=cs_min)
    call disu%create_mga(bbi, bbx, cs_min)
    if (.false.) then
      fp = trim(debug_d)//ta([this%gid])//'_l'
      call disu%grid_mga_nod%write_vrt(fp)
    end if
    !
    ! convert the mga to x
    call disu%mga_to_x()
    !
    return
  end subroutine tQuad_grid_init
  
  subroutine tQuad_mf6_get_data(this, dmdat, i4a, r8a, r8x)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    type(tDataModelData), pointer, intent(in) :: dmdat
    integer(I4B), dimension(:), allocatable, intent(inout) :: i4a
    real(R8B), dimension(:), allocatable, intent(inout) :: r8a
    real(R8B), dimension(:,:), allocatable, intent(inout) :: r8x
    !
    ! -- local
    type(tBbx) :: bbx
    type(tMF6Disu), pointer :: disu
    real(R8B) :: cs_min_rea
! ------------------------------------------------------------------------------
    disu => this%disu
    call this%get_bb(child_bbx=bbx)
    call this%get_prop_csv(key='cs_min_rea', r8v=cs_min_rea)
    !
    if (allocated(i4a)) deallocate(i4a)
    if (allocated(r8a)) deallocate(r8a)
    if (allocated(r8x)) deallocate(r8x)
    !
    select case(dmdat%i_in_file_type)
    case(i_vrt)
      select case(dmdat%i_type)
      case(i_type_list)      
        call dmdat%mf6_get_data_vrt_list(disu, bbx, cs_min_rea, i4a, r8x)
      case(i_type_array)
        call dmdat%mf6_get_data_vrt_array(disu, bbx, cs_min_rea, i4a, r8a)
      end select
    case(i_vrt_array)
      call dmdat%mf6_get_data_vrt_array_array(disu, bbx, r8a)
    case(i_csv)
      call dmdat%mf6_get_data_csv_list(disu, bbx, i4a, r8a)
    end select
    !
    return
  end subroutine tQuad_mf6_get_data

  subroutine mf6_data_write(wbd, id, i_out_file_type, i4a, r8a, r8x, ndat)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    type(tMf6Wbd), pointer, intent(in) :: wbd
    character(len=*), intent(in) :: id
    integer(I4B), intent(in) :: i_out_file_type
    integer(I4B), dimension(:),   allocatable, intent(in), optional :: i4a
    real(R8B),    dimension(:),   allocatable, intent(in), optional :: r8a
    real(R8B),    dimension(:,:), allocatable, intent(in), optional :: r8x
    integer(I4B), intent(in), optional :: ndat
    ! -- local
    logical :: li4a, lr8a, lr8x, lwrite_list_1d, lwrite_list_nd, lwrite_array
    character(len=MXSLEN) :: d
    integer(I4B) :: ir, ic, nc
    integer(I4B), dimension(:),   allocatable :: i4a_wrk
    real(R8B),    dimension(:),   allocatable :: r8a_wrk
    real(R8B),    dimension(:,:), allocatable :: r8x_wrk
! ------------------------------------------------------------------------------
    !
    ! checks
    li4a = .false.; lr8a = .false.; lr8x = .false.
    if (present(i4a)) then
      if (allocated(i4a)) then
        li4a = .true.
        if (present(ndat)) then
          allocate(i4a_wrk(ndat))
          do ir = 1, ndat
            i4a_wrk(ir) = i4a(ir)
          end do
        else
          allocate(i4a_wrk, source=i4a)
        end if
      end if
    end if
    if (present(r8a)) then
      if (allocated(r8a)) then
        lr8a = .true.
        if (present(ndat)) then
          allocate(r8a_wrk(ndat))
          do ir = 1, ndat
            r8a_wrk(ir) = r8a(ir)
          end do
        else
          allocate(r8a_wrk, source=r8a)
        end if
      end if
    end if
    if (present(r8x)) then
      if (allocated(r8x)) then
        lr8x = .true.
        if (present(ndat)) then
          nc = size(r8x,1)
          allocate(r8x_wrk(nc,ndat))
          do ir = 1, ndat; do ic = 1, nc
            r8x_wrk(ic,ir) = r8x(ic,ir)
          end do; end do
        else
          allocate(r8x_wrk, source=r8x)
        end if
      end if
    end if
    !
    if ((.not.li4a).and.(.not.lr8a).and.(.not.lr8x)) then
      call logmsg('No writing of MODFLOW 6 data...')
      return
    end if
    !
    lwrite_list_1d = .false.
    lwrite_list_nd = .false.
    lwrite_array   = .false.
    if (li4a) then
      if ((.not.lr8a).and.(.not.lr8x)) then
        call errmsg('mf6_data_write: program error.')
      end if
      if (lr8a.and.lr8x) then
        call errmsg('mf6_data_write: program error.')
      end if
      if (lr8a) then
        lwrite_list_1d = .true.
      else
        lwrite_list_nd = .true.
      end if
    else
      lwrite_array = .true.
      if (.not.lr8a) then
        call errmsg('tQuad_mf6_data_write: program error.')
      end if
    end if
    !
    select case(i_out_file_type)
    case(i_binpos)
      if (lwrite_list_1d) then
        call wbd%write_list(id=id, i4a1=i4a_wrk, r8a1=r8a_wrk)
      elseif(lwrite_list_nd) then
        call wbd%write_list(id=id, i4a1=i4a_wrk, r8x=r8x_wrk)
      else
        call wbd%write_array(id=id, r8a=r8a_wrk)
      end if
    case(i_bin)
      d = get_dir(wbd%csv%file)
      if (lwrite_list_1d) then
        call wbd%write_list(id=id, i4a1=i4a_wrk, r8a1=r8a_wrk, &
          f_bin=trim(d)//trim(id)//'.bin')
      elseif(lwrite_list_nd) then
        call wbd%write_list(id=id, i4a1=i4a_wrk, r8x=r8x_wrk, &
          f_bin=trim(d)//trim(id)//'.bin')
      else
        call wbd%write_array(id=id, r8a=r8a_wrk, &
          f_bin=trim(d)//trim(id)//'.bin')
      end if
    case(i_asc)
      d = get_dir(wbd%csv%file)
      if (lwrite_list_1d) then
        call wbd%write_list(id=id, i4a1=i4a_wrk, r8a1=r8a_wrk, &
          f_asc=trim(d)//trim(id)//'.asc')
      elseif(lwrite_list_nd) then
        call wbd%write_list(id=id, i4a1=i4a_wrk, r8x=r8x_wrk, &
          f_asc=trim(d)//trim(id)//'.asc')
      else
        call wbd%write_array(id=id, r8a=r8a_wrk, &
          f_asc=trim(d)//trim(id)//'.asc')
      end if
    end select
    !
    if (allocated(i4a_wrk)) deallocate(i4a_wrk)
    if (allocated(r8a_wrk)) deallocate(r8a_wrk)
    if (allocated(r8x_wrk)) deallocate(r8x_wrk)
    !
    return
  end subroutine mf6_data_write
  
  subroutine tDataModelData_mf6_get_data_vrt_list(this, disu, bbx, cs_min_rea, &
    i4a, r8x)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModelData) :: this
    type(tMF6Disu), intent(in), pointer :: disu
    type(tBbx), intent(in) :: bbx
    real(R8B), intent(in) :: cs_min_rea
    integer(I4B), dimension(:), allocatable, intent(inout) :: i4a
    real(R8B), dimension(:,:), allocatable, intent(inout) :: r8x
    !
    ! -- local
    type(tMultiGridArray), allocatable :: mga_read
    type(tMultiGrid), pointer :: mg      => null()
    type(tMultiGrid), pointer :: mg_top  => null()
    type(tMultiGrid), pointer :: mg_bot  => null()
    type(tMultiGrid), pointer :: mg_dist => null()
    type(tGrid),      pointer :: g       => null()
    type(tData), pointer :: dat => null()
    type(tMf6Wbd), pointer :: wbd => null()
    !
    logical :: ltopisbot
    !
    real(R8B), dimension(:), allocatable :: cs_read
    real(R8B) :: cs_min_tgt
    !
    real(R4B), dimension(:), allocatable :: f1d
    real(R4B) :: mvr4, top, bot, dist, gtop, gbot, dz, f, r4v
    !
    integer(I4B) :: i, ierr, il, jl, ig, ir, ic, n, ncell, nr8a, m
    !
    integer(I1B), dimension(:), allocatable :: top_nodes 
! ------------------------------------------------------------------------------
    !
    if (allocated(i4a)) deallocate(i4a)
    if (allocated(r8x)) deallocate(r8x)
    !
    call disu%x_to_top_nodes(top_nodes)
    ! 
    allocate(mga_read)
    call mga_read%init(this%ndat)
    !
    cs_min_tgt = min(bbx%cs, cs_min_rea)
    allocate(cs_read(disu%ngrid))
    cs_read(disu%ngrid) = cs_min_tgt
    do ig = disu%ngrid-1, 1, -1
      cs_read(ig) = cs_read(ig+1)*2.d0
    end do
    !
    ! step 1: read the data for *ALL* the ngrid levels
    do i = 1, this%ndat
      dat => this%dat(i); mg => mga_read%mga(i)
      !
      if (dat%i_type == i_vrt) then
        call vrt_read_extent_mg(vrt=dat%vrt, bbx=bbx, csa=cs_read, &
          i_uscl=dat%i_uscl, i_dscl=dat%i_dscl, mg=mg, mvr4=mvr4)
      end if
    end do
    !
    allocate(f1d(disu%nodes)); f1d = R4ZERO
    !
    select case(this%i_assign)
    case(i_assign_intersect)
      ierr = 0
      if ((this%itop  <= 0).or.(this%itop  > this%ndat)) ierr = 1
      if ((this%ibot  <= 0).or.(this%ibot  > this%ndat)) ierr = 1
      if ((this%idist <= 0).or.(this%idist > this%ndat)) ierr = 1
      if (ierr == 1) call errmsg('Invalid itop, ibot or idist.')
      !
      !if ((dmdat%ndat == 2).and.(dmdat%itop == dmdat%ibot)) then
      !  ltopisbot = .true.
      !else
      !  ltopisbot = .false.
      !end if
      !
      mg_top  => mga_read%mga(this%itop)
      mg_bot  => mga_read%mga(this%ibot)
      mg_dist => mga_read%mga(this%idist)
      !
      do jl = 1, disu%nlay_act
         mg => disu%grid_mga_nod%get_mg(jl)
         do ig = 1, mg%ngrid
           g => mg%grid(ig)
           do ir = 1, g%nr; do ic = 1, g%nc
             n = g%xi4(ic,ir)
             if (n /= g%mvi4) then
               top  = mg_top%grid(ig)%xr4(ic,ir)
               bot  = mg_bot%grid(ig)%xr4(ic,ir)
               dist = mg_dist%grid(ig)%xr4(ic,ir)
               !
               if ((top  /= mg_top%grid(ig)%mvr4).and. &
                   (bot  /= mg_bot%grid(ig)%mvr4).and. &
                   (dist /= mg_dist%grid(ig)%mvr4)) then
                 !
                 ! make consistent
                 if (bot > top) then
                   call logmsg('Warning: inconsistent top/bot found, setting bot equal to top.')
                   bot = top
                   mg_bot%grid(ig)%xr4(ic,ir) = bot
                 end if
                 !
                 if (top == bot) then
                   ltopisbot = .true.
                 else
                   ltopisbot = .false.
                 end if
                 !
                 gtop = disu%grid_mga_top%mga(jl)%grid(ig)%xr4(ic,ir)
                 gbot = disu%grid_mga_bot%mga(jl)%grid(ig)%xr4(ic,ir)
                 !
                 if (ltopisbot) then
                   if (top_nodes(n) == 1) then
                     if (top >= gbot) then
                       f1d(n) = R4ONE
                     end if
                   else
                     if ((top >= gbot).and.(top < gtop)) then
                       f1d(n) = R4ONE
                     end if
                   end if
                 else
                   dz = top - bot
                   if (dz >= R4ZERO) then
                     f1d(n) = calc_frac(top, bot, top_nodes(n), gtop, gbot)
                     !
                     ! overrule in case of top node
                     if (top_nodes(n) == 1) then
                       if (bot >= gtop) then
                         f1d(n) = R4ONE
                       end if
                     end if
                   end if
                 end if
               end if
             end if
           end do; end do
         end do
      end do
    case(i_assign_exact_layer)
      !
      do jl = 1, disu%nlay_act
        il = disu%lay_act(jl)
         mg => disu%grid_mga_nod%get_mg(jl)
         do ig = 1, mg%ngrid
           g => mg%grid(ig)
           do ir = 1, g%nr; do ic = 1, g%nc
             n = g%xi4(ic,ir)
             if (n /= g%mvi4) then
               if (il == this%ilay) then
                 f1d(n) = R4ONE
               end if
             end if
           end do; end do
         end do
      end do
    case(i_assign_first_layer)
      !
      do jl = 1, disu%nlay_act
         mg => disu%grid_mga_nod%get_mg(jl)
         do ig = 1, mg%ngrid
           g => mg%grid(ig)
           do ir = 1, g%nr; do ic = 1, g%nc
             n = g%xi4(ic,ir)
             if (n /= g%mvi4) then
               if (top_nodes(n) == 1) then
                 f1d(n) = R4ONE
               end if
             end if
           end do; end do
         end do
      end do
    end select
    !
    ! set fraction to zero when dist is zero
    if (.true.) then
      ncell = 0
      do jl = 1, disu%nlay_act
         mg => disu%grid_mga_nod%get_mg(jl)
         do ig = 1, mg%ngrid
           g => mg%grid(ig)
           do ir = 1, g%nr; do ic = 1, g%nc
             n = g%xi4(ic,ir)
             if (n /= g%mvi4) then
               f = f1d(n)
               if (f > R4ZERO) then
                 i = this%idist
                 r4v = mga_read%mga(i)%grid(ig)%xr4(ic,ir)
                 if (r4v == R4ZERO) then
                   f1d(n) = R4ZERO
                   ncell = ncell + 1
                 end if
               end if
             end if
          end do; end do
        end do
      end do
      if (ncell > 0) then
        call logmsg('Skipping '//ta([ncell])//' cells with zero value...')
      end if
    end if
    !
    ! count using the fractions and allocate
    ncell = 0
    do i = 1, disu%nodes
      if (f1d(i) > R4ZERO) ncell = ncell + 1
    end do
    if (ncell == 0) then
      call logmsg('No data found...')
      ! clean up
      if (allocated(mga_read)) then
        call mga_read%clean()
        deallocate(mga_read)
      end if
      if (allocated(i4a)) deallocate(i4a)
      if (allocated(i4a)) deallocate(i4a)
      
      if (allocated(r8x)) deallocate(r8x)
      if (allocated(f1d)) deallocate(f1d)
      if (allocated(cs_read)) deallocate(cs_read)
      if (allocated(top_nodes)) deallocate(top_nodes)
      return
    end if
    nr8a = this%ndat
    allocate(i4a(ncell), r8x(nr8a,ncell))
    !
    ! buffer the data
    m = 0
    do jl = 1, disu%nlay_act
       mg => disu%grid_mga_nod%get_mg(jl)
       do ig = 1, mg%ngrid
         g => mg%grid(ig)
         do ir = 1, g%nr; do ic = 1, g%nc
           n = g%xi4(ic,ir)
           if (n /= g%mvi4) then
             f = f1d(n)
             if (f > R4ZERO) then
               m = m + 1
               i4a(m) = n
               do i = 1, nr8a
                 r4v = mga_read%mga(i)%grid(ig)%xr4(ic,ir)
                 if (i == this%idist) then
                   r4v = r4v*f
                 end if
                 r8x(i,m) = real(r4v, R8B)
               end do
             end if
           end if
         end do; end do
       end do
    end do
    !
    ! check the consistency top/bot
    if (this%itop /= this%ibot) then
      do m = 1, ncell
        if (r8x(this%ibot,m) > r8x(this%itop,m)) then
          call errmsg('Inconsistent data: bot > top.')
        end if
      end do
    end if
    !
    ! clean up
    if (allocated(mga_read)) then
      call mga_read%clean()
      deallocate(mga_read)
    end if
    if (allocated(f1d)) deallocate(f1d)
    if (allocated(cs_read)) deallocate(cs_read)
    if (allocated(top_nodes)) deallocate(top_nodes)
    !
    return 
  end subroutine tDataModelData_mf6_get_data_vrt_list
  
  subroutine tDataModelData_mf6_get_data_vrt_array(this, disu, bbx, &
    cs_min_rea, i4a, r8a)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModelData) :: this
    type(tMF6Disu), intent(in), pointer :: disu
    type(tBbx), intent(in) :: bbx
    real(R8B), intent(in) :: cs_min_rea
    integer(I4B), dimension(:), allocatable, intent(inout) :: i4a
    real(R8B), dimension(:), allocatable, intent(inout) :: r8a
    !
    ! -- local
    type(tMultiGridArray), allocatable :: mga_read
    type(tMultiGrid), pointer :: mg      => null()
    type(tGrid),      pointer :: g       => null()
    type(tData), pointer :: dat => null()
    type(tMf6Wbd), pointer :: wbd => null()
    !
    logical :: ltopisbot
    !
    real(R8B), dimension(:), allocatable :: cs_read
    real(R8B) :: cs_min_tgt
    !
    real(R4B) :: mvr4, r4v
    !
    integer(I4B) :: i, jl, ig, ic, ir, n
    !
    integer(I1B), dimension(:), allocatable :: top_nodes 
! ------------------------------------------------------------------------------
    !
    if (allocated(i4a)) deallocate(i4a)
    if (allocated(r8a)) deallocate(r8a)
    !
    call disu%x_to_top_nodes(top_nodes)
    ! 
    allocate(mga_read)
    call mga_read%init(this%ndat)
    !
    cs_min_tgt = min(bbx%cs, cs_min_rea)
    allocate(cs_read(disu%ngrid))
    cs_read(disu%ngrid) = cs_min_tgt
    do ig = disu%ngrid-1, 1, -1
      cs_read(ig) = cs_read(ig+1)*2.d0
    end do
    
    if (this%i_assign == i_assign_first_layer) then
      !
      ! step 1: read the data for *ALL* the ngrid levels
      do i = 1, this%ndat
        dat => this%dat(i); mg => mga_read%mga(i)
        !
        if (dat%i_type == i_vrt) then
          call vrt_read_extent_mg(vrt=dat%vrt, bbx=bbx, csa=cs_read, &
            i_uscl=dat%i_uscl, i_dscl=dat%i_dscl, mg=mg, mvr4=mvr4)
        end if
      end do
      !
      allocate(r8a(disu%nodes)); r8a = this%r8const
      !
      do jl = 1, disu%nlay_act
         mg => disu%grid_mga_nod%get_mg(jl)
         do ig = 1, mg%ngrid
           g => mg%grid(ig)
           do ir = 1, g%nr; do ic = 1, g%nc
             n = g%xi4(ic,ir)
             if (n /= g%mvi4) then
               if (top_nodes(n) == 1) then
                 r4v = mga_read%mga(1)%grid(ig)%xr4(ic,ir)
                 r8a(n) = real(r4v,R8B)
               end if
             end if
           end do; end do
         end do
      end do
      !
    else
      call errmsg('Array output not supported for vrt input.')
    end if
    !
    ! clean up
    call mga_read%clean(); deallocate(mga_read)
    if (allocated(cs_read)) deallocate(cs_read)
    if (allocated(top_nodes)) deallocate(top_nodes)
    !
    return 
  end subroutine tDataModelData_mf6_get_data_vrt_array

  subroutine tDataModelData_mf6_get_data_vrt_array_array(this, disu, bbx, r8a)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModelData) :: this
    type(tMF6Disu), intent(in), pointer :: disu
    type(tBbx), intent(in) :: bbx
    real(R8B), dimension(:), allocatable, intent(inout) :: r8a
    ! -- local
    type(tMultiGrid), pointer :: mg   => null()
    type(tGrid),      pointer :: g    => null()
    type(tData),      pointer :: dat  => null()
    type(tVrtArray),  pointer :: vrta => null()
    type(tVrt),       pointer :: vrt  => null()
    !
    logical :: lconst
    !
    integer(I4B) :: m, il, jl, ig, ir, ic, n
    !
    real(R4B), dimension(:,:), allocatable :: r4x
    real(R4B) :: mvr4, r4const, r4v
! ------------------------------------------------------------------------------
  
    dat => this%dat(1); vrta => dat%vrta
    allocate(r8a(disu%nodes)); r8a = R8ZERO
    m = 0
    do jl = 1, disu%nlay_act
      il = disu%lay_act(jl)
      !
      mg => disu%grid_mga_nod%get_mg(jl)
      !
      ! read grid in the original resolution and buffer
      call vrta%check_constant(il, lconst, r4const)
      if (.not.lconst) then
        vrt => vrta%get_vrt(il)
        call vrt%read_extent_native_cs(bbx=bbx)
        call vrt%buffer_data()
      end if
      !
      do ig = 1, mg%ngrid
        if (mg%grid_count(ig) > 0) then
          g => mg%grid(ig)
          if (.not.lconst) then
            call vrt%read_extent(xr4=r4x, mvr4=mvr4, bbx=g%bbx, &
              i_uscl=dat%i_uscl, i_dscl=dat%i_dscl, clean_tile=.false.)
          end if
          do ir = 1, g%nr; do ic = 1, g%nc
            n = g%xi4(ic,ir)
            if (n > 0) then
              if (.not.lconst) then
                r4v = r4x(ic,ir)
              else
                r4v = r4const
              end if
              if (r4v == mvr4) then
                call errmsg('tQuad_mf6_write_data: array.')
              else
                r8a(n) = real(r4v,R8B)
                m = m + 1
              end if
            end if
          end do; end do
        end if
      end do
      !
      ! clean the vrt arrays
      if (.not.lconst) then
        call vrt%clean_x()
      end if
    end do
    !
    if (m /= disu%nodes) then
      call errmsg('tQuad_mf6_write_data: array.')
    end if
    !
    ! clean up
    if (allocated(r4x)) deallocate(r4x)
    !
    return 
  end subroutine tDataModelData_mf6_get_data_vrt_array_array
  
  subroutine tDataModelData_mf6_get_data_csv_list(this, disu, bbx, i4a, r8a)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tDataModelData) :: this
    type(tMF6Disu), intent(in), pointer :: disu
    type(tBbx), intent(in) :: bbx
    integer(I4B), dimension(:), allocatable, intent(inout) :: i4a
    real(R8B), dimension(:), allocatable, intent(inout) :: r8a
    !
    ! -- local
    type(tData),      pointer :: dat  => null()
    type(tCsv),       pointer :: csv  => null()
    type(tMultiGrid), pointer :: mg   => null()
    type(tGrid),      pointer :: g    => null()
    type(tVrtArray),  pointer :: vrta => null()
    type(tVrt),       pointer :: vrt  => null()
    !
    logical :: lfound, lconst
    !
    real(R8B), dimension(:), allocatable :: csv_x, csv_y
    !
    real(R4B), dimension(:,:), allocatable :: r4x, f2d
    real(R4B), dimension(:), allocatable :: csv_top, csv_bot, csv_v, fnod
    real(R4B) :: top, bot, gtop, gbot, f, r4const, mvr4, r4v, ftot
    !
    integer(I4B), dimension(:), allocatable :: csv_ilay, mapping
    integer(I4B), dimension(:,:), allocatable :: nod2d
    integer(I4B) :: ierr, nc, nr, nl, np, iact, ic, ir, il, jl, mp, n, m, ip, ig, i
    !
    integer(I1B), dimension(:), allocatable :: top_nodes
! ------------------------------------------------------------------------------
    if (allocated(i4a)) deallocate(i4a)
    if (allocated(r8a)) deallocate(r8a)
    !
    call disu%x_to_top_nodes(top_nodes)
    ! 
    dat => this%dat(1); csv => dat%csv
    ierr = 0
    if (this%idist <= 0) ierr = 1
    if (ierr == 1) call errmsg('Invalid idist.')
    !
    call csv%get_column(ic=1, r8a=csv_x)
    call csv%get_column(ic=2, r8a=csv_y)
    call csv%get_column(ic=this%idist, r4a=csv_v)
    !
    nc = size(disu%grid_x_nod,1)
    nr = size(disu%grid_x_nod,2)
    nl = size(disu%grid_x_nod,3)
    !
    select case(this%i_assign)
    case(i_assign_exact_layer)
      call csv%get_column(ic=this%ilay,  i4a=csv_ilay)
      !
      np = size(csv_x)
      !
      if (np == 0) then
        call errmsg('No CSV input found.')
      end if
      !
      do iact = 1, 2
        mp = 0
        do ip = 1, np
          if (point_in_bb(reshape([csv_x(ip),csv_y(ip)],[2,1]), bbx)) then
            call get_icr(ic, ir, csv_x(ip), csv_y(ip), bbx%xll, bbx%yur, disu%cs_min)
            if (valid_icr(ic, ir, nc, nr)) then
              do jl = 1, nl
                il = disu%lay_act(jl)
                n = disu%grid_x_nod(ic, ir, jl)
                if ((n > 0).and.(csv_ilay(ip) == il)) then
                  mp = mp + 1
                  if (iact == 2) then
                    i4a(mp) = n
                    r8a(mp) = real(csv_v(ip),R8B)
                  end if
                end if
              end do
            end if
          end if
        end do
        if (iact == 1) then
          if (allocated(i4a)) deallocate(i4a)
          if (allocated(r8a)) deallocate(r8a)
          if (mp > 0) then
            allocate(i4a(mp)); i4a = I4ZERO
            allocate(r8a(mp)); r8a = R8ZERO
          end if
        end if
      end do
    !
    case(i_assign_intersect)
      !
      allocate(fnod(disu%nodes))
      fnod = R4ONE
      !
      if (this%itop  <= 0) ierr = 1
      if (this%ibot  <= 0) ierr = 1
      if (ierr == 1) call errmsg('Invalid itop, ibot.')

      call csv%get_column(ic=this%itop, r4a=csv_top)
      call csv%get_column(ic=this%ibot, r4a=csv_bot)
      !
      np = size(csv_x)
      !
      if (np == 0) then
        call errmsg('No CSV input found.')
      end if
      if (allocated(i4a)) deallocate(i4a)
      allocate(mapping(np)); mapping = 0
      !
      do iact = 1, 2
        mp = 0
        do ip = 1, np
          if (point_in_bb(reshape([csv_x(ip),csv_y(ip)],[2,1]), bbx)) then
            call get_icr(ic, ir, csv_x(ip), csv_y(ip), bbx%xll, bbx%yur, disu%cs_min)
            if (valid_icr(ic, ir, nc, nr)) then
              lfound = .false.
              do il = 1, nl
                n = disu%grid_x_nod(ic, ir, il)
                if (n > 0) then
                  if (.not.lfound) then
                    mp = mp + 1
                    mapping(mp) = ip ! mapping
                    lfound = .true.
                  end if
                  if (iact == 2) then
                    top = csv_top(ip); bot = csv_bot(ip)
                    gtop = disu%top(n); gbot = disu%bot(n)
                    f = calc_frac(top, bot, top_nodes(n), gtop, gbot)
                    f2d(mp,il) = f
                    if (this%ilay0_sfac) then
                      f2d(mp,il) = f2d(mp,il)*fnod(n)*(gtop-gbot) !fraction KD
                    end if
                    nod2d(mp,il) = n
                  end if
                end if
              end do
            end if
          end if
        end do
        if (iact == 1) then
          if (allocated(f2d)) deallocate(f2d)
          if (allocated(nod2d)) deallocate(nod2d)
          if (mp > 0) then
            allocate(f2d(mp,nl)); f2d = R4ZERO
            allocate(nod2d(mp,nl)); nod2d = 0
            !
            if (this%ilay0_sfac) then
              dat => this%dat(1); vrta => dat%vrta
              call vrta%init()
              m = 0
              do jl = 1, disu%nlay_act
                il = disu%lay_act(jl)
                !
                mg => disu%grid_mga_nod%get_mg(jl)
                !
                ! read grid in the original resolution and buffer
                call vrta%check_constant(il, lconst, r4const)
                if (.not.lconst) then
                  vrt => vrta%get_vrt(il)
                  call vrt%read_extent_native_cs(bbx=bbx)
                  call vrt%buffer_data()
                end if
                !
                do ig = 1, mg%ngrid
                  if (mg%grid_count(ig) > 0) then
                    g => mg%grid(ig)
                    if (.not.lconst) then
                      call vrt%read_extent(xr4=r4x, mvr4=mvr4, bbx=g%bbx, &
                        i_uscl=dat%i_uscl, i_dscl=dat%i_dscl, clean_tile=.false.)
                    end if
                    do ir = 1, g%nr; do ic = 1, g%nc
                      n = g%xi4(ic,ir)
                      if (n > 0) then
                        if (.not.lconst) then
                          r4v = r4x(ic,ir)
                        else
                          r4v = r4const
                        end if
                        if (r4v == mvr4) then
                          call errmsg('tQuad_mf6_write_data: array.')
                        else
                          fnod(n) = r4v
                          m = m + 1
                        end if
                      end if
                    end do; end do
                  end if
                end do
                !
                ! clean the vrt arrays
                if (.not.lconst) then
                  call vrt%clean_x()
                end if
              end do
              !
              if (m /= disu%nodes) then
                call errmsg('tQuad_mf6_write_data: array.')
              end if
              !
            end if
          end if
        end if
      end do
      !
      ! scale the fractions to 1
      do i = 1, mp
        ftot = sum(f2d(i,:))
        do il = 1, nl
          f2d(i,il) = f2d(i,il)/(ftot+R4TINY)
        end do
      end do
      if (.false.) then
        do i = 1, mp
          write(*,*) 'sum frac: ',sum(f2d(i,:))
        end do
      end if
      !
      ! count and buffer
      do iact = 1, 2
        n = 0
        do il = 1, nl
          do i = 1, mp
            f = f2d(i,il)
            if (f > R4ZERO) then
              ip = mapping(i)
              if (csv_v(ip) /= R4ZERO) then
                n = n + 1
                if (iact == 2) then
                  i4a(n) = nod2d(i,il)
                  ip = mapping(i)
                  r8a(n) = f*csv_v(ip)
                end if
              end if
            end if
          end do
        end do
        if (iact == 1) then
          if (n > 0) then
            if (allocated(i4a)) deallocate(i4a)
            if (allocated(r8a)) deallocate(r8a)
            allocate(i4a(n), r8a(n))
          end if
        end if
      end do
      !
    case(i_assign_first_layer)
      call errmsg('Assigned to first layer is not supported for CSV input.')
    end select
    !
    ! clean up
    if (allocated(top_nodes)) deallocate(top_nodes)
    if (allocated(csv_x))     deallocate(csv_x)
    if (allocated(csv_y))     deallocate(csv_y)
    if (allocated(csv_v))     deallocate(csv_v)
    if (allocated(csv_ilay))  deallocate(csv_ilay)
    if (allocated(fnod))      deallocate(fnod)
    if (allocated(csv_top))   deallocate(csv_top)
    if (allocated(csv_bot))   deallocate(csv_bot)
    if (allocated(mapping))   deallocate(mapping)
    if (allocated(f2d))       deallocate(f2d)
    if (allocated(nod2d))     deallocate(nod2d)
    if (allocated(r4x))       deallocate(r4x)
    !
    return
  end subroutine tDataModelData_mf6_get_data_csv_list
  
  subroutine tQuad_mf6_write_data(this)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    ! -- local
    type(tBB) :: bbi
    type(tBbx) :: bbx
    type(tMultiGridArray), allocatable :: mga_read
    type(tMultiGrid), pointer :: mg      => null()
    type(tMultiGrid), pointer :: mg_top  => null()
    type(tMultiGrid), pointer :: mg_bot  => null()
    type(tMultiGrid), pointer :: mg_dist => null()
    type(tGrid),      pointer :: g       => null()
    type(tDataModelData), pointer :: dmdat => null()
    type(tData), pointer :: dat => null()
    type(tMF6disu), pointer :: disu => null()
    type(tCsv), pointer :: csv => null()
    type(tMf6Wbd), pointer :: wbd => null()
    type(tVrtArray), pointer :: vrta => null()
    type(tVrt), pointer :: vrt => null()
    !
    logical :: lconst
    !
    integer(I4B), dimension(:), allocatable :: csv_ilay
    real(R8B), dimension(:), allocatable :: csv_x, csv_y
    real(R4B), dimension(:), allocatable :: csv_top, csv_bot, csv_v
    !
    real(R8B), dimension(:,:), allocatable :: r8x
    real(R8B), dimension(:), allocatable :: r8a, cs, cs_read
    real(R8B) :: cs_min_tgt, cs_min_rea
    !
    real(R4B), dimension(:,:), allocatable :: r4x, f2d
    real(R4B), dimension(:), allocatable :: f1d, fnod
    real(R4B) :: mvr4, top, bot, dist, gtop, gbot, dz, z1, z2, r4v, f, ftot
    real(R4B) :: r4const
    !
    integer(I4B), dimension(:), allocatable :: nod1d, i4a, mapping
    integer(I4B), dimension(:,:), allocatable :: nod2d
    integer(I4B) :: ndat, idat, ilay, ig, i, ip, ierr
    integer(I4B) :: il, jl, ir, ic, np, mp, nl, nr, nc
    integer(I4B) :: n, m, ncell, nr8a, iact
    !
    integer(I1B), dimension(:), allocatable :: top_nodes
    !
    logical :: lfound, lwrite, ltopisbot
! ------------------------------------------------------------------------------
    call this%grid_init(); disu => this%disu; wbd => disu%wbd
    call disu%x_to_top_nodes(top_nodes)
    !
    call this%get_bb(child_bbi=bbi, child_bbx=bbx)
    call this%get_prop_csv(key='cs_min_rea', r8v=cs_min_rea) !25-08
    cs_min_tgt = min(bbx%cs, cs_min_rea)
    !
    allocate(cs_read(disu%ngrid))
    cs_read(disu%ngrid) = cs_min_tgt
    do ig = disu%ngrid-1, 1, -1
      cs_read(ig) = cs_read(ig+1)*2.d0
    end do
    !
    do idat = 1, this%dat_mod%get_ndat()
      call this%dat_mod%get_dat(idat, dmdat)
      !
      select case(dmdat%i_in_file_type)
      case(i_vrt)
        !
        select case(dmdat%i_type)
        case(i_type_list)
          !
          allocate(mga_read)
          call mga_read%init(dmdat%ndat)
          !
          ! step 1: read the data for *ALL* the ngrid levels
          do i = 1, dmdat%ndat
            dat => dmdat%dat(i); mg => mga_read%mga(i)
            !
            if (dat%i_type == i_vrt) then
              call vrt_read_extent_mg(vrt=dat%vrt, bbx=bbx, csa=cs_read, &
                i_uscl=dat%i_uscl, i_dscl=dat%i_dscl, mg=mg, mvr4=mvr4)
            end if
          end do
          !
          allocate(f1d(disu%nodes)); f1d = R4ZERO
          !
          select case(dmdat%i_assign)
          case(i_assign_intersect)
            ierr = 0
            if ((dmdat%itop  <= 0).or.(dmdat%itop  > dmdat%ndat)) ierr = 1
            if ((dmdat%ibot  <= 0).or.(dmdat%ibot  > dmdat%ndat)) ierr = 1
            if ((dmdat%idist <= 0).or.(dmdat%idist > dmdat%ndat)) ierr = 1
            if (ierr == 1) call errmsg('Invalid itop, ibot or idist.')
            !
            !if ((dmdat%ndat == 2).and.(dmdat%itop == dmdat%ibot)) then
            !  ltopisbot = .true.
            !else
            !  ltopisbot = .false.
            !end if
            !
            mg_top  => mga_read%mga(dmdat%itop)
            mg_bot  => mga_read%mga(dmdat%ibot)
            mg_dist => mga_read%mga(dmdat%idist)
            !
            do jl = 1, disu%nlay_act
               mg => disu%grid_mga_nod%get_mg(jl)
               do ig = 1, mg%ngrid
                 g => mg%grid(ig)
                 do ir = 1, g%nr; do ic = 1, g%nc
                   n = g%xi4(ic,ir)
                   if (n /= g%mvi4) then
                     top  = mg_top%grid(ig)%xr4(ic,ir)
                     bot  = mg_bot%grid(ig)%xr4(ic,ir)
                     dist = mg_dist%grid(ig)%xr4(ic,ir)
                     !
                     if ((top  /= mg_top%grid(ig)%mvr4).and. &
                         (bot  /= mg_bot%grid(ig)%mvr4).and. &
                         (dist /= mg_dist%grid(ig)%mvr4)) then
                       !
                       ! make consistent
                       if (bot > top) then
                         call logmsg('Warning: inconsistent top/bot found, setting bot equal to top.')
                         bot = top
                         mg_bot%grid(ig)%xr4(ic,ir) = bot
                       end if
                       !
                       if (top == bot) then
                         ltopisbot = .true.
                       else
                         ltopisbot = .false.
                       end if
                       !
                       gtop = disu%grid_mga_top%mga(jl)%grid(ig)%xr4(ic,ir)
                       gbot = disu%grid_mga_bot%mga(jl)%grid(ig)%xr4(ic,ir)
                       !
                       if (ltopisbot) then
                         if (top_nodes(n) == 1) then
                           if (top >= gbot) then
                             f1d(n) = R4ONE
                           end if
                         else
                           if ((top >= gbot).and.(top < gtop)) then
                             f1d(n) = R4ONE
                           end if
                         end if
                       else
                         dz = top - bot
                         if (dz >= R4ZERO) then
                           f1d(n) = calc_frac(top, bot, top_nodes(n), gtop, gbot)
                           !
                           ! overrule in case of top node
                           if (top_nodes(n) == 1) then
                             if (bot >= gtop) then
                               f1d(n) = R4ONE
                             end if
                           end if
                         end if
                       end if
                     end if
                   end if
                 end do; end do
               end do
            end do
          case(i_assign_exact_layer)
            !
            do jl = 1, disu%nlay_act
              il = disu%lay_act(jl)
               mg => disu%grid_mga_nod%get_mg(jl)
               do ig = 1, mg%ngrid
                 g => mg%grid(ig)
                 do ir = 1, g%nr; do ic = 1, g%nc
                   n = g%xi4(ic,ir)
                   if (n /= g%mvi4) then
                     if (il == dmdat%ilay) then
                       f1d(n) = R4ONE
                     end if
                   end if
                 end do; end do
               end do
            end do
          case(i_assign_first_layer)
            !
            do jl = 1, disu%nlay_act
               mg => disu%grid_mga_nod%get_mg(jl)
               do ig = 1, mg%ngrid
                 g => mg%grid(ig)
                 do ir = 1, g%nr; do ic = 1, g%nc
                   n = g%xi4(ic,ir)
                   if (n /= g%mvi4) then
                     if (top_nodes(n) == 1) then
                       f1d(n) = R4ONE
                     end if
                   end if
                 end do; end do
               end do
            end do
          end select
          !
          ! set fraction to zero when dist is zero
          if (.true.) then
            ncell = 0
            do jl = 1, disu%nlay_act
               mg => disu%grid_mga_nod%get_mg(jl)
               do ig = 1, mg%ngrid
                 g => mg%grid(ig)
                 do ir = 1, g%nr; do ic = 1, g%nc
                   n = g%xi4(ic,ir)
                   if (n /= g%mvi4) then
                     f = f1d(n)
                     if (f > R4ZERO) then
                       i = dmdat%idist
                       r4v = mga_read%mga(i)%grid(ig)%xr4(ic,ir)
                       if (r4v == R4ZERO) then
                         f1d(n) = R4ZERO
                         ncell = ncell + 1
                       end if
                     end if
                   end if
                end do; end do
              end do
            end do
            if (ncell > 0) then
              call logmsg('Skipping '//ta([ncell])//' cells with zero value...')
            end if
          end if
          !
          ! count using the fractions and allocate
          ncell = 0
          do i = 1, disu%nodes
            if (f1d(i) > R4ZERO) ncell = ncell + 1
          end do
          if (ncell == 0) then
            call logmsg('No data found...')
            ! clean up
            if (allocated(mga_read)) then
              call mga_read%clean()
              deallocate(mga_read)
            end if
            if (allocated(i4a)) deallocate(i4a)
            if (allocated(i4a)) deallocate(i4a)
            
            if (allocated(r8x)) deallocate(r8x)
            if (allocated(f1d)) deallocate(f1d)
            cycle
          end if
          nr8a = dmdat%ndat
          allocate(i4a(ncell), r8x(nr8a,ncell))
          !
          ! buffer the data
          m = 0
          do jl = 1, disu%nlay_act
             mg => disu%grid_mga_nod%get_mg(jl)
             do ig = 1, mg%ngrid
               g => mg%grid(ig)
               do ir = 1, g%nr; do ic = 1, g%nc
                 n = g%xi4(ic,ir)
                 if (n /= g%mvi4) then
                   f = f1d(n)
                   if (f > R4ZERO) then
                     m = m + 1
                     i4a(m) = n
                     do i = 1, nr8a
                       r4v = mga_read%mga(i)%grid(ig)%xr4(ic,ir)
                       if (i == dmdat%idist) then
                         r4v = r4v*f
                       end if
                       r8x(i,m) = real(r4v, R8B)
                     end do
                   end if
                 end if
               end do; end do
             end do
          end do
          !
          ! check the consistency top/bot
          if (dmdat%itop /= dmdat%ibot) then
            do m = 1, ncell
              if (r8x(dmdat%ibot,m) > r8x(dmdat%itop,m)) then
                call errmsg('Inconsistent data: bot > top.')
              end if
            end do
          end if
          !
          select case(dmdat%i_out_file_type)
          case(i_binpos)
            call wbd%write_list(id=dmdat%id, i4a1=i4a, r8x=r8x)
          case(i_bin)
            call wbd%write_list(id=dmdat%id, i4a1=i4a, r8x=r8x, &
              f_bin=trim(this%mod_dir)//'\'//trim(dmdat%id)//'.bin')
          case(i_asc)
            call wbd%write_list(id=dmdat%id, i4a1=i4a, r8x=r8x, &
              f_asc=trim(this%mod_dir)//'\'//trim(dmdat%id)//'.asc')
          end select
        case(i_type_array)
          if (dmdat%i_assign == i_assign_first_layer) then
            !
            allocate(mga_read)
            call mga_read%init(dmdat%ndat)
            !
            ! step 1: read the data for *ALL* the ngrid levels
            do i = 1, dmdat%ndat
              dat => dmdat%dat(i); mg => mga_read%mga(i)
              !
              if (dat%i_type == i_vrt) then
                call vrt_read_extent_mg(vrt=dat%vrt, bbx=bbx, csa=cs_read, &
                  i_uscl=dat%i_uscl, i_dscl=dat%i_dscl, mg=mg, mvr4=mvr4)
              end if
            end do
            !
            allocate(r8a(disu%nodes)); r8a = dmdat%r8const
            !
            do jl = 1, disu%nlay_act
               mg => disu%grid_mga_nod%get_mg(jl)
               do ig = 1, mg%ngrid
                 g => mg%grid(ig)
                 do ir = 1, g%nr; do ic = 1, g%nc
                   n = g%xi4(ic,ir)
                   if (n /= g%mvi4) then
                     if (top_nodes(n) == 1) then
                       r4v = mga_read%mga(1)%grid(ig)%xr4(ic,ir)
                       r8a(n) = real(r4v,R8B)
                     end if
                   end if
                 end do; end do
               end do
            end do
            !
            select case(dmdat%i_out_file_type)
            case(i_binpos)
              call wbd%write_array(id=dmdat%id, r8a=r8a)
            case(i_bin)
              call wbd%write_array(id=dmdat%id, r8a=r8a, &
                f_bin=trim(this%mod_dir)//'\'//trim(dmdat%id)//'.bin')
            case(i_asc)
              call wbd%write_array(id=dmdat%id, r8a=r8a, &
                f_asc=trim(this%mod_dir)//'\'//trim(dmdat%id)//'.asc')
            end select
          else
            call errmsg('Array output not supported for vrt input.')
          end if
        end select
      case(i_vrt_array)
        select case(dmdat%i_type)
        case(i_type_list)
          call errmsg('List output not supported for vrt array input.')
        case(i_type_array)
          dat => dmdat%dat(1); vrta => dat%vrta
          allocate(r8a(disu%nodes)); r8a = R8ZERO
          m = 0
          do jl = 1, disu%nlay_act
            il = disu%lay_act(jl)
            !
            mg => disu%grid_mga_nod%get_mg(jl)
            !
            ! read grid in the original resolution and buffer
            call vrta%check_constant(il, lconst, r4const)
            if (.not.lconst) then
              vrt => vrta%get_vrt(il)
              call vrt%read_extent_native_cs(bbx=bbx)
              call vrt%buffer_data()
            end if
            !
            do ig = 1, mg%ngrid
              if (mg%grid_count(ig) > 0) then
                g => mg%grid(ig)
                if (.not.lconst) then
                  call vrt%read_extent(xr4=r4x, mvr4=mvr4, bbx=g%bbx, &
                    i_uscl=dat%i_uscl, i_dscl=dat%i_dscl, clean_tile=.false.)
                end if
                do ir = 1, g%nr; do ic = 1, g%nc
                  n = g%xi4(ic,ir)
                  if (n > 0) then
                    if (.not.lconst) then
                      r4v = r4x(ic,ir)
                    else
                      r4v = r4const
                    end if
                    if (r4v == mvr4) then
                      call errmsg('tQuad_mf6_write_data: array.')
                    else
                      r8a(n) = real(r4v,R8B)
                      m = m + 1
                    end if
                  end if
                end do; end do
              end if
            end do
            !
            ! clean the vrt arrays
            if (.not.lconst) then
              call vrt%clean_x()
            end if
          end do
          !
          if (m /= disu%nodes) then
            call errmsg('tQuad_mf6_write_data: array.')
          end if
          !
          select case(dmdat%i_out_file_type)
          case(i_binpos)
            call wbd%write_array(id=dmdat%id, r8a=r8a)
          case(i_bin)
            call wbd%write_array(id=dmdat%id, r8a=r8a, &
              f_bin=trim(this%mod_dir)//'\'//trim(dmdat%id)//'.bin')
          case(i_asc)
            call wbd%write_array(id=dmdat%id, r8a=r8a, &
              f_asc=trim(this%mod_dir)//'\'//trim(dmdat%id)//'.asc')
          end select
          !
        end select
      case(i_csv)
        select case(dmdat%i_type)
        case(i_type_list)
          dat => dmdat%dat(1); csv => dat%csv
          ierr = 0
          if (dmdat%idist <= 0) ierr = 1
          if (ierr == 1) call errmsg('Invalid idist.')
          !
          call csv%get_column(ic=1, r8a=csv_x)
          call csv%get_column(ic=2, r8a=csv_y)
          call csv%get_column(ic=dmdat%idist, r4a=csv_v)
          !
          nc = size(disu%grid_x_nod,1)
          nr = size(disu%grid_x_nod,2)
          nl = size(disu%grid_x_nod,3)
          !
          select case(dmdat%i_assign)
          case(i_assign_exact_layer)
            call csv%get_column(ic=dmdat%ilay,  i4a=csv_ilay)
            !
            np = size(csv_x)
            !
            if (np == 0) then
              call errmsg('No CSV input found.')
            end if
            !
            do iact = 1, 2
              mp = 0
              do ip = 1, np
                if (point_in_bb(reshape([csv_x(ip),csv_y(ip)],[2,1]), bbx)) then
                  call get_icr(ic, ir, csv_x(ip), csv_y(ip), bbx%xll, bbx%yur, disu%cs_min)
                  if (valid_icr(ic, ir, nc, nr)) then
                    do jl = 1, nl
                      il = disu%lay_act(jl)
                      n = disu%grid_x_nod(ic, ir, jl)
                      if ((n > 0).and.(csv_ilay(ip) == il)) then
                        mp = mp + 1
                        if (iact == 2) then
                          i4a(mp) = n
                          r8a(mp) = real(csv_v(ip),R8B)
                        end if
                      end if
                    end do
                  end if
                end if
              end do
              if (iact == 1) then
                if (allocated(i4a)) deallocate(i4a)
                if (allocated(r8a)) deallocate(r8a)
                if (mp > 0) then
                  allocate(i4a(mp)); i4a = I4ZERO
                  allocate(r8a(mp)); r8a = R8ZERO
                end if
              end if
            end do
          !
          case(i_assign_intersect)
            !
            allocate(fnod(disu%nodes))
            fnod = R4ONE
            !
            if (dmdat%itop  <= 0) ierr = 1
            if (dmdat%ibot  <= 0) ierr = 1
            if (ierr == 1) call errmsg('Invalid itop, ibot.')

            call csv%get_column(ic=dmdat%itop, r4a=csv_top)
            call csv%get_column(ic=dmdat%ibot, r4a=csv_bot)
            !
            np = size(csv_x)
            !
            if (np == 0) then
              call errmsg('No CSV input found.')
            end if
            if (allocated(i4a)) deallocate(i4a)
            allocate(mapping(np)); mapping = 0
            !
            do iact = 1, 2
              mp = 0
              do ip = 1, np
                if (point_in_bb(reshape([csv_x(ip),csv_y(ip)],[2,1]), bbx)) then
                  call get_icr(ic, ir, csv_x(ip), csv_y(ip), bbx%xll, bbx%yur, disu%cs_min)
                  if (valid_icr(ic, ir, nc, nr)) then
                    lfound = .false.
                    do il = 1, nl
                      n = disu%grid_x_nod(ic, ir, il)
                      if (n > 0) then
                        if (.not.lfound) then
                          mp = mp + 1
                          mapping(mp) = ip ! mapping
                          lfound = .true.
                        end if
                        if (iact == 2) then
                          top = csv_top(ip); bot = csv_bot(ip)
                          gtop = disu%top(n); gbot = disu%bot(n)
                          f = calc_frac(top, bot, top_nodes(n), gtop, gbot)
                          f2d(mp,il) = f
                          if (dmdat%ilay0_sfac) then
                            f2d(mp,il) = f2d(mp,il)*fnod(n)*(gtop-gbot) !fraction KD
                          end if
                          nod2d(mp,il) = n
                        end if
                      end if
                    end do
                  end if
                end if
              end do
              if (iact == 1) then
                if (allocated(f2d)) deallocate(f2d)
                if (allocated(nod2d)) deallocate(nod2d)
                if (mp > 0) then
                  allocate(f2d(mp,nl)); f2d = R4ZERO
                  allocate(nod2d(mp,nl)); nod2d = 0
                  !
                  if (dmdat%ilay0_sfac) then
                    dat => dmdat%dat(1); vrta => dat%vrta
                    call vrta%init()
                    m = 0
                    do jl = 1, disu%nlay_act
                      il = disu%lay_act(jl)
                      !
                      mg => disu%grid_mga_nod%get_mg(jl)
                      !
                      ! read grid in the original resolution and buffer
                      call vrta%check_constant(il, lconst, r4const)
                      if (.not.lconst) then
                        vrt => vrta%get_vrt(il)
                        call vrt%read_extent_native_cs(bbx=bbx)
                        call vrt%buffer_data()
                      end if
                      !
                      do ig = 1, mg%ngrid
                        if (mg%grid_count(ig) > 0) then
                          g => mg%grid(ig)
                          if (.not.lconst) then
                            call vrt%read_extent(xr4=r4x, mvr4=mvr4, bbx=g%bbx, &
                              i_uscl=dat%i_uscl, i_dscl=dat%i_dscl, clean_tile=.false.)
                          end if
                          do ir = 1, g%nr; do ic = 1, g%nc
                            n = g%xi4(ic,ir)
                            if (n > 0) then
                              if (.not.lconst) then
                                r4v = r4x(ic,ir)
                              else
                                r4v = r4const
                              end if
                              if (r4v == mvr4) then
                                call errmsg('tQuad_mf6_write_data: array.')
                              else
                                fnod(n) = r4v
                                m = m + 1
                              end if
                            end if
                          end do; end do
                        end if
                      end do
                      !
                      ! clean the vrt arrays
                      if (.not.lconst) then
                        call vrt%clean_x()
                      end if
                    end do
                    !
                    if (m /= disu%nodes) then
                      call errmsg('tQuad_mf6_write_data: array.')
                    end if
                    !
                  end if
                end if
              end if
            end do
            !
            ! scale the fractions to 1
            do i = 1, mp
              ftot = sum(f2d(i,:))
              do il = 1, nl
                f2d(i,il) = f2d(i,il)/(ftot+R4TINY)
              end do
            end do
            if (.false.) then
              do i = 1, mp
                write(*,*) 'sum frac: ',sum(f2d(i,:))
              end do
            end if
            !
            ! count and buffer
            do iact = 1, 2
              n = 0
              do il = 1, nl
                do i = 1, mp
                  f = f2d(i,il)
                  if (f > R4ZERO) then
                    ip = mapping(i)
                    if (csv_v(ip) /= R4ZERO) then
                      n = n + 1
                      if (iact == 2) then
                        i4a(n) = nod2d(i,il)
                        ip = mapping(i)
                        r8a(n) = f*csv_v(ip)
                      end if
                    end if
                  end if
                end do
              end do
              if (iact == 1) then
                if (n > 0) then
                  if (allocated(i4a)) deallocate(i4a)
                  if (allocated(r8a)) deallocate(r8a)
                  allocate(i4a(n), r8a(n))
                end if
              end if
            end do
            !
          case(i_assign_first_layer)
            call errmsg('Assigned to first layer is not supported for CSV input.')
          end select
          !
          if (allocated(r8a)) then
            lwrite = .true.
          else
            lwrite = .false.
          end if
          !
          if (lwrite) then
            select case(dmdat%i_out_file_type)
            case(i_binpos)
              call wbd%write_list(id=dmdat%id, i4a1=i4a, r8a1=r8a)
            case(i_bin)
              call wbd%write_list(id=dmdat%id, i4a1=i4a, r8a1=r8a, &
                f_bin=trim(this%mod_dir)//'\'//trim(dmdat%id)//'.bin')
            case(i_asc)
              call wbd%write_list(id=dmdat%id, i4a1=i4a, r8a1=r8a, &
                f_asc=trim(this%mod_dir)//'\'//trim(dmdat%id)//'.asc')
            end select
          end if
          !
        case(i_type_array)
          call errmsg('Array output is not supported for CSV input.')
        end select
      case default
        call errmsg('tQuad_mf6_write_data')
      end select
      !
      ! clean up
      if (allocated(mga_read)) then
        call mga_read%clean()
        deallocate(mga_read)
      end if
      if (allocated(i4a))     deallocate(i4a)
      if (allocated(r4x))     deallocate(r4x)
      if (allocated(r8a))     deallocate(r8a)
      if (allocated(r8x))     deallocate(r8x)
      if (allocated(mapping)) deallocate(mapping)
      if (allocated(nod1d))   deallocate(nod1d)
      if (allocated(nod2d))   deallocate(nod2d)
      if (allocated(csv_x))   deallocate(csv_x)
      if (allocated(csv_y))   deallocate(csv_y)
      if (allocated(csv_top)) deallocate(csv_top)
      if (allocated(csv_bot)) deallocate(csv_bot)
      if (allocated(csv_v))   deallocate(csv_v)
      if (allocated(f1d))     deallocate(f1d)
      if (allocated(f2d))     deallocate(f2d)
      if (allocated(fnod))    deallocate(fnod)
    end do
    !
    ! write the csv-file
    call wbd%write_csv()
    !
    ! clean up
    call disu%clean()
    deallocate(this%disu); this%disu => null()
    if (allocated(top_nodes)) deallocate(top_nodes)
    if (allocated(cs_read)) deallocate(cs_read)
    !
    return
  end subroutine tQuad_mf6_write_data
  
  function calc_frac(t1, b1, tn1, t2, b2) result(f)
! ******************************************************************************
!
!    SPECIFICATIONS:
!    Calculate the fraction of 1 in 2.
! ------------------------------------------------------------------------------
    ! -- dummy
    real(R4B), intent(in)    :: t1
    real(R4B), intent(in)    :: b1
    integer(I1B), intent(in) :: tn1
    real(R4B), intent(in)    :: t2
    real(R4B), intent(in)    :: b2
    real(R4B) :: f
    ! -- local
    real(R4B) :: dz1, dz2
! ------------------------------------------------------------------------------
    !
    dz1 = t1 - b1
    !
    if (b1 >= t2) then ! I
      f = R4ZERO
      return
    end if
    if (t1 <= b2) then ! II
      f = R4ZERO
      return
    end if
    if ((t1 >= t2).and.(b1 <= b2)) then ! III
      if (tn1 == 1) then
        dz2 = t1 - b2
      else
        dz2 = t2 - b2
      end if
      f = dz2/dz1
      return
    end if
    if ((t1 <= t2).and.(b1 >= b2)) then ! IV
      f = R4ONE
      return
    end if
    if ((t1 > t2).and.(b1 < t2).and.(b1 > b2)) then ! V
      if (tn1 == 1) then
        dz2 = t1 - b1
      else
        dz2 = t2 - b1
      end if
      f = dz2/dz1
      return
    end if
    if ((b1 < b2).and.(t1 < t2).and.(t1 > b2)) then ! VI
      dz2 = t1 - b2
      f = dz2/dz1
      return
    end if
    !
    return
  end function calc_frac
  
  subroutine vrt_read_extent_mg(vrt, bbx, csa, i_uscl, i_dscl, mg, &
    mvi1, mvi2, mvi4, mvi8, mvr4, mvr8)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    type(tVrt), intent(in) :: vrt
    !
    type(tBbx), intent(in) :: bbx
    real(R8B), dimension(:) :: csa
    integer(I4B), intent(in), optional :: i_uscl
    integer(I4B), intent(in), optional :: i_dscl
    type(tMultiGrid), intent(inout) :: mg
    !
    integer(I1B), intent(out), optional :: mvi1
    integer(I2B), intent(out), optional :: mvi2
    integer(I4B), intent(out), optional :: mvi4
    integer(I8B), intent(out), optional :: mvi8
    real(R4B),    intent(out), optional :: mvr4
    real(R8B),    intent(out), optional :: mvr8
    !
    ! -- local
    type(tBbx) :: bbx_read
    type(tGrid), pointer :: g   => null()
    
    integer(I4B) :: ngrid, n, ig
! ------------------------------------------------------------------------------
    !
    ngrid = size(csa)
    call mg%init(ngrid, xll=bbx%xll, xur=bbx%xur, &
      yll=bbx%yll, yur=bbx%yur, csa=csa, &
      mvi1=mvi1, mvi2=mvi2, mvi4=mvi4, mvi8=mvi8, mvr4=mvr4, mvr8=mvr8)
    !
    ! read the original grid only
    bbx_read = bbx; bbx_read%cs = vrt%dst_bbx%cs
    call vrt%read_extent(bbx=bbx_read, &
       i_uscl=i_uscl_nodata, i_dscl=i_dscl_nodata, clean_tile=.false.)
    !
    ! buffer the data read
    call vrt%buffer_data()
    !
    ! up- or downscaling only
    do ig = 1, ngrid
      g => mg%grid(ig)
      if ((g%i_data_type /= i_r4).and.(g%i_data_type /= i_r8)) then
        call errmsg('vrt_read_extent_mg: unsupported data type.')
      end if
      if (g%i_data_type == i_r4) then
        call vrt%read_extent(xr4=g%xr4, mvr4=g%mvr4, bbx=g%bbx, &
          i_uscl=i_uscl, i_dscl=i_dscl, clean_tile=.false.)
        mvr4 = g%mvr4
      else
        call vrt%read_extent(xr8=g%xr8, mvr8=g%mvr8, bbx=g%bbx, &
          i_uscl=i_uscl, i_dscl=i_dscl, clean_tile=.false.)
        mvr8 = g%mvr8
      end if
    end do
    !
    ! set the grid count
    call mg%set_grid_count()
    !
    ! clean
    call vrt%clean()
    !
    return
  end subroutine vrt_read_extent_mg
  
  subroutine tQuad_set_flag(this, mapped_gids, active, work_done, remove)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    logical, intent(in), optional :: mapped_gids
    logical, intent(in), optional :: active
    logical, intent(in), optional :: work_done
    logical, intent(in), optional :: remove
    ! -- local
! ------------------------------------------------------------------------------
    if (present(mapped_gids)) this%mapped_gids = mapped_gids
    if (present(active))      this%active      = active
    if (present(work_done))   this%work_done   = work_done
    !
    return
  end subroutine tQuad_set_flag
  
  function tQuad_get_flag(this, mapped_gids, active, work_done, &
    generated, has_bnd) result(flag)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    logical, intent(in), optional :: mapped_gids
    logical, intent(in), optional :: active
    logical, intent(in), optional :: work_done
    logical, intent(in), optional :: generated
    logical, intent(in), optional :: has_bnd
    logical :: flag
    ! -- local
    integer :: i
    type(tIntf), pointer :: intf => null()
! ------------------------------------------------------------------------------
    i = 0
    if (present(mapped_gids)) then
      flag = this%mapped_gids; i = i + 1
    end if
    if (present(active)) then
      flag = this%active; i = i + 1
    end if
    if (present(work_done)) then
      flag = this%work_done; i = i + 1
    end if
    if (present(generated)) then
      flag = this%generated; i = i + 1
    end if
    if (present(has_bnd)) then
      intf => this%intf
      if (.not.associated(intf)) then
        call errmsg('tQuad_get_flag: flag has_bnd could not be set.')
      end if
      if (intf%n_mv > 0) then
        flag = .true.
      else
        flag = .false.
      end if
      i = i + 1
    end if
    !
    if ((i == 0).or.(i > 1)) then
      call errmsg('tQuad_get_flag: flag could not be set.')
    end if
    !
    return
  end function tQuad_get_flag
    
  subroutine tQuad_get_bb(this, prent_bbi, prent_bbx, child_bbi, child_bbx)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    type(tBb),  intent(out), optional :: prent_bbi
    type(tBbX), intent(out), optional :: prent_bbx
    type(tBb),  intent(out), optional :: child_bbi
    type(tBbX), intent(out), optional :: child_bbx
    ! -- local
! ------------------------------------------------------------------------------
    if (present(prent_bbi)) then
      prent_bbi = this%bbo%prent_bbi
    end if
    if (present(prent_bbx)) then
      prent_bbx = this%bbo%prent_bbx
    end if
    if (present(child_bbi)) then
      child_bbi = this%bbo%child_bbi
    end if
    if (present(child_bbx)) then
      child_bbx = this%bbo%child_bbx
    end if
    !
    return
  end subroutine tQuad_get_bb

!  function tQuad_get_grid_count(this) result(grid_count)
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- dummy
!    class(tQuad) :: this
!    integer(I4B), dimension(:), allocatable :: grid_count
!    ! -- local
!! ------------------------------------------------------------------------------
!    !
!    allocate(grid_count, source=this%grid_count)
!    !
!    return
!  end function tQuad_get_grid_count
  
  function tQuad_get_mo_nbr(this) result(id)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(tQuad) :: this
    integer(I4B) :: id
    ! -- local
! ------------------------------------------------------------------------------
    id = this%intf%mo_nbr
    !
    return
  end function tQuad_get_mo_nbr

! ==============================================================================
! ==============================================================================
! GENERIC MODULE ROUTINES
! ==============================================================================
! ==============================================================================

  function get_number_of_levels(cs_src, cs_dst) result(nlev)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    real(R8B), intent(in) :: cs_src
    real(R8B), intent(in) :: cs_dst
    integer(I4B) :: nlev
! ------------------------------------------------------------------------------
  
    if (cs_dst > cs_src) then
      nlev = log(cs_dst/cs_src)/log(2.d0) + 1
      nlev = -nlev
    else
      nlev = log(cs_src/cs_dst)/log(2.d0) + 1
    end if
    !
    return
  end function get_number_of_levels
  
  function get_refinement_level(nlev) result(ref)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I4B), intent(in) :: nlev
    real(R8B) :: ref
! ------------------------------------------------------------------------------
    ref = 2**(real(abs(nlev),R8B)-1)
    if (nlev > 0) then
      ref = R8ONE/ref
    end if
    !
    return
  end function get_refinement_level
  
  subroutine label_int_bnd_i1(x, x_mv, x_src, x_tgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer(I1B), dimension(:,:), intent(inout) :: x
    integer(I1B), intent(in) :: x_mv, x_src, x_tgt
    ! -- local
    logical :: lbnd
    integer(I4B) :: nc, nr, ir, ic, jr, jc, kr, kc, ist, n
! ------------------------------------------------------------------------------
    !
    !return
    nc = size(x,1); nr = size(x,2)
    n = 0
    !
    do ir = 1, nr
      do ic = 1, nc
        ! get 9-point stencil
        sticir_np = reshape([ic, ir, ic+1, ir, ic-1, ir, ic, ir+1, ic, ir-1, &
        ic+1, ir-1, ic-1, ir-1, ic+1, ir+1, ic-1, ir+1], shape(sticir_np))
        !
        ! check if stencil exceeds bounds
        stl_np = .true.
        do ist = 1, nst_np
          if ((sticir_np(1,ist) < 1).or.(sticir_np(1,ist) > nc)) stl_np(ist) = .false.
          if ((sticir_np(2,ist) < 1).or.(sticir_np(2,ist) > nr)) stl_np(ist) = .false.
        end do
        ! apply stencil
        if (stl_np(i_p)) then
          jc = sticir_np(1,i_p); jr = sticir_np(2,i_p)
          if (x(jc,jr) == x_src) then
            lbnd = .false.
            do ist = 2, nst_np
              if (stl_np(ist)) then
                kc = sticir_np(1,ist); kr = sticir_np(2,ist)
                if ((x(kc,kr) == x_mv)) lbnd = .true.
              end if
            end do
            if (lbnd) then
              x(jc,jr) = x_tgt
              n = n + 1
            end if
          end if
        end if
      end do
    end do
    !
    call logmsg('# boundary cells found: '//ta([n]))
    !
    return
  end subroutine label_int_bnd_i1
  
end module quad2dModule

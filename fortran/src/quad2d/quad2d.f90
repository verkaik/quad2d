module main_module
  ! -- modules
  use utilsmod, only: MXSLEN, MXSLENLONG, I1B, I4B, I8B, R4B, R8B, tBb, tBbX, tBbObj, ta, logmsg, &
    errmsg, I4ZERO, R4ZERO, I4MINONE, I8ZERO, I8ONE, DHALF, quicksort_r, bbi_intersect, node_to_icrl, &
    open_file, get_xy, get_icr, R8ONE,R8ZERO, get_args, tIni, tCSV, fill_with_nearest, tUnp, &
    calc_unique, change_case, get_compiler, split_str, get_ext, read_line, parse_line, fileexist, fillgap, &
    get_unique, grid_load_imbalance, readidf, get_dir_files, strip_ext, get_slash, create_dir, &
    quicksort_r, I_I1, I_I2, I_I4, I_I8, I_R4, I_R8, I_C, tGrid, get_neighbors, replace_token, &
    get_elapsed_time
  use vrt_module, only: tVrt, i_SourceFilename
  !
  use hdrModule, only: tHdr, tHdrHdr, writeflt, &
    i_dscl_nointp, i_dscl_intp, i_uscl_arith
  !
  use quad2dModule, only: tQuads, tQuad, tIntf, tNbrIntf, &
    tLayerModels, tLayerModel, tDataModels, tDataModel, tDataModelData, &
    tProps, get_number_of_levels, get_refinement_level, &
    i_lid, i_gid, i_i_graph, i_lay_mod,  &
    i_tgt_cs_min, i_tgt_cs_min_lay_beg, i_tgt_cs_min_lay_end, i_tgt_cs_min_dz_top, &
    i_tgt_cs_max, i_head, i_merge, n_prop_field, &
    tMF6Disu, mf6_data_write, tMF6Exchange, valid_icir
  !
  use multigrid_module, only: tMultiGrid, tMultiGridArray
  use mf6_wbd_mod, only: tMf6Wbd, MAX_NR_CSV
  use metis_module, only: tMetis
  use chaco_module, only: tChaco
  
  implicit none
  !
! -- locals
  integer(I4B), parameter :: i_nodes             = 1
  integer(I4B), parameter :: i_neighbor_lid      = 2
  integer(I4B), parameter :: i_neighbor_nr_cells = 3
  integer(I4B), parameter :: n_part_field = i_neighbor_nr_cells
  character(len=MXSLEN), dimension(n_part_field) :: part_fields
  !
  type(tVrt), pointer :: vrt => null()
  type(tHdr), pointer :: hdrg => null()
  type(tHdr), pointer :: hdrg_gid => null(), hdrg_mask => null()
  type(tQuad), pointer :: q => null()
  type(tQuads), pointer :: xq => null()
  type(tBb), pointer :: bb => null()
  type(tBb), dimension(:), pointer :: bb_gid => null()
  type(tBbObj) :: bbo
  type(tUnp), dimension(:), allocatable :: regbb
  !
  type(tBb)  :: bbi
  type(tBbX) :: bbx
  type(tBbX) :: bbx_t
  !
  character(len=MXSLEN), dimension(:), allocatable :: args, csv_val
  character(len=MXSLEN) :: f_vrt, f_tile, d_out, f_in_csv, f_in_csv_merge, f_out_csv, f_part_csv
  character(len=MXSLEN) :: f_in_csv_pref, f_in_csv_post
  character(len=MXSLEN) :: f_mod_def_inp, f_lay_mod_csv, fp
  character(len=MXSLEN) :: d_in, uuid_in, uuid_out
  character(len=MXSLEN) :: f_gid_in, f_gid_out, f_gid_mask
  character(len=MXSLEN) :: mod_root_dir, xch_root_dir, mod_sub_dir_fields, xch_id_field
  character(len=MXSLEN) :: sel_field, sel_val, sel_npart, sel_nodes
  character(len=MXSLEN) :: chd_lid
  character(len=MXSLEN) :: f_hiera_in, f_hiera_field, f_weight
  character(len=MXSLEN) :: gid_exclude, gid_separate, gid_first_num
  character(len=MXSLEN) :: csv_field, csv_fill_field
  character(len=MXSLEN) :: f_in_idf
  character(len=MXSLEN) :: d_log
  character(len=MXSLEN) :: f_in_flt, f_in_vrt_1, f_in_vrt_2, f_out_vrt, post_fix
  character(len=MXSLEN) :: f_exe
  character(len=MXSLEN) :: elapsed_line
  character(len=MXSLEN) :: vtk_lid, vtk_csv
  !
  ! fields
  character(len=MXSLEN), dimension(n_prop_field) :: fields
  !
  logical :: lrenumber, lwrite, ljoin, lsplit, lremove, lwrite_props
  logical :: loverwrite_props, lwrite_asc, luse_uuid, lwrite_disu
  logical, parameter :: LDUM = .true.
  logical :: write_nod_map, write_chd, write_hiera
  logical :: luse_chaco, lwrite_ximbal
  
  integer(I4B) :: mask_mv, xid_mv
  integer(I4B), dimension(:,:), allocatable :: mask
  !
  integer(I4B) :: nblk, nid, nlid, gid_max, ib, lid, lid_min, lid_max, ljd, gid, nlev, n, m, ntry
  integer(I4B) :: icb, irb, ic, ir, ic0, ic1, ir0, ir1, gnr, gnc, nr, nc, gic, gir, refi4
  integer(I4B) :: icb0, icb1, irb0, irb1, jc, jr, kc, kr, jc0, jc1, jr0, jr1, ntile, it, it0, it1, iu, ilev
  integer(I4B) :: mvxid, gid_mv, area_min, area_max, n_inact, area, i, idum, ireg, jreg, nreg
  integer(I4B) :: run_opt, nzlev_max, gid_mask
  integer(I4B) :: npart, max_iter, max_np, n_csv_val
  integer(I4B), dimension(:), allocatable :: gids, l2gid, g2lid
  integer(I4B), dimension(:,:), allocatable :: xid, i4w2d, regun
  integer(I4B) :: np_beg, np_step, np_end
  integer(I4B), dimension(8) :: ibdt, iedt
  !
  real(R4B), dimension(:,:), allocatable :: r4w2d
  !
  real(R4B) :: r4mv, tgt_imbal
  !
  real(R8B) :: x0, x1, y0, y1
  real(R8B) :: refr8, xll, yll, yur, cs_gid, cs_max
  real(R8B) :: weight_fac
  real(R8B) :: wgt_fac, wgt_per
  real(R8B) :: elsec
  !
  integer(I4B) :: kper_beg, kper_end, tile_nc, tile_nr
  character(len=MXSLEN) :: head_layers
  !
  save
  
contains

subroutine quad_settings()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  !
  ! -- local
  type(tIni) :: ini
  character(len=MXSLEN) :: sect, cdate, cversion, s
! ------------------------------------------------------------------------------
  args = get_args()
  !
  if (size(args) < 2) then
    call errmsg('Invalid number of program arguments.')
  end if
  lid_min = 0
  lid_max = 0 
  if (size(args) > 2) then
    read(args(3),*) lid_min
  end if
  if (size(args) > 3) then
    read(args(4),*) lid_max
  end if
  !
  call logmsg('====================================================================')
  call get_compiler(cdate, cversion)
  call logmsg('Program QUAD2D compiled '//trim(cdate)//' using '//trim(cversion))
  !
  s = args(1); sect = change_case(trim(adjustl(s)),'l')
  !
  ! read the options required for each run option
  call ini%init(args(2))
  !
  ! variable used for ALL options
  call ini%get_val('_general', 'gid_max', i4v=gid_max, i4v_def=5000000)
  !
  run_opt = -1
  select case(sect)
  case('filter')
    !=========!
    run_opt = 0
    !=========!
    call ini%get_val(sect, 'f_gid_in',  cv=f_gid_in)
    call ini%get_val(sect, 'f_gid_out',  cv=f_gid_out)
    call ini%get_val(sect, 'f_gid_mask', cv=f_gid_mask)
    call ini%get_val(sect, 'area_min', i4v=area_min)
    call ini%get_val(sect, 'area_max', i4v=area_max, i4v_def=0) ! not yet used
  case('init')
    !=========!
    run_opt = 1
    !=========!
    call ini%get_val(sect, 'f_vrt',    cv=f_vrt)
    call ini%get_val(sect, 'd_out',    cv=d_out)
    call ini%get_val(sect, 'renumber', l4v=lrenumber, l4v_def=.false.)
    call ini%get_val(sect, 'split',    l4v=lsplit, l4v_def=.false.)
    call ini%get_val(sect, 'join',     l4v=ljoin, l4v_def=.false.)
    call ini%get_val(sect, 'area_min', i4v=area_min, i4v_def=0)
    call ini%get_val(sect, 'area_max', i4v=area_max, i4v_def=huge(I4ZERO))
    call ini%get_val(sect, 'cs_max',   r8v=cs_max)
    call ini%get_val(sect, 'gid_mv',   i4v=mvxid, i4v_def=0)
    call ini%get_val(sect, 'use_uuid', l4v=luse_uuid, l4v_def=.false.)
    call ini%get_val(sect, 'uuid_out', cv=uuid_out, cv_def='quad2d')
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv, cv_def='')
  case('balance')
    !=========!
    run_opt = 2
    !=========!
    call ini%get_val(sect, 'write_props', l4v=lwrite_props, l4v_def=.true.)
    !
    call ini%get_val(sect, 'd_in',    cv=d_in)
    call ini%get_val(sect, 'uuid_in', cv=uuid_in, cv_def='quad2d')
    !
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv)
  case('interface_gen')
    !=========!
    run_opt = 3
    !=========!
    call ini%get_val(sect, 'write_props', l4v=lwrite_props, l4v_def=.false.)
    !
    call ini%get_val(sect, 'd_in',    cv=d_in)
    call ini%get_val(sect, 'uuid_in', cv=uuid_in, cv_def='quad2d')
    !
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv, cv_def='')
    !
    call ini%get_val(sect, 'f_mod_def_inp', cv=f_mod_def_inp)
    !
    call ini%get_val(sect, 'lid_field',       cv=fields(i_lid),     cv_def='lid')
    call ini%get_val(sect, 'gid_field',       cv=fields(i_gid),     cv_def='gid')
    call ini%get_val(sect, 'i_i_graph_field', cv=fields(i_i_graph), cv_def='i_i_graph')
    !
    call ini%get_val(sect, 'lay_mod_field',            cv=fields(i_lay_mod))
  case('grid_gen','grid_gen_merge')
    !=========!
    if (sect == 'grid_gen') then
      run_opt = 4
    else
      run_opt = 17
    end if
    !=========!
    call ini%get_val(sect, 'write_props', l4v=lwrite_props, l4v_def=.true.)
    call ini%get_val(sect, 'overwrite_props', l4v=loverwrite_props, l4v_def=.false.)
    call ini%get_val(sect, 'write_asc', l4v=lwrite_asc, l4v_def=.false.)
    !
    call ini%get_val(sect, 'd_in',    cv=d_in)
    call ini%get_val(sect, 'uuid_in', cv=uuid_in, cv_def='quad2d')
    !
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv, cv_def='')
    if (size(args) > 4) then
      f_out_csv = trim(strip_ext(f_out_csv))//trim(args(5))//'.csv'
    end if
    !
    call ini%get_val(sect, 'f_mod_def_inp', cv=f_mod_def_inp)
    !
    call ini%get_val(sect, 'lid_field',       cv=fields(i_lid),     cv_def='lid')
    call ini%get_val(sect, 'gid_field',       cv=fields(i_gid),     cv_def='gid')
    call ini%get_val(sect, 'i_i_graph_field', cv=fields(i_i_graph), cv_def='i_i_graph')
    !
    call ini%get_val(sect, 'lay_mod_field',            cv=fields(i_lay_mod))
    call ini%get_val(sect, 'tgt_cs_min_field',         cv=fields(i_tgt_cs_min))
    call ini%get_val(sect, 'tgt_cs_max_field',         cv=fields(i_tgt_cs_max))
    call ini%get_val(sect, 'tgt_cs_min_lay_beg_field', cv=fields(i_tgt_cs_min_lay_beg), cv_def='')
    call ini%get_val(sect, 'tgt_cs_min_lay_end_field', cv=fields(i_tgt_cs_min_lay_end), cv_def='')
    call ini%get_val(sect, 'tgt_cs_min_dz_top_field',  cv=fields(i_tgt_cs_min_dz_top), cv_def='')
    !
    call ini%get_val(sect, 'mod_root_dir', cv=mod_root_dir)
    call ini%get_val(sect, 'mod_sub_dir_fields', cv=mod_sub_dir_fields, cv_def='gid')
    !
    call ini%get_val(sect, 'write_disu', l4v=lwrite_disu, l4v_def=.true.)
    !
    if (run_opt == 17) then
      call ini%get_val(sect, 'i_merge', cv=fields(i_merge))
    end if
    !
  case('mf6_xch_write','mf6_xch_write_merge')
    if (sect == 'mf6_xch_write') then
      run_opt = 5
    else
      run_opt = 19
    end if
    !=========!
    call ini%get_val(sect, 'd_in',    cv=d_in)
    call ini%get_val(sect, 'uuid_in', cv=uuid_in, cv_def='quad2d')
    !
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    if (run_opt == 19) then
      call ini%get_val(sect, 'f_in_csv_merge', cv=f_in_csv_merge)
    end if
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv)
    !
    call ini%get_val(sect, 'xch_root_dir', cv=xch_root_dir)
    call ini%get_val(sect, 'xch_id_field', cv=xch_id_field, cv_def='lid')
    
    call ini%get_val(sect, 'lid_field',       cv=fields(i_lid),     cv_def='lid')
    call ini%get_val(sect, 'gid_field',       cv=fields(i_gid),     cv_def='gid')
    call ini%get_val(sect, 'i_i_graph_field', cv=fields(i_i_graph), cv_def='i_i_graph')
    !
    call ini%get_val(sect, 'lay_mod_field',            cv=fields(i_lay_mod))
  case('mf6_data_write', 'mf6_data_write_merge')
    if (sect == 'mf6_data_write') then
      run_opt = 6
    else
      run_opt = 18
    end if
    call ini%get_val(sect, 'd_in',    cv=d_in)
    call ini%get_val(sect, 'd_log',   cv=d_log, cv_def='')
    call ini%get_val(sect, 'uuid_in', cv=uuid_in, cv_def='quad2d')
    !
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    if (run_opt == 18) then
      call ini%get_val(sect, 'f_in_csv_merge', cv=f_in_csv_merge)
    end if
    !
    call ini%get_val(sect, 'f_mod_def_inp', cv=f_mod_def_inp)
    !
    call ini%get_val(sect, 'lid_field',        cv=fields(i_lid),        cv_def='lid')
    call ini%get_val(sect, 'gid_field',        cv=fields(i_gid),        cv_def='gid')
    call ini%get_val(sect, 'tgt_cs_min_field', cv=fields(i_tgt_cs_min), cv_def='tgt_cs_min')
    call ini%get_val(sect, 'lay_mod_field',    cv=fields(i_lay_mod))
    !
    call ini%get_val(sect, 'mod_root_dir', cv=mod_root_dir)
    call ini%get_val(sect, 'mod_sub_dir_fields', cv=mod_sub_dir_fields, cv_def='')
  case('partition')
    !=========!
    run_opt = 7
    !=========!
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv)
    call ini%get_val(sect, 'f_part_csv', cv=f_part_csv, cv_def='')
    !
    call ini%get_val(sect, 'lid_field',  cv=fields(i_lid), cv_def='lid')
    !
    call ini%get_val(sect, 'nodes_field',              cv=part_fields(i_nodes), &
      cv_def='nodes')
    call ini%get_val(sect, 'neighbor_lid_field',       cv=part_fields(i_neighbor_lid), &
      cv_def='neighbor_lid')
    call ini%get_val(sect, 'neighbor_nr_cells_field',  cv=part_fields(i_neighbor_nr_cells), &
      cv_def='neighbor_nr_cells')
    
    call ini%get_val(sect, 'sel_field', cv=sel_field)
    call ini%get_val(sect, 'sel_val',   cv=sel_val)
    call ini%get_val(sect, 'sel_npart', cv=sel_npart, cv_def='')
    call ini%get_val(sect, 'sel_nodes', cv=sel_nodes, cv_def='')
  case('mf6_post')
    !=========!
    run_opt = 8
    !=========!
    call ini%get_val(sect, 'd_in',    cv=d_in)
    call ini%get_val(sect, 'uuid_in', cv=uuid_in, cv_def='quad2d')
    !
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    call ini%get_val(sect, 'f_in_csv_merge', cv=f_in_csv_merge, cv_def='')
    !
    call ini%get_val(sect, 'f_mod_def_inp', cv=f_mod_def_inp)
    call ini%get_val(sect, 'lay_mod_field',            cv=fields(i_lay_mod))
    !
    call ini%get_val(sect, 'lid_field',        cv=fields(i_lid),        cv_def='lid')
    call ini%get_val(sect, 'gid_field',        cv=fields(i_gid),        cv_def='gid')
    call ini%get_val(sect, 'tgt_cs_min_field', cv=fields(i_tgt_cs_min), cv_def='tgt_cs_min')
    call ini%get_val(sect, 'head_field',       cv=fields(i_head))
    call ini%get_val(sect, 'head_layers',      cv=head_layers)
    !
    call ini%get_val(sect, 'kper_beg', i4v=kper_beg, i4v_def=1)
    call ini%get_val(sect, 'kper_end', i4v=kper_end, i4v_def=1)
    call ini%get_val(sect, 'tile_nc', i4v=tile_nc, i4v_def=huge(I4ZERO))
    call ini%get_val(sect, 'tile_nr', i4v=tile_nr, i4v_def=huge(I4ZERO))
    call ini%get_val(sect, 'f_out_vrt_pref', cv=f_vrt)
    call ini%get_val(sect, 'write_nod_map', l4v=write_nod_map, l4v_def=.false.)
    !
    call ini%get_val(sect, 'write_chd', l4v=write_chd, l4v_def=.false.)
    call ini%get_val(sect, 'chd_lid',cv=chd_lid, cv_def='')
    !
    if (write_chd) then
      call ini%get_val(sect, 'mod_root_dir', cv=mod_root_dir)
      call ini%get_val(sect, 'mod_sub_dir_fields', cv=mod_sub_dir_fields)
    end if
    !
    call ini%get_val(sect, 'vtk_lid',cv=vtk_lid, cv_def='')
    call ini%get_val(sect, 'vtk_csv',cv=vtk_csv, cv_def='')
  case('mf6_post_wtd')
    !=========!
    run_opt = 20
    !=========!
    call ini%get_val(sect, 'f_in_flt_ahn',   cv=f_in_flt)
    call ini%get_val(sect, 'f_in_vrt_head',  cv=f_in_vrt_1)
    call ini%get_val(sect, 'f_in_vrt_nod',   cv=f_in_vrt_2)
    call ini%get_val(sect, 'f_out_vrt_wtd',  cv=f_out_vrt)
    call ini%get_val(sect, 'f_out_flt_post', cv=post_fix)
  case('fill_gap')
    !=========!
    run_opt = 9
    !=========!
    call ini%get_val(sect, 'f_gid_mask', cv=f_gid_mask)
    call ini%get_val(sect, 'gid_mask', i4v=gid_mask)
    call ini%get_val(sect, 'f_gid_in',  cv=f_gid_in)
    call ini%get_val(sect, 'f_gid_out',  cv=f_gid_out)
  case('unique')
    !=========!
    run_opt = 10
    !=========!
    call ini%get_val(sect, 'f_gid_in',  cv=f_gid_in)
    call ini%get_val(sect, 'f_gid_out',  cv=f_gid_out)
  case('sub_grid_part')
    !=========!
    run_opt = 11
    !=========!
    call ini%get_val(sect, 'f_gid_in',      cv =f_gid_in)
    call ini%get_val(sect, 'f_weight',      cv =f_weight,       cv_def='')
    call ini%get_val(sect, 'f_gid_out',     cv =f_gid_out)
    call ini%get_val(sect, 'gid_separate',  cv =gid_separate,  cv_def='')
    call ini%get_val(sect, 'max_np',        i4v =max_np)
    call ini%get_val(sect, 'chaco_exe',     cv  =f_exe, cv_def='')
    call ini%get_val(sect, 'write_ximbal',  l4v=lwrite_ximbal, l4v_def=.false.)
    call ini%get_val(sect, 'use_chaco',     l4v=luse_chaco, l4v_def=.false.)
    call ini%get_val(sect, 'wgt_fac',       r8v=wgt_fac, r8v_def=0.7d0)
    call ini%get_val(sect, 'wgt_per',       r8v=wgt_per, r8v_def=0.05d0)
   
  case('check_sub_grid_part')
    !=========!
    run_opt = 21
    !=========!
    call ini%get_val(sect, 'f_gid_in',      cv =f_gid_in)
    call ini%get_val(sect, 'f_weight',      cv =f_weight,       cv_def='')
    call ini%get_val(sect, 'np_beg',        i4v =np_beg)
    call ini%get_val(sect, 'np_step',       i4v =np_step, i4v_def=1)
    call ini%get_val(sect, 'np_end',        i4v =np_end)
    call ini%get_val(sect, 'chaco_exe',     cv  =f_exe, cv_def='')
    call ini%get_val(sect, 'use_chaco',     l4v=luse_chaco, l4v_def=.false.)
    
  !case('sub_grid_part')
  !  !=========!
  !  run_opt = 11
  !  !=========!
  !  call ini%get_val(sect, 'f_gid_in',      cv=f_gid_in)
  !  call ini%get_val(sect, 'f_weight',      cv=f_weight, cv_def='')
  !  call ini%get_val(sect, 'f_hiera_in',    cv=f_hiera_in, cv_def='')
  !  call ini%get_val(sect, 'npart',        i4v=npart)
  !  call ini%get_val(sect, 'f_gid_out',     cv=f_gid_out)
  !  call ini%get_val(sect, 'f_hiera_field', cv=f_hiera_field, cv_def='')
  !  call ini%get_val(sect, 'f_out_csv',     cv=f_out_csv)
  !  call ini%get_val(sect, 'write_hiera',   l4v=write_hiera, l4v_def=.false.)
  !  call ini%get_val(sect, 'gid_separate',  cv=gid_separate, cv_def='')
  !  call ini%get_val(sect, 'weight_fac',    r8v=weight_fac, r8v_def=R8ONE)
  case('full_grid_part')
    !=========!
    run_opt = 12
    !=========!
    call ini%get_val(sect, 'f_gid_mask',  cv=f_gid_mask)
    call ini%get_val(sect, 'f_weight',    cv=f_weight, cv_def='')
    call ini%get_val(sect, 'npart',       i4v=npart)
    call ini%get_val(sect, 'f_gid_out',   cv=f_gid_out, cv_def='')
    call ini%get_val(sect, 'gid_exclude', cv=gid_exclude, cv_def='')
    call ini%get_val(sect, 'gid_separate', cv=gid_separate, cv_def='')
  case('grid_balancing')
    !=========!
    run_opt = 13
    !=========!
    call ini%get_val(sect, 'f_gid_in',      cv =f_gid_in)
    call ini%get_val(sect, 'f_weight',      cv =f_weight,       cv_def='')
    call ini%get_val(sect, 'tgt_imbal',     r4v=tgt_imbal,     r4v_def=1.1)
    call ini%get_val(sect, 'f_gid_out',     cv =f_gid_out)
    call ini%get_val(sect, 'gid_first_num', cv =gid_first_num,  cv_def='')
    call ini%get_val(sect, 'max_iter',      i4v =max_iter,     i4v_def=100)
    call ini%get_val(sect, 'max_np',        i4v =max_np,       i4v_def=10000)
    call ini%get_val(sect, 'chaco_exe',     cv  =f_exe, cv_def='')
  case('idf_to_flt')
    !=========!
    run_opt = 15
    !=========!
    call ini%get_val(sect, 'f_in_idf', cv=f_in_idf)
    call ini%get_val(sect, 'mv', r4v=r4mv)
  end select
  !
  if (sect(1:9) == 'merge_csv') then
    !=========!
    run_opt = 16
    !=========!
    call ini%get_val(sect, 'f_in_csv_pref', cv=f_in_csv_pref)
    call ini%get_val(sect, 'f_in_csv_post', cv=f_in_csv_post)
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv)
  end if
  if (sect(1:13) == 'csv_add_field') then
    !=========!
    run_opt = 14
    !=========!
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv)
    call ini%get_val(sect, 'csv_field', cv=csv_field)
    call ini%get_val(sect, 'csv_fill_field', cv=csv_fill_field, cv_def='')
    call ini%get_val(sect, 'n_csv_val', i4v =n_csv_val, i4v_def=1)
    allocate(csv_val(n_csv_val)); csv_val = ''
    if (n_csv_val > 1) then
      do i = 1, n_csv_val
        call ini%get_val(sect, 'csv_val_'//trim(ta([i])), cv=csv_val(i))
      end do
    else
      call ini%get_val(sect, 'csv_val', cv=csv_val(1))
    end if
  end if
  !
  if (run_opt == -1) then
    call errmsg('Invalid run option: '//trim(sect))
  end if
  !
  call logmsg('***** Run option: '//trim(sect)//' *****')
  call logmsg('====================================================================')
  !
  ! clean up
  call ini%clean()
  !
  return
end subroutine quad_settings
!
subroutine quad_filter()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  allocate(hdrg_gid, hdrg_mask)
  !
  call hdrg_mask%read_full_grid(f_gid_mask)
  call hdrg_mask%get_grid(xi4=mask, mvi4=mask_mv)
  call hdrg_mask%clean(); deallocate(hdrg_mask); hdrg_mask => null()
  !
  call hdrg_gid%read_full_grid(f_gid_in)
  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
  !
  nc = size(xid,1); nr = size(xid,2)
  if ((nc /= size(mask,1)).or.(nr /= size(mask,2))) then
    call errmsg('quad_filter: invalid dimensions')
  end if
  !
  allocate(r4w2d(nc,nr))
  !
  ! apply the mask

  gid_max = maxval(xid)
  !
  allocate(bb_gid(gid_max), gids(gid_max))
  gids = 0
  !
  do ir = 1, nr
    do ic = 1, nc
      gid = xid(ic,ir)
      if (gid /= xid_mv) then
        gids(gid) = 1
        bb => bb_gid(gid)
        bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
        bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
        bb%nrow = bb%ir1 - bb%ir0 + 1
        bb%ncol = bb%ic1 - bb%ic0 + 1
      end if
    end do
  end do
  !
  do ir = 1, nr
    do ic = 1, nc
      r4w2d(ic,ir) = real(xid(ic,ir),R4B)
      if (mask(ic,ir) /= mask_mv) then
        if (xid(ic,ir) == xid_mv) then
          r4w2d(ic,ir) = -1.
        end if
      end if
    end do
  end do
  !
  do gid = 1, gid_max
    if (gids(gid) == 1) then
      bb => bb_gid(gid)
      !
      if(allocated(i4w2d)) deallocate(i4w2d)
      allocate(i4w2d(bb%ncol, bb%nrow))
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
        if (xid(ic,ir) == gid) then
          i4w2d(jc,jr) = 1
        else
          i4w2d(jc,jr) = 0
        end if
      end do; end do
      !
      call calc_unique(i4w2d, 5, regun, regbb, nreg, idum, 0., 0., 0.)
      !
      ! determine the largest region
      n = 0
      do ireg = 1, nreg
        if (regbb(ireg)%n > n) then
          jreg = ireg; n = regbb(ireg)%n
        end if
      end do
      !
      do ireg = 1, nreg
        lremove = .false.
        if (ireg /= jreg) then
          lremove = .true.
        else
          ! check area
          if (regbb(ireg)%n < area_min) then
            lremove = .true.
          end if
          !
          if (.not. lremove) then
            n = 0
            do ir = 1, bb%nrow
              m = 0
              do ic = 1, bb%ncol
                if (regun(ic,ir) == ireg) then
                  m = m + 1
                end if
              end do
              n = max(n,m)
            end do
            if (n <= 5) then
              lremove = .true.
            else
              n = 0
              do ic = 1, bb%ncol
                m = 0
                do ir = 1, bb%nrow
                  if (regun(ic,ir) == ireg) then
                    m = m + 1
                  end if
                end do
                n = max(n,m)
              end do
            end if
            if (n <= 5) then
              lremove = .true.
            end if
          end if
          !
        end if
        if (lremove) then
          do ir = 1, bb%nrow; do ic = 1, bb%ncol
            if (regun(ic,ir) == ireg) then
              jr = bb%ir0 + ir - 1; jc = bb%ic0 + ic - 1
              r4w2d(jc,jr) = -1.
            end if
          end do; end do
        end if
      end do
      !
    end if
  end do
  !
  call fill_with_nearest(r4w2d, real(xid_mv,R4B), -1.)
  !
  xid = int(r4w2d, I4B)
  call hdrg_gid%set_grid(xid, xid_mv)
  call hdrg_gid%write(f_gid_out)
  !
  call hdrg_gid%clean(); deallocate(hdrg_gid); hdrg_gid => null()
  deallocate(r4w2d, bb_gid, gids, i4w2d)
  !
  return
end subroutine quad_filter

subroutine quad_unique()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- local
  integer(I4B) :: n, i
  integer(I4B), dimension(:,:), allocatable :: xid_new
! ------------------------------------------------------------------------------
  !
  allocate(hdrg_gid)
  call hdrg_gid%read_full_grid(f_gid_in)
  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
  !
  nc = size(xid,1); nr = size(xid,2)
  !
  if (allocated(i4w2d)) deallocate(i4w2d)
  allocate(xid_new(nc,nr)); xid_new = xid_mv
  !
  gid_max = maxval(xid)
  !
  allocate(bb_gid(gid_max), gids(gid_max))
  gids = 0
  !
  do ir = 1, nr
    do ic = 1, nc
      gid = xid(ic,ir)
      if (gid /= xid_mv) then
        gids(gid) = 1
        bb => bb_gid(gid)
        bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
        bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
        bb%nrow = bb%ir1 - bb%ir0 + 1
        bb%ncol = bb%ic1 - bb%ic0 + 1
      end if
    end do
  end do
  !
  n = 0
  do gid = 1, gid_max
    if (gids(gid) == 1) then
      bb => bb_gid(gid)
      !
      if(allocated(i4w2d)) deallocate(i4w2d)
      allocate(i4w2d(bb%ncol, bb%nrow))
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
        if (xid(ic,ir) == gid) then
          i4w2d(jc,jr) = 1
        else
          i4w2d(jc,jr) = 0
        end if
      end do; end do
      !
      call calc_unique(i4w2d, 5, regun, regbb, nreg, idum, 0., 0., 0.)
      !
      do ireg = 1, nreg
        n = n + 1
        do ir = 1, bb%nrow; do ic = 1, bb%ncol
          if (regun(ic,ir) == ireg) then
            jr = bb%ir0 + ir - 1; jc = bb%ic0 + ic - 1
            if (xid(jc,jr) == gid) then
              xid_new(jc,jr) = n
            end if
          end if
        end do; end do
      end do
      !
    end if
  end do
  !
  call hdrg_gid%replace_grid(xi4=xid_new, mvi4=xid_mv)
  call hdrg_gid%write(f_gid_out)
  !
  ! clean-up
  if (allocated(xid)) deallocate(xid)
  if (allocated(xid_new)) deallocate(xid_new)
  if (allocated(i4w2d)) deallocate(i4w2d)
  if (associated(hdrg_gid)) then
    call hdrg_gid%clean(); deallocate(hdrg_gid); hdrg_gid => null()
  end if
  !
  return
end subroutine quad_unique

subroutine quad_full_grid_part()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- local
  type(tMetis), pointer :: met => null()
  character(len=MXSLEN) :: f_out
  real(R4B) :: wtot, wloc
  real(R8B) :: imbal
  integer(I4B), dimension(:), allocatable :: ids
  integer(I4B), dimension(:,:), allocatable :: weight, weight_loc
  integer(I4B) :: weight_mv, id, npart_loc, npart_tot, ip, ip_offset, nid
  integer(I4B) :: n_disjoint
! ------------------------------------------------------------------------------
  !
  if (len_trim(f_gid_out) > 0) then
    f_out= f_gid_out
  else
    f_out = trim(strip_ext(f_gid_mask))//'_fp_np'//ta([npart],'(i3.3)')
  end if
  !
  allocate(hdrg_gid)
  call hdrg_gid%read_full_grid(f_gid_mask)
  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
  nc =  size(xid,1); nr = size(xid,2)
  !
  ! replace missing values
  if (xid_mv /= 0) then
    do ir = 1, nr; do ic = 1, nc
      if (xid(ic,ir) == xid_mv) then
        xid(ic,ir) = 0
     end if
    end do; end do
    xid_mv = 0
  end if
  !
  if (len_trim(gid_exclude) > 0) then
    call parse_line(s=gid_exclude, i4a=ids, token_in=',')
    do i = 1, size(ids)
      do ir = 1, nr; do ic = 1, nc
        if (ids(i) == xid(ic,ir)) then
          xid(ic,ir) = xid_mv
        end if
      end do; end do
    end do
  end if
  !
  ! read the weights
  if (len_trim(f_weight) > 0) then
    allocate(hdrg)
    call hdrg%read_full_grid(f_weight)
    call hdrg%get_grid(xi4=weight, mvi4=weight_mv)
    call hdrg%clean(); deallocate(hdrg); hdrg => null()
  else
    allocate(weight(nc,nr)); weight = 1
  end if
  !
  ! filter for xid
  do ir = 1, nr; do ic = 1, nc
   if (xid(ic,ir) > 0) then
     if (weight(ic,ir) == 0) then
       call logmsg('**** weight zero found! ****')
     end if
     weight(ic,ir) = max(1,weight(ic,ir))
   else
     weight(ic,ir) = 0
    end if
  end do; end do
  !
  if (len_trim(gid_separate) > 0) then
    wtot = real(sum(weight),R4B)
    npart_tot = npart
    call parse_line(s=gid_separate, i4a=ids, token_in=',')
    ip_offset = 0
    nid = size(ids)
    do i = 1, nid + 1
      if (allocated(weight_loc)) deallocate(weight_loc)
      allocate(weight_loc(nc,nr)); weight_loc = 0
      wloc = 0
      if (i <= nid) then
        id = ids(i)
        call logmsg('Applying METIS to ID='//ta([id])//'...')
        do ir = 1, nr; do ic = 1, nc
          if (xid(ic,ir) /= xid_mv) then
            if (xid(ic,ir) == id) then
              wloc = wloc + real(weight(ic,ir),R4B)
              weight_loc(ic,ir) = weight(ic,ir)
            end if
          end if
        end do; end do
      else ! remainder
        do ir = 1, nr; do ic = 1, nc
          if (xid(ic,ir) /= xid_mv) then
            if (xid(ic,ir) > 0) then
              wloc = wloc + real(weight(ic,ir),R4B)
              weight_loc(ic,ir) = weight(ic,ir)
            end if
          end if
        end do; end do
      end if
      !
      npart_loc = max(1,nint(npart*wloc/wtot))
      npart_loc = min(npart_loc, npart_tot)
      if (npart_loc > 1) then
        allocate(met)
        call met%init(weight_loc, npart_loc)
        call met%set_opts(niter_in=100)
        call met%recur()
        call met%set_ids(weight_loc, id_offset_in=ip_offset)
        call met%clean(); deallocate(met); met => null()
        do ir = 1, nr; do ic = 1, nc
          ip = weight_loc(ic,ir)
          if (ip > 0) then
            xid(ic,ir) = -ip
          end if
        end do; end do
        deallocate(weight_loc)
      else
        if (i <= nid) then
          call logmsg('***** No METIS partitioning for ID='//ta([id])//' (1 parts) *****')
        else
          call logmsg('***** No METIS partitioning for remainder (1 parts) *****')
        end if
        do ir = 1, nr; do ic = 1, nc
          ip = weight_loc(ic,ir)
          if (ip > 0) then
            xid(ic,ir) = -(ip_offset+1)
          end if
        end do; end do
      end if
      ip_offset = ip_offset + npart_loc
      npart_tot = max(0, npart_tot - npart_loc)
    end do
    !
    xid = abs(xid)
    n_disjoint = quad_count_disjoint(xid, xid_mv)
    call logmsg('# disjoint: '//ta([n_disjoint]))
    call grid_load_imbalance(xid, xid_mv, weight, imbal, npart_loc)
    call logmsg('Overall load imbalance for '//ta([npart_loc])//' parts: '//ta([imbal]))
    !f_out = trim(f_out)//'_kib'//adjustl(ta([1000.d0*imbal],'(f10.0)'))
  else
    allocate(met)
    call met%init(weight, npart)
    call met%set_opts(niter_in=100)
    call met%recur()
    call met%set_ids(xid)
    call met%clean(); deallocate(met); met => null()
  end if
  !
  ! write the new ids
  call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
  call hdrg_gid%write(f_out)
  
  ! clean-up
  if (allocated(xid)) deallocate(xid)
  if (associated(hdrg_gid)) then
    call hdrg_gid%clean(); deallocate(hdrg_gid); hdrg_gid => null()
  end if
  !
  return
end subroutine quad_full_grid_part

!subroutine quad_grid_balancing()
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------  
!! -- local
!  !
!  logical, parameter :: lwriteximbal = .true.
!  !
!  type(tQuad), pointer :: q_new => null()
!  type(tBb)  :: bbi_gid, bb_new
!  type(tBbX) :: bbx_gid, bbx_hid, bbx_q
!  type(tHdr), pointer :: hdrg_hid => null()
!  type(tHdrHdr), pointer :: hdr => null()
!  type(tCSV), pointer    :: csv => null()
!  type(tMetis), pointer :: met => null()
!  !
!  logical :: lgidfirst, ldone, lmetis, lfirst
!  character(len=MXSLEN) :: f
!  integer(I4B) :: weight_mv, gid_max_loc, iter, j
!  integer(I4B) :: lid_new, gid_new, iact, nact, np_eval, d_np, nmet, np_met
!  integer(I4B) :: np_strt, np_tgt, np, np_full, np_diff, ip, p
!  integer(I4B), dimension(:,:), allocatable :: mask, part, xid_org
!  integer(I4B), dimension(:,:), allocatable :: weight, qweight
!  integer(I4B), dimension(:), allocatable :: gid_first
!  integer(I4B), dimension(:), allocatable :: gids_new
!  real(R4B), dimension(:), allocatable :: wgt_sort
!  real(R4B), dimension(:,:), allocatable :: ximbal
!  real(R8B) :: wgt_avg, wgt_tgt, wgt_tot, wgt_q, imbal
!! ------------------------------------------------------------------------------
!  !
!  allocate(hdrg_gid)
!  call hdrg_gid%read_full_grid(f_gid_in)
!  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
!  allocate(xid_org, source=xid)
!  nc =  size(xid,1); nr = size(xid,2)
!  !
!  ! read the weights
!  if (len_trim(f_weight) > 0) then
!    allocate(hdrg)
!    call hdrg%read_full_grid(f_weight)
!    call hdrg%get_grid(xi4=weight, mvi4=weight_mv)
!    call hdrg%clean(); deallocate(hdrg); hdrg => null()
!  else
!    allocate(weight(nc,nr)); weight = 1
!  end if
!  !
!  if (len_trim(gid_first_num) > 0) then
!    lgidfirst = .true.
!    call parse_line(s=gid_first_num, i4a=gid_first, token_in=',')
!  else
!    lgidfirst = .false.
!  end if
!  !
!  ! set bbi and bbx for gid
!  bbi_gid%ic0 = 1;   bbi_gid%ic1 = nc
!  bbi_gid%ir0 = 1;   bbi_gid%ir1 = nr
!  bbi_gid%ncol = nc; bbi_gid%nrow = nr
!  hdr => hdrg_gid%hdr; call hdr%get_bbx(bbx_gid)
!  xll = bbx_gid%xll; yll = bbx_gid%yll; cs_gid = bbx_gid%cs
!  call bbo%set(prent_bbi=bbi_gid, prent_bbx=bbx_gid)
!  !
!  ! determine the global ids and bounding box
!  allocate(bb_gid(gid_max), gids(gid_max), g2lid(gid_max))
!  gids = 0; g2lid = 0
!  do gid = 1, gid_max
!    bb => bb_gid(gid)
!    call bb%init()
!  end do
!  do ir = 1, nr; do ic = 1, nc
!    gid = xid(ic,ir)
!    if (gid /= xid_mv) then
!      gids(gid) = 1
!      bb => bb_gid(gid)
!      bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
!      bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
!      bb%ncol = bb%ic1 - bb%ic0 + 1; bb%nrow = bb%ir1 - bb%ir0 + 1
!    end if
!  end do; enddo
!  !
!  gid_max_loc = maxval(xid)
!  nlid = sum(gids)
!  if (nlid > gid_max) then
!    call errmsg('Increase gid_max.')
!  end if
!  !
!  d_np = 1000; lfirst = .true.; iter = 0
!  do iter = 1, max_iter
!    !
!    gid_new = gid_max_loc
!    if (allocated(xid)) deallocate(xid)
!    allocate(xid, source=xid_org)
!    allocate(xq)
!    call xq%init(nlid, gid_max)
!    !
!    lid = 0; wgt_tot = R8ZERO
!    do gid = 1, gid_max
!      if (gids(gid) == 1) then
!        lid = lid + 1
!        g2lid(gid) = lid
!        q => xq%get_quad(lid)
!        bb => bb_gid(gid) ! local bounding box
!        bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb%ir1)*cs_gid
!        bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
!        bbx%cs = cs_gid
!        call bbo%set(child_bbi=bb, child_bbx=bbx)
!        call q%init(gid=gid, lid=lid, bbo=bbo)
!        call get_mask(xid, gid, bb, mask)
!        if (allocated(qweight)) deallocate(qweight)
!        allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
!        do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!          jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!          if (mask(jc,jr) > 0) then
!            qweight(jc,jr) = weight(ic,ir)
!          end if
!        end do; end do;
!        call q%calc_prop(mask=mask, weight=qweight)
!        call q%get_prop(weight=wgt_q)
!        wgt_tot = wgt_tot + wgt_q
!      end if
!    end do
!    !
!    np_eval = xq%get_number_active()
!    np_eval = np_eval + (iter-1)*d_np
!    wgt_avg = wgt_tot/np_eval; wgt_tgt = tgt_imbal * wgt_avg
!    !
!    ldone = .true.
!    n = xq%n; lid_new = n; nmet = 0; np_met = 0
!    do lid = 1, n
!      q => xq%get_quad(lid)
!      if (q%get_flag(active=LDUM)) then
!        call q%get_prop(weight=wgt_q)
!        !
!        lmetis = .false.
!        if (wgt_q > wgt_tgt) then
!          np_full = nint(wgt_q/wgt_tgt)
!          np_full = max(2,np_full)
!          if (np_full > 1) then
!            lmetis = .true.
!          end if
!        end if
!        if (lmetis) then
!          ldone = .false.
!          ! number of parts
!          nmet = nmet + 1
!          np_met = np_met + np_full - 1
!          !
!          call q%get_bb(child_bbi=bb)
!          call get_mask(xid, q%gid, bb, mask)
!          if (allocated(qweight)) deallocate(qweight)
!          allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
!          do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!            jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!            if (mask(jc,jr) > 0) then
!              qweight(jc,jr) = weight(ic,ir)
!            end if
!          end do; end do;
!          !
!          ! full grid METIS
!          allocate(met)
!          call met%init(qweight, np_full)
!          call met%set_opts()
!          call met%recur(verbose=.true.)
!          if (allocated(part)) deallocate(part)
!          allocate(part,source=mask)
!          call met%set_ids(part)
!          call met%clean(); deallocate(met); met => null()
!          !
!          do ip = 1, np_full
!            lid_new = lid_new + 1; gid_new = gid_new + 1
!            q_new => xq%get_quad(lid_new)
!            call bb_new%init()
!            do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!              jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!              p = part(jc,jr)
!              if (p == ip) then
!               ! determine the new bounding box
!                bb_new%ic0 = min(bb_new%ic0,ic); bb_new%ic1 = max(bb_new%ic1,ic)
!                bb_new%ir0 = min(bb_new%ir0,ir); bb_new%ir1 = max(bb_new%ir1,ir)
!                bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1;  bb_new%nrow = bb_new%ir1 - bb%ir0 + 1
!                if (xid(ic,ir) /= q%gid) then
!                  call errmsg('quad_sub_grid_part: program error.')
!                else
!                  xid(ic,ir) = gid_new
!                end if
!              end if 
!            end do; end do
!            !
!            bbx%xll = xll + (bb_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb_new%ir1)*cs_gid
!            bbx%xur = bbx%xll + bb_new%ncol*cs_gid; bbx%yur = bbx%yll + bb_new%nrow*cs_gid
!            bbx%cs = cs_gid
!            call bbo%set(child_bbi=bb_new, child_bbx=bbx)
!            call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
!            if (q%gid_prent > 0) then
!              call q_new%set(hlev=q%hlev, gid_prent=q%gid_prent, lid_prent=q%lid_prent)
!            else
!              call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid)
!            end if
!            call get_mask(xid, gid_new, bb_new, mask)
!            if (allocated(qweight)) deallocate(qweight)
!            allocate(qweight(bb_new%ncol,bb_new%nrow)); qweight = I4ZERO
!            do ir = bb_new%ir0, bb_new%ir1; do ic = bb_new%ic0, bb_new%ic1
!              jr = ir - bb_new%ir0 + 1; jc = ic - bb_new%ic0 + 1
!              if (mask(jc,jr) > 0) then
!                qweight(jc,jr) = weight(ic,ir)
!              end if
!            end do; end do;
!            call q_new%calc_prop(mask=mask, weight=qweight)
!            xq%n = lid_new
!          end do !ip
!          !
!          call q%set_flag(active=.false.)
!        end if
!      end if
!    end do
!    ! 
!    call logmsg('***** iteration '//ta([iter])//' *****')
!    call logmsg('Number of times METIS was called: '//ta([nmet]))
!    call logmsg('Number of new METIS parts:        '//ta([np_met]))
!    call grid_load_imbalance(xid, xid_mv, weight, imbal, np)
!    call logmsg('Overall load imbalance for '//ta([np])//' parts: '//ta([imbal]))
!    if (ldone) then
!      call logmsg('Nothing to be done...')
!      exit
!    end if
!    if (imbal <= tgt_imbal) then
!      call logmsg('Target load imbalance reached...')
!      exit
!    end if
!    if (iter >= max_iter) then
!      call logmsg('No convergence, maximum iterations reached: '//ta([max_iter]))
!      exit
!    end if
!    if (np >= max_np) then
!      call logmsg('Maximum number of partitions reached: '//ta([max_np]))
!      exit
!    end if
!    !
!    ! clean up
!    call xq%clean(); deallocate(xq)
!  end do
!  
!  ! renumber
!  allocate(gids_new(xq%n)); gids_new = 0
!  n = 0
!  !
!  do j = 1, size(gid_first)
!     gid = gid_first(j)
!     call logmsg('***** gid = '//ta([gid])//' beg: '//ta([n+1]))
!     do lid = 1, xq%n
!       q => xq%get_quad(lid)
!       if (q%get_flag(active=LDUM)) then
!         if (q%gid_prent == gid) then
!           n = n + 1
!           gids_new(q%gid) = n
!         end if
!       end if
!     end do
!     call logmsg('***** gid = '//ta([gid])//' end: '//ta([n]))
!  end do
!  !
!  do lid = 1, xq%n
!    q => xq%get_quad(lid)
!    if (q%get_flag(active=LDUM)) then
!      if (gids_new(q%gid) == 0) then
!        n = n + 1
!        gids_new(q%gid) = n
!      end if
!    end if
!  end do
!  !
!  ! replace the global id
!  do lid = 1, xq%n
!    q => xq%get_quad(lid)
!    if (q%get_flag(active=LDUM)) then
!      gid = gids_new(q%gid)
!      call q%get_bb(child_bbi=bb)
!      call get_mask(xid, q%gid, bb, mask)
!      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!        if (mask(jc,jr) == 1) then
!          xid(ic,ir) = -gid
!        end if
!      end do; end do
!    end if
!  end do
!  xid = abs(xid)
!  !
!  call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
!  call hdrg_gid%write(f_gid_out)
!  if (lwriteximbal) then
!    call grid_load_imbalance(xid, xid_mv, weight, imbal, np, ximbal)
!    f = trim(f_gid_out)//'_imbal'
!    call writeflt(f, ximbal, nc, nr, xll, yll, cs_gid, R4ZERO)
!  end if
!  !
!  return
!end subroutine quad_grid_balancing

subroutine quad_grid_balancing()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- local
  !
  logical, parameter :: lwriteximbal = .true.
  !
  type(tQuad), pointer :: q_new => null(), q_nbr => null()
  type(tBb), pointer :: bbip => null()
  type(tBb)  :: bbi_gid, bbi_new, bb_nbr
  type(tBbX) :: bbx_gid, bbx_hid, bbx_q
  type(tHdr), pointer :: hdrg_hid => null()
  type(tHdrHdr), pointer :: hdr => null()
  type(tCSV), pointer    :: csv => null()
  type(tMetis), pointer :: met => null()
  type(tChaco), pointer :: chac => null()
  type(tBb), dimension(:), pointer :: bbi_arr
  !
  logical :: lgidfirst, ldone, lmetis
  character(len=MXSLEN) :: f
  integer(I4B) :: weight_mv, gid_max_loc, iter, j, iopt, ipart, n_disjoint
  integer(I4B) :: lid_new, gid_new, iact, nact, nmet, np_met, lid_smallest, lid_nbr
  integer(I4B) :: np_strt, np_tgt, np, np_full, np_diff, ip, p, n_wgt
  integer(I4B), dimension(:,:), allocatable :: mask, part, xid_split
  integer(I4B), dimension(:,:), allocatable :: weight, qweight
  integer(I4B), dimension(:), allocatable :: gid_first, wgt_sort_lid
  integer(I4B), dimension(:), allocatable :: gids_new, gid2part, xid_nbr
  real(R4B), dimension(:), allocatable :: wgt_sort
  real(R4B), dimension(:,:), allocatable :: ximbal
  real(R8B) :: wgt_avg, wgt_tgt, wgt_max, wgt_tot, wgt_q, wgt_qnbr, imbal
  real(R8B) :: wgt_min, wgt_new
! ------------------------------------------------------------------------------
  !
  allocate(hdrg_gid)
  call hdrg_gid%read_full_grid(f_gid_in)
  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
  nc =  size(xid,1); nr = size(xid,2)
  if (.false.) then
    do ir = 1, nr; do ic = 1, nc
      if (xid(ic,ir) /= xid_mv) then
        if (xid(ic,ir) == 6) xid(ic,ir) = xid_mv
       end if
    end do; end do
  end if
  !
  ! read the weights
  if (len_trim(f_weight) > 0) then
    allocate(hdrg)
    call hdrg%read_full_grid(f_weight)
    call hdrg%get_grid(xi4=weight, mvi4=weight_mv)
    call hdrg%clean(); deallocate(hdrg); hdrg => null()
  else
    allocate(weight(nc,nr)); weight = 1
  end if
  !
  if (len_trim(gid_first_num) > 0) then
    lgidfirst = .true.
    call parse_line(s=gid_first_num, i4a=gid_first, token_in=',')
  else
    lgidfirst = .false.
  end if
  !
  ! set bbi and bbx for gid
  bbi_gid%ic0 = 1;   bbi_gid%ic1 = nc
  bbi_gid%ir0 = 1;   bbi_gid%ir1 = nr
  bbi_gid%ncol = nc; bbi_gid%nrow = nr
  hdr => hdrg_gid%hdr; call hdr%get_bbx(bbx_gid)
  xll = bbx_gid%xll; yll = bbx_gid%yll; cs_gid = bbx_gid%cs
  call bbo%set(prent_bbi=bbi_gid, prent_bbx=bbx_gid)
  !
  ! determine the global ids and bounding box
  allocate(bb_gid(gid_max), gids(gid_max), g2lid(gid_max))
  gids = 0; g2lid = 0
  do gid = 1, gid_max
    bb => bb_gid(gid)
    call bb%init()
  end do
  do ir = 1, nr; do ic = 1, nc
    gid = xid(ic,ir)
    if (gid /= xid_mv) then
      gids(gid) = 1
      bb => bb_gid(gid)
      bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
      bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
      bb%ncol = bb%ic1 - bb%ic0 + 1; bb%nrow = bb%ir1 - bb%ir0 + 1
    end if
  end do; enddo
  !
  gid_max_loc = maxval(xid)
  nlid = sum(gids)
  if (nlid > gid_max) then
    call errmsg('Increase gid_max.')
  end if
  allocate(xq)
  call xq%init(nlid, gid_max)
  !
  ! allocate the weight arrays
  allocate(wgt_sort(gid_max), wgt_sort_lid(gid_max))
  !
  lid = 0; wgt_tot = R8ZERO
  do gid = 1, gid_max
    if (gids(gid) == 1) then
      lid = lid + 1
      g2lid(gid) = lid
      q => xq%get_quad(lid)
      bb => bb_gid(gid) ! local bounding box
      bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb%ir1)*cs_gid
      bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
      bbx%cs = cs_gid
      call bbo%set(child_bbi=bb, child_bbx=bbx)
      call q%init(gid=gid, lid=lid, bbo=bbo)
      call get_mask(xid, gid, bb, mask)
      if (allocated(qweight)) deallocate(qweight)
      allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
        if (mask(jc,jr) > 0) then
          qweight(jc,jr) = weight(ic,ir)
        end if
      end do; end do;
      call q%calc_prop(mask=mask, weight=qweight)
      call q%get_prop(weight=wgt_q)
      wgt_tot = wgt_tot + wgt_q
    end if
  end do
  !
  np_tgt = xq%get_number_active(); np_strt = np_tgt
  call grid_load_imbalance(xid, xid_mv, weight, imbal, np_strt)
  call logmsg('Start with load imbalance for '//ta([np_strt])//' parts: '//ta([imbal]))
  !
  gid_new = gid_max_loc; iter = 0
  allocate(xid_split(nc,nr)); xid_split = 0
  do while(.true.)
    iter = iter + 1
    !
    ! fill and sort the weights
    do lid = 1, xq%n
      q => xq%get_quad(lid)
      wgt_sort_lid(lid) = lid
      if (q%get_flag(active=LDUM)) then
        call q%get_prop(weight=wgt_q)
        wgt_sort(lid) = wgt_q 
      else
        wgt_sort(lid) = R4ZERO
      end if
    end do
    !wgt_sort = -wgt_sort
    call quicksort_r(wgt_sort, wgt_sort_lid, xq%n)
    !wgt_sort = -wgt_sort
    ldone = .true.
    n = xq%n; lid_new = n; nmet = 0; np_met = 0
    do i = 1, n
      lid = wgt_sort_lid(i)
      q => xq%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        call q%get_prop(weight=wgt_q)
        ! check
        if (wgt_q /= wgt_sort(i)) then
          call errmsg('Quad_grid_balancing: program error.')
        end if
        !
        wgt_avg = wgt_tot/np_tgt
        wgt_tgt = tgt_imbal * wgt_avg
        !call logmsg('New average: '//ta([wgt_avg]))
        lmetis = .false.
        if (wgt_q > wgt_tgt) then
          np_full = nint(wgt_q/wgt_tgt)
          np_full = max(2,np_full)
          if (np_full > 1) then
            lmetis = .true.
          end if
          np_full = min(np_full, max_np - np_tgt + 1)
          if (np_full <= 1) then
            lmetis = .false.
          end if
        end if
        if (lmetis) then
          ldone = .false.
          ! number of parts
          nmet = nmet + 1
          np_met = np_met + np_full - 1
          np_tgt = np_tgt + np_full - 1 !!! OPTION 1 !!!
          !
          call q%get_bb(child_bbi=bb)
          call get_mask(xid, q%gid, bb, mask)
          if (allocated(qweight)) deallocate(qweight)
          allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
          do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
            jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
            if (mask(jc,jr) > 0) then
              xid_split(ic,ir) = np_full
              qweight(jc,jr) = weight(ic,ir)
            end if
          end do; end do;
          !
          ! full grid METIS
          allocate(met)
          call met%init(qweight, np_full)
          call met%set_opts()
          call met%recur(verbose=.true.)
          if (allocated(part)) deallocate(part)
          allocate(part,source=mask)
          call met%set_ids(part)
          call met%clean(); deallocate(met); met => null()
          !
          do ip = 1, np_full
            lid_new = lid_new + 1; gid_new = gid_new + 1
            q_new => xq%get_quad(lid_new)
            call bbi_new%init()
            do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
              jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
              p = part(jc,jr)
              if (p == ip) then
               ! determine the new bounding box
                bbi_new%ic0 = min(bbi_new%ic0,ic); bbi_new%ic1 = max(bbi_new%ic1,ic)
                bbi_new%ir0 = min(bbi_new%ir0,ir); bbi_new%ir1 = max(bbi_new%ir1,ir)
                bbi_new%ncol = bbi_new%ic1 - bbi_new%ic0 + 1;  bbi_new%nrow = bbi_new%ir1 - bbi_new%ir0 + 1
                if (xid(ic,ir) /= q%gid) then
                  call errmsg('quad_sub_grid_part: program error.')
                else
                  xid(ic,ir) = gid_new
                end if
              end if 
            end do; end do
            !
            bbx%xll = xll + (bbi_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bbi_new%ir1)*cs_gid
            bbx%xur = bbx%xll + bbi_new%ncol*cs_gid; bbx%yur = bbx%yll + bbi_new%nrow*cs_gid
            bbx%cs = cs_gid
            call bbo%set(child_bbi=bbi_new, child_bbx=bbx)
            call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
            if (q%gid_prent > 0) then
              call q_new%set(hlev=q%hlev, gid_prent=q%gid_prent, lid_prent=q%lid_prent)
            else
              call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid)
            end if
            call get_mask(xid, gid_new, bbi_new, mask)
            if (allocated(qweight)) deallocate(qweight)
            allocate(qweight(bbi_new%ncol,bbi_new%nrow)); qweight = I4ZERO
            do ir = bbi_new%ir0, bbi_new%ir1; do ic = bbi_new%ic0, bbi_new%ic1
              jr = ir - bbi_new%ir0 + 1; jc = ic - bbi_new%ic0 + 1
              if (mask(jc,jr) > 0) then
                qweight(jc,jr) = weight(ic,ir)
              end if
            end do; end do;
            call q_new%calc_prop(mask=mask, weight=qweight)
            g2lid(gid_new) = lid_new
            xq%n = lid_new
          end do !ip
          !
          call q%set_flag(active=.false.)
        end if
      end if
    end do
    ! 
    !np_tgt = np_tgt + np_met !!! OPTION 2 !!!
    
    call logmsg('***** iteration '//ta([iter])//' *****')
    call logmsg('Number of times METIS was called: '//ta([nmet]))
    call logmsg('Number of new METIS parts:        '//ta([np_met]))
    call grid_load_imbalance(xid, xid_mv, weight, imbal, np)
    call logmsg('Overall load imbalance for '//ta([np])//' parts: '//ta([imbal]))
    if (ldone) then
      call logmsg('Nothing to be done...')
      exit
    end if
    if (imbal <= tgt_imbal) then
      call logmsg('Target load imbalance reached...')
      exit
    end if
    if (iter >= max_iter) then
      call logmsg('No convergence, maximum iterations reached: '//ta([max_iter]))
      exit
    end if
    if (iopt == 1) then
      if (np >= max_np) then
        call logmsg('Maximum number of partitions reached: '//ta([max_np]))
        exit
      end if
    end if
  end do
  !
  np_diff = np-np_strt
  call logmsg('Done, generated '//ta([np_diff])//' parts ('// &
    ta([100.*real(np_diff,R4B)/real(np,R4B)],'(f7.2)')//' % increment)')
  !
  ! renumber
  allocate(gids_new(gid_max)); gids_new = 0
  n = 0
  !
  do j = 1, size(gid_first)
     gid = gid_first(j)
     call logmsg('***** gid = '//ta([gid])//' beg: '//ta([n+1]))
     do lid = 1, xq%n
       q => xq%get_quad(lid)
       if (q%get_flag(active=LDUM)) then
         if (q%gid_prent == gid) then
           n = n + 1
           gids_new(q%gid) = n
         end if
       end if
     end do
     call logmsg('***** gid = '//ta([gid])//' end: '//ta([n]))
  end do
  !
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) then
      if (gids_new(q%gid) == 0) then
        n = n + 1
        gids_new(q%gid) = n
      end if
    end if
  end do
  !
  ! replace the global id
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) then
      gid = gids_new(q%gid)
      call q%get_bb(child_bbi=bb)
      call get_mask(xid, q%gid, bb, mask)
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
        if (mask(jc,jr) == 1) then
          xid(ic,ir) = -gid
        end if
      end do; end do
    end if
  end do
  xid = abs(xid)
  !
  call grid_load_imbalance(xid, xid_mv, weight, imbal, np, wgt_max=wgt_max)
  call logmsg('Overall load imbalance for '//ta([np])// &
    ' parts: '//ta([imbal],'(f6.2)')//', max weight: '//ta([wgt_max]))
  
  call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
  call hdrg_gid%write(trim(f_gid_out)//'_np'//ta([np]))
  if (lwriteximbal) then
    call grid_load_imbalance(xid, xid_mv, weight, imbal, np, ximbal)
    f = trim(f_gid_out)//'_imbal_np'//ta([np])
    call writeflt(f, ximbal, nc, nr, xll, yll, cs_gid, R4ZERO)
  end if
  call hdrg_gid%replace_grid(xi4=xid_split, mvi4=0)
  call hdrg_gid%write(trim(f_gid_out)//'_split_np'//ta([np]))
  !
  return
end subroutine quad_grid_balancing

!subroutine quad_grid_balancing()
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------  
!! -- local
!  !
!  logical, parameter :: lwriteximbal = .true.
!  !
!  type(tQuad), pointer :: q_new => null(), q_nbr => null()
!  type(tBb), pointer :: bbip => null()
!  type(tBb)  :: bbi_gid, bbi_new, bb_nbr
!  type(tBbX) :: bbx_gid, bbx_hid, bbx_q
!  type(tHdr), pointer :: hdrg_hid => null()
!  type(tHdrHdr), pointer :: hdr => null()
!  type(tCSV), pointer    :: csv => null()
!  type(tMetis), pointer :: met => null()
!  type(tChaco), pointer :: chac => null()
!  type(tBb), dimension(:), pointer :: bbi_arr
!  !
!  !real(R8B), parameter :: wgt_fac = 0.7d0
!  !real(R8B), parameter :: wgt_per = 0.05d0
!  real(R8B), parameter :: wgt_fac = 0.7d0
!  real(R8B), parameter :: wgt_per = 0.05d0
!  
!  logical :: lgidfirst, ldone, lmetis, luse_chaco, lcluster
!  character(len=MXSLEN) :: f
!  integer(I4B) :: weight_mv, gid_max_loc, iter, j, iopt, ipart, n_disjoint
!  integer(I4B) :: lid_new, gid_new, iact, nact, nmet, np_met, lid_smallest, lid_nbr
!  integer(I4B) :: np_strt, np_tgt, np, np_full, np_diff, ip, p, area, n_wgt
!  integer(I4B), dimension(:,:), allocatable :: mask, part, xid_split
!  integer(I4B), dimension(:,:), allocatable :: weight, qweight, ini_weight, unique
!  integer(I4B), dimension(:), allocatable :: gid_first, wgt_sort_lid
!  integer(I4B), dimension(:), allocatable :: gids_new, gid2part, xid_nbr
!  real(R4B), dimension(:), allocatable :: wgt_sort
!  real(R4B), dimension(:,:), allocatable :: ximbal
!  real(R8B) :: wgt_avg, wgt_tgt, wgt_max, wgt_tot, wgt_q, wgt_qnbr, imbal, area_q, area_max
!  real(R8B) :: wgt_min, wgt_new
!! ------------------------------------------------------------------------------
!  
!  area_max = 12.5d0 !km2
!  area_max = area_max * 1000.d0 * 1000.d0 ! m2
!  !
!  iopt = 2
!  luse_chaco = .false.
!  lcluster = .false.
!  !
!  allocate(hdrg_gid)
!  call hdrg_gid%read_full_grid(f_gid_in)
!  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
!  nc =  size(xid,1); nr = size(xid,2)
!  if (.false.) then
!    do ir = 1, nr; do ic = 1, nc
!      if (xid(ic,ir) /= xid_mv) then
!        if (xid(ic,ir) == 6) xid(ic,ir) = xid_mv
!       end if
!    end do; end do
!  end if
!  !
!  ! read the weights
!  if (len_trim(f_weight) > 0) then
!    allocate(hdrg)
!    call hdrg%read_full_grid(f_weight)
!    call hdrg%get_grid(xi4=weight, mvi4=weight_mv)
!    call hdrg%clean(); deallocate(hdrg); hdrg => null()
!  else
!    allocate(weight(nc,nr)); weight = 1
!  end if
!  !
!  if (len_trim(gid_first_num) > 0) then
!    lgidfirst = .true.
!    call parse_line(s=gid_first_num, i4a=gid_first, token_in=',')
!  else
!    lgidfirst = .false.
!  end if
!  !
!  ! set bbi and bbx for gid
!  bbi_gid%ic0 = 1;   bbi_gid%ic1 = nc
!  bbi_gid%ir0 = 1;   bbi_gid%ir1 = nr
!  bbi_gid%ncol = nc; bbi_gid%nrow = nr
!  hdr => hdrg_gid%hdr; call hdr%get_bbx(bbx_gid)
!  xll = bbx_gid%xll; yll = bbx_gid%yll; cs_gid = bbx_gid%cs
!  call bbo%set(prent_bbi=bbi_gid, prent_bbx=bbx_gid)
!  !
!  ! determine the global ids and bounding box
!  allocate(bb_gid(gid_max), gids(gid_max), g2lid(gid_max))
!  gids = 0; g2lid = 0
!  do gid = 1, gid_max
!    bb => bb_gid(gid)
!    call bb%init()
!  end do
!  do ir = 1, nr; do ic = 1, nc
!    gid = xid(ic,ir)
!    if (gid /= xid_mv) then
!      gids(gid) = 1
!      bb => bb_gid(gid)
!      bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
!      bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
!      bb%ncol = bb%ic1 - bb%ic0 + 1; bb%nrow = bb%ir1 - bb%ir0 + 1
!    end if
!  end do; enddo
!  !
!  gid_max_loc = maxval(xid)
!  nlid = sum(gids)
!  if (nlid > gid_max) then
!    call errmsg('Increase gid_max.')
!  end if
!  allocate(xq)
!  call xq%init(nlid, gid_max)
!  !
!  ! allocate the weight arrays
!  allocate(wgt_sort(gid_max), wgt_sort_lid(gid_max))
!  !
!  lid = 0; wgt_tot = R8ZERO
!  do gid = 1, gid_max
!    if (gids(gid) == 1) then
!      lid = lid + 1
!      g2lid(gid) = lid
!      q => xq%get_quad(lid)
!      bb => bb_gid(gid) ! local bounding box
!      bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb%ir1)*cs_gid
!      bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
!      bbx%cs = cs_gid
!      call bbo%set(child_bbi=bb, child_bbx=bbx)
!      call q%init(gid=gid, lid=lid, bbo=bbo)
!      call get_mask(xid, gid, bb, mask)
!      if (allocated(qweight)) deallocate(qweight)
!      allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
!      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!        if (mask(jc,jr) > 0) then
!          qweight(jc,jr) = weight(ic,ir)
!        end if
!      end do; end do;
!      call q%calc_prop(mask=mask, weight=qweight)
!      call q%get_prop(weight=wgt_q)
!      wgt_tot = wgt_tot + wgt_q
!    end if
!  end do
!  !
!  if (.true.) then
!    allocate(ini_weight(nc,nr)); ini_weight = 0
!    do lid = 1, xq%n
!      q => xq%get_quad(lid)
!      if (q%get_flag(active=LDUM)) then
!        call q%get_bb(child_bbi=bb)
!        call get_mask(xid, q%gid, bb, mask)
!        call q%get_prop(weight=wgt_q)
!        do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!          jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!          if (mask(jc,jr) == 1) then
!            ini_weight(ic,ir) = wgt_q
!          end if
!        end do; end do
!      end if
!    end do
!    f = trim(f_gid_out)//'_ini_weight'
!    call writeflt(f, ini_weight, nc, nr, xll, yll, cs_gid, I4ZERO)
!    deallocate(ini_weight)
!  end if
!  
!  np_tgt = xq%get_number_active(); np_strt = np_tgt
!  call grid_load_imbalance(xid, xid_mv, weight, imbal, np_strt)
!  call logmsg('Start with load imbalance for '//ta([np_strt])//' parts: '//ta([imbal]))
!  !
!  gid_new = gid_max_loc; iter = 0
!  allocate(xid_split(nc,nr)); xid_split = 0
!  do while(.true.)
!    iter = iter + 1
!    !
!    ! fill and sort the weights
!    do lid = 1, xq%n
!      q => xq%get_quad(lid)
!      wgt_sort_lid(lid) = lid
!      if (q%get_flag(active=LDUM)) then
!        call q%get_prop(weight=wgt_q)
!        wgt_sort(lid) = wgt_q 
!      else
!        wgt_sort(lid) = R4ZERO
!      end if
!    end do
!    if (iopt == 2) then
!      wgt_sort = -wgt_sort
!      call quicksort_r(wgt_sort, wgt_sort_lid, xq%n)
!      wgt_sort = -wgt_sort
!      i = nint(wgt_per *real(xq%n))
!      wgt_tgt = wgt_sort(i)
!      call logmsg('Target weight: '//ta([wgt_tgt]))
!      !max_np = wgt_tot/wgt_tgt
!      !call logmsg('Target number of parts: '//ta([max_np])// &
!      !  ', max weight: '//ta([wgt_tgt]))
!    elseif(iopt == 3) then
!      wgt_tgt = wgt_tot / max_np
!    else
!      !wgt_sort = -wgt_sort
!      call quicksort_r(wgt_sort, wgt_sort_lid, xq%n)
!      !wgt_sort = -wgt_sort
!    end if
!    ldone = .true.
!    n = xq%n; lid_new = n; nmet = 0; np_met = 0
!    do i = 1, n
!      lid = wgt_sort_lid(i)
!      q => xq%get_quad(lid)
!      if (q%get_flag(active=LDUM)) then
!        call q%get_prop(weight=wgt_q)
!        ! check
!        if (wgt_q /= wgt_sort(i)) then
!          call errmsg('Quad_grid_balancing: program error.')
!        end if
!        !
!        if (iopt == 1) then
!          wgt_avg = wgt_tot/np_tgt
!          wgt_tgt = tgt_imbal * wgt_avg
!          !call logmsg('New average: '//ta([wgt_avg]))
!          lmetis = .false.
!          if (wgt_q > wgt_tgt) then
!            np_full = nint(wgt_q/wgt_tgt)
!            np_full = max(2,np_full)
!            if (np_full > 1) then
!              lmetis = .true.
!            end if
!            np_full = min(np_full, max_np - np_tgt + 1)
!            if (np_full <= 1) then
!              lmetis = .false.
!            end if
!          end if
!        else
!          if (wgt_q > wgt_tgt) then
!            lmetis = .true.
!            np_full = nint(wgt_q/(wgt_fac*wgt_tgt))
!            np_full = max(2,np_full)
!          else
!            lmetis = .false.
!          end if
!          
!          !call q%get_bb(child_bbx=bbx_q)
!          !call q%get_prop(area=area)
!          !area_q = area*bbx_q%cs**2
!          !if (area_q > area_max) then
!          !   lmetis = .true.
!          !   wgt_tgt = wgt_fac*wgt_tot/max_np
!          !   np_full = nint(wgt_q/wgt_tgt)
!          !   if (np_full == 1) then
!          !     lmetis = .false.
!          !   end if
!          !else
!          !  lmetis = .false.
!          !end if
!        end if
!        if (lmetis) then
!          ldone = .false.
!          ! number of parts
!          nmet = nmet + 1
!          np_met = np_met + np_full - 1
!          np_tgt = np_tgt + np_full - 1 !!! OPTION 1 !!!
!          !
!          call q%get_bb(child_bbi=bb)
!          call get_mask(xid, q%gid, bb, mask)
!          if (allocated(qweight)) deallocate(qweight)
!          allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
!          do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!            jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!            if (mask(jc,jr) > 0) then
!              xid_split(ic,ir) = np_full
!              qweight(jc,jr) = weight(ic,ir)
!            end if
!          end do; end do;
!          !
!          ! full grid METIS
!          allocate(met)
!          call met%init(qweight, np_full)
!          call met%set_opts()
!          call met%recur(verbose=.true.)
!          if (allocated(part)) deallocate(part)
!          allocate(part,source=mask)
!          call met%set_ids(part)
!          call met%clean(); deallocate(met); met => null()
!          !
!          do ip = 1, np_full
!            lid_new = lid_new + 1; gid_new = gid_new + 1
!            q_new => xq%get_quad(lid_new)
!            call bbi_new%init()
!            do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!              jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!              p = part(jc,jr)
!              if (p == ip) then
!               ! determine the new bounding box
!                bbi_new%ic0 = min(bbi_new%ic0,ic); bbi_new%ic1 = max(bbi_new%ic1,ic)
!                bbi_new%ir0 = min(bbi_new%ir0,ir); bbi_new%ir1 = max(bbi_new%ir1,ir)
!                bbi_new%ncol = bbi_new%ic1 - bbi_new%ic0 + 1;  bbi_new%nrow = bbi_new%ir1 - bbi_new%ir0 + 1
!                if (xid(ic,ir) /= q%gid) then
!                  call errmsg('quad_sub_grid_part: program error.')
!                else
!                  xid(ic,ir) = gid_new
!                end if
!              end if 
!            end do; end do
!            !
!            bbx%xll = xll + (bbi_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bbi_new%ir1)*cs_gid
!            bbx%xur = bbx%xll + bbi_new%ncol*cs_gid; bbx%yur = bbx%yll + bbi_new%nrow*cs_gid
!            bbx%cs = cs_gid
!            call bbo%set(child_bbi=bbi_new, child_bbx=bbx)
!            call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
!            if (q%gid_prent > 0) then
!              call q_new%set(hlev=q%hlev, gid_prent=q%gid_prent, lid_prent=q%lid_prent)
!            else
!              call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid)
!            end if
!            call get_mask(xid, gid_new, bbi_new, mask)
!            if (allocated(qweight)) deallocate(qweight)
!            allocate(qweight(bbi_new%ncol,bbi_new%nrow)); qweight = I4ZERO
!            do ir = bbi_new%ir0, bbi_new%ir1; do ic = bbi_new%ic0, bbi_new%ic1
!              jr = ir - bbi_new%ir0 + 1; jc = ic - bbi_new%ic0 + 1
!              if (mask(jc,jr) > 0) then
!                qweight(jc,jr) = weight(ic,ir)
!              end if
!            end do; end do;
!            call q_new%calc_prop(mask=mask, weight=qweight)
!            g2lid(gid_new) = lid_new
!            xq%n = lid_new
!          end do !ip
!          !
!          call q%set_flag(active=.false.)
!        end if
!      end if
!    end do
!    ! 
!    !np_tgt = np_tgt + np_met !!! OPTION 2 !!!
!    
!    call logmsg('***** iteration '//ta([iter])//' *****')
!    call logmsg('Number of times METIS was called: '//ta([nmet]))
!    call logmsg('Number of new METIS parts:        '//ta([np_met]))
!    call grid_load_imbalance(xid, xid_mv, weight, imbal, np)
!    call logmsg('Overall load imbalance for '//ta([np])//' parts: '//ta([imbal]))
!    if (iopt >= 2) then
!      ldone = .true.
!    end if
!    if (ldone) then
!      call logmsg('Nothing to be done...')
!      exit
!    end if
!    if (imbal <= tgt_imbal) then
!      call logmsg('Target load imbalance reached...')
!      exit
!    end if
!    if (iter >= max_iter) then
!      call logmsg('No convergence, maximum iterations reached: '//ta([max_iter]))
!      exit
!    end if
!    if (iopt == 1) then
!      if (np >= max_np) then
!        call logmsg('Maximum number of partitions reached: '//ta([max_np]))
!        exit
!      end if
!    end if
!  end do
!  !
!  np_diff = np-np_strt
!  call logmsg('Done, generated '//ta([np_diff])//' parts ('// &
!    ta([100.*real(np_diff,R4B)/real(np,R4B)],'(f7.2)')//' % increment)')
!  !
!  ! renumber
!  allocate(gids_new(gid_max)); gids_new = 0
!  n = 0
!  !
!  do j = 1, size(gid_first)
!     gid = gid_first(j)
!     call logmsg('***** gid = '//ta([gid])//' beg: '//ta([n+1]))
!     do lid = 1, xq%n
!       q => xq%get_quad(lid)
!       if (q%get_flag(active=LDUM)) then
!         if (q%gid_prent == gid) then
!           n = n + 1
!           gids_new(q%gid) = n
!         end if
!       end if
!     end do
!     call logmsg('***** gid = '//ta([gid])//' end: '//ta([n]))
!  end do
!  !
!  ! balancing: remove the small weights
!  if ((iopt >= 2).and.(lcluster)) then
!    call logmsg('Start clustering smallest weights.')
!    do while(.true.)
!      ! fill and sort the weights
!      n_wgt = 0
!      do lid = 1, xq%n
!        q => xq%get_quad(lid)
!        if (q%get_flag(active=LDUM)) then
!          n_wgt = n_wgt + 1
!          wgt_sort_lid(n_wgt) = lid
!          call q%get_prop(weight=wgt_q)
!          wgt_sort(n_wgt) = wgt_q 
!        end if
!      end do
!      !call logmsg('Number of weights: '//ta([n_wgt]))
!      !
!      ! sort from small to large
!      call quicksort_r(wgt_sort, wgt_sort_lid, n_wgt)
!      !
!      ldone = .true.
!      do j = 1, n_wgt
!        lid_smallest = wgt_sort_lid(j)
!        q => xq%get_quad(lid_smallest)
!        call q%get_prop(weight=wgt_q)
!        call q%get_bb(child_bbi=bb)
!        call get_mask(xid, q%gid, bb, mask)
!        !
!        ! find the smallest neightbor
!        call get_neighbors(xid, xid_mv, q%gid, xid_nbr)
!        !
!        wgt_min = huge(0.d0)
!        do i = 1, size(xid_nbr)
!          gid = xid_nbr(i)
!          lid = g2lid(gid)
!          if (lid == 0) then
!            call errmsg('Program error.')
!          else
!            q_nbr => xq%get_quad(lid)
!            call q_nbr%get_prop(weight=wgt_qnbr)
!            if (wgt_qnbr < wgt_min) then
!              wgt_min = wgt_qnbr
!              lid_nbr = lid
!            end if
!          end if
!        end do
!        !
!        ! add the load to the neighbor
!        call q%set_flag(active=.false.)
!        q_nbr => xq%get_quad(lid_nbr)
!        call q_nbr%get_prop(weight=wgt_qnbr)
!        call q_nbr%get_bb(child_bbi=bb_nbr)
!        bb_nbr%ic0 = min(bb_nbr%ic0,bb%ic0); bb_nbr%ic1 = max(bb_nbr%ic1,bb%ic1)
!        bb_nbr%ir0 = min(bb_nbr%ir0,bb%ir0); bb_nbr%ir1 = max(bb_nbr%ir1,bb%ir1)
!        bb_nbr%ncol =  bb_nbr%ic1 -  bb_nbr%ic0 + 1
!        bb_nbr%nrow =  bb_nbr%ir1 -  bb_nbr%ir0 + 1
!        wgt_new = wgt_qnbr + wgt_q
!        if (wgt_new <= wgt_sort(n_wgt)) then
!          !call logmsg('Weight gid '//ta([q_nbr%gid])//': '//ta([wgt_qnbr])// ' -> '//ta([wgt_new]))
!          call q_nbr%set(weight=wgt_qnbr+wgt_q)
!          call q_nbr%set(child_bbi=bb_nbr)
!          !
!          ! overwrite the gid
!          do ir = bb%ir0, bb%ir1
!            do ic = bb%ic0, bb%ic1
!              if (xid(ic,ir) /= xid_mv) then
!                if (xid(ic,ir) == q%gid) then
!                  xid(ic,ir) = q_nbr%gid
!                end if
!              end if
!           end do
!          end do
!          g2lid(q%gid) = 0
!          ldone = .false.; exit
!        end if
!      end do
!      if (ldone) then
!        write(*,*) '**** Done ****'
!        exit
!      end if
!    end do
!    call grid_load_imbalance(xid, xid_mv, weight, imbal, np)
!    call logmsg('Overall load imbalance after clustering '//ta([np])//' parts: '//ta([imbal]))
!    call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
!    call hdrg_gid%write(trim(f_gid_out)//'_tmp_np'//ta([max_np]))
!    !stop
!  end if
!  !
!  do lid = 1, xq%n
!    q => xq%get_quad(lid)
!    if (q%get_flag(active=LDUM)) then
!      if (gids_new(q%gid) == 0) then
!        n = n + 1
!        gids_new(q%gid) = n
!      end if
!    end if
!  end do
!  !
!  ! replace the global id
!  do lid = 1, xq%n
!    q => xq%get_quad(lid)
!    if (q%get_flag(active=LDUM)) then
!      gid = gids_new(q%gid)
!      call q%get_bb(child_bbi=bb)
!      call get_mask(xid, q%gid, bb, mask)
!      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!        if (mask(jc,jr) == 1) then
!          xid(ic,ir) = -gid
!        end if
!      end do; end do
!    end if
!  end do
!  xid = abs(xid)
!  !
!  ! lumped METIS
!  if (iopt >= 2) then
!    !max_np = 200
!    allocate(met)
!    call met%init_lump(ids=xid, nparts=max_np, verbose=.false., weight=weight)
!    if (luse_chaco) then ! chaco
!      allocate(chac)
!      call chac%init_from_metis(met)
!      !call chac%run(f_exe)
!      call chac%set_metis(met)
!      call chac%clean(); deallocate(chac); chac => null()
!    else
!      call met%set_opts() !(niter_in=1000)
!      call met%recur()
!      !call met%kway()
!    end if
!    !
!    allocate(gid2part(gid_max)); gid2part = 0
!    do i = 1, met%nvtxs
!      gid = met%idmapinv(i)
!      gid2part(gid) = met%part(i) + 1
!    end do
!    call met%clean(); deallocate(met); met => null()
!    !
!    do ir = 1, nr; do ic = 1, nc
!      gid = xid(ic,ir)
!      if (gid /= xid_mv) then
!        ipart = gid2part(gid)
!        xid(ic,ir) = ipart
!      end if
!    end do; end do
!  end if ! iopt 2
!  !
!  call quad_repair_disjoint(xid, xid_mv, weight)
!  
!  ! compute disjoint parts
!  allocate(bbi_arr(gid_max))
!!  do i = 1, max_np
!  do i = 1, gid_max
!    call bbi_arr(i)%init()
!  end do
!  do ir = 1, nr; do ic = 1, nc
!    gid = xid(ic,ir)
!    if (gid /= xid_mv) then
!      !ipart = gid2part(gid)
!      !bbi_arr(ipart)%ic0 = min(bbi_arr(ipart)%ic0, ic)
!      !bbi_arr(ipart)%ic1 = max(bbi_arr(ipart)%ic1, ic)
!      !bbi_arr(ipart)%ir0 = min(bbi_arr(ipart)%ir0, ir)
!      !bbi_arr(ipart)%ir1 = max(bbi_arr(ipart)%ir1, ir)
!      !bbi_arr(ipart)%ncol = bbi_arr(ipart)%ic1 - bbi_arr(ipart)%ic0 + 1
!      !bbi_arr(ipart)%nrow = bbi_arr(ipart)%ir1 - bbi_arr(ipart)%ir0 + 1
!      bbi_arr(gid)%ic0 = min(bbi_arr(gid)%ic0, ic)
!      bbi_arr(gid)%ic1 = max(bbi_arr(gid)%ic1, ic)
!      bbi_arr(gid)%ir0 = min(bbi_arr(gid)%ir0, ir)
!      bbi_arr(gid)%ir1 = max(bbi_arr(gid)%ir1, ir)
!      bbi_arr(gid)%ncol = bbi_arr(gid)%ic1 - bbi_arr(gid)%ic0 + 1
!      bbi_arr(gid)%nrow = bbi_arr(gid)%ir1 - bbi_arr(gid)%ir0 + 1
!    end if
!  end do; end do
!  
!  n_disjoint = 0
!  !do ipart = 1, max_np
!  do ipart = 1, gid_max
!    if (g2lid(ipart) == 0) cycle
!    bbip => bbi_arr(ipart)
!    if (allocated(unique)) deallocate(unique)
!    allocate(unique(bbip%ncol,bbip%nrow))
!    unique = 0
!    do ir = bbip%ir0, bbip%ir1; do ic = bbip%ic0, bbip%ic1
!      gid = xid(ic,ir)
!      if (gid /= xid_mv) then
!        if (gid == ipart) then
!          jc = ic - bbip%ic0 + 1; jr = ir - bbip%ir0 + 1
!          unique(jc,jr) = 1
!        end if
!      end if
!    end do; end do
!    call calc_unique(unique, 5, regun, regbb, nreg, idum, 0., 0., 0.)
!    if (nreg > 1) then
!      n_disjoint = n_disjoint + 1
!      do ireg = 1, nreg
!        do jr = regbb(ireg)%ir0, regbb(ireg)%ir1
!          do jc = regbb(ireg)%ic0, regbb(ireg)%ic1
!            if (unique(jc,jr) == 1) then
!              ic = jc + bbip%ic0 - 1
!              ir = jr + bbip%ir0 - 1
!              xid(ic,ir) = -abs(xid(ic,ir))
!            end if
!          end do
!        end do
!      end do
!    end if
!  end do
!  !
!  call logmsg('# disjoint: '//ta([n_disjoint])//'/'//ta([max_np]))
!  !
!  call grid_load_imbalance(xid, xid_mv, weight, imbal, np, wgt_max=wgt_max)
!  call logmsg('Overall load imbalance for '//ta([np])// &
!    ' parts: '//ta([imbal],'(f6.2)')//', max weight: '//ta([wgt_max]))
!  
!  call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
!  call hdrg_gid%write(trim(f_gid_out)//'_np'//ta([np]))
!  if (lwriteximbal) then
!    call grid_load_imbalance(xid, xid_mv, weight, imbal, np, ximbal)
!    f = trim(f_gid_out)//'_imbal_np'//ta([np])
!    call writeflt(f, ximbal, nc, nr, xll, yll, cs_gid, R4ZERO)
!  end if
!  call hdrg_gid%replace_grid(xi4=xid_split, mvi4=0)
!  call hdrg_gid%write(trim(f_gid_out)//'_split_np'//ta([np]))
!  !
!  return
!end subroutine quad_grid_balancing

subroutine quad_get_bb(xid, xid_mv, bba)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- dummy
  integer(I4B), dimension(:,:), intent(in) :: xid
  integer(I4B), intent(in) :: xid_mv
  type(tBB), dimension(:), pointer, intent(inout) :: bba
  ! -- local
  type(tBB), pointer :: bb => null()
  integer(I4B) :: nic, ir, ic, id, id_max
! ------------------------------------------------------------------------------
  if (associated(bba)) deallocate(bba); bba => null()
  !
  ! determine the maximum id
  id_max = 0
  do ir = 1, nr; do ic = 1, nc
    id = xid(ic,ir)
    if (id /= xid_mv) then
      id_max = max(id_max, id)
    end if
  end do; end do
  !
  if (id_max == 0) then
    call errmsg('quad_count_disjoint: maximum id  = 0')
  end if
  allocate(bba(id_max))
  do id = 1, id_max
    bb => bba(id)
    call bb%init()
  end do
  !
  ! determine the bounding box
  do ir = 1, nr; do ic = 1, nc
    id = xid(ic,ir)
    if (id /= xid_mv) then
      bb => bba(id)
      bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
      bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
      bb%ncol = bb%ic1 - bb%ic0 + 1
      bb%nrow = bb%ir1 - bb%ir0 + 1
    end if
  end do; end do
  !
  return
end subroutine quad_get_bb

function quad_count_disjoint(xid, xid_mv) result(n_child)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- dummy
  integer(I4B), dimension(:,:), intent(inout) :: xid
  integer(I4B), intent(in) :: xid_mv
  integer(I4B) :: n_child ! result
  ! -- local
  type(tBB), pointer :: bb => null()
  type(tBB), dimension(:), pointer :: bba => null()
  integer(I4B), dimension(:,:), allocatable :: unique
  integer(I4B) :: id, id_max
! ------------------------------------------------------------------------------
  n_child = 0
  call quad_get_bb(xid, xid_mv, bba)
  id_max = size(bba)
  !
  do id = 1, id_max
    if (allocated(unique)) deallocate(unique)
    bb => bba(id)
    if (bb%defined()) then
      allocate(unique(bb%ncol,bb%nrow))
      unique = 0
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        if (xid(ic,ir) /= xid_mv) then
          if (xid(ic,ir) == id) then
            jc = ic - bb%ic0 + 1; jr = ir - bb%ir0 + 1
            unique(jc,jr) = 1
          end if
        end if
      end do; end do
      call calc_unique(unique, 5, regun, regbb, nreg, idum, 0., 0., 0.)
      n_child = n_child + nreg - 1
    end if
  end do
  !
  ! clean up
  if (allocated(unique)) deallocate(unique)
  if (associated(bba)) deallocate(bba); bba => null()
  if (allocated(unique)) deallocate(unique)
  if (allocated(regun)) deallocate(regun)
  if (allocated(regbb)) deallocate(regbb)
  !
  return
end function quad_count_disjoint

subroutine quad_get_disjoint(xid, xid_mv, weight, renumber, &
  n_child, bba_child, id_child, id_parent, wgt_child, wgt_child_idx)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- dummy
  integer(I4B), dimension(:,:),              intent(inout) :: xid
  integer(I4B),                              intent(in)    :: xid_mv
  integer(I4B), dimension(:,:),              intent(in)    :: weight
  logical,                                   intent(in)    :: renumber
  integer(I4B),                              intent(out)   :: n_child
  type(tBB),    dimension(:),   pointer,     intent(inout) :: bba_child
  integer(I4B), dimension(:),   allocatable, intent(inout) :: id_child
  integer(I4B), dimension(:),   allocatable, intent(inout) :: id_parent
  integer(I4B), dimension(:),   allocatable, intent(inout) :: wgt_child
  integer(I4B), dimension(:),   allocatable, intent(inout) :: wgt_child_idx
  !
  ! -- local
  type(tBB), pointer :: bb => null(), bb_child
  type(tBB), dimension(:), pointer :: bba => null()
  real(R4B), dimension(:), allocatable :: reg_wgt, r4wk
  integer(I4B), dimension(:,:), allocatable :: unique
  integer(I4B), dimension(:), allocatable :: reg_idx
  integer(I4B) :: id_max, id, ir, ic, jr, jc, id_new, nreg, ireg, jreg
! ------------------------------------------------------------------------------
  if (associated(bba_child))    deallocate(bba_child); bba_child => null()
  if (allocated(id_child))      deallocate(id_child)
  if (allocated(id_parent))     deallocate(id_parent)
  if (allocated(wgt_child))     deallocate(wgt_child)
  if (allocated(wgt_child_idx)) deallocate(wgt_child_idx)
  !
  ! get the bounding box
  call quad_get_bb(xid, xid_mv, bba)
  id_max = size(bba)
  !
  n_child = quad_count_disjoint(xid, xid_mv)
  if (n_child == 0) then
    call logmsg('quad_get_disjoint: no disjoint parts found.')
    return
  end if
  !
  ! allocate
  allocate(bba_child(n_child))
  allocate(id_child(n_child)); id_child = 0
  allocate(id_parent(n_child)); id_parent = 0
  allocate(wgt_child(n_child)); wgt_child = 0
  allocate(wgt_child_idx(n_child)); wgt_child_idx = 0
  !
  ! determine the
  id_new = id_max
  n_child = 0
  do id = 1, id_max
    if (allocated(unique)) deallocate(unique)
    bb => bba(id)
    if (.not.bb%defined()) cycle
    !
    allocate(unique(bb%ncol,bb%nrow))
    unique = 0
    do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
      if (xid(ic,ir) /= xid_mv) then
        if (xid(ic,ir) == id) then
          jc = ic - bb%ic0 + 1; jr = ir - bb%ir0 + 1
          unique(jc,jr) = 1
        end if
      end if
    end do; end do
    call calc_unique(unique, 5, regun, regbb, nreg, idum, 0., 0., 0.)
    if (nreg > 1) then
      if (allocated(reg_wgt)) deallocate(reg_wgt)
      if (allocated(reg_idx)) deallocate(reg_idx)
      allocate(reg_wgt(nreg), reg_idx(nreg))
      reg_wgt = R4ZERO
      do ireg = 1, nreg
        reg_idx(ireg) = ireg
        do jr = regbb(ireg)%ir0, regbb(ireg)%ir1
          do jc = regbb(ireg)%ic0, regbb(ireg)%ic1
            if (regun(jc,jr) == ireg) then
              ir = bb%ir0 + jr - 1; ic = bb%ic0 + jc - 1
              reg_wgt(ireg) = reg_wgt(ireg) + real(weight(ic,ir),R4B)
            end if
          end do
        end do
      end do
      !
      ! sort from large to small
      reg_wgt = -reg_wgt
      call quicksort_r(reg_wgt, reg_idx, nreg)
      reg_wgt = -reg_wgt
      !
      ! loop over the largest
      do jreg = 2, nreg
        n_child = n_child + 1
        ireg = reg_idx(jreg)
        if (renumber) then
          id_new = id_new + 1
          !
          id_child(n_child)  = id_new
          id_parent(n_child) = id
          !
          do jr = regbb(ireg)%ir0, regbb(ireg)%ir1
            do jc = regbb(ireg)%ic0, regbb(ireg)%ic1
            if (regun(jc,jr) == ireg) then
              ir = bb%ir0 + jr - 1; ic = bb%ic0 + jc - 1
              xid(ic,ir) = id_new
            end if
          end do; end do
        end if ! renumber
        !
        ! store the weight
        wgt_child(n_child) = reg_wgt(jreg)
        !
        ! store the bounding box
        bb_child => bba_child(n_child)
        bb_child%ic0 = regbb(ireg)%ic0 + bb%ic0 - 1 ! correct
        bb_child%ic1 = regbb(ireg)%ic1 + bb%ic0 - 1 ! correct
        bb_child%ir0 = regbb(ireg)%ir0 + bb%ir0 - 1 ! correct
        bb_child%ir1 = regbb(ireg)%ir1 + bb%ir0 - 1 ! correct
        bb_child%ncol = regbb(ireg)%ncol
        bb_child%nrow = regbb(ireg)%nrow
        !
      end do
    end if
  end do
  !
  ! sort from large to small
  allocate(r4wk(n_child))
  do i = 1, n_child
    wgt_child_idx(i) = i
    r4wk(i) = real(wgt_child(i),R4B)
  end do
  r4wk = -r4wk
  call quicksort_r(r4wk, wgt_child_idx, n_child)
  r4wk = -r4wk
  deallocate(r4wk)
  !
  ! clean up
  if (associated(bba)) deallocate(bba); bba => null()
  if (allocated(unique)) deallocate(unique)
  if (allocated(regun)) deallocate(regun)
  if (allocated(regbb)) deallocate(regbb)
  if (allocated(reg_wgt)) deallocate(reg_wgt)
  if (allocated(reg_idx)) deallocate(reg_idx)
  !
  return
end subroutine quad_get_disjoint

subroutine quad_get_id_weight(xid, xid_mv, weight, id_wgt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- dummy
  integer(I4B), dimension(:,:),              intent(inout) :: xid
  integer(I4B),                              intent(in)    :: xid_mv
  integer(I4B), dimension(:,:),              intent(in)    :: weight
  integer(I4B), dimension(:),   allocatable, intent(inout) :: id_wgt
  !
  ! -- local
  type(tBB), dimension(:), pointer :: bba => null()
  integer(I4B) :: id_max, id, ir, ic
! ------------------------------------------------------------------------------
  if (allocated(id_wgt)) deallocate(id_wgt)
  !
  ! get the bounding box
  call quad_get_bb(xid, xid_mv, bba)
  id_max = size(bba)
  !
  allocate(id_wgt(id_max)); id_wgt = 0
  do id = 1, id_max
    bb => bba(id)
    if (bb%defined()) then
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        if (xid(ic,ir) /= xid_mv) then
          if ((xid(ic,ir) == id)) then
            id_wgt(id) = id_wgt(id) + weight(ic,ir)
          end if
        end if
      end do; end do
    end if
  end do
  !
  ! clean up
  deallocate(bba); bba => null()
  !
  return
end subroutine quad_get_id_weight
  
subroutine quad_repair_disjoint(xid, xid_mv, weight)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- dummy
  integer(I4B), dimension(:,:), intent(inout) :: xid
  integer(I4B), intent(in) :: xid_mv
  integer(I4B), dimension(:,:), intent(in) :: weight
  ! -- local
  type(tBB), pointer :: bb => null()
  type(tBB), dimension(:), pointer :: bba_child => null()
  real(R4B), dimension(:), allocatable :: nbr_wgt
  logical :: found
  integer(I4B), dimension(:), allocatable :: id_child, id_parent, wgt_child
  integer(I4B), dimension(:), allocatable :: wgt_child_idx
  integer(I4B), dimension(:), allocatable :: id_nbr, nbr_idx
  integer(I4B), dimension(:), allocatable :: wgt_id
  integer(I4B) :: n_child, i, j, k, nbr, id_max, wgt, idc, idp, idn, n
! ------------------------------------------------------------------------------
  !
  ! get the weights for the ids
  call quad_get_id_weight(xid, xid_mv, weight, wgt_id)
  id_max = size(wgt_id)
  !
  ! get the disjoint regions
  call quad_get_disjoint(xid, xid_mv, weight, .true., &
    n_child, bba_child, id_child, id_parent, wgt_child, wgt_child_idx)
  !
  if (n_child == 0) then
    call logmsg('No disjoint parts found, returning...')
  else
    call logmsg('Merging '//ta([n_child])//' areas from large to small...')
  end if
  !
  ! loop over the regions from large to small
  do k = 1, n_child
    i = wgt_child_idx(k)
    idc = id_child(i)
    idp = id_parent(i)
    bb => bba_child(i)
    ! 
    ! get the neighbors
    call get_neighbors(xid, xid_mv, id_child(i), id_nbr)
    !
    nbr = size(id_nbr)
    if (allocated(nbr_wgt)) deallocate(nbr_wgt)
    if (allocated(nbr_idx)) deallocate(nbr_idx)
    allocate(nbr_wgt(nbr), nbr_idx(nbr))
    do j = 1, nbr
      idn = id_nbr(j)
      nbr_idx(j) = j
      if (idn > id_max) then
        nbr_wgt(j) = huge(nbr_wgt(j))
      else
        nbr_wgt(j) = wgt_id(idn)
      end if
    end do
    !
    ! sort the neighbors from small to large
    call quicksort_r(nbr_wgt, nbr_idx, nbr)
    !
    ! get smallest neighbor
    found = .false.
    do j = 1, nbr
      idn = id_nbr(nbr_idx(j))
      if (idn <= id_max) then
        found = .true.
        exit
      else
        call logmsg('Searching for non-child neighbor...')
      end if
    end do
    !
    if (.not.found) then
      call errmsg('quad_repair_disjoint: only child neighbors found.')
    end if
    !
    ! merge the smallest regions
    do ir = bb%ir0, bb%ir1
      do ic = bb%ic0, bb%ic1
        if (xid(ic,ir) /= xid_mv) then
          if (xid(ic,ir) == idc) then
            xid(ic,ir) = idn
          end if
        end if
      end do
    end do
    !
    ! update the weights
    wgt_id(idn)  = wgt_id(idn) + wgt_child(i) ! neighbor
    wgt_id(idp)  = wgt_id(idp) - wgt_child(i) ! parent
  end do
  !
  ! clean up
  if (allocated(wgt_id))        deallocate(wgt_id)
  if (associated(bba_child))    deallocate(bba_child); bba_child => null()
  if (allocated(id_child))      deallocate(id_child)
  if (allocated(id_parent))     deallocate(id_parent)
  if (allocated(wgt_child))     deallocate(wgt_child)
  if (allocated(wgt_child_idx)) deallocate(wgt_child_idx)
  if (allocated(nbr_wgt))       deallocate(nbr_wgt)
  if (allocated(nbr_idx))       deallocate(nbr_idx)
  !
  return
end subroutine quad_repair_disjoint
  
!subroutine quad_repair_disjoint(xid, xid_mv, weight)
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!  ! -- dummy
!  integer(I4B), dimension(:,:), intent(inout) :: xid
!  integer(I4B), intent(in) :: xid_mv
!  integer(I4B), dimension(:,:), intent(in) :: weight
!  ! -- local
!  type(tBb), pointer :: bb => null()
!  type(tBb), dimension(:), pointer :: bba => null()
!  real(R4B), dimension(:), allocatable :: reg_wgt, nbr_wgt
!  integer(I4B), dimension(:,:), allocatable :: unique
!  integer(I4B), dimension(:), allocatable :: id_wgt, reg_idx, id_nbr, nbr_idx
!  integer(I4B) :: id_max, nc, nr, ic, ir, jc, jr, id, i, j, ireg, jreg, nbr
!! ------------------------------------------------------------------------------
!  nc = size(xid,1); nr = size(xid,2)
!  id_max = 0
!  do ir = 1, nr; do ic = 1, nc
!    id = xid(ic,ir)
!    if (id /= xid_mv) then
!      id_max = max(id_max, id)
!    end if
!  end do; end do
!  !
!  allocate(bba(id_max), id_wgt(id_max))
!  do i = 1, id_max
!    call bba(i)%init()
!  end do
!  !
!  id_wgt = 0
!  do ir = 1, nr; do ic = 1, nc
!    id = xid(ic,ir)
!    if (id /= xid_mv) then
!      !
!      ! bounding box
!      bba(id)%ic0 = min(bba(id)%ic0, ic)
!      bba(id)%ic1 = max(bba(id)%ic1, ic)
!      bba(id)%ir0 = min(bba(id)%ir0, ir)
!      bba(id)%ir1 = max(bba(id)%ir1, ir)
!      bba(id)%ncol = bba(id)%ic1 - bba(id)%ic0 + 1
!      bba(id)%nrow = bba(id)%ir1 - bba(id)%ir0 + 1
!      !
!      ! weight
!      id_wgt(id) = id_wgt(id) + weight(ic,ir)
!    end if
!  end do; end do
!  !
!  do i = 1, id_max
!    if (allocated(unique)) deallocate(unique)
!    bb => bba(i)
!    !
!    allocate(unique(bb%ncol,bb%nrow))
!    unique = 0
!    do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!      id = xid(ic,ir)
!      if (id /= xid_mv) then
!        if (id == i) then
!          jc = ic - bb%ic0 + 1; jr = ir - bb%ir0 + 1
!          unique(jc,jr) = 1
!        end if
!      end if
!    end do; end do
!    call calc_unique(unique, 5, regun, regbb, nreg, idum, 0., 0., 0.)
!    if (nreg > 1) then
!      if (allocated(reg_wgt)) deallocate(reg_wgt)
!      if (allocated(reg_idx)) deallocate(reg_idx)
!      allocate(reg_wgt(nreg), reg_idx(nreg))
!      reg_wgt = R4ZERO
!      do ireg = 1, nreg
!        reg_idx(ireg) = ireg
!        do jr = regbb(ireg)%ir0, regbb(ireg)%ir1
!          do jc = regbb(ireg)%ic0, regbb(ireg)%ic1
!            if (regun(jc,jr) == ireg) then
!              ir = bb%ir0 + jr - 1; ic = bb%ic0 + jc - 1
!              reg_wgt(ireg) = reg_wgt(ireg) + real(weight(ic,ir),R4B)
!            end if
!          end do
!        end do
!      end do
!      !
!      ! sort from large to small
!      reg_wgt = -reg_wgt
!      call quicksort_r(reg_wgt, reg_idx, nreg)
!      reg_wgt = -reg_wgt
!      !
!      ! first, label the smaller region to id_max + ireg
!      do jreg = 2, nreg
!        ireg = reg_idx(jreg)
!        do jr = regbb(ireg)%ir0, regbb(ireg)%ir1
!          do jc = regbb(ireg)%ic0, regbb(ireg)%ic1
!            if (regun(jc,jr) == ireg) then
!              ir = bb%ir0 + jr - 1; ic = bb%ic0 + jc - 1
!              xid(ic,ir) = id_max + ireg
!            end if
!          end do
!        end do
!      end do
!      !
!      ! loop over the smallest regions
!      do jreg = 2, nreg
!        ireg = reg_idx(jreg)
!        
!        ! get the neighbors
!        call get_neighbors(xid, xid_mv, id_max + ireg, id_nbr)
!        !
!        nbr = size(id_nbr)
!        if (allocated(nbr_wgt)) deallocate(nbr_wgt)
!        if (allocated(nbr_idx)) deallocate(nbr_idx)
!        allocate(nbr_wgt(nbr), nbr_idx(nbr))
!        do j = 1, nbr
!          id = id_nbr(j)
!          nbr_idx(j) = j
!          nbr_wgt(j) = id_wgt(id)
!        end do
!        call quicksort_r(nbr_wgt, nbr_idx, nbr)
!        !
!        ! smallest neighbor
!        id = id_nbr(nbr_idx(1))
!        !
!        ! merge the smallest regions
!        do jr = regbb(ireg)%ir0, regbb(ireg)%ir1
!          do jc = regbb(ireg)%ic0, regbb(ireg)%ic1
!            if (regun(jc,jr) == ireg) then
!              ir = bb%ir0 + jr - 1; ic = bb%ic0 + jc - 1
!              xid(ic,ir) = id
!            end if
!          end do
!        end do
!        !
!        ! update the weights
!        id_wgt(id) = id_wgt(id) + reg_wgt(jreg)
!        id_wgt(i)  = id_wgt(i)  - reg_wgt(jreg)
!      end do
!    end if
!  end do
!  !
!  return
!end subroutine quad_repair_disjoint

subroutine quad_check_sub_grid_part()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  type(tMetis), pointer :: met => null()
  type(tChaco), pointer :: chac => null()
  integer(I4B), dimension(:,:), allocatable :: weight
  integer(I4B) :: weight_mv, nc, nr, np
! ------------------------------------------------------------------------------
  !
  allocate(hdrg_gid)
  call hdrg_gid%read_full_grid(f_gid_in)
  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
  call hdrg_gid%clean(); deallocate(hdrg_gid); hdrg_gid => null()
  nc =  size(xid,1); nr = size(xid,2)

  ! read the weights
  if (len_trim(f_weight) > 0) then
    allocate(hdrg)
    call hdrg%read_full_grid(f_weight)
    call hdrg%get_grid(xi4=weight, mvi4=weight_mv)
    call hdrg%clean(); deallocate(hdrg); hdrg => null()
  else
    allocate(weight(nc,nr)); weight = 1
  end if
  !
  do np = np_beg, np_end, np_step
    allocate(met)
    call met%init_lump(ids=xid, nparts=np, verbose=.false., weight=weight)
    if (luse_chaco) then ! chaco
      allocate(chac)
      call chac%init_from_metis(met)
      call chac%run(f_exe)
      call chac%set_metis(met)
      call chac%clean(); deallocate(chac); chac => null()
    else
      call met%set_opts() !(niter_in=1000)
      call met%recur()
      !call met%kway()
    end if
    call met%clean(); deallocate(met); met => null()
  end do
  !
  ! clean up
  if (allocated(xid)) deallocate(xid)
  if (allocated(weight)) deallocate(weight)
  !
  return
end subroutine quad_check_sub_grid_part

subroutine quad_sub_grid_part()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
  !
  logical, parameter :: ldebug = .false.
  !
  type(tQuad), pointer :: q_new => null(), q_nbr => null()
  type(tBb), pointer :: bbip => null()
  type(tBb)  :: bbi_gid, bbi_new, bb_nbr
  type(tBbX) :: bbx_gid, bbx_hid, bbx_q
  type(tHdrHdr), pointer :: hdr => null()
  type(tCSV), pointer    :: csv => null()
  type(tMetis), pointer :: met => null()
  type(tChaco), pointer :: chac => null()
  !
  logical :: lgidsep, lmetis
  character(len=MXSLEN) :: f
  real(R4B), dimension(:), allocatable :: wgt_sort
  real(R4B), dimension(:,:), allocatable :: ximbal
  real(R8B) :: wgt_avg, wgt_tgt, wgt_max, wgt_tot, wgt_q, wgt_qnbr, imbal
  real(R8B) :: wgt_min, wgt_new, wgt_sep, wgt_lump
  integer(I4B), dimension(:,:), allocatable :: mask, part, xid_split, xid_sep, mask_sep
  integer(I4B), dimension(:,:), allocatable :: weight, qweight, ini_weight
  integer(I4B), dimension(:), allocatable :: gid_sep, wgt_sort_lid
  integer(I4B), dimension(:), allocatable :: gids_new, gid2part, xid_nbr
  integer(I4B), dimension(1) :: loc
  integer(I4B) :: weight_mv, gid_max_loc, j, ipart, n_disjoint
  integer(I4B) :: lid_new, gid_new, gid_offset
  integer(I4B) :: np, np_full, nmet, np_met, np_diff, ip, p, n_wgt, np_strt
  integer(I4B) :: np_lump, np_tot, offset
! ------------------------------------------------------------------------------
  !
  call logmsg('***** wgt_fac: '//ta([wgt_fac]))
  call logmsg('***** wgt_per: '//ta([wgt_per]))
  !
  allocate(hdrg_gid)
  call hdrg_gid%read_full_grid(f_gid_in)
  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
  nc =  size(xid,1); nr = size(xid,2)
  !
  ! read the weights
  if (len_trim(f_weight) > 0) then
    allocate(hdrg)
    call hdrg%read_full_grid(f_weight)
    call hdrg%get_grid(xi4=weight, mvi4=weight_mv)
    call hdrg%clean(); deallocate(hdrg); hdrg => null()
  else
    allocate(weight(nc,nr)); weight = 1
  end if
  !
  if (len_trim(gid_separate) > 0) then
    lgidsep = .true.
    call parse_line(s=gid_separate, i4a=gid_sep, token_in=',')
  else
    lgidsep = .false.
  end if
  !
  ! set bbi and bbx for gid
  bbi_gid%ic0 = 1;   bbi_gid%ic1 = nc
  bbi_gid%ir0 = 1;   bbi_gid%ir1 = nr
  bbi_gid%ncol = nc; bbi_gid%nrow = nr
  hdr => hdrg_gid%hdr; call hdr%get_bbx(bbx_gid)
  xll = bbx_gid%xll; yll = bbx_gid%yll; cs_gid = bbx_gid%cs
  call bbo%set(prent_bbi=bbi_gid, prent_bbx=bbx_gid)
  !
  ! determine the global ids and bounding box
  allocate(bb_gid(gid_max), gids(gid_max), g2lid(gid_max))
  gids = 0; g2lid = 0
  do gid = 1, gid_max
    bb => bb_gid(gid)
    call bb%init()
  end do
  do ir = 1, nr; do ic = 1, nc
    gid = xid(ic,ir)
    if (gid /= xid_mv) then
      gids(gid) = 1
      bb => bb_gid(gid)
      bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
      bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
      bb%ncol = bb%ic1 - bb%ic0 + 1; bb%nrow = bb%ir1 - bb%ir0 + 1
    end if
  end do; enddo
  !
  gid_max_loc = maxval(xid)
  nlid = sum(gids)
  if (nlid > gid_max) then
    call errmsg('Increase gid_max.')
  end if
  allocate(xq)
  call xq%init(nlid, gid_max)
  !
  ! allocate the weight arrays
  allocate(wgt_sort(gid_max), wgt_sort_lid(gid_max))
  !
  lid = 0; wgt_tot = R8ZERO
  do gid = 1, gid_max
    if (gids(gid) == 1) then
      lid = lid + 1
      g2lid(gid) = lid
      q => xq%get_quad(lid)
      bb => bb_gid(gid) ! local bounding box
      bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb%ir1)*cs_gid
      bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
      bbx%cs = cs_gid
      call bbo%set(child_bbi=bb, child_bbx=bbx)
      call q%init(gid=gid, lid=lid, bbo=bbo)
      call get_mask(xid, gid, bb, mask)
      if (allocated(qweight)) deallocate(qweight)
      allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
        if (mask(jc,jr) > 0) then
          qweight(jc,jr) = weight(ic,ir)
        end if
      end do; end do;
      call q%calc_prop(mask=mask, weight=qweight)
      call q%get_prop(weight=wgt_q)
      wgt_tot = wgt_tot + wgt_q
    end if
  end do
  !
  call logmsg('Total weight: '//ta([int(wgt_tot)]))
  !
  gid_new = gid_max_loc
  if (lgidsep) then
    allocate(xid_split(nc,nr)); xid_split = 0
    wgt_sep = R8ZERO
    allocate(xid_sep(nc,nr)); xid_sep = xid_mv
    allocate(mask_sep(nc,nr)); mask_sep = 0
    np_tot = 0
    do i = 1, size(gid_sep)
      gid = gid_sep(i); lid = g2lid(gid)
      q => xq%get_quad(lid)
      call q%get_prop(weight=wgt_q)
      !
      np_full =  max(nint(real(max_np,R8B)*wgt_q/wgt_tot),1)
      !
      wgt_sep = wgt_sep + wgt_q
      call q%get_bb(child_bbi=bb)
      call get_mask(xid, q%gid, bb, mask)
      if (allocated(qweight)) deallocate(qweight)
      allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
        if (mask(jc,jr) > 0) then
          xid_split(ic,ir) = np_full
          qweight(jc,jr) = weight(ic,ir)
        end if
      end do; end do;
      !
      ! full grid METIS
      if (np_full == 1) then
        call logmsg('***** Skipping METIS full partitioning! *****')
        np_tot = np_tot + 1
        do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
          jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
          if (mask(jc,jr) > 0) then
            xid_sep(ic,ir) = np_tot
          end if
        end do; end do
      else
        allocate(met)
        call met%init(qweight, np_full)
        call met%set_opts()
        call met%recur()
        if (allocated(part)) deallocate(part)
        allocate(part,source=mask)
        call met%set_ids(part)
        call met%clean(); deallocate(met); met => null()
        !
        do ip = 1, np_full
          np_tot = np_tot + 1
          gid_new = gid_new + 1
          do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
            jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
            p = part(jc,jr)
            if (p == ip) then
              if (xid(ic,ir) /= gid) then
                call errmsg('quad_sub_grid_part: program error.')
              else
                xid_sep(ic,ir) = np_tot
              end if
            end if 
          end do; end do
        end do
      end if
      !
      ! deactive the quad and remove the ids's, store mask
      call q%set_flag(active=.false.)
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        if (xid(ic,ir) == gid) then
          mask_sep(ic,ir) = gid
          xid(ic,ir) = xid_mv
        end if
      end do; end do
    end do
    wgt_lump = wgt_tot - wgt_sep
    call logmsg('Weight for sub-optimal partitioning: '//ta([int(wgt_lump)]))
    call logmsg('Estimate # lumped partitions: '//ta([nint(real(max_np)*wgt_lump/wgt_tot)]))
    np_lump = max_np - np_tot
    call logmsg('Realized # lumped partitions: '//ta([np_lump]))
  else
    np_lump = max_np
  end if
  !
  if (ldebug) then
    allocate(ini_weight(nc,nr)); ini_weight = 0
    do lid = 1, xq%n
      q => xq%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        call q%get_bb(child_bbi=bb)
        call get_mask(xid, q%gid, bb, mask)
        call q%get_prop(weight=wgt_q)
        do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
          jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
          if (mask(jc,jr) == 1) then
            ini_weight(ic,ir) = wgt_q
          end if
        end do; end do
      end if
    end do
    f = trim(f_gid_out)//'_ini_weight'
    call writeflt(f, ini_weight, nc, nr, xll, yll, cs_gid, I4ZERO)
    deallocate(ini_weight)
  end if
  
  np_strt = xq%get_number_active()
  call grid_load_imbalance(xid, xid_mv, weight, imbal, np_strt)
  call logmsg('Start with load imbalance for '//ta([np_strt])//' parts: '//ta([imbal]))
  !
  gid_new = gid_max_loc
  !
  ! fill and sort the weights
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    wgt_sort_lid(lid) = lid
    if (q%get_flag(active=LDUM)) then
      call q%get_prop(weight=wgt_q)
      wgt_sort(lid) = wgt_q 
    else
      wgt_sort(lid) = R4ZERO
    end if
  end do
  wgt_sort = -wgt_sort
  call quicksort_r(wgt_sort, wgt_sort_lid, xq%n)
  wgt_sort = -wgt_sort
  i = nint(wgt_per *real(xq%n))
  wgt_tgt = wgt_sort(i)
  call logmsg('Target weight: '//ta([wgt_tgt]))
  
  n = xq%n; lid_new = n; nmet = 0; np_met = 0
  do i = 1, n
    lid = wgt_sort_lid(i)
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) then
      call q%get_prop(weight=wgt_q)
      ! check
      if (wgt_q /= wgt_sort(i)) then
        call errmsg('Quad_grid_balancing: program error.')
      end if
      !
      if (wgt_q > wgt_tgt) then
        lmetis = .true.
        np_full = nint(wgt_q/(wgt_fac*wgt_tgt))
        np_full = max(2,np_full)
      else
        lmetis = .false.
      end if
      if (lmetis) then
        ! number of parts
        !
        nmet = nmet + 1
        np_met = np_met + np_full - 1
        !
        call q%get_bb(child_bbi=bb)
        call get_mask(xid, q%gid, bb, mask)
        if (allocated(qweight)) deallocate(qweight)
        allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
        do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
          jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
          if (mask(jc,jr) > 0) then
            xid_split(ic,ir) = np_full
            qweight(jc,jr) = weight(ic,ir)
          end if
        end do; end do;
        !
        ! full grid METIS
        allocate(met)
        call met%init(qweight, np_full)
        call met%set_opts()
        call met%recur(verbose=.true.)
        if (allocated(part)) deallocate(part)
        allocate(part,source=mask)
        call met%set_ids(part)
        call met%clean(); deallocate(met); met => null()
        !
        do ip = 1, np_full
          lid_new = lid_new + 1; gid_new = gid_new + 1
          q_new => xq%get_quad(lid_new)
          call bbi_new%init()
          do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
            jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
            p = part(jc,jr)
            if (p == ip) then
             ! determine the new bounding box
              bbi_new%ic0 = min(bbi_new%ic0,ic); bbi_new%ic1 = max(bbi_new%ic1,ic)
              bbi_new%ir0 = min(bbi_new%ir0,ir); bbi_new%ir1 = max(bbi_new%ir1,ir)
              bbi_new%ncol = bbi_new%ic1 - bbi_new%ic0 + 1;  bbi_new%nrow = bbi_new%ir1 - bbi_new%ir0 + 1
              if (xid(ic,ir) /= q%gid) then
                call errmsg('quad_sub_grid_part: program error.')
              else
                xid(ic,ir) = gid_new
              end if
            end if 
          end do; end do
          !
          bbx%xll = xll + (bbi_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bbi_new%ir1)*cs_gid
          bbx%xur = bbx%xll + bbi_new%ncol*cs_gid; bbx%yur = bbx%yll + bbi_new%nrow*cs_gid
          bbx%cs = cs_gid
          call bbo%set(child_bbi=bbi_new, child_bbx=bbx)
          call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
          if (q%gid_prent > 0) then
            call q_new%set(hlev=q%hlev, gid_prent=q%gid_prent, lid_prent=q%lid_prent)
          else
            call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid)
          end if
          call get_mask(xid, gid_new, bbi_new, mask)
          if (allocated(qweight)) deallocate(qweight)
          allocate(qweight(bbi_new%ncol,bbi_new%nrow)); qweight = I4ZERO
          do ir = bbi_new%ir0, bbi_new%ir1; do ic = bbi_new%ic0, bbi_new%ic1
            jr = ir - bbi_new%ir0 + 1; jc = ic - bbi_new%ic0 + 1
            if (mask(jc,jr) > 0) then
              qweight(jc,jr) = weight(ic,ir)
            end if
          end do; end do;
          call q_new%calc_prop(mask=mask, weight=qweight)
          g2lid(gid_new) = lid_new
          xq%n = lid_new
        end do !ip
        !
        call q%set_flag(active=.false.)
      end if
    end if
  end do
  ! 
  call logmsg('Number of times METIS was called: '//ta([nmet]))
  call logmsg('Number of new METIS parts:        '//ta([np_met]))
  call grid_load_imbalance(xid, xid_mv, weight, imbal, np)
  call logmsg('Overall load imbalance for '//ta([np])//' parts: '//ta([imbal],'(f6.2)'))
  !
  np_diff = np-np_strt
  call logmsg('Done, generated '//ta([np_diff])//' parts ('// &
    ta([100.*real(np_diff,R4B)/real(np,R4B)],'(f7.2)')//' % increment)')
  !
  ! renumber
  allocate(gids_new(gid_max)); gids_new = 0
  n = 0
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) then
      if (gids_new(q%gid) == 0) then
        n = n + 1
        gids_new(q%gid) = n
      end if
    end if
  end do
  !
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) then
      gid = gids_new(q%gid)
      call q%get_bb(child_bbi=bb)
      call get_mask(xid, q%gid, bb, mask)
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
        if (mask(jc,jr) == 1) then
          xid(ic,ir) = gid
        end if
      end do; end do
    end if
  end do
  !
  ! lumped METIS
  allocate(met)
  call met%init_lump(ids=xid, nparts=np_lump, verbose=.false., weight=weight)
  if (luse_chaco) then ! chaco
    allocate(chac)
    call chac%init_from_metis(met)
    call chac%run(f_exe)
    call chac%set_metis(met)
    call chac%clean(); deallocate(chac); chac => null()
  else
    call met%set_opts() !(niter_in=1000)
    call met%recur()
    !call met%kway()
  end if
  !
  allocate(gid2part(gid_max)); gid2part = 0
  do i = 1, met%nvtxs
    gid = met%idmapinv(i)
    gid2part(gid) = met%part(i) + 1
  end do
  call met%clean(); deallocate(met); met => null()
  !
  do ir = 1, nr; do ic = 1, nc
    gid = xid(ic,ir)
    if (gid /= xid_mv) then
      ipart = gid2part(gid)
      xid(ic,ir) = ipart
    end if
  end do; end do
  !
  ! repair for disjoint partitions
  call quad_repair_disjoint(xid, xid_mv, weight)
  !
  ! extra check
  n_disjoint = quad_count_disjoint(xid, xid_mv)
  call logmsg('# disjoint: '//ta([n_disjoint])//'/'//ta([np_lump]))
  !
  if (ldebug) then
    call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
    call hdrg_gid%write(trim(f_gid_out)//'_debug'//ta([np]))
  end if
  
  if (lgidsep) then
    !
    ! first, renumber
    do ir = 1, nr; do ic = 1, nc
      gid = xid(ic,ir)
      if (gid /= xid_mv) then
        xid(ic,ir) = gid + np_tot
      end if
    end do; end do
    !
    ! repair for disjoint partitions
    call quad_repair_disjoint(xid_sep, xid_mv, weight)
    !
    ! second, add the separate partitions
    do i = 1, size(gid_sep)
      gid = gid_sep(i); lid = g2lid(gid)
      q => xq%get_quad(lid)
      call q%get_bb(child_bbi=bb)
      !
      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
        if (mask_sep(ic,ir) == gid) then
          if (xid_sep(ic,ir) == xid_mv) then
            call errmsg('quad_sub_grid_part: program error.')
          else
            xid(ic,ir) = xid_sep(ic,ir)
          end if
        end if
      end do; end do
    end do
  end if
  !
  call grid_load_imbalance(xid, xid_mv, weight, imbal, np, wgt_max=wgt_max)
  call logmsg('Overall load imbalance for '//ta([np])// &
    ' parts: '//ta([imbal],'(f6.2)')//', max weight: '//ta([wgt_max]))
  
  call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
  call hdrg_gid%write(trim(f_gid_out)//'_np'//ta([np]))
  if (lwrite_ximbal) then
    call grid_load_imbalance(xid, xid_mv, weight, imbal, np, ximbal)
    f = trim(f_gid_out)//'_imbal_np'//ta([np])
    call writeflt(f, ximbal, nc, nr, xll, yll, cs_gid, R4ZERO)
  end if
  if (ldebug) then
    call hdrg_gid%replace_grid(xi4=xid_split, mvi4=0)
    call hdrg_gid%write(trim(f_gid_out)//'_split_np'//ta([np]))
  end if
  !
  ! clean up
  call hdrg_gid%clean(); deallocate(hdrg); hdrg => null()
  if (allocated(xid)) deallocate(xid)
  if (allocated(xid_sep)) deallocate(xid_sep)
  if (allocated(mask)) deallocate(mask)
  if (allocated(mask_sep)) deallocate(mask_sep)
  
  return
end subroutine quad_sub_grid_part

!subroutine quad_sub_grid_part()
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------  
!! -- local
!  !
!  logical, parameter :: lwriteximbal = .true.
!  !
!  type(tQuad), pointer :: q_new => null()
!  type(tBb)  :: bbi_gid, bb_new
!  type(tBbX) :: bbx_gid, bbx_hid, bbx_q
!  type(tHdr), pointer :: hdrg_hid => null()
!  type(tHdrHdr), pointer :: hdr => null()
!  type(tCSV), pointer    :: csv => null()
!  type(tMetis), pointer :: met => null()
!  !
!  logical :: lfound, lfirst, llast, lhiera, lgidsep, lok
!  character(len=MXSLEN) :: f
!  character(len=MXSLEN), dimension(:), allocatable :: sa
!  integer(I4B) :: ngid, gid_max_loc, nh, ih, xidh_mv, weight_mv, id, nhid, j, n_old
!  integer(I4B) :: np_rem, np, np_full, np_lump, np_lump_loc, ip, p, lid_new, gid_new, iact
!  integer(I4B), dimension(:,:), allocatable :: mask, part, xidh, hmap
!  integer(I4B), dimension(:,:), allocatable :: weight, qweight, lweight
!  integer(I4B), dimension(:), allocatable :: lev_id, uplev_id, hid, hlev, gid_sep, gids_new
!  real(R4B), dimension(:,:), allocatable :: ximbal
!  real(R8B), dimension(:), allocatable :: wgth
!  real(R8B) :: wgt, wgt_metis, wgt_q, xm, ym, wgt_tot, wgt_tgt, imbal
!! ------------------------------------------------------------------------------
!  !
!  allocate(hdrg_gid)
!  call hdrg_gid%read_full_grid(f_gid_in)
!  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
!  nc =  size(xid,1); nr = size(xid,2)
!  !
!  ! read the weights
!  if (len_trim(f_weight) > 0) then
!    allocate(hdrg)
!    call hdrg%read_full_grid(f_weight)
!    call hdrg%get_grid(xi4=weight, mvi4=weight_mv)
!    call hdrg%clean(); deallocate(hdrg); hdrg => null()
!  else
!    allocate(weight(nc,nr)); weight = 1
!  end if
!  !
!  if (len_trim(gid_separate) > 0) then
!    lgidsep = .true.
!    call parse_line(s=gid_separate, i4a=gid_sep, token_in=',')
!  else
!    lgidsep = .false.
!  end if
!  !
!  ! set bbi and bbx for gid
!  bbi_gid%ic0 = 1;   bbi_gid%ic1 = nc
!  bbi_gid%ir0 = 1;   bbi_gid%ir1 = nr
!  bbi_gid%ncol = nc; bbi_gid%nrow = nr
!  hdr => hdrg_gid%hdr; call hdr%get_bbx(bbx_gid)
!  xll = bbx_gid%xll; yll = bbx_gid%yll; cs_gid = bbx_gid%cs
!  call bbo%set(prent_bbi=bbi_gid, prent_bbx=bbx_gid)
!  !
!  ! determine the global ids and bounding box
!  allocate(bb_gid(gid_max), gids(gid_max), g2lid(gid_max))
!  gids = 0; g2lid = 0
!  do gid = 1, gid_max
!    bb => bb_gid(gid)
!    call bb%init()
!  end do
!  do ir = 1, nr; do ic = 1, nc
!    gid = xid(ic,ir)
!    if (gid /= xid_mv) then
!      gids(gid) = 1
!      bb => bb_gid(gid)
!      bb%ic0 = min(bb%ic0, ic); bb%ic1 = max(bb%ic1, ic)
!      bb%ir0 = min(bb%ir0, ir); bb%ir1 = max(bb%ir1, ir)
!      bb%ncol = bb%ic1 - bb%ic0 + 1; bb%nrow = bb%ir1 - bb%ir0 + 1
!    end if
!  end do; enddo
!  !
!  gid_max_loc = maxval(xid)
!  nlid = sum(gids)
!  if (nlid > gid_max) then
!    call errmsg('Increase gid_max.')
!  end if
!  allocate(xq)
!  call xq%init(nlid, gid_max)
!  !
!  lid = 0; wgt_tot = R8ZERO
!  do gid = 1, gid_max
!    if (gids(gid) == 1) then
!      lid = lid + 1
!      g2lid(gid) = lid
!      q => xq%get_quad(lid)
!      bb => bb_gid(gid) ! local bounding box
!      bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb%ir1)*cs_gid
!      bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
!      bbx%cs = cs_gid
!      call bbo%set(child_bbi=bb, child_bbx=bbx)
!      call q%init(gid=gid, lid=lid, bbo=bbo)
!      call get_mask(xid, gid, bb, mask)
!      if (allocated(qweight)) deallocate(qweight)
!      allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
!      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!        if (mask(jc,jr) > 0) then
!          qweight(jc,jr) = weight(ic,ir)
!        end if
!      end do; end do;
!      call q%calc_prop(mask=mask, weight=qweight)
!      call q%get_prop(weight=wgt_q)
!      wgt_tot = wgt_tot + wgt_q
!    end if
!  end do
!  !
!  ! create the hierarchy
!  if (len_trim(f_hiera_in) > 0) then
!    lhiera = .true.
!    call parse_line(f_hiera_in, sa, token_in=',')
!    nh = size(sa)
!    allocate(hmap(nh,nlid), hlev(nh)); hmap = 0
!  else
!    lhiera = .false.
!    nh = 1
!    allocate(hmap(nh,nlid), hlev(nh)); hmap = 1
!  end if
!  !
!  if (lhiera) then
!    do ih = 1, nh
!      allocate(hdrg_hid)
!      call hdrg_hid%read_full_grid(sa(ih))
!      call hdrg_hid%get_grid(xi4=xidh, mvi4=xidh_mv)
!      hdr => hdrg_hid%hdr; call hdr%get_bbx(bbx_hid)
!      !
!      do lid = 1, xq%n
!        q => xq%get_quad(lid)
!        call q%get_prop(xm=xm, ym=ym)
!        call get_icr(ic, ir, xm, ym, bbx_hid%xll, bbx_hid%yur, bbx_hid%cs)
!        if ((ic > 0).and.(ir> 0)) then
!          id = xidh(ic,ir)
!          if (id /= xidh_mv) then
!            hmap(ih,lid) = id
!          end if
!        else
!          call errmsg('quad_sub_grid_part: could not sample point.')
!        end if
!      end do
!      call hdrg_hid%clean(); deallocate(hdrg_hid); hdrg_hid => null()
!    end do
!    !
!    ! check
!    if (minval(hmap(1,:)) == 0) then
!      call errmsg('quad_sub_grid_part: not all areas are classified.')
!    end if
!    !
!  end if  
!  !
!  ! set the levels
!  do lid = 1, xq%n
!    q => xq%get_quad(lid)
!    call q%set(hlev=hmap(:,lid))
!  end do
!  !
!  wgt_tgt = wgt_tot / npart
!  np_rem = npart
!  !area_tgt = area_tot / 5
!  !
!  do ih = 1, nh
!    ! get the IDs
!    call get_unique(hmap(ih,:), hid)
!    nhid = size(hid)
!    allocate(wgth(nhid)); wgth = R8ZERO
!    !
!    do i = 1, nhid
!      id = hid(i)
!      if (id == 0) cycle
!      do lid = 1, xq%n
!        q => xq%get_quad(lid)
!        if (q%get_flag(active=LDUM)) then
!          if (q%get_hlev(ih) == id) then
!            call q%get_prop(weight=wgt_q)
!            wgth(i) = wgth(i) + wgt_q
!          end if
!        end if
!      end do
!    end do
!    !
!    gid_new = gid_max_loc
!    !
!    ! check the target area
!    do i = 1, nhid
!      id = hid(i)
!      if (wgth(i) == R8ZERO) cycle
!      if (wgth(i) < wgt_tgt) then ! merge the quads for this level
!        lid_new = xq%n + 1; gid_new = gid_new + 1
!        call bb_new%init()
!        ! merge the quads
!        lfound = .false.; lfirst = .true.
!        do lid = 1, xq%n
!          q => xq%get_quad(lid)
!          if (q%get_flag(active=LDUM)) then
!            if (q%get_hlev(ih) == id) then
!              lfound = .true.
!              call q%set_flag(active=.false.)
!              call q%get_bb(child_bbi=bb)
!              ! change the id
!              do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                if (xid(ic,ir) == q%gid) then
!                  xid(ic,ir) = gid_new
!                end if
!              end do; end do
!              !
!              ! determine the new bounding box
!              bb_new%ic0 = min(bb_new%ic0,bb%ic0); bb_new%ic1 = max(bb_new%ic1,bb%ic1)
!              bb_new%ir0 = min(bb_new%ir0,bb%ir0); bb_new%ir1 = max(bb_new%ir1,bb%ir1)
!              bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1;  bb_new%nrow = bb_new%ir1 - bb_new%ir0 + 1
!              !
!              if (lfirst) then
!                hlev = 0; hlev(1:ih) = q%hlev(1:ih)
!                lfirst = .false.
!              end if
!            end if
!          end if
!        end do
!        if (lfound) then
!          bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb%ir1)*cs_gid
!          bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
!          bbx%cs = cs_gid
!          call bbo%set(child_bbi=bb_new, child_bbx=bbx)
!          !
!          q_new => xq%get_quad(lid_new)
!          call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
!          call q_new%set(hlev=hlev)
!          call q_new%get_bb(child_bbi=bb)
!          call get_mask(xid, q_new%gid, bb, mask)
!          call q_new%calc_prop(mask=mask)
!          xq%n = lid_new
!        end if
!      else
!        llast = .false.
!        if (ih == nh) then
!          llast = .true.
!        else
!          llast = .true.
!          do lid = 1, xq%n
!            q => xq%get_quad(lid)
!            if (q%get_flag(active=LDUM)) then
!              if (q%get_hlev(ih) == id) then
!                if (q%get_hlev(ih+1) /= 0) then
!                  llast = .false.
!                end if
!              end if
!            end if
!          end do
!        end if
!        !
!        if (llast) then
!          ! set 0: count weight for seperate IDS subject to full METIS partitioning
!          wgt_metis = R8ZERO
!          if (lgidsep) then
!            do j = 1, size(gid_sep)
!              gid = gid_sep(j); lid = g2lid(gid)
!              q => xq%get_quad(lid)
!              if (q%get_flag(active=LDUM)) then
!                if (q%get_hlev(ih) == id) then
!                  call q%get_prop(weight=wgt_q)
!                  if (wgt_q > wgt_tgt) then
!                    wgt_metis = wgt_metis + wgt_q
!                    call q%set_flag(active=.false.)
!                  end if
!                end if
!              end if
!            end do
!          end if
!          !
!          ! step 1: apply full METIS for quads with weight > weight_fac*weight_tgt
!          n = xq%n; lid_new = n
!          do lid = 1, n
!            q => xq%get_quad(lid)
!            if (q%get_flag(active=LDUM)) then
!              if (q%get_hlev(ih) == id) then
!                call q%get_prop(weight=wgt_q)
!                if (wgt_q > weight_fac*wgt_tgt) then
!                  ! number of parts
!                  np_full = nint(wgt_q/(weight_fac*wgt_tgt))
!                  np_full = max(2,np_full)
!                  call q%get_bb(child_bbi=bb)
!                  call get_mask(xid, q%gid, bb, mask)
!                  if (allocated(qweight)) deallocate(qweight)
!                  allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
!                  do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                    jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!                    if (mask(jc,jr) > 0) then
!                      qweight(jc,jr) = weight(ic,ir)
!                    end if
!                  end do; end do;
!                  !
!                  ! full grid METIS
!                  allocate(met)
!                  call met%init(qweight, np_full)
!                  call met%set_opts()
!                  call met%recur(verbose=.true.)
!                  if (allocated(part)) deallocate(part)
!                  allocate(part,source=mask)
!                  call met%set_ids(part)
!                  call met%clean(); deallocate(met); met => null()
!                  !
!                  do ip = 1, np_full
!                    lid_new = lid_new + 1; gid_new = gid_new + 1
!                    q_new => xq%get_quad(lid_new)
!                    call bb_new%init()
!                    do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                      jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!                      p = part(jc,jr)
!                      if (p == ip) then
!                       ! determine the new bounding box
!                        bb_new%ic0 = min(bb_new%ic0,ic); bb_new%ic1 = max(bb_new%ic1,ic)
!                        bb_new%ir0 = min(bb_new%ir0,ir); bb_new%ir1 = max(bb_new%ir1,ir)
!                        bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1;  bb_new%nrow = bb_new%ir1 - bb%ir0 + 1
!                        if (xid(ic,ir) /= q%gid) then
!                          call errmsg('quad_sub_grid_part: program error.')
!                        else
!                          xid(ic,ir) = gid_new
!                        end if
!                      end if 
!                    end do; end do
!                    !
!                    bbx%xll = xll + (bb_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb_new%ir1)*cs_gid
!                    bbx%xur = bbx%xll + bb_new%ncol*cs_gid; bbx%yur = bbx%yll + bb_new%nrow*cs_gid
!                    bbx%cs = cs_gid
!                    call bbo%set(child_bbi=bb_new, child_bbx=bbx)
!                    call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
!                    call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid)
!                    call get_mask(xid, gid_new, bb_new, mask)
!                    if (allocated(qweight)) deallocate(qweight)
!                    allocate(qweight(bb_new%ncol,bb_new%nrow)); qweight = I4ZERO
!                    do ir = bb_new%ir0, bb_new%ir1; do ic = bb_new%ic0, bb_new%ic1
!                      jr = ir - bb_new%ir0 + 1; jc = ic - bb_new%ic0 + 1
!                      if (mask(jc,jr) > 0) then
!                        qweight(jc,jr) = weight(ic,ir)
!                      end if
!                    end do; end do;
!                    call q_new%calc_prop(mask=mask, weight=qweight)
!                    xq%n = lid_new
!                  end do !ip
!                  !
!                  call q%set_flag(active=.false.)
!                end if
!              end if
!            end if
!          end do
!          !
!          ! set: apply lumped METIS to the other
!          n_old = xq%n
!          wgt = wgth(i) - wgt_metis
!          if (wgt < R8ZERO) then 
!            call errmsg('quad_sub_grid_part: program error')
!          end if
!          np_lump = 0
!          if (wgt > wgt_tgt) then
!            np_lump = max(2,nint(wgt/wgt_tgt))
!            !
!            ! determine the bounding box and work array
!            call bb_new%init(); lfirst = .true.
!            do iact = 1, 2
!              do lid = 1, xq%n
!                q => xq%get_quad(lid)
!                if (q%get_flag(active=LDUM)) then
!                  if (q%get_hlev(ih) == id) then
!                    call q%get_bb(child_bbi=bb)
!                    if (iact == 1) then
!                      bb_new%ic0 = min(bb_new%ic0,bb%ic0); bb_new%ic1 = max(bb_new%ic1,bb%ic1)
!                      bb_new%ir0 = min(bb_new%ir0,bb%ir0); bb_new%ir1 = max(bb_new%ir1,bb%ir1)
!                      bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1
!                      bb_new%nrow = bb_new%ir1 - bb_new%ir0 + 1
!                    else
!                      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                        if (xid(ic,ir) == q%gid) then
!                          jr = ir - bb_new%ir0 + 1; jc = ic - bb_new%ic0 + 1
!                          i4w2d(jc,jr) = q%lid
!                          lweight(jc,jr) = weight(ic,ir)
!                        end if
!                      end do; end do
!                    end if
!                    if (lfirst) then
!                      hlev = 0; hlev(1:ih) = q%hlev(1:ih)
!                      lfirst = .false.
!                    end if
!                  end if
!                end if
!              end do
!              if (iact == 1) then
!                if (allocated(i4w2d)) deallocate(i4w2d)
!                if (allocated(lweight)) deallocate(lweight)
!                allocate(i4w2d(bb_new%ncol,bb_new%nrow))
!                allocate(lweight(bb_new%ncol,bb_new%nrow))
!                i4w2d = I4ZERO; lweight = I4ZERO
!              end if
!            end do
!            !
!            np_lump_loc = np_lump
!            do while(.true.)
!              if (associated(met)) then
!                call met%clean(); deallocate(met); met => null()
!              end if
!              allocate(met)
!              call logmsg('**** Lumped METIS for '//ta([np_lump_loc])//' parts *****')
!              call met%init_lump(i4w2d, np_lump_loc, verbose=.true., weight=lweight)
!              call logmsg('**** Number of graph vertices: '//ta([met%nvtxs]))
!              call met%set_opts()
!              call met%recur(ok=lok)
!              !call met%set_opts(contig_in=0)
!              !call met%kway(ok=lok)
!              if (.not.lok) then
!                call logmsg('**** Warning, lumped METIS failed for '//ta([np_lump_loc])//' parts *****')
!                np_lump_loc = np_lump_loc - 1
!              else
!                exit
!              end if
!            end do
!            if (allocated(part)) deallocate(part)
!            !
!            do ip = 1, np_lump_loc
!              lid_new = lid_new + 1; gid_new = gid_new + 1
!              q_new => xq%get_quad(lid_new)
!              call bb_new%init()
!              do j = 1, met%nvtxs
!                if (met%part(j) == ip - 1) then
!                  lid = met%idmapinv(j)
!                  q => xq%get_quad(lid)
!                  if ((q%get_hlev(ih) /= id).or.(.not.q%get_flag(active=LDUM))) then
!                    call errmsg('quad_sub_grid_part: program error')
!                  end if
!                  call q%set_flag(active=.false.)
!                  call q%get_bb(child_bbi=bb)
!                  bb_new%ic0 = min(bb_new%ic0,bb%ic0); bb_new%ic1 = max(bb_new%ic1,bb%ic1)
!                  bb_new%ir0 = min(bb_new%ir0,bb%ir0); bb_new%ir1 = max(bb_new%ir1,bb%ir1)
!                  bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1
!                  bb_new%nrow = bb_new%ir1 - bb_new%ir0 + 1
!                  !
!                  do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                    if (xid(ic,ir) == q%gid) then
!                      xid(ic,ir) = gid_new
!                    end if
!                  end do; end do
!                end if
!              end do
!              bbx%xll = xll + (bb_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb_new%ir1)*cs_gid
!              bbx%xur = bbx%xll + bb_new%ncol*cs_gid; bbx%yur = bbx%yll + bb_new%nrow*cs_gid
!              bbx%cs = cs_gid
!              call bbo%set(child_bbi=bb_new, child_bbx=bbx)
!              call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo) !TODO: set parent LIST
!              call q_new%set(hlev=hlev) 
!              call q_new%get_bb(child_bbi=bb)
!              call get_mask(xid, q_new%gid, bb, mask)
!              call q_new%calc_prop(mask=mask)
!              xq%n = lid_new
!            end do
!            call met%clean(); deallocate(met); met => null()
!          end if
!          !
!          ! check
!          do lid = 1, n_old
!            q => xq%get_quad(lid)
!            if (q%get_flag(active=LDUM)) then
!              call errmsg('quad_sub_grid_part: program error.')
!            end if
!          end do
!          !
!          n = xq%get_number_active()
!          call logmsg('# unique IDs BEFORE splitting: '//ta([n]))
!          !
!          ! split in to non-unique
!          n = xq%n; lid_new = n
!          do lid = 1, n
!            q => xq%get_quad(lid)
!            if (q%get_flag(active=LDUM)) then
!              call q%get_bb(child_bbi=bb)
!              call get_mask(xid, q%gid, bb, i4w2d)
!              call calc_unique(i4w2d, 5, regun, regbb, nreg, idum, 0., 0., 0.)
!              if (nreg > 1) then
!                do ireg = 1, nreg
!                  lid_new = lid_new + 1; gid_new = gid_new + 1
!                  q_new => xq%get_quad(lid_new)
!                  call bb_new%init()
!                  do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                    jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!                    if (regun(jc,jr) == ireg) then
!                      bb_new%ic0 = min(bb_new%ic0,bb%ic0); bb_new%ic1 = max(bb_new%ic1,bb%ic1)
!                      bb_new%ir0 = min(bb_new%ir0,bb%ir0); bb_new%ir1 = max(bb_new%ir1,bb%ir1)
!                      bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1
!                      bb_new%nrow = bb_new%ir1 - bb_new%ir0 + 1
!                      xid(ic,ir) = gid_new
!                    end if
!                  end do; end do
!                  bbx%xll = xll + (bb_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb_new%ir1)*cs_gid
!                  bbx%xur = bbx%xll + bb_new%ncol*cs_gid; bbx%yur = bbx%yll + bb_new%nrow*cs_gid
!                  bbx%cs = cs_gid
!                  call bbo%set(child_bbi=bb_new, child_bbx=bbx)
!                  call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
!                  call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid) 
!                  call q_new%get_bb(child_bbi=bb)
!                  call get_mask(xid, q_new%gid, bb, mask)
!                  call q_new%calc_prop(mask=mask)
!                  xq%n = lid_new
!                end do
!                call q%set_flag(active=.false.)
!              end if
!            end if
!          end do
!          !
!          n = xq%get_number_active()
!          call logmsg('# unique IDs AFTER splitting: '//ta([n]))
!          !
!          if (lgidsep) then
!            n = xq%get_number_active()
!            np_rem = npart - n
!            wgt_tgt = wgt_metis / np_rem
!            !
!            n = xq%n; lid_new = n
!            do j = 1, size(gid_sep)
!              gid = gid_sep(j); lid = g2lid(gid)
!              q => xq%get_quad(lid)
!              call q%set_flag(active=.true.)
!              if (q%get_flag(active=LDUM)) then
!                if (q%get_hlev(ih) == id) then
!                  call q%get_prop(weight=wgt_q)
!                  if (wgt_q > wgt_tgt) then
!                    ! number of parts
!                    np_full = max(2,nint(wgt_q/wgt_tgt))
!                    np_full = min(np_rem, np_full)
!                    np_rem = np_rem - np_full
!                    !
!                    call q%get_bb(child_bbi=bb)
!                    call get_mask(xid, q%gid, bb, mask)
!                    if (allocated(qweight)) deallocate(qweight)
!                    allocate(qweight(bb%ncol,bb%nrow)); qweight = I4ZERO
!                    do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                      jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!                      if (mask(jc,jr) > 0) then
!                        qweight(jc,jr) = weight(ic,ir)
!                      end if
!                    end do; end do;
!                    !
!                    ! full grid METIS
!                    allocate(met)
!                    call met%init(qweight, np_full)
!                    call met%set_opts()
!                    call logmsg('***** Full METIS partitioning for GID '//ta([gid])// &
!                      ' into '//ta([np_full])//' parts *****')
!                    call met%recur()
!                    if (allocated(part)) deallocate(part)
!                    allocate(part,source=mask)
!                    call met%set_ids(part)
!                    call met%clean(); deallocate(met); met => null()
!                    !
!                    do ip = 1, np_full
!                      lid_new = lid_new + 1; gid_new = gid_new + 1
!                      q_new => xq%get_quad(lid_new)
!                      call bb_new%init()
!                      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!                        p = part(jc,jr)
!                        if (p == ip) then
!                         ! determine the new bounding box
!                          bb_new%ic0 = min(bb_new%ic0,ic); bb_new%ic1 = max(bb_new%ic1,ic)
!                          bb_new%ir0 = min(bb_new%ir0,ir); bb_new%ir1 = max(bb_new%ir1,ir)
!                          bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1;  bb_new%nrow = bb_new%ir1 - bb%ir0 + 1
!                          if (xid(ic,ir) /= q%gid) then
!                            call errmsg('quad_sub_grid_part: program error.')
!                          else
!                            xid(ic,ir) = gid_new
!                          end if
!                        end if 
!                      end do; end do
!                      !
!                      bbx%xll = xll + (bb_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb_new%ir1)*cs_gid
!                      bbx%xur = bbx%xll + bb_new%ncol*cs_gid; bbx%yur = bbx%yll + bb_new%nrow*cs_gid
!                      bbx%cs = cs_gid
!                      call bbo%set(child_bbi=bb_new, child_bbx=bbx)
!                      call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
!                      call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid)
!                      call get_mask(xid, gid_new, bb_new, mask)
!                      if (allocated(qweight)) deallocate(qweight)
!                      allocate(qweight(bb_new%ncol,bb_new%nrow)); qweight = I4ZERO
!                      do ir = bb_new%ir0, bb_new%ir1; do ic = bb_new%ic0, bb_new%ic1
!                        jr = ir - bb_new%ir0 + 1; jc = ic - bb_new%ic0 + 1
!                        if (mask(jc,jr) > 0) then
!                          qweight(jc,jr) = weight(ic,ir)
!                        end if
!                      end do; end do;
!                      call q_new%calc_prop(mask=mask, weight=qweight)
!                      call q_new%get_prop(weight=wgt_q)
!                      wgt_metis = wgt_metis + wgt_q
!                      xq%n = lid_new
!                    end do !ip
!                    !
!                    call q%set_flag(active=.false.)
!                  end if
!                end if
!              end if
!            end do
!          end if
!        end if ! ih = nh
!      end if
!    end do
!    !
!    deallocate(hid, wgth)
!  end do
!  !
!  !
!  ! renumber
!  allocate(gids_new(xq%n)); gids_new = 0
!  n = 0
!  !
!  do j = 1, size(gid_sep)
!     gid = gid_sep(j)
!     do lid = 1, xq%n
!       q => xq%get_quad(lid)
!       if (q%get_flag(active=LDUM)) then
!         if (q%gid_prent == gid) then
!           n = n + 1
!           gids_new(q%gid) = n
!         end if
!       end if
!     end do
!  end do
!  !
!  do lid = 1, xq%n
!    q => xq%get_quad(lid)
!    if (q%get_flag(active=LDUM)) then
!      if (gids_new(q%gid) == 0) then
!        n = n + 1
!        gids_new(q%gid) = n
!      end if
!    end if
!  end do
!  !
!  ! replace the global id
!  do lid = 1, xq%n
!    q => xq%get_quad(lid)
!    if (q%get_flag(active=LDUM)) then
!      gid = gids_new(q%gid)
!      call q%get_bb(child_bbi=bb)
!      call get_mask(xid, q%gid, bb, mask)
!      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!        if (mask(jc,jr) == 1) then
!          xid(ic,ir) = -gid
!        end if
!      end do; end do
!    end if
!  end do
!  xid = abs(xid)
!  !
!  ! write the new ids
!  call grid_load_imbalance(xid, xid_mv, weight, imbal, np, ximbal)
!  call logmsg('Overall load imbalance for '//ta([np])//' parts: '//ta([imbal]))
!  !
!  call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
!  call hdrg_gid%write(f_gid_out)
!  if (lwriteximbal) then
!    f = trim(f_gid_out)//'_imbal'
!    call writeflt(f, ximbal, nc, nr, xll, yll, cs_gid, R4ZERO)
!  end if
!  !
!  ! write the hierarchical grids
!  if (write_hiera) then
!    do ih = 1, nh
!      ! get the IDs
!      call get_unique(hmap(ih,:), hid)
!      nhid = size(hid)
!      !
!      if (allocated(i4w2d)) deallocate(i4w2d)
!      allocate(i4w2d(nc,nr)); i4w2d = 0
!      !
!      do i = 1, nhid
!        id = hid(i)
!        if (id == 0) cycle
!        do lid = 1, xq%n
!          q => xq%get_quad(lid)
!          if (q%get_flag(active=LDUM)) then
!            if (q%get_hlev(ih) == id) then
!              call q%get_bb(child_bbi=bb)
!              call get_mask(xid, q%gid, bb, mask)
!              do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
!                jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
!                if (mask(jc,jr) > 0) then
!                  i4w2d(ic,ir) = id
!                end if
!              end do; end do
!            end if
!          end if
!        end do
!      end do
!      !
!      fp = trim(sa(ih))//'_calc'
!      call writeflt(fp, i4w2d, nc, nr, xll, yll, cs_gid, 0)
!    end do
!  end if
!  
!  ! write the csv
!  !allocate(csv)
!  !call csv%init(file=f_out_csv, &
!  !    hdr_keys=['lid', 'gid', 'lid_prent', 'gid_prent', &
!  !      'xm', 'ym', 'area'],&
!  !    nr=xq%n, hdr_i_type=[i_i4, i_i4, i_i4, i_i4, &
!  !      i_r8, i_r8, i_r8,])
!  
!  ! clean-up
!  if (allocated(hmap)) deallocate(hmap)
!  if (allocated(hlev)) deallocate(hlev)
!  if (allocated(xid)) deallocate(xid)
!  if (allocated(i4w2d)) deallocate(i4w2d)
!  if (associated(hdrg_gid)) then
!    call hdrg_gid%clean(); deallocate(hdrg_gid); hdrg_gid => null()
!  end if
!  !
!  return
!end subroutine quad_sub_grid_part

subroutine quad_csv_add_field()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- local
  type(tCSV), pointer :: csv => null()
! -- local
  logical :: lfill
  character(len=MXSLEN) :: s
  character(len=MXSLEN), dimension(:), allocatable :: hdr, hdr_add, val_add
  character(len=MXSLEN), dimension(:), allocatable :: cfill
  integer(I4B) :: nc, nc_add, nc_max, ir, ic, jc, jc0
! ------------------------------------------------------------------------------
  if (len_trim(csv_fill_field) > 0) then
    lfill = .true.
  else
    lfill = .false.
  end if
  !
  allocate(csv)
  call csv%read_hdr(f_in_csv, hdr)
  nc = size(hdr)
  !
  call parse_line(csv_field, hdr_add, token_in=',')
  nc_add = size(hdr_add)
  nc_max = nc + nc_add
  !
  call csv%read(f_in_csv, nc_max=nc_max)
  !
  ! get the fill columns
  if (lfill) then
    call csv%get_column(key=csv_fill_field, ca=cfill)
  end if
  !
  ! add the new data
  call csv%add_hdr(hdr_add)
  !
  if (n_csv_val > 1) then
    jc0 = 2
  else
    jc0 = 0
  end if
  !
  do i = 1, n_csv_val
    call parse_line(csv_val(i), val_add, token_in=',')
    if (nc_add /= (size(val_add)-jc0)) then
      call errmsg('quad_csv_add_field: size csv_field differs from csv_val.')
    end if
    if (n_csv_val > 1) then
      read(val_add(1),*) ir0
      read(val_add(2),*) ir1
      if ((ir0 < 1).or.(ir0 > csv%nr).or. &
          (ir1 < 1).or.(ir1 > csv%nr).or.(ir0 > ir1)) then
        call errmsg('quad_csv_add_field: invalid row, '//trim(csv_val(i)))
      end if
    else
      ir0 = 1
      ir1 = csv%nr
    end if
    !
    do ir = ir0, ir1
      do jc = 1, nc_add
        ic = nc + jc
        if (lfill) then
          s = replace_token(val_add(jc+jc0),'$',cfill(ir))
        else
          s = val_add(jc+jc0)
        end if
        call csv%set_val(ic=ic, ir=ir, cv=trim(s))
      end do
    end do
    !
  end do
  !
  csv%file = trim(f_out_csv); call csv%write()
  !
  ! clean up
  call csv%clean(); deallocate(csv); csv => null()
  if (allocated(hdr)) deallocate(hdr)
  if (allocated(hdr_add)) deallocate(hdr_add)
  if (allocated(val_add)) deallocate(val_add)
  if (allocated(cfill)) deallocate(cfill)
  !
  return
end subroutine quad_csv_add_field

subroutine quad_idf_to_flt()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- local
  character(len=MXSLEN), dimension(:), allocatable :: f_idf_list
  character(len=MXSLEN) :: f_idf, f_flt
  integer(I4B) :: n, nc, nr, ic, ir, i
  real(R4B), dimension(:,:), allocatable :: x
  real(R4B) :: r4mv_read, r4xll, r4yll, r4cs
! ------------------------------------------------------------------------------
  !
  f_idf_list = get_dir_files(f_in_idf)
  n = size(f_idf_list)
  !
  call logmsg('Converting '//ta([n])//' IDF files...')
  do i = 1, n
    call logmsg('***** '//ta([i])//'/'//ta([n])//' *****')
    f_idf = f_idf_list(i)
    f_flt = strip_ext(f_idf)
    call readidf( f_idf, x, nc, nr, r4xll, r4yll, r4cs, r4mv)
    call writeflt(f_flt, x, nc, nr, r4xll, r4yll, r4cs, r4mv)
  end do
  !
  ! clean up
  deallocate(f_idf_list, x)
  !
  return
end subroutine quad_idf_to_flt

subroutine quad_merge_csv()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- local
  type(tCSV), pointer :: csv => null(), csv_tot => null()
  character(len=MXSLEN), dimension(:), allocatable :: sa
  character(len=MXSLEN) :: f
  integer(I4B), dimension(:), allocatable :: nr_max_arr
  integer(I4B) :: n_csv, nr_max, i, ir_offset
! ------------------------------------------------------------------------------
  !
  call parse_line(f_in_csv_post, sa, token_in=',')
  n_csv = size(sa)
  allocate(nr_max_arr(n_csv)); nr_max_arr = 0
  !
  ! first, count the total number
  nr_max = 0
  do i = 1, n_csv
    f = trim(f_in_csv_pref)//trim(sa(i))//'.csv'
    allocate(csv); call csv%read(f)
    nr_max_arr(i) = csv%get_nr()
    call csv%clean(); deallocate(csv); csv => null()
  end do
  nr_max = sum(nr_max_arr)
  !
  ! fill the merged csv
  allocate(csv_tot)
  f = trim(f_in_csv_pref)//trim(sa(1))//'.csv'
  call csv_tot%read(file=f, nr_max=nr_max)
  !
  ir_offset = nr_max_arr(1)
  do i = 2, n_csv
    f = trim(f_in_csv_pref)//trim(sa(i))//'.csv'
    call csv_tot%read(file=f, ir_offset=ir_offset)
    ir_offset = ir_offset + nr_max_arr(i)
  end do
  !
  ! write the merged csv-file
  csv_tot%nr = nr_max; csv_tot%file = trim(f_out_csv); call csv_tot%write()
  !
  ! clean up
  call csv_tot%clean(); deallocate(csv_tot); csv_tot => null()
  if (allocated(nr_max_arr)) deallocate(nr_max_arr)
  !
  return
end subroutine quad_merge_csv

subroutine get_mask(xid, id, bb, mask)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- arguments
  integer(I4B), dimension(:,:), intent(in) :: xid
  integer(I4B), intent(in) :: id
  type(tBB), intent(in) :: bb
  integer(I4B), dimension(:,:), allocatable, intent(inout) :: mask
! -- local
  integer(I4B) :: jr, jc, n
! ------------------------------------------------------------------------------
  !
  if (allocated(mask)) deallocate(mask)
  allocate(mask(bb%ncol,bb%nrow)); mask = 0
  n = 0
  do ir = bb%ir0, bb%ir1
    do ic = bb%ic0, bb%ic1
      if (xid(ic,ir) == id) then
        n = n + 1
        jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
        mask(jc,jr) = 1
      end if
    end do
  end do
  !
  if (n == 0) then
    call errmsg('get_mask: no id found.')
  end if
  !
  return
end subroutine get_mask

subroutine quad_fill_gap()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  !
  ! -- local
  character(len=MXSLEN) :: fi, fo
  character(len=MXSLEN), dimension(:), allocatable :: fi_list, fo_list
  integer(I4B) :: n, i
! ------------------------------------------------------------------------------

  allocate(hdrg_gid, hdrg_mask)
  !
  call hdrg_mask%read_full_grid(f_gid_mask)
  call hdrg_mask%get_grid(xi4=mask, mvi4=mask_mv)
  call hdrg_mask%clean(); deallocate(hdrg_mask); hdrg_mask => null()
  
  call parse_line(f_gid_in,  fi_list, token_in=',')
  n = size(fi_list)
  call parse_line(f_gid_out, fo_list, token_in=',')
  if (n /= size(fo_list)) then
    call errmsg('Inconsistent input/output.')
  end if
  !
  allocate(r4w2d(nc,nr))
  !
  do i = 1, n
    fi = fi_list(i); fo = fo_list(i)
    
    call hdrg_gid%read_full_grid(fi)
    call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
    !
    nc = size(xid,1); nr = size(xid,2)
    if ((nc /= size(mask,1)).or.(nr /= size(mask,2))) then
      call errmsg('quad_filter: invalid dimensions')
    end if
    !
    do ir = 1, nr; do ic = 1, nc
      if (mask(ic,ir) == gid_mask) then
        if (xid(ic,ir) == xid_mv) then
          xid(ic,ir) = -1
        end if
      else
        xid(ic,ir) = xid_mv
      end if
    end do; end do
    !call fill_with_nearest(xid, xid_mv, -1)
    r4w2d = real(xid,R4B)
    call fillgap(r4w2d, real(xid_mv,R4B), -1., nine_point=.false.)
    xid = int(r4w2d,I4B)
    !
    call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
    call hdrg_gid%write(fo)
    !
    call hdrg_gid%clean()
  end do
  !
  ! clean-up
  if (allocated(xid)) deallocate(xid)
  if (allocated(r4w2d)) deallocate(r4w2d)
  if (allocated(fi_list)) deallocate(fi_list)
  if (allocated(fo_list)) deallocate(fo_list)
  if (associated(hdrg_mask)) then
    call hdrg_mask%clean(); deallocate(hdrg_mask); hdrg_mask => null()
  end if
  if (associated(hdrg_gid)) then
    call hdrg_gid%clean(); deallocate(hdrg_gid); hdrg_gid => null()
  end if
  !
  return
end subroutine quad_fill_gap


subroutine quad_init()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  !
  allocate(vrt); call vrt%init(f_vrt);
  ntile = vrt%get_ntiles(); call logmsg('# tiles: '//ta([ntile]))
  !
  it0 = 1; it1 = ntile
  !it0 = 1; it1 = 2
  
  call vrt%get_bb(dst_bbi=bbi, dst_bbx=bbx)
  gnc = bbi%ncol; gnr = bbi%nrow
  call logmsg('Size of '//trim(f_vrt)//': '//ta([gnc])//' cols x '// &
    ta([gnr])//' rows')
  cs_gid = bbx%cs; xll = bbx%xll; yll = bbx%yll; yur =  bbx%yur
  !
  nlev = get_number_of_levels(cs_gid, cs_max)
  refr8 = get_refinement_level(nlev)
  refi4 = max(0, int(refr8, I4B))
  !
  ! loop over the blocks and determine the bounding boxes
  allocate(bb_gid(gid_max), gids(gid_max), g2lid(gid_max))
  gids = 0; nid = 0
  do it = it0, it1
    if (lrenumber) then
      f_tile = vrt%get_tile_file(itile=it)
      call vrt%read_full_tile(itile=it, xi4=xid, mvi4=mvxid, &
        renum=.true., nid=nid, f_csv=trim(f_tile)//'.renum.csv')
    else
      call vrt%read_full_tile(itile=it, xi4=xid, mvi4=mvxid)
    end if
    !
    call vrt%get_bb(itile=it, dst_bbi=bbi)
    do ir = 1, size(xid,2)
      do ic = 1, size(xid,1)
        gid = xid(ic,ir)
        if (gid /= mvxid) then
          if (gid > gid_max) then
            call errmsg('Increase gid_max > '//ta([gid]))
          end if
          gids(gid) = 1
          gir = bbi%ir0 + ir - 1; gic = bbi%ic0 + ic - 1
          bb => bb_gid(gid)
          bb%ic0 = min(bb%ic0, gic); bb%ic1 = max(bb%ic1, gic)
          bb%ir0 = min(bb%ir0, gir); bb%ir1 = max(bb%ir1, gir)
        end if
      end do
    end do
  end do
  !
  ! resize the bounding boxes when necessary
  nlid = 0
  do gid = 1, gid_max
    if (gids(gid) == 1) then
      nlid = nlid + 1
      !
      bb => bb_gid(gid)
      !
      ! make bounding box slightly larger
      bb%ic0 = bb%ic0 - 1; bb%ic0 = max(bb%ic0,1)
      bb%ic1 = bb%ic1 + 1; bb%ic1 = min(bb%ic1,gnc)
      bb%ir0 = bb%ir0 - 1; bb%ir0 = max(bb%ir0,1)
      bb%ir1 = bb%ir1 + 1; bb%ir1 = min(bb%ir1,gnr)
      !
      ! match grids
      if (nlev < 0) then
        call get_xy(x0, y0, bb%ic0, bb%ir0, xll, yur, cs_gid)
        call get_xy(x1, y1, bb%ic1, bb%ir1, xll, yur, cs_gid)
        !
        call get_icr(ic0, ir0, x0, y0, xll, yur, cs_max)
        call get_icr(ic1, ir1, x1, y1, xll, yur, cs_max)
        !
        bb%ic0 = (ic0-1)*refi4 + 1; bb%ic1 = ic1*refi4
        bb%ir0 = (ir0-1)*refi4 + 1; bb%ir1 = ir1*refi4
      end if
      !
      bb%ncol = bb%ic1 - bb%ic0 + 1
      bb%nrow = bb%ir1 - bb%ir0 + 1
    end if
  end do
  call logmsg('# ids found: '//ta([nlid]))
  !
  ! first, determine the interior & exterior boundaries of each quad (ids!)
  allocate(xq)
  call xq%init(nlid, gid_max, dir=d_out)
  call xq%generate_uuid(luse_uuid, uuid_out)
  !
  allocate(l2gid(gid_max))
  lid = 0; g2lid = 0
  bbi%ic0 = 1; bbi%ic1 = gnc; bbi%ir0 = 1; bbi%ir1 = gnr; bbi%ncol = gnc; bbi%nrow = gnr
  bbx%xll = xll; bbx%yll = yll; bbx%cs = cs_gid
  bbx%xur = bbx%xll + bbi%ncol*bbx%cs; bbx%yur = bbx%yll + bbi%nrow*bbx%cs
  call bbo%set(prent_bbi=bbi, prent_bbx=bbx)
  do gid = 1, gid_max
    if (gids(gid) == 1) then
      lid = lid + 1
      q => xq%get_quad(lid)
      bb => bb_gid(gid) ! local bounding box
      bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (gnr-bb%ir1)*cs_gid
      bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
      call bbo%set(child_bbi=bb, child_bbx=bbx)
      call q%init(gid=gid, lid=lid, bbo=bbo, ntile=ntile, cs_max=cs_max)
      l2gid(lid) = gid
      g2lid(gid) = lid
    end if
  end do
  !
  ! set the itile flag
  do it = 1, ntile
    call vrt%get_bb(itile=it, dst_bbi=bbi)
    do lid = 1, nlid
      q => xq%get_quad(lid)
      call q%set(itile=it, tile_bbi=bbi)
    end do
  end do
  !
  ! map the block data to the local global ids
  do lid = 1, nlid
    q => xq%get_quad(lid)
    call q%init_map_gids(mv=I4ZERO, mv_map=I4MINONE)
  end do
  !
  ! determine the interfaces and properties
  nid = 0
  do it = it0, it1
    if (lrenumber) then
      call logmsg('# IDs begin: '//ta((/nid/)))
      call vrt%read_full_tile(itile=it, hdrg=hdrg, renum=.true., nid=nid, &
        clean_hdrg=.false.)
      call logmsg('# IDs end  : '//ta((/nid/)))
    else
      call vrt%read_full_tile(itile=it, hdrg=hdrg, clean_hdrg=.false.)
    end if
    !
    do lid = 1, nlid
      !
      q => xq%get_quad(lid)
      !
      if (.not.q%get_flag(mapped_gids=LDUM)) then
        call q%map_gids(hdrg, it)
      end if
      if (      q%get_flag(mapped_gids=LDUM).and.&
          (.not.q%get_flag(  work_done=LDUM))) then
        !
        call vrt%get_bb(itile=it, dst_bbx=bbx_t)
        !
        !if (q%gid /= 600240) cycle
        !
        call q%calc_interface(g2lid)
        if (.not.q%check_gid_mask()) then
          call logmsg('Tried to repair inconsistent mask....')
        end if
        if (.not.q%check_gid_mask()) then
          call errmsg('Inconsistent mask.')
        end if
        call q%calc_prop(itile=it, tile_bbx=bbx_t)
        call q%clean_gids()
        call q%set_flag(work_done=.true.)
      end if
    end do
    call hdrg%clean()
  end do
  !
  ! connect the interface cell numbers 1-on-1
  call xq%reorder_intf()
  !
  ! deactivate quads that are not read
  n_inact = 0
  do lid = 1, nlid
    q => xq%get_quad(lid)
    if (.not.q%get_flag(mapped_gids=LDUM)) then
      n_inact = n_inact + 1
      call q%set_flag(active=.false.)
    end if
  end do
  call logmsg('# inactive quads: '//ta([n_inact])//'...')
  !
  return
end subroutine quad_init
!
subroutine quad_graph()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------

  ! construct the graph
  call xq%construct_graph()
  call xq%disconnect_graph()
  !
  return
end subroutine quad_graph

subroutine quad_split()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  call xq%split(area_max, g2lid)
  !
  return
end subroutine quad_split

subroutine quad_join()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  call xq%join(area_min)
  !
  return
end subroutine quad_join

subroutine quad_write()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  call xq%write_bb()
  call xq%write_intf('.intf.bin')
  call xq%write_graphs()
  if (len_trim(f_out_csv) == 0) then
    call xq%write_props()
  else
    call xq%write_props(f_out_csv)
  end if
  !
  return
end subroutine quad_write
!
subroutine quad_clean()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  if (associated(xq)) then
    call xq%clean(); deallocate(xq); xq => null()
  end if
  !
  return
end subroutine quad_clean

subroutine quad_read(fp_intf, read_vintf)
! ******************************************************************************
!
!    SPECIFICATIONS:
  ! -- local
  character(len=*), intent(in), optional :: fp_intf
  logical, intent(in), optional :: read_vintf
! ------------------------------------------------------------------------------
  allocate(xq)
  call xq%init_select_bb(gid_max, f_in_csv, uuid_in, d_in, fields)
  if (present(fp_intf)) then
    call xq%init_select_intf(fp_intf, read_vintf)
    call xq%init_select_graph()
  end if
  !call xq%init_select(gid_max, f_in_csv, uuid_in, d_in, fields, fp_intf, &
  !  read_vintf)
  !
  return
end subroutine quad_read

subroutine quad_balancing()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  call xq%balance_graphs(lwrite_props, f_out_csv)
  !
  return
end subroutine quad_balancing

subroutine quad_init_models(set_dat_mod)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- dummy
  logical, intent(in), optional :: set_dat_mod
  ! -- local
  character(len=1), dimension(2), parameter :: comment = [';','#']
  !
  type(tLayerModel),  pointer :: lay_mod  => null()
  type(tLayerModels), pointer :: lay_mods => null()
  type(tDataModel),   pointer :: dat_mod  => null()
  type(tDataModels),  pointer :: dat_mods => null()
  logical :: lfound, lread, lcompress, set_dat_mod_loc, ldone_zp
  character(len=MXSLEN) :: f, s, s1, slc, ext
  character(len=MXSLEN), dimension(:), allocatable :: keys, fz, sa, fdat
  integer(I4B), dimension(:), allocatable :: nlay
  integer(I4B) :: ios, n_inp_mod, n_dat, ilm, ju
! ------------------------------------------------------------------------------
  !
  if (present(set_dat_mod)) then
    set_dat_mod_loc = set_dat_mod
  else
    set_dat_mod_loc = .true.
  end if
  !
  call open_file(f_mod_def_inp, iu, 'r')
  call read_line(iu, f, ios, comment)
  !
  allocate(xq%lay_mods)
  lay_mods => xq%lay_mods
  call lay_mods%init(f)
  call lay_mods%get(lookup_keys=keys, nlay=nlay, n_inp_mod=n_inp_mod)
  !
  if (set_dat_mod_loc) then
    allocate(xq%dat_mods)
    dat_mods => xq%dat_mods
    call dat_mods%init(n_inp_mod)
  end if
  !
  allocate(fdat(2))
  do i = 1, n_inp_mod
    if (allocated(fz)) deallocate(fz)
    allocate(fz(nlay(i)+1))
    !
    ! read the model ID
    call read_line(iu, s1, ios, comment)
    if (ios /= 0) call errmsg('Could not read '//trim(f_mod_def_inp)//'.')
    !
    slc = change_case(adjustl(s1),'l')
    if (trim(slc) == trim(keys(i))) then
      lfound = .true.
    else 
      call errmsg('quad_init_layer_models: could not find data'//&
        ' for layer model '//trim(keys(i))//'.')
    end if
    !
    ! read the zp definition file
    call read_line(iu, f, ios, comment)
    if (ios /= 0) call errmsg('Could not read '//trim(f_mod_def_inp)//'.')
    ext = get_ext(f)
    select case(ext)
    case('.txt')
      lcompress = .false.
      call open_file(f, ju, 'r')
      n = 0
      do while(.true.)
        call read_line(ju, f, ios, comment)
        if (ios /= 0) exit
        if (len_trim(f) > 0) then
          n = n + 1
          fz(n) = f
        end if
      end do
      close(ju)
      !
      if (n /= (nlay(i)+1)) then
        call errmsg('Error reading: '//trim(f)//'.')
      end if
    case('.vrt')
      lcompress = .true.
      fz(1) = f
    case default
      call errmsg('Unrecognized filetype: '//trim(f))
    end select
    !
    call read_line(iu, fdat(1), ios, comment)
    if (ios /= 0) call errmsg('Could not read '//trim(f_mod_def_inp)//'.')
    call read_line(iu, fdat(2), ios, comment)
    if (ios /= 0) call errmsg('Could not read '//trim(f_mod_def_inp)//'.')
    
    ! set the data 
    lay_mod => lay_mods%lay_mods(i)
    call lay_mod%init(keys(i), nlay(i), fz, lcompress)
    !
    if (set_dat_mod_loc) then
      dat_mod => dat_mods%dat_mods(i)
      call dat_mod%init(keys(i), fdat(1), fdat(2), n_inp_mod)
    end if
  end do
  !
  close(iu)
  !
  ! determine the unique mappings
  if (set_dat_mod_loc) then
    call xq%dat_mods%create_uni_map()
  end if
  !
  !allocate(fdat(2))
  !do i = 1, n_inp_mod
  !  rewind(iu); lread = .false.; lfound = .false.; lcompress = .false.
  !  !
  !  if (allocated(fz)) deallocate(fz)
  !  allocate(fz(nlay(i)+1)); n = 0
  !  !
  !  do while(.true.)
  !    read(unit=iu,iostat=ios,fmt='(a)') s
  !    if (ios /= 0) exit
  !    if (len_trim(s) == 0) cycle
  !    !
  !    call split_str(s, ' ', sa)
  !    s1 = sa(1)
  !    !
  !    if (s1(1:1) == '#') cycle
  !    !
  !    if (lread) then
  !      if (n_dat > 0) then
  !        fdat(n_dat) = s
  !      else
  !        n = n + 1; fz(n) = s
  !      end if
  !      ldone_zp = .false.
  !      if (lcompress) ldone_zp = .true.
  !      if (n == nlay(i)+1) ldone_zp = .true.
  !      if (ldone_zp) then
  !        n_dat = n_dat + 1
  !      end if
  !      if (n_dat > 2) lread = .false.
  !    end if
  !    !
  !    slc = change_case(adjustl(s1),'l')
  !    if (trim(slc) == trim(keys(i))) then
  !      lread = .true.; lfound = .true.; lcompress = .false.
  !      n_dat = 0
  !      if (size(sa) >  1) then
  !        slc = change_case(adjustl(sa(2)),'l')
  !        if (trim(slc) == 'compressed') then
  !          lcompress = .true.
  !        end if
  !      end if
  !    end if
  !  end do
  !  !
  !  if (.not.lfound) then
  !    call errmsg('quad_init_layer_models: could not find data'//&
  !      ' for layer model '//trim(keys(i))//'.')
  !  end if
  !  !
  !  ! set the data 
  !  lay_mod => lay_mods%lay_mods(i)
  !  call lay_mod%init(keys(i), nlay(i), fz, lcompress)
  !  !
  !  if (set_dat_mod_loc) then
  !    dat_mod => dat_mods%dat_mods(i)
  !    call dat_mod%init(keys(i), fdat(1), fdat(2))
  !  end if
  !end do
  !!
  !close(iu)
  !
  do i = 1, xq%n
    q => xq%get_quad(i)
    if (q%get_flag(active=LDUM)) then
      q%lay_mods => lay_mods
    end if
  end do
  
  if (set_dat_mod_loc) then
    do i = 1, xq%n
      q => xq%get_quad(i)
      if (q%get_flag(active=LDUM)) then
        call q%get_prop_csv(ikey=i_lay_mod, i4v=ilm)
        q%dat_mod => dat_mods%dat_mods(ilm)
      end if
    end do
  end if
  !
  if (allocated(fz)) deallocate(fz)
  if (allocated(fdat)) deallocate(fdat)
  !call lay_mods%clean()
  !
  return
end subroutine quad_init_models

subroutine quad_layer_models_interface()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  if (lwrite_props) then
    call xq%add_lm_intf(f_out_csv)
  else
    call xq%add_lm_intf()
  end if
  !
  return
end subroutine quad_layer_models_interface

subroutine quad_mf6_xch_write()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  call xq%write_mf6_xch_intf(xch_root_dir, xch_id_field, f_out_csv)
  !
  return
end subroutine quad_mf6_xch_write

subroutine quad_mf6_xch_write_merge()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- local
  type(tQuad), pointer :: q_m1 => null(), q_m2 => null()
  type(tCsv), pointer :: csv => null()
  type(tIntf), pointer :: intf => null(), intf_m1 => null(), intf_m2 => null()
  type(tNbrIntf), pointer :: nintf => null(), nintf_m1 => null(), nintf_m2 => null()
  type(tMf6Wbd), pointer :: wbd => null()
  type(tMF6Exchange), pointer :: xch_m1 => null(), xch_m2 => null(), xch_im => null()
  !
  logical :: lfound, lskip, lwrite
  !
  character(len=1) :: slash
  character(len=MXSLEN) :: d, s, id, f
  character(len=MXSLENLONG) :: s_long
  !
  integer(I4B) :: lid, im, jm, km, nim, inbr, jnbr, lid_nbr, iact, nlid_nbr
  integer(I4B) :: ihc_min, ihc_max, n1, n2, iexg, imm, jmm
  integer(I4B), dimension(:), allocatable :: im_arr, lid2im_arr, lid_arr, lid_arr_inv
  integer(I4B), dimension(:), allocatable :: lid_arr_nbr, lid_arr_nbr_inv, nbr_im
  integer(I4B), dimension(:), allocatable :: nodes_arr, lid2im_lidx, lid2jm_lidx
  integer(I4B), dimension(:), allocatable :: im_nodes_offset, jm_nodes_offset
  integer(I4B), dimension(:), allocatable :: i4w1, i4w2
  !
  real(R8B) :: cl1_min, cl1_max, cl2_min, cl2_max, hwva_min, hwva_max
! ------------------------------------------------------------------------------
  !
  ! create the directory for output
  d = xch_root_dir
  call create_dir(d, .true.)
  !
  allocate(wbd)
  call wbd%init(f_out_csv)
  !
  ! get the nodes array for all quads
  call xq%props%csv%get_column(key='lid', i4a=i4w1)
  call xq%props%csv%get_column(key='nodes', i4a=i4w2)
  allocate(nodes_arr(xq%n)); nodes_arr = 0
  do i = 1, size(i4w1)
    lid = i4w1(i)
    nodes_arr(lid) = i4w2(i)
  end do
  deallocate(i4w1, i4w2)
  !
  allocate(csv)
  call csv%read(f_in_csv_merge)
  call csv%get_column(key='lid', i4a=im_arr)
  nim = size(im_arr)
  !
  allocate(lid2im_arr(gid_max)); lid2im_arr = 0
  !
  ! first, flag the quads and set im
  do im = 1, nim
    ! get the local lids
    call csv%get_val(ir=im, ic=csv%get_col('lid_merged'), cv=s_long)
    call parse_line(s=s_long, i4a=lid_arr, token_in=';'); nlid = size(lid_arr)
    do i = 1, nlid
      lid = lid_arr(i)
      lid2im_arr(lid) = im
    end do
  end do
  !
  allocate(nbr_im(nim))
  do im = 1, nim
    ! create the mappings
    call csv%get_val(ir=im, ic=csv%get_col('lid_merged'), cv=s_long)
    call parse_line(s=s_long, i4a=lid_arr, token_in=';'); nlid = size(lid_arr)
    allocate(im_nodes_offset(nlid), lid_arr_inv(xq%n)); lid_arr_inv = 0
    n = 0
    do i = 1, nlid
      lid = lid_arr(i)
      lid_arr_inv(lid) = i
      im_nodes_offset(i) = n
      n = n + nodes_arr(lid)
    end do
    !
    ! determine the neighbor
    nbr_im = 0
    do i = 1, nlid
      lid = lid_arr(i)
      q => xq%get_quad(lid); intf => q%intf
      do inbr = 1, intf%n_nbr ! loop over interfaces
        nintf => intf%nbr_intf(inbr)
        lid_nbr = nintf%nbr_lid
        jm = lid2im_arr(lid_nbr)
        if ((jm > 0).and.(jm > im)) then ! jm > im: symmetric part
          nbr_im(jm) = 1
        end if
      end do
    end do
    !
    ! loop over the neighbors
    slash = get_slash()
    do jm = 1, nim
      if (nbr_im(jm) == 1) then
        !
        ! create the mappings for the neighbor
        call csv%get_val(ir=jm, ic=csv%get_col('lid_merged'), cv=s_long)
        call parse_line(s=s_long, i4a=lid_arr_nbr, token_in=';'); nlid_nbr = size(lid_arr_nbr)
        allocate(jm_nodes_offset(nlid_nbr), lid_arr_nbr_inv(xq%n)); lid_arr_nbr_inv = 0
        n = 0
        do i = 1, nlid_nbr
          lid = lid_arr_nbr(i)
          lid_arr_nbr_inv(lid) = i
          jm_nodes_offset(i) = n
          n = n + nodes_arr(lid)
        end do
        !
        ihc_min  = huge(ihc_min);   ihc_max = -huge(ihc_max)
        cl1_min  = huge(cl1_min);   cl1_max = -huge(cl1_max)
        cl2_min  = huge(cl2_min);   cl2_max = -huge(cl2_max)
        hwva_min = huge(hwva_min); hwva_max = -huge(hwva_max)
        !
        ! create the node offsets
        allocate(xch_im)
        call xch_im%init()
        !
        do iact = 1, 2
          xch_im%nexg = 0
          do i = 1, nlid
            lid = lid_arr(i)
            q_m1 => xq%get_quad(lid); intf_m1 => q_m1%intf
            do inbr = 1, intf_m1%n_nbr ! loop over interfaces
              nintf_m1 => intf_m1%nbr_intf(inbr)
              lid_nbr = nintf_m1%nbr_lid
              km = lid2im_arr(lid_nbr)
              if (km == jm) then ! take the non-interior interface
                if (nintf_m1%xch_active == 1) then
                  xch_m1 => nintf_m1%xch
                  ihc_min = min(ihc_min,xch_m1%ihc); ihc_max = max(ihc_max,xch_m1%ihc)
                  cl1_min = min(cl1_min,xch_m1%cl1); cl1_max = max(cl1_max,xch_m1%cl1)
                  cl2_min = min(cl2_min,xch_m1%cl2); cl2_max = max(cl2_max,xch_m1%cl2)
                  hwva_min = min(hwva_min,xch_m1%hwva); hwva_max = max(hwva_max,xch_m1%hwva)
                  if (iact == 1) then
                    xch_im%nexg = xch_im%nexg + xch_m1%nexg
                  else
                    imm = lid_arr_inv(lid); jmm = lid_arr_nbr_inv(lid_nbr)
                    if ((imm == 0).or.(jmm == 0)) then
                      call errmsg('quad_mf6_xch_write_merge: program error.')
                    end if
                    do iexg = 1, xch_m1%nexg
                      n1 = xch_m1%cellidm1(iexg)
                      n2 = xch_m1%cellidm2(iexg)
                      !
                      xch_im%nexg = xch_im%nexg + 1
                      xch_im%cellidm1(xch_im%nexg) = n1 + im_nodes_offset(imm)
                      xch_im%cellidm2(xch_im%nexg) = n2 + jm_nodes_offset(jmm)
                    end do
                  end if
                else
                  ! find the quad containing the exchange data
                  lfound = .false.; lskip = .false.
                  q_m2 => xq%get_quad(lid_nbr); intf_m2 => q_m2%intf
                  do jnbr = 1, intf_m2%n_nbr
                    nintf_m2 => intf_m2%nbr_intf(jnbr)
                    if (nintf_m2%nbr_lid == lid) then
                      if (nintf_m2%xch_active == 0) then
                        lskip = .true.
                        if (iact == 2) then
                          call logmsg('Warning: no connection found for lid '// &
                            trim(ta([lid]))//' -> '//trim(ta([lid_nbr]))//'!')
                        end if
                      end if
                      xch_m2 => nintf_m2%xch
                      lfound = .true.
                      exit
                    end if
                  end do
                  !
                  if (.not.lfound) then
                    call errmsg('quad_mf6_xch_write_merge: program error 2.')
                  end if
                  if (.not.lskip) then
                    ihc_min = min(ihc_min,xch_m2%ihc); ihc_max = max(ihc_max,xch_m2%ihc)
                    cl1_min = min(cl1_min,xch_m2%cl1); cl1_max = max(cl1_max,xch_m2%cl1)
                    cl2_min = min(cl2_min,xch_m2%cl2); cl2_max = max(cl2_max,xch_m2%cl2)
                    hwva_min = min(hwva_min,xch_m2%hwva); hwva_max = max(hwva_max,xch_m2%hwva)
                    if (iact == 1) then
                      xch_im%nexg = xch_im%nexg + xch_m2%nexg
                    else
                      imm = lid_arr_inv(lid); jmm = lid_arr_nbr_inv(lid_nbr)
                      if ((imm == 0).or.(jmm == 0)) then
                        call errmsg('quad_mf6_xch_write_merge: program error.')
                      end if
                      do iexg = 1, xch_m2%nexg
                        n2 = xch_m2%cellidm1(iexg)
                        n1 = xch_m2%cellidm2(iexg)
                        !
                        xch_im%nexg = xch_im%nexg + 1
                        xch_im%cellidm1(xch_im%nexg) = n1 + im_nodes_offset(imm)
                        xch_im%cellidm2(xch_im%nexg) = n2 + jm_nodes_offset(jmm)
                      end do
                    end if
                  end if !lskip
                end if
              end if
            end do !inbr
          end do
          lwrite = .true.
          if (iact == 1) then
            if (xch_im%nexg > 0) then
              if ((ihc_min /= ihc_max).or.(cl1_min /= cl1_max).or.&
                  (cl2_min /= cl2_max).or.(hwva_min /= hwva_max)) then
                call errmsg('quad_mf6_xch_write_merge: program error 4.')
              else
                xch_im%ihc  = ihc_min
                xch_im%cl1  = cl1_min
                xch_im%cl2  = cl2_min
                xch_im%hwva = hwva_min
              end if
              allocate(xch_im%cellidm1(xch_im%nexg)); xch_im%cellidm1 = 0
              allocate(xch_im%cellidm2(xch_im%nexg)); xch_im%cellidm2 = 0
            else
              lwrite = .false.
            end if
          end if
          if (.not.lwrite) then
            call logmsg('WARNING: no exchanges found for part '// &
              trim(ta([im]))//' -> '//trim(ta([jm]))//'!')
            exit
          end if
        end do !iact
        !
        ! write the exchange
        if (lwrite) then
          id = trim(ta([im]))//'-'//trim(ta([jm]))
          f = trim(xch_root_dir)//slash//'exchangedata_'//trim(id)//'.asc'
          call xch_im%write(f, id, wbd)
        end if
        !
        ! clean up
        if (associated(xch_im)) then
          call xch_im%clean(); deallocate(xch_im); xch_im => null()
        end if
        if (allocated(jm_nodes_offset)) deallocate(jm_nodes_offset)
        if (allocated(lid_arr_nbr_inv)) deallocate(lid_arr_nbr_inv)
      end if
    end do ! jm-loop
    !
    ! clean up
    if (allocated(im_nodes_offset)) deallocate(im_nodes_offset)
    if (allocated(lid_arr_inv)) deallocate(lid_arr_inv)
  end do ! im-loop
  !
  ! write the csv-file
  call wbd%write_csv()
  !
  ! clean up
  call wbd%clean(); deallocate(wbd); wbd => null()
  call csv%clean(); deallocate(csv); csv => null()
  if (allocated(nodes_arr))  deallocate(nodes_arr)
  if (allocated(im_arr))     deallocate(im_arr)
  if (allocated(lid2im_arr)) deallocate(lid2im_arr)
  if (allocated(nbr_im))     deallocate(nbr_im)
  !
  return
end subroutine quad_mf6_xch_write_merge

subroutine quad_mf6_data_write_old()
! ******************************************************************************
!
!    SPECIFICATIONS:
  ! -- local
  character(len=MXSLEN) :: d, f_log
  logical :: writelog
  integer(I4B) :: lid0, lid1, n_act, iu
! ------------------------------------------------------------------------------
  ! set the ranges
  if (lid_min > 0) then
     lid0 = max(1,lid_min)
  else
    lid0 = 1
  end if
  if (lid_max > 0) then
     lid1 = min(xq%n,lid_max)
  else
    lid1 = xq%n
  end if
  if (lid0 > lid1) then
    call errmsg('quad_mf6_data_write: lid_min > lid_max.')
  end if
  !
  ! determine total number of active quads
  n_act = 0
  do lid = lid0, lid1
    if (q%get_flag(active=LDUM)) n_act = n_act + 1
  end do
  !
  if (len_trim(d_log) > 0) then
    writelog = .true.
    d = trim(d_log)
    call create_dir(d, .true.)
  else
    writelog = .false.
  end if
  !
  n = 0
  do lid = lid0, lid1
    q => xq%get_quad(lid)
    !
    ! DEBUG:
    !if (q%gid /= 637452) cycle
    !
    if (q%get_flag(active=LDUM)) then
      call logmsg('***** Processing quad '//ta([q%gid])//' *****')
      call q%mf6_write_data()
      n = n + 1
      call logmsg('***** '//ta([100.*real(n,R4B)/n_act],'(f6.2)')//' % *****')
      if (writelog) then
        f_log = trim(d)//'done_lid_'//ta([lid],'(i10.10)')//'.txt'
        call open_file(f_log, iu, 'w'); close(iu)
      end if
    end if
    !
  end do
  !
  return
end subroutine quad_mf6_data_write_old
  
subroutine quad_mf6_data_write()
! ******************************************************************************
!
!    SPECIFICATIONS:
  ! -- local
  type(tDataModelData), pointer :: dmdat => null()
  !
  character(len=MXSLEN) :: d, f_log
  logical :: writelog
  integer(I4B) :: lid0, lid1, n_act, iu, idat
  !
  integer(I4B), dimension(:), allocatable :: i4a
  real(R8B), dimension(:), allocatable :: r8a
  real(R8B), dimension(:,:), allocatable :: r8x
! ------------------------------------------------------------------------------
  ! set the ranges
  if (lid_min > 0) then
     lid0 = max(1,lid_min)
  else
    lid0 = 1
  end if
  if (lid_max > 0) then
     lid1 = min(xq%n,lid_max)
  else
    lid1 = xq%n
  end if
  if (lid0 > lid1) then
    call errmsg('quad_mf6_data_write: lid_min > lid_max.')
  end if
  !
  ! determine total number of active quads
  n_act = 0
  do lid = lid0, lid1
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) n_act = n_act + 1
  end do
  !
  if (len_trim(d_log) > 0) then
    writelog = .true.
    d = trim(d_log)
    call create_dir(d, .true.)
  else
    writelog = .false.
  end if
  !
  n = 0
  do lid = lid0, lid1
    q => xq%get_quad(lid)
    !
    ! DEBUG:
    !if (q%gid /= 637452) cycle
    !
    if (q%get_flag(active=LDUM)) then
      call logmsg('***** Processing quad '//ta([q%gid])//' *****')
      !
      call q%grid_init()
      do idat = 1, q%dat_mod%get_ndat()
        call q%dat_mod%get_dat(idat, dmdat)
        call q%mf6_get_data(dmdat, i4a, r8a, r8x)
        call mf6_data_write(q%disu%wbd, dmdat%id, dmdat%i_out_file_type, i4a, r8a, r8x)
      end do
      call q%disu%wbd%write_csv()
      !
      ! clean up disu
      call q%disu%clean(); deallocate(q%disu); q%disu => null()
      !
      n = n + 1
      call logmsg('***** '//ta([100.*real(n,R4B)/n_act],'(f6.2)')//' % *****')
      if (writelog) then
        f_log = trim(d)//'done_lid_'//ta([lid],'(i10.10)')//'.txt'
        call open_file(f_log, iu, 'w'); close(iu)
      end if
    end if
    !
  end do
  !
  return
end subroutine quad_mf6_data_write

subroutine quad_mf6_data_write_merge()
! ******************************************************************************
!
!    SPECIFICATIONS:
  ! -- local
  type(tCsv), pointer :: csv => null()
  type(tMF6Disu), pointer :: disu => null()
  type(tDataModelData), pointer :: dmdat => null()
  type(tMf6Wbd), pointer :: wbd => null()
  !
  logical :: li4a, lr8a, lr8x, lfirst, lwritelog
  !
  character(len=MXSLEN) :: d, f_csv_dat, f_binpos, f_log, id_pref, id
  character(len=MXSLENLONG) :: s_long
  !
  integer(I4B), dimension(:), allocatable :: im_arr, lid_arr, qi4a, i4a
  integer(I4B), dimension(:), allocatable ::  nodesa, nodesa_offset
  integer(I4B) :: im, im0, im1, nim, nlid, idat, i, j, k, ilm, nodes_sum
  integer(I4B) :: nodes, nodes_im, n, m, ndat_q, ndat, ic, ir, nc
  integer(I4B) :: i_out_file_type, iu, nr_max
  !
  real(R8B), dimension(:), allocatable :: qr8a, r8a
  real(R8B), dimension(:,:), allocatable :: qr8x, r8x
! ------------------------------------------------------------------------------
  !
  if (len_trim(d_log) > 0) then
    lwritelog = .true.
    d = trim(d_log)
    call create_dir(d, .true.)
  else
    lwritelog = .false.
  end if
  !
  allocate(csv)
  call csv%read(f_in_csv_merge)
  call csv%get_column(key='lid', i4a=im_arr)
  nim = size(im_arr)
  !
  ! set the ranges
  if (lid_min > 0) then
    im0 = max(1,lid_min)
  else
    im0 = 1
  end if
  if (lid_max > 0) then
    im1 = min(nim,lid_max)
  else
    im1 = nim
  end if
  if (im0 > im1) then
    call errmsg('quad_mf6_data_write_merge: im_min > im_max.')
  end if
  !
  do im = im0, im1
    ! get the local lids
    call csv%get_val(ir=im, ic=csv%get_col('lid_merged'), cv=s_long)
    call parse_line(s=s_long, i4a=lid_arr, token_in=';'); nlid = size(lid_arr)
    !
    ! read the nodes array
    if (allocated(nodesa)) deallocate(nodesa)
    allocate(nodesa(nlid))
    do j = 1, nlid
      lid = lid_arr(j)
      q => xq%get_quad(lid)
      call q%get_prop_csv(key='nodes', i4v=nodesa(j))
    end do
    !
    ! check
    call csv%get_val(ir=im, ic=csv%get_col('nodes'), i4v=nodes_im)
    if (nodes_im /= sum(nodesa)) then
      call errmsg('quad_mf6_data_write_merge: invalid nodes.')
    end if
    !
    ! determine the node offsets
    if (allocated(nodesa_offset)) deallocate(nodesa_offset)
    allocate(nodesa_offset(nlid))
    nodes = 0
    do i = 1, nlid
      nodesa_offset(i) = nodes
      nodes = nodes + nodesa(i)
    end do
    !
    ! read the disu
    nr_max = 9 + 13*nlid + xq%dat_mods%n_uni_map
    call csv%get_val(ir=im, ic=csv%get_col('csv_dat'), cv=f_csv_dat)
    do j = 1, nlid
      lid = lid_arr(j)
      q => xq%get_quad(lid)
      if (.not.q%get_flag(active=LDUM)) then
        call errmsg('quad_mf6_data_write_merge: inactive quad.')
      end if
      id_pref = '_'//ta([lid_arr(j)])
      call logmsg('Initializing grid for lid '//ta([lid_arr(j)])// &
        ' ('//ta([j])//'/'//ta([nlid])//')...')
      call q%grid_init(f_csv_dat, id_pref, nr_max)
    end do
    !
    ! set the wbd
    f_binpos = trim(strip_ext(f_csv_dat))//'.binpos'
    allocate(wbd)
    call wbd%read_csv(f_csv_dat, nr_max=MAX_NR_CSV)
    wbd%f_binpos = f_binpos
    !
    ! loop over the unique mapping ids
    do i = 1, xq%dat_mods%n_uni_map
    !do i = 46, 46
      ndat = 0
      id = xq%dat_mods%uni_map_id(i)
      lfirst = .true.
      !
      ! set the array
      if (allocated(i4a)) deallocate(i4a)
      if (allocated(r8a)) deallocate(r8a)
      if (allocated(r8x)) deallocate(r8x)
      !
      ! get the data and merge
      do j = 1, nlid
        lid = lid_arr(j)
        q => xq%get_quad(lid)
        call q%get_prop_csv(ikey=i_lay_mod, i4v=ilm)
        idat = xq%dat_mods%uni_map_ir(ilm, i)
        if (idat > 0) then
          !
          call q%dat_mod%get_dat(idat, dmdat)
          if (id /= dmdat%id) then
            call errmsg('quad_mf6_data_write_merge: program error.')
          end if
          !
          ! ----- read the data -----
          !do k = 1, 100
          !  call q%dat_mod%get_dat(idat, dmdat)
          ! ----- read the data -----
          call q%mf6_get_data(dmdat, qi4a, qr8a, qr8x)
          !end do
          ! ----- read the data -----
          !
          if (lfirst) then
            i_out_file_type = dmdat%i_out_file_type
            lfirst = .false.
          end if
          !
          ! set the data
          li4a = .false.; lr8a = .false.; lr8x = .false.
          ndat_q = 0
          if (allocated(qi4a)) then
            li4a = .true.
            ndat_q = size(qi4a)
          end if
          if (allocated(qr8a)) then
            lr8a = .true.
            if (ndat_q > 0) then
              if (size(qr8a) /= ndat_q) then
                call errmsg('quad_mf6_data_write_merge: invalid data size.')
              end if
            else
              ndat_q = size(qr8a)
            end if
          end if
          if (allocated(qr8x)) then
            lr8x = .true.
            if (ndat_q > 0) then
              if (size(qr8x,2) /= ndat_q) then
                call errmsg('quad_mf6_data_write_merge: invalid data size.')
              end if
            else
              ndat_q = size(qr8a)
            end if
          end if
          !
          if ((.not.li4a).and.(.not.lr8a).and.(.not.lr8x)) then
            call logmsg('WARNING: nothing to do for '//trim(id)//'!')
            cycle
          end if
          !
          if (li4a) then
            if (.not.allocated(i4a)) then
              allocate(i4a(nodes_im))
            end if
            do ir = 1, ndat_q
              n = qi4a(ir) + nodesa_offset(j)
              i4a(ir+ndat) = n
            end do
          end if
          if (lr8a) then
            if (.not.allocated(r8a)) then
              allocate(r8a(nodes_im))
            end if
            do ir = 1, ndat_q
              r8a(ir+ndat) = qr8a(ir)
            end do
          end if
          if (lr8x) then
            nc = size(qr8x,1)
            if (.not.allocated(r8x)) then
              allocate(r8x(nc,nodes_im))
            end if
            do ir = 1, ndat_q
              do ic = 1, nc
                r8x(ic,ir+ndat) = qr8x(ic,ir)
              end do
            end do
          end if
          ndat = ndat + ndat_q
        end if !idat
      end do
      !
      ! write the data
      if (ndat > 0) then
        call logmsg('Writing MODFLOW 6 data for: '//trim(id)//'...')
        call mf6_data_write(wbd, id, i_out_file_type, i4a, r8a, r8x, ndat)
        call wbd%write_csv()
      end if
      !
    end do
    !
    ! clean up
    if (allocated(qi4a)) deallocate(qi4a)
    if (allocated(qr8a)) deallocate(qr8a)
    if (allocated(qr8x)) deallocate(qr8x)
    if (allocated(i4a)) deallocate(i4a)
    if (allocated(r8a)) deallocate(r8a)
    if (allocated(r8x)) deallocate(r8x)
    call wbd%clean(); deallocate(wbd)
    do j = 1, nlid
      lid = lid_arr(j)
      q => xq%get_quad(lid)
      call q%disu%clean()
    end do
    !
    if (lwritelog) then
      f_log = trim(d)//'done_im_'//ta([im],'(i10.10)')//'.txt'
      call open_file(f_log, iu, 'w'); close(iu)
    end if
  end do
  !
  call csv%clean(); deallocate(csv); csv => null()
  !
  return
end subroutine quad_mf6_data_write_merge
  
subroutine quad_set_mod_dir()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  call logmsg('Setting the model directories...')
  call xq%set_mod_dir(mod_root_dir, mod_sub_dir_fields)
  !
  return
end subroutine quad_set_mod_dir

subroutine quad_grid_gen()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- local
  type(tMF6Disu), pointer :: disu => null()
  type(tCsv), pointer :: csv => null(), csv_old => null()
  !
  integer(I4B), parameter :: perc_intv = 10
  !
  character(len=1) :: slash
  character(len=MXSLEN) :: f_csv_dat, f_binpos, s_lay_act
  logical :: lactive, lskip, lwrite
  integer(I4B), dimension(:), allocatable :: lay, lid_map, lid_arr
  integer(I4B) :: ncell_tot, ncell, nja, nlay_act, ngrid, write_csv_delta
  integer(I4B) :: lid0, lid1, n_act, i, n
  real(R8B) :: cs_min_rea, cs_max_rea
! ------------------------------------------------------------------------------
  ncell_tot = 0; n = 0
  write_csv_delta = max(1,int(real(xq%n_act,R4B)/perc_intv,I4B))
  !
  if (lwrite_props) then
    if (len_trim(f_out_csv) == 0) then
      call errmsg('f_out_csv not set.')
    end if
    csv => xq%props%csv
    call csv%add_key('nodes')
    call csv%add_key('nja')
    call csv%add_key('nlay_act')
    call csv%add_key('ngrid_lev')
    call csv%add_key('lay_act')
    call csv%add_key('cs_min_rea')
    call csv%add_key('cs_max_rea')
    if (lwrite_disu) call csv%add_key('csv_dat')
  end if
  if (loverwrite_props) then
    if (fileexist(f_out_csv)) then
      allocate(csv_old)
      call csv_old%read(f_out_csv)
      call csv_old%get_column(key=fields(i_lid),i4a=lid_arr)
      n = maxval(lid_arr); allocate(lid_map(n)); lid_map = 0
      do i = 1, size(lid_arr)
        lid = lid_arr(i)
        lid_map(lid) = i
      end do
    end if
  end if
  !
  ! set the ranges
  if (lid_min > 0) then
     lid0 = max(1,lid_min)
  else
    lid0 = 1
  end if
  if (lid_max > 0) then
     lid1 = min(xq%n,lid_max)
  else
    lid1 = xq%n
  end if
  if (lid0 > lid1) then
    call errmsg('lid_min > lid_max.')
  end if
  !
  ! determine total number of active quads
  n_act = 0
  do lid = lid0, lid1
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) n_act = n_act + 1
  end do
  !
  n = 0
  slash = get_slash()
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    lactive = q%get_flag(active=LDUM)
    !
    lwrite = .false.
    if (lactive) then
      if (loverwrite_props) then
        lwrite = .true.
      end if
    end if
    !
    ! DEBUG:
    !if (lid == 24728) then
    !  write(*,*) '@@@@ 24728'
    !end if
    !
    if (lactive) then
      !
      lskip = .false.
      if ((lid >= lid0).and.(lid <= lid1)) then
        call logmsg('***** Processing quad '//ta([q%gid])//' *****')
        lwrite = .true.
        f_csv_dat = trim(q%mod_dir)//slash//'dat.csv'
        
        ! create the grid and store the arrays
        allocate(disu)
        call disu%init()
        if(.not.lwrite_asc) then
          f_binpos = trim(strip_ext(f_csv_dat))//'.binpos'
          call disu%init_csv(f_csv_dat, f_binpos)
        else
          call disu%init_csv(f_csv_dat)
        end if
        !
        ! call the main grid generation subroutine
        call q%grid_gen(nlay_act, lay, lskip, disu)
        !
        if (lwrite_disu) then
          if (lwrite_asc) then
            call disu%write(write_opt=1, write_asc=.true., d_out=trim(q%mod_dir)//'\')
            call disu%write(write_opt=2, write_asc=.true., d_out=trim(q%mod_dir)//'\')
          else
            call disu%write(write_opt=1)
            call disu%write(write_opt=2)
          end if
        end if
        s_lay_act = ta(lay,sep_in=';')
        !
        ncell = disu%nodes
        nja   = disu%nja 
        ngrid = disu%ngrid
        cs_min_rea = disu%cs_min_rea
        cs_max_rea = disu%cs_max_rea
      else
        if (associated(csv_old)) then
          call logmsg('***** Using data from old csv for quad '//ta([q%gid])//' *****')
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('nodes'), i4v=ncell)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('nja'), i4v=nja)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('nlay_act'), i4v=nlay_act)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('ngrid_lev'), i4v=ngrid)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('csv_dat'), cv=f_csv_dat)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('lay_act'), cv=s_lay_act)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('cs_min_rea'), r8v=cs_min_rea)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('cs_max_rea'), r8v=cs_max_rea)
        else
          cycle
        end if
      end if
      !
      n = n + 1; ncell_tot = ncell_tot + ncell
      call logmsg('***** '//ta([100.*real(n,R4B)/n_act],'(f6.2)')//' %: '// &
        trim(adjustl(ta([real(ncell_tot,R4B)/1000000.],'(f10.2)')))//' M cells *****')
      if (lskip) then
        call logmsg('********** SKIPPING QUAD *********')
        call q%set_flag(active=.false.)
      end if
      if (lwrite) then
        call q%set_prop_csv(key='nodes', i4v=ncell)
        call q%set_prop_csv(key='nja', i4v=nja)
        call q%set_prop_csv(key='nlay_act', i4v=nlay_act)
        call q%set_prop_csv(key='ngrid_lev', i4v=ngrid)
        call q%set_prop_csv(key='lay_act', cv=trim(s_lay_act))
        call q%set_prop_csv(key='cs_min_rea', r8v=cs_min_rea)
        call q%set_prop_csv(key='cs_max_rea', r8v=cs_max_rea)
        if (lwrite_disu) then
          call q%set_prop_csv(key='csv_dat', cv=trim(f_csv_dat))
        end if
        if ((n == xq%n_act).or.(mod(n,write_csv_delta) == 0).or.(lid == lid1)) then
          csv%file = trim(f_out_csv)
          if (loverwrite_props) then
            call csv%write()
          else
            call csv%write(ir0=lid0,ir1=lid1)
          end if
        end if
      end if
      call disu%clean(); deallocate(disu)
    end if
    !
    !if (lwrite_props.and.((n == xq%n_act).or.(mod(n,write_csv_delta) == 0))) then
    !  csv%file = trim(f_out_csv); call csv%write()
    !end if
  end do
  !
  ! clean-up
  if (associated(csv_old)) then
    call csv_old%clean(); deallocate(csv_old); csv_old => null()
  end if
  if (allocated(lid_arr)) deallocate(lid_arr)
  if (allocated(lid_map)) deallocate(lid_map)
  !
  return
end subroutine quad_grid_gen

subroutine quad_grid_gen_merge()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- local
  type(tCsv), pointer :: csv => null()
  type(tMF6Disu), dimension(:), pointer :: disu_arr => null()
  type(tMF6Disu), pointer :: disu => null(), disu_merge => null()
  character(len=1) :: slash
  character(len=MXSLEN) :: d, f_csv_dat, f_binpos, s_lay_act, id_pref
  logical :: lactive, lskip
  integer(I4B), dimension(:), allocatable :: lay, lid_arr, imerge, lid2imerge
  integer(I4B) :: lid, nim, im, im0, im1, jm, km, nlid, nlay_act, n, nodes_tot
! ------------------------------------------------------------------------------
  if (lwrite_props) then
    if (len_trim(f_out_csv) == 0) then
      call errmsg('f_out_csv not set.')
    end if
    allocate(csv)
    call csv%init(file=f_out_csv, &
      hdr_keys=['lid', 'nlid_merged', 'lid_merged', 'nodes', 'nja', 'csv_dat'], &
      nr_max=xq%n, hdr_i_type=[i_i4, i_i4, i_c, i_i4, i_i4, i_c])
  end if
  !
  ! read the parts to be merged
  call xq%props%csv%get_column(key=fields(i_merge),i4a=imerge)
  if (size(imerge) /= xq%n_act) then
    call errmsg('quad_grid_gen_merge: invalid size.')
  end if
  allocate(lid2imerge(xq%n)); lid2imerge = 0
  call xq%props%csv%get_column(key=fields(i_lid),i4a=lid_arr)
  do i = 1, size(lid_arr)
    lid = lid_arr(i)
    lid2imerge(lid) = i
  end do
  deallocate(lid_arr)
  !
  nim = maxval(imerge)
  !
  ! set the ranges
  if (lid_min > 0) then
    im0 = max(1,lid_min)
  else
    im0 = 1
  end if
  if (lid_max > 0) then
    im1 = min(nim,lid_max)
  else
    im1 = nim
  end if
  if (im0 > im1) then
    call errmsg('quad_grid_gen_merge: im_min > im_max.')
  end if
  !
  slash = get_slash()
  n = 0; nodes_tot = 0
  do im = im0, im1
!  do im = 1, 1
    d = trim(mod_root_dir)//slash//ta([im]); call create_dir(d, .true.)
    f_csv_dat = trim(d)//slash//'dat.csv'
    if(.not.lwrite_asc) then
      f_binpos = trim(strip_ext(f_csv_dat))//'.binpos'
    end if
    !
    ! count
    nlid = 0
    do lid = 1, xq%n
      q => xq%get_quad(lid)
      lactive = q%get_flag(active=LDUM)
      if (lactive) then
        jm = imerge(lid2imerge(lid))
        if (jm == im) then
          nlid = nlid + 1
        end if
      end if
    end do
    if (nlid == 0) then
      call errmsg('quad_grid_gen_merge: nothing to merge')
    end if
    if (allocated(lid_arr)) deallocate(lid_arr)
    allocate(lid_arr(nlid))
    nlid = 0
    do lid = 1, xq%n
      q => xq%get_quad(lid)
      lactive = q%get_flag(active=LDUM)
      if (lactive) then
        jm = imerge(lid2imerge(lid))
        if (jm == im) then
          nlid = nlid + 1
          lid_arr(nlid) = lid
        end if
      end if
    end do
    !
    ! allocate the disuu array
    allocate(disu_arr(nlid))
    !
    ! fill disu
    km = 0
    do lid = 1, xq%n
      q => xq%get_quad(lid)
      lactive = q%get_flag(active=LDUM)
      if (lactive) then
        jm = imerge(lid2imerge(lid))
        if (jm == im) then
          km = km + 1
          !
          disu => disu_arr(km)
          call disu%init()
          !
          ! generate the seperate grid
          call q%grid_gen(nlay_act, lay, lskip, disu)
        end if
      end if
    end do
    !
    ! initialize the merged disu
    allocate(disu_merge)
    call disu_merge%init()
    if(.not.lwrite_asc) then
      call disu_merge%init_csv(f_csv_dat, f_binpos)
    else
      call disu_merge%init_csv(f_csv_dat)
    end if
    !
    ! merge the separate disu
    call xq%merge_disu(lid_arr, disu_arr, disu_merge)
    !
    ! write the merged disu data
    if (lwrite_disu) then
      call disu_merge%write(write_opt=1)
    end if
    !
    ! write all seperate disu data; and clean
    do km = 1, nlid
      disu => disu_arr(km)
      !
      if(.not.lwrite_asc) then
        call disu%init_csv(f_csv_dat, f_binpos, reuse=.true.)
      else
        call disu%init_csv(f_csv_dat, reuse=.true.)
      end if
      !
      id_pref = '_'//ta([lid_arr(km)])
      if (lwrite_disu) then
        call disu%write(write_opt=1, id_pref=id_pref)
        call disu%write(write_opt=2, id_pref=id_pref)
      end if
      call disu%clean()
    end do
    !
    n = n + 1; nodes_tot = nodes_tot + disu_merge%nodes
    call logmsg('***** '//ta([100.*real(n,R4B)/(im1-im0+1)],'(f6.2)')//' %: '// &
        trim(adjustl(ta([real(nodes_tot,R4B)/1000000.],'(f10.2)')))//' M nodes *****')
    !
    ! write to the csv
    if (lwrite_props) then
      call csv%increase_nr()
      call csv%set_val_by_key(key='lid', i4v=im)
      call csv%set_val_by_key(key='nlid_merged', i4v=nlid)
      call csv%set_val_by_key(key='lid_merged', cv=ta(lid_arr, sep_in=';'))
      call csv%set_val_by_key(key='nodes', i4v=disu_merge%nodes)
      call csv%set_val_by_key(key='nja', i4v=disu_merge%nja)
      call csv%set_val_by_key(key='csv_dat', cv=f_csv_dat)
    end if
    !
    ! clean up
    deallocate(disu_arr); disu_arr => null()
    call disu_merge%clean(); deallocate(disu_merge); disu_merge => null()
  end do
  !
  ! write the csv file and clean up
  if (lwrite_props) then
    call csv%write()
    call csv%clean(); deallocate(csv)
  end if
  !
  return
end subroutine quad_grid_gen_merge

subroutine quad_partition()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- modules
  use metis_module, only: tMetis
  !
  ! -- local
  real(R8B), parameter :: tolimbal = 1.00001d0  ! load imbalance tolerance
  !
  type(tCsv), pointer :: csv => null(), csv_sel => null(), csv_part => null()
  type(tMetis), pointer :: met => null()
  logical :: lnpart, lwpartinfo
  character(len=MXSLEN) :: s, key
  character(len=MXSLEN), dimension(:), allocatable :: sel_val_arr, key_arr
  character(len=MXSLENLONG), dimension(:), allocatable :: nbr_sa, nbr_n_sa
  real(R4B) :: r
  integer(I4B), dimension(:), allocatable :: sel_npart_arr, sel_nodes_arr
  integer(I4B), dimension(:), allocatable :: sel_rows, nbr, nbr_n
  integer(I4B), dimension(:), allocatable :: lid_arr, nodes_arr, lid_map
  integer(I4B), dimension(:), allocatable :: i_type_arr
  integer(I4B) :: i_sel, n_sel, lid_max, i, j, n, nvtxs, nja, lid, iact
  integer(I4B) :: nodes, nodes_tot, nodes_avg, nodes_min, nodes_max
  integer(I4B) :: nparts, part, nparts_max
! ------------------------------------------------------------------------------
  if (len_trim(f_part_csv) > 0) then
    lwpartinfo = .true.
  else
    lwpartinfo = .false.
  end if
  !
  allocate(csv)
  call csv%read(f_in_csv)
  !
  ! check presence of the fields
  do i = 1, n_part_field
    part_fields(i) = change_case(part_fields(i), 'l')
    if (.not.csv%exist_col(part_fields(i))) then
      call errmsg('Column '//trim(part_fields(i))//' does not exist.')
    end if
  end do
  sel_field = change_case(sel_field, 'l')
  if (.not.csv%exist_col(sel_field)) then
    call errmsg('Column '//trim(sel_field)//' does not exist.')
  end if
  !
  if (len_trim(sel_npart) > 0) then
    lnpart = .true.
    call parse_line(s=sel_npart, i4a=sel_npart_arr, token_in=',')
    n = size(sel_npart_arr)
  else
    lnpart = .false.
    call parse_line(s=sel_nodes, i4a=sel_nodes_arr, token_in=',')
    n = size(sel_nodes_arr)
  end if
  !
  call parse_line(s=sel_val, sa_short=sel_val_arr, token_in=',')
  if (n /= size(sel_val_arr)) then
    call errmsg('key sel_npart and/or sel_nodes and/or sel_val.')
  else
    n_sel = n
  end if
  !
  ! first, determine the number of partitions
  if (.not.lnpart) then
    allocate(sel_npart_arr(n_sel)); sel_npart_arr = 0
    do i_sel = 1, n_sel
      ! get the selected csv
      call csv%get_selection(sel_field, sel_val_arr(i_sel), csv_sel, sel_rows)
      !
      if (associated(csv_sel)) then
        call csv_sel%get_column(key=part_fields(i_nodes), i4a=nodes_arr)
        nodes_tot = sum(nodes_arr)
        nodes = sel_nodes_arr(i_sel)
        nparts = nint(real(nodes_tot,R4B)/real(nodes,r4B),I4B)
        sel_npart_arr(i_sel) = max(nparts,1)
      end if
    end do
  end if
  !
  if (lwpartinfo) then
    nparts_max = maxval(sel_npart_arr)
    allocate(csv_part)
    allocate(key_arr(n_sel*2+1), i_type_arr(n_sel*2+1))
    key_arr(1) = 'part'
    i_type_arr(1) = I_I4
    do i_sel = 1, n_sel
      i = 2*(i_sel-1)+2
      key_arr(i)   = 'load_'//ta([i_sel])
      i_type_arr(i) = I_R8
      key_arr(i+1) = 'load_imbalance_'//ta([i_sel])
      i_type_arr(i+1) = I_R8
    end do
    call csv_part%init(file=f_part_csv, hdr_keys=key_arr, &
      nr=nparts_max, hdr_i_type=i_type_arr)
    ic = csv_part%get_col('part')
    do ir = 1, nparts_max
      call csv_part%set_val(ic=ic, ir=ir, i4v=ir, create_s=.true.)
    end do
  end if
  !
  csv_sel => null()
  allocate(met)
  !
  do i_sel = 1, n_sel
    ! get the selected csv
    call csv%get_selection(sel_field, sel_val_arr(i_sel), csv_sel, sel_rows)
    !
    if (associated(csv_sel)) then
      !
      ! get all necessary arrays
      call csv_sel%get_column(key=fields(i_lid),i4a=lid_arr)
      call csv_sel%get_column(key=part_fields(i_nodes), i4a=nodes_arr)
      call csv_sel%get_column(key=part_fields(i_neighbor_lid), ca=nbr_sa)
      call csv_sel%get_column(key=part_fields(i_neighbor_nr_cells), ca=nbr_n_sa)
      !
      nodes_tot = sum(nodes_arr)
      nodes_avg = int(real(nodes_tot,R4B)/size(nodes_arr),I4B)
      !nodes_min = 1
      !nodes_max = huge(nodes_max) !TODO
      nvtxs = csv_sel%get_nr()
      nparts = sel_npart_arr(i_sel)
      !
      ! allocate for METIS
      call met%clean()
      allocate(met%nparts);                      met%nparts = int(nparts,I8B)
      allocate(met%ncon);                        met%ncon   = I8ONE
      allocate(met%tpwgts(met%nparts*met%ncon)); met%tpwgts = R8ONE/real(met%nparts,R8B)
      allocate(met%ubvec(met%ncon));             met%ubvec  = tolimbal
      allocate(met%nvtxs);                       met%nvtxs  = int(nvtxs,I8B)
      allocate(met%vwgt(met%nvtxs));             met%vwgt   = I8ZERO
      allocate(met%vsize(met%nvtxs));            met%vsize  = I8ONE
      allocate(met%part(met%nvtxs));             met%part   = I8ZERO
      allocate(met%xadj(met%nvtxs+1));           met%xadj   = I8ZERO
      call met%set_opts()
      !
      lid_max= maxval(lid_arr)
      if (allocated(lid_map)) deallocate(lid_map)
      allocate(lid_map(lid_max)); lid_map = 0
      do i = 1, nvtxs
        lid = lid_arr(i)
        lid_map(lid) = i
      end do
      !
      ! fill the array
      do iact = 1, 2
        nja = 0
        do i = 1, nvtxs
          call parse_line(nbr_sa(i), i4a=nbr, token_in=';');
          !
          if (iact == 2) then
            call parse_line(nbr_n_sa(i), i4a=nbr_n, token_in=';')
            if (size(nbr) /= size(nbr_n)) then
              call errmsg('Neigboring data for partitioning is inconsistent.')
            end if
            n = nodes_arr(i)
            !n = max(n,nodes_min); n = min(n,nodes_max)
            met%vwgt(i) = int(n,I8B)
            met%xadj(i+1) = met%xadj(i)
          end if
          do j = 1, size(nbr)
            lid = nbr(j)
            if (lid <= lid_max) then
              n = lid_map(lid)
            else
              n = 0
            end if
            if (n > 0) then
              nja = nja + 1
              if (iact == 2) then
                met%adjncy(nja) = int(n-1,I8B)
                met%adjwgt(nja) = int(nbr_n(j),I8B)
                !met%adjwgt(nja) = 1 !DEBUG
                met%xadj(i+1)   = met%xadj(i+1) + 1
              end if
            end if
          end do
        end do
        if (iact == 1) then
          allocate(met%adjncy(max(nja,1)))
          allocate(met%adjwgt(max(nja,1)))
        end if
      end do
      !
      ! Perform METIS partitioning
      if (met%nparts > 1) then
        if (nvtxs > met%nparts) then 
          call met%recur()
          !call met%kway()
        else
          n = 0
          do i = 1, nvtxs
            met%part(i) = int(n,I8B)
            n = n + 1 
          end do
          call met%calc_imbal()
        end if
      end if
      !
      ! set the partition number
      key = 'part_'//ta([i_sel]); call csv%add_key(key)
      ic = csv%get_col(key)
      do i = 1, nvtxs
        ir = sel_rows(i)
        part = int(met%part(i),I4B) + 1
        call csv%set_val(ic=ic, ir=ir, i4v=part, create_s=.true.)
      end do
      !
      ! store the partition information
      if (lwpartinfo) then
        if (met%nparts > 1) then
          ic = csv_part%get_col('load_'//ta([i_sel]))
          do ir = 1, met%nparts
            call csv_part%set_val(ic=ic, ir=ir, r8v=real(met%load(ir),R8B), create_s=.true.)
          end do
          do ir = met%nparts + 1, nparts_max
            call csv_part%set_val(ic=ic, ir=ir, r8v=R8ZERO, create_s=.true.)
          end do
          ic = csv_part%get_col('load_imbalance_'//ta([i_sel]))
          do ir = 1, met%nparts
            call csv_part%set_val(ic=ic, ir=ir, r8v=met%loadimbal(ir), create_s=.true.)
          end do
          do ir = met%nparts + 1, nparts_max
            call csv_part%set_val(ic=ic, ir=ir, r8v=R8ZERO, create_s=.true.)
          end do
        else
          ic = csv_part%get_col('load_'//ta([i_sel]))
          call csv_part%set_val(ic=ic, ir=1, r8v=real(sum(nodes_arr),R8B), create_s=.true.)
          ic = csv_part%get_col('load_imbalance_'//ta([i_sel]))
          call csv_part%set_val(ic=ic, ir=1, r8v=R8ZERO, create_s=.true.)
        end if
      end if
    else
      call logmsg('No data for selection '//ta([i_sel])//'...')
    end if
  end do
  !
  ! write the csv
  csv%file = trim(f_out_csv); call csv%write()
  if (lwpartinfo) then
    call csv_part%write(r8mv=R8ZERO)
  end if
  !
  ! clean up
  if (allocated(sel_val_arr))   deallocate(sel_val_arr)
  if (allocated(nbr_sa))        deallocate(nbr_sa)
  if (allocated(nbr_n_sa))      deallocate(nbr_n_sa)
  if (allocated(sel_npart_arr)) deallocate(sel_npart_arr)
  if (allocated(sel_rows))      deallocate(sel_rows)
  if (allocated(nbr))           deallocate(nbr)
  if (allocated(nbr_n))         deallocate(nbr_n)
  if (allocated(lid_arr))       deallocate(lid_arr)
  if (allocated(nodes_arr))     deallocate(nodes_arr)
  if (allocated(lid_map))       deallocate(lid_map)
  if (associated(met)) then
    call met%clean(); deallocate(met); met => null()
  end if
  if (associated(csv)) then
    call csv%clean(); deallocate(csv); csv => null()
  end if
  if (associated(csv_sel)) then
    call csv_sel%clean(); deallocate(csv_sel); csv_sel => null()
  end if
  if (associated(csv_part)) then
    call csv_part%clean(); deallocate(csv_part); csv_part => null()
  end if
  !
  return
end subroutine quad_partition

subroutine quad_mf6_write_heads(lmerge)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- local
  integer(I4B) :: lid0, lid1
  logical, intent(in) :: lmerge
! ------------------------------------------------------------------------------

  ! set the ranges
  if (lid_min > 0) then
     lid0 = max(1,lid_min)
  else
    lid0 = 1
  end if
  if (lid_max > 0) then
     lid1 = min(xq%n,lid_max)
  else
    lid1 = xq%n
  end if
  if (lid0 > lid1) then
    call errmsg('quad_mf6_write_heads: lid_min > lid_max.')
  end if
  if (lmerge) then
    call xq%write_mf6_heads(lid0, lid1, head_layers, kper_beg, kper_end, &
      tile_nc, tile_nr, f_vrt, write_nod_map, vtk_lid, vtk_csv, f_in_csv_merge)
  else
    call xq%write_mf6_heads(lid0, lid1, head_layers, kper_beg, kper_end, &
      tile_nc, tile_nr, f_vrt, write_nod_map, vtk_lid, vtk_csv)
  end if
  !
  return
end subroutine quad_mf6_write_heads

subroutine quad_mf6_write_chd()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- local
  type(tQuad), pointer :: q => null()
  integer(I4B), dimension(:), allocatable :: lid_act, lids
  integer(I4B) :: i, lid
! ------------------------------------------------------------------------------
  call logmsg('Processing constant head boundary heads.')
  !
  ! filter for selected quads
  if (len_trim(chd_lid) > 0) then
    allocate(lid_act(xq%n)); lid_act = 0
    call parse_line(s=chd_lid, i4a=lids, token_in=',')
    do i = 1, size(lids)
      lid = lids(i)
      lid_act(lid) = 1
    end do
    do lid = 1, xq%n
      q => xq%get_quad(lid)
      if (q%get_flag(active=LDUM)) then
        if (lid_act(lid) == 0) then
          call q%set_flag(active=.false.)
          xq%n_act = xq%n_act - 1
        end if
      end if
    end do
    deallocate(lid_act, lids)
  end if
  !
  call xq%write_mf6_bnd_heads(kper_beg, kper_end)

  return
end subroutine quad_mf6_write_chd

subroutine quad_mf6_write_wtd()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- local
  type(tVrt), pointer       :: vrt_head => null(), vrt_node => null()
  type(tHdr), pointer       :: hdr_ahn => null()
  type(tMultiGrid), pointer :: mg_ahn => null()
  type(tGrid), pointer      :: g_ahn => null()
  type(tBB)  :: tile_bbi
  type(tBB), dimension(:), allocatable :: bbi_arr
  type(tBBx) :: tile_bbx
  !
  logical :: found
  character(len=MXSLEN) :: f
  real(R8B), dimension(:), allocatable :: cs_read
  real(R8B) :: cs_min, cs, xc, yc
  real(R4B), dimension(:,:), allocatable :: head, ahn, wtd
  real(R4B) :: mvhead, mvahn, head_val, ahn_val
  integer(I4B), dimension(:,:), allocatable :: node
  integer(I4B) :: itile, ntiles, mvnode, node_max, node_min, n, max_nc, max_nr
  integer(I4B) :: nlev, ilev, bs, jc, jr
! ------------------------------------------------------------------------------
  !
  allocate(vrt_head, vrt_node)
  call vrt_head%init(f_in_vrt_1)
  call vrt_node%init(f_in_vrt_2)
  f_in_flt = strip_ext(f_in_flt)
  !
  ntiles = vrt_head%get_ntiles()
  if (ntiles /= vrt_node%get_ntiles()) then
    call errmsg('Invalid number of tiles found.')
  end if
  !
  allocate(mg_ahn)
  !
  ! read the tiles
  do itile = 1, ntiles
    !if (itile /= 66) cycle !DEBUG
    call logmsg('***** Processing tile '//ta([itile])//'/'//ta([ntiles])//' *****')
    call vrt_head%read_full_tile(itile=itile, xr4=head, mvr4=mvhead)
    call vrt_node%read_full_tile(itile=itile, xi4=node, mvi4=mvnode)
    call vrt_head%get_bb(itile=itile, src_bbi=tile_bbi, src_bbx=tile_bbx)
    !
    ! check
    if ((size(head,1) /= size(node,1)).or.(size(head,2) /= size(node,2))) then
      call errmsg('Non matching tile sizes.')
    else
      nc = size(head,1); nr = size(head,2)
    end if
    !
    ! determine minimum/maximum node
    found = .false.
    node_min = huge(node_min); node_max = -huge(node_max)
    do ir = 1, nr; do ic = 1, nc
      if (node(ic,ir) /= mvnode) then
        found = .true.
        node_min = min(node_min,node(ic,ir))
        node_max = max(node_max,node(ic,ir))
      end if
    end do; end do
    !
    ! intitialize the wtd
    if (allocated(wtd)) deallocate(wtd)
    allocate(wtd(tile_bbi%ncol, tile_bbi%nrow))
    wtd = mvhead
    !
    if (.not.found) then
      call logmsg('***** No nodes found for tile '//ta([itile])//'! *****')
    else
      !
      n = node_max - node_min + 1
      if (allocated(bbi_arr)) deallocate(bbi_arr)
      allocate(bbi_arr(n))
      do i = 1, n
        call bbi_arr(i)%init()
      end do
      max_nc = 0; max_nr = 0
      do ir = 1, nr; do ic = 1, nc
         if (node(ic,ir) /= mvnode) then
           i = node(ic,ir) - node_min + 1
           bbi_arr(i)%ir0 = min(bbi_arr(i)%ir0,ir)
           bbi_arr(i)%ir1 = max(bbi_arr(i)%ir1,ir)
           bbi_arr(i)%ic0 = min(bbi_arr(i)%ic0,ic)
           bbi_arr(i)%ic1 = max(bbi_arr(i)%ic1,ic)
           bbi_arr(i)%ncol = bbi_arr(i)%ic1 - bbi_arr(i)%ic0 + 1
           bbi_arr(i)%nrow = bbi_arr(i)%ir1 - bbi_arr(i)%ir0 + 1
           max_nc = max(max_nc, bbi_arr(i)%ncol)
           max_nr = max(max_nr, bbi_arr(i)%nrow)
         end if
      end do; end do
      if (max_nc /= max_nr) then
        call errmsg('Invalid block size.')
      end if
      !
      cs_min = tile_bbx%cs
      !
      nlev = int(log(real(max_nc,R8B))/log(2.d0)) + 1
      if (allocated(cs_read)) deallocate(cs_read)
      allocate(cs_read(nlev))
      !
      do ilev = 1, nlev
        cs_read(ilev) = cs_min * 2**(ilev-1)
      end do
      !
      call mg_ahn%clean()
      !
      do ilev = 1, nlev
        tile_bbx%cs = cs_read(ilev)
        !
        allocate(hdr_ahn)
        call hdr_ahn%read_extent(f_in_flt, tile_bbx, i_uscl_arith, i_dscl_nointp)
        call hdr_ahn%get_grid(xr4=ahn, mvr4=mvahn)
        call hdr_ahn%clean(); deallocate(hdr_ahn); hdr_ahn => null()
        !
        if (.false.) then
          f = 'e:\data\lhm-flex\snellius\200m_var_01_merge\ahn_l'//ta([ilev])
          call writeflt(f, ahn, size(ahn,1), size(ahn,2), &
            tile_bbx%xll, tile_bbx%yll, tile_bbx%cs, mvahn)
        end if
        !
        if (ilev == 1) then
          call mg_ahn%init(nlev, tile_bbx%xll, tile_bbx%xur, &
            tile_bbx%yll, tile_bbx%yur, cs_read, mvr4=mvahn)
        end if
        g_ahn => mg_ahn%grid(ilev); call g_ahn%clean_xi()
        call g_ahn%set_arr(xr4=ahn)
      end do
      !
      ! set the water table depth
      do ir = 1, nr; do ic = 1, nc
        n = node(ic,ir)
        if (n /= mvnode) then
          i = n - node_min + 1
          bs = bbi_arr(i)%ncol
          ilev = int(log(real(bs,R8B))/log(2.d0)) + 1
          g_ahn => mg_ahn%grid(ilev)
          call  get_xy(xc, yc, ic, ir, tile_bbx%xll, tile_bbx%yur, cs_min)
          call get_icr(jc, jr, xc, yc, tile_bbx%xll, tile_bbx%yur, cs_read(ilev))
          if (valid_icir(jc, jr, g_ahn%nc, g_ahn%nr)) then
            head_val = head(ic,ir)
            ahn_val  = g_ahn%xr4(jc,jr)
            if ((head_val /= mvhead).and.(ahn_val /= mvahn)) then
              wtd(ic,ir) = ahn_val - head_val
            end if
          end if
        end if
      end do; end do
    end if
    !
    ! overwrite the head vrt, and write the wtd flt
    f = vrt_head%raw%raw_mid(itile)%raw(i_SourceFilename)%s(4)
    f = trim(strip_ext(f))//trim(post_fix)//'.flt'
    vrt_head%raw%raw_mid(itile)%raw(i_SourceFilename)%s(4) = f
    call writeflt(strip_ext(f), wtd, tile_bbi%ncol, tile_bbi%nrow, &
      tile_bbx%xll, tile_bbx%yll, cs_min, mvhead)
  end do
  !
  ! write the vrt
  vrt_head%f = f_out_vrt
  call vrt_head%write()
  !
  ! clean-up
  call vrt_head%clean(); deallocate(vrt_head); vrt_head => null()
  call vrt_node%clean(); deallocate(vrt_node); vrt_node => null()
  if (allocated(head)) deallocate(head)
  if (allocated(node)) deallocate(node)
  if (allocated(wtd)) deallocate(wtd)
  !
  return
end subroutine

end module main_module

! ==============================================================================
! ==============================================================================
program quad2d
  use main_module
  !
  ! start timing
  call date_and_time(values=ibdt)
  !
  call quad_settings()
  !
  if (run_opt == 0) then
    call quad_filter()
  end if
  if (run_opt == 9) then
    call quad_fill_gap()
  end if
  if (run_opt == 10) then
    call quad_unique()
  end if
  if (run_opt == 11) then
    call quad_sub_grid_part()
  end if
  if (run_opt == 21) then
    call quad_check_sub_grid_part()
  end if
  if (run_opt == 12) then
    call quad_full_grid_part()
  end if
  if (run_opt == 13) then
    call quad_grid_balancing()
  end if
  if (run_opt == 14) then
    call quad_csv_add_field()
  end if
  if (run_opt == 15) then
    call quad_idf_to_flt()
  end if
  if (run_opt == 16) then
    call quad_merge_csv()
  end if
  !
  if (run_opt == 1) then
    call quad_init()
    call quad_graph()
    if (lsplit) call quad_split()
    if (ljoin)  call quad_join()
    if (lsplit.or.ljoin) call quad_graph()
    call quad_write()
  end if
  !
  ! balancing
  if (run_opt == 2) then
    call quad_clean()
    call quad_read('.intf.bin')
    call quad_balancing()
  end if
  !
  ! interface_gen
  if (run_opt == 3) then
    call quad_clean()
    call quad_read('.intf.bin')
    call quad_init_models()
    call quad_layer_models_interface()
  end if
  !
  ! grid_gen
  if (run_opt == 4) then
    call quad_clean()
    call quad_read('.intf.lm.bin', read_vintf=.true.)
    if (lwrite_disu) then
      call quad_set_mod_dir()
    end if
    call quad_init_models()
    call quad_grid_gen()
  end if
  !
  ! grid_gen_merge
  if (run_opt == 17) then
    call quad_clean()
    call quad_read('.intf.lm.bin', read_vintf=.true.)
    call quad_init_models()
    call quad_grid_gen_merge()
  end if
  !
  ! mf6_xch_write
  if (run_opt == 5) then
    call quad_clean()
    call quad_read(fp_intf='.intf.lm.bin', read_vintf=.false.)
    call quad_mf6_xch_write()
  end if
  !
  ! mf6_xch_write_merge
  if (run_opt == 19) then
    call quad_clean()
    call quad_read(fp_intf='.intf.lm.bin', read_vintf=.false.)
    call quad_mf6_xch_write_merge()
  end if
  !
  ! mf6_data_write
  if (run_opt == 6) then
    call quad_clean()
    call quad_read()
    call quad_set_mod_dir()
    call quad_init_models(set_dat_mod=.true.)
    !call quad_mf6_data_write_old()
    call quad_mf6_data_write()
  end if
  !
  ! mf6_data_write_merge
  if (run_opt == 18) then
    call quad_clean()
    call quad_read()
    call quad_init_models(set_dat_mod=.true.)
    call quad_mf6_data_write_merge()
  end if
  !
  ! partition
  if (run_opt == 7) then
    call quad_partition()
  end if
  !
  ! mf6_post
  if (run_opt == 8) then
    call quad_clean()
    call quad_read('.intf.lm.bin', read_vintf=.false.)
    if (write_chd) then
      call quad_set_mod_dir()
    end if
    call quad_init_models(set_dat_mod=.false.)
    if (write_chd) then
      call quad_mf6_write_chd()
    else
      if (len_trim(f_in_csv) > 0) then
        call quad_mf6_write_heads(lmerge=.true.)
      else
        call quad_mf6_write_heads(lmerge=.false.)
      end if
    end if
  end if
  !
  ! mf6_post_wtd
  if (run_opt == 20) then
    call quad_mf6_write_wtd()
  end if
  !
  ! clean up
  call quad_clean()
  !
  ! write the timing information
  call date_and_time(values=iedt)
  call get_elapsed_time(ibdt, iedt, elsec, elapsed_line)
  call logmsg('')
  call logmsg(trim(elapsed_line))
  call logmsg(' Total number of seconds: '// adjustl(ta([elsec],'(f12.3)')))
  !
end program quad2d

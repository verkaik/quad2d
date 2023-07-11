module main_module
  ! -- modules
  use utilsmod, only: MXSLEN, MXSLENLONG, I1B, I4B, I8B, R4B, R8B, tBb, tBbX, tBbObj, ta, logmsg, &
    errmsg, I4ZERO, I4MINONE, I8ZERO, I8ONE, DHALF, quicksort_r, bbi_intersect, node_to_icrl, &
    open_file, get_xy, get_icr, R8ONE,R8ZERO, get_args, tIni, tCSV, fill_with_nearest, tUnp, &
    calc_unique, change_case, get_compiler, split_str, get_ext, read_line, parse_line, fileexist, fillgap, &
    get_unique, grid_load_imbalance, readidf, get_dir_files, strip_ext, get_slash, create_dir, &
    I_I1, I_I2, I_I4, I_I8, I_R4, I_R8, I_C
  use vrt_module, only: tVrt
  !
  use hdrModule, only: tHdr, tHdrHdr, writeflt, &
    i_dscl_nointp, i_dscl_intp
  !
  use quad2dModule, only: tQuads, tQuad, tIntf, tNbrIntf, &
    tLayerModels, tLayerModel, tDataModels, tDataModel, &
    tProps, get_number_of_levels, get_refinement_level, &
    i_lid, i_gid, i_i_graph, i_lay_mod,  &
    i_tgt_cs_min, i_tgt_cs_min_lay_beg, i_tgt_cs_min_lay_end, i_tgt_cs_min_dz_top, &
    i_tgt_cs_max, i_head, n_prop_field
  !
  use metis_module, only: tMetis
  
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
  character(len=MXSLEN), dimension(:), allocatable :: args
  character(len=MXSLEN) :: f_vrt, f_tile, d_out, f_in_csv, f_out_csv, f_part_csv
  character(len=MXSLEN) :: f_mod_def_inp, f_lay_mod_csv, fp
  character(len=MXSLEN) :: d_in, uuid_in
  character(len=MXSLEN) :: f_gid_in, f_gid_out, f_gid_mask
  character(len=MXSLEN) :: mod_root_dir, xch_root_dir, mod_sub_dir_fields, xch_id_field
  character(len=MXSLEN) :: sel_field, sel_val, sel_npart, sel_nodes
  character(len=MXSLEN) :: chd_lid
  character(len=MXSLEN) :: f_hiera_in, f_hiera_field, f_weight
  character(len=MXSLEN) :: gid_exclude, gid_separate
  character(len=MXSLEN) :: csv_field, csv_val
  character(len=MXSLEN) :: f_in_idf
  character(len=MXSLEN) :: d_log
  !
  ! fields
  character(len=MXSLEN), dimension(n_prop_field) :: fields
  !
  logical :: lrenumber, lwrite, ljoin, lsplit, lremove, lwrite_props, loverwrite_props, lwrite_asc, luse_uuid, lwrite_disu
  logical, parameter :: LDUM = .true.
  logical :: write_nod_map, write_chd, write_hiera
  
  integer(I4B) :: mask_mv, xid_mv
  integer(I4B), dimension(:,:), allocatable :: mask
  !
  integer(I4B) :: nblk, nid, nlid, gid_max, ib, lid, lid_min, lid_max, ljd, gid, nlev, n, m, ntry
  integer(I4B) :: icb, irb, ic, ir, ic0, ic1, ir0, ir1, gnr, gnc, nr, nc, gic, gir, refi4
  integer(I4B) :: icb0, icb1, irb0, irb1, jc, jr, kc, kr, jc0, jc1, jr0, jr1, ntile, it, it0, it1, iu, ilev
  integer(I4B) :: mvxid, gid_mv, area_min, area_max, n_inact, area, i, idum, ireg, jreg, nreg
  integer(I4B) :: run_opt, nzlev_max, gid_mask
  integer(I4B) :: npart
  integer(I4B), dimension(:), allocatable :: gids, l2gid, g2lid
  integer(I4B), dimension(:,:), allocatable :: xid, i4w2d, regun
  !
  real(R4B), dimension(:,:), allocatable :: r4w2d
  !
  real(R4B) :: r4mv
  !
  real(R8B) :: x0, x1, y0, y1
  real(R8B) :: refr8, xll, yll, yur, cs_gid, cs_max
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
  call logmsg('==================================')
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
  case('grid_gen')
    !=========!
    run_opt = 4
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
    call ini%get_val(sect, 'tgt_cs_min_dz_top_field', cv=fields(i_tgt_cs_min_dz_top), cv_def='')
    !
    call ini%get_val(sect, 'mod_root_dir', cv=mod_root_dir)
    call ini%get_val(sect, 'mod_sub_dir_fields', cv=mod_sub_dir_fields, cv_def='gid')
    !
    call ini%get_val(sect, 'write_disu', l4v=lwrite_disu, l4v_def=.true.)

  case('mf6_xch_write')
    !=========!
    run_opt = 5
    !=========!
    call ini%get_val(sect, 'd_in',    cv=d_in)
    call ini%get_val(sect, 'uuid_in', cv=uuid_in, cv_def='quad2d')
    !
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv)
    !
    call ini%get_val(sect, 'xch_root_dir', cv=xch_root_dir)
    call ini%get_val(sect, 'xch_id_field', cv=xch_id_field, cv_def='lid')
    
    call ini%get_val(sect, 'lid_field',       cv=fields(i_lid),     cv_def='lid')
    call ini%get_val(sect, 'gid_field',       cv=fields(i_gid),     cv_def='gid')
    call ini%get_val(sect, 'i_i_graph_field', cv=fields(i_i_graph), cv_def='i_i_graph')
    !
    call ini%get_val(sect, 'lay_mod_field',            cv=fields(i_lay_mod))
  case('mf6_data_write')
    !=========!
    run_opt = 6
    !=========!
    call ini%get_val(sect, 'd_in',    cv=d_in)
    call ini%get_val(sect, 'd_log',   cv=d_log, cv_def='')
    call ini%get_val(sect, 'uuid_in', cv=uuid_in, cv_def='quad2d')
    !
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    !
    call ini%get_val(sect, 'f_mod_def_inp', cv=f_mod_def_inp)
    !
    call ini%get_val(sect, 'lid_field',        cv=fields(i_lid),        cv_def='lid')
    call ini%get_val(sect, 'gid_field',        cv=fields(i_gid),        cv_def='gid')
    call ini%get_val(sect, 'tgt_cs_min_field', cv=fields(i_tgt_cs_min), cv_def='tgt_cs_min')
    !
    call ini%get_val(sect, 'lay_mod_field',            cv=fields(i_lay_mod))
    
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
  case('hiera_grid_part')
    !=========!
    run_opt = 11
    !=========!
    call ini%get_val(sect, 'f_gid_in',      cv=f_gid_in)
    call ini%get_val(sect, 'f_hiera_in',    cv=f_hiera_in)
    call ini%get_val(sect, 'npart',        i4v=npart)
    call ini%get_val(sect, 'f_gid_out',     cv=f_gid_out)
    call ini%get_val(sect, 'f_hiera_field', cv=f_hiera_field)
    call ini%get_val(sect, 'f_out_csv',     cv=f_out_csv)
    call ini%get_val(sect, 'write_hiera',   l4v=write_hiera, l4v_def=.false.)
  case('full_grid_part')
    !=========!
    run_opt = 12
    !=========!
    call ini%get_val(sect, 'f_gid_mask',  cv=f_gid_mask)
    call ini%get_val(sect, 'f_weight',    cv=f_weight, cv_def='')
    call ini%get_val(sect, 'npart',       i4v=npart)
    call ini%get_val(sect, 'f_gid_out',   cv=f_gid_out)
    call ini%get_val(sect, 'gid_exclude', cv=gid_exclude, cv_def='')
    call ini%get_val(sect, 'gid_separate', cv=gid_separate, cv_def='')
  case('lump_grid_part')
    !=========!
    run_opt = 13
    !=========!
    call ini%get_val(sect, 'f_gid_in',    cv=f_gid_in)
    call ini%get_val(sect, 'f_weight',    cv=f_weight, cv_def='')
    call ini%get_val(sect, 'npart',       i4v=npart)
    call ini%get_val(sect, 'f_gid_out',   cv=f_gid_out)
  case('csv_add_field')
    !=========!
    run_opt = 14
    !=========!
    call ini%get_val(sect, 'f_in_csv', cv=f_in_csv)
    call ini%get_val(sect, 'f_out_csv', cv=f_out_csv)
    call ini%get_val(sect, 'csv_field', cv=csv_field)
    call ini%get_val(sect, 'csv_val', cv=csv_val)
  case('idf_to_flt')
    !=========!
    run_opt = 15
    !=========!
    call ini%get_val(sect, 'f_in_idf', cv=f_in_idf)
    call ini%get_val(sect, 'mv', r4v=r4mv)
  case default
    call errmsg('Invalid run option: '//trim(sect))
  end select
  !
  call logmsg('***** Run option: '//trim(sect)//' *****')
  call logmsg('==================================')
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
  integer(I4B) :: weight_mv, id, npart_loc, npart_tot, ip, ip_offset, nid
  real(R4B) :: wtot, wloc
  real(R8B) :: imbal
  integer(I4B), dimension(:), allocatable :: ids
  integer(I4B), dimension(:,:), allocatable :: weight, weight_loc
! ------------------------------------------------------------------------------
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
        ip_offset = ip_offset + npart_loc
        npart_tot = max(0, npart_tot - npart_loc)
        deallocate(weight_loc)
      end if
    end do
    !
    xid = abs(xid)
    call grid_load_imbalance(xid, xid_mv, weight, imbal, npart_loc)
    call logmsg('Overall load imbalance for '//ta([npart_loc])//' parts: '//ta([imbal]))
  else
    allocate(met)
    call met%init(weight, npart)
    call met%set_opts(niter_in=100)
    call met%recur()
    call met%set_ids(xid)
    call met%clean(); deallocate(met); met => null()
  end if
  
  ! write the new ids
  call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
  call hdrg_gid%write(f_gid_out)
  
  ! clean-up
  if (allocated(xid)) deallocate(xid)
  if (associated(hdrg_gid)) then
    call hdrg_gid%clean(); deallocate(hdrg_gid); hdrg_gid => null()
  end if
  !
  return
end subroutine quad_full_grid_part

subroutine quad_hiera_grid_part()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- local
  type(tQuad), pointer :: q_new => null()
  type(tBb)  :: bbi_gid, bb_new
  type(tBbX) :: bbx_gid, bbx_hid, bbx_q
  type(tHdr), pointer :: hdrg_hid => null()
  type(tHdrHdr), pointer :: hdr => null()
  type(tCSV), pointer    :: csv => null()
  type(tMetis), pointer :: met => null()
  !
  logical :: lfound, lfirst, llast
  character(len=MXSLEN), dimension(:), allocatable :: sa
  integer(I4B) :: ngid, gid_max_loc, nh, ih, xidh_mv, id, area_n, nhid, j
  integer(I4B) :: np_rem, np, np_full, np_lump, ip, p, lid_new, gid_new, iact
  integer(I4B), dimension(:,:), allocatable :: mask, part, xidh, hmap
  integer(I4B), dimension(:), allocatable :: lev_id, uplev_id, hid, hlev
  real(R8B), dimension(:), allocatable :: areah
  real(R8B) :: area, area_metis, area_q, xm, ym, area_tot, area_tgt
! ------------------------------------------------------------------------------
  !
  allocate(hdrg_gid)
  call hdrg_gid%read_full_grid(f_gid_in)
  call hdrg_gid%get_grid(xi4=xid, mvi4=xid_mv)
  nc =  size(xid,1); nr = size(xid,2)
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
  lid = 0; area_tot = R8ZERO
  do gid = 1, gid_max
    if (gids(gid) == 1) then
      lid = lid + 1
      g2lid(lid) = gid
      q => xq%get_quad(lid)
      bb => bb_gid(gid) ! local bounding box
      bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb%ir1)*cs_gid
      bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
      bbx%cs = cs_gid
      call bbo%set(child_bbi=bb, child_bbx=bbx)
      call q%init(gid=gid, lid=lid, bbo=bbo)
      call get_mask(xid, gid, bb, mask)
      call q%calc_prop(mask=mask)
      call q%get_prop(area=area_n)
      area_tot = area_tot + area_n*cs_gid*cs_gid
    end if
  end do
  !
  ! create the hierarchy
  call parse_line(f_hiera_in, sa, token_in=',')
  nh = size(sa)
  allocate(hmap(nh,nlid), hlev(nh)); hmap = 0
  !
  do ih = 1, nh
    allocate(hdrg_hid)
    call hdrg_hid%read_full_grid(sa(ih))
    call hdrg_hid%get_grid(xi4=xidh, mvi4=xidh_mv)
    hdr => hdrg_hid%hdr; call hdr%get_bbx(bbx_hid)
    !
    do lid = 1, xq%n
      q => xq%get_quad(lid)
      call q%get_prop(xm=xm, ym=ym)
      call get_icr(ic, ir, xm, ym, bbx_hid%xll, bbx_hid%yur, bbx_hid%cs)
      if ((ic > 0).and.(ir> 0)) then
        id = xidh(ic,ir)
        if (id /= xidh_mv) then
          hmap(ih,lid) = id
        end if
      else
        call errmsg('quad_hiera_grid_part: could not sample point.')
      end if
    end do
    call hdrg_hid%clean(); deallocate(hdrg_hid); hdrg_hid => null()
  end do
  !
  ! check
  if (minval(hmap(1,:)) == 0) then
    call errmsg('quad_hiera_grid_part: not all areas are classified.')
  end if
  !
  ! set the levels
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    call q%set(hlev=hmap(:,lid))
  end do
  !
  area_tgt = area_tot / npart
  np_rem = npart
  !area_tgt = area_tot / 5
  !
  do ih = 1, nh
    ! get the IDs
    call get_unique(hmap(ih,:), hid)
    nhid = size(hid)
    allocate(areah(nhid)); areah = R8ZERO
    !
    do i = 1, nhid
      id = hid(i)
      if (id == 0) cycle
      do lid = 1, xq%n
        q => xq%get_quad(lid)
        if (q%get_flag(active=LDUM)) then
          if (q%get_hlev(ih) == id) then
            call q%get_prop(area=area_n)
            areah(i) = areah(i) + area_n*cs_gid*cs_gid
          end if
        end if
      end do
    end do
    !
    gid_new = gid_max_loc
    !
    ! check the target area
    do i = 1, nhid
      id = hid(i)
      if (areah(i) == R8ZERO) cycle
      if (areah(i) < area_tgt) then ! merge the quads for this level
        lid_new = xq%n + 1; gid_new = gid_new + 1
        if (gid_new == 25134) then
          write(*,*) '@@@@'
        end if
        call bb_new%init()
        ! merge the quads
        lfound = .false.; lfirst = .true.
        do lid = 1, xq%n
          q => xq%get_quad(lid)
          if (q%get_flag(active=LDUM)) then
            if (q%get_hlev(ih) == id) then
              lfound = .true.
              call q%set_flag(active=.false.)
              call q%get_bb(child_bbi=bb)
              ! change the id
              do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
                if (xid(ic,ir) == q%gid) then
                  xid(ic,ir) = gid_new
                end if
              end do; end do
              !
              ! determine the new bounding box
              bb_new%ic0 = min(bb_new%ic0,bb%ic0); bb_new%ic1 = max(bb_new%ic1,bb%ic1)
              bb_new%ir0 = min(bb_new%ir0,bb%ir0); bb_new%ir1 = max(bb_new%ir1,bb%ir1)
              bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1;  bb_new%nrow = bb_new%ir1 - bb_new%ir0 + 1
              !
              if (lfirst) then
                hlev = 0; hlev(1:ih) = q%hlev(1:ih)
                lfirst = .false.
              end if
            end if
          end if
        end do
        if (lfound) then
          bbx%xll = xll + (bb%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb%ir1)*cs_gid
          bbx%xur = bbx%xll + bb%ncol*cs_gid; bbx%yur = bbx%yll + bb%nrow*cs_gid
          bbx%cs = cs_gid
          call bbo%set(child_bbi=bb_new, child_bbx=bbx)
          !
          q_new => xq%get_quad(lid_new)
          call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
          call q_new%set(hlev=hlev)
          call q_new%get_bb(child_bbi=bb)
          call get_mask(xid, q_new%gid, bb, mask)
          call q_new%calc_prop(mask=mask)
          xq%n = lid_new
        end if
      else
        llast = .false.
        if (ih == nh) then
          llast = .true.
        else
          llast = .true.
          do lid = 1, xq%n
            q => xq%get_quad(lid)
            if (q%get_flag(active=LDUM)) then
              if (q%get_hlev(ih) == id) then
                if (q%get_hlev(ih+1) /= 0) then
                  llast = .false.
                end if
              end if
            end if
          end do
        end if
        !
        if (llast) then
          ! step 1: apply full METIS for quads with area > area_tgt
          area_metis = R8ZERO
          n = xq%n; lid_new = n
          do lid = 1, n
            q => xq%get_quad(lid)
            if (q%get_flag(active=LDUM)) then
              if (q%get_hlev(ih) == id) then
                call q%get_prop(area=area_n)
                area_q = area_n*cs_gid*cs_gid
                if (area_q > area_tgt) then
                  ! number of parts
                  np_full = max(2,nint(area_q/area_tgt))
                  !
                  call q%get_bb(child_bbi=bb)
                  call get_mask(xid, q%gid, bb, mask)
                  call q%get_prop(area=area_n) !DEBUG
                  area = area_n*cs_gid*cs_gid !DEBUG
                  !
                  ! full grid METIS
                  allocate(met)
                  call met%init(mask, np_full)
                  call met%set_opts()
                  call met%recur(verbose=.true.)
                  if (allocated(part)) deallocate(part)
                  allocate(part,source=mask)
                  call met%set_ids(part)
                  !write(*,*) '@@',maxval(part)
                  call met%clean(); deallocate(met); met => null()
                  !
                  do ip = 1, np_full
                    lid_new = lid_new + 1; gid_new = gid_new + 1
                    !if (gid_new == 25134) then
                    !  write(*,*) '@@@@'
                    !end if
                    q_new => xq%get_quad(lid_new)
                    call bb_new%init()
                    do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
                      jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
                      p = part(jc,jr)
                      if (p == ip) then
                       ! determine the new bounding box
                        bb_new%ic0 = min(bb_new%ic0,ic); bb_new%ic1 = max(bb_new%ic1,ic)
                        bb_new%ir0 = min(bb_new%ir0,ir); bb_new%ir1 = max(bb_new%ir1,ir)
                        bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1;  bb_new%nrow = bb_new%ir1 - bb%ir0 + 1
                        if (xid(ic,ir) /= q%gid) then
                          call errmsg('quad_hiera_grid_part: program error.')
                        else
                          xid(ic,ir) = gid_new
                        end if
                      end if 
                    end do; end do
                    !
                    bbx%xll = xll + (bb_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb_new%ir1)*cs_gid
                    bbx%xur = bbx%xll + bb_new%ncol*cs_gid; bbx%yur = bbx%yll + bb_new%nrow*cs_gid
                    bbx%cs = cs_gid
                    call bbo%set(child_bbi=bb_new, child_bbx=bbx)
                    call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
                    call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid)
                    call get_mask(xid, gid_new, bb_new, mask)
                    call q_new%calc_prop(mask=mask)
                    call q_new%get_prop(area=area_n)
                    area_metis = area_metis + area_n*cs_gid*cs_gid
                    xq%n = lid_new
                  end do !ip
                  !
                  call q%set_flag(active=.false.)
                end if
              end if
            end if
          end do
          !
          ! set: apply lumped METIS to the other 
          area = areah(i) - area_metis
          if (area < R8ZERO) then 
            call errmsg('quad_hiera_grid_part: program error')
          end if
          np_lump = 0
          if (area > area_tgt) then
            np_lump = max(2,nint(area/area_tgt))
            !
            ! determine the bounding box and work array
            call bb_new%init(); lfirst = .true.
            do iact = 1, 2
              do lid = 1, xq%n
                q => xq%get_quad(lid)
                if ((q%get_flag(active=LDUM)).and.(q%gid_prent == 0)) then
                  if (q%get_hlev(ih) == id) then
                    call q%get_bb(child_bbi=bb)
                    if (iact == 1) then
                      bb_new%ic0 = min(bb_new%ic0,bb%ic0); bb_new%ic1 = max(bb_new%ic1,bb%ic1)
                      bb_new%ir0 = min(bb_new%ir0,bb%ir0); bb_new%ir1 = max(bb_new%ir1,bb%ir1)
                      bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1
                      bb_new%nrow = bb_new%ir1 - bb_new%ir0 + 1
                    else
                      do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
                        if (xid(ic,ir) == q%gid) then
                          jr = ir - bb_new%ir0 + 1; jc = ic - bb_new%ic0 + 1
                          i4w2d(jc,jr) = q%lid
                        end if
                      end do; end do
                    end if
                    if (lfirst) then
                      hlev = 0; hlev(1:ih) = q%hlev(1:ih)
                      lfirst = .false.
                    end if
                  end if
                end if
              end do
              if (iact == 1) then
                if (allocated(i4w2d)) deallocate(i4w2d)
                allocate(i4w2d(bb_new%ncol,bb_new%nrow))
                i4w2d = 0
              end if
            end do
            !
            allocate(met)
            call met%init_lump(i4w2d, np_lump, verbose=.true.)
            call met%set_opts()
            call met%recur(verbose=.true.)
            if (allocated(part)) deallocate(part)
            !
            do ip = 1, np_lump
              lid_new = lid_new + 1; gid_new = gid_new + 1
              q_new => xq%get_quad(lid_new)
              call bb_new%init()
              do j = 1, met%nvtxs
                if (met%part(j) == ip - 1) then
                  lid = met%idmapinv(j)
                  q => xq%get_quad(lid)
                  if ((q%get_hlev(ih) /= id).or.(.not.q%get_flag(active=LDUM))) then
                    call errmsg('quad_hiera_grid_part: program error')
                  end if
                  call q%set_flag(active=.false.)
                  call q%get_bb(child_bbi=bb)
                  bb_new%ic0 = min(bb_new%ic0,bb%ic0); bb_new%ic1 = max(bb_new%ic1,bb%ic1)
                  bb_new%ir0 = min(bb_new%ir0,bb%ir0); bb_new%ir1 = max(bb_new%ir1,bb%ir1)
                  bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1
                  bb_new%nrow = bb_new%ir1 - bb_new%ir0 + 1
                  !
                  do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
                    if (xid(ic,ir) == q%gid) then
                      xid(ic,ir) = gid_new
                    end if
                  end do; end do
                end if
              end do
              bbx%xll = xll + (bb_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb_new%ir1)*cs_gid
              bbx%xur = bbx%xll + bb_new%ncol*cs_gid; bbx%yur = bbx%yll + bb_new%nrow*cs_gid
              bbx%cs = cs_gid
              call bbo%set(child_bbi=bb_new, child_bbx=bbx)
              call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo) !TODO: set parent LIST
              call q_new%set(hlev=hlev) 
              call q_new%get_bb(child_bbi=bb)
              call get_mask(xid, q_new%gid, bb, mask)
              call q_new%calc_prop(mask=mask)
              xq%n = lid_new
            end do
            call met%clean(); deallocate(met); met => null()
          end if
        end if ! ih = nh
      end if
    end do
    !
    deallocate(hid, areah)
  end do
  !
  n = xq%get_number_active()
  call logmsg('# unique IDs BEFORE splitting: '//ta([n]))
  !
  ! check
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) then
      call q%get_bb(child_bbi=bb)
      call get_mask(xid, q%gid, bb, mask)
    end if
  end do
  !
  ! split non-unique parts
  n = xq%n; lid_new = n
  do lid = 1, n
    if (lid == 19775) then
      write(*,*) '@@@@'
    end if
    q => xq%get_quad(lid)
    if (q%get_flag(active=LDUM)) then
      call q%get_bb(child_bbi=bb)
      call get_mask(xid, q%gid, bb, i4w2d)
      call calc_unique(i4w2d, 5, regun, regbb, nreg, idum, 0., 0., 0.)
      if (nreg > 1) then
        do ireg = 1, nreg
          lid_new = lid_new + 1; gid_new = gid_new + 1
          call bb_new%init()
          do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
            jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
            if (regun(jc,jr) == ireg) then
              bb_new%ic0 = min(bb_new%ic0,bb%ic0); bb_new%ic1 = max(bb_new%ic1,bb%ic1)
              bb_new%ir0 = min(bb_new%ir0,bb%ir0); bb_new%ir1 = max(bb_new%ir1,bb%ir1)
              bb_new%ncol = bb_new%ic1 - bb_new%ic0 + 1
              bb_new%nrow = bb_new%ir1 - bb_new%ir0 + 1
              if ((ic==314).and.(ir==993).and.(gid_new==25134)) then
                write(*,*) '@@@@'
              end if
              xid(ic,ir) = gid_new
            end if
          end do; end do
          bbx%xll = xll + (bb_new%ic0-1)*cs_gid; bbx%yll = yll + (nr-bb_new%ir1)*cs_gid
          bbx%xur = bbx%xll + bb_new%ncol*cs_gid; bbx%yur = bbx%yll + bb_new%nrow*cs_gid
          bbx%cs = cs_gid
          call bbo%set(child_bbi=bb_new, child_bbx=bbx)
          call q_new%init(gid=gid_new, lid=lid_new, bbo=bbo)
          call q_new%set(hlev=q%hlev, gid_prent=q%gid, lid_prent=q%lid) 
          call q_new%get_bb(child_bbi=bb)
          call get_mask(xid, q_new%gid, bb, mask)
          call q_new%calc_prop(mask=mask)
          xq%n = lid_new
        end do
        call q%set_flag(active=.false.)
      end if
    end if
  end do
  !
  n = xq%get_number_active()
  call logmsg('# unique IDs AFTER splitting: '//ta([n]))
  
  ! write the new ids
  call hdrg_gid%replace_grid(xi4=xid, mvi4=xid_mv)
  call hdrg_gid%write(f_gid_out)
  !
  ! write the hierarchical grids
  if (write_hiera) then
    do ih = 1, nh
      ! get the IDs
      call get_unique(hmap(ih,:), hid)
      nhid = size(hid)
      !
      if (allocated(i4w2d)) deallocate(i4w2d)
      allocate(i4w2d(nc,nr)); i4w2d = 0
      !
      do i = 1, nhid
        id = hid(i)
        if (id == 0) cycle
        do lid = 1, xq%n
          q => xq%get_quad(lid)
          if (q%get_flag(active=LDUM)) then
            if (q%get_hlev(ih) == id) then
              call q%get_bb(child_bbi=bb)
              call get_mask(xid, q%gid, bb, mask)
              do ir = bb%ir0, bb%ir1; do ic = bb%ic0, bb%ic1
                jr = ir - bb%ir0 + 1; jc = ic - bb%ic0 + 1
                if (mask(jc,jr) > 0) then
                  i4w2d(ic,ir) = id
                end if
              end do; end do
            end if
          end if
        end do
      end do
      !
      fp = trim(sa(ih))//'_calc'
      call writeflt(fp, i4w2d, nc, nr, xll, yll, cs_gid, 0)
    end do
  end if
  
  ! write the csv
  !allocate(csv)
  !call csv%init(file=f_out_csv, &
  !    hdr_keys=['lid', 'gid', 'lid_prent', 'gid_prent', &
  !      'xm', 'ym', 'area'],&
  !    nr=xq%n, hdr_i_type=[i_i4, i_i4, i_i4, i_i4, &
  !      i_r8, i_r8, i_r8,])
  
  ! clean-up
  if (allocated(hmap)) deallocate(hmap)
  if (allocated(hlev)) deallocate(hlev)
  if (allocated(xid)) deallocate(xid)
  if (allocated(i4w2d)) deallocate(i4w2d)
  if (associated(hdrg_gid)) then
    call hdrg_gid%clean(); deallocate(hdrg_gid); hdrg_gid => null()
  end if

  !
  return
end subroutine quad_hiera_grid_part

subroutine quad_csv_add_field()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- local
  type(tCSV), pointer    :: csv => null()
! -- local
  character(len=MXSLEN), dimension(:), allocatable :: hdr, hdr_add, val_add
  integer(I4B) :: nc, nc_add, nc_max, ir, ic, jc
! ------------------------------------------------------------------------------
  allocate(csv)
  call csv%read_hdr(f_in_csv, hdr)
  nc = size(hdr)
  !
  call parse_line(csv_field, hdr_add, token_in=',')
  call parse_line(csv_val,   val_add, token_in=',')
  nc_add = size(hdr_add)
  if (nc_add /= size(val_add)) then
    call errmsg('Size csv_field differs from csv_val.')
  end if
  nc_max = nc + nc_add
  !
  call csv%read(f_in_csv, nc_max=nc_max)
  !
  ! add the new data
  call csv%add_hdr(hdr_add)
  do ir = 1, csv%nr
    do jc = 1, nc_add
      ic = nc + jc
      call csv%set_val(ic=ic, ir=ir, cv=val_add(jc))
    end do
  end do
  !
  csv%file = trim(f_out_csv); call csv%write()
  call csv%clean(); deallocate(csv)
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

subroutine get_mask(xid, id, bb, mask)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------  
! -- arguments
  integer(I4B), dimension(:,:), intent(in) :: xid
  integer(I4B), intent(in) :: id
  type(tBB) :: bb
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
  call xq%generate_uuid(luse_uuid)
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
      call dat_mod%init(keys(i), fdat(1), fdat(2))
    end if
  end do
  !
  close(iu)
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

subroutine quad_mf6_data_write()
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
    !if (q%gid /= 500005) cycle
    !if (q%gid /= 636924) cycle
    !if (q%gid /= 990378) cycle
    !if (q%gid /= 603528) cycle
    !if (q%gid /= 603486) cycle
    !if (q%gid /= 600157) cycle
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
    !if (.false.) then
    !  allocate(wbd)
    !  !
    !  ! array
    !  call wbd%init(trim(this%mod_dir)//'\dat.csv', trim(this%mod_dir)//'\dat.binpos')
    !  call wbd%write_array(id='v1', nod_ptr=[0_I1B,1_I1B,1_I1B], i4a=[3, 4, 5])
    !  call wbd%write_array(id='v2', i4a=[3, 4, 5])
    !  call wbd%write_array(id='v3', r8a=[3.d0, 4.d0, 5.d0])
    !  call wbd%write_array(id='v4', nod_ptr=[0_I1B,1_I1B,1_I1B], i4a=[3, 4, 5], f_bin=trim(this%mod_dir)//'\dat.bin')
    !  call wbd%write_array(id='v5', nod_ptr=[0_I1B,1_I1B,1_I1B], r8a=[3.d0, 4.d0, 5.d0])
    !  call wbd%write_array(id='v6', nod_ptr=[0_I1B,1_I1B,1_I1B], r8a=[3.d0, 4.d0, 5.d0], f_asc=trim(this%mod_dir)//'\dat.asc')
    !  !
    !  ! lists:
    !  call wbd%write_list(id='riv', period_beg=1, period_end=10, nod_ptr=[0_I1B,1_I1B,1_I1B], i4a1=[1, 2, 3], &
    !    r8a1=[1.d0, 1.d0, 1.d0],r8a2=[2.d0, 2.d0, 2.d0],r8a3=[3.d0, 3.d0, 3.d0], f_asc=trim(this%mod_dir)//'\riv.asc')
    !  call wbd%write_list(id='riv', period_beg=1, period_end=10, nod_ptr=[0_I1B,1_I1B,1_I1B], i4a1=[1, 2, 3], &
    !    r8a1=[1.d0, 1.d0, 1.d0],r8a2=[2.d0, 2.d0, 2.d0],r8a3=[3.d0, 3.d0, 3.d0], f_bin=trim(this%mod_dir)//'\riv.bin')
    !  call wbd%write_list(id='exg', &
    !    i4a1=[1, 2, 3], i4a2=[4, 5, 6], i4a3=[7, 8, 9], &
    !    r8a1=[1.d0, 2.d0, 3.d0],r8a2=[4.d0, 5.d0, 6.d0],r8a3=[7.d0, 8.d0, 9.d0], f_bin=trim(this%mod_dir)//'\exg.bin')
    !  call wbd%write_list(id='exg', &
    !    i4a1=[1, 2, 3], i4a2=[4, 5, 6], i4a3=[7, 8, 9], &
    !    r8a1=[1.d0, 2.d0, 3.d0],r8a2=[4.d0, 5.d0, 6.d0],r8a3=[7.d0, 8.d0, 9.d0])
    !  
    !  call wbd%write_csv()
    !  stop
    !end if
  
  return
end subroutine quad_mf6_data_write

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
  type(tCsv), pointer :: csv => null(), csv_old => null()
  !
  integer(I4B), parameter :: perc_intv = 10
  !
  character(len=1) :: slash
  character(len=MXSLEN) :: f_csv_dat, s_lay_act
  integer(I4B), dimension(:), allocatable :: lay, lid_map, lid_arr
  integer(I4B) :: ncell_tot, ncell, nja, nlay_act, ngrid, write_csv_delta
  integer(I4B) :: lid0, lid1, n_act, i, n
  logical :: lskip
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
    if (q%get_flag(active=LDUM)) n_act = n_act + 1
  end do
  !
  n = 0
  slash = get_slash()
  do lid = 1, xq%n
    q => xq%get_quad(lid)
    !
    ! DEBUG:
    !if (lid == 24728) then
    !  write(*,*) '@@@@ 24728'
    !end if
    !
    if (q%get_flag(active=LDUM)) then
      call logmsg('***** Processing quad '//ta([q%gid])//' *****')
      !
      lskip = .false.
      if ((lid >= lid0).and.(lid <= lid1)) then
        f_csv_dat = trim(q%mod_dir)//slash//'dat.csv'
        call q%grid_gen(f_csv_dat, lwrite_disu, lwrite_asc, &
          ncell, nja, nlay_act, lay, ngrid, lskip)
        s_lay_act = ta(lay,sep_in=';')
      else
        if (associated(csv_old)) then
          call logmsg('***** Using data from old csv for quad '//ta([q%gid])//' *****')
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('nodes'), i4v=ncell)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('nja'), i4v=nja)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('nlay_act'), i4v=nlay_act)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('ngrid_lev'), i4v=ngrid)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('csv_dat'), cv=f_csv_dat)
          call csv_old%get_val(ir=lid_map(lid), ic=csv_old%get_col('lay_act'), cv=s_lay_act)
        else
          cycle
        end if
      end if
      !
      n = n + 1; ncell_tot = ncell_tot + ncell
      call logmsg('***** '//ta([100.*real(n,R4B)/n_act],'(f6.2)')//' %: '// &
        trim(adjustl(ta([real(ncell_tot,R4B)/1000000.],'(f10.2)')))//' M cells *****')
      if (lskip) then
        call logmsg('********** SKIPPING: NON-EXISTING LAYER MODEL *********')
        call q%set_flag(active=.false.)
      end if
      if (lwrite_props) then
        call q%set_prop_csv(key='nodes', i4v=ncell)
        call q%set_prop_csv(key='nja', i4v=nja)
        call q%set_prop_csv(key='nlay_act', i4v=nlay_act)
        call q%set_prop_csv(key='ngrid_lev', i4v=ngrid)
        call q%set_prop_csv(key='lay_act', cv=trim(s_lay_act))
        if (lwrite_disu) then
          call q%set_prop_csv(key='csv_dat', cv=trim(f_csv_dat))
        end if
        if ((n == xq%n_act).or.(mod(n,write_csv_delta) == 0)) then
          csv%file = trim(f_out_csv)
          call csv%write()
        end if
      end if
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
  call parse_line(s=sel_val, sa=sel_val_arr, token_in=',')
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

subroutine quad_mf6_write_heads()
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- local
  integer(I4B) :: lid0, lid1
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

  call xq%write_mf6_heads(lid0, lid1, head_layers, kper_beg, kper_end, &
    tile_nc, tile_nr, f_vrt, write_nod_map)
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

end module main_module

! ==============================================================================
! ==============================================================================
program quad2d
  use main_module
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
    call quad_hiera_grid_part()
  end if
  if (run_opt == 12) then
    call quad_full_grid_part()
  end if
  if (run_opt == 13) then
    !call quad_lump_grid_part()
  end if
  if (run_opt == 14) then
    call quad_csv_add_field()
  end if
  if (run_opt == 15) then
    call quad_idf_to_flt()
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
  ! mf6_xch_write
  if (run_opt == 5) then
    call quad_clean()
    call quad_read(fp_intf='.intf.lm.bin', read_vintf=.false.)
    call quad_mf6_xch_write()
  end if
  
  ! mf6_data_write
  if (run_opt == 6) then
    call quad_clean()
    call quad_read()
    call quad_set_mod_dir()
    call quad_init_models(set_dat_mod=.true.)
    call quad_mf6_data_write()
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
      call quad_mf6_write_heads()
    end if
  end if
  !
  ! clean up
  call quad_clean()
end program quad2d

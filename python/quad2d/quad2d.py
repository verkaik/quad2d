"""
Author: Jarno Verkaik
Program: quad2d
Date: may 2023
"""

# import the modules
import os
import sys
from jinja2 import FileSystemLoader, Environment, meta
import configparser
import csv
import logging as log
#import logging
import datetime
import time
import argparse
from subprocess import Popen, PIPE
from collections import OrderedDict as od
from pathlib import Path

# get the starting time
start_time = datetime.datetime.now().strftime('%d-%m-%Y %H:%M:%S')

class DeltaTimeFormatter(log.Formatter):
    def format(self, record):
        duration = datetime.datetime.utcfromtimestamp(record.relativeCreated / 1000)
        record.delta = duration.strftime("%H:%M:%S")
        return super().format(record)

# logger
handler_1 = log.FileHandler("quad2d.log")
handler_2 = log.StreamHandler(sys.stdout)
LOGFORMAT = '+%(delta)s %(levelname)s: %(message)s'
fmt = DeltaTimeFormatter(LOGFORMAT)
handler_1.setFormatter(fmt)
handler_2.setFormatter(fmt)
log.basicConfig(handlers=[handler_1,handler_2],
                format=LOGFORMAT,
                level=log.DEBUG)
log.info(f'Start of program: {start_time}')

# command line parser
clp = argparse.ArgumentParser()
clp.add_argument('-mpi', '--mpi', type=str, \
                default=r'c:\Program Files\Microsoft MPI\Bin\mpiexec.exe', \
                help='Message Passing Interface program.')
clp.add_argument('-np', '--np', type=int, default=1,
                help='Number of MPI processes.')
clp.add_argument('-mf6', '--mf6', type=str,
                default=r'c:\Users\verkaik_jo\data\codes\git\modflow6-parallel-fork-18-03-21\bin\mf6.exe',
                help='MODFLOW 6 executable.')
clp.add_argument('-mfsim', '--mfsim', type=str,
                default=r'mfsim.nam',
                help='MODFLOW 6 simulation name file.')
clp.add_argument('-templates', '--templates', type=str, \
                default=r'.\mf6_templates/', \
                help='Template directory for MODFLOW 6 files.')
clp.add_argument('-pre', '--pre', action='store_true', default=False, \
                help='Flag for pre-processing MODFLOW 6')
clp.add_argument('-run', '--run', action='store_true', default=False, \
                help='Flag for running MODFLOW 6.')
clp.add_argument('-post', '--post', action='store_true', default=False, \
                help='Flag for post-processing MODFLOW 6')
clp.add_argument('-coupled', '--coupled', action='store_true', default=True, \
                help='Flag using coupled using exchanges.')
clp.add_argument('-parallel', '--parallel', action='store_true', default=False, \
                help='Flag for enabling parallel computing .')
clp.add_argument('-serial_decoupled', '--serial_decoupled', action='store_true', default=False, \
                help='Flag for enabling running with multiple decoupled models .')
clp.add_argument('-cgc', '--cgc', action='store_true', default=False,
                help='Coarse grid correction option.')
clp.add_argument('-cgc_solver', '--cgc_solver', type=int, default=1,
                help='Coarse grid correction solver (1: LU; 2: ILU(0)).')
clp.add_argument('-ini', '--ini', type=str, help='INI-file.')
cla = clp.parse_args().__dict__
log.info(10*'='+'BEGIN command line arguments'+10*'=')
for key in cla:
    log.info('%s: %s'%(key,str(cla[key])))
log.info(10*'='+'END command line arguments'+10*'=')

i_const  = 0
i_asc    = 1
i_bin    = 2
i_binpos = 3

# {%- by itself means current line should have no empty lines between current and previous line
# -%} by itself means current line should have a single empty line above it
# {%- and -%} means current line should be flush with previous line

#############################################################################
def get_cla_key(key):
#############################################################################
    if key not in cla:
        raise Exception('Key %s not found' % key)
    return cla[key]

#############################################################################
def read_csv(f, filter_ids=[]):
#############################################################################
    log.info('Reading %s...'%f)
    reader = csv.reader(open(f, 'r'))
    i = 0; d = {}
    for row in reader:
        if (i == 0):
            hdr = row[1:]
        else:
            id = row[0]
            #
            if not filter_ids:
                add_id = True
            else:
                if '-' in id:
                   id_list = id.split('-')
                else:
                   id_list = [id]
                #
                add_id = True
                for jd in id_list:
                     if eval(jd) not in filter_ids:
                         add_id = False
            if add_id:
                #print(id)
                #print(row)
                d[id] = {hdr[i]:row[i+1] for i in range(len(hdr))}
        i += 1
    return d

#############################################################################
class mf6_template():
#############################################################################
    def __init__(self, f):
        templates_dir = str(Path(get_cla_key('templates')))
        loader = FileSystemLoader(templates_dir)
        env = Environment(loader=loader, keep_trailing_newline=True)
        self.env = env
        self.f = f+'.j2'
        self.template = env.get_template(self.f)

    def set(self, d, valid_keys=None):
        # first, skip the booleans
        for key in d:
            v = d[key]
            if isinstance(v, str):
                if ((v.lower() == 't') or (v.lower() == 'true')):
                    d[key] = True
                if ((v.lower() == 'f') or (v.lower() == 'false')):
                   d.pop(key)

        if valid_keys:
            for key in d:
                if key not in valid_keys:
                    v = d[key]
                    if (type(v) != type(True)):
                        print(valid_keys)
                        log_error(f'Invalid key found for processing {self.f}: {key}')
        self.d = d

    def render(self, f):
        log.info(f'Writing {f}...')
        f = open(f,'w')
        f.write(self.template.render(self.d))
        f.close()

    def get_variables(self):
        template_source = self.env.loader.get_source(self.env, self.f)[0]
        parsed_content = self.env.parse(template_source)
        return list(meta.find_undeclared_variables(parsed_content))

#############################################################################
def log_error(msg):
#############################################################################
    log.error(msg)
    sys.exit(1)

#############################################################################
def read_ini(f):
#############################################################################
    log.info(f'Reading {f}...')
    config = configparser.ConfigParser()
    config.read(f)
    d1 = config._sections
    d = od()
    for k1 in d1:
        k1_lc = k1.lower()
        d[k1_lc] = od()
        #
        d2 = d1[k1]
        for k2 in d2:
            k2_lc = k2.lower()
            v2 = d2[k2]
            v = v2
            if ('periods' in k2):
                try:
                    v_tmp = eval(v2.replace('\n',''))
                    if not isinstance(v_tmp, dict):
                        log_error(f'Periods in [{k1}] should be a dictionary')
                    v = v_tmp
                    if set(map(type, v)) == {range}: # range as dictionary key
                        v_tmp = {}
                        for k in v.keys():
                            range_list = list(list(k))
                            for i in range_list:
                                if i in v_tmp:
                                    log_error(f'Overlapping data for periods in [{k1}]')
                                v_tmp[i] = v[k]
                    v = od(sorted(v_tmp.items()))
                except:
                    log_error(f'Could not read periods in [{k1}]')
            elif('perioddata' in k2):
                try:
                    v_tmp = eval(v2.replace('\n',''))
                    v = v_tmp
                except:
                    log_error(f'Could not read perioddata in [{k1}]')
            d[k1_lc][k2_lc] = v
    return d

#############################################################################
def check_main_ini(d, sect_list):
#############################################################################
    for sect in sect_list:
        if (sect not in d.keys()):
            log_error(f'Section {sect} not found.')

#############################################################################
def check_obl_keys(d, section, obl_keys):
#############################################################################
    for key in obl_keys:
       if key not in d[section]:
           log_error(f'Key {key} not found in [{section}].')

#############################################################################
def key_present(d, s, k):
#############################################################################
    if (s not in d):
        return False
    if (k not in d[s]):
        return False
    return True

#############################################################################
def get_key(d, s, k, eval_k=False):
#############################################################################
    if not key_present(d, s, k):
        log_error(f'Could not find [{s}] {k}.')

    if eval_k:
        try:
            return eval(d[s][k])
        except:
            log_error(f'Keyword value could not be evaluated {d[s][k]}.')

    else:
        return d[s][k]

#############################################################################
def get_dat_nr(d, k):
#############################################################################
    return get_key(d, k, 'nr', eval_k=True)

#############################################################################
def get_dat_tuple(d, k):
#############################################################################
    file      = get_key(d, k, 'file')
    file_type = get_key(d, k, 'file_type')
    if file_type == 'binpos':
        pos_beg = get_key(d, k, 'binpos_beg', eval_k=True)
        pos_end = get_key(d, k, 'binpos_end', eval_k=True)
        t = (i_binpos, file, pos_beg, pos_end)
    elif file_type == 'bin':
        t = (i_bin, file)
    elif file_type == 'asc':
        t = (i_asc, file)
    return t

#############################################################################
def read_mod_ini(d_mod_ini, section, obl_key_list):
#############################################################################
    d = {}; d_map = od()
    for key in obl_key_list:
        if key_present(d_mod_ini, section, key):
            v = get_key(d_mod_ini, section, key)
            if v.isnumeric():
                t = (i_const, eval(v))
                #print(f'key = {key}, value = {t}')
                d[key] = t
            else:
                d_map[key] = v
        else:
            d_map[key] = key

    if section in d_mod_ini:
        for k, v in d_mod_ini[section].items():
            if ((k not in d) and (k not in d_map)):
                d[k] = v
    return d, d_map

#############################################################################
def get_csv_id(id, d_map):
#############################################################################
    if id not in d_map:
        return id
    else:
        return d_map[id]

#############################################################################
def determine_periods(section, d_mod_ini, d_mod_csv):
#############################################################################
    d = []; ids_csv = []
    if section not in d_mod_ini:
        return d
    ds = d_mod_ini[section]

    if 'periods' not in ds:
        log_error(f'Key periods not found in [{section}].')
    #
    # set the merge_binpos flag
    if 'merge_binpos' not in ds:
        merge_binpos = False
    else:
        try:
            merge_binpos = eval(ds['merge_binpos'])
            if (type(merge_binpos) != type(True)):
                log_error(f'Invalid value for merge_binpos in [{section}].')
        except:
            merge_binpos = False
    #
    # create the mapping
    d_map = {}
    for key in ds:
        if (key != 'merge_binpos' and key != 'periods'):
            d_map[key] = ds[key]

    # get the dictionary
    dsp_ini = ds['periods']

    # check if anything is to do for this submodel
    data_found = False
    for iper in dsp_ini.keys():
         for id in dsp_ini[iper]:
             id_csv = get_csv_id(id, d_map)
             if (id in d_mod_csv):
                  data_found = True
    if not data_found:
        return d, ids_csv

    # check if merging is possible
    if merge_binpos:
        merge = True
        for iper in dsp_ini.keys():
            d_tmp = od()
            for id in dsp_ini[iper]:
                id_csv = get_csv_id(id, d_map)
                if (id in d_mod_csv):
                    d_tmp[id] = get_dat_tuple(d_mod_csv, id_csv)
            file_type_list = []; file_list = []
            for id in d_tmp:
                file_type_list.append(d_tmp[id][0])
                file_list.append(d_tmp[id][1])

            if len(set(file_type_list)) != 1:
                merge = False
            else:
                if list(set(file_type_list))[0] != i_binpos:
                    merge = False
            if len(set(file_list)) != 1:
                merge = False
    else:
        merge = True
        #
        # merge in case there is only 1 input per stress period
        for iper in dsp_ini.keys():
            if (len(dsp_ini[iper]) > 1):
                merge = False
    #
    if merge:
        d.append({})
        d[0] = {'maxbound': 0, 'periods': od()}
        for iper in dsp_ini.keys():
            maxbound = 0
            d_tmp = od()
            for id in dsp_ini[iper]:
                id_csv = get_csv_id(id, d_map)
                if (id_csv in d_mod_csv):
                    maxbound += get_dat_nr(d_mod_csv, id_csv)
                    d_tmp[id] = get_dat_tuple(d_mod_csv, id_csv)

            first_id = list(d_tmp.keys())[0]
            last_id  = list(d_tmp.keys())[-1]
            d[0]['periods'][iper] = (i_binpos, d_tmp[first_id][1], \
                                          d_tmp[first_id][2], d_tmp[last_id][3])
            d[0]['maxbound'] = max(d[0]['maxbound'], maxbound)
    else:
        # get the unique ids and set as labels
        id_map = {}
        for iper in dsp_ini.keys():
            for id in dsp_ini[iper]:
                id_csv = get_csv_id(id, d_map)
                if (id_csv in d_mod_csv):
                    ids_csv.append(id_csv)
                    if id_csv not in id_map:
                        id_map[id_csv] = id
        ids_csv = list(set(ids_csv)); ids_csv.sort()
        #
        i = 0
        for id_csv in ids_csv:
            id = id_map[id_csv]
            d_tmp =  {'maxbound': 0, 'periods': od()}
            for iper in dsp_ini.keys():
               if (id in dsp_ini[iper]):
                   maxbound = get_dat_nr(d_mod_csv, id_csv)
                   d_tmp['maxbound'] = max(d_tmp['maxbound'], maxbound)
                   d_tmp['periods'][iper] = get_dat_tuple(d_mod_csv, id_csv)
            d.append(d_tmp)
    return d, ids_csv

#############################################################################
def write_exchanges(d_ini, d_xch, d_template, d_mf6_mod):
#############################################################################
    d_xch_files = {}
    if get_cla_key('serial_decoupled'):
        return d_xch_files

    log.info('Writing exchange files...')
    #
    d_xch_type = {}
    d_xch_type['exg-gwfgwf'] = 'gwf6-gwf6'
    d_xch_type['exg-gwtgwt'] = 'gwt6-gwt6'
    d_xch_type['exg-gwfgwt'] = 'gwf6-gwt6'
    #
    for mt in d_mf6_mod:
        section = d_mf6_mod[mt]['int_xch']
        if section is not None:
            xch_type = d_xch_type[section]
            d_xch_files[xch_type] = []
            if section is not None:
                fdir = get_key(d_ini, section, 'fdir')
                d_sect = {}
                for k, v in d_ini[section].items():
                    d_sect[k] = v
                d_sect.pop('fdir')
                #
                p = Path(fdir)
                if not p.exists():
                    log.info(f'Creating folder {p}')
                    p.mkdir()
                for xch in d_xch:
                    f = p.joinpath(xch + '.' + section)
                    xch_type_a = xch_type.split('-')[0]
                    xch_type_b = xch_type.split('-')[1]
                    exgmnamea = xch_type_a + '-' + xch.split('-')[0]
                    exgmnameb = xch_type_b + '-' + xch.split('-')[1]
                    d_xch_files[xch_type].append((f, exgmnamea, exgmnameb))
                    d = d_sect
                    d['nexg'] = d_xch[xch]['nr']
                    d['path'] = d_xch[xch]['file']
                    template = mf6_template(section)
                    template.set(d, valid_keys=d_template[section])
                    template.render(f)
    return d_xch_files

#############################################################################
def write_simulation(d_ini, d_xch_files, d_template, d_mf6_mod):
#############################################################################

    log.info('Writing simulation files...')

    parallel         = get_cla_key('parallel')
    nrproc = 1
    cgc              = get_cla_key('cgc')
    serial_decoupled = get_cla_key('serial_decoupled')
    if serial_decoupled:
        parallel = False

    ############
    # SIM-TDIS #
    ############
    section = 'sim-tdis'
    fname_tdis = get_key(d_ini, 'sim-tdis', 'fname')
    d = {}
    for k, v in d_ini[section].items():
        d[k] = v
    d.pop('fname')
    template = mf6_template(section)
    template.set(d, valid_keys=d_template[section])
    template.render(fname_tdis)

    f_base = get_key(d_ini, 'sim-nam', 'fname')

    #########################
    # SIM-NAM: models ascii #
    #########################
    fname_mod = None
    if not serial_decoupled:
        fname_mod = Path(f_base + '.mod.asc')
        tp_name = 'sim-nam-models'
        d = {}
        nrproc = 0
        for mt in d_mf6_mod:
            nrproc = max(nrproc, d_mf6_mod['gwf6']['nrproc'])
        if nrproc > 1:
            d['parallel'] = True
        d['models'] = []
        for mt in d_mf6_mod:
            d['models'].extend(d_mf6_mod[mt]['models'])
        template = mf6_template(tp_name)
        template.set(d, valid_keys=d_template[tp_name])
        template.render(fname_mod)

    ############################
    # SIM-NAM: exchanges ascii #
    ############################
    fname_xch = None
    if d_xch_files != {}:
        tp_name = 'sim-nam-exchanges'
        d = {}
        d['exchanges'] = []
        fname_xch = Path(f_base + '.xch.asc')
        for xch_type in d_xch_files:
            for xch in d_xch_files[xch_type]:
                exgtype   = xch_type
                exgfile   = xch[0]
                exgmnamea = xch[1]
                exgmnameb = xch[2]
                d['exchanges'].append((exgtype, exgfile, exgmnamea, exgmnameb))
            template = mf6_template(tp_name)
            template.set(d, valid_keys=d_template[tp_name])
            template.render(fname_xch)

    ###########
    # SIM-IMS #
    ###########
    d_ims = {}
    obl_keys = ['outer_dvclose', 'outer_maximum', \
                'inner_maximum', 'inner_dvclose', 'inner_rclose', \
                'linear_acceleration']
    for mt in d_mf6_mod:
        tp_name = 'sln-ims'
        section = mt + '-' + tp_name
        check_main_ini(d_ini, [section])
        check_obl_keys(d_ini, section, obl_keys)
        fname = get_key(d_ini, section, 'fname')
        d_ims[mt] = fname
        d = {}
        for k, v in d_ini[section].items():
            d[k] = v
        d.pop('fname')
        template = mf6_template(tp_name)
        template.set(d, valid_keys=d_template[tp_name])
        template.render(fname)

    ###################################
    # SIM-NAM: solution models ascii #
    ###################################
    d_smo = {}
    if not serial_decoupled:
        tp_name = 'sim-nam-solution-models'
        if cgc:
            d['coarse_grid_correction'] = True
        for mt in d_mf6_mod:
            fname =  Path(f_base + '.smo.'+ mt +'.asc')
            d_smo[mt] = fname
            d = {}
            d['solutionmodels'] = []
            i_cgc = 0
            for mod in d_mf6_mod[mt]['models']:
                id = mod[2]
                i_cgc += 1
                if cgc:
                    d['solutionmodels'].append((id, i_cgc))
                else:
                    d['solutionmodels'].append((id))
            template = mf6_template(tp_name)
            template.set(d, valid_keys=d_template[tp_name])
            template.render(fname)

    ##########################################
    # SIM-NAM: solution models ascii wrapper #
    ##########################################
    d_smw = {}
    if not serial_decoupled:
        tp_name = 'sim-nam-solution-models-wrapper'
        for mt in d_mf6_mod:
            fname = Path(f_base + '.smw.asc')
            d_smw[mt] = fname
            d = {}
            d['solmodels'] = d_smo[mt]
            template = mf6_template(tp_name)
            template.set(d, valid_keys=d_template[tp_name])
            template.render(fname)

    ###########
    # SIM-NAM #
    ###########
    fname_nam = Path(f_base + '.nam')
    mfsim_list = write_mfsim(d_ini, d_template, d_mf6_mod, d_ims, d_smw, \
        fname_nam, fname_tdis, fname_mod, fname_xch, \
        nrproc)
    return mfsim_list


#############################################################################
def write_mfsim(d_ini, d_template, d_mf6_mod, d_ims, d_smw, \
    fname_nam, fname_tdis, fname_mod, fname_xch, \
    nrproc):
#############################################################################
    serial_decoupled = get_cla_key('serial_decoupled')
    section = 'sim-nam'
    d = {}
    for k, v in d_ini[section].items():
        d[k] = v
    d.pop('fname')
    d['tdis6'] = fname_tdis
    #
    mfsim_list = []
    if not serial_decoupled:
        if nrproc > 1:
            d['domain_decomposition'] = nrproc
        d['models'] = fname_mod
        d['models_is_file'] = True
        d['filein'] = True
        if fname_xch is not None:
            d['exchanges'] = fname_xch
        lst = []
        for mt in d_mf6_mod:
            lst.append(('ims6', d_ims[mt], d_smw[mt]))
        d['solutiongroups'] = []
        d['solutiongroups'].append(lst)
        template = mf6_template(section)
        template.set(d, valid_keys=d_template[section])
        template.render(fname_nam)
        mfsim_list.append(fname_nam)
    else:
        for mt in d_mf6_mod:
            for mod in d_mf6_mod[mt]['models']:
                id = mod[2]
                d['models'] = " ".join(mod)
                d['solutiongroups'] = []
                p = Path(fname_nam)
                fname = p.parents[0].joinpath(p.stem + '.' + str(id) + '.nam')
                lst = [('ims6', d_ims[mt], [id])]
                d['solutiongroups'].append(lst)
                template = mf6_template(section)
                template.set(d, valid_keys=d_template[section])
                template.render(fname)
                mfsim_list.append(fname)

    #d['solutiongroups'] = [[('ims6', 'test.gwf.ims', 'sol1.asc'), \
    #                        ('ims6', 'test.gwt.ims', 'sol1.asc')], \
    #                       [('ims6', 'test.gwf.ims', 'sol1.asc')]]

    return mfsim_list

#############################################################################
def write_gwf_model(id, nodes, nja, d_mod_ini, d_template, mod_dir, d_mod_csv):
#############################################################################
    # first, check if the nam file is present
    if not 'gwf-nam' in d_mod_ini:
        return None

    ###############
    # GWF modules #
    ###############
    d_mf6_m = od()
    d_mf6_m['gwf-disu'] = {'ftype':'disu6','fname':str(), \
        'obl_keys': ['top','bot','area','idomain','iac','ja','ihc','cl12','hwva']}
    d_mf6_m['gwf-npf']  = {'ftype':'npf6' ,'fname':str(), \
        'obl_keys': ['icelltype','k','k33']}
    d_mf6_m['gwf-sto']  = {'ftype':'sto6' ,'fname':str(), \
        'obl_keys': ['iconvert','ss','sy']}
    d_mf6_m['gwf-ic']   = {'ftype':'ic6'  ,'fname':str(), \
        'obl_keys': ['strt']}
    for module in d_mf6_m:
        log.info(f'Writing file for {module}...')
        template = mf6_template(module)
        d, d_map = read_mod_ini(d_mod_ini, module, d_mf6_m[module]['obl_keys'])
        if (module == 'gwf-disu'):
            d['nodes'] = nodes
            d['nja']   = nja
        if module == 'gwf-ic':
            for key in d_map:
                if key == 'strt':
                    v = d_map[key]
                    if v not in d_mod_csv:
                        fname = Path(v.replace('{{model_id}}',str(id)))
                        if not fname.is_file():
                           log_error(f'Invalid key value for {key} in {module}: {v}')
                        else:
                            d[key] = (i_bin, fname)
                    else:
                        d[key] = get_dat_tuple(d_mod_csv, d_map[key])
                else:
                    d[key] = get_dat_tuple(d_mod_csv, d_map[key])
        else:
            for key in d_map:
                d[key] = get_dat_tuple(d_mod_csv, d_map[key])
        d_mf6_m[module]['fname'] = mod_dir / f"{id}.gwf.{d_mf6_m[module]['ftype']}"
        template.set(d, valid_keys=d_template[module])
        template.render(d_mf6_m[module]['fname'])

    ###############
    # GWF-OC      #
    ###############
    l_mf6_p = []
    section = 'gwf-oc'
    if section in d_mod_ini:
        d = {}
        for k, v in d_mod_ini[section].items():
            d[k] = v
        ftype = 'oc6'
        fname = mod_dir / f"{id}.gwf.{ftype}"
        replace_file = ['budgetfile','headfile','concentrationfile']
        for f in replace_file:
            if f in d:
                d[f] = Path(d[f].replace('{{model_id}}',str(id)))
        template = mf6_template(section)
        template.set(d, valid_keys=d_template[section])
        template.render(fname)
        #
        d_package = {}
        d_package['ftype'] = ftype
        d_package['fname'] = fname
        l_mf6_p.append(d_package)

    #######################
    # GWF stress-packages #
    #######################
    packages = ['gwf-riv','gwf-drn','gwf-wel','gwf-ghb','gwf-chd','gwf-rch']
    for package in packages:
        if package in d_mod_ini:
            periods_list, ids_csv = determine_periods(package, d_mod_ini, d_mod_csv)
            d = {}
            if len(periods_list) == 1:
                log.info(f'***** Writing single file for {package}...')
                d = periods_list[0]
                ftype =  package.split('-')[1]+'6'
                fname = mod_dir / f"{id}.gwf.{ftype}"
                #
                section = f'{package}-options'
                if section in d_mod_ini:
                    for k, v in d_mod_ini[section].items():
                        d[k] = v
                template = mf6_template(package)
                template.set(d, valid_keys=d_template[package])
                template.render(fname)
                #
                d_package = {}
                d_package['ftype'] = ftype
                d_package['fname'] = fname
                l_mf6_p.append(d_package)
            elif len(periods_list) > 1:
                log.info(f'***** Writing multiple file for {package}...')
                i = 0
                for id_csv in ids_csv:
                    d = periods_list[i]; i += 1
                    ftype =  package.split('-')[1]+'6'
                    fname = mod_dir / f"{id}.gwf.{id_csv}.{ftype}"
                    section = f'{package}-options'
                    if section in d_mod_ini:
                        for k, v in d_mod_ini[section].items():
                            d[k] = v
                    template = mf6_template(package)
                    template.set(d, valid_keys=d_template[package])
                    template.render(fname)
                    #
                    d_package = {}
                    d_package['ftype'] = ftype
                    d_package['fname'] = fname
                    l_mf6_p.append(d_package)

    ###########
    # GWF-NAM #
    ###########
    section = 'gwf-nam'
    template = mf6_template(section)
    d, d_map = read_mod_ini(d_mod_ini, section, [])
    d['packages'] = []
    for package in d_mf6_m:
        d['packages'].append((d_mf6_m[package]['ftype'], d_mf6_m[package]['fname'], ''))
    for package in l_mf6_p:
        d['packages'].append((package['ftype'], package['fname'], ''))

    if 'listing_file' in d:
        s = d['listing_file']
        d['listing_file'] = Path(s.replace('{{model_id}}',str(id)))
    template.set(d, valid_keys=d_template[section])
    fname = mod_dir / f"{id}.gwf.nam"
    template.render(fname)

    return fname

#############################################################################
def write_gwt_model(id, d_mod_ini, d_template, mod_dir, d_mod_csv): #TODO
#############################################################################
    # first, check if the nam file is present
    if not 'gwt-nam' in d_mod_ini:
        return None
    else:
        log_error('GWT is not yet supported.')

#############################################################################
def mf6_model_admin(d_ini, d_xch, mf6_mod_lst):
#############################################################################
    parallel = get_cla_key('parallel')

    model_type = [];
    for i in range(len(mf6_mod_lst)):
        mtype, mfname, name, part = mf6_mod_lst[i]
        model_type.append(mtype)
    model_type = list(set(model_type))

    d_mf6_mod = {}
    for mt in model_type:
       # internal: submodel <--> submodel
       d_mf6_mod[mt] = {'int_xch': None, 'ext_xch': None, \
        'nrproc': 1, 'nr_models': 0, 'models': []}

    # set the number of models
    for mt in model_type:
        for i in range(len(mf6_mod_lst)):
            mtype, mfname, name, part = mf6_mod_lst[i]
            if (mtype == mt):
                d_mf6_mod[mt]['nr_models'] += 1

   # if parallel, set nrproc, set renumber partitions
    if parallel:
        for mt in model_type:
            partitions = []
            for i in range(len(mf6_mod_lst)):
                mtype, mfname, name, part = mf6_mod_lst[i]
                if (mtype == mt):
                    partitions.append(part)
            partitions = list(set(partitions))
            nrproc = len(partitions)
            d_mf6_mod[mt]['nrproc'] = nrproc
            if nrproc > 1:
                d_map = {}
                i = 0
                for part in partitions:
                    i += 1
                    d_map[part] = i
                for i in range(len(mf6_mod_lst)):
                    mtype, mfname, name, part = mf6_mod_lst[i]
                    if (mtype == mt):
                        mf6_mod_lst[i][3] = d_map[part]

    # the main trigger for coupling is the keyword exg-gwfgwf
    if 'exg-gwfgwf' in d_ini:
        if 'gwf6' in d_mf6_mod:
            d_mf6_mod['gwf6']['int_xch'] = 'exg-gwfgwf'
    if 'exg-gwtgwt' in d_ini:
        if 'gwt6' in d_mf6_mod:
            log.error('exg-gwtgwt is not yet supported')
            d_mf6_mod['gwt6']['int_xch'] = 'exg-gwtgwt'
    if 'exg-gwfgwt' in d_ini:
        if (('gwf6' in d_mf6_mod) and ('gwt6' in d_mf6_mod)):
            log.error('exg-gwfgwt is not yet supported')
            d_mf6_mod['gwf6']['ext_xch'] = 'exg-gwfgwt'
            d_mf6_mod['gwt6']['ext_xch'] = 'exg-gwfgwt'

    # then, check for nrproc and numerical connection do disable internal coupling
    for mt in model_type:
        if d_xch == {}:
            d_mf6_mod[mt]['int_xch'] = None
        #if d_mf6_mod[mt]['nrproc'] == 1:
        #    raise
        #    d_mf6_mod[mt]['int_xch'] = None
        #else:
        #    if d_xch == {}:
        #        d_mf6_mod[mt]['int_xch'] = None

    if get_cla_key('serial_decoupled'):
        for mt in model_type:
            d_mf6_mod[mt]['int_xch'] = None
            d_mf6_mod[mt]['ext_xch'] = None

    for mt in model_type:
        for i in range(len(mf6_mod_lst)):
            mtype, mfname, name, part = mf6_mod_lst[i]
            if d_mf6_mod[mt]['nrproc'] > 1:
                d_mf6_mod[mt]['models'].append((mtype, str(mfname), name, part))
            else:
                d_mf6_mod[mt]['models'].append((mtype, str(mfname), name))
    return d_mf6_mod


#############################################################################
def write_model(id, d_ini, d_props, d_mod_ini_list, d_template):
#############################################################################
    log.info(f'Writing MODFLOW 6 model files for id={id}...')
    #
    f_csv_dat = get_key(d_props, str(id), 'csv_dat')
    nodes     = int(get_key(d_props, str(id), 'nodes'))
    nja       = int(get_key(d_props, str(id), 'nja'))
    #
    lay_mod = get_key(d_props, str(id), 'lay_mod', eval_k=True)
    d_mod_ini = d_mod_ini_list[lay_mod-1]
    #
    # read the file for csv_dat
    p = Path(f_csv_dat)
    mod_dir = p.parents[0]
    if not p.is_file():
        log_error(f'File {f_csv_dat} not found.')
    d_mod_csv = read_csv(str(p))

    # models
    gwf_mfname = write_gwf_model(id, nodes, nja, d_mod_ini, d_template, mod_dir, d_mod_csv)
    gwt_mfname = write_gwt_model(id, d_mod_ini, d_template, mod_dir, d_mod_csv)

    # get partition number
    parallel = get_cla_key('parallel')
    if parallel:
        field = get_key(d_ini, '_general', 'props_part_field')
        part = get_key(d_props, str(id), field, eval_k=True)
    else:
        part = 0

    mf6_mod_lst = []
    if gwf_mfname is not None:
        mf6_mod_lst.append(['gwf6', gwf_mfname, 'gwf6-'+str(id), part])
    if gwt_mfname is not None:
        mf6_mod_lst.append(['gwt6', gwt_mfname, 'gwt6-'+str(id), part])

    return mf6_mod_lst

#############################################################################
def run_model(mf6, mfsim, np=1, mpi=None):
#############################################################################

    p_sim = Path(mfsim)
    moddir  = p_sim.parents[0]
    namfile = p_sim.name
    #
    if np > 1:
        cmd = f'{mpi} -np {np} {mf6} -s {namfile}'
    else:
        cmd = f'{mf6} -s {namfile}'
    log.info(f'Running {cmd}...')

    os.chdir(moddir)
    p = Popen(cmd, cwd=moddir, stdout=PIPE, stderr=PIPE)

    stdout = []
    while True:
        line = p.stdout.readline().decode('ascii')
        if not line:
            break
        else:
            print(' '+line.rstrip())
            stdout.append(line)
    out, err = p.communicate()
    ierr = p.returncode
    log.info(f'Done running with errorcode {ierr}...')

    return ierr, stdout

#############################################################################
def pre():
#############################################################################
    # read the input ini-file
    f_ini = get_cla_key('ini')

    d_ini = read_ini(f_ini)
    check_main_ini(d_ini,['_general', 'sim-tdis', 'sim-nam'])
    #
    # get the list of model ids
    if key_present(d_ini, '_general', 'model_id'):
        # read the model ids
        id_list = get_key(d_ini, '_general', 'model_id', eval_k=True)
    else:
        # ALL ids
        id_list = []

    # read the properties
    f = get_key(d_ini, '_general', 'props_csv')
    d_props = read_csv(f, filter_ids=id_list)

    # read the exchanges
    f = get_key(d_ini, '_general', 'exchanges_csv')
    d_xch = read_csv(f, filter_ids=id_list)

    # set the id list
    if not id_list:
        id_list = list(d_props.keys())

    # read the model definitions
    mod_def = get_key(d_ini, '_general', 'model_definition', eval_k=True)
    d_mod_ini_list = []
    for f in mod_def:
        d_mod_ini_list.append(read_ini(f))


    # get keys for supported templates
    supp_templates = ['gwf-disu','gwf-npf','gwf-ic', 'gwf-oc', 'gwf-nam', \
                      'gwf-riv', 'gwf-drn', 'gwf-wel', 'gwf-ghb', 'gwf-chd', 'gwf-rch', \
                      'sim-nam', 'sim-nam-models', 'sim-nam-exchanges', \
                      'sim-nam-solution-models', 'sim-nam-solution-models-wrapper', 'sim-tdis',
                      'exg-gwfgwf','sln-ims']
    d_template = {}
    for section in supp_templates:
        log.info(f'{section}')
        template = mf6_template(section)
        d_template[section] = template.get_variables()

    # write the model files
    n = len(id_list); i = 0
    mf6_mod_lst = []
    mfsim_list = []
    for id in id_list:
        i += 1
        log.info(10*'='+f' Processing {i:06d}/{n:06d} '+10*'=')
        lst = write_model(id, d_ini, d_props, d_mod_ini_list, d_template)
        mf6_mod_lst.extend(lst)

    # set up the administration
    d_mf6_mod = mf6_model_admin(d_ini, d_xch, mf6_mod_lst)

    # write the exchanges
    d_xch_files = write_exchanges(d_ini, d_xch, d_template, d_mf6_mod)

    # write simulation files
    mfsim_list = write_simulation(d_ini, d_xch_files, d_template, d_mf6_mod)

    return mfsim_list

#############################################################################
def run(mfsim_list):
#############################################################################
    mf6 = get_cla_key('mf6')
    mpi = get_cla_key('mpi')
    np  = get_cla_key('np')

    for mfsim in mfsim_list:
        if np > 1:
            ierr, stdout = run_model(mf6, mfsim, np, mpi)
        else:
            ierr, stdout = run_model(mf6, mfsim)

        if ierr != 0:
           log.info('***** MODFLOW 6 RUN FAILED: ' + str(mfsim))
           log.info(20*'=' + ' BEGIN STDOUT MODFLOW 6 ' + 20*'=')
           for line in stdout:
               log.info(line.rstrip())
           log.info(20*'='+ ' END STDOUT MODFLOW 6 ' + 20*'=')

#############################################################################
def post():
#############################################################################
    pass

#############################################################################
def main():
#############################################################################
    if get_cla_key('pre'):
        mfsim_list = pre()
    else:
        mfsim_list = [get_cla_key('mfsim')]
    if get_cla_key('run'):
        run(mfsim_list)
    if get_cla_key('post'):
        post()

if __name__ == '__main__':
    main()


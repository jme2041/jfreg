#!/usr/bin/env python3
# jfreg: Joint fMRI Registration
#
# Copyright (c) 2019, Jeffrey M. Engelmann
#
# jfreg is released under the revised (3-clause) BSD license.
# For details, see LICENSE.txt
#

import sys
import os
import shutil
import re
import argparse
import datetime
import tempfile
import subprocess

class FuncDset:
    """Store information about a functional dataset"""
    def __init__(self, f, outdir):
        assert os.path.isfile(f)
        assert os.path.isdir(outdir)
        self.outdir = outdir
        self.f = f

        self.mc       = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_mc.nii')
        self.par      = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_mc.par')
        self.mc_mat   = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_mc.mat')
        self.base     = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_base.nii')
        self.base_t1  = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_base_t1.nii')
        self.base2t1  = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_base_to_t1.mat')
        self.base2mni = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_base_to_mni.mat')
        self.func2mni = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_to_mni.mat')
        self.mc_mni   = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_mc_mni.nii')
        self.extents  = os.path.join(outdir, os.path.basename(self.f)[:-4] + '_mc_mni_extents.nii')

        self.files = [
            self.mc,
            self.par,
            self.base,
            self.base_t1,
            self.base2t1,
            self.base2mni,
            self.mc_mni,
            self.extents ]

        self.dirs = [ self.mc_mat, self.func2mni ]

        cp = cmd('fslnvols', self.f, echo=False, capture_output=True)
        self.nvols = int(cp.stdout.decode().splitlines()[0])

    def check_if_exists(self, overwrite):

        for f in self.files:
            if os.path.isfile(f):
                if overwrite:
                    os.remove(f)
                else:
                    raise IOError('File already exists: %s' % f)

        for d in self.dirs:
            if os.path.isdir(d):
                if overwrite:
                    shutil.rmtree(d)
                else:
                    raise IOError('Directory already exists: %s' % d)

    def delete_all(self):
        for f in self.files:
            if os.path.isfile(f):
                os.remove(f)

        for d in self.dirs:
            if os.path.isdir(d):
                shutil.rmtree(d)

def check_nifti_extension(f):
#    """Check NIfTI extension (jfreg requires single-file unzipped)"""
    if not os.path.isfile(f):
        raise IOError('Could not find file: %s' % f)
    if not f.endswith('.nii'):
        raise IOError(
            'Files must be non-gzipped single-file NIfTI datasets (.nii)')

def cmd(*cmd, env=None, echo=True, capture_output=False):
    """Run a shell command (echoing it to stdout if requested)"""
    # This needs to be compatible with Python 3.6, so use stdout= and stderr=
    # instead of capture_output=
    if(echo):
        print(' '.join(cmd))
    if capture_output:
        return subprocess.run(cmd, check=True, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        return subprocess.run(cmd, check=True, env=env)

def main(argv):
    """Main jfreg routine"""
    version = '0.1.0'

    parser = argparse.ArgumentParser(
        prog='jfreg',
        add_help=False,
        allow_abbrev=False,
        description='Joint fMRI registration',
        usage='%(prog)s [options] -head HEAD -brain BRAIN FUNC [FUNC ...]')

    g1 = parser.add_argument_group('required parameters')
    g1.add_argument('-head',
        required=True,
        metavar='HEAD',
        help='T1-weighted dataset in subject space (with skull)')
    g1.add_argument('-brain',
        required=True,
        metavar='BRAIN',
        help='Brain-extracted T1-weighted dataset in subject space')

    g2 = parser.add_argument_group('options')
    g2.add_argument('-basevol',
        type=int,
        default=0,
        metavar='VOL',
        help='Functional base volume number (0-indexed; default: %(default)d)')
    g2.add_argument('-ores',
        type=float,
        default=3.0,
        metavar='MM',
        help='''Output resolution for functional dataset (mm;
        default: %(default)1.1f)''')
    g2.add_argument('-bbr',
        metavar='WMSEG',
        help='Use brain-based registration (WMSEG: white-matter mask)')
    g2.add_argument('-search',
        type=int,
        default=180,
        choices=[0, 90, 180],
        metavar='DEG',
        help='''Search angles in degrees (%(choices)s;
        default: %(default)d)''')
    g2.add_argument('-outdir',
        default=os.getcwd(),
        help='Output directory (default: current working directory)')
    g2.add_argument('-overwrite',
        action='store_true',
        help='Overwrite output datasets that already exist')
    g2.add_argument('-version',
        action='version',
        version=version,
        help='Show version number and exit')
    g2.add_argument('-help',
        action='help',
        help='Show this help message and exit')

    parser.add_argument('func',
        nargs='+',
        metavar='FUNC',
        help='Functional dataset(s)')

    # Parse arguments
    opts = parser.parse_args(argv)

    # Set environment variables
    jfreg_env = os.environ
    jfreg_env['FSLOUTPUTTYPE'] = 'NIFTI'

    # Check FLIRT version
    FLIRT_MIN = (6, 0)
    cp = cmd('flirt', '-version', env=jfreg_env, echo=False, capture_output=True)
    flirt_verstr = cp.stdout.decode().splitlines()[0]
    pattern = r'^FLIRT\s+version\s+(\d+)\.(\d+).*$'
    match = re.search(pattern, flirt_verstr)
    assert match is not None
    flirt_ver = tuple(map(int, match.groups()))
    if flirt_ver < FLIRT_MIN:
        raise RuntimeError('FLIRT %d.%d or newer is required' % FLIRT_MIN)

    # Get FSL MNI template and FLIRT schedule
    # This verifies that FSL is configured correctly
    mni_template = os.path.abspath(os.path.join(os.environ['FSLDIR'],
        'data', 'standard', 'MNI152_T1_2mm_brain.nii.gz'))
    if not os.path.isfile(mni_template):
        raise IOError('Could not find MNI template: %s' % mni_template)

    bbr_sch = os.path.abspath(os.path.join(os.environ['FSLDIR'],
        'etc', 'flirtsch', 'bbr.sch'))
    if not os.path.isfile(bbr_sch):
        raise IOError('Could not find BBR FLIRT schedule')

    # Check the range of numeric parameters
        if opts.ores <= 0.0:
            raise ValueError('Output resolution must be greater than 0.0')
        if opts.basevol < 0:
            raise ValueError('Base volume cannot be negative')

    # Check that input datasets exist and are .nii files
    head = os.path.abspath(opts.head)
    brain = os.path.abspath(opts.brain)
    func = [os.path.abspath(f) for f in opts.func]
    wmseg = os.path.abspath(opts.bbr) if opts.bbr else None
    check_nifti_extension(head)
    check_nifti_extension(brain)
    [check_nifti_extension(f) for f in func]
    if wmseg:
        check_nifti_extension(wmseg)

    # Check output directory
    outdir = os.path.abspath(opts.outdir)
    if not os.path.isdir(opts.outdir):
        raise IOError('Could not find output directory: %s' % outdir)

    # Check if structural output already exists. Delete if overwriting.
    brain2mni   = os.path.join(outdir, os.path.basename(brain)[:-4] + '_to_mni.mat')
    brain_mni   = os.path.join(outdir, os.path.basename(brain)[:-4] + '_mni.nii')
    head_mni    = os.path.join(outdir, os.path.basename(head)[:-4]  + '_mni.nii')

    structs = [ brain2mni, brain_mni, head_mni ]
    for f in structs:
        if os.path.isfile(f):
            if opts.overwrite:
                os.remove(f)
            else:
                raise IOError('File already exists: %s' % f)

    # Loop on functional datasets.
    # Check if output already exists. Delete if overwriting.
    funcs = [ FuncDset(f, outdir) for f in func ]
    for f in funcs:
        f.check_if_exists(opts.overwrite)

    func_mni_template = None

    try:
        print('jfreg version %s begins (%s)' % (version,
            datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

        # ################################################
        # SPATIAL NORMALIZATION OF THE T1-WEIGHTED DATASET
        # ################################################

        # Register the brain-extracted T1-weighted dataset to the MNI template.
        # Use 12-parameter affine transformation.
        cmd('flirt',
            '-in'       , brain,
            '-ref'      , mni_template,
            '-out'      , brain_mni[:-4],
            '-omat'     , brain2mni,
            '-cost'     , 'corratio',
            '-dof'      , str(12),
            '-searchrx' , str(-opts.search), str(opts.search),
            '-searchry' , str(-opts.search), str(opts.search),
            '-searchrz' , str(-opts.search), str(opts.search),
            '-interp'   , 'trilinear',
            env=jfreg_env)
        assert os.path.isfile(brain_mni)
        assert os.path.isfile(brain2mni)

        # Apply the resulting transformation to the non-brain extracted
        # T1-weighted dataset.
        cmd('flirt',
            '-in'       , head,
            '-ref'      , mni_template,
            '-out'      , head_mni[:-4],
            '-init'     , brain2mni,
            '-applyxfm' ,
            '-interp'   , 'trilinear',
            env=jfreg_env)
        assert os.path.isfile(head_mni)

        # Create a template dataset in MNI space at the desired output
        # resolution. This is a temporary file that is only used to set the
        # grid. Base the name of this file on the brain dataset to avoid
        # conflict with existing files
        func_mni_template = os.path.join(outdir, os.path.basename(brain)[:-4] + '_ftemp.nii')
        assert not os.path.isfile(func_mni_template)
        cmd('flirt',
            '-in'          , mni_template,
            '-ref'         , mni_template,
            '-out'         , func_mni_template[:-4],
            '-applyisoxfm' , str(opts.ores),
            '-interp'      , 'trilinear',
            env=jfreg_env)
        assert os.path.isfile(func_mni_template)

        # ###########################
        # LOOP ON FUNCTIONAL DATASETS
        # ###########################

        for func in funcs:

            # #################
            # MOTION CORRECTION
            # #################

            cmd('mcflirt',
                '-in'           , func.f,
                '-out'          , func.mc[:-4],
                '-refvol'       , str(opts.basevol),
                '-mats'         ,
                '-plots'        ,
                '-spline_final' ,
                env=jfreg_env)
            assert os.path.isfile(func.mc)
            assert os.path.isfile(func.par)
            assert os.path.isdir(func.mc_mat)

            # ####################################
            # FUNCTIONAL-STRUCTURAL COREGISTRATION
            # ####################################

            # Extract functional base volume
            cmd('fslroi',
                func.f,
                func.base[:-4],
                str(opts.basevol),
                str(1),
                env=jfreg_env)
            assert os.path.isfile(func.base)

            # Temporary transformation matrix file
            init_mat = os.path.join(outdir, func.base2t1)[:-4] + '_init.mat'
            if(os.path.exists(init_mat)):
                os.remove(init_mat)

            try:
                # Begin with 6-DOF registration
                cmd('flirt',
                    '-in'       , func.base,
                    '-ref'      , brain,
                    '-out'      , func.base_t1[:-4],
                    '-omat'     , init_mat,
                    '-cost'     , 'corratio',
                    '-dof'      , str(6),
                    '-searchrx' , str(-opts.search), str(opts.search),
                    '-searchry' , str(-opts.search), str(opts.search),
                    '-searchrz' , str(-opts.search), str(opts.search),
                    '-interp'   , 'trilinear',
                    env=jfreg_env)
                assert os.path.isfile(func.base_t1)
                assert os.path.isfile(init_mat)

                if opts.bbr:
                    # BBR (if requested)
                    os.remove(func.base_t1)
                    cmd('flirt',
                        '-in'       , func.base,
                        '-ref'      , head,
                        '-out'      , func.base_t1[:-4],
                        '-omat'     , func.base2t1,
                        '-cost'     , 'bbr',
                        '-dof'      , str(6),
                        '-wmseg'    , wmseg,
                        '-init'     , init_mat,
                        '-schedule' , bbr_sch,
                        '-searchrx' , str(-opts.search), str(opts.search),
                        '-searchry' , str(-opts.search), str(opts.search),
                        '-searchrz' , str(-opts.search), str(opts.search),
                        '-interp'   , 'trilinear',
                        env=jfreg_env)
                    assert os.path.isfile(func.base_t1)
                    assert os.path.isfile(func.base2t1)

                else:
                    assert not os.path.isfile(func.base2t1)
                    shutil.move(init_mat, func.base2t1)

            finally:
                if os.path.isfile(init_mat):
                    os.remove(init_mat)

            # Final output with spline interpolation
            os.remove(func.base_t1)
            cmd('flirt',
                '-in'       , func.base,
                '-ref'      , head,
                '-out'      , func.base_t1[:-4],
                '-init'     , func.base2t1,
                '-applyxfm' ,
                '-interp'   , 'spline',
                env=jfreg_env)
            assert os.path.isfile(func.base_t1)

            # ####################################
            # WARP FUNCTIONAL DATASET TO MNI SPACE
            # ####################################

            ones = None
            ones_mni = None

            try:
                # Create a temporary, all-1 dataset for extents masking
                ones = os.path.join(outdir, os.path.basename(func.f[:-4]) + '_ones.nii')
                cmd('fslmaths', func.f, '-abs', '-add', str(1), '-bin', ones,
                    '-odt', 'char', env=jfreg_env)
                assert(os.path.isfile(ones))

                # Concatenate T1 to MNI and functional base to T1 transforms to
                # get the functional base to MNI transform
                cmd('convert_xfm', '-omat', func.base2mni,
                    '-concat', brain2mni, func.base2t1, env=jfreg_env)
                assert os.path.isfile(func.base2mni)

                # Create a directory for storing the functional-to-MNI
                # transformation matrices.
                os.makedirs(func.func2mni)

                with tempfile.TemporaryDirectory(dir=outdir) as td:
                    # Split the functional time series into individual volumes
                    # This is necessary because applyxfm4D does not provide
                    # control over interpolation method (nearest neighbor is
                    # needed for extents masking)
                    func_vol_base = os.path.join(td, 'func_vol')
                    func_mni_base = os.path.join(td, 'func_mni')
                    cmd('fslsplit', func.f, func_vol_base, '-t', env=jfreg_env)

                    # Split the all-ones dataset into individual volumes
                    ones_vol_base = os.path.join(td, 'ones_vol')
                    ones_mni_base = os.path.join(td, 'ones_mni')
                    cmd('fslsplit', ones, ones_vol_base, '-t', env=jfreg_env)

                    # Loop on volumes
                    for v in range (0, func.nvols):
                        func_in  = os.path.join(td, '%s%04d.nii' % (func_vol_base, v))
                        func_out = os.path.join(td, '%s%04d.nii' % (func_mni_base, v))
                        ones_in  = os.path.join(td, '%s%04d.nii' % (ones_vol_base, v))
                        ones_out = os.path.join(td, '%s%04d.nii' % (ones_mni_base, v))

                        assert(os.path.isfile(func_in))
                        assert(os.path.isfile(ones_in))
                        assert(not os.path.isfile(func_out))
                        assert(not os.path.isfile(ones_out))

                        # Concatenate the functional to base transform for each
                        # volume with the base to MNI transform to get the
                        # functional to MNI transform for this volume

                        func2base = os.path.join(func.mc_mat  , 'MAT_%04d' % v)
                        func2mni  = os.path.join(func.func2mni, 'F2MNI_%04d' % v)

                        assert(os.path.isfile(func2base))
                        assert(not os.path.isfile(func2mni))

                        cmd('convert_xfm', '-omat', func2mni,
                            '-concat', func.base2mni, func2base,
                            env=jfreg_env)

                        assert os.path.isfile(func2mni)

                        # Apply the functional to MNI transform to get the volume
                        # in standard space
                        # Use sync (Blackman) interpolation to match applyxfm4D
                        cmd('flirt',
                            '-in'         , func_in,
                            '-ref'        , func_mni_template,
                            '-out'        , func_out,
                            '-applyxfm'   ,
                            '-init'       , func2mni,
                            '-interp'     , 'sinc',
                            '-sincwindow' , 'blackman',
                            env=jfreg_env)

                        # Apply the functional to MNI transform to the ones
                        # dataset to get the extents of the transform for each
                        # volume. Use nearest neighbor interpolation.
                        cmd('flirt',
                            '-in'       , ones_in,
                            '-ref'      , func_mni_template,
                            '-out'      , ones_out,
                            '-applyxfm' ,
                            '-init'     , func2mni,
                            '-interp'   , 'nearestneighbour')

                    # Merge the functional data
                    func_mni_files = ['%s%04d.nii' % (func_mni_base, v) for v in range(0, func.nvols)]
                    cmd('fslmerge', '-t', func.mc_mni, *func_mni_files, env=jfreg_env)
                    assert os.path.isfile(func.mc_mni)

                    # Merge the ones dataset
                    ones_mni_files = ['%s%04d.nii' % (ones_mni_base, v) for v in range(0, func.nvols)]
                    ones_mni = ones[:-4] + '_mni.nii'
                    cmd('fslmerge', '-t', ones_mni, *ones_mni_files, env=jfreg_env)

                    # Create the extents mask
                    cmd('fslmaths', ones_mni, '-Tmin', func.extents,
                        '-odt', 'char', env=jfreg_env)

            finally:
                if os.path.isfile(ones):
                    os.remove(ones)
                if os.path.isfile(ones_mni):
                    os.remove(ones_mni)

    except:
        for f in structs:
            if os.path.isfile(f):
                os.remove(f)

        for f in funcs:
            f.delete_all()

        raise

    finally:
        if os.path.isfile(func_mni_template):
            os.remove(func_mni_template)

        print('jfreg ends (%s)' %
            datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

if __name__ == '__main__':
    main(sys.argv[1:])

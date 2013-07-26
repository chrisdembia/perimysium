"""

To be used in Jython.

"""
import datetime
import os
import xml.etree.ElementTree as etree

import dataman

def experiment(cmc_setup_fpath, parent_dir, name, description, fcn,
        minimal=True, run_command=False, overwrite=False):
    """

    cmc_input is 'setup', 'model', 'control_constraints', 'tasks', 'actuators'

    Parameters
    ----------
    minimal : str, optional (default: True)
        If True, copy over only the modified files; rely otherwise on the
        original files that CMC setup file referred to. If False, all
        dependencies are copied over.
    run_command : str, optional
        If CMC is to be run after all the experiment files are written, this is
        the command to run OpenSim CMC (e.g., 'cmc').
    overwrite : bool, optional (default: False)
        If True, files will be written even if experiment directory already
        exists.

    """
    # TODO optional, doc into a pytable.
    # TODO deal with potentially overwriting changes that the user also makes
    # to the setup file (force set files).
    # TODO make it easy for the user to place additiona files in this
    # directory.

    # If the file is changed, save the new version of it in the experiment
    # directory. Otherwise, depend on files in the cmc_input directory.
    # 1: must give filenames
    # 2: then the user must have a way of writing the new file somewhere.

    destination = os.path.join(parent_dir, name)

    if os.path.exists(destination):
        if not overwrite:
            raise Exception('Directory %s exists; not overwriting.' %
                    destination)
    else:
        os.makedirs(destination)

    # Copy over all setup files to destination.
    # -----------------------------------------
    # The user can modify these.  This way, we are not giving the user an
    # opportunity to edit the original setup files.
    if minimal:
        do_not_copy = ['desired_kinematics', 'force_plates',
                'extload_kinematics']
    else:
        do_not_copy = None
    orig_fpaths, exp_fpaths = dataman.copy_cmc_inputs(cmc_setup_fpath,
            destination, do_not_copy=do_not_copy)

    # Let the user change the files.
    # ------------------------------
    cmc_input = {
            'setup': exp_fpaths['setup'],
            'model': exp_fpaths['model'],
            'control_constraints': exp_fpaths['control_constraints'],
            'tasks': exp_fpaths['tasks'],
            'actuators': exp_fpaths['actuators']
            }
    # Get last-time-modified time stamps.
    mtimes = dict()
    for key, val in cmc_input.items():
        if val != None:
            mtimes[key] = os.path.getmtime(val)

    # Give the user a chance to edit these files.
    fcn(cmc_input)

    # Delete the files that didn't change, and properly update the setup files.
    # -------------------------------------------------------------------------
    # Necessarily need a setup file, and an external_loads file.
    if minimal:
        setup = etree.parse(cmc_input['setup'])
        tags = {'model': 'model_file', 'tasks': 'task_set_file', 'actuators':
                'force_set_files', 'control_constraints': 'constraints_file'}
    
        for key in ['model', 'control_constraints', 'tasks', 'actuators']:
            if (cmc_input[key] != None and os.path.getmtime(cmc_input[key]) <=
                    mtimes[key]):
                os.remove(cmc_input[key])
                setup.findall('.//%s' % tags[key])[0].text = \
                        os.path.relpath(orig_fpaths[key], destination)
        setup.write(exp_fpaths['setup'])

    # Write a README file.
    # --------------------
    readme = open(os.path.join(destination, 'README.txt'), 'a')
    readme.writelines("OpenSim CMC experiment '%s': %s.\n" % (name, description))
    readme.writelines("These files were generated on/at %s." %
            datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%MZ'))
    readme.close()

    if run_command:
        curdir = os.path.abspath(os.curdir())
        os.chdir(destination)
        try:
            os.system('%s -S %s' % (run_command, exp_cmc_setup_fpath))
        except:
            os.chdir(curdir)
            raise
        os.chdir(curdir)



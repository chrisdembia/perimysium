"""

To be used in Jython.

"""
import datetime
# TODO import difflib
import filecmp
import os
import xml.etree.ElementTree as etree

import dataman


def experiment(cmc_setup_fpath, parent_dir, name, description, fcn,
        minimal=True, run_command=False, overwrite=False):
    """

    cmc_input is 'setup', 'model', 'control_constraints', 'tasks', 'actuators'
    TODO if additional force set files are provided, update the correct field of cmc setup file; we'll take care of making sure the original actuators.xml is placed back in the field.

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
    # TODO clean up 

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

    # Write a README file.
    # --------------------
    readme = open(os.path.join(destination, 'README.txt'), 'a')
    readme.write("OpenSim CMC experiment '%s': %s\n" % (name, description))
    readme.write("These files were generated on/at %s.\n" %
            datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%MZ'))
    readme.write("This simulation is based on %s.\n" % cmc_setup_fpath)

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
    setup_fsf = etree.parse(exp_fpaths['setup']).findall(
            './/force_set_files')[0].text.strip()
    cmc_input = {
            'setup': exp_fpaths['setup'],
            'model': exp_fpaths['model'],
            'control_constraints': exp_fpaths['control_constraints'],
            'tasks': exp_fpaths['tasks'],
            'actuators': exp_fpaths['actuators']
            }

    # Give the user a chance to edit these files.
    cmc_input = fcn(cmc_input)

    print cmc_input

    # Delete the files that didn't change, and properly update the setup files.
    # -------------------------------------------------------------------------
    # Necessarily need a setup file, and an external_loads file.
    if minimal:
        readme.write("\n\nThe files that have changed:\n")
        setup = etree.parse(cmc_input['setup'])
        tags = {'model': 'model_file', 'tasks': 'task_set_file', 'actuators':
                'force_set_files', 'control_constraints': 'constraints_file'}

        setup_fields_to_edit = ['model', 'control_constraints', 'tasks']

        if type(cmc_input['actuators']) == str:
            # There's only one actuators file specified.
            setup_fields_to_edit += ['actuators']
        else:
            # Deal with the fact that there are multiple files specified.
            # Assume files are not renamed or deleted.
            newval = ''
            for i in range(len(orig_fpaths['actuators'])):
                if filecmp.cmp(orig_fpaths['actuators'][i],
                        cmc_input['actuators'][i]):
                    os.remove(cmc_input['actuators'][i])
                    newval += ' ' + os.path.relpath(orig_fpaths['actuators'][i],
                            destination)
                else:
                    newval += ' ' + cmc_input['actuators'][i]
            if len(orig_fpaths['actuators']) < len(cmc_input['actuators']):
                n_added = (len(cmc_input['actuators']) -
                        len(orig_fpaths['actuators']))
                n_prev = len(orig_fpaths['actuators'])
                for i in range(n_prev, n_prev + n_added):
                    newval += ' ' + os.path.relpath(cmc_input['actuators'][i],
                            destination)
            setup.findall('.//%s' % tags['actuators'])[0].text = newval
            del newval


        for key in setup_fields_to_edit:
            if cmc_input[key] != None:
                if (orig_fpaths[key] and
                        filecmp.cmp(orig_fpaths[key], cmc_input[key])):
                    # File unchanged. Delete the copy we made.
                    os.remove(cmc_input[key])
                    # Point the setup file to the original input file.
                    newval = os.path.relpath(orig_fpaths[key], destination)
                    if key == 'actuators':
                        if (setup.findall('.//force_set_files')[0].text.strip()
                                != setup_fsf):
                            # The `fcn` edited this field.
                            # In case other files are already specified, don't
                            # overwrite.
                            newval += setup.findall('.//%s' % tags[key])[0].text
                    setup.findall('.//%s' % tags[key])[0].text = newval
                else:
                    readme.write("%s --> %s\n" % (orig_fpaths[key],
                        cmc_input[key]))
                    # Put diff in the README.
                    # TODO
                    #readme.writelines(difflib.context_diff(
                    #    open(orig_fpaths[key]).readlines(),
                    #    open(cmc_input[key]).readlines(),
                    #    fromfile=orig_fpaths[key],
                    #    tofile=cmc_input[key]))
                    #readme.write('\n')
        setup.write(exp_fpaths['setup'])

    readme.write('\n')
    readme.close()

    if run_command:
        curdir = os.path.abspath(os.curdir)
        os.chdir(destination)
        try:
            os.system('%s -S %s' % (run_command, 
                os.path.basename(cmc_input['setup'])))
        except:
            os.chdir(curdir)
            raise
        os.chdir(curdir)



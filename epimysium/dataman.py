"""Manages the movement and use of data files."""

# TODO must edit cmc setup after it's moved, as well as grf.
# TODO put a comment in all files giving their original Hamner location? Not
# really necessary unless I share this code.
# TODO why is the model file so nested? there's definitely not that many
# changes that need to be made to the model file.
# TODO allow specifying which cycles to manage.

import abc
import csv
import difflib
import filecmp
import os
import shutil
import sys
import xml.etree.ElementTree as etree

import tables
import numpy as np

from cherithon import log
import swirl

def storage2numpy(storage_file):
    """Returns the data from a storage file in a numpy format. Skips all lines
    up to and including the line that says 'endheader'.

    Parameters
    ----------
    storage_file : str
        Path to an OpenSim Storage (.sto) file.

    Returns
    -------
    data : np.ndarry (or numpy structure array or something?)
        Contains all columns from the storage file, indexable by column name.

    Examples
    --------
    Columns from the storage file can be obtained as follows:
    
        >>> data = storage2numpy('<filename>')
        >>> data['ground_force_vy']

    """
    # What's the line number of the line containing 'endheader'?
    f = open(storage_file, 'r')

    for i, line in enumerate(f):
        if line.count('endheader') != 0:
            line_number_of_line_containing_endheader = i + 1
            break
    f.close()

    # With this information, go get the data.
    data = np.genfromtxt(storage_file, names=True,
            skip_header=line_number_of_line_containing_endheader)

    return data

def _splitall(path):
    """Splits a path into a list of the directories in the path. Copied from http://my.safaribooksonline.com/book/programming/python/0596001673/files/pythoncook-chp-4-sect-16.

    Parameters
    ----------
    path : str
        The path to split up.

    Returns
    -------
    allparts : list of str's
        One entry for each directory in the provided path; kept in the correct
        order.

    """
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts

def cmc_input_fpaths():
    """Assumes the input files refer to valid file paths (which copy_cmc_inputs
    does not assume).

    """
    pass

def copy_cmc_inputs(cmc_setup_fpath, destination, replace=None, **kwargs):
    """Given a CMC setup file, copies all files necessary to run CMC over to
    `destination`. All files necessary to run CMC are stored in the same
    directory. The CMC setup file and the external loads files are edited so
    that they refer to the correct copied files.

    Also sets the results_directory to be 'results'.

    TODO allow placing files in another location (not a flat structure); but
    then update the path in the file so that it's relative to the file.

    Parameters
    ----------
    cmc_setup_fpath : str
        Path to a CMC setup file.
    destination : str
        Directory in which to place the setup files.
    replace : dict, optional
        In case the paths in the files may be invalid, replace parts of the
        paths with the strings given in this dict. Keys are strings to look
        for, values are what to replace the key with.
    setup : str, optional
        A new filename for the cmc setup file (the first argument to this
        method).
    model : str, optional
        A new filename for this file.
    tasks : str, optional
        A new filename for this file.
    actuators : str, optional
        A new filename for this file.
    control_constraints : str, optional
        A new filename for this file.
    desired_kinematics : str, optional
        A new filename for this file.
    external_loads : str, optional
        A new filename for this file.
    force_plates : str, optional
        A new filename for this file.
    extload_kinematics : str, optional
        A new filename for this file.

    Returns
    -------
    new_fpaths : dict
        A valid filepath to the all the new files that were just copied over.
        The keys are as follows:
        - setup
        - model
        - tasks
        - actuators
        - control_constraints
        - desired_kinematics
        - external_loads
        - force_plates
        - extload_kinematics

    """
    fname = cmc_setup_fpath

    setup = etree.parse(fname)

    old = dict()

    def valid_path(file_containing_path, xml, tag):

        path = xml.findall('.//%s' % tag)[0].text
        if path == None or path.lstrip() == '': return None
        if os.path.exists(path): return path

        if replace:
            for key, val in replace.items():
                path = path.replace(key, val)
        if os.path.exists(path): return path

        path = path.replace('\\', '/')
        if os.path.exists(path): return path

        path2 = os.path.normpath(
                os.path.join(os.path.split(file_containing_path)[0], path))
        if os.path.exists(path2): return path2

        path2 = os.path.normpath(
                os.path.join(os.path.split(file_containing_path)[0], path.lstrip()))
        if os.path.exists(path2): return path2

        path2 = os.path.normpath(
                os.path.join(os.path.split(file_containing_path)[0], path.rstrip()))
        if os.path.exists(path2): return path2

        raise Exception('Path %s does not exist.' % path)

    # Get file names.
    # ---------------
    # Settings / parameters.
    old['model'] = valid_path(fname, setup, 'model_file')
    old['tasks'] = valid_path(fname, setup, 'task_set_file')
    old['actuators'] = valid_path(fname, setup, 'force_set_files')
    old['control_constraints'] = valid_path(fname, setup, 'constraints_file')

    # Data.
    old['desired_kinematics'] = valid_path(fname, setup,
            'desired_kinematics_file')
    old['external_loads'] = valid_path(fname, setup, 'external_loads_file')

    # Try to open the external loads file.
    extloads = etree.parse(old['external_loads'])
    old['force_plates'] = valid_path(old['external_loads'], extloads,
            'datafile')
    old['extload_kinematics'] = valid_path(old['external_loads'],
            extloads, 'external_loads_model_kinematics_file')

    # Copy files over.
    # ----------------
    # We'll store the location of the copies.
    new_fpaths = dict()
    new_fpaths['setup'] = None
    new_fpaths['model'] = None
    new_fpaths['tasks'] = None
    new_fpaths['actuators'] = None
    new_fpaths['control_constraints'] = None
    new_fpaths['desired_kinematics'] = None
    new_fpaths['external_loads'] = None
    new_fpaths['force_plates'] = None
    new_fpaths['extload_kinematics'] = None

    if not os.path.exists(destination): os.makedirs(destination)
    for key, val in old.items():
        if val and key != 'external_loads':
            if key in kwargs:
                new_fpath = os.path.join(destination, kwargs[key])
                shutil.copy(val, new_fpath)
                new_fpaths[key] = new_fpath
            else:
                shutil.copy(val, destination)
                new_fpaths[key] = os.path.join(destination, os.path.basename(val))

    # Edit the names of the files in the setup files.
    # -----------------------------------------------
    def edit_field(xml, tag, key):
        if old[key]:
            if key in kwargs:
                newvalue = kwargs[key]
            else:
                newvalue = os.path.basename(old[key])
            xml.findall('.//%s' % tag)[0].text = newvalue

    setup.findall('.//results_directory')[0].text = 'results'
    edit_field(setup, 'model_file', 'model')
    edit_field(setup, 'task_set_file', 'tasks')
    edit_field(setup, 'force_set_files', 'actuators')
    edit_field(setup, 'constraints_file', 'control_constraints')
    edit_field(setup, 'desired_kinematics_file', 'desired_kinematics')
    edit_field(setup, 'external_loads_file', 'external_loads')

    edit_field(extloads, 'datafile', 'force_plates')
    edit_field(extloads, 'external_loads_model_kinematics_file',
            'extload_kinematics')

    if 'setup' in kwargs:
        setup_new_fpath = os.path.join(destination, kwargs['setup'])
    else:
        setup_new_fpath = os.path.basename(cmc_setup_fpath)
    setup.write(setup_new_fpath)
    new_fpaths['setup'] = setup_new_fpath

    if 'external_loads' in kwargs:
        extloads_new_fpath = os.path.join(destination, kwargs['external_loads'])
    else:
        extloads_new_fpath = os.path.basename(cmc_setup_fpath)
    extloads.write(extloads_new_fpath)
    new_fpaths['external_loads'] = extloads_new_fpath

    return new_fpaths


def dock_output_in_pytable(h5file, output_path, group_path, allow_one=False,
        title=''):
    """Docks an OpenSim output, via a table for each STO file, in a pyTable
    file.  It's assumed the tables don't already exist in the last group
    specified.

    Parameters
    ----------
    h5file : tables.File
        pyTables File object, opened using tables.openFile(...). Does NOT close
        the file.
    output_path : str
        File path in which the OpenSim output is located (e.g. .STO files).
        Only .STO files are loaded, and it is assumed that all .STO files in
        this directory are from one run. That is, they have the same prefix,
        which is the name of the run.
    group_path : str, or list of str's
        The group tree hierarchy specifying where the output is to be docked in
        the h5file; as a path or as list of each directory's name (e.g.:
        'path/to/file' or ['path', 'to', 'file'])
    allow_one : bool, optional (default: False)
        Allows the loading of just one STO file. Otherwise, an exception is
        thrown. It is common that if only one STO file exists, it is a partial
        states file and means that the simulation did not complete.
    title : str, optional
        Title, in the pyTables file, for this group.

    Returns
    -------
    current_group : tables.Group
        The pyTables group in which the output has been stored.

    """
    # If output_path doesn't exist, can't do anything.
    if not os.path.exists(output_path):
        raise Exception("Output path {0:r} doesn't exist.".format(output_path))

    # Convert group_path to list of str's, if necessary.
    if type(group_path) == str:
        group_path = _splitall(group_path)

    # -- Make all necessary groups to get to where we're going.
    current_group = _blaze_group_trail(h5file, group_path, title)

    # -- Determine which files we want to use to create tables.

    # Make a list of all files in this directory ending is 'sto'.
    storage_files = [f for f in os.listdir(output_path) if f.endswith('.sto')]

    # If there are no storage files, the user probably gave a bad path.
    if len(storage_files) == 0:
        raise Exception("No .STO files found in {0}.".format(output_path))

    # If there's only one, usually the states file, forget about this output.
    if (not allow_one) and len(storage_files) == 1:
        raise Exception("Only one .STO file found: {0}.".format(
            storage_files[0]))


    # Get the length of the common prefix of these files.
    n_shared = _length_of_shared_prefix(storage_files)

    # -- Add tables in the current group.

    # Loop through all storage files.
    for f in storage_files:

        # Path to the data file.
        filepath = os.path.join(output_path, f)

        # Get name of the table: after the run name and before the file ext.
        if len(storage_files) == 1:
            table_name = os.path.splitext(f)[0]
        else:
            table_name = os.path.splitext(f)[0][n_shared:]

        # Create and populate a table with the data from this file.
        _populate_table(h5file, current_group, table_name, filepath)

    return current_group

def _blaze_group_trail(h5file, group_path, title=''):

    # Start at the root.
    current_group = h5file.root

    # Loop through each group in the path the user provided.
    for next_group_name in group_path:

        # This group has not been created yet.
        if not hasattr(current_group, next_group_name):

            # Create the group.
            h5file.createGroup(current_group, next_group_name, title=title)

        # Set this group as the current group.
        current_group = getattr(current_group, next_group_name)

    # Return this leaf group.
    return current_group


def _populate_table(h5file, group, table_name, filepath):
    """Populates a pyTables file with a table, using data from CSV file at
    filepath.

    Parameters
    ----------
    h5file : tables.File
        The file to which the table is to be added.
    group : tables.Group
        The group in the file to which the table is to be added.
    table_name : str
        Name of the table to be added.
    filepath : str
        Path to the data file containing data to put in this table.

    Returns
    -------
    table : tables.Table
        The table that has just been created.

    """
    # Open data file.
    csvread = csv.reader(open(filepath, 'r'), delimiter='\t',
            skipinitialspace=True)

    # - Parse the data as a CSV file.
    do_parse = False
    take_header = False

    # For each row in the CSV file.
    for csvrow in csvread:

        # If this is the header line.
        if take_header:

            # Save this row for later; we'll need it.
            title_row = csvrow

            # Can't have periods in table column names in pyTables.
            for i in range(len(title_row)):
                title_row[i] = title_row[i].replace('.', '_')

            take_header = False

            # Grab table column names.
            table_cols = dict()
            for col in title_row: 
                # Checking if the column is empty. This is a
                # once-in-a-blue-moon bug fix as a result of inconsistency in
                # Hamner's files. See CMC results for subject 2, speed 2 m/s,
                # cycle 1, states_OG.sto file.
                if col != '':
                    table_cols[col] = tables.Float32Col()

            # Create pyTables table.
            table = h5file.createTable(group, table_name, table_cols,
                    'Output file {0}'.format(table_name))

        # If this is a data row.
        elif do_parse:

            # For each column in the data file in this row.
            for i in range(len(table_cols.keys())):

                # Append the data into the table.
                table.row[title_row[i]] = csvrow[i]

            # Tell pyTables to append this data to the table.
            table.row.append()

        # The header is over.
        if csvrow == ['endheader']:
            take_header = True
            do_parse = True

    # Save (?).
    table.flush()

    # Give access to it.
    return table

def _length_of_shared_prefix(strings):
    """Determines, from a list of strings, the length of the string that is
    shared at the beginning of all strings.

    Parameters
    ----------
    strings : list of str's
        List of strings to compare.

    Returns
    -------
    n_shared : int
        The number of characters shared at the beginning of all strings.

    """
    # Initialize the number of shared caracters to something equal to or
    # greater than what it will finally be.
    n_shared = len(strings[0])
    for i_string in range(1, len(strings)):

        # Compare the 0th and 1st strings.
        diff = difflib.SequenceMatcher(a=strings[0], b=strings[i_string])

        # The 3rd element of the first element of matching blocks tells
        # how many characters these strings share.
        n_shared = min(n_shared, diff.get_matching_blocks()[0][2])

    return n_shared

class Log(log.Log):
    """ TODO """

    def __init__(self, fname):
        """ TODO """
        self.fname = fname

    def entry_new(self, priority, cond, obj=None, func_name=""):
        fid = open(self.fname, 'a')
        if obj != None: fid.write("'" + obj.name + "'.")
        if func_name != "": fid.write(func_name + "():")
        fid.write(cond + '\n')
        fid.close()


class Dataman(object):
    """Abstract base class for managing data."""

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass


class ChandMLwalk(Dataman):
    """Mediolateral analysis of walking."""
    pass


class HamnerXX(Dataman):
    """Manages Hamner's 20-subject running data. XX is for 20 in Roman
    numerals.

    """
    expdir = os.path.join("hpl_data", "hamner_test")
    osimdir = "opensimresults"
    cmc_actuators_from = "gait2392_cmc_actuators.xml"
    cmc_tasks_from = "gait2392_cmc_tasks.xml"
    cmc_actuators_to = "cmc_actuators.xml"
    cmc_tasks_to = "cmc_tasks.xml"
    cop_from = os.path.join("session 2", "exporteddata",
            "run_%i00 02_newcop3_v24.mot")
    cop_to = "cop.mot"
    # TODO move from-names from Cycle to here.

    grf_to = "grf.xml"
    concon_to = "controlcontraints.xml"
    cmc_setup_to = "cmc_setup.xml"
    kin_to = "rra_kinematics_q.xml"
    model_to = "rra_model.osim"

    n_total_subjects = 20

    rra_output_tables = ['actuation_force',
            'actuation_power',
            'actuation_speed',
            'controls',
            'kinematics_dudt',
            'kinematics_q',
            'kinematics_u',
            'perr',
            'states']

    cmc_output_tables = ['Actuation_force',
            'Actuation_power',
            'Actuation_speed',
            'Kinematics_dudt',
            'Kinematics_q',
            'Kinematics_u',
            'controls',
            'pErr',
            'states']


    class CopyMode:
        MINIMAL=1
        CONSERVATIVE=2

    def __init__(self, cmc_exec, from_path, h5fname, subject_idxs, speed_idxs):
        """
        Parameters
        ----------
        cmc_exec : str
            Path to cmc executable.
        from_path : str
            Something that probably ends in 'Hamner_Data'.
        h5fname : str
            Name of the pyTables/HDF5 file that data is stored in.
        n_subjects : int
            The number of subjects to copy over. Max is
            HamnerXX.n_total_subjects.

        """
        self.log = Log("hamnerxxlog.txt")
        self.cmc_exec = cmc_exec
        self._copy_mode = self.CopyMode.CONSERVATIVE
        self.from_path = from_path

        # Create h5 file if it doesn't exist.
        self.h5fname = h5fname
        h5 = tables.openFile(h5fname, mode='a', title='HamnerXX')
        self.n_subjects = len(subject_idxs)
        if max(subject_idxs) > self.n_total_subjects:
            raise Exception("Greatest subject is {0}; less than {1}.".format(
                self.n_total_subjects, self.n_subjects))
        self._subjects = []
        for idx in subject_idxs:
            self._subjects += [Subject(self, idx, speed_idxs)]
        h5.close()

    @property
    def copy_mode(self): return self._copy_mode

    @copy_mode.setter
    def copy_mode(self, mode): self._copy_mode = mode

    def _actuallycopy(self, pairs):
        """ TODO """
        # Loop through the pairs of files, checking for errors and moving
        # them.
        for key, val in pairs.iteritems():

            # Local copy for the file to copy, and where to copy it to.
            fro = val[0]
            to = val[1]

            try:

                # If 'to' is there, compare it to 'from'.
                if (os.path.exists(to) and filecmp.cmp(fro, to)):
                    self.log.entry_new(Log.ERROR(), "[" + fro +
                            "] and [" + to +
                            "] are not the same; overwriting.")

                    # Try copying the file.
                shutil.copy2(fro, to)

            except (IOError, OSError) as e:
                self.log.entry_new(Log.ERROR(), e.strerror + " " + e.filename)

    def working_copy_cmc_new(self, to_path):
        """Copies and renames files necessary for CMC. Files are stored in a
        hierarchy that reveals the extent to which files are shared. Files are
        also modified so that a CMC run can actually be performed after all the
        copying. The options int he CMC setup files are not modified;
        generating a setup from scratch wouldn't be the responsibility of a
        Dataman.

        The files necessary for CMC are:

        Given organization:

            0.  osim/subject/cmc/[cmc setup template xml].
            1.  osim/subject/cmc/[cmc actuators xml].
            2.  osim/subject/cmc/[cmc tasks xml].
            3.  osim/subject/rra/speed/[cycle - model file osim].
            4.  osim/subject/rra/speed/cycle/[kinematics sto].
            5.  osim/subject/cmc/speed/[cycle - cmc setup xml].
            6.  osim/subject/cmc/speed/[cycle - grf xml].
            7.  osim/subject/cmc/speed/[cycle - control constraints xml].
            8.  exp/subject/session/exported/[speed - cop mot].

        Minimal organization:

            0.  subject/speed/cycle/[cmc setup xml].
            1.  subject/speed/cycle/[model file osim].
            2.  [cmc actuators xml].
            3.  [cmc tasks xml].
            4.  subject/speed/cycle/[kinematics sto].
            5.  subject/speed/cycle/[grf xml].
            6.  subject/speed/cycle/[control constraints xml].
            7.  subject/speed/[cop mot].

        Conservative organization:

            0.  subject/speed/cycle/simtemplate/[cmc setup xml].
            1.  subject/speed/cycle/simtemplate/[model file osim].
            2.  subject/speed/cycle/simtemplate/[cmc actuators xml].
            3.  subject/speed/cycle/simtemplate/[cmc tasks xml].
            4.  subject/speed/cycle/simtemplate/[kinematics sto].
            5.  subject/speed/cycle/simtemplate/[grf xml].
            6.  subject/speed/cycle/simtemplate/[control constraints xml].
            7.  subject/speed/cycle/simtemplate/[cop mot].

        This necessitates calls to Subject, Speed, and Cycle instances. The
        code for a specific file to be copied is located in the most-nested
        class necessary, whether that's from the copying-from location or the
        copying-to location.

        Parameters
        ----------
        to_path : str
            Probably some local directory on a machine (probably not in a
            repository or Dropbox) where the simulations are to be run.

        """
        # If to_path doesn't exist, make it!
        if os.path.exists(to_path):
            raise Exception("Path {0!r} already exists.".format(to_path))
        os.mkdir(to_path)

        # Only copies files over.
        for subj in self._subjects:
            subj._copy_cmc_new(to_path)

        # Edits the cmc setup files to account for the new file names.
        for subj in self._subjects:
            subj._modify_cmc_setups(to_path)

    def cmc_run_new(self, to_path, name):
        """Run CMC for all subjects/speeds/cycles specified. Stores results in
        a pyTables/HDF5 file.

        Parameters
        ----------
        to_path : str
            The to_path to which a working copy was placed.
        name : str
            The name of the folder, placed within each 'cycle' folder, in which
            to place simulation results.

        """
        h5 = tables.openFile(self.h5fname, mode='a')
        invest_group = h5.createGroup(h5.root, name,
                'Investigation {0}'.format(name))
        for subj in self._subjects:
            subj._cmc_run_new(to_path, name, h5, invest_group)
        h5.close()

    def store_rra_output(self):
        """Stores the output of Hamner's RRA simulations into the pyTables/HDF5
        file. No copying necessary.

        """
        h5 = tables.openFile(self.h5fname, mode='a')
        try:
            base_group = h5.createGroup(h5.root, 'base', "Hamner's original data.")
        except tables.NodeError as e:
            base_group = h5.getNode(h5.root, 'base')
        for subj in self._subjects:
            subj._store_rra_output(h5, base_group)
        h5.close()

    def torque_comparison_report(self, run_name, fname, stiffness):
        """Creates a PDF report comparing, for each cyle, hip moment to
        associated spring torque.

        Parameters
        ----------
        run_name : str
            The name of the cmc_run from which to extract results.
        fname : str
            Name/path of the PDf to generate.
        stiffness : float
            Stiffness of torsional spring, in units of N-m/rad.

        """
        report = swirl.FigureReport(fname)

        report.front_matter("Torque comparison report\\\\"
                "Hip spring assistive device\\\\"
                "Neuromuscular Biomechanics Laboratory\\\\"
                "Stanford University", "Chris Dembia, Carmichael Ong")

        for subj in self._subjects:
            subj._torque_comparison_report(run_name, report, stiffness)

        report.next_is("\\end{document}")

        pathparts = os.path.split(fname)
        report.write(pathparts[0], pathparts[1])

class Subject(object):

    name_pre = "subject"
    cmc_dirname = "cmc_multiplesteps_v24_tendon_040_vmax_15_passive_10_2x/"
    rra_dirname = "rra_multiplesteps"
    speed_map = {1: 2, 2: 3, 3: 4, 4: 5}

    def __init__(self, hamner, index, speed_idxs):
        """ TODO """
        self.hamner = hamner
        self.name = self.name_pre + "{:02g}".format(index)
        self.from_dir = os.path.join(hamner.osimdir, self.name)
        self.from_cmc_dir = os.path.join(self.from_dir, self.cmc_dirname)
        self.from_rra_dir = os.path.join(self.from_dir, self.rra_dirname)
        self._speeds = []
        for idx in speed_idxs:
            self._speeds += [Speed(hamner, self, idx, self.speed_map[idx])]

    def _copy_cmc_new(self, to_path):
        """ TODO """

        # For convenience.
        from_path = self.hamner.from_path
        # If path doesn't exist, make it!
        subj_dir = os.path.join(to_path, self.name)
        if not os.path.exists(subj_dir): os.mkdir(subj_dir)

        if self.hamner.copy_mode == self.hamner.CopyMode.MINIMAL:

            # TODO do diff's to see if files actually differ, precluding the
            # minimal copying.

            # Prepare filenames for movement.
            actuators_from = os.path.join(from_path, self.hamner.osimdir,
                    self.name, self.cmc_dirname,
                    self.hamner.cmc_actuators_from)
            actuators_to = os.path.join(to_path, self.hamner.cmc_actuators_to)
            tasks_from = os.path.join(from_path, self.hamner.osimdir,
                    self.name, self.cmc_dirname, self.hamner.cmc_tasks_from)
            tasks_to = os.path.join(to_path, self.hamner.cmc_tasks_to)

            # Package for looping.
            pairs = {
                    'actuators': (actuators_from, actuators_to),
                    'tasks': (tasks_from, tasks_to)
                    }

            self.hamner._actuallycopy(pairs)

        for speed in self._speeds:
            speed._copy_cmc_new(to_path)

    def _modify_cmc_setups(self, to_path):
        """Modifies pre-existing CMC setup files. This action depends on
        whether copy mode is set to be MINIMAl or CONSERVATIVE.

        """
        for speed in self._speeds:
            speed._modify_cmc_setups(to_path)

    def _cmc_run_new(self, to_path, name, h5, invest_group):
        subj_group = h5.createGroup(invest_group, self.name, self.name)
        for speed in self._speeds:
            speed._cmc_run_new(to_path, name, h5, subj_group)

    def _store_rra_output(self, h5, base_group):
        try:
            subj_group = h5.createGroup(base_group, self.name, self.name)
        except tables.NodeError as e:
            subj_group = h5.getNode(base_group, self.name)
        for speed in self._speeds:
            speed._store_rra_output(h5, subj_group)

    def _torque_comparison_report(self, run_name, report, stiffness):
        report.next_is("\\section{%s}" % self.name)
        for speed in self._speeds:
            speed._torque_comparison_report(run_name, report, stiffness)


class Speed(object):

    name_pre = "speed"
    speed_cmc_dirname = "cmc_results_v240_run_%i0002"
    speed_rra_dirname = "rra_results_v191_run_%i0002"
    cycle_cmc_from_dir = "cmc_results_v240_run_%i0002_cycle%i"
    cycle_rra_from_dir = "rra_results_v191_run_%i0002_cycle%i"

    def __init__(self, hamner, subject, index, speed):
        self.hamner = hamner
        self.subject = subject
        self.index = index
        self.name = self.name_pre + "{:02g}".format(index)
        self.speed = speed
        self.cmc_from_dir = os.path.join(subject.from_cmc_dir, 
                self.speed_cmc_dirname % speed)
        self.rra_from_dir = os.path.join(subject.from_rra_dir, 
                self.speed_rra_dirname % speed)
        # Get the number of cycles.
        self._cycles = []
        keepGoing = True
        n_cycles = 0
        while keepGoing:
            if os.path.exists(os.path.join(hamner.from_path, self.cmc_from_dir,
                self.cycle_cmc_from_dir % (self.speed, n_cycles+1))):
                n_cycles += 1
                self._cycles += [Cycle(hamner, subject, self, n_cycles)]
            else:
                keepGoing = False

    def _copy_cmc_new(self, to_path):
        """ TODO """

        # For convenience.
        from_path = self.hamner.from_path
        # If path doesn't exist, make it!
        speed_dir = os.path.join(to_path, self.subject.name, self.name)
        if not os.path.exists(speed_dir): os.mkdir(speed_dir)

        if self.hamner.copy_mode == self.hamner.CopyMode.MINIMAL:

            # Prepare filenames for movement.
            cop_from = os.path.join(from_path, self.hamner.expdir,
                    self.subject.name, self.hamner.cop_from % self.speed)
            cop_to = os.path.join(to_path, self.subject.name, self.name,
                    self.hamner.cop_to)

            # Package for looping.
            pairs = {
                    'cop': (cop_from, cop_to),
                    }

            self.hamner._actuallycopy(pairs)

        for cycle in self._cycles:
            cycle._copy_cmc_new(to_path)

    def _modify_cmc_setups(self, to_path):
        """Modifies pre-existing CMC setup files. This action depends on
        whether copy mode is set to be MINIMAl or CONSERVATIVE.

        """
        for cycle in self._cycles:
            cycle._modify_cmc_setups(to_path)

    def _cmc_run_new(self, to_path, name, h5, subj_group):
        speed_group = h5.createGroup(subj_group, self.name, self.name)
        for cycle in self._cycles:
            cycle._cmc_run_new(to_path, name, h5, speed_group)

    def _store_rra_output(self, h5, subj_group):
        try:
            speed_group = h5.createGroup(subj_group, self.name, self.name)
        except tables.NodeError as e:
            speed_group = h5.getNode(subj_group, self.name)
        for cycle in self._cycles:
            cycle._store_rra_output(h5, speed_group)

    def _torque_comparison_report(self, run_name, report, stiffness):
        report.next_is("\\subsection{Speed %i m/s}" % self.speed)
        h5 = tables.openFile(self.hamner.h5fname)
        try:

            # Initialize.
            act_time = np.array([])
            hip_flexion_moment_left = np.array([])
            hip_flexion_moment_right = np.array([])

            # Prep for calculating spring torque.
            kin_time = np.array([])
            hip_flexion_left = np.array([])
            hip_flexion_right = np.array([])

            for cycle in self._cycles:
                # Group with the necessary data.
                rra = h5.getNode(h5.getNode(h5.getNode(h5.root.base,
                    self.subject.name), self.name), cycle.name).rra

                # Actual hip flexion moments.
                # Data is so dense, we only need every 10th point.
                act_time = np.concatenate((act_time, 
                    rra.actuation_force.cols.time[::10]))
                hip_flexion_moment_left = np.concatenate((hip_flexion_moment_left,
                    rra.actuation_force.cols.hip_flexion_l[::10]))
                hip_flexion_moment_right = np.concatenate((hip_flexion_moment_right,
                    rra.actuation_force.cols.hip_flexion_r[::10]))

                # Prep for calculating spring torque.
                kin_time = np.concatenate((kin_time,
                    rra.kinematics_q.cols.time[::10]))
                hip_flexion_left = np.concatenate((hip_flexion_left,
                    rra.kinematics_q.cols.hip_flexion_l[::10]))
                hip_flexion_right = np.concatenate((hip_flexion_right,
                    rra.kinematics_q.cols.hip_flexion_r[::10]))

            # Actual calculation.
            stiffness_deg = stiffness * np.pi / 180.0
            spring_torque = -stiffness_deg * (hip_flexion_right - hip_flexion_left)

            # Create the plot.
            report.simple_data_plot(
                    'xlabel=time (s), ylabel=moment (N-m), no markers, '
                    'legend pos=outer north east, '
                    'width=6in, '
                    'height=2in',
                [
                    'left hip flexion',
                    'right hip flexion',
                    'associated hip spring'
                    ],
                [
                    act_time, act_time,
                    kin_time
                    ],
                [
                    hip_flexion_moment_left,
                    hip_flexion_moment_right,
                    spring_torque
                    ])
        except tables.NoSuchNodeError as e:
            report.next_is("Data not available.")

        h5.close()
        for cycle in self._cycles:
            cycle._torque_comparison_report(run_name, report, stiffness)


class Cycle(object):

    name_pre = "cycle"

    def __init__(self, hamner, subject, speed, index):
        self.hamner = hamner
        self.subject = subject
        self.speed = speed
        self.cycle_idx = index
        self.name = self.name_pre + "{:02g}".format(index)

        self.grf_from = "%s_run_%i0002_cycle%i_grf_v240.xml" % (
                self.subject.name, self.speed.speed, self.cycle_idx)
        self.concon_from = "%s_run_%i0002_cycle%i_controlconstraints.xml" % (
                self.subject.name, self.speed.speed, self.cycle_idx)
        self.cmc_setup_from = "%s_setup_cmc_run_%i0002_cycle%i_v240.xml" % (
                self.subject.name, self.speed.speed, self.cycle_idx)
        self.kin_from = "%s_run_%i0002_cycle%i_kinematics_q.sto" % (
                self.subject.name, self.speed.speed, self.cycle_idx)
        self.model_from = ("%s_rra_adjusted_run_%i0002_cycle%i_v191_" +
                "tendon_040_vmax_15_passive_10_2x.osim") % (
                    self.subject.name, self.speed.speed, self.cycle_idx)

    def _copy_cmc_new(self, to_path):
        """ TODO """

        # For convenience.
        from_path = self.hamner.from_path
        # If path doesn't exist, make it!
        cycle_dir = os.path.join(to_path, self.subject.name, self.speed.name,
                self.name)
        if not os.path.exists(cycle_dir): os.mkdir(cycle_dir)

        # Prepare filenames for movement.
        grf_from = os.path.join(from_path, self.speed.cmc_from_dir,
                self.grf_from)
        concon_from = os.path.join(from_path, self.speed.cmc_from_dir,
                self.concon_from)
        cmc_setup_from = os.path.join(from_path,
                self.speed.cmc_from_dir, self.cmc_setup_from)
        kin_from = os.path.join(from_path, self.speed.rra_from_dir,
                self.speed.cycle_rra_from_dir % (
                    self.speed.speed, self.cycle_idx),
                self.kin_from)
        model_from = os.path.join(from_path, self.speed.rra_from_dir,
                self.model_from)

        grf_to = os.path.join(to_path, self.subject.name, self.speed.name,
                self.name, self.hamner.grf_to)
        concon_to = os.path.join(to_path, self.subject.name,
                self.speed.name, self.name, self.hamner.concon_to)
        cmc_setup_to = os.path.join(to_path, self.subject.name,
                self.speed.name, self.name, self.hamner.cmc_setup_to)
        kin_to = os.path.join(to_path, self.subject.name,
                self.speed.name, self.name, self.hamner.kin_to)
        model_to = os.path.join(to_path, self.subject.name,
                self.speed.name, self.name, self.hamner.model_to)

        # Package for looping.
        pairs = {
                'grf': (grf_from, grf_to),
                'concon': (concon_from, concon_to),
                'cmc_setup': (cmc_setup_from, cmc_setup_to),
                'kin': (kin_from, kin_to),
                'model': (model_from, model_to),
                }

        self.hamner._actuallycopy(pairs);

        if self.hamner.copy_mode == self.hamner.CopyMode.CONSERVATIVE:
            # These are the files that need to be duplicated.
            # Subject-level: actuators, tasks.
            # Speed-level: cop

            actuators_from = os.path.join(from_path, self.hamner.osimdir,
                    self.subject.name, self.subject.cmc_dirname,
                    self.hamner.cmc_actuators_from)
            tasks_from = os.path.join(from_path, self.hamner.osimdir,
                    self.subject.name, self.subject.cmc_dirname,
                    self.hamner.cmc_tasks_from)
            cop_from = os.path.join(from_path, self.hamner.expdir,
                    self.subject.name, self.hamner.cop_from % self.speed.speed)

            actuators_to = os.path.join(to_path, self.subject.name,
                    self.speed.name, self.name, self.hamner.cmc_actuators_to)
            tasks_to = os.path.join(to_path, self.subject.name,
                    self.speed.name, self.name, self.hamner.cmc_tasks_to)
            cop_to = os.path.join(to_path, self.subject.name, self.speed.name,
                    self.name, self.hamner.cop_to)

            # Package for looping.
            pairs = {
                    'actuators': (actuators_from, actuators_to),
                    'tasks': (tasks_from, tasks_to),
                    'cop': (cop_from, cop_to)
                    }

            self.hamner._actuallycopy(pairs)

    def _modify_cmc_setups(self, to_path):

        # --- Modify the CMC file.
        fname = os.path.join(to_path, self.subject.name, self.speed.name,
                self.name, self.hamner.cmc_setup_to)
        self._modify_cmc_setup_impl(fname)

        # --- Modify the GRF file.
        grf_fname = os.path.join(to_path, self.subject.name,
                self.speed.name, self.name, self.hamner.grf_to)
        self._modify_grf_impl(grf_fname)

    def _modify_cmc_setup_impl(self, fname, run_name=None, prepend_path='.',
            results_dir='results'):
        """Modifies pre-existing CMC setup files. This action depends on
        whether copy mode is set to be MINIMAl or CONSERVATIVE.

        """
        # Read into an ElementTree for parsing and modification.
        setup = etree.parse(fname)

        # -- name in output files.
        name = setup.findall('.//CMCTool')
        # Error check.
        if len(name) != 1: self._xml_find_raise('CMCTool', fname)
        # Change the entry.
        newname = '{0}_{1}_{2}'.format(self.subject.name, self.speed.name,
                self.name)
        if run_name != None:
            newname = '{0}_{1}'.format(run_name, newname)
        name[0].attrib['name'] = newname

        # -- model_file.
        mf = setup.findall('.//model_file')
        # TODO consider logging instead.
        if len(mf) != 1: self._xml_find_raise('model_file', fname)
        # Change the entry.
        mf[0].text = os.path.join(prepend_path, self.hamner.model_to)

        # -- results_directory
        rd = setup.findall('.//results_directory')
        # Error check.
        if len(rd) != 1: self._xml_find_raise('results_directory', fname)
        # Change the entry: place holder.
        rd[0].text = os.path.join(prepend_path, results_dir)

        # -- external_loads_file
        ext = setup.findall('.//external_loads_file')
        # Error check.
        if len(ext) != 1: self._xml_find_raise('external_loads_file', fname)
        # Change the entry.
        ext[0].text = os.path.join(self.hamner.grf_to)
        # TODO inconsistent naming btn external loads and grf.

        # -- desired_kinematics_file
        kin = setup.findall('.//desired_kinematics_file')
        # Error check.
        if len(kin) != 1: self._xml_find_raise('desired_kinematics_file',
                fname)
        # Change the entry.
        kin[0].text = os.path.join(prepend_path, self.hamner.kin_to)

        # -- constraints_file
        concon = setup.findall('.//constraints_file')
        # Error check.
        if len(concon) != 1: self._xml_find_raise('constraints_file', fname)
        # Change the entry.
        concon[0].text = os.path.join(prepend_path, self.hamner.concon_to)

        if self.hamner.copy_mode == self.hamner.CopyMode.MINIMAL:
            # -- force_set_files
            act = setup.findall('.//force_set_files')
            # Error-check.
            if len(act) != 1: self._xml_find_raise('force_set_files', fname)
            # Change the entry: file located in to_path
            act[0].text = os.path.join(prepend_path, '..', '..', '..',
                    self.hamner.cmc_actuators_to)

            # -- task_set_file
            task = setup.findall('.//task_set_file')
            # Error check.
            if len(task) != 1: self._xml_find_raise('task_set_file', fname)
            # Change the entry: file located in to_path
            task[0].text = os.path.join(prepend_path, '..', '..', '..',
                    self.hamner.cmc_tasks_to)

        elif self.hamner.copy_mode == self.hamner.CopyMode.CONSERVATIVE:
            # -- force_set_files
            act = setup.findall('.//force_set_files')
            # Error-check.
            if len(act) != 1: self._xml_find_raise('force_set_files', fname)
            # Change the entry: file located in to_path
            act[0].text = os.path.join(prepend_path,
                    self.hamner.cmc_actuators_to)

            # -- task_set_file
            task = setup.findall('.//task_set_file')
            # Error check.
            if len(task) != 1: self._xml_find_raise('task_set_file', fname)
            # Change the entry: file located in to_path
            task[0].text = os.path.join(prepend_path, self.hamner.cmc_tasks_to)

        # Save the modified file.
        setup.write(fname)

    def _modify_grf_impl(self, fname, prepend_path='.'):
        """Modifies a pre-existing GRF file. This action depends on
        whether copy mode is set to be MINIMAl or CONSERVATIVE.

        """
        # Read into an ElementTree for parsing and modification.
        grf_tree = etree.parse(fname)

        if self.hamner.copy_mode == self.hamner.CopyMode.MINIMAL:
            # -- datafile
            dat = grf_tree.findall('.//datafile')
            # Error check.
            if len(dat) != 1: self._xml_find_raise('datafile', grf_fname)
            # Change the entry.
            dat[0].text = os.path.join(prepend_path, '..', self.hamner.cop_to)

        elif self.hamner.copy_mode == self.hamner.CopyMode.CONSERVATIVE:
            # -- datafile
            dat = grf_tree.findall('.//datafile')
            # Error check.
            if len(dat) != 1: self._xml_find_raise('datafile', grf_fname)
            # Change the entry.
            dat[0].text = os.path.join(prepend_path, self.hamner.cop_to)

        # Save the file.
        grf_tree.write(fname)

    @staticmethod
    def _xml_find_raise(tag, fname):
        raise Exception("XML tag '{0}' doesn't exist or occurs multiple " +
                "times in file {1}".format(tag, fname))

    def _cmc_run_new(self, to_path, name, h5, speed_group):
        """Does the actual execution."""
        # TODO belongs in the hipspring package.

        # This cycle folder.
        cycle_path = os.path.join(to_path, self.subject.name, self.speed.name,
                self.name)

        # --- Manage inputs.
        input_dirname = "input_{0}".format(name)
        input_path = os.path.join(cycle_path, input_dirname)
        if os.path.exists(input_path):
            pass
            #raise Exception("A run with name '{0}' may already exist. " +
            #        "Not overwriting.".format(name))
        else:
            # Create the input directory.
            os.mkdir(input_path)

        # Copy over the default setup file to the input directory.
        cmc_setup_path = os.path.join(input_path, self.hamner.cmc_setup_to)
        shutil.copy2(os.path.join(cycle_path, self.hamner.cmc_setup_to),
                cmc_setup_path)

        # TODO also need grf.xml in the current directory.
        shutil.copy2(os.path.join(cycle_path, self.hamner.grf_to), input_path)

        # Edit the setup file.
        results_dirname = "output_{0}".format(name)
        results_path = os.path.join(cycle_path, results_dirname)
        if os.path.exists(results_path):
            pass
            #raise Exception(("A run with name '{0}' may already exist. " +
            #        "Not overwriting.").format(name))

        # Reflect file movements in the cmc setup file.
        #self._modify_cmc_setup_impl(cmc_setup_path, prepend_path='..',
        #        results_dir=results_dirname)
        self._modify_cmc_setup_impl(cmc_setup_path, name,
                prepend_path=os.path.join(input_path, '..'),
                results_dir=results_dirname)

        # --- Modifying GRF file's knowledge of where cop.xml is.
        # TODO also need grf.xml in the current directory.
        self._modify_grf_impl(os.path.join(input_path, self.hamner.grf_to),
                prepend_path='..')

        # Create the results directory.
        #os.mkdir(results_path)

        # Run CMC.
        os.chdir(results_path)
        #os.system("{0} -S {1}".format(self.hamner.cmc_exec,
        #    os.path.join('..', input_dirname, self.hamner.cmc_setup_to)))

        # Parse output and save an HDF5 file.
        cycle_group = h5.createGroup(speed_group, self.name, self.name)
        
        for tablename in self.hamner.cmc_output_tables:
            filepath = os.path.join(results_path,
                    '{0}_{1}_{2}_{3}_{4}.sto'.format(name, self.subject.name,
                        self.speed.name, self.name, tablename))
            csvread = csv.reader(open(filepath, 'r'), delimiter='\t',
                    skipinitialspace=True)

            do_parse = False
            take_header = False
            for csvrow in csvread:
                if take_header:
                    title_row = csvrow
                    take_header = False
                    table_cols = dict()
                    for col in title_row:
                        table_cols[col.replace('.', '_')] = tables.Float32Col()
                    table = h5.createTable(cycle_group, tablename, table_cols,
                            'Output file {0}'.format(tablename))
                elif do_parse:
                    tablerow = table.row
                    for i in range(len(table_cols.keys())):
                        tablerow[title_row[i].replace('.', '_')] = csvrow[i]
                    tablerow.append()
                if csvrow == ['endheader']:
                    take_header = True
                    do_parse = True
            table.flush()

    def _store_rra_output(self, h5, speed_group):
        try:
            cycle_group = h5.createGroup(speed_group, self.name, self.name)
        except tables.NodeError as e:
            cycle_group = h5.getNode(speed_group, self.name)
        try:
            rra_group = h5.createGroup(cycle_group, 'rra', "Hamner's RRA output")
        except tables.NodeError as e:
            rra_group = h5.getNode(cycle_group, 'rra')

        output_templ = "%s_run_%i0002_cycle%i_{0}.sto" % (
                self.subject.name, self.speed.speed, self.cycle_idx)
        for tablename in self.hamner.rra_output_tables:
            filepath = os.path.join(self.hamner.from_path,
                    self.speed.rra_from_dir,
                    self.speed.cycle_rra_from_dir % (self.speed.speed,
                        self.cycle_idx),
                    output_templ.format(tablename))
            csvread = csv.reader(open(filepath, 'r'), delimiter='\t',
                    skipinitialspace=True)

            try:
                do_parse = False
                take_header = False
                for csvrow in csvread:
                    if take_header:
                        title_row = csvrow
                        take_header = False
                        table_cols = dict()
                        for col in title_row:
                            table_cols[col.replace('.', '_')] = tables.Float32Col()
                            table = h5.createTable(rra_group, tablename, table_cols,
                                    'Output file {0}'.format(tablename))
                    elif do_parse:
                        tablerow = table.row
                        for i in range(len(table_cols.keys())):
                            tablerow[title_row[i].replace('.', '_')] = csvrow[i]
                        tablerow.append()
                    if csvrow == ['endheader']:
                        take_header = True
                        do_parse = True
                table.flush()
            except tables.NodeError as e:
                # TODO If the table already exists, don't overwrite it.
                pass

    def _torque_comparison_report(self, run_name, report, stiffness):
        h5 = tables.openFile(self.hamner.h5fname)
        report.next_is("\\subsubsection{Cycle %i}" % self.cycle_idx)
        try:

            # Group with the necessary data.
            rra = h5.getNode(h5.getNode(h5.getNode(h5.root.base,
                self.subject.name), self.speed.name), self.name).rra

            # Actual hip flexion moments.
            # Data is so dense, we only need every 10th point.
            act_time = rra.actuation_force.cols.time[::10]
            hip_flexion_moment_left = rra.actuation_force.cols.hip_flexion_l[::10]
            hip_flexion_moment_right = rra.actuation_force.cols.hip_flexion_r[::10]

            # Prep for calculating spring torque.
            kin_time = rra.kinematics_q.cols.time[::10]
            hip_flexion_left = rra.kinematics_q.cols.hip_flexion_l[::10]
            hip_flexion_right = rra.kinematics_q.cols.hip_flexion_r[::10]

            # Actual calculation.
            stiffness_deg = stiffness * np.pi / 180.0
            spring_torque = -stiffness_deg * (hip_flexion_right - hip_flexion_left)

            # Create the plot.
            report.simple_data_plot(
                    'xlabel=time (s), ylabel=moment (N-m), no markers, '
                    'legend pos=outer north east',
                [
                    'left hip flexion',
                    'right hip flexion',
                    'associated hip spring'
                    ],
                [
                    act_time, act_time,
                    kin_time
                    ],
                [
                    hip_flexion_moment_left,
                    hip_flexion_moment_right,
                    spring_torque
                    ])
        except tables.NoSuchNodeError as e:
            report.next_is("Data not available.")
        h5.close()




















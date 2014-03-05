"""Manages the movement and use of data files."""

# TODO must edit cmc setup after it's moved, as well as grf.
# TODO put a comment in all files giving their original Hamner location? Not
# really necessary unless I share this code.
# TODO why is the model file so nested? there's definitely not that many
# changes that need to be made to the model file.
# TODO allow specifying which cycles to manage.

import csv
import difflib
import filecmp
import os
import shutil
import sys
import re
import xml.etree.ElementTree as etree

try: import tables
except ImportError, e: print e.message

try: import numpy as np
except ImportError, e: print e.message

if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    # Taken from /usr/lib/python2.7/posixpath.py
    # This method does not exist prior to python2.6.
    def relpath(path, start=os.path.curdir):
        """Return a relative version of a path"""
    
        if not path:
            raise ValueError("no path specified")
    
        start_list = [x for x in os.path.abspath(start).split(os.path.sep) if x]
        path_list = [x for x in os.path.abspath(path).split(os.path.sep) if x]
    
        # Work out how much of the filepath is shared by start and path.
        i = len(os.path.commonprefix([start_list, path_list]))
    
        rel_list = [os.path.pardir] * (len(start_list)-i) + path_list[i:]
        if not rel_list:
            return curdir
        return os.path.join(*rel_list)
    os.path.relpath = relpath

class GaitLandmarks(object):
    def __init__(self,
            primary_leg=None,
            cycle_start=None,
            cycle_end=None,
            left_strike=None,
            left_toeoff=None,
            right_strike=None,
            right_toeoff=None):
        self.primary_leg  = primary_leg
        self.cycle_start  = cycle_start
        self.cycle_end    = cycle_end
        self.left_strike  = left_strike
        self.left_toeoff  = left_toeoff
        self.right_strike = right_strike
        self.right_toeoff = right_toeoff

    def cycle_duration(self):
        return self.cycle_end - self.cycle_start


class ANCFile(object):
    """A plain-text file format for storing analog data from Motion Analysis
    Realtime. They have a file extension '.anc'. The file extension '.anb' is
    for binary files.

    The metadata for the file is stored in attributes of this object.

    This class is based off of similar code written by Amy Silder.

    """
    def __init__(self, fpath):
        """
        Parameters
        ----------
        fpath : str
            Valid file path to an ANC (.anc) file.

        """
        with open(fpath) as f:
            line1 = f.readline()
            line1list = line1.split('\t')
            self.file_type = line1list[1].strip()
            self.generation = line1list[3].strip()

            line2 = f.readline()
            line2list = line2.split('\t')
            self.board_type = line2list[1].strip()
            self.polarity = line2list[3].strip()

            line3 = f.readline()
            line3list = line3.split('\t')
            self.trial_name = line3list[1]
            self.trial_num = int(line3list[3])
            self.duration = float(line3list[5])
            self.num_channels = int(line3list[7])

            line4 = f.readline()
            line4list = line4.split('\t')
            self.bit_depth = int(line4list[1])
            self.precise_rate = float(line4list[3])

            line = f.readline()
            iline = 5
            while line.strip() == '':
                # There will most likely be a few empty lines.
                line = f.readline()
                iline += 1

            # Metadata for each column.
            header_row = line
            self.names = header_row.split()[1:]
            rate_row = f.readline()
            iline += 1
            self.rates = {self.names[i]: float(v) for i, v in
                    enumerate(rate_row.split()[1:])}
            range_row = f.readline()
            iline += 1
            self.ranges = {self.names[i]: float(v) for i, v in
                    enumerate(range_row.split()[1:])}

        dtype = {'names': ['time'] + self.names,
                'formats': (len(self.names) + 1) * ['float64']}
        self.data = np.loadtxt(fpath, delimiter='\t', skiprows=iline,
                    dtype=dtype)
        self.time = self.data['time']

    def __getitem__(self, name):
        """See `column()`.

        """
        return self.column(name)

    def __setitem__(self, name, val):
        """self.data[name] = val

        """
        self.data[name] = val

    def column(self, name):
        """
        Parameters
        ----------
        name : str
            Name of a column in the file (e.g., 'F1X'). For the 'time', column,
            just get the 'time' attribute.

        Returns
        -------
        col : np.array
            The data you are looking for.

        """
        return self.data[name]



class TRCFile(object):
    """A plain-text file format for storing motion capture marker trajectories.
    TRC stands for Track Row Column.

    The metadata for the file is stored in attributes of this object.

    See
    http://simtk-confluence.stanford.edu:8080/display/OpenSim/Marker+(.trc)+Files
    for more information.

    """
    def __init__(self, fpath):
        """
        Parameters
        ----------
        fpath : str
            Valid file path to a TRC (.trc) file.

        """
        # Read the header lines / metadata.
        # ---------------------------------
        # Split by any whitespace.
        # TODO may cause issues with paths that have spaces in them.
        f = open(fpath)
        # These are lists of each entry on the first few lines.
        first_line = f.readline().split()
        # Skip the 2nd line.
        f.readline()
        third_line = f.readline().split()
        fourth_line = f.readline().split()
        f.close()

        # First line.
        self.path = first_line[3]

        # Third line.
        self.data_rate = float(third_line[0])
        self.camera_rate = float(third_line[1])
        self.num_frames = int(third_line[2])
        self.num_markers = int(third_line[3])
        self.units = third_line[4]
        self.orig_data_rate = float(third_line[5])
        self.orig_data_start_frame = int(third_line[6])
        self.orig_num_frames = int(third_line[7])

        # Marker names.
        # The first and second column names are 'Frame#' and 'Time'.
        self.marker_names = fourth_line[2:]

        len_marker_names = len(self.marker_names)
        if len_marker_names != self.num_markers:
            warnings.warn('Header entry NumMarkers, %i, does not '
                    'match actual number of markers, %i. Changing '
                    'NumMarkers to match actual number.' % (
                        self.num_markers, len_marker_names))
            self.num_markers = len_marker_names

        # Load the actual data.
        # ---------------------
        col_names = ['frame_num', 'time']
        # This naming convention comes from OpenSim's Inverse Kinematics tool,
        # when it writes model marker locations.
        for mark in self.marker_names:
            col_names += [mark + '_tx', mark + '_ty', mark + '_tz']
        dtype = {'names': col_names,
                'formats': ['int'] + ['float64'] * (3 * self.num_markers + 1)}
        self.data = np.loadtxt(fpath, delimiter='\t', skiprows=6, dtype=dtype)
        self.time = self.data['time']

        # Check the number of rows.
        n_rows = self.time.shape[0]
        if n_rows != self.num_frames:
            warnings.warn('Header entry NumFrames, %i, does not '
                    'match actual number of frames, %i, Changing '
                    'NumFrames to match actual number.' % (
                        self.num_frames, n_rows))
            self.num_frames = n_rows

    def __getitem__(self, key):
        """See `marker()`.

        """
        return self.marker(key)

    def marker(self, name):
        """The trajectory of marker `name`, given as a `self.num_frames` x 3
        array. The order of the columns is x, y, z.

        """
        this_dat = np.empty((self.num_frames, 3))
        this_dat[:, 0] = self.data[name + '_tx']
        this_dat[:, 1] = self.data[name + '_ty']
        this_dat[:, 2] = self.data[name + '_tz']
        return this_dat

    def add_marker(self, name, x, y, z):
        """Add a marker, with name `name` to the TRCFile.

        Parameters
        ----------
        name : str
            Name of the marker; e.g., 'R.Hip'.
        x, y, z: array_like
            Coordinates of the marker trajectory. All 3 must have the same
            length.

        """
        if (len(x) != self.num_frames or len(y) != self.num_frames or len(z) !=
                self.num_frames):
            raise Exception('Length of data (%i, %i, %i) is not '
                    'NumFrames (%i).', len(x), len(y), len(z), self.num_frames)
        self.marker_names += [name]
        self.num_markers += 1
        self.data[name] = np.concatenate(x, y, z, axis=1)

    def marker_exists(self, name):
        """
        Returns
        -------
        exists : bool
            Is the marker in the TRCFile?

        """
        return name in self.marker_names

    def write(self, fpath):
        """Write this TRCFile object to a TRC file.

        Parameters
        ----------
        fpath : str
            Valid file path to which this TRCFile is saved.

        """
        f = open(fpath, 'w')

        # Line 1.
        f.write('PathFileType  4\t(X/Y/Z) %s\n' % os.path.split(fpath)[0])

        # Line 2.
        f.write('DataRate\tCameraRate\NumFrames\tNumMarkers\t'
                'Units\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n')

        # Line 3.
        f.write('%.1f\t%.1f\t%i\t%i\t%s\t%.1f\t%i\t%i\n' % (
            self.data_rate, self.camera_rate, self.num_frames,
            self.num_markers, self.units, self.orig_data_rate,
            self.orig_data_start_frame, self.orig_num_frames))

        # Line 4.
        f.write('Frame#\tTime\t')
        for imark in self.num_markers:
            f.write('%s\t\t\t' % self.marker_names[imark])
        f.write('\n')

        # Line 5.
        f.write('\t\t')
        for imark in self.num_markers:
            f.write('X%i\tY%s\tZ%s\t' % (imark, imark, imark))
        f.write('\n')

        # Line 6.
        f.write('\n')

        # Data.
        for iframe in self.num_frames:
            f.write('%i' % iframe)
            f.write('\t%.5f', self.time[iframe])
            for mark in self.marker_names:
                idxs = [mark + '_tx', mark + '_ty', mark + '_tz']
                f.write('\t%.3f\t%.3f\t%.3f' % (self.data[idxs]))
            f.write('\n')

        f.close()


def ndarray2storage(ndarray, storage_fpath, name=None, in_degrees=False):
    """Saves an ndarray, with named dtypes, to an OpenSim Storage file.

    Parameters
    ----------
    ndarray : numpy.ndarray
    storage_fpath : str
    in_degrees : bool, optional
    name : str
        Name of Storage object.

    """
    n_rows = ndarray.shape[0]
    n_cols = len(ndarray.dtype.names)

    f = open(storage_fpath, 'w')
    f.write('%s\n' % (name if name else storage_fpath,))
    f.write('version=1\n')
    f.write('nRows=%i\n' % n_rows)
    f.write('nColumns=%i\n' % n_cols)
    f.write('inDegrees=%s\n' % ('yes' if in_degrees else 'no',))
    f.write('endheader\n')
    for line_num, col in enumerate(ndarray.dtype.names):
        if line_num != 0:
            f.write('\t')
        f.write('%s' % col)
    f.write('\n')

    for i_row in range(n_rows):
        for line_num, col in enumerate(ndarray.dtype.names):
            if line_num != 0:
                f.write('\t')
            f.write('%f' % ndarray[col][i_row])
        f.write('\n')

    f.close()

def storage2numpy(storage_file, excess_header_entries=0):
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
    excess_header_entries : int, optional
        If the header row has more names in it than there are data columns.
        We'll ignore this many header row entries from the end of the header
        row. This argument allows for a hacky fix to an issue that arises from
        Static Optimization '.sto' outputs.

    Examples
    --------
    Columns from the storage file can be obtained as follows:
    
        >>> data = storage2numpy('<filename>')
        >>> data['ground_force_vy']

    """
    # What's the line number of the line containing 'endheader'?
    f = open(storage_file, 'r')

    header_line = False
    for i, line in enumerate(f):
        if header_line:
            column_names = line.split()
            break
        if line.count('endheader') != 0:
            line_number_of_line_containing_endheader = i + 1
            header_line = True
    f.close()

    # With this information, go get the data.
    if excess_header_entries == 0:
        names = True
        skip_header = line_number_of_line_containing_endheader
    else:
        names = column_names[:-excess_header_entries]
        skip_header = line_number_of_line_containing_endheader + 1
    data = np.genfromtxt(storage_file, names=names,
            skip_header=skip_header)

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

def cmc_input_fpaths(cmc_setup_fpath, replace=None):
    """Given a CMC setup file, returns the paths to all the files that the cmc
    setup file depends on.

    Parameters
    ----------
    cmc_setup_fpath : str
        Path to a CMC setup file.
    replace : dict, optional
        In case the paths in the files may be invalid, replace parts of the
        paths with the strings given in this dict. Keys are strings to look
        for, values are what to replace the key with.

    Returns
    -------
    inputs : dict
        A valid filepath to all the original CMC input files related to the
        provided CMC setup file.
        - setup
        - model
        - tasks
        - actuators (str if 1 file, list of str's if multiple files.
        - control_constraints
        - desired_kinematics
        - external_loads
        - force_plates
        - extload_kinematics

    """
    fname = cmc_setup_fpath

    setup = etree.parse(fname)

    inputs = dict()

    def valid_path_helper(file_containing_path, xml, tag):

        path = xml.findall('.//%s' % tag)[0].text
        return valid_path(path, file_containing_path)

    def valid_path(path, file_containing_path):
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
                os.path.join(os.path.split(file_containing_path)[0],
                    path.lstrip()))
        if os.path.exists(path2): return path2

        path2 = os.path.normpath(
                os.path.join(os.path.split(file_containing_path)[0],
                    path.rstrip()))
        if os.path.exists(path2): return path2

        raise Exception("Paths '%s' and '%s' do not exist." % (path, path2))

    # Get file names.
    # ---------------
    # Settings / parameters.
    inputs['model'] = valid_path_helper(fname, setup, 'model_file')
    inputs['tasks'] = valid_path_helper(fname, setup, 'task_set_file')

    # This is a list of files, not just 1 file.
    actu = list()
    for path in setup.findall('.//force_set_files')[0].text.split():
        actu.append(valid_path(path, fname))
    
    inputs['control_constraints'] = valid_path_helper(fname, setup,
            'constraints_file')

    # Data.
    inputs['desired_kinematics'] = valid_path_helper(fname, setup,
            'desired_kinematics_file')
    inputs['external_loads'] = valid_path_helper(fname, setup, 'external_loads_file')

    # Try to open the external loads file.
    extloads = etree.parse(inputs['external_loads'])
    inputs['force_plates'] = valid_path_helper(inputs['external_loads'],
            extloads, 'datafile')
    inputs['extload_kinematics'] = valid_path_helper(inputs['external_loads'],
            extloads, 'external_loads_model_kinematics_file')

    if len(actu) == 1:
        inputs['actuators'] = actu[0]
    else:
        inputs['actuators'] = actu

    return inputs


def copy_cmc_inputs(cmc_setup_fpath, destination, replace=None,
        do_not_copy=None, **kwargs):
    """Given a CMC setup file, copies all files necessary to run CMC over to
    `destination`. All files necessary to run CMC are stored in the same
    directory. The CMC setup file and the external loads files are edited so
    that they refer to the correct copied files.

    Also sets the results_directory to be 'results'.

    TODO allow placing files in another location (not a flat structure); but
    then update the path in the file so that it's relative to the file.

    TODO Assumes the path to each ForceSet file specified in force_set_files
    are separate by spaces, and that the paths themselves contain no spaces.
    Putting quotes around a path is not sufficient to allow a space in a path.

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
    do_not_copy : list of str's, optional
        Names of keys ('model', 'tasks', etc.; see the remaining parameters)
        for which files should not be copied over. The corresponding tags in
        the files will be updated so they refer to these original files, even
        if the setup files have moved. This takes precedence over the
        specification of new filenames for the remaining optional arguments.
        'setup' and 'external_loads' are necessarily copied over. 'actuators'
        are treated as a group: all the files in this field are copied, or none
        of them are copied.
    setup : str, optional
        A new filename for the cmc setup file (the first argument to this
        method).
    model : str, optional
        A new filename for this file.
    tasks : str, optional
        A new filename for this file.
    actuators : str, optional
        A new filename for this file. TODO NO LONGER WORKS.
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
    old_fpaths : dict
        A valid filepath to all the original CMC input files related to the
        provided CMC setup file.
        - setup
        - model
        - tasks
        - actuators (str if 1 file, list of str's if multiple files.
        - control_constraints
        - desired_kinematics
        - external_loads
        - force_plates
        - extload_kinematics
    new_fpaths : dict
        A valid filepath to the all the new files that were just copied over.
        The keys are as above for `old_fpaths`.
        NOTE/TODO: For backwards compatibility, new_fpaths contains an entry
        for 'actuators' if there is only one actuators (force set) file listed
        in the cmc setup file. If there is more than one 'actuators' file, then
        new_fpaths['actuators'] is a list of the force set files.

    """
    if do_not_copy != None:
        if 'setup' in do_not_copy or 'external_loads' in do_not_copy:
            raise Exception('`do_not_copy` cannot contain `setup` or '
                    '`external_loads`.')
    fname = cmc_setup_fpath

    setup = etree.parse(fname)

    old = dict()

    def valid_path_helper(file_containing_path, xml, tag):

        path = xml.findall('.//%s' % tag)[0].text
        return valid_path(path, file_containing_path)

    def valid_path(path, file_containing_path):
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
                os.path.join(os.path.split(file_containing_path)[0],
                    path.lstrip()))
        if os.path.exists(path2): return path2

        path2 = os.path.normpath(
                os.path.join(os.path.split(file_containing_path)[0],
                    path.rstrip()))
        if os.path.exists(path2): return path2

        raise Exception("Paths '%s' and '%s' do not exist." % (path, path2))

    # Get file names.
    # ---------------
    # Settings / parameters.
    old['model'] = valid_path_helper(fname, setup, 'model_file')
    old['tasks'] = valid_path_helper(fname, setup, 'task_set_file')

    # This is a list of files, not just 1 file.
    old_actu = list()
    for path in setup.findall('.//force_set_files')[0].text.split():
        old_actu.append(valid_path(path, fname))
    
    old['control_constraints'] = valid_path_helper(fname, setup, 'constraints_file')

    # Data.
    old['desired_kinematics'] = valid_path_helper(fname, setup,
            'desired_kinematics_file')
    old['external_loads'] = valid_path_helper(fname, setup, 'external_loads_file')

    # Try to open the external loads file.
    extloads = etree.parse(old['external_loads'])
    old['force_plates'] = valid_path_helper(old['external_loads'], extloads,
            'datafile')
    old['extload_kinematics'] = valid_path_helper(old['external_loads'],
            extloads, 'external_loads_model_kinematics_file')

    # Copy files over.
    # ----------------
    # We'll store the location of the copies.
    new_fpaths = dict()
    new_fpaths['setup'] = None
    new_fpaths['model'] = None
    new_fpaths['tasks'] = None
    new_actu_fpaths = list()
    new_fpaths['actuators'] = None
    new_fpaths['control_constraints'] = None
    new_fpaths['desired_kinematics'] = None
    new_fpaths['external_loads'] = None
    new_fpaths['force_plates'] = None
    new_fpaths['extload_kinematics'] = None

    if not os.path.exists(destination): os.makedirs(destination)
    for key, val in old.items():
        if val and key != 'external_loads' and (do_not_copy == None or key not
                in do_not_copy):
            if key in kwargs:
                new_fpath = os.path.join(destination, kwargs[key])
                shutil.copy(val, new_fpath)
                new_fpaths[key] = new_fpath
            else:
                shutil.copy(val, destination)
                new_fpaths[key] = os.path.join(destination,
                        os.path.basename(val))
    if do_not_copy == None or 'actuators' not in do_not_copy:
        for actu in old_actu:
            shutil.copy(actu, destination)
            new_actu_fpaths.append(
                    os.path.join(destination, os.path.basename(actu)))

    # Edit the names of the files in the setup files.
    # -----------------------------------------------
    def edit_field(xml, tag, key):
        if old[key]:
            if do_not_copy != None and key in do_not_copy:
                xml.findall('.//%s' % tag)[0].text = \
                        os.path.relpath(old[key], destination)
            else:
                if key in kwargs:
                    newvalue = kwargs[key]
                else:
                    newvalue = os.path.basename(old[key])
                xml.findall('.//%s' % tag)[0].text = newvalue

    setup.findall('.//results_directory')[0].text = 'results'
    edit_field(setup, 'model_file', 'model')
    edit_field(setup, 'task_set_file', 'tasks')
    #edit_field(setup, 'force_set_files', 'actuators')
    edit_field(setup, 'constraints_file', 'control_constraints')
    edit_field(setup, 'desired_kinematics_file', 'desired_kinematics')
    edit_field(setup, 'external_loads_file', 'external_loads')

    # We cannot use edit_field() to edit the force_set_files field, because
    # it's more complicated. So we have these next bunch of lines to do that.
    new_actu_value = ''
    if do_not_copy != None and 'actuators' in do_not_copy:
        for this_path in old_actu:
            new_actu_value += ' %s' % os.path.relpath(this_path, destination)
    else:
        for this_path in new_actu_fpaths:
            new_actu_value += ' %s' % os.path.basename(this_path)
    setup.findall('.//%s' % 'force_set_files')[0].text = new_actu_value

    # Give the user the new actuator file paths.
    if len(new_actu_fpaths) > 0:
        if len(new_actu_fpaths) == 1:
            # Return just that one file path.
            new_fpaths['actuators'] = new_actu_fpaths[0]
        else:
            # Return a list of file paths.
            new_fpaths['actuators'] = new_actu_fpaths

    edit_field(extloads, 'datafile', 'force_plates')
    edit_field(extloads, 'external_loads_model_kinematics_file',
            'extload_kinematics')

    if 'setup' in kwargs:
        setup_new_fpath = os.path.join(destination, kwargs['setup'])
    else:
        setup_new_fpath = os.path.join(destination,
                os.path.basename(cmc_setup_fpath))
    setup.write(setup_new_fpath)
    new_fpaths['setup'] = setup_new_fpath

    if 'external_loads' in kwargs:
        extloads_new_fpath = os.path.join(destination, kwargs['external_loads'])
    else:
        extloads_new_fpath = os.path.join(destination,
                os.path.basename(old['external_loads']))
    extloads.write(extloads_new_fpath)
    new_fpaths['external_loads'] = extloads_new_fpath

    # Also, now we can also give the user the list of old actuator paths.
    if len(old_actu) == 1:
        old['actuators'] = old_actu[0]
    else:
        old['actuators'] = old_actu

    return old, new_fpaths


def dock_simulation_tree_in_pytable(h5fname, study_root, h5_root, verbose=True,
        **kwargs):
    """Docks all simulations in the tree into the h5file at the desired
    location. The directory structure used in the h5 file is the same that is
    used for the simulation output. All leaves in the tree MUST be simulation
    outputs (e.g. contain STO files).

    Parameters
    ----------
    h5fname : str
        Name of the pyTables/HDF5 file in which to dock the output.
    study_root : str
        Path to the directory tree containing simulation output to be docked.
    h5_root : list of str's
        The group in the pyTables/HDF5 file in which these simulations will be
        docked.
    verbose : bool, optional (default: True)
        Prints a notice for every output loaded into the pyTables file.
    kwargs : passed onto `dock_output_in_pytable`.

    """
    # TODO overwrite : bool, optional (default : False)
    # TODO     If a group already exists, delete it and rewrite it with the
    # TODO     newly-found data.  Otherwise, the group is skipped.
    # Open the pyTables file.
    h5file = tables.open_file(h5fname, mode='a')

    # Report the number of exceptions we get.
    exception_count = 0

    # Walk the entire directory structure.
    for (path, dirs, files) in os.walk(study_root):

        # Only want leaves: there are no more directories in the directory.
        # This path also needs to have files to be considered.
        # TODO an alternative check is os.path.split(path)[1] == 'output'
        if len(dirs) == 0 and len(files) != 0:

            if verbose:
                print "Loading {0}.".format(path)

            # Create a list describing the path, excluding 'output' at the end.
            # Also make sure we do this for a path that's relative to the study
            # root.
            path_parts = _splitall(os.path.relpath(os.path.split(path)[0],
                study_root))

            # Prepend the user's desired root for this study in the h5 file.
            group_path = h5_root + path_parts

            # Dock this specific simulation output.
            try:
                this_table = dock_output_in_pytable(h5file, path,
                        group_path, **kwargs)
            except Exception, e:
                print "Exception at path {0}: {1}".format(path, e.message)
                exception_count += 1

    print "Number of exceptions: %i" % exception_count

    # Close the pyTables file.
    h5file.close()

def dock_trc_in_pytable(h5file, trc_fpath, table_name, group_path, title='',
        overwrite_if_newer=False):
    """Write the contents of a TRC file to the database.

    Parameters
    ----------
    h5file : tables.File
        pyTables File object, opened using tables.open_file(...). Does NOT
        close the file.
    trc_fpath : str
        Valid path to a TRC file. We load this into memory using
        `dataman.TRCFile`.
    table_name : str
        Name of the table in the database that'll hold this data.
    group_path : str, or list of str's
        The group tree hierarchy specifying where the output is to be docked in
        the h5file; as a path or as list of each directory's name (e.g.:
        'path/to/file' or ['path', 'to', 'file'])
    title : str, optional
        Title, in the pyTables file, for this group.
    overwrite_if_newer : bool, optional (default: False)
        By default, the tables cannot already exist in the group specified.
        However, if this is set to True, and the output_path's mtime is greater
        than the mtime stored in this group, we will first delete the tables so
        that the newer tables can be written to the database. This is done on a
        per-table basis. Skip popluation of a table if the table already exists
        AND the data isn't any newer.

    Returns
    -------
    current_group : tables.Group
        The pyTables group in which the output has been stored.

    """
    # If trc_fpath doesn't exist, can't do anything.
    if not os.path.exists(trc_fpath):
        raise Exception("TRC file {0:r} doesn't exist.".format(trc_fpath))

    # Convert group_path to list of str's, if necessary.
    if type(group_path) == str:
        group_path = _splitall(group_path)

    # -- Make all necessary groups to get to where we're going.
    current_group = _blaze_group_trail(h5file, group_path, title)

    # If we are considering overwriting and we SHOULD overwrite (is_newer),
    # then remove the existing group.
    if overwrite_if_newer:
        if hasattr(current_group, table_name):
            if not hasattr(getattr(current_group, table_name).attrs, 'mtime'):
                raise Exception("Table exists, but is not labeled with an "
                        "mtime. Not overwriting.")
            if (getattr(current_group, table_name).attrs.mtime <
                    os.path.getmtime(trc_fpath)):
                getattr(current_group, table_name)._f_remove(True) 
                _populate_table_with_trc(h5file, current_group, table_name,
                        trc_fpath)
            else:
                # Table exists and isn't newer; skip writing.
                pass
        else:
            # Table doesn't exist; write it!
            _populate_table_with_trc(h5file, current_group, table_name,
                    trc_fpath)
    else:
        # Try to populate indiscriminantly.
        _populate_table_with_trc(h5file, current_group, table_name, trc_fpath)

    # Update the attribute for when this group was last updated.
    # This must go after _populate_table, because otherwise the table is
    # not necessarily created yet.
    getattr(current_group, table_name).attrs.mtime = \
            os.path.getmtime(trc_fpath)

    return current_group

def _populate_table_with_trc(h5file, group, table_name, trc_fpath):
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
    trc_fpath : str
        Valid path to a TRC file. We load this into memory using
        `dataman.TRCFile`.

    Returns
    -------
    table : tables.Table
        The table that has just been created.

    """
    # Open data file.
    trcf = TRCFile(trc_fpath)

    # Create the table.
    # -----------------
    # Save this row for later; we'll need it.
    orig_col_names = trcf.data.dtype.names
    col_names = list(trcf.data.dtype.names)

    # Can't have periods in table column names in pyTables.
    for i in range(len(col_names)):
        col_names[i] = col_names[i].replace('.', '_')

    # Create pyTables table columns.
    table_cols = dict()
    for col in col_names:
        # Checking if the column is empty. This is a
        # once-in-a-blue-moon bug fix as a result of inconsistency in
        # Hamner's files. See CMC results for subject 2, speed 2 m/s,
        # cycle 1, states_OG.sto file.
        if col != '':
            table_cols[col] = tables.Float32Col()

    # Create pyTables table.
    table = h5file.create_table(group, table_name, table_cols,
            'Output file {0}'.format(trc_fpath))

    # Add data to the table.
    # ----------------------
    # TODO could probably do this more efficiently than adding rows one by one.
    # For each column in the data file in this row.
    for it in range(trcf.num_frames):
        for i in range(len(table_cols.keys())):

            # Append the data into the table.
            table.row[col_names[i]] = trcf.data[orig_col_names[i]][it]

        # Tell pyTables to append this data to the table.
        table.row.append()

    # Save (?).
    table.flush()

    # Give access to it.
    return table

def dock_output_in_pytable(h5file, output_path, group_path, allow_one=False,
        title='', ext='.sto', overwrite_if_newer=False,
        remove_shared_name=True, table_name_repl={}, **kwargs):
    """Docks an OpenSim output, via a table for each STO (see `ext`) file, in a
    pyTable file.

    Parameters
    ----------
    h5file : tables.File
        pyTables File object, opened using tables.open_file(...). Does NOT
        close the file.
    output_path : str
        File path in which the OpenSim output is located (e.g. .STO files).
        Only .STO files are loaded, and it is assumed that all .STO (see `ext`)
        files in this directory are from one run. That is, they have the same
        prefix, which is the name of the run.
    group_path : str, or list of str's
        The group tree hierarchy specifying where the output is to be docked in
        the h5file; as a path or as list of each directory's name (e.g.:
        'path/to/file' or ['path', 'to', 'file'])
    allow_one : bool, optional (default: False)
        Allows the loading of just one Storage file. Otherwise, an exception is
        thrown. It is common that if only one Storage file exists, it is a
        partial states file and means that the simulation did not complete.
    title : str, optional
        Title, in the pyTables file, for this group.
    ext : str, optional
        The filename ending of the Storage files. Default is '.sto', but the
        user may want to use '.mot' also.
    overwrite_if_newer : bool, optional (default: False)
        By default, the tables cannot already exist in the group specified.
        However, if this is set to True, and the output_path's mtime is greater
        than the mtime stored in this group, we will first delete the tables so
        that the newer tables can be written to the database. This is done on a
        per-table basis. Skip popluation of a table if the table already exists
        AND the data isn't any newer.
    remove_shared_name : bool, optional
        When creating names of tables, remove the shared portion of the Storage
        file names.
    table_name_repl : dict, optional
        Key-value pairs of regular expressions and corresponding replacements
        for table names.
    **kwargs : optional
        Passed onto _populate_table. May want to use the kwarg 'replacements'.

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
    storage_files = [f for f in os.listdir(output_path) if f.endswith(ext)]

    # If there are no storage files, the user probably gave a bad path.
    if len(storage_files) == 0:
        raise Exception("No {0} files found in {1}.".format(ext, output_path))

    # If there's only one, usually the states file, forget about this output.
    if (not allow_one) and len(storage_files) == 1:
        raise Exception("Only one {0} file found: {1}.".format(ext,
            storage_files[0]))

    # Get the length of the common prefix of these files.
    if remove_shared_name:
        n_shared = _length_of_shared_prefix(storage_files)

    # -- Add tables in the current group.

    # Loop through all storage files.
    for f in storage_files:

        # Path to the data file.
        filepath = os.path.join(output_path, f)

        # Get name of the table: after the run name and before the file ext.
        if len(storage_files) == 1 or not remove_shared_name:
            table_name = os.path.splitext(f)[0]
        else:
            table_name = os.path.splitext(f)[0][n_shared:]
        for find, rep in table_name_repl.items():
            regex = re.compile(find)
            table_name = regex.sub(rep, table_name)

        # If we are considering overwriting and we SHOULD overwrite (is_newer),
        # then remove the existing group.
        if overwrite_if_newer:
            if hasattr(current_group, table_name):
                # Table exists.
                table = getattr(current_group, table_name)
                if ((not hasattr(table.attrs, 'mtime')) or 
                        (table.attrs.mtime < os.path.getmtime(filepath))):
                    # Table exists and ((there is newer data available) OR (we
                    # don't know how old the stored table is)).
                    getattr(current_group, table_name)._f_remove(True) 
                    _populate_table(h5file, current_group, table_name,
                            filepath, **kwargs)
                else:
                    # Table exists and isn't newer; skip writing.
                    pass
            else:
                # Table doesn't exist; write it!
                _populate_table(h5file, current_group, table_name, filepath,
                        **kwargs)
        else:
            # Try to populate indiscriminantly.
            _populate_table(h5file, current_group, table_name, filepath,
                    **kwargs)

        # Update the attribute for when this group was last updated.
        # This must go after _populate_table, because otherwise the table is
        # not necessarily created yet.
        if hasattr(current_group, table_name):
            getattr(current_group, table_name).attrs.mtime = \
                    os.path.getmtime(filepath)

    return current_group

def _blaze_group_trail(h5file, group_path, title=''):

    # Start at the root.
    current_group = h5file.root

    # Loop through each group in the path the user provided.
    for next_group_name in group_path:

        # This group has not been created yet.
        if not hasattr(current_group, next_group_name):

            # Create the group.
            h5file.create_group(current_group, next_group_name, title=title)

        # Set this group as the current group.
        current_group = getattr(current_group, next_group_name)

    # Return this leaf group.
    return current_group


def _populate_table(h5file, group, table_name, filepath, replacements={}):
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
    replacements : dict, optional
        In Storage file column names, replace each instance of key with the
        corresponding value. This is useful, for example, in converting the
        Storage column names 1_ground_force_x to l_ground_force_x. The keys can
        be regular expressions.

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
                for find, rep in replacements.items():
                    regex = re.compile(find)
                    title_row[i] = regex.sub(rep, title_row[i])

            # TODO Hack to deal with bug in static optimization.
            if 'Right_GRF' in title_row: title_row.remove('Right_GRF')
            if 'Left_GRF' in title_row: title_row.remove('Left_GRF')

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
            table = h5file.create_table(group, table_name, table_cols,
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

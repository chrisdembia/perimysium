"""Manages the movement and use of data files."""

import shutil
import os
import abc

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
    expdir = 'HPL_Data/Hamner_Test/'
    osimdir = 'OpenSimResults/'
    cmc_actuators_from = 'gait2392_CMC_Actuators.xml'
    cmc_tasks_from = 'gait2392_CMC_Tasks.xml'
    cmc_actuators_to = 'cmc_actuators.xml'
    cmc_tasks_to = 'cmc_tasks.xml'

    n_total_subjects = 20

    class CopyMode:
        MINIMAL=1
        CONSERVATIVE=1

    def __init__(self, n_subjects):
        """
        Parameters
        ----------
        n_subjects : int
            The number of subjects to copy over. Max is
            HamnerXX.n_total_subjects.
        """
        if n_subjects > self.n_total_subjects:
            raise Exception("Max subjects is {0}.".format(
                self.n_total_subjects))
        self.n_subjects = n_subjects
        self._subjects = []
        for i_subj in range(n_subjects):
            self._subjects += [Subject(self, i_subj+1)]

    @property
    def copy_mode(self): return self._copy_mode

    @copy_mode.setter
    def copy_mode(self, mode): self._copy_mode = mode

    def copy_new(self, from_path, to_path):
        """Copies and renames files necessary for CMC. Files are stored in a
        hierarchy that reveals the extent to which files are shared. Files are
        NOT modified.
        
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

        Parameters
        ----------
        from_path : str
            Something that probably ends in 'Hamner_Data'.
        to_path : str
            Probably some local directory on a machine (probably not in a
            repository or Dropbox) where the simulations are to be run.

        """

        # Append separation if the user did not provide it.
        if from_path[len(from_path)-1] != os.sep: from_path += os.sep
        if to_path[len(to_path)-1] != os.sep: to_path += os.sep

        # If to_path doesn't exist, make it!
        if not os.path.isdir(to_path): os.mkdir(to_path)

        for subj in self._subjects:
            subj.copy_new(self, from_path, to_path)


class Subject(object):

    name_pre = 'subject'
    cmc_dirname = 'cmc_multipleSteps_v24_Tendon_040_Vmax_15_Passive_10_2X/'

    def __init__(self, hamner, number):
        """ TODO """
        self.hamner = hamner
        self.name = self.name_pre + str(number)
        self._speeds = []

    def copy_new(self, from_path, to_path):


        if hamner.copy_mode == Hamner.CopyMode.MINIMAL:

            # TODO do diff's to see if files actually differ, precluding the
            # minimal copying.
            
            # cmc_actuators
            shutil.copy2(os.path.join(from_path, hamner.osimdir, self.name,
                                      self.cmc_dirname,
                                      hamner.cmc_actuators_from),
                         os.path.join(to_path, hamner.cmc_actuators_to))

        if hamner.copy_mode == Hamner.CopyMode.CONSERVATIVE:
            raise Exception("TODO")


class Cycle(object):

    def __init__(self):
        pass

    def copy_new(self, from_path, to_path):
        raise Exception("TODO")




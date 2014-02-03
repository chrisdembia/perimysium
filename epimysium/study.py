"""
s = Study('loadedwalking')
s01 = s.subject_new(1)
s01nl = s01.condition_new('noload')
s01nls = s01nl.condition_new('static')
s01nlf = s01nl.condition_new('free')
s01nlf.trial_new(1)
s01nlf.static_trial_new(1)
s01nlf.motion_capture_fname_is(', 'motion_capture.trc', 'ground_reaction.mot')

s01nlf01.base_sim
s01nlf01.predict_sim.

shelve
jsonpickle
json
"""

class Object(object):
    def __init__(self, name):
        self._name = name

    @property
    def name(self):
        return self._name

class Study(Object):
    """
    Attributes
    ----------
    subject : dict
        of Subject's.

    """
    def __init__(self, name):
        super(Study, self).__init__(name)
        self.subjects

    def subject(self, key):
        return self._subjects[key]

    def subject_new(self, key):


class Subject(Object):
    """

    Attributes
    ----------
    number : int

    """
    def __init__(self, number):
        super(Subject, self).__init__('subject%02i' % number)
        self.number = number

class Condition(Object):
    """A node in the Study tree.

    Attributes
    ----------
    conditions : dict
        of subconditions.
    trial : list
        of Trial's.

    """
    def __init__(self, name):
        super(Trial, self).__init__(name)
        pass

class Trial(Object):
    """

    Attributes
    ----------
    number : int

    """
    def __init__(self, number, condition):
        super(Trial, self).__init__('trial%02i' % number)
        self.number = number
        self._condition = condition

class TrackingSimulation(Object):
    """

    """
    def __init__(self, name):
        super(TrackingSimulation, self).__init__(name)

    def average_metabolic_expenditure_rate(self):
        raise NotImplementedError()

class CMC(Object):
    """
    """
    def __init__(self, trial):


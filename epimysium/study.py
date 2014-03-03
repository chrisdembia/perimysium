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

import yaml

class Object(yaml.YAMLObject):
    yaml_tag = u'Object'

    def __init__(self, name):
        self.name = name

#    @classmethod
#    def to_yaml(cls, dumper, data):
#        return dumper.represent_mapping(data.yaml_tag, {'name': data.name})
#
#    @classmethod
#    def from_yaml(cls, loader, node):
#        print 'DEBUGGG'
#        print dir(node)
#        print node.value
#        print dir(node.value)
#        print node.value[0]
#        print dir(node.value[0])
#        print dir(node.value[0][0])
#        print node.value[0][0]
#        print node.value[0][0].value
#        value = loader.construct_scalar(node)
#        return Object(node.value)

#    def __repr__(self):
#        return "%s(name=%r)" % (self.__class__.__name__, self.name)

    def save(self, fpath):
        with open(fpath, 'w') as f:
            f.write(yaml.dump(self, default_flow_style=False, indent=4))

    @classmethod
    def load(cls, fpath):
        with open(fpath) as f:
            return yaml.load(f.read())

class Study(Object):
    """
    Attributes
    ----------
    subjects : dict
        of Subject's.

    """

    yaml_tag = u'!Study'

    def __init__(self, name):
        super(Study, self).__init__(name)
        self.subjects = []

class Subject(Object):
    """

    Attributes
    ----------
    number : int
    conditions : dict
        of subconditions.

    """

    yaml_tag = u'!Subject'

    def __init__(self, number):
        super(Subject, self).__init__('subject%02i' % number)
        self.number = number
        self.conditions = dict()

    def condition_new(self, cond):
        self.conditions[cond.name] = cond

    def condition(self, name):
        return self.conditions[name]

class Condition(Object):
    """A node in the Study tree.

    Attributes
    ----------
    conditions : dict
        of subconditions.
    trial : list
        of Trial's.

    """

    yaml_tag = u'!Condition'

    def __init__(self, name):
        super(Condition, self).__init__(name)
        self.conditions = dict()
        self.trials = []

    def condition_new(self, cond):
        self.conditions[cond.name] = cond

    def condition(self, name):
        return self.conditions[name]

class Trial(Object):
    """

    Attributes
    ----------
    number : int

    """

    yaml_tag = u'!Trial'

    def __init__(self, number):
        super(Trial, self).__init__('trial%02i' % number)
        self.number = number

# class TrackingSimulation(Object):
#     """
# 
#     """
#     def __init__(self, name):
#         super(TrackingSimulation, self).__init__(name)
# 
#     def average_metabolic_expenditure_rate(self):
#         raise NotImplementedError()
# 
# class CMC(Object):
#     """
#     """
#     def __init__(self, trial):
# 

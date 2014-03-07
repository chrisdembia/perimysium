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

    yaml_tag = u'!Object'

    def __init__(self, name):
        self.name = name

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'name': data.name,
            })
#
#    @classmethod
#    def from_yaml(cls, loader, node):
#        value = loader.construct_mapping(node)
#        return Object(node.value)

#    def __repr__(self):
#        return '%s(name=%r)' % (self.__class__.__name__, self.name)

    def save(self, fpath):
        # TODO move to Study.
        with open(fpath, 'w') as f:
            Dumper = yaml.Dumper
            Dumper.ignore_aliases = lambda self, data: True
            yaml.dump(self, stream=f, default_flow_style=False, indent=4,
                    Dumper=Dumper)

    def __str__(self):
        return self.name

class Study(Object):
    """
    Attributes
    ----------
    subjects : dict
        of Subject's.

    """

    yaml_tag = u'!Study'

    def __init__(self, name, pytable_fpath):
        super(Study, self).__init__(name)
        self.pytable_fpath = pytable_fpath
        self.subjects = dict()

    def subject_new(self, subject):
        self.subjects[subject.number] = subject
        # TODO subject.study = self

    def subject(self, number):
        return self.subjects[number]

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'name': data.name,
            'pytable_fpath': data.pytable_fpath,
            'subjects': data.subjects,
            })

    @classmethod
    def from_yaml(cls, loader, node):
        study_dict = loader.construct_mapping(node)
        study = Study(study_dict['name'], study_dict['pytable_fpath'])
        for study_node in node.value:
            key = study_node[0].value
            if key == 'subjects':
                subject_dict_mapping = study_node[1]
                subject_dict = loader.construct_mapping(subject_dict_mapping)
                for k, v in subject_dict.items():
                    study.subject_new(v)
                    v.study = study

        return study

    @classmethod
    def load(cls, fpath):
        with open(fpath) as f:
            return yaml.load(f.read())

#    def __repr__(self):
#        print self.subjects
#        return '%s(name=%r, pytable_fpath=%r, subjects=%r)' % (
#                self.__class__.__name__,
#                self.name,
#                self.pytable_fpath,
#                self.subjects)

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
        self.study = None

    def condition_new(self, cond):
        self.conditions[cond.name] = cond

    def condition(self, name):
        return self.conditions[name]

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'number': data.number,
            'conditions': data.conditions,
            })

    @classmethod
    def from_yaml(cls, loader, node):
        subject_dict = loader.construct_mapping(node)
        subject = Subject(subject_dict['number'])
        for subject_node in node.value:
            key = subject_node[0].value
            if key == 'conditions':
                condition_dict_mapping = subject_node[1]
                condition_dict = loader.construct_mapping(
                        condition_dict_mapping)
                for k, v in condition_dict.items():
                    subject.condition_new(v)
                    v.subject = subject
        return subject

#    def __repr__(self):
#        return '%s(number=%i, conditions=%r)' % (
#                self.__class__.__name__,
#                self.number,
#                self.conditions,
#                )

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
        self.trials = dict()
        self.subject = None
        self.parent = None

    def condition_new(self, cond):
        self.conditions[cond.name] = cond

    def condition(self, name):
        return self.conditions[name]

    def trial_new(self, trial):
        self.trials[trial.number] = trial

    def trial(self, number):
        return self.trials[number]

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'name': data.name,
            'conditions': data.conditions,
            'trials': data.trials,
            })

    @classmethod
    def from_yaml(cls, loader, node):
        condition_dict = loader.construct_mapping(node)
        condition = Condition(condition_dict['name'])
        for condition_node in node.value:
            key = condition_node[0].value
            if key == 'conditions':
                subcond_dict_mapping = condition_node[1]
                subcond_dict = loader.construct_mapping(subcond_dict_mapping)
                for k, v in subcond_dict.items():
                    condition.condition_new(v)
                    v.parent = condition
                    v.subject = condition.subject
            if key == 'trials':
                trial_dict_mapping = condition_node[1]
                trial_dict = loader.construct_mapping(trial_dict_mapping)
                for k, v in trial_dict.items():
                    condition.trial_new(v)
                    v.condition = condition
        return condition

#    def __repr__(self):
#        if len(self.conditions) > 1:
#            cond_repr = {k: repr(v) for k, v in self.conditions.items()}
#        else:
#            cond_repr = 'dict()'
#        return '%s(name=%r, conditions=%s, trials=%r)' % (
#                self.__class__.__name__,
#                self.name,
#                cond_repr,
#                self.trials,
#                )

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
        self.condition = None

    def specific_metabolic_cost(self):
        pass

    def average_over_gait_cycle(self, h5_table_path, column_name):
        pass

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'number': data.number,
            })

    @classmethod
    def from_yaml(cls, loader, node):
        trial_dict = loader.construct_mapping(node)
        return Trial(trial_dict['number'])

#    def __repr__(self):
#        return '%s(number=%i)' % (self.__class__.__name__, self.number)

# class TrackingSimulation(Object):
#     """
#     # TODO knows where input files are.
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

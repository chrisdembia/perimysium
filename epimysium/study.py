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

#    def __repr__(self):
#        return "%s(name=%r)" % (self.__class__.__name__, self.name)

    def save(self, fpath):
        with open(fpath, 'w') as f:
            # YAML will automatically create anchors and aliases for classes it
            # thinks are the same and used in multiple places. We don't want
            # that.
            # http://pyyaml.org/ticket/91
            Dumper = yaml.SafeDumper
            Dumper.ignore_aliases = lambda self, data: True
            yaml.dump(self, stream=f, default_flow_style=False, indent=4,
                    Dumper=Dumper)

    @classmethod
    def load(cls, fpath):
        with open(fpath) as f:
            return yaml.load(f.read())

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'name': data.name,
            })

    @classmethod
    def from_yaml(cls, loader, node):
        value = loader.construct_mapping(node)
        return Object(node.value)

class Study(Object):
    """
    Attributes
    ----------
    subjects : dict
        of Subject's.

    """

    yaml_tag = u'!Study'

    def __init__(self, name, subjects=dict()):
        super(Study, self).__init__(name)
        self.subjects = subjects

    def subject_new(self, subject):
        self.subjects[subject.number] = subject

    def subject(self, number):
        return self.subjects[number]

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'name': data.name,
            'subjects': data.subjects,
            })

    @classmethod
    def from_yaml(cls, loader, node):
        study_dict = loader.construct_mgpping(node)
        print study_dict
        study = Study(study_dict['name'])
        for study_node in node.value:
            key = study_node[0].value
            if key == 'subjects':
                subject_sequence = study_node[1]
                print 'HELLO'
                print type(study_node[1])
                print type(subject_sequence)
                #for subject_mapping in subject_sequence.value:
                #    print ''
                #    print subject_mapping
                #    #subject = loader.construct_scalar(subject_mapping)
                #    #print subject
                #    #study.subject_new(subject)

        #            #print 'hi', v
                subject_objects = loader.construct_sequence(subject_sequence)
                print subject_objects
                print subject_objects[0]
                #for isubj, subj in enumerate(subject_objects):
                #    study.subject_new(subj)

        #            this_subject_mapping = subject_sequence.value[isubj]
        #            for subject_node in this_subject_mapping.value:
        #                key = subject_node[0].value
        #                if key == 'conditions':
        #                    cond_sequence = subject_node[1]
        #                    cond_objects = loader.construct_sequence(cond_sequence)
        #                    for icond, cond in enumerate(cond_objects):
        #                        subj.condition_new(cond)

        #        #import pdb; pdb.set_trace()
        #        #print value

        return study

class Subject(Object):
    """

    Attributes
    ----------
    number : int
    conditions : dict
        of subconditions.

    """

    yaml_tag = u'!Subject'

    def __init__(self, number, conditions=dict()):
        super(Subject, self).__init__('subject%02i' % number)
        self.number = number
        self.conditions = conditions

    def condition_new(self, cond):
        self.conditions[cond.name] = cond

    def condition(self, name):
        return self.conditions[name]

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'number': data.number,
            'conditions': data.conditions.values(),
            })

#    @classmethod
#    def from_yaml(cls, loader, node):
#        value = loader.construct_mapping(node)
#        return Subject(value['number'])

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

    def __init__(self, name, conditions=dict(), trials=dict()):
        super(Condition, self).__init__(name)
        self.conditions = conditions
        self.trials = trials

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
            'conditions': data.conditions.values(),
            'trials': data.trials.values(),
            })

#    @classmethod
#    def from_yaml(cls, loader, node):
#        value = loader.construct_mapping(node)
#        return Condition(value['name'])


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

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_mapping(data.yaml_tag, {
            'number': data.number,
            })

#    @classmethod
#    def from_yaml(cls, loader, node):
#        value = loader.construct_mapping(node)
#        return Trial(value['number'])

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

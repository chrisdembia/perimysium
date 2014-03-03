import os

from nose.tools import assert_equal
import yaml

from epimysium.study import Study, Subject, Condition, Trial

parentdir = os.path.abspath(os.path.dirname(__file__))

def fpath(fname):
    return os.path.join(parentdir, fname)

def test_serialization_dump():
    s = Study('walking')
    s.subjects.append(Subject(1))
    s.subjects.append(Subject(2))
    s.subjects[0].condition_new(Condition('noload'))
    s.subjects[0].condition_new(Condition('loaded'))
    s.subjects[0].condition('noload').condition_new(Condition('free'))
    s.subjects[0].condition('noload').condition('free').trials.append(Trial(1))
    
    outfpath = fpath('walking_study.yml')
    outfpath_des = fpath('walking_study_des.yml')
    s.save(outfpath)

    assert_equal(open(outfpath).read(), open(outfpath_des).read())

    os.remove(outfpath)

def test_serialization_load():
    s = Study.load(fpath('walking_study_des.yml'))
    assert_equal(s.name, 'walking')
    assert_equal(s.subjects[0].name, 'subject01')

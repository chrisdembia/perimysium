import os

from nose.tools import assert_equal
import yaml

from epimysium.study import Study, Subject, Condition, Trial

parentdir = os.path.abspath(os.path.dirname(__file__))

def fpath(fname):
    return os.path.join(parentdir, fname)

def test_serialization_dump():
    s = Study('walking')
    s.subject_new(Subject(1))
    s.subject_new(Subject(2))
    s.subject(1).condition_new(Condition('noload'))
    s.subject(1).condition_new(Condition('loaded'))
    s.subject(1).condition('noload').condition_new(Condition('free'))
    s.subject(1).condition('noload').condition('free').trial_new(Trial(1))
    
    outfpath = fpath('walking_study_dump.yml')
    outfpath_des = fpath('walking_study_des.yml')
    s.save(outfpath)

    assert_equal(open(outfpath).read(), open(outfpath_des).read())

    os.remove(outfpath)

def test_serialization_load():
    outfpath = fpath('walking_study_load.yml')
    outfpath_des = fpath('walking_study_des.yml')

    s = Study.load(outfpath_des)
    s.save(outfpath)

    assert_equal(open(outfpath).read(), open(outfpath_des).read())

    os.remove(outfpath)

    assert_equal(s.name, 'walking')
    assert_equal(s.subject(1).name, 'subject01')

if __name__ == '__main__':
    test_serialization_load()

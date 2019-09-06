from nose.tools import ok_, eq_

from hommod.models.range import SequenceRange


def test_left_right_from():

    r1 = SequenceRange(1, 11, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(9, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    ok_(r1.is_left_from(r2))
    ok_(r2.is_right_from(r1))


def test_distance_from():

    r1 = SequenceRange(1, 9, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(11, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    eq_(r1.get_distance_from(r2), 2)


def test_intersection():

    r1 = SequenceRange(1, 15, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(10, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    i = r1.get_intersection(r2)

    eq_(i.start, 10)
    eq_(i.end, 15)


def test_subtraction():

    r = SequenceRange(5, 15, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    r -= 5

    eq_(r.start, 0)
    eq_(r.end, 10)


def test_eq():

    r1 = SequenceRange(0, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(0, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    ok_(r1 == r2)

def test_ne():

    r1 = SequenceRange(10, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(0, 10, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    ok_(r1 != r2)


def test_length():

    r = SequenceRange(5, 15, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    eq_(r.get_length(), 10)


def test_sub_sequence():

    r = SequenceRange(5, 15, "AAAAAAAVAVAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    eq_(r.get_sub_sequence(), "AAVAVAAAAA")


def test_includes_residue():

    r = SequenceRange(5, 15, "AAAAAAAVAVAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    ok_(r.includes_residue(6))
    ok_(r.includes_residue(10))
    ok_(r.includes_residue(15))


def test_percentage_overlap():

    r1 = SequenceRange(0, 10, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(5, 15, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    eq_(r1.get_percentage_overlap(r2), 50.0)


def test_overlaps_with():

    r1 = SequenceRange(0, 10, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(5, 15, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    ok_(r1.overlaps_with(r2))


def test_not_overlaps_with():
    r1 = SequenceRange(0, 10, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(10, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    ok_(not r1.overlaps_with(r2))


def test_encloses():

    r1 = SequenceRange(0, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(5, 15, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    ok_(r1.encloses(r2))


def test_merge():

    r1 = SequenceRange(0, 11, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(5, 22, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    m = r1.merge_with(r2)

    eq_(m.start, 0)
    eq_(m.end, 22)


def test_hashable():

    r1 = SequenceRange(0, 20, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    r2 = SequenceRange(5, 15, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

    d = {
        r1: 1,
        r2: 2
    }

    eq_(d[r1], 1)
    eq_(d[r2], 2)

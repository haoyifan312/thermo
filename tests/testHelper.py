DEL = 1.0e-6
import numpy

def verify_rel(answer, expected, delta=DEL):
    if expected == 0.0:
        return answer == 0.0
    else:
        if abs((answer-expected)/expected) < delta:
            return True
        else:
            return False


def verify_abs(answer, expected, delta=DEL):
    if expected == 0.0:
        return answer == 0.0
    else:
        if abs(answer-expected) < delta:
            return True
        else:
            return False


def verify_vec_abs(answer, expected, delta=DEL):
    if len(answer) != len(expected):
        return False

    for i, element in enumerate(answer):
        if type(element) != type(expected[i]):
            return False

        if isinstance(element, numpy.ndarray):
            return verify_vec_abs(element, expected[i], delta)
        else:
            return verify_abs(element, expected[i], delta)


def verify_vec_rel(answer, expected, delta=DEL):
    if len(answer) != len(expected):
        return False

    for i, element in enumerate(answer):
        if type(element) != type(expected[i]):
            return False

        if isinstance(element, numpy.ndarray):
            return verify_vec_rel(element, expected[i], delta)
        else:
            return verify_rel(element, expected[i], delta)
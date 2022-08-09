import numpy as np


def random_unit():
    """ Returns a random vector with unit length. """
    r = np.random.rand(3) - 0.5
    return r / np.linalg.norm(r)


def calculate_distance(*coordinate):
    if len(coordinate) == 1:
        assert isinstance(coordinate[0], (list, np.ndarray)) and len(coordinate[0]) == 3
        return np.linalg.norm(coordinate[0])
    elif len(coordinate) == 2:
        assert isinstance(coordinate[0], (list, np.ndarray)) and isinstance(coordinate[1], (list, np.ndarray))
        assert len(coordinate[0]) == len(coordinate[1]) and len(coordinate[1]) == 3
        d1 = np.array(coordinate[0])
        d2 = np.array(coordinate[1])
        return np.linalg.norm(d1 - d2)
    else:
        raise ValueError('too much input!')


def calculate_angle(co1, co2, co3):
    assert len(co1) == 3 and len(co2) == 3 and len(co3) == 3
    co1 = np.array(co1)
    co2 = np.array(co2)
    co3 = np.array(co3)
    return angle_between_arrays(co1 - co2, co3 - co2) * 180 / np.pi


def calculate_dihedral(co1, co2, co3, co4):
    assert len(co1) == 3 and len(co2) == 3 and len(co3) == 3 and len(co4) == 3
    co1 = np.array(co1)
    co2 = np.array(co2)
    co3 = np.array(co3)
    co4 = np.array(co4)
    v1 = co2 - co1
    v2 = co3 - co2
    v3 = co4 - co3
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    sig = np.sign(np.dot(np.cross(n2, n1), v2))
    return angle_between_arrays(n1, n2) * 180 / np.pi * sig


def calculate_dihedral_no_sign(co1, co2, co3, co4):
    assert len(co1) == 3 and len(co2) == 3 and len(co3) == 3 and len(co4) == 3
    co1 = np.array(co1)
    co2 = np.array(co2)
    co3 = np.array(co3)
    co4 = np.array(co4)
    v1 = co2 - co1
    v2 = co3 - co2
    v3 = co4 - co3
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    return angle_between_arrays(n1, n2) * 180 / np.pi


def calculate_square(dif):
    return np.sum(np.square(dif))


def normalization(array):
    return np.array(array) / np.linalg.norm(array)


def angle_between_arrays(array1, array2):
    c = np.dot(normalization(array1), normalization(array2))
    if c > 1 or c < -1 :
        print("unreasonable cos value", c)
        c = np.round(c)
    return np.arccos(c)


def rotation_dot(dot_xyz, centroid, array1, array2):
    # rotate the dot around the centroid, from array1 to array2
    n_vector = normalization(np.cross(normalization(array1), normalization(array2)))
    theta = angle_between_arrays(array1, array2)
    rotate_vector = dot_xyz - centroid
    new_dot = centroid + rotate_vector * np.cos(theta) + (np.cross(n_vector, rotate_vector)) * np.sin(
        theta) + n_vector * (
                  np.dot(n_vector, rotate_vector)) * (1. - np.cos(theta))
    return new_dot


def cal_n_vector_theta(array1, array2):
    n_vector = normalization(np.cross(normalization(array1), normalization(array2)))
    theta = angle_between_arrays(array1, array2)
    return n_vector, theta


def rotate_n_theta(old_dot, centroid, n_vector, theta):
    rotate_vector = old_dot - centroid
    new_dot = centroid + rotate_vector * np.cos(theta) + (np.cross(n_vector, rotate_vector)) * np.sin(
        theta) + n_vector * (
                  np.dot(n_vector, rotate_vector)) * (1. - np.cos(theta))
    return new_dot


def matrix_distance(list1, list2, cut_off_distance=15., black_list=None):
    m = len(list1)
    n = len(list2)
    black_matrix = np.mat(np.ones((m, n)))
    if black_list:
        for black in black_list:
            black_matrix[black] = 0.
    matrix1_list = [np.mat(np.zeros(shape=(m, n))), np.mat(np.zeros(shape=(m, n))), np.mat(np.zeros(shape=(m, n)))]
    matrix2_list = [np.mat(np.zeros(shape=(m, n))), np.mat(np.zeros(shape=(m, n))), np.mat(np.zeros(shape=(m, n)))]
    for i in range(3):
        for index, point in enumerate(list1):
            matrix1_list[i][index, :] = point[i]
    for i in range(3):
        for index, point in enumerate(list2):
            matrix2_list[i][:, index] = point[i]
    difference = [matrix1_list[i] - matrix2_list[i] for i in range(3)]
    abs_dif = [np.abs(difference[i]) for i in range(3)]
    b = [np.where(abs_dif[i] <= cut_off_distance, 1, 0) for i in range(3)]
    judge = np.multiply((b[0] * b[1] * b[2]), black_matrix)
    sum_dif = np.square(np.multiply(difference[0], judge)) + np.square(np.multiply(difference[1], judge)) + np.square(
        np.multiply(difference[2], judge))
    return sum_dif


def matrix_distance_range(list1, list2, cut_off_range=(0., 15.), black_list=None):
    m = len(list1)
    n = len(list2)
    black_matrix = np.mat(np.ones((m, n)))
    if black_list:
        for black in black_list:
            black_matrix[black] = 0.
    matrix1_list = [np.mat(np.zeros(shape=(m, n))), np.mat(np.zeros(shape=(m, n))), np.mat(np.zeros(shape=(m, n)))]
    matrix2_list = [np.mat(np.zeros(shape=(m, n))), np.mat(np.zeros(shape=(m, n))), np.mat(np.zeros(shape=(m, n)))]
    # print("called1")
    for i in range(3):
        for index, point in enumerate(list1):
            matrix1_list[i][index, :] = point[i]
    for i in range(3):
        for index, point in enumerate(list2):
            matrix2_list[i][:, index] = point[i]
    # print("called2")
    difference = [matrix1_list[i] - matrix2_list[i] for i in range(3)]
    abs_dif = [np.abs(difference[i]) for i in range(3)]
    b = [np.where(abs_dif[i] <= cut_off_range[1], 1, 0) for i in range(3)]
    judge = np.multiply((b[0] * b[1] * b[2]), black_matrix)
    # print("called3")

    sum_dif = np.square(np.multiply(difference[0], judge)) + np.square(np.multiply(difference[1], judge)) + np.square(
        np.multiply(difference[2], judge))
    # print("called4")

    return sum_dif


def get_longest_path(all_path):
    assert isinstance(all_path, dict)
    longest_path = set()
    chain_length = 0
    for each_path in all_path.values():
        assert isinstance(each_path, list)
        if len(each_path) > chain_length:
            chain_length = len(each_path)
            longest_path.clear()
            longest_path.add(tuple(each_path))
        elif len(each_path) == chain_length:
            longest_path.add(tuple(each_path))
    return longest_path


def d_2(co1, co2):
    _co = co1 - co2
    return np.inner(_co, _co)

# def angle_between_arrays(array1, array2):
#     return np.arccos(np.dot(unitizaion(array1), unitizaion(array2)))
#
#
# def rotation_dot(dot_xyz, centroid, array1, array2):
#     # rotate the dot around the centroid, from array1 to array2
#     n_vector = np.cross(unitizaion(array1), unitizaion(array2))
#     theta = angle_between_arrays(array1, array2)
#     rotate_vector = dot_xyz - centroid
#     new_dot = centroid + rotate_vector * np.cos(theta) + (np.cross(n_vector, rotate_vector)) * np.sin(theta) + n_vector * (
#         np.dot(n_vector, rotate_vector)) * (1. - np.cos(theta))
#     return new_dot

if __name__ == "__main__":
    # a = [np.array([1., 5., 1.]), np.array([1., 2., 1.])]
    # b = [np.array([1., 1., 1.]), np.array([1., 15., 3.]), np.array([1., 2., 1.])]
    # sum = matrix_distance(a, b)
    # print(sum)
    # import MDWorkUp_old.parameter_calculation as mc
    #
    # count = pm.count_points(sum, 10, 10)
    # print(count)
    import random


    def test():
        # ro = np.array([random.random() for _ in range(1000)])
        # r1 = np.array([random.random() for _ in range(1000)])
        ro = np.arange(1000)
        r1 = np.arange(1000)
        # r1 = np.array([2.23453, -0.3452452, 1.23452345]*1000)
        ro = np.sqrt(np.square(ro) + np.square(r1))
        pass


    import timeit

    t = timeit.timeit(stmt='test()', setup='from __main__ import test', number=10000)
    print(t)

    # print(calculate_distance(ro))
    # cent = np.array([0, 0, 0])
    # a = np.array([2.54139987, -0.96383034, 1.45959034])
    # b = np.array([1, 0, 0])
    # print(rotation_dot(ro, cent, b, a))
    # pass

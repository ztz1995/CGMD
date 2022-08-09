import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import special as spc
from pyMD import cartesian_operation as cso
from scipy import integrate


def find_maximum(y, num=2, cal_range=25, sort=False):
    max_idx_old = -1
    max_idx = list()
    count = 0
    for i in range(len(y) - cal_range + 1):
        piece = y[i:i + cal_range]
        if isinstance(piece, np.ndarray):
            max_idx_new = piece.argmax() + i
        else:
            max_idx_new = piece.idxmax()
        # print(max_idx_new, count)
        if max_idx_new == max_idx_old:
            count += 1
            if count == cal_range // 2:
                max_idx.append(max_idx_new)
        else:
            max_idx_old = max_idx_new
            count = 0
    if sort:
        max_idx.sort(key=lambda idx: abs(y[idx]), reverse=True)
    return max_idx[:num]


def harmonic(x, *parameters):
    k, r = parameters
    return k * np.square(x - r)


def harmonic_d(x, *parameters):
    k, d = parameters
    theta = x / 180 * np.pi
    return k * (1 + d * np.cos(theta))


def harmonic_p(x, k):
    theta = x / 180 * np.pi
    return k * (1 + np.cos(theta))


def harmonic_n(x, k):
    theta = x / 180 * np.pi
    return k * (1 - np.cos(theta))


def multi_harmonic(x, *parameters):
    theta = x / 180 * np.pi
    ret = np.zeros(x.shape)
    for i in range(5):
        ret += parameters[i] * np.power(np.cos(theta), i)
    return ret


def fene_expand(x, *parameters):
    k, delta, r0 = parameters
    return -0.5 * k * r0 * r0 * np.log(1 - np.square((x - delta) / r0))


def fit_bond(x, E, style, p0=None, d=0., k=1.5):
    if style == "harmonic":
        func = harmonic
        if p0 is None:
            p = [10, x[np.argmin(E)]]
        else:
            p = [p0[0], p0[1] + d / 2.]
        bound_l = [0., p[1] - 0.2]
        bound_r = [np.inf, p[1] + 0.2]
    elif style == "fene/expand":
        func = fene_expand
        if p0 is None:
            p = [10., x[np.argmin(E)], x[np.argmin(E)] * k]
        else:
            p = [p0[0], p0[1] + d / 2, (p0[1] + d / 2) * k]
        bound_l = [0., p[1] * 0.98, p[1] * 0.95 * k]
        bound_r = [np.inf, p[1] * 1.02, p[1] * 1.05 * k]
    else:
        raise KeyError("No such style!")
    # print(p)
    # print(bound_l, bound_r)
    try:
        popt, _ = curve_fit(func, x, E, p0=p, bounds=(bound_l, bound_r), maxfev=int(1e5))
    except ValueError:
        print(p, bound_l, bound_r)
        raise KeyboardInterrupt

    return popt


def fit_angle(x, E, style, p0=None, d=0.):
    xr = x / 180 * np.pi
    if style == "harmonic":
        func = harmonic
        if p0 is None:
            p = [10, xr[np.argmin(E)]]
        else:
            p = [p0[0], (p0[1] + d / 2.) / 180. * np.pi]
        bound_l = [0., p[1] * 0.9]
        bound_r = [np.inf, p[1] * 1.1]
    else:
        raise KeyError("No such style!")
    # print(p)
    # plt.plot(x, E)
    # plt.show()
    popt, _ = curve_fit(func, xr, E, p0=p, bounds=(bound_l, bound_r), maxfev=int(1e5))
    return [popt[0], popt[1] / np.pi * 180]


def fit_dihedral(x, E, style, p0=None, d=0.):
    if style == "multi/harmonic":
        func = multi_harmonic
        if p0 is None:
            p = [1., 1., 1., 1., 1.]
        else:
            p = p0
        bound_l = [-np.inf, -np.inf, -np.inf, -np.inf, -np.inf]
        bound_r = [np.inf, np.inf, np.inf, np.inf, np.inf]

    elif style == "harmonic_d":
        if x[np.argmin(E)] > 90:
            func = harmonic_p
        else:
            func = harmonic_n
        if p0 is None:
            p = [10]
        else:
            p = [p0[0]]
        bound_l = [0.]
        bound_r = [np.inf]
    else:
        raise KeyError("No such style!")
    popt, _ = curve_fit(func, x, E, p0=p, bounds=(bound_l, bound_r), maxfev=int(1e5))
    if func == harmonic_p:
        if p0 is None:
            return [popt[0], 1]
        else:
            return [p0[0] * d, 1]
    elif func == harmonic_n:
        if p0 is None:
            return [popt[0], -1]
        else:
            return [p0[0] * d, -1]
    else:
        return popt


def pow2_9(x, *parameters):
    k1, r, b, k2 = parameters
    return k1 * pow(x - r, 2) + b + k2 * np.power(x / r, -9)

def pow2(x, *params):
    k, r = params
    judge = np.where(x < r, 1, 0)
    return k * np.square(x/r-1) * judge


def pow9_6_1(x, *parameters):
    k1, k2, r = parameters
    return k1 / x + k2 * (2 * np.power(x / r, -9) - 3 * np.power(x / r, -6))


def fit_non_bond(x, E, l, aa=False):
    if aa:
        fit_func = pow9_6_1
        if min(E) < -5.:
            p0 = [-5., 1., 5.]
        else:
            p0 = [5, 1., 5.]
        bound_l = [-np.inf, 0.001, 0.1]
        bound_r = [np.inf, 10., 10.]
        popt, _ = curve_fit(fit_func, x[150:500], E[150:500], p0=p0, bounds=(bound_l, bound_r), maxfev=int(1e5))
        fit_e = fit_func(x, *popt)
        fit_e = np.append(fit_e[:300] - fit_e[299] + E[299], E[300:])
        #
        print(popt)
        plt.plot(x[l + 10:], fit_e[l + 10:])
        plt.plot(x[l + 10:], E[l + 10:])
        plt.show()
    else:
        fit_func = pow2_9
        p0 = [0.5, 4.5, 0., 1]
        bound_l = [0.1, 0.5, -20, 0.1]
        bound_r = [np.inf, 6., 10, 10]
        popt, _ = curve_fit(fit_func, x[l + 10:l + 50], E[l + 10:l + 50], p0=p0,
                            bounds=(bound_l, bound_r), maxfev=int(1e5))
        fit_e = fit_func(x[:l + 25], *popt)
        fit_e = np.append(fit_e - fit_e[-1] + E[l + 24], E[l + 25:])
    return fit_e


def cal_derivative(x, y):
    y2 = np.append(y[1:], y[-1])
    y1 = np.append(y[0], y[:-1])
    x2 = np.append(x[1:], x[-1])
    x1 = np.append(x[0], x[:-1])
    return (y2 - y1) / (x2 - x1)


def error_l2(present_rdf, target_rdf):
    numerator = sum(np.square(present_rdf - target_rdf))
    denominator = sum(np.square(target_rdf))
    return np.sqrt(numerator / denominator)


def error_l1(present_rdf, target_rdf):
    numerator = sum(abs(present_rdf - target_rdf))
    denominator = sum(abs(present_rdf) + abs(target_rdf))
    return numerator / denominator


def linear_func(x, *parameters):
    k, b = parameters
    return k * x + b


def smooth_step(points, x_l, x_r, reverse=False):
    judge_l = np.where(points <= x_l, 1, 0)
    judge_r = np.where(points > x_r, 1, 0)
    judge_m = np.ones(points.shape) - judge_l - judge_r
    if not reverse:
        y = (x_r - points) / (x_r - x_l)
        return judge_l + y * y * (3 - 2 * y) * judge_m
    else:
        y = (points - x_l) / (x_r - x_l)
        return judge_r + y * y * (3 - 2 * y) * judge_m


def smooth_charmm(x, x_l, x_r):
    judge_l = np.where(x <= x_l, 1, 0)
    judge_r = np.where(x > x_r, 1, 0)
    judge_m = np.ones(x.shape) - judge_l - judge_r
    xr2 = x_r * x_r
    xl2 = x_l * x_l
    x2 = x * x
    y = np.square(xr2 - x2) * (xr2 + 2 * x2 - 3 * xl2) / np.power(xr2 - xl2, 3)
    return judge_l + y * judge_m


def lj96(_r, e1, e2, s1, s2):
    e = (2 * np.sqrt(e1 * e2) * np.power(s1 * s2, 3)) / (np.power(s1, 6) + np.power(s2, 6))
    s = np.power((np.power(s1, 6) + np.power(s2, 6)) / 2, 1 / 6)
    return e * (2 * np.power(s / _r, 9) - 3 * np.power(s / _r, 6))


def lj_amber(_r, e1, e2, s1, s2):
    e = np.sqrt(e1 * e2)
    s = 0.5 * (s1 + s2) * 2
    return 4 * e * (np.power(s / _r, 12) - np.power(s / _r, 6))


def coulomb(_r, _q1, _q2):
    return 332.0647 * _q1 * _q2 / _r


def coulomb_dsf(_r, _q1, _q2, rc, a):
    judge = np.where(_r <= rc, 1., 0.)
    return (332.0647 * _q1 * _q2 * (spc.erfc(a * _r) / _r - spc.erfc(a * rc) / rc + (
            spc.erfc(a * rc) / rc / rc + 2 * a / np.sqrt(np.pi) * np.exp(-a * a * rc * rc) / rc) * (_r - rc))) * judge


def cal_aa_force(points, type_tuple, cut_off=12., q=False, dielectric=1.0, dsf=False, a=0.01):
    # COMPASS
    # q_dict = {"n1": -0.462, "h1": 0.351, "c3": 0.532, "o1": -0.5}
    # q_dict = {"n1": -0.367, "h1": 0.351, "c3": 0.532, "o1": -0.5}
    # epsilon_dict = {"n1": 0.2, "h1": 0.01, "c3": 0.064, "o1": 0.192}
    # sigma_dict = {"n1": 3.7, "h1": 1.45, "c3": 3.9, "o1": 3.43}

    # GAFF urea
    # q_dict = {"n1": -0.888, "h1": 0.388, "c3": 0.884, "o1": -0.66}
    # q_dict = {"n1": -0.595, "h1": 0.388, "c3": 0.884, "o1": -0.66}
    # q_dict = {"n1": -0.5, "h1": 0.388, "c3": 0.884, "o1": -0.66}
    # epsilon_dict = {"n1": 0.17, "h1": 0.0157, "c3": 0.086, "o1": 0.21}
    # sigma_dict = {"n1": 1.824, "h1": 0.6, "c3": 1.9080, "o1": 1.6612}

    # q_dict = {"n1": -0.595, "h1": 0.388, "c1": 0.884, "o1": -0.66}
    # epsilon_dict = {"n1": 0.17, "h1": 0.0157, "c1": 0.086, "o1": 0.21}
    # sigma_dict = {"n1": 1.824, "h1": 0.6, "c1": 1.9080, "o1": 1.6612}

    # new_type
    epsilon_dict = {"n1": 0.2, "h1": 0.01, "c1": 0.064, "o1": 0.192}
    sigma_dict = {"n1": 3.7, "h1": 1.45, "c1": 3.9, "o1": 3.43}
    # q_dict = {"n1": -0.595, "h1": 0.388, "c1": 0.884, "o1": -0.66}
    q_dict = {"n1": -0.462, "h1": 0.351, "c1": 0.532, "o1": -0.5}

    t1 = type_tuple[0]
    t2 = type_tuple[1]
    vdw_f = lj96(points, epsilon_dict[t1], epsilon_dict[t2], sigma_dict[t1], sigma_dict[t2])
    # vdw_f = lj_amber(points, epsilon_dict[t1], epsilon_dict[t2], sigma_dict[t1], sigma_dict[t2])
    ret = vdw_f * smooth_step(points, cut_off - 2, cut_off)
    if q:
        if dsf:
            c_f = coulomb_dsf(points, q_dict[t1], q_dict[t2], cut_off, a) / dielectric
        else:
            c_f = coulomb(points, q_dict[t1], q_dict[t2]) / dielectric
        return c_f + ret
    else:
        return ret


def pdf_uniform(_x):
    return 1


def pdf_gaussian(_x, mu, s_2):
    _y = np.exp(-np.square(_x - mu) / s_2) / np.sqrt(2 * np.pi * s_2)
    return _y


def pdf_cauchy(_x, x0, r):
    return 1 / np.pi * (r / (r * r + np.square(_x - x0)))


def single_gaussian(_x, k, mu, s2, b):
    return k / np.sqrt(2 * np.pi * s2) * np.exp(-np.square(_x - mu) / s2) + b


def random_gen(fx, x_min, x_max, seed=None):
    if seed is not None:
        np.random.seed(seed)
    while True:
        a, b = np.random.rand(2)
        a = a * (x_max - x_min) + x_min
        if b < fx(a):
            return a


def create_fix_chain(r1, r2, b, atom_num):
    def cal_p(co):
        ks = s * (N - s) / N
        rs = r1 * (N - s) / N + r2 * s / N
        p = np.power(3 / 2 / np.pi / ks / b / b, 1.5) * np.exp(-3 / 2 / ks / b / b * cso.d_2(co, rs))
        return p

    r1 = np.asarray(r1)
    r2 = np.asarray(r2)
    N = atom_num + 1
    last_co = r1
    co_list = list()
    for s in range(1, atom_num + 1):
        new_co = last_co + b * cso.random_unit()
        # print(cal_p(new_co))
        while cal_p(new_co) < np.random.rand():
            new_co = last_co + b * cso.random_unit()
            # print(cal_p(new_co))

        co_list.append(new_co)
        last_co = new_co
    return co_list


def free_joint(r1, r2, co, s, N, b):
    ks = s * (N - s) / N
    rs = r1 * (N - s) / N + r2 * s / N
    p = np.power(3 / 2 / np.pi / ks / b / b, 1.5) * np.exp(-3 / 2 / ks / b / b * cso.d_2(co, rs))
    return p


def create_fix_end_rotate_chain(p_func, r1, r2, b, atom_num, angle=109.5):
    def new_co_gen(_l_co, _ll_co=None):
        if _ll_co is None:
            return _l_co + b * cso.random_unit()
        else:
            l_vec = _l_co - _ll_co
            ran_vec = cso.random_unit()
            n_vec = np.cross(l_vec, ran_vec)
            n_vec /= np.linalg.norm(n_vec)
            _new_co = l_vec * cost + n_vec * b * sint + _l_co
            # print(np.inner(l_vec, _new_co - _l_co)/b/b)
            return _new_co

    cost = np.sin(angle / 180 * np.pi)
    sint = -np.cos(angle / 180 * np.pi)
    r1 = np.asarray(r1)
    r2 = np.asarray(r2)
    N = atom_num + 1
    l_co = r1
    ll_co = None
    co_list = list()
    for s in range(1, atom_num + 1):
        print(s)
        new_co = new_co_gen(l_co, ll_co)
        while p_func(r1, r2, new_co, s, N, b) < np.random.rand():
            new_co = new_co_gen(l_co, ll_co)
            # print(p_func(r1, new_co, s, b) * p_func(r2, new_co, N-s, b))
            # print(p_func(r1, r2, new_co, s, N, b))
            # print(p_func(r2, new_co, N-s, b))
        co_list.append(new_co)
        ll_co = l_co
        l_co = new_co
    return co_list


class ChainCreator:

    # def __init__(self, r1, r2, atom_num, b, angle=None, model=0):
    def __init__(self, max_N, b, angle=None, model=0):
        # self.r1 = np.asarray(r1)
        # self.r2 = np.asarray(r2)
        # self.atom_num = atom_num
        self.b = b
        if angle is None and model == 1:
            self.angle = 109.5
        else:
            self.angle = angle
        if self.angle is not None:
            self.cost = np.sin(self.angle / 180 * np.pi)
            self.sint = -np.cos(self.angle / 180 * np.pi)

        self.model = model
        self.ro_params = dict()
        self.p_rotate_r = dict()
        if model == 0:
            self.max_r = self.b * max_N
        elif model == 1:
            self.max_r = self.b * max_N * np.cos(self.angle / 2 / 180 * np.pi)
        # print(self.max_r)
        self.r_list = np.arange(0, self.max_r, 0.005) + 1E-8

    def p_free_joint_chain(self, org, r, N):
        p = np.power(3 / 2 / np.pi / N / self.b / self.b, 1.5) * np.exp(-3 / 2 / N / self.b / self.b * cso.d_2(org, r))
        return p

    def cal_ro_param(self, n):
        if n not in self.ro_params:
            n3 = np.power(3, n)
            l2 = self.b * self.b
            R2 = (2 * n - 1.5 + 1.5 / n3) * l2
            R4 = (20 / 3 * n * n - 64 / 3 * n + 57 / 2 - 10 * n / n3 - 57 / 2 / n3) * l2 * l2
            R6 = (280 / 9 * n * n * n - 686 / 3 * n * n + 7240 / 9 * n - 21667 / 18 - 4 / 9 *
                  np.power(-1, n) + 70 / n3 * n * n + 518 / n3 * n + 7225 / 6 / n3) * l2 * l2 * l2
            self.ro_params[n] = (R2, R4, R6)
        return self.ro_params[n]

    def p_free_rotate_chain(self, org, r, n):
        def integrand(k, R):
            S = np.exp(-R2 / 6 * np.power(k, 2) - (R2 ** 2 / 72 - R4 / 120) * np.power(k, 4) - (
                    R2 ** 3 / 648 - R2 * R4 / 720 + R6 / 5040) * np.power(k, 6))
            f = S * np.sin(k * R) / R * k * 4 * np.pi
            return f

        R2, R4, R6 = self.cal_ro_param(n)
        if n not in self.p_rotate_r:
            p_list = list()
            print("called")
            for ri in self.r_list:
                # print("called")
                p_list.append(max(integrate.quad(integrand, 0, np.inf, (ri,))[0], 0))
            p_list = np.asarray(p_list)
            p_list /= np.sum(p_list) * (self.r_list[1] - self.r_list[0])
            self.p_rotate_r[n] = p_list

        R = cso.calculate_distance(org, r)
        return np.interp(R, self.r_list, self.p_rotate_r[n])

    def new_co_free_joint(self, l_co, ll_co=None):
        return l_co + self.b * cso.random_unit()

    def new_co_free_rotate(self, l_co, ll_co=None):
        if ll_co is None:
            return l_co + self.b * cso.random_unit()
        else:
            l_vec = l_co - ll_co
            ran_vec = cso.random_unit()
            n_vec = np.cross(l_vec, ran_vec)
            n_vec /= np.linalg.norm(n_vec)
            _new_co = l_vec * self.cost + n_vec * self.b * self.sint + l_co
            # print(np.inner(l_vec, _new_co - _l_co)/b/b)
            return _new_co

    def evaluate(self, p_func, gen_func, l_co, ll_co, s, N, r1, r2):
        p_max = -1
        co = None
        for i in range(100):
            new_co = gen_func(l_co, ll_co)
            p = p_func(r1, new_co, s) * p_func(r2, new_co, N - s)
            # print(i)
            # print(new_co)
            # print(p_func(r1, new_co, s))
            # print(p_func(r2, new_co, N - s))
            if p > p_max:
                p_max = p
                co = new_co
        return p_max, co

    def gen_new_co(self, p_func, gen_func, l_co, ll_co, s, N, p_max, r1, r2):
        new_co = gen_func(l_co, ll_co)
        p = p_func(r1, new_co, s) * p_func(r2, new_co, N - s)
        while p < p_max * np.random.rand():
            new_co = gen_func(l_co, ll_co)
            p = p_func(r1, new_co, s) * p_func(r2, new_co, N - s)
        return new_co

    def create_chain(self, r1, r2, atom_num):
        if self.model == 0:
            p_func = self.p_free_joint_chain
            gen_func = self.new_co_free_joint
        elif self.model == 1:
            p_func = self.p_free_rotate_chain
            # p_func = self.p_free_joint_chain
            gen_func = self.new_co_free_rotate
            # gen_func = self.new_co_free_joint

        else:
            raise KeyError("No such model!")
        r1 = np.asarray(r1)
        r2 = np.asarray(r2)
        l_co = r1
        ll_co = None
        co_list = list()
        for s in range(1, atom_num + 1):
            p_max, new_co = self.evaluate(p_func, gen_func, l_co, ll_co, s, atom_num + 1, r1, r2)
            # print(p_max)
            new_co = self.gen_new_co(p_func, gen_func, l_co, ll_co, s, atom_num + 1, p_max * 10, r1, r2)
            # print(new_co)
            co_list.append(new_co)
            ll_co = l_co
            l_co = new_co
        return co_list


def sin_shift(x, A, w, fi):
    return A * np.sin((w * x + fi))


if __name__ == '__main__':
    # d = 1
    # begin = 1.2
    # end = 12.
    # di = 1.
    # tt1 = ("h1", "o1")
    # tt2 = ("o1", "n1")
    # x1 = np.arange(begin, end, 0.001)
    # # y1 = cal_aa_force(x1, tt2, True, q=True, dielectric=di)
    # y1 = cal_aa_force(x1, tt1, True, q=True, dielectric=di, dsf=True, rc=12., a=0.01)
    #
    # x2 = np.arange(begin + d, end + d, 0.001)
    # y2 = cal_aa_force(x2, tt2, True, q=True, dielectric=di, dsf=True, rc=12., a=0.01)
    #
    # # x2 = np.arange(begin + d, end + d, 0.001)
    # # y2 = cal_aa_force(x2, ("o1", "n1"), True, q=True, dielectric=di)
    # import matplotlib.pyplot as plt
    #
    # plt.plot(x1, y1 + y2)
    # plt.show()
    # x = np.arange(1, 15, 0.1)
    # y = smooth_charmm(x, 10., 12.)
    # y2 = smooth_step(x, 10., 12.)
    # plt.plot(x, y)
    # plt.plot(x, y2)
    # plt.show()

    # l = create_fix_end_rotate_chain(free_joint, [0, 0, 0], [30, 0, 0], 1.54, 20)
    # print(l)

    atom_creator = ChainCreator(50, 1.54, model=1)
    l_l = atom_creator.create_chain([0, 0, 0], [10, 0, 0], 50)
    print(l_l)

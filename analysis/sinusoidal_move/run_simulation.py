import os
import numpy as np
from multiprocessing import Process
import time
from pyMD import file_parser as fps


def prepare_path(_root_path, _cal_range):
    os.makedirs(_root_path + "cal/", exist_ok=True)
    os.system("cp %sinit/* %scal/" % (_root_path, _root_path))
    with open(_root_path + "cal/in.template", "r") as f:
        in_lines = f.readlines()

    for tp in _cal_range:
        v_lines = ["# Variables\n", "variable        tp equal %d\n" % tp, "variable        loops equal 40\n"]
        with open(_root_path + "cal/in.%d" % tp, "w") as f:
            f.writelines(v_lines)
            f.writelines(in_lines)


def invoke_lmp_process(file_path, in_name, num_parallel, screen=False):
    p_lmp = Process(target=fps.invoke_lmp, args=(file_path, in_name, num_parallel, screen))
    p_lmp.start()
    start = time.time()
    while True:
        if p_lmp.is_alive():
            time.sleep(1)
        else:
            if not os.path.exists(file_path + "finish." + in_name):
                print("%s, failed" % in_name)
                return False
            else:
                end = time.time()
                dt = int(end - start)
                print("%s, succeeded, time: %dh %dm %ds" % (in_name, dt // 3600, dt // 60 % 60, dt % 60))
                return True


def simulate(_root_path, _cal_range, num_parallel=24):
    for tp in _cal_range:
        invoke_lmp_process(_root_path + "cal/", "in.%d" % tp, num_parallel, screen=True)


if __name__ == '__main__':
    work_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/sin_range/"
    cal_range = np.logspace(4.1, 5, num=9, dtype=int, endpoint=False)
    prepare_path(work_path, cal_range)
    simulate(work_path, cal_range, num_parallel=12)

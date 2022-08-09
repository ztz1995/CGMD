if __name__ == '__main__':

    fp = "/home/centos/ztz/SH6S/cgu_1blk_2000chain/1blk_1E-6_test/dump.param"
    start = 0
    nevery1 = 10
    n1 = 10
    nevery2 = 50000
    total = 4000000
    dump_ids = list()
    for i in range(start, total+nevery2, nevery2):
        for j in range(n1):
            dump_ids.append(str(int(i + nevery1*j))+"\n")

    with open(fp, "w") as file:
        file.writelines(dump_ids[1:])

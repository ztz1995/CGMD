
def cal_segment_para(mole_num, lines):
    mole_id = list()
    segment_id = list()
    for i in range(mole_num):
        mole_id.append(i+1)
        s_id = int(i // (mole_num // lines)) + 1
        segment_id.append(s_id)
    return mole_id, segment_id


def output(fp, mole_id, segment_id):
    file_lines = ["# mole id -> segment id\n"]
    for i in range(len(mole_id)):
        file_lines.append("%d %d\n" % (mole_id[i], segment_id[i]))
    with open(fp, "w") as file:
        file.writelines(file_lines)
    pass


if __name__ == '__main__':

    name = "para2"
    fp = "./segment.%s" % name
    mole_num = 40
    lines = 2
    mole_id, segment_id = cal_segment_para(mole_num, lines)
    output(fp, mole_id, segment_id)
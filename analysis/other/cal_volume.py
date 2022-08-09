import pandas as pd

if __name__ == '__main__':

    fp = "log.csv"
    data = pd.read_csv(fp, sep="\\s+")
    print(pow(data.Volume.mean(), 1/3))
    # 43.986
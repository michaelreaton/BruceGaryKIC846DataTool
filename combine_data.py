import csv
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.dates as mdates
import datetime


def main():
    combined = read_all_csvs()
    combined = clean_bad_data(combined)
    combined = sorted(combined,key=getKey)
    print('#observations:',len(combined))
    combined = [['MJD','V-mag','SE','air mass']] + combined
    write_combined_csv(combined)
    raw_plot(combined)

def raw_plot(combined):
    x = []
    y = []
    for i in range(1,len(combined)):
        x.append(float(combined[i][0]))
        y.append(float(combined[i][1]))
    x = mjddates_to_gregoriandates(x)
    plt.xlabel("Date")
    plt.ylabel("V-Mag")
    start_date = x[0]
    end_date = x[len(x)-1]
    title = "Bruce Gary Raw Data (" + str(start_date) + " to " + str(end_date) + ")"
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.gcf().autofmt_xdate()
    plt.grid()
    plt.scatter(x,y,s=1,marker='s')
    plt.savefig("SavedPlots/scatter.png")
    plt.show()

def mjddates_to_gregoriandates(mjd):
    dates = [mjd_to_gregorian(mjd_val) for mjd_val in mjd]
    return dates

def mjd_to_gregorian(mjd_val):
    start_date = datetime.date(1858, 11, 17)
    return start_date + datetime.timedelta(mjd_val)

def clean_bad_data(combined):
    clean = []
    bad_data_point_counter = 0
    for i in range(len(combined)):
        if(len(combined[i]) == 4 and is_float(combined[i][0]) and is_float(combined[i][1])):
            clean.append(combined[i])
        else:
            bad_data_point_counter+=1
    print("bad data points discarded:", bad_data_point_counter)
    return clean

def write_combined_csv(list):
    with open("Combined/bruce_gary_raw_data_combined.csv", "w") as data_file:
        writer = csv.writer(data_file,lineterminator='\n')
        writer.writerows(list)

def read_all_csvs():
    files = glob.glob('RawData/*')
    combined = []
    for file in files:
        combined = combined + read_csv(file)
    return combined

def getKey(item):
    return item[0]

def read_csv(path):
    with open(path, 'r') as data_file:
        # skip the header rows
        input_csv = csv.reader(data_file, delimiter='\t')
        for hr in range(1,13):
            next(input_csv)
        return list(input_csv)

def is_float(value):
  try:
    float(value)
    return True
  except:
    return False

if __name__ == "__main__":
    main()

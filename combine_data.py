import csv
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.dates as mdates
import datetime


def main():
    header = [['MJD','V-mag','SE','air mass']]
    combined = read_all_csvs()
    combined = clean_bad_data(combined)
    combined = sorted(combined,key=getKey)
    combined_good_air = filter_by_air_mass(combined)
    print('#observations:',len(combined))
    print('#observations air mass <= 2.5:',len(combined_good_air))

    write_csv(header+combined,"bruce_gary_raw_data_combined")
    scatter_plot(combined,plot_name="scatter", plot_title="Bruce Gary Raw Data Unmodified")
    write_csv(header+combined_good_air,"bruce_gary_raw_data_good_air_combined")
    scatter_plot(combined_good_air,plot_name="scatter_good_air", plot_title="Bruce Gary Raw Data air mass <= 2.5")


def filter_by_air_mass(data):
    filtered = []
    for i in range(1,len(data)):
        if(float(data[i][3]) <= 2.5):
            filtered.append(data[i])
    return filtered

def scatter_plot(combined, plot_name, plot_title):
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
    title = plot_title + "\n(" + str(start_date) + " to " + str(end_date) + ")"
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.gcf().autofmt_xdate()
    plt.grid()
    plt.scatter(x,y,s=1,marker='s')
    plt.savefig("SavedPlots/" + plot_name + ".png")
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
        if(len(combined[i]) == 4 and is_float(combined[i][0]) and is_float(combined[i][1]) and is_float(combined[i][2]) and is_float(combined[i][3])):
            clean.append(combined[i])
        else:
            bad_data_point_counter+=1
    print("bad data points discarded:", bad_data_point_counter)
    return clean

def write_csv(list, name):
    with open("Combined/" + name + ".csv", "w") as data_file:
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

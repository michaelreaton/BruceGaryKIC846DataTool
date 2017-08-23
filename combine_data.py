import csv
import glob
import matplotlib.pyplot as plt
import math
import datetime


def main():
    header = [['MJD','V-mag','SE','air mass']]
    header2 = [['Avg MJD','Avg V-mag','Avg SE','Avg Air Mass','Observation Count','Uncertainty']]

    combined = read_all_csvs()
    combined = clean_bad_data(combined)
    combined = sorted(combined,key=getKey)
    combined_good_air = filter_by_air_mass(combined)
    print('#observations:',len(combined))
    print('#observations air mass <= 2.0:',len(combined_good_air))

    write_csv(header+combined,"bruce_gary_raw_data_combined")
    scatter_plot(combined,plot_name="scatter", plot_title="Bruce Gary Raw Data Unmodified",marker_size=1)
    write_csv(header+combined_good_air,"bruce_gary_raw_data_good_air_combined")
    scatter_plot(combined_good_air,plot_name="scatter_good_air", plot_title="Bruce Gary Raw Data Air Mass <= 2.0",marker_size=1)
    daily_bins = get_daily_binned_data(combined_good_air)
    write_csv(header2+daily_bins,"bruce_gary_raw_data_good_air_daily_bins")
    scatter_plot(daily_bins,plot_name="scatter_good_air_daily_bins", plot_title="Bruce Gary Daily Bins Air Mass <= 2.0",marker_size=16)
    print(daily_bins)

def get_daily_binned_data(data):
    days = []
    day_time_sum = []
    day_mag_sum = []
    day_se_sum = []
    day_airmass_sum = []
    day_counts = []
    for i in range(len(data)):
        x = int(float(data[i][0]))
        if(x not in days):
            days.append(x)
            day_time_sum.append(0)
            day_mag_sum.append(0)
            day_se_sum.append(0)
            day_airmass_sum.append(0)
            day_counts.append(0)
        index = days.index(x)
        day_time_sum[index] += float(data[i][0])
        day_mag_sum[index] += float(data[i][1])
        day_se_sum[index] += float(data[i][2])
        day_airmass_sum[index] += float(data[i][3])
        day_counts[index] += 1
    binned = []
    for i in range(len(days)):
        time_avg = day_time_sum[i]/day_counts[i]
        mag_avg = day_mag_sum[i]/day_counts[i]
        se_avg = day_se_sum[i]/day_counts[i]
        airmass_avg = day_airmass_sum[i]/day_counts[i]
        uncertainty = se_avg/math.sqrt(day_counts[i])
        binned.append([time_avg,mag_avg,se_avg,airmass_avg,day_counts[i],uncertainty])
    return binned

def filter_by_air_mass(data):
    filtered = []
    for i in range(1,len(data)):
        if(float(data[i][3]) <= 2.0):
            filtered.append(data[i])
    return filtered

def scatter_plot(combined, plot_name, plot_title, marker_size):
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
    # default scale is 1 in your original case, scales with other cases:

    plt.title(title)
    plt.gca().invert_yaxis()
    plt.gcf().autofmt_xdate()
    plt.grid()
    plt.scatter(x,y,s=marker_size,marker='o')
    plt.tight_layout()
    plt.rcParams["figure.figsize"] = (8,4.5)
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

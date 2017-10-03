import csv
import glob
import matplotlib.pyplot as plt
import math
import datetime
import numpy as np
import scipy.optimize as opt


start_mag = 11.905
dip_mjd_ranges = [[57875.,57908.],[57910.,57930.],[57969.,57988.],[58000.,58012.]]
elsie_range = [[57875.,57908.]]
celeste_range = [[57910.,57930.]]
skara_brae_range = [[57969.,57988.]]
angkor_range = [[58000.,58012.]]

def main():
    header = [['MJD','V-mag','SE','air mass']]
    header2 = [['Avg MJD','Avg V-mag','Avg SE','Avg Air Mass','Observation Count','Uncertainty']]

    vbandCombined = read_all_csvs('RawData/vband/*')
    vbandCombined = clean_bad_data(vbandCombined)
    vbandCombined = sorted(vbandCombined,key=getKey)
    vband_combined_good_air = filter_by_air_mass(vbandCombined)
    print('#observations vband:',len(vbandCombined))
    print('#observations air mass <= 2.0 vband:',len(vband_combined_good_air))

    gbandCombined = read_all_csvs("RawData/g'band/*")
    gbandCombined = clean_bad_data(gbandCombined)
    gbandCombined = sorted(gbandCombined,key=getKey)
    gband_combined_good_air = filter_by_air_mass(gbandCombined)
    print('#observations gband:',len(gbandCombined))
    print('#observations air mass <= 2.0 gband:',len(gband_combined_good_air))

    write_csv(header+vbandCombined,"vband_bruce_gary_raw_data_combined")
    scatter_plot(vbandCombined,plot_name="vband_scatter", plot_title="Bruce Gary V Band Raw Data Unmodified",marker_size=1)
    write_csv(header+vband_combined_good_air,"vband_bruce_gary_raw_data_good_air_combined")
    scatter_plot(vband_combined_good_air,plot_name="vband_scatter_good_air", plot_title="Bruce Gary V Band Raw Data Air Mass <= 2.0",marker_size=1)
    vband_daily_bins = get_bins(vband_combined_good_air,60*60*24)
    write_csv(header2+vband_daily_bins,"vband_bruce_gary_raw_data_good_air_daily_bins")
    scatter_plot(vband_daily_bins,plot_name="vband_scatter_good_air_daily_bins", plot_title="Bruce Gary V Band Daily Bins Air Mass <= 2.0",marker_size=16)
    vband_hourly_bins = get_bins(vband_combined_good_air,60*60)
    write_csv(header2+vband_hourly_bins,"vband_bruce_gary_raw_data_good_air_hourly_bins")
    scatter_plot(vband_hourly_bins,plot_name="vband_scatter_good_air_hourly_bins", plot_title="Bruce Gary V Band Hourly Bins Air Mass <= 2.0",marker_size=16)

    write_csv(header+gbandCombined,"g'band_bruce_gary_raw_data_combined")
    scatter_plot(vbandCombined,plot_name="g'band_scatter", plot_title="Bruce Gary g' Band Raw Data Unmodified",marker_size=1)
    write_csv(header+gband_combined_good_air,"g'band_bruce_gary_raw_data_good_air_combined")
    scatter_plot(gband_combined_good_air,plot_name="g'band_scatter_good_air", plot_title="Bruce Gary g' Band Raw Data Air Mass <= 2.0",marker_size=1)
    gband_daily_bins = get_bins(gband_combined_good_air,60*60*24)
    write_csv(header2+gband_daily_bins,"g'band_bruce_gary_raw_data_good_air_daily_bins")
    scatter_plot(gband_daily_bins,plot_name="g'band_scatter_good_air_daily_bins", plot_title="Bruce Gary Daily g' Band Bins Air Mass <= 2.0",marker_size=16)
    gband_hourly_bins = get_bins(gband_combined_good_air,60*60)
    write_csv(header2+gband_hourly_bins,"g'band_bruce_gary_raw_data_good_air_hourly_bins")
    scatter_plot(gband_hourly_bins,plot_name="g'band_scatter_good_air_hourly_bins", plot_title="Bruce Gary g' Band Hourly Bins Air Mass <= 2.0",marker_size=16)

    vband_elsie_hourly_bins = include_data(vband_hourly_bins,elsie_range)
    scatter_plot(vband_elsie_hourly_bins,plot_name="vband_scatter_good_air_hourly_bins_elsie", plot_title="Bruce Gary V Band Hourly Bins Elsie Air Mass <= 2.0",marker_size=16)

    vband_celeste_hourly_bins = include_data(vband_hourly_bins,celeste_range)
    scatter_plot(vband_celeste_hourly_bins,plot_name="vband_scatter_good_air_hourly_bins_celeste", plot_title="Bruce Gary V Band Hourly Bins Celeste Air Mass <= 2.0",marker_size=16)

    vband_skara_brae_hourly_bins = include_data(vband_hourly_bins,skara_brae_range)
    scatter_plot(vband_skara_brae_hourly_bins,plot_name="vband_scatter_good_air_hourly_bins_skara_brae", plot_title="Bruce Gary V Band Hourly Bins Skara Brae Air Mass <= 2.0",marker_size=16)

    angkor_hourly_bins = include_data(vband_hourly_bins,angkor_range)
    scatter_plot(angkor_hourly_bins,plot_name="vband_scatter_good_air_hourly_bins_angkor", plot_title="Bruce Gary V Band Hourly Bins Angkor Air Mass <= 2.0",marker_size=16)

    #print(daily_bins)
    dips_excluded = exclude_data(vband_combined_good_air, dip_mjd_ranges)
    daily_bins_dips_excluded = get_bins(dips_excluded,60*60*24)
    write_csv(header2+daily_bins_dips_excluded,"vband_bruce_gary_raw_data_good_air_daily_bins_dips_excluded")
    scatter_plot(daily_bins_dips_excluded,plot_name="vband_scatter_good_air_daily_bins_dips_excluded", plot_title="Bruce Gary V Band Daily Bins (Dips Excluded) Air Mass <= 2.0",marker_size=16,fit_type="gaussian")

def get_bins(data_list,bin_seconds):
    #data_list = filter_to_mjd(data_list,mjd)
    bin_increment = (1./(24.*60.*60.))*bin_seconds
    start = float(data_list[0][0])
    end = float(data_list[len(data_list)-1][0])
    if bin_seconds == 60*60*24:
        start = int(start)
        end = int(end)
        bin_increment = 1
    bin_count = int((start-end)/bin_increment)+1
    bins = []
    bin_time_sum = []
    bin_mag_sum = []
    bin_se_sum = []
    bin_airmass_sum = []
    bin_counts = []
    for i in range(len(data_list)):
        relative_mjd = float(data_list[i][0])-float(start)
        if bin_seconds == 60*60*24:
            bin = int(relative_mjd)
        else:
            bin = int(relative_mjd/bin_increment)
        if bin not in bins:
            bins.append(bin)
            bin_time_sum.append(0)
            bin_mag_sum.append(0)
            bin_se_sum.append(0)
            bin_airmass_sum.append(0)
            bin_counts.append(0)
        index = bins.index(bin)
        bin_time_sum[index] += float(data_list[i][0])
        bin_mag_sum[index] += float(data_list[i][1])
        bin_se_sum[index] += float(data_list[i][2])
        bin_airmass_sum[index] += float(data_list[i][3])
        bin_counts[index]+=1
    binned = []
    for i in range(len(bins)):
        time_avg = bin_time_sum[i]/bin_counts[i]
        mag_avg = bin_mag_sum[i]/bin_counts[i]
        se_avg = bin_se_sum[i]/bin_counts[i]
        airmass_avg = bin_airmass_sum[i]/bin_counts[i]
        uncertainty = se_avg/math.sqrt(bin_counts[i])
        binned.append([time_avg,mag_avg,se_avg,airmass_avg,bin_counts[i],uncertainty])
    return binned

def filter_by_air_mass(data):
    filtered = []
    for i in range(1,len(data)):
        if(float(data[i][3]) <= 2.0):
            filtered.append(data[i])
    return filtered

def lin_func(x, a):
    return start_mag+(a*x)

def gaussian_func(x, a, mu, sigma):
    return a * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

def scatter_plot(combined, plot_name, plot_title, marker_size, fit_type="none"):
    x = []
    y = []
    error = []
    for i in range(len(combined)):
        x.append(float(combined[i][0]))
        y.append(float(combined[i][1]))
        if len(combined[i]) > 5:
            error.append(combined[i][5])
        else:
            error.append(0)
    x = np.array(x).astype("float")
    y = np.array(y).astype("float")
    x_float = np.array(x).astype("float")
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
    if len(combined[0]) < 6:
        plt.scatter(x, y, s=marker_size, marker='o')
    else:
        plt.errorbar(x=x,y=y, yerr=error, fmt='o',markersize=4.0)
    plt.tight_layout()
    plt.rcParams["figure.figsize"] = (8, 4.5)

    #x_inc, y_inc = exclude_dips(x_float,y)
    # if fit_type=="linear":
    #     optimizedParameters, pcov = opt.curve_fit(lin_func, x_float-x_float[0], y)
    #     print(*optimizedParameters)
    #     plt.plot(x, lin_func(x_float-x_float[0], *optimizedParameters), label="linear fit")
    #
    # if fit_type=="gaussian":
    #     X = np.arange(len(x_float))
    #     mu_guess = np.sum(X * y) / np.sum(y)
    #     sigma_guess = np.sqrt(np.abs(np.sum((X - x_float) ** 2 * y) / np.sum(y)))
    #
    #     a_guess = y.max()
    #     optimizedParameters, pcov = opt.curve_fit(gaussian_func,x_float-x_float[0] ,y,p0=[a_guess, mu_guess, sigma_guess],maxfev=10000)
    #     print(*optimizedParameters)
    #     plt.plot(x, gaussian_func(x_float-x_float[0] , *optimizedParameters), label="gaussian fit")
    plt.savefig("BruceGarySavedPlots/" + plot_name + ".png")
    plt.show()

def exclude_data(data, excluded_ranges):
    data_inc = []
    for i in range(1,len(data)):
        exclude = 0
        for j in range(len(excluded_ranges)):
            if float(data[i][0]) >= excluded_ranges[j][0] and float(data[i][0]) <= excluded_ranges[j][1]:
                exclude = 1
        if exclude == 0:
            data_inc.append([data[i][0],data[i][1],data[i][2],data[i][3]])
    data_inc = np.array(data_inc)
    return data_inc

def include_data(data, included_ranges):
    data_inc = []
    for i in range(1,len(data)):
        include = 0
        for j in range(len(included_ranges)):
            if float(data[i][0]) >= included_ranges[j][0] and float(data[i][0]) <= included_ranges[j][1]:
                include = 1
        if include == 1:
            data_inc.append([data[i][0],data[i][1],data[i][2],data[i][3]])
    data_inc = np.array(data_inc)
    return data_inc

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

def read_all_csvs(path):
    files = glob.glob(path)
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

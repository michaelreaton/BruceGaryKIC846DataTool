import csv
import matplotlib.pyplot as plt
import math
import datetime
import numpy as np
import scipy.optimize as opt
dip_jd_ranges = [[2457891.,2457897.],[2457915.,2457932.],[2457964.,2457981.]]

def main():
    header = [['Avg JD','Avg V-mag','Avg UNC','Avg Air Mass','Observation Count','Uncertainty']]
    parsed = read_all_band_data()
    v_band = filter_by_band(parsed,"V")
    v_band = exclude_data(v_band,dip_jd_ranges)
    v_band_weekly_bins = get_bins(v_band,60*60*24*7,.001)
    scatter_plot(v_band_weekly_bins,plot_name="aavso_vband_1_week_bins_scatter", plot_title="AAVSO 1 weeks bins, air mass <= 2, uncertainty < .001",marker_size=16)
    write_csv(header+v_band_weekly_bins,"aavso_vband_data_good_air_weekly_bins_dips_excluded")

    #print(v_band_daily_bins)

def read_all_band_data():
    csv_file = open('AAVSORawData/aavsodata_20170824.txt', 'rt')
    reader = csv.reader(csv_file)
    parsed_data = []
    next(reader)
    for row in reader:
        if (is_float(row[0]) and is_float(row[1]) and is_float(row[2]) and is_float(row[12])):
            data_row = [float(row[0]),float(row[1]),float(row[2]), float(row[12]), row[4]]
            parsed_data.append(data_row)
            #print(data_row)
    return parsed_data

def filter_by_band(data, band):
    filtered_data = []
    for i in range(len(data)):
        if(data[i][4] == band and data[i][3] <= 2.0):
            filtered_data.append(data[i])
    filtered_data = sorted(filtered_data,key=getKey)
    return filtered_data

def get_bins(data_list,bin_seconds,min_uncertainty):
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
    bin_unc_sum = []
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
            bin_unc_sum.append(0)
            bin_airmass_sum.append(0)
            bin_counts.append(0)
        index = bins.index(bin)
        bin_time_sum[index] += float(data_list[i][0])
        bin_mag_sum[index] += float(data_list[i][1])
        bin_unc_sum[index] += float(data_list[i][2])
        bin_airmass_sum[index] += float(data_list[i][3])
        bin_counts[index]+=1
    binned = []
    for i in range(len(bins)):
        time_avg = bin_time_sum[i]/bin_counts[i]
        mag_avg = bin_mag_sum[i]/bin_counts[i]
        unc_avg = bin_unc_sum[i]/bin_counts[i]
        airmass_avg = bin_airmass_sum[i]/bin_counts[i]
        uncertainty = unc_avg/math.sqrt(bin_counts[i])
        if (uncertainty <= min_uncertainty):
            binned.append([time_avg,mag_avg,unc_avg,airmass_avg,bin_counts[i],uncertainty])
    return binned

def scatter_plot(combined, plot_name, plot_title, marker_size, fit_type="none"):
    x = []
    y = []
    x_doy_2016 = []
    error = []
    for i in range(len(combined)):
        x.append(float(combined[i][0]))
        x_doy_2016.append(jd_to_doy_2016(x[i]))
        y.append(float(combined[i][1]))
        if len(combined[i]) > 5:
            error.append(combined[i][5])
        else:
            error.append(0)

    x = np.array(x).astype("float")
    y = np.array(y).astype("float")
    x_float = np.array(x).astype("float")

    x_doy_2016 = np.array(x_doy_2016).astype("float")
    x = jddates_to_gregoriandates(x)
    plt.xlabel("Date")
    plt.ylabel("V-Mag")
    start_date = x[0]
    end_date = x[len(x)-1]
    title = plot_title + "\nCurve: 11.8400+ .0880*exp(-(doy_2016 - 1083.) ^ 2 / (617 ^ 2))\n(" + str(start_date) + " to " + str(end_date) + ")"
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
    if fit_type=="linear":
        optimizedParameters, pcov = opt.curve_fit(lin_func, x_float-x_float[0], y)
        print(*optimizedParameters)
        plt.plot(x, lin_func(x_float-x_float[0], *optimizedParameters), label="linear fit")

    if fit_type=="gaussian":
        X = np.arange(len(x_float))
        mu_guess = np.sum(X * y) / np.sum(y)
        sigma_guess = np.sqrt(np.abs(np.sum((X - x_float) ** 2 * y) / np.sum(y)))

        a_guess = y.max()
        optimizedParameters, pcov = opt.curve_fit(gaussian_func,x_float-x_float[0] ,y,p0=[a_guess, mu_guess, sigma_guess],maxfev=10000)
        print(*optimizedParameters)
    plt.plot(x, aavso_gaussian_func(x_doy_2016), label="Bruce Gary gaussian fit")
    plt.savefig("AAVSOSavedPlots/" + plot_name + ".png")
    plt.show()

def jddates_to_gregoriandates(jd):
    dates = [jd_to_gregorian(jd_val) for jd_val in jd]
    return dates

def jd_to_gregorian(jd_val):
    start_date = datetime.date(1858, 11, 17)
    return start_date + datetime.timedelta(jd_val - 2400000.5)

def jddates_to_doy_2016(jd):
    dates = [jd_to_doy_2016(jd_val) for jd_val in jd]
    return dates

def jd_to_doy_2016(jd_val):
    start_date = datetime.date(1858, 11, 17)
    gregorian_date = start_date + datetime.timedelta(jd_val - 2400000.5)
    year2016 = datetime.date(2016,1,1)
    return (gregorian_date - year2016).total_seconds()/float(86400)

def lin_func(x, a):
    return 11.80 +(a*x)

def bruce_gary_gaussian_func(x):
    return 11.8934+ .0334* np.exp(-(x - 1083.) ** 2 / (617 ** 2))

def aavso_gaussian_func(x):
    return 11.84+ .0880* np.exp(-(x - 1083.) ** 2 / (617 ** 2))

def gaussian_func(x, a, mu, sigma):
    return a * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

def getKey(item):
    return item[0]

def write_csv(list, name):
    with open("AAVSOCombined/" + name + ".csv", "w") as data_file:
        writer = csv.writer(data_file,lineterminator='\n')
        writer.writerows(list)

def exclude_data(data, excluded_ranges):
    data_inc = []
    for i in range(1,len(data)):
        exclude = 0
        for j in range(len(excluded_ranges)):
            if float(data[i][0]) >= excluded_ranges[j][0] and float(data[i][0]) <= excluded_ranges[j][1]:
                exclude = 1
        if exclude == 0:
            data_inc.append([data[i][0],data[i][1],data[i][2],data[i][3],data[i][4]])
    data_inc = np.array(data_inc)
    return data_inc

def is_float(value):
  try:
    float(value)
    return True
  except:
    return False

if __name__ == "__main__":
    main()
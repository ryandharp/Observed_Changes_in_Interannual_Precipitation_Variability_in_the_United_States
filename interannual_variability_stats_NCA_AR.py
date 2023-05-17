#%% Importing Modules
import pandas as pd
import numpy as np
import pymannkendall as mk
import time
import sys
import multiprocessing as mp
import pymannkendall_rdh as mk_rdh



#%% Collecting Input Values
py_script_name, begin_year, window_length = sys.argv
begin_year = int(begin_year)
window_length = int(window_length)

# begin_year = 1950  # first year for trend analysis
percent_qualifying_years = 90  # percentage of years qualifying necessary for an acceptable station record
yearly_pct_threshold = 10  # percent of missing values allowed
precip_min_cutoff = 1  # mm recording needed for a 'rainy' day

# window_length = 11  # number of years for smoothing
window_buffer = int((window_length-1)/2)
threshold = .75  # threshold for percentage of years needed to qualify for a smoothing
zero_year_threshold = 0  # number of zero precip years allowed


#%% Defining Functions
# streamlined mann-kendall test
def mannkendall(time_series):
    trend, h, p, z, Tau, s, var_s, slope, intercept = mk.original_test(time_series)
    return p, slope

# running mean
def window_calculations(time_series, window_length, year_list, threshold):
    half_window = window_length // 2
    windowed_data_mean = np.zeros(np.shape(year_list[half_window:-half_window])); windowed_data_mean[:] = np.nan
    windowed_data_std = np.zeros(np.shape(year_list[half_window:-half_window])); windowed_data_std[:] = np.nan
    windowed_data_cov = np.zeros(np.shape(year_list[half_window:-half_window])); windowed_data_cov[:] = np.nan
    c = 0
    for year in year_list[half_window:-half_window]:  # looping through years to use as window centers
        windowed_data = time_series[(time_series.index >= (year - half_window)) & (time_series.index <= (year + half_window))]
        if windowed_data.count() / window_length > threshold:
            windowed_data_mean[c] = windowed_data.mean()
            windowed_data_std[c] = windowed_data.std()
            windowed_data_cov[c] = windowed_data_std[c]/windowed_data_mean[c]
        c += 1
    return windowed_data_mean, windowed_data_std, windowed_data_cov

def window_percentiles(time_series, window_length, year_list):
    half_window = window_length // 2
    windowed_data_10 = np.zeros(np.shape(year_list[half_window:-half_window])); windowed_data_10[:] = np.nan
    windowed_data_25 = np.zeros(np.shape(year_list[half_window:-half_window])); windowed_data_25[:] = np.nan
    c = 0
    for year in year_list[half_window:-half_window]:  # looping through years to use as window centers
        windowed_data = time_series[(year_list >= (year - half_window)) & (year_list <= (year + half_window))]
        windowed_data_mean[c] = windowed_data.mean()
        windowed_data_std[c] = windowed_data.std()
        c += 1
    return windowed_data_mean, windowed_data_std

# data threshold analysis
def load_station(fname):
    df = pd.read_csv(fname, dtype={'DATE': object, 'PRCP': float, 'TMAX': float, 'TMIN': float, 'PRCP_ATTRIBUTES': str,
                                   'TMAX_ATTRIBUTES': str, 'TMIN_ATTRIBUTES': str}, low_memory=False)
    data = df.filter(['DATE', 'PRCP', 'TMAX', 'TMIN', 'PRCP_ATTRIBUTES', 'TMAX_ATTRIBUTES', 'TMIN_ATTRIBUTES'])
    data['DATE'] = pd.to_datetime(data['DATE'])
    del df
    return data

def explode_flags(compact_data):
    min_ind = compact_data.index[compact_data['PRCP_ATTRIBUTES'].notnull()].min()
    if compact_data['PRCP_ATTRIBUTES'][min_ind].count(',') == 3:
        compact_data[['PRCP_MEAS_FLAG', 'PRCP_QUAL_FLAG', 'PRCP_SOURCE_CODE', 'PRCP_TIME']] = compact_data['PRCP_ATTRIBUTES'].str.split(',', expand=True)
    elif compact_data['PRCP_ATTRIBUTES'][min_ind].count(',') == 2:
        compact_data[['PRCP_MEAS_FLAG', 'PRCP_QUAL_FLAG', 'PRCP_SOURCE_CODE']] = compact_data['PRCP_ATTRIBUTES'].str.split(',', expand=True)
    flag_mask = (compact_data['PRCP_QUAL_FLAG'].isna()) | (compact_data['PRCP_QUAL_FLAG'].isin(qual_flags))
    compact_data['PRCP'] = compact_data['PRCP'].where(~flag_mask, np.nan)
    return compact_data

def get_days_in_year2(years_array, first_year, last_year):
    days_in_year_array = np.zeros(np.shape(years_array))
    i = 0
    for yr in years_array:
        if yr % 4 == 0:
            days_in_year_array[i] = 366
        else:
            days_in_year_array[i] = 365
        i += 1
    total_years_array = np.arange(first_year, last_year+1)
    total_days_array = np.zeros(np.shape(total_years_array))
    i = 0
    for yr in total_years_array:
        if yr % 4 == 0:
            total_days_array[i] = 366
        else:
            total_days_array[i] = 365
        i += 1
    total_days = sum(total_days_array)
    total_days_series = pd.Series(total_days_array, index=total_years_array)
    return days_in_year_array, total_days, total_days_series, total_years_array

def get_pct_missing(station_data, days_in_year, total_days, years_array, total_days_series, start_year, last_year):
    date_mask = (station_data['DATE'].dt.year >= start_year) & (station_data['DATE'].dt.year <= last_year)
    num_days_in_data_year = station_data['PRCP'].groupby(station_data['DATE'][date_mask].dt.year).count()
    num_days_in_data_year_series = pd.Series(num_days_in_data_year, years_array)
    # num_missing_flagged = station_data['PRCP'].isnull().groupby(station_data['DATE'][date_mask].dt.year).sum().astype(int)
    # years_mask = np.isin(years, num_days_in_data_year.index)
    yearly_pct_missing  = 100 - (num_days_in_data_year_series.divide(total_days_series, fill_value=0)*100)
    total_pct_missing = 100 - sum(num_days_in_data_year)/total_days*100
    # yearly_pct_missing = 100 - (days_in_year[years_mask] - num_missing_flagged)/days_in_year[years_mask]*100
    # total_pct_missing = 100 - (days_in_year.sum() - num_missing_flagged.sum())/days_in_year.sum()*100
    return yearly_pct_missing, total_pct_missing

def region_mannkendall(regional_df):
    regional_ts = np.transpose(np.array(regional_df))
    regional_mk_results = mk.regional_test(np.transpose(regional_ts[1:]))
    region_var_p = regional_mk_results[2]
    region_var_slope = regional_mk_results[7]
    return region_var_p, region_var_slope

def region_mannkendall_rdh(regional_df):
    x = np.array(regional_df)
    x_old = x[window_buffer:-window_buffer, :]
    x_old = x_old[:, np.nansum(x_old, axis=0) > 0]
    regional_mk_results = mk_rdh.regional_test(x_old)
    region_var_p = regional_mk_results[2]
    region_var_slope = regional_mk_results[7]
    return region_var_p, region_var_slope



#%% Loading File List
ghcn_stations = pd.read_csv('/home/rdh0715/US_stations_with_NCA.csv')

qual_stations = ghcn_stations[(ghcn_stations['len_years'] >= ((2020-begin_year+1) * percent_qualifying_years / 100))]
# qual_stations = qual_stations[qual_stations['NCA_region'] == domain]

# qual_stations = qual_stations[qual_stations['station_id'].str.slice(0,2) == 'US']
qual_station_id = qual_stations['station_id']
qual_lat = qual_stations['latitude']
qual_lon = qual_stations['longitude']
qual_year_length = qual_stations['num_qual_years']
nca_region = qual_stations['NCA_region']


#%% Establishing Arrays
array_size = np.shape(qual_station_id)
station_id_list = []
station_lat = np.zeros(array_size); station_lat[:] = np.nan
station_lon = np.zeros(array_size); station_lon[:] = np.nan
start_year_array = np.zeros(array_size); start_year_array[:] = np.nan
last_year_array = np.zeros(array_size); last_year_array[:] = np.nan
station_qual_year_length = np.zeros(array_size); station_qual_year_length[:] = np.nan

smoothed_annual_precip_missing_percentage_array = np.zeros(array_size); smoothed_annual_precip_missing_percentage_array[:] = np.nan
smoothed_annual_precip_freq_missing_percentage_array = np.zeros(array_size); smoothed_annual_precip_freq_missing_percentage_array[:] = np.nan
smoothed_mean_dry_missing_percentage_array = np.zeros(array_size); smoothed_mean_dry_missing_percentage_array[:] = np.nan

del array_size


#%%
# preparing for analysis
qual_flags = ['D', 'G', 'I', 'K', 'L', 'M', 'N', 'O', 'R', 'S', 'T', 'W', 'X', 'Z']

len_qual_stations = len(qual_station_id)
print(len_qual_stations)

num_regions = len(np.unique(nca_region)) + 1
region_annual_mean_p = np.zeros(num_regions); region_annual_mean_p[:] = np.nan
region_annual_mean_slope = np.zeros(num_regions); region_annual_mean_slope[:] = np.nan
region_annual_mean = np.zeros(num_regions); region_annual_mean[:] = np.nan
region_annual_var_p = np.zeros(num_regions); region_annual_var_p[:] = np.nan
region_annual_var_slope = np.zeros(num_regions); region_annual_var_slope[:] = np.nan
region_annual_var = np.zeros(num_regions); region_annual_var[:] = np.nan
region_annual_cov_p = np.zeros(num_regions); region_annual_cov_p[:] = np.nan
region_annual_cov_slope = np.zeros(num_regions); region_annual_cov_slope[:] = np.nan
region_annual_cov = np.zeros(num_regions); region_annual_cov[:] = np.nan
region_annual_freq_p = np.zeros(num_regions); region_annual_freq_p[:] = np.nan
region_annual_freq_slope = np.zeros(num_regions); region_annual_freq_slope[:] = np.nan
region_annual_freq = np.zeros(num_regions); region_annual_freq[:] = np.nan
domain_num_stations = np.zeros(num_regions); domain_num_stations[:] = np.nan


start = time.time()
def wrapper(domain):
# for domain in np.unique(nca_region):
    # domain = 5  # using Great Lakes as an example

    domain_stations = qual_stations[nca_region == domain]

    domain = int(domain)
    print(domain)
    print(domain_stations)
    print(time.time()-start)

    # domain_year_array = np.arange(qual_stations['start_year'].min(), qual_stations['last_year'].max(), 1)
    domain_year_array = np.arange(begin_year, qual_stations['last_year'].max(), 1)

    domain_annual_mean = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_mean[:] = np.nan
    domain_annual_precip_var = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_precip_var[:] = np.nan
    domain_annual_precip_cov = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_precip_cov[:] = np.nan
    domain_annual_freq = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_freq[:] = np.nan
    domain_annual_mean_df = pd.DataFrame(np.transpose(domain_annual_mean)); domain_annual_mean_df.index = domain_year_array
    domain_annual_precip_var_df = pd.DataFrame(np.transpose(domain_annual_precip_var)); domain_annual_precip_var_df.index = domain_year_array
    domain_annual_precip_cov_df = pd.DataFrame(np.transpose(domain_annual_precip_cov)); domain_annual_precip_cov_df.index = domain_year_array
    domain_annual_freq_df = pd.DataFrame(np.transpose(domain_annual_freq)); domain_annual_freq_df.index = domain_year_array

    j = 0
    num_domain_stations = 0
    for i, station in domain_stations['station_id'].iteritems():
        fname = station + '.csv'
        print(fname)

        station_id_list.append(station)
        station_lat[j] = qual_lat.loc[i]
        station_lon[j] = qual_lon.loc[i]
        station_qual_year_length[j] = qual_year_length.loc[i]


        #%% Loading Data
        # fname = 'USW00094728.csv'
        # df = pd.read_csv('~/Documents/great_lakes_precip_variability/raw_data/'+fname)
        df = pd.read_csv('/projects/b1045/GHCN/'+fname, low_memory=False)

        # temp_data = load_station('~/Documents/great_lakes_precip_variability/raw_data/'+fname)
        temp_data = load_station('/projects/b1045/GHCN/'+fname)
        station_data = explode_flags(temp_data)

        start_date = station_data['DATE'].min()
        start_year = start_date.date().year
        start_year_array[j] = start_year
        last_date = station_data['DATE'].max()
        last_year = last_date.date().year
        last_year_array[j] = last_year
        years_array = station_data['DATE'].dt.year.unique()
        # years = np.arange(start_year, last_year)
        years = np.arange(begin_year, last_year)

        days_in_year, total_days, total_days_series, total_years_array = get_days_in_year2(years_array, start_year, last_year)
        yearly_pct_missing, total_pct_missing = get_pct_missing(station_data, days_in_year, total_days, years_array, total_days_series, start_year, last_year)
        qual_years = total_years_array[yearly_pct_missing <= yearly_pct_threshold]
        qual_years = qual_years[qual_years >= begin_year]
        pct_qual_years = len(qual_years) / (2020-begin_year+1) * 100
        if pct_qual_years < percent_qualifying_years:
            j += 1
            continue

        data = df.filter(['DATE', 'PRCP', 'TMAX', 'TMIN'])
        data['DATE'] = pd.to_datetime(data['DATE'])

        del df, temp_data, station_data, start_date, last_date, years_array
        del days_in_year, total_days, total_days_series, total_years_array, yearly_pct_missing, total_pct_missing


        #%% Data Cleaning
        # converting mm and removing all events less than 1/3 mm
        data['PRCP'] = data['PRCP']/10  # converting to mm from tenths-mm
        data['PRCP_DATE'] = data['PRCP'] >= precip_min_cutoff  # could use 3 mm cutoff instead
        data['PRCP_QUAL'] = data['PRCP'][data['PRCP_DATE']]
        # year_list = data['DATE'].dt.year.unique()
        full_precip = data['PRCP_QUAL'][data['DATE'].dt.year.isin(qual_years)]
        zero_years = qual_years[full_precip.groupby(data['DATE'].dt.year).sum() == 0]  # testing for weird "valid" records with complete years of zero precipitation
        if len(zero_years) > zero_year_threshold:
            j += 1
            continue
        del zero_years


        #%% Starting with Annual Statistics
        # annual precipitation
        year_bool = data['DATE'].dt.year.isin(qual_years)  # boolean of qualifying years

        annual_precip_sum = data['PRCP_QUAL'][year_bool].groupby(data['DATE'].dt.year).sum()
        domain_annual_mean_df[j] = annual_precip_sum
        annual_precip_smoothed, annual_precip_smoothed_var, annual_precip_smoothed_cov = window_calculations(annual_precip_sum, window_length, years, threshold)
        annual_precip_smoothed_var = pd.Series(annual_precip_smoothed_var)
        annual_precip_smoothed_var.index = years[window_buffer:-window_buffer]
        domain_annual_precip_var_df[j] = annual_precip_smoothed_var
        annual_precip_smoothed_cov = pd.Series(annual_precip_smoothed_cov)
        annual_precip_smoothed_cov.index = years[window_buffer:-window_buffer]
        domain_annual_precip_cov_df[j] = annual_precip_smoothed_cov
        smoothed_annual_precip_missing_percentage_array[j] = np.sum(pd.isna(annual_precip_smoothed))/len(annual_precip_smoothed)*100

        del annual_precip_smoothed, annual_precip_smoothed_var, annual_precip_smoothed_cov


        #%% Precipitation Frequency
        annual_precip_freq = data['PRCP_QUAL'][year_bool].groupby(data['DATE'].dt.year).count()
        domain_annual_freq_df[j] = annual_precip_freq

        del annual_precip_freq

        j += 1
        num_domain_stations += 1


#%% Printing Test Results

    # annual mean
    region_annual_mean_p[domain-1], region_annual_mean_slope[domain-1] = region_mannkendall_rdh(domain_annual_mean_df)
    region_annual_mean[domain-1] = domain_annual_mean_df.mean(axis=0).mean()

    # interannual variability
    region_annual_var_p[domain-1], region_annual_var_slope[domain-1] = region_mannkendall_rdh(domain_annual_precip_var_df)
    region_annual_var[domain-1] = domain_annual_precip_var_df.mean(axis=0).mean()

    # relative interannual variability
    region_annual_cov_p[domain-1], region_annual_cov_slope[domain-1] = region_mannkendall_rdh(domain_annual_precip_cov_df)
    region_annual_cov[domain-1] = domain_annual_precip_cov_df.mean(axis=0).mean()

    # annual frequency
    region_annual_freq_p[domain-1], region_annual_freq_slope[domain-1] = region_mannkendall_rdh(domain_annual_freq_df)
    region_annual_freq[domain-1] = domain_annual_freq_df.mean(axis=0).mean()
    domain_num_stations[domain-1] = num_domain_stations

    results_df = pd.DataFrame(columns = ['annual_mean_p', 'annual_mean_slope', 'annual_mean', 'annual_var_p',
                                         'annual_var_slope', 'annual_var', 'annual_cov_p', 'annual_cov_slope', 'annual_cov',
                                         'annual_freq_p', 'annual_freq_slope', 'annual_freq', 'num_domain_stations'])

    results_df['annual_mean_p'] = region_annual_mean_p
    results_df['annual_mean_slope'] = region_annual_mean_slope
    results_df['annual_mean'] = region_annual_mean
    results_df['annual_var_p'] = region_annual_var_p
    results_df['annual_var_slope'] = region_annual_var_slope
    results_df['annual_var'] = region_annual_var
    results_df['annual_cov_p'] = region_annual_cov_p
    results_df['annual_cov_slope'] = region_annual_cov_slope
    results_df['annual_cov'] = region_annual_cov
    results_df['annual_freq_p'] = region_annual_freq_p
    results_df['annual_freq_slope'] = region_annual_freq_slope
    results_df['annual_freq'] = region_annual_freq
    results_df['num_domain_stations'] = domain_num_stations

    results_df.to_csv('/home/rdh0715/interannual_variability_results/AR_results/NCA_domain_precip_stats_start_year_' + str(begin_year) + '_window_length_' + str(window_length) + '_region_' + str(domain) + '_AR2.csv')

if __name__ == '__main__':
    pool = mp.Pool(processes=11)
    pool.map(wrapper, np.array([1, 3, 4, 5, 6, 7, 8, 9, 10]))


#%% Old Code
# region_per_10_p = np.zeros(num_regions); region_per_10_p[:] = np.nan
# region_per_10_slope = np.zeros(num_regions); region_per_10_slope[:] = np.nan
# region_per_25_p = np.zeros(num_regions); region_per_25_p[:] = np.nan
# region_per_25_slope = np.zeros(num_regions); region_per_25_slope[:] = np.nan
# region_per_50_p = np.zeros(num_regions); region_per_50_p[:] = np.nan
# region_per_50_slope = np.zeros(num_regions); region_per_50_slope[:] = np.nan
# region_per_75_p = np.zeros(num_regions); region_per_75_p[:] = np.nan
# region_per_75_slope = np.zeros(num_regions); region_per_75_slope[:] = np.nan
# region_per_90_p = np.zeros(num_regions); region_per_90_p[:] = np.nan
# region_per_90_slope = np.zeros(num_regions); region_per_90_slope[:] = np.nan
# region_per_95_p = np.zeros(num_regions); region_per_95_p[:] = np.nan
# region_per_95_slope = np.zeros(num_regions); region_per_95_slope[:] = np.nan
# region_per_99_p = np.zeros(num_regions); region_per_99_p[:] = np.nan
# region_per_99_slope = np.zeros(num_regions); region_per_99_slope[:] = np.nan
# region_per_99p5_p = np.zeros(num_regions); region_per_99p5_p[:] = np.nan
# region_per_99p5_slope = np.zeros(num_regions); region_per_99p5_slope[:] = np.nan


    # domain_annual_10_per_exceedances = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_10_per_exceedances[:] = np.nan
    # domain_annual_25_per_exceedances = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_25_per_exceedances[:] = np.nan
    # domain_annual_50_per_exceedances = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_50_per_exceedances[:] = np.nan
    # domain_annual_75_per_exceedances = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_75_per_exceedances[:] = np.nan
    # domain_annual_90_per_exceedances = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_90_per_exceedances[:] = np.nan
    # domain_annual_95_per_exceedances = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_95_per_exceedances[:] = np.nan
    # domain_annual_99_per_exceedances = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_99_per_exceedances[:] = np.nan
    # domain_annual_99p5_per_exceedances = np.zeros([np.shape(domain_stations)[0], len(domain_year_array)]); domain_annual_99p5_per_exceedances[:] = np.nan
    # domain_annual_10_per_exceedances_df = pd.DataFrame(np.transpose(domain_annual_10_per_exceedances)); domain_annual_10_per_exceedances_df.index = domain_year_array
    # domain_annual_25_per_exceedances_df = pd.DataFrame(np.transpose(domain_annual_25_per_exceedances)); domain_annual_25_per_exceedances_df.index = domain_year_array
    # domain_annual_50_per_exceedances_df = pd.DataFrame(np.transpose(domain_annual_50_per_exceedances)); domain_annual_50_per_exceedances_df.index = domain_year_array
    # domain_annual_75_per_exceedances_df = pd.DataFrame(np.transpose(domain_annual_75_per_exceedances)); domain_annual_75_per_exceedances_df.index = domain_year_array
    # domain_annual_90_per_exceedances_df = pd.DataFrame(np.transpose(domain_annual_90_per_exceedances)); domain_annual_90_per_exceedances_df.index = domain_year_array
    # domain_annual_95_per_exceedances_df = pd.DataFrame(np.transpose(domain_annual_95_per_exceedances)); domain_annual_95_per_exceedances_df.index = domain_year_array
    # domain_annual_99_per_exceedances_df = pd.DataFrame(np.transpose(domain_annual_99_per_exceedances)); domain_annual_99_per_exceedances_df.index = domain_year_array
    # domain_annual_99p5_per_exceedances_df = pd.DataFrame(np.transpose(domain_annual_99p5_per_exceedances)); domain_annual_99p5_per_exceedances_df.index = domain_year_array


# %% Exceedence of percentile thresholds
# if (start_year <= first_half_year_min) & (last_year >= first_half_year_max):
#     percentile_threshold_10 = first_half_precip.quantile(q=0.10, interpolation='linear')
#     domain_annual_10_per_exceedances_df[j] = 1000*(full_precip > percentile_threshold_10).groupby(data['DATE'].dt.year).sum()
#
#     percentile_threshold_25 = first_half_precip.quantile(q=0.25, interpolation='linear')
#     domain_annual_25_per_exceedances_df[j] = 1000*(full_precip > percentile_threshold_25).groupby(data['DATE'].dt.year).sum()
#
#     percentile_threshold_50 = first_half_precip.quantile(q=0.50, interpolation='linear')
#     domain_annual_50_per_exceedances_df[j] = 1000*(full_precip > percentile_threshold_50).groupby(data['DATE'].dt.year).sum()
#
#     percentile_threshold_75 = first_half_precip.quantile(q=0.75, interpolation='linear')
#     domain_annual_75_per_exceedances_df[j] = 1000*(full_precip > percentile_threshold_75).groupby(data['DATE'].dt.year).sum()
#
#     percentile_threshold_90 = first_half_precip.quantile(q=0.90, interpolation='linear')
#     domain_annual_90_per_exceedances_df[j] = 1000*(full_precip > percentile_threshold_90).groupby(data['DATE'].dt.year).sum()
#
#     percentile_threshold_95 = first_half_precip.quantile(q=0.95, interpolation='linear')
#     domain_annual_95_per_exceedances_df[j] = 1000*(full_precip > percentile_threshold_95).groupby(data['DATE'].dt.year).sum()
#
#     percentile_threshold_99 = first_half_precip.quantile(q=0.99, interpolation='linear')
#     domain_annual_99_per_exceedances_df[j] = 1000*(full_precip > percentile_threshold_99).groupby(data['DATE'].dt.year).sum()
#
#     percentile_threshold_99p5 = first_half_precip.quantile(q=0.995, interpolation='linear')
#     domain_annual_99p5_per_exceedances_df[j] = 1000*(full_precip > percentile_threshold_99p5).groupby(data['DATE'].dt.year).sum()
#
#     del percentile_threshold_10, percentile_threshold_25, percentile_threshold_50, percentile_threshold_75
#     del percentile_threshold_90, percentile_threshold_95, percentile_threshold_99, percentile_threshold_99p5

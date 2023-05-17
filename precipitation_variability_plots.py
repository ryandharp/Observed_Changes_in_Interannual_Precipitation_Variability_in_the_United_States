import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import scipy.stats as ss
from scipy.stats import skewnorm




#%% Figure: Region Choropleth Maps with GHCN Station Scatter

NEON_domains = gpd.read_file('/Users/ryanharp/Documents/great_lakes_precip_variability/supplemental_files/NEONDomains_0/NEON_Domains.shp')
NEON_domains = NEON_domains.drop([1, 5], axis = 0)  # dropping 'extra' domain portions for centroid graphing purposes
NEON_domains.index = np.arange(1, 21, 1)
NEON_domains = NEON_domains.drop([18, 19, 20], axis = 0)

# domain_trends = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/NEON_domain_results/NEON_domain_precip_stats_v2.csv')
domain_trends = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/AR_results2/merged/NEON_domain_precip_stats_start_year_1950_window_length_11_merged.csv')
domain_trends = domain_trends.drop(labels = 'Unnamed: 0', axis = 1)  # dropping repeat index column
domain_trends.index = np.arange(1, 21, 1)  # setting index to be equal to domain ID
domain_trends = domain_trends.drop([18, 19, 20], axis = 0)

# results = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/precip_stats_50_years_1mm_threshold_v2.csv')
# results = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/smoothed_precip_stats_50_years_1mm_threshold.csv')
results = domain_trends

sig_value = 0.05
var = 'annual_mean'
cmap = 'bwr_r'  # 'PRGn_r'
results_sig = results[var + '_p']
sig_results = results[results[var + '_p'] < sig_value]
sig_results_p = sig_results[var + '_p']
sig_results_slope = sig_results[var + '_slope']
non_sig_results = results[results[var + '_p'] >= sig_value]
non_sig_results_p = non_sig_results[var + '_p']
non_sig_results_slope = non_sig_results[var + '_slope']

crs_new = ccrs.PlateCarree()
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')

test = domain_trends[var + '_slope']*10
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(14, 8))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray')
c_low = np.nanquantile(test, q=0.05)
c_high = np.nanquantile(test, q=0.95)
c_val = np.max([abs(c_low), abs(c_high)])
NEON_domains.plot(column = domain_trends[var + '_slope']*10, ax = ax, cmap=cmap, edgecolor='k', vmin=-c_val, vmax=c_val)
domain_trends[var + '_slope'][domain_trends[var + '_p'] > sig_value] = np.nan
domain_trends[var + '_slope'][domain_trends[var + '_slope'] == 0] = np.nan
NEON_domains.plot(column = domain_trends[var + '_slope']*10, ax = ax, cmap=cmap, edgecolor='k', vmin=-c_val, vmax=c_val, missing_kwds={'color':'none', 'edgecolor':'grey', 'hatch':'/'})
domain_trends[var + '_slope'] = np.nan
NEON_domains.plot(column = domain_trends[var + '_slope'], ax = ax, cmap=cmap, edgecolor='k', vmin=-c_val, vmax=c_val, missing_kwds={'color':'none', 'edgecolor':'k'})

plt.set_cmap(cmap)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-c_val, c_val))
cbar = fig.colorbar(sm, fraction=0.0275, pad=0.04)
cbar.set_label('Change in Annual Mean (mm/decade)', fontsize=20)
cbar.ax.tick_params(labelsize=16)
# cbar.set_label('Change in Relative Variability ($\mathregular{decade^{ -1}}$)', fontsize=20)
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=18)
plt.yticks(np.arange(30, 60, 10), fontsize=18)
plt.title('Trend in Annual Mean Precipitation', fontsize=24)
plt.savefig('../figures/interannual_variability_revised_plots/figs5a.svg')
plt.show()




#%% Same as above but for NCA regions

NCA_domains = gpd.read_file('/Users/ryanharp/Documents/great_lakes_precip_variability/supplemental_files/NCARegions/cb_2020_us_state_500k_ncaregions.shp')

NCA_domains = NCA_domains.drop([1], axis = 0)  # dropping 'extra' domain portions for centroid graphing purposes
NCA_domains.index = [1, 3, 4, 5, 6, 7, 8, 9, 10]
NCA_domains = NCA_domains.drop([1, 3], axis = 0)

domain_trends = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/AR_results2/merged/NCA_domain_precip_stats_start_year_1950_window_length_11_merged.csv')
domain_trends = domain_trends.drop(labels = 'Unnamed: 0', axis = 1)  # dropping repeat index column
domain_trends.index = [1, 3, 4, 5, 6, 7, 8, 9, 10]  # setting index to be equal to domain ID
domain_trends = domain_trends.drop([1, 3], axis = 0)
results = domain_trends

sig_value = 0.05
var = 'annual_freq'
cmap =  'bwr_r'  # 'PRGn_r'
results_sig = results[var + '_p']
sig_results = results[results[var + '_p'] < sig_value]
sig_results_p = sig_results[var + '_p']
sig_results_slope = sig_results[var + '_slope']
non_sig_results = results[results[var + '_p'] >= sig_value]
non_sig_results_p = non_sig_results[var + '_p']
non_sig_results_slope = non_sig_results[var + '_slope']

crs_new = ccrs.PlateCarree()
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')

test = domain_trends[var + '_slope']*10
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(14, 8))
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(states_provinces, edgecolor='gray')
c_low = np.nanquantile(test, q=0.05)
c_high = np.nanquantile(test, q=0.95)
c_val = np.max([abs(c_low), abs(c_high)])
NCA_domains.plot(column = domain_trends[var + '_slope']*10, ax = ax, cmap=cmap, edgecolor='k', vmin=-c_val, vmax=c_val)
domain_trends[var + '_slope'][domain_trends[var + '_p'] > sig_value] = np.nan
domain_trends[var + '_slope'][domain_trends[var + '_slope'] == 0] = np.nan
NCA_domains.plot(column = domain_trends[var + '_slope']*10, ax = ax, cmap=cmap, edgecolor='k', vmin=-c_val, vmax=c_val, missing_kwds={'color':'none', 'edgecolor':'grey', 'hatch':'/'})
domain_trends[var + '_slope'] = np.nan
NCA_domains.plot(column = domain_trends[var + '_slope'], ax = ax, cmap=cmap, edgecolor='k', vmin=-c_val, vmax=c_val, missing_kwds={'color':'none', 'edgecolor':'k'})

plt.set_cmap(cmap)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-c_val, c_val))
cbar = fig.colorbar(sm, fraction=0.0275, pad=0.04)
cbar.set_label('Change in Wet Day Frequency (days/decade)', fontsize=20)
cbar.ax.tick_params(labelsize=16)
# cbar.set_label('Change in Relative Variability ($\mathregular{decade^{ -1}}$)', fontsize=20)
plt.xlim([-129, -62])
plt.ylim([23, 52])
plt.xticks(np.arange(-125, -55, 10), fontsize=18)
plt.yticks(np.arange(30, 60, 10), fontsize=18)
plt.title('Trend in Annual Wet Day Frequency', fontsize=24)
plt.savefig('../figures/interannual_variability_revised_plots/figs7b.svg')
plt.show()





#%% Figure: Scatter plot with all info on it
sig_value = 0.05
# domain_trends = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/NEON_domain_results/NEON_domain_precip_stats_v2.csv')
domain_trends = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/AR_results2/merged/NEON_domain_precip_stats_start_year_1970_window_length_11_merged.csv')

domain_trends = domain_trends.drop(labels = 'Unnamed: 0', axis = 1)  # dropping repeat index column
domain_trends.index = np.arange(1, 21, 1)  # setting index to be equal to domain ID
domain_trends = domain_trends.drop([18, 19, 20], axis = 0)
domain_trends_sig = domain_trends.loc[:, ['annual_mean_p', 'annual_freq_p', 'annual_var_p', 'annual_cov_p']]
domain_trends_sig = domain_trends_sig < sig_value
domain_trends_slopes = domain_trends.loc[:, ['annual_mean_slope', 'annual_freq_slope', 'annual_var_slope', 'annual_cov_slope']]
domain_trends_means = domain_trends.loc[:, ['annual_mean', 'annual_freq', 'annual_var', 'annual_cov']]
domain_trends_norm = pd.DataFrame(np.array(domain_trends_slopes)/np.array(domain_trends_means))*1000  # *1000 to transform into %/decade
domain_trends_norm.columns = ['annual_mean', 'annual_freq', 'annual_var', 'annual_cov']
domain_trends_norm.index = np.arange(1, 18, 1)  # setting index to be equal to domain ID

fig, ax = plt.subplots(figsize=(8, 14))
plt.grid()
plt.axvline(color='k', linewidth=1)
# for i in np.arange(1, 18):
#     plt.axhline(y=i, color='grey', linewidth=0.5)
plt.scatter(domain_trends_norm['annual_mean'], domain_trends_norm['annual_mean'].index, facecolor='#0077BB', alpha=0.85, s=650)
plt.scatter(domain_trends_norm['annual_mean'][~domain_trends_sig['annual_mean_p']], domain_trends_norm['annual_mean'][~domain_trends_sig['annual_mean_p']].index, facecolor='w', alpha=0.85, s=150)
plt.scatter(domain_trends_norm['annual_freq'], domain_trends_norm['annual_freq'].index, facecolor='#88CCEE', alpha=0.85, s=650)
plt.scatter(domain_trends_norm['annual_freq'][~domain_trends_sig['annual_freq_p']], domain_trends_norm['annual_freq'][~domain_trends_sig['annual_freq_p']].index, facecolor='w', alpha=0.85, s=150)
plt.scatter(domain_trends_norm['annual_freq'][domain_trends_norm['annual_freq']==0], domain_trends_norm['annual_freq'][domain_trends_norm['annual_freq']==0].index, facecolor='w', alpha=0.85, s=150)
plt.scatter(domain_trends_norm['annual_var'], domain_trends_norm['annual_var'].index, facecolor='#117733', alpha=0.85, s=650)
plt.scatter(domain_trends_norm['annual_var'][~domain_trends_sig['annual_var_p']], domain_trends_norm['annual_var'][~domain_trends_sig['annual_var_p']].index, facecolor='w', alpha=0.85, s=150)
plt.scatter(domain_trends_norm['annual_cov'], domain_trends_norm['annual_cov'].index, facecolor='#44BB99', alpha=0.85, s=650)
plt.scatter(domain_trends_norm['annual_cov'][~domain_trends_sig['annual_cov_p']], domain_trends_norm['annual_cov'][~domain_trends_sig['annual_cov_p']].index, facecolor='w', alpha=0.85, s=150)
plt.xlim([-5, 5])
plt.xticks([-4, -2, 0, 2, 4], fontsize=16)
plt.xlabel('Trend (%/decade)', fontsize=20)
plt.yticks(np.arange(1,18), [])
ax.invert_yaxis()
plt.title('Change in Precipitation Metrics', fontsize=24)
plt.savefig('../figures/interannual_variability_revised_plots/fig2.svg')
plt.show()




#%% Same as above but for NCA
sig_value = 0.05

domain_trends = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/AR_results2/merged/NCA_domain_precip_stats_start_year_1970_window_length_11_merged.csv')
domain_trends = domain_trends.drop(labels = 'Unnamed: 0', axis = 1)  # dropping repeat index column
domain_trends.index = [1, 3, 4, 5, 6, 7, 8, 9, 10]  # setting index to be equal to domain ID
domain_trends = domain_trends.drop([1, 3], axis = 0)

domain_trends_sig = domain_trends.loc[:, ['annual_mean_p', 'annual_freq_p', 'annual_var_p', 'annual_cov_p']]
domain_trends_sig = domain_trends_sig < sig_value
domain_trends_slopes = domain_trends.loc[:, ['annual_mean_slope', 'annual_freq_slope', 'annual_var_slope', 'annual_cov_slope']]
domain_trends_means = domain_trends.loc[:, ['annual_mean', 'annual_freq', 'annual_var', 'annual_cov']]
domain_trends_norm = pd.DataFrame(np.array(domain_trends_slopes)/np.array(domain_trends_means))*1000  # *1000 to transform into %/decade
domain_trends_norm.columns = ['annual_mean', 'annual_freq', 'annual_var', 'annual_cov']
domain_trends_norm.index = [4, 5, 6, 7, 8, 9, 10]  # setting index to be equal to domain ID

fig, ax = plt.subplots(figsize=(8, 8))
plt.grid()
plt.axvline(color='k', linewidth=1)
plt.scatter(domain_trends_norm['annual_mean'], domain_trends_norm['annual_mean'].index, facecolor='#0077BB', alpha=0.85, s=650)
plt.scatter(domain_trends_norm['annual_mean'][~domain_trends_sig['annual_mean_p']], domain_trends_norm['annual_mean'][~domain_trends_sig['annual_mean_p']].index, facecolor='w', alpha=0.85, s=150)
plt.scatter(domain_trends_norm['annual_freq'], domain_trends_norm['annual_freq'].index, facecolor='#88CCEE', alpha=0.85, s=650)
plt.scatter(domain_trends_norm['annual_freq'][~domain_trends_sig['annual_freq_p']], domain_trends_norm['annual_freq'][~domain_trends_sig['annual_freq_p']].index, facecolor='w', alpha=0.85, s=150)
plt.scatter(domain_trends_norm['annual_freq'][domain_trends_norm['annual_freq']==0], domain_trends_norm['annual_freq'][domain_trends_norm['annual_freq']==0].index, facecolor='w', alpha=0.85, s=150)
plt.scatter(domain_trends_norm['annual_var'], domain_trends_norm['annual_var'].index, facecolor='#117733', alpha=0.85, s=650)
plt.scatter(domain_trends_norm['annual_var'][~domain_trends_sig['annual_var_p']], domain_trends_norm['annual_var'][~domain_trends_sig['annual_var_p']].index, facecolor='w', alpha=0.85, s=150)
plt.scatter(domain_trends_norm['annual_cov'], domain_trends_norm['annual_cov'].index, facecolor='#44BB99', alpha=0.85, s=650)
plt.scatter(domain_trends_norm['annual_cov'][~domain_trends_sig['annual_cov_p']], domain_trends_norm['annual_cov'][~domain_trends_sig['annual_cov_p']].index, facecolor='w', alpha=0.85, s=150)
plt.xlim([-4, 4])
plt.xticks([-4, -2, 0, 2, 4], fontsize=16)
plt.xlabel('Trend (%/decade)', fontsize=20)
plt.yticks(np.arange(4, 11), [])
ax.invert_yaxis()
plt.title('Change in Precipitation Metrics', fontsize=24)
plt.savefig('../figures/interannual_variability_revised_plots/figs2.svg')
plt.show()





# #%% Figure 4: Histograms of p-values and Slopes from Variability and Coefficient of Variation
#
# # loading p-value data
# results = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/smoothed_precip_stats_50_years_1mm_threshold.csv')
# US_stations = results[results['station_id'].str.slice(0,2) == 'US']
# US_stations_sig = US_stations[US_stations['annual_var_p'] <= 0.05]
#
# # np.sum(US_stations['annual_var_slope'][US_stations['annual_var_p']<=0.05] >0)  # 184 < 0, 339 > 0
# # np.sum(US_stations['annual_cov_slope'][US_stations['annual_cov_p']<=0.05] <0)  # 248 < 0, 260 > 0
#
# # p-value histogram of variability
# plt.hist(US_stations['annual_var_p'], bins=np.arange(0, 1, 0.01), facecolor = matplotlib.colors.to_rgba('royalblue', 0.67), edgecolor='black', density=True)
# plt.axhline(5, 0, 1, color='black', alpha=0.67)
# plt.xlabel('p-value')
# plt.ylabel('density (%)')
# plt.title('P-values for Annual Variability')
# plt.show()
#
# # histogram of variability slope
# plt.hist(US_stations_sig['annual_var_slope'], bins=np.arange(-6, 6, 0.20), facecolor = matplotlib.colors.to_rgba('royalblue', 0.67), edgecolor='black', density=True, stacked=True)
# plt.xlabel('mm/year')
# plt.ylabel('density (%)')
# plt.title('Trend in Annual Variability')
# plt.show()
#
# # p-value histogram of coefficient of variation
# plt.hist(US_stations['annual_cov_p'], bins=np.arange(0, 1, 0.01), facecolor = matplotlib.colors.to_rgba('royalblue', 0.67), edgecolor='black', density=True)
# plt.axhline(5, 0, 1, color='black', alpha=0.67)
# plt.xlabel('p-value')
# plt.ylabel('density (%)')
# plt.title('P-values for Annual Coefficient of Variation')
# plt.show()
#
# # histogram of CoV slope
# plt.hist(US_stations['annual_cov_slope'][US_stations['annual_cov_p']<=0.05], bins=np.arange(-0.005, 0.005, 0.0001), facecolor = matplotlib.colors.to_rgba('royalblue', 0.67), edgecolor='black', density=True, stacked=True)
# plt.xlabel('mm/year')
# plt.ylabel('density (%)')
# plt.title('Trend in Annual Coefficient of Variation')
# plt.show()




#%% Figure: variance/coefficient of variance dependency example
# attempting with early/late daily precip from intensity analysis
domain = 5
pdf_1 = np.load('/Users/ryanharp/Documents/great_lakes_precip_variability/cleaned_data/domain_' + str(
    domain) + '_daily_both_first_half_pdf.npy')
pdf_1 = np.sort(pdf_1[~np.isnan(pdf_1)])
if domain in [1, 5, 7, 10]:
    pdf_1 = np.sort(pdf_1)[:-1]
pdf_2 = np.load('/Users/ryanharp/Documents/great_lakes_precip_variability/cleaned_data/domain_' + str(
    domain) + '_daily_both_second_half_pdf.npy')
pdf_2 = np.sort(pdf_2[~np.isnan(pdf_2)])
pdf_2 = pdf_1*1.1
# pdf_1_20 = len(pdf_1)/5
# pdf_2 = np.concatenate((pdf_1[:np.int(np.round(pdf_1_20, 0))]*1,
#          pdf_1[np.int(np.round(pdf_1_20, 0)):np.int(np.round(pdf_1_20, 0))*2]*1.01,
#          pdf_1[np.int(np.round(pdf_1_20, 0))*2:np.int(np.round(pdf_1_20, 0))*3]*1.025,
#          pdf_1[np.int(np.round(pdf_1_20, 0))*3:np.int(np.round(pdf_1_20, 0))*4]*1.055,
#          pdf_1[np.int(np.round(pdf_1_20, 0))*4:np.int(np.round(pdf_1_20, 0))*5]*1.15))
# pdf_2 = np.concatenate((pdf_1[:np.int(np.round(pdf_1_20, 0))]*1.35,
#          pdf_1[np.int(np.round(pdf_1_20, 0)):np.int(np.round(pdf_1_20, 0))*2]*1.35,
#          pdf_1[np.int(np.round(pdf_1_20, 0))*2:np.int(np.round(pdf_1_20, 0))*3]*1.25,
#          pdf_1[np.int(np.round(pdf_1_20, 0))*3:np.int(np.round(pdf_1_20, 0))*4]*1.17,
#          pdf_1[np.int(np.round(pdf_1_20, 0))*4:np.int(np.round(pdf_1_20, 0))*5]*1))

# actual bootstrapping example starts here
bootstrap_num = 1000
len_pdf = len(pdf_1)
low_freq = np.zeros([100,bootstrap_num])
high_freq = np.zeros([100,bootstrap_num])
low_freq_mean = np.zeros([bootstrap_num])
high_freq_mean = np.zeros([bootstrap_num])
low_freq_std = np.zeros([bootstrap_num])
high_freq_std = np.zeros([bootstrap_num])
low_freq_cov = np.zeros([bootstrap_num])
high_freq_cov = np.zeros([bootstrap_num])

for j in np.arange(1,bootstrap_num+1,1):
    for i in np.arange(1, 101, 1):
        low_freq[i-1, j-1] = np.random.choice(pdf_1, size = 100, replace = True).sum()
        high_freq[i-1, j-1] = np.random.choice(pdf_1, size = 110, replace = True).sum()
    low_freq_mean[j-1] = low_freq[:,j-1].mean()
    high_freq_mean[j-1] = high_freq[:,j-1].mean()
    low_freq_std[j-1] = low_freq[:,j-1].std()
    high_freq_std[j-1] = high_freq[:,j-1].std()
    low_freq_cov[j-1] = low_freq[:,j-1].std()/low_freq[:,j-1].mean()
    high_freq_cov[j-1] = high_freq[:,j-1].std()/high_freq[:,j-1].mean()

# std_diff_per = (high_freq_std-low_freq_std)/low_freq_std.mean()*100
# x_max = np.max([np.abs(np.max(std_diff_per)), np.abs(np.min(std_diff_per))])
# plt.hist(std_diff_per, alpha=0.7, density=True, bins=50)
# plt.axvline(np.mean(std_diff_per), linewidth=1.5)
# plt.axvline(x=0, c='k', linewidth=0.5)
# # plt.xlim(-x_max*1.1, x_max*1.1)
# plt.xlim(-40, 40)
# plt.show()
#
# cov_diff_per = (high_freq_cov-low_freq_cov)/low_freq_cov.mean()*100
# x_max = np.max([np.abs(np.max(cov_diff_per)), np.abs(np.min(cov_diff_per))])
# plt.hist(cov_diff_per, alpha=0.7, density=True, bins=50)
# plt.axvline(np.mean(cov_diff_per), linewidth=1.5)
# plt.axvline(x=0, c='k', linewidth=0.5)
# # plt.xlim(-x_max*1.1, x_max*1.1)
# plt.xlim(-40, 40)
# plt.show()


low_flat = low_freq.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq.flatten()
high_params = skewnorm.fit(high_flat)

x = np.linspace(500, 1400, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3a.svg')
plt.show()


low_flat = low_freq_std.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq_std.flatten()
high_params = skewnorm.fit(high_flat)

x = np.linspace(70, 140, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3b.svg')
plt.show()


low_flat = low_freq_cov.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq_cov.flatten()
high_params = skewnorm.fit(high_flat)

x = np.linspace(0.085, 0.145, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xlim([0.085, 0.145])
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3c.svg')
plt.show()








bootstrap_num = 1000
len_pdf = len(pdf_1)
low_freq = np.zeros([100,bootstrap_num])
high_freq = np.zeros([100,bootstrap_num])
low_freq_mean = np.zeros([bootstrap_num])
high_freq_mean = np.zeros([bootstrap_num])
low_freq_std = np.zeros([bootstrap_num])
high_freq_std = np.zeros([bootstrap_num])
low_freq_cov = np.zeros([bootstrap_num])
high_freq_cov = np.zeros([bootstrap_num])

for j in np.arange(1,bootstrap_num+1,1):
    for i in np.arange(1,101,1):
        low_freq[i-1, j-1] = np.random.choice(pdf_1, size = 100, replace = True).sum()
        high_freq[i-1, j-1] = np.random.choice(pdf_2, size = 100, replace = True).sum()
    low_freq_mean[j-1] = low_freq[:,j-1].mean()
    high_freq_mean[j-1] = high_freq[:,j-1].mean()
    low_freq_std[j-1] = low_freq[:,j-1].std()
    high_freq_std[j-1] = high_freq[:,j-1].std()
    low_freq_cov[j-1] = low_freq[:,j-1].std()/low_freq[:,j-1].mean()
    high_freq_cov[j-1] = high_freq[:,j-1].std()/high_freq[:,j-1].mean()

# std_diff_per = (high_freq_std-low_freq_std)/low_freq_std.mean()*100
# x_max = np.max([np.abs(np.max(std_diff_per)), np.abs(np.min(std_diff_per))])
# plt.hist(std_diff_per, alpha=0.7, density=True, bins=25)
# plt.axvline(np.mean(std_diff_per), linewidth=1.5)
# plt.axvline(x=0, c='k', linewidth=0.5)
# # plt.xlim(-x_max*1.1, x_max*1.1)
# plt.xlim(-40, 40)
# plt.show()
#
# cov_diff_per = (high_freq_cov-low_freq_cov)/low_freq_cov.mean()*100
# x_max = np.max([np.abs(np.max(cov_diff_per)), np.abs(np.min(cov_diff_per))])
# plt.hist(cov_diff_per, alpha=0.7, density=True, bins=25)
# plt.axvline(np.mean(cov_diff_per), linewidth=1.5)
# plt.axvline(x=0, c='k', linewidth=0.5)
# # plt.xlim(-x_max*1.1, x_max*1.1)
# plt.xlim(-40, 40)
# plt.show()


low_flat = low_freq.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq.flatten()
high_params = skewnorm.fit(high_flat)

x = np.linspace(500, 1400, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3d.svg')
plt.show()


low_flat = low_freq_std.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq_std.flatten()
high_params = skewnorm.fit(high_flat)

x = np.linspace(70, 140, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3e.svg')
plt.show()


low_flat = low_freq_cov.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq_cov.flatten()
high_params = skewnorm.fit(high_flat)

x = np.linspace(0.085, 0.145, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xlim([0.085, 0.145])
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3f.svg')
plt.show()






bootstrap_num = 1000
len_pdf = len(pdf_1)
low_freq = np.zeros([100,bootstrap_num])
high_freq = np.zeros([100,bootstrap_num])
low_freq_mean = np.zeros([bootstrap_num])
high_freq_mean = np.zeros([bootstrap_num])
low_freq_std = np.zeros([bootstrap_num])
high_freq_std = np.zeros([bootstrap_num])
low_freq_cov = np.zeros([bootstrap_num])
high_freq_cov = np.zeros([bootstrap_num])

for j in np.arange(1,bootstrap_num+1,1):
    for i in np.arange(1,101,1):
        low_freq[i-1, j-1] = np.random.choice(pdf_1, size = 100, replace = True).sum()
        high_freq[i-1, j-1] = np.random.choice(pdf_2, size = 110, replace = True).sum()
    low_freq_mean[j-1] = low_freq[:,j-1].mean()
    high_freq_mean[j-1] = high_freq[:,j-1].mean()
    low_freq_std[j-1] = low_freq[:,j-1].std()
    high_freq_std[j-1] = high_freq[:,j-1].std()
    low_freq_cov[j-1] = low_freq[:,j-1].std()/low_freq[:,j-1].mean()
    high_freq_cov[j-1] = high_freq[:,j-1].std()/high_freq[:,j-1].mean()

# std_diff_per = (high_freq_std-low_freq_std)/low_freq_std.mean()*100
# x_max = np.max([np.abs(np.max(std_diff_per)), np.abs(np.min(std_diff_per))])
# plt.hist(std_diff_per, alpha=0.7, density=True, bins=25)
# plt.axvline(np.mean(std_diff_per), linewidth=1.5)
# plt.axvline(x=0, c='k', linewidth=0.5)
# # plt.xlim(-x_max*1.1, x_max*1.1)
# plt.xlim(-40, 40)
# plt.show()
#
# cov_diff_per = (high_freq_cov-low_freq_cov)/low_freq_cov.mean()*100
# x_max = np.max([np.abs(np.max(cov_diff_per)), np.abs(np.min(cov_diff_per))])
# plt.hist(cov_diff_per, alpha=0.7, density=True, bins=25)
# plt.axvline(np.mean(cov_diff_per), linewidth=1.5)
# plt.axvline(x=0, c='k', linewidth=0.5)
# # plt.xlim(-x_max*1.1, x_max*1.1)
# plt.xlim(-40, 40)
# plt.show()


low_flat = low_freq.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq.flatten()
high_params = skewnorm.fit(high_flat)

# x = np.linspace(np.quantile(low_flat, q=0.0001), np.quantile(high_flat, q=0.9999))
x = np.linspace(500, 1400, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3g.svg')
plt.show()


low_flat = low_freq_std.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq_std.flatten()
high_params = skewnorm.fit(high_flat)

# x = np.linspace(np.quantile(low_flat, q=0.001), np.quantile(high_flat, q=0.999))
x = np.linspace(70, 140, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3h.svg')
plt.show()


low_flat = low_freq_cov.flatten()
low_params = skewnorm.fit(low_flat)

high_flat = high_freq_cov.flatten()
high_params = skewnorm.fit(high_flat)

# x = np.linspace(np.quantile(low_flat, q=0.001), np.quantile(high_flat, q=0.999), 1000)
x = np.linspace(0.085, 0.145, 1000)
low_p = skewnorm.pdf(x, low_params[0], low_params[1], low_params[2])
high_p = skewnorm.pdf(x, high_params[0], high_params[1], high_params[2])

plt.plot(x, low_p, c = 'slategrey', alpha=0.5, linewidth=2.5)
plt.plot(x, high_p, c = 'mediumblue', linewidth=2.5)
plt.xlim([0.085, 0.145])
plt.xticks([])
plt.yticks([])
plt.savefig('../figures/interannual_variability_revised_plots/fig3i.svg')
plt.show()




np.mean(pdf_2)
np.std(pdf_2)
ss.skew(pdf_2)
ss.kurtosis(pdf_2)

np.mean(pdf_3)
np.std(pdf_3)
ss.skew(pdf_3)
ss.kurtosis(pdf_3)

np.mean(pdf_4)
np.std(pdf_4)
ss.skew(pdf_4)
ss.kurtosis(pdf_4)



#%% Figure S4: GHCN Station Scatter
#
# NEON_domains = gpd.read_file('/Users/ryanharp/Documents/great_lakes_precip_variability/supplemental_files/NEONDomains_0/NEON_Domains.shp')
# NEON_domains = NEON_domains.drop([1, 5], axis = 0)  # dropping 'extra' domain portions for centroid graphing purposes
# NEON_domains.index = np.arange(1, 21, 1)
# NEON_domains = NEON_domains.drop([18, 19, 20], axis = 0)
#
# domain_trends = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/NEON_domain_results/NEON_domain_precip_stats_v2.csv')
# domain_trends = domain_trends.drop(labels = 'Unnamed: 0', axis = 1)  # dropping repeat index column
# domain_trends.index = np.arange(1, 21, 1)  # setting index to be equal to domain ID
# domain_trends = domain_trends.drop([18, 19, 20], axis = 0)
#
# stations = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/precip_stats_50_years_1mm_threshold_v2.csv')
# # stations = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/results/smoothed_precip_stats_50_years_1mm_threshold.csv')
#
# results = stations[stations['station_id'].str.slice(0,2) == 'US']  # US only
#
#
# sig_value = 1
# var = 'annual_freq'
# cmap = 'bwr_r'
# # cmap = 'PRGn_r'
# results_sig = results[var + '_p']
# sig_results = results[results[var + '_p'] < sig_value]
# sig_results_p = sig_results[var + '_p']
# sig_results_slope = sig_results[var + '_slope']*10
# non_sig_results = results[results[var + '_p'] >= sig_value]
# non_sig_results_p = non_sig_results[var + '_p']
# non_sig_results_slope = non_sig_results[var + '_slope']
#
# crs_new = ccrs.PlateCarree()
# states_provinces = cfeature.NaturalEarthFeature(
#     category='cultural',
#     name='admin_1_states_provinces_lines',
#     scale='50m',
#     facecolor='none')
#
# fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(14, 8))
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(states_provinces, edgecolor='gray')
# NEON_domains.plot(ax=ax, color='white', edgecolor='black')  # plotting neon boundaries
# c_low = np.nanquantile(sig_results_slope, q=0.10)
# c_high = np.nanquantile(sig_results_slope, q=0.90)
# c_val = np.max([abs(c_low), abs(c_high)])
#
# # adding in non-significant station locations
# plt.scatter(
#     x=non_sig_results['longitude'],
#     y=non_sig_results['latitude'],
#     color='grey',
#     s=10,
#     alpha=0.25,
#     transform=ccrs.PlateCarree()
# )
# # modify the plot by adding a scatterplot over the map
# plt.scatter(
#     x=sig_results['longitude'],
#     y=sig_results['latitude'],
#     c=sig_results_slope,
#     s=20,
#     alpha=1,
#     # edgecolors='black',
#     linewidths=0.5,
#     transform=ccrs.PlateCarree()
# )
#
# plt.set_cmap(cmap)
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-c_val, c_val))
# # sm = plt.cm.ScalarMappable(cmap=cmap)
# plt.clim(-c_val, c_val)
# cbar = fig.colorbar(sm, fraction=0.0275, pad=0.04)
# cbar.set_label('Change in Precipitation Frequency (days/decade)', fontsize=20)
# cbar.ax.tick_params(labelsize=16)
# # cbar.set_label('Change in Relative Variability ($\mathregular{decade^{ -1}}$)', fontsize=20)
# plt.xlim([-129, -62])
# plt.ylim([23, 52])
# plt.xticks(np.arange(-125, -55, 10), fontsize=18)
# plt.yticks(np.arange(30, 60, 10), fontsize=18)
# plt.title('Trend in Annual Precipitation Frequency', fontsize=24)
# plt.show()





#%% Old Code
#%% Figure 1: Boundaries of the NEON Domains with Stations and Histogram

# # loading NEON domain shapefile
# neon_gdf = gpd.read_file('/Users/ryanharp/Documents/great_lakes_precip_variability/supplemental_files/NEONDomains_0/NEON_Domains.shp')
#
# # loading and filtering station data
# num_qual_years = 50
# percent_qualifying_years = 90
# ghcn_stations = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/ghcn_summary.csv')
# US_ghcn_stations = ghcn_stations[ghcn_stations['station_id'].str.slice(0,2) == 'US']  # US only
# qual_stations = US_ghcn_stations[(US_ghcn_stations['num_qual_years'] > num_qual_years) & (ghcn_stations['pct_qual_years'] >= percent_qualifying_years)]  # checking for completeness of record
#
# # prepping plotting boundaries
# states_provinces = cfeature.NaturalEarthFeature(
#     category='cultural',
#     name='admin_1_states_provinces_lines',
#     scale='50m',
#     facecolor='none')
# crs_new = ccrs.PlateCarree()
#
# # plotting country, state boundaries
# fig, ax = plt.subplots(1, 1, subplot_kw={'projection': crs_new}, figsize=(17, 11))
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.5)
# ax.set_aspect('equal')
# # plotting NEON domain boundaries
# neon_gdf.plot(ax=ax, color='white', edgecolor='black')  # plotting neon boundaries
# # plotting location of qualifying stations
# plt.scatter(
#     x = qual_stations['longitude'],
#     y = qual_stations['latitude'],
#     color='cornflowerblue',
#     marker = '.',
#     s = 8,
#     alpha = 0.75,
#     transform=ccrs.PlateCarree()
# )
# # framing around US
# plt.xlim([-129, -62])
# plt.ylim([23, 52])
# plt.xticks(np.arange(-125, -55, 10), fontsize=18)
# plt.yticks(np.arange(30, 60, 10), fontsize=18)
# plt.title('NEON Ecoregions and GHCN Station Locations', fontsize=24)
# plt.show()
#
#
# # plotting histogram of length of records for qualifying stations
# fig, ax = plt.subplots(figsize = (10, 12))
# plt.hist(qual_stations['num_qual_years'], bins = np.arange(50, 150, 5), facecolor=matplotlib.colors.to_rgba('royalblue', 0.8), edgecolor = 'black')
# plt.xticks(np.arange(60, 150, 20), fontsize=20)
# plt.xlabel('Number of Qualifying Years', fontsize=20)
# plt.yticks(np.arange(0, 170, 20), fontsize=20)
# plt.ylabel('Number of Stations', fontsize=20)
# plt.title('Length of Station Records', fontsize=24)
# plt.show()
#
#
# # plotting histogram of number of records within each domain
# stations_with_domain = pd.read_csv('/Users/ryanharp/Documents/great_lakes_precip_variability/cleaned_data/NEON_domain_daily_precip_stats.csv')
# domain_station_counts = stations_with_domain.groupby('station_domain').count()
# domain_station_counts = domain_station_counts['station_id']
#
# fig, ax = plt.subplots(figsize=(10, 12))
# ax.barh(np.arange(1, 18, 1), domain_station_counts[1:17],
#         facecolor=matplotlib.colors.to_rgba('royalblue', 0.8), edgecolor='black')
# ax.invert_yaxis()
# # ax.invert_xaxis()
# plt.yticks(np.arange(1, 18, 1))
# # plt.yticks([])
# ax.axes.yaxis.set_ticklabels([])
# plt.xticks(np.arange(0, 150, 25), fontsize=20)
# ax.axvline(color='k')
# plt.xlabel('Number of Stations', fontsize=20)
# plt.title('Stations per Region', fontsize=24)
# plt.show()
#
#

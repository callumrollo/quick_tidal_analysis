T = readtable('../sample_data/omura_temp_200m.csv');
time_yday = T.time_yday;
temperature = T.temp_200m;

[xtrend,xamp,xph,xts,xr,stats] = quick_tidal_analysis(temperature,time_yday);
figure
hold on

plot(time_yday, temperature-mean(temperature))

plot(time_yday, xr)
plot(time_yday, xts)
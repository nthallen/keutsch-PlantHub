function [ t_avg, data_avg , means_avg ] = binavg_plant ( t, data, t_res )

% Constraints:
% (1) Bin size of each array must equal t_res since raw data is 1 Hz
% (2) Each bin must have less than 50% NaNs in order to be averaged

t0=min(t);
tend=max(t);
j=1;
for i = t0:t_res:tend
    precheck = data( find( t >= i & t < i+t_res ) ,: );
    if size(precheck,1) == t_res && mean(isnan(precheck)) < 0.65
            means_avg(j) = mean(isnan(precheck));
            data_avg(j,:)= nanmean( precheck ,1 );
            t_avg(j)=i+t_res/2;
            j=j+1;
    else
            means_avg(j)= NaN;
            data_avg(j,:)= NaN;
            t_avg(j)=i+t_res/2;
            j=j+1;
    end
end

figure,plot(t_res*means_avg,'.','MarkerSize',15)
ylabel('Number of NaNs in Bin')

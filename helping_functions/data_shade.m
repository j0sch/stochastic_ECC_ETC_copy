function [mean_data,confidence]=data_shade(data,dim)

mean_data=nanmean(data,dim);
stds=nanstd(data,0,dim);

confidence=[mean_data-1.96*stds; flip(mean_data+1.96*stds)];

end


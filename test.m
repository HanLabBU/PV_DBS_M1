pv_data_path = ['~/Projects' f 'PV DBS Project' f 'PV_Data' f];

ses = dir([pv_data_path '*.mat']);
all_matfiles = {ses.name};

for i=1:length(all_matfiles)
    data = load([pv_data_path all_matfiles{i}]);
    try
        if size(data.align.trial{1}.detrend_traces, 2) > 1
            disp(all_matfiles{i});
        end
    catch
        try
            if size(data.align.trial{2}.detrend_traces, 2) > 1
                disp(all_matfiles{i});
            end
        catch
            if size(data.align.trial{3}.detrend_traces, 2) > 1
                disp(all_matfiles{i});
            end
        end
    end
end

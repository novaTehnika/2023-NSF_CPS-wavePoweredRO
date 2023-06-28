function [] = startParPool(n_cpu)
    flag = 1;
    while(flag)
        try
            p = gcp('nocreate');
            if isempty(p); parpool('local',n_cpu); end
            flag = 0;
        catch
            flag = 1;
        end
    
    end
end
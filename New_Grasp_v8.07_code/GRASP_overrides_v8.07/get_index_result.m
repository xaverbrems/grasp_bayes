function foreimage = get_index_result(index)
% index = [nbr dpth]
% or for backward compatibility [idx nbr depth]
%ATH updated 05/04/13 to avoid correcting bayesian results
%ATH updated 05/04/13 to GRASP v6.70
%could replace get_selector_result with no options

%***** Worksheet types *****
% 1 = sample scattering
% 2 = sample background
% 3 = sample cadmium
% 4 = sample transmission
% 5 = sample empty transmission
% 6 = sample empty beam transmission
% 7 = sample mask
% 8 = I0 Beam Intensity
% 99 = detector efficiency map
% 10 = FR Trans Check
% 11,12,13,14 = PA Sample Sampleing ++ -+ -- +-
% 19,20,21 Bayesian results

global status_flags
global inst_params
global grasp_env

if length(index)==2
    index = [1 index];
end

    %Back up selector values to be overwritten
    selector.fw = status_flags.selector.fw;
    selector.fn = status_flags.selector.fn;
    selector.fd = status_flags.selector.fd;

    selector.bw = status_flags.selector.bw;
    selector.bn = status_flags.selector.bn;
    selector.bd = status_flags.selector.bd;
    % override states flags - to be reset at end  
    % always get Sample and Background
    status_flags.selector.fw = index(1);
    status_flags.selector.fn = index(2);
    status_flags.selector.fd = index(3);
    status_flags.selector.bw = index(1)+1;
    status_flags.selector.bn = index(2);
    status_flags.selector.bd = index(3);


    foreimage = get_selector_result;
    

%*****  reset status flags

    status_flags.selector.fw = selector.fw;
    status_flags.selector.fn = selector.fn;
    status_flags.selector.fd = selector.fd;
    status_flags.selector.bw = selector.bw;
    status_flags.selector.bn = selector.bn;
    status_flags.selector.bd = selector.bd;
end

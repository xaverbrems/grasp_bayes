%% run_OO.m
% Object oriented version of Bayesian rock.  Allows rocking width and angle
% offsets to be optimised.  Comment or condition out lines to skip this step if
% you want.
% Requires following .m files:
% Bayesian_Rock_new.m
% minFunc_012
% 

global inst_params
global status_flags
global grasp_handles
status_flags.command_window.display_params=0;
status_flags.display.refresh = 0;

% check whether GUI bayes window is open and take bayes_input params from
% there

if ishandle(grasp_handles.window_modules.bayes.window)
    
    disp('GUI bayes window is open: Data taken from GUI')
    reset = status_flags.user_modules.bayes.reset % reads in data again - set to 1 if you have new data or have modified code, 0 if you are using same data again.
    input_index = status_flags.user_modules.bayes.input_index ; %data location  %nb. will be overwritten if set below
    output_index = status_flags.user_modules.bayes.output_index; %output location
    calc_cumulative = status_flags.user_modules.bayes.calc_cumulative;  %show cumulative results
    intensity_sd_multiplier = status_flags.user_modules.bayes.intensity_sd_multiplier; %scales default intensity s.d.
    nonsensefactor = status_flags.user_modules.bayes.nonsensefactor;  %mask unmeasured points, low value is strict, 1 is no masking.
    shape = status_flags.user_modules.bayes.shape; %'l' Lorentzian, 'g' Gaussian
    norm = status_flags.user_modules.bayes.norm; %'h' or 'a' for for height or area normalised scaling factors
    inputs = {strcat('GUI: ', status_flags.user_modules.bayes.input_name)}; %Cell array for multiple measurments
    informative_prior = status_flags.user_modules.bayes.informative_prior; %use previous posterior for prior
    pixel_prior = status_flags.user_modules.bayes.pixel_prior ; % uses individual pixel errors for prior.
    masktype = status_flags.user_modules.bayes.boxing_type;% only use sectors/sector boxes for fit, recommended. Syntax {0,'sectors' or 'sector_boxes'}
    sanoffset = status_flags.user_modules.bayes.sanoffset;
    phioffset = status_flags.user_modules.bayes.phioffset;%n.b. fminunc is best, but needs optimization toolbox
    spot = [status_flags.user_modules.bayes.spot_x status_flags.user_modules.bayes.spot_y];
    eta0 = status_flags.user_modules.bayes.eta0;
    rock_type = status_flags.user_modules.bayes.rock_type;
    
    
    if strcmp(shape,'l')
     shape = 'lorentz'
    else
	   shape = 'gauss'
    end
    A.shape = [shape norm]
    
    fit=1;  %fit fwhm and/or offsets. 
    fixed = [0 0 0];  %for fit: fixed = 1: [rocking_width sanoffset phioffset].  Can be set individually below.
    calcerrors = 1; %calculate errors on fit if not provided (i.e. no optimisation toolbox - takes extra time)
    fitmethod = 'check'; %'check','fminsearch','fminlbfgs','fminunc','minFunc'; 'check' uses 'fminunc' if available, otherwise 'minFunc'
else


reset = 1 % reads in data again - set to 1 if you have new data or have modified code, 0 if you are using same data again.
input_index = 1; %data location  %nb. will be overwritten if set below
output_index = 1; %output location
calc_cumulative = 1;  %show cumulative results
intensity_sd_multiplier =1; %scales default intensity s.d.
nonsensefactor = 1;  %mask unmeasured points, low value is strict, 1 is no masking.
shape = 'g'; %'l' Lorentzian, 'g' Gaussian
norm = 'a'; %'h' or 'a' for for height or area normalised scaling factors
inputs = {'Heidi omega 175mT'} %Cell array for multiple measurments
informative_prior = 0; %use previous posterior for prior
pixel_prior = 0; % uses individual pixel errors for prior.

    if strcmp(shape,'l')
     shape = 'lorentz'
    else
	   shape = 'gauss'
    end
    A.shape = [shape norm]

end

for i = 1:length(inputs)
    
    
    %first check whether GUI is active
    if ~ishandle(grasp_handles.window_modules.bayes.window)
        %general options set below, can be overridden for each input (e.g. to
        %fix sanoffset and phioffset after first fit).
        %% set fitting options
        fit=1;  %fit fwhm and/or offsets. 
        fixed = [0 0 0];  %for fit: fixed = 1: [rocking_width sanoffset phioffset].  Can be set individually below.
        calcerrors = 1; %calculate errors on fit if not provided (i.e. no optimisation toolbox - takes extra time)
        masktype = 0;% only use sectors/sector boxes for fit, recommended. Syntax {0,'sectors' or 'sector_boxes'}
        fitmethod = 'check'; %'check','fminsearch','fminlbfgs','fminunc','minFunc'; 'check' uses 'fminunc' if available, otherwise 'minFunc'
        sanoffset = 0;
        phioffset = 0;%n.b. fminunc is best, but needs optimization toolbox
        %fminlbfgs sometimes causes trouble with 3 free params if you are
        %already nearly at the right answer.
    
        %% keep record of previous data in 'input_record.m'.
    
        [input_index, output_index, eta0, spot, rock_type, sanoffset, phioffset] = input_record(inputs{i})
        
    end

    %% set up object
    if exist('A') && informative_prior
        posterior = A.posterior;
    end
    
    if reset || ~exist('A')
        clear A;        % clear previous object (needed if code changes)
        A = Bayesian_Rock_new;  %create object
        A = GetGraspData(A,input_index,1,'last');             %read in data from Grasp
        A = CalcErrors(A,'2d');             %CalcErrors(obj,errormode,consterr) create error matrices '3d' '2d' or 'const' with 
        A = InitialiseGraspOutput(A, input_index, output_index,calc_cumulative); %propagate params into output space (obj, input_index, output_index)
    end
    
    %% create mask
    if masktype ~= 0
        A = CreateUserMask(A,masktype);
    else
        A = CreateUserMask(A);
    end
    %% set priors
    if informative_prior
        A.prior = posterior;
    else
        A.prior.intensity.mean = 0;
        A.prior.rocking_fwhm.mean = rockwidthtoq(A,eta0,spot,rock_type);
        A.prior.rocking_fwhm.sd = A.prior.rocking_fwhm.mean*2;
        A.prior.sanoffset.mean = sanoffset*pi/180;         %n.b. these are in radians
        A.prior.sanoffset.sd = 10*pi/180;
        A.prior.phioffset.mean = phioffset*pi/180;
        A.prior.phioffset.sd = 10*pi/180;
        
        % Calculate prior sd from maximum spot intensity
        A_params_temp = A.rocking_data.params(:,1)
        if strcmp(status_flags.normalization.status,'none') %i.e. abs counts
            mon=1; stmon=1;
        elseif strcmp(status_flags.normalization.status,'mon') %i.e. normalise data to standard monitor
            divider = cell2mat(cellfun(@(s)s.monitor,A_params_temp,'uni',0)); divider_standard = status_flags.normalization.standard_monitor;
        elseif strcmp(status_flags.normalization.status,'mon2') %i.e. normalise data to standard monitor
            divider = cell2mat(cellfun(@(s)s.monitor2,A_params_temp,'uni',0)); divider_standard = status_flags.normalization.standard_monitor;        
        end
        if pixel_prior
            priormethod = 'pixelmax'
        else
           priormethod = 'max'
        end
        %A.prior.intensity.sd = intensity_sd_multiplier*0.5*pi*A.prior.rocking_fwhm.mean*(0.5+rawmax-0.5*sqrt(1+4*rawmax))*divider_standard/divider;   % use max value - error for bayes sigma
       
        A = setPriorIntensitysd(A,priormethod,intensity_sd_multiplier)
        %A.prior.intensity.sd = intensity_sd_multiplier.*rawmax
%         if strcmp(shape,'la')
%             A.prior.intensity.sd=A.prior.intensity.sd*A.prior.rocking_fwhm.mean*pi/2
%         end
%         
    end
    
    %% Fit eta0, sanoffset, phioffset if requested
    if fit
        A = fitParams(A, fixed,fitmethod,calcerrors);% no method = fminunc
        
        optimalparams = ...
             [A.prior.intensity.mean;...
                 A.prior.intensity.sd;...
                 A.posterior.rocking_fwhm.mean;...
                 A.posterior.sanoffset.mean;...
                 A.posterior.phioffset.mean;];
    else
        A.posterior.rocking_fwhm.mean = A.prior.rocking_fwhm.mean;  % if not fitting
        A.posterior.rocking_fwhm.sd = A.prior.rocking_fwhm.sd;  % if not fitting
        A.posterior.sanoffset.mean =  A.prior.sanoffset.mean;
        A.posterior.sanoffset.sd =  A.prior.sanoffset.sd;
        A.posterior.phioffset.mean = A.prior.phioffset.mean;
        A.posterior.phioffset.sd = A.prior.phioffset.sd;
        optimalparams = A.prior;
        %          [A.prior.intensity.mean;...
        %          A.prior.intensity.sd;...
        %          A.prior.rocking_fwhm.mean;...
        %          A.prior.sanoffset.mean;...
        %          A.prior.phioffset.mean;];
    end
    %% use fit result to calculate intensities
    if masktype ~= 0
        A = CreateUserMask(A);
    end
    
    A = calcIntensities(A,optimalparams,calc_cumulative);%calcIntensities(obj,params,shape,calc_cumulative)
    
    eta0_post = qtorockwidth(A,A.posterior.rocking_fwhm.mean,spot,rock_type)
    
    if fit
        f = fopen('./results.txt','a');
        fprintf(f,'%s\t',inputs{i});
        fprintf(f,'%g\t',A.posterior.rocking_fwhm.mean);
        fprintf(f,'%g\t',A.posterior.rocking_fwhm.sd);
        fprintf(f,'%g\t',eta0_post);
        eta0_post_err = abs((eta0_post-qtorockwidth(A,(A.posterior.rocking_fwhm.mean+A.posterior.rocking_fwhm.sd),spot,rock_type)))
        disp(['sanoffset = ' num2str(A.posterior.sanoffset.mean*180/pi)]);
        disp(['sanoffseterr = ' num2str(A.posterior.sanoffset.sd*180/pi)]);
        disp(['phioffset = ' num2str(A.posterior.phioffset.mean*180/pi)]);
        disp(['phioffseterr = ' num2str(A.posterior.phioffset.sd*180/pi)]);
        fprintf(f,'%g\t',eta0_post_err);
        fprintf(f,'%g\t',A.posterior.sanoffset.mean*180/pi);
        fprintf(f,'%g\t',A.posterior.sanoffset.sd*180/pi);
        fprintf(f,'%g\t',A.posterior.phioffset.mean*180/pi);
        fprintf(f,'%g\n',A.posterior.phioffset.sd*180/pi);
        fclose(f);
    end
    updateGrasp(A);
    
    informative_prior = 1; %use previous posterior for prior
end

%% hide unrocked regions
A = hidenonsense(A,nonsensefactor,0);
%A.posterior.intensity.mean = A.posterior.intensity.mean.*exp(-64*(A.posterior.intensity.sd./A.prior.intensity.sd).^2)
%A.posterior.intensity.mean = A.posterior.intensity.mean.*exp(-2*(A.posterior.intensity.sd./A.posterior.intensity.mean).^2)

%% update Grasp display and show results in degrees
status_flags.command_window.display_params=1;
status_flags.display.refresh = 1;
updateGrasp(A);

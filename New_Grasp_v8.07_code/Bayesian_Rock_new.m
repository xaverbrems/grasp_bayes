%% classdef Bayesian_Rock_new
classdef Bayesian_Rock_new
    %Bayesian_Rock  - Combines SANS data to give a weighted rocking image
    %   Detailed explanation goes here
    
    properties
        rocking_data
        mask
        usermask
        %prior = Bayes_Rocking_Params
        %posterior = Bayes_Rocking_Params
        prior
        posterior
        posterior_degrees %added so integrated intensity can also be displayed in counts * degrees
        cumulative = []
        cumulative_degrees = [] %added so cumujlative integrated intensity can also be displayed in counts * degrees
        weights
        shape = 'lorentza' %'lorentza','lorentzh','gaussa','gaussh'
        factor = 1
        priorfactor % correction to prior when using integrated intensity normalised scaling factors
        start_depth
        end_depth
        step
        final_result_index
        cumulative_result_index
        weighted_data_index
        weights_index
        pixels
        frames
        san
        phi
        %    sanoffset = 0
        %    phioffset = 0
        qx
        qy
        modq
        ewald
        max %maximum value over unmasked detector
        pixelmax %maximum value in each pixel
        axes = [1 0; 0 1] %rocking axes enabled and direction.  If set to 0 ignore param, but allow offset.
        
        %to be calculated:
        cumulative_posterior = []
        deltaq
        
    end % properties
    
    methods
        %% GetGraspData
        function obj = GetGraspData(obj, input_index, start_depth, end_depth)
            % GetGraspData(obj, input_index, start_depth, end_depth)
            % extras measurement data and relevant parameters from GRASP
            % global variables, and populates the following properties:
            % rocking_data, start_depth,end_depth, step, pixels, frames,
            % san, phi, qx, qy
            global grasp_data
            global inst_params
            global grasp_env
            
            %Patch error on zero counts to 1
            forezeros = grasp_data(1).data1{input_index}==0;
            backzeros = grasp_data(2).data1{input_index}==0;
            grasp_data(1).error1{input_index}(forezeros)=1;
            grasp_data(2).error1{input_index}(backzeros)=1;
            
            if nargin < 4
                end_depth = 'last';
            end
            if nargin < 3
                start_depth = 1;
            end
            
            pixels = circshift(inst_params.detector1.pixels,[0 1]);obj.pixels = pixels;  %get x and y correct and in correct format
            frames = 0;  %initialise frame count
            
            for j = 1:length(input_index)
                depth=grasp_data(1).dpth{input_index(j)};
                
                if start_depth == 'last'
                    obj.start_depth(j) = grasp_data(1).dpth{input_index(j)};
                else
                    obj.start_depth(j) = start_depth;
                end
                if end_depth == 'last'
                    obj.end_depth(j) = grasp_data(1).dpth{input_index(j)};
                else
                    obj.end_depth(j) = end_depth;
                end
                if end_depth(j) > start_depth(j)
                    obj.step(j)=1;
                else
                    obj.step(j)=-1;
                end
                frames=frames + abs(obj.end_depth(j)-obj.start_depth(j)+1);
            end
            
            obj.frames = frames;
            %initialise destination arrays
            obj.rocking_data.data = zeros([pixels frames]);
            obj.rocking_data.errors = zeros([pixels frames]);
            obj.rocking_data.params = cell(frames);
            obj.rocking_data.qmatrix = zeros([pixels 21]); %12/07/16 updated for graspv7.11
            
            j=1;
            for i = obj.start_depth(j):obj.step(j):obj.end_depth(j)
                corrected_data=get_index_result([1 input_index i+1]);
                corrected_data.data1(corrected_data.mask1~=1)=NaN;  %set masked points to NaN
                corrected_data.error1(corrected_data.mask1~=1)=NaN;  %set masked points to NaN
                obj.rocking_data.data(:,:,i) = corrected_data.data1;
                obj.rocking_data.errors(:,:,i) = corrected_data.error1;
                obj.rocking_data.params{i} = corrected_data.params1;
                obj.rocking_data.qmatrix = obj.rocking_data.qmatrix + corrected_data.qmatrix1;
            end
            obj.rocking_data.qmatrix = obj.rocking_data.qmatrix./frames;  %average qmatrix (main difference is in wavelength)
            obj.qx = obj.rocking_data.qmatrix(:,:,3);
            obj.qy = obj.rocking_data.qmatrix(:,:,4);
            obj.modq = obj.rocking_data.qmatrix(:,:,5);
            obj.mask = corrected_data.mask1;
            
            % set up axis directions and enabling for individual
            % instruments
            if (strcmp(grasp_env.inst,'ORNL_cg2') && strcmp(grasp_env.inst_option,'old_detector192y'));
                obj.axes = [-1 0; 0 1]
            elseif strcmp(grasp_env.inst,'ORNL_cg2');
                obj.axes = [-1 0;0 1]
            elseif strcmp(grasp_env.inst,'ILL_d33')
                obj.axes = [1 0;0 -1]
            elseif strcmp(grasp_env.inst,'ILL_d11')
                obj.axes = [1 0; 0 0]
            elseif strcmp(grasp_env.inst,'ILL_d22')
                obj.axes = [-1 0; 0 -1]
            elseif strcmp(grasp_env.inst,'HZB_V4')
                obj.axes = [1 0; 0 -1]
            elseif strcmp(grasp_env.inst,'FRM2_SANS_I')
                obj.axes = [-1 0; 0 -1]
            end
            
            %this version doesn't adapt to different angles san, phi, chi
            %etc --> need to check
%             if strcmp(grasp_env.inst,'HZB_V4')
%                 sanparam = inst_params.vectors.san;   %%%to update for grasp 8.07
%                 phiparam = inst_params.vectors.chi;
%             elseif strcmp(grasp_env.inst,'FRM2_SANS_I')
%                 sanparam = inst_params.vectors.omega_2b;
%                 phiparam = inst_params.vectors.chi_2b;
%             else
%                 sanparam = inst_params.vectors.san;
%                 phiparam = inst_params.vectors.phi;
%             end
            
            % convert the cell array into normal array to extract san and phi
           
            obj_params_temp = obj.rocking_data.params(:,1) % (frames x 1) cell array with struct param in it
            san_temp = cell2mat(cellfun(@(s)s.san,obj_params_temp,'uni',0)) %extract san into a normal array
            phi_temp = cell2mat(cellfun(@(s)s.phi,obj_params_temp,'uni',0)) %extract phi into a normal array
            
            obj.san = pi/180.*(obj.axes(1,1).*san_temp...%%%to update for grasp 8.07
                +obj.axes(1,2).*phi_temp);
            
            obj.phi = pi/180.*(obj.axes(2,1).*san_temp...
                +obj.axes(2,2).*phi_temp);
            
            obj.max = max(max(max(obj.rocking_data.data.*repmat(corrected_data.mask1, [1 1 depth])))); %find maximum spot intensity
            obj.pixelmax = max(obj.rocking_data.data.*repmat(corrected_data.mask1, [1 1 depth]),[],3);
            %set up Ewald sphere
            lambda = obj.rocking_data.params{1,1}.wav;%%%to update for grasp 8.14 only one single lambda value --> take from first scan
            k = (2*pi)/lambda;
            obj.ewald = k-(k^2-obj.modq.^2).^0.5;
            
            %patch zero errors back for correct box sums etc.
            grasp_data(1).error1{input_index}(forezeros) = 0;
            grasp_data(2).error1{input_index}(backzeros) = 0;
            clear forezeros backzeros
            
        end
        %% InitialiseGraspOutput
        function obj = InitialiseGraspOutput(obj, input_index, output_index,calc_cumulative) %needs to be done after GetGraspData
            %InitialiseGraspOutput(obj, input_index, output_index)
            %initialises and clears GRASP worksheets for destination of
            %calculations.  Params are copied from input range.
            %
            if nargin < 4
                calc_cumulative = 0;
            end
            global alex_grasp_extras %used to tell grasp not to correct calculated data.
            global inst_params
            if isempty(obj.rocking_data.params);
                error('InitGraspOutput:nodata','Please run GetGraspData first');
            end
            %             populate destination with empty input data to fill up fields.
            %updated ATH for v7.11
            final_result_type = 41;
            cumulative_result_type = 42;
            weights_type = 43
            
            final_result_index = [data_index(final_result_type) output_index];obj.final_result_index=final_result_index;  %12/07/16 updated for v7.11
            cumulative_result_index = [data_index(cumulative_result_type) output_index];obj.cumulative_result_index=cumulative_result_index; %12/07/16 updated for v7.11
            %            weighted_data_index = [1 output_index+2];obj.weighted_data_index=weighted_data_index;
            weights_index = [data_index(weights_type) output_index];obj.weights_index=weights_index;%12/07/16 updated for v7.11
            
            start_depth=obj.start_depth;
            step = obj.step;
            end_depth = obj.end_depth;
            
            initialise_data_arrays(final_result_type,final_result_index(2));  %updated ATH for v7.11
            %            initialise_data_arrays(final_result_index(1), final_result_index(2));
            if calc_cumulative % intermediate data
                initialise_data_arrays(cumulative_result_type,cumulative_result_index(2)); %updated ATH for v7.11
                %                initialise_data_arrays(, cumulative_result_index(2));
                %                initialise_data_arrays(weighted_data_index(1),weighted_data_index(2));
                initialise_data_arrays(weights_type,weights_index(2));  %updated ATH for v7.11
            end
            
            %single depth
            % adapt the data format so it can be processed with
            % place_data() which has the params in a struct which is placed in  the first cell of an
            % array have to loop thriough all detectors
            
            template=retrieve_data([1 input_index 2]);
                
            template.n_frames = 1;
            template.file_type = 'single frame';
            dets = inst_params.detectors;
            for det = 1:dets
                template_params_temp = template.(['params' num2str(det)]);
                template.(['params' num2str(det)]) = cell(1);
                template.(['params' num2str(det)]){1,1} = template_params_temp;
                %template.data1 = zeros(size(template.data1));
                template.(['data' num2str(det)])=zeros(size(template.(['data' num2str(det)])));
            end
            place_data(template, final_result_type,final_result_index(2),1); %updated ATH for v7.11
            
            %set "alex_grasp_extras.show_direct" flag so GRASP knows not to do further processing
            %to be used in get_selector_result.m
            %alex_grasp_extras.showdirect(final_result_index(1),final_result_index(2))=1;
            
            %multiple depths
            if calc_cumulative % intermediate data
                for i = start_depth:step:end_depth
                    
                    template=retrieve_data([1 input_index i+1]);
                    for det = 1:dets
                        %store params temporary into variable
                        template_params_temp = template.(['params' num2str(det)]);
                        template.(['params' num2str(det)]) = cell(1);
                        template.(['params' num2str(det)]){1,1} = template_params_temp;
                    
                        template.(['data' num2str(det)])=zeros(size(template.(['data' num2str(det)])));
                    end
                    template.n_frames=1;
                    template.file_type = 'single frame';
                    %                   place_data(template, weighted_data_index(1),weighted_data_index(2),i);
                    place_data(template, weights_type,weights_index(2),i);
                    place_data(template, cumulative_result_type,cumulative_result_index(2),i);
                    %       alex_grasp_extras.showdirect(cumulative_result_index(1),cumulative_result_index(2)) = 1;
                    %                    alex_grasp_extras.showdirect(weighted_data_index(1),weighted_data_index(2)) = 1;
                    %      alex_grasp_extras.showdirect(weights_index(1),weights_index(2)) = 1;
                end
            end
        end %InitialiseGraspOutput
        %%      CreateUserMask
        function obj = CreateUserMask(obj,mode)
            %Combines instrument, user and sector box masks to create logical index mask
            %for treating subset of data
            if nargin < 2
                mode = 'basic'
            end
            switch mode
                case 'basic'
                    obj.usermask = repmat(obj.mask, [1 1 obj.frames]);      %only use instrument and beam masks
                case 'sectors'
                    smask = sector_callbacks('build_sector_mask');
                    mask = obj.mask & smask.det1;
                    obj.usermask = repmat(mask, [1 1 obj.frames]);
                    clear smask mask;
                case 'sector_boxes'
                    smask =  sector_box_callbacks('build_the_masks');
                    mask = obj.mask & smask.sum_mask.det1;
                    obj.usermask = repmat(mask, [1 1 obj.frames]);
                    clear smask mask;
            end
            obj.usermask = logical(obj.usermask);
        end
        %%      GetDeltaQ
        function obj = GetDeltaQ(obj,sanoffset,phioffset)
            %sanoffset, phioffset in radians
            if sum(obj.axes(1,:))==0
                san = obj.san-sanoffset;
            else
                san = obj.san-sum(obj.axes(1,:)).*sanoffset;
            end
            if sum(obj.axes(2,:))==0
                phi = obj.phi-phioffset;
            else
                phi = obj.phi-sum(obj.axes(2,:)).*phioffset;
            end
            
            obj.deltaq = obj.deltaqcalc(obj.ewald,obj.qx,obj.qy,san,phi);
        end %GetDeltaQ
        %%      CalcErrors
        function obj = CalcErrors(obj,errormode,consterr)
            switch errormode
                case '3d'
                    return
                case '2d'
                    meanerrors = mean(obj.rocking_data.errors,3);
                    obj.rocking_data.errors = repmat(meanerrors,[1 1 obj.frames]);
                    clear meanerrors;
                case 'const'
                    ind=find(repmat(obj.mask(:,:,1),[1 1 obj.frames]));
                    if nargin <3
                        obj.rocking_data.errors =  mean(obj.rocking_data.errors(ind))*ones([obj.pixels obj.frames]);
                    else
                        obj.rocking_data.errors =  consterr*ones([obj.pixels obj.frames]);
                    end
            end
        end %CalcErrors
        %%      rockwidthtoq
        function fwhm = rockwidthtoq(obj,eta0,spot_coords,axis)
            %rockwidthtoq(eta0,spot_coords)
            %eta0 is fwhm in degrees for rocking measurement.
            %spot_coords = (x0,y0)
            %axis = 'san' (rock about vertical axis) or
            %       'phi' (about horizonal axis)
            if strcmp(axis,'san')
                fwhm = 2*abs(obj.qx(spot_coords(2),spot_coords(1))*tan(eta0*pi/360));
            elseif strcmp(axis,'phi')
                fwhm = 2*abs(obj.qy(spot_coords(2),spot_coords(1))*tan(eta0*pi/360));
            end
            
        end %rockwidthtoq
        %%      qtorockwidth
        function eta0 = qtorockwidth(obj,qfwhm,spot_coords,axis)
            %qtorockwidth(eta0,spot_coords)
            %eta0 is fwhm in degrees for rocking measurement.
            %spot_coords = (x0,y0)
            %axis = 'san' (rock about vertical axis) or
            %       'phi' (about horizonal axis)
            if strcmp(axis,'san')
                eta0 = abs(2*atan(0.5*qfwhm/obj.qx(spot_coords(2),spot_coords(1)))*180/pi);
            elseif strcmp(axis,'phi')
                eta0 = abs(2*atan(0.5*qfwhm/obj.qy(spot_coords(2),spot_coords(1)))*180/pi);
            end
        end %qtorockwidth
        %%      setPriorIntensity
        function obj = setPriorIntensitysd(obj,method,value)
            %function obj = setPriorIntensitysd(obj,method,value)
            %method is 'max' or 'pixelmax'
            %anything else uses 'value' in data units (i.e. will be
            %multiplied by fwhm to give integrated intensity if necessary)

            if strcmpi(method,'max')
                obj.prior.intensity.sd = value.*obj.max;
            elseif strcmpi(method,'pixelmax')
                obj.prior.intensity.sd = value.*obj.pixelmax;
            elseif strcmpi(method,'value')%just use value
                obj.prior.intensity.sd = value;
            else
                error('''value'' should be ''max'', ''pixelmax'', or ''value''')
            end
            if strcmp(obj.shape,'lorentza')
                 obj.prior.intensity.sd = obj.prior.intensity.sd.*obj.prior.rocking_fwhm.mean*pi*0.5;
            elseif strcmp(obj.shape,'gaussa')
                obj.prior.intensity.sd = obj.prior.intensity.sd.*obj.prior.rocking_fwhm.mean*0.424661*(2*pi)^0.5;
            end
        end
        %%      calcIntensities
        function obj = calcIntensities(obj,params,calc_cumulative)
            % function obj = calcIntensities(obj,params,calc_cumulative)
            % params is either numeric array [priormean priorsd fwhm ...
            % sanoffset phioffset]
            % or Bayes_Rocking_Params class to use previous posterior as
            % prior
            if nargin < 3
                calc_cumulative = 0;  % by default don't calculate cumulative results.
            end
            if isnumeric(params)
                priormean =  params(1);
                priorsd   =  params(2);
                fwhm      =  params(3);
                sanoffset =  params(4);
                phioffset =  params(5);
            else
                priormean =  params.intensity.mean;
                priorsd   =  params.intensity.sd;
                fwhm      =  params.rocking_fwhm.mean;
                sanoffset =  params.sanoffset.mean;
                phioffset =  params.phioffset.mean;
            end
            
            
            obj.posterior.intensity.mean = NaN(obj.pixels);
            obj.posterior.intensity.sd = NaN(obj.pixels);
            obj.factor = NaN([obj.pixels obj.frames]);
            obj.weights = NaN([obj.pixels obj.frames]);
            weights_denom3d = NaN([obj.pixels obj.frames]);
            
            %calculate deltaq
            obj = GetDeltaQ(obj,sanoffset,phioffset);
            mask = obj.usermask;
            %calculate weighting factors
            if strcmp(obj.shape,'lorentzh')
                obj.factor(mask) = obj.lorentzh(obj.deltaq(mask),fwhm*0.5);
            elseif strcmp(obj.shape,'lorentza')
                obj.factor(mask) = obj.lorentza(obj.deltaq(mask),fwhm*0.5);
                %priorsd=priorsd*0.5.*pi.*fwhm; scaling done in setPriorIntensitysd
            elseif strcmp(obj.shape,'gaussa')
                obj.factor(mask) = obj.gaussa(obj.deltaq(mask), fwhm*0.424661);
            elseif strcmp(obj.shape,'gaussh')
                obj.factor(mask) = obj.gaussh(obj.deltaq(mask), fwhm*0.424661);
            end
            
            obj.weights(mask) = obj.factor(mask)./obj.rocking_data.errors(mask).^2;
            weights_denom = sum((obj.factor.*obj.weights),3) + priorsd.^-2;
            weights_denom3d = repmat(weights_denom, [1 1 obj.frames]);
            obj.posterior.intensity.mean = (priormean.*priorsd.^-2+...
                sum(obj.weights.*obj.rocking_data.data,3))./weights_denom;
            obj.posterior.intensity.sd(mask(:,:,1)) = weights_denom(mask(:,:,1)).^-0.5;
            if calc_cumulative
                cumulative_num=priormean.*priorsd.^-2;
                cumulative_denom = priorsd.^-2;
                for i = obj.start_depth:obj.step:obj.end_depth
                    cumulative_num = cumulative_num + obj.weights(:,:,i).*obj.rocking_data.data(:,:,i);
                    cumulative_denom = cumulative_denom + obj.factor(:,:,i).*obj.weights(:,:,i);
                    obj.cumulative.mean(:,:,i) = cumulative_num./cumulative_denom;
                    obj.cumulative.sd(:,:,i) = cumulative_denom.^-0.5;
                end
            end
            obj.weights(mask) = obj.weights(mask)./weights_denom3d(mask);
            
            %calculate the posterior mean and sd in counts/degrees by just multiplying
            %it by (180/(qx^2 + qy^2)^(1/2)
            obj.posterior_degrees.intensity.mean = obj.posterior.intensity.mean .*(180./(obj.modq*pi));
            obj.posterior_degrees.intensity.sd = obj.posterior.intensity.sd .*(180./(obj.modq*pi));
            
            % do the same for the cumulative
            if calc_cumulative
                for i = obj.start_depth:obj.step:obj.end_depth
                obj.cumulative_degrees.mean(:,:,i) = obj.cumulative.mean(:,:,i).*(180./(obj.modq*pi));
                obj.cumulative_degrees.sd(:,:,i) = obj.cumulative.sd(:,:,i).*(180./(obj.modq*pi));
                end
            end
            %if status_flags.user_modules.bayes.qztodegrees == 1
            %obj.posterior.intensity.mean
            %end
            
        end %calcIntensities
        %%      LogPosterior
        function LogPosterior = LogPosterior(obj,params,use_mask)
            % LogPosterior = LogPosterior(obj,params,use_mask)
            % calculates (negative) logarithm of unnormalised posterior probability
            % given a set of parameters.  Optionally uses sector masks to
            % only consider certain pixels
            if nargin < 3, use_mask = 1; end
            
            fwhm = params(1);
            sanoffset = params(2);
            phioffset = params(3);

            if length(obj.prior.intensity.mean)==1
                intensityparams =...
                    [obj.prior.intensity.mean;...
                    obj.prior.intensity.sd;...
                    fwhm; sanoffset; phioffset;];
            else
                intensityparams = obj.prior;
                intensityparams.rocking_fwhm.mean = fwhm;
                intensityparams.sanoffset.mean = sanoffset;
                intensityparams.phioffset.mean = phioffset;
                
            end
            
            if (strcmp(obj.shape,'lorentza')||strcmp(obj.shape,'gaussa'))
                priorfactor = 1;
            elseif strcmp(obj.shape,'lorentzh')
                priorfactor = 0.5*pi*fwhm;
            elseif strcmp(obj.shape,'gaussh')
                priorfactor = 0.424661*(2*pi)^0.5*fwhm
            end
            
            if use_mask
                obj = obj.calcIntensities(intensityparams,1);
            else
                obj = obj.calcIntensities(intensityparams);
            end
            intensities = repmat(obj.posterior.intensity.mean,[1 1 obj.frames]);
            
            likelihood = 0.5*(obj.rocking_data.data - obj.factor.*intensities).^2 ...
                ./obj.rocking_data.errors.^2;
            priorintensity = 0.5*(obj.posterior.intensity.mean-obj.prior.intensity.mean).^2 ...
                .*(priorfactor./obj.prior.intensity.sd).^2;
            LogPosterior = sum(likelihood(~isnan(likelihood)))+sum(priorintensity(~isnan(priorintensity)))...           %remove NaNs
                +0.5*(sanoffset-obj.prior.sanoffset.mean)^2/obj.prior.sanoffset.sd ...
                +0.5*(phioffset-obj.prior.phioffset.mean)^2/obj.prior.phioffset.sd ...
                +0.5*(fwhm-obj.prior.rocking_fwhm.mean).^2./obj.prior.rocking_fwhm.sd.^2;
            %-log(fwhm);  %use Jeffreys prior for fwhm.
        end
        %% derivatives (not implemented)
        %         function obj = dLogPosterior(obj)
        %         end
        % %%
        %         function obj = ddLogPosterior(obj)
        %         end
        %%  fitParams
        function obj = fitParams(obj,fixed,method,calcerrors)
            if nargin < 3
                method = 'check';  %use fminunc if licence allows
            elseif nargin < 4
                calcerrors = 0;
            end
            
            initparams=[obj.prior.rocking_fwhm.mean ...
                obj.prior.sanoffset.mean ...
                obj.prior.phioffset.mean]
            
            scalefactor = [1e-2 1 1];
            scaledinitparams=initparams./scalefactor
            
            f=@(x)LogPosterior(obj,obj.interlace(x.*scalefactor(~fixed),initparams.*scalefactor,fixed));
            
            clear hessian
            switch method
                case 'check'
                    if license('test','optimization_toolbox')
                        options = optimset('Display','iter-detailed','TolX',1e-8,'Tolfun',1e-8,'DiffMaxChange',0.1);
                        [scaledparamsout,fval,exitflag,output,grad,hessianmat] = fminunc(f,scaledinitparams(~fixed),options);
                    else %use minFunc
                        %[scaledparamsout,fval,exitflag,output] = fminsearch(f,scaledinitparams(~fixed),options);
                        options = []
                        options.Method = 'lbfgs';
                        options.numDiff = 1;
                        options.useMex = 0;
                        f=@(x)LogPosterior(obj,(obj.interlace(x'.*scalefactor(~fixed),initparams.*scalefactor,fixed))');
                        [scaledparamsout,fval,exitflag,output] = minFunc(f,(scaledinitparams(~fixed))',options);
                        
                    end
                case 'fminunc'
                    [scaledparamsout,fval,exitflag,output,grad,hessianmat] = fminunc(f,scaledinitparams(~fixed),options);
                case 'fminsearch'
                    %[scaledparamsout,fval,exitflag,output] = fminsearch(f,scaledinitparams(~fixed),options);
                    scaledparamsout=scaledinitparams(~fixed);exitflag=1;
                case 'fminlbfgs'
                    options.GoalsExactAchieve = 0;
                    options.Display = 'iter';%use 'plot' to show line searches
                    [scaledparamsout,fval,exitflag,output,grad]=fminlbfgs(f,scaledinitparams(~fixed),options);
                case 'minFunc'
                    options = []
                    options.Method = 'lbfgs';
                    options.numDiff = 1;
                    options.useMex = 1;
                    f=@(x)LogPosterior(obj,(obj.interlace(x'.*scalefactor(~fixed),initparams.*scalefactor,fixed))');
                    
                    [scaledparamsout,fval,exitflag,output] = minFunc(f,(scaledinitparams(~fixed))',options);
            end
            exitflag
            
            allscaledparamsout=obj.interlace(scaledparamsout,scaledinitparams,fixed);
            obj.posterior.rocking_fwhm.sd=NaN;
            obj.posterior.sanoffset.sd=NaN;
            obj.posterior.phioffset.sd=NaN;
            if exitflag <0
                return;
            end
            obj.posterior.rocking_fwhm.mean=allscaledparamsout(1)*scalefactor(1);
            obj.posterior.sanoffset.mean=allscaledparamsout(2)*scalefactor(2);
            obj.posterior.phioffset.mean=allscaledparamsout(3)*scalefactor(3);
            
            
            if exist('hessianmat')  %calculate errors if hessian available
                covariance = inv(hessianmat);
            elseif calcerrors %calculate hessian if instructed
                disp('Calculating Hessian matrix')
                [hessianmat, err] = hessian(f,scaledparamsout);
                covariance = inv(hessianmat);
            end
            if  exist('covariance')
                for i = 1:size(covariance,1)
                    diagcovar(i)=covariance(i,i);
                end
                alldiag = obj.interlace(diagcovar,nan(size(fixed)),fixed);
                
                obj.posterior.rocking_fwhm.sd=alldiag(1)^0.5*scalefactor(1);
                obj.posterior.sanoffset.sd=alldiag(2)^0.5*scalefactor(2);
                obj.posterior.phioffset.sd=alldiag(3)^0.5*scalefactor(3);
            end
        end
        


        %% hidecrap
        function obj = hidenonsense(obj,nonsensefactor,substvalue)
            global status_flags
            global inst_params
            
            %             if max(size(obj.prior.intensity.sd))==1
            %                priorsd = obj.prior.intensity.sd
            %             else
            %                 if strcmp(status_flags.normalization.status,'none') %i.e. abs counts
            %                     mon=1; stmon=1;
            %                 else
            %                     mon = obj.rocking_data.params(inst_params.vectors.monitor); stmon = status_flags.normalization.standard_monitor;
            %                 end
            %                 rawmax = obj.max*mon/stmon;
            %                 priorsd = 0.5*pi*obj.prior.rocking_fwhm.mean*(0.5+rawmax-0.5*sqrt(1+4*rawmax))*stmon/mon;
            %             end
            %obj.posterior.intensity.mean = obj.posterior_degrees.intensity.mean
            obj.posterior.intensity.mean(obj.posterior.intensity.sd > obj.prior.intensity.sd.*nonsensefactor) = substvalue;
            obj.posterior.intensity.mean(~obj.mask)=0;
            %do the same for the degree posterior
            obj.posterior_degrees.intensity.mean(obj.posterior.intensity.sd > obj.prior.intensity.sd.*nonsensefactor) = substvalue;
            obj.posterior_degrees.mean(~obj.mask)=0;
            
        end
        
        
        %% qztodegrees
        function obj = qztodegrees(obj)
            %replaces the intensities (both integrated and cumulative) and
            %its errors with the corresponding intensities in
            %counts * degree
            obj.posterior.intensity = obj.posterior_degrees.intensity
            obj.cumulative = obj.cumulative_degrees
        end
        
        %% updateGrasp
        function updateGrasp(obj)
            global grasp_data
            %remove NaN's
            obj.weights(isnan(obj.weights))=0;
            obj.posterior.intensity.mean(isnan(obj.posterior.intensity.mean))=0;
            obj.posterior.intensity.sd(isnan(obj.posterior.intensity.sd))=0;
            obj.cumulative.mean(isnan(obj.cumulative.mean))=0;
            obj.cumulative.sd(isnan(obj.cumulative.sd))=0;
            
            
            
            grasp_data(obj.final_result_index(1)).data1{obj.final_result_index(2)} = obj.posterior.intensity.mean; %sample data, set number 2
            grasp_data(obj.final_result_index(1)).error1{obj.final_result_index(2)} = obj.posterior.intensity.sd; %sample data, set number 2
            grasp_data(obj.final_result_index(1)).attenuation= 1; %
            
            if isstruct(obj.cumulative) % intermediate data here
                obj.cumulative.mean(isnan(obj.cumulative.mean))=0;
                grasp_data(obj.weights_index(1)).data1{obj.weights_index(2)} = obj.weights; %weights => sample data, set number 3
                grasp_data(obj.cumulative_result_index(1)).data1{obj.cumulative_result_index(2)} = obj.cumulative.mean; %weights => sample data, set number 3
                grasp_data(obj.cumulative_result_index(1)).error1{obj.cumulative_result_index(2)} = obj.cumulative.sd; %weights => sample data, set number 3
                grasp_data(obj.weights_index(1)).attenuation= 1; %
                
            end
            grasp_update
        end
    end %methods
    %% static methods
    methods (Static = true)
        function output = interlace(x,a,fixed)
            %function paramscombined = interlace(x,a,fixed)
            %produces output with numel(a) elements
            %where 'fixed' is nonzero or true, elements are replaced
            %sequentially by elements of x
            fixed = logical(fixed);
            if any(fixed)
                a+fixed; % gives error if sizes don't match
                output = zeros(size(fixed));
                output(~fixed) = x;
                output(fixed) = a(fixed);
            else
                output = x;
            end
            %disp(output) %show optimization steps
        end
        
        function deltaq = deltaqcalc(ewald,qx,qy,san,phi)
            %function deltaq = deltaqcalc(ewald,qx,qy,san,phi)
            %calculates perpendicular distance from lattice plane to Ewald
            %sphere
            %
            %ewald,qx,qy are 2D arrays
            %san, phi are either scalars or column vectors
            %san is rocking angle about vertical axis
            %phi is rocking angle horizontal vertical axis
            %angles are assumed to be small enough that rotations commute
            pixels = size(ewald);
            san = sin(san);
            phi = sin(phi);
            if isvector(san)== 0 || isvector(phi) == 0
                error('deltaq: san and phi must be scalars or 1d vectors');
            elseif (length(san) > 1) && (length(phi) > 1) && any(size(san) ~= size(phi))
                error('deltaq: san and phi must be vectors of equal length');
            end
            
            if isscalar(san) && isscalar(phi)
                deltaq = ewald - qx.*san-qy.*phi;
                return
            elseif ~(isscalar(san)) && ~(isscalar(phi))
                frames = length(san);
                san = reshape(san,1,1,frames);
                phi = reshape(phi,1,1,frames);
                san = repmat(san,pixels);
                phi = repmat(phi,pixels);
                
            elseif ~(isscalar(san)) && isscalar(phi)
                frames = length(san);
                san = reshape(san,1,1,frames);
                san = repmat(san,pixels);
                
            elseif isscalar(san) && ~(isscalar(phi))
                frames = length(phi);
                phi = reshape(phi,1,1,frames);
                phi = repmat(phi,pixels);
            end
            ewald = repmat(ewald, [1 1 frames]);
            qx = repmat(qx, [1 1 frames]);
            qy = repmat(qy, [1 1 frames]);
            
            deltaq = ewald - qx.*san+qy.*phi;
            clear qx qy san phi
            
        end %deltaqcalc
        function factor = lorentzh(x,hwhm)
            %zero-centred lorentzian with peak = 1
            factor = hwhm^2./(hwhm^2+x.^2);
        end
        function factor = lorentza(x,hwhm)
            %zero-centred lorentzian with area = 1
            factor = abs(hwhm)./(pi*(hwhm^2+x.^2));
        end
        
        function factor = gaussh(x,sigma)
            factor = exp(-0.5*x.^2*sigma^-2);
        end
        function factor = gaussa(x,sigma)
            factor = (2*pi)^-0.5*sigma^-1*exp(-0.5*x.^2*sigma^-2);
        end
        
        
    end
end


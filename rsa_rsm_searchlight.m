function varargout = rsa_rsm_searchlight(args)
% Computes representational similarity analysis (RSA) by searchlight.
% This function implements RSA in a partial correlation approach, as
% described in Viviani (2021). Various terms can be specified for
% partialling out in the correlation, including the estimated covariance of
% the beta coefficients involved in the analysis ('BCov') or an estimate of
% these covariance from the masked volume ('SCov') or simply the
% cross-products of these coefficients ('BB'), or any user-specified term; 
% see Viviani (2021) for details. This function is implemented as an add-on
% to the SPM12 package (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/),
% which is a requirement for its use. By using this function, you agree to
% the terms for its use as specified in the license file enclosed to this
% add-on functions.
%
% The flavour of rsa implemented here, beyond the partial correlation
% features, is characterized by using correlation to assess concordance
% (Pearson correlation, instead of rank correlation as in the original rsa
% software by Nili et al.) for the brain signal and, implementationally, 
% that it works on similarities rather than  dissimilarities (to match the 
% equations of Viviani 2021). No multidimensional scaling or Procrustes 
% rotations of the brain signal are computed. Correspondances are assessed 
% by (partial) correlation of the off-diagonal terms of the similarity 
% matrices of the brain signal and of the stimuli, as in Nili et al. 2014. 
% For these parameters, the algorithm is expected to work identically to a 
% dissimilarity implementation, given that correlations are invariant to 
% scale and location of data.
%
% The input to the analysis is specified as a list of directories. Each
% directory should contain: 1. a set of Nifti images containing estimated
% coefficients of a first-level model or contrasts, the same in each
% directory; 3. the SPM.mat files computed by SPM (when BCov is specified 
% as a confound) 2. a saved associative array (struct) with the similarity 
% or dissimilarity maps to be included in the model as targets of the
% correlation or as confounding terms to be partialled out. As of writing, 
% this software can read the files with dissimilarity maps produced by Nili 
% et al. (2014) or by the helper funtion rsa_RSM.m. BCov, SCov, BB, if 
% required, are not specified explicitly but are computed internally.
%
% The correlation maps are written into the input directories or in 
% args.OutputDir, if specified.
%
% Usage:
% args = rsa_rsm_searchlight(ask_args); %collect arguments in args, no execution
% str = rsa_rsm_searchlight(args);      %execute
%
% str = rsa_rsm_searchlight;    %specify arguments interactively and execute
%
% Arguments:
% args.Directories      a char array containing the input directories
% args.Method           (opt) one of 'Pearson', 'Spearman', 'regression',
%                       'slsize', or empty. If empty defaults to 'Pearson'.
% args.BrainMapType     (opt) one of 'sscp', 'cov', 'cor', or empty. If
%                       empty, defaults to 'sscp'. Ignored if method is
%                       'regression'.
% args.MapFile          the name of the file containing the similarity or
%                       dissimilarity maps in each directory. If this file
%                       was produced by the Nili et al. software, it will
%                       contain dissimilarity maps that will be converted
%                       into similarity maps. If this name is a full path
%                       file name, the same maps will be used in the whole
%                       input set. If it is a simple name (without path),
%                       it will direct the function to select this file to 
%                       read in the similarity input maps in each input 
%                       directory.
% args.MapSel           (opt) a cell array of names specifying which maps
%                       from the map file should be used in the
%                       correlation. Defaults to all maps (if empty). May
%                       also contain 'BCov', 'SCov', 'BB'.
% args.MapSelPcorr      (opt) a cell array of names specifying which
%                       maps from the map file should be partialled out in
%                       the correlation. Specify 'BCov', 'SCov', 'BB' to 
%                       partial out these terms. Defaults to none (if 
%                       empty). If args.Method is 'regression', the terms
%                       specified here are simply added to the model.
% args.BetaFiles        (opt) If empty (default), directs the function to
%                       look for the beta images in the directories
%                       specified in args.Directories as input. If a simple
%                       string, it specifies a prefix to be prepended to
%                       the names of the beta images in these directories
%                       (for example, in case the images were preliminary
%                       smoothed). If a cell of a string, a regular 
%                       expression to select the nifti volumes as input to
%                       the rsa.
% args.BetaIdx          (opt) If an array of logical values or values
%                       consisting of 0 and 1, or numerical indices to
%                       select a subset of files as directed by
%                       args.BetaFiles. If empty, (default), maps to the 
%                       first beta images (or the first images as specified 
%                       by args.BetaFiles) that match the number of stimuli 
%                       in the representation maps specified in args.RsmFile.
% args.MaskFile         (opt) mask file. If empty, the file mask.nii is
%                       used in each of the input directories (this is an 
%                       SPM-generated file).
% args.MaskConfound     (opts) the mask used to compute confoundd terms
%                       such as BB or SCov.
% args.SearchlightDef   (opt) one of 'sphere', 'box'. If empty, defaults
%                       to 'sphere'.
% args.SearchlightSize  (opt) searchlight size (in mm.). If empty defaults
%                       to 8. Must be a positive value at least vox size*2.
% args.OffDiagOffset    (opt) off-diagonal offset. Defaults to zero. Zero
%                       selects all off-diagonal terms of the mx (but not 
%                       the diagonal). Numbers larger than zero exclude the
%                       band of elements near the diagonal. Must be a
%                       positive value or zero.
% args.ModelName        (opt) a model name to be included in the output. If
%                       empty, a model name will be automatically
%                       generated.
% args.OutputDir        (opt) the directory where the correlation volumes
%                       are written. If empty, the output is written in the
%                       same directories as the input files.
%
% Output:
% str.ModelName         the name of the analysis
% str.Rsms              the similarity structures of the stimuli and
%                       confounds used in the analysis
% str.Mn                the mean volume correlations of the stimuli
%                       similarity maps (study average)
% str.Mns               the mean volume correlations of the stimuli
%                       similarity maps (per subject)
% str.Cmx               correlation of predictors
% str.Output            the file where this structure was saved
% str.Args              the arguments used in the analysis
% str.disp()            a function printing the defining features of the
%                       model
% str.diagn()           a function printing the summary volume correlation
%                       diagnostics for each stimulus map
% str.boxplot()         a function displaying a boxplot of the mean volume
%                       correlations of the stimuli similarity maps
%                       (str.Mns)
% str.getmap('name')    returns the similarity map of 'name' (average
%                       across study). The individual maps are in str.rsms
% str.mapcorr('map1', 'map2') returns the correlation between the
%                       similarity maps, averaged
%
% References: 
%
% Nili, H., Wingfield, C., Walther, A., Su, L., Marslen-Wilson, W., 
% Kriegeskorte, N., 2014. A toolbox for representational similarity 
% analysis. PLoS Comp. Biol. 10, e1003553.
%
% Viviani, R., 2021. Overcoming bias in representational similarity
% analysis. arXiv, 2102.08931.
%
% 2020-21 Roberto Viviani - December 2020
% Institute of Psychology, University of Innsbruck
% Last modified: January 23rd 2021

    %-Check parameters--------------------------------------------------
    args.MFile = mfilename;
    args.Version = '1.0';
    args.SPMVersion = spm('Ver');
    args.MatlabVersion = version;
    args.Timestamp = datestr(now);
    args.Code = regexp(fileread([mfilename('fullpath'),'.m']), '\n', 'split')';
    
    %-Ask user if any parameter is missing.-----------------------------
    if nargin > 0
        args = Ui_(args, inputname(1));
    else
        args = Ui_(args, '');
    end
    if isfield(args, 'Args') 
        args = rmfield(args, 'Args'); 
        if isfield(args, 'Quit')
            args = rmfield(args, 'Quit'); end
        varargout{1} = args;
        return; 
    end
    if isfield(args, 'Quit') && args.Quit
        varargout{1} = rmfield(args, 'Quit'); return; end
    
    %-Do the real job---------------------------------------------------
    dirnames = args.Directories;
    csubj = size(dirnames, 1);
    
    %load inputs
    istr = load_input_(args);
    if isempty(istr), varargout{1} = []; return; end

    %mail loop thru subjects
    rsa = rsa_rsmio();
    mns = []; rsms = rsa.create_rsm();
    rsmcovs = rsa.create_rsm(); cmx = cell(csubj,1);
    figure(spm_figure('GetWin', 'Interactive'));
    hf = figure(spm_figure('GetWin', 'Interactive'));
    for i = 1 : csubj
        fprintf('Subject %d of %d', i, csubj);
        set(hf, 'Name', sprintf('Subject %d of %d', i, csubj));
        str_ = srchlght_(istr(i), args);
        if isempty(str_), continue; end
        mns = [mns; str_.mns];
        rsms = [rsms; str_.rsms];
        rsmcovs = [rsmcovs; str_.rsmcovs];
        cmx{i} = str_.cmx;
    end
    set(hf, 'Name', 'rsa add-on');
    if isempty(mns), varargout{1} = []; return; end
    maps = arrayfun(@(x) x.name, rsms(1,:), 'UniformOutput', false);
    mapcovs = {};
    if ~isempty(rsmcovs)
        mapcovs = arrayfun(@(x) x.name, rsmcovs(1,:), 'UniformOutput', false);
    end
    
    function meancorr_
        disp('Summary of mean volume correlations');
        mn = mean(mns, 1);
        sd = std(mns, [], 1);
        if size(mns, 1) > 1
            for i_ = 1 : size(rsms, 2)
                fprintf('%s, mean volume: %0.4f, sd: %0.4f, sign test: p=%0.4f\n', ...
                    rsms(1,i_).name, mn(i_), sd(i_), signtest(mns(:,i_)));
            end
        else
            for i_ = 1 : size(rsms, 2)
                fprintf('%s, mean volume: %0.4f\n', ...
                    rsms(1,i_).name, mn(i_));
            end
        end
    end
    meancorr_();
    
    modelname = get_modelname_(args);
    str.ModelName = modelname;
    str.MapNames = maps;
    str.Rsms = [rsms, rsmcovs];
    str.Mn = mean(mns, 1);
    str.Mns = mns;
    mx = zeros(size(cmx{1}));
    for m = 1:length(cmx), mx = mx + cmx{m}; end
    str.Cmx = mx ./ length(cmx);

    if isempty(args.OutputDir)
        utils = rsa_utils();
        args.OutputDir = utils.get_commonpath(char(istr(:).directory)); 
    end
    str.Output = fullfile(args.OutputDir, [modelname, '.mat']);

    str.Args = args;
    str.diagn = @meancorr_;
    str.getmap = @(mapname) cov_(str, mapname);
    str.disp = @() print_names_(str, maps, mapcovs);
    str.boxplot = @(varargin) boxplot_(str, varargin);
    str.mapcorr = @(map1, map2) mapcorr_(str, map1, map2);
    
    save(str.Output, 'str');
    if nargout > 0, varargout{1} = str; end
end

% Functions for output structure-------------------------------------------
function cov = cov_(str, mapname)
% Returns the similarity map for mapname. Mapname can also be 'BCov',
% 'SCov', or 'BB' if these were specified in the analysis.

    names = cellstr(char(str.Rsms(1,:).name));
    if ~ischar(mapname) || isempty(mapname)
        disp('Required input: one of');
        cellfun(@disp, names);
        cov = [];
        return;
    end
    idx = strcmp(mapname, names);
    if any(idx)
        idx = find(idx);
        cov = zeros(size(str.Rsms(1,idx).RSM));
        for i_ = 1 : size(str.Rsms, 1)
            cov = cov + str.Rsms(i_,idx).RSM;
        end
        cov = cov ./ size(str.Rsms, 1);
    else
        disp(['No map named ', mapname]);
        disp('Valid map names are');
        cellfun(@disp, names);
        cov = [];
        return;
    end
end

function boxplot_(str, idx)
% Draws a boxplot of the average volume correlations of the similarity
% maps. idx indexes the maps one wants to plot (optional)

    if size(str.Mns, 1) < 2
        disp('Only one dataset, cannot draw boxplot');
        return;
    end

    figure;
    if isempty(idx), 
        Mns = str.Mns;  MapNames = str.MapNames;
    else
        Mns = str.Mns(:,idx{1});  MapNames = str.MapNames(idx{1});
    end
    boxplot(Mns);
    ha = gca;
    h = findobj(ha, 'Color', [0 0 1]);
    set(h, 'Color', [0 0 0.6], 'LineWidth', 1.5);
    hold on;
    h = plot(kron(ones(size(Mns,1),1), get(ha, 'XTick')), Mns, 'bo');
    set(h, 'Color', [0.6 0.6 1], 'LineWidth', 0.5);
    plot(get(ha, 'XLim'), [0 0], 'k:')
    set(ha, 'XTickLabel', MapNames, 'FontSize', 12);
    xlabel('');
    ylabel('mean volume correlation'); 
    set(gcf, 'Position', [850   540   400   365]);
    hold off;
end

function r = mapcorr_(str, mapname1, mapname2)
% Computes the concordance of mapname1 and mapname2.

    args = str.Args;
    args.maxslsize = Inf;
    map1 = str.getmap(mapname1);
    if isempty(map1), return; end
    map2 = str.getmap(mapname2);
    if isempty(map2), return; end
    
    rsa = rsa_algorithms();
    ccond = size(str.Rsms(1).RSM, 1);
    args.maxslsize = ccond;
    utils = rsa_utils();
    trmsk = utils.triangmsk(ccond, args.OffDiagOffset);
    f = rsa.corrfun(map2(trmsk), [], trmsk, args);
    [L, D] = svd(map1, 0);
    r = f(L * sqrt(D), []);
end

function print_names_(str, maps, mapcovs)
% Prints information about the specification of the rsa on the console.

    disp('*Modelled similarity maps:')
    cellfun(@disp, maps);
    if ~isempty(mapcovs)
        disp('*Confound maps:')
        cellfun(@disp, mapcovs);
    end
    if str.Args.OffDiagOffset > 0
        fprintf('*Off diagonal offset: %i\n', str.Args.OffDiagOffset);
    end
    fprintf('*Method: %s\n*Searchlight: %s %i\n', str.Args.Method, ...
        str.Args.SearchlightDef, str.Args.SearchlightSize);
    fprintf('*Brain map type: %s\n', str.Args.BrainMapType);
    if ~isempty(str.Args.MaskConfound)
        disp(['Mask confound maps: ', str.Args.MaskConfound]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = srchlght_(istr, args)
% Compute the searchlight correlations for the input specified by the istr
% structure. The correlations are saved to disk. Returns a struct str
% containing the used representational maps in the model and the average
% correlations over the whole volume.

    str = [];
    rsa = rsa_rsmio();
    defs = rsa_defaults();
    utils = rsa_utils();

    %Select desired RDMs for model
    names = cellstr(char(istr.rsm(:).name));
    if isempty(args.MapSel), args.MapSel = names; end
    rsm_idx = ismember(names, args.MapSel);
    if ~any(rsm_idx)
        disp('Invalid RSM selection indices');
        return;
    end
    RSMs = istr.rsm(rsm_idx);
    
    %Select, if given, RSMs for partialling out in correlation
    RSMcovs = rsa.create_rsm();
    if isfield(args, 'MapSelPcorr') && ~isempty(args.MapSelPcorr)
        rsm_idx = ismember(names, args.MapSelPcorr);
        if any(rsm_idx)
            RSMcovs = istr.rsm(rsm_idx);
        end
    end
    
    %Identify beta images
    try
        SPM.xY.VY = spm_vol(istr.files);
    catch LE
        disp(['Could not load file headers from ', istr.directory]);
        disp(LE.message);
        return;
    end
    [flag, msg] = spm_check_orientations(SPM.xY.VY);
    if ~flag
        disp(msg);
        return;
    end
    
    %Effective searchlight size checks.
    voxsize = utils.voxsize(SPM.xY.VY(1));
    slsize = args.SearchlightSize / mean(voxsize);
    if strcmpi(args.SearchlightDef, 'box')
        if slsize < 2
            fprintf(['\nSpecified searchlight size is %1.1f, but voxel ',...
                'size is %1.1f.\n'], args.SearchlightSize, slsize);
            disp('This size is too small for a box searchlight');
            return;
        else
            maxslsize = 27;
        end
    elseif strcmpi(args.SearchlightDef, 'sphere')
        if slsize < 1
            fprintf(['\nSpecified searchlight size is %1.1f, but voxel ',...
                'size is %1.1f.\n'], args.SearchlightSize, slsize);
            disp('This size is too small for a sphere searchlight');
            return;
        else
            if slsize < 1.5,   maxslsize = 7;
            elseif slsize >= 2, maxslsize = 27;
            else maxslsize = 14;
            end
        end
    end
    %the naming of this variable is misleading. It is used to compute the
    %minimal searchlight size relative to what the maximal size can be
    args.maxslsize = maxslsize;  %the full searchlight size
                
    %mask
    try
        SPM.VM = spm_vol(istr.mskfile);
    catch LE
        disp(['Could not load mask header: ', istr.mskfile]);
        disp(LE.message);
        return;
    end
    [flag, msg] = spm_check_orientations([SPM.xY.VY(1), SPM.VM]);
    if ~flag
        disp('Mask is not in the same space as the input volumes');
        disp(msg);
        return;
    end
    
    %additional inputs to test the effect of Bcov & co. At present, the
    %input spec is checked to prevent selecting Bcov when the input is not
    %given by beta images.
    Bcov = [];
    if any(strcmp({'BCov'}, [args.MapSel, args.MapSelPcorr]))
        %we recompute Bcov to improve on rounding, instead of relying on
        %the Bcov stored in the SPM struct
        try
            spmstr = load(istr.spmfile);
            pKX = spmstr.SPM.xX.pKX;
        catch LE
            disp('Could not parse spm.mat file: wrong/damaged file?');
            disp(LE.message);
            return;
        end
        Bcov = pKX * pKX';
        if isempty(Bcov)
            disp('Invalid spm.mat file: model yet to be estimated?');
            return;
        end
        clear pKX spmstr
        % Select columns/rows of Bcov according to the loaded betas.
        if ~strcmp(istr.seltype, 'regexp')
            Bcov = Bcov(istr.idx, istr.idx);
        else
            Bcov = [];%Bcov cannot be selected in args if seltype is regexp
        end
    end
    betas = [];
    if any(ismember({'BB', 'SCov'}, [args.MapSel, args.MapSelPcorr]))
        mskhdr = [];
        if ~isempty(istr.mskconf)
            try
                mskhdr = spm_vol(istr.mskconf);
                spm_check_orientations([mskhdr, SPM.VM]);
            catch LE
                msg{1} = 'The mask for BB/SCov could not be loaded.';
                msg{2} = LE.message;
                msg{3} = 'Choose ''continue'' to ignore mask';
                if utils.input_ui('Critical issue', msg, 'b', ...
                        'Continue|Abort', [0 1], 2)
                    error('Operation terminated by user');
                else
                    mskhdr = [];
                end
            end
        end
        betas = utils.load_vols(SPM.xY.VY, [SPM.VM, mskhdr]);
        betas(:,isnan(sum(betas))) = [];
    end
    corflag = strcmpi(args.BrainMapType, 'cor') || ...
        strcmpi(args.Method, 'regression');
    for m = {'BCov', 'BB', 'SCov'}
        if any(strcmp(m, args.MapSel))
            RSMs(end+1) = rsa.create_rsm(char(m), ...
                getcov_(char(m), betas, Bcov, corflag));
        end
    end
    if isfield(args, 'MapSelPcorr') && ~isempty(args.MapSelPcorr)
        for m = {'BCov', 'BB', 'SCov'}
            if any(strcmp(m, args.MapSelPcorr))
                RSMcovs(end+1) = rsa.create_rsm(char(m), ...
                    getcov_(char(m), betas, Bcov, corflag));
            end
        end
    end
    
    %inputs
    ccond = size(istr.rsm(1).RSM, 1);
    trmsk = utils.triangmsk(ccond, args.OffDiagOffset);
    if sum(trmsk(:)) < 10
        error('The off-diagonal area is too small: %i', sum(trmsk(:)));
    end
    trmsk = logical(trmsk(:));

    covv(:,1) = RSMs(1).RSM(trmsk);
    for i = 2 : numel(RSMs)
        covv(:,end + 1) = RSMs(i).RSM(trmsk);
    end
    for i = 1 : size(covv, 2)  %check there are different similarities
        if 1 == length(unique(covv(:,i)))
            error(['Off-diagonal elements of map ', RSMs(i).name, ...
                ' are identical']);
        end
        if max(abs(covv(:,i))) < defs.offdiagtol || ...
                std(covv(:,i)) / abs(mean(covv(:,i))) < defs.offdiagtol
            msg{1} = ['Off-diagonal elements of map ', RSMs(i).name];
            msg{2} = 'appear to carry no information (low or similar';
            msg{3} = 'values). Recommended action: abort';
            msg{4} = 'or change offdiagtol in rsa_defaults to suppress';
            msg{5} = 'this warning.';
            if utils.input_ui('Critical issue', msg, 'b', ...
                    'Continue|Abort', [0 1], 2)
                error('Operation terminated by user');
            end
            clear msg;
        end
    end
    ccovv = size(covv, 2);
    
    pcovv = [];
    if ~isempty(RSMcovs)
        pcovv(:,1) = RSMcovs(1).RSM(trmsk);
        for i = 2 : numel(RSMcovs)
            pcovv(:,end + 1) = RSMcovs(i).RSM(trmsk);
        end
    end
    for i = 1 : size(pcovv, 2)  %check there are different similarities
        if 1 == length(unique(pcovv(:,i)))
            error(['Off-diagonal elements of map ', RSMcovs(i).name, ...
                ' are identical']);
        end
        if max(abs(pcovv(:,i))) < defs.offdiagtol
            msg{1} = ['Off-diagonal elements of map ', RSMcovs(i).name];
            msg{2} = 'appear to carry no information (low values).';
            msg{3} = 'Recommended action: abort';
            msg{4} = 'or change offdiagtol in rsa_defaults to suppress';
            msg{5} = 'this warning.';
            if utils.input_ui('Critical issue', msg, 'b', ...
                    'Continue|Abort', [0 1], 2)
                error('Operation terminated by user');
            end
            clear msg;
        end
        %you could scale confounds for stability here. Instead, at present
        %a diagnostic is provided below. The idea is that if the scales of
        %the maps are wildly different, something is wrong. If the said
        %diagnostic is switched off, scales are indeed scaled.
        if defs.predscaletol < 0
            pcovv(:,i) = pcovv(:,i) ./ max(abs(pcovv(:,i)));
        end
    end
    
    S = svd([covv, pcovv]);
    if defs.ranktol > 0, tol = defs.ranktol;
    else, tol = size(covv,1) * eps(max(S));
    end
    if sum(S > tol) < (size(covv,2)+size(pcovv,2))
        msg{1} = 'The selected maps are largely collinear';
        msg{2} = 'Recommended action: abort';
        msg{3} = 'or change ranktol in rsa_defaults to suppress this warning';
        if utils.input_ui('Critical issue', msg, 'b', ...
                'Continue|Abort', [0 1], 2)
            error('Operation terminated by user');
        end
        clear msg;
    end
    if defs.predscaletol > 0 && max(S) / min(S) > defs.predscaletol
        msg{1} = 'The maps are on very different scales (very large values';
        msg{2} = 'in beta images/source files?)';
        msg{3} = 'Recommended action: rescale input data if appropriate,';
        msg{4} = 'or change predscaletol in rsa_defaults to suppress this';
        msg{5} = 'warning';
        msg{6} = '(This warning may also arise with FAST when autoregression';
        msg{7} = 'coefficients are high, creating large covariance estimates)';
        msg{8} = 'Absolute peak values in maps:';
        extrema = max(abs([covv, pcovv]), [], 1);
        mapps = [args.MapSel, args.MapSelPcorr];
        for i = 1 : numel(extrema)
            msg{end+1} = sprintf('%s:   %2.3f', mapps{i}, extrema(i));
        end
        if utils.input_ui('Critical issue', msg, 'b', ...
                'Continue|Abort', [0 1], 2)
            error('Operation terminated by user');
        end
        clear msg extrema mapps;
    end
    clear S tol;
    
    %fill-in opts struct for spm's searchlight function.
    opts.def = args.SearchlightDef;
    opts.spec = args.SearchlightSize;
    opts.model = get_modelname_(args);
    opts.maps = cellstr(char(RSMs(:).name));
    if strcmpi(args.Method, 'regression')
        opts.maps = [opts.maps; cellstr(char(RSMcovs(:).name))];
    end
    opts.outputdir = istr.outputdir;

    %searchlight computation
    rsaalgo = rsa_algorithms();
    corrf = rsaalgo.corrfun(covv, pcovv, trmsk, args);
    tic
    R = spm12_searchlight(SPM, opts, corrf);
    toc
    if ~iscell(R), P{1} = R; R = P; clear P; end
    if isempty(R)
        disp('Analysis could not be computed'); 
        disp('Invalid/anomalous maps or map types?');
        return;
    end
    
    %mean rsa correlation values over the volume
    msk = ~isnan(R{1});
    mns = [];
    disp('Volume avg rsa correlations');
    for i = 1 : numel(R)
        mns(i) = mean(R{i}(msk));
        fprintf('%s\t%0.4f\n', opts.maps{i}, mns(i));
    end
    disp(' ');
    
    str.mns = mns;
    str.rsms = RSMs;
    str.rsmcovs = RSMcovs;
    str.corrfun = corrf;
    
    pred = [covv, pcovv, ones(size(covv,1),1)];
    str.cmx = utils.cov2cor(pred' * pred);
end

function cov = getcov_(spec, betas, Bcov, corrflag)
% Create the confound maps BB, SCov, and BCov. corrflag transforms the
% confound maps into correlation matrices.

    switch spec
        case 'BB'
            cov = betas * betas' ./ size(betas, 2);
        case 'SCov'
            bb = bsxfun(@minus, betas, mean(betas));
            cov = bb * bb' ./ size(bb, 2);
        case 'BCov'
            cov = Bcov;
        otherwise
            error(spec);
    end
    cov = (cov + cov') ./ 2;
    if corrflag
        utils = rsa_utils();
        cov = utils.cov2cor(cov); 
    end
end

function model = get_modelname_(args)
% Form a reasonable model name to identify the model. The name is based on
% args.ModelName, if nonempty; otherwise a name based on the selected
% features will be formed automatically.

    if ~isfield(args, 'ModelName'), args.ModelName = []; end
    model = args.ModelName;
    if strcmpi(args.Method, 'slsize')
        if isempty(model)
            model = sprintf('%s%i-slsize', args.SearchlightDef, ...
                args.SearchlightSize);
        else
            model = sprintf('%s%i-%s-slsize', args.SearchlightDef, ...
                args.SearchlightSize, model);
        end
        return;
    end
    if isempty(model)
        %do not encode maps, as they are expressed in the individual files
        %anyway
        %model = strjoin(encode_(args.MapSel), '');
        if ~isempty(args.MapSelPcorr)
            model = [model, 'x', ...
                strjoin(encode_(args.MapSelPcorr), 'x')];
        end
    end
    if args.OffDiagOffset > 0
        model = ['ODO', int2str(args.OffDiagOffset), '-', model];
    end
    if strcmpi(args.Method, 'regression')
        model = ['lm-', model];
        args.BrainMapType = '';  %avoid recording this in name
    elseif ~strcmpi(args.Method, 'Pearson')
        model = [char(encode_({args.Method})), '-', model];
    end
    model = sprintf('%s%i-%s-%s', args.SearchlightDef, ...
        args.SearchlightSize, args.BrainMapType, model);
    model = strrep(model, '--', '-');
    if model(end) == '-', model = model(1:end-1); end
    
    function codes = encode_(strcell)
        codes = cellfun(@(x) regexprep(x, '[aeiouAEIOU]', ''), strcell, ...
            'UniformOutput', false);
        for i_ = 1 : length(codes)
            if length(codes{i_}) > 3, codes{i_} = codes{i_}(1:3); end
            if length(codes{i_}) > 2
                codes{i_} = [upper(codes{i_}(1)), codes{i_}(2:end)];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utilities
function istrout = load_input_(args)    
% Load and validate the input for all subjects. Returns empty for invalid
% input. Diagnostic directly printed on console.

    istrout = [];
    rsa = rsa_rsmio();
    utils = rsa_utils();
    
    dirnames = args.Directories;
    cmnpth = utils.get_commonpath(dirnames);
    csubj = size(dirnames, 1);
    
    %initialize istr struct array.
    istr = cellfun(@(x) struct(...
            'directory', x,...
            'rsm', '', ...
            'files', '', ...
            'seltype', '', ...
            'idx', [], ...
            'mskfile', '', ...
            'mskconf', '', ...
            'spmfile', fullfile(x, 'spm.mat'), ...
            'outputdir', ''), cellstr(dirnames));
    
    %rsm files
    flag = false;
    if any(strfind(args.MapFile, filesep))
        rsm = rsa.load_rsm(args.MapFile);
        if isempty(rsm), flag = true;
        else, [istr.rsm] = deal(rsm); 
        end
    else
        for i = 1 : csubj
            rsm = rsa.load_rsm(fullfile(deblank(dirnames(i,:)), ...
                args.MapFile));
            if isempty(rsm), flag = true; continue; end
            names{i} = cellstr(char(rsm(:).name))'; %#ok<*AGROW>
            istr(i).rsm = rsm(ismember(names{i}, [args.MapSel, args.MapSelPcorr]));
        end
        idx = cellfun(@(x) ~isempty(setxor(names{1}, x)), names(2:end));
        if any(idx)
            disp('The rsm files do not contain the same maps:');
            for i = (find(idx)' + 1), disp(deblank(dirnames(i,:))); end
            flag = true;
        end
        clear names
    end
    if flag, return; end
    names = cellstr(char(istr(1).rsm(:).name))'; 
    idx = ~ismember(args.MapSel, [names, {'BCov', 'BB', 'SCov'}]);
    if any(idx)
        disp('Specified maps are not in the rsm files:');
        cellfun(@disp, args.MapSel(idx));
        disp(' ');
        disp('Available maps in the rsm files:');
        cellfun(@disp, names);
        flag = true;
    end
    if ~isempty(args.MapSelPcorr)
        idx = ~ismember(args.MapSelPcorr, [names, {'BCov', 'BB', 'SCov'}]);
        if any(idx)
            disp('Specified maps for partial correlation are not in the map files:');
            cellfun(@disp, args.MapSelPcorr(idx));
            disp(' ');
            disp('Available maps for partialling out:');
            cellfun(@disp, setdiff([names, {'BCov', 'BB', 'SCov'}], ...
                args.MapSel));
            flag = true;
        end
    end
    if flag, return; end
    
    %selected maps have same numer of dimensions. If idx specified, will be
    %validated to select this number of beta images.
    idx = ismember(names, union(args.MapSel, args.MapSelPcorr));
    rsm = istr(1).rsm(idx);
    ccond = size(rsm(1).RSM, 1);
    for i = 1 : csubj
        idxx(i) = any(arrayfun(@(x) ccond ~= size(x.RSM, 1), ...
            istr(i).rsm(idx)));
    end
    if all(idxx)
        disp('Not all selected maps have the same number of dimensions.');
        disp('Dimensions of selected maps in subject 1:');
        arrayfun(@(x) fprintf('%s\t%i\n', x.name, size(x.RSM, 1)), ...
            istr(1).rsm(idx));
        return;
    elseif any(idxx)
        disp('Selected maps have differing dimensions in these subjects:')
        for i = find(idxx), disp(dirnames(i,:)); end
    end
    clear idxx
    
    % Beta images for input
    function [files, seltype, b_idx] = files_input_(dirname, beta_idx)
        files = []; seltype = ''; b_idx = [];
        if ~isempty(args.BetaFiles) && iscell(args.BetaFiles)
            %input files specified by regular expression
            seltype = 'regexp'; b_idx = [];
            files = spm_select('FPList', dirname, args.BetaFiles{1});
            if isempty(files)
                disp(['No files found in ', dirname, ' using ' ...
                    args.BetaFiles{1}]);
                return;
            end
            files = sortrows(files);
            if isempty(beta_idx),  beta_idx = 1:ccond; end
            beta_idx = utils.validate_idx(beta_idx);
            if isempty(beta_idx), return; end
            try
                files = files(beta_idx,:);
            catch AE
                disp('Invalid BetaIdx specification:');
                disp(AE.message);
                return;
            end
            if size(files, 1) ~= ccond
                disp(['Invalid BetaIdx field value: wrong number of ', ...
                    'selected files']);
                disp(['BetaIdx value directs selection of ', ...
                    int2str(size(files,1)), ' files in ', dirname]);
                disp(['but there are ', int2str(ccond), ' dimensions in maps']);
                files = []; return;
            end
            disp(' ');
            disp(['Selected files in ', dirname, ':']);
            cellfun(@disp, cellstr(files));
            b_idx = beta_idx;
        else
            %input files are beta images in input directories
            prefix = args.BetaFiles;
            if isempty(beta_idx)
                %here, we assume that the RDM dimension corresponds to the
                %first beta images we load, thus excluding realign pars etc.
                seltype = 'beta'; b_idx = 1:ccond;
            else
                %BetaIdx is an array indicating which images to load
                seltype = 'idx';
                b_idx = utils.validate_idx(beta_idx);
                if isempty(b_idx), return; end
                if ccond ~= length(b_idx)
                    disp(['Invalid BetaIdx field value: wrong number ', ...
                        'of selected files']);
                    disp(['BetaIdx value directs selection of ', ...
                        int2str(length(b_idx)), ' files in ', dirname]);
                    disp(['but there are ', int2str(ccond), ...
                        ' dimensions in maps']);
                    return;
                end
            end
            try
                files = arrayfun(@(x) fullfile(dirname, ...
                    sprintf('%sbeta_%04i.nii', prefix, x)), b_idx, ...
                    'UniformOutput', false);
            catch AE
                disp(dirname);
                disp('Invalid BetaIdx field argument');
                disp(AE.message);
                return;
            end
            idx_ = cellfun(@(x) ~exist(x, 'file'), files);
            if any(idx_)
                disp('Following files are required as input');
                cellfun(@disp, files(idx_));
                disp(['but were not found in ', dirname]);
                files = []; return;
            end
            files = char(files);
            if strcmp(seltype, 'idx')
                disp(' ');
                disp(['Selected files in ', dirname, ':']);
                cellfun(@disp, cellstr(files));
            end                
        end
    end  %function
    flag = false;
    if any(arrayfun(@(x) ~isempty(x.idx), istr(1).rsm)) && ~isempty(args.BetaIdx)
        disp('Error: beta idx spec in arguments and in rsm structure');
        return;
    end
    for i = 1 : csubj
        if numel(istr(i).rsm) > 1
            %this will work also if rsm.idx is empty
            if ~all(arrayfun(@(x) all(eq(istr(i).rsm(1).idx, x.idx)), istr(i).rsm))
                disp('The idx fields of maps are not the same across maps');
                return;
            end
        end
        if ~isempty(istr(i).rsm(1).idx), beta_idx = istr(i).rsm(1).idx;
        else beta_idx = args.BetaIdx; end
        [files, seltype, bidx] = files_input_(istr(i).directory, beta_idx);
        if isempty(files), flag = true;
        else
            istr(i).files = files;
            istr(i).seltype = seltype;
            istr(i).idx = bidx;
            if ~isempty(args.OutputDir)
                istr(i).outputdir = strrep(istr(i).directory, cmnpth, ...
                    args.OutputDir);
            else
                istr(i).outputdir = istr(i).directory;
            end
        end
    end
    if flag, return; end
    
    % masks
    if isempty(args.MaskFile)
        %Seek mask.nii in each input directory
        idx = cellfun(@(x) ~exist(fullfile(x, 'mask.nii'), 'file'), ...
            cellstr(dirnames));
        if all(idx)
            disp('Could not find mask files in input directories.');
            return;
        elseif any(idx)
            disp('Mask files missing in the following input directories:');
            for d = find(idx)', disp(dirnames(d,:)); end
            return;
        end
        for i = 1 : csubj
            istr(i).mskfile = fullfile(deblank(dirnames(i,:)), 'mask.nii');
        end
    else
        if any(strfind(args.MaskFile, filesep))
            if ~exist(args.MaskFile, 'file')
                disp(['Specified mask ''', args.MaskFile]);
                disp('does not exist');
                return;
            end
            [istr.mskfile] = deal(args.MaskFile);
        else
            for i = 1 : csubj
                istr(i).mskfile = fullfile(deblank(dirnames(i,:)), ...
                    args.MaskFile);
            end
            idx = arrayfun(@(x) ~exist(x.mskfile, 'file'), istr);
            if all(idx)
                disp('Could not find mask files in input directories.');
                return;
            elseif any(idx)
                disp('Mask files missing in the following input directories:');
                for d = find(idx)', disp(dirnames(d,:)); end
                return;
            end
        end
    end
    
    if any(ismember({'BB', 'SCov'}, [args.MapSel, args.MapSelPcorr])) && ...
            (isfield(args, 'MaskConfound') && ~isempty(args.MaskConfound))
        if any(strfind(args.MaskConfound, filesep))
            if ~exist(args.MaskConfound, 'file')
                disp(['Specified mask ''', args.MaskConfound]);
                disp('does not exist');
                return;
            end
            [istr.mskconf] = deal(args.MaskConfound);
        else
            idx = cellfun(@(x) ~exist(fullfile(x, args.MaskConfound), ...
                'file'), cellstr(dirnames));
            if all(idx)
                disp('Could not find confound mask files in input directories.');
                return;
            elseif any(idx)
                disp('Mask file missing in the following input directories:');
                for d = find(idx)', disp(dirnames(d,:)); end
                return;
            end
            for i = 1 : csubj
                istr(i).mskconf = fullfile(deblank(dirnames(i,:)), ...
                    args.MaskConfound);
            end
        end
    end
    
    % using BCov as confound map requires spm.mat
    if ismember({'BCov'}, [args.MapSel, args.MapSelPcorr])
        spmmats = cellfun(@(x) fullfile(x, 'spm.mat'), ...
            cellstr(dirnames), 'UniformOutput', false);
        idx = cellfun(@(x) ~exist(x, 'file'), spmmats);
        if any(idx)
            if all(idx)
                disp('Directories do not contain spm.mat, but');
                disp('spm.mat is required to read BCov');
            else
                disp('The following directories did not contain spm.mat:');
                for i = find(idx)'
                    disp(deblank(dirnames(i,:)));
                end
                disp('however, spm.mat is required to read BCov, as');
                disp('specified in other directive');
            end
            return;
        end
    end
    
    if ~isempty(args.OutputDir)
        if ~exist(args.OutputDir, 'dir')
            mkdir(args.OutputDir);
        end
    end
    
    %empty return for abort
    istrout = istr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input interface
function args = Ui_(args, iname)
    figure(spm_figure('GetWin', 'Interactive'));
    spm_progress_bar('Clear');
    spm('Pointer');
       
    rsa = rsa_rsmio();
    utils = rsa_utils();
    %utility to print out calling var name
    fdisp = @(s) fprintf([strrep(s, '\', '\\'), '\n'], iname);
    
    if ~isfield(args, 'Directories') || isempty(args.Directories)
        args.Directories = spm_select(Inf, 'dir', ...
            'Select first-level directories');
    end
    if isempty(args.Directories), args.Quit = 1; return; end
    if ~ischar(args.Directories)
        fdisp('Invalid argument in %s.Directories: must be a string');
        args.Quit = 1;
        return;
    end
    idx = cellfun(@(x) ~exist(x, 'dir'), cellstr(args.Directories));
    if any(idx)
        args.Quit = 1;
        if all(idx), fdisp('%s.Directories: these directories do not exist.');
        else
            fdisp('%s.Directories: The following directories could not be found:');
            for i = find(idx)', disp(deblank(args.Directories(i,:))); end
        end
        return;
    end
    
    %Method
    if ~isfield(args, 'Method') || isempty(args.Method)
        args.Method = 'Pearson';
    end
    if ~ischar(args.Method)
        fdisp('Invalid argument in %s.Method: must be a string');
        disp('Method must be one of ''Pearson'', ''Spearman'', ''regression''');
        args.Quit = 1;
        return;
    end
    if ~ismember(args.Method, {'Pearson', 'Spearman', 'regression', ...
            'Regression', 'slsize'})
        fdisp(['%s.Method: Invalid method: ', args.Method]);
        disp('Method must be one of ''Pearson'', ''Spearman'', ''regression''');
        args.Quit = 1;
        return;
    end
    
    %Type of similarity from bran signal
    if ~isfield(args, 'BrainMapType') || isempty(args.BrainMapType)
        %args.BrainMapType = 'sscp';
        msg{1} = 'Specify similarity map for brain signal.';
        msg{2} = 'sscp: sum of squares and cross products;';
        msg{3} = 'covariance: mean-adjusted sscp;';
        msg{4} = 'correlation: variance-adjusted covariance.';
        msg{5} = '''covariance'' and ''correlation'' eliminate';
        msg{6} = 'effect of mean searchlight activity.';
        args.BrainMapType = utils.input_ui('Brain similarity map:', ...
            msg, 'm', 'sscp|covariance|correlation', ...
            char('sscp', 'cov', 'cor'), 1);
        clear msg;
    end
    if ~ischar(args.BrainMapType)
        fdisp('Invalid argument in %s.BrainMapType: must be a string')
        disp('with content ''sscp'', ''cov'', or ''cor''');
        args.Quit = 1;
        return;
    end
    if size(args.BrainMapType, 1) ~= 1
        fdisp('%s.BrainMapType must be a single string,');
        disp('''sscp'', ''cov'', or ''cor''');
        args.Quit = 1;
        return;
    end
    args.BrainMapType = lower(args.BrainMapType);
    distfs = {'euclidean', 'squaredeuclidean', 'seuclidean', 'cosine', ...
        'correlation', 'spearman', 'mds'};
    if ~any(strcmp(args.BrainMapType, [{'sscp', 'cov', 'cor'}, distfs]))
        fdisp(['%s.BrainMapType: invalid value: ', args.BrainMapType]);
        disp('Possible methods: ''sscp'', ''cov'', ''cor''');
        args.Quit = 1;
        return;
    end
    
    %Representational similarity maps
    if ~isfield(args, 'MapFile') || isempty(args.MapFile)
        msg{1} = 'Map files contain the representational similarity maps.';
        msg{2} = 'One map file can be selected if the maps are the same';
        msg{3} = 'in the whole experiment, or multiple map files can be';
        msg{4} = 'loaded from the input directories.';
        flag = utils.input_ui('Which map file?', msg, 'm', ...
            'One file for all data|One file in each dir', [0, 1]);
        clear msg;
        if flag
            args.MapFile = utils.input_ui('Name of map file:', '+1', 's');
        else
            args.MapFile = spm_select(1, '.*mat$', 'Select file');
        end
    end
    if ~ischar(args.MapFile)
        fdisp('Invalid argument in %s.MapFile: must be a string');
        args.Quit = 1;
        return;
    end
    if any(strfind(args.MapFile, filesep))
        if ~exist(args.MapFile, 'file')
            fdisp(['%s.MapFile: File ''', args.MapFile, '''does not exist.']);
            args.Quit = 1;
            return;
        end
    else
        idx = cellfun(@(x) ~exist(fullfile(x, args.MapFile), 'file'), ...
            cellstr(args.Directories));
        if any(idx)
            args.Quit = 1;
            if all(idx)
                fdisp(['%s.MapFile: Files ''', args.MapFile, ...
                    ''' not found in the input directories.']);
            else
                fdisp('%s.MapFile: These map files do not exist: ');
                for i = find(idx)'
                    disp(fullfile(deblank(args.Directories(i,:)), ...
                        args.MapFile));
                end
            end
            return;
        end
    end
    
    %We load a map file here to get the names of available maps, which will
    %be used in several places below.
    %Load a map file.
    if any(strfind(args.MapFile, filesep)), rsm = args.MapFile;
    else
        rsm = fullfile(deblank(args.Directories(1,:)), args.MapFile);
    end
    rsm = rsa.load_rsm(rsm);
    if isempty(rsm)
        args.Quit = 1;
        return;
    end
    %Get names
    mapnames = arrayfun(@(x) x.name, rsm, 'UniformOutput', false);
    
    %Map selection
    DONE_CODE = 'DONE (exit menu)';
    CONF_CODE = 'Confound maps...';
    if ~isfield(args, 'MapSel')
        %Input
        names = mapnames;
        msg{1} = 'Choose the maps for the analysis. At least one map';
        msg{2} = ['must be selected. When finished, choose ''', DONE_CODE, ''''];
        msg{3} = 'to exit.';
        args.MapSel(1) = utils.input_ui('Choose input maps:', msg, 'm', ...
            names, names); clear msg;
        %revised: do not offer confounds at this stage in the interface
        %menu (manual selection always allowed)
        %names = [names, {CONF_CODE, DONE_CODE}];
        names = [names, {DONE_CODE}];
        while true
            names(ismember(names, args.MapSel)) = [];
            if numel(names) < 2, break; end  %only 'DONE' left
            sel = utils.input_ui('Choose input maps:', '+1', 'm', names, names);
            if strcmp(char(sel), DONE_CODE), break;
            elseif strcmp(char(sel), CONF_CODE)
                names(strcmp(CONF_CODE, names)) = [];
                names = [names(1:end-1), {'BCov', 'SCov', 'BB', DONE_CODE}];
                continue;
            else, args.MapSel(end+1) = sel; 
            end
        end
    end
    if isempty(args.MapSel) %include all names
        if any(strfind(args.MapFile, filesep)), rsm = args.MapFile;
        else
            rsm = fullfile(deblank(args.Directories(1,:)), args.MapFile);
        end
        rsm = rsa.load_rsm(rsm);
        if isempty(rsm)
            args.Quit = 1;
            return;
        end
        %Get all names
        args.MapSel = mapnames;
    end
    if ~iscell(args.MapSel)
        fdisp('Invalid argument in %s.MapSel: must be a cell array of strings');
        args.Quit = 1;
        return;
    end
    if any(~cellfun(@ischar, args.MapSel))
        fdisp(['Invalid argument in %s.MapSel: the cell array must ', ...
            'containg strings']);
        args.Quit = 1;
        return;
    end
    if all(ismember([mapnames, {'BCov', 'SCov', 'BB'}], args.MapSel))
        if isfield(args, 'MapSelPcorr') && ~isempty(args.MapSelPcorr)
            fdisp('Invalid %s.MapSelPcorr: ignored');
        end
        args.MapSelPcorr = {};
    end
    
    %give info about options for adjustment
    msg{1}='Select the options to ajdust for non-orthogonality of design.';
    msg{2}='''partial corr. maps'' allows selecting one of BCov, Scov, or';
    msg{3}='BB (see the manual) to partial out the effect of design.';
    msg{4}='''diagonal offset'' leaves out from the correlation elements';
    msg{5}='of the representational matrix close to the diagonal.';
    msg{6}='Importantly, the diagonal offset only makes sense if the ';
    msg{7}='columns of the design matrix are individual trials arranged';
    msg{8}='in temporal order.';
    msg{9}='The combination of both is meaningful if the TR is short and/or';
    msg{10}='large autocorrelation models cannot be otherwise corrected.';
    if ~any(isfield(args, {'MapSelPcorr', 'OffDiagOffset'}))
        flag = utils.input_ui('Bias adjustment?', msg, 'm', ...
            'partial corr. maps|diagonal offset|both|none', ...
            [1 2 3 0]);  clear msg;
        switch flag
            case 1
                args.OffDiagOffset = 0;
            case 2
                args.MapSelPcorr = [];
            case 0
                args.OffDiagOffset = 0;
                args.MapSelPcorr = [];
        end
    end
    
    if ~isfield(args, 'MapSelPcorr')
        names = [{'BCov', 'SCov', 'BB'}, mapnames, {'NONE (exit menu)'}];
        names(ismember(names, args.MapSel)) = [];
        %Input
        sel = utils.input_ui('Choose confound maps:', '0', 'm', names, names);
        if strcmp(char(sel), 'NONE (exit menu)')
            args.MapSelPcorr = [];
        else
            args.MapSelPcorr(1) = sel;
            names{end} = DONE_CODE;
            while true
                names(ismember(names, args.MapSelPcorr)) = [];
                if any(ismember(args.MapSelPcorr, {'BB', 'SCov'}))
                    names(ismember(names, {'BB', 'SCov'})) = [];
                end
                if numel(names) < 2, break; end
                sel = utils.input_ui('Choose confound maps:', '+1', 'm', ...
                    names, names);
                if strcmp(char(sel), DONE_CODE), break;
                else, args.MapSelPcorr(end+1) = sel; end
            end
        end
    end    
    if isempty(args.MapSelPcorr) && ischar(args.MapSelPcorr)
        args.MapSelPcorr = [];
    end
    if ~isempty(args.MapSelPcorr) && ...
            (~iscell(args.MapSelPcorr) || any(~cellfun(@ischar, args.MapSelPcorr)))
        fdisp('Invalid argument in %s.MapSelPcorr: must be a cell string');
        args.Quit = 1;
        return;
    end
    if ~isempty(intersect(args.MapSel, args.MapSelPcorr))
        disp('Predictor maps cannot be selected twice:');
        for i = intersect(args.MapSel, args.MapSelPcorr), disp(i); end
        fdisp('(note that empty %s.MapSel selects all maps).');
        args.Quit = 1;
        return;
    end
    
    if all(ismember({'BB', 'SCov'}, union(args.MapSel, args.MapSelPcorr)))
        if utils.input_ui('BB and SCov are likely collinear.', '+1', 'b', ...
                'continue?|stop', [0 1], 2)
            args.Quit = 1;
            return;
        end
    end
    
    %Off-diagonal offset
    if ~isfield(args, 'OffDiagOffset')
        args.OffDiagOffset = utils.input_ui('Off diagonal offset?', '+1', ...
            'i');
    elseif isempty(args.OffDiagOffset)
        args.OffDiagOffset = 0;
    end
    if ~isnumeric(args.OffDiagOffset)
        fdisp(['Invalid %s.OffDiagOffset argument: must be zero or ', ...
            'positive integer scalar']);
        args.Quit = 1;
        return;
    end
    if ~isscalar(args.OffDiagOffset) || 0 ~= mod(args.OffDiagOffset, 1)
        fdisp('%s.OffDiagOffset must be an integer scalar');
        args.Quit = 1;
    end
    if args.OffDiagOffset < 0
        fdisp('%s.OffDiagOffset: Off diagonal offset cannot be negative.');
        args.Quit = 1;
        return;
    elseif args.OffDiagOffset > 7
        if utils.input_ui('Unusually large offset.', '+1', 'b', ...
                'continue?|stop', [0 1], 2)
            args = rmfield(args, 'OffDiagOffset');
            args.Quit = 1;
            return;
        end
    end

    %Beta or con images
    if ~isfield(args, 'BetaFiles')
        msg{1} = 'The option ''beta img in directories'' selects the beta';
        msg{2} = 'images produced by SPM. Use ''beta img with prefix'' to';
        msg{3} = 'select preprocessed beta images (for example with an ''s''-';
        msg{4} = 'prefix).';
        if ismember({'BCov'}, union(args.MapSel, args.MapSelPcorr))
            flag = utils.input_ui('Which input volumes?', msg, 'm', ...
                'beta img in directories|beta img with prefix', ...
                [0 1], 1);
        else
            msg{end+1} = 'Choose ''regular expr. files'' to select other inputs';
            msg{end+1} = 'through a regular expression.';
            flag = utils.input_ui('Which input volumes?', msg, 'm', ...
         'beta img in directories|beta img with prefix|regular expr. files', ...
                [0 1 2], 1);
        end
        clear msg;
        switch flag
            case 0
                args.BetaFiles = [];
            case 1
                args.BetaFiles = utils.input_ui('Specify prefix', '+1', 's', 's');
            case 2
                args.BetaFiles = cellstr(utils.input_ui('Specify reg expr.', ...
                    '+1', 's', '^beta_\d{4}\.nii$'));
        end
    end
    if ~isempty(args.BetaFiles)
        if iscell(args.BetaFiles)
            if ~ischar(args.BetaFiles{1})
                fdisp('Not a regular expression in %s.BetaFiles');
                args.Quit = 1;
                return;
            end
            if isempty(args.BetaFiles{1})
                fdisp('%s.BetaFiles: Regular expression expected, empty value given');
                args.Quit = 1;
                return;
            end
            if 1 ~= numel(args.BetaFiles)
                fdisp('Only one regular expression in %s.BetaFiles');
                args.Quit = 1;
                return;
            end
            if ismember({'BCov'}, union(args.MapSel, args.MapSelPcorr))
                fdisp('Regular expression given in %s.BetaFiles');
                disp('while also using BCov as a predictor or confound.');
                disp('Invalid directive: do not know how to compute BCov');
                disp('from spm.mat for arbitrary input files.');
                disp('Use ''BB'' or ''SCov'' instead of ''BCov''.');
                args.Quit = 1;
                return;
            end
        else
            if ~ischar(args.BetaFiles)
                fdisp('%s.BetaFiles should be empty (to select beta images)');
                disp('or contain prefix string to select preprocessed beta images');
                if ~ismember({'BCov'}, union(args.MapSel, args.MapSelPcorr))
                    disp('or contain a regular expression to select volumes');
                end
                args.Quit = 1;
                return;
            end
            if 1 ~= size(args.BetaFiles, 1)
                fdisp('Char matrix in %s.BetaFiles:');
                disp('Are you trying to select input volumes directly?');
                disp('Use a cell with a regular expression to select volumes');
                disp('or one prefix string to select preprocessed beta images.');
                args.Quit = 1;
                return;
            end
        end
    end
    %empty, string, or singleton cell of string
    if ~isempty(args.BetaFiles) && ~ischar(args.BetaFiles) && ...
        ~iscell(args.BetaFiles)
        fdisp('Invalid %s.BetaFiles argument.');
        fdisp('%s.BetaFiles should be empty (to select the beta_00xx.nii');
        disp('files of the model), a prefix string (to select preprocessed');
        disp('volumes), or a cell string (a regular expression to select');
        disp('the input volumes).');
        args.Quit = 1;
        return;
    elseif iscell(args.BetaFiles)
        if 1 ~= numel(args.BetaFiles) || ~ischar(args.BetaFiles{1})
            fdisp('Invalid %s.BetaFiles argument.');
            disp('Cell must contain one string element');
            fdisp('%s.BetaFiles should be empty (to select the beta_00xx.nii');
            disp('files of the model), a prefix string (to select preprocessed');
            disp('volumes), or a cell string (a regular expression to select');
            disp('the input volumes).');
            args.Quit = 1;
            return;
        end
    end
    
    if ~isfield(args, 'BetaIdx')
        args.BetaIdx = [];
    end
    if ~isempty(args.BetaIdx)
        if ~(isnumeric(args.BetaIdx) || islogical(args.BetaIdx))
            fdisp('Invalid %s.BetaIdx value: array of numbers or logical');
            args.Quit = 1;
            return;
        end
        if min(size(args.BetaIdx)) > 1
            fdisp('Invalid %s.BetaIdx value: must be vector, not matrix');
            args.Quit = 1;
            return;
        end
        args.BetaIdx = args.BetaIdx(:)';
    end
    
    %Mask
    if ~isfield(args, 'MaskFile')
        idx = cellfun(@(x) exist(fullfile(x, 'mask.nii'), 'file'), ...
            cellstr(args.Directories));
        if all(idx)
            msg{1} = 'Choosing ''masks in input dirs'' directs the programme';
            msg{2} = 'to use the ''mask.nii'' files produced by SPM as the';
            msg{3} = 'mask of each input directory.'; 
            flag = utils.input_ui('Mask selection', msg, 'm', ...
                'masks in input dirs|mask I select now', [0 1]);
            clear msg;
            if flag
                args.MaskFile = spm_select(1, '.*(nii|img)$', 'Select mask');
            else
                args.MaskFile = [];
            end
        else
            args.MaskFile = spm_select(1, '.*(nii|img)$', 'Select mask');
        end
    else
        if isempty(args.MaskFile)
            idx = cellfun(@(x) ~exist(fullfile(x, 'mask.nii'), 'file'), ...
                cellstr(args.Directories));
            if all(idx)
                fdisp('Empty %s.MaskFile directs selection of ''mask.nii'' files');
                disp('in the input directories. However, there are no ''mask.nii''');
                disp('files in the input directories, or some mask is missing.');
                args.Quit = 1;
                return;
            elseif any(idx)
                fdisp('Empty %s.MaskFile directs selection of ''mask.nii'' files');
                disp('in the input directories. However, some ''mask.nii'' files');
                disp('are missing:');
                cellfun(@disp, cellstr(args.Directories(idx,:)));
            end
        else
            if ~ischar(args.MaskFile)
                fdisp('Invalid argument %s.MaskFile: must be a string');
                args.Quit = 1;
                return;
            end
            if size(args.MaskFile, 1) > 1
                fdisp('Invalid argument %s.MaskFile: only one mask allowed');
                args.Quit = 1;
                return; 
            end
            % check for existence in input checks
        end
    end
    
    if any(ismember({'BB', 'SCov'}, union(args.MapSel, args.MapSelPcorr)))
        msg{1} = 'This is an optional mask that may be specified to';
        msg{2} = 'compute the BB or the SCov terms. The mask will be combined';
        msg{3} = 'with an AND operation with the volume mask.';
        if ~isfield(args, 'MaskConfound') && ...
            utils.input_ui('Specific mask for BB/Scov?', msg, 'y/n', [1 0])
            clear msg;
            flag = utils.input_ui('Mask selection', '+1', 'm', ...
                'masks in input dirs|mask I select now', [0 1]);
            if flag
                args.MaskConfound = spm_select(1, '.*(nii|img)$', ...
                    'Select mask');
            else
                args.MaskConfound = utils.input_ui('Mask name:', '+1', 's');
            end
        else
            args.MaskConfound = [];
        end
        if ~isempty(args.MaskConfound)
            if ~ischar(args.MaskConfound)
                fdisp('Invalid argument %s.MaskConfound: must be a string');
                args.Quit = 1;
                return;
            end
            if size(args.MaskConfound, 1) > 1
                fdisp('Invalid argument %s.MaskConfound: only one mask allowed');
                args.Quit = 1;
                return; 
            end
            if any(strfind(args.MaskConfound, filesep)) && ...
                    ~exist(args.MaskConfound, 'file')
                fdisp(['%s.MaskConfound: ''', args.MaskConfound, ...
                    ''' does not exist']);
                args.Quit = 1;
                return;
            end
        end
    end
        
    %Searchlight definition
    if ~isfield(args, 'SearchlightDef') || isempty(args.SearchlightDef)
        args.SearchlightDef = 'sphere';
    else
        if ~ischar(args.SearchlightDef)
            fdisp('Invalid %s.SearchlightDef argument: must be a string');
            disp('Possible vales: ''sphere'', ''box''');
            args.Quit = 1;
            return;
        end
        if ~ismember(args.SearchlightDef, {'sphere', 'box'})
            disp(['Invalid searchlight: ', args.SearchlightDef]);
            disp('Possible vales: ''sphere'', ''box''');
            args.Quit = 1;
            return;
        end
    end
    if ~isfield(args, 'SearchlightSize')
        args.SearchlightSize = utils.input_ui('Searchlight size (mm):', ...
            '+1', 'n', 8);
    end
    if isempty(args.SearchlightSize), args.SearchlightSize = 8; end
    if ~isnumeric(args.SearchlightSize)
        fdisp('Invalid %s.SearchlightSize argument: must be numeric');
        args.Quit = 1;
        return;
    end
    if ~isscalar(args.SearchlightSize) || 0 ~= mod(args.SearchlightSize, 1)
        fdisp('%s.SearchlightSize must be an integer scalar');
        args.Quit = 1;
    end
    if args.SearchlightSize < 1
        fdisp('Invalid %s.SearchlightSize argument: must be positive');
        args.Quit = 1;
        return;
    end
    if args.SearchlightSize <= 2
        if utils.input_ui('Unusually small searchlight size.', ...
                '+1', 'b', 'continue?|stop', [0 1], 2)
            args = rmfield(args, 'SearchlightSize');
            args.Quit = 1;
            return;
        end
    elseif args.SearchlightSize >= 20
        if utils.input_ui('Unusually large searchlight size.', ...
                '+1', 'b', 'continue?|stop', [0 1], 2)
            args = rmfield(args, 'SearchlightSize');
            args.Quit = 1;
            return;
        end
    end
    
    %Model name
    if ~isfield(args, 'ModelName')
        msg{1} = 'Enter here a model name to identify the output files.';
        msg{2} = 'If you leave the model name to ''AUTO'', a model name';
        msg{3} = 'will be formed automatically. The automatic model name';
        msg{4} = 'for the current directives is';
        msg{5} = strrep(get_modelname_(args), '_', '\_');
        args.ModelName = utils.input_ui('Model name', msg, 's', 'AUTO');
        if strcmp(args.ModelName, 'AUTO'), args.ModelName = []; end
    end
    if ~isempty(args.ModelName) 
        if ~ischar(args.ModelName)
            fdisp('Invalid %s.ModelName argument: must be a string');
            args.Quit = 1;
            return;
        end
        if size(args.ModelName, 1) ~= 1
            fdisp('%s.ModelName must be a single string');
            args.Quit = 1;
            return;
        end
    end
    
    %Output directory
    if ~isfield(args, 'OutputDir')
        args.OutputDir = spm_select([0,1], 'dir', ...
            'Output directory (or none)');
    end
    if ~isempty(args.OutputDir) 
        if ~ischar(args.OutputDir)
            fdisp('Invalid %s.OutputDir argument: must be string');
            args.Quit = 1;
            return;
        end
        if size(args.OutputDir, 1) ~= 1
            fdisp('%s.OutputDir must be a single string');
            args.Quit = 1;
            return;
        end
        if args.OutputDir(end) ~= filesep
            args.OutputDir = [args.OutputDir, filesep];
        end
    end
    
    %Check fields
    flds = {            
       'MFile'          
       'Version'        
       'SPMVersion'    
       'MatlabVersion'
       'Timestamp'      
       'Code'
       'Directories'   
       'BrainMapType'
       'Method'         
       'MapFile'        
       'MapSel'         
       'MapSelPcorr'    
       'BetaFiles'      
       'BetaIdx'       
       'MaskFile'
       'MaskConfound'
       'MDSDims'
       'SearchlightDef' 
       'SearchlightSize'
       'OffDiagOffset'
       'ModelName'
       'OutputDir'
       'Args'
       'Quit'
       }; 
    if any(~ismember(fields(args), flds))
        argsflds = fields(args);
        idx = ~ismember(argsflds, flds);
        fdisp('Invalid arguments: the following fields of %s are unkown');
        for i = find(idx)',     disp(argsflds{i});      end
        args.Quit = 1;
    end
end

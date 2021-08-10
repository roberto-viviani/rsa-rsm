function rsa = rsa_rsmio()
% A library of i/o functions for the rsa add-on functions. This library
% centralizes loading and validating rsm structures from rsm files. rsm
% files are saved manually to disk after being created with rsa_create_rsm.
%
% 2020-21 Roberto Viviani - December 2020
% Institute of Psychology, University of Innsbruck
% Last modified: April 3rd 2021

    rsa.create_rsm   = @create_rsm_;
    rsa.load_rsm     = @load_rsm_;
    rsa.validate_rsm = @validate_rsm_;
end

function rsm = create_rsm_(name, mapmx)
% This function creates an rsm structure. This structure has the following
% fields:
%       name    the name of the feature
%       type    one of 'similarity', 'converted', 'dissimilarity'. For
%               documentation purposes only -- always treated as similarity
%       rsm     the representation similarity map (a matrix). Only the
%               lower triangular part is used.
%       idx     (optional) a logical array to select the beta images to
%               which the rsm applies.
% This function only creates the function wrapper. No checks on the
% arguments (use validate_rsm_).

    if 0 == nargin
        rsm = struct('name', {}, 'type', {}, 'RSM', {}, 'idx', {});
        return;
    end
    
    %structure to be converted to rsm -- already validated (fixes order of
    %fields to compose arrays)
    if 1 == nargin && isstruct(name)
        rsm.name = name.name;
        rsm.type = name.type;
        rsm.RSM = name.RSM;
        rsm.idx = name.idx;
        return;
    end
    
    %values to initialize rsm
    if nargin < 2, error('create_rsm_: two arguments required'); end
    rsm = struct('name', name, 'type', 'similarity', ...
        'RSM', mapmx, 'idx', []);
end

function flag = validate_rsm_(rsm)
% Validate rsm as an rsm structure. See create_rsm_. Returns false on
% invalid arguments, printing in the console informative message.
%
    utils = rsa_utils();
    flag = false;
    if isempty(rsm)
        disp('Invalid rsm structure: empty value');
        return;
    end
    if ~isstruct(rsm)
        disp('Invalid rsm structure: not a structure');
        return;
    end
    flds = fields(create_rsm_());
    idx = ~ismember(flds, fields(rsm));
    if any(idx)
        disp('Invalid rsm structure: following fields are missing');
        cellfun(@disp, flds(idx));
        return;
    end
    if ~ischar(rsm.name)
        disp('Invalid rsm structure: invalid name field');
        return;
    end
    if ~ismember({rsm.type}, {'similarity', 'converted', 'dissimilarity'})
        disp('Invalid rsm type: should be one of ''similarity'', ''dissimilarity''');
        return;
    end
    if ~isnumeric(rsm.RSM)
        disp('Invalid rsm map: should be a numeric type');
        return;
    end
    sz = size(rsm.RSM);
    if 2 ~= numel(sz) || sz(1) ~= sz(2)
        disp('Similarity map should be a square matrix');
        return;
    end
    if any(1 ~= diag(rsm.RSM))
        disp('Invalid similarity map (not a similarity matrix?)');
        return;
    end
    msk = utils.triangmsk(size(rsm.RSM,1));
    if any(isnan(rsm.RSM(msk)))
        disp('Invalid similarity map: NaN elements');
        return;
    end
    if 1 == length(unique(rsm.RSM(msk)))
        disp('Invalid similarity map: off diagonal elements are identical');
        return;
    end
    if ~isempty(rsm.idx)
        rsm.idx = utils.validate_idx(rsm.idx);
        if isempty(rsm.idx), return; end
        if islogical(rsm.idx),  cidx = sum(rsm.idx);
        else                    cidx = numel(rsm.idx);
        end
        if cidx ~= size(rsm.RSM, 1)
            disp('An invalid idx field in the rsm structure was provided');
            disp(['The field should select ' int2str(size(rsm.RSM, 1)) ...
                ' elements, it selects ' int2str(sum(rsm.idx)) '.']);
            return;
        end
        rsm.idx = rsm.idx(:)';
    end
    flag = true;
end

function rsmout = load_rsm_(rsmfile)
% load an rsm file, keeping in mind possible different formats. We return
% here an array of similarity maps. This function handles maps generated 
% with the original Nili et al. sofwtare (dissimilarity maps to be 
% converted into similarity), as well as rsm's with an additional field 
% 'type' with content 'similarity' (no conversion). This function will 
% return an empty value on errors at parsing the rsm file.
%
    rsmout = create_rsm_();
    rsm = load(rsmfile);
    if isempty(rsm) || ~isstruct(rsm)
        disp('Invalid map file');
        return;
    end
    if numel(fields(rsm(1))) == 1 && numel(rsm) == 1
        rsm = rsm.(char(fields(rsm)));
    end
    if isempty(rsm)
        disp('Empty map file');
        return;
    end
    if ~isfield(rsm(1), 'name')
        disp(['Could not parse ' rsmfile]);
        return;
    end

    for i = 1 : length(rsm)
        rsm_ = rsm(i);
        if isfield(rsm_, 'RDM')   %Nili et al's RSM struct
            rsm_.type = 'converted';
            rsm_.RSM = 1 - rsm_.RDM ./ max(rsm_.RDM(:));
            rsm_ = rmfield(rsm_, 'RDM');
        end
        if ~isfield(rsm_, 'idx'), rsm_.idx = []; end
        if ~validate_rsm_(rsm_)
            fprintf('-> %s - %s: map not loaded\n', rsmfile, ...
                rsm_.name);
            continue;
        end
        if ~any(strcmp(rsm_.name, cellstr(char(rsmout(:).name))))
            rsmout(end+1) = create_rsm_(rsm_); %#ok<*AGROW>
        else
            disp(['Duplicate map name: ', rsm_.name]);
            rsmout = [];
            return;
        end
    end
end



function str = rsa_utils()
% RSA_UTILS Utility functions for the rsa package, mainly input/output to
% user and mass storage.
%
% 2020-21 Roberto Viviani - December 2020
% Institute of Psychology, University of Innsbruck
% Last modified: March 21th 2021

    str.load_vols = @load_vols_;
    str.voxsize = @voxdim_;
    str.input_ui = @input_ui_;
    str.get_commonpath = @get_commonpath_;
    str.cov2cor = @cov2cor;
    str.validate_idx = @validate_idx_;
    str.triangmsk = @triangmsk_;
end

%Load vols%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vols = load_vols_(hdrvols, hdrmsk)
% Load volumes into a nvols x nvox matrix.

    flag = spm_check_orientations([hdrvols(1), hdrmsk]);
    if flag,  f = @spm_read_vols;
    else
        error('Volume and mask must have the same voxel size and orientation');
    end
    msk = logical(round(f(hdrmsk(1))));
    for i = 2 : length(hdrmsk)
        msk = msk & logical(round(f(hdrmsk(i))));
    end
    cvoxels = sum(msk(:));
    fcount = length(hdrvols);
    vols = zeros([fcount, cvoxels]);
    spm_progress_bar('Init', fcount, 'Loading vols');
    for f = 1 : fcount
        buf = spm_read_vols(hdrvols(f));
        vols(f,:) = buf(msk)';
        spm_progress_bar('Set', f);
    end
    spm_progress_bar('Clear');
end

function voxdim = voxdim_(hdr)
% Compute voxel size from header, a row array

    vprm = spm_imatrix(hdr.mat);
    voxdim = abs(vprm(7:9));
end

%Input ui%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reply = input_ui_(varargin)
% A wrapper function to spm_input, allowing to display an explanatory
% message in the input window.
% See spm_input for credits.

    if 0 == nargin
        figure(spm_figure('GetWin', 'Interactive'));
        return; 
    end
    
    hf = spm_figure('GetWin', 'Interactive');
    set(hf, 'DefaultUicontrolFontSize', 10);
    
    flag = test_comment__(varargin);
    pos = get(gcf, 'Position');
    FONT_SIZE = floor(pos(3) / 38);
    if flag
        comment = varargin{2};
        figure(spm_figure('GetWin', 'Interactive'));
        spm_progress_bar('Clear');
        %write comment
        reset(gca);
        for c = 1 : length(comment)
            if isempty(comment{c}), comment{c} = ''; end
        end
        comment = cellfun(@(x) strrep(x, '_', '\_'), comment, ...
            'UniformOutput', false);
        text(0, .55, char(comment), 'FontSize', FONT_SIZE)
        set(gca, 'Visible', 'off');
        varargin{2} = '1';
    end
    spm('Pointer');
    ver_ = spm('ver');
    if strcmp(ver_, 'SPM99')
        cmd = 'spm_input_ui(';
    elseif strcmp(ver_, 'SPM2')
        cmd = 'spm_input(';
    elseif strcmp(ver_, 'SPM5')
        cmd = 'spm_input(';
    elseif strcmp(ver_, 'SPM8')
        cmd = 'spm_input(';
    elseif strcmp(ver_(1:5), 'SPM12')
        cmd = 'spm_input(';
    else
        error(['Unrecognized SPM version: ', spm('ver')]);
    end
    for a = 1 : nargin
        cmd = [cmd, ' varargin{', int2str(a), '},']; %#ok<*AGROW>
    end
    cmd = cmd(1:length(cmd)-1);
    cmd = [cmd, ');'];
    reply = eval(cmd);
    if flag, cla(gca); end
end

function flag = test_comment__(vrarg)
    flag = 0;
    if length(vrarg) < 2, return; end
    if ~iscell(vrarg{2}), return; end
    if ~ischar(vrarg{2}{1}), return; end
    flag = 1;
end

%Common path%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pth = get_commonpath_(fname)
% This function returns the minimal common path of a character array of
% file names in fname. If no common path is found, an empty value is
% returned. If any of the rows of fname is blank, an empty value is
% returned.

    fname = char(fname);
    if 1 == size(fname,1)
        pth = fileparts(fname);
        return;
    end
    for f = 1 : size(fname,1)
        pth = fileparts(fname(f,:));
        buff{f} = [pth, filesep];
    end
    buff = char(buff);
    indx = double(buff);
    %column_subtract_(indx, double(buff(1,:)));
    subtrh = double(buff(1,:));
    for i = 1 : size(indx, 2)
        indx(:,i) = indx(:,i) - subtrh(i);
    end
    indx = sum(indx);
    pos = 0;
    for p = 1 : length(indx)
        if indx(p) == 0, pos = p;
        else, break;
        end
    end
    if 0 == pos, pth = []; return; end
    %this checks that commonality was on a whole
    %directory name, not just part of it.
    fpos = strfind(buff(1,:), filesep);
    indx = fpos <= pos;
    indx(indx == 0) = [];
    if pos > fpos(length(indx)), pos = fpos(length(indx)); end
    pth = buff(1,1:pos);
end

function cor = cov2cor(cov)
% Covariance to correlation.

    dg = 1 ./ sqrt(diag(cov));
    cor = cov;
    for i = 1 : size(cov, 2)
        cor(:,i) = cor(:,i) .* dg(i);
        cor(i,:) = cor(i,:) .* dg(i);
    end
end

function idx = validate_idx_(idx)
% An index is valid if logical or composed only by 0 and 1, or unique
% integers. If valid, returns the index in the format find(idx), otherwise
% and empty vector.
    if isempty(idx), return; end
    if numel(idx) ~= length(idx)
        disp('Invalid index: should be a vector');
        idx = [];
        return;
    end
    if islogical(idx) || isempty(setdiff(idx, [0 1]))
        idx = find(idx);
    else
        if max(idx - round(idx)) > 0
            disp('Index should be composed of integers');
            idx = [];
            return;
        elseif min(idx) < 1
            disp('Index should not contain zero or negative values');
            idx = [];
            return;
        elseif numel(unique(idx)) ~= numel(idx)
            disp('Index should not contain repeated values');
            idx = [];
            return;
        end
    end
end
        
function msk = triangmsk_(order, offset)
    % Returns a binary mask for a matrix of size order x order, selecting
    % the triangle at offset offset. The lower triangular area is selected
    % (to maintain compatibility with squareform).
    
    if nargin < 2, offset = 0; end
    msk = tril(true(order, order), -(1 + offset));
end
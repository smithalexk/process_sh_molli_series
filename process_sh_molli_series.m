% Process (sh)MOLLI data
%
% Usage 1:
% T1MAP = process_sh_molli_series(directory,(plot),(user selection));
%
% Usage 2:
% process_sh_molli_series(T1MAP);

function T1MAP = process_sh_molli_series(varargin)

% Check user has Image Processing Toolbox
box = 'Image Processing Toolbox';
v = ver;
if ~(any(strcmp(box, {v.Name})))
    disp([box ' required']);
    return;
end

% Check number of arguments
if numel(varargin) < 1
    % No arguemnts - use current directory with user selection
    directory = uigetdir;
    if directory == 0
        directory = '.';
    end
    T1MAP = read_and_process(directory,true);
    interactive_plot(T1MAP);
elseif isa(varargin{1},'char')
    % First argument is a string - assume this is a valid path. No checking
    % currently
    directory = varargin{1};
    % Has the user specified whether they want to select the ROI?
    if numel(varargin) < 3
        user_selection = true;
    else
        user_selection = varargin{3};
    end
    % Run the processing here
    T1MAP = read_and_process(directory,user_selection);
    % Has the user asked to do some plotting?
    if numel(varargin) < 2
        plot_flag = 1;
    else
        plot_flag = varargin{2};
    end
    % Do the plotting here
    if plot_flag
        interactive_plot(T1MAP);
    end
elseif isa(varargin{1},'struct')
    % If a struct is passed, assume that we want to plot the data. No
    % checking is done
    T1MAP = varargin{1};
    interactive_plot(T1MAP);
end

end

function T1MAP = read_and_process(directory,user_selection)

% Read the directory contents
D = dir(directory);
D = D(3:end);

% Some empty places to store the data
im = [];
time = [];
inv_time = [];

% Read all the DICOM images and meta data
for i = 1:length(D)
    info = dicominfo([directory filesep D(i).name]);
    image = dicomread([directory filesep D(i).name]);
    im(:,:,:,i) = image; % Image data
    try
        time(i,1) = info.TriggerTime; % TriggerTime tag
    catch
        %disp('No TriggerTime');
    end
    try
        inv_time(i,1) = info.InversionTime; % InversionTime tag <- this is the one we want
    catch
        %disp('No InversionTime');
    end
end

% Optionally, ask user which to use
if user_selection
    prompt = ['Select which DICOM tag to use:\n 1. TriggerTime   - ' num2str(time') '\n 2. InversionTime - ' num2str(inv_time') '\n   > '];
    reply = input(prompt,'s');
else
    % If the use isn't going to make a choice, do some guess work
    if isempty(time) || ~isempty(inv_time)
        reply = '2';
    elseif isempty(inv_time)
        reply = '1';
    else
        if range(time) > range(inv_time)
            reply = '1';
        else
            reply = '2';
        end
    end
end
if isempty(reply)
    reply = '2';
end
if reply == '1'
    inv_time = time;
end

% Re-order (sh)MOLLI data by InversionTime tag
[inv_time,n] = sort(inv_time);
im = im(:,:,:,n);

% Reshape the data into a column. One row per pixel. One column per time
% point of the data
im_column = reshape(im, size(im,1)*size(im,2), size(im,4));

% Create some masks to flip 0-3 of the first points negative
mask0 = ones(size(inv_time));
mask1 = mask0; mask1(1) = -1;
mask2 = mask0; mask2(1:2) = -1;
mask3 = mask0; mask3(1:3) = -1;
masks = [mask0 mask1 mask2 mask3];

% Choose between python and MATLAB fitting
py_fit = 0;
if ~isempty(pyversion)
    py_fit = 1;
    disp('python fitting to be used');
    dummy_cfit = fit([0,0,0,0]',[1,1,1,1]',fittype('(Axy - Bxy * exp(-tinv/tonestar))','dependent',{'y'},...
        'independent',{'tinv'},'coefficients',{'Axy','Bxy','tonestar'},'options',fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0],'StartPoint',[1000,1000,1000])));
else
    disp('MATLAB fitting to be used (slower)');
end
% Create the fit types - Currently don't use the 'abs' approach
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0],'StartPoint',[1000,1000,1000]);
%abs_molli = fittype('abs(Axy - Bxy * exp(-tinv/tonestar))','dependent',{'y'},...
%    'independent',{'tinv'},'coefficients',{'Axy','Bxy','tonestar'},'options',fo);
molli = fittype('Axy - Bxy * exp(-tinv/tonestar)','dependent',{'y'},...
    'independent',{'tinv'},'coefficients',{'Axy','Bxy','tonestar'},'options',fo);

% Display the first image and let the user choose a poly ROI
if user_selection
    disp('Select ROI to process.\nPress ESC to process the whole image.');
    figure, imagesc(im(:,:,1,1)); axis image; hold on; drawnow;
    h = impoly;
    position = wait(h);
else
    position = [];
end

% Pressing ESC will analyse the whole image
if numel(position) == 0
    position = [1 1 size(im,1) size(im,1); 1 size(im,2) size(im,2) 1]';
end

% Create coordinate list for x and y
[X,Y] = meshgrid(1:size(im,1),1:size(im,2));
X = X(:); Y = Y(:);

% Find all the points inside the polygon and their indicies
ROI_mask = inpolygon(X,Y,position(:,1), position(:,2));% & (range(im_column,2) > 333) & (range(im_column,2) < 4000);
ROI_inds = find(ROI_mask == 1);

% Create some blank variables to store the results in
ALL_FITS = cell(size(ROI_inds));
ALL_SEE = nan(length(ROI_inds), 4);
tic
parfor_progress(numel(ROI_inds));
for I = 1:numel(ROI_inds)
    i = ROI_inds(I); % This is the image index
    this_y = im_column(i,:)';
    
    % Do the fit with each mask
    %tic
    if py_fit
        
        % Import python module
        if count(py.sys.path,'') == 0
            insert(py.sys.path,int32(0),'');
        end
        FunPath = strsplit(mfilename('fullpath'),filesep);
        ModPath = strjoin(FunPath(1:end-1),filesep);
        
        P = py.sys.path;
        if count(P, ModPath) == 0
            insert(P, int32(0), ModPath);
        end
        mod = py.importlib.import_module('sh_molli_fit');
        py.importlib.reload(mod);
        x_fit = py.list(inv_time'); % Only need to do this once
        
        warning('off', 'curvefit:cfit:subsasgn:coeffsClearingConfBounds');
        % Create python lists
        y_fit0 = py.list((this_y.*mask0)');
        y_fit1 = py.list((this_y.*mask1)');
        y_fit2 = py.list((this_y.*mask2)');
        y_fit3 = py.list((this_y.*mask3)');
        % Do python fits
        try py_fit0 = mod.my_fit(x_fit,y_fit0); catch, py_fit0 = {0,0,0, Inf}; end% caught = caught + 1; end
        %py_fit0 = mod.my_fit(x_fit,y_fit0);
        try py_fit1 = mod.my_fit(x_fit,y_fit1); catch, py_fit1 = {0,0,0, Inf}; end% caught = caught + 1; end
        try py_fit2 = mod.my_fit(x_fit,y_fit2); catch, py_fit2 = {0,0,0, Inf}; end% caught = caught + 1; end
        try py_fit3 = mod.my_fit(x_fit,y_fit3); catch, py_fit3 = {0,0,0, Inf}; end% caught = caught + 1; end
        % Hack together some fit objects so I don't have to rewrite the
        % rest of the code
        fit0 = dummy_cfit;
        fit0.Axy = py_fit0{1};
        fit0.Bxy = py_fit0{2};
        fit0.tonestar = py_fit0{3};
        fit1 = dummy_cfit;
        fit1.Axy = py_fit1{1};
        fit1.Bxy = py_fit1{2};
        fit1.tonestar = py_fit1{3};
        fit2 = dummy_cfit;
        fit2.Axy = py_fit2{1};
        fit2.Bxy = py_fit2{2};
        fit2.tonestar = py_fit2{3};
        fit3 = dummy_cfit;
        fit3.Axy = py_fit3{1};
        fit3.Bxy = py_fit3{2};
        fit3.tonestar = py_fit3{3};
        % Save the SSE
        G0.sse = py_fit0{4};
        G1.sse = py_fit1{4};
        G2.sse = py_fit2{4};
        G3.sse = py_fit3{4};
        warning('on', 'curvefit:cfit:subsasgn:coeffsClearingConfBounds');
    else
        [fit0,G0] = fit(inv_time, this_y.*mask0, molli);
        [fit1,G1] = fit(inv_time, this_y.*mask1, molli);
        [fit2,G2] = fit(inv_time, this_y.*mask1, molli);
        [fit3,G3] = fit(inv_time, this_y.*mask3, molli);
        %toc
    end
    % Save all the results
    ALL_FITS{I} = {fit0, fit1, fit2, fit3};
    ALL_SEE(I,:) = [G0.sse G1.sse G2.sse G3.sse];
    
    
    parfor_progress;
end
parfor_progress(0);
toc

% Find the best fits from each fit
BEST_FITS = cell(size(ROI_inds));
BEST_MASKS = [];
BEST_SSE = [];
BEST_IND = nan(size(ROI_inds));
for I = 1:numel(ROI_inds)
    ind = find(ALL_SEE(I,:) == min(ALL_SEE(I,:)));
    ind = ind(1);
    BEST_FITS(I) = ALL_FITS{I}(ind);
    BEST_MASKS(I,:) = masks(:,ind)';
    BEST_IND(I) = ind - 1; % mask0 is at ind1
    BEST_SSE(I) = ALL_SEE(I,ind);
end

% Create an output image
output = zeros(size(im,1), size(im,2));
output_LL = zeros(size(im,1), size(im,2));
FITS = cell(size(im,1), size(im,2));
MASKS = cell(size(im,1), size(im,2));
INDS = nan(size(im,1), size(im,2));
SSE = zeros(size(im,1), size(im,2));
for i = 1:length(ROI_inds)
    try
        % t1vec(i)= Tonestar*(coeffvals(2)/coeffvals(1)-1); % LL correction
        output(ROI_inds(i)) = BEST_FITS{i}.tonestar;%
        output_LL(ROI_inds(i)) = BEST_FITS{i}.tonestar * ((BEST_FITS{i}.Bxy / BEST_FITS{i}.Axy) - 1);
        FITS{ROI_inds(i)} = BEST_FITS{i};
        MASKS{ROI_inds(i)} = BEST_MASKS(i,:);
        INDS(ROI_inds(i)) = BEST_IND(i);
        SSE(ROI_inds(i)) = BEST_SSE(i);
    catch
        i
    end
end

% Any T1s that are unrealistic, set to 0
output_LL(output_LL > 2000 | output_LL < 0) = 0;

% Create the output structure
T1MAP.directory = directory;
T1MAP.input = im;
T1MAP.InversionTime = inv_time;
T1MAP.output = output;
T1MAP.output_LL_corrected = output_LL;
T1MAP.FITS = FITS;
T1MAP.INDS = INDS;
T1MAP.ROI_mask = ROI_mask;
T1MAP.MASKS = MASKS;
T1MAP.im_column = im_column;
T1MAP.SSE = SSE;

end

function interactive_plot(T1MAP)

% Plot the output data together
subplot(2,2,1); imagesc(T1MAP.input(:,:,1,1)); axis image
subplot(2,2,2); imagesc(T1MAP.SSE); axis image;
subplot(2,2,3); imagesc(T1MAP.output_LL_corrected); axis image;

% Activate interactive graph
while true
    [y,x] = ginput(1);
    if round(x) < 1 ...
            || round(x) > size(T1MAP.output,1) ...
            || round(y) < 1 ...
            || round(y) > size(T1MAP.output,2)
        break
    else
        i = sub2ind(size(T1MAP.output), round(x), round(y));
        subplot(2,2,4);
        hold off;
        plot(0);
        hold on;
        if T1MAP.ROI_mask(i)
            this_mask = T1MAP.MASKS{i};
            %plot(FITS{round(x),round(y)});
            plot(T1MAP.InversionTime, T1MAP.im_column(i,:), 'bx'); drawnow;
            plot(T1MAP.InversionTime, T1MAP.im_column(i,:).*this_mask,'rx');
            plot(T1MAP.FITS{round(x),round(y)});
        else
            plot(T1MAP.InversionTime, T1MAP.im_column(i,:), 'bx'); drawnow;
        end
    end
end


end

function y = fit_func(x,a,b,t)
y = (a-b.*exp(-x./t));
end

function percent = parfor_progress(N)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%a
%   PARFOR_PROGRESS updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress;
%      end
%      parfor_progress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/

narginchk(0, 1);

if nargin < 1
    N = -1;
end

percent = 0;
w = 50; % Width of progress bar
tn = now; % Time now

if N > 0
    
    % If no GCP, create it now. This will stop the message being displayed
    % after the parfor display has started
    if isempty(gcp('nocreate'))
        gcp;
    end
    
    f = fopen('parfor_progress.txt', 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%5.8f\n%5.8f\n', N, tn); % Save N and t at the top of progress.txt
    fclose(f);
    
    time_fmt = 'HH:MM:SS';
    t_start = datestr(now, time_fmt);
    t_now = datestr(tn, time_fmt);
    
    %waitbar(0, 'Progress');
    if nargout == 0
        disp(['  0%[>', repmat(' ', 1, w), ']', ' ', t_start, ' -> ', t_now, ' -> ', time_fmt]);
    end
elseif N == 0
    f = fopen('parfor_progress.txt', 'r');
    progress = fscanf(f, '%f');
    fclose(f);
    try
        delete('parfor_progress.txt');
    catch
    end
    percent = 100;
    
    time_fmt = 'HH:MM:SS';
    t_start = datestr(progress(2), time_fmt);
    t_now = datestr(tn, time_fmt);
    t_end = datestr(now, time_fmt);
    
    if nargout == 0
        disp([repmat(char(8), 1, (w+9+33)), newline, '100%[', repmat('=', 1, w+1), ']', ' ', t_start, ' -> ', t_now, ' -> ', t_end]);
    end
    %h = findall(0,'Tag','TMWWaitbar');
    %close(h);
else
    try
        % Much quicker to just try this and catch the error
        % It was taking ~33% of the time on this if statement
        %if ~exist('parfor_progress.txt', 'file')
        %error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
        %end
        
        f = fopen('parfor_progress.txt', 'a');
        fprintf(f, '%5.8f\n', tn);
        fclose(f);
        
        f = fopen('parfor_progress.txt', 'r');
        progress = fscanf(f, '%f');
        fclose(f);
        percent = (length(progress)-2)/progress(1)*100;
        
        % Caclaulated elased time
        t1 = progress(2);
        time_fmt = 'HH:MM:SS';
        t_start = datestr(progress(2), time_fmt);
        t_now = datestr(tn, time_fmt);
        t_estend = datestr((((tn-t1)/percent) * 100) + t1, time_fmt);
        
        if nargout == 0
            perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
            % This is probably quicker as sprintf, but it is currnely an
            % insignificant amout of time compared to the file operations,
            % so may as well keep it as it is.
            disp([repmat(char(8), 1, w+9+33), newline, perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']', ' ', t_start, ' -> ', t_now, ' -> ', t_estend]);
        end
        %h = findall(0,'Tag','TMWWaitbar');
        %waitbar(percent/100, h);
    catch
        error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
    end
end
end
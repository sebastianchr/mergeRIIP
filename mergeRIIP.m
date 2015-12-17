% mergeRIIP
function mergeRIIP(varargin)
if nargin==1
    %Single parameter must be the name of a configuration file
    settings0=ascii2struct(varargin{1});
else
    settings0=struct(varargin{:});
end

default.filenames={};
default.lambda=[];  %Å at P08

default.preview=''; %only the data
default.tth_scaling='standard';   %options: overlap, standard
default.peak_finding='manual'; %options: automatic, manual
default.cell=[];
default.centering ='';  %options: P, F, I, A, B, C
default.save_result='merged.dat';
default.save_settings='settings.ini';

% As standard the default values are used
settings=default;
% If a parameter is defined in the settings-file (as loaded into
% settings0), then this value is used instead
if (isfield(settings0,'preview'));            settings.preview=settings0.preview;  end
if (isfield(settings0,'filenames'));            settings.filenames=settings0.filenames;  end
if (isfield(settings0,'lambda'));               settings.lambda=settings0.lambda;  end
if (isfield(settings0,'tth_scaling'));   settings.tth_scaling=settings0.tth_scaling;  end
if (isfield(settings0,'peak_finding'));  settings.peak_finding=settings0.peak_finding;  end
if (isfield(settings0,'cell'));                 settings.cell=settings0.cell;  end
if (isfield(settings0,'centering'));            settings.centering=settings0.centering;  end
if (isfield(settings0,'save_result'));          settings.save_result=settings0.save_result;  end
if (isfield(settings0,'save_settings'));        settings.save_settings=settings0.save_settings;  end


if ~isempty(settings.preview)
    filebase=settings.preview;
    fdir=dir([filebase '*']);
    filenames={fdir.name};
    number_of_files=length(filenames);
    [tth, int, esd]=readdata(filenames);
    
    figure;
    hold on
    for i=1:number_of_files
        plot(tth{i},int{i});
    end
    
    return
    
end


%Define some general limits for the data. 'master_min' cuts away the
%000-reflection
master_min=2; %degrees
master_max=170;%degrees

if isempty(settings.filenames)
    settings.filenames=get_filenames;
end

if strcmpi(settings.tth_scaling,'standard')
    if isempty(settings.cell)
        settings.cell=get_cell;
    end
    if isempty(settings.centering)
        settings.centering=get_centering;
    end
    if isempty(settings.lambda)
        settings.lambda=get_lambda;
    end
end


%Read in the raw data
[tth0, int0, esd0]=readdata(settings.filenames);
number_of_datasets=length(settings.filenames);

global fig_temp verbose
verbose = 1;
if verbose>1
    fig_temp=figure;
end


%% Settings for auto_peak_find
sampling_rate=20; %distance between points used to determine background
signal_over_background=5; %how intense should a peak be bfore it is
smooth_factor=10;
%considered a peak
%%

%groom dataset

for i=1:number_of_datasets;
    %  clean for NaN
    ind_nan=isnan(int0{i});
    tth0{i}(ind_nan)=[];
    int0{i}(ind_nan)=[];
    esd0{i}(ind_nan)=[];
    %remove stuff outside the master limits
    ind_limits=master_min<tth0{i} & tth0{i}<master_max;
    tth0{i}=tth0{i}(ind_limits);
    int0{i}=int0{i}(ind_limits);
    esd0{i}=esd0{i}(ind_limits);
end

out={};
switch settings.tth_scaling
    case 'overlap'
        %If the image plates are cut at an angle and the overlapping region
        %contain one or more peaks
       
        
        %The final merged data set is collected in these variables.
        % The 2theta scale of the first detector is assumed to be correct
        %the subsequent datasets are scaled with respect to this one
        tth0_scaled{1}=tth0{1};
        off_set=zeros(number_of_datasets,1);
        
        for i=1:number_of_datasets-1
            x1=tth0{i}; x2=tth0{i+1};
            y1=int0{i}; y2=int0{i+1};
            e1=esd0{i}; e2=esd0{i+1};
            
            %Find the overlapping region
            lower_bound=min(x2); upper_bound=max(x1);
            ind1=lower_bound<x1 & x1 < upper_bound;
            ind2=lower_bound<x2 & x2 < upper_bound;
            
            if sum(ind1)==0
                error(['There is no overlap between datasets in 2theta. '...
                    ' Use ''standard'' for ''tth_scaling''.'])
            end
            
            x1o=x1(ind1); x2o=x2(ind2);
            y1o=y1(ind1); y2o=y2(ind2);
            e1o=e1(ind1); e2o=e2(ind2);
            
            
            fig1=figure;
            plot(x1o,y1o,x2o,y2o)
            legend('Click first','Click last')
            switch settings.peak_finding
                case 'automatic'
                    tth_select1=auto_find_peaks(x1o,y1o,sampling_rate,...
signal_over_background,...
smooth_factor);
                    tth_select2=auto_find_peaks(x2o,y2o,sampling_rate,...
signal_over_background,...
smooth_factor);
                case 'manual'
                    questdlg(['Click on the peaks with the mouse. ' ...
                        'End with right click.'],'Select peaks','OK','OK');
                    [tth_select, ~] = getpts(fig1);
                    tth_select1=tth_select(1:2:end-1); %user must end with double-click or right-click
                    tth_select2=tth_select(2:2:end-1); %user must end with double-click or right-click
                    
            end
            n_peaks=length(tth_select1);
            [tth_fit1, esd_fit1]=fit_n_peaks(x1o,y1o,e1o,tth_select1,0.1);
            [tth_fit2, esd_fit2]=fit_n_peaks(x2o,y2o,e2o,tth_select2,0.1);
            
            %make sure the reflections are matched
            d=repmat(tth_fit1,1,n_peaks)-repmat(tth_fit2',n_peaks,1);
            [~, ind]=min(abs(d),[],2);
            
            tth_fit2=tth_fit2(ind);
            off_set(i+1)=off_set(i)+mean((tth_fit2-tth_fit1)); %include esd
            
            %correct the tth-scale and merge data
            tth0_scaled{i+1}=tth0{i+1}-off_set(i+1);
        end
        
        %monitor these shifts
        
        
    case 'standard'
        %get wavelength and unitcell
        lambda=settings.lambda;
        
        A=get_unit_cell(settings.cell);
        %define list of reflection
        %generate list of unique hkl-values
        h=0:15; k=0:15; l=0:15;
        [H,K,L]=ndgrid(h,k,l);
        hkl_all=[H(:), K(:),L(:)];
        
        %Remove systematically absent reflections due to centering
        hkl_all=centering_systematic_absence(hkl_all,settings.centering);
        
        tth_calc_all=hkl2tth(A,lambda,hkl_all);
        % sort accoding to tth
        [~,ind_sort] =sort(tth_calc_all);
        tth_calc_all=tth_calc_all(ind_sort);
        hkl_all=hkl_all(ind_sort,:);
        %hkl's are considered identical if their respective 2theta are
        %separated by less than min_diff
        min_diff=0.0001; %degrees
        d_tth=diff(tth_calc_all);
        ind_remove=find(d_tth<min_diff)+1;
        hkl_all(ind_remove,:)=[];
        tth_calc_all(ind_remove)=[];
        
        
        %loop over datasets
        for i=1:number_of_datasets;
            
            %calculate expected peak positions
            fig1=figure; hold on
            plot(tth0{i}, int0{i})
            
            
            %select those appropriate for each databank
            tth_min=min(tth0{i}); tth_max=max(tth0{i});
            if tth_min<master_min; tth_min=master_min; end
            if tth_max>master_max; tth_max=master_max; end
            
            in_range=tth_min-0.5<=tth_calc_all & tth_calc_all<=tth_max+0.5;
            tth_calc=tth_calc_all(in_range);
            hkl_calc=hkl_all(in_range,:);
            n_hkl_calc=size(hkl_calc,1);
            
            figure(fig1);
            bottom_mark=max(int0{i}); top_mark=bottom_mark*1.05;
            plot(repmat(tth_calc',2,1),...
                repmat([bottom_mark top_mark]',1, n_hkl_calc),'r')
            
            switch settings.peak_finding
                case 'automatic'
                    %automatic peak finding
                    tth_select=auto_find_peaks(tth0{i},int0{i},sampling_rate,...
signal_over_background,...
smooth_factor);
                case 'manual'
                    questdlg(['Click on the peaks with the mouse. ' ...
                        'End with right click.'],'Select peaks','OK','OK');
                    %manual peak selection
                    [tth_select, ~] = getpts(fig1);
                    %user must end with double-click or right-click
                    tth_select=tth_select(1:end-1);
            end
            
            n_select=length(tth_select);
            figure(fig1);
            plot(repmat(tth_select',2,1),...
                repmat([bottom_mark-top_mark 0]',1, n_select),'k')
            
            %determine peakpositions in the raw data
            [tth_fit_raw, tth_fit_esd]=fit_n_peaks(tth0{i},int0{i},...
                esd0{i},tth_select,0.1);
            
            %match observed and calculated reflections
            n_calc=length(tth_calc);
            %Determine difference between all observed and calculated peaks
            d=repmat(tth_fit_raw,1,n_calc)-repmat(tth_calc',n_select,1);
            %For each observed reflection, find the calculated reflection
            %that is nearest
            [d0, ind]=min(abs(d),[],2);
            
            tth_calc_matched=tth_calc(ind);
            hkl_calc_matched=hkl_calc(ind,:);
            
            figure(fig1);
            plot(repmat(tth_calc_matched',2,1),...
                repmat([bottom_mark top_mark]',1, n_select),'b')
            d_obs_calc=tth_fit_raw-tth_calc_matched;
            
            [p1 p1_esd]=fit_line(tth_calc_matched,d_obs_calc,tth_fit_esd);
            
            tth_fit_scaled=(tth_fit_raw-p1(2))./(p1(1)+1);
            tth0_scaled{i}=(tth0{i}-p1(2))./(p1(1)+1);
            out{i}=[hkl_calc_matched tth_calc_matched tth_fit_raw tth_fit_scaled];
            
        end
        
        
        
        
        fig_results=figure; hold on
        for i=1:number_of_datasets
            plot(out{i}(:,4),out{i}(:,5)-out{i}(:,4))
            plot(out{i}(:,4),out{i}(:,6)-out{i}(:,4))
        end
        
        %collect all scaled dataset
        
        
        
        xlabel('2theta calc (degrees)')
        ylabel('obs-calc (degrees)')
        legend(settings.filenames,'interpreter','none')
        
        
    otherwise
        error('Unknown scaling_method')
        
end

% merge the datasets with scaled 2theta-axis
for i=1:number_of_datasets
    if i==1
        tth_merge=tth0_scaled{i};
        int_merge=int0{i};
        esd_merge=esd0{i};
    else
        x2=tth0_scaled{i};
        y2=int0{i};
        e2=esd0{i};
        %check for overlapping regions
        lower_bound=min(x2); upper_bound=max(tth_merge);
        ind1=lower_bound<=tth_merge & tth_merge <= upper_bound;
        ind2=lower_bound<=x2 & x2 <= upper_bound;
        if lower_bound<upper_bound %data IS overlapping
            x1o=tth_merge(ind1); x2o=x2(ind2);
            y1o=int_merge(ind1); y2o=y2(ind2);
            e1o=esd_merge(ind1); e2o=e2(ind2);
            
            y2o1=interp1(x2o,y2o,x1o);
            e2o1=interp1(x2o,e2o,x1o);
            
            x_overlap=x1o;
            %From Statistik - viden fra data page 115
            y_overlap=(y1o./e1o.^2+y2o1./e2o1.^2)./(1./e1o.^2+1./e2o1.^2);
            e_overlap=1./(1./e1o.^2+1./e2o1.^2);
            
            tth_merge=[tth_merge(~ind1); x_overlap; x2(~ind2)];
            int_merge=[int_merge(~ind1); y_overlap; y2(~ind2)];
            esd_merge=[esd_merge(~ind1); e_overlap; e2(~ind2)];
        else %data IS NOT overlapping
            tth_merge=[tth_merge; x2];  %#ok
            int_merge=[int_merge; y2];  %#ok
            esd_merge=[esd_merge; e2];  %#ok
        end
    end
end

ind_remove=isnan(int_merge);
tth_merge(ind_remove)=[];
int_merge(ind_remove)=[];
esd_merge(ind_remove)=[];

if verbose > 0
    figure('Name','Merged data');
    plot(tth_merge,int_merge)
    xlabel('2\theta (degrees)')
    ylabel('Intensity (arb. unit)')
end


save_result_to_file(settings.save_result,[tth_merge,int_merge esd_merge])
%save settings
if settings.save_settings
    struct2ascii(settings,'settings.ini')
end
end

function save_result_to_file(filename,out)
% save a diffraction pattern to a ascii file with a single line header
%input:
%filename: string
%out: n x 3 matrix
fid=fopen(filename,'w');
fprintf(fid,'%16s%16s%16s\r\n','2theta(degree)','intensity','esd');
fprintf(fid,'%16.5f%16.6f%16.6f\r\n',out');
fclose(fid);
end

function xi=auto_find_peaks(x0,y0,...
sampling_rate,...
signal_over_background,...
smooth_factor)
%automatically find peaks in a diffraction pattern. First the background is
%determined. Peaks are defined as values exceeding "signal_over_bckground"
%times the background. One maximum is dertermined for each connected region
%exceeding this threshold.

%sampling_rate=20; %distance between points used to determine background
%signal_over_background=5; %how intense should a peak be bfore it is
%                           %considered a peak
x=x0(1:sampling_rate:end); y=y0(1:sampling_rate:end);
%smooth and ignore 'outliers' i.e. the bragg-peaks
ybg = smooth(x,y,smooth_factor,'rloess');

int0_bg=interp1(x,ybg,x0);
above=y0>int0_bg*signal_over_background;
% above(i)=0 and above(i+1)=1 means a peaks is beginning
% above(i)=1 and above(i+1)=0 means a peaks is endning
above([1 end])=0; %ensure that a peak is always ending/starting at the borders
% find peak ends and beginnings
ind_peak_start=find((above(1:end-1)-above(2:end)==-1));
ind_peak_end=find((above(1:end-1)-above(2:end)==1));

if ~(length(ind_peak_start)==length(ind_peak_end))
    %Checks that the number of peak-beginnings and -endings are equal. I
    %have not foreseen a situation where this is not fulfilled.
    error('Number of peak beginnings and endings is not equal')
else
    n_peaks=length(ind_peak_start);
end
% For each peak region find the maximum value. This is considered the peak
% position
xi=zeros(n_peaks,1);
for k=1:n_peaks
    x_peak=x0(ind_peak_start(k):ind_peak_end(k));
    y_peak=y0(ind_peak_start(k):ind_peak_end(k));
    [~, ind]=max(y_peak);
    xi(k)=x_peak(ind);
end
end

function par_guess=make_peak_guess(tth,int)
%Estimate the intensity, position, width and background for a single
%gaussian on a constant backgroun
[position_guess, int_guess, ~]=get_peak_centre(tth,int);
FWHM_guess=get_peak_width(tth,int,0.5);
width_guess=FWHM_guess/2.3548;
%Background is estimated as the average of the endpoints
bg_guess=mean(int([1 end]));
par_guess=[int_guess, position_guess, width_guess, bg_guess];

    function FW=get_peak_width(x,y,fraction)
        [~,yc, indc]=get_peak_centre(x,y);
        %left side of peak
        xl=x(1:indc-1); yl=y(1:indc-1);
        y0=yc*fraction;
        [~,indn]=minn(abs(yl-y0),2);
        xn=xl(indn);
        yn=yl(indn);
        [xn, ind]=sort(xn);
        yn=yn(ind);
        xlf=(y0-yn(1))/(yn(2)-yn(1))*(xn(2)-xn(1))+xn(1);
        
        
        
        %right side of peak
        xr=x(indc+1:end); yr=y(indc+1:end);
        [~,indn]=minn(abs(yr-y0),2);
        xn=xr(indn);
        yn=yr(indn);
        [xn, ind]=sort(xn);
        yn=yn(ind);
        xrf=(y0-yn(1))/(yn(2)-yn(1))*(xn(2)-xn(1))+xn(1);
        FW=xrf-xlf;
    end
    function [xc, yc, ind]=get_peak_centre(x,y)
        [yc, ind]=max(y);
        xc=x(ind);
    end
    function [xm, ind]=minn(x,n)
        %get the n lowest numbers
        x=x(:);
        if length(x)<n
            error(['Number of requested values exceeds the number ' ...
                'of elements in input'])
        end
        xm=zeros(n,1);
        ind=zeros(n,1);
        for i=1:n
            [xm(i), ind(i)]=min(x);
            x(ind(i))=NaN;
        end
    end
end

function [par, par_esd]=fit_peak(x,y,esd,startpoints)
% Fit a single gaussian to the
% w=1./esd; %NOTE NOTE this is not correct but if square the BG is weighted to high
w=ones(size(y));
[xData, yData, weights] = prepareCurveData( x, y, w );


% Set up fittype and options.
ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0];
opts.StartPoint = startpoints(1:4);
% opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

par=coeffvalues(fitresult);
par_esd=coefferror(fitresult,gof);

% Plot fit with data.
global fig_temp verbose

if verbose>1
    figure(fig_temp);
    h = plot( fitresult, xData, yData );
    legend( h, 'Fit', 'Data', 'Location', 'NorthEast' );
    % Label axes
    xlabel('2\theta (degrees)')
    ylabel('Intensity')
    waitforbuttonpress
end


end

function [xi_fit, ei_fit]=fit_n_peaks(x,y,e,xi,fw)
n=length(xi);
xi_fit=zeros(n,1);
ei_fit=zeros(n,1);
for j=1:n
    xmin=xi(j)-fw/2;
    xmax=xi(j)+fw/2;
    ind_fit=xmin<x & x<xmax;
    x_fit=x(ind_fit);
    y_fit=y(ind_fit);
    e_fit=e(ind_fit);
    
    starting_par=make_peak_guess(x_fit,y_fit);
    [par_fit, esd_fit]=fit_peak(x_fit,y_fit,e_fit,starting_par);
    xi_fit(j,:)=par_fit(2);
    ei_fit(j)= esd_fit(2);
end
end

function [par, par_esd]=fit_line(x,y,esd)

% w=1./esd.^2;
w=ones(size(y));
[xData, yData, weights] = prepareCurveData( x, y, w );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

par=coeffvalues(fitresult);
par_esd=coefferror(fitresult,gof);

end

function [hkl_sa, ind]=centering_systematic_absence(hkl,centering)
switch centering
    case 'F'
        ind=(mod(sum(hkl(:,[1 2]),2),2)==0 & ...
            mod(sum(hkl(:,[2 3]),2),2)==0  & ...
            mod(sum(hkl(:,[1 3]),2),2)==0 );
    case 'I'
        ind=(mod(sum(hkl,2),2)==0);
    case 'P'
        ind=ones(size(hkl,1),1);
    case 'A'
        ind=(mod(sum(hkl(:,[2 3]),2),2)==0);
    case 'B'
        ind=(mod(sum(hkl(:,[1 3]),2),2)==0);
    case 'C'
        ind=(mod(sum(hkl(:,[1 2]),2),2)==0);
    otherwise
        ind=ones(size(hkl,1),1);
end

hkl_sa=hkl(ind,:);


end

function error=coefferror(c,gof)

con=confint(c);
width=abs(con(2,:)-con(1,:));
if nargin==2
    degrees_of_freedom=gof.dfe;
else
    degrees_of_freedom=Inf;
    disp('Assumes infinite degrees of freedom')
end
error=width/(2*tinv(0.975,degrees_of_freedom));

end

function files=get_filenames()

%Open GUI to select muliple files
[tempFileName,tempPathName] = uigetfile({'*.dat'},...
    'Select Image Plate File (tif/dat)','multiselect','on');

% End script execution cleanly if the user press cancel
if ~iscell(tempFileName)
    if tempFileName == 0
        error('You pressed cancel. No file selected.')
    end
end

number_of_files=length(tempFileName);
files=cell(number_of_files,1);
for i=1:number_of_files
    files{i}=[tempPathName tempFileName{i}];
end
end

function uc=get_cell
prompt = {'a (Å)' 'b (Å)' 'c (Å)' '\alpha' '\beta' '\gamma'};
defaultanswer = {'5' '5' '5' '90' '90' '90'};
answer=inputdlg(prompt,'Cell parameters',1,...
    defaultanswer);
if isempty(answer)
    error('User did not define mandatory settings. Program ended.')
else
    uc=zeros(1,6);
    for i = 1:6
        uc(i) = str2double(answer{i});
    end
end
end

function centering=get_centering
str = {'P','F', 'I','A','B','C'};
str_code = str;
[s,v] = listdlg('PromptString',...
    'Centering of unit cell',...
    'SelectionMode','single',...
    'ListString',str,...
    'Listsize',[160 80]);
if v
    centering=str_code{s};
else
    error('Did not choose a shape for integration path')
end
end

function lambda=get_lambda
prompt = {'\lambda (Å)'};
defaultanswer = {'0.495837'};
answer=inputdlg(prompt,'Provide wavelength:',1,...
    defaultanswer);
if isempty(answer)
    error('User did not define mandatory settings. Program ended.')
else
    lambda = str2double(answer{1});
end
end

function [tth, int, esd]=readdata(files)

number_of_files=length(files);
tth=cell(1,number_of_files);
int=cell(1,number_of_files);
esd=cell(1,number_of_files);


for i=1:number_of_files
    try
        dat=importdata(files{i});
    catch err
        if strcmp(err.identifier,'MATLAB:importdata:FileNotFound')
            error('Could not find the file: %s', files{i})
        else
            rethrow(err)
        end
        
    end
    if isstruct(dat)
        tth{i}=dat.data(:,1);
        int{i}=dat.data(:,2);
        esd{i}=dat.data(:,3);
    else
        tth{i}=dat(:,1);
        int{i}=dat(:,2);
        esd{i}=dat(:,3);
    end
end
end

function tth0=hkl2tth(A,lambda0,hkl0)
% hkl is n x 3 matrix
% A is 3 x 3 matrix

M=2*pi*(A^-1);
qnorm_aa=sqrt(sum((M*hkl0').^2,1))';

qnorm2tth = @(q_aa,lambda) 2*asind(q_aa*lambda/(4*pi));
tth0=qnorm2tth(qnorm_aa,lambda0);

end
function A=get_unit_cell(unit_cell)

if isequal(size(unit_cell),[3,3])
    A=unit_cell;
elseif isequal(size(unit_cell),[6,1]) || isequal(size(unit_cell),[1,6])
    %from Giacovazzo s 75
    a=unit_cell(1); b=unit_cell(2); c=unit_cell(3);
    alpha=unit_cell(4); beta=unit_cell(5); gamma=unit_cell(6);
    cosalphastar=(cosd(beta)*cosd(gamma)-cosd(alpha))/(sind(beta)*sind(gamma));
    V=a*b*c*(1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma)).^.5;
    cstar=a*b*sind(gamma)/V;
    A=[a 0 0; b*cosd(gamma) b*sind(gamma) 0 ; c*cosd(beta) -c*sind(beta)*cosalphastar 1/cstar];
else
    error('The format of ''unit_cell'' was not recognized')
end
end

function settings=ascii2struct(inputfile)
fid=fopen(inputfile);
if fid==-1
    error('Settings file (%s) not found',inputfile)
end
settings_str=fread(fid,'*char')';

settings=[];
set=regexp(settings_str,'(\w+)[ ]*\=(.+?)[\r\n^]+','tokens');

for i=  1:length(set)
    field_name=set{i}{1};
    field_str=set{i}{2};
    
    %str2num(field_str) will return empty if field str does not contain a
    %string which can be interpreted as a number; however e.g. 'line' and
    %'rectangle' are names of MATLAB are matlab objects and they will be
    %interpreted as their numeric values
    %The work around below therefore first checks for letters (except e).
    a=regexp(field_str,'({.*})','tokens');
    if ~isempty(a) %is cell-array
        try
            eval(['field_value = ' field_str ';'])
        catch err
            error(fprintf('Failed to read keyword: %s',field_str))
            
        end
    else
        %does it contain letters
        if isempty(regexp(field_str,'[a-zA-Z]','once'));
            %could it be a number
            if ~isempty(regexp('1e234','^[0-9.,+ ]*[eE][0-9+-]*$','match'));
                %converting to number
                field_value=str2num(field_str);  %#ok
            else
                %saving as string
                field_value=regexprep(field_str,'["'' ]','');
            end
        else
            field_value=str2num(field_str);  %#ok
            if isempty(field_value) %could not convert to number
                %saving as string
                field_value=regexprep(field_str,'["'' ]','');
            end
        end
    end
    
    settings.(field_name)=field_value;
end
end

function struct2ascii(settings,outputfile)


fid=fopen(outputfile,'w');

all_fields=fieldnames(settings);
for i=1:length(all_fields)
    name =all_fields{i};
    value=settings.(name);
    outstr=ex_func(value,6);
    fprintf(fid,'%s = %s \r\n',name,outstr);
    
end

fclose(fid);

    function [ out ] = ex_func( in, prec )
        
        in_datatype = class(in);
        
        switch in_datatype
            case 'char'
                out =  in ;
            case {'single','double','logical'}
                out=matrix2str(double(in),prec);
                out=out(2:end-1); %remove the surrounding square brackets
            case 'cell'
                out=cell2str(in,prec);
            otherwise
                error('Unknown type');
        end
        
        function str=cell2str(in_cell,prec0)
            str='{';
            for j=1:length(in_cell)
                switch class(in_cell{j})
                    case 'double'
                        str=[str  matrix2str(in_cell{j},prec0) ' ,'];   %#ok
                    case 'char'
                        str=[str  '''' in_cell{j} '''' ' ,'];   %#ok
                    case 'cell'
                        str=[str  cell2str(in_cell{j},prec0) ' ,'];  %#ok
                end
            end
            
            str=str(1:end-1); %remove the last comma
            str=[str '}']; %close with curly bracket
            
        end
        
        function str=matrix2str(in_mat,prec0)
            str='[';
            for j=1:size(in_mat,1)
                str =[str num2str(in_mat(j,:),prec0) ';'];   %#ok
            end
            if length(str) > 1
                str=str(1:end-1); %remove the last ";"
            end
            str=[str ']']; %close
        end
        
        
    end
end
function varargout = iCSD_0_1_1(CSDMode,ePot,varargin)

% implementation of iCSD toolbox to use as functions
%
% [...] = iCSD_0_1_1(CSDMode,ePot, ... )
%
% ePot .... evoked potential
%
% CSDMode
%
% CSD = iCSD_0_1_1('STANDARD', ... ,el_pos,b0,b1,cond,VakninFlag)
% el_pos ....... electrode positions in mm
% b0 ........... threepoint filter: center
% b1 ........... threepoint filter: neighbor
% exCond ....... extracellular cortical conductivity parameter S/m
% VakninFlag ... repeat first and last channel
%
% CSD = iCSD_0_1_1('DELTA', ... ,el_pos,b0,b1,exCond,topCond,diam)
% el_pos ....... electrode positions in mm
% b0 ........... threepoint filter: center
% b1 ........... threepoint filter: neighbor
% exCond ....... extracellular cortical conductivity parameter S/m
% topCond ...... topical (above cortex) cortical conductivity parameter S/m
% diam ......... 
%
% [CSD,z] = iCSD_0_1_1('STEP', ... ,el_pos,gauss_sigma,exCond,topCond,diam)
% el_pos ........ electrode positions in mm
% gauss_sigma ... gaussian filter: std. dev.
% exCond ........ extracellular cortical conductivity parameter S/m
% topCond ....... topical (above cortex) cortical conductivity parameter S/m
% diam .......... 
%
% [CSD,z] = iCSD_0_1_1('SPLINE', ... ,el_pos,gauss_sigma,exCond,topCond,diam)
% NOT IMPLEMENTED YET
%
%
% from CSDplotter-0.1.1 by Klas H. Pettersen

switch upper(CSDMode)
    case 'STANDARD'
         CSD = standardCSD(ePot,varargin{:});
         varargout{1} = CSD;
    case 'DELTA'
        CSD = deltaCSD(ePot,varargin{:});
         varargout{1} = CSD;
    case 'STEP'
        [CSD,zs] = stepCSD(ePot,varargin{:});
         varargout{1} = CSD;
         varargout{2} = zs;
    case 'SPLINE'
        [CSD,zs] = splineCSD(ePot,varargin{:});
         varargout{1} = CSD;
         varargout{2} = zs;
    otherwise
        CSD = standardCSD(ePot,el_pos,0,0,0.3,false);
         varargout{1} = CSD;

end
        
function [CSD,el_pos_plot] = standardCSD(ePot,el_pos,b0,b1,cond,VakninFlag)
% cond ......... cortical parameter S/m
% b0 ........... threepoint filter: center
% b1 ........... threepoint filter: neighbor
% el_pos ....... electrode positions in mm
% VakninFlag ... repeat first and last channel

% filter parameters:
if b0+2*b1 == 0 && b1~=0
    errordlg('Singularity: b0+2*b1 cannot equal zero.');
    return
end;

% electrical parameters:
if cond<=0
    errordlg('ex. cond. has to be a positive number');
    return
end;

% size, potential (m1 has to equal number of electrode contacts)
[m1,m2] = size(ePot);

% electrode parameters:
el_pos = el_pos.*1e-3; % mm -> m
el_pos_plot = el_pos(2:length(el_pos)-1); % if not Vaknin electrodes
N = length(el_pos);
h = mean(diff(el_pos));
pot = ePot;
if m1~=N
    errordlg(['Number of electrode contacts has to equal number of rows in potential matrix. Currently there are ',...
        num2str(N),' electrodes contacts, while the potential matrix has ',num2str(m1),' rows.']) 
    return
end;

% compute standard CSD with vaknin el.
if VakninFlag
  el_pos_plot = el_pos;
  pot(1,:) = ePot(1,:);
  pot(2:m1+1,:)=ePot;
  pot(m1+2,:)=ePot(m1,:);
end;

CSD = -cond*D1(length(pot(:,1)),h)*pot;

if b1~=0 %filter iCSD (does not change size of CSD matrix)
  [n1,n2]=size(CSD);            
  CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
  CSD_add(n1+2,:)=zeros(1,n2);
  CSD_add(2:n1+1,:)=CSD;        %CSD_add has n1+2 rows
  CSD = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows
end;

function [CSD] = deltaCSD(ePot,el_pos,b0,b1,exCond,topCond,diam)
% exCond ......... cortical parameter S/m
% topCond .......
% b0 ........... threepoint filter: center
% b1 ........... threepoint filter: neighbor
% el_pos ....... electrode positions in mm
% diam ... repeat first and last channel

% filter parameters:
if b0+2*b1 == 0 && b1~=0
    errordlg('Singularity: b0+2*b1 cannot equal zero.');
    return
end;

% electrical parameters:
if exCond<=0
    errordlg('Ex. cond. has to be a positive number');
    return
end;

% size, potential (m1 has to equal number of electrode contacts)
[m1,m2] = size(ePot);

el_pos = el_pos.*1e-3; % mm -> m

% geometrical parameters:
diam = diam*1e-3; %diameter in [m]
if diam<=0
    errordlg('Diameter has to be a positive number.');
    return
end;

if topCond~=exCond && (el_pos~=abs(el_pos) || length(el_pos)~=length(nonzeros(el_pos)))
    errordlg('Electrode contact positions must be positive when top cond. is different from ex. cond.')
    return;
end;
if m1~=length(el_pos)
    errordlg(['Number of electrode contacts has to equal number of rows in potential matrix. Currently there are ',...
        num2str(length(el_pos)),' electrodes contacts, while the potential matrix has ',num2str(m1),' rows.']) 
    return
end;

% compute delta iCSD:
CSD = F_delta(el_pos,diam,exCond,topCond)^-1*ePot;

if b1~=0 %filter iCSD
  [n1,n2]=size(CSD);            
  CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
  CSD_add(n1+2,:)=zeros(1,n2);
  CSD_add(2:n1+1,:)=CSD;        %CSD_add has n1+2 rows
  CSD = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows
end;

function [new_CSD_matrix,zs] = stepCSD(ePot,el_pos,gauss_sigma,exCond,topCond,diam)
% filter parameters:
gauss_sigma = gauss_sigma*1e-3;
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
if gauss_sigma<0
    errordlg('The gaussian filter width cannot be negative.')
    return
end;

% electrical parameters:
if exCond<=0
    errordlg('Ex. cond. has to be a positive number');
    return
end;

% size, potential (m1 has to equal number of electrode contacts)
[m1,m2] = size(ePot);

el_pos = el_pos.*1e-3; % mm -> m

% geometrical parameters:
diam = diam*1e-3; %diameter in [m]
if diam<=0
    errordlg('Diameter has to be a positive number.');
    return
end;

if topCond~=exCond & (el_pos~=abs(el_pos) | length(el_pos)~=length(nonzeros(el_pos)))
    errordlg('Electrode contact positions must be positive when top cond. is different from ex. cond.')
    return;
end;
if m1~=length(el_pos)
    errordlg(['Number of electrode contacts has to equal number of rows in potential matrix. Currently there are ',...
        num2str(length(el_pos)),' electrodes contacts, while the potential matrix has ',num2str(m1),' rows.']) 
    return
end;

% compute step iCSD:
CSD = F_const(el_pos,diam,exCond,topCond)^-1*ePot;

% make CSD continous (~200 points):
le = length(el_pos);
h = el_pos(2)-el_pos(1);
first_z = el_pos(1)-h/2; %plot starts at z1-h/2;
mfactor = ceil(200/le);
minizs = 0:h/mfactor:(mfactor-1 )*h/mfactor;
for i=1:size(CSD,1) % all rows
    zs((1:mfactor)+mfactor*(i-1)) = first_z+(i-1)*h+minizs;
    new_CSD_matrix((1:mfactor)+mfactor*(i-1),:) = repmat(CSD(i,:),mfactor,1);
end;


if gauss_sigma~=0 %filter iCSD
  [zs,new_CSD_matrix]=gaussian_filtering(zs,new_CSD_matrix,gauss_sigma,filter_range);
%  [new_positions,gfiltered_spline_CSD]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
end;

function [CSD_cs,zs] = splineCSD(ePot,el_pos,gauss_sigma,exCond,topCond,diam)

% filter parameters:
gauss_sigma = gauss_sigma*1e-3; % mm -> m
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
if gauss_sigma<0
    errordlg('The gaussian filter width cannot be negative.')
    return
end;

% electrical parameters:
if exCond<=0
    errordlg('Ex. cond. has to be a positive number');
    return
end;

% size, potential (m1 has to equal number of electrode contacts)
[m1,m2] = size(ePot);

% geometrical parameters:
diam = diam*1e-3; %diameter in [m]
if diam<=0
    errordlg('Diameter has to be a positive number.');
    return
end;

el_pos = el_pos*1e-3;
if topCond~=exCond && (el_pos~=abs(el_pos) || length(el_pos)~=length(nonzeros(el_pos)))
    errordlg('Electrode contact positions must be positive when top cond. is different from ex. cond.')
    return;
end;
if m1~=length(el_pos)
    errordlg(['Number of electrode contacts has to equal number of rows in potential matrix. Currently there are ',...
        num2str(length(el_pos)),' electrodes contacts, while the potential matrix has ',num2str(m1),' rows.']) 
    return
end;

% compute spline iCSD:
Fcs = F_cubic_spline(el_pos,diam,exCond,topCond);
[zs,CSD_cs] = make_cubic_splines(el_pos,ePot,Fcs);
%[pos1,my_CSD_spline]=new_CSD_range(zs,CSD_cs,0,2.4e-3);
if gauss_sigma~=0 %filter iCSD
  [zs,CSD_cs]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
%  [new_positions,gfiltered_spline_CSD]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
end;













function out = D1(N,h)
%function out = D1(N,h)
%
%The matrix form of the standard double derivative formula, called D1 in
%Freeman and Nicholson (1975).
%
%N: number of electrodes
%h: inter-contact distance.
%
%out is a (N-2)x(N) matrix.

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

if nargin < 2, h = 100e-6 ;end;
if nargin < 1, N = 20 ;end;

for i=1:N-2
    for j=1:N
        if (i == j-1)
            out(i,j) = -2/h^2;
        elseif (abs(i-j+1) == 1)
            out(i,j) = 1/h^2;
        else
            out(i,j) = 0;
        end;
    end;
end;

function out = S_general(N,b0,b1)
%S = S_general(N,b0,b1)
%This is the three point filter matrix.
%Returns matrix of size (N-2)x(N),
%which represents a three point "spatial noise" filter with mid
%coeffescient b0 (will be normalized by the function) and neighbouring
%coeffescients b1 (will also be normalized).
%Default filter has b0 = 2 and b1 = 1 and number_of_electrodes = 20.
%
%The Hamming-filter has b0 = 0.54 and b1 = 0.23.

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

if nargin < 1, N = 20; end
if nargin < 3, b1 = 1; b0 = 2; end;

c = b0 + 2*b1;

out = zeros(N-2,N);
for i=1:N-2
    for j=1:N
        if (i == j-1)
            out(i,j) = b0/c;
        elseif (abs(i-j+1) == 1)
            out(i,j) = b1/c;
        end;
    end;
end;

function out = F_delta(el_pos,d,cond,cond_top)
%function out = F_delta(el_pos,d,cond,cond_top)
%
%Computes the F-matrix from infinitesimally thin current source density
%sheets with diameter d and homogenous activity throughout the sheet.
%
%el_pos:    the z-positions of the electrode contacts, default:
%100e-6:100e-6:2300e-6 
%d:         activity diameter, default: 500e-6
%cond:      cortical conductivity, default: 0.3
%cond_top:  conductivity on top of cortex, default: cond
%
%out is a (number_of_electrodes)x(number_of_electrodes) matrix. 

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

%addpath('fCSD_matrixes');

if nargin < 1, el_pos = 100e-6:100e-6:2300e-6; end
if nargin < 2, d = 500e-6; end;
if nargin < 3, cond = 0.3; end;
if nargin < 4, cond_top = cond; end;

%DEFINE FILENAME
N = length(el_pos);
r_off = 0;
z1 = el_pos(1);
h = el_pos(2)-z1;

% compute matrix
out = zeros(N);
for j=1:N                     %zj is position of CSD-plane
    zj = z1 + (j-1)*h;
    for i=1:N                   %zi is position of electrode
        zi = z1 + (i-1)*h;
        out(j,i) = h/(2*cond)*((sqrt((zj-zi)^2+(d/2)^2)-abs(zj-zi))+ ...
            (cond-cond_top)/(cond+cond_top)*(sqrt((zj+zi)^2+(d/2)^2)-abs(zj+zi)));
    end; %for i
end; %for j

% handle save file
% 
% full_filename = [matrix_folder() filesep 'Fd' make_filename(d,r_off,N,h,z1,cond,cond_top) '.mat'];
% 
% try,
%   load(full_filename,'Fd');
%   out = Fd;
% catch,
%   msgstr = lasterr;
%   
%   
%   Fd = out;
%   save(full_filename, 'Fd');
% end;  %catch

function out_string = make_filename(diameter,off,N,h,z1,cond,cond_top)
%out_string = make_filename(diameter,off,N,h,z1,cond,cond_top)
%
%makes filename (without ending) from parameters used in simulation
%diameter: diameter
%off: off center position of electrode [m]
%N: number of electrodes
%z1: position of first electrode [m]
%cond: extracellular conductivity [S/m]
%cond_top: conductivity above cortex [S/m]

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html
tdiameter=['_d' num2str(diameter*1e6)];
if off=='0'; toff='_off0'; elseif off==0; toff=''; else; toff = ['_off' num2str(off/(diameter/2))];end;
if N==23; tN=''; else; tN=['_N' num2str(N)];end;
if h==100e-6; th=''; else; th=['_h' num2str(h*1e6)];end;
if z1==100e-6; tf=''; else; tf=['_f' num2str(z1*1e6)];end;
if cond==0.05; tcond='';else;tcond=['_s' num2str(cond)];end;
if cond_top==cond; tcond_top='';else;tcond_top = ['_t' num2str(cond_top)];end;

out_string = [tdiameter toff tN th tf tcond tcond_top];

function out = matrix_folder()
%Creates the string (path) of the matrix folder which is used to store the
%inversion matrixes for the iCSD methods: this_folder/methods/saved.
%Creates this folder if it does not exist.

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

[pathstr, name, ext] = fileparts(mfilename('fullpath'));
out = [pathstr filesep 'saved']; %folder in which matrixes are saved
if exist(out,'dir')==0; mkdir(out); end; %create folder if it does not exist

function Fc = F_const(el_pos,d, cond, cond_top,this_tol)
%function Fc = F_const(el_pos,d, cond, cond_top,this_tol)
%
%Gives the transformation matrix F from CSDs to potential for the constant
%iCSD method.
%
%el_pos: electrode positions, default: 100e-6:100e-6:2300e-6
%d: activity diameter, default: 500e-6
%cond: extracellular conductivity, default: 0.3
%cond_top: conductivity above cortex, default: cond
%this_tol: tolerance of the integral, default: 1e-6


%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

%DEFAULT VALUES
if nargin<1, el_pos = 100e-6:100e-6:2300e-6;end;
if nargin<2, d = 500e-6; end;
if nargin<3, cond = 0.3; end;
if nargin<4, cond_top = cond; end;
if nargin<5, this_tol = 1e-6; end;

%DEFINE FILENAME
N = length(el_pos);
r_off = 0;
z1 = el_pos(1);
h = el_pos(2)-z1;
full_filename = [matrix_folder() filesep 'Fc' make_filename(d,r_off,N,h,z1,cond,cond_top) '.mat'];

%GET/COMUTE Fc
try, %see if Fc exists
  load(full_filename,'Fc','tol');
  if tol>this_tol
    Fc = compute_F_constant(full_filename,el_pos,d,cond,cond_top,this_tol); %local
  end;
catch, %compute Fc
  msgstr_Fc = lasterr;
  Fc = compute_F_constant(full_filename,el_pos,d,cond,cond_top,this_tol);
end;

function out = compute_F_constant(full_filename,el_pos,d,cond,cond_top,tol);
  N = length(el_pos);
  h = el_pos(2)-el_pos(1);
  out = zeros(N);        %define matrix
  for j = 1:N            %rows
    zj = el_pos(j);   %electrode positions
    for i = 1:N %columns
        if i~=1 %define lower integral limit
          lower_int = el_pos(i)-(el_pos(i)-el_pos(i-1))/2;
        else
          lower_int = max(0,el_pos(i)-h/2);
        end;
        if i~=N %define upper integral limit
          upper_int = el_pos(i)+(el_pos(i+1)-el_pos(i))/2;
        else
          upper_int = el_pos(i)+h/2;
        end;
  %      zi = el_pos(i);   %mid CSD position    
        out(j,i) = quad(@f_cylinder,lower_int,upper_int,tol,[],zj,d,cond) ...
            +(cond-cond_top)/(cond+cond_top) ...
            *quad(@f_cylinder,lower_int,upper_int,tol,[],-zj,d,cond);
    end; %for i
  end; %for j
  Fc = out;
  save(full_filename, 'Fc','tol');
return;

function out1 = f_cylinder(zeta,z,diam,cond)
  out1 = 1./(2.*cond).*(sqrt((diam/2)^2+((z-zeta)).^2)-abs(z-zeta));
return;

function [new_positions,gfiltered_CSD] = gaussian_filtering(positions,unfiltered_CSD,gauss_sigma,filter_range)
%function [new_positions, gfiltered_CSD]= ...
%gaussian_filtering(positions,unfiltered_CSD,gauss_sigma,filter_range)
%
%This function filters the CSD using a gaussian filter.
%
%positions:     The CSD positions
%unfiltered_CSD: The unfiltered CSD matrix
%gauss_sigma:   standard deviation of the gaussian
%filter_range:  the filter width, default: 5*gauss_sigma

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

if nargin<4; filter_range = 5*gauss_sigma; end;

step = positions(2)-positions(1);
filter_positions = -filter_range/2:step:filter_range/2;
gaussian_filter = 1/(gauss_sigma*sqrt(2*pi))*exp(-filter_positions.^2/(2*gauss_sigma^2));
filter_length = length(gaussian_filter);
[m,n]=size(unfiltered_CSD);
temp_CSD=zeros(m+2*filter_length,n);
temp_CSD(filter_length+1:filter_length+m,:)=unfiltered_CSD(:,:); % one filter length of zeros on each side
scaling_factor = sum(gaussian_filter);
temp_CSD = filter(gaussian_filter/scaling_factor,1,temp_CSD); % filter works such that the first filter_length positions are crap
gfiltered_CSD=temp_CSD(round(1.5*filter_length)+1:round(1.5*filter_length)+m,:); % first filter_length is crap, next 0.5 filter length corresponds to positions smaller than the original positions
new_positions = positions;

function F = F_cubic_spline(el_pos,d,cond,cond_top,this_tol)
%function F = F_cubic_spline(el_pos,d,cond,cond_top,this_tol)
%
%Modified from general_spline_method_v3
%
%Creates the F matrix of the cubic spline method.
%
%el_pos:    the z-positions of the electrode contacts, default:
%100e-6:100e-6:2300e-6 
%d:         activity diameter, default: 500e-6
%cond:      cortical conductivity, default: 0.3
%cond_top:  conductivity on top of cortex, default: cond
%this_tol: tolerance of integral, default: 1e-6

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

if nargin<1; el_pos = 100e-6:100e-6:2300e-6;end;
if nargin<2; d = 500e-6; end;
if nargin<3; cond = 0.3; end;
if nargin<4; cond_top = cond; end;
if nargin<5; this_tol=1e-6; end;

tfilename = make_filename(d,0,length(el_pos),el_pos(2)-el_pos(1),el_pos(1),cond,cond_top); %part of Fd filename
full_filename = [matrix_folder() filesep 'Fcs' tfilename '.mat'];

try,
  load(full_filename,'Fcs','tol');
  if tol<=this_tol
    F = Fcs;
  else
    F = compute_F_cubic_spline(full_filename,el_pos,d,cond,cond_top,this_tol);
  end;
catch,
  msgstr = lasterr;
  F = compute_F_cubic_spline(full_filename,el_pos,d,cond,cond_top,this_tol);
end;

function F = compute_F_cubic_spline(full_filename,el_pos,d,cond,cond_top,tol)
  % Define positions and constants
  N = length(el_pos);     %number of electrodes
  z_js = zeros(1,N+2);            %declare electrode positions included
  z_js(1,2:N+1) = el_pos(1,:);    %two imaginary, first electrode in z = 0
  h_av = sum(diff(el_pos))/(N-1); %average inter-contact distance
  z_js(1,N+2) = z_js(1,N+1)+h_av; %last imaginary electrode position

  [E0,E1,E2,E3] = compute_Ematrixes(el_pos);

  %   Define integration matrixes
  F0  = zeros(N,N+1);
  F1  = zeros(N,N+1);
  F2  = zeros(N,N+1);
  F3  = zeros(N,N+1);

  for j = 1:N            %rows
 %   progress = [num2str(j) ' of ' num2str(N) ': int. tolereance: ' num2str(tol)]
    for i = 1:N+1           %columns
      F0(j,i) = quad(@f0,z_js(i),z_js(i+1),tol,[],z_js(j+1),d,cond);
      F1(j,i) = quad(@f1,z_js(i),z_js(i+1),tol,[],z_js(j+1),z_js(i),d,cond);
      F2(j,i) = quad(@f2,z_js(i),z_js(i+1),tol,[],z_js(j+1),z_js(i),d,cond);
      F3(j,i) = quad(@f3,z_js(i),z_js(i+1),tol,[],z_js(j+1),z_js(i),d,cond);
      if cond ~= cond_top     %image technique if conductivity not constant
        F0(j,i) = F0(j,i) + (cond-cond_top)/(cond+cond_top) ...
            *quad(@f0,z_js(i),z_js(i+1),tol,[],-z_js(j+1),d,cond);
        F1(j,i) = F1(j,i) + (cond-cond_top)/(cond+cond_top) ...
            *quad(@f1,z_js(i),z_js(i+1),tol,[],-z_js(j+1),z_js(i),d,cond);
        F2(j,i) = F2(j,i) + (cond-cond_top)/(cond+cond_top) ...
            *quad(@f2,z_js(i),z_js(i+1),tol,[],-z_js(j+1),z_js(i),d,cond);
        F3(j,i) = F3(j,i) + (cond-cond_top)/(cond+cond_top) ...
            *quad(@f3,z_js(i),z_js(i+1),tol,[],-z_js(j+1),z_js(i),d,cond);
      end;
    end;
  end;

  temp_F = F0*E0+F1*E1+F2*E2+F3*E3;  %the F matrix, (N x N+2)

  %   Convert to (N+2xN+2) matrixes by applying the boundary I_0 = I_N+1 = 0.
  F = zeros(N+2);
  F(2:N+1,:) = temp_F(:,:);
  F(1,1) = 1;          %implies I_N+1 = Phi_N+1
  F(N+2,N+2) = 1;      %implies I_N+2 = Phi_N+2
  Fcs = F;
  save(full_filename, 'Fcs','tol');
return;

%Potential functions:
function out0 = f0(zeta,zj,diam,cond)
  out0 = 1./(2.*cond).*(sqrt((diam/2)^2+((zj-zeta)).^2)-abs(zj-zeta));
return;

function out1 = f1(zeta,zj,zi,diam,cond)
  out1 = (zeta-zi).*f0(zeta,zj,diam,cond);
return;

function out2 = f2(zeta,zj,zi,diam,cond)
  out2 = (zeta-zi).^2.*f0(zeta,zj,diam,cond);
return;

function out3 = f3(zeta,zj,zi,diam,cond)
  out3 = (zeta-zi).^3.*f0(zeta,zj,diam,cond);
return;

function [out_zs,my_CSD] = make_cubic_splines(contact_positions,pot,Fcs,my_E0,my_E1,my_E2,my_E3,num_out_zs)
%[out_zs,my_CSD] = make_cubic_splines(contact_positions,pot,Fcs,...
%my_E0,my_E1,my_E2,my_E3,num_out_zs)
%
%Makes the cubic spline function(s).
%
%contact_positions: contact positions
%pot: measured potentials
%Fcs: the cubic spline transformation matrix corresponding to the given
%contact positions
%my_E0,...,my_E3: the matrixes containing the "recursive rules", see paper
%appendix
%num_out_zs: number of out-parameters
%
%Arguments 4 to 8 are optional:
% if nargin<8; num_out_zs = 200; end;
% if nargin<4; %compute E matrixes
%     [my_E0,my_E1,my_E2,my_E3] = compute_Ematrixes(contact_positions);
% end;

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

if nargin<8; num_out_zs = 200; end;
if nargin<4; %compute E matrixes
    [my_E0,my_E1,my_E2,my_E3] = compute_Ematrixes(contact_positions);
end;

%Cubic spline method
[N,num_of_timesteps] = size(pot);
cs_pot = zeros(N+2,num_of_timesteps);
cs_pot(2:N+1,:) = pot(:,:);    % Phi_1 = Phi_N+2 = 0

CSD_coeff = Fcs^(-1)*cs_pot;

%The cubic spline polynomial coeffescients
A0 = my_E0*CSD_coeff;
A1 = my_E1*CSD_coeff;
A2 = my_E2*CSD_coeff;
A3 = my_E3*CSD_coeff;

h = mean(diff(contact_positions));
el_pos_with_ends = zeros(1,length(contact_positions));
el_pos_with_ends(1,1) = 0;
el_pos_with_ends(2:N+1) = contact_positions(1:N);
el_pos_with_ends(1,N+2) = contact_positions(N)+h;

out_zs = el_pos_with_ends(1):(el_pos_with_ends(N+2)-el_pos_with_ends(1))/(num_out_zs-1):el_pos_with_ends(N+2);
i = 1;
for j=1:length(out_zs)
  if out_zs(j)>el_pos_with_ends(i+1)
      i=i+1;
  end;
   my_CSD(j,:) = A0(i,:) + A1(i,:).*(out_zs(j)-el_pos_with_ends(i)) + ...
        A2(i,:)*(out_zs(j)-el_pos_with_ends(i)).^2 + A3(i,:)*(out_zs(j)-el_pos_with_ends(i)).^3;
end;

function [E0,E1,E2,E3] = compute_Ematrixes(el_pos)
%function [E0,E1,E2,E3] = compute_Ematrixes(el_pos)
%
%Computes the E0, E1, E2 and E3 matrixes used in the cubic spline iCSD
%method. These matrixes contains the recursive formulas for finding the F
%matrix (see paper appendix).

%Copyright 2005 Klas H. Pettersen under the General Public License,
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%See: http://www.gnu.org/copyleft/gpl.html

  N = length(el_pos);
  z_js = zeros(1,N+2);            %declare electrode positions included ...
  z_js(1,2:N+1) = el_pos(1,:);    %two imaginary, first electrode in z = 0.
  h_av = sum(diff(el_pos))/(N-1); %average inter-contact distance
  z_js(1,N+2) = z_js(1,N+1)+h_av; %last imaginary electrode position

  C_vec = 1./diff(z_js);          %length: N+1
  % Define transformation matrixes
  C_jm1 = zeros(N+2);
  C_j0 = zeros(N+2);
  C_jall = zeros(N+2);
  C_mat3 = zeros(N+1);

  for i=1:N+1
    for j=1:N+1
      if i == j
        C_jm1(i+1,j+1) = C_vec(i);
        C_j0(i,j) = C_jm1(i+1,j+1);
        C_mat3(i,j) = C_vec(i);
      end;
    end;
  end;
  C_jm1(N+2,N+2) = 0;

  C_jall = C_j0;
  C_jall(1,1) = 1;
  C_jall(N+2,N+2) = 1;

  C_j0(1,1) = 0;

  Tjp1 = zeros(N+2);         %converting an element k_j to k_j+1
  Tjm1 = zeros(N+2);         %converting an element k_j to k_j-1
  Tj0  = eye(N+2);
  Tj0(1,1) = 0;
  Tj0(N+2,N+2) = 0;

  %C to K
  for i=2:N+2
    for j=1:N+2
      if i==j-1
        Tjp1(i,j) = 1;
      end;
      if i==j+1
        Tjm1(i,j) = 1;
      end;
    end;
  end;


  % C to K transformation matrix
  K = (C_jm1*Tjm1+2*C_jm1*Tj0+2*C_jall+C_j0*Tjp1)^(-1)*3*...
    (C_jm1^2*Tj0-C_jm1^2*Tjm1+C_j0^2*Tjp1-C_j0^2*Tj0);

  %   Define matrixes for C to A transformation
  Tja  = zeros(N+1,N+2);      %identity matrix except that it cuts off last elenent
  Tjp1a  = zeros(N+1,N+2);    %converting k_j to k_j+1 and cutting off last element


  %C to A
  for i=1:N+1
    for j=1:N+2
      if i==j-1
        Tjp1a(i,j) = 1;
      end;
      if i==j
        Tja(i,j) = 1;
      end;
    end;
  end;


  %   Define spline coeffiscients
  E0  = Tja;    
  E1  = Tja*K; 
  E2  = 3*C_mat3^2*(Tjp1a-Tja)-C_mat3*(Tjp1a+2*Tja)*K;
  E3  = 2*C_mat3^3*(Tja-Tjp1a)+C_mat3^2*(Tjp1a+Tja)*K;

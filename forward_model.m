% This script creates magnetic data using positions from a field survey.
% The magnetic moment and position of the anomalies is given in the
% Pos_ENV table.
% Errors for sensor positions and magnetic noise can be set in line 33/34.

%% load data
% The files are located in "PROCESSED DATA MATLAB"
load('survey_final.mat')
load('BASEMAG\basemag.mat')

%% Preallocation of variables
points=height(survey.mag_final); % number of data points
XYZ1=[survey.mag_final.N1 survey.mag_final.E1 -survey.mag_final.height1]; % Sensor 1 position matrix
XYZ2=[survey.mag_final.N2 survey.mag_final.E2 -survey.mag_final.height2]; % Sensor 2 position matrix
% GNSS data are recorded with the z axis pointing upward, whereas for
% magnetic data. the z acis points downward. Therefore, the sign of the z
% component in the position matrix is changed.

%% base magnetometer data & anomaly position data preparation
basemag.MAG_POSIXtime_korr=(basemag.MAG_POSIXtime(1):0.01:(basemag.MAG_POSIXtime(1)+(height(basemag)-1)*0.01))'; % correct timestamps of basemag data
Bext=movmean(interp1(basemag.MAG_POSIXtime_korr,basemag.B1,survey.mag_final.MAG_POSIXtime,'linear'),100,1); % interpolate basemag data onto survey timestamps

% rotate external data onto reoriented survey magnetic data
RM=RotMatVecMin(basemag.B1(1,:)/norm(basemag.B1(1,:)),survey.mag_final.B1_Reor(1,:)/norm(survey.mag_final.B1_Reor(1,:)));
Bext=Bext*RM';

% Set coordinates to relative. Line 3 is the South West corner of the
% survey area.
Pos_ENV.Height_abs=Pos_ENV.Height_abs-Pos_ENV.Height_abs(3);
Pos_ENV.North_abs=Pos_ENV.North_abs-Pos_ENV.North_abs(3);
Pos_ENV.East_abs=Pos_ENV.East_abs-Pos_ENV.East_abs(3);
%% ENTER POSITION ERROR AND MAGNETIC NOISE HERE
Pos_err=0.00; % Position error in m
B_err=0; % Magnetic field error in nT

Err_XYZ=rand(size(XYZ1)); % random values for position error
Err_B1=rand(size(XYZ1)); % random values for B1 data 
Err_B2=rand(size(XYZ1)); % random values for B2 data

% The errors in sensor position correspond to error patterns due to
% uncertainties in the GNSS position and the calculation of the sensor
% positions with the help of IMU data.
% The magnetic field error corresponds to noise originating from the UAV,
% the induction effect and vibrations of the boom.
PosErr=movmean(Err_XYZ.*2.*Pos_err-Pos_err,200)*10;
B1err=movmean(Err_B1.*2.*B_err-B_err,200)*10;
B2err=movmean(Err_B2.*2.*B_err-B_err,200)*10;

%% Calculate magnetic field
% First step: external field
B1=Bext;
B2=Bext;

% Second step: calculation of anomaly fields. The resulting field is the
% sum of external field and all anomaly fields
for i=13:height(Pos_ENV)
	B1=B1+MagFeld_Dipol_remanent(Pos_ENV.M(i,:),1,0,[Pos_ENV.North_abs(i) Pos_ENV.East_abs(i) -Pos_ENV.Height_abs(i)],XYZ1+PosErr);
	B2=B2+MagFeld_Dipol_remanent(Pos_ENV.M(i,:),1,0,[Pos_ENV.North_abs(i) Pos_ENV.East_abs(i) -Pos_ENV.Height_abs(i)],XYZ2+PosErr);
end

B1=B1+B1err; % add magnetic noise if present
B2=B2+B2err; % add magnetic noise if present

%% prepare and export data
% The structure of the exported data is equel to the original mag_final.mat
survey.mag_final.B1_Reor=B1;
survey.mag_final.B2_Reor_korr=B2;
survey.mag_final.BT1=sqrt(sum(B1.^2,2));
survey.mag_final.BT2=sqrt(sum(B2.^2,2));
save('sim_survey_final.mat','survey','Pos_ENV','h','lat','lon','Basis_UTM')

%% functions
function B=MagFeld_Dipol_remanent(Mrem,V,chi,Pos,XYZ)
	% 
	% Input:
    	% Mrem : remanent magnetization (x,y,z) [A/m]
    	% V	   : volume [m^3]
    	% Pos  : dipole position
    	% XYZ  : matrix defining the coordinates where the anomaly field is
		%			calculated
	% 
	% 
	% Output:
    	% H : magnetic anomaly [nT]
	%     
	%     
    	
	if ~exist('V','var')||isempty(V)
		V=1; % Set Volume to 1 if no input is given (Mrem is equal to the magnetic moment)
	end
	if ~exist('chi','var')||isempty(chi)
		chi=0; % ignore demagnetization if no susceptibility is given
	end
	if ~exist('Pos','var')||isempty(Pos)
		Pos=[0 0 0]; % set the dipole position to the origin of the coordinate system if no input is given
	end
	
	% check formatting of position vector
	if size(Pos,1)>size(Pos,2)
		Pos=Pos';
	end
	
	% change to dipole-centred coordinate system
	XYZ=XYZ-Pos;
	
	% 1/4*pi(3r(m*r)-m * r2)/r5
	
	Q=V/(4*pi);
	m=Mrem./(1+chi/3); % magnetic moment [Am^2]
	mr=(m*XYZ')';
	rmr=[XYZ(:,1).*mr XYZ(:,2).*mr XYZ(:,3).*mr];
	r_abs=vecnorm(XYZ')';
	mr2=(m'*(r_abs.^2)')';
	r_abs5=r_abs.^5;
	B=Q*(3*rmr-mr2)./r_abs5.*4e2.*pi; % anomaly field [nT]
end
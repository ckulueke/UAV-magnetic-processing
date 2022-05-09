% This script is used to reproduce the plots shown in the article. The
% processing steps performed on the model data are identical to the
% processing steps of the survey data.
%
% In the article, only a part of the survey is shown. To plot the complete
% survey, the marked lines need to be commented out in each plot.

%% load and prepare survey data
% The files are located in "PROCESSED DATA MATLAB"
load('survey_final.mat');

% extract data from struct (cf. "FINAL MATLAB STRUCTURE README.txt")
mag=survey.mag_final;

%% load and prepare model data
sim=load('simsurvey_final.mat', 'survey');
% extract data from struct
mag_sim=sim.survey.mag_final;

mag_sim=mag_sim(~isoutlier(mag_sim.height1),:); % drops in flight altitude are marked as outliers and correspondig data are eliminated

% calculate gradients from reoriented magnetic data
mag_sim.gradB=(mag_sim.B1_Reor-mag_sim.B2_Reor_korr)./0.5...
.*sign(atan2d(mag_sim.B1(:,1),mag_sim.B1(:,2)));
mag_sim.gradBT=(mag_sim.BT1-mag_sim.BT2)./0.5...
.*sign(atan2d(mag_sim.B1(:,1),mag_sim.B1(:,2)));

%% interpolate model data
disp('Starting interpolation...')

% Interpolate gradients to regular grid. A moving median is subtracted from
% the data to account for large-scale and regional spatial variations
% (substitution for tie-line corrections, wavelengths approximately equal
% to profile length)
interp.gradT_sim=scatteredInterpolant([mag_sim.E_Grad mag_sim.N_Grad],...
mag_sim.gradBT-movmedian(mag_sim.gradBT,1000),'natural','none');
interp.gradX_sim=scatteredInterpolant([mag_sim.E_Grad mag_sim.N_Grad],...
mag_sim.gradB(:,1)-movmedian(mag_sim.gradB(:,1),1000),'natural','none');
interp.gradY_sim=scatteredInterpolant([mag_sim.E_Grad mag_sim.N_Grad],...
mag_sim.gradB(:,2)-movmedian(mag_sim.gradB(:,2),1000),'natural','none');
interp.gradZ_sim=scatteredInterpolant([mag_sim.E_Grad mag_sim.N_Grad],...
mag_sim.gradB(:,3)-movmedian(mag_sim.gradB(:,3),1000),'natural','none');

[X_grad,Y_grad]=meshgrid(min(mag.E_Grad):0.1:max(mag.E_Grad),...
		min(mag.N_Grad):0.1:max(mag.N_Grad)); % regular mesh grid with 10 cm spacing

gradT_sim=interp.gradT_sim(X_grad,Y_grad);
gradX_sim=interp.gradX_sim(X_grad,Y_grad);
gradY_sim=interp.gradY_sim(X_grad,Y_grad);
gradZ_sim=interp.gradZ_sim(X_grad,Y_grad);
disp('Gradient data interpolated')


% Interpolate TMI and component data to regular grid. A moving median is 
% subtracted from the data to account for large-scale and regional spatial 
% variations (substitution for tie-line corrections, wavelengths
% approximately equal to profile length)

% Convert table to matrix. For TMI and component data, there is no
% distinction between the sensors, therefore the data are concatenated
% magdata matrix columns: E, N, BE, BN, BV, BT
magdata=[mag_sim.E1 mag_sim.N1 mag_sim.B1_Reor-movmedian(sgolayfilt(mag_sim.B1_Reor,3,51),1000,'omitnan') mag_sim.BT1-movmedian(sgolayfilt(mag_sim.BT1,3,51),1000,'omitnan');...
mag_sim.E2 mag_sim.N2 mag_sim.B2_Reor_korr-movmedian(sgolayfilt(mag_sim.B2_Reor_korr,3,51),1000,'omitnan') mag_sim.BT2-movmedian(sgolayfilt(mag_sim.BT2,3,51),1000,'omitnan')];

interp.TMI_sim=scatteredInterpolant(magdata(:,1:2),magdata(:,6),'natural','none');
interp.X_sim=scatteredInterpolant(magdata(:,1:2),magdata(:,3),'natural','none');
interp.Y_sim=scatteredInterpolant(magdata(:,1:2),magdata(:,4),'natural','none');
interp.Z_sim=scatteredInterpolant(magdata(:,1:2),magdata(:,5),'natural','none');

[X,Y]=meshgrid(min(magdata(:,1)):0.1:max(magdata(:,1)),...
		min(magdata(:,2)):0.1:max(magdata(:,2))); % regular mesh grid with 10 cm spacing

Bx_sim=interp.X_sim(X,Y);
By_sim=interp.Y_sim(X,Y);
Bz_sim=interp.Z_sim(X,Y);
TMI_sim=interp.TMI_sim(X,Y);
disp('TMI and component data interpolated')


% Change anomaly positions to relative. Line 3 denotes the South-West
% corner of the survey area.
Pos_ENV.North_abs=Pos_ENV.North_abs-Pos_ENV.North_abs(3);
Pos_ENV.East_abs=Pos_ENV.East_abs-Pos_ENV.East_abs(3);
Pos_ENV.Height_abs=Pos_ENV.Height_abs-Pos_ENV.Height_abs(3);
Pos_ENV=Pos_ENV(13:end,:); % Lines 1 - 12 denote the base station and corners of the survey area.
% The anomalies shown in the article are anomalies 1 - 6 (lines 13 - 18).
disp('PROGRAM FINISHED')

%% plot interpolated gradient maps
% limits used in the article. Only anomalies 1 - 6 are shown
limx=[24 54];
limy=[-3 16];
lims(1)=duration(0,36,37);
lims(2)=duration(0,36,45);

% Interpolated TMI gradient map
f(1)=figure(1);clf;
	f(1).Name='Interpolated TMI gradient map';
	f(1).Units='centimeters';
	f(1).Position=[0 3 16 22];
	s1(1)=subplot(2,1,1); % survey
		pcolor(X_grad,Y_grad,asinh(survey.gradT))
		hold on
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		clim=caxis;
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated TMI gradient map (data)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'a)')
	s1(2)=subplot(2,1,2); % model
		pcolor(X_grad,Y_grad,asinh(gradT_sim))
		hold on
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		% -------- optional: plot positions of data points --------------
		% 		scatter(mag_sim.E1,mag_sim.N1,10,'k.')
		% 		scatter(mag_sim.E2,mag_sim.N2,10,'k.')
		% ---------------------------------------------------------------
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		h=asinhcolorbar('\partial_x|B| [nT/m]');
 		h.Location='southoutside';
		caxis(clim);
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated TMI gradient map (model)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'b)')
	s1(1).Position=[0.13,0.58,0.775,0.35];
	s1(1).InnerPosition=[0.13,0.58,0.775,0.35];
	s1(1).OuterPosition=[0 0.53,1,0.42];
	s1(2).Position=[0.13,0.13,0.775,0.35];
	s1(2).InnerPosition=[0.13,0.13,0.775,0.35];
	s1(2).OuterPosition=[0 0.09,1,0.42];

% Interpolated B_N gradient map
f(2)=figure(2);clf;
	f(2).Name='Interpolated B_N gradient map';
	f(2).Units='centimeters';
	f(2).Position=[0 3 16 22];
	s2(1)=subplot(2,1,1); % survey
		pcolor(X_grad,Y_grad,asinh(survey.gradX))
		hold on
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		clim=caxis;
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated B_N gradient map (data)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'a)')
	s2(2)=subplot(2,1,2); % model
		pcolor(X_grad,Y_grad,asinh(gradX_sim))
		hold on
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		h=asinhcolorbar('\partial_xB_N [nT/m]');
 		h.Location='southoutside';
		caxis(clim);
		xlabel('Easting in m')
		ylabel('Northing in m')
		title('Interpolated B_N gradient map (model)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'b)')
	s2(1).Position=[0.13,0.58,0.775,0.35];
	s2(1).InnerPosition=[0.13,0.58,0.775,0.35];
	s2(1).OuterPosition=[0 0.53,1,0.42];
	s2(2).Position=[0.13,0.13,0.775,0.35];
	s2(2).InnerPosition=[0.13,0.13,0.775,0.35];
	s2(2).OuterPosition=[0 0.09,1,0.42];

% Interpolated B_E gradient map
f(3)=figure(3);clf;
	f(3).Name='Interpolated B_E gradient map';
	f(3).Units='centimeters';
	f(3).Position=[0 3 16 22];
	s3(1)=subplot(2,1,1); % survey
		pcolor(X_grad,Y_grad,asinh(survey.gradY))
		hold on
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		clim=caxis;
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated B_E gradient map (data)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'a)')
	s3(2)=subplot(2,1,2); % model
		pcolor(X_grad,Y_grad,asinh(gradY_sim))
		hold on
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		xlim(limx);
		ylim(limy);
		h=asinhcolorbar('\partial_xB_N [nT/m]');
 		h.Location='southoutside';
		caxis(clim);
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated B_E gradient map (model)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'b)')
	s3(1).Position=[0.13,0.58,0.775,0.35];
	s3(1).InnerPosition=[0.13,0.58,0.775,0.35];
	s3(1).OuterPosition=[0 0.53,1,0.42];
	s3(2).Position=[0.13,0.13,0.775,0.35];
	s3(2).InnerPosition=[0.13,0.13,0.775,0.35];
	s3(2).OuterPosition=[0 0.09,1,0.42];

% Interpolated B_V gradient map
f(4)=figure(4);clf;
	f(4).Name='Interpolated B_V gradient map';
	f(4).Units='centimeters';
	f(4).Position=[0 3 16 22];
	s4(1)=subplot(2,1,1); % survey
		pcolor(X_grad,Y_grad,asinh(survey.gradZ))
		hold on
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		clim=caxis;
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated B_V gradient map (data)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'a)')
	s4(2)=subplot(2,1,2); % model
		pcolor(X_grad,Y_grad,asinh(gradZ_sim))
		hold on
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		h=asinhcolorbar('\partial_xB_V [nT/m]');
 		h.Location='southoutside';
		caxis(clim);
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated B_V gradient map (model)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'b)')
	s4(1).Position=[0.13,0.58,0.775,0.35];
	s4(1).InnerPosition=[0.13,0.58,0.775,0.35];
	s4(1).OuterPosition=[0 0.53,1,0.42];
	s4(2).Position=[0.13,0.13,0.775,0.35];
	s4(2).InnerPosition=[0.13,0.13,0.775,0.35];
	s4(2).OuterPosition=[0 0.09,1,0.42];

%% plot interpolated TMI anomaly map
f(5)=figure(5);clf;
	f(5).Name='Interpolated TMI map';
	f(5).Units='centimeters';
	f(5).Position=[0 3 16 22];
	s5(1)=subplot(2,1,1); % survey
		pcolor(X,Y,asinh(survey.TMI))
		hold on
		scatter(mag.E1,mag.N1,10,'w.')
		scatter(mag.E2,mag.N2,10,'w.')
		scatter(mag.E1(mag.time-mag.time(1)>lims(1) & mag.time-mag.time(1)<lims(2)) ,mag.N1(mag.time-mag.time(1)>lims(1) & mag.time-mag.time(1)<lims(2)),20,'k.')
		scatter(mag.E2(mag.time-mag.time(1)>lims(1) & mag.time-mag.time(1)<lims(2)),mag.N2(mag.time-mag.time(1)>lims(1) & mag.time-mag.time(1)<lims(2)),20,'k.')
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		clim=caxis;
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated TMI map (data)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'a)')
	s5(2)=subplot(2,1,2); % model
		pcolor(X,Y,asinh(TMI_sim))
		hold on
		scatter(mag.E1,mag.N1,10,'w.')
		scatter(mag.E2,mag.N2,10,'w.')
		scatter(Pos_ENV.East_abs(1:6),Pos_ENV.North_abs(1:6),'k*')
		text(Pos_ENV.East_abs(1:6)+1,Pos_ENV.North_abs(1:6)+1,num2str(Pos_ENV.AnomalyLabel(1:6)))
		axis equal
		shading interp
		colormap jet
		% ------- to plot the whole survey, comment out these lines -------
		xlim(limx);
		ylim(limy);
		% -----------------------------------------------------------------
		h=asinhcolorbar('|B| [nT]');
 		h.Location='southoutside';
		caxis(clim);
		xlabel('Easting [m]')
		ylabel('Northing [m]')
		title('Interpolated TMI map (model)')
		text(limx(1)+(limx(2)-limx(1))*0.02,limy(1)+(limy(2)-limy(1))*0.05,'b)')
	s5(1).Position=[0.13,0.58,0.775,0.35];
	s5(1).InnerPosition=[0.13,0.58,0.775,0.35];
	s5(1).OuterPosition=[0 0.53,1,0.42];
	s5(2).Position=[0.13,0.13,0.775,0.35];
	s5(2).InnerPosition=[0.13,0.13,0.775,0.35];
	s5(2).OuterPosition=[0 0.09,1,0.42];

%% plot gradient timeseries
f(10)=figure(10);clf;
	s10(1)=subplot(4,1,1);
		plot(mag.time-mag.time(1),mag.gradB(:,1)-mag_sim.gradB(:,1))
		grid minor
		ylabel('\partial_x B_N [nT/m]')
		xlabel('time [hh:mm:ss]')
		title('Gradient residuals (data-model)')
	s10(2)=subplot(4,1,2);
		plot(mag.time-mag.time(1),mag.gradB(:,2)-mag_sim.gradB(:,2))
		grid minor
		ylabel('\partial_x B_E [nT/m]')
		xlabel('time [hh:mm:ss]')
	s10(3)=subplot(4,1,3);
		plot(mag.time-mag.time(1),mag.gradB(:,3)-mag_sim.gradB(:,3))
		grid minor
		ylabel('\partial_x B_V [nT/m]')
		xlabel('time [hh:mm:ss]')
	s10(4)=subplot(4,1,4);
		plot(mag.time-mag.time(1),mag.gradBT-mag_sim.gradBT)
		grid minor
		ylabel('\partial_x BT [nT/m]')
		xlabel('time [hh:mm:ss]')
	linkaxes(s10,'x')
f(11)=figure(11);clf;
	s11(1)=subplot(3,1,1);
		plot(mag.time-mag.time(1),mag.gradB(:,1),'k','LineWidth',1.5)
		hold on
		plot(mag.time-mag.time(1),mag_sim.gradB(:,1),'g','LineWidth',1.5)
		grid minor
		ylabel('\partial_x B_N [nT/m]')
		xlabel('time [hh:mm:ss]')
		title('Gradient timeseries')
		legend('measured','model')
	s11(2)=subplot(3,1,2);
		plot(mag.time-mag.time(1),mag.gradB(:,2),'k','LineWidth',1.5)
		hold on
		plot(mag.time-mag.time(1),mag_sim.gradB(:,2),'g','LineWidth',1.5)
		grid minor
		ylabel('\partial_x B_E [nT/m]')
		xlabel('time [hh:mm:ss]')
	s11(3)=subplot(3,1,3);
		plot(mag.time-mag.time(1),mag.gradB(:,3),'k','LineWidth',1.5)
		hold on
		plot(mag.time-mag.time(1),mag_sim.gradB(:,3),'g','LineWidth',1.5)
		grid minor
		ylabel('\partial_x B_V [nT/m]')
		xlabel('time [hh:mm:ss]')
	linkaxes(s11,'x')
	xlim(lims)
f(11).Units='centimeters';
f(11).PaperUnits='centimeters';
f(11).Position=[0 0 16 9];
f(11).PaperPosition=[0 0 16 9];

%% plot TMI timeseries
f(12)=figure(12);clf;
	s12(1)=subplot(2,1,1);
		plot(mag.time-mag.time(1),mag.BT1,'k','LineWidth',1.5)
		hold on
		plot(mag.time-mag.time(1),mag_sim.BT1,'g','LineWidth',1.5)
		grid minor
		ylabel('BT1 [nT]')
		xlabel('time [hh:mm:ss]')
		title('TMI timeseries')
		legend('measured','model')
		ax=gca;
		ax.YAxis.Exponent=0;
	s12(2)=subplot(2,1,2);
		plot(mag.time-mag.time(1),mag.BT2,'k','LineWidth',1.5)
		hold on
		plot(mag.time-mag.time(1),mag_sim.BT2,'g','LineWidth',1.5)
		grid minor
		ylabel('BT2 [nT]')
		xlabel('time [hh:mm:ss]')
		ax=gca;
		ax.YAxis.Exponent=0;
	linkaxes(s12,'x')
	xlim(lims)
f(12).Units='centimeters';
f(12).PaperUnits='centimeters';
f(12).Position=[0 0 16 9];
f(12).PaperPosition=[0 0 16 9];
%% load data
% Script for reducing vibration-related (high-frequency) noise. Can be 
% performed on reoriented and non-reoriented data. If performed on
% non-reoriented data, exchange lines 16/17 and 27/28 

load("survey_final.mat",'survey') % load data from matlab workscape file
mag=survey.mag_final; % extracting magnetic data table from struct
%% Preallocation
% Preallocation of Roll, Yaw and Pitch vectors denoting the relative angles
% between the two sensors
Roll=zeros(height(mag),1);
Yaw=Roll;
Pitch=Roll;
%% Calculate relative angles
for i=1:height(mag)
	[Roll(i),Yaw(i),Pitch(i)]=RYP(RotMatVecMin(mag.B2_Reor(i,:),mag.B1_Reor(i,:))); % select this line when performing the correction on reoriented data
% 	[Roll(i),Yaw(i),Pitch(i)]=RYP(RotMatVecMin(mag.B2_Reor(i,:),mag.B1(i,:))); % select this line when performing the correction on non-reoriented data
end
%% highpass filter
% angles are high-pass filtered by removing low-pass filtered
% angles from calculated angles 
Roll=Roll-sgolayfilt(Roll,3,51);
Yaw=Yaw-sgolayfilt(Yaw,3,51);
Pitch=Pitch-sgolayfilt(Pitch,3,51);
%% correction of Sensor 2 reoriented data
for i=1:height(mag)
	mag.B2_Reor_korr(i,:)=(RotMatRYP(Roll(i),Yaw(i),Pitch(i))*mag.B2_Reor(i,:)')'; % select this line when performing the correction on reoriented data
% 	mag.B2_korr(i,:)=(RotMatRYP(Roll(i),Yaw(i),Pitch(i))*mag.B2(i,:)')'; % select this line when performing the correction on non-reoriented data
end
% END OF PROGRAM
%% functions
function [Roll,Yaw,Pitch]=RYP(A)
	% Calculate Roll, Yaw and Pitch angles from rotation matrix A
	% Input: rotation matrix A
	% Output: Roll, Yaw, Pitch angles (in degrees)
	
	Pitch=atan2d(-A(3,1),sqrt(A(1,1)^2+A(2,1)^2));
	Yaw=atan2d(A(2,1)/cos(Pitch),A(1,1)/cos(Pitch));
	Roll=atan2d(A(3,2)/cos(Pitch),A(3,3)/cos(Pitch));
end

function [RotMatVecMin,a,phi]=RotMatVecMin(v,w)
	% Calculate ratation matrix that rotates vector v onto vector w using
	% the smallest possible rotation angle
	
	% v,w: vectors
	% a: unit vector denoting the rotation axis
	% phi: rotation angle in degrees
	
	if v./norm(v)==w./norm(w) % no rotation needed if vectors are parallel
    	a=[0. 0. 0.];
    	phi=0.;
    	RotMatVecMin=eye(3);
    	
	else
    	a=cross(v,w)/norm(cross(v,w));
    	phi=acosd(dot(v,w)/norm(v)/norm(w));
	
    	RotMatVecMin=[[cosd(phi)+a(1)^2*(1-cosd(phi)) a(1)*a(2)*(1-cosd(phi))-a(3)*sind(phi) a(1)*a(3)*(1-cosd(phi))+a(2)*sind(phi)];...
               	[a(1)*a(2)*(1-cosd(phi))+a(3)*sind(phi) cosd(phi)+a(2)^2*(1-cosd(phi)) a(2)*a(3)*(1-cosd(phi))-a(1)*sind(phi)];...
               	[a(1)*a(3)*(1-cosd(phi))-a(2)*sind(phi) a(2)*a(3)*(1-cosd(phi))+a(1)*sind(phi) cosd(phi)+a(3)^2*(1-cosd(phi))]];
	end
end
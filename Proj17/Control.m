% This script sets the controller parameters for the SLS 3-D Printer

% ================
% CONTROLLER GAINS
% ================

% Enter optimized PID values here.
% No more than 3 significant figures per gain value.

Kd0 = round(0.245, 3, 'significant');
Kp0 = round(0.9248*Kd0, 3, 'significant');
Ki0 = round(45.9*Kd0, 3, 'significant');

Kd1 = round(0.00321, 3, 'significant');
Kp1 = round(46.69*Kd1, 3, 'significant');
Ki1 = round(0, 3, 'significant');

PID0 = [Kp0 Ki0 Kd0]; %[Kp Ki Kd]
PID1 = [Kp1 Ki1 Kd1]; %[Kp Ki Kd]

% Enter feedback sensor values here.

FB0 = 1/Sens0;
FB1 = 1/Sens1;


% =====================
% Set-Point Time Vector
% =====================

% The Time Vector is stored in a variable called "Time".
% It's initial value is equally spaced for all samples and is
% set in TRAJECTORY.M
%
% Redefine this vector here to optimize the build time of the part.
% You can define it analytically or type in the elements manually but
% you must not change the length because it must contain one value for
% each Xd/Yd position pair.
% In the Matlab window, enter "length(Time)" to see how big it is.

% The Time vector must range from 0 to TotalTime

%Time       = 0:SampleTime:TotalTime; % DO NOT CHANGE TotalTime


%%
%Calculating transfer function
%Integrator (For use later)
intr = tf(1, [1 0]);

%For Q0
%Existing transfer functions
AmpTF0 = tf(Amp0n, Amp0d);
ElecTF0 = tf(Elec0n, Elec0d);
MechTF0 = tf(Mech0n, Mech0d);

%Feedback loop in middle
motorG0 = ElecTF0*TConst0*MechTF0;
motorH0 = BackEMF0;
motorTF0 = feedback(motorG0, motorH0);

G0 = AmpTF0*motorTF0*intr;
G0 = zpk(minreal(G0)); %For ease of viewing
H0 = Sens0*FB0;

%Kd0 = round(1, 3, 'significant');
%Kp0 = round(0.9248*Kd0, 3, 'significant');
%Ki0 = round(45.9*Kd0, 3, 'significant');
pid0 = pid(Kp0, Ki0, Kd0); %tfpid0 = pid(PID0, [1 0]);
olPID0 = pid0*G0;
tf0 = feedback(G0, H0); %Transfer function of whole system, no PID
tfpid0 = feedback(G0*pid0, H0); %Transfer function of whole system, PID

%rlocus(olPID1)
%axis([-60 0 -10 10])
%step(tfpid0)
%stepinfo(tfpid0)

%For Q1
%Existing transfer functions
AmpTF1 = tf(Amp1n, Amp1d);
ElecTF1 = tf(Elec1n, Elec1d);
MechTF1 = tf(Mech1n, Mech1d);

%Feedback loop in middle
motorG1 = ElecTF1*TConst1*MechTF1;
motorH1 = BackEMF1;
motorTF1 = feedback(motorG1, motorH1);

G1 = AmpTF1*motorTF1*intr;
G1 = zpk(minreal(G1)); %For ease of viewing
H1 = Sens1*FB1;

%Kd1 = round(0.003, 3, 'significant');
%Kp1 = round(46.69*Kd1, 3, 'significant');
%Ki1 = round(0, 3, 'significant');
pid1 = pid(Kp1, Ki1, Kd1); %tfpid0 = pid(PID0, [1 0]);
olPID1 = pid1*G1;
tf1 = feedback(G1, H1); %Transfer function of whole system, no PID
tfpid1 = feedback(G1*pid1, H1); %Transfer function of whole system, PID

%rlocus(olPID1)
%axis([-60 0 -10 10])
%step(tfpid1)
%stepinfo(tfpid1)
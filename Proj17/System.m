% This script sets the model parameters for the SLS 3-D Printer

% Example: Specifying a Dynamics Block
% n = [1 2 3];
% d = [4 5 6];
% Transfer Function = (s^2 + 2s + 3) / (4s^2 + 5s + 6)

% ========================
% PHYSICAL UNIT CONVERSION
% ========================
% Example: if you decide to work in (Kg), all masses must be represented
%          in (Kg) but the spec sheet may provide masses in (g)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values you should over-write in Control.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ================
% CONTROLLER GAINS
% ================
% PID Controller
PID0 = [1 0 0];
PID1 = [1 0 0];

% Controller Feedback Gain
FB0  = 0;
FB1  = 0;

% Assign the Set-Point Trajectory & Time Vector
TRAJECTORY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values you should over-write in System.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==========================
% Choose Motors
% ==========================
AMAX22_5W_SB;                % Default Maxon motor
Q0 = MotorParam;
Q1 = MotorParam;

% ==========================
% Motor Parameters
% ==========================
NomI0   = 1;                 % Max average current
StallI0 = 1;                 % Max peak current
NomI1   = 1;                 % Max average current
StallI1 = 1;                 % Max peak current


% ==========================
% Q0 - Rotation about y-axis
% ==========================

% Amplifier Dynamics
Amp0n   = [49.1383];               % Numerator
Amp0d   = [1 49.1667];               % Denominator
AmpSat0 = Big;

% Electrical Motor Dynamics
% I/V = 1/Z = 1/(sL+R) = 1/ (0.362L+5.39)
Elec0n  = [1];               % Numerator
Elec0d  = [0.000362 5.39];      % Denominator

% Torque Const & Back EMF
TConst0  = 1;
BackEMF0 = 0;

% Mechanical Motor Dynamics
Mech0n  = [1];               % Numerator
Mech0d  = [1];               % Denominator
JntSat0 =  Big;

% Sensor Dynamics
Sens0    =  0;
SensSat0 =  Big;

% Static Friction
StFric0 = 0;


% ==========================
% Q1 - Rotation about x-axis
% ==========================

% Amplifier Dynamics
Amp1n   = [1];               % Numerator
Amp1d   = [1];               % Denominator
AmpSat1 = Big;

% Electrical Motor Dynamics
Elec1n = [1];                % Numerator
Elec1d = [1];                % Denominator

% Torque Const & Back EMF
TConst1  = 1;
BackEMF1 = 0;

% Mechanical Motor Dynamics
Mech1n = [1];                % Numerator
Mech1d = [1];                % Denominator
JntSat1 =  Big;

% Sensor Dynamics
Sens1    =  0;
SensSat1 =  Big;

% Static Friction
StFric1 = 0;


% ==================
% TRANSFER FUNCTIONS
% ==================
% Compute transfer functions from above values and perform system analysis
% You may prefer to put this section in a separate .m file

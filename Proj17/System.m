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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Over-write the default values from DEFAULT.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==========================
% Choose Motors
% ==========================

% Motor Unit Conversions
% ----------------------


% ==========================
% Motor Parameters
% ==========================

% Maximum Current
% ---------------


% =============================
% Q0 : Rotation about y-axis
% =============================

% Amplifier Dynamics
% ------------------

% Electrical Motor Dynamics
% -------------------------

% Torque Const & Back EMF
% -----------------------

% Mechanical Motor Dynamics
% -------------------------

% Sensor Dynamics
% ---------------

% Static Friction
% ---------------


% =============================
% Q1 : Rotation about x-axis
% =============================

% Amplifier Dynamics
% ------------------

% Electrical Motor Dynamics
% -------------------------

% Torque Const & Back EMF
% -----------------------

% Mechanical Motor Dynamics
% -------------------------

% Sensor Dynamics
% ---------------

% Static Friction
% ---------------


% ==================
% TRANSFER FUNCTIONS
% ==================
% Compute transfer functions from above values and perform system analysis
% You may prefer to put this section in a separate .m file

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

%Ring Mass Calculations

%Converting units to SI units
LinkR1M = LinkR1/1000; %m
LinkR2M = LinkR2/1000; %m
LinkDM = LinkD/1000; %m
LinkOffM = LinkOff/1000; %m

%Calculating inner and outer volumes, then total ring for finding mass
RingOuterV = pi*(LinkR2M^2)*(LinkDM); %m^3
RingInnerV = pi*(LinkR1M^2)*(LinkDM); %m^3
RingTotalV = RingOuterV-RingInnerV; %m^3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Over-write the default values from DEFAULT.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==========================
% Choose Motors
% ==========================

%Default motors for both
AMAX22_6W_SB;
Q0 = MotorParam;
AMAX16_2W_SB;
Q1 = MotorParam;

% Motor Unit Conversions
% ----------------------

%Converting everything to SI units
MotorMass0KG = Q0(Weight)/1000; %kg
MotorMass1KG = Q1(Weight)/1000; %kg
MotorRadius0M = Q0(OuterDiam)/2000; %m
MotorRadius1M = Q1(OuterDiam)/2000; %m
MotorLength0M = Q0(Length)/1000; %m
MotorLength1M = Q1(Length)/1000; %m

% ==========================
% Motor Parameters
% ==========================

% Maximum Current
% ---------------

NomI0   = Q0(NomCurr);        % Max average current, A
StallI0 = Q0(StallCurr);      % Max peak current, A
NomI1   = Q1(NomCurr);        % Max average current, A
StallI1 = Q1(StallCurr);      % Max peak current, A

% =============================
% Q0 : Rotation about y-axis
% =============================

% Amplifier Dynamics
% ------------------

%Transfer function obtained using MNA

%All amplifier values are converted to prefix-normalized SI units L, C, R
Amp0n = [(C*10^(-6))*(R1*10^6)*R2-L/1000]; %Putting constants into MNA result
Amp0d = [(L/1000)*R1*(10^6)*(C*10^(-6)) (C*10^(-6))*R2*R1*(10^6)]; %Putting constants into transfer function
AmpSat0 = Q0(NomV); %Amplifier specified to not exceed motor nominal voltage, V

% Electrical Motor Dynamics
% -------------------------

Elec0n = [1];
Elec0d = [Q0(TermL)/1000 Q0(TermR)]; %Electrical side is inductance (H)*s + Resistance (ohms)

% Torque Const & Back EMF
% -----------------------

TConst0 = Q0(TorqueK)/1000; %Nm/A
BackEMF0 = 1/(Q0(SpdK)*RadPSecPerRPM); %Using speed constant, find voltage of back EMF (V/(rad/s))

% Mechanical Motor Dynamics
% -------------------------

%JTotal = JMotor + JQ1 + JCounterWeight + JRing

JM = Q0(RotJ)/10000000; %Inertia of the motor in kg*m^2 (Given)

%Finding the moment of inertia for Q1 using superposition
%Radius of motor is the same for superposition, find density in terms of kg/m
%instead of full density for simplicity
MotorDensity1 = MotorMass1KG/MotorLength1M; %kg/m (radius constant, variable length)

%Finding Q1 and counterweight inertia using superposition
%Finding equivalent mass of both 'pieces' for superposition
OuterMass = MotorDensity1*(MotorLength1M+LinkOffM); %kg
InnerMassM = MotorDensity1*(LinkOffM); %kg
%Using equations for moments of inertia, calculate the moment of
%inertia for the larger and smaller piece
JMout = (1/3)*(OuterMass*((MotorLength1M+LinkOffM)^2)); %kgm^2
JMin = (1/3)*InnerMassM*(LinkOffM^2); %kgm^2
JQ1 = JMout-JMin; %Superposition, kgm^2
JCW = JQ1; %Assuming motor and counterweight are physically the same, kgm^2

%Finding ring inertia using superposition
%First, finding mass using volume
%Calculations done modeling ring as a solid aluminum cylinder cut into a
%ring, simplifying portions that hold Q1 and counterweight
OuterMassR = RingOuterV*RhoAl*1000; %kg, for finding inertia
InnerMass = RingInnerV*RhoAl*1000; %kg, for finding inertia
%Then use these masses to calculate the total inertia of the ring
Jout = (3*OuterMassR*(LinkR2M^2)+OuterMassR*(LinkDM^2))/12; %kgm^2
Jin = (3*InnerMass*(LinkR1M^2)+InnerMass*(LinkDM^2))/12; %kgm^2
JR = Jout-Jin; %kgm^2

JTotal = JM + JQ1 + JCW + JR; %Summing up moments of inertia, kgm^2

%The internal friction of the motor shaft is represented by the B value at
%the no load speed. B = 1/((Speed Torque Gradient)*1000*((rad/s)/rpm))
%The bearing opposite Q0 is represented by multiplying this B value by 2.
DynamF = 2*(((Q0(NoLoadCurr)*Q0(TorqueK))/(1000000*Q0(NoLoadSpd)*RadPSecPerRPM))); %Nm/(rad/s)

K = SpringK/(1000*2*pi); %Given in project document, Nm/rad (corrrecting to SI units)

Mech0n = [1 0]; %Note that term in denominator puts a zero in the numerator
Mech0d = [JTotal DynamF K]; %All units in SI units as above

JntSat0 = Big; %Range of motion of sensor is unlimited for Q0 as per project document

% Sensor Dynamics
% ---------------

%Derived from information given in project document
%Converting sensor range to degrees
SensRange = SensAng*RadPerDeg; %Rad
Sens0 = SensV/SensRange; %V/Rad
%Sensor outputs -5V at -180 degrees, 5V at 180 degrees, linear in between
SensSat0 = SensV; %Sensor saturates at max given voltage, V       

% Static Friction
% ---------------

%Finding mass then weight of everything Q0 supports
%Mass is sum of Q0, equivalent counterweight, and ring
RingMass = (RingTotalV*RhoAl)*1000; %kg
Q0HeldMass0 = Q0(Weight)/1000+Q0(Weight)/1000+RingMass; %kg
Weight0 = G*Q0HeldMass0; %Calculating weight in N using gravity coefficient

StFric0 = (Weight0*uSF/1000000); %N

% =============================
% Q1 : Rotation about x-axis
% =============================

% Amplifier Dynamics
% ------------------

%Transfer function obtained using MNA

%All amplifier values are converted to prefix-normalized SI units L, C, R
Amp1n = [(C*10^(-6))*(R1*10^6)*R2-L/1000];
Amp1d = [(L/1000)*R1*(10^6)*(C*10^(-6)) (C*10^(-6))*R2*R1*(10^6)]; %Putting in values for transfer function
AmpSat1 = Q0(NomV); %Amplifier specified to not exceed motor nominal voltage

% Electrical Motor Dynamics
% -------------------------

Elec1n = [1];
Elec1d = [Q1(TermL)/1000 Q1(TermR)]; %Electrical side is inductance (H)*s + Resistance (ohms)

% Torque Const & Back EMF
% -----------------------

TConst1 = Q1(TorqueK)/1000; %Nm/A
BackEMF1 = 1/(Q1(SpdK)*RadPSecPerRPM); %Using speed constant, find voltage of back EMF (V/(rad/s))

% Mechanical Motor Dynamics
% -------------------------

%JTotal = J1 as outlined in project document

J1 = Q1(RotJ)/10000000; %Q1 has no inertia besides rotor friction, kg*m^2

%The internal friction of the motor shaft is represented by the B value at
%the no load speed. B = 1/((Speed Torque Gradient)*1000*((rad/s)/rpm))
%No bearing is mentioned, so only the motor dynamic friction is considered.
B1 = (((Q1(NoLoadCurr)*Q1(TorqueK))/(1000000*Q1(NoLoadSpd)*RadPSecPerRPM))); %Nm/(rad/s)

Mech1n = [1];
Mech1d = [J1 B1]; %Notable that Q1 has no return spring, therefore no K and 1st order denominator

JntSat1 = JntLim*RadPerDeg; %As per controller, the jnt sat checker takes in radians not degrees, rad
%Joint is considered saturated when it reaches max angle

% Sensor Dynamics
% ---------------

%Converting sensor range to degrees
SensRange = SensAng*RadPerDeg; %rad
Sens1 = SensV/SensRange; %V/rad
%Derived from the info given on sensor
%Sensor outputs -5V at -180 degrees, 5V at 180 degrees, linear in between
SensSat1 = SensV; %Sensor saturates at max given voltage, V

% Static Friction
% ---------------

StFric1 = 0; %Q1 has no significant static friction as laser is weightless

% ==================
% TRANSFER FUNCTIONS
% ==================
% Compute transfer functions from above values and perform system analysis
% You may prefer to put this section in a separate .m file

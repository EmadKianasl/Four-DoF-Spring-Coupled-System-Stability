function [DATA,AxisRange]=demoparm(system)
%DEMOPARM    Parameters for demo systems
%            The demo systems include:
%             1) Logistic map
%             2) Henon map
%             3) Duffing's equation
%             4) Lorenz equation
%             5) Rossler equation
%             6) Van Der Pol equation
%             7) Stewart-McCumber model

%            by Steve W. K. SIU, July 5, 1998.

%-------Common parameters--------
output=0;                   %Don't check "Output File": 1="check", 0="uncheck"
LEout=0;                    %Don't check "Lyapunov Exponents"
ODEout=0;                   %Don't check "Lyapunov Dimension"
LEprecision=1;              %Precision of output values of the 
ODEprecision=1;             %       Lyapunov exponents and dimension
                            %       1="%.4f", 2="%.6f', ..., 5=".12f"
%Line Colors
Blue=1; Black=2; Green=3; Red=4; Yellow=5; Magenta=6; Cyan=7;
LineColor=Blue;              %  line color: Blue


switch system

case 'parallelbars'
   
   %Parameters for Rossler-hyperchaos
   IntMethod=2;             %Integration method: 1=Discrete map, 2=ODE45, 3=ODE23
                            % 4=ODE113, 5=ODE23S, 6=ODE15S
   InitialTime=0;           %Initial time: 0
   FinalTime=10;          %Final Time: 10 sec 
   TimeStep=0.001;           %Time step: 0.01 sec
   RelTol=1e-5;             %Relative tolerance
   AbsTol=1e-6;             %Absolute tolerance
   IC=[0 0 0 0 0 0 0 0];              %Initial coniditions
   LODEnum=64;              %No. of linearized ODEs
   
   %PLOTTING OPTIONS: 	Only one of them can be set "on" (i.e. 1)
   plot1=0;                 %Plot immediately
   plot2=1;                 %Plot every  ItrNum iterations
   ItrNum=5;				
   
   
   Discard=50;             %Transient iterations to be discarded: 
                            %   50 Iterations = 50*10 time steps = 20 sec
   UpdateSteps=5;          %Update the LEs every 5 time steps
   
   %Axis range for plotting
   AxisRange=[InitialTime,FinalTime,-20,20];
   

otherwise
   error('Invalid system!')
end

%Save the parameters in a matrix 
DATA=[   output,       LEout, LEprecision, ODEout,  ODEprecision, ...
      IntMethod, InitialTime,   FinalTime, TimeStep,      RelTol, ...
         AbsTol,       plot1,       plot2, ItrNum,     LineColor, ...
        Discard, UpdateSteps,     LODEnum,     IC];
    
   


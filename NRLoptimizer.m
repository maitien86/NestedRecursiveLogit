%   MAIN PROGRAM
%   ---------------------------------------------------
Credits;
globalVar;

global resultsTXT; 
file_linkIncidence = './Input/linkIncidence.txt';
file_AttEstimatedtime = './Input/TravelTime.txt';
file_turnAngles = './Input/TurnAngles.txt';
file_observations = './Input/observationsAll.txt';

%% Set estimation type
global isFixedMu;
isFixedMu = 0;
SampleObs  = [];

%% Initialize the optimizer structure

isLinkSizeInclusive = false;
isFixedUturn = false;
loadData;
Op = Op_structure;
initialize_optimization_structure();

if isLinkSizeInclusive ==  true
    Op.x = [-2.139;-0.748;-0.224;-3.301;-0.155;0.341;-0.581;-0.092];
else
    Op.x = [-2.139;-0.748;-0.224;-3.301;0.341;-0.581;-0.092];
end
%% Relax Att
getAtt();
for i = Op.m+1: Op.n
        u = sparse((size(incidenceFull)));
        Atts(i).value = (u);
end

Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
Op.Hessian_approx = OptimizeConstant.BHHH;
Gradient = zeros(nbobs,Op.n);

%% Optimizing ... 
fprintf('\n Estimating ... \n');
tic;
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','GradObj','on');
[x,fval,exitflag,output,grad] = fminunc(@LL,Op.x,options);
Op.value = fval;
Op.x = x;
Op.grad = grad;

getCov;
%Finishing ...
ElapsedTtime = toc
resultsTXT = [resultsTXT sprintf('\n Number of function evaluation %d \n', Op.nFev)];
resultsTXT = [resultsTXT sprintf('\n Estimated time %d \n', ElapsedTtime)];
%% Send email nitification 
try
   notifyMail('send', resultsTXT);
catch exection
   fprintf('\n Can not send email notification !!! \n');
end


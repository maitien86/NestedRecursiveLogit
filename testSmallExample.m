....
%   Linked based network route choice model with unrestricted choice set
%   Optmization algorithm
%   MAI ANH TIEN - DIRO
%   29 - July - 2013
%   MAIN PROGRAM
%   ---------------------------------------------------
%%
Credits;
% Initialize email notification
globalVar;
global resultsTXT; 
notifyMail('set','amyeuphich@gmail.com','sntal2908');


%------------------------------------------------------

file_linkIncidence = './ExampleNested/incidence.txt';
file_AttEstimatedtime = './ExampleNested/time.txt';
file_turnAngles = './ExampleNested/turnAngle.txt';
file_observations = './ExampleNested/obs.txt';
%file_observations = './simulatedData/ObservationsAll_NoLS.txt';

% Initialize the optimizer structure
global isFixedMu;
isFixedMu = 1;
global mu;
mu = ones(10,1);
mu(2) = 0.5;
mu(3) = 0.8;
isLinkSizeInclusive = false;
isFixedUturn = false;
loadData;
Gobs = Obs;
Gnbobs = size(Obs,1);
nbobs = 1;
Op = Op_structure;
initialize_optimization_structure();
Gradient = zeros(nbobs,Op.n);
Op.x = [-1,0,0,0]';
%Op.x = [-1,0,0,0]';
vt = zeros(Gnbobs,1);
for i = 1:Gnbobs
    Obs = Gobs(i,:);
    [Op.value, Op.grad ] = getLL2_nested();
    fprintf(' Path: ');
    disp(full(Obs(2:end)));
    fprintf('     Probability: %f \n', exp(-Op.value));
    vt(i) = exp(-Op.value);
end

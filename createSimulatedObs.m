%   Generate obs
%%
global Op;
global Obs;

disp('Observation generating ....')
Op.n = 4;
x0 = [-1.8 -1.0 -1.5 -20.0 ]';
ODpairs = Obs(:,1:2);
nbobsOD = 1;
filename = './simulatedData/Obs_for_PSL.txt';
generateObs(filename, x0, ODpairs, nbobsOD);

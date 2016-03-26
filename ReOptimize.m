....
%%  Return optimizing

Op.x = Op.x - Op.step;
Op.delta = 0.005;
while (true)    
  Op.k = Op.k + 1;
  if strcmp(Op.Optim_Method,OptimizeConstant.LINE_SEARCH_METHOD);
    ok = line_search_iterate();
    if ok == true
        PrintOut(Op);
    else
        disp(' Unsuccessful process ...')
        break;
    end
  else
    ok = btr_interate();
    PrintOut(Op);
  end
  [isStop, Stoppingtype, isSuccess] = CheckStopping(Op);  
  %----------------------------------------
  % Check stopping criteria
  if(isStop == true)
      isSuccess
      fprintf('The algorithm stops, due to %s \n', Stoppingtype);
      resultsTXT = [resultsTXT sprintf('The algorithm stops, due to %s \n', Stoppingtype)];
      break;
  end
end
%% Compute variance - Covariance matrix
%PrintOut(Op);
%disp(' Calculating VAR-COV ...');
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


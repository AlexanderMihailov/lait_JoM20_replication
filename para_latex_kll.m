% Create matrix of paramater names - text/LateX strings
% Last modified: AM150805 -> AM160518 _kll appended and last three
% parameters below commented out

paraname = '\beta';  % AM150805, added missing in THIS LINE ONLY str2mat() - AND DELETED: caused an error message 
paraname = str2mat(paraname,'\alpha');
paraname = str2mat(paraname,'h');
paraname = str2mat(paraname,'\sigma');
paraname = str2mat(paraname,'\phi');
paraname = str2mat(paraname,'\eta');
paraname = str2mat(paraname,'\delta_{H}');
paraname = str2mat(paraname,'\delta_{F}');
paraname = str2mat(paraname,'\theta_{H}');
paraname = str2mat(paraname,'\theta_{F}');
paraname = str2mat(paraname,'a_{1}');
paraname = str2mat(paraname,'a_{2}');
paraname = str2mat(paraname,'a_{3}');
paraname = str2mat(paraname,'b_{1}');
paraname = str2mat(paraname,'b_{2}');
paraname = str2mat(paraname,'b_{3}');
paraname = str2mat(paraname,'c_{1}');
paraname = str2mat(paraname,'c_{2}');
paraname = str2mat(paraname,'c_{3}');
paraname = str2mat(paraname,'\rho_{H}');
paraname = str2mat(paraname,'\rho_{F}');
paraname = str2mat(paraname,'\rho_{a}');
paraname = str2mat(paraname,'\rho_{q}');
paraname = str2mat(paraname,'\rho_{s}');
paraname = str2mat(paraname,'\rho_{r}');
paraname = str2mat(paraname,'\mu_{q}');
paraname = str2mat(paraname,'\mu_{y}');
paraname = str2mat(paraname,'\mu_{r}');
paraname = str2mat(paraname,'\sigma_{H}');
paraname = str2mat(paraname,'\sigma_{F}');
paraname = str2mat(paraname,'\sigma_{a}');
paraname = str2mat(paraname,'\sigma_{q}');
paraname = str2mat(paraname,'\sigma_{s}');
paraname = str2mat(paraname,'\sigma_{\pi*}');
paraname = str2mat(paraname,'\sigma_{y*}');
paraname = str2mat(paraname,'\sigma_{r*}');
paraname = str2mat(paraname,'\sigma_{r}');
%paraname = str2mat(paraname,'\chi');    %Added by Antonio 23-09-2014 -> AM160518: commented out
%paraname = str2mat(paraname,'\rho_{g}');    %Added by Antonio 23-09-2014 -> AM160518: commented out
%paraname = str2mat(paraname,'\sigma_{g}');    %Added by Antonio 23-09-2014 -> AM160518: commented out

% paraname = str2mat(paraname,'const_{\pi_{F}}');
% paraname = str2mat(paraname,'const_{\pi}');
% paraname = str2mat(paraname,'const_{r}  ');
% paraname = str2mat(paraname,'const_{\pi^{\ast}}');
% paraname = str2mat(paraname,'const_{r^{\ast}}');
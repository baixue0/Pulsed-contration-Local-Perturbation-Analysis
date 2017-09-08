classdef DefaultValues
    
    properties(Constant)
        STARTDATA = struct(  ...
            'iterationN' , 1, ...
            'ADnumber' , 7, ...
            'amplitude' , 1e-03, ...
            'orbitpoints' , 100, ...
            'increment' , 1e-05, ...
            'showmultipliers' , true);
        
        CONTINUERDATA = struct(    ...
            'InitStepSize' , 0.01,    ...
            'MinStepSize' , 1e-5,    ...
            'MaxStepSize' , 0.1,    ...
            'MaxNewtonIters' , 3,    ...
            'MaxCorrIters' , 10,    ...
            'TestIters' , 10,    ...
            'VarTolerance' , 1e-06,    ...
            'FunTolerance' , 1e-06,    ...
            'TestTolerance' , 1e-05,    ...
            'Adapt' , 3,    ...
            'MaxNumPoints' , 300,    ...
            'ClosedCurve' , 50 );
        
        OPTIONS = struct( ...
            ... %Suspend Computation
            'suspend' , 1 , ...%  1 = at special points, 2 = each point , 3 = never.
            ... % Maximum number of untitled curves of a particular type'
            'archive' , 4 , ...
            ... % Plot after X points
            'plotoutput' , 1);
        
        CURVECOLORCFG = struct( ...
            'LP' , {{'Color' , 'green'}} , ...
            'PD' ,  {{'Color' , 'blue'}}, ...
            'FP' ,  {{'Color' , 'black'}}, ...
            'NS' ,  {{'Color' , 'magenta'}}, ...
            'HO' ,  {{'Color' , 'black'}}, ...
            'HE' ,  {{'Color' , 'black'}}, ...
            'HomT', {{'Color' ,  'red'}}, ...
            'HetT', {{'Color' , 'red'}} , ...
            'O', {{'Color' , [0.8,0.8,0.8] , 'linestyle' , '--'}} );  

        OTHERCOLORCFG = struct( ...
            'LABEL' , {{'Color' ,  'black' , 'FontSize' , 10}} , ...
            'SING' , {{'Marker','*','Color','red'}} , ...
            'USERSING' , {{'Marker','*','Color','magenta'}} , ...
            'ORBITPOINT' , {{'Marker','.','Color','black'}} ...
        );
      
        
        POINTNAMES = struct( ...
            'PD' , 'Period Doubling', ...
            'BP' , 'Branch Point',  ...
            'FP' , 'Fixed Point', ...
            'LP' , 'Limit Point', ...
            'NS' , 'Neimark-Sacker', ...
            'CP' , 'Cusp', ...
            'GPD' , 'Generalized Period Doubling', ...
            'CH' , 'Chenciner', ...
            'R1' , 'Resonance 1:1', ...
            'R2' , 'Resonance 1:2', ...
            'R3' , 'Resonance 1:3', ...
            'R4' , 'Resonance 1:4', ...
            'LPPD' , 'Fold-Flip', ...
            'LPNS' , 'Fold-Neimark-Sacker', ...
            'PDNS' , 'Flip-Neimark-Sacker', ...
            'NSNS' , 'Double Neimark-Sacker', ...
            'P'    , 'Point' , ...
            'CO'    , 'Connecting Orbit' , ...
            'LP_HO'    , 'Limit Point (connecting orbit)' , ...
            'LP_HE'    , 'Limit Point (connecting orbit)' , ...
            'O'    , 'Orbit');
        
        
    end
    
    
    
    
    
    
    
end
